

function phi = PWF(Simg, Rimg, L1, L2, imgPix, photonE,varargin)

%% parsing
p = inputParser();

%%% FOVcrop paramters
addParameter(p, 'bgWingRatio', 0.05, @isscalar);

%%% preprocessing parameters
addParameter(p, 'flatten_background', true, @islogical);

%%% iteration parameters
addParameter(p, 'eta', [1,1], @isnumeric);
addParameter(p, 'iterMax', 5000, @isnumeric);
addParameter(p, 'corrCrit', 1e-5, @isnumeric);

%%% regularization parameter
addParameter(p, 'gamma', 1, @isnumeric);
addParameter(p, 'regImageCoeff', [0.1,0.01], @isnumeric);

%%% show parameter
addParameter(p, 'figVisPeriod',100, @isscalar);
addParameter(p, 'coor_comp_frequency',10, @isscalar);

%%%
parse(p, varargin{:});
bgWingRatio  = p.Results.bgWingRatio;

flatten_background    = p.Results.flatten_background;
eta    = p.Results.eta;
iterMax = p.Results.iterMax;
corrCrit = p.Results.corrCrit;

gamma=p.Results.gamma;
regImageCoeff=p.Results.regImageCoeff;

figVisPeriod=p.Results.figVisPeriod;
coor_comp_frequency=p.Results.coor_comp_frequency;


%% setup
if isempty(figVisPeriod)
    figVisPeriod = iterMax;
end

%%% gpuInd
if gpuDeviceCount("available")
    gpuInd = 1;    
end

XrayWl  = Etowl(photonE);

%% FOV
Yrange = [1, size(Rimg,1)];
bgWingPix0 = ceil( bgWingRatio * size(Rimg,2) );
Xrange = [1+bgWingPix0, size(Rimg,2)-bgWingPix0];

%%% add bgWings to Xrange
bgWingPix0 = ceil( bgWingRatio * (Xrange(2)-Xrange(1)+1) ); % left and right

%%% prevent the bgWing exceeds original FOV
bgWingPix(1) = Xrange(1) - max(1, Xrange(1)-bgWingPix0);
bgWingPix(2) = min(size(Rimg,2), Xrange(2)+bgWingPix0) - Xrange(2);
if any(bgWingPix < bgWingPix0)
    warning( 'Backround pixels may not sufficient: [%d, %d]', bgWingPix(1), bgWingPix(2) )
end

%%% update Xrange
Xrange_wBG = [Xrange(1)-bgWingPix(1), Xrange(2)+bgWingPix(2)];
cropSize = [Yrange(2)-Yrange(1)+1, Xrange_wBG(2)-Xrange_wBG(1)+1];

%%% def. indices
Yind = Yrange(1):Yrange(2);
Xind = Xrange_wBG(1):Xrange_wBG(2);

%%% crop reference speckle
Rimg = Rimg(Yind,Xind,:);


%% spkRes and samRes def.
figure(221141),
cohLength = sqrt(pi)*speckleGrainSizeFitXY( Rimg ) * imgPix; % nm, 
fieldBWsigma = 1/sqrt(pi)./cohLength;     % nm^-1 
intensityBWsigma = sqrt(2/pi)./cohLength; % nm^-1 

%%% def. diffAngle and xpadWingSize 
diffAngle = XrayWl.*fieldBWsigma/2; % rad
minAddWingRatio = ((L1+L2)*diffAngle/imgPix) ./ cropSize * 2; % last factor is for safety
xpadWingRatio = max([minAddWingRatio, bgWingRatio])*3; %  % last factor is for safety


%% flaten_background
%%% to remove vertical stripe that happened in the experiements due to the vertical instability of the focusing optics
bgXind1 = 1:bgWingPix(1);
bgXind2 = (cropSize(2)-bgWingPix(2)+1):cropSize(2);
refMeanBG1 = Rimg(:,bgXind1,:);
refMeanBG2 = Rimg(:,bgXind2,:);
xvec = cast(1:cropSize(2),'like',Rimg);
xvec = (xvec - mean(bgXind1))/(mean(bgXind2)-mean(bgXind1));


%% padding
%%% x-space: cropSize --> xpadSize
xpadSize = cropSize + 2*round(xpadWingRatio.*cropSize);
refISpk  = borderPadding(Rimg, xpadWingRatio,'uniform');

%%% imgPix update
imgPix2  = min(1./(3*intensityBWsigma))/2; % +-3*sqrt(2) sigma range
pixRatio = imgPix2/imgPix; % considering intensity bandwidth doubling
padSize  = ceil( xpadSize/pixRatio );
pixRatio = mean(xpadSize./padSize);
imgPix2  = imgPix * pixRatio;

%%% x-space: xpadSize --> padSize
refISpk = fftshift(ifft2(ifftshift( msize(fftshift(fft2(ifftshift(refISpk))), padSize)),'symmetric'));
refISpk = refISpk * prod(padSize./xpadSize);

%%% def. camMask
camMask = borderPadding(true(cropSize), xpadWingRatio,'tukey');
camMask = fftshift(ifft2(ifftshift( msize(fftshift(fft2(ifftshift(camMask))), padSize)),'symmetric'));
camMask = camMask * prod(padSize./xpadSize);
if gpuInd
    camMask =gpuArray(single(camMask));
end
% camMask = camMask >0.99;
camMasksq = camMask.^2; %MODIFIED

%%% value normalization
paddedCamSize = round(cropSize/pixRatio);
meanRefSpeckle = mean(mcrop(refISpk, paddedCamSize),[1 2]);
refISpk = refISpk./meanRefSpeckle;

%% Q1, Q2, and incohFilter
%%% constants
type_precise=@(x) single(x);
numerical_error_tweaking=0.00001;
kgrid = ndgrid_matSizeIn(padSize, 1,'centerZero_ifftshift');
kgrid = cat(3, kgrid{1}, kgrid{2})/imgPix2; % nm^-1
if gpuInd
    kgrid = gpuArray(kgrid);
end
krhosq = sum(kgrid.^2,3); % nm^-2

%%% Q1 and Q2
Q1 = single(exp( -1i*pi* krhosq * L1 * XrayWl ));
Q2 = single(exp( -1i*pi* krhosq * L2 * XrayWl ));

%%% incohFilter = exp(-1/2*(pi*lc*kho)^2)
incohFilter = type_precise(exp( -sum( (kgrid .* reshape(1./intensityBWsigma,1,1,[])).^2, 3) ));

%% create the illumination field
%%% Wiener deconvolution with NSR map
FourierNoiseMask = (sqrt(krhosq) > 1/imgPix2/2) & all( abs(kgrid) < (1/imgPix2/2), 3);
% figure, imagesc(FourierNoiseMask)

Sf = abs(fft2(refISpk)).^2;
Sf = reshape(Sf, [], 1);
% figure, imagesc(log10(abs(Sf(:,:,1)))); axis image

bgMean = mean( Sf(FourierNoiseMask,:),  'all');
bgStd  =  std( Sf(FourierNoiseMask,:),0,'all');

% def. NSR
NSR = bgStd./(mean(Sf,2) - bgMean); % Sf->= SNR
NSR(NSR<0) = inf;
NSR = reshape(NSR, padSize);

% Wiener deconv
incohFilter_deconv = incohFilter./ (abs(incohFilter).^2 + NSR); 
illFieldCamPlane = sqrt(ifft2(fft2(refISpk) .* incohFilter_deconv));

%%% backpropagation: camera plane --> diffuser plane
illFieldObjPlane = single(ifft2(fft2(illFieldCamPlane) .* conj(Q2)));


%% Preconditioner and regularizations
%%% Preconditioner
imagPart = krhosq.*imgPix2^2;
precond_kernel = single(1./imagPart);
precond_kernel(1,1) = 0;

%%% regularizations
grad_reg_kernel = 1-exp(-gamma*sum( (kgrid./reshape(fieldBWsigma,1,1,[])).^2, 3));
grad_reg_kernel = single(grad_reg_kernel);

%% Sample speckle preprossessing
%%% gpu
if gpuInd
    samISpk = gpuArray(Simg);
end
samISpk = samISpk(Yind,Xind,:,:);

%%% flaten background
if flatten_background
    samISpk = flaten_background_func(samISpk, xvec, bgXind1, bgXind2, refMeanBG1, refMeanBG2);
end

%%% preprocessing
samISpk = borderPadding(samISpk, xpadWingRatio,'uniform');
samISpk = fftshift(ifft2(ifftshift( msize(fftshift(fft2(ifftshift(samISpk))), padSize)),'symmetric'));
samISpk = samISpk * prod(padSize./xpadSize);
samISpk = samISpk./meanRefSpeckle;
samISpk_deconv = single(ifft2(fft2(double(samISpk))...
    .*double(incohFilter)./ (abs(double(incohFilter)).^2 + numerical_error_tweaking ),'symmetric'));



%% iteration
%%% NAGcoeff
NAGcoeff = single(zeros(iterMax,1));
t0 = 1;
for ii = 1:iterMax
    t1 = (1+sqrt(1+4*t0^2))/2;
    NAGcoeff(ii) = (t0-1)/t1;
    t0=t1;
end

%%% intialization
showPhiShiftParams = struct([]);
etaConst_prob = max( sum(abs(illFieldObjPlane).^2, 3), [], 'all' )^2;
if gpuInd
    probOut = complex(zeros(padSize,'single','gpuArray'));
    prob_GD = complex(zeros(padSize,'single','gpuArray'));
    errorList = zeros(iterMax,1,'single','gpuArray');
    probCorr  = zeros(iterMax,1,'single','gpuArray');
else
    probOut = complex(zeros(padSize,'single'));
    prob_GD = complex(zeros(padSize,'single'));
    errorList = zeros(iterMax,1,'single');
    probCorr  = zeros(iterMax,1,'single');
end

%%% iteration
breakFlag = false;
for ii = 1:iterMax
    probImg = exp(probOut);
    probImg2 = ifft2( fft2(probImg) .* Q1);
    Yiter = ifft2( fft2(illFieldObjPlane.* probImg2) .* Q2);

    %%% ERROR TERM
    target = single(ifft2(fft2(type_precise(abs(Yiter).^2) - samISpk_deconv) .* incohFilter) ).* camMasksq;%MODIFIED FAST AND SAME RESULT
    if mod(ii,coor_comp_frequency) == 0
        errorList(ii-coor_comp_frequency+1:ii) = rms(target,'all');
    end

    %%% BACKWARD
    gPsi = ifft2(fft2(target).* incohFilter);%MODIFIED
    gPsi = ifft2(fft2(gPsi.* Yiter) .* conj(Q2)); % [yy, xx, Nscan, 2]
    gPsi = gPsi .* conj(illFieldObjPlane);
    gPsi = ifft2( fft2(gPsi) .* conj(Q1));
    gPsi = gPsi .* conj(probImg);
    g_prob = mean(gPsi,3);

    % PRECONDITIONING
    gRyrov_imag = imag(g_prob);
    gRyrov_imag = ifft2(fft2(gRyrov_imag) .* precond_kernel,'symmetric');
    g_prob = real(g_prob) + 1i*gRyrov_imag;

    % regularization:
    g_prob = g_prob/etaConst_prob; 
    reg_comp = ifft2( fft2(probOut).*grad_reg_kernel);
    g_prob = g_prob + regImageCoeff(1).*real(reg_comp)+1i.*regImageCoeff(2).*imag(reg_comp);% + bgCorrectCoeff*probOut.*bgMask;

    %%% NAG1 step
    prob_GD_prev = prob_GD;
    prob_GD = probOut - (eta(1).*real(g_prob) + 1i*eta(2).*imag(g_prob));
    probBefore = probOut;
    probOut = prob_GD + NAGcoeff(ii).*(prob_GD - prob_GD_prev);

    %%% error calc
    if mod(ii,coor_comp_frequency) == 0
        probCorr(ii-coor_comp_frequency+1:ii) = -log10( imgCorrCalc(exp(probBefore), exp(probOut), camMask > 0.99) ); % correlation calculation with last iteration.

        %%% breakFlag
        if (probCorr(ii) < corrCrit) && ii > 1
            breakFlag = true;
        end
    end


    %%% show
    if (mod(ii,figVisPeriod) == 0) || (breakFlag && ~isinf(figVisPeriod))
        figure(54356245);
        showProbOut = mcrop(probOut, paddedCamSize);

        subplot(221)
        imagesc(-real(showProbOut)); axis image; colorbar; title('attenuation coeff.'); %colormap gray;
        colormap gray;

        subplot(222)
        [imagShowProbOut, showPhiShiftParams] = PhiShift_XRAY(imag(showProbOut), round(bgWingPix/pixRatio), showPhiShiftParams);
        imagesc(imagShowProbOut); axis image;colorbar; title('phase')
        colormap gray;

        subplot(223)
        yyaxis left;
        semilogy(1:coor_comp_frequency:ii, errorList(1:coor_comp_frequency:ii)); ylabel('RMSE'); axis square;
        yyaxis right; semilogy(1:coor_comp_frequency:ii, probCorr(1:coor_comp_frequency:ii)); ylabel('probCorr');
        yline(corrCrit,'b-','corrCrit')
        title(sprintf('RMSE:%.2e, corr:%.2e',errorList(ii),probCorr(ii)))

        subplot(224)
        imagesc(abs(target)); axis image; colorbar; title('RMSE')
        colormap gray;
    
        drawnow;
    end

    %%% break
    if breakFlag
        break;
    end

    if ii == iterMax
        warning("iteration reaches iterMax..")
    end

end  

%%% undo padding and save
phiTemp = fftshift(ifft2(ifftshift( msize(fftshift(fft2(ifftshift(probOut))), xpadSize) )));
phiTemp = phiTemp * prod(xpadSize./padSize); % xpadSize
phiTemp = mcrop(phiTemp, cropSize);
imagPhiTemp = PhiShift_XRAY(imag(phiTemp), bgWingPix, struct([]) );
phi =  gather(real(phiTemp)+1i*imagPhiTemp);

end



%% subfuntions
function spkIin = flaten_background_func(spkIin, xvec, bgXind1, bgXind2, refMeanBG1, refMeanBG2)

profile_L = mean( refMeanBG1 ./ spkIin(:,bgXind1,:), 2);
profile_R = mean( refMeanBG2 ./ spkIin(:,bgXind2,:), 2);
spkIin = spkIin.*( (profile_R - profile_L).*xvec + profile_L );

end

function corrOut = imgCorrCalc(img1, img2, mask)
corrOut = abs((img1(mask)'*img2(mask)/norm(img1(mask))/norm(img2(mask))));


% Fimg1 = fft2(img1/norm(img1(mask)).*mask );
% Fimg2 = fft2(img2/norm(img2(mask)).*mask );
% acMap = abs(ifft2(conj(Fimg1).*Fimg2));
% corrOut = max(acMap,[],'all');


end

function [phaseImg, PhiShiftParams] = PhiShift_XRAY(phaseImg, bgWingPix, PhiShiftParams )

if isempty(PhiShiftParams)
    %Uimg=Timg(1:200,1:240);
    [yy, xx] = size(phaseImg);

    %%% def. bgMask
    bgMask = false(yy,xx);
    bgInds = [1:bgWingPix(1), (xx-bgWingPix(2)+1):xx];
    bgMask(:,bgInds) = true;

    %%% def. meshes
    xgrids = ndgrid_matSizeIn([yy,xx],0,'centerZero');
    YY = xgrids{1};
    XX = xgrids{2};
    A = [YY(bgMask), XX(bgMask)];
    A = [A, ones(size(A,1),1)]; % [N x 3]

    %%% gpuArray
    if isgpuarray(phaseImg)
        A = gpuArray(A);
        YY = gpuArray(YY);
        XX = gpuArray(XX);
        bgMask = gpuArray(bgMask);
    end

    %%% save
    PhiShiftParams = struct;
    PhiShiftParams.A = A;
    PhiShiftParams.YY = YY;
    PhiShiftParams.XX = XX;
    PhiShiftParams.bgMask = bgMask;
else

    %%% load
    A = PhiShiftParams.A ;
    YY = PhiShiftParams.YY;
    XX = PhiShiftParams.XX;
    bgMask = PhiShiftParams.bgMask;

end



b = phaseImg(bgMask); % [N x 1]
vecOut = A\b; % solution of Ax = b;
phaseImg = phaseImg - (YY*vecOut(1) +XX*vecOut(2) +vecOut(3));

% figure(100),imagesc(phaseImg)
end





