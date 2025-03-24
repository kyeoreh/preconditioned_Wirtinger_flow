function [cohLength, fitResult2D] = speckleGrainSizeFitXY(imgStack, lowPassWingSize)
if nargin < 2
    lowPassWingSize = [1,1];
end

%%% x-direction
yy = size(imgStack,1);
xx = size(imgStack,2);
mc = mcoor([yy,xx]);
ky = ((1:yy) - mc(1)).'/yy;
kx = ((1:xx) - mc(2)).'/xx;

Fimgx = abs(fftshift(fft(ifftshift( imgStack./mean(imgStack,2)),[],2)));
Fimgy = abs(fftshift(fft(ifftshift( imgStack./mean(imgStack,1)),[],1)));
Fimgx = mean(Fimgx,3);
Fimgy = mean(Fimgy,3);

Fimgx(:,mc(2)-lowPassWingSize(2):mc(2)+lowPassWingSize(2)) = 0; % high-pass filter to erase flat field convolution
Fimgy(mc(1)-lowPassWingSize(1):mc(1)+lowPassWingSize(1),:) = 0;
FIx = mean(Fimgx, 1);
FIy = mean(Fimgy, 2);
FIx = FIx - (FIx(1) + FIx(end))/2;
FIy = FIy - (FIy(1) + FIy(end))/2;

cohLength = zeros(1,2);
fitResult2D  = ones(yy,xx);

for pp = 1:2
    switch pp
        case 1
            test = double( gather( FIy ) );
            ktest = ky;
        case 2
            test = double( gather( FIx.') );
            ktest = kx;
    end
    [maxVal, maxInd] = max(test);
    validInd = true(size(test));
    validInd(maxInd:(2*mc(pp) - maxInd)) = false;

    %%% gauss fit
    fitFunc = @(a,c,x) a.*exp(-1/2*(pi*x*c).^2);
    fitOp = fitoptions(fitFunc);
    fitOp.StartPoint = [maxVal,  1/4];
    f1 = fit(ktest(validInd), test(validInd), fitFunc, fitOp);
    
    switch pp
        case 1
            fitOut = f1(ktest);
            plot(ktest, fitOut,'k', ktest(validInd), test(validInd).','r.')
            xline(1/f1.c,'-.r','1/( coherence length )')

        case 2
       
            fitOut = f1(ktest);               
            hold on; 
            plot(ktest, fitOut,'k', ktest(validInd), test(validInd).','b.'); hold off;
            xline(1/f1.c,'-.b','1/( coherence length )')
            hold off;
            fitOut = fitOut.';
    end
    cohLength(pp) = f1.c; 
    fitResult2D = fitResult2D .* fitOut;

end


title(sprintf('coherence length = [%.2f, %.2f] pixels', cohLength))
drawnow;

end