addpath(genpath('subfunctions'));

%% setup parameters
%%% source
photonE = 10; % keV

%%% lengths (nm)
L1 = 3e6;  % sample-diffuser
L2 = 20e6; % diffuser-detector

%%% microscope LMPLFLN10x
mag = 10;

%%% camera: pco.edge 5.5
camPix = 6.5e3;
imgPix = camPix/mag;

%% data load
dfFiles = dir("data\DF_*.tif");
ffFiles = dir("data\FF_*.tif");
refFiles = dir("data\dif_FF_*.tif");
samFiles = dir("data\TOMO_*.tif");

%%% load: dark field
dfimg = single(0);
for ii = 1:length(dfFiles)
    dfimg = dfimg + single(imread(fullfile(dfFiles(ii).folder, dfFiles(ii).name)));    
end

%%% load: flat field
ffimg = single(0);
for ii = 1:length(dfFiles)
    ffimg = ffimg + single(imread(fullfile(ffFiles(ii).folder, ffFiles(ii).name)));    
end

%%% load: reference speckle
Rimg = single(0);
for ii = 1:length(dfFiles)
    Rimg = Rimg + single(imread(fullfile(refFiles(ii).folder, refFiles(ii).name)));    
end

%%% load: reference speckle
Simg = single(imread(fullfile(samFiles(1).folder, samFiles(1).name)));


%%% flat field correction
Simg = (double(Simg)-dfimg)./(ffimg-dfimg);
Rimg = (double(Rimg)-dfimg)./(ffimg-dfimg);

figure,
subplot(121), imagesc(Simg); axis image; title('sample speckle'); colorbar;
subplot(122), imagesc(Rimg); axis image; title('reference speckle'); colorbar;


%% PWF
phi = PWF(Simg, Rimg, L1, L2, imgPix, photonE);
figure,
subplot(121), imagesc(-real(phi)); axis image; title('attenuation coeff.'); colorbar;
subplot(122), imagesc(imag(phi)); axis image; title('phase'); colorbar;
colormap gray;
