function [img2,windowOut] = borderPadding(img, wingRatio,style)

if isscalar(wingRatio)
    wingRatio = [wingRatio,wingRatio];
end

sz = size(img,[1,2]);
Nd = ndims(img);
assert(Nd >= 2)

mm = round(wingRatio.*sz);
padSz = sz + 2*mm;
% img = img0;

img2 = img;
dimVec = 1:Nd;
dimVec(1:2) = [2,1];
for dd = 1:2
    if dd == 2 % transpose        
        img2 = permute(img2, dimVec);
    end    
    sztemp = size(img2);
    img2 = reshape(img2,sztemp(1),[]);
    if strcmpi(style,'uniform')
        img2 = cat(1, repmat(img2(1,:),[mm(dd),1]), img2(:,:), repmat(img2(end,:),[mm(dd),1]) );
    else
        img2 = cat(1, flipud(img2(1:mm(dd),:)), img2(:,:), flipud(img2(end-mm(dd)+1:end,:)));
    end
    img2 = reshape(img2,[size(img2,1), sztemp(2:end)]);

    if dd == 2 % transpose        
        img2 = ipermute(img2, dimVec);
    end
end

%%% tukey windowing
if strcmp(style,'tukey')
    w1 = tukeywin(padSz(1),2*wingRatio(1));
    w2 = tukeywin(padSz(2),2*wingRatio(2));

    windowOut = w1.*w2';
    img2 = img2.*windowOut;
end

end