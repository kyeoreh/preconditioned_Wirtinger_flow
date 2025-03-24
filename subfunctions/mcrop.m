% cimg = mcrop(img,[yy,xx,zz])

function img = mcrop(img, outputSize)

inputSize = size(img);
inputDim = length(inputSize);
outputDim = length(outputSize);

if inputDim < outputDim
    error('Dimensions of input images must be larger than given dimesion')
       
elseif inputDim > outputDim    
    conservingDim = (outputDim+1):inputDim;
    outputSize(conservingDim) = inputSize(conservingDim);
end
if  any(outputSize > inputSize)
    error('Final image size must be smaller than initial image size')
end
if all(outputSize == inputSize)
    return;
end

dimInds = cell(inputDim, 1);
for kk = 1:1:inputDim
    mp = floor(inputSize(kk)/2)+1;
    mside = round((outputSize(kk)-1)/2);
    dimInds{kk} = (mp-mside) : (mp+outputSize(kk)-mside-1);
end

img = img(dimInds{:});
