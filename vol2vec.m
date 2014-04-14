
function [vec, dim] = vol2vec(vol)
dim = size(vol);
vec = zeros(1, dim(1)*dim(2)*dim(3));
area = dim(1)*dim(2);

for ind3 = 1:dim(3)
    for ind2 = 1:dim(2)
        for ind1 = 1:dim(1)
    vec((ind3 -1)*area + (dim(1)*(ind2 - 1)) + ind1) = vol(ind1, ind2, ind3);
        end
    end
end



