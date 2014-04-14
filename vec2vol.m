function [vol] = vec2vol(vec, dim)
%ensure dim is a vector with 3 elements
%ensure length(vec) = dim(1)*dim(2)*dim(3)
area = dim(1)*dim(2);
vol = zeros(dim);

for i = 0:length(vec)-1
         ind3 = floor(i / area) + 1;
     ind2 = floor(rem(i, area) / dim(1)) +1;
    ind1 = rem(rem(i, area), dim(1))+1;
%     ind3 = idivide(i, area) + 1;
%     ind2 = idivide(rem(i, area), dim(1)) +1;
%     ind1 = rem(rem(i, area), dim1);
    
    vol(ind1, ind2, ind3) = vec(i+1);
end