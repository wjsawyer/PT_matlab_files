function [I] = sym_grad_map(A,B)
%apply symmetrical gradient map centered on the median value
%where I is a 16 bit 2d image and B is a binary image, covering the general
%area of interst
if size(A) ~=size(B)
    error('size mismatch');
end

A = im2uint16(A);

%alternate way of adjusting range

A = stretch(A);


%alternate way of adjusting range, imadjust doesn't want to work, not sure why...
% [vec, dim] = vol2vec(A);
% vec = imadjust(vec); %%IMADJUST DOESN'T WORK ON VECTOR? ALL ONES...
% A = vec2vol(vec, dim);


B = im2uint16(B);
B (B ~= 0) = 1;
area = sum(sum(sum(B)));

C = A .* B;
average = sum(sum(sum(C)))/area;


%apply the gradient map
%sets the average as the new max and goes to zero on either side with same
%rate defined by distance to closest end. 
I = A; % normal
if average < intmax(class(A))/2
    
     I(I > 1.5*average) = .5*average;
     I(I < .5*average) = .5*average;
     I(I > average) = 2*average - I(I > average); %set higher values symmetrical

else %avg is greater
    I = intmax(class(I)) - I; %flip around so same as above
   I(I > 1.5*average) = .5*average;
     I(I < .5*average) = .5*average;
     I(I > average) = 2*average - I(I > average); %set higher values symmetrical

end    

I = stretch(I);

% figure;
% imagesc(A(:,:,149));
% figure;
% imagesc(B(:,:,149));
% figure;
% imagesc(I(:,:,149));
% 
% use graythresh function
% 

end


function [OUT] = stretch(IN)
mx = max(max(max(IN)));
mn = min(min(min(IN)));
range = mx-mn;
OUT = IN-mn;
scalar = double(intmax(class(IN))) / double(range);
OUT = OUT*scalar;
end

% 
% %%/
%         aux = load('B.mat');
%         %aux = load(strcat(directory, '\B.mat'));
%         names=fieldnames(aux);
%         Avizo_B_mat = aux.(names{1});
%         [~,ydim, xdim, zdim] = size(Avizo_B_mat);
%         B = reshape(Avizo_B_mat(1,:,:,:),ydim,xdim,zdim);
%         disp('using avizo label field');
%         
%            aux = load('A.mat');
%     %aux = load(strcat(directory, '\A.mat'));
%     names=fieldnames(aux);
%     Avizo_A_mat = aux.(names{1});
%     A = reshape(Avizo_A_mat(1,:,:,:),ydim,xdim,zdim);

%%