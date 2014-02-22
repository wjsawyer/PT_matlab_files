function [epvector, pieces] = ctdir(bin, ct)

eq = imggrad(bin, ct);


bin = bwareaopen(bin, 2); %removes single points useful direction script

%define matrix with endpoints
endpoints = im2double(bwmorph(bin, 'endpoints'));
%remove points based on CT image brightness
% bavg = 29000; %edge brightness values should fall within bavg +/- brange
% brange = 2000;
% endpoints(abs(A(:,:,z)- bavg)< brange) =0;

%define matrix of segments, each with a different value
%define pieces just once per layer, or every step?
cc = bwconncomp(bin);
pieces = labelmatrix(cc);

% set the value of each endpoint to that of the segment it is on
% endpoints = pieces.*endpoints;




%%%%%%%convert matrix of endpoints to a vecotr of their x,y coordinates
epvector = zeros(3,sum(sum(endpoints)));
%columns 1 and 2 are coordinates, column 3 is the CT gradient direction
i = 1;
   [xdim, ydim] = size(bin);
for x = 1 : xdim 
    for y = 1 : ydim
        if endpoints(x,y) == 1
            epvector(1:2,i) = [x;y];
            epvector(3,i) = eq(x,y);
            i = i + 1;
        end
    end
end
