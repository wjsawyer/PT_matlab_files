%%script to take MULTIEDGE result and clean up: 
    %this includes connecting edges and filling in 
    %to be run after multiedge while variables still exist
    %eventually append to end of multiedge
   
    
    %load in binedges 
if (exist('binedges','var') == 1) && (exist('A','var') == 1)
binedges = evalin('base', 'binedges'); %main result from multiedge
A = evalin('base', 'A');               %greyscale image volume
else
    error('no A or binedges available. run multiedge.m');
end   

[xdim, ydim, zdim] = size(binedges);
connect = binedges;
counter = zeros(zdim, 2);
%a vector to store the number of objects in each layer when done

%check to see if running with gradient already done
if exist('equi','var') == 1
    equi = evalin('base','equi');
else
    error('no gradient available. run imggrad.m');
end    
  

for z = 7%1:zdim
z;
%collect each point in a set of endpoints to its nearest neighbor  not
%including endpoints from the same segment

layer = bwmorph(binedges(:,:,z),'skel', Inf); %skeletonize image
%layer = bwmorph(layer, 'bridge');

isum = 1;
fsum = 2;



while isum ~=fsum && counter(z,1)<10   %set maximum number of loops per z layer
    isum = sum(sum(layer));
    counter(z,1) = counter(z,1) + 1;
 layer = bwareaopen(layer, 2);   

%define matrix with endpoints
endpoints = im2double(bwmorph(layer, 'endpoints'));
%remove points based on CT image brightness
% bavg = 29000; %edge brightness values should fall within bavg +/- brange
% brange = 2000;
% endpoints(abs(A(:,:,z)- bavg)< brange) =0;

%define matrix of segments, each with a different value
%define pieces just once per layer, or every step?
cc = bwconncomp(layer);
pieces = im2double(labelmatrix(cc));

% set the value of each endpoint to that of the segment it is on
% endpoints = pieces.*endpoints;


%%%%%%%convert matrix of endpoints to a vecotr of their x,y coordinates
epvector = zeros(3,sum(sum(endpoints)));
%columns 1 and 2 are coordinates, column 3 is the CT gradient direction
i = 1;
for x = 1 : xdim 
    for y = 1 : ydim
        if endpoints(x,y) == 1
            epvector(1:2,i) = [x;y];
            epvector(3,i) = equi(x,y,z);
            i = i + 1;
        end
    end
end

epclose = zeros(4,size(epvector, 2));
%columns 1 and 2 are coordinates, column 3 is the ct gradient direction at
%the point, column 4 is the distance between points
epclose(4,:) = sqrt(xdim^2 + ydim^2); %% or an arb. large value

for m = 1:size(epvector, 2) %select point to find nearest neighbor to

    for j = 1:size(epvector, 2) %run through all points
           if pieces(epvector(1,m),epvector(2,m)) ~= pieces(epvector(1,j),epvector(2,j)) %%%|| (epvector(2,m)==xdim && epclose(2,m) == xdim) %if second point is on a different segment...or if they are both on the righthand edge
               dist = sqrt((epvector(1,m)-epvector(1,j))^2+(epvector(2,m)-epvector(2,j))^2); %calculate distance
               if dist < epclose(4,m) %if dist is smaller than previous dist values....
                   epclose(1:3,m) = epvector(1:3,j); %...then record the point and grad direction
                   epclose(4,m) = dist; %and then the distance between points
                   
               end
           end
    end
end
%epclose(:,j) stores the coordinates of the closest endpoint to j and the
%distance between them





%now loop through each point and draw a line to its nearest neighbor.
%manypoint's will be their nearestneighbor's nearestneighbor, so like will
%be drawn twice
%also, have a upper limit for distance 
%think of drawing a line from point to closest point
dthresh = 29; %pixel distance cut off
rng = 30; %+/- range in degrees for comparing dir with equi
bin2 = zeros([xdim, ydim]);
bin3 = zeros([xdim, ydim]);

for p = 1:size(epvector, 2)
   if epclose(4,p) <=dthresh
       hdist = epclose(2,p) - epvector(2,p); %horizontal distance between points
       hdist = hdist - sign(hdist); %reduce length by 1, regardless of sign
       vdist = epclose(1,p) - epvector(1,p); %vertical distance between points
       vdist = vdist - sign(vdist); %reduce length by 1, regardless of sign
       dir = atand(vdist/hdist); %in degrees
     if (((dir-rng)<epvector(3,p) && epvector(3,p)<(dir+rng)) ||....
        ((dir-rng)<(epvector(3,p)+180) && (epvector(3,p)+180)<(dir+rng)) ||....
        ((dir-rng)<(epvector(3,p)-180) && (epvector(3,p)-180)<(dir+rng))) &&...
        (((dir-rng)<epclose(4,p) && epclose(4,p)<(dir+rng)) ||....
        ((dir-rng)<(epclose(4,p)+180) && (epclose(4,p)+180)<(dir+rng)) ||....
        ((dir-rng)<(epclose(4,p)-180) && (epclose(4,p)-180)<(dir+rng)))
    
      if hdist >= vdist 
       for r = 1: abs(hdist)% for each pixel column between the two (horizontal distance)
           %one pixel changed per column (for horizontalish lines)
           %calc height value for that position, round to a whole number and set point = 1
           bin2(epvector(1,p) + round(vdist*(r/abs(hdist))), epvector(2,p)+sign(hdist)*r) = 1; 
       end
      else 
       for r = 1: abs (vdist)    
           %one pixel changed per row (for verticalish lines)
           %calc height value for that position, round to a whole number and set point = 1
           bin3(epvector(1,p)+sign(vdist)*r, epvector(2,p) + round(hdist*(r/abs(vdist)))) = 1; 
       end
      end
      end 
   end
end
%combine layers and skeletonize
layer = layer + bin2 + bin3;
layer = bwmorph(layer,'skel', Inf);


fsum = sum(sum(layer));
%layer = bwmorph(layer,'spur'); %spur removes pixel regardless, so if it happens before the fsum check fsum will decrease and never reach a s.s.

end

% layer = bwmorph(layer, 'diag');
% layer = imfill(layer, 'holes');
connect(:,:,z) = layer;

cc = bwconncomp(layer);
counter(z,2) = cc.NumObjects;
end

figure;
imagesc(connect(:,:,z));
title('layer from final connect data');


% %save connect matrix to .mat file for avizo processing
% cd(directory);
% save('./matlab_out/edgecon.mat', 'connect');
% cd(olddir);

%%figures
% figure; 
% imagesc(bin2 + 2*endpoints);
% title('edges added  one pixel per column');
% figure;
% imagesc(bin3  + 2*endpoints);
% title('edges added  one pixel per row');

