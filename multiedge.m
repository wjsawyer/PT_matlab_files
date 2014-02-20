% user defined paramaters
    % Voxel size in microns
    CELLSIZE=36.04;

    % thresholds 
    othresh = 5; % threshold for finding edges outward from the surface
    ithresh = 2; % threshold for finding edges inward from the surface
    mt = 5; % outward distance at which mineral phases are identified
    binthresh = 0; %gradient value above which is rounded to 1 when converting gradedge to binedge
    
    
%     % Mineral ID grayscale threshold values
%     bio = 105; % 8-bit minimum value for biotite
%     qf = 90; % 8-bit maximum value for undiferentiated quartz/feldspar
 %   bio = 30000; % 16-bit minimum value for biotite
  %  qf = 26000; % 16-bit maximum value for undiferentiated quartz/feldspar

    
% read in data volumes and reformat to 3D from 4D

    % load data
    olddir = pwd;
    directory = uigetdir; %use this line for gui to pick directory
    %directory = 'F:\pseudo_harddrive\Avizo\WAG12-1A_Avizo\will\matlab\A_median_filtered\65-256-interp';
    
   
    cd(directory);
    
     % binary label map volume
% %     if exist('autolabel','file')==2
% %         aux = load('autolabel.mat');
% %         names=fieldnames(aux);
% %         B = aux.(names{1});
% %         [ydim, xdim, zdim] = size(B);
% %         disp('using label field output from autolabel');
% %     else
        aux = load('B.mat');
        %aux = load(strcat(directory, '\B.mat'));
        names=fieldnames(aux);
        Avizo_B_mat = aux.(names{1});
        [~,ydim, xdim, zdim] = size(Avizo_B_mat);
        B = reshape(Avizo_B_mat(1,:,:,:),ydim,xdim,zdim);
        disp('using avizo label field');
% %     end
    
    
    BB = zeros((ydim+(2*othresh)), (xdim+(2*othresh)), (zdim+(2*othresh)));
    for x = 1:ydim %(othresh + 1):(ydim + othresh)
       for y = 1:xdim %(othresh + 1):(xdim + othresh)
           for z = 1:zdim %(othresh + 1):(zdim + othresh)
               BB(x+othresh,y+othresh,z+othresh) = B(x,y,z);
           end
       end
    end
    for x = 1:othresh
    BB(x,:,:) = BB(othresh+1,:,:);
    BB(((ydim + (2 * othresh)) + 1 - x),:,:) = BB((ydim + othresh),:,:);
    end
    for y = 1:othresh
    BB(:,y,:) = BB(:,othresh+1,:);
    BB(:,((xdim + (2 * othresh)) + 1 - y),:) = BB(:,(xdim + othresh),:);
    end
    for z = 1:othresh
    BB(:,:,z) = BB(:,:,othresh+1);
    BB(:,:,((zdim + (2 * othresh)) + 1 - z)) = BB(:,:,(zdim + othresh));
    end

    % grayscale image volume
    aux = load('A.mat');
    %aux = load(strcat(directory, '\A.mat'));
    names=fieldnames(aux);
    Avizo_A_mat = aux.(names{1});
    A = reshape(Avizo_A_mat(1,:,:,:),ydim,xdim,zdim);
    
        AA = zeros((ydim+(2*othresh)), (xdim+(2*othresh)), (zdim+(2*othresh)));
    for x = 1:ydim %(othresh + 1):(ydim + othresh)
       for y = 1:xdim %(othresh + 1):(xdim + othresh)
           for z = 1:zdim %(othresh + 1):(zdim + othresh)
               AA(x+othresh,y+othresh,z+othresh) = A(x,y,z);
           end
       end
    end
    for x = 1:othresh
    AA(x,:,:) = AA(othresh+1,:,:);
    AA(((ydim + (2 * othresh)) + 1 - x),:,:) = AA((ydim + othresh),:,:);
    end
    for y = 1:othresh
    AA(:,y,:) = AA(:,othresh+1,:);
    AA(:,((xdim + (2 * othresh)) + 1 - y),:) = AA(:,(xdim + othresh),:);
    end
    for z = 1:othresh
    AA(:,:,z) = AA(:,:,othresh+1);
    AA(:,:,((zdim + (2 * othresh)) + 1 - z)) = AA(:,:,(zdim + othresh));
    end
   
    
    clear aux names
    cd(olddir);
    

    % find size of 3D volume
    
%     y = 1:ydim; % define vector of rows
%     x = 1:xdim; % define vector of columns

    % convert 4D variables to 3D
    
    

    %3d array to store the overlapping edges with gradient values
    %gradient value = distance from hand trace (value of n/thresh means
    %voxel is thresh-n away from hand traced line (+ up)).
    
    gradedges = zeros(size(AA)); %stores initial gradient
    ygrad = gradedges;
    xgrad = gradedges; 
    binedges = gradedges;

    canny = zeros(size(AA));

for z=othresh + 1:zdim + othresh
    B1 = BB(:,:,z);
    A1 = AA(:,:,z);
    canny(:,:,z) = edge(A1, 'canny');
    
    %below for loop creates gradient adge of the label field
    %changed all gradients to be positive and added to the current value so
    %that gradients wouldn't overwrite or negate eachother
    for i=(othresh+1):(xdim+othresh)
        for j=(othresh+1):(ydim+othresh) 
          if B1(j,i)==1 %if we are on colored surface   
              %%%apply gradient in the x direction
              if B1(j-1,i) ==0 % if on 'top side' of a section
                  for m=0:othresh %apply outside gradient above
                    xgrad(j-m,i,z) = xgrad(j-m,i,z) + (1-(m/othresh));
                  end
                   for m=1:ithresh %apply inside gradient below (-)
                    xgrad(j+m,i,z) = xgrad(j+m,i,z) + (1-(m/ithresh));
                   end
              end
              if B1(j+1,i) == 0 %if on bottom side of a section
                  for m=0:othresh %apply outside gradient below
                    xgrad(j+m,i,z) = xgrad(j+m,i,z) + (1-(m/othresh));
                  end
                for m=1:ithresh %apply inside gradient above (-)
                    xgrad(j-m,i,z) = xgrad(j-m,i,z) + (1-(m/ithresh));
                end
              end
              
              
              %%%apply gradient in the y direction
              if B1(j,i-1) ==0 % if on 'left side' of a section
                  for m=0:othresh %apply outside gradient to left
                    ygrad(j,i-m,z) = ygrad(j,i-m,z) + (1-(m/othresh));
                  end
                   for m=1:ithresh %apply inside gradient to right (-)
                    ygrad(j,i+m,z) = ygrad(j,i+m,z) + (1-(m/ithresh));
                   end
              end
              if B1(j,i+1) == 0 %if on right side of a section
                  for m=0:othresh %apply outside gradient below to right
                    ygrad(j,i+m,z) = ygrad(j,i+m,z) + (1-(m/othresh));
                  end
                for m=1:ithresh %apply inside gradient to left (-)
                    ygrad(j,i-m,z) = ygrad(j,i-m,z) + (1-(m/ithresh));
                end
              end
              
              
          end
        end
    end
end
gradedges = (xgrad + ygrad) .*canny;


%this makes issues with edgecon, should be implemented later in the process
% %set borders equal to 1 if they were in B
% binedges(1+othresh,:,:) = binedges(1+othresh,:,:) + BB(1+othresh,:,:);
% binedges(ydim+othresh,:,:) = binedges(ydim+othresh,:,:) + BB(ydim+othresh,:,:);
% binedges(:,1+othresh,:) = binedges(:,1+othresh,:) + BB(:,1+othresh,:);
% binedges(:,xdim+othresh,:) = binedges(:,ydim+othresh,:) + BB(:,ydim+othresh,:);

%convert to binary
binedges(gradedges>binthresh) = 1;

%get rid of excess areas around perimeter
binedges = binedges(othresh+1:ydim+othresh, othresh+1:xdim+othresh, othresh+1:zdim+othresh);







    
    %k-nearest neighbor searching
    
 