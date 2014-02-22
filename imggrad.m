%image gradient
%mostly copied from autolabel on 2/7/14
%produces an image of gradient magnitude and one of gradient direction
%processes entire volume but calculations are done in 2d, hopefully 3d
%support

function equi = imggrad(B, A)
%     olddir = pwd;
%     if exist('directory','var')~=1
%     directory = uigetdir;
%     end
%     cd(directory);
%     aux = load('A.mat');
%     names=fieldnames(aux);
%     Avizo_A_mat = aux.(names{1});
%     clear aux names
%     % find size of 3D volume
%     [~,xdim, ydim, zdim] = size(Avizo_A_mat);
%     A = reshape(Avizo_A_mat(1,:,:,:),xdim,ydim,zdim);
%     
%     
%     
%     
%     aux = load('B.mat');
%     names=fieldnames(aux);
%     Avizo_B_mat = aux.(names{1});
%     B = reshape(Avizo_B_mat(1,:,:,:),xdim,ydim,zdim);
%     clear aux names
%     cd(olddir);
%     
    A(A>50000) = 50000;
    SE = strel('rectangle', [12 12]);
    B = imdilate(B, SE);

%% plots a filled in circle for testing
%      A = zeros(100,100 ,100);
%      for i = 1:100
%         for j = 1:100
%              for k = 1:100
%                 if sqrt((i-50)^2 + (j-50)^2 + (k-50)^2) <=40
%                     A(i,j,k) = 1;
%                 end
%              end
%         end
%      end
     
%      for f = 1:100
%      filt=fspecial('gaussian',7,1);
% A(:,:,f)=conv2(A(:,:,f),filt,'same');
%      end
     
%% set up space and define constants
     [xdim, ydim] = size(A);
    dx = zeros(size(A));
    dy = dx;
    
    %values for checking change in brightness - defined distance to check
    %in given direction. 
    xx = 2;
    yy = xx;

    %%define a matrix of equal values to A with a 'buffer region' on all
    %%edges of width defined above so that calculations of the change in
    %%brightness can be done for all voxels in initial space. 
    AA = zeros((xdim+(2*xx)), (ydim+(2*yy)));
    for x = 1:xdim %(xx + 1):(xdim + xx)
       for y = 1:ydim %(yy + 1):(ydim + yy)
           AA(x+xx,y+yy) = A(x,y);
       end
    end
    for x = 1:xx
    AA(x,:) = AA(xx+1,:);
    AA(((xdim + (2 * xx)) + 1 - x),:) = AA((xdim + xx),:);
    end
    for y = 1:yy
    AA(:,y) = AA(:,yy+1);
    AA(:,((ydim + (2 * yy)) + 1 - y)) = AA(:,(ydim + yy));
    end
%     for z = 1:zz
%     AA(:,:,z) = AA(:,:,zz+1);
%     AA(:,:,((zdim + (2 * zz)) + 1 - z)) = AA(:,:,(zdim + zz));
%     end
%% 
    %calculation of gradient magnitude in each direction
    filter = fspecial('gaussian');
    for x = (xx + 1):(xdim + xx)
       for y = (yy + 1):(ydim + yy)
          % for z = (zz + 1):(zdim + zz)
               if B(x-xx, y-yy) == 1 %only calculate gradient where the label field is ==1
                     for q = 1: xx
                   dx(x-xx,y-yy) = dx(x-xx,y-yy) + (AA(x-q,y) - AA(x+q,y))/2; %gradient in x direction
                   dy(x-xx,y-yy) = dy(x-xx,y-yy) + (AA(x,y+q) - AA(x,y-q))/2; %gradient in y direction
                      %dz(x-xx,y-yy,z-zz) = (AA(x,y,z+zz) - AA(x,y,z-zz))/2; %gradient in z direction
                     end
               end
           %end
       end
    end
    
   % mag = sqrt(dx.*dx + dy.*dy); %for now dz = 0
    
    %IS THIS DX/DY OR DY/DX????
    dr = atand(dx./dy); %calculates the direction in degrees (for radians change to atan
    
    
    equi = dr + 90; %values perpendicular to the gradient, ie. eqipotential lines
    %for j = 1:zdim
       equi = conv2(equi, filter, 'same'); 
    %end
    
%     %% plotting
%     
%     slice = 7;
%     figure;
%     subplot(2,2,1);
%     imagesc(dr(:,:,slice));
%     title('direction of gradient in xy plane');
%     
%     subplot(2,2,2);
%     imagesc(equi(:,:,slice));
%     title('direction of equi in xy plane'); 
%     
%     
%     subplot(2,2,3);
%     imagesc(A(:,:,slice));
%     title('A'); 
%     
%     subplot(2,2,4);
%     imagesc(mag(:,:,slice));
%     title('magnitude of gradient in xy plane'); 