

%threshold limits for PT color
llimit = 27500;
ulimit = 31000;

    olddir = pwd;
    if exist('directory','var')~=1
    directory = uigetdir;
    end
    cd(directory);
    aux = load('A.mat');
    names=fieldnames(aux);
    Avizo_A_mat = aux.(names{1});
    clear aux names
    % find size of 3D volume
    [~,xdim, ydim, zdim] = size(Avizo_A_mat);

    A = reshape(Avizo_A_mat(1,:,:,:),xdim,ydim,zdim);
    
    %this will produce false positives along edges between lighter and
    %darker materials.
    label = zeros(size(A));
    label(A<ulimit & A>llimit) = 1;
    edge = zeros(size(A));
    edgethresh = 16; % if log(edge)>this, them will be set to 0
    size = 2000; %minimum size of piece to keep. 
    
    
    %values for checking change in brightness - defined distance to check
    %in given direction. 
    xx = 3;
    yy = 3;
    zz = 3;
   
    %%define a matrix of equal values to A with a 'buffer region' on all
    %%edges of width defined above so that calculations of the change in
    %%brightness can be done for all voxels in initial space. 
    AA = zeros((xdim+(2*xx)), (ydim+(2*yy)), (zdim+(2*zz)));
    for x = 1:xdim %(xx + 1):(xdim + xx)
       for y = 1:ydim %(yy + 1):(ydim + yy)
           for z = 1:zdim %(zz + 1):(zdim + zz)
               AA(x+xx,y+yy,z+zz) = A(x,y,z);
           end
       end
    end
    for x = 1:xx
    AA(x,:,:) = AA(xx+1,:,:);
    AA(((xdim + (2 * xx)) + 1 - x),:,:) = AA((xdim + xx),:,:);
    end
    for y = 1:yy
    AA(:,y,:) = AA(:,yy+1,:);
    AA(:,((ydim + (2 * yy)) + 1 - y),:) = AA(:,(ydim + yy),:);
    end
    for z = 1:zz
    AA(:,:,z) = AA(:,:,zz+1);
    AA(:,:,((zdim + (2 * zz)) + 1 - z)) = AA(:,:,(zdim + zz));
    end
    
    %fill in values for slope using essentially a simple edge detection
    %method. this provides a criteria for identifying areas which may have
    %been found as false positives during the label generation
    for x = (xx + 1):(xdim + xx)
       for y = (yy + 1):(ydim + yy)
           for z = (zz + 1):(zdim + zz)
     edge(x-xx,y-xx,z-zz) = abs(AA(x-xx,y,z) - AA(x+xx,y,z))^2 +....
                             abs(AA(x,y-yy,z) - AA(x,y+yy,z))^2 +....
                             abs(AA(x,y,z+zz) - AA(x,y,z+zz))^2;

           end
       end
    end
    
    %label2 takes raw label and subtracts areas where edges are
    label2 = label;
    label2(log(edge) > edgethresh) = 0;
    
    
    %label3 takes raw label and applies an erosion filter
    SE = strel('disk',3);
    label3 = imopen(label, SE);
  
    
    %label4 deletes groups of voxels smaller than a thresh size
    %should be used after some processing so that everything isn't
    %connected
      cc = bwconncomp(label2);
      numobj = cc.NumObjects;
label4 = bwareaopen(label3, size);
    
    
    save('./autolabel.mat','label2')
    
    cd(olddir);
    
    %DIFFERENT APPROACH
    %find edges that have an average value of x on one side, where x is the
    %value of PT
 