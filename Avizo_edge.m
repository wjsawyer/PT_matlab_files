% This script is designed to read in a CT image volume as well as a
% preliminary interpretation (label map) volume.  The script is written to
% read these volumes from AVIZO matlab export files, but could be modified
% to work with other volume interpretation/segmentation packages.  The
% script uses Canny edge detection to find edges within the grayscale
% volume and then finds the nearest edges to the preliminary interpretation
% (within user defined thresholds).
%
% Input files are grayscale image volume (A.mat) and binary label map
% volume (B.mat) with ones for material of interest.
%
% User input parameters: cellsize of the image volumes, 
% +/- distance thresholds (othresh/ithresh) to 
% determine how far from the intital interpretation the refined
% interpretation can stray, a mineral threshold (mt) that determines
% how far from the edge mineral gray scales will be analyzed and grayscale
% thresholds (bio, qf) for the mineral phases to be identified.
%
% Outputs include digital elevation models of the fault surfaces, label map
% volumes of these surfaces for QC work, interpreted mineralogy for the
% surfaces, and a thickness map.
% 
% note that matrix axes are transposed from image axes

% user defined paramaters
    % Voxel size in microns
    CELLSIZE=36.04;

    % thresholds (all values in voxels)
    othresh = 5; % threshold for finding edges outward from the surface
    ithresh = 2; % threshold for finding edges inward from the surface
    mt = 5; % outward distance at which mineral phases are identified
    
%     % Mineral ID grayscale threshold values
%     bio = 105; % 8-bit minimum value for biotite
%     qf = 90; % 8-bit maximum value for undiferentiated quartz/feldspar
    bio = 30000; % 16-bit minimum value for biotite
    qf = 26000; % 16-bit maximum value for undiferentiated quartz/feldspar

    
% read in data volumes and reformat to 3D from 4D

    % load data
    olddir = pwd;
    directory = uigetdir;
    cd(directory);
    % grayscale image volume
    aux = load('A.mat');
    %aux = load(strcat(directory, '\A.mat'));
    names=fieldnames(aux);
    Avizo_A_mat = aux.(names{1});
    % binary label map volume
    aux = load('B.mat');
    %aux = load(strcat(directory, '\B.mat'));
    names=fieldnames(aux);
    Avizo_B_mat = aux.(names{1});
    cd(olddir);
    clear aux names

    % find size of 3D volume
    [~,ydim, xdim, zdim] = size(Avizo_A_mat);
    y = 1:ydim; % define vector of rows
    x = 1:xdim; % define vector of columns

    % create output surface arrays
    topsurf = zeros(zdim,xdim);
    botsurf=topsurf;
    topmin=topsurf;
    botmin=topsurf;

    % convert 4D variables to 3D
    A = reshape(Avizo_A_mat(1,:,:,:),ydim,xdim,zdim);
    B = reshape(Avizo_B_mat(1,:,:,:),ydim,xdim,zdim);

    % create a matrix filled with row values (for distance calculation)
    row = repmat(y',1,xdim);

    
% loop through xy slices using edge detection to refine interpretation
    
    for j=1:zdim
        img=reshape(A(:,:,j),ydim,xdim);
        lbl=reshape(B(:,:,j),ydim,xdim);
        edg=edge(img,'canny');
        edg=double(edg);
        edg(edg==0)=NaN;

        % find initial guess row indices by column

            % top surface
            % use the sum of the cumulative product to calculate the number
            % of cells before the first zero
            idxt = mod(sum(cumprod(double(lbl ~= 1),1),1)+1,ydim+1);
            idxt(idxt==0) = NaN;

            % bottom
            idxb = ydim-mod(sum(cumprod(double(flipud(lbl ~= 1)),1),1)...
                +1,ydim+1)+1;
            idxb = idxb.*(idxb<=ydim); % change values > ydim to 0
            idxb(idxb==0) = NaN;

        % find nearest edge pixels
        
            % first find distance between all edge pixels and initial
            % interpretation
            
            % create matrices of top and bottom row indices repeating in
            % columns
            top = repmat(idxt, ydim, 1);
            bot = repmat(idxb, ydim, 1);

            % subtract initial surface location from rows to get in-row
            % distance
            top = row-top;
            bot = row-bot;

            % find distance to edge pixels
            top = top.*edg;
            bot = bot.*edg;
            
        % find nearest edges to intial interpretation within user-defined
        % thresholds.  Start with top surface.
            
            % Replace all edge values outside of thresholds with NaN
            top(top<-othresh | top>ithresh)=NaN;
            % find nearest remaining edge by column
            [~, I] = min(abs(top),[],1);
            % create a temporary surface from the minimum values
            temp = top(sub2ind([ydim,xdim],I,x));
            % add distance (correction) to initial surface elevations
            tempsurf = temp+idxt;
            % place row elevations into DEM matrix
            topsurf(j,:) = tempsurf;
            
        % identify mineral phases for each edge voxel
            
            mintest=NaN(size(tempsurf)); % initialize vector for results
            % find locations with data and get the gray values from
            % locations offset by the user identified threshold
%             mintest(~isnan(tempsurf)) = img(sub2ind([ydim,xdim],...
%                 tempsurf(~isnan(tempsurf))-mt,x(~isnan(tempsurf))));

            % with check for out of bounds values with mt
            mintest(tempsurf>mt) = img(sub2ind([ydim,xdim],...
                tempsurf(tempsurf>mt)-mt,x(tempsurf>mt)));
            
            % set values of temoprary vector  for each mineral phase
            mineral=zeros(1,xdim);
            mineral(mintest>bio)=2;
            mineral(mintest<qf)=1;
            % place row mineral values into mineral matrix
            topmin(j,:) = mineral;
            
        % find nearest edges to intial interpretation within user-defined
        % thresholds.  Repeat for bottom surface.

            % Replace all edge values outside of thresholds with NaN
            bot(bot<-ithresh | bot>othresh)=NaN;
            % find nearest remaining edge by column
            [~, I] = min(abs(bot),[],1);
            % create a temporary surface from the minimum values
            temp = bot(sub2ind([ydim,xdim],I,x));
            % add distance (correction) to initial surface elevations
            tempsurf = temp+idxb;
            % place row elevations into DEM matrix
            botsurf(j,:)=tempsurf;
            
        % identify mineral phases for each edge voxel            
            
            mintest=NaN(size(tempsurf)); % initialize vector for results
            % find locations with data and get teh gray values from
            % locations offset by the user identified threshold
            mintest(tempsurf<(ydim-mt)) = img(sub2ind([ydim,xdim],...
                tempsurf(tempsurf<(ydim-mt))+mt,x(tempsurf<(ydim-mt))));
            % set values of temoprary vector  for each mineral phase
            mineral=zeros(1,xdim);
            mineral(mintest>bio)=2;
            mineral(mintest<qf)=1;
            % place row mineral values into mineral matrix
            botmin(j,:) = mineral;
    end 

%replace NaN values with -9999 for Arc and save min surface
topmin=fliplr(topmin);
topmin(isnan(topmin))=-9999;
botmin(isnan(botmin))=-9999;

% axis equal
% save('./L12_06_MATLAB/matlab_out/L12_06_topmin.txt','topmin','-ascii')
% save('./L12_06_MATLAB/matlab_out/L12_06_botmin.txt','botmin','-ascii')
cd(directory);
save('./matlab_out/topmin.txt','topmin','-ascii')
save('./matlab_out/botmin.txt','botmin','-ascii')


% create new label map by creating 3D arrays of surface x coordinates and
% comparing to x coordinate matrix.
C=uint8(zeros(size(B)));
Ct=C;
Cb=C;

temp=topsurf';
DEM1 = temp(:,:,ones(ydim,1));
DEM1 = shiftdim(DEM1,2);

temp=botsurf';
DEM2 = temp(:,:,ones(ydim,1));
DEM2 = shiftdim(DEM2,2);

[~,Y,~] = meshgrid(x,y,[1:zdim]);
C(Y>=DEM1 & Y<=DEM2) = 1;
Ct(Y==DEM1) = 1;
Cb(Y==DEM2) = 1;

save('./matlab_out/C.mat','C')
save('./matlab_out/Ct.mat','Ct')
save('./matlab_out/Cb.mat','Cb')

%save raw version of topsurf/botsurf
rawtopsurf=topsurf;
rawtopsurf=rawtopsurf.*CELLSIZE;
rawtopsurf(isnan(rawtopsurf))=-9999;
rawbotsurf=botsurf;
rawbotsurf=rawbotsurf.*CELLSIZE;
rawbotsurf(isnan(rawbotsurf))=-9999;

% save('./L12_06_MATLAB/matlab_out/L12_06_rawt.txt','rawtopsurf','-ascii')
% save('./L12_06_MATLAB/matlab_out/L12_06_rawb.txt','rawbotsurf','-ascii')
save('./matlab_out/rawt.txt','rawtopsurf','-ascii')
save('./matlab_out/rawb.txt','rawbotsurf','-ascii')

%thickness, raw
th= botsurf-topsurf;
th=th.*CELLSIZE;
th(isnan(th))=-9999;
% save('./L12_06_MATLAB/matlab_out/L12_06_rawSN.txt','th','-ascii')
save('./matlab_out/rawth.txt','th','-ascii')
cd(olddir);

% subtract best fit plane and calculate dip angle
[Zm, ~, YCoeff] = fitplane_mtrx(topsurf);
topsurf = topsurf-round(Zm);
topdip = atan(YCoeff)*180/pi;
[Zm, ~, YCoeff] = fitplane_mtrx(botsurf);
botsurf = botsurf-round(Zm);
botdip = atan(YCoeff)*180/pi;

topsurf=fliplr(topsurf);
botsurf=-1*botsurf;

%replace NaN values with -9999 for Arc and save surface
topsurf=topsurf.*CELLSIZE;
botsurf=botsurf.*CELLSIZE;
topsurf(isnan(topsurf))=-9999;
botsurf(isnan(botsurf))=-9999;
cd(directory);
% save('./L12_06_MATLAB/matlab_out/L12_06_topsurf.txt', 'topsurf','-ascii')
% save('./L12_06_MATLAB/matlab_out/L12_06_botsurf.txt', 'botsurf', '-ascii')
save('./matlab_out/topsurf.txt', 'topsurf','-ascii')
save('./matlab_out/botsurf.txt', 'botsurf', '-ascii')

cd(olddir);

% header values for esri ascii import
% NCOLS 256 
% NROWS 256
% XLLCORNER 0
% YLLCORNER 0
% CELLSIZE 36.04
% NODATA_VALUE -9999

        