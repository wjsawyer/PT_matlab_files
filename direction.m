%calculates the angle of the line segment near an endpoint
function epvector = direction(bin)
%where bin is a 2d matrix
%epvector is a 3 by n matrix where n is the number of endpoints
[xdim, ydim] = size(bin); % gets dimensions of array
bin = bwareaopen(bin, 2); % removes single points (can't evaluate a slope)


%set up pieces array and endpoints array
cc = bwconncomp(bin);
pieces = (labelmatrix(cc));
endpoints = (bwmorph(bin, 'endpoints'));


%init epvector with length equal to number of endpoints
numendpoints = sum(sum(endpoints));
epvector = zeros(3,numendpoints);

%columns 1 and 2 are coordinates, column 3 is the gradient direction
i = 1;
for x = 1 : xdim 
    for y = 1 : ydim
        if endpoints(y,x) == 1
            epvector(1:2,i) = [x;y];
            epvector(3,i) = pieces(y,x);
            i = i + 1;
        end
    end
end

r = 5; %radius of square around endpoint to include in slope calculation.
            %only involves points on segment
for q = 1:numendpoints
    index = cc.PixelIdxList{1,(epvector(3,q))}(:, 1);
    [X, Y] = ind2vec(index, xdim); %should this be y dim?
   
    %done in 3 steps to ensure the (x,y) pairs keep the same col number
    % set to zero if outside range
    X(X<(epvector(1,q) - r) | X>(epvector(1,q) + r)) = 0;    
    Y(Y<(epvector(2,q) - r) | Y>(epvector(2,q) + r)) = 0;
    %set to zero if the other is outside range
    X(Y == 0) = 0;
    Y(X == 0) = 0;
    %remove zeros
    X(X == 0) = [];
    Y(Y == 0) = [];
    
    %fit a 1st order polynomial (line) to the points
    %the 1st coefficiant is slope, 2nd is intercept
    %throws warning b/c not all x may be distinct
    warning('off', 'MATLAB:polyfit:RepeatedPointsOrRescale');
    P  = polyfit(X,Y,1);
    warning('on', 'MATLAB:polyfit:RepeatedPointsOrRescale');

    slope = -(P(1)); %polyfit finds slope when plotting (where y increases up)
                        %want slope in array, where Y increases down
     %convert to angle                   
    epvector(3,q) = atand(slope);
end


%makes pieces
%label each piece with the int corrresponding to its spot in pixelidlist
%make endpoints

%for each endpoint, find what piece it is on
%call piece in index form from cc.PixelIdxList{1,N}(:, 1);

%turn index into x and y

%remove points not close enough to the endpoint

%do polyfit on the remainer and extract a slope

%convert slope to angle
%output epvector(3, length) where first two columns are coordinates, 3rd in the direction

    

