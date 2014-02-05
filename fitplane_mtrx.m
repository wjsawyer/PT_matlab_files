 function [Zm, XCoeff, YCoeff] = fitplane_mtrx(Z)
% function [Zm] = fitplane_mtrx(Z)
% Generates a best fit plane through a matrix of elevation data
% assumes inf for no data values

% find data points
Z(isnan(Z))=0;
[I,J,K]=find(Z);
[ydim,xdim]=size(Z);

Xcolv = J(:); % Make X a column vector
Ycolv = I(:); % Make Y a column vector
Zcolv = K(:); % Make Z a column vector
Const = ones(size(Xcolv)); % Vector of ones for constant term

Coefficients = [Xcolv Ycolv Const]\Zcolv; % Find the coefficients
XCoeff = Coefficients(1); % X coefficient
YCoeff = Coefficients(2); % X coefficient
CCoeff = Coefficients(3); % constant term
% Using the above variables, z = XCoeff * x + YCoeff * y + CCoeff

% L=plot3(x,y,z,'r.'); % Plot the original data points
% %set(L,'Markersize',2*get(L,'Markersize')) % Making the circle markers larger
% %set(L,'Markerfacecolor','r') % Filling in the markers
% hold on
[xm, ym]=meshgrid(1:xdim,1:ydim); % Generating a regular grid for plotting
Zm = XCoeff * xm + YCoeff * ym + CCoeff;
% surf(xx,yy,zz) % Plotting the surface
% title(sprintf('Plotting plane z=(%f)*x+(%f)*y+(%f)',XCoeff, YCoeff, CCoeff))
% By rotating the surface, you can see that the points lie on the plane
% Also, if you multiply both sides of the equation in the title by 4,
% you get the equation in the comment on the third line of this example