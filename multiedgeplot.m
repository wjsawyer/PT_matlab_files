function z = multiedgeplot(n)
%define which z-slice should be plotted
B = evalin('base', 'B');
gradedges = evalin('base', 'gradedges');
canny = evalin('base', 'canny');
binedges = evalin('base', 'binedges');
xgrad = evalin('base', 'xgrad');
ygrad = evalin('base', 'ygrad');

z = n;
figure(1);
%plot of labelfield
subplot(2,3,1); 
imagesc(B(:,:,z));
title('initial label field');
colorbar;

%plot of labelfield
subplot(2,3,2); 
imagesc(canny(:,:,z));
title('canny edge detection');
colorbar;

%plot of labelfield
subplot(2,3,3); 
imagesc(binedges(:,:,z));
title('final binary edges');
colorbar;

%plot of labelfield
subplot(2,3,4); 
imagesc(xgrad(:,:,z));
title('gradient produced in x direction from label field');
colorbar;

%plot of labelfield
subplot(2,3,5); 
imagesc(ygrad(:,:,z));
title('gradient produced in y direction from label field');
colorbar;

%plot of canny 
subplot(2,3,6); 
imagesc(gradedges(:,:,z));
title('sum of x and y gradients');
colorbar;







end
