function z = autolabelplot(n)
%define which z-slice should be plotted
A = evalin('base', 'A');
edge = evalin('base', 'edge');
label = evalin('base', 'label');
label2 = evalin('base', 'label2');
label3 = evalin('base', 'label3');
llimit = evalin('base', 'llimit');
ulimit = evalin('base', 'ulimit');
edgethresh = evalin('base', 'edgethresh');
xx = evalin('base', 'xx');
yy = evalin('base', 'yy');
zz = evalin('base', 'zz');

z = n;
figure;
%plot of labelfield
subplot(2,3,1); 
imagesc(A(:,:,z));
title(['grayscale image from slice ' num2str(z)]);
colorbar;

%plot of labelfield
subplot(2,3,2); 
imagesc(edge(:,:,z));
title(['edge intensities where dist in x is' num2str(xx) ' and y is '  num2str(yy) ' and z is ' num2str(zz) ' from slice ' num2str(z)]);
colorbar;

%plot of labelfield
subplot(2,3,3); 
imagesc(label(:,:,z));
title(['labelfield: guess PT. brightness between ' num2str(llimit) ' and ' num2str(ulimit) ' from slice ' num2str(z)]);
colorbar;

%plot of labelfield
subplot(2,3,4); 
imagesc(label2(:,:,z));
title(['labelfield minus points where log(edge intensity)> ' num2str(edgethresh) ' from slice ' num2str(z)]);
colorbar;


%plot of labelfield
subplot(2,3,5); 
imagesc(label3(:,:,z));
title(['eroded labelfield from slice ' num2str(z)]);
colorbar;








end
