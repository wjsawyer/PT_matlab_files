% Copyright (C) 2012 Quan Wang <wangq10@rpi.edu>, 
% Signal Analysis and Machine Perception Laboratory, 
% Department of Electrical, Computer, and Systems Engineering, 
% Rensselaer Polytechnic Institute, Troy, NY 12180, USA

% This package implements the Gradient Vector Flow (GVF) in C++/MEX, which
% is significantly faster than the Matlab implementation. 
% Please run this script as a simple demo. 

% This implementation is used by the following paper: 
% Quan Wang, Kim L. Boyer, 
% The active geometric shape model: A new robust deformable shape model and its applications, 
% Computer Vision and Image Understanding, Volume 116, Issue 12, December 2012, Pages 1178-1194, 
% ISSN 1077-3142, 10.1016/j.cviu.2012.08.004. 

clear;clc;close all;

%% compile
mex GVF.cpp;

%% load image
I0=imread('image.png');
I=double(I0);

%% blur image
s=5;
I = gaussianBlur(I,s);

%% parameters
alpha=0.2;
mu=0.2;
iter=20;

%% run GVF
tic;
[u,v] = GVF(I, alpha, mu, iter);
t=toc;

fprintf('Computing GVF uses %f seconds \n',t);

%% visualize result
imshow(I0);
hold on;
quiver(u,v);
axis ij off;

