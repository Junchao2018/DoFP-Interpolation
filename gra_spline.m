function [img1 img2 img3 img4] = gra_spline(img)

% =========================================================================
% Image interpolation for DoFP, Version 1.0
% Copyright(c) 2016 Junchao Zhang, Haibo Luo, Bin Hui and Zheng Chang
% All Rights Reserved.
%
% ----------------------------------------------------------------------
% Permission to use, copy, or modify this software and its documentation
% for educational and research purposes only and without fee is here
% granted, provided that this copyright notice and the original authors'
% names appear on all copies and supporting documentation. This program
% shall not be used, rewritten, or adapted as the basis of a commercial
% software or hardware product without first obtaining permission of the
% authors. The authors make no representations about the suitability of
% this software for any purpose. It is provided "as is" without express
% or implied warranty.
%----------------------------------------------------------------------
%
% This is an implementation of the algorithm for image interpolation
% 
% Please cite the following paper if you use this code:
%
% Junchao Zhang, Haibo Luo, Bin Hui and Zheng Chang, "Image interpolation 
% for division of focal plane polarimeters with intensity correlation," 
% Opt. Express 24, 20799-20807 (2016)
% 
%--------------------------------------------------------------------------
[m n]=size(img);
img1=zeros(m,n);
img2=zeros(m,n);
img3=zeros(m,n);
img4=zeros(m,n);
img1(1:2:m,1:2:n)=img(1:2:m,1:2:n);
img2(1:2:m,2:2:n)=img(1:2:m,2:2:n);
img3(2:2:m,2:2:n)=img(2:2:m,2:2:n);
img4(2:2:m,1:2:n)=img(2:2:m,1:2:n);
%%
x=1:2:n;
y=img1(1:2:m,1:2:n);
pp=csapi(x,y);
img1(1:2:m,1:n-1)=fnval(pp,1:n-1);
x1=1:2:m;
y1=img1(1:2:m,1:n-1);
y1=y1';
pp=csapi(x1,y1);
img1(1:m-1,1:n-1)=(fnval(pp,1:m-1))';
img1(m,:)=img1(m-1,:);
img1(:,n)=img1(:,n-1);
%%
x=2:2:n;
y=img2(1:2:m,2:2:n);
pp=csapi(x,y);
img2(1:2:m,2:n)=fnval(pp,2:n);
x1=1:2:m;
y1=img2(1:2:m,2:n);
y1=y1';
pp=csapi(x1,y1);
img2(1:m-1,2:n)=(fnval(pp,1:m-1))';
img2(m,:)=img2(m-1,:);
img2(:,1)=img2(:,2);
%%
x=2:2:n;
y=img3(2:2:m,2:2:n);
pp=csapi(x,y);
img3(2:2:m,2:n)=fnval(pp,2:n);
x1=2:2:m;
y1=img3(2:2:m,2:n);
y1=y1';
pp=csapi(x1,y1);
img3(2:m,2:n)=(fnval(pp,2:m))';
img3(1,:)=img3(2,:);
img3(:,1)=img3(:,2);
%%
x=1:2:n;
y=img4(2:2:m,1:2:n);
pp=csapi(x,y);
img4(2:2:m,1:n-1)=fnval(pp,1:n-1);
x1=2:2:m;
y1=img4(2:2:m,1:n);
y1=y1';
pp=csapi(x1,y1);
img4(2:m,1:n)=(fnval(pp,2:m))';
img4(1,:)=img4(2,:);
img4(:,n)=img4(:,n-1);
end

