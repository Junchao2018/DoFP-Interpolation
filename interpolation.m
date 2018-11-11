function [img45, img135, img0, img90] = interpolation(img)

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



%%%%Micro-polarizator pattern
%%45 135
%%90  0
[rows, cols]=size(img);
img_inv=zeros(rows,cols);
img_errorh=zeros(rows,cols);
img_errorv=zeros(rows,cols);

H=[ 1   0   1;
    0   0   0;
    -1  0  -1];
V=[ 1   0  -1;
    0   0   0;
    1   0  -1];
img_errorh=abs(imfilter(img,H,'symmetric','same'));
img_errorv=abs(imfilter(img,V,'symmetric','same'));
%%
img_errorh(1,:)=0;
img_errorh(end,:)=0;
img_errorh(:,1)=0;
img_errorh(:,end)=0;
img_errorv(1,:)=0;
img_errorv(end,:)=0;
img_errorv(:,1)=0;
img_errorv(:,end)=0;
%%
[img45 img135 img0 img90]=gra_spline(img);
TB=[1/4 0 1/4;0 0 0;1/4 0 1/4];
img_inv=imfilter(img,TB,'same');
img_inv(1,:)=img_inv(3,:);
img_inv(end,:)=img_inv(end-2,:);
img_inv(:,1)=img_inv(:,3);
img_inv(:,end)=img_inv(:,end-2);
%%%%%%%%%%Interpolation in diagonal direction
aa=1;%1.1343;
for i=4:rows-3
    for j=4:cols-3
        d1=abs(img(i-1,j+1)-img(i-3,j+3))+abs(img(i+1,j-1)-img(i-1,j+1))+abs(img(i+3,j-3)-img(i+1,j-1))+...
            abs(img(i-1,j-1)-img(i-3,j+1))+abs(img(i+1,j-3)-img(i-1,j-1))+abs(img(i-1,j-3)-img(i-3,j-1))+...
            abs(img(i+1,j+1)-img(i-1,j+3))+abs(img(i+3,j-1)-img(i+1,j+1))+abs(img(i+3,j+1)-img(i+1,j+3))+...
            abs(img(i+2,j-2)+img(i-2,j+2)-2*img(i,j));%45
        d2=abs(img(i-1,j-1)-img(i-3,j-3))+abs(img(i+1,j+1)-img(i-1,j-1))+abs(img(i+3,j+3)-img(i+1,j+1))+...
            abs(img(i-3,j-1)-img(i-1,j+1))+abs(img(i+1,j+3)-img(i-1,j+1))+abs(img(i-3,j+1)-img(i-1,j+3))+...
            abs(img(i-1,j-3)-img(i+1,j-1))+abs(img(i+3,j+1)-img(i+1,j-1))+abs(img(i+3,j-1)-img(i+1,j-3))+...
            abs(img(i-2,j-2)+img(i+2,j+2)-2*img(i,j));%135
        if d1>aa*d2
            img_inv(i,j)=-(img(i-3,j-3)+img(i+3,j+3))/16+(img(i-1,j-1)+img(i+1,j+1))*9/16;
        else if aa*d1<d2
            img_inv(i,j)=-(img(i-3,j+3)+img(i+3,j-3))/16+(img(i-1,j+1)+img(i+1,j-1))*9/16;
            end
        end
    end
end

coeff=1;
for i=4:rows-3
    for j=4:cols-3
        verror=img_errorh(i,j)+img_errorh(i,j-2)+img_errorh(i,j+2)+...
                   img_errorh(i-2,j)+img_errorh(i-2,j-2)+img_errorh(i-2,j+2)+...
                   img_errorh(i+2,j)+img_errorh(i+2,j-2)+img_errorh(i+2,j+2);
        herror=img_errorv(i,j)+img_errorv(i,j-2)+img_errorv(i,j+2)+...
               img_errorv(i-2,j)+img_errorv(i-2,j-2)+img_errorv(i-2,j+2)+...
               img_errorv(i+2,j)+img_errorv(i+2,j-2)+img_errorv(i+2,j+2);
        if mod(i,2)~=0&&mod(j,2)~=0%45
            img45(i,j)=img(i,j);
            img0(i,j)=img_inv(i,j);
            if (verror/herror>coeff)%Horizonal
                img135(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                img90(i,j)=-3*(img_inv(i,j-3)+img_inv(i,j+3))/40+(img_inv(i,j+1)+img_inv(i,j-1))*23/40;
            else if (herror/verror>coeff)%Vertical
                    img135(i,j)=-3*(img_inv(i-3,j)+img_inv(i+3,j))/40+(img_inv(i+1,j)+img_inv(i-1,j))*23/40;
                    img90(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                else
                    img135(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                    img90(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                end
            end
        else if mod(i,2)~=0&&mod(j,2)==0%135
                img135(i,j)=img(i,j);
                img90(i,j)=img_inv(i,j);
                if (verror/herror>coeff)%Horizonal
                    img45(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                    img0(i,j)=-3*(img_inv(i,j-3)+img_inv(i,j+3))/40+(img_inv(i,j+1)+img_inv(i,j-1))*23/40;
                else if (herror/verror>coeff)%Vertical
                        img45(i,j)=-3*(img_inv(i-3,j)+img_inv(i+3,j))/40+(img_inv(i+1,j)+img_inv(i-1,j))*23/40;
                        img0(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                    else
                        img45(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                        img0(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                    end
                end
            else if mod(i,2)==0&&mod(j,2)==0%0
                    img0(i,j)=img(i,j);
                    img45(i,j)=img_inv(i,j);
                    if(verror/herror>coeff)%Horizonal
                        img135(i,j)=-3*(img_inv(i,j-3)+img_inv(i,j+3))/40+(img_inv(i,j+1)+img_inv(i,j-1))*23/40;
                        img90(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                    else if (herror/verror>coeff)%Vertical
                            img135(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                            img90(i,j)=-3*(img_inv(i-3,j)+img_inv(i+3,j))/40+(img_inv(i+1,j)+img_inv(i-1,j))*23/40;
                        else
                            img135(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                            img90(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                        end
                    end
                else%%%%90
                    img90(i,j)=img(i,j);
                    img135(i,j)=img_inv(i,j);
                    if (verror/herror>coeff)%Horizonal
                        img45(i,j)=-3*(img_inv(i,j-3)+img_inv(i,j+3))/40+(img_inv(i,j+1)+img_inv(i,j-1))*23/40;
                        img0(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                    else if (herror/verror>coeff)%Vertical
                            img45(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                            img0(i,j)=-3*(img_inv(i-3,j)+img_inv(i+3,j))/40+(img_inv(i+1,j)+img_inv(i-1,j))*23/40;
                        else
                            img45(i,j)=-3*(img(i-3,j)+img(i+3,j))/40+(img(i+1,j)+img(i-1,j))*23/40;
                            img0(i,j)=-3*(img(i,j-3)+img(i,j+3))/40+(img(i,j+1)+img(i,j-1))*23/40;
                        end
                    end
                end
            end
        end
    end
end
end

