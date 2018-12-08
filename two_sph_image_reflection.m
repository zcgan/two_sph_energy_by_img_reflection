
function [energy, img_num]=two_sph_image_reflection(Q1,Q2,a1,a2,d1,d2,eps_1,eps_2)
%two dielectric spheres, with central charges. 
%Image reflection method, based on the analytical image charge solution for a single dielectric sphere.
%The line integral of the image charge solution is discretized using the
%Gauss-Jacobi Quadrature rule.
%the detailed reflection parameters can be modified at lines 54-60;


%suppose we have two spheres with
% radius a1 and a2, central charge Q1 and Q2, locating at (0,0,d1) and (0,0,-d2)
%the contact distance is delta=d1+d2-a1-a2. The outside dielectric eps_o=1.0. while the
%inside dielectric eps_i is a variable.
%note that for simplicity, please set d1>a1, d2>a2.

%%%%input variables%%%%%%%%%%%
%Q1 central charge of sphere #1
%Q2 central charge of sphere #2
% 
% a1 radius of sphere #1
% a2 radius of sphere #2
% 
% d1 z-coordiante of sphere #1, set d1>a1
% -d2 z-coordiante of sphere #2 (d2>0), set d2>a2
% 
% eps_i dielectric constant inside the sphere #i


%%%%%%output variables%%%%%%%%%
%energy: the interaction energy between the two dielectric spheres with
%central charges
%img_num: total number of image charges used.
energy=0.0;
img_num=0;

%first check the overlap condition
if(d1<a1|d2<a2)
    disp('Exiting...Please set d1>a1 and d2>a2 to make sure the spheres do not overlap!');
    return;
end

global img; %vector storing the location and charge for all the images
global img_counter; % the counter for image charges while generating them
format long;
eps_o=1.0; %dielectric outside For simplicity. Do NOT change it. Rescale your result for other value of eps_o

%%%%%%%set the source and image charge array%%%%%%%

src=zeros(2,4); %2 charges, each has its x,y,z coordinates and charge q
img_max=10000; % the maximum number of images allowed.

img=zeros(img_max,4); %the array storing the image location and charges

img_counter=1; %integer for counting the number of images

%%%initialize src array%%%
src(1,1)=0;src(1,2)=0;src(1,3)=d1;src(1,4)=Q1;
src(2,1)=0;src(2,2)=0;src(2,3)=-d2;src(2,4)=Q2;

%%%initialize the image generation process%%%%
 nlevel=5; %image series is reflected nlevel times (higher levels are neglected)
 order=zeros(1,nlevel); %order is the number of discretized Gauss-Jacobi quadrature points in each level of reflection
 order(1)=5; %first level, set order =5, i.e., 5+1 Kelvin image charges for each source.
 order(2)=4; %second level, set order =4
 for i=3:nlevel
     order(i)=3; %higher levels, set order =3 
 end
 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%generate the first level images due to the two point sources%%%
    for i=1:2
        pos=src(i,:);
        j=3-i;
        generate_img(pos,order(1),Q1,Q2,a1,a2,d1,d2,eps_1,eps_2)
    end
    
  
%%%generate the higher level images due to reflection%%%%%%%%%%%%%
%%%image reflection up to nlevel%%%%%%%%%
 img_ind_start=1;
 img_ind_end=img_counter-1;
for i=2:nlevel
    for j=img_ind_start:img_ind_end
        pos=img(j,:);
        generate_img(pos,order(i),Q1,Q2,a1,a2,d1,d2,eps_1,eps_2)
    end
    img_ind_start=img_ind_end+1;
    img_ind_end=img_counter-1;
end
  
%%%%%%%%%%%%%%output%%%%%%%%%%%%%%%%%%

img_num=img_counter-1; %total number of images
%img(1:22,:)

%%%%energy based on the definition that E=0.5*\sum q_i\phi_i

%%%Coulomb part%%%%%%
dx=0.0-0.0;
dy=0.0-0.0;
dz=d1-(-d2);
dr=sqrt(dx*dx+dy*dy+dz*dz);
energy=energy+Q1*Q2/dr;
%%%Source-image interaction (i.e. the polarization part)%%%
%%%for first sphere%%%
for j=1:img_num
dx=img(j,1)-src(1,1);
dy=img(j,2)-src(1,2);
dz=img(j,3)-src(1,3);
dr=sqrt(dx*dx+dy*dy+dz*dz);
if(dr>a1) %only needs the image-source contribution due to the other spheres
energy=energy+0.5*Q1*img(j,4)/dr;
end 
end

%%%for second sphere%%%
for j=1:img_num
dx=img(j,1)-src(2,1);
dy=img(j,2)-src(2,2);
dz=img(j,3)-src(2,3);
dr=sqrt(dx*dx+dy*dy+dz*dz);
if(dr>a2)
energy=energy+0.5*Q2*img(j,4)/dr;
end
end


end

function generate_img(pos,order,Q1,Q2,a1,a2,d1,d2,eps_1,eps_2)
global img;
global img_counter;
%constants related to eps_i, assuming eps_o=1
eps_i=zeros(1,2); lamda=zeros(1,2);beta=zeros(1,2);
gamma=zeros(1,2);

eps_i(1)=eps_1;
eps_i(2)=eps_2;

lamda=1.0./(1.0+eps_i);
beta=lamda-1.0;
gamma=(eps_i-1.0)./(1.0+eps_i);

%%%image charge in the first sphere%%%
dx=pos(1)-0.0;
dy=pos(2)-0.0;
dz=pos(3)-d1;
dr=sqrt(dx*dx+dy*dy+dz*dz);

if(dr>a1)
%Kelvin image first%
rk=a1*a1/dr;
rim=rk/dr;
img(img_counter,1)=rim*dx+0.0;
img(img_counter,2)=rim*dy+0.0;
img(img_counter,3)=rim*dz+d1;
img(img_counter,4)=-gamma(1)*a1*pos(4)/dr; 		

img_counter=img_counter+1;
%%%%%Get Gauss quadrature weights and locations%%%%%%%%
%order=2;%number of Gauss pts
alpha=0.0;%always zero
a=0;b=rk; %a is always 0
x=zeros(order,1);w=zeros(order,1);
[x,w]=jacobi_rule_xw (order, alpha, beta(1), a, b);
for i=1:order
img(img_counter,1)=0+x(i)*dx/dr;
img(img_counter,2)=0+x(i)*dy/dr;
img(img_counter,3)=d1+x(i)*dz/dr;
img(img_counter,4)=w(i)*gamma(1)*lamda(1)*rk^(1-lamda(1))*pos(4)/a1;
img_counter=img_counter+1;
end						
end

%%%image charge in the second sphere%%%
dx=pos(1)-0.0;
dy=pos(2)-0.0;
dz=pos(3)-(-d2);
dr=sqrt(dx*dx+dy*dy+dz*dz);

if(dr>a2)
%Kelvin image first%
rk=a2*a2/dr;
rim=rk/dr;
img(img_counter,1)=rim*dx+0.0;
img(img_counter,2)=rim*dy+0.0;
img(img_counter,3)=rim*dz+(-d2);
img(img_counter,4)=-gamma(2)*a2*pos(4)/dr;	

img_counter=img_counter+1;
%%%%%Get Gauss quadrature weights and locations%%%%%%%%
%order=2;%number of Gauss pts
alpha=0.0;%always zero
a=0;b=rk; %a is always 0
x=zeros(order,1);w=zeros(order,1);               
[x,w]=jacobi_rule_xw (order, alpha, beta(2), a, b); 
for i=1:order
img(img_counter,1)=0+x(i)*dx/dr;
img(img_counter,2)=0+x(i)*dy/dr;
img(img_counter,3)=-d2+x(i)*dz/dr;
img(img_counter,4)=w(i)*gamma(2)*lamda(2)*rk^(1-lamda(2))*pos(4)/a2;
img_counter=img_counter+1;
end						
end

end
