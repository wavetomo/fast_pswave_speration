
%%Fast P/S wavefield decomposition for vertically transversely isotropic media based on wavefront phase direction
clear;clc;

%% example 1 : homogeneous model
 
% nx=600;nz=600;dx=10;dz=10; pml=0;
% 
% vxfile='homogeneous_model_snapshot_600_600_0700ms.vx';
% vzfile='homogeneous_model_snapshot_600_600_0700ms.vz';
% 
% 
% epsilon=ones(nz,nx)*0.3;
% delta=ones(nz,nx)*0.25;
% vp=ones(nz,nx)*3000;
% vs=ones(nz,nx)*1765;
% rho=ones(nz,nx)*1000;

%% example 2 : layered model

% nx=860;nz=460;  dx=10;dz=10; pml=30;
% 
% vxfile='layered_model_snapshot_460_860_1100ms.vx';
% vzfile='layered_model_snapshot_460_860_1100ms.vz';
% 
% vp=read_matrix('.\layer model\vp_400x800_10m.bin',400,800);
% vp=imgaussfilt(vp,3);
% vs=read_matrix('.\layer model\vs_400x800_10m.bin',400,800);
% vs=imgaussfilt(vs,3);
% rho=read_matrix('.\layer model\rho_400x800_10m.bin',400,800);
% rho=imgaussfilt(rho,3);
% epsilon=read_matrix('.\layer model\epsilon_400x800_10m.bin',400,800);
% epsilon=imgaussfilt(epsilon,3);
% delta=read_matrix('.\layer model\delta_400x800_10m.bin',400,800);
% delta=imgaussfilt(delta,3);


%% example 3 : Hess model

nx=1560;nz=860; pml=30; dx=10;dz=10; 
vxfile='Hess_model_snapshot_860_1560_2600ms.vx';
vzfile='Hess_model_snapshot_860_1560_2600ms.vz';

vp=read_matrix('.\Hess model\hess_vp_10m.bin',915,2206);
vp=vp(1:800,1:1500);
vp=imgaussfilt(vp,3);
vs=read_matrix('.\Hess model\hess_vs_10m.bin',915,2206);
vs=vs(1:800,1:1500);
vs=imgaussfilt(vs,3);
rho=read_matrix('.\Hess model\hess_rho_10m.bin',915,2206);
rho=rho(1:800,1:1500);
rho=imgaussfilt(rho,3);
epsilon=read_matrix('.\Hess model\hess_epsilon_10m.bin',915,2206);
epsilon=epsilon(1:800,1:1500);
epsilon=imgaussfilt(epsilon,3);
delta=read_matrix('.\Hess model\hess_delta_10m.bin',915,2206);
delta=delta(1:800,1:1500);
delta=imgaussfilt(delta,3);



vp=modpad2d(vp,pml,nz,nx);
vs=modpad2d(vs,pml,nz,nx);
rho=modpad2d(rho,pml,nz,nx);
epsilon=modpad2d(epsilon,pml,nz,nx);
delta=modpad2d(delta,pml,nz,nx);



temp_vx=read_matrix(vxfile,nz,nx);
temp_vz=read_matrix(vzfile,nz,nx);

% %-shift wavefield---
[kz,kx]=generate_wavenumber(nz,nx,dz,dx);
shift_x=exp(1i*(dx/2).*kx);
shift_z=exp(1i*(-dz/2).*kz);
if mod(nx,2)==0
    shift_x(:,nx/2+1) = real(shift_x(:,nx/2+1));
end
if mod(nz,2)==0
    shift_z(nz/2+1,:) = real(shift_z(nz/2+1,:));
end

jkx=1i.*kx;
jkz=1i.*kz;
if mod(nx,2)==0
    jkx(:,nx/2+1) = real(jkx(:,nx/2+1));
end
if mod(nz,2)==0
    jkz(nz/2+1,:) = real(jkz(nz/2+1,:));
end

vx=temp_vx;

vz=ifft(shift_z.*fft(temp_vz,[],1),[],1);
vz=ifft(shift_x.*fft(vz,[],2),[],2);  

tic;
%-----------------
r1=(1+2*epsilon).*vp.*vp - vs.*vs;
r2=sqrt(((1+2*delta).*vp.*vp - vs.*vs).*(vp.*vp - vs.*vs));
r3=vp.*vp - vs.*vs;
r4=2*(delta-epsilon).*vp.*vp.*(vp.*vp - vs.*vs);
r5=vp.*vp./(vp.*vp - vs.*vs);
% compute direction for vx
gz_vx=derivate1_fd8(vx,1,dz); % d vx/dz
gx_vx=derivate1_fd8(vx,2,dx); % d vx/dx

nxs_vx = zeros(nz,nx);
nzs_vx = zeros(nz,nx);
rx_vx = zeros(nz,nx);
rz = zeros(nz,nx);
rr_vx = gx_vx.*gx_vx + gz_vx.*gz_vx;
for ix=1:nx
    for iz=1:nz        
        if rr_vx(iz,ix)==0.0
            nxs_vx(iz,ix) = 0.0;
            nzs_vx(iz,ix) = 0.0;
        else
            nxs_vx(iz,ix) = gx_vx(iz,ix)*gx_vx(iz,ix)/(gx_vx(iz,ix)*gx_vx(iz,ix) + gz_vx(iz,ix)*gz_vx(iz,ix)); % nx*nx
            nzs_vx(iz,ix) = gz_vx(iz,ix)*gz_vx(iz,ix)/(gx_vx(iz,ix)*gx_vx(iz,ix) + gz_vx(iz,ix)*gz_vx(iz,ix)); % nz*nz
        end
    end
end
nxs_vx = imgaussfilt(nxs_vx,5);
nzs_vx = imgaussfilt(nzs_vx,5);


%% compute direction for vz
gz_vz=derivate1_fd8(vz,1,dz); % d vz/dz
gx_vz=derivate1_fd8(vz,2,dx); % d vz/dx

nxs_vz = zeros(nz,nx);
nzs_vz = zeros(nz,nx);
rx_vz = zeros(nz,nx);
rr_vz = gx_vz.*gx_vz + gz_vz.*gz_vz;
for ix=1:nx
    for iz=1:nz        
        if rr_vz(iz,ix)==0.0
            nxs_vz(iz,ix) = 0.0;
            nzs_vz(iz,ix) = 0.0;
        else
            nxs_vz(iz,ix) = gx_vz(iz,ix)*gx_vz(iz,ix)/(gx_vz(iz,ix)*gx_vz(iz,ix) + gz_vz(iz,ix)*gz_vz(iz,ix)); % nx*nx
            nzs_vz(iz,ix) = gz_vz(iz,ix)*gz_vz(iz,ix)/(gx_vz(iz,ix)*gx_vz(iz,ix) + gz_vz(iz,ix)*gz_vz(iz,ix)); % nz*nz
        end
    end
end
nxs_vz = imgaussfilt(nxs_vz,5);
nzs_vz = imgaussfilt(nzs_vz,5);

%% Solve the auxiliary wavefield w

Ux=fft2(vx);
Uz=fft2(vz);

Ux_1=Ux./(kx.^2+kz.^2);
Uz_1=Uz./(kx.^2+kz.^2);
Ux_1(1,1)=0; Uz_1(1,1)=0;

Ux_2=kx.^2.*kz.^2.*Ux./(kx.^2+kz.^2).^3;
Uz_2=kx.^2.*kz.^2.*Uz./(kx.^2+kz.^2).^3;
Ux_2(1,1)=0; Uz_2(1,1)=0;


Ux_3=(kx.^2.*kz.^2-kz.^4).*Ux./(kx.^2+kz.^2).^3;
Uz_3=(kx.^2.*kz.^2-kz.^4).*Uz./(kx.^2+kz.^2).^3;
Ux_3(1,1)=0; Uz_3(1,1)=0;



ux1=ifft2(Ux_1,'symmetric');  uz1=ifft2(Uz_1,'symmetric');
ux2=ifft2(Ux_2,'symmetric');  uz2=ifft2(Uz_2,'symmetric');
ux3=ifft2(Ux_3,'symmetric');  uz3=ifft2(Uz_3,'symmetric');


alpha=1e-20;

% Hemholtz decomposition

r_vx=r2./(r1+r4.*nzs_vx./(r1.*nxs_vx + r3.*nzs_vx + alpha));
r_vz=r2./(r1+r4.*nzs_vz./(r1.*nxs_vz + r3.*nzs_vz + alpha));

ux1dxx = derivate2_fd8(ux1,2,dx);
uz1dxx = derivate2_fd8(uz1,2,dx);

ux1dzz = derivate2_fd8(ux1,1,dz);
uz1dzz = derivate2_fd8(uz1,1,dz);

ux1dx = derivate1_fd8(ux1,2,dx);
ux1dxz = derivate1_fd8(ux1dx,1,dz);

uz1dx = derivate1_fd8(uz1,2,dx);
uz1dxz = derivate1_fd8(uz1dx,1,dz);

vxp1 = ux1dxx + r_vz.*uz1dxz;

vzp1 = r_vx.*ux1dxz + r_vz.^2.*uz1dzz;

vxs1 = r_vx.^2.*ux1dzz - r_vz.*uz1dxz;

vzs1 = uz1dxx - r_vx.*ux1dxz;



ux2dxx = derivate2_fd8(ux2,2,dx);
uz2dxx = derivate2_fd8(uz2,2,dx);

ux2dzz = derivate2_fd8(ux2,1,dz);
uz2dzz = derivate2_fd8(uz2,1,dz);

ux2dx = derivate1_fd8(ux2,2,dx);
ux2dxz = derivate1_fd8(ux2dx,1,dz);

uz2dx = derivate1_fd8(uz2,2,dx);
uz2dxz = derivate1_fd8(uz2dx,1,dz);

vxp2 = ux2dxx + r_vz.*uz2dxz;

vzp2 = r_vx.*ux2dxz + r_vz.^2.*uz2dzz;

vxs2 = r_vx.^2.*ux2dzz - r_vz.*uz2dxz;

vzs2 = uz2dxx - r_vx.*ux2dxz;



ux3dxx = derivate2_fd8(ux3,2,dx);
uz3dxx = derivate2_fd8(uz3,2,dx);

ux3dzz = derivate2_fd8(ux3,1,dz);
uz3dzz = derivate2_fd8(uz3,1,dz);

ux3dx = derivate1_fd8(ux3,2,dx);
ux3dxz = derivate1_fd8(ux3dx,1,dz);

uz3dx = derivate1_fd8(uz3,2,dx);
uz3dxz = derivate1_fd8(uz3dx,1,dz);

vxp3 = ux3dxx + r_vz.*uz3dxz;

vzp3 = r_vx.*ux3dxz + r_vz.^2.*uz3dzz;

vxs3 = r_vx.^2.*ux3dzz - r_vz.*uz3dxz;

vzs3 =uz3dxx - r_vx.*ux3dxz;


vxp = -(vxp1 + 4*epsilon.*r5.*vxp2 - 2*delta.*r5.*vxp3);
vzp = -(vzp1 + 4*epsilon.*r5.*vzp2 - 2*delta.*r5.*vzp3);
vxs = -(vxs1 + 4*epsilon.*r5.*vxs2 - 2*delta.*r5.*vxs3);
vzs = -(vzs1 + 4*epsilon.*r5.*vzs2 - 2*delta.*r5.*vzs3);




toc;

vx=vx(pml+1:nz-pml,pml+1:nx-pml);
vz=vz(pml+1:nz-pml,pml+1:nx-pml);
vxp=vxp(pml+1:nz-pml,pml+1:nx-pml);
vxs=vxs(pml+1:nz-pml,pml+1:nx-pml);
vzp=vzp(pml+1:nz-pml,pml+1:nx-pml);
vzs=vzs(pml+1:nz-pml,pml+1:nx-pml);



[colormax0,colormin0] = perclip([vx;vz;vxs;vxp;vzs;vzp],0.004);  
x=(0:nx-2*pml)*dx/1000;
z=(0:nz-2*pml)*dz/1000;

figure(1);imagesc(x,z,vx);shading interp;sucolormap(gca,0);
clim([colormin0, colormax0]);
set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
title('Horizontal','FontName','Times New Roman','FontSize',11);

figure(2);imagesc(x,z,vxp);shading interp;sucolormap(gca,0);
clim([colormin0, colormax0]);
set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
title('P Horizontal','FontName','Times New Roman','FontSize',11);


figure(3);imagesc(x,z,vxs);shading interp;sucolormap(gca,0);
clim([colormin0, colormax0]);
set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
title('S Horizontal','FontName','Times New Roman','FontSize',11);

figure(4);imagesc(x,z,vz);shading interp;sucolormap(gca,0);
clim([colormin0, colormax0]);
set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
title('Vertical','FontName','Times New Roman','FontSize',11);

figure(5);imagesc(x,z,vzp);shading interp;sucolormap(gca,0);
clim([colormin0, colormax0]);
set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
title('P Vertical','FontName','Times New Roman','FontSize',11);


figure(6);imagesc(x,z,vzs);shading interp;sucolormap(gca,0);
clim([colormin0, colormax0]);
set(gca,'TickDir','out');set(gca,'Ydir','reverse');set(gca,'Units','centimeters');
set(gca,'FontName','Times New Roman','FontSize',10);
xlabel('Distance (km)','FontName','Times New Roman','FontSize',11);
ylabel('Depth (km)','FontName','Times New Roman','FontSize',11);
title('S Vertical','FontName','Times New Roman','FontSize',11);
