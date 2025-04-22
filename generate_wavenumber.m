function [kz,kx]=generate_wavenumber(nz,nx,dz,dx)
 %====================================================================
 %    [kz,kx]=generate_wavenumber(nz,nx,dz,dx)
 %    create wavenumber kx & kz
 %    nx: number of element along x (2nd) direction
 %    nz: number of element along z (1st) direction
 %    dx: interval along x (2nd) direction 
 %    dz: interval along z (1st) direction 
 %    return value:
 %    kx: kx array with size of 'nz x nx'  
 %    kz: kz array with size of 'nz x nx'  
 %====================================================================
    % kx
    deltF=1/(nx*dx);
    %==========
    wx=2*pi*deltF;
    index_max_wavenumber=fix(nx/2);
    k=zeros(nz,nx);
    k(1,1:index_max_wavenumber+1)=0:index_max_wavenumber;%求0：N/2的波数
    for j=index_max_wavenumber+2:nx
       k(1,j)=-(nx-(j-1));%N/2+1:N-1的波数
    end
    for iz=2:nz
        k(iz,:)=k(1,:);
    end
    kx=wx*k;
    % kz
    deltF=1/(nz*dz);
    %==========
    wz=2*pi*deltF;
    index_max_wavenumber=fix(nz/2);
    k=zeros(nz,nx);
    k(1:index_max_wavenumber+1,1)=(0:index_max_wavenumber)';%求0：N/2的波数
    for j=index_max_wavenumber+2:nz
       k(j,1)=-(nz-(j-1));%N/2+1:N-1的波数
    end
    for ix=2:nx
        k(:,ix)=k(:,1);
    end
    kz=wz*k;
end