function out=derivate1_fd8(in,direction,interval)
% the 8th order accuracy finite difference to calculate the first derivative
% input parameters:
% 		in: input 2D array
% 		out: out 2D array
% 		direction: the direction for calculating directive, =2: calculate derivative along x (horizontal);
%                   =1 : calcualte derivative along z (vertical)
% 		interval: the interval of the calculation direction
% return:
%       out: out 2D array
% 		¡¾NB: make sure the input pointer, in and out, are independent, otherwise the vectorization may produce wrong result¡¿
% 		¡¾NB: this function has been tested regorously¡¿

% 	//float C84={0.8,-0.20,0.038095238095238,-0.003571428571429};//classic coefficient
% 	//float C42={0.666666666666667,-0.083333333333333};//classic coefficient
% 	//float C84={0.839725908387582,-0.243216342043984,0.059588302317624,-0.008074531266178};//my optimized coefficient
% 	//float C42={0.678723260800114,-0.089577669496057};//my optimized coefficient
inv_interval1=1/interval;
inv_interval2=1/(2*interval);

C80=0.839725908387582*inv_interval1;
C81=-0.243216342043984*inv_interval1;
C82=0.059588302317624*inv_interval1;
C83=-0.008074531266178*inv_interval1;

C40=0.678723260800114*inv_interval1;
C41=-0.089577669496057*inv_interval1;

[nz,nx]=size(in);
out=zeros(nz,nx);
if direction==1
    for ix=1:nx
        out(1,ix)=(in(2,ix)-in(1,ix))*inv_interval1;
        out(2,ix)=(in(3,ix)-in(1,ix))*inv_interval2;
        out(3,ix)=C40*(in(4,ix)-in(2,ix))+C41*(in(5,ix)-in(1,ix));
        out(4,ix)=C40*(in(5,ix)-in(3,ix))+C41*(in(6,ix)-in(2,ix));
        
        for iz=5:nz-4
            out(iz,ix)=C80*(in(iz+1,ix)-in(iz-1,ix))+C81*(in(iz+2,ix)-in(iz-2,ix))...
                + C82*(in(iz+3,ix)-in(iz-3,ix))+C83*(in(iz+4,ix)-in(iz-4,ix));
        end
        
        out(nz-3,ix)=C40*(in(nz-2,ix)-in(nz-4,ix))+C41*(in(nz-1,ix)-in(nz-5,ix));
        out(nz-2,ix)=C40*(in(nz-1,ix)-in(nz-3,ix))+C41*(in(nz,ix)-in(nz-4,ix));
        
        out(nz-1,ix)=(in(nz,ix)-in(nz-2,ix))*inv_interval2;
        out(nz,ix)=(in(nz,ix)-in(nz-1,ix))*inv_interval1;
    end
end

if direction==2
    for iz=1:nz
        out(iz,1)=(in(iz,2)-in(iz,1))*inv_interval1;
        out(iz,2)=(in(iz,3)-in(iz,1))*inv_interval2;
        out(iz,3)=C40*(in(iz,4)-in(iz,2))+C41*(in(iz,5)-in(iz,1));
        out(iz,4)=C40*(in(iz,5)-in(iz,3))+C41*(in(iz,6)-in(iz,2));
        
        for ix=5:nx-4
            out(iz,ix)=C80*(in(iz,ix+1)-in(iz,ix-1))+C81*(in(iz,ix+2)-in(iz,ix-2))...
                + C82*(in(iz,ix+3)-in(iz,ix-3))+C83*(in(iz,ix+4)-in(iz,ix-4));
        end
        
        out(iz,nx-3)=C40*(in(iz,nx-2)-in(iz,nx-4))+C41*(in(iz,nx-1)-in(iz,nx-5));
        out(iz,nx-2)=C40*(in(iz,nx-1)-in(iz,nx-3))+C41*(in(iz,nx)-in(iz,nx-4));
        
        out(iz,nx-1)=(in(iz,nx)-in(iz,nx-2))*inv_interval2;
        out(iz,nx)=(in(iz,nx)-in(iz,nx-1))*inv_interval1;
    end    
end
end