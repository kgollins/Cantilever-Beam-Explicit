%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Dynamic Impact on a Cantelever Beam
%
% Kenneth Gollins
% Finite Element Methods
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
close all
clc

%material properties
E = 70e9;%GPa
v = 0.25;
t = .01; %m
rho = 2500; %kg/m^3
z1 = .01; %D.C. for mode 1
z2 = .02; %D.C. for mode 2

%coordinates
x=[-5 0 5 5 5 0 -5 -5]*.01;
y=[-.5 -.5 -.5 0 .5 .5 .5 0]*.01;

%-------------------Compute J,B and M------------------------------------
% Compute the Jacobian Matrix [J], the strain interpolation maxtrix [B]
% and the Mass Matrix [M]
% Input: X and Y coordinates
% Output: J,B,M

[J,B,M,Ke] = jacB(x,y,rho,t,E,v);
    %Ke = stiffness(J,B,E,v,t) %Gauss 3x3 (weird results???)
%-------------------Apply BC---------------------------------------------
% U1=V1=U7=V7=U8=V8=0
Kg_r=double(Ke(3:12,3:12)); %reduced Stiffness Matrix
Mg_r=double(M(3:12,3:12)); %reduced Mass Matrix

%------------------Eigenvalue problem-----------------------------------
[phi2,omega] = eig(Kg_r,Mg_r);
phi = phi2(:,1:2);%reduced phi for the first 2 modes

%------------------Solve for C------------------------------------------

alpha=[1 omega(1,1);1 omega(2,2)];
beta=[2*sqrt(omega(1,1))*z1;2*sqrt(omega(2,2))*z2];

X=inv(alpha)*beta;

Cg_r=X(1)*Mg_r + X(2)*Kg_r;

%-------------------Diagonalize-----------------------------------------
Ks=phi'*Kg_r*phi;
Ms=phi'*Mg_r*phi;
Cs=phi'*Cg_r*phi;

%-------------------solve for displacement using CDM-------------------
    [v,time] = cdm(Ks,Ms,Cs,phi);
    u1=zeros(1,length(v));

for x=1:length(v)
    v1(x)=v{x}(1);%mode 1 amplitudes
end
for x=1:length(v)
    v2(x)=v{x}(2);%mode 2 amplitudes
end
    vv=[v1' v2'];

    disp = phi*vv';
    max(abs(disp(8,:)))

%--------------Plot Dispacement at node 5---------------
figure(1)
plot(time*1e3,disp(8,1:end-1)*1e3)
xlabel('Time (ms)')
ylabel('Displacement (mm)')

%--------------Plot Mode Shapes-------------------------

A=t^2;
m=rho*A;
I=t*t^3/12;
w1=sqrt(omega(1,1))/(2*pi); %hz
w2=sqrt(omega(2,2))/(2*pi); %hz  

% % Mode 1
% %     gamma1=sqrt((w1/(sqrt(E*I/(m*(0.1^4))))));
% %     i=1e-3;
% %     x=0;
% %     for n=1:100
% %         y=(sin(gamma1)-sinh(gamma1))*(sin(x)-sinh(x))+ (cos(gamma1)+cosh(gamma1))...
% %             *(cos(x)-cosh(x));
% %         x=x+(gamma1/0.1)*i;
% %         mode_shape(n)=y;
% %     end
% %     figure(2)
% %     plot(length,mode_shape/max(abs(mode_shape)));
% %     ylabel('Mode Amplitude (dimensionless)')
% %     xlabel('Beam Length (m)')
% % 
% % Mode 2
% %     gamma1=sqrt((w2/(sqrt(E*I/(m*(0.1^4))))));
% %     i=1e-3;
% %     x=0;
% %     for n=1:100
% %         y=(sin(gamma1)-sinh(gamma1))*(sin(x)-sinh(x))+ (cos(gamma1)+cosh(gamma1))...
% %             *(cos(x)-cosh(x));
% %         x=x+(gamma1/0.1)*i;
% %         mode_shape(n)=y;
% %     end
% %     figure(3)
% %     plot(linspace(0,.1),mode_shape/max(abs(mode_shape)));
% %     ylabel('Mode Amplitude (dimensionless)')
% %     xlabel('Beam Length (m)')


















