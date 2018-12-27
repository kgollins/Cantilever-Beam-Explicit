% Compute the Jacobian Matrix [J], the strain interpolation maxtrix [B],
% and the Mass Matrix [M]
%
% Kenneth Gollins
% Assignment 7
% 12/2/13

function [J,B,M,Ke] = jacB(x,y,rho,t,E,v)
syms r s

N1 = .25*(1-r)*(1-s)*(-r-s-1);
N2 = .5*(1-s)*(1+r)*(1-r);
N3 = .25*(1+r)*(1-s)*(r-s-1);
N4 = .5*(1+r)*(1+s)*(1-s);
N5 = .25*(1+r)*(1+s)*(r+s-1);
N6 = .5*(1+s)*(1+r)*(1-r);
N7 = .25*(1-r)*(1+s)*(-r+s-1);
N8 = .5*(1-r)*(1+s)*(1-s);

N=[N1 N2 N3 N4 N5 N6 N7 N8];

X=N*x';
Y=N*y';

dNxr = diff(X,r);
dNyr = diff(Y,r);
dNxs = diff(X,s);
dNys = diff(Y,s);

Jsym=[dNxr dNyr; dNxs dNys];
Js=inv(Jsym);
J=matlabFunction(Jsym);

%Compute the Mass Matrix
Nr = [N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8 0;...
    0 N1 0 N2 0 N3 0 N4 0 N5 0 N6 0 N7 0 N8];
inside = int(Nr'*Nr*det(Jsym),r,-1,1);
M = rho*t *int(inside,s,-1,1);

dN1r = diff(N1,r);
dN1s = diff(N1,s);
dN2r = diff(N2,r);
dN2s = diff(N2,s);
dN3r = diff(N3,r);
dN3s = diff(N3,s);
dN4r = diff(N4,r);
dN4s = diff(N4,s);
dN5r = diff(N5,r);
dN5s = diff(N5,s);
dN6r = diff(N6,r);
dN6s = diff(N6,s);
dN7r = diff(N7,r);
dN7s = diff(N7,s);
dN8r = diff(N8,r);
dN8s = diff(N8,s);


B1= [Js(1,1)*dN1r+Js(1,2)*dN1s 0 Js(1,1)*dN2r+Js(1,2)*dN2s 0 Js(1,1)*dN3r+Js(1,2)*dN3s 0 Js(1,1)*dN4r+Js(1,2)*dN4s 0,...
    Js(1,1)*dN5r+Js(1,2)*dN5s 0 Js(1,1)*dN6r+Js(1,2)*dN6s 0 Js(1,1)*dN7r+Js(1,2)*dN7s 0 Js(1,1)*dN8r+Js(1,2)*dN8s 0];

B2= [0 Js(2,1)*dN1r+Js(2,2)*dN1s 0 Js(2,1)*dN2r+Js(2,2)*dN2s 0 Js(2,1)*dN3r+Js(2,2)*dN3s 0 Js(2,1)*dN4r+Js(2,2)*dN4s,...
    0 Js(2,1)*dN5r+Js(2,2)*dN5s 0 Js(2,1)*dN6r+Js(2,2)*dN6s 0 Js(2,1)*dN7r+Js(2,2)*dN7s 0 Js(2,1)*dN8r+Js(2,2)*dN8s];

B3= [Js(2,1)*dN1r+Js(2,2)*dN1s Js(1,1)*dN1r+Js(1,2)*dN1s Js(2,1)*dN2r+Js(2,2)*dN2s Js(1,1)*dN2r+Js(1,2)*dN2s Js(2,1)*dN3r+Js(2,2)*dN3s Js(1,1)*dN3r+Js(1,2)*dN3s Js(2,1)*dN4r+Js(2,2)*dN4s Js(1,1)*dN4r+Js(1,2)*dN4s,...
    Js(2,1)*dN5r+Js(2,2)*dN5s Js(1,1)*dN5r+Js(1,2)*dN5s Js(2,1)*dN6r+Js(2,2)*dN6s Js(1,1)*dN6r+Js(1,2)*dN6s Js(2,1)*dN7r+Js(2,2)*dN7s Js(1,1)*dN7r+Js(1,2)*dN7s Js(2,1)*dN8r+Js(2,2)*dN8s Js(1,1)*dN8r+Js(1,2)*dN8s];

Bsym=[B1;B2;B3];
B=matlabFunction(Bsym);

D=(E*t)/(1-v^2)*[1 v 0; v 1 0; 0 0 (1-v)/2];
% Integral
    inside = int(Bsym'*D*Bsym*det(Jsym),r,-1,1);
    Ke=int(inside,s,-1,1);

% pretty(Jsym)
% pretty(Bsym)
% 
% latex(Bsym(3,1))


end