%Central Difference Method

function [v,t] = cdm(k,m,c,phi)

dt=1e-6;
time=100e-4;
step=time/dt;

v{1}=zeros(length(k),1);%t-1
v{2}=zeros(length(k),1);%t
M_hat = (m/dt^2)+(c/(2*dt));

for n=2:step
t(n)=n*dt;

    if n<=100
        Fr=[0;0;0;0;0;0;0;800*t(n);0;0];
        Fs=phi'*-Fr;
    elseif n>=101 && n<=200
        Fr=[0;0;0;0;0;0;0;(800-800*(t(n)-1e-4));0;0];
        Fs=phi'*-Fr;
    else
        Fr=[0;0;0;0;0;0;0;0;0;0];
        Fs=phi'*Fr;
    end
        R_hat = Fs - (k-(2*m)/dt^2)*v{n}- ((m/dt^2)-(c/(2*dt)))*v{n-1};

    v{n+1}=M_hat\R_hat;

end

end