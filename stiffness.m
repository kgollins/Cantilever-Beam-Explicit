% Kenneth Gollins

function Ke = stiffness(J,B,E,v,t)

D=(E*t)/(1-v^2)*[1 v 0; v 1 0; 0 0 (1-v)/2];

%Gauss 3x3
Ke=0;
w=[8/9 5/9 5/9];
r=[0 sqrt(3/5) -sqrt(3/5)];
s=[0 sqrt(3/5) -sqrt(3/5)];

for i=1:3
    for j=1:3
        Ke=Ke + w(i)*w(j)*B(r(i),s(j))'*D*B(r(i),s(j))*det(J(r(i),s(j)));
    end
end

end







