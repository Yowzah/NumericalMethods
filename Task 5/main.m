% phi=sin(4pi*x)
% c_3=1/sqrt2, psi_3=sin(4pi*x)
% u(x,t)=exp(-9*pi^2*t)*sin(4pi*x)
syms x t;
%clc
u(x,t)=exp(-9*pi^2*t)*sin(4*pi*x);
Uex=zeros(6,6);
U5=zeros(6,6);
U10=zeros(6,6);
U20=zeros(6,6);

for i=0:0.2:1
    for j=0:0.02:0.1
        Uex(j*50+1, i*5+1)=u(i,j);
    end;
end

fprintf("uf: \n");
Uex
fprintf("udsf in N=5: \n");
U5=udsf(5)
fprintf("udsf in N=10: \n");
U10=udsf(10)
fprintf("udsf in N=20: \n");
U20=udsf(20)
norm(solweight(0,10,20)-Uex, inf);
C=zeros(4,4);
solweight(0,5,5);
solweight(0,5,5)-Uex;
C(1,1)=max(max(abs(solweight(0,5,5)-Uex)));
C(2,1)=max(max(abs(solweight(1,5,5)-Uex)));
C(3,1)=max(max(abs(solweight(1/2,5,5)-Uex)));
C(4,1)=max(max(abs(solweight(1/2-0.04/(12*0.02),5,5)-Uex)));

C(1,2)=max(max(abs(solweight(0,10,20)-Uex)));
C(2,2)=max(max(abs(solweight(1,10,20)-Uex)));
C(3,2)=max(max(abs(solweight(1/2,10,20)-Uex)));
C(4,2)=max(max(abs(solweight(1/2-0.1*0.1/(12*0.005),10,20)-Uex)));

C(1,3)=max(max(abs(solweight(0,20,80)-Uex)));
C(2,3)=max(max(abs(solweight(1,20,80)-Uex)));
C(3,3)=max(max(abs(solweight(1/2,20,80)-Uex)));
C(4,3)=max(max(abs(solweight(1/2-0.05*0.05/(12*0.00125),20,80)-Uex)));

C(1,4)=max(max(abs(solweight(0,20,20)-Uex)));
C(2,4)=max(max(abs(solweight(1,20,20)-Uex)));
C(3,4)=max(max(abs(solweight(1/2,20,20)-Uex)));
C(4,4)=max(max(abs(solweight(1/2-0.05*0.05/(12*0.005),20,20)-Uex)));

fprintf("sigma=0             ");
for i=1:4
    fprintf("%7.0e ", C(1,i));
end
fprintf("\n");
fprintf("sigma=1             ");
for i=1:4 
    fprintf("%7.0e ", C(2,i));
end
fprintf("\n");
fprintf("sigma=1/2           ");
for i=1:4
    fprintf("%7.0e ", C(3,i));
end
fprintf("\n");
fprintf("sigma=1/2-h^2/12tau ");
for i=1:4
    fprintf("%7.0e ", C(4,i));
end
fprintf("\n");


