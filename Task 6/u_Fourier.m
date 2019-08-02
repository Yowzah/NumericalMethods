function U = u_Fourier()
format long;
syms x;
syms t;
p=2;
     C = 1/sqrt(2);
     T = C*exp(-(pi*p)^2*t);
     Psis = sin(pi*p*x);

u = sqrt(2)*T*Psis;

     xx = 0:0.2:1;
     tt = 0:0.1/5:0.1;
     for i=1:6
          for j=1:6
               U(i,j) = vpa(subs(subs(u,t,tt(i)),x,xx(j)));
          end;
     end;

end