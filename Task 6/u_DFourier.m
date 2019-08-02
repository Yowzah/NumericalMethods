function U = u_DFourier(N)
syms x;
syms t;
h = 1/N;
xx = linspace(0,1, N+1);

for p=1:N-1
     for i=1:N-1
          Phii(i) = fun_phi(i*h);
          Psii(i) = sin(pi*p*i*h);
     end;
     C(p) = sqrt(2)*h*(Phii*Psii');
     
     temp(p) = exp(-(pi*p)^2*t)*sin(pi*p*x);
end;

u = sqrt(2)*(C*temp');

     xx = 0:0.2:1;
     tt = 0:0.1/5:0.1;
     for i=1:6
          for j=1:6
               U(i,j) = vpa(subs(subs(u,t,tt(i)),x,xx(j)));
          end;
     end;
end