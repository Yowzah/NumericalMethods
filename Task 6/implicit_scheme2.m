function u = implicit_scheme2(N,M,sigma)
syms x;
syms t;
h = 1/N;
tau = 0.1/M;

%if sigma>=0 & sigma<0.5
%     if tau<=h^2/(2*(1-2*sigma))
%          disp('Условие устойчивости выполнено.');
%     else
%          disp('Условие устойчивости НЕ выполнено.');
%     end;
%end;
%if sigma>=0.5 & sigma<=1
%     disp('Условие устойчивости выполнено.');
%end;


for p=1:N-1
     for i=1:N-1
          temp(p,i) = sin(p*pi*i*h);
     end
     phis(p) = fun_phi(p*h);
end

for k=0:5
     u(k+1,0+1) = 0;
     u(k+1,5+1) = 0;
     for i=1:4
          for p=1:N-1
               C(p) = h*sqrt(2)*(phis*temp(p,:)');
               lambda(p) = (1-4*tau*(1-sigma)/(h^2)*(sin(p*pi*h/2))^2)/(1+4*sigma*tau/(h^2)*(sin(p*pi*h/2))^2);
          end;
          j = fix(1/5*i/h);
          n = fix(0.1/5*k/tau);
     u(k+1,i+1) = sqrt(2)*(C.*(lambda.^n)*temp(:,j));
     end;
end;

end