function u = explicit_scheme2(N, M)
syms x;
syms t;
h = 1/N;
tau = 0.1/M;

%if tau/h^2<0.5
%     disp('Условие устойчивости выполнено.');
%else
%     disp('Условие устойчивости НЕ выполнено.');
%end

for p=1:N-1
     for i=1:N-1
          temp(p,i) = sin(p*pi*i*h);
     end
     phis(p) = fun_phi(p*h);
end



k=0;
u(k+1,0+1) = 0;
u(k+1,5+1) = 0; 
for i=1:4
     for p=1:N-1
          C(p) = h*sqrt(2)*(phis*temp(p,:)');
     end;
     j = fix(1/5*i/h);
     u(k+1,i+1) = sqrt(2)*(C*temp(:,j));
end;

for k=1:5
     u(k+1,0+1) = 0;
     u(k+1,5+1) = 0;
     for i=1:4
          for p=1:N-1
               C(p) = h*sqrt(2)*(phis*temp(p,:)');
               lambda(p) = 1-4*tau/(h^2)*(sin(p*pi*h/2))^2;
          end;
          j = fix(1/5*i/h);
          
          n = fix(0.1/5*k/tau);
          u(k+1,i+1) = sqrt(2)*(C.*(lambda.^n)*temp(:,j));
     end;
end;

end