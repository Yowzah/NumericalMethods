function A = solweight(weight,N,M)
A=zeros(6,6);
h=1/N;
h2=1/5;
tau=0.1/M;
cp=zeros(1,N-1);
for i=1:(N-1)
    for j=1:(N-1)
        cp(i)=cp(i)+sin(i*j*pi*h)*sin(4*pi*j*h); %20
    end
    cp(i)=cp(i)*h*sqrt(2);
end
cp(3);
for i=1:6
    A(1,i)=sin(4*pi*h2*(i-1));
end
lyambda3=(1-4*(1-weight)*tau*sin(4*pi*h/2)*sin(4*pi*h/2)/(h^2))/(1+3*tau*weight*sin(4*pi*h/2)*sin(4*pi*h/2)/(h^2)); %21
for i=2:6
    for j=2:5
        A(i,j)=sqrt(2)*cp(3)*sin(4*pi*(j-1)*h2)*((lyambda3)^((i-1)*M/5)); %21.5
    end
end
end

