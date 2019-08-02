function A = udsf(N)
A=zeros(6,6);
h=1/N;
cp=zeros(1,N-1);
for i=1:(N-1)
    for j=1:(N-1)
        cp(i)=cp(i)+sin(i*j*pi*h)*sin(4*pi*j*h);
    end
    cp(i)=cp(i)*h*sqrt(2);
end

for i=1:6
    for j=1:6
        for p=1:(N-1)
            A(j,i)=A(j,i)+cp(p)*sin(p*pi*(i-1)/5)*exp(-pi*pi*p^2*(j-1)/50);
        end
        A(j,i)=A(j,i)*sqrt(2);
    end
end
end

