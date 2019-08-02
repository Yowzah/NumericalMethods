function FredholmKF(n)
    syms x y h u f;
    
    h(x,y) = 0.5*tanh(x*y);
    
    f(x) = x - 0.5;
    
    X = 0:1/(2*n):1;
    A = zeros(2*n+1,1);
    
    for i = 2:2:2*n
        A(i,1) = 4;
        A(i+1,1) = 2;
    end
    A(1,1)=1;
    A(2*n+1,1)=1;
    A = A/(6*n);
    
    G = vpa(eye(2*n+1) - reshape(subs(A'.*subs(formula(h),y,X),x,X),[2*n+1,2*n+1])); % (16)
    disp(G);
    F = f(X');
    disp(F);
    z = linsolve(G,F);
    u(x) = subs(h,y,X)*(A.*z) + f(x);
    a = u(0);
    b = u(1/2);
    c = u(1);
    
    disp(vpa([a; b; c], 6));
end

