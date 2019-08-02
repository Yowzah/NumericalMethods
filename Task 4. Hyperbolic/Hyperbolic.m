function Hyperbolic %вариант 12 (гиперболический тип)
    clc
    syms x t fun;
    p(x) = x+2 ;
    u(x,t) = x + t^2;
    q = diff(u,x);
    fun = symfun(diff(u,t,2) - diff(p(x)*q(x,t),x),[t,x]);
    phi(x) = subs(u,t,0);
    psi(x) = subs(diff(u,t),t,0);
    a(t) = subs(u,x,0);
    b(t) = subs(u + diff(u,x),x,1);
    
    N = 5;
    J = zeros(N+1,3);
    T = 0:0.2:1;
    X = 0:0.2:1;
    Y = formula(subs(subs(u,t,T),x,X'));
    V = zeros(6,6);
    H = zeros(3,1);
    Q = zeros(3,1);
    G = zeros(3,1);
    l = 1;
    while N <= 20
        disp('Step');
        M = 2*N;
        h = 1/N;
        H(l,1) = h;
        r = 1/M;
        X = 0:h:1;
        T = 0:r:1;
        W = formula(subs(subs(u,t,T),x,X')); %page 53 whr
        disp('Matrix Whr:');
        disp(vpa(W(1:0.2*N:end,1:0.2*M:end),6));
        
        U = zeros(N+1,M+1);
        L = zeros(N,M);
        U(:,1) = phi(X)';
        U(:,2) = U(:,1) + r*psi(X)';
        
        for k = 2:M
            for i = 2:N
                L(i,k) = p(X(i)+h/2)*(U(i+1,k)-U(i,k))/(h^2) - ...
                    p(X(i)-h/2)*(U(i,k)-U(i-1,k))/(h^2);
                U(i,k+1) = 2*U(i,k)-U(i,k-1)+(r^2)*(L(i,k)+fun(X(i),T(k)));
            end
            U(1,k+1) = a(T(k+1));
            U(N+1,k+1) = (b(T(k+1))+(4*U(N,k+1)-U(N-1,k+1))/(2*h))/(1+3/(2*h));
        end 
        disp('U matrix after aproximation in points: ');
        disp(U(1:0.2*N:end,1:0.2*M:end));
        fprintf('h = %.6f\n', h);
        
        Q(l,1) = max(max(Y-U(1:0.2*N:end,1:0.2*M:end)));
        fprintf('discrepancy = %.6f\n', vpa(max(max(Y-U(1:0.2*N:end,1:0.2*M:end))),6));

        if h ~= 0.2
            G(l,1) = max(max(V-U(1:0.2*N:end,1:0.2*M:end)));
            %disp(vpa(max(max(V-U(1:0.2*N:end,1:0.2*M:end))),6));
        end
        V = U(1:0.2*N:end,1:0.2*M:end);
        N = 2*N;
        l = l+1;
        fprintf('\n');
    end
    disp('Aproximation ended');
    Tabl = table(H(:,1),Q(:,1),G(:,1));
    disp(Tabl);
end

