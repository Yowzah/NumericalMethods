    function FredholmKer
    syms x y H h f R a b u;
    f(x) = x - 0.5;
    
    a(x) = [x/2, -x^3/6, x^5/15, -17*x^7/630];
    b(x) = [x x^3 x^5 x^7];

    G = eye(4) - int(b'*a, 0, 1); % A from (6)
    F = int(b'*f, 0, 1); % B from (6)
    
    form_a = formula(a);
    form_b = formula(b);
    
    C = G(1:3,1:3)\F(1:3); % C = A/B
    
    u(x) = f(x) + form_a(1:3)*C;
    
    %vpa(expand(u),6);
    start_1 = vpa(u(0));
    middle_1 = vpa(u(0.5));
    end_1 = vpa(u(1));
    
    
    

    R(x,y) = form_a(1:3)*inv(G(1:3,1:3))*subs(form_b(1:3).',x,y); % G^~
    nu = 17/630;
    B = R(1,1);
    
    mark = vpa((1+B)*nu/(1-(1+B)*nu)*subs(abs(u), fminbnd(-abs(u),0,1)));
    
    C = G\F;
    u(x) = f(x) + a*C;
    %vpa(expand(u),6)
    start_2 = vpa(u(0));
    middle_2 = vpa(u(0.5));
    end_2 = vpa(u(1));
    
    X = {'u^3'; 'u^4'};
    U_3 = [start_1;middle_1;end_1];
    U_4 = [start_2; middle_2;end_2];
    
    disp(U_3)
    disp(U_4)
    disp(U_3-U_4)

    disp(mark)
    
    FredholmKF(3)
end

