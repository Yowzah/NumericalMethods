function [G_l] = Ritz()
    syms W p q x;
    k = 1.57894;
    l = 8.59453;
    p(x) = k*x + l;
    q(x) = k^2*(1/(k*x + l) - k*x);
    fprintf("Метод Ритца\n");
    n = 7;
    for i = 0:(n-1) %"выбор координатной системы"
        d = 1/4*sqrt((i+3)*(i+4)*(2*i+5)/(2*(i+1)*(i+2)));
        W(i+1,1) = (1 - x^2)*d*jacobiP(i,2,2,x);
    end

    G = zeros(n);
    G_l = zeros(n);
    for i = 1:n
        for j = 1:i
            G(i,j) = vpaintegral(W(i,1)*W(j,1), [-1 1]);
            G(j,i) = G(i,j);
            G_l(i,j) = vpaintegral(p*diff(W(i,1))*diff(W(j,1)) + q*W(i,1)*W(j,1), [-1 1]);
            G_l(j,i) = G_l(i,j);
        end
    end
    e = sort(eig(G_l));
    disp(vpa(e(1)));
    disp(vpa(e(2)));
end

