function [] = scalarMethodComposition(G_l)
    disp('Метод скалярных произведений');
    epsilon = input('epsilon: ');
    y = ones(size(G_l, 1), 1);
    z = G_l\y;
    lambda0 = 0;
    lambda1 = (z'*y)/(y'*y);
    while (abs(lambda1 - lambda0) > epsilon)
        y = z;
        [~,I] = max(abs(y));
        y = y / y(I,1);
        z = inv(G_l)*y;
        lambda0 = lambda1;
        lambda1 = (z'*y)/(y'*y);
        disp(vpa(1/lambda1));
    end
    lambda1 = 1/lambda1;
    disp(vpa(lambda1));
end

