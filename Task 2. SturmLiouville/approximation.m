function [] = approximation()
    syms x y_1 y_2 p q;
    k = 1.57894;
    l = 8.59453;
    p(x) = k*x + l;
    q(x) = k^2*(1/(k*x + l) - k*x);
    p_min = p(-1);
    p_max = p(1);   
    q_min = q(1);
    q_max = q(-1);
    
    disp("Получить верхние и нижние оценки для первых двух собственных чисел");

    y_1(x) = cos(pi/2 * x);
    disp('cos(pi/2 * x)');

    lymbda_min_1 = vpa((pi/2)^2 * p_min + q_min, 8);
    nu_1 = vpa(subs(-diff(p_min*diff(cos(pi/2 *x))) + q_min*cos(pi/2 *x) - ((pi/2)^2 * p_min + q_min)*cos(pi/2 *x), 0), 13); %невязка используя формулу 1
    lymbda_max_1 = vpa((pi/2)^2 * p_max + q_max, 6);
    nu_2 = vpa(subs(-diff(p_max*diff(cos(pi/2 *x))) + q_max*cos(pi/2 *x) - ((pi/2)^2 * p_max + q_max)*cos(pi/2 *x), 0), 13);

    y_2(x) = sin(pi*x);
    disp('sin(pi*x)');

    lymbda_min_2 = vpa(pi^2 * p_min + q_min, 6);
    nu_3 = vpa(subs(-diff(p_min*diff(sin(pi*x))) + q_min*sin(pi*x) - ((pi)^2 * p_min + q_min)*sin(pi*x), 0), 13);
    lymbda_max_2 = vpa(pi^2 * p_max + q_max, 6);
    nu_4 = vpa(subs(-diff(p_max*diff(sin(pi*x))) + q_max*sin(pi*x) - ((pi)^2 * p_max + q_max)*sin(pi*x), 0), 13);
    
    str = {'min';'max'};
    lymbda_1 = [double(lymbda_min_1); double(lymbda_max_1)];
    e_1 = [double(nu_2); double(nu_1)];
    e_2 = [double(nu_4); double(nu_3)];
    lymbda_2 = [double(lymbda_min_2); double(lymbda_max_2)];
    T = table(str, lymbda_1, e_1, lymbda_2, e_2);
    disp(T);
    
    disp('Получить приближения для первых двух собственных чисел:');
    int_1 = vpaintegral(p(x)*((diff(cos(pi/2 *x)))^2) + q(x)*(cos(pi/2 *x))^2, x, -1, 1);
    int_2 = vpaintegral((cos(pi/2 *x))^2, -1, 1);
    lymbda_1 = int_1/int_2;
    disp(lymbda_1);

    int_1 = vpaintegral(p(x)*((diff(sin(pi*x)))^2) + q(x)*(sin(pi*x))^2, x, -1, 1);
    int_2 = vpaintegral((sin(pi*x))^2, x, -1, 1);
    lymbda_2 = int_1/int_2;
    disp(lymbda_2);
end

