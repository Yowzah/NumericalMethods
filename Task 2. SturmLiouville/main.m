format long;
clc;
fprintf("Проблема собственных значений в задаче Штурма-Лиувилля\n");
fprintf("Вариант 7\n");
syms x W y_1 y_2 phi p q;
k = 1.57894;
l = 8.59453;
disp(k);
disp(l);
fprintf(" u(-1) = u(1) = 0\n");

approximation();

G_l = Ritz();

scalarMethodComposition(G_l);

