%clear all;
format long;
syms x;

uf = u_Fourier();

udsf_5 = u_DFourier(5);
udsf_10 = u_DFourier(10);
udsf_20 = u_DFourier(20);

uf

udsf_5
udsf_10
udsf_20

sigma = 0;
uimpl5_5 = explicit_scheme2(5,5);
uimpl10_20 = explicit_scheme2(10,20);
uimpl20_80 = explicit_scheme2(20,80);
uimpl20_20 = explicit_scheme2(20,20);

error(1,:) = [max(max(abs(uf-uimpl5_5))), max(max(abs(uf-uimpl10_20))), max(max(abs(uf-uimpl20_80))), max(max(abs(uf-uimpl20_20)))];

sigma = 1;
uimpl5_5 = implicit_scheme2(5,5, sigma);
uimpl10_20 = implicit_scheme2(10,20,sigma);
uimpl20_80 = implicit_scheme2(20,80,sigma);
uimpl20_20 = implicit_scheme2(20,20,sigma);

error(2,:) = [max(max(abs(uf-uimpl5_5))), max(max(abs(uf-uimpl10_20))), max(max(abs(uf-uimpl20_80))), max(max(abs(uf-uimpl20_20)))];

sigma = 0.5;
uimpl5_5 = implicit_scheme2(5,5, sigma);
uimpl10_20 = implicit_scheme2(10,20,sigma);
uimpl20_80 = implicit_scheme2(20,80,sigma);
uimpl20_20 = implicit_scheme2(20,20,sigma);

error(3,:) = [max(max(abs(uf-uimpl5_5))), max(max(abs(uf-uimpl10_20))), max(max(abs(uf-uimpl20_80))), max(max(abs(uf-uimpl20_20)))];

sigma = 0.5 - 0.2^2/(12*0.02);
uimpl5_5 = implicit_scheme2(5,5, sigma);
sigma = 0.5 - 0.1^2/(12*0.005);
uimpl10_20 = implicit_scheme2(10,20,sigma);
sigma = 0.5 - 0.05^2/(12*0.00125);
uimpl20_80 = implicit_scheme2(20,80,sigma);
sigma = 0.5 - 0.05^2/(12*0.005);
uimpl20_20 = implicit_scheme2(20,20,sigma);
error(4,:) = [max(max(abs(uf-uimpl5_5))), max(max(abs(uf-uimpl10_20))), max(max(abs(uf-uimpl20_80))), max(max(abs(uf-uimpl20_20)))];

double(uf)

temp(1) = max(max(abs(uf-udsf_5)));
temp(2) = max(max(abs(uf-udsf_10)));
temp(3) = max(max(abs(uf-udsf_20)));
double(temp)'

text = [0; 1; 0.5; 42];
cols3 = [text, double(error)];
cnames = {'sigma', '(0.2,0.02)', '(0.1,0.005)','(0.05,0.00125)', '(0.05,0.005'};
uitable('Parent', figure('Name', 'Table2', 'Position', [500 250 500 150]), 'Position', [50 50 380 100], 'Data', cols3, 'ColumnName', cnames, 'RowName', ([]));