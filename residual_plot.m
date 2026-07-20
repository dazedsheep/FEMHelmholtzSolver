lin = load("residue_linear_case.mat");
nonlin = load("residue.mat");

col = 10;
iters = 1:20;
% J
figure, plot(iters, nonlin.residue_base_case(iters,col),'-*');
hold on
plot(iters, (lin.residue(iters,col)),'-.');
plot(iters, 0.6.^(iters));
xlabel('Iteration (n)')
ylabel('J_n(x,x_n)');
yscale log
ylim([10^-5,5*10^(-3)])
legend("nonlinear base case", "linear base case");
set(gca,'fontname','Arial');  % Set it to times
set(gca, 'FontWeight', 'bold');

hold off

%%
% boundary residual
res_nonlinear = sqrt(sum(nonlin.residue_base_case(iters,1:3),2));
res_linear = sqrt(sum(lin.residue(iters,1:3),2));

figure, plot(iters, (res_nonlinear),'-*');
hold on
plot(iters, (res_linear),'-.');
xlabel('Iteration (n)')
ylabel('||F(x_n) - h||_L^2(\partial \Omega)');
%ylim([0,5*10^(-3)])
legend("nonlinear base case", "linear base case");
set(gca,'fontname','Arial');  % Set it to times
set(gca, 'FontWeight', 'bold');
yscale log
