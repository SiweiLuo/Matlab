% laplace1d.m
%
% Solve:
%
%   -u''(x) +10*u(x) = f(x) = (4*pi^2+10)sin(2*pi*x),  x_min < x < x_max
%   u(x_min) = u(x_max) = 0.
%
% DPN 4/12/16

clear all; close all;

% For RunNumber=2 see Johnson 1.19 to explain machine precision accuracy
% (up to conditioning of A).
RunNumber = 1;

x_min = 0.0;
x_max = 1.0;
M = 20;

x = linspace(x_min,x_max,M)';

% Nodes (M of them)

Z = zeros(1,M);
Z(1,:) = x;

% Triangles (N of them)

N = M-1;
T = zeros(2,N);
T(1,:) = [1:M-1];
T(2,:) = [2:M];

% Right-hand side for the exact solution u = sin(2*pi*x)

if(RunNumber==1)
  f = ((2*pi)^2+10.)*sin(2*pi*x);
elseif(RunNumber==2)
  f = 2*ones(M,1);
end

% Assembly

A = zeros(M,M);
b = zeros(M,1);

for n=1:N
  x_1 = Z(1,T(1,n));
  x_2 = Z(1,T(2,n));
  f_1 = f(T(1,n));
  f_2 = f(T(2,n));
  h = x_2 - x_1;
  a_n = zeros(2,2);
  a_n(1,1) = 1.0/h;
  a_n(1,2) = -1.0/h;
  a_n(2,1) = -1.0/h;
  a_n(2,2) = 1.0/h;
  q_n = zeros(2,2);
  q_n(1,1) = h/0.3;
  q_n(1,2) = h/0.6;
  q_n(2,1) = h/0.6;
  q_n(2,2) = h/0.3;
  b_n = zeros(2,1);
  b_n(1) = f_1*h/2.0;
  b_n(2) = f_2*h/2.0;
  for alpha=1:2
    for beta=1:2   
      A(T(alpha,n),T(beta,n)) = A(T(alpha,n),T(beta,n))+ a_n(alpha,beta)+q_n(alpha,beta);
    end
    b(T(alpha,n)) = b(T(alpha,n)) + b_n(alpha);
  end
end

A = A(2:(M-1),2:(M-1));
b = b(2:(M-1));

x_plot = x(2:(M-1));
if(RunNumber==1)
  u_exact = sin(2*pi*x_plot);
elseif(RunNumber==2)
  u_exact = (x_plot-x_min).*(x_max-x_plot);
end

% Solution

xi=A\b;

% Plot the solution

figure(1);
hh = gca;

plot(x_plot,xi,'g-*',x_plot,u_exact,'b-o');
xlabel('$x$','interpreter','latex');
ylabel('$y$','interpreter','latex');
title('Plot of Approximation and Exact Solution','interpreter','latex');
legend('Approximation','Exact');
set(hh,'fontsize',16);

%saveas(hh,'laplace1d','epsc');

fprintf('|xi-u_exact| = %g\n',norm(xi-u_exact,inf));