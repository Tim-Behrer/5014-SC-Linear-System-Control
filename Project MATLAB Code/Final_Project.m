%% Header
% ASEN 5014, Linear Control Systems
% Final Project
% Author:
% Date Last Modified: 05/01/24 [mm/dd/yy]
% NOTE: 
%   - This final project is the analysis and design for the spacraft flying
%   formation.
%% Housekeeping 
clc;
clear all;
close all;
%% Simulation Stuff
% initial condition
x0 = [10;10;10;0.2;0.2;0.2];
% time stuff
t0 = 0;
dt = 0.5;
tf = 35*60;
t = (t0:dt:tf)';
%% Defining System
n = sqrt((398600/(6778^3))); % mean motion of the primary spacecraft [1/s]
% state-space system
%   - x = [x,y,z,dx,dy,dz]'
%   - y = [x,y,z]'
%   - u = [ux,uy,uz]'
A = [0,0,0,1,0,0;0,0,0,0,1,0;0,0,0,0,0,1;(3*(n^2)),0,0,0,(2*n),0;0,0,0,(-2*n),0,0;0,0,-(n^2),0,0,0];
B = [0,0,0;0,0,0;0,0,0;1,0,0;0,1,0;0,0,1];
C = [1,0,0,0,0,0;0,1,0,0,0,0;0,0,1,0,0,0];
D = zeros(3,3);
% characterizing open-loop dynamics
[e_vec_OL,e_val_OL] = eig(A);
dis_val = eig(A);
m = [2;2;2];
q = [size(null(A-(dis_val(1)*eye(size(A)))),2);size(null(A-(dis_val(3)*eye(size(A)))),2);size(null(A-(dis_val(4)*eye(size(A)))),2)];
fprintf('---|OPEN-LOOP BEHAVIOR|---\n');
fprintf('Open-Loop Eigenvalues:\n')
fprintf('   Eigenvalue 1: %0.2f | m = %d | q = %d\n',dis_val(1),m(1),q(1));
fprintf('   Eigenvalue 2: %0.2f%+0.5fi | m = %d | q = %d\n',real(dis_val(3)),imag(dis_val(3)),m(1),q(2));
fprintf('   Eigenvalue 3: %0.2f%+0.5fi | m = %d | q = %d\n',real(dis_val(4)),imag(dis_val(4)),m(1),q(3));
% simulating open loop response 
r = [0;0;0]; % desired refrence point 
u = r*ones(1,length(t)); % making our desired refrence a point in time
% simulating the system with only u=Fr
sys = ss(A,B,C,D);
y_ol = lsim(sys,u,t,x0);
% plot of initial system response 
figure
sgtitle('Open-Loop System Response'); 
subplot(1,3,1);
hold on;
grid minor;
plot(t,y_ol(:,1),'-b','LineWidth',1.2);
xlabel('Time [s]')
ylabel('x')
str = sprintf('x_0=%0.2f | dx_0/dt=%0.2f',x0(1),x0(4));
title(str)
subtitle('NO CONTROL INPUT');
subplot(1,3,2);
hold on;
grid minor;
plot(t,y_ol(:,2),'-b','LineWidth',1.2);
xlabel('Time [s]')
ylabel('y')
str = sprintf('y_0=%0.2f | dy_0/dt=%0.2f',x0(1),x0(4));
title(str)
subtitle('NO CONTROL INPUT');
subplot(1,3,3);
hold on;
grid minor;
plot(t,y_ol(:,3),'-b','LineWidth',1.2);
xlabel('Time [s]')
ylabel('z')
str = sprintf('z_0=%0.2f | dz_0/dt=%0.2f',x0(1),x0(4));
title(str)
subtitle('NO CONTROL INPUT');

%% System Controllability and Observability 
fprintf('---|CONTROLLABILITY AND OBSERVABILITY|---\n');
% controllability
Control = [B,A*B,(A^2)*B,(A^3)*B,(A^4)*B,(A^5)*B]; % controllabilty matrix
N = size(A,2);
control_check = rank(Control); % check condition
if control_check==N
    fprintf('Open-Loop system is controllable.\n');
    fprintf('   Rank(C) = %d = n\n',control_check);
else
    fprintf('Open-Loop system is not controllable.\n');
    fprintf('   Rank(C) = %d =/= n\n',control_check);
end
% oberservability
sigma = [C;C*A;C*(A^2);C*(A^3);C*(A^4);C*(A^5)]; % oberservability matrix
observe_check = rank(sigma); % check condition
if observe_check==N
    fprintf('Open-Loop system is observable.\n');
    fprintf('   Rank(sigma) = %d = n\n',control_check);
else
    fprintf('Open-Loop system is not observable.\n');
    fprintf('   Rank(sigma) = %d =/= n\n',control_check);
end
fprintf('System is both controllable and observable and is therefore minimal.\n')
%% Closed-Loop System controller 
fprintf('---|CLOSED-LOOP CONTROLLER|---\n');
% CONTROL LAW: U=-Kx+Fr
% Desired System Parameters
Mp = 0.05; % 5% max overshoot
ts = 15*60; % setteling time (15 min in seconds) 
xi = sqrt((log(Mp)^2)/((pi^2)+(log(Mp)^2))); % damping ratio
w_n = -log(0.05)/(xi*ts); % natural frequency
% finding our dominant 2nd order poles (states x,z)
desired_poles_1 = -(xi*w_n)+(w_n*sqrt(1-(xi^2))*1i);
desired_poles_2 = -(xi*w_n)-(w_n*sqrt(1-(xi^2))*1i);
scaling = 1.1; % other poles are 10% faster 
faster_pole_real = scaling*real(desired_poles_1);
% designing K 
eig_val_cl = [desired_poles_1;(faster_pole_real+(imag(desired_poles_2)*1i));desired_poles_2;(faster_pole_real+(imag(desired_poles_1)*1i));(faster_pole_real+(imag(desired_poles_2)*1i));(faster_pole_real+(imag(desired_poles_1)*1i))]; % desired closed-loop poles
K = place(A,B,eig_val_cl); % gain matrix based on desired poles
% designing F
F = (C*((-A+(B*K))^-1)*B)^-1; % 'feed-forward' matrix
% defining closed loop system
A_cl = A-(B*K);
B_cl = B*F;
C_cl = C;
D_cl = D;
r = [5;5;5]; % desired refrence point 
u = r*ones(1,length(t)); % making our desired refrence a point in time
u(:,1:floor(length(t)/10)) = zeros(size(u(:,1:floor(length(t)/10))));
% simulating the system with only u=Fr
F_mod = F; %= (C*((-A+(B*zeros(3,6)))^-1)*B)^-1;
B_cl_F = B*F_mod;
sys2 = ss(A,B_cl_F,C,D);
y_F = lsim(sys2,u,t,x0);
create_plots(y_F,u',t,3,x0,r,'System Response u=Fr');
% sumulating the system with full u=-Kx+Fr
sys3 = ss(A_cl,B_cl,C_cl,D_cl);
y_F = lsim(sys3,u,t,x0);
create_plots(y_F,u',t,3,x0,r,'System Response u=-Kx+Fx');
fprintf('Control Input u=Fr is NOT STABLE\n');
fprintf('Control Input u=-Kx+Fr is STABLE\n');
%% Desiging qnd Evaluating the Observer
des_freq = 10*-w_n;
q_obs = [des_freq,des_freq,des_freq-(1i*imag(desired_poles_1)),des_freq-(1i*imag(desired_poles_1)),des_freq+(1i*imag(desired_poles_1)),des_freq+(1i*imag(desired_poles_1))]; % desired poles
L = place(A',C',q_obs)'; 
A_obs = A-(L*C);
B_obs = zeros(size(B));
C_obs = eye(size(A_obs));
D_obs = zeros(size(B));
obs_sys = ss(A_obs,B_obs,C_obs,D_obs);
einit = [1;1;1;1;1;1];
[error,tOut] = initial(obs_sys,einit,400);
for ii=1:6
    check(ii) = find(abs(error(:,ii))<=0.05,1);
    t_con_e(ii) = tOut(check(ii));
end
% outputting to command line
fprintf('---|OBSERVER|---\n');
fprintf('Average time to achieve 5%% positional error %0.3f s \n',mean(t_con_e(1:3)));
fprintf('Average time to achieve 5%% derivative error %0.3f s \n',mean(t_con_e(4:6)));
% plot
figure
sgtitle('Observer Error vs. Time')
% positions
subplot(2,3,1);
hold on;
plot(tOut,error(:,1),'-b','LineWidth',1.2);
xline(tOut(check(1)),'-m','LineWidth',1);
grid on;
legend('','5% Error Time');
legend('Location','best');
xlabel('time (s)');
ylabel('$e_{x}$','interpreter','latex')
subplot(2,3,2);
hold on;
plot(tOut,error(:,2),'-r','LineWidth',1.2);
xline(tOut(check(2)),'-m','LineWidth',1);
grid on;
legend('','5% Error Time');
legend('Location','best');
xlabel('time (s)');
ylabel('$e_{y}$','interpreter','latex')
subplot(2,3,3);
hold on;
plot(tOut,error(:,3),'-k','LineWidth',1.2);
xline(tOut(check(3)),'-m','LineWidth',1);
grid on;
legend('','5% Error Time');
legend('Location','best');
xlabel('time (s)');
ylabel('$e_{z}$','interpreter','latex')
% velocities 
subplot(2,3,4);
hold on;
plot(tOut,error(:,4),'--b','LineWidth',1.2);
xline(tOut(check(4)),'-m','LineWidth',1);
grid on;
legend('','5% Error Time');
legend('Location','best');
xlabel('time (s)');
ylabel('$\dot{e_{x}}$','interpreter','latex')
subplot(2,3,5);
hold on;
plot(tOut,error(:,5),'--r','LineWidth',1.2);
xline(tOut(check(5)),'-m','LineWidth',1);
grid on;
legend('','5% Error Time');
legend('Location','best');
xlabel('time (s)');
ylabel('$\dot{e_{y}}$','interpreter','latex')
subplot(2,3,6);
hold on;
plot(tOut,error(:,6),'--k','LineWidth',1.2);
xline(tOut(check(6)),'-m','LineWidth',1);
grid on;
legend('','5% Error Time');
legend('Location','best');
xlabel('time (s)');
ylabel('$\dot{e_{z}}$','interpreter','latex')
%% LQR System
% ---[simulating system without observer]--- 
% determining ideal system 
Q = diag([3,5,3,1,1,1]);
R = [30,0,0;0,50,0;0,0,30];
N = 0; % cross-matrix terms 
[K_lqr,S_lqr,P_lqr] = lqr(sys,Q,R,N);
F_lqr = (C*((-A+(B*K_lqr))^-1)*B)^-1; % 'feed-forward' matrix
% defining new closed loop system
A_lqr = A-(B*K_lqr);
B_lqr = B*F_lqr;
C_lqr = C;
D_lqr = D;
% simulating optimal System
t_lqr = (0:0.1:60)';
u = r*ones(1,length(t_lqr));
sys4 = ss(A_lqr,B_lqr,C_lqr,D_lqr);
y_lqr = lsim(sys4,u,t_lqr,x0);
create_plots(y_lqr,u',t_lqr,3,x0,r,'Optimal LQR System Response u=-Kx+Fx');
% ---[simulating full system with observer]--- 
% redefining out observer
des_freq = 100*-w_n;
q_obs_new = [des_freq,des_freq,des_freq-(1i*imag(desired_poles_1)),des_freq-(1i*imag(desired_poles_1)),des_freq+(1i*imag(desired_poles_1)),des_freq+(1i*imag(desired_poles_1))]; % desired poles
L_new = place(A',C',q_obs_new)'; 
A_obs_new = A-(L_new*C);
% define system 
A_big = [A_lqr,B*K_lqr;zeros(6,6),A_obs_new]; 
B_big = [B_lqr;zeros(6,3)];
C_big = [1,0,0,0,0,0,0,0,0,0,0,0;
         0,1,0,0,0,0,0,0,0,0,0,0;
         0,0,1,0,0,0,0,0,0,0,0,0;
         0,0,0,0,0,0,1,0,0,0,0,0;
         0,0,0,0,0,0,0,1,0,0,0,0;
         0,0,0,0,0,0,0,0,1,0,0,0]; 
D_big = zeros(6,3);
% simulation 
x0_full = [10;10;10;0.2;0.2;0.2;0.6;2;0.01;0;0;0];
sys5 = ss(A_big,B_big,C_big,D_big);
y_lqr_full = lsim(sys5,u,t_lqr,x0_full);
% plotting entire system
figure;
sgtitle('Full LQR System Response')
    % state plots
        subplot(2,3,1);
        hold on;
        plot(t_lqr,y_lqr_full(:,1),'-b','LineWidth',1.2);
        plot(t_lqr,u,'--r','LineWidth',1);
        grid on;
        xlabel('time (s)');
        ylabel('$x$','interpreter','latex')
        str = sprintf('x_0=%0.2f | dx_0/dt=%0.2f',x0_full(1),x0_full(3));
        title(str)
        str = sprintf('r_x=%0.2f',r(2));
        subtitle(str);
        subplot(2,3,2);
        hold on;
        plot(t_lqr,y_lqr_full(:,2),'-r','LineWidth',1.2);
        plot(t_lqr,u,'--r','LineWidth',1);
        grid on;
        xlabel('time (s)');
        ylabel('$y$','interpreter','latex')
        str = sprintf('y_0=%0.2f | dy_0/dt=%0.2f',x0_full(2),x0_full(4));
        title(str)
        str = sprintf('r_y=%0.2f',r(2));
        subtitle(str);
        subplot(2,3,3);
        hold on;
        plot(t_lqr,y_lqr_full(:,3),'-k','LineWidth',1.2);
        plot(t_lqr,u,'--r','LineWidth',1);
        grid on;
        xlabel('time (s)');
        ylabel('$z$','interpreter','latex')
        str = sprintf('z_0=%0.2f | dz_0/dt=%0.2f',x0_full(3),x0_full(6));
        title(str)
        str = sprintf('r_z=%0.2f',r(2));
        subtitle(str);
    % error plots
        subplot(2,3,4);
        hold on;
        plot(t_lqr,y_lqr_full(:,4),'-b','LineWidth',1.2);
        grid on;
        xlabel('time (s)');
        ylabel('$e_{x}$','interpreter','latex')
        str = sprintf('e_{x0}=%0.2f | de_{x0}/dt=%0.2f',x0_full(7),x0_full(10));
        title(str)
        subplot(2,3,5);
        hold on;
        plot(t_lqr,y_lqr_full(:,5),'-r','LineWidth',1.2);
        grid on;
        xlabel('time (s)');
        ylabel('$e_{y}$','interpreter','latex')
        str = sprintf('e_{y0}=%0.2f | de_{y0}/dt=%0.2f',x0_full(8),x0_full(11));
        title(str)
        subplot(2,3,6);
        hold on;
        plot(t_lqr,y_lqr_full(:,6),'-k','LineWidth',1.2);
        grid on;
        xlabel('time (s)');
        ylabel('$e_{z}$','interpreter','latex')
        str = sprintf('e_{z0}=%0.2f | de_{z0}/dt=%0.2f',x0_full(9),x0_full(12));
        title(str)
