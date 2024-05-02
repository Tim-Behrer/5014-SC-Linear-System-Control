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
clear;
close all;
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
%% Uncontrolled system response
sys_uncontrolled = ss(A,B,C,D);
simulate_response(sys_uncontrolled,[],1)
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
%% Closed Loop system response
sys_closedLoop = ss(A_cl,B_cl,C_cl,D_cl);

simulate_response(sys_closedLoop,[],1) %Simulate just closed loop
simulate_response(sys_uncontrolled,sys_closedLoop,2) %Simulate uncontrolled and closed loop
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


%% Cost Function Development
Q = eye(6).*[3 5 3 1 1 1];
R = 10*Q(1:3,1:3);

[K_opt,~,eig_opt] = lqr(sys_uncontrolled,Q,R,0);
F_opt = (C*((-A+(B*K_opt))^-1)*B)^-1; % 'feed-forward' matrix
% defining optimal loop system
A_opt = A-(B*K_opt);
B_opt = B*F_opt;
C_opt = C;
D_opt = D;
sys_opt = ss(A_opt,B_opt,C_opt,D_opt);
[SI1, ~] = simulate_response(sys_opt,[],1);
[SI1, SI2] = simulate_response(sys_opt,sys_closedLoop,2);
