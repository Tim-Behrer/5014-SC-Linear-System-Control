function [] = create_plots(y,u,t,n,x0,r,plot_title)
% Creates a nice set of plots to show our system response 
% INPUTS:
%   y: system output (length(time),n)
%   u: system control signal
%   t: time vector [s]
%   n: number of states 
%   x0: initial condition
%   r: desired position
%   plot_title: string for the plot title
% OUTPUTS:
%   1,n set of plots
%% Code
figure
sgtitle(plot_title); 
subplot(1,n,1);
hold on;
grid minor;
plot(t,y(:,1),'-b','LineWidth',1.2);
plot(t,u(:,1),'--r','LineWidth',1);
xlabel('Time [s]')
ylabel('x')
str = sprintf('x_0=%0.2f | dx_0/dt=%0.2f',x0(1),x0(4));
title(str)
str = sprintf('r_x=%0.2f',r(1));
subtitle(str);
legend('Output','Control Signal')
subplot(1,n,2);
hold on;
grid minor;
plot(t,y(:,2),'-b','LineWidth',1.2);
plot(t,u(:,2),'--r','LineWidth',1);
xlabel('Time [s]')
ylabel('y')
legend('Output','Control Signal')
str = sprintf('y_0=%0.2f | dy_0/dt=%0.2f',x0(2),x0(5));
title(str)
str = sprintf('r_y=%0.2f',r(2));
subtitle(str);
subplot(1,n,3);
hold on;
grid minor;
plot(t,y(:,3),'-b','LineWidth',1.2);
plot(t,u(:,3),'--r','LineWidth',1);
xlabel('Time [s]')
ylabel('z')
legend('Output','Control Signal')
str = sprintf('z_0=%0.2f | dz_0/dt=%0.2f',x0(3),x0(6));
title(str)
str = sprintf('r_z=%0.2f',r(3));
subtitle(str);

end

