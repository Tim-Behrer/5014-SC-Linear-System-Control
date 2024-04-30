function [] = simulate_response(sys1,sys2,n)
%Given a state space model, plot the response of the system

%% Plotting
if n == 1
    figure();
    grid on
    SR = stepplot(sys1);
    SR.InputName = {'Radial Gas Thruster','Along-Track Gas Thruster','Cross-Track Gas Thruster'};
    SR.OutputName = {'Radial Gas Position','Along-Track Position','Cross-Track Position'};
    figure();
    IR = impulseplot(sys1);
    IR.InputName = {'Radial Gas Thruster','Along-Track Gas Thruster','Cross-Track Gas Thruster'};
    IR.OutputName = {'Radial Gas Position','Along-Track Position','Cross-Track Position'};
elseif n == 2
    figure();
    grid on
    SR = stepplot(sys1,sys2);
    setoptions(SR,'Normalize','on')
    SR.InputName(:) = {'Radial Gas Thruster','Along-Track Gas Thruster','Cross-Track Gas Thruster'};
    SR.OutputName(:) = {'Radial Gas Position','Along-Track Position','Cross-Track Position'};
    figure();
    IR = impulseplot(sys1,sys2);
    setoptions(IR,'Normalize','on')
    IR.InputName(:) = {'Radial Gas Thruster','Along-Track Gas Thruster','Cross-Track Gas Thruster'};
    IR.OutputName(:) = {'Radial Gas Position','Along-Track Position','Cross-Track Position'};
else
    fprintf("Error: Number of comparisons > 2")
end

end