function [RRinterval, RRtimes] = RRinterval_calc(simOut,trange)
% Calculate RR intervals from Simulink simulation output (simOut)
% Inputs: 
    % simOut:  Simulink simulation output 
    % trange: 1x2 vector with simulation start and stop time for calculating RR interval
% Outputs: 
    % RRinterval: row vector of RR intervals
    % RRtimes: row vector of times corresponding to RR intervals

% Michelle Gee, 2/16/24
%% Simulation input
simTime         = simOut.time;

% convert intervals to model indices
tbaseIdx        = find(simTime>= trange(1)  & simTime <= trange(2));

time = simOut.time(tbaseIdx);
phi = simOut.Phi(tbaseIdx);

% RR interval calculation for baseline interval
% Determine indices for when phi=0 (beginning and end of heart beat)
for j = length(phi):-1:2
    if phi(j) == phi(j-1)
        phi(j-1) = [];
        time(j-1) = [];
    end
end

RRidx = find(~phi);
RRtimes = time(RRidx);
RRinterval = diff(RRtimes);
end
