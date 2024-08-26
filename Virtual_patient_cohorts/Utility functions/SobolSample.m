function [sampleSet] = SobolSample(minvals, maxvals, d, n)
% Michelle Gee
% June 29, 2021
% Sobol sample generation to fit parameter space baroreceptor reflex model
% with ICN. Generates n random sets of parameter values for each
% transfer function
% Inputs
    % d = dimension of sample. Each dimension is a parameter
    % n = # of sobol points
    % minvals = minimum value of possible parameter range
    % maxvals = vector of maximum values of possible parameter range
% Output
    % Parameter values to sample as a matrix with each column the 
    % values for a single parameter 
    
if length(minvals) ~= d || length(maxvals) ~= d
    disp('Check minval and maxval vectors')
    return 
end
delta = zeros(d);
sampleSet = zeros(n,d);

for i = 1:d
    delta(i) = maxvals(i) - minvals(i);
end

p = sobolset(d); % generate sobol set
X0 = net(p,n); % display first n points

for i = 1:d
    sampleSet(:,i) = minvals(i) + X0(:,i) * delta(i);
end

%plot(X0(:,1),X0(:,2),'o')
end