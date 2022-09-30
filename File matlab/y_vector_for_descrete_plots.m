function v = y_vector_for_descrete_plots(input_vector)

% This function is used in main snake matlab file, and in particular in the
% post processing section. It is used in alpha, omega, beta and gamma plots
% in order to plot them as piecewise constant signals.

n = length(input_vector);
v = zeros(2*n,1)';

for i = 1:n
    v((i-1)*2+1) = input_vector(i);
    v((i-1)*2+2)= v((i-1)*2+1);
end
