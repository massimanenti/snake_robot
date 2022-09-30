function v = time_vector_for_descrete_plots(T)

% This function is used in main snake matlab file, and in particular in the
% post processing section. It is used in alpha, omega, beta and gamma plots
% in order to plot them as piecewise constant signals.

v = zeros((T+1-2)*2+2,1)';
v(1) = 0;
v(2*T) = T;
for i = 2:2:2*T-1
    v(i)=i/2;
    v(i+1)= v(i);
end
