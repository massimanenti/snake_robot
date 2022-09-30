% This script contains the main code of the Snake Robot OCP

clear all
close all
clc

fstar_best = realmax;               % initialization of the best value of the cost function found (realmax returns the largest finite floating-point number in IEEE® double precision)
initial_conditions_guesses = 1;     % counter, used to count the number of optimization problems solved (and consequently the number of initial guesses done)

for initial_conditions_index = 1:1:1
    
    %% Model parameters
    % In this section, some parameters of the mathematical model of the snake are defined, together with
    % some secondary parameter (related to the procedure instead of the snake itself). 
    % For semplicity, a struct called "parameters" is defined.
    % This struct is passed, as one of the inputs, in the needed functions.
    
    parameters.n        =   10;             % number of links (1)
    parameters.l        =   0.07;           % half length of a link (m)
    parameters.m        =   1;              % mass of each link (kg)
    parameters.J        =   0.0065;         % moment of inertia of each link (kg*m^2)
    parameters.ct       =   0.5;            % tangential viscous friction (N*s/m)
    parameters.cn       =   5;              % normal viscous friction (N*s/m)
    
    parameters.n_vary   =   5;              % n_vary stands for number of variations (1). The entire optimization is based on a parametrization of the reference trajectory of the snake joint angles. 
    % The parameters of the reference trajectory will serve as optimization variables. 
    % The reference trajectory will be characterized by 4 time varying parameters: alpha, omega, beta and gamma. 
    % The number of times they can change in a single simulation will be n_vary.
    % The number of optimization variables will therefore be 4*n_vary. 
    
    % Both n_vary and n are redefined here, so that from now on, in this script 
    % there is no need to write parameter.xxx to use them
    n_vary = parameters.n_vary;
    n = parameters.n;
    
    
    %% TARGET POINT COORDINATES
    
    % In this section, the coordinates of the goal point are defined. The user can
    % freely change them, remembering that the goal corresponds to the
    % desired final position of the center of mass of the robot.
    
    goal     =   [3 ; 5];   % x and y coordinates of the goal, both in (m)
    xgoal    = goal(1,1);
    ygoal    = goal(2,1);
    
    %% Time interval parameters
    
    % In this section, every parameter that concern sampling time or time intervals is
    % defined.
    
    Ts             =   0.1;                % Sampling time (s)
    Tend           =   11;                 % Time horizon (s)
    Nsim           =   Tend/Ts;            % Number of simulation steps (1)
    
    % Given that the parameters alpha, beta, omega and gamma will change at most
    % n_vary times during the simulation, it follows that there will be time
    % intervals when the parameters will remain constant. These time intervals will be a
    % designer choice in order to keep low the number of optimization
    % variables. However, the general idea is to keep the length of these
    % intervals as constant as possible, even if the first and last intervals
    % will be exceptions because they have some additional purposes during the
    % optimization.
    
    first_interval = floor(Nsim/n_vary/2);  % time step in which the first interval ends (1)
    last_interval = Nsim-floor(Nsim/n_vary/2);  % time step in which the last interval starts (1)
    
    % last_interval and T are then inserted into the parameters struct
    parameters.last_interval = last_interval;
    parameters.T = zeros(1, n_vary + 1);
    
    % The sequence of time instants is entirely represented by the vector T
    parameters.T(1,1) = 1;                                  % the first time istant is t = 1
    parameters.T(1, 2) = first_interval;                    % the first interval ends with t = first_interval
    parameters.T(1, n_vary) = last_interval;     % the last interval begins with t = last_interval
    parameters.T(1, n_vary + 1) = Nsim;          % the last interval ends with t = Nsim
    
    % In the case that n_vary is bigger than 3 (i.e. we have 4 or more
    % intervals), the duration of each interval (with the exception of the
    % first and last interval) is the same (in reality they could differ of
    % a few time steps due to round-off to the floor value). This behaviour
    % is enforced in the following lines of code:
    if n_vary > 3
        for ind_interval=3:n_vary-1
            parameters.T(1, ind_interval) = parameters.T(1, ind_interval - 1) + floor( (Nsim - (first_interval + Nsim - last_interval)) / (n_vary - 2));
        end 
    end
    % note that one could change the length of each interval by working
    % with parameter.T():
    % parameter.T(1, interval_number) = initial_time_step
    % where "initial_time_step" is the time step in which the
    % "interval_number" interval starts.
    
    %% Physical limits
    
    % A serpenoid curve is a mathematical description of lateral undulation 
    % that real snakes perform in order to achieve locomotion.
    %
    % The gait pattern lateral undulation is achieved by moving the joints
    % of a planar snake robot according to:
    %
    % phi_ref(ind_joint,1)= alpha*sin(omega*Ts*t+(ind_joint-1)*beta)+gamma;
    %
    % where ind_joint ∈ {1, ... , n−1}, alpha and omega are the amplitude and
    % angular frequency, respecitively, of the sinusoidal joint motion.
    % beta determines the phase shift between the joints and gamma is a joint
    % offset.
    % Ts*t is the time (we are working in discrete time)

    % In this section, the maximum and minimum limits for robot speed and
    % optimization variables are defined. Moreover, the weight used in the cost
    % function are also initialized and the vector "weights" is then created.
    % For each weight, a normalization term is computed.
    
    v_max = 1;         % speed limit of our robot, based on considerations about real applications and results of our simulation (m/s). Note that it is used to scale a cost function term, so it does not have to be the real maximum velocity value
    
    alpha_min = 0;      % min amplitude of the serpenoid curve (rad)
    alpha_max = 0.5;    % max amplitude of the serpenoid curve (rad)

    omega_min = 0;                  % min frequency of the serpenoid curve (rad/s)
    omega_max = min([30 2*pi/Ts]);  % max frequency of the serpenoid curve (rad/s)
    
    beta_min = pi/12;           % min offset between two consecutive joints. If beta is equal to 0, there is no forward propelling. (rad)
    beta_max = 2*pi-0.01;       % max offset between two consecutive joints. (rad)
    
    gamma_max = 3;    % not enforced constraint (rad)
    gamma_min = -3;   % not enforced constraint (rad)

    %% Cost function weights 
    
    % Each weight is composed of a normalization term (e.g. Q_norm) and a
    % tuning factor (e.g. q)
    
    Q_norm = sqrt((Nsim-1)*(goal(1,1)^2 + goal(2,1)^2));
    q           = sqrt(0.1);
    Q           =   q*1/Q_norm*diag(ones(2*Nsim-2,1));      % distance from goal weight
    
    Qf_norm = sqrt((goal(1,1)^2 + goal(2,1)^2));
    qf          = sqrt(0.89);
    Qf          =  qf*1/Qf_norm*eye(2);                     % distance from goal (terminal) weight
    
    Qv_norm = sqrt((Nsim-last_interval)*(2*v_max^2));
    qv          = sqrt(0.01);        % 350
    Qv          =  qv/Qv_norm*eye(2*(Nsim-last_interval));  % weights on last interval velocities
    
    
    weights.Q = Q;
    weights.Qf= Qf;
    weights.Qv= Qv;
    
    %% PD parameters
    
    % Here the parameters of the PD tracking controller are defined 
    kp = 43.1279;
    kd = 6.7288;
    
    %% Initial guess
    
    % In this section, the initial guesses for the optimization algorithm
    % are defined. 
    % Therefore, the vector of initial condition,
    % serpenoid_curve_parameters_0,
    % is computed, which is a column vector with 4*n_vary rows.
    % Note that for the first initial guess (initial_conditions_guesses == 1), 
    % a user defined initial guess is set,
    % whereas for all the other guesses a random choice is made
    % (still guaranteeing the initial feasibility for the linear constraints).

    serpenoid_curve_parameters_0 = zeros(4*n_vary, 1);   % vector initialization
        
    if initial_conditions_guesses == 1  % first initial guess
        alpha0 = pi/12;
        omega0 = 10;
        beta0 = pi/16;
        gamma0 = 0;
        for ind_vary=1:4:n_vary*4-3
            serpenoid_curve_parameters_0(ind_vary, 1) = alpha0;
            serpenoid_curve_parameters_0(ind_vary + 1, 1) = omega0;
            serpenoid_curve_parameters_0(ind_vary + 2, 1) = beta0;
            serpenoid_curve_parameters_0(ind_vary + 3, 1) = gamma0;
        end
    else    % other initial guesses
        for ind_vary=1:4:n_vary*4-3
            serpenoid_curve_parameters_0(ind_vary, 1) = alpha_min + (alpha_max - alpha_min)* rand;
            serpenoid_curve_parameters_0(ind_vary + 1, 1) = omega_min + (omega_max - omega_min)* rand;
            serpenoid_curve_parameters_0(ind_vary + 2, 1) = beta_min + (beta_max - beta_min)* rand;
            serpenoid_curve_parameters_0(ind_vary + 3, 1) = gamma_min + (gamma_max - gamma_min)* rand;
        end
    end
                                
    %% Constraints 
    % In this section, the main constraints are set.
    
    %  phi saturation constraints
    parameters.phi_min = -(pi-pi*(n-2)/n);
    parameters.phi_max =  (pi-pi*(n-2)/n);
    
    % alpha, beta, omega, gamma inequality saturation constraints.
    % The inequality constraints take the form C x >= d, so the matrix C and
    % the vector d must be defined.
    
    C0 = eye(length(serpenoid_curve_parameters_0));
    
    % Since gamma is not subjected to any constraint, the corresponding element
    % of matrix C must be 0.
    
    for ind_gamma = 4:4:4*n_vary
        C0(ind_gamma, ind_gamma) = 0;
    end    
    
    C = [ C0;
         -C0];
     
    % In the following lines, the vector d is built
     
    d = zeros(8*n_vary,1);
    
    for ind_d = 1:4:4*(n_vary-1)+1
        d(ind_d:ind_d+3) = [alpha_min; omega_min; beta_min; 0]; 
    end  
    
    for ind_d = 4*n_vary+1:4:4*(2*n_vary-1)+1
        d(ind_d:ind_d+3) = -[alpha_max; omega_max; beta_max; 0]; 
    end  
    
    
    % Number of nonlinear inequality constraints 
    p = 0;
    
    
    % Simulation time is divided in 3 time slots: interval_phi = [1, floor(Nsim/3), floor(2*Nsim/3), Nsim]'
    % A "step_phi" value is associated to each interval: step_phi = [3, 3, 3]'
    % The idea is to constrain the value of all the phis of the robot to be
    % inside an interval (given by phi_min and phi_max). In order to avoid
    % bounding all the phis for all the simulation time, the constaint is given
    % every "step_phi" time steps to every phi of the robot.
    
    interval_phi = [1, floor(Nsim/3), floor(2*Nsim/3), Nsim]';
    step_phi = [3, 3, 3]';
    parameters.interval_phi = interval_phi;
    parameters.step_phi = step_phi;
    
    % number of inequality constraints
    q = 2*(n-1)*( ceil((interval_phi(4,1)-(interval_phi(3,1)))/step_phi(1,1)) + ceil((interval_phi(3,1)-(interval_phi(2,1)))/step_phi(2,1)) + ceil((interval_phi(2,1)-(interval_phi(1,1)-1))/step_phi(3,1)) );
                                 
    %% INITIALIZATION OF STATE
    % In this section, each variable is initialized. Moreover some variables and
    % matrices needed to develop the model are defined. 
    
    % In order to define the initial position of the links (i.e. theta0(i)), we consider the
    % serpenoid curve, with given parameters.
    %
    % theta0(i) represents the angle that the i-th link of the snake robot
    % forms with the global x axis (positive CCW), at the initial time.
    
    alpha00 = pi/12;       % amplitude of serpenoid curve [rad]
    omega00 = 10;         % angular frequency of serpenoid curve [rad/s]
    beta00 = pi/12;        % phase shift in the serpenoid curve, between a link and its successor [rad]
    gamma00 = 0;          % joint offset, parameter of serpenoid curve [rad]
    
    theta0 = zeros(n,1);
    for ind=1:n
        theta0(ind) = alpha00*sin(omega00*Ts+(ind-1)*beta00) + gamma00; 
    end 

    % phi0(i) is the angle between two joints, i.e. theta0(i)-theta0(i+1), at the initial time step.
    % NOTE: phi(i) with i ∈ {1, ... , n−1} and theta(n) will be directly in
    % the state vector zsim of the simulation (with other variables, e.g. speeds).
    %
    % phi_bar0 is the column vector containing all the phi(i) (from 1 to n-1) as
    % first elements and theta(n) as the last one, at the initial time, i.e.:
    % phi_bar0 = [phi0(1); 
    %             phi0(2); 
    %              ... ; 
    %             phi0(n-1); 
    %             theta0(n)]
    %
    % H is the matrix that relates theta0 and phi_bar0, phi_bar0 =
    % H^(-1)*theta0
    
    H = zeros(n);           
    for ind_joint=1:n
        for column=1:n
            if column >= ind_joint
                H(ind_joint,column)=1;  
            end 
        end    
    end
    
    % From theta0 and H, phi_bar0 is computed
    phi_bar0 = H^-1*theta0;
    
    % Since we consider the robot still at the initial time step, the time derivatives should be equal to zero.
    % position of center of mass is in the zero of the global reference frame, velocity of center of mass is equal to 0.
    % For this reason the remaining variables are initialized with all
    % zeros:
    px0 = zeros(1, 1);
    py0 = zeros(1, 1);
    phi_dot0 = zeros(n-1, 1);
    theta_n_dot0 = zeros(1, 1);
    px_dot0 = zeros(1, 1);
    py_dot0 = zeros(1, 1);
    
    % All the previous variables are then grouped together in the initial state vector
    z0 = zeros(2*n+4, 1); 
    z0(:,:) = [phi_bar0; px0; py0; phi_dot0; theta_n_dot0; px_dot0; py_dot0];
    % NOTE: here we can appreciate the variables that compose 
    % the state vector zsim of the simulation:
    % zsim will be a sequence of 2*n+4 elements column vectors (each column represent a different time instant):
    % zsim(:, 1) = [phi0(1);        % first phi angle
    %               phi0(2);        ...
    %                ...            ...
    %               phi0(n-1);      % n-1-th phi angle
    %               theta0(n);      % orientation of the n-th link
    %               px0;            % x coordinate of center of mass
    %               py0;            % y coordinate of center of mass
    %               phi_dot0(1);    % angular velocity related to the first phi angle
    %               phi_dot0(2);    ...
    %                ...            ...
    %               phi_dot0(n-1);  % angular velocity related to the n-1-th phi angle
    %               theta_n_dot0;   % angular velocity of the n-th link
    %               px_dot0;        % x velocity of center of mass
    %               py_dot0]        % y velocity of center of mass
    
    %% Optimization parameters
    % In this section, some of the parameters of myoptions are set.
    
    myoptions               =   myoptimset;
    myoptions.ls_beta       = 0.3;        
    myoptions.ls_c          = 0.1;
    myoptions.gradmethod    = 'CD';     
    myoptions.graddx        = eps^(1/3);
    myoptions.nitermax      = 300;
    myoptions.Hessmethod    = 'GN';
    % myoptions.GN_funF       =
    % @(serpenoid_curve_parameters)cost_function_con_snake_GN(serpenoid_curve_parameters,z0,Ts, Tend,parameters,goal,kp,kd, weights);       
    myoptions.GN_funF       = @(serpenoid_curve_parameters)cost_function_con_snake_GN_mex(serpenoid_curve_parameters,z0,Ts, Tend,parameters,goal,kp,kd, weights);    % if using compiled version
    
    %% Generating code (mex files) for faster computation (to be run only once whenever Tend changes)
    
    codegen cost_function_con_snake -args {serpenoid_curve_parameters_0, z0, Ts, Tend, parameters, goal, kp, kd, weights}
    codegen cost_function_con_snake_GN -args {serpenoid_curve_parameters_0,z0,Ts, Tend, parameters, goal, kp, kd, weights}
                                     
    %% Optimization routine
    % In this section, the optimization routine is performed, calling the
    % function con_optimization_routine
    
    % The following line does not use the generated mex code
    % [serpenoid_curve_parameters_star, fxstar, k, exitflag, xsequence] = con_optimization_routine(@(serpenoid_curve_parameters)cost_function_con_snake(serpenoid_curve_parameters, z0, Ts, Tend, parameters, goal, kp, kd, weights), serpenoid_curve_parameters_0, [], [], C, d, p, q,  myoptions);
    
    % The following line does use the generated mex code
    [serpenoid_curve_parameters_star, fxstar, k, exitflag, xsequence] = con_optimization_routine(@(serpenoid_curve_parameters)cost_function_con_snake_mex(serpenoid_curve_parameters, z0, Ts, Tend, parameters, goal, kp, kd, weights), serpenoid_curve_parameters_0,[], [], C, d, p, q, myoptions);
        
    % The cost function is computed one last time, with the optimal values of
    % the optimization variables
    v = cost_function_con_snake(serpenoid_curve_parameters_star,z0,Ts, Tend, parameters, goal, kp, kd, weights);
    fstar = v(1,1);

    if fstar < fstar_best
        serpenoid_curve_parameters_star_best = serpenoid_curve_parameters_star;
        fstar_best = fstar;
    end

    initial_conditions_guesses = initial_conditions_guesses + 1;
end


%% Post-processing the results and plots
% In this section, all the post processing is performed. 
% To set up the animation:
animation_on = 1;   % set it to 1 to see the animation
time_of_pause_animation = 0.1;   % change this parameter will show the animation slower (higher value) or faster (lower value)

% Simulation performed with the optimal variables in order to obtain all
% the informations needed for the plots
[zstar, Xj, Yj, phi_ref_plot, u_plot] = snake_trajectory_with_reference(serpenoid_curve_parameters_star_best,z0,Ts, Tend, parameters, goal, kp, kd, animation_on, time_of_pause_animation); 

% All the variables are then extracted from the optimal state vector zstar
phi(1:n-1, :) = zstar(1:n-1, :);
theta_n(1, :) = zstar(n, :);
px(1, :) = zstar(n+1, :);
py(1, :) = zstar(n+2, :);
phi_dot(1:n-1, :) = zstar(n+3:2*n+1, :);
theta_n_dot(1, :) = zstar(2*n+2, :);
px_dot(1, :) = zstar(2*n+3, :);
py_dot(1, :) = zstar(2*n+4, :);

vect_t_total = Ts * (1:Nsim);                                          % time vector used for the plots
start_iteration_plot = 1;                                       % selects first sample to plot (= 1 runs it from the beginning)
end_iteration_plot =Nsim;                                      % selects last sample to plot (= Nsim runs it untill the end)
vect_t_selected = Ts* (start_iteration_plot:end_iteration_plot);

if animation_on == 1
    % The single components of the cost function are computed here:


    F1(:, 1)  = Q*[ones(Nsim-1,1)*xgoal - px(1,1:Nsim-1)';      % distance from goal (x axis)
                   ones(Nsim-1,1)*ygoal - py(1,1:Nsim-1)'];     % distance from goal (y axis)
      
      
    F2(:, 1) = Qf*[xgoal - px(1,Nsim)';        % distance from goal (x axis), last time instant
                   ygoal - py(1,Nsim)'];          % distance from goal (y axis), last time instant
                         
    F3(:, 1) = Qv*[ -px_dot(1,last_interval+1:Nsim)';  % last interval x velocity
                    -py_dot(1,last_interval+1:Nsim)']; % last interval y velocity



    
    f1 = F1'*F1
    f2 = F2'*F2
    f3 = F3'*F3

    f = f1 + f2 + f3



%     Q_special           =   1/Q_norm*diag(ones(2*Nsim-2,1));
%     Qf_special          =   1/Qf_norm*eye(2); 
%     Qv_special          =   1/Qv_norm*eye(2*(Nsim-last_interval));
% 
%     F1_special(:, 1)  = Q_special*[ones(Nsim-1,1)*xgoal - px(1,1:Nsim-1)';      % distance from goal (x axis)
%                    ones(Nsim-1,1)*ygoal - py(1,1:Nsim-1)'];     % distance from goal (y axis)
%       
%       
%     F2_special(:, 1) = Qf_special*[xgoal - px(1,Nsim)';        % distance from goal (x axis), last time instant
%                    ygoal - py(1,Nsim)'];          % distance from goal (y axis), last time instant
%                          
%     F3_special(:, 1) = Qv_special*[ -px_dot(1,last_interval+1:Nsim)';  % last interval x velocity
%                     -py_dot(1,last_interval+1:Nsim)']; % last interval y velocity
% 
%     f1_special = F1_special'*F1_special
%     f2_special = F2_special'*F2_special
%     f3_special = F3_special'*F3_special
% 
%     f_special = f1_special + f2_special + f3_special
end


%----------------------------------------------------------------------------------------------------------------------------------
% Plot setup: it is possible to choose wich plot to visualize by setting
% the corresponding variable to "1"
plot_of_com_traj = 1;                                           % 1 to have the plot of the trajectory of the center of mass 

plot_of_all_phi = 0;                                            % 1 to have the plot of all the phi angles (and their references) with respect to time 

plot_of_chosen_phi = 1;                                         % 1 to have the plot of some selected phi angles (and their references) that are chosen in the "phi_to_be_plotted" vector
phi_to_be_plotted = [1,2,3];                                        % in a row vector you can put which phi to plot, for example: phi_to_be_plotted = [1, 3, 5, 6];

plot_of_px_dot = 1;                                             % 1 to have the plot of px_dot with respect to time 
plot_of_py_dot = 1;                                             % 1 to have the plot of py_dot with respect to time 
plot_of_mod_speed = 1;                                          % 1 to have the plot of the velocity of the center of mass of the robot with respect to time 

plot_of_alpha = 1;                                              % 1 to have the plot of alpha with respect to time and the corresponding limits
plot_of_omega = 1;                                              % 1 to have the plot of omega with respect to time and the corresponding limits
plot_of_beta  = 1;                                              % 1 to have the plot of beta with respect to time and the corresponding limits
plot_of_gamma = 1;                                              % 1 to have the plot of gamma with respect to time and the corresponding limits

plot_of_chosen_u = 1;
u_to_be_plotted = [1,2,3,4,5];
%-----------------------------------------------------------------------------------------------------------------------------------

% 1. plot of the center of mass in the 2D plane
if plot_of_com_traj == 1
    figure,
    plot(px(1,:),py(1,:)),
    hold on
    for interval=1:n_vary-1
        
        plot(px(1,parameters.T(1, interval+1)), py(parameters.T(1, interval+1)),'ok','markersize',5,'markerfacecolor','g')
        hold on
    end
    plot(xgoal,ygoal,'ok','markersize',2,'markerfacecolor','b')
    hold on 

    grid on, 
    xlabel('X (m)'),
    ylabel('Y (m)'),
    title('C.O.M. Trajectory')
end

% 2. plot of ALL joint angles compared to their reference (with respect to time)
if plot_of_all_phi == 1
    for i=1:n-1
        
        figure 
        hold on
        plot(vect_t_selected, phi(i,start_iteration_plot:end_iteration_plot),'DisplayName',strcat('phi ',num2str(i))) 
        plot(vect_t_selected, phi_ref_plot(i,start_iteration_plot:end_iteration_plot), 'DisplayName',strcat('phi ',num2str(i), ' ref'))
        plot(vect_t_selected, parameters.phi_min*ones(size(start_iteration_plot:end_iteration_plot)),'DisplayName','lower limit')
        plot(vect_t_selected, parameters.phi_max*ones(size(start_iteration_plot:end_iteration_plot)),'DisplayName', 'upper limit')        
        grid on, 
        hold off
        xlabel('time (sec)'),
        
        ytxt = ['phi', num2str(i), '(t) (rad)'];
        ylabel(ytxt),
        
        titletxt = ['Phi ' num2str(i), ' VS Phi ', num2str(i), ' ref'];
        title(titletxt),
        
        legend
        
    end
end

% 3. plot of CHOSEN joint angles compared to its reference (with respect to time)
if plot_of_chosen_phi == 1
    for index_plot_chosen_phi = phi_to_be_plotted
        
        figure 
        hold on
        plot(vect_t_selected,zstar(index_plot_chosen_phi,start_iteration_plot:end_iteration_plot),'DisplayName',strcat('phi ',num2str(index_plot_chosen_phi)))
        plot(vect_t_selected, phi_ref_plot(index_plot_chosen_phi,start_iteration_plot:end_iteration_plot),'DisplayName',strcat('phi ',num2str(index_plot_chosen_phi),' ref')),
        plot(vect_t_selected, parameters.phi_min*ones(size(start_iteration_plot:end_iteration_plot)),'DisplayName','lower limit')
        plot(vect_t_selected, parameters.phi_max*ones(size(start_iteration_plot:end_iteration_plot)),'DisplayName','upper limit')
        grid on,
        
        xlabel('time (sec)'),
        
        ytxt = ['phi', num2str(index_plot_chosen_phi), '(t) (rad)'];
        ylabel(ytxt),
        
        titletxt = ['Phi ', num2str(index_plot_chosen_phi), ' VS Phi ', num2str(index_plot_chosen_phi), ' ref'];
        title(titletxt),
        
        legend
        
    end
end    

% 4. plot of px_dot (with respect to time)
if plot_of_px_dot == 1
    figure
    plot (vect_t_total, px_dot(1,:))
    grid on, 
    xlabel('time (sec)'),
    ylabel('px dot (m/s)')
    title('Px dot')
end
    
% 5. plot of py_dot (with respect to time)
if plot_of_py_dot == 1
    figure
    plot (vect_t_total,py_dot(1,:)) 
    grid on, 
    xlabel('time (sec)'),
    ylabel('py dot (m/s)')
    title('Py dot')
end

% 6. plot of the module of total velocity (with respect to time)
if plot_of_mod_speed == 1
    for ind_mod_vel= 1:length(px)
        mod_vel(1,ind_mod_vel) = sqrt(px_dot(1,ind_mod_vel)^2 + py_dot(1,ind_mod_vel)^2);
    end
    figure
    plot (vect_t_total,mod_vel(1,:)) 
    grid on, 
    xlabel('time (sec)'),
    ylabel('|vel| (m/s)')
    title('Module of the total velocity')
end

% 7. plot of alpha (with respect to time)
if plot_of_alpha == 1
    
    figure
    alpha_vector_star = serpenoid_curve_parameters_star_best(1:4:4*(n_vary-1)+1);
    
    hold on
    plot (time_vector_for_descrete_plots(n_vary), y_vector_for_descrete_plots(alpha_vector_star))
    plot (0:n_vary, alpha_max*ones(n_vary+1,1))
    plot (0:n_vary, alpha_min*ones(n_vary+1,1))
    
    grid on, 
    xlabel('time intervals'),
    ylabel('alpha (rad)')
    title('Alpha wrt time')
      
    legend('alpha wrt time','alpha max','alpha min')
end

% 8. plot of omega (with respect to time)
if plot_of_omega == 1
    
    figure
    omega_vector_star = serpenoid_curve_parameters_star_best(2:4:4*(n_vary-1)+2);
    
    hold on
    plot (time_vector_for_descrete_plots(n_vary), y_vector_for_descrete_plots(omega_vector_star)) 
    plot (0:n_vary, omega_max*ones(n_vary+1,1))
    plot (0:n_vary, omega_min*ones(n_vary+1,1))
    
    grid on, 
    xlabel('time interval'),
    ylabel('omega (rad/s)')
    title('Omega wrt time')
    
    legend('omega wrt time','omega max','omega min')
end

% 9. plot of beta (with respect to time)
if plot_of_beta == 1
    
    figure
    beta_vector_star = serpenoid_curve_parameters_star_best(3:4:4*(n_vary-1)+3);
    
    hold on
    plot (time_vector_for_descrete_plots(n_vary), y_vector_for_descrete_plots(beta_vector_star)) 
    plot (0:n_vary, beta_max*ones(n_vary+1,1))
    plot (0:n_vary, beta_min*ones(n_vary+1,1))
    
    grid on, 
    xlabel('time interval'),
    ylabel('beta (rad)')
    title('Beta wrt time')

    legend('beta wrt time','beta max','beta min')
end

% 10. plot of gamma (with respect to time)
if plot_of_gamma == 1
    
    figure
    gamma_vector_star = serpenoid_curve_parameters_star_best(4:4:4*n_vary);
    
    hold on
    plot (time_vector_for_descrete_plots(n_vary), y_vector_for_descrete_plots(gamma_vector_star)) 
    plot (0:n_vary, gamma_max*ones(n_vary+1,1))
    plot (0:n_vary, gamma_min*ones(n_vary+1,1))
    
    grid on, 
    xlabel('time interval'),
    ylabel('gamma (rad)')
    title('Gamma wrt time')
    
    legend('gamma wrt time','gamma max','gamma min')
end

% 11 plot of u 
if plot_of_chosen_u == 1
    for index_plot_chosen_u = u_to_be_plotted
        
        figure 
        hold on
        plot(vect_t_total,u_plot(index_plot_chosen_u,:),'DisplayName',strcat('u ',num2str(index_plot_chosen_u)))
        grid on,
        
        xlabel('time (sec)'),
        
        ytxt = ['u', num2str(index_plot_chosen_u), '(t) (Nm)'];
        ylabel(ytxt),
        
        titletxt = ['u', num2str(index_plot_chosen_u)];
        title(titletxt),
        
        legend
        
    end
end    






