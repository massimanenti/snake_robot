function  v = cost_function_snake(serpenoid_curve_parameters, z0, Ts, Tend, parameters, goal, kp, kd, weights)

% function that simulates the snake robot, computes the cost function and the constraints.
% 
% OUTPUT
% v           =   [f;h];
% f represents the cost function value (that should be minimized)
% h represents the nonlinear inequality constraints
% 
% INPUT
% serpenoid_curve_parameters, optimization variables (alpha, omega, beta, gamma in the time intervals)
% z0, initial condition
% Ts, sampling time
% Tend, time of the simulation
% parameters, structure where snake robot parameters are stored
% goal, x and y coordinate of the goal
% kp, closed loop control constant
% kd, closed loop control constant
% weights, weights of the cost function

% Taking the parameters from the data structure
n = parameters.n;   % number of links
l = parameters.l;   % length of the single link
n_vary = parameters.n_vary; % number of intervals in which time is divided
last_interval = parameters.last_interval;   % starting time instant of the last interval
phi_min = parameters.phi_min;   % constraint on phi_min
phi_max = parameters.phi_max;   % constraint on phi_max

% Taking the weight matrices from the data structure
Q = weights.Q;      % distance from goal weight
Qf = weights.Qf;    % distance from goal (terminal point)
Qv = weights.Qv;    % weights on last interval velocities


% Taking remaining parameters
xgoal = goal(1,1);  % x coordinate of the goal
ygoal = goal(2,1);  % y coordinate of the goal

px0 = z0(n+1,1);    % x position of the snake robot (initial condition)
py0 = z0(n+2,1);    % y position of the snake robot (initial condition)

Nsim = Tend/Ts;

%% COST FUNCTION INITIALIZATION

% F is used to compute the cost function (f = F'*F).
% It is a column vector, where different costs are progressively updated

F = [zeros(Nsim-1,1);               % distance from goal (x axis)
     zeros(Nsim-1,1);               % distance from goal (y axis)
     zeros(1,1);                    % distance from goal (x axis), last time instant
     zeros(1,1);                    % distance from goal (y axis), last time instant
     zeros(Nsim-last_interval,1);   % last interval x velocity
     zeros(Nsim-last_interval,1)   % last interval y velocity
     ];          

%% INITIALIATION OF REFERENCE SIGNAL FOR THE CONTROLLER AND CONTROLLER PARAMETERS

phi_ref = zeros(n-1,1);
phi_ref_dot = zeros(n-1,1);
phi_ref_ddot = zeros(n-1,1);

phi = zeros(n-1,Nsim);
phi_dot = zeros(n-1,Nsim);
theta_n= zeros(1,Nsim);
px = zeros(1, Nsim);
py = zeros(1, Nsim);
theta_n_dot = zeros(1, Nsim);
px_dot = zeros(1, Nsim);    
py_dot = zeros(1, Nsim);
mod_vel = zeros(1, Nsim);   % absolute value of the velocity of the center of mass of the robot

u = zeros(n-1,1);       % the vector of inputs (actuator torques)
                        % in each joint of the snake robot there is a motor
                        % capable of generating a torque.

                        
% The state of our model will be the vector zsim.
% zsim is composed of 2*n+4 elements column vectors, where each column is referred to a time instant:
% zsim(:, generic_time) = [phi(1);      % phi(i) is the i-th phi angle
%                          phi(2);
%                           ...
%                          phi(n-1);
%                          theta(n);    % n-th theta angle
%                          px;          % x coordinate of center of mass
%                          py;          % y coordinate of center of mass
%                          phi_dot(1);  % the last elements are the velocities of the first elements
%                          phi_dot(2);
%                           ...
%                          phi_dot(n-1);
%                          theta_dot(n);
%                          px_dot;
%                          py_dot]
zsim = zeros(2*n+4,Nsim);   % states of the snake robot (rows represent variables, columns the evolution in time)
zsim(:,1)=z0;

T = parameters.T;   % vector where the initial times of the intervals are stored (last element is the last time instant of the simulation)


%% SIMULATION

% n_vary stands for number of variations (1). The entire optimization is based on a parametrization of the reference trajectory of the snake joint angles. 
% The parameters of the reference trajectory will serve as optimization variables. 
% The reference trajectory will be characterized by 4 time varying parameters: alpha, omega, beta and gamma. 
% The number of times they can change in a single simulation will be n_vary.
% The number of optimization variables will therefore be 4*n_vary.

for ind_vary=1:n_vary       % ind_vary: index to select the time intervals
        
    % optimization variables 8used in the definition of the reference signal):
    alpha = serpenoid_curve_parameters(4*(ind_vary-1)+1,1);
    omega = serpenoid_curve_parameters(4*(ind_vary-1)+2,1);
    beta  = serpenoid_curve_parameters(4*(ind_vary-1)+3,1);
    gamma = serpenoid_curve_parameters(4*(ind_vary-1)+4,1);
    
    for t=T(1,ind_vary)+1:T(1,ind_vary+1)   % t: index to select the time istants

        for ind_joint=1:n-1  % ind_joint: index to select the joints         

            % COMPUTATION OF PHI_REF and its derivatives at each time instant.
            % The equations of these references come from the literature
            % that studied snake's movement.
            % They are single column vectors.

            phi_ref(ind_joint,1)= alpha*sin(omega*Ts*(t-1)+(ind_joint-1)*beta)+gamma;
            phi_ref_dot(ind_joint,1)= alpha*omega*cos(omega*Ts*(t-1)+(ind_joint-1)*beta);
            phi_ref_ddot(ind_joint,1)= -alpha*omega^2*sin(omega*Ts*(t-1)+(ind_joint-1)*beta);

        end
        
        % SIMULATION    
        % now the vector of state is computed in descrete time with Euler Forward Method
        % Xj, Yj contains the coordinates of the joints, they are not
        % needed here        
        [zdot, Xj, Yj, mult, add] = model_simulation_snake(0,zsim(:,t-1), u, 0, parameters);
        zsim(:,t) = zsim(:,t-1) + Ts*zdot;
        
        phi(1:n-1, t) = zsim(1:n-1, t);
        theta_n(1, t) = zsim(n, t);
        px(1, t) = zsim(n+1, t);
        py(1, t) = zsim(n+2, t);
        phi_dot(1:n-1, t) = zsim(n+3:2*n+1, t);
        theta_n_dot(1, t) = zsim(2*n+2, t);
        px_dot(1, t) = zsim(2*n+3, t);
        py_dot(1, t) = zsim(2*n+4, t);

        % closed loop control action:
        u = phi_ref_ddot + kd*(phi_ref_dot - phi_dot(1:n-1,t)) + kp*(phi_ref - phi(1:n-1,t));
        mod_vel(1,t) = sqrt( px_dot(1,t)^2 + (py_dot(1,t))^2 ); % absolute value of the velocity of the center of mass of the robot
        
    end
end

%% Constraint

% Constraint on straightness
% Simulation time is divided in 3 time slots: interval_phi = [1, floor(Nsim/3), floor(2*Nsim/3), Nsim]'
% these 3 time solts correspond to the 3 "for" cycles below.
% A "step_phi" value is associated to each interval: step_phi = [3, 3, 3]'
% The idea is to constrain the value of all the phis of the robot to be
% inside an interval (given by phi_min and phi_max). In order to avoid
% bounding all the phis for all the simulation time, the constaint is given
% every "step_phi" time steps to every phi of the robot.

ind_h = 1;
interval_phi = parameters.interval_phi;
step_phi = parameters.step_phi;

h = zeros(2*(n-1)*( ceil((interval_phi(4,1)-(interval_phi(3,1)))/step_phi(1,1)) + ceil((interval_phi(3,1)-(interval_phi(2,1)))/step_phi(2,1)) + ceil((interval_phi(2,1)-(interval_phi(1,1)-1))/step_phi(3,1)) ), 1);

for colonna=interval_phi(1,1):step_phi(1,1):interval_phi(2,1)
 
    h(ind_h:ind_h+2*(n-1)-1,1) = [ phi(1:n-1, colonna) - phi_min * ones(n-1, 1);
                                  -phi(1:n-1,colonna) + phi_max * ones(n-1, 1)];
 
    ind_h = ind_h + 2*(n-1);
end   

for colonna=interval_phi(2,1)+1:step_phi(2,1):interval_phi(3,1)
 
    h(ind_h:ind_h+2*(n-1)-1,1) = [ phi(1:n-1, colonna) - phi_min * ones(n-1, 1);
                                  -phi(1:n-1,colonna) + phi_max * ones(n-1, 1)];
 
    ind_h = ind_h + 2*(n-1);
end   

for colonna=interval_phi(3,1)+1:step_phi(3,1):interval_phi(4,1)
 
    h(ind_h:ind_h+2*(n-1)-1,1) = [ phi(1:n-1, colonna) - phi_min * ones(n-1, 1);
                                  -phi(1:n-1,colonna) + phi_max * ones(n-1, 1)];
 
    ind_h = ind_h + 2*(n-1);
end   


%% UPDATE COST FUNCTION 

F(1:2*Nsim-2, 1)  = Q*[ones(Nsim-1,1)*xgoal - px(1,1:Nsim-1)';      % distance from goal (x axis)
                       ones(Nsim-1,1)*ygoal - py(1,1:Nsim-1)'];     % distance from goal (y axis)
  
  
F(2*Nsim-1:2*Nsim, 1) = Qf*[xgoal - px(1,Nsim)';        % distance from goal (x axis), last time instant
                         ygoal - py(1,Nsim)'];          % distance from goal (y axis), last time instant
                     
F(2*Nsim+1:4*Nsim-2*last_interval, 1) = Qv*[ -px_dot(1,last_interval+1:Nsim)';  % last interval x velocity
                                             -py_dot(1,last_interval+1:Nsim)']; % last interval y velocity
                                        
  
  f = F'*F;

%% Stack cost and constraints

v           =   [f;h];
