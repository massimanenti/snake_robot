function  [zsim, Xj, Yj, phi_ref_plot, u_plot] = snake_trajectory_with_reference(serpenoid_curve_parameters, z0, Ts, Tend, parameters, goal, kp, kd, animation_on, time_of_pause_animation)

% This function is called in the matlab file main_snake, and in particular
% in the post processing section. It is used to return some useful values
% for the plot as well as the entire evolution of the state. These useful
% variables are the coordinates of the joints (Xj and Yj) which are then
% used to perform the animation, and phi_ref_plot which is a variable
% containing the reference for the joint angles computed using the optimal
% serpenoid curve parameters. These reference phi are used afterwards in
% the plots.
% This function uses the model_simulation_snake function, model of the
% snake robot.

%% PARAMETERS EXTRACTION

n       = parameters.n;
n_vary  = parameters.n_vary;
last_interval = parameters.last_interval;
xgoal   = goal(1,1);
ygoal   = goal(2,1);
px0     = z0(n+1,1);
py0     = z0(n+2,1);
Nsim    = Tend/Ts;

%% ANIMATION SETUP       
% animation of the snake robot on the 2D plane.
% In the animation, the center of mass, position of the joints and links are plotted.

if animation_on == 1
    handle_fig = figure;
    hold all
    xlim([-10, 10])
    ylim([-10, 10])
    xlabel(' X [m] ')
    ylabel(' Y [m] ')
    title('Animation of snake robot')
    pbaspect([1 1 1])
    grid on
    for i = 1:n+1
        if i<=n
            animation(i) = plot([-1,1],[-1,1],'k','linewidth',1);
        end
        animation(i+n) = plot (1,1,'ok','markersize',5);
            
    end 
    animation(2*n+3) = plot (0,0,'ok','markersize',20, 'markerfacecolor','r');
end

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

mod_vel = zeros(1, Nsim);

u = zeros(n-1,1);
u_plot = zeros(n-1,Nsim);
zsim = zeros(2*n+4,Nsim);
zsim(:,1)=z0;

T = parameters.T;

%% SIMULATION

for ind_vary=1:n_vary    % ind_vary is the index to select the intervals
    
    alpha = serpenoid_curve_parameters(4*(ind_vary-1)+1,1);
    omega = serpenoid_curve_parameters(4*(ind_vary-1)+2,1);
    beta  = serpenoid_curve_parameters(4*(ind_vary-1)+3,1);
    gamma = serpenoid_curve_parameters(4*(ind_vary-1)+4,1);
    
    for t=T(1,ind_vary)+1:T(1,ind_vary+1) % t is the index to select the time istants

        for ind_joint=1:n-1    % ind_joint is the index to select the joints       

            % COMPUTATION OF PHI_REF and its derivatives at each time instant.
            % They are single column vectors.

            phi_ref(ind_joint,1)= alpha*sin(omega*Ts*(t-1)+(ind_joint-1)*beta)+gamma;
            phi_ref_dot(ind_joint,1)= alpha*omega*cos(omega*Ts*(t-1)+(ind_joint-1)*beta);
            phi_ref_ddot(ind_joint,1)= -alpha*omega^2*sin(omega*Ts*(t-1)+(ind_joint-1)*beta);
            
            phi_ref_plot(ind_joint,t) =  phi_ref(ind_joint,1);              % phi_ref used in the plot,
                                                                            % it must became a matrix in order to keep track of its time evolution

        end


        
        % SIMULATION    
        % now the vector of state in descrete time is computed with Euler Forward Method
        % Xj contains the x coordinates of the joints
        % Yj contains the y coordinates of the joints
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
        u_plot(:,t) = mult*u+add;
        
        % SNAKE ROBOT MOVEMENT ANIMATION
        if animation_on == 1
            if ishandle(handle_fig)
                hold on
                plot (xgoal,ygoal,'ok','markersize',5,'markerfacecolor','b');
                plot (px0,py0,'ok','markersize',5,'markerfacecolor','g');
                delete(animation)       % in order to delete the old position of the snake on the plot
                for i = 1:n+1           % for every i-th link of the snake robot, it plots a line connecting the extremities of the joints 
                    if i<(n+1)
                        animation(i) = plot([Xj(i,1),Xj(i+1,1)],[Yj(i,1),Yj(i+1,1)],'k','linewidth',2); 
                    end
                    animation(i+n) = plot (Xj(i,1),Yj(i,1),'ok','markersize',5,'markerfacecolor','w');                    
                end 
                animation(2*n+3) = plot (px(1,t),py(1,t),'ok','markersize',5, 'markerfacecolor','r');
                pause(time_of_pause_animation)
            end

        end
    end
end


