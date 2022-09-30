function  [zdot,Xj,Yj,mult,add] = model_simulation_snake(t,z,u,d,parameters)

% This function simulates the model of the snake robot.
%
% OUTPUT
% zdot, derivative of the state
% Xj, vector containing the x coordinates of the joints
% Yj, vector containing the y coordinates of the joints
%
% INPUT
% t, time (not needed since the model is time invariant)
% z, state vector
% u, input vector
% d, disturbances
% parameters, stuct containing the parameters of the model
%
% z, the state vector, has the following form:
%
% -  is a 2*n+4 elements column vector:
% -  z = [phi(1);       % first phi angle
%         phi(2);       ...
%          ...          ...
%         phi(n-1);     % n-1-th phi angle
%         theta(n);     % orientation of the n-th link
%         px;           % x coordinate of center of mass
%         py;           % y coordinate of center of mass
%         phi_dot(1);   % angular velocity related to the first phi angle
%         phi_dot(2);   ...
%          ...          ...
%         phi_dot(n-1); % angular velocity related to the n-1-th phi angle
%         theta_dot(n); % angular velocity of the n-th link
%         px_dot;       % x velocity of center of mass
%         py_dot]       % y velocity of center of mass
%

% The values contained in the struct 'parameters',
% which is an input of this function, is assigned to each
% parameter (n,l,m.. )
n = parameters.n;   % number of links of the robot [1]
l = parameters.l;   % half length of a link [m]
m = parameters.m;   % mass of a link [kg]
J = parameters.J;   % moment of inertia of a link [kg*m^2]
ct = parameters.ct; % tangential viscous friction coefficient (longitudinal direction of link i) [N*s/m]
cn = parameters.cn; % normal viscous friction coefficient (lateral direction of link i) [N*s/m]

% The vectors phi and phi_bar are created starting from the state vector z, input of the function. 
% According to the model derivation:
% phi vector contains phi(i), i=1..n-1 , which are the first n-1 elements of
% the state vector z;
% phi_bar vector contains phi(i), i=1,..n-1 and theta(n), which are the first n
% elements of the state vector z
phi = z(1:n-1,1);
phi_bar = z(1:n,1);

% Then, the vectors of derivatives phidot and phidot_bar
% are created. In particular, given the following structure for z:
% z = [phi(1);           (1)
%      phi(2);           (2)
%       ...
%      phi(n-1);         (n-1)
%      theta(n);         (n)

%      px;               (n+1)
%      py;               (n+2)

%      phi_dot(1);       (n+3)
%      phi_dot(2);       
%       ...
%      phi_dot(n-1);     (2*n+1)
%      theta_dot(n);     (2*n+2)
%      px_dot;           (2*n+3)
%      py_dot]           (2*n+4)
%
% it is clear that phidot can be created extracting elements from index n+3
% to 2*n+1, while phidot_bar from n+3 to 2*n+2

phidot = z(n+3:2*n+1,1);
phidot_bar = z(n+3:2*n+2,1);

% The same reasoning is used to create the vectors px, py and their
% derivatives, starting from z.
px=z(n+1,1);
py=z(n+2,1);
pxdot=z(2*n+3,1);
pydot=z(2*n+4,1);

% Next, the model of the snake robot must be computed. To do so, first
% the following matrices are defined:

% Matrix H
% The matrix H is used to perform a change of coordinates, from phi_bar to
% theta, so as to pass from absolute link angles to relative joint angles.
% This will be useful for performing partial feedback linearization

% H has the following form:
%   H = [1 1 1 1.. 1
%        0 1 1 1.. 1
%        .         .
%        .         .
%        .         .
%        0 0 0 0.. 0]   ,   square matrix (R^(n*n))

H = zeros(n);

for row=1:n
    for column=1:n
        if column >= row
            H(row,column)=1;  
        end 
    end    
end

% Vectors theta and thetadot
% Once H is computed, the change of coordinates is performed, obtaining
% theta, vector of n joint angles, and thetadot, vector of its derivatives
theta = H*phi_bar;             
thetadot=H*phidot_bar;

% Stheta and Ctheta
% First, the vectors sintheta, costheta are built.
% They are column vectors containing the cosine and sine of each element of
% theta. Then, the two diagonal matrices Stheta and Ctheta are built.
% Inside the 'for' cycle, also the vector that contains the square of the
% derivative of each theta(i) is computed.

sintheta=zeros(n,1);  
costheta=zeros(n,1);
thetadot_sq=zeros(n,1);

for ind=1:n
    sintheta(ind,1) = sin(theta(ind,1));
    costheta(ind,1) = cos(theta(ind,1));
    thetadot_sq(ind,1) = thetadot(ind,1)^2;       
end

Stheta = diag(sintheta);
Ctheta = diag(costheta);

% Martices A, D
% The matrices A and D represent, respectively, an addition and a difference matrix,
% which will be used, respectively, for adding and subtracting pairs of adjacent elements of a vector. 
% Furthermore, the vector e represents a summation vector, which
% will be used for adding all elements of an N-dimensional vector.

% A is a rectangular matrix (n-1xn) with all zeros except all 1 in the
% main diagonal and the upper diagonal  --> e.g. [1 1 0; 0 1 1]
A = zeros(n-1,n);
for ind=1:n-1
    A(ind,ind)=1;
    A(ind,ind+1)=1;
end

% D is a rectangular matrix (n-1xn) with all zeros except all 1 in the
% main diagonal and -1 in the upper diagonal  --> e.g. [1 -1 0; 0 1 -1]
D = zeros(n-1,n);
for ind=1:n-1
    D(ind,ind)=1;
    D(ind,ind+1)=-1;
end


e = ones(n,1);              % vector of ones, dimension n, needed for matrix E and other computations
E = [e zeros(n,1); 
     zeros(n,1) e];         % matrix of two columns and 2*n raws

% Matrices Cn and Ct
% these matrices are diagonal matrices. Each element on the diagonal is,
% correspondigly, the normal and tangential viscous friction coefficient
Cn = cn*eye(n);
Ct = ct*eye(n); 

V=A'*(D*D')^-1*A;
K=A'*(D*D')^-1*D;
W=m*l^2*Stheta*V*Ctheta-m*l^2*Ctheta*V*Stheta;                  
Mtheta = J*eye(n)+m*l^2*Stheta*V*Stheta+m*l^2*Ctheta*V*Ctheta;

% In order to compute the model of the robot, the following idea is
% considered:
% In order to write the model in a simpler form, partial feedback
% linearisation of underactuated systems is performed. This consists of
% linearising the dynamics corresponding to the actuated degrees of freedom
% of the system. To do so, first the model must be partitioned into two
% parts representing the actuated and unactuated degrees of freedom,
% respectively.

M_bar = [H'*Mtheta*H, zeros(n,2);zeros(2,n), n*m*eye(2)]; 
W_bar = [H'*W*diag(thetadot)*thetadot; zeros(2,1)];     
G_bar = [-l*H'*Stheta*K l*H'*Ctheta*K; -e' zeros(1,n); zeros(1,n) -e'];
B_bar = [eye(n-1); zeros(3,n-1)];                                

% the dynamics of the relative joint angles of the snake robot are the 
% actuated degrees of freedom of the snake robot, while the dynamics of 
% the absolute orientation and position of the snake robot are the 
% unactuated degrees of freedom. 
% The model may therefore be partitioned as:

M11 = M_bar(1:n-1,1:n-1);
M12 = M_bar(1:n-1,n:n+2);
M21 = M_bar(n:n+2,1:n-1);
M22 = M_bar(n:n+2,n:n+2);

W1 = W_bar(1:n-1,1);
W2 = W_bar(n:n+2,1);

G1 = G_bar(1:n-1, :);
G2 = G_bar(n:n+2, :);

% X and Y
% Vectors of positions and velocities of center of mass of each link are
% computed
X    = -l*K'*costheta+e*px;
Y    = -l*K'*sintheta+e*py;
Xdot = l*K'*Stheta*thetadot+e*pxdot;
Ydot = -l*K'*Ctheta*thetadot+e*pydot;

% Viscous friction forces 
% We define the viscous friction force on link i in the local link frame,
% f:    f(i) = - [ct    0
%               0   cn] * v(i),  
% where v(i) is the link velocity in the local link frame 
% Then, it must be expressed in the global frame.
% Doing it in a matrix form, the following is obtained:
frx= -[ct*Ctheta^2+cn*Stheta^2 (ct-cn)*Stheta*Ctheta]*[Xdot; Ydot];
fry= -[(ct-cn)*Stheta*Ctheta ct*Stheta^2+cn*Ctheta^2]*[Xdot; Ydot];
fr = [frx; fry];

% Based on the partitioned model, we are now ready to transform the model 
% of the snake robot to a simpler form through partial feedback linearisation 
% by introducing an input transformation which linearises the dynamics of 
% the actuated degrees of freedom.

% The idea is to obtain the following state equations: 
% zdot = [actuated dof
%         unactuated dof
%         derivatives of actuated dof
%         derivatives of unactuated dof]
%     
% where, calling them zdot1, zdot2, zdot3, zdot4, the following relation holds:
% 
% zdot3 = u
% zdot4 = A_italic+B_italic*u 

A_italic = -M22^-1*(W2+G2*fr);
B_italic = -M22^-1*M21;


%% OUTPUT COMPUTATION

% vector initialization
zdot = zeros(2*n+4,1);     
Xj = zeros(n+1,1);          
Yj = zeros(n+1,1);          

zdot(1:n-1,1)=z(n+3:2*n+1,1);       % z actuated: phi from 1 to n-1
zdot(n:n+2,1)=z(2*n+2:2*n+4,1);     % z unactuated: theta_n, px py
zdot(n+3:2*n+1,1)=u;                % u, corresponds to phi_ddot
zdot(2*n+2:2*n+4,1)= A_italic+B_italic*u;

Xj(1:n,1) = X - l*costheta;
Xj(n+1,1) = X(n,1) + l*costheta(n,1);

Yj(1:n,1) = Y - l*sintheta;
Yj(n+1,1) = Y(n,1) + l*sintheta(n,1);


mult = (M11 - M12*M22^-1*M21);
add = W1 + G1*fr -  M12*M22^-1*(W2 + G2*fr);












