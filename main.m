% Main Script
clc;
%clear all;
%close all;
global m
global I
global F M
global g

%% Load Aircraft Data from CSV and assigning them to variables
filename_density_L = 'ourflightcondation.csv'; 
aircraft_data = readmatrix(filename_density_L, 'Range', 'B2:B61'); 

% Time vector parameters
dt = aircraft_data(1);    
tfinal = aircraft_data(2); 
time_V = (0:dt:tfinal)';

% Initial conditions for states
s0 = aircraft_data(4:15);
sdot0 = zeros(12,1);  % Initial derivative of states
Vto = sqrt(s0(1)^2 + s0(2)^2 + s0(3)^2);  % Compute initial total velocity (Vto)
initialstatemodified=[s0(1:3);0;s0(4:6)];

% Control actions values
dc = [ aircraft_data(57:59) * pi/180 ; aircraft_data(60) ];  % Control inputs in radians

% Gravity, mass, and inertia values from aircraft data
m = aircraft_data(51);
g = aircraft_data(52);
Ixx = aircraft_data(53);
Iyy = aircraft_data(54);
Izz = aircraft_data(55);
Ixz = aircraft_data(56);
Ixy = 0;  Iyz = 0;

% Inertia matrix and inverse inertia matrix
I = [Ixx, -Ixy, -Ixz; -Ixy, Iyy, -Iyz; -Ixz, -Iyz, Izz];
invI = inv(I);

% Angle rad to deg and vice versa
D2R=pi/180;
R2D=180/pi;

% Stability derivatives for longitudinal motion
SD_Long = aircraft_data(21:36);
SD_Long_final = SD_Long;
templong  =num2cell(SD_Long);
[Xu,Zu,Mu,Xw,Zw,Mw,Zwd,Zq,Mwd,Mq,Xde,Zde,Mde,Xth,Zdth,Mdth] = deal(templong{:});
clear templong;

% Stability derivatives for lateral motion
SD_Lat_dash = aircraft_data(37:50);
SD_Lat = LateralFunc(SD_Lat_dash, Ixz, Izz, Ixx, Vto);
SD_Lat_final = SD_Lat;
templatrel = num2cell(SD_Lat);
[Yv,Yb,Lb,Nb,Lp,Np,Lr,Nr,Yda,Ydr,Lda,Nda,Ldr,Ndr] = deal(templatrel{:});
clear templatrel;

Lv=Lb/Vto;
Nv=Nb/Vto;
Yp=0; % not found in Nasa data
Yr=0; % not found in nasa data 

% Initial gravity force in body frame
mg0 = m * g * [sin(s0(8)); -cos(s0(8)) * sin(s0(7)); -cos(s0(8)) * cos(s0(7))];
initial_F_M = m * g * [sin(s0(8)); -cos(s0(8)) * sin(s0(7)); -cos(s0(8)) * cos(s0(7));0;0;0]; % Simulink Input
M_0 = [0;0;0];

%Initial velocities and rates
u_0 = s0(1);
v_0 = s0(2);
w_0 = s0(3);
p_0 = s0(4);
q_0 = s0(5);
r_0 = s0(6);
phi_0 = s0(7);
theta_0 = s0(8);
psi_0 = s0(9);
x_0 = s0(10);
y_0 = s0(11);
z_0 = s0(12);
wdot_0 = 0;
wdot = 0;
mass = m;

%% Matlab Nonlinear Model
for i = 1:16
    SD(i) = SD_Long_final(i);      
end
for j = 1:14
    SD(j+16) = SD_Lat_final(j);  % All stability derivatives in one array
end

n = 1000;
velocity_0=s0(1:3);
angular_v0=s0(4:6);
angles_0=s0(7:9);
position_0=s0(10:12);
Initial_state=[velocity_0;angular_v0;angles_0;position_0];

% solve and get the state and data  
[F,M]=F_M_Cal(SD,s0,Vto,dc,Initial_state,0);
[t, states] = RK4_solver(@getstates,Initial_state , 0, tfinal,n, SD,Vto, dc);

%% Simulink Nonlinear Model

State_Matrix = [ Xu  0   Xw  0   0   0   0;
                 0   Yv  0   Yp  0   Yr  0;
                 Zu  0   Zw  0   Zq  Lr  Zwd;
                 0   Lv  0   Lp  0   Lr  0;
                 Mu  0   Mw  0   Mq  0  Mwd;
                 0   Nv  0   Np  0   Nr 0];

Control_Matrix = [0 Xde Xth 0
                  Yda 0 0 Ydr 
                  0  Zde Zdth 0
                  Lda 0  0  Ldr
                  0  Mde Mdth 0 
                  Nda 0  0 Ndr];

% intial struc state 
so_struct.x = s0(10);
so_struct.y = s0(11);
so_struct.H = -1*s0(12);
so_struct.phi = s0(7);
so_struct.theta = s0(8);
so_struct.psi = s0(9);
so_struct.p = s0(4);
so_struct.q = s0(5);
so_struct.r = s0(6);
so_struct.u = s0(1);
so_struct.v = s0(2);
so_struct.w = s0(3);
so_struct.alpha = s0(8); 
so_struct.Beta = 0; 
so_struct.V_total = Vto; 
so_struct.W_dot = 0; 
so_struct.gamma=0;

% % run the model and get data
% modelName = 'Final_Task3_Simulink';
% load_system(modelName);
% set_param(modelName, 'StartTime', '0', 'StopTime', num2str(tfinal));
% simOut = sim(modelName, 'ReturnWorkspaceOutputs', 'on');
% 
% t2 = simOut.tout;   
% y2 = simOut.simout.Data;    
% y2 = y2';  

%% Linear Model 
x0=[u_0;w_0;q_0;theta_0]; %input
%% Full Longitudinal System
% Define constants
theta_0 = s0(8); wo = s0(3);u_0 = s0(1);wdot=0;

% Define constants and matrices for the longitudinal state-space model
u = [dc(3)* ones(size(t)); (dc(4)*ones(size(t)))];

A_longitudinal = [Xu, Xw,-w_0, -g*cos(theta_0);
                 Zu/(1-Zwd), Zw/(1-Zwd), (Zq + u_0)/(1-Zwd), -g*sin(theta_0)/(1-Zwd);
                  Mu + (Zu * Mwd) / (1 - Zwd), Mw + (Zw * Mwd) / (1 - Zwd), Mq + ((Zq + u_0) * Mwd) / (1 - Zwd), - (g * Mwd * sin(theta_0)) / (1 - Zwd);
                 0, 0, 1, 0];
B_longitudinal = [Xde, Xth;
                 Zde/(1-Zwd), Zdth/(1-Zwd);
                 Mde + Mwd*Zde/(1-Zwd), Mdth + Mwd*Zdth/(1-Zwd);
                 0, 0];
C_longitudinal = eye(4); % Observing all states
D_longitudinal = zeros(4, 2); % No direct feedthrough
sys_longitudinal = ss(A_longitudinal, B_longitudinal, C_longitudinal , D_longitudinal );
tf_full_longitudinal = tf(sys_longitudinal);

[y_full, t, x] = lsim(sys_longitudinal,u, t);

% Transfer function from elevator deflection (De) to each state output
tf_full_longitudinal_De_with_stateu = minreal(tf_full_longitudinal(1,1)); % u state
tf_full_longitudinal_De_with_statew = minreal(tf_full_longitudinal(2,1)); % w state
tf_full_longitudinal_De_with_stateq = minreal(tf_full_longitudinal(3,1)); % q state
tf_full_longitudinal_De_with_state_theta = minreal(tf_full_longitudinal(4,1)); % theta state
% Transfer function from throttle (Th) to each state output
tf_full_longitudinal_Th_with_stateu = minreal(tf_full_longitudinal(1,2)); % u state
tf_full_longitudinal_Th_with_statew = minreal(tf_full_longitudinal(2,2)); % w state
tf_full_longitudinal_Th_with_stateq = minreal(tf_full_longitudinal(3,2)); % q state
tf_full_longitudinal_Th_with_state_theta = minreal(tf_full_longitudinal(4,2)); % theta state

%% Long Mode For Longitudinal 
A_longitudinal_long = [Xu, -g*cos(theta_0); -Zu/(Zq+u_0), g*sin(theta_0)/(Zq+u_0)];
B_longitudinal_long = [Xde, Xth; -Zde/(Zq+u_0), -Zdth/(Zq+u_0)];

% Assuming the output directly maps to the state (no C and D matrices provided, use identity matrix for direct mapping)
C_longitudinal_long  = eye(2); % Direct mapping of state to output
D_longitudinal_long  = zeros(2, 2); % No direct feedthrough

% Define state-space model for the long-period longitudinal system
sys_longitudinal_long = ss(A_longitudinal_long, B_longitudinal_long, C_longitudinal_long, D_longitudinal_long);
tf_long_period = tf(sys_longitudinal_long);

% Run the simulation for the long-period mode with external input 'u'
[y_long, t, x] = lsim(sys_longitudinal_long, u, t);

% Transfer function from elevator deflection (De) to each state output
tf_long_period_De_with_state_u = minreal(tf_long_period(1,1)); % u state
tf_long_period_De_with_state_theta = minreal(tf_long_period(2,1)); % theta state

% Transfer function from throttle (Th) to each state output
tf_long_period_Th_with_state_u = minreal(tf_long_period(1,2)); % u state
tf_long_period_Th_with_state_theta = minreal(tf_long_period(2,2)); % theta state

%% Short Mode For Longitudinal 
A_longitudinal_short = [
    Zw / (1 - Zwd), (Zq + u_0) / (1 - Zwd);
    Mw + (Mwd * (Zw / (1 - Zwd))), Mq + (Mwd * ((Zq + u_0) / (1 - Zwd)))
];
B_longitudinal_short = [
    Zde / (1 - Zwd), Zdth / (1 - Zwd);
    Mde + (Mwd * Zde / (1 - Zwd)), Mdth + (Mwd*Zdth / (1 - Zwd))
];
C_longitudinal_short = eye(2); % Direct mapping of state to output
D_longitudinal_short = zeros(2, 2); % No direct feedthrough

% Define state-space model for the short-period longitudinal system
sys_longitudinal_short = ss(A_longitudinal_short, B_longitudinal_short, C_longitudinal_short, D_longitudinal_short);
tf_short_period = tf(sys_longitudinal_short);

% Initial state vector for short-period mode ()
x0_short = [w_0; q_0]; % Corrected to match the 2-dimensional system

% Run the simulation for the short-period mode
[y_short, t, x] = lsim(sys_longitudinal_short, u, t);

% Extract and simplify transfer functions with corresponding labels
% Transfer function from elevator deflection (De) to each state output
tf_short_period_De_with_state_w = minreal(tf_short_period(1,1)); % w state
tf_short_period_De_with_state_q = minreal(tf_short_period(2,1)); % q state

% Transfer function from throttle (Th) to each state output
tf_short_period_Th_with_state_w = minreal(tf_short_period(1,2)); % w state
tf_short_period_Th_with_state_q = minreal(tf_short_period(2,2)); % q state

%% Full Lateral Mode System
x0_lateral = [v_0; p_0; r_0; phi_0; psi_0];


A_lateral = [Yv, Yp+wo, Yr-u_0, g*cos(theta_0), 0;
            Lv, Lp, Lr, 0, 0;
            Nv, Np, Nr, 0, 0;
            0, 1, tan(theta_0), 0, 0;
            0, 0, 1/cos(theta_0), 0, 0];

B_lateral = [Yda, Ydr;
            Lda, Ldr;
            Nda, Ndr;
            0, 0;
            0, 0];

C_lateral = eye(5); % Observing all states
D_lateral = zeros(5, 2); % No direct feedthrough

u_lat = [dc(1)* ones(size(t)),1*ones(size(t))]';

% Create the state-space model for the lateral dynamics
sys_lateral = ss(A_lateral, B_lateral, C_lateral, D_lateral);
tf_full_lateral = tf(sys_lateral);

% Initial state vector for the lateral dynamics

% Run the simulation for the lateral dynamics with the specified input
[y_lateral, t, x_lateral] = lsim(sys_lateral, u_lat, t, x0_lateral);

% Transfer function from aileron deflection (Da) to each state output
tf_full_lateral_Da_with_state_v = minreal(tf_full_lateral(1,1));
tf_full_lateral_Da_with_state_p = minreal(tf_full_lateral(2,1));
tf_full_lateral_Da_with_state_r = minreal(tf_full_lateral(3,1)); 
tf_full_lateral_Da_with_state_phi = minreal(tf_full_lateral(4,1));
tf_full_lateral_Da_with_state_epsi = minreal(tf_full_lateral(5,1)); 

% Transfer function from rudder deflection (Dr) to each state output
tf_full_lateral_Dr_with_state_v = minreal(tf_full_lateral(1,2)); 
tf_full_lateral_Dr_with_state_p = minreal(tf_full_lateral(2,2)); 
tf_full_lateral_Dr_with_state_r = minreal(tf_full_lateral(3,2)); 
tf_full_lateral_Dr_with_state_phi = minreal(tf_full_lateral(4,2));
tf_full_lateral_Dr_with_state_epsi = minreal(tf_full_lateral(5,2)); 

%% 3-DOF Approximations for Spiral lateral Mode
A_lateral_3DOF_spiral = [Lp, Lr, 0;
                         Np, Nr, 0;
                         1, tan(theta_0), 0];
% B_lateral_3DOF_spiral = [Ldr;Ndr;0];
B_lateral_3DOF_spiral = [Lda, Ldr;
                         Nda, Ndr;
                         0, 0];

% Define the output matrix to observe all states (p, r, and v)
C_lateral_3DOF_spiral = eye(3);  % Observing all three states
D_lateral_3DOF_spiral = zeros(3, 2);  % No direct feedthrough from input to output

% Create the state-space model for the spiral mode approximation
sys_lateral_3DOF_spiral = ss(A_lateral_3DOF_spiral, B_lateral_3DOF_spiral, C_lateral_3DOF_spiral, D_lateral_3DOF_spiral);
tf_3DOF_spiral = tf(sys_lateral_3DOF_spiral);

% Initial state vector for the spiral mode with 3 DOF (initial conditions for p, r, and v)
x0_spiral_3DOF = [p_0; r_0; phi_0];  % Assuming zero initial conditions

% Run the simulation for the spiral mode with 3 DOF using the specified input
[y_spiral_3DOF, t, x_spiral_3DOF] = lsim(sys_lateral_3DOF_spiral, u_lat, t, x0_spiral_3DOF);

% Transfer functions from rudder deflection (Dr) to each state output
tf_3DOF_spiral_Dr_with_state_p = minreal(tf_3DOF_spiral(1,1)); % Roll Rate (p)
tf_3DOF_spiral_Dr_with_state_r = minreal(tf_3DOF_spiral(2,1)); % Yaw Rate (r)
tf_3DOF_spiral_Dr_with_state_phi = minreal(tf_3DOF_spiral(3,1)); % Roll Angle (phi)

%% 3-DOF Approximations for Dutch roll lateral Mode
A_lateral_3DOF_dutch = [Yv,  Yp+wo, Yr-u_0;
                        Lv, Lp, 0;
                        Nv, 0, Nr];
B_lateral_3DOF_dutch = [Yda, Ydr;
                        Lda, Ldr;
                        Nda, Ndr];
% Output matrix to observe all states (v, p, and r)
C_lateral_3DOF_dutch = eye(3);  % Observing all three states
D_lateral_3DOF_dutch = zeros(3, 2);  % No direct feedthrough from input to output

% Create the state-space model for the Dutch Roll mode approximation
sys_lateral_3DOF_dutch = ss(A_lateral_3DOF_dutch, B_lateral_3DOF_dutch, C_lateral_3DOF_dutch, D_lateral_3DOF_dutch);
tf_3DOF_dutch = tf(sys_lateral_3DOF_dutch);

% Initial state vector for the Dutch Roll mode with 3 DOF (initial conditions for v, p, and r)
x0_dutch_3DOF = [v_0; p_0; r_0];  % Assuming zero initial conditions

% Run the simulation for the Dutch Roll mode with 3 DOF using the specified input
[y_dutch_3DOF, t, x_dutch_3DOF] = lsim(sys_lateral_3DOF_dutch, u_lat, t, x0_dutch_3DOF);

% Transfer functions from aileron deflection (Da) to each state output
tf_3DOF_dutch_Da_with_state_v = minreal(tf_3DOF_dutch(1,1)); % Side velocity (v)
tf_3DOF_dutch_Da_with_state_p = minreal(tf_3DOF_dutch(2,1)); % Roll angle (phi)
tf_3DOF_dutch_Da_with_state_r = minreal(tf_3DOF_dutch(3,1)); % Yaw rate (r)

% Transfer functions from rudder deflection (Dr) to each state output
tf_3DOF_dutch_Dr_with_state_v = minreal(tf_3DOF_dutch(1,2));
tf_3DOF_dutch_Dr_with_state_p = minreal(tf_3DOF_dutch(2,2));
tf_3DOF_dutch_Dr_with_state_r = minreal(tf_3DOF_dutch(3,2));

%% 2-DOF Approximations for the Dutch roll lateral Mode
A_lateral_2DOF_dutch = [Yv, Yr-u_0;
                        Nv, Nr];
B_lateral_2DOF_dutch = [Yda, Ydr;
                        Nda, Ndr];

% Output matrix to observe both states (v and r)
C_lateral_2DOF_dutch = eye(2);  % Observing both states
D_lateral_2DOF_dutch = zeros(2, 2);  % No direct feedthrough from input to output

% Create the state-space model for the 2-DOF Dutch Roll mode
sys_lateral_2DOF_dutch = ss(A_lateral_2DOF_dutch, B_lateral_2DOF_dutch, C_lateral_2DOF_dutch, D_lateral_2DOF_dutch);
tf_2DOF_dutch = tf(sys_lateral_2DOF_dutch);

% Initial state vector for the 2-DOF Dutch Roll mode (initial conditions for v and r)
x0_dutch_2DOF = [v_0;r_0];  % Assuming zero initial conditions

% Run the simulation for the Dutch Roll mode with 2 DOF using the specified input
[y_dutch_2DOF, t, x_dutch_2DOF] = lsim(sys_lateral_2DOF_dutch, u_lat, t, x0_dutch_2DOF);

% Transfer functions from aileron deflection (Da) to each state output
tf_2DOF_dutch_Da_with_state_v = minreal(tf_2DOF_dutch(1,1)); % Side velocity (v)
tf_2DOF_dutch_Da_with_state_r = minreal(tf_2DOF_dutch(2,1)); % Yaw rate (r)

% Transfer functions from rudder deflection (Dr) to each state output
tf_2DOF_dutch_Dr_with_state_v = minreal(tf_2DOF_dutch(1,2));
tf_2DOF_dutch_Dr_with_state_r = minreal(tf_2DOF_dutch(2,2));

%% 1-DOF Approximations for the Roll lateral Mode
A_lateral_1DOF = Lp; % Only the roll damping term Lp
B_lateral_1DOF = Lda; % Aileron control input affecting roll rate
C_lateral_1DOF = 1; % Observe roll rate (p) directly
D_lateral_1DOF = 0; % No direct feedthrough

sys_lateral_1DOF = ss(A_lateral_1DOF, B_lateral_1DOF, C_lateral_1DOF, D_lateral_1DOF);
% Calculate the transfer function from aileron deflection (Da) to roll rate (p)
tf_1DOF_roll_Da_with_state_p = minreal(tf(sys_lateral_1DOF));

% Define the input signal (u), representing aileron deflection
u_1dof = ones(size(t));  % Constant aileron deflection for demonstration

% Initial state vector for roll rate
x0_roll_1DOF = p_0;  % Assuming zero initial roll rate

% Run the simulation for the roll mode with the specified input
[y_roll_1DOF, t, x_roll_1DOF] = lsim(sys_lateral_1DOF, u_1dof, t, x0_roll_1DOF);


%% Longitudinal Controller Design 

% Transfer function from elevator deflection (De) to each state output
De_u = minreal(tf_full_longitudinal(1,1)); % u state
De_w = minreal(tf_full_longitudinal(2,1)); % w state
De_q = minreal(tf_full_longitudinal(3,1)); % q state
De_theta = minreal(tf_full_longitudinal(4,1)); % theta state

% Transfer function from throttle (Th) to each state output
Th_u = minreal(tf_full_longitudinal(1,2)); % u state
Th_w = minreal(tf_full_longitudinal(2,2)); % w state
Th_q = minreal(tf_full_longitudinal(3,2)); % q state
Th_theta = minreal(tf_full_longitudinal(4,2)); % theta state
u_dT= tf_full_longitudinal(1,2);

% Transfer function used in building the controller
servo= tf(10,[1 10]);
integrator=tf(1,[1 0]);
diffrentiator=tf([1 0],1);
engine_timelag=tf(0.1,[1 0 1]);

De_Wd=De_w*diffrentiator;
De_az=De_Wd-u_0*De_q;
Theta_az=De_az/De_theta;
De_alpha=De_w/u_0;
De_gamma=De_theta-De_alpha;

%pitch control
ol_theta_thetacomm=-servo*De_theta;

% Velocity Control
OL_u_ucom=servo*engine_timelag*u_dT;

load('Velocity_Controller_tf.mat')
load('eltaw2am_c.mat')
load('Pitch_Design2.mat')
load('Altitude_Design.mat')

%% Lateral Controller Design 

% Transfer function from aileron deflection 
v_da=minreal(tf_full_lateral(1,1));
p_da=minreal(tf_full_lateral(2,1));
r_da=minreal(tf_full_lateral(3,1));
phi_da=minreal(tf_full_lateral(4,1));
psi_da=minreal(tf_full_lateral(5,1));
beta_da=v_da/u_0;

% Transfer function from rudder deflection 
v_dr=minreal(tf_full_lateral(1,2));
p_dr=minreal(tf_full_lateral(2,2));
r_dr=minreal(tf_full_lateral(3,2));
phi_dr=minreal(tf_full_lateral(4,2));
psi_dr=minreal(tf_full_lateral(5,2));
beta_dr=v_da/u_0;

% yaw controller
ol_r_rcomm=servo*r_dr;
load('Yaw_Damper.mat')

Lat_YawDamper_ss_series=feedback(series(append(1,servo),sys_lateral,[1 2],[1 2]),HPF_YawDamper,2,3,1);
Lat_YawDamper_tf_series=tf(Lat_YawDamper_ss_series);
phi_da_AfterYawDamper=Lat_YawDamper_tf_series(4,1);
% phi controller
ol_phi_phicomm=minreal(servo*phi_da_AfterYawDamper);


% %% Ploting
% %% MATLAB Plotting
% rad_to_deg = 180 / pi;  % Conversion from radians to degrees
% 
% % Plot results (e.g., velocities, angles, positions)
% figure;
% sgtitle('Matlab Code Plots');
% 
% % Plot forward, side, and vertical velocities (u, v, w)
% subplot(4,3,1);
% plot(t, states(1,:), 'b');  % u (forward velocity)
% title('Forward Velocity (u)');
% xlabel('Time (s)');
% ylabel('Velocity (ft/s)');
% 
% % subplot(4,3,2);
% % plot(t2, y2(2,:), 'r');  % v (side velocity)
% % title('Side Velocity (v)');
% % xlabel('Time (s)');
% % ylabel('Velocity (ft/s)');
% % 
% % subplot(4,3,3);
% % plot(t2, y2(3,:), 'g');  % w (vertical velocity)
% % title('Vertical Velocity (w)');
% % xlabel('Time (s)');
% % ylabel('Velocity (ft/s)');
% 
% subplot(4,3,2);
% plot(t, asin(states(2,:)/Vto)*rad_to_deg, 'r');  % bita (side velocity)
% title('bita');
% xlabel('Time (s)');
% ylabel('Velocity (ft/s)');
% 
% subplot(4,3,3);
% plot(t, atan(states(3,:)./states(1,:))*rad_to_deg, 'g');  % alpha (vertical velocity)
% title('Alpa');
% xlabel('Time (s)');
% ylabel('Velocity (ft/s)');
% 
% % Angular velocities (p, q, r)
% 
% subplot(4,3,4);
% plot(t, states(4,:), 'b');  % p (roll rate)
% title('Roll Rate (p)');
% xlabel('Time (s)');
% ylabel('Angular Velocity (rad/s)');
% 
% subplot(4,3,5);
% plot(t, states(5,:), 'r');  % q (pitch rate)
% title('Pitch Rate (q)');
% xlabel('Time (s)');
% ylabel('Angular Velocity (rad/s)');
% 
% subplot(4,3,6);
% plot(t, states(6,:), 'g');  % r (yaw rate)
% title('Yaw Rate (r)');
% xlabel('Time (s)');
% ylabel('Angular Velocity (rad/s)');
% 
% % Euler angles (phi, theta, psi) in degrees
% 
% subplot(4,3,7);
% plot(t, states(7,:) * rad_to_deg, 'b');  % phi (roll angle in degrees)
% title('Roll Angle (phi) in degrees');
% xlabel('Time (s)');
% ylabel('Angle (deg)');
% 
% subplot(4,3,8);
% plot(t, states(8,:) * rad_to_deg, 'r');  % theta (pitch angle in degrees)
% title('Pitch Angle (theta) in degrees');
% xlabel('Time (s)');
% ylabel('Angle (deg)');
% 
% subplot(4,3,9);
% plot(t, states(9,:) * rad_to_deg, 'g');  % psi (yaw angle in degrees)
% title('Yaw Angle (psi) in degrees');
% xlabel('Time (s)');
% ylabel('Angle (deg)');
% 
% % Position (x, y, z) in feet
% 
% subplot(4,3,10);
% plot(t, states(10,:), 'b');  
% title('Position X (ft)');
% xlabel('Time (s)');
% ylabel('Position (ft)');
% 
% subplot(4,3,11);
% plot(t, states(11,:), 'r'); 
% title('Position Y (ft)');
% xlabel('Time (s)');
% ylabel('Position (ft)');
% 
% subplot(4,3,12);
% plot(t, round(states(12,:),3), 'g');  
% title('Position Z (ft)');
% xlabel('Time (s)');
% ylabel('Position (ft)');
% 
% 
% % Trajectory (x vs y vs z in feet)
% figure;
% plot3(states(10,:), states(11,:),  round(states(12,:),3), 'b');  % 3D trajectory plot
% title('Matlab Code Trajectory (ft)', 'FontWeight', 'bold', 'FontSize', 14);
% xlabel('X (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% ylabel('Y (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% zlabel('Z (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% grid on;
% 
% 
% %% Simulink Plot
% % Plot results (e.g., velocities, angles, positions)
% figure;
% sgtitle('Simulink Code Plots');
% % Plot forward, side, and vertical velocities (u, v, w)
% subplot(4,3,1);
% plot(t2, y2(1,:), 'b');  % u (forward velocity)
% title('Forward Velocity (u)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Velocity (ft/s)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% % subplot(4,3,2);
% % plot(t2, y2(2,:), 'r');  % v (side velocity)
% % title('Side Velocity (v)', 'FontWeight', 'bold', 'FontSize', 12);
% % xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% % ylabel('Velocity (ft/s)', 'FontWeight', 'bold', 'FontSize', 10);
% % 
% % subplot(4,3,3);
% % plot(t2, y2(3,:), 'g');  % w (vertical velocity)
% % title('Vertical Velocity (w)', 'FontWeight', 'bold', 'FontSize', 12);
% % xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% % ylabel('Velocity (ft/s)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% 
% % plots of alpha and beta
% subplot(4,3,2);
% plot(t2, y2(13,:)*rad_to_deg, 'r');  % v (side angle beta)
% title('Side angle (Beta)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Beta (deg)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% subplot(4,3,3);
% plot(t2, y2(14,:)*rad_to_deg, 'g');  % w (alpha )
% title('angle of attach (Alpha)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('alpha (deg)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% 
% % Angular velocities (p, q, r)
% subplot(4,3,4);
% plot(t2, y2(4,:)* rad_to_deg, 'b');  % p (roll rate)
% title('Roll Rate (p)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Angular Velocity (rad/s)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% subplot(4,3,5);
% plot(t2, y2(5,:)*rad_to_deg, 'r');  % q (pitch rate)
% title('Pitch Rate (q)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Angular Velocity (rad/s)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% subplot(4,3,6);
% plot(t2, y2(6,:)*rad_to_deg, 'g');  % r (yaw rate)
% title('Yaw Rate (r)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Angular Velocity (rad/s)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% % Euler angles (phi, theta, psi) in degrees
% subplot(4,3,7);
% plot(t2, y2(7,:) * rad_to_deg, 'b');  % phi (roll angle in degrees)
% title('Roll Angle (phi) in degrees', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Angle (deg)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% subplot(4,3,8);
% plot(t2, y2(8,:) * rad_to_deg, 'r');  % theta (pitch angle in degrees)
% title('Pitch Angle (theta) in degrees', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Angle (deg)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% subplot(4,3,9);
% plot(t2, unwrap(y2(9,:)) * rad_to_deg, 'g');  % psi (yaw angle in degrees)
% title('Yaw Angle (psi) in degrees', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Angle (deg)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% % Position (x, y, z) in feet
% subplot(4,3,10);
% plot(t2, y2(10,:), 'b');  
% title('Position X (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Position (ft)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% subplot(4,3,11);
% plot(t2, y2(11,:), 'r'); 
% title('Position Y (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Position (ft)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% subplot(4,3,12);
% plot(t2, round(y2(12,:), 3), 'g');  
% title('Position Z (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% xlabel('Time (s)', 'FontWeight', 'bold', 'FontSize', 10);
% ylabel('Position (ft)', 'FontWeight', 'bold', 'FontSize', 10);
% 
% % Trajectory (x vs y vs z in feet)
% figure;
% plot3(y2(10,:), y2(11,:), round(y2(12,:), 3), 'b');  % 3D trajectory plot
% title('Simulink Code Trajectory (ft)', 'FontWeight', 'bold', 'FontSize', 14);
% xlabel('X (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% ylabel('Y (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% zlabel('Z (ft)', 'FontWeight', 'bold', 'FontSize', 12);
% grid on;
% 
% %% Step response plots
% 
% %Longitudinal
% 
% % Compare 'u' (velocity) from full and long-period models
% subplot(4,1,1); 
% plot(t, y_full(:,1)+x0(1), 'r', 'DisplayName', 'Fully linearized Model'); hold on;
% plot(t, y_long(:,1)+x0(1), 'b--', 'DisplayName', 'Long-Period Model');
% plot(t2, y2(1,:), 'g--', 'DisplayName', 'Non-linear Model');
% title('Response of u (velocity)'); xlabel('Time (s)'); ylabel('u (m/s)');
% legend;
% hold off;
% 
% % Compare 'w' (vertical speed) from full and short-period models
% subplot(4,1,2); 
% plot(t, y_full(:,2)+x0(2), 'r', 'DisplayName', 'Fully linearized Model'); hold on;
% plot(t, y_short(:,1)+x0(2), 'b--', 'DisplayName', 'Short-Period Model');
% plot(t, states(3,:), 'g--', 'DisplayName', 'Non-Linear Model');
% title('Response of w (vertical speed)'); xlabel('Time (s)'); ylabel('w (m/s)');
% legend;
% hold off;
% 
% % Compare 'q' (pitch rate) from full and short-period models
% subplot(4,1,3); 
% plot(t, y_full(:,3)+x0(3), 'r', 'DisplayName', 'Fully linearized Model'); hold on;
% plot(t, y_short(:,2)+x0(3), 'b--', 'DisplayName', 'Short-Period Model');
% plot(t2, y2(5,:), 'g--', 'DisplayName', 'Non-Linear Model');
% title('Response of q (pitch rate)'); xlabel('Time (s)'); ylabel('q (rad/s)');
% legend;
% hold off;
% 
% % Compare 'theta' (pitch angle) from full and long-period models
% subplot(4,1,4); 
% plot(t,y_full(:,4)+x0(4), 'r', 'DisplayName', 'Fully linearized Model'); hold on;
% plot(t, y_long(:,2)+x0(4), 'b--', 'DisplayName', 'Long-Period Model');
% plot(t2, y2(8,:), 'g--', 'DisplayName', 'Non-Linear Model');
% title('Response of \theta (pitch angle)'); xlabel('Time (s)'); ylabel('\theta (rad)');
% legend;
% hold off;
% 
% 
% 
% %% Ploting Rotloucs and Bodeplot For Longitudinal 
% % Call the function for each transfer function comparison between full and short modes
% 
% % Full longitudinal vs short-period
% plotAndSaveSubplots(tf_full_longitudinal_De_with_statew, tf_short_period_De_with_state_w, 'Short', 'De', 'w');
% plotAndSaveSubplots(tf_full_longitudinal_De_with_stateq, tf_short_period_De_with_state_q, 'Short', 'De', 'q');
% plotAndSaveSubplots(tf_full_longitudinal_Th_with_statew, tf_short_period_Th_with_state_w, 'Short', 'Th', 'w');
% plotAndSaveSubplots(tf_full_longitudinal_Th_with_stateq, tf_short_period_Th_with_state_q, 'Short', 'Th', 'q');
% 
% % Full longitudinal vs long-period
% plotAndSaveSubplots(tf_full_longitudinal_De_with_stateu, tf_long_period_De_with_state_u, 'Long', 'De', 'u');
% plotAndSaveSubplots(tf_full_longitudinal_De_with_state_theta, tf_long_period_De_with_state_theta, 'Long', 'De', 'theta');
% plotAndSaveSubplots(tf_full_longitudinal_Th_with_stateu, tf_long_period_Th_with_state_u, 'Long', 'Th', 'u');
% plotAndSaveSubplots(tf_full_longitudinal_Th_with_state_theta, tf_long_period_Th_with_state_theta, 'Long', 'Th', 'theta');
% 
% 
% %% Plot and save Root Locus and Bode plots for each mode
% % Define the transfer functions and labels for comparisons
% 
% % Aileron Deflection (Da) Comparisons
% Da_v = {tf_full_lateral_Da_with_state_v, tf_3DOF_dutch_Da_with_state_v, tf_2DOF_dutch_Da_with_state_v};
% Da_p = {tf_full_lateral_Da_with_state_p, tf_3DOF_dutch_Da_with_state_p,tf_1DOF_roll_Da_with_state_p};
% Da_r = {tf_full_lateral_Da_with_state_r, tf_3DOF_dutch_Da_with_state_r, tf_2DOF_dutch_Da_with_state_r};
% Da_phi = {tf_full_lateral_Da_with_state_phi};
% Da_epsi = {tf_full_lateral_Da_with_state_epsi};
% 
% % Rudder Deflection (Dr) Comparisons
% Dr_v = {tf_full_lateral_Dr_with_state_v, tf_3DOF_dutch_Dr_with_state_v, tf_2DOF_dutch_Dr_with_state_v};
% Dr_p = {tf_full_lateral_Dr_with_state_p, tf_3DOF_spiral_Dr_with_state_p,tf_3DOF_dutch_Dr_with_state_p};
% Dr_r = {tf_full_lateral_Dr_with_state_r, tf_3DOF_spiral_Dr_with_state_r,tf_3DOF_dutch_Dr_with_state_r, tf_2DOF_dutch_Dr_with_state_r};
% Dr_phi = {tf_full_lateral_Dr_with_state_phi,tf_3DOF_spiral_Dr_with_state_phi};
% Dr_epsi = {tf_full_lateral_Dr_with_state_epsi};
% 
% % Mode labels for Aileron Deflection (Da) Comparisons
% mode_labels_Da_v = {'Full Mode', '3-DOF Dutch Roll', '2-DOF Dutch Roll'};
% mode_labels_Da_p = {'Full Mode','3-DOF Dutch Roll','1-DOF Roll'};
% mode_labels_Da_r = {'Full Mode', '3-DOF Dutch Roll', '2-DOF Dutch Roll'};
% mode_labels_Da_phi = {'Full Mode'};
% mode_labels_Da_epsi = {'Full Mode'};
% 
% % Mode labels for Rudder Deflection (Dr) Comparisons
% mode_labels_Dr_v = {'Full Mode', '3-DOF Dutch Roll', '2-DOF Dutch Roll'};
% mode_labels_Dr_p = {'Full Mode', '3-DOF Spiral', '3-DOF Dutch Roll'};
% mode_labels_Dr_r = {'Full Mode', '3-DOF Spiral', '3-DOF Dutch Roll', '2-DOF Dutch Roll'};
% mode_labels_Dr_phi = {'Full Mode', '3-DOF Spiral'};
% mode_labels_Dr_epsi = {'Full Mode'};
% 
% %Call the function to compare each mode for specific input-output pairs
% 
% % Aileron Deflection (Da) Comparisons
% plotAndCompareModesWithSubplots(Da_v, mode_labels_Da_v, 'Da', 'v');
% plotAndCompareModesWithSubplots(Da_p, mode_labels_Da_p, 'Da', 'p');
% plotAndCompareModesWithSubplots(Da_r, mode_labels_Da_r, 'Da', 'r');
% plotAndCompareModesWithSubplots(Da_phi, mode_labels_Da_phi, 'Da', 'phi');
% plotAndCompareModesWithSubplots(Da_epsi, mode_labels_Da_epsi, 'Da', 'epsi');
% 
% % Rudder Deflection (Dr) Comparisons
% plotAndCompareModesWithSubplots(Dr_v, mode_labels_Dr_v, 'Dr', 'v');
% plotAndCompareModesWithSubplots(Dr_p, mode_labels_Dr_p, 'Dr', 'p');
% plotAndCompareModesWithSubplots(Dr_r, mode_labels_Dr_r, 'Dr', 'r');
% plotAndCompareModesWithSubplots(Dr_phi, mode_labels_Dr_phi, 'Dr', 'phi');
% plotAndCompareModesWithSubplots(Dr_epsi, mode_labels_Dr_epsi, 'Dr', 'epsi');
% 
% %% Display Tranfer Function For Longitudinal 
% % Call the function below for each input-output pair comparison
% 
% % Full longitudinal vs short-period
% displayTransferFunctionAsFraction(tf_full_longitudinal_De_with_statew, tf_short_period_De_with_state_w, 'Short', 'De', 'w');
% displayTransferFunctionAsFraction(tf_full_longitudinal_De_with_stateq, tf_short_period_De_with_state_q, 'Short', 'De', 'q');
% displayTransferFunctionAsFraction(tf_full_longitudinal_Th_with_statew, tf_short_period_De_with_state_w, 'Short', 'Th', 'w');
% displayTransferFunctionAsFraction(tf_full_longitudinal_Th_with_stateq, tf_short_period_De_with_state_q, 'Short', 'Th', 'q');
% 
% % Full longitudinal vs long-period
% displayTransferFunctionAsFraction(tf_full_longitudinal_De_with_stateu, tf_long_period_De_with_state_u, 'Long', 'De', 'u');
% displayTransferFunctionAsFraction(tf_full_longitudinal_De_with_state_theta, tf_long_period_De_with_state_theta, 'Long', 'De', 'theta');
% displayTransferFunctionAsFraction(tf_full_longitudinal_Th_with_stateu, tf_long_period_Th_with_state_u, 'Long', 'Th', 'u');
% displayTransferFunctionAsFraction(tf_full_longitudinal_Th_with_state_theta, tf_long_period_Th_with_state_theta, 'Long', 'Th', 'theta');
% 
% %% Display the transfer functions in fraction format
% % Call the function below for each group of transfer functions
% 
% % Aileron Deflection (Da) Comparisons
% displayTransferFunctionGroup(Da_v, mode_labels_Da_v, 'Da', 'v');
% displayTransferFunctionGroup(Da_p, mode_labels_Da_p, 'Da', 'p');
% displayTransferFunctionGroup(Da_r, mode_labels_Da_r, 'Da', 'r');
% displayTransferFunctionGroup(Da_phi, mode_labels_Da_phi, 'Da', 'phi');
% displayTransferFunctionGroup(Da_epsi, mode_labels_Da_epsi, 'Da', 'epsi');
% 
% % Rudder Deflection (Dr) Comparisons
% displayTransferFunctionGroup(Dr_v, mode_labels_Dr_v, 'Dr', 'v');
% displayTransferFunctionGroup(Dr_p, mode_labels_Dr_p, 'Dr', 'p');
% displayTransferFunctionGroup(Dr_r, mode_labels_Dr_r, 'Dr', 'r');
% displayTransferFunctionGroup(Dr_phi, mode_labels_Dr_phi, 'Dr', 'phi');
% displayTransferFunctionGroup(Dr_epsi, mode_labels_Dr_epsi, 'Dr', 'epsi');
% 
% %% Functions
% 
% % Define a function to display the transfer function expressions as a large fraction
% function displayTransferFunctionAsFraction(tf_full, tf_mode, mode_name, input_name, output_name)
%     % Extract numerator and denominator for the full mode transfer function
%     [num_full, den_full] = tfdata(tf_full, 'v'); % Get numerator and denominator as vectors
%     % Convert numerator and denominator to strings
%     num_full_str = poly2str(num_full, 's');
%     den_full_str = poly2str(den_full, 's');
% 
%     % Extract numerator and denominator for the mode transfer function (short/long)
%     [num_mode, den_mode] = tfdata(tf_mode, 'v'); % Get numerator and denominator as vectors
%     % Convert numerator and denominator to strings
%     num_mode_str = poly2str(num_mode, 's');
%     den_mode_str = poly2str(den_mode, 's');
% 
%     % Display the formatted transfer functions as a large fraction
%     fprintf('Comparison for Transfer Function (%s to %s):\n', input_name, output_name);
%     fprintf('------------------------------------------------------\n');
% 
%     % Display the full mode transfer function in S domain format as a large fraction
%     fprintf('Full Mode Transfer Function:\n');
%     fprintf('    %s\n', num_full_str);  % Numerator
%     fprintf('    %s\n', repmat('-', max(length(num_full_str), length(den_full_str)), 1)); % Divider line
%     fprintf('    %s\n\n', den_full_str); % Denominator
% 
%     % Display the alternative mode (short or long) transfer function as a large fraction
%     fprintf('%s Mode Transfer Function:\n', mode_name);
%     fprintf('    %s\n', num_mode_str);  % Numerator
%     fprintf('    %s\n', repmat('-', max(length(num_mode_str), length(den_mode_str)), 1)); % Divider line
%     fprintf('    %s\n\n', den_mode_str); % Denominator
% 
%     % Add a line break for readability between comparisons
%     fprintf('\n\n');
% end
% 
% % Define function to display the transfer functions for a group in fraction format
% function displayTransferFunctionGroup(tf_list, mode_labels, input_name, output_name)
%     num_modes = length(tf_list);
%     fprintf('Transfer Function Comparisons for %s to %s:\n', input_name, output_name);
%     fprintf('------------------------------------------------------\n');
% 
%     for k = 1:num_modes
%         % Extract numerator and denominator
%         [num, den] = tfdata(tf_list{k}, 'v');
%         num_str = poly2str(num, 's');
%         den_str = poly2str(den, 's');
% 
%         % Display transfer function for each mode
%         fprintf('%s Mode Transfer Function:\n', mode_labels{k});
%         fprintf('    %s\n', num_str);  % Numerator
%         fprintf('    %s\n', repmat('-', max(length(num_str), length(den_str)), 1)); % Divider line
%         fprintf('    %s\n\n', den_str); % Denominator
%     end
%     fprintf('\n\n'); % Add extra line spacing between groups for readability
% end
% 
% % Define function to plot and save Root Locus and Bode plots with subplots for each mode comparison
% function plotAndCompareModesWithSubplots(tf_list, mode_labels, input_name, output_name)
%     % Ensure 'graphs' directory exists
%     if ~exist('graphs', 'dir')
%         mkdir('graphs');
%     end
% 
%     num_modes = length(tf_list);
% 
%     % Plot Root Locus with each mode in its own subplot
%     figure;
%     for k = 1:num_modes
%         subplot(1, num_modes, k); % Create a subplot for each mode
%         rlocus(tf_list{k});
%         title([mode_labels{k}]);
%         xlabel('Real Axis');
%         ylabel('Imaginary Axis');
%     end
%     sgtitle(['Root Locus - ', input_name, ' to ', output_name]);
%     saveas(gcf, fullfile('graphs', ['RootLocus_Comparison_', input_name, '_to_', output_name, '.png']));
% 
%     % Plot Bode Plot with each mode in the same plot for overlay comparison
%     figure;
%     hold on;
%     for k = 1:num_modes
%         bode(tf_list{k});
%     end
%     hold off;
%     title(['Bode Plot Comparison - ', input_name, ' to ', output_name]);
%     legend(mode_labels);
%     saveas(gcf, fullfile('graphs', ['BodePlot_Comparison_', input_name, '_to_', output_name, '.png']));
% end
% 
% % Define a function to plot and save root locus and Bode plots with subplots for each transfer function
% function plotAndSaveSubplots(tf_full, tf_short, mode_name, input_name, output_name)
%     % Ensure 'graphs' directory exists
%     if ~exist('graphs', 'dir')
%         mkdir('graphs');
%     end
% 
%     % Create figure for Root Locus with two subplots
%     figure;
% 
%     % Full Mode Root Locus in the first subplot
%     subplot(1, 2, 1);
%     rlocus(tf_full);
%     title(['Root Locus - Full Mode (', input_name, ' to ', output_name, ')']);
%     xlabel('Real Axis');
%     ylabel('Imaginary Axis');
% 
%     % Short Mode Root Locus in the second subplot
%     subplot(1, 2, 2);
%     rlocus(tf_short);
%     title(['Root Locus - ', mode_name, ' Mode (', input_name, ' to ', output_name, ')']);
%     xlabel('Real Axis');
%     ylabel('Imaginary Axis');
% 
%     % Save the combined root locus figure
%     saveas(gcf, fullfile('graphs', ['RootLocus_', input_name, '_to_', output_name, '_Full_vs_', mode_name, '.png']));
% 
%     % Bode plot (Separate for Full and Short Mode)
%     figure;
%     bode(tf_full);
%     hold on;
%     bode(tf_short);
%     hold off;
%     title(['Bode Plot - ', input_name, ' to ', output_name, ' (Full vs ', mode_name, ' Mode)']);
%     legend('Full Mode', [mode_name, ' Mode']);
%     saveas(gcf, fullfile('graphs', ['BodePlot_', input_name, '_to_', output_name, '_', mode_name, '.png']));
% end
airport_name ="Hilo International (PHTO)";


% %% Test_1
% % 
% % Get the folder of the current MATLAB script
% current_folder = fileparts(mfilename('fullpath'));
% output_folder = fullfile(current_folder, 'Task7_result_graphs'); % Define the folder as 'Task7_result_graphs'
% 
% % Create the folder if it does not exist
% if ~exist(output_folder, 'dir')
%     mkdir(output_folder);
% end
% 
% % Figure 1: Pitching Angle Output for Pitch Controller
% figure(1)
% plot(out.tout, out.theta_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.theta_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.theta_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$Pitching$ $Angle$ $Response$ $for$ $Pitch$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\theta$ $command$', '$\theta$ $non$ $linear$', '$\theta$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\theta$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'test1_Pitching_Angle_Response_for_Pitch_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Figure 2: U Output for Pitch Controller
% figure(2)
% plot(out.tout, out.u_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.u_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.u_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$U$ $Response$ $for$ $Pitch$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$U$ $command$', '$U$ $non$ $linear$', '$U$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); 
% ylabel('$U$ $(ft/s)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'test1_U_Response_for_Pitch_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Figure 3: Flight Path Angle for Pitch Controller
% figure(3)
% plot(out.tout, out.gamma_states.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.gamma_states.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\gamma$ $Response$ $for$ $Pitch$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\gamma$ $non$ $linear$', '$\gamma$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\gamma$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'test1_Flight_Path_Angle_for_Pitch_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Figure 4: H Output for Pitch Controller
% figure(4)
% plot(out.tout, out.H_states.Data(:, 2), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.H_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$H$ $Response$ $for$ $Pitch$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$H$ $non$ $linear$', '$H$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$H$ $(ft)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'test1_H_Response_for_Pitch_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% % Plot δe (elevator control action)
% figure(5)
% plot(out.tout, out.elevator_controller.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.elevator_controller.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\delta_e$ $for$ $Pitch$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_e$ $non$ $linear$', '$\delta_e$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_e$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test1_delta_e_for_Pitch_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);


% 
% % %% Test 2
%%
% % Plot θ (pitch angle)
% figure(6)
% plot(out.tout, out.theta_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.theta_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.theta_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$\theta$ $Response$ $for$ $Pitch$ $\&$ $Velocity$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\theta$ $command$', '$\theta$ $non$ $linear$', '$\theta$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\theta$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test2_theta_Response_for_Pitch_Velocity_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot u (forward velocity)
% figure(7)
% plot(out.tout, out.u_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.u_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.u_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$U$ $Response$ $for$ $Pitch$ $\&$ $Velocity$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$U$ $command$', '$U$ $non$ $linear$', '$U$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$U$ $(ft/s)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test2_u_Response_for_Pitch_Velocity_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot γ (flight path angle)
% figure(8)
% plot(out.tout, out.gamma_states.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.gamma_states.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\gamma$ $Response$ $for$ $Pitch$ $\&$ $Velocity$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\gamma$ $non$ $linear$', '$\gamma$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\gamma$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test2_gamma_Response_for_Pitch_Velocity_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot h (altitude)
% figure(9)
% %plot(out.tout, out.H_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.H_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on 
% plot(out.tout, out.H_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$H$ $Response$ $for$ $Pitch$ $\&$ $Velocity$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$H$ $non$ $linear$', '$H$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$H$ $(ft)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test2_h_Response_for_Pitch_Velocity_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δe (elevator control action)
% figure(10)
% plot(out.tout, out.elevator_controller.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.elevator_controller.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\delta_e$ $for$ $Pitch$ $\&$ $Velocity$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_e$ $non$ $linear$', '$\delta_e$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_e$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test2_delta_e_for_Pitch_Velocity_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δth (throttle control action)
% figure(11)
% plot(out.tout, out.Trust_controller.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.Trust_controller.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\delta_{th}$ $for$ $Pitch$ $\&$ $Velocity$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_{th}$ $non$ $linear$', '$\delta_{th}$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_{th}$ ', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test2_delta_th_for_Pitch_Velocity_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);


% Test 3

% figure(12)
% plot(out.tout, out.theta_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.theta_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.theta_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$\theta$ $Response$ $for$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\theta$ $command$', '$\theta$ $non$ $linear$', '$\theta$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$\theta$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test3_theta_Response_for_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% figure(13)
% plot(out.tout, out.u_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.u_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.u_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$U$ $Response$ $for$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$U$ $command$', '$U$ $non$ $linear$', '$U$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$U$ $(ft/s)$', 'interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test3_u_Response_for_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% figure(14)
% plot(out.tout, out.gamma_states.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.gamma_states.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\gamma$ $Response$ $for$ $Altitude$ $Hold$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\gamma$ $non$ $linear$', '$\gamma$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$\gamma$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test3_gamma_Response_for_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% figure(15)
% plot(out.tout, out.H_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.H_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.H_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$H$ $Response$ $for$ $Altitude$ $Hold$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$H$ $command$', '$h$ $non$ $linear$', '$h$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$H$ $(ft)$', 'interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test3_h_Response_for_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% figure(16)
% plot(out.tout, out.elevator_controller.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.elevator_controller.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\delta_e$ $for$ $Altitude$ $Hold$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_e$ $non$ $linear$', '$\delta_e$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$\delta_e$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test3_delta_e_for_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% figure(17)
% plot(out.tout, out.Trust_controller.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.Trust_controller.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\delta_{th}$ $for$ $Altitude$ $Hold$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_{th}$ $non$ $linear$', '$\delta_{th}$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$\delta_{th}$ ', 'interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test3_delta_th_for_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);

% %% Test 4
% 
% % Plot β (sideslip angle)
% figure(18)
% plot(out.tout, out.beta_states.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.beta_states.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\beta$ $Response$ $for$ $Lateral$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\beta$ $non$ $linear$', '$\beta$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$\beta$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test4_beta_Response_for_Lateral_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot φ (bank angle)
% figure(19)
% plot(out.tout, out.phi_states.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.phi_states.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\phi$ $Response$ $for$ $Lateral$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\phi$ $non$ $linear$', '$\phi$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12);
% ylabel('$\phi$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test4_phi_Response_for_Lateral_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot ψ (heading angle)
% figure(20)
% plot(out.tout, out.psi_states.Data(:, 1), '--k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.psi_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on
% plot(out.tout, out.psi_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$\psi$ $Response$ $for$ $Lateral$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\psi$ $command$', '$\psi$ $non$ $linear$', '$\psi$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\psi$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test4_psi_Response_for_Lateral_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot h (altitude)
% figure(21)
% plot(out.tout, out.H_states.Data(:, 1), '-.k', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.H_states.Data(:, 2), 'r', 'LineWidth', 2);
% hold on 
% plot(out.tout, out.H_states.Data(:, 3), '--b', 'LineWidth', 2);
% title('$H$ $Response$ $for$ $Lateral$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$H$ $command$', '$h$ $non$ $linear$', '$h$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$H$ $(ft)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test4_h_Response_for_Lateral_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δr (rudder control action)
% figure(22)
% plot(out.tout, out.rudder_controller.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.rudder_controller.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\delta_r$ $for$ $Lateral$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_r$ $non$ $linear$', '$\delta_r$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_r$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test4_delta_r_for_Lateral_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δa (aileron control action)
% figure(23)
% plot(out.tout, out.aileron_controller.Data(:, 1), 'r', 'LineWidth', 2);
% grid on
% grid minor
% hold on
% plot(out.tout, out.aileron_controller.Data(:, 2), '--b', 'LineWidth', 2);
% title('$\delta_a$ $for$ $Lateral$ $Controller$', ...
%     'interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_a$ $non$ $linear$', '$\delta_a$ $linear$', ...
%     'interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_a$ $(degrees)$', 'interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% save_path = fullfile(output_folder, 'Test4_delta_a_for_Lateral_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% %% Test 5
% % Plot β (sideslip angle)
% figure(24)
% plot(out.tout, out.beta_states.Data(:, 1), 'r', 'LineWidth', 2); % Nonlinear
% grid on
% grid minor
% hold on
% plot(out.tout, out.beta_states.Data(:, 2), '--b', 'LineWidth', 2); % Linear
% title('$\beta$ $Response$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\beta$ $non$ $linear$', '$\beta$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\beta$ $(degrees)$', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_beta_Response_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot φ (bank angle)
% figure(25)
% plot(out.tout, out.phi_states.Data(:, 1), 'r', 'LineWidth', 2); % Nonlinear
% grid on
% grid minor
% hold on
% plot(out.tout, out.phi_states.Data(:, 2), '--b', 'LineWidth', 2); % Linear
% title('$\phi$ $Response$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\phi$ $non$ $linear$', '$\phi$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\phi$ $(degrees)$', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_phi_Response_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot ψ (heading angle)
% figure(26)
% plot(out.tout, out.psi_states.Data(:, 1), '-.k', 'LineWidth', 2); % Command
% grid on
% grid minor
% hold on
% plot(out.tout, out.psi_states.Data(:, 2), 'r', 'LineWidth', 2); % Nonlinear
% hold on
% plot(out.tout, out.psi_states.Data(:, 3), '--b', 'LineWidth', 2); % Linear
% title('$\psi$ $Response$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\psi$ $command$', '$\psi$ $non$ $linear$', '$\psi$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\psi$ $(degrees)$', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_psi_Response_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot h (altitude)
% figure(27)
% plot(out.tout, out.H_states.Data(:, 1), '-.k', 'LineWidth', 2); % Command
% grid on
% grid minor
% hold on
% plot(out.tout, out.H_states.Data(:, 2), 'r', 'LineWidth', 2); % Nonlinear
% hold on
% plot(out.tout, out.H_states.Data(:, 3), '--b', 'LineWidth', 2); % Linear
% title('$H$ $Response$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$H$ $command$', '$h$ $non$ $linear$', '$h$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$H$ $(ft)$', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_h_Response_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δr (rudder control action)
% figure(28)
% plot(out.tout, out.rudder_controller.Data(:, 1), 'r', 'LineWidth', 2); % Nonlinear
% grid on
% grid minor
% hold on
% plot(out.tout, out.rudder_controller.Data(:, 2), '--b', 'LineWidth', 2); % Linear
% title('$\delta_r$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_r$ $non$ $linear$', '$\delta_r$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_r$ $(degrees)$', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_delta_r_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δa (aileron control action)
% figure(29)
% plot(out.tout, out.aileron_controller.Data(:, 1), 'r', 'LineWidth', 2); % Nonlinear
% grid on
% grid minor
% hold on
% plot(out.tout, out.aileron_controller.Data(:, 2), '--b', 'LineWidth', 2); % Linear
% title('$\delta_a$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_a$ $non$ $linear$', '$\delta_a$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_a$ $(degrees)$', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_delta_a_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δe (elevator control action)
% figure(30)
% plot(out.tout, out.elevator_controller.Data(:, 1), 'r', 'LineWidth', 2); % Nonlinear
% grid on
% grid minor
% hold on
% plot(out.tout, out.elevator_controller.Data(:, 2), '--b', 'LineWidth', 2); % Linear
% title('$\delta_e$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_e$ $non$ $linear$', '$\delta_e$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_e$ $(degrees)$', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_delta_e_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% 
% % Plot δth (throttle control action)
% figure(31)
% plot(out.tout, out.Trust_controller.Data(:, 1), 'r', 'LineWidth', 2); % Nonlinear
% grid on
% grid minor
% hold on
% plot(out.tout, out.Trust_controller.Data(:, 2), '--b', 'LineWidth', 2); % Linear
% title('$\delta_{th}$ $for$ $Lateral$ $\&$ $Altitude$ $Hold$ $Controller$', ...
%     'Interpreter', 'latex', 'FontSize', 15);
% legend('$\delta_{th}$ $non$ $linear$', '$\delta_{th}$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\delta_{th}$ ', 'Interpreter', 'latex', 'FontSize', 12);
% save_path = fullfile(output_folder, 'Test5_delta_th_for_Lateral_Altitude_Hold_Controller.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);

% % Test 6
% %%
% % Plot ψ (heading angle) for Mission
% figure(32)
% plot(out.tout, out.psi_states.Data(:, 1), '-.k', 'LineWidth', 2); % Command
% grid on
% grid minor
% hold on
% plot(out.tout, out.psi_states.Data(:, 2), 'r', 'LineWidth', 1); % Nonlinear
% hold on
% plot(out.tout, out.psi_states.Data(:, 3), '--b', 'LineWidth', 2); % Linear
% title('$\psi$ $Mission$', 'Interpreter', 'latex', 'FontSize', 15);
% legend('$\psi$ $command$', '$\psi$ $non$ $linear$', '$\psi$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$\psi$ $(degrees)$', 'Interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% % % save_path = fullfile(output_folder, 'Test6_Mission_psi.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);
% % Plot H (altitude) for Mission
% figure(33)
% plot(out.tout, out.H_states.Data(:, 1), '-.k', 'LineWidth', 2); % Command
% grid on
% grid minor
% hold on
% plot(out.tout, out.H_states.Data(:, 2), 'r', 'LineWidth', 1); % Nonlinear
% hold on
% plot(out.tout, out.H_states.Data(:, 3), '--b', 'LineWidth', 2); % Linear
% title('$H$ $Mission$', 'Interpreter', 'latex', 'FontSize', 15);
% legend('$H$ $command$', '$H$ $non$ $linear$', '$H$ $linear$', ...
%     'Interpreter', 'latex', 'FontSize', 10);
% xlabel('Time (s)', 'Interpreter', 'latex', 'FontSize', 12); % Added x-axis label
% ylabel('$H$ $(ft)$', 'Interpreter', 'latex', 'FontSize', 12);
% % Save Figure
% % save_path = fullfile(output_folder, 'Test6_Mission_H.png');
% exportgraphics(gcf, save_path, 'Resolution', 300);