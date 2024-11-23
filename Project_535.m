% Francois-Xavier Duclos 261050648
% 535 Project

% ------- Step 0: Variable identification -------

gamma = 1.4;
R_constant = 287;

r_h = 0.45;              % Hub radius (minimum r value)
r_s = 0.5;               % Shroud radius (maximum r value)
C_x_inlet = 136;
rho_inlet = 1.5;
T_inlet_static = 25 + 273.15;
A_inlet = pi * (r_s^2 - r_h^2);                     
c_p = 1005;
alpha = 30;                         % Angle for rC_theta calculation in degrees
tan_alpha = tan(deg2rad(alpha));    % Convert angle to radians and calculate tangent

N = 6000 * (2*pi)/60; % [rad/s]

delta_rC_theta_h = 82.3; % [m^2/s]
delta_rC_theta_s = 84.4; % [m^2/s]

% ------- Step 1: Grid -------

% Define parameters for the grid
num_points_x = 15;      % Number of grid points in the x direction
num_points_r = 15;      % Number of grid points in the r direction
blade_width = 0.1;      % Blade width (arbitrary value, adjust based on system)

% Define the x and r range based on blade width
x_min = 0;               
x_max = 5 * blade_width; 

% Generate grid points in x and r directions
x_values = linspace(x_min, x_max, num_points_x);
r_values = linspace(r_h, r_s, num_points_r);

% Create mesh grid for x and r coordinates
[X, R] = meshgrid(x_values, r_values);

% Initialize Psi_values as a 2D array to store Psi(x, r) for each grid point
Psi_values = zeros(num_points_r, num_points_x);

% Nested loop to iterate through each (x, r) point and calculate Psi
for i = 1:num_points_r
    for j = 1:num_points_x
        Psi_values(i, j) = calculatePsi(R(i, j), r_h, r_s);
    end
end

function Psi = calculatePsi(r, r_h, r_s)
    % This function calculates Psi based on the radial position r,
    % the hub radius r_h, and the shroud radius r_s.
    Psi = (r^2 - r_h^2) / (r_s^2 - r_h^2);

    % Ensure Psi is zero when r is exactly equal to r_h
    if r == r_h
        Psi = 0;
    end
end

% ------- Step 2: Variables initiatization -------

% Variables Initialization
P_inlet_static = rho_inlet * R * T_inlet_static;            % Inlet Static Pressure
m_dot = rho_inlet * C_x_inlet * A_inlet;                    % Mass flow rate
H0_inlet = c_p * T_inlet_static + (C_x_inlet^2) / 2;        % Total enthalpy at inlet

% Calculating Velocities

% Initialize matrices to store velocities
C_x = C_x_inlet * ones(num_points_r, num_points_x);
C_r = zeros(num_points_r, num_points_x);
rC_theta = zeros(num_points_r, num_points_x);

dx = x_values(2) - x_values(1); % Spacing in the x direction
dr = r_values(2) - r_values(1); % Spacing in the r direction

% Calculate C_x and C_r 
for i = 2:num_points_r-1        % Exclude boundary nodes for r
    for j = 2:num_points_x-1    % Exclude boundary nodes for x

        % Calculate C_x
        C_x(i, j) = m_dot / (2 * pi * rho_inlet * R(i, j)) * ...
                    (Psi_values(i, j + 1) - Psi_values(i, j - 1)) / (2 * dr);
        
        % Calculate C_r
        C_r(i, j) = -m_dot / (2 * pi * rho_inlet * R(i, j)) * ...
                    (Psi_values(i + 1, j) - Psi_values(i - 1, j)) / (2 * dx);
    end
end

% Initialize C_m
C_m = C_x * ones(num_points_r, num_points_x);

for i = 2:num_points_r-1        % Exclude boundary nodes for r
    for j = 2:num_points_x-1    % Exclude boundary nodes for x
        
        C_m(i,j) = sqrt(C_r(i, j)^2 + C_x(i, j)^2);
    end
end

% Initialize rC_theta
for i = 1:num_points_r
    for j = 1:num_points_x
        rC_theta(i, j) = R(i, j) * C_x_inlet * tan_alpha;
    end
end

% Initialize Loss Factor

% Define the x positions of the leading edge (LE) and trailing edge (TE)
x_LE = blade_width;           % x position for LE
x_TE = 2 * blade_width;       % x position for TE

% Find the indices in x_values closest to x_LE and x_TE
[~, i_LE] = min(abs(x_values - x_LE));
[~, i_TE] = min(abs(x_values - x_TE));

% Initialize the loss factor matrix w as zeros
w = zeros(num_points_r, num_points_x);

% Calculate the uniform loss factor for the blade section
loss_factor = 0.5 / (i_TE - i_LE);

% Apply the loss factor uniformly between i = LE + 1 and i = TE
for i = (i_LE + 1):i_TE
    w(:, i) = loss_factor;  
end

% ------- Step 3: Start the rotor -------
 
% Define the parameters
M = num_points_r;  % Number of radial points from hub to shroud

% Initialize rC_theta and calculate the swirl increase at TE
rC_theta = zeros(num_points_r, num_points_x);  % Initialize rC_theta to zero

% Calculate the prescribed swirl increase at the trailing edge for each radial station
delta_rC_theta_TE = zeros(1, M);

for j = 1:M
    delta_rC_theta_TE(j) = delta_rC_theta_h + (delta_rC_theta_s - delta_rC_theta_h) * (j - 1) / (M - 1);
    delta_rC_theta_TE = delta_rC_theta_TE.';
end

% Distribute rC_theta incrementally from i_LE to i_TE along the axial direction
for j = 1:M  % Loop over radial points from hub to shroud
    for i = (i_LE + 1):i_TE  % Loop over axial points from L.E. to T.E.
        % Linearly interpolate in the x-direction between i_LE and i_TE
        rC_theta(i, j) = delta_rC_theta_TE(j) * (i - i_LE) / (i_TE - i_LE);
    end
end

rC_theta = rC_theta.';
% Calculate C_theta
C_theta = rC_theta ./ R; 

% ------- Step 3.5: Calculation Block -------

% Initialize U as zeros across the entire grid
U = zeros(num_points_r, num_points_x);

% Assign non-zero values of U only within the blade region
for i = i_LE:i_TE
    for j = 1:num_points_r
        U(j, i) = N * R(j); 
    end
end

% C_theta_2 
C_theta_2 = C_theta(:,2);

% V_theta_2
V_theta_2 = zeros(num_points_r, 1);

for j = 1:num_points_r
    V_theta_2(j) = N * R(j,2) - C_theta_2(j);
end

% Beta_2 in rad
Beta_2 = zeros(num_points_r, 1); 
Beta_2_deg = zeros(num_points_r, 1); 

for j = 1:num_points_r
    Beta_2(j) = atan(V_theta_2(j) / C_m(j,2));
    Beta_2_deg(j) = rad2deg(Beta_2(j));
end

% V_2
V_2 = zeros(num_points_r, 1);

for j = 1:num_points_r
    V_2(j) = sqrt(C_m(j,2)^2 + V_theta_2(j)^2);
end

% Initialize T_o1 array to store total temperature for each radial position j
T_o1 = zeros(num_points_r, 1);
h_o1 = zeros(num_points_r, 1);

% Calculate T_o1 for each radial position j
for j = 1:num_points_r
    T_o1(j) = T_inlet_static + (C_x_inlet^2 + C_theta(j, 1)^2) / (2 * c_p);
    h_o1(j) = c_p * T_o1(j);
end

% Initialize T_o1_rel array to store total temperature for each radial position j
T_o1_rel = zeros(num_points_r, 1);
h_o1_rel = zeros(num_points_r, 1);

% Calculate T_o1_rel for each radial position j
for j = 1:num_points_r
    T_o1_rel(j) = T_o1(j) + (U(j, 1)^2 - 2 * U(j, 1) * C_theta(j, 1)) / (2 * c_p);
    h_o1_rel(j) = c_p * T_o1_rel(j);
end

% Initialize P_o1_rel array to store total temperature for each radial position j
P_o1_rel = zeros(num_points_r, 1);

% Calculate P_o1_rel for each radial position j
for j = 1:num_points_r
    P_o1_rel(j) = P_inlet_static(j,1) * (T_o1_rel(j) / T_o1(j))^(gamma / (gamma - 1));
end

S_1 = zeros(num_points_r, 1);   

% Initialize T_o1_rel array to store total temperature for each radial position j
h_o2_rel = zeros(num_points_r, 1);

% Calculate T_o1_rel for each radial position j
for j = 1:num_points_r
    h_o2_rel(j) = h_o1_rel(j) - (U(j, 1)^2) / 2 + (U(j, 2)^2) / 2;
end

% h_2
h_2 = zeros(num_points_r, 1);

for j = 1:num_points_r
    h_2(j) = h_o2_rel(j) - V_2(j)^2 / 2;
end

% h_o2
h_o2 = zeros(num_points_r, 1);

for j = 1:num_points_r
    h_o2(j) = h_2(j) - (C_m(j,2)^2 + C_theta(j,2)^2) / 2;
end

% Initialize P_o2_rel_i array to store total temperature for each radial position j
P_o2_rel_i = zeros(num_points_r, 1);

% Calculate P_o2_rel_i for each radial position j
for j = 1:num_points_r
    P_o2_rel_i(j) = P_o1_rel(j) * (h_o2_rel(j) / h_o1_rel(j))^(gamma / (gamma - 1));
end

% Initialize P_o2_rel array to store total temperature for each radial position j
P_o2_rel = zeros(num_points_r, 1);

% Calculate P_o2_rel for each radial position j
for j = 1:num_points_r
    P_o2_rel(j) = P_o2_rel_i(j) - w(j, 2) * (P_o1_rel(j) - P_inlet_static(j,1));
end

% Initialize P_o2 array to store total temperature for each radial position j
P_o2 = zeros(num_points_r, 1);

% Calculate P_o2 for each radial position j
for j = 1:num_points_r
    P_o2(j) = P_o2_rel(j) * (h_o2(j) / h_o2_rel(j))^(gamma / (gamma - 1));
end

% Initialize P_2 array to store total temperature for each radial position j
P_2 = zeros(num_points_r, 1);

% Calculate P_2 for each radial position j
for j = 1:num_points_r
    P_2(j) = P_o2(j) * (h_2(j) / h_o2(j))^(gamma / (gamma - 1));
end

% Initialize rho_2 array to store total temperature for each radial position j
rho_2 = zeros(num_points_r, 1);

% Calculate rho_2 for each radial position j
for j = 1:num_points_r
    rho_2(j) = c_p *  P_2(j) / (R_constant * h_2(j));
end

% Initialize S_2 array to store total temperature for each radial position j
S_2 = zeros(num_points_r, 1);

% Calculate S_2 for each radial position j
for j = 1:num_points_r
    S_2(j) = S_1(j) + c_p *  log(h_o2_rel(j) / h_o1_rel(j)) - R_constant * log(P_o2_rel(j) / P_o1_rel(j));
end

% ------- Step 4: Calculate Vorticity -------

omega = 9;