% Francois-Xavier Duclos 261050648
% 535 Project

%% ------- Step 0: Variable identification -------

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

%% ------- Step 1: Grid -------

% Define parameters for the grid
num_points_x = 15;      % Number of grid points in the x direction
num_points_r = 15;      % Number of grid points in the r direction
blade_width = 0.1;      % Blade width (Specified for this system)

% Define the x and r range based on blade width
x_min = 0;               
x_max = 2 * blade_width; 

% Generate grid points in x and r directions
x_values = linspace(x_min, x_max, num_points_x);
r_values = linspace(r_h, r_s, num_points_r);

% Create mesh grid for x and r coordinates
[X, R] = meshgrid(x_values, r_values);

% Initialize Psi_values as a 2D array to store Psi(x, r) for each grid point
Psi_values = zeros(num_points_r, num_points_x);

% Nested loop to iterate through each (x, r) point and calculate Psi
for r = 1:num_points_r
    for x = 1:num_points_x
        Psi_values(r, x) = calculatePsi(R(r, x), r_h, r_s);
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

%% ------- Step 2: Variables initiatization -------

% Variables Initialization
P_inlet_static = rho_inlet * R * T_inlet_static;            % Inlet Static Pressure
m_dot = rho_inlet * C_x_inlet * A_inlet;                    % Mass flow rate
H0_inlet = c_p * T_inlet_static + (C_x_inlet^2) / 2;        % Total enthalpy at inlet

% Calculating Velocities

% Initialize matrices to store velocities
C_x = zeros(num_points_r, num_points_x);
C_r = zeros(num_points_r, num_points_x);
rC_theta = zeros(num_points_r, num_points_x);

dx = x_values(2) - x_values(1); % Spacing in the x direction
dr = r_values(2) - r_values(1); % Spacing in the r direction

% Calculate C_x and C_r 
for r = 1:num_points_r       % Exclude boundary nodes for r
    for x = 1:num_points_x    % Exclude boundary nodes for x
        if r==1 && x==1 % hub-LE
            % Calculate C_x
            C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                        (Psi_values(r + 1, x) - Psi_values(r, x)) / (dr);
            
            % Calculate C_r
            C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                        (Psi_values(r, x + 1) - Psi_values(r, x)) / (dx);

        elseif r==1 && x==num_points_x % hub-TE
                % Calculate C_x
                C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r + 1, x) - Psi_values(r, x)) / (dr);
                
                % Calculate C_r
                C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x) - Psi_values(r, x - 1)) / (dx);

        elseif r==num_points_r && x==num_points_x % shroud-TE
                % Calculate C_x
                C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x) - Psi_values(r - 1, x)) / (dr);
                
                % Calculate C_r
                C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x) - Psi_values(r, x - 1)) / (dx);

        elseif r==num_points_r && x==1 % shroud-LE
                % Calculate C_x
                C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x) - Psi_values(r - 1, x)) / (dr);
                
                % Calculate C_r
                C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x + 1) - Psi_values(r, x)) / (dx);

        elseif r==1 % hub
            % Calculate C_x
            C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                        (Psi_values(r + 1, x) - Psi_values(r, x)) / (dr);
            
            % Calculate C_r
            C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                        (Psi_values(r, x + 1) - Psi_values(r, x-1)) / (2*dx);

        elseif x==1 % LE
                % Calculate C_x
                C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r + 1, x) - Psi_values(r - 1, x)) / (2 * dr);
                
                % Calculate C_r
                C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x + 1) - Psi_values(r, x)) / (dx);

        elseif r==num_points_r % shroud
                % Calculate C_x
                C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x) - Psi_values(r - 1, x)) / (dr);
                
                % Calculate C_r
                C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x + 1) - Psi_values(r, x - 1)) / (2*dx);

        elseif x==num_points_x % TE
                % Calculate C_x
                C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r + 1, x) - Psi_values(r - 1, x)) / (2*dr);
                
                % Calculate C_r
                C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                            (Psi_values(r, x) - Psi_values(r, x - 1)) / (dx);
        else
            % Calculate C_x
            C_x(r, x) = m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                        (Psi_values(r + 1, x) - Psi_values(r - 1, x)) / (2 * dr);
            
            % Calculate C_r
            C_r(r, x) = -m_dot / (2 * pi * rho_inlet * R(r, x)) * ...
                        (Psi_values(r, x + 1) - Psi_values(r, x - 1)) / (2 * dx);
        end
    end
end

% Initialize C_m
C_m = C_x * ones(num_points_r, num_points_x);

for r = 1:num_points_r      % Exclude boundary nodes for r
    for x = 1:num_points_x   % Exclude boundary nodes for x
        
        C_m(r,x) = sqrt(C_r(r, x)^2 + C_x(r, x)^2);
    end
end

% Initialize rC_theta
for r = 1:num_points_r
    for x = 1:num_points_x
        rC_theta(r, x) = R(r, x) * C_x_inlet * tan_alpha;
    end
end

% Initialize Loss Factor

% Define the x positions of the leading edge (LE) and trailing edge (TE)
x_LE = x_max/(num_points_x-1)*4;           % x position for LE
x_TE = x_LE + blade_width;       % x position for TE

% Find the indices in x_values closest to x_LE and x_TE
[~, i_LE] = min(abs(x_values - x_LE));
[~, i_TE] = min(abs(x_values - x_TE));

% Initialize the loss factor matrix w as zeros
w = zeros(num_points_r, num_points_x);

% Calculate the uniform loss factor for the blade section
loss_factor = 0.5 / (i_TE - i_LE);

% Apply the loss factor uniformly between i = LE + 1 and i = TE
for r = (i_LE + 1):i_TE
    w(:, r) = loss_factor;  
end

% ------- Step 3: Start the rotor -------
 
% Define the parameters
M = num_points_r;  % Number of radial points from hub to shroud

% Initialize rC_theta and calculate the swirl increase at TE
rC_theta = zeros(num_points_r, num_points_x);  % Initialize rC_theta to zero

% Calculate the prescribed swirl increase at the trailing edge for each radial station
delta_rC_theta_TE = zeros(1, M);

for x = 1:M
    delta_rC_theta_TE(x) = delta_rC_theta_h + (delta_rC_theta_s - delta_rC_theta_h) * (x - 1) / (M - 1);
    delta_rC_theta_TE = delta_rC_theta_TE.';
end

% Distribute rC_theta incrementally from i_LE to i_TE along the axial direction
for x = 1:M  % Loop over radial points from hub to shroud
    for r = (i_LE + 1):i_TE  % Loop over axial points from L.E. to T.E.
        % Linearly interpolate in the x-direction between i_LE and i_TE
        rC_theta(r, x) = delta_rC_theta_TE(x) * (r - i_LE) / (i_TE - i_LE);
    end
end

rC_theta = rC_theta.';
% Calculate C_theta
C_theta = rC_theta ./ R; 

% Initialize U as zeros across the entire grid
U = zeros(num_points_r, num_points_x);

% Assign non-zero values of U only within the blade region
for r = i_LE:i_TE
    for x = 1:num_points_r
        U(x, r) = N * R(x); 
    end
end

S_1 = zeros(num_points_r, 1); 

%% ------- Step 3.5: Calculation Block -------

% Define the radial positions and initialize parameters
M = num_points_r;
N = num_points_x;

% Initialize 2D matrices to store properties at each station (each column is a station)
T_o = zeros(M, N);         % Total temperature at each station
h_o = zeros(M, N);         % Total enthalpy at each station
T_o_rel = zeros(M, N);     % Relative total temperature at each station
h_o_rel = zeros(M, N);     % Relative total enthalpy at each station
P_o_rel = zeros(M, N);     % Relative total pressure at each station
S = zeros(M, N);           % Entropy at each station
C_theta_global = zeros(M, N); % Tangential velocity component at each station
V_theta_global = zeros(M, N); % Relative tangential velocity component at each station
Beta_global = zeros(M, N);    % Flow angle at each station
V_global = zeros(M, N);       % Magnitude of velocity at each station

% New global matrices for static temperature and static pressure
T_static_global = zeros(M, N); % Static temperature at each station
P_static_global = zeros(M, N); % Static pressure at each station
rho_global = zeros(M, N);      % Density at each station

% Initialize the first column (station 1) with inlet values

% Initialize T_o1 array to store total temperature for each radial position j
T_o1 = zeros(num_points_r, 1);
h_o1 = zeros(num_points_r, 1);

% Calculate T_o1 for each radial position j
for x = 1:num_points_r
    T_o1(x) = T_inlet_static + (C_x_inlet^2 + C_theta(x, 1)^2) / (2 * c_p);
    h_o1(x) = c_p * T_o1(x);
end

T_o(:, 1) = T_o1;
h_o(:, 1) = h_o1;

% Initialize T_o1_rel array to store total temperature for each radial position j
T_o1_rel = zeros(num_points_r, 1);
h_o1_rel = zeros(num_points_r, 1);

% Calculate T_o1_rel for each radial position j
for x = 1:num_points_r
    T_o1_rel(x) = T_o1(x) + (U(x, 1)^2 - 2 * U(x, 1) * C_theta(x, 1)) / (2 * c_p);
    h_o1_rel(x) = c_p * T_o1_rel(x);
end

T_o_rel(:, 1) = T_o1_rel;
h_o_rel(:, 1) = h_o1_rel;

% Initialize P_o1_rel array to store total pressure for each radial position j
P_o1_rel = zeros(num_points_r, 1);

% Calculate P_o1_rel for each radial position j
for x = 1:num_points_r
    P_o1_rel(x) = P_inlet_static(x,1) * (T_o1_rel(x) / T_o1(x))^(gamma / (gamma - 1));
end

P_o_rel(:, 1) = P_o1_rel;

S(:, 1) = 0;  
C_theta_global(:, 1) = C_theta(:, 1);

% Initialize static temperature and pressure at station 1
T_static_global(:, 1) = T_inlet_static(:,1);  % Set initial static temperature to inlet static temperature
P_static_global(:, 1) = P_inlet_static(:,1);  % Set initial static pressure to inlet static pressure
rho_global(:, 1) = P_static_global(:, 1) ./ (R_constant .* T_static_global(:, 1));


% Loop over each axial station (starting from station 2)
for r = 2:N
    % Compute V_theta and other parameters at each radial position
    for x = 1:M
        % Calculate V_theta and Beta at station i
        V_theta_global(x, r) = N * R(x, r) - C_theta_global(x, r - 1);
        Beta_global(x, r) = atan(V_theta_global(x, r) / C_m(x, r));
        
        % Calculate the magnitude of velocity
        V_global(x, r) = sqrt(C_m(x, r)^2 + V_theta_global(x, r)^2);
        
        % Calculate total temperature and enthalpy at station i
        T_o(x, r) = T_o(x, r - 1) + (C_x(x, r - 1)^2 + C_theta_global(x, r - 1)^2) / (2 * c_p);
        h_o(x, r) = c_p * T_o(x, r);
        
        % Calculate relative total temperature and enthalpy at station i
        T_o_rel(x, r) = T_o_rel(x, r - 1) + (U(x, r)^2 - 2 * U(x, r) * C_theta_global(x, r - 1)) / (2 * c_p);
        h_o_rel(x, r) = c_p * T_o_rel(x, r);
        
        % Calculate relative total pressure at station i
        P_o_rel(x, r) = P_o_rel(x, r - 1) * (h_o_rel(x, r) / h_o_rel(x, r - 1))^(gamma / (gamma - 1));
        
        % Update entropy at station i
        S(x, r) = S(x, r - 1) + c_p * log(h_o_rel(x, r) / h_o_rel(x, r - 1)) - R_constant * log(P_o_rel(x, r) / P_o_rel(x, r - 1));

        % Static temperature and pressure calculations (assuming isentropic relations)
        T_static_global(x, r) = T_o(x, r) - (V_global(x, r)^2 / (2 * c_p));
        P_static_global(x, r) = P_o_rel(x, r) * (T_static_global(x, r) / T_o(x, r))^(gamma / (gamma - 1));
        % Calculate density at station i
        rho_global(x, r) = P_static_global(x, r) ./ (R_constant .* T_static_global(x, r));
    end
    
    % Update C_theta values for station i based on the last calculated station 1 values
    C_theta_global(:, r) = C_theta(:, r);  % Assuming C_theta is defined for each x position in original data
end

% ------- Step 4: Calculate Vorticity -------

omega = zeros(num_points_x, num_points_r);

% Loop through each (i, j) node, excluding boundaries
for r = 2:num_points_x - 1
    for x = 2:num_points_r - 1
        % Term 1: C_theta(i, j) / (R(i, j) * delta_r * m_dot * C_x(i, j))
        term1 = C_theta(r, x) / R(r, x) * (rC_theta(r+1, x) - rC_theta(r-1, x));
        
        % Term 2: T(i, j) / (delta_r * m_dot * C_x(i, j)) * (S(i, j+1) - S(i, j-1))
        term2 = T_static_global(r, x) * (S(r+1, x) - S(r-1, x));
        
        % Term 3: - (1 / (delta_r * m_dot * C_x(i, j))) * (H_o(i, j+1) - H_o(i, j-1))
        term3 = -1 * (h_o(r+1, x) - h_o(r-1, x));
        
        % Calculate omega_i,j
        omega(r, x) = pi / (dr * m_dot*C_x(r,x)) .* (term1 + term2 + term3);
    end
end

% S(:,:) = 0;

% ------- Step 5: Update Psi -------

