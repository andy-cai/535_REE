% Francois-Xavier Duclos 261050648
% 535 Project

%% ------- Step 0: Variable identification -------

gamma = 1.4;
k_gamma = (gamma / (gamma - 1));
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
        Psi_values(r, x) = calculateInitPsi(R(r, x), r_h, r_s);
    end
end

%% ------- Step 2: Variables initiatization -------

% Variables Initialization
P_inlet_static = rho_inlet * R(1:num_points_r, 1) * T_inlet_static;            % Inlet Static Pressure
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

% hub-LE
C_x(1, 1) = m_dot / ((2 * pi * rho_inlet) .* R(1, 1)) * (Psi_values(2, 1) - Psi_values(1, 1)) / (dr);
C_r(1, 1) = -m_dot / (2 * pi * rho_inlet * R(1, 1)) * (Psi_values(1, 2) - Psi_values(1, 1)) / (dx);

% hub-TE
C_x(1, num_points_x) = m_dot / (2 * pi * rho_inlet * R(1, num_points_x)) * (Psi_values(2, num_points_x) - Psi_values(1, num_points_x)) / (dr);
C_r(1, num_points_x) = -m_dot / (2 * pi * rho_inlet * R(1, num_points_x)) * (Psi_values(1, num_points_x) - Psi_values(1, num_points_x - 1)) / (dx);

% shroud-LE
C_x(num_points_r, 1) = m_dot / (2 * pi * rho_inlet * R(num_points_r, 1)) * (Psi_values(num_points_r, 1) - Psi_values(num_points_r - 1, 1)) / (dr);
C_r(num_points_r, 1) = -m_dot / (2 * pi * rho_inlet * R(num_points_r, 1)) * (Psi_values(num_points_r, 2) - Psi_values(num_points_r, 1)) / (dx);

% shroud-TE
C_x(num_points_r, num_points_x) = m_dot / (2 * pi * rho_inlet * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r - 1, num_points_x)) / (dr);
C_r(num_points_r, num_points_x) = -m_dot / (2 * pi * rho_inlet * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r, num_points_x - 1)) / (dx);


% hub
C_x(1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(1, 2:num_points_x-1)) .* (Psi_values(2, 2:num_points_x-1) - Psi_values(1, 2:num_points_x-1)) ./ (dr);
C_r(1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(1, 2:num_points_x-1)) .* (Psi_values(1, 3:num_points_x) - Psi_values(1, 1:num_points_x-2)) ./ (2*dx);

% LE
C_x(2:num_points_r-1, 1) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 1)) .* (Psi_values(3:num_points_r, 1) - Psi_values(1:num_points_r-2, 1)) ./ (2 * dr);
C_r(2:num_points_r-1, 1) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 1)) .* (Psi_values(2:num_points_r-1, 2) - Psi_values(2:num_points_r-1, 1)) ./ (dx);

% shroud
C_x(num_points_r, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 2:num_points_x-1) - Psi_values(num_points_r - 1, 2:num_points_x-1)) ./ (dr);
C_r(num_points_r, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 3:num_points_x) - Psi_values(r, 1:num_points_x-2)) ./ (2*dx);

% TE
C_x(2:num_points_r-1, num_points_x) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(3:num_points_r, num_points_x) - Psi_values(1:num_points_r - 2, num_points_x)) ./ (2*dr);
C_r(2:num_points_r-1, num_points_x) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(2:num_points_r-1, num_points_x) - Psi_values(2:num_points_r-1, num_points_x - 1)) ./ (dx);

% Interior of Blade
C_x(2:num_points_r-1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(3:num_points_r, 2:num_points_x-1) - Psi_values(1:num_points_r-2, 2:num_points_x-1)) ./ (2 * dr);
C_r(2:num_points_r-1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(2:num_points_r-1, 3:num_points_x) - Psi_values(2:num_points_r-1, 1:num_points_x-2)) ./ (2 * dx);


% Initialize C_m
C_m = sqrt(C_r.^2 + C_x.^2);

% Initialize rC_theta 
%% Assignment says 'initialize everywhere' does that include in front of LE?
%% if we want whirl to be = 0 before LE, need to only do after LE
rC_theta = R .* C_x_inlet * tan_alpha;

% Initialize Loss Factor

% Define the x positions of the leading edge (LE) and trailing edge (TE)
n_cells_before_LE = 4;
x_LE = x_max/(num_points_x-1) * n_cells_before_LE;   % x position for LE
x_TE = x_LE + blade_width;                           % x position for TE

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

% Calculate the prescribed swirl increase at the trailing edge for each radial station
delta_rC_theta_TE = zeros(M, 1);

for r = 1:M
    delta_rC_theta_TE(r) = delta_rC_theta_h + (delta_rC_theta_s - delta_rC_theta_h) * (r - 1) / (M - 1);
end
%% NEED TO CONFIRM ORDER IS CORRECT, WE WANT TO DISPLAY HUB at r=1, SHROUD at r = M

%% Setting whirl before LE to zero, can include if needed, don't think this is correct (in free vortex)
% for x = 1:i_LE  % Loop over axial points before LE
%     % Linearly interpolate in the x-direction between i_LE and i_TE
%     rC_theta(1:M, x) = zeros(M, 1);
% end

for x = (i_LE + 1):i_TE  % Loop over axial points from L.E. to T.E.
    % Linearly interpolate in the x-direction between i_LE and i_TE
    rC_theta(1:M, x) = rC_theta(1:M, x - 1) + delta_rC_theta_TE .* ((x - i_LE) / (i_TE - i_LE));
end

% Initialize U as zeros across the entire grid
U = zeros(num_points_r, num_points_x);

% Assign non-zero values of U only within the blade region
U(1:num_points_r, i_LE:i_TE) = N .* R(1:num_points_r, i_LE:i_TE); 

% Initialize entropy to zero
S_1 = zeros(num_points_r, 1);

%% ------- Step 3.5a: Calculation Block: Initialization -------

% Define the radial positions and initialize parameters
M = num_points_r;
N = num_points_x;

% Initialize 2D matrices to store properties at each station (each column is a station)
T_o = zeros(M, N);         % Total temperature at each station
h_o = zeros(M, N);         % Total enthalpy at each station
T_o_rel = zeros(M, N);     % Relative total temperature at each station
h_o_rel = zeros(M, N);     % Relative total enthalpy at each station
P_o_rel = zeros(M, N);     % Relative total pressure at each station
P_o_rel_ideal = zeros(M, N); % ideal (isentropic total pressure at each station)
S = zeros(M, N);           % Entropy at each station
% C_theta_global = zeros(M, N); % Tangential velocity component at each station (already defined ? )
V_theta_global = zeros(M, N); % Relative tangential velocity component at each station
Beta_global = zeros(M, N);    % Flow angle at each station
V_global = zeros(M, N);       % Magnitude of velocity at each station

% New global matrices for static temperature and static pressure
T_static_global = zeros(M, N); % Static temperature at each station
P_static_global = zeros(M, N); % Static pressure at each station
rho_global = zeros(M, N);      % Density at each station

loss_global = zeros(M, N);


% Initialize the first column (station 1) with inlet values
% % Station 1: Before leading edge at inlet: x = 1
% % Station 2: After trailing edge at exit: x = N
%% someone confirm this isn't LE, TE

%% ------- Step 3.5b: Calculation Block: Calculations -------

% Calculate C_theta
C_theta = rC_theta ./ R; 

%% Inlet
% Calculate T_o1 for each radial position j at entry
T_o1 = T_inlet_static + (C_x(1:M, 1).^2 + C_theta(1:M, 1).^2) / (2 * c_p);
h_o1 = c_p .* T_o1;

% populate initial column of static temperature and enthalpy
T_o(1:M, 1) = T_o1;
h_o(1:M, 1) = h_o1;

% Initialize T_o1_rel array to store relative total temperature for each radial position j
% % only valid over entire domain as U initialized as 0 outside rotor
T_o1_rel = T_o1 + (U(1:M, 1) .* (U(1:M, 1) - 2.*C_theta(1:M, 1)) / (2 * c_p));
h_o1_rel = c_p .* T_o1_rel;

% populate initial column of relative temp and enthalpy
T_o_rel(:, 1) = T_o1_rel;
h_o_rel(:, 1) = h_o1_rel;

% Calculate P_o1_rel for each radial position j
P_o1_rel = P_inlet_static(1:M) .* (T_o1_rel ./ T_o1).^k_gamma;

% populate initial column of relative pressure
P_o_rel(:, 1) = P_o1_rel;

P_o_rel_ideal(:, 1) = P_o1_rel; %% ** confirm **

% initialize entropy at inlet to be zero
S(:, 1) = 0;  

% initialize loss factor
loss_global(1:M, (i_LE + 1):i_TE) = loss_factor;

% Initialize static temperature and pressure at station 1
T_static_global(:, 1) = T_inlet_static(:,1);  % Set initial static temperature to inlet static temperature
P_static_global(:, 1) = P_inlet_static(:);  % Set initial static pressure to inlet static pressure
rho_global(:, 1) = P_static_global(:, 1) ./ (R_constant .* T_static_global(:, 1));

% calculate I at each radial location (inlet)
I_inlet = h_o1_rel - (0.5 * U(1:M, 1).^2);

%% * confirm *
% note that I = H0_rel - 0.5*U^2 is constant throughout (Conservation of Rothalpy)

%% ERROR IN S: NEGATIVE VALUES CALCULATED


%% Outlet
% these calcs are not required, but retained to match notes for clarity
V_theta_global(1:M, N) = U(1:M, N) - C_theta(1:M, N);
Beta_global(1:M, N) = atan(V_theta_global(1:M, N) ./ C_m(1:M, N));
V_global(1:M, N) = sqrt(C_m(1:M, N).^2 + V_theta_global(1:M, N).^2);


% Move from 1 -> N updating using common calc block
% note that I and velocities are already calculated

% relative quantities from Conservation of Rothalpy (assuming no c_r for first iteration)
h_o_rel(1:M, 2:N) = I_inlet + 0.5 * U(1:M, 2:N).^2;
T_o_rel(1:M, 2:N) = h_o_rel(1:M, 2:N) ./ c_p;

% stagnation temperatures from velocities and T_o_rel
T_o(1:M, 2:N) = T_o_rel(1:M, 2:N) - (U(1:M, 2:N) .* (U(1:M, 2:N) - 2.*C_theta(1:M, 2:N)) ./ (2 * c_p));

% constant T_o after TE (no work)
T_o(1:M, (i_TE+1):N) = repmat(T_o(1:M, i_TE), 1, N - i_TE);

% constant T_o_rel = T_o after TE (no work, no U)
T_o(1:M, (i_TE+1):N) = repmat(T_o(1:M, i_TE), 1, N - i_TE);

% constant T_o_rel = T_o after TE (no work, no U)
T_o_rel(1:M, (i_TE+1):N) = T_o(1:M, (i_TE+1):N);

h_o(1:M, 2:N) = c_p .* T_o(1:M, 2:N);
h_o_rel(1:M, 2:N) = c_p .* T_o_rel(1:M, 2:N);

% static properties from stagnation and velocities
T_static_global(1:M, 2:N) = T_o(1:M, 2:N) - (V_global(1:M, 2:N).^2 ./ (2 * c_p));


% finding pressures
for x = 2:N
    % isentropic relative pressure
    P_o_rel_ideal(1:M, x) = P_o_rel_ideal(1:M, x-1) .* ((T_o_rel(1:M, x) ./ T_o_rel(1:M, x-1)).^k_gamma);
    % relative pressure using loss coeff
    P_o_rel(1:M, x) = P_o_rel_ideal(1:M, x) - loss_global(1:M, x) .* (P_o_rel(1:M, x-1) - P_static_global(1:M, x-1));

    % static pressure from relative and stagnation properties
    P_static_global(1:M, x) = P_o_rel(1:M, x) .* (T_o(1:M, x) ./ T_o_rel(1:M, x)).^(k_gamma);
end

% density throughout
rho_global(1:M, 2:N) = P_static_global(1:M, 2:N) ./ (R_constant .* T_static_global(1:M, 2:N));

% update entropy
for x = 2:N    
    S(1:M, x) = S(1:M, x - 1) + c_p .* log(h_o_rel(1:M, x) ./ h_o_rel(1:M, x - 1)) - R_constant * log(P_o_rel(1:M, x) ./ P_o_rel(1:M, x - 1));
end


%% Iteration Loop
stop_condition = 1;
conv = [];

min_iter = 50;
max_iter = 1000;
iteration = 1;

while (stop_condition && iteration < max_iter) || iteration <= min_iter
        
    % ------- Step 4: Calculate Vorticity -------
    % exclude hub, shroud, LE = (2:M-1, 2:N)
    
    % Term 1: C_theta(i, j) / (R(i, j) * (rC(i, j+1) - rC(i, j-1))
    term1 = C_theta(2:M-1, 2:N) ./ R(2:M-1, 2:N) .* (rC_theta(3:M, 2:N) - rC_theta(1:M-2, 2:N));
    
    % Term 2: T(i, j) * (S(i, j+1) - S(i, j-1))
    term2 = T_static_global(2:M-1, 2:N) .* (S(3:M, 2:N) - S(1:M-2, 2:N));
    
    % Term 3: (H_o(i, j+1) - H_o(i, j-1))
    term3 = (h_o(3:M, 2:N) -  h_o(1:M-2, 2:N));
    
    % Denom: delta_r * m_dot * c_x / pi
    denom = (dr .* m_dot .* C_x(2:M-1, 2:N)) ./ pi;
    
    % Note: doesn't include inlet, hub, shroud
    omega = (term1 + term2 - term3) ./ denom;
    
    clear term1 term2 term3 denom
    
    % set entropy for iteration to zero
    S_old = S; % store last iteration for convergence
    Psi_old = Psi_values;
    S(:, :) = 0;
    
    % ------- Step 5: Update Psi -------
    
    %% NOTE: I do not know how to compute the (rho*r)1/2 -> is there a formula or is this the average of adjacent cells?
    %% -> currently assuming the average of adjacent cells ***** CHECK ******
    
    % computing adjacent averages of rho*r
    % % technically plus and minus r are reversed, as indices reversed
    rho_r_plus_x = (rho_global(2:M-1, 3:N).*R(2:M-1, 3:N) + rho_global(2:M-1, 2:N-1).* R(2:M-1, 2:N-1)) / 2; % average of rho*r at cell and at next cell in x
    rho_r_minus_x = (rho_global(2:M-1, 2:N-1).*R(2:M-1, 2:N-1) + rho_global(2:M-1, 1:N-2).* R(2:M-1, 1:N-2)) / 2; % average of rho*r at cell and at prior cell in x
    rho_r_plus_r = (rho_global(3:M, 2:N-1).*R(3:M, 2:N-1) + rho_global(2:M-1, 2:N-1).* R(2:M-1, 2:N-1)) / 2; % average of rho*r at cell and at next cell in r (towards shroud)
    rho_r_minus_r = (rho_global(2:M-1, 2:N-1).*R(2:M-1, 2:N-1) + rho_global(1:M-2, 2:N-1).* R(1:M-2, 2:N-1)) / 2; % average of rho*r at cell and at prior cell in r (towards hub)
    
    % coefficients
    A = (rho_r_plus_x.^(-1) + rho_r_minus_x.^(-1) + rho_r_plus_r.^(-1) + rho_r_minus_r.^(-1)).^(-1);
    
    B = (Psi_values(2:M-1, 3:N) ./ rho_r_plus_x) + (Psi_values(2:M-1, 1:N-2) ./ rho_r_minus_x) + (Psi_values(3:M, 2:N-1) ./ rho_r_plus_r) + (Psi_values(1:M-2, 2:N-1) ./ rho_r_minus_r);
    
    % update psi
    Psi_values(2:M-1, 2:N-1) = A .* (B + dx^2.*omega(:, 2:N-1));
    
    % update exit
    A_exit = (2*(rho_r_minus_x(:, end)).^(-1) + rho_r_plus_r(:, end).^(-1) + rho_r_minus_r(:, end).^(-1)).^(-1);
    B_exit = (2.*Psi_values(2:M-1, N-1) ./ rho_r_minus_x(:, end)) + (Psi_values(3:M, N) ./ rho_r_plus_r(:, end)) + (Psi_values(1:M-2, N) ./ rho_r_minus_r(:, end));
    
    Psi_values(2:M-1, N) = A_exit .* (B_exit + dx^2.*omega(:, N-1));
    % note that omega doesn't include leading edge, so N-1 = N (exit)
    
    
    % ------- Step 6: Calculate velocities -------
        % (inlet remains constant)
        % repeat of prior calcs
    
    % hub-TE
    C_x(1, num_points_x) = m_dot / (2 * pi * rho_inlet * R(1, num_points_x)) * (Psi_values(2, num_points_x) - Psi_values(1, num_points_x)) / (dr);
    C_r(1, num_points_x) = -m_dot / (2 * pi * rho_inlet * R(1, num_points_x)) * (Psi_values(1, num_points_x) - Psi_values(1, num_points_x - 1)) / (dx);
    
    % shroud-TE
    C_x(num_points_r, num_points_x) = m_dot / (2 * pi * rho_inlet * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r - 1, num_points_x)) / (dr);
    C_r(num_points_r, num_points_x) = -m_dot / (2 * pi * rho_inlet * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r, num_points_x - 1)) / (dx);
    
    % hub
    C_x(1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(1, 2:num_points_x-1)) .* (Psi_values(2, 2:num_points_x-1) - Psi_values(1, 2:num_points_x-1)) ./ (dr);
    C_r(1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(1, 2:num_points_x-1)) .* (Psi_values(1, 3:num_points_x) - Psi_values(1, 1:num_points_x-2)) ./ (2*dx);
    
    % shroud
    C_x(num_points_r, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 2:num_points_x-1) - Psi_values(num_points_r - 1, 2:num_points_x-1)) ./ (dr);
    C_r(num_points_r, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 3:num_points_x) - Psi_values(r, 1:num_points_x-2)) ./ (2*dx);
    
    % TE
    C_x(2:num_points_r-1, num_points_x) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(3:num_points_r, num_points_x) - Psi_values(1:num_points_r - 2, num_points_x)) ./ (2*dr);
    C_r(2:num_points_r-1, num_points_x) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(2:num_points_r-1, num_points_x) - Psi_values(2:num_points_r-1, num_points_x - 1)) ./ (dx);
    
    % Interior of Blade
    C_x(2:num_points_r-1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(3:num_points_r, 2:num_points_x-1) - Psi_values(1:num_points_r-2, 2:num_points_x-1)) ./ (2 * dr);
    C_r(2:num_points_r-1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(2:num_points_r-1, 3:num_points_x) - Psi_values(2:num_points_r-1, 1:num_points_x-2)) ./ (2 * dx);
    
    % calculate C_m
    C_m = sqrt(C_r.^2 + C_x.^2);
    
    % ------- Step 7: Calculation Block: Calculations -------> Start at TE and trace backwards to LE
    % swirl remains the same, imposed by work distribution
    
    % interpolation factors % [(M - 2) x (N - 1)]
    % % easier as psi decreases from shroud -> hub
    % r = 0 at hub, r = M at shroud
    
    %% Interpolation - Tracing streamline backwards
    % streamlines at hub and shroud not calculated, as straight lines
    a_towards_shroud = (Psi_values(2:M-1, 2:N) - Psi_values(1:M-2, 1:N-1)) ./ (Psi_values(2:M-1, 1:N-1) - Psi_values(1:M-2, 1:N-1));
    shroud_side = ~(a_towards_shroud > 1);
    a_towards_shroud(~shroud_side) = 0; % removing invalid interpolations (streamline goes other direction)
    b_towards_shroud = 1 - a_towards_shroud;
    
    b_towards_hub = (Psi_values(2:M-1, 2:N) - Psi_values(2:M-1, 1:N-1)) ./ (Psi_values(3:M, 1:N-1) - Psi_values(2:M-1, 1:N-1));
    hub_side = ~(b_towards_hub < 0);
    b_towards_hub(~hub_side) = 1; % removing invalid interpolations (streamline goes other direction)
    a_towards_hub = 1 - b_towards_hub;
    
    
    %% Note: remember to swap stag / rel properties at LE and TE depending on the cells (slide 23)
    
    %% Note: There has got to be a cleaner way to do this -_-
    
    for x = 2:N % from inlet to exit (inlet already defined)
        for r = 2:M-1
    
            if x <= i_LE || x > i_TE % for values before LE or after TE, use stagnation enthalpy
                
                if shroud_side(r-1, x-1) % if streamline goes towards shroud
                    h_o_streamline = h_o(r-1, x-1).* a_towards_shroud(r-1, x-1) + h_o(r, x-1).* b_towards_shroud(r-1, x-1); % h_o at origin (i-1)
                    
                    % calculate other properties at streamline origin
                    R_streamline = R(r-1, x-1).* a_towards_shroud(r-1, x-1) + R(r, x-1).* b_towards_shroud(r-1, x-1); % R at origin (i-1)
                    % U_streamline = U(r-1, x-1).* a_towards_shroud(r-1, x-1) + U(r, x-1).* b_towards_shroud(r-1, x-1); % U at origin ** is this zero outside of blade region? **
    
    
                else % if streamline goes towards hub
                    h_o_streamline = h_o(r, x-1).* a_towards_hub(r-1, x-1) + h_o(r+1, x-1).* b_towards_hub(r-1, x-1); % h_o at origin (i-1)
                    
                    % calculate other properties at streamline origin
                    R_streamline = R(r, x-1).* a_towards_hub(r-1, x-1) + R(r+1, x-1).* b_towards_hub(r-1, x-1); % R at origin (i-1)
                    % U_streamline = U(r, x-1).* a_towards_hub(r-1, x-1) + U(r+1, x-1).* b_towards_hub(r-1, x-1); % U at origin (i-1) ** is this zero outside of blade region? **
    
                end
                    U_streamline = 0; % depends on prof's answer
    
                    % trace forward to define enthalpy (i)
                    h_o(r, x) = h_o_streamline - (0.5 * U_streamline^2) + (0.5 * U(r, x)^2); % defining same point, just from different cell depending on sttreamline
    
                    % update relative enthalpy
                    h_o_rel(r, x) = h_o(r, x) + U(r, x) * (U(r, x) - 2*C_theta(r, x)) / 2;
        
            else % for values inside blade, use relative enthalpy
    
                if shroud_side(r-1, x-1) % if streamline goes towards shroud
                    h_o_rel_streamline = h_o_rel(r-1, x-1).* a_towards_shroud(r-1, x-1) + h_o_rel(r, x-1).* b_towards_shroud(r-1, x-1); % h_o_rel at origin (i-1)
                    
                 % calculate other properties at streamline origin
                    R_streamline = R(r-1, x-1).* a_towards_shroud(r-1, x-1) + R(r, x-1).* b_towards_shroud(r-1, x-1); % R at origin (i-1)
                    U_streamline = U(r-1, x-1).* a_towards_shroud(r-1, x-1) + U(r, x-1).* b_towards_shroud(r-1, x-1); % U at origin (non-zero)
    
    
                else % if streamline goes towards hub
                    h_o_rel_streamline = h_o_rel(r, x-1).* a_towards_hub(r-1, x-1) + h_o_rel(r+1, x-1).* b_towards_hub(r-1, x-1); % h_o_rel at origin (i-1)
                    
                    % calculate other properties at streamline origin
                    R_streamline = R(r, x-1).* a_towards_hub(r-1, x-1) + R(r+1, x-1).* b_towards_hub(r-1, x-1); % R at origin (i-1)
                    U_streamline = U(r, x-1).* a_towards_hub(r-1, x-1) + U(r+1, x-1).* b_towards_hub(r-1, x-1); % U at origin (i-1) (non-zero)
    
                end
    
                % trace forward to define relative enthalpy (i)
                h_o_rel(r, x) = h_o_rel_streamline - (0.5 * U_streamline^2) + (0.5 * U(r, x)^2); % defining same point, just from different nodes depending on streamline
    
                % update enthalpy
                h_o(r, x) = h_o_rel(r, x) - U(r, x) * (U(r, x) - 2*C_theta(r, x)) / 2;
            end
        end
    end
    
    %% Common calculation block from h_o / h_o_rel & velocities
    
    % V_theta_global is constant (U & C_theta don't change)
    Beta_global(1:M, N) = atan(V_theta_global(1:M, N) ./ C_m(1:M, N));
    V_global(1:M, N) = sqrt(C_m(1:M, N).^2 + V_theta_global(1:M, N).^2);
    
    % stagnation (rel & abs) temperatures from h_o & h_o_rel -> from CoRothalpy
    T_o_rel(1:M, 2:N) = h_o_rel(1:M, 2:N) ./ c_p;
    T_o(1:M, 2:N) = h_o(1:M, 2:N) ./ c_p;
    
    % static properties from stagnation and velocities
    T_static_global(1:M, 2:N) = T_o(1:M, 2:N) - (V_global(1:M, 2:N).^2 ./ (2 * c_p));
    
    % sweeping through domain to update pressures
    for x = 2:N
        % isentropic relative pressure
        P_o_rel_ideal(1:M, x) = P_o_rel_ideal(1:M, x-1) .* ((h_o_rel(1:M, x) ./ h_o_rel(1:M, x-1)).^k_gamma);
    
        % relative pressure using loss coeff
        P_o_rel(1:M, x) = P_o_rel_ideal(1:M, x) - loss_global(1:M, x) .* (P_o_rel(1:M, x-1) - P_static_global(1:M, x-1));
    
        % static pressure from relative and stagnation properties
        P_static_global(1:M, x) = P_o_rel(1:M, x).*(T_o(1:M, x) ./ T_o_rel(1:M, x)).^(k_gamma);
    end
    
    % density throughout
    rho_global(1:M, 2:N) = P_static_global(1:M, 2:N) ./ (R_constant .* T_static_global(1:M, 2:N));
    
    % update entropy
    for x = 2:N
        S(1:M, x) = S(1:M, x - 1) + c_p .* log(h_o_rel(1:M, x) ./ h_o_rel(1:M, x - 1)) - R_constant * log(P_o_rel(1:M, x) ./ P_o_rel(1:M, x - 1));
    end
    
    % update convergence parameters
    eps = norm(Psi_values - Psi_old);
    
    conv = [conv eps];

    stop_condition = eps > 1e-5;
    iteration = iteration + 1;


%% end of loop
end 

semilogy(conv)

%% Questions for Prof:
% 1. for the TE, do you set U = 0 when outside the blade and U =/= 0
% when inside? Currently the prior value (at the TE) will have a U value,
% and therefore the streamline has a U value, although h_o is assumed to be
% constant.
% 2. Confirm (rho*r)1/2 is the average of adjacent cells?




% ------- Step X: Functions ---------

function Psi = calculateInitPsi(r, r_h, r_s)
    % This function calculates Psi based on the radial position r,
    % the hub radius r_h, and the shroud radius r_s.
    Psi = (r^2 - r_h^2) / (r_s^2 - r_h^2);

    % Ensure Psi is zero when r is exactly equal to r_h
    if r == r_h
        Psi = 0;
    end
end


