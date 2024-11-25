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

delta_rC_theta_hub = 82.3; % [m^2/s]
delta_rC_theta_hub = 0;
delta_rC_theta_shroud = 84.4; % [m^2/s]
delta_rC_theta_shroud = 0;

%% ------- Step 1: Grid -------

blade_width = 0.1;      % Blade width (Specified for this system)

% Define the x and r range based on blade width
widths_before_rotor = 4;     % Ensure 4 blade widths before rotor
widths_after_rotor = 4;      % Ensure 4 blade widths before rotor

% Define size of the grid
x_min = 0;               
x_max = (widths_before_rotor + widths_after_rotor + 1) * blade_width;
num_points_rotor = 10;        % Number of grid points within blade
num_points_x = 1+num_points_rotor*(widths_before_rotor+widths_after_rotor+1); % Ensure LE, TE sit exactly on grid point

x_values = linspace(x_min, x_max, num_points_x);
dx = x_values(2) - x_values(1); % Spacing in the x direction
dr = dx; % Spacing in the r direction
num_points_r = 1 + blade_width/dr;


% Generate grid points in x and r directions
r_values = linspace(r_h, r_s, num_points_r);


% Create mesh grid for x and r coordinates
[X, R] = meshgrid(x_values, r_values);

% 2D plot to visualize grid
%{
figure;
hold on;

plot(X, R, 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1);    % Horizontal grid lines
plot(X', R', 'Color', [0.7, 0.7, 0.7], 'LineWidth', 1);  % Vertical grid lines

% Add vertical lines indicating rotor LE and TE
xline(widths_before_rotor*blade_width, 'r--', 'LineWidth', 2);
xline((widths_before_rotor+1)*blade_width, 'r--', 'LineWidth', 2);

% Set labels and title
xlabel('Width');
ylabel('Radius');
title('Grid for Finite Difference Calculation');

% Adjust grid and axis
grid on;
axis equal;
%}


%% ------- Step 2: Variables initiatization -------

% Variables Initialization
P_inlet_static = rho_inlet * R_constant * T_inlet_static;            % Inlet Static Pressure
m_dot = rho_inlet * C_x_inlet * A_inlet;                    % Mass flow rate
H0_inlet = c_p * T_inlet_static + (C_x_inlet^2) / 2;        % Total enthalpy at inlet

% Initialize Psi_values as a 2D array to store Psi(x, r) for each grid point
Psi_values = zeros(num_points_r, num_points_x);
% Nested loop to iterate through each (x, r) point and calculate Psi
for r = 1:num_points_r
    for x = 1:num_points_x
        Psi_values(r, x) = calculateInitPsi(R(r, x), r_h, r_s);
    end
end
clear r

% ------- Calculating Velocities -------
% Initialize matrices to store velocities
C_x = zeros(num_points_r, num_points_x);
C_r = zeros(num_points_r, num_points_x);
% rC_theta = zeros(num_points_r, num_points_x); Initialized later


% -----Calculate C_x and C_r -----
% Treat each corner and perimeter row/column separately...

% Hub, inlet. Modified as there are no points below, or left of, (1, 1)
C_x(1, 1) = m_dot / ((2 * pi * rho_inlet) .* R(1, 1)) * (Psi_values(2, 1) - Psi_values(1, 1)) / (dr);
C_r(1, 1) = -m_dot / (2 * pi * rho_inlet * R(1, 1)) * (Psi_values(1, 2) - Psi_values(1, 1)) / (dx);

% Hub, outlet. Modified as there are no points below, or right of (1, num_points_x)
C_x(1, num_points_x) = m_dot / (2 * pi * rho_inlet * R(1, num_points_x)) * (Psi_values(2, num_points_x) - Psi_values(1, num_points_x)) / (dr);
C_r(1, num_points_x) = -m_dot / (2 * pi * rho_inlet * R(1, num_points_x)) * (Psi_values(1, num_points_x) - Psi_values(1, num_points_x - 1)) / (dx);

% Shroud, inlet. Modified as there are no points above, or left of, (num_points_r, 1)
C_x(num_points_r, 1) = m_dot / (2 * pi * rho_inlet * R(num_points_r, 1)) * (Psi_values(num_points_r, 1) - Psi_values(num_points_r - 1, 1)) / (dr);
C_r(num_points_r, 1) = -m_dot / (2 * pi * rho_inlet * R(num_points_r, 1)) * (Psi_values(num_points_r, 2) - Psi_values(num_points_r, 1)) / (dx);

% Shroud, outlet. Modified as there are no points above, or right of, (num_points_r, num_points_x)
C_x(num_points_r, num_points_x) = m_dot / (2 * pi * rho_inlet * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r - 1, num_points_x)) / (dr); 
C_r(num_points_r, num_points_x) = -m_dot / (2 * pi * rho_inlet * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r, num_points_x - 1)) / (dx);

% Hub. Modified as there are no points below the row
C_x(1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(1, 2:num_points_x-1)) .* (Psi_values(2, 2:num_points_x-1) - Psi_values(1, 2:num_points_x-1)) ./ (dr);
C_r(1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(1, 2:num_points_x-1)) .* (Psi_values(1, 3:num_points_x) - Psi_values(1, 1:num_points_x-2)) ./ (2*dx); 

% Shroud. Modified as there are no points above the row
C_x(num_points_r, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 2:num_points_x-1) - Psi_values(num_points_r - 1, 2:num_points_x-1)) ./ (dr);
C_r(num_points_r, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 3:num_points_x) - Psi_values(num_points_r, 1:num_points_x-2)) ./ (2*dx);

% Inlet. Modified as there are no points left of the column
C_x(2:num_points_r-1, 1) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 1)) .* (Psi_values(3:num_points_r, 1) - Psi_values(1:num_points_r-2, 1)) ./ (2 * dr);
C_r(2:num_points_r-1, 1) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 1)) .* (Psi_values(2:num_points_r-1, 2) - Psi_values(2:num_points_r-1, 1)) ./ (dx);

% Outlet. Modified as there are no points right of the column
C_x(2:num_points_r-1, num_points_x) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(3:num_points_r, num_points_x) - Psi_values(1:num_points_r - 2, num_points_x)) ./ (2*dr);
C_r(2:num_points_r-1, num_points_x) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(2:num_points_r-1, num_points_x) - Psi_values(2:num_points_r-1, num_points_x - 1)) ./ (dx);

% Everywhere else
C_x(2:num_points_r-1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(3:num_points_r, 2:num_points_x-1) - Psi_values(1:num_points_r-2, 2:num_points_x-1)) ./ (2 * dr);
C_r(2:num_points_r-1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_inlet .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(2:num_points_r-1, 3:num_points_x) - Psi_values(2:num_points_r-1, 1:num_points_x-2)) ./ (2 * dx);

% Initialize C_m
C_m = sqrt(C_r.^2 + C_x.^2);

% Initialize rC_theta
rC_theta = zeros(num_points_r, num_points_x);
rC_theta(1:num_points_r, :) = R(1:num_points_r, :) .* C_x_inlet * tan_alpha;

% ----- Initialize Loss Factor -----
% Define the x positions of the leading edge (LE) and trailing edge (TE)
n_cells_before_LE = widths_before_rotor*num_points_rotor;
x_LE = x_max/(num_points_x-1) * n_cells_before_LE;   % x position for LE
x_TE = x_LE + blade_width;                           % x position for TE

% Find the indices in x_values closest to x_LE and x_TE
[~, i_LE] = min(abs(x_values - x_LE));
[~, i_TE] = min(abs(x_values - x_TE));

% Initialize the loss factor matrix w as zeros
loss_global = zeros(num_points_r, num_points_x);

% Calculate the uniform loss factor for the blade section
loss_factor = 0.05 / (i_TE - i_LE);

% Apply the loss factor uniformly between i = LE + 1 and i = TE
loss_global(:, (i_LE + 1):i_TE) = loss_factor;

%% ------- Step 3: Start the rotor -------
 
% Define the parameters
M = num_points_r;  % Number of radial points from hub to shroud
% Calculate the prescribed swirl increase at the trailing edge for each radial station
delta_rC_theta_TE = zeros(M, 1);
for r = 1:M
    delta_rC_theta_TE(r) = delta_rC_theta_hub + (delta_rC_theta_shroud - delta_rC_theta_hub) * (r - 1) / (M - 1);
end

%% Apply increase in whirl and distribute evenly across blade width
for x = (i_LE + 1):i_TE  % Loop over axial points from L.E. to T.E.
    % Linearly interpolate in the x-direction between i_LE and i_TE
    % % Note; delta_rC_theta_LE = 0 based on this interpolation and
    % therefore is not added
    rC_theta(1:M, x) = rC_theta(1:M, x) + delta_rC_theta_TE .* ((x - i_LE) / (i_TE - i_LE));
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
h = zeros(M, N);           % Static enthalpy at each station
T_o_rel = zeros(M, N);     % Relative total temperature at each station
h_o_rel = zeros(M, N);     % Relative total enthalpy at each station
P_o_rel = zeros(M, N);     % Relative total pressure at each station
P_o = zeros(M, N);         % Total pressure at each station
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

%% ------- Step 3.5b: Calculation Block: Calculations -------

% Calculate C_theta and other non-iterated quantities
C_theta = rC_theta ./ R; 

% These seem to sit nicer in the code here than below, but commenting out in case it breaks
V_theta_global(1:M, 1:N) = U(1:M, 1:N) - C_theta(1:M, 1:N);
Beta_global(1:M, 1:N) = atan(V_theta_global(1:M, 1:N) ./ C_m(1:M, 1:N));
V_global(1:M, 1:N) = sqrt(C_m(1:M, 1:N).^2 + V_theta_global(1:M, 1:N).^2);


%% Inlet Calculations, following in-class example
% Calculate T_o1 for each radial position j at entry
T_o1 = T_inlet_static + (C_x(1:M, 1).^2 + C_theta(1:M, 1).^2) / (2 * c_p);
h_o1 = c_p .* T_o1;

% Populate initial column of stagnation temperature and enthalpy
T_o(1:M, 1) = T_o1;
h_o(1:M, 1) = h_o1;

% Initialize T_o1_rel array to store relative total temperature for each radial position j
% % only valid over entire domain as U initialized as 0 outside rotor
T_o1_rel = T_o1 + (U(1:M, 1) .* (U(1:M, 1) - 2.*C_theta(1:M, 1).*U(1:M,1)) / (2 * c_p));
h_o1_rel = c_p .* T_o1_rel;

% populate initial column of relative temp and enthalpy
T_o_rel(:, 1) = T_o1_rel;
h_o_rel(:, 1) = h_o1_rel;

% Calculate P_o1_rel for each radial position j
P_o1_rel = P_inlet_static .* (T_o1_rel ./ T_o1).^k_gamma;

% populate initial column of relative pressure
P_o_rel(:, 1) = P_o1_rel;

P_o_rel_ideal(:, 1) = P_o1_rel; 

% initialize entropy at inlet to be zero
S(:, 1) = 0;  

% Initialize static temperature and pressure at inlet
T_static_global(:, 1) = T_inlet_static;  % Set initial static temperature to inlet static temperature
P_static_global(:, 1) = P_inlet_static;  % Set initial static pressure to inlet static pressure
P_o(:, 1) = P_inlet_static .* (T_o1 ./ T_inlet_static).^(k_gamma);

rho_global(:, 1) = P_static_global(:, 1) ./ (R_constant .* T_static_global(:, 1));


% note that I = H0_rel - 0.5*U^2 is constant throughout (Conservation of Rothalpy)

%% Middle Chunk of Calculations
% Move from 1 -> N updating using common calc block
% note that I and velocities are already calculated

% relative quantities from Conservation of Rothalpy (assuming no c_r for first iteration)
h_o_rel(1:M, 2:N) = h_o1_rel + 0.5 * U(1:M, 2:N).^2; % Conservation of Rothalpy, U(1, 1) = 0

T_o_rel(1:M, 2:N) = h_o_rel(1:M, 2:N) ./ c_p;

% stagnation temperatures from velocities and T_o_rel
T_o(1:M, 2:N) = T_o_rel(1:M, 2:N) - (U(1:M, 2:N) .* (U(1:M, 2:N) - 2.*C_theta(1:M, 2:N)) ./ (2 * c_p));

% constant T_o after TE (no work)
T_o(1:M, (i_TE+1):N) = repmat(T_o(1:M, i_TE), 1, N - i_TE);

% constant T_o = T_o after TE (no work, no U)
% T_o(1:M, (i_TE+1):N) = repmat(T_o(1:M, i_TE), 1, N - i_TE);

% constant T_o_rel = T_o after TE (no work, no U)
T_o_rel(1:M, (i_TE+1):N) = T_o(1:M, (i_TE+1):N);

h_o(1:M, 2:N) = c_p .* T_o(1:M, 2:N);
h_o_rel(1:M, 2:N) = c_p .* T_o_rel(1:M, 2:N);

h(1:M, 1:N) = h_o_rel(1:M, 1:N) - 0.5 .* V_global(1:M, 1:N).^2;

% static properties from stagnation and velocities
T_static_global(1:M, 2:N) = T_o(1:M, 2:N) - (V_global(1:M, 2:N).^2 ./ (2 * c_p));


% finding pressures
for x = 2:N
    % isentropic relative pressure
    P_o_rel_ideal(1:M, x) = P_o_rel_ideal(1:M, x-1) .* ((T_o_rel(1:M, x) ./ T_o_rel(1:M, x-1)).^k_gamma);
    % relative pressure using loss coeff
    P_o_rel(1:M, x) = P_o_rel_ideal(1:M, x) - loss_global(1:M, x) .* (P_o_rel(1:M, x-1) - P_static_global(1:M, x-1));
end

% stagnation and static static pressure from relative and stagnation properties
P_o(1:M, 2:N) = P_o_rel(1:M, 2:N) .* (T_o(1:M, 2:N) ./ T_o_rel(1:M, 2:N)).^(k_gamma);
P_static_global(1:M, 2:N) = P_o(1:M, 2:N) .* (T_static_global(1:M, 2:N) ./ T_o(1:M, 2:N)).^(k_gamma);

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

%{
 % Streamlines
figure;
contourf(X, R, Psi_values, 20);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Streamline Functions');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
% xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

% Entropy plot
figure;
contourf(X, R, S, 20);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Entropy');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

figure;
contourf(X, R, h_o, 20);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Stagnation Enthalpy');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

figure;
contourf(X, R, rho_global, 20);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Density');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bard to indicate value scale

% Static pressure plot
figure;
contourf(X, R, P_static_global, 20);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Static Pressure');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale
%}

while (stop_condition && iteration < max_iter) || iteration <= min_iter
    % if iteration == min_iter
    %     % Streamlines
    %     figure;
    %     contourf(X, R, Psi_values, 10);  % Replace T with the desired output, 20 is the number of contour levels
    %     xlabel('x (Axial Coordinate)');
    %     ylabel('r (Radial Coordinate)');
    %     title('Streamline Functions');
    %     x_lower = (0.5 + widths_before_rotor)*blade_width;
    %     x_upper = x_lower + 2*blade_width;
    %     % xlim([x_lower x_upper]);
    %     colorbar;  % Add a color bar to indicate value scale
    % 
    %     % Stagnation enthalpy plot
    %     figure;
    %     contourf(X, R, h_o, 10);  % Replace T with the desired output, 20 is the number of contour levels
    %     xlabel('x (Axial Coordinate)');
    %     ylabel('r (Radial Coordinate)');
    %     title('Stagnation Enthalpy');
    %     x_lower = (0.5 + widths_before_rotor)*blade_width;
    %     x_upper = x_lower + 2*blade_width;
    %     %xlim([x_lower x_upper]);
    %     colorbar;  % Add a color bar to indicate value scale
    % 
    %     % Static pressure plot
    %     figure;
    %     contourf(X, R, P_static_global, 10);  % Replace T with the desired output, 20 is the number of contour levels
    %     xlabel('x (Axial Coordinate)');
    %     ylabel('r (Radial Coordinate)');
    %     title('Static Pressure');
    %     x_lower = (0.5 + widths_before_rotor)*blade_width;
    %     x_upper = x_lower + 2*blade_width;
    %     %xlim([x_lower x_upper]);
    %     colorbar;  % Add a color bar to indicate value scale
    % 
    %     figure;
    %     % semilogy(conv)
    % 
    % end
        
    % ------- Step 4: Calculate Vorticity -------
    % exclude hub, shroud, inlet = (2:M-1, 2:N)
    
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
    
    %% -> currently assuming the average of adjacent cells ***** CHECK ******. This alignes with instructions
    
    % computing adjacent averages of rho*r
    rho_r_plus_x = (rho_global(2:M-1, 3:N).*R(2:M-1, 3:N) + rho_global(2:M-1, 2:N-1).* R(2:M-1, 2:N-1)) / 2; % average of rho*r at cell and at next cell in x
    rho_r_minus_x = (rho_global(2:M-1, 2:N-1).*R(2:M-1, 2:N-1) + rho_global(2:M-1, 1:N-2).* R(2:M-1, 1:N-2)) / 2; % average of rho*r at cell and at prior cell in x
    rho_r_plus_r = (rho_global(3:M, 2:N-1).*R(3:M, 2:N-1) + rho_global(2:M-1, 2:N-1).* R(2:M-1, 2:N-1)) / 2; % average of rho*r at cell and at next cell in r (towards shroud)
    rho_r_minus_r = (rho_global(2:M-1, 2:N-1).*R(2:M-1, 2:N-1) + rho_global(1:M-2, 2:N-1).* R(1:M-2, 2:N-1)) / 2; % average of rho*r at cell and at prior cell in r (towards hub)
    
    % coefficients
    A = (rho_r_plus_x.^(-1) + rho_r_minus_x.^(-1) + rho_r_plus_r.^(-1) + rho_r_minus_r.^(-1)).^(-1);
    
    B = (Psi_values(2:M-1, 3:N) ./ rho_r_plus_x) + (Psi_values(2:M-1, 1:N-2) ./ rho_r_minus_x) + (Psi_values(3:M, 2:N-1) ./ rho_r_plus_r) + (Psi_values(1:M-2, 2:N-1) ./ rho_r_minus_r);
    
    % update psi
    Psi_values(2:M-1, 2:N-1) = A .* (B + dx^2.*omega(:, 1:N-2));
    
    % update exit
    A_exit = (2 ./ (rho_r_minus_x(:, end)) + rho_r_plus_r(:, end).^(-1) + rho_r_minus_r(:, end).^(-1)).^(-1);
    B_exit = (2.*Psi_values(2:M-1, N-1) ./ rho_r_minus_x(:, end)) + (Psi_values(3:M, N) ./ rho_r_plus_r(:, end)) + (Psi_values(1:M-2, N) ./ rho_r_minus_r(:, end));
    
    Psi_values(2:M-1, N) = A_exit .* (B_exit + dx^2.*omega(:, N-1));
    % note that omega doesn't include leading edge, so N-1 = N (exit)
    
    
    % ------- Step 6: Calculate velocities -------
        % (inlet remains constant)
        % repeat of prior calcs
        
    % Hub, outlet. Modified as there are no points below, or right of, (1, num_points_x)
    C_x(1, num_points_x) = m_dot / (2 * pi * rho_global(1, num_points_x) * R(1, num_points_x)) * (Psi_values(2, num_points_x) - Psi_values(1, num_points_x)) / (dr);
    % Should the denominator be 0.5dr or dr?
    % Why is there a '+' here until end... Suspicious is something is missing?
    C_r(1, num_points_x) = -m_dot / (2 * pi * rho_global(1, num_points_x) * R(1, num_points_x)) * (Psi_values(1, num_points_x) - Psi_values(1, num_points_x - 1)) / (dx);
    
    % Shroud, outlet. Modified as there are no points above, or right of, (num_points_r, num_points_x)
    C_x(num_points_r, num_points_x) = m_dot / (2 * pi * rho_global(num_points_r, num_points_x) * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r - 1, num_points_x)) / (dr);
    C_r(num_points_r, num_points_x) = -m_dot / (2 * pi * rho_global(num_points_r, num_points_x) * R(num_points_r, num_points_x)) * (Psi_values(num_points_r, num_points_x) - Psi_values(num_points_r, num_points_x - 1)) / (dx);
    
    % Hub. Modified as there are no points below the row
    C_x(1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_global(1, 2:num_points_x-1) .* R(1, 2:num_points_x-1)) .* (Psi_values(2, 2:num_points_x-1) - Psi_values(1, 2:num_points_x-1)) ./ (dr);
    C_r(1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_global(1, 2:num_points_x-1) .* R(1, 2:num_points_x-1)) .* (Psi_values(1, 3:num_points_x) - Psi_values(1, 1:num_points_x-2)) ./ (2*dx); 
    
    % Shroud. Modified as there are no points above the row
    C_x(num_points_r, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_global(num_points_r, 2:num_points_x-1) .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 2:num_points_x-1) - Psi_values(num_points_r - 1, 2:num_points_x-1)) ./ (dr);
    C_r(num_points_r, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_global(num_points_r, 2:num_points_x-1) .* R(num_points_r, 2:num_points_x-1)) .* (Psi_values(num_points_r, 3:num_points_x) - Psi_values(num_points_r, 1:num_points_x-2)) ./ (2*dx);
    
    % Outlet. Modified as there are no points right of the column
    C_x(2:num_points_r-1, num_points_x) = m_dot ./ (2 * pi * rho_global(2:num_points_r-1, num_points_x) .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(3:num_points_r, num_points_x) - Psi_values(1:num_points_r - 2, num_points_x)) ./ (2*dr);
    C_r(2:num_points_r-1, num_points_x) = -m_dot ./ (2 * pi * rho_global(2:num_points_r-1, num_points_x) .* R(2:num_points_r-1, num_points_x)) .* (Psi_values(2:num_points_r-1, num_points_x) - Psi_values(2:num_points_r-1, num_points_x - 1)) ./ (dx);
    
    % Everywhere else
    C_x(2:num_points_r-1, 2:num_points_x-1) = m_dot ./ (2 * pi * rho_global(2:num_points_r-1, 2:num_points_x-1) .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(3:num_points_r, 2:num_points_x-1) - Psi_values(1:num_points_r-2, 2:num_points_x-1)) ./ (2 * dr);
    C_r(2:num_points_r-1, 2:num_points_x-1) = -m_dot ./ (2 * pi * rho_global(2:num_points_r-1, 2:num_points_x-1) .* R(2:num_points_r-1, 2:num_points_x-1)) .* (Psi_values(2:num_points_r-1, 3:num_points_x) - Psi_values(2:num_points_r-1, 1:num_points_x-2)) ./ (2 * dx);

    % calculate C_m
    C_m = sqrt(C_r.^2 + C_x.^2);
    
    % ------- Step 7: Calculation Block: Calculations -------> Start at TE and trace backwards to LE
    % swirl remains the same, imposed by work distribution
    
    % interpolation factors % [(M - 2) x (N - 1)]
    % % easier as psi decreases from shroud -> hub
    % r = 0 at hub, r = M at shroud
    
    %% Interpolation - Tracing streamline backwards
    % streamlines at hub and shroud not calculated, as straight lines
        % a is interpolation factor closer to shroud
        % b is interpolation factor closer to hub

    a_towards_shroud = (Psi_values(2:M-1, 2:N) - Psi_values(1:M-2, 1:N-1)) ./ (Psi_values(2:M-1, 1:N-1) - Psi_values(1:M-2, 1:N-1));
    shroud_side = ~(a_towards_shroud > 1);
    a_towards_shroud(~shroud_side) = NaN; % removing invalid interpolations (streamline goes other direction)
    b_towards_shroud = 1 - a_towards_shroud;
    
    a_towards_hub = (Psi_values(2:M-1, 2:N) - Psi_values(2:M-1, 1:N-1)) ./ (Psi_values(3:M, 1:N-1) - Psi_values(2:M-1, 1:N-1));
    hub_side = ~(a_towards_hub <= 0);
    a_towards_hub(~hub_side) = NaN; % removing invalid interpolations (streamline goes other direction)
    b_towards_hub = 1 - a_towards_hub;

    
    
    %% Note: remember to swap stag / rel properties at LE and TE depending on the cells (slide 23)
    
    %% Note: There has got to be a cleaner way to do this -_-

    %% HUB / SHROUD
    for x = 2:N % from inlet to exit (inlet already defined)
        if x <= i_LE || x > i_TE % for values before LE or after TE, use stagnation enthalpy

            % hub
            h_o(1, x) = h_o(1, x-1); % h_o at origin (i-1)
            h_o_rel(1, x) = h_o(1, x) + U(1, x) * (U(1, x) - 2*C_theta(1, x)) / 2;

            % shroud
            h_o(M, x) = h_o(M, x-1); % h_o at origin (i-1)
            h_o_rel(M, x) = h_o(M, x) + U(M, x) * (U(M, x) - 2*C_theta(M, x)) / 2;


        else % relative enthalpies
            % hub
            h_o_rel(1, x) = h_o_rel(1, x-1);
            h_o(1, x) = h_o_rel(1, x) - U(1, x) * (U(1, x) - 2*C_theta(1, x)) / 2;

            % shroud
            h_o_rel(M, x) = h_o_rel(M, x-1);
            h_o(M, x) = h_o_rel(M, x) - U(M, x) * (U(M, x) - 2*C_theta(M, x)) / 2;

        end
    end


    
    for x = 2:N % from inlet to exit (inlet already defined)
        for r = 2:M-1
    
            if x <= i_LE || x > i_TE % for values before LE or after TE, use stagnation enthalpy
                
                if shroud_side(r-1, x-1) % if streamline goes towards shroud
                    h_o_streamline = h_o(r, x-1).* a_towards_shroud(r-1, x-1) + h_o(r-1, x-1).* b_towards_shroud(r-1, x-1); % h_o at origin (i-1)
                      
                else % if streamline goes towards hub
                    h_o_streamline = h_o(r+1, x-1).* a_towards_hub(r-1, x-1) + h_o(r, x-1).* b_towards_hub(r-1, x-1); % h_o at origin (i-1)

                end

                U_streamline = 0;

                % trace forward to define enthalpy (i)
                h_o(r, x) = h_o_streamline; % - (0.5 * U_streamline^2) + (0.5 * U(r, x)^2); % defining same point, just from different cell depending on streamline

                % update relative enthalpy
                h_o_rel(r, x) = h_o(r, x) + U(r, x) * (U(r, x) - 2*C_theta(r, x)) / 2;
        
            else % for values inside blade, use relative enthalpy
    
                if shroud_side(r-1, x-1) % if streamline goes towards shroud
                    h_o_rel_streamline = h_o_rel(r, x-1).* a_towards_shroud(r-1, x-1) + h_o_rel(r-1, x-1).* b_towards_shroud(r-1, x-1); % h_o_rel at origin (i-1)
                    
                 % calculate other properties at streamline origin
                    U_streamline = U(r, x-1).* a_towards_shroud(r-1, x-1) + U(r-1, x-1).* b_towards_shroud(r-1, x-1); % U at origin (non-zero)
                 
                else % if streamline goes towards hub
                    h_o_rel_streamline = h_o_rel(r+1, x-1).* a_towards_hub(r-1, x-1) + h_o_rel(r, x-1).* b_towards_hub(r-1, x-1); % h_o_rel at origin (i-1)
                    
                    % calculate other properties at streamline origin
                    U_streamline = U(r+1, x-1).* a_towards_hub(r-1, x-1) + U(r, x-1).* b_towards_hub(r-1, x-1); % U at origin (i-1) (non-zero)
    
                end
    
                % trace forward to define relative enthalpy (i)
                h_o_rel(r, x) = h_o_rel_streamline - (0.5 * U_streamline^2) + (0.5 * U(r, x)^2); % defining same point, just from different nodes depending on streamline
    
                % update stagnation enthalpy
                h_o(r, x) = h_o_rel(r, x) - U(r, x) * (U(r, x) - 2*C_theta(r, x)) / 2;
            end
        end
    end
    
    %% Common calculation block from h_o / h_o_rel & velocities
    
    % V_theta_global is constant (U & C_theta don't change)
    Beta_global(1:M, 2:N) = atan(V_theta_global(1:M, 2:N) ./ C_m(1:M, 2:N));
    V_global(1:M, 2:N) = sqrt(C_m(1:M, 2:N).^2 + V_theta_global(1:M, 2:N).^2);
    
    % stagnation (rel & abs) temperatures from h_o & h_o_rel
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
    
    end

    % stagnation and static static pressure from relative and stagnation properties
    P_o(1:M, 2:N) = P_o_rel(1:M, 2:N) .* (T_o(1:M, 2:N) ./ T_o_rel(1:M, 2:N)).^(k_gamma);
    P_static_global(1:M, 2:N) = P_o(1:M, 2:N) .* (T_static_global(1:M, 2:N) ./ T_o(1:M, 2:N)).^(k_gamma);
    
    % density throughout
    rho_global(1:M, 2:N) = P_static_global(1:M, 2:N) ./ (R_constant .* T_static_global(1:M, 2:N));
    
    % update entropy
    for x = 2:N
        S(1:M, x) = S(1:M, x - 1) + c_p .* log(h_o_rel(1:M, x) ./ h_o_rel(1:M, x - 1)) - R_constant * log(P_o_rel(1:M, x) ./ P_o_rel(1:M, x - 1));
    end
    
    % update convergence parameters
    eps = norm(Psi_values - Psi_old);
    eps = norm(S - S_old);
    
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


%{
[X, R] = meshgrid(x_values, r_values);  % Create grid of x and r values
% Surface plot
figure;% Contour plot
figure;
contourf(X, R, T, 20);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Temperature Contour Plot');
colorbar;  % Add a color bar to indicate value scale
surf(X, R, V_global);  % Replace T with the desired output (temperature, pressure, etc.)
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
zlabel('Temperature');  % Replace with the variable you are plotting
title('Temperature Distribution');
shading interp;  % Smooth the color transition
colorbar;  % Add a color bar to indicate value scale
%}

%% ----- Plots -----
% Streamlines
figure;
contourf(X, R, Psi_values, 10);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Streamline Functions');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

% Stagnation enthalpy plot
figure;
contourf(X, R, h_o, 10);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Stagnation Enthalpy');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

% Static pressure plot
figure;
contourf(X, R, P_static_global, 10);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Static Pressure');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

% C_x
figure;
contourf(X, R, C_x, 10);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('C_x');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

% C_r
figure;
contourf(X, R, C_r, 10);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('C_r');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale


figure;
contourf(X, R, C_m, 10);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('C_m');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale

% Entropy plot
figure;
contourf(X, R, S, 20);  % Replace T with the desired output, 20 is the number of contour levels
xlabel('x (Axial Coordinate)');
ylabel('r (Radial Coordinate)');
title('Entropy');
x_lower = (0.5 + widths_before_rotor)*blade_width;
x_upper = x_lower + 2*blade_width;
%xlim([x_lower x_upper]);
colorbar;  % Add a color bar to indicate value scale


% ------- Step X: Functions ---------

function plot_contour(name, X, R, f)
    figure;
    n = 10;
    contourf(X, R, f, n);  % Replace T with the desired output, n is the number of contour levels
    xlabel('x (Axial Coordinate)');
    ylabel('r (Radial Coordinate)');
    title(name);
    colorbar;

end

function Psi = calculateInitPsi(r, r_h, r_s)
    % This function calculates Psi based on the radial position r,
    % the hub radius r_h, and the shroud radius r_s.
    Psi = (r^2 - r_h^2) / (r_s^2 - r_h^2);

    % Ensure Psi is zero when r is exactly equal to r_h
    if r == r_h
        Psi = 0;
    end
end


