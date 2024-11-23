% MATLAB script for solving flow in an axial compressor using the Finite Difference Method (FDM)

% Step 1: Define parameters and grid
N_x = 10; % Number of grid points in axial direction (including LE and TE)
N_r = 10; % Number of grid points in radial direction (hub to shroud)

x = linspace(0, 1, N_x); % Axial coordinate grid (normalized from LE to TE)
r = linspace(0.9, 1, N_r); % Radial coordinate grid (hub to shroud)

% Step 2: Initialize variables
psi = zeros(N_r, N_x); % Stream function
rho = ones(N_r, N_x) * 1.5; % Density initialized to 1.5 kg/m^3
H0 = ones(N_r, N_x) * (25 + 273.15); % Stagnation enthalpy initialized based on inlet temp (K)

% Initial velocities
Cx_const = 136; % Constant axial velocity (m/s)
Cr = zeros(N_r, N_x); % Initialize radial velocity to zero
Ctheta = zeros(N_r, N_x); % Initialize tangential velocity to zero
rCtheta = zeros(N_r, N_x); % Initialize rCtheta to zero

% Step 3: Impose initial swirl conditions
rCtheta_TE = linspace(82.3, 84.4, N_r); % Swirl imposed at trailing edge (hub to shroud)
for j = 1:N_r
    for i = 2:N_x % Inside blade (from LE+1 to TE)
        rCtheta(i, j) = rCtheta(i-1, j) + (rCtheta_TE(j) - 82.3) * (i-1) / (N_x-1);
    end
end

% Step 4: Iterate and solve for vorticity and stream function
omega = zeros(N_r, N_x); % Initialize vorticity
epsilon_psi = 1e-5; % Convergence criterion for stream function
max_iterations = 100;

for iter = 1:max_iterations
    % Update vorticity (omega) using finite differences
    for j = 2:N_r-1
        for i = 2:N_x-1
            omega(i, j) = (rho(i, j) * Cx_const * (psi(i, j+1) - psi(i, j-1)) / (2 * (r(j+1) - r(j-1))) ...
                         + (H0(i, j+1) - H0(i, j-1)) / (2 * (r(j+1) - r(j-1))));
        end
    end
    
    % Update stream function (psi) node-by-node
    for j = 2:N_r-1
        for i = 2:N_x-1
            A = 1 / ((r(i+1) - r(i)) / 2) + 1 / ((r(i-1) - r(i)) / 2);
            B = (psi(i+1, j) + psi(i-1, j) + psi(i, j+1) + psi(i, j-1)) / 4;
            psi(i, j) = A * B;
        end
    end
    
    % Check convergence
    if max(abs(psi(:) - omega(:))) < epsilon_psi
        break;
    end
end

% Step 5: Calculate velocities based on updated stream function
for j = 1:N_r
    for i = 2:N_x-1
        Cr(i, j) = - (psi(i+1, j) - psi(i-1, j)) / (2 * (x(i+1) - x(i-1)));
        Ctheta(i, j) = rCtheta(i, j) / r(j);
    end
end

% Step 6: Display results
fprintf('Number of iterations for convergence: %d\n', iter);
mesh(x, r, psi);
title('Stream Function Distribution');
xlabel('Axial Coordinate (x)');
ylabel('Radial Coordinate (r)');
zlabel('Stream Function (Psi)');
