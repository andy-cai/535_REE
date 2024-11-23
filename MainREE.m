%{
Main
Authors:
Description: This script applies the finite difference method to
numerically computer derivatives in the Radial Equilibrium Equation.

Inputs: 

Outputs:


%}

% Define the geometry of the axial flow compressor stage. All units in SI
rotor_width = 0.1;          % x-dimension 
d_hub = 0.9;                % limit of y-dimension
d_tip = 1;                  % limit of y-dimension
r_hub = d_hub/2;
r_tip = d_tip/2;
N = 6000*pi/30;
c_x_in = 136;
rho_in = 1.5;
T_in = 25 + 273.15;
loss_coeff = 0.05;
c_p = 1005; % Assume constant. This might be illegal
p_in = 100000;
gas_const = 287;
A_in = pi*r_tip^2-r_hub^2;

alpha = 30*pi/180;   %Unclear which angle this corresponds to tbh

% beta_1_h = 55.15*pi/180;
% alpha_1_t = 30*pi/180;
% beta_1_t = 60*pi/180;



%% 1. Generate the grid
% Create square matrix to store properties at different geometries
n = 14; % Number of divisions in square
x = linspace(0, rotor_width, n+1);
r = linspace(r_hub, r_tip, n+1);
[X, Y] = meshgrid(x,r);     %ndgrid() might end up being better


%% Plots
% 2D plot
figure;
plot(X, Y, 'k.', 'MarkerSize', 10);  % Plot each point with black dots
hold on;
plot(X', Y', 'k.', 'MarkerSize', 10); % Plot the transpose to show vertical lines
xlabel('Width');
ylabel('Radius');
title('Grid for Finite Difference Calculation');
grid on;
axis equal;
plot(X, Y, 'k');    % Horizontal grid lines
plot(X', Y', 'k');  % Vertical grid lines

% % Convert to cylindrical coordinates to create an annular 3D shape
% theta = X * (2 * pi / rotor_width); % Mapping the width to an angular range from 0 to 2*pi
% R = Y;  % Radial distance stays as Y values
% % Convert to Cartesian coordinates for 3D plotting
% [X_cart, Y_cart] = pol2cart(theta, R);
% % Plot the phi distribution in 3D as an annular surface
% figure;
% surf(X_cart, Y_cart, phi, 'EdgeColor', 'none');  % Create the surface plot, colored by phi values
% colorbar;  % Add a color bar to show the scale of phi
% xlabel('X-axis (Cartesian)');
% ylabel('Y-axis (Cartesian)');
% zlabel('\phi (Phi)');
% title('3D Annular Plot of \phi Distribution');
% axis equal;
% view(3);  % Set the view to 3D


%% Initialize all variables
denom = r_tip^2 - r_hub^2;
% Note that using Y. instead of r. is necessary. Y. yield a matrix, while r. gives a vector
psi = (Y.^2 - r_hub^2)/(denom);
rho_in = p_in / gas_const / T_in;
m = rho_in * c_x_in * A_in;
Stag_enthalp_in = c_p*T_in + c_x_in^2/2;

% Create matrices for c_x, c_r 
c_x = zeros(n,n);
c_r = zeros(n,n);
% Populate velocities for when x = 0
for i = 1:n
    c_x(i, 1) = c_x_in;
end

del_r = 0.1/n;
%for i = 1:n-1
 %   c_x(i+1, 1) = 1/(r(i)*del_r)*psi()
%end

%Graph psi
figure;
surf(X, Y, psi);
xlabel('X-axis');
ylabel('R');
zlabel('\phi (Phi)');
title('Initial \psi Distribution');