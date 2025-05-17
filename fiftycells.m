% Define parameters
L = 0.5; % Length of the rod (m)
N = 50; % Number of cells
dx = L / N; % Cell size (m)
k = 1; % Thermal conductivity (W/m.K)
q = 0; % Heat generation (W/m^3)

% Boundary conditions
T_left = 100; % Left boundary temperature (C)
T_right = 300; % Right boundary temperature (C)

% Initialize temperature array
T = zeros(N+2, 1);

% Set boundary conditions
T(1) = T_left;
T(end) = T_right;

% Coefficient matrix and RHS vector
A = zeros(N, N);
b = zeros(N, 1);

% Fill the matrix A and vector b
for i = 1:N
    if i == 1
        A(i,i) = 2 * k / dx^2;
        A(i,i+1) = -k / dx^2;
        b(i) = q + T_left * k / dx^2;
    elseif i == N
        A(i,i-1) = -k / dx^2;
        A(i,i) = 2 * k / dx^2;
        b(i) = q + T_right * k / dx^2;
    else
        A(i,i-1) = -k / dx^2;
        A(i,i) = 2 * k / dx^2;
        A(i,i+1) = -k / dx^2;
        b(i) = q;
    end
end

% Solve the linear system
T_interior = A \ b;

% Combine the interior and boundary temperatures
T(2:end-1) = T_interior;
% Print numerical solutions
disp('Numerical Solution:');
for i = 1:10:N+2
    fprintf('T(%d) = %f °C\n', i, T(i));
end
% Plot the results
x = linspace(0, L, N+2);
plot(x, T, '-o');
xlabel('Position (m)');
ylabel('Temperature (°C)');
title('Temperature Distribution in the Plate');
legend('Numerical Solution');
grid on;

