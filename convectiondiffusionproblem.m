% Define parameters
L = 1.0;         % Length of domain
rho = 1.0;       % Density
Gamma = 0.1;     % Diffusion coefficient

% Define grid
n_cells = 5;     % Number of cells
n_points = n_cells + 1;  % Number of grid points
dx = L / n_cells; % Grid spacing
x = linspace(0, L, n_points); % Grid points

% Boundary conditions
phi_0 = 1;       % Boundary condition at x = 0
phi_L = 0;       % Boundary condition at x = L

% Define cases
cases = [0.1, 2.5]; % velocities

% Analytical solution function
function phi_analytical = analytical_solution(u, x, L, rho, Gamma, phi_0, phi_L)
    phi_analytical = phi_0 + (phi_L - phi_0) * (exp(rho * u * x / Gamma) - 1) / (exp(rho * u * L / Gamma) - 1);
endfunction

% Function to solve the problem using upwind scheme
function phi = solve_convection_diffusion_upwind(u, n_points, dx, rho, Gamma, phi_0, phi_L)
    % Create matrix A and vector B
    A = zeros(n_points - 2, n_points - 2);
    B = zeros(n_points - 2, 1);

    % Upwind scheme
    if u > 0
        % For positive velocities (upwind from the left)
        for i = 1:(n_points - 2)
            A(i, i) = -2 * (Gamma / dx^2) - u / dx;
            if i > 1
                A(i, i-1) = Gamma / dx^2 + u / dx;
            end
            if i < (n_points - 2)
                A(i, i+1) = Gamma / dx^2;
            end
        end
    else
        % For negative velocities (upwind from the right)
        for i = 1:(n_points - 2)
            A(i, i) = -2 * (Gamma / dx^2) + u / dx;
            if i > 1
                A(i, i-1) = Gamma / dx^2;
            end
            if i < (n_points - 2)
                A(i, i+1) = Gamma / dx^2 - u / dx;
            end
        end
    end

    % Boundary conditions
    B(1) = -phi_0 * (Gamma / dx^2 - u / dx);
    B(end) = -phi_L * (Gamma / dx^2 + u / dx);

    % Solve the system
    phi_interior = A \ B;

    % Add boundary conditions
    phi = [phi_0; phi_interior; phi_L];
endfunction

% Function to solve the problem using central differencing scheme
function phi = solve_convection_diffusion_central(u, n_points, dx, rho, Gamma, phi_0, phi_L)
    % Create matrix A and vector B
    A = zeros(n_points - 2, n_points - 2);
    B = zeros(n_points - 2, 1);

    % Central differencing scheme
    for i = 1:(n_points - 2)
        A(i, i) = -2 * (Gamma / dx^2) - u / (2 * dx);
        if i > 1
            A(i, i-1) = Gamma / dx^2 + u / (2 * dx);
        end
        if i < (n_points - 2)
            A(i, i+1) = Gamma / dx^2 - u / (2 * dx);
        end
    end

    % Boundary conditions
    B(1) = -phi_0 * (Gamma / dx^2 - u / (2 * dx));
    B(end) = -phi_L * (Gamma / dx^2 + u / (2 * dx));

    % Solve the system
    phi_interior = A \ B;

    % Add boundary conditions
    phi = [phi_0; phi_interior; phi_L];
endfunction

% Loop through each case and calculate
for i = 1:length(cases)
    u = cases(i);

    % Solve for 5 nodes using upwind scheme
    phi_5_nodes_upwind = solve_convection_diffusion_upwind(u, n_points, dx, rho, Gamma, phi_0, phi_L);

    % Analytical solution for 5 nodes
    phi_analytical_5_nodes = analytical_solution(u, x, L, rho, Gamma, phi_0, phi_L);

    % Solve for 5 nodes using central differencing scheme
    phi_5_nodes_central = solve_convection_diffusion_central(u, n_points, dx, rho, Gamma, phi_0, phi_L);

    % Display results
    disp(['Case ', num2str(i), ' (u = ', num2str(u), ' m/s)']);
    disp('Numerical solution with 5 nodes (Upwind scheme):');
    disp(phi_5_nodes_upwind);
    disp('Numerical solution with 5 nodes (Central differencing scheme):');
    disp(phi_5_nodes_central);
    disp('Analytical solution with 5 nodes:');
    disp(phi_analytical_5_nodes);

    % Plotting
    figure;
    plot(x, phi_5_nodes_upwind, 'ro-', x, phi_5_nodes_central, 'g*-', x, phi_analytical_5_nodes, 'b*-');
    title(['Comparison for u = ', num2str(u), ' m/s']);
    legend('Upwind Scheme', 'Central Differencing', 'Analytical');
    xlabel('x (m)');
    ylabel('φ');
    grid on;
end

% Recalculate with 20 nodes for Case 2
n_cells_20 = 20;
n_points_20 = n_cells_20 + 1;
dx_20 = L / n_cells_20;
x_20 = linspace(0, L, n_points_20);

% Solve for 20 nodes using upwind scheme
phi_20_nodes_upwind = solve_convection_diffusion_upwind(cases(2), n_points_20, dx_20, rho, Gamma, phi_0, phi_L);

% Analytical solution for 20 nodes
phi_analytical_20_nodes = analytical_solution(cases(2), x_20, L, rho, Gamma, phi_0, phi_L);

% Plotting for 20 nodes
figure;
plot(x_20, phi_20_nodes_upwind, 'ro-', x_20, phi_analytical_20_nodes, 'b*-');
title('Comparison for u = 2.5 m/s with 20 nodes');
legend('Numerical (Upwind)', 'Analytical');
xlabel('x (m)');
ylabel('φ');
grid on;

