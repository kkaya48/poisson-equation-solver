%getting started, problem parameters
L = 1;
N = input('Enter the value of N:');
sigma = input('Enter the value of sigma: ');

dx = L / (N - 1);
dy = L / (N - 1);

x = linspace(0, L, N);
y = linspace(0, L, N);

[X, Y] = meshgrid(x, y); % meshgrid creates the 2D coordinate matrices X and Y from the x and y grid vectors

inside = (X + Y) <= L; % inside: logical mask for grid points inside the triangular domain (true/false)

disp(inside);

figure;
plot(X(inside), Y(inside), 'bo');
axis equal;
grid on;
xlabel('x');
ylabel('y');
title('Points inside the triangular domain');

bottom = inside & (Y == 0);
left = inside & (X == 0);
diagonal = inside & ((X + Y) == L);
interior = inside & ~bottom & ~left & ~diagonal;

% testing: visualize boundary and interior point classificationne
figure;
hold on;

plot(X(bottom), Y(bottom), 'ro');
plot(X(left), Y(left), 'go');
plot(X(diagonal), Y(diagonal), 'mo');
plot(X(interior), Y(interior), 'bo');

axis equal;
grid on;
xlabel('x');
ylabel('y');
title('Boundary and interior points');
legend('bottom boundary', 'left boundary', 'diagonal boundary', 'interior');

hold off;

% index formula from the assignment
mIndex = zeros (N, N);

for i = 1:N
    for j= 1:(N - i + 1)
        mIndex(j,i) = j + (i -1) * (N - i/2 + 1);
    end
end

disp(mIndex);

% testing the computation time
tic

% Build the linear system A*u = b for f(x,y) = sigma
M = N * (N + 1) / 2;
A = zeros(M, M);
b = zeros(M, 1);

for i = 1:N
    for j = 1:(N - i + 1)
        m = mIndex(j, i);

        % Boundary conditions
        if j == N - i + 1
            % diagonal boundary: x + y = L
            A(m, m) = 1;
            b(m) = 0;
            
        elseif i == 1
            % left boundary: x = 0
            A(m, m) = 1;
            b(m) = -1;
        
        elseif j == 1
            % bottom boundary: y = 0
            A(m, m) = 1;
            b(m) = 1;
        else
            % interior point: finite difference equation
            m_right = mIndex(j, i + 1);
            m_left = mIndex(j, i - 1);
            m_up = mIndex(j+1, i);
            m_down = mIndex(j - 1, i);

            A(m, m) = -2 / dx^2 - 2 / dy^2;
            A(m, m_right) = 1 / dx^2;
            A(m, m_left)  = 1 / dx^2;
            A(m, m_up)    = 1 / dy^2;
            A(m, m_down)  = 1 / dy^2;

            b(m) = -sigma;
        end
    end
end

u = A \ b;

U = NaN(N, N);

for i = 1:N
    for j = 1:(N - i + 1)
        m = mIndex(j, i);
        U(j, i) = u(m);
    end
end

elapsedTime = toc;
disp(['Computation time: ', num2str(elapsedTime), ' seconds'])

figure;
surf(X, Y, U);
xlabel('x');
ylabel('y');
zlabel('u(x,y)');
title(['Numerical Solution of the Poisson Problem for N = ', num2str(N),', sigma = ',num2str(sigma),', comp time = ',num2str(elapsedTime),' s']);
grid on;