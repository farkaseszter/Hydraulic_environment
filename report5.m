clear all
% Given Parameters
Qp    = 128.0;    % Flow rate (m^3/s)
B     = 10.0;     % Channel width (m)
fi    = 1.0/250;  % Slope
fn    = 0.025;    % Manning's n
arfa  = 1.1;      % Empirical coefficient
error = 0.0001;   % Convergence tolerance

% Constants
C1 = fn * Qp / (fi^0.5 * B);
C  = C1^1.5;
bc = 2.0 * C / B;
g  = 9.8;

% Initial guess for normal depth
h0 = C1 * 0.6;

% Newton-Raphson Method to find normal depth
while true
    f    = h0^2.5 - bc * h0 - C;
    dfdh = 2.5 * h0^1.5 - bc;   % Derivative of f with respect to h0
    h0new = h0 - f / dfdh;
    e = abs((h0new - h0) / h0);
    
    fprintf('Iter: h0=%.6f, f=%.6f, dfdh=%.6f, h0new=%.6f\n', h0, f, dfdh, h0new);
    
    h0 = h0new;
    
    if e < error
        break;
    end
end

% Velocity calculation
R = h0 / (1.0 + 2.0 * h0 / B);
v = R^(2/3) * sqrt(fi) / fn;

% Critical depth calculation
hc = (arfa * Qp^2 / (g * B^2))^(1/3);

% Output results
fprintf('\nResults:\n');
fprintf('Normal depth h0 = %.6f m\n', h0);
fprintf('Hydraulic radius R = %.6f m\n', R);
fprintf('Velocity v = %.6f m/s\n', v);
fprintf('Critical depth hc = %.6f m\n', hc);
