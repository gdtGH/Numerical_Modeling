% GodunovFlux.m
% Computes the numerical flux at a single interface x_{i+1/2} for the
% 1D nonlinear Shallow Water Equations using the Roe approximate
% Riemann solver.
%
% The SWE in conservative form:
%   U_t + F(U)_x = 0
%   U = [h; q],   F(U) = [q; q^2/h + 0.5*g*h^2]
%
% INPUT:
%   U_L  (2x1)  left  state  [h_L; q_L]
%   U_R  (2x1)  right state  [h_R; q_R]
%   g    (scalar) gravitational acceleration
%
% OUTPUT:
%   Fnum (2x1)  numerical flux at the interface

function Fnum = GodunovFlux(U_L, U_R, g)

%--------------------------------------------------------------------------
% Extract left and right states — guard h > 0
%--------------------------------------------------------------------------
h_L = max(U_L(1), eps);
q_L = U_L(2);
h_R = max(U_R(1), eps);
q_R = U_R(2);

u_L = q_L / h_L;   % velocity left
u_R = q_R / h_R;   % velocity right

%--------------------------------------------------------------------------
% Physical fluxes F(U_L) and F(U_R)
%--------------------------------------------------------------------------
F_L = [q_L;
       q_L^2 / h_L + 0.5 * g * h_L^2];

F_R = [q_R;
       q_R^2 / h_R + 0.5 * g * h_R^2];

%--------------------------------------------------------------------------
% Roe-averaged quantities
%--------------------------------------------------------------------------
sL = sqrt(h_L);
sR = sqrt(h_R);

% Roe-averaged velocity
u_hat = (sL * u_L + sR * u_R) / (sL + sR);

% Roe-averaged celerity
c_hat = sqrt(0.5 * g * (h_L + h_R));

%--------------------------------------------------------------------------
% Eigenvalues of the Roe matrix
%--------------------------------------------------------------------------
lambda_1 = u_hat - c_hat;
lambda_2 = u_hat + c_hat;

%--------------------------------------------------------------------------
% Right eigenvectors of the Roe matrix
%   R_1 = [1;  u_hat - c_hat]
%   R_2 = [1;  u_hat + c_hat]
%--------------------------------------------------------------------------
R_1 = [1; lambda_1];
R_2 = [1; lambda_2];

%--------------------------------------------------------------------------
% Wave strengths (alpha) — projection of dU = U_R - U_L onto eigenvectors
%   Solve: alpha_1*R_1 + alpha_2*R_2 = dU
%   dh = alpha_1 + alpha_2
%   dq = alpha_1*lambda_1 + alpha_2*lambda_2
%--------------------------------------------------------------------------
dh = h_R - h_L;
dq = q_R - q_L;

% From the 2x2 system:
%   [1,        1      ] [alpha_1]   [dh]
%   [lambda_1, lambda_2] [alpha_2] = [dq]
denom = lambda_2 - lambda_1;

% Avoid degenerate case (should not occur for valid h > 0)
if abs(denom) < eps
    denom = eps;
end

alpha_1 = ( lambda_2 * dh - dq) / denom;
alpha_2 = (-lambda_1 * dh + dq) / denom;

%--------------------------------------------------------------------------
% Roe flux: F_roe = 0.5*(F_L + F_R) - 0.5 * sum( |lambda_k| * alpha_k * R_k )
%--------------------------------------------------------------------------
Fnum = 0.5 * (F_L + F_R) ...
     - 0.5 * (abs(lambda_1) * alpha_1 * R_1 ...
            + abs(lambda_2) * alpha_2 * R_2);

end
