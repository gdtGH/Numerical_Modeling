function verify_patch
%% verify_patch — unit checks for the consistent-mass / quadrature patch
%==========================================================================
% TEST A : local P1 mass matrix assembled with Quadrature(3) equals the
%          consistent matrix  (h/6)*[2 1; 1 2]  to within 1e-12.
% TEST B : Cavalieri-Simpson [Quadrature(3)] integrates x^2 on [0,1]
%          (BJ = 1) exactly to 1/3.
%
% These reproduce exactly the assembly path used by Assembly/Matrix1D.m
% (ShapeBasis -> Quadrature -> EvalShapeBasis -> GetJacobian -> Mass).
%==========================================================================

clc;
addpath FESpace
addpath Assembly

tol = 1e-12;
fprintf('============================================================\n');
fprintf(' verify_patch — consistent mass / quadrature checks\n');
fprintf('============================================================\n\n');

%% ---------- shared FE objects ----------
basis            = ShapeBasis;                  % P1 shape functions
[nodes_1D, w_1D] = Quadrature(3);               % Cavalieri-Simpson
[Phi, ~]         = EvalShapeBasis(basis, nodes_1D);
nln              = 2;

fprintf('Quadrature(3) nodes   = [%g %g %g]\n', nodes_1D);
fprintf('Quadrature(3) weights = [%g %g %g]  (sum = %.15g)\n\n', ...
        w_1D(1), w_1D(2), w_1D(3), sum(w_1D));

%% ================= TEST A : local mass matrix =================
hList  = [0.37, 1.0, 2.5];     % a few element lengths
passA  = true;
maxErrA = 0;

fprintf('--- TEST A : local mass matrix M_loc = (h/6)*[2 1; 1 2] ---\n');
for h = hList
    coord      = [0; h];                         % 1 element, length h
    [BJ, ~]    = GetJacobian(coord, nodes_1D);   % BJ = h
    M_loc      = Mass(Phi, w_1D, nln, BJ, ones(numel(w_1D),1));
    M_expected = (h/6) * [2 1; 1 2];
    errH       = max(abs(M_loc(:) - M_expected(:)));
    maxErrA    = max(maxErrA, errH);
    passA      = passA && (errH < tol);
    fprintf('  h = %5.3f :  max|M_loc - (h/6)[2 1;1 2]| = %.3e\n', h, errH);
end
fprintf('  M_loc (h = 1) =\n');
M1 = Mass(Phi, w_1D, nln, 1.0, ones(numel(w_1D),1));
disp(M1);
fprintf('  TEST A : %s   (max error = %.3e)\n\n', ...
        ternary(passA, 'PASS', 'FAIL'), maxErrA);

%% ================= TEST B : integrate x^2 on [0,1] =================
BJ   = 1;                       % GetJacobian uses full element length
fq   = nodes_1D.^2;             % f(x) = x^2 at quadrature nodes
Iq   = sum(w_1D(:) .* BJ .* fq(:));
errB = abs(Iq - 1/3);
passB = errB < tol;

fprintf('--- TEST B : Quadrature(3) integral of x^2 on [0,1] ---\n');
fprintf('  computed = %.15g   exact = %.15g   err = %.3e\n', Iq, 1/3, errB);
fprintf('  TEST B : %s\n\n', ternary(passB, 'PASS', 'FAIL'));

%% ---------- summary ----------
fprintf('============================================================\n');
fprintf(' SUMMARY :  A = %s   |   B = %s\n', ...
        ternary(passA,'PASS','FAIL'), ternary(passB,'PASS','FAIL'));
fprintf('============================================================\n');

if ~(passA && passB)
    error('verify_patch:FAIL', 'One or more unit checks failed.');
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
