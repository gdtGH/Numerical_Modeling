% DataTest.m
% Defines the three Riemann test cases for HW3 — Nonlinear SWE with FVM.
% Returns a cell array of Data structs, one per test case.
%
% Usage:
%   Tests = DataTest();
%   Data  = Tests{1};   % Dam Break
%   Data  = Tests{2};   % Reverse Dam Break
%   Data  = Tests{3};   % Stationary Shock

function Tests = DataTest()

%==========================================================================
% SHARED PARAMETERS
%==========================================================================
L   = 1;       % domain length
g   = 9.81;    % gravitational acceleration [m/s^2]
M   = 200;     % number of FVM cells
CFL = 0.9;     % CFL number
T   = 1.0;     % final time [s]

%==========================================================================
% TEST a — Dam Break
%   h(x,0) = 1.0  if x < L/2
%   h(x,0) = 0.1  if x >= L/2
%   q(x,0) = 0    everywhere
%   Expected: rarefaction leftward + shock rightward
%==========================================================================
Da.name   = 'HW3_DamBreak';
Da.domain = [0, L];
Da.T      = T;
Da.M      = M;
Da.CFL    = CFL;
Da.g      = g;
Da.h0     = @(x)  1.0 .* (x <  L/2) + 0.1 .* (x >= L/2);
Da.q0     = @(x)  0 .* x;

%==========================================================================
% TEST b — Reverse Dam Break
%   h(x,0) = 0.1  if x < L/2
%   h(x,0) = 1.0  if x >= L/2
%   q(x,0) = 0    everywhere
%   Expected: specular mirror of Test a
%==========================================================================
Db.name   = 'HW3_RevDamBreak';
Db.domain = [0, L];
Db.T      = T;
Db.M      = M;
Db.CFL    = CFL;
Db.g      = g;
Db.h0     = @(x)  0.1 .* (x <  L/2) + 1.0 .* (x >= L/2);
Db.q0     = @(x)  0 .* x;

%==========================================================================
% TEST c — Stationary Shock (hydraulic jump)
%   h(x,0) = 1.0  if x < L/2
%   h(x,0) = 0.5  if x >= L/2
%   q(x,0) = 0    everywhere
%   Expected: standing shock wave at the interface
%==========================================================================
Dc.name   = 'HW3_Shock';
Dc.domain = [0, L];
Dc.T      = T;
Dc.M      = M;
Dc.CFL    = CFL;
Dc.g      = g;
Dc.h0     = @(x)  1.0 .* (x <  L/2) + 0.5 .* (x >= L/2);
Dc.q0     = @(x)  0 .* x;

%==========================================================================
% ASSEMBLE OUTPUT
%==========================================================================
Tests = {Da, Db, Dc};

end
