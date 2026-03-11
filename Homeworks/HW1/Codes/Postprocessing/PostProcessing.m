function [solutions] = PostProcessing(Data, femregion, uh)
    x    = femregion.dof(:,1);
    u_ex = Data.uex(x, Data.omega, Data.ro, Data.vel);
    solutions = struct('u_ex', u_ex, 'uh', uh);
    fprintf('============================================================\n')
end