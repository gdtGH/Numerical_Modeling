clear; close all; clc;

addpath Assembly
addpath BoundaryConditions
addpath Errors
addpath MeshGeneration
addpath FESpace
addpath Postprocessing

%% =========================================================
%GLOBAL STYLE — modifica solo qui
%=========================================================
fontsize= 25;% fontsize assi e legenda
fontsizeT = 50;% fontsize titoli
lw= 2; % linewidth
c1= '#0072BD';
c2= '#D95319';
c3= '#77AC30';

%% =========================================================
%POINT 2 — Verification and convergence
%=========================================================

Data = DataTest('HW1_P2');
Data.calc_errors= 1;
Data.eig_analysis = 0;
Data.plot_eigvct= 0;
Data.visual_graph = 0;

nElVec = [10, 20, 40, 80];
err = cell(4,1); sol = cell(4,1); fem = cell(4,1);

for i = 1:4
    [err{i}, sol{i}, fem{i}, ~] = MainHMZ(Data, nElVec(i));
end

hVec = cellfun(@(f) f.h,fem);
eVecL2 = cellfun(@(e) e.L2, err);

%% Plot — FEM vs Exact, coarse mesh (N=10)

x_h= fem{1}.dof(:,1);
u_h= real(sol{1}.uh);
u_ex = real(sol{1}.u_ex);

figure('Name','FEM vs Exact — coarse','NumberTitle','off');
plot(x_h, u_ex, '--', 'Color',c2, 'LineWidth',lw, 'DisplayName','Re(u_{ex})'); hold on;
plot(x_h, u_h,'-','Color',c1, 'LineWidth',lw, 'DisplayName','Re(u_h)');
xlabel('x [m]','FontSize',fontsize);
ylabel('Pressure [Pa]','FontSize',fontsize);
title(sprintf('FEM vs Exact Solution — N = %d (coarse mesh)', nElVec(1)),...
'FontSize',fontsizeT,'FontWeight','bold');
legend('Location','best','FontSize',fontsize);
grid on; box on; ax=gca; ax.FontSize=fontsize;
xlim([0 1]);

%% Plot — FEM vs Exact, fine mesh (N=80)

x_h= fem{4}.dof(:,1);
u_h= real(sol{4}.uh);
u_ex = real(sol{4}.u_ex);

figure('Name','FEM vs Exact — fine','NumberTitle','off');
plot(x_h, u_ex, '--', 'Color',c2, 'LineWidth',lw, 'DisplayName','Re(u_{ex})'); hold on;
plot(x_h, u_h,'-','Color',c1, 'LineWidth',lw, 'DisplayName','Re(u_h)');
xlabel('x [m]','FontSize',fontsize);
ylabel('Pressure [Pa]','FontSize',fontsize);
title(sprintf('FEM vs Exact Solution — N = %d (fine mesh)', nElVec(4)),...
'FontSize',fontsizeT,'FontWeight','bold');
legend('Location','best','FontSize',fontsize);
grid on; box on; ax=gca; ax.FontSize=fontsize;
xlim([0 1]);

%% Plot — Convergence L2

figure('Name','Convergence L2','NumberTitle','off');

ref = hVec.^2 * (eVecL2(end)/hVec(end)^2);
loglog(hVec, ref, '--k', 'LineWidth', 1.5, 'DisplayName','ReferenceO(h^2)');
hold on;
loglog(hVec, eVecL2, '-o', 'Color',c2, 'LineWidth',lw, ...
    'MarkerSize',8, 'MarkerFaceColor',c2, 'DisplayName','||u - u_h||_{L^2}');

rates = log(eVecL2(1:end-1)./eVecL2(2:end)) ./ log(hVec(1:end-1)./hVec(2:end));
for i = 1:length(rates)
    xm = sqrt(hVec(i)*hVec(i+1));
    ym = sqrt(eVecL2(i)*eVecL2(i+1)) * 1.8;
    text(xm, ym, sprintf('p = %.2f', rates(i)), ...
        'FontSize',fontsize,'HorizontalAlignment','center',...
        'Color',c2,'FontWeight','bold');
end

xlabel('h(mesh size)','FontSize',fontsize);
ylabel('||u - u_h||_{L^2(0,1)}','FontSize',fontsize);
title('Convergence Analysis — L^2 Norm','FontSize',fontsizeT,'FontWeight','bold');
legend('Location','best','FontSize',fontsize);
grid on; box on; ax=gca; ax.FontSize=fontsize;

fprintf('\n--- L2 Convergence Rates ---\n');
for i = 1:length(rates)
    fprintf('%.4f -> %.4f|p = %.4f\n', hVec(i), hVec(i+1), rates(i));
end

%% =========================================================
%POINT 3 — Wave reflection at impedance jump
%=========================================================

Data3 = DataTest('HW1_P3');
Data3.visual_graph = 0;
Data3.calc_errors= 0;
Data3.eig_analysis = 0;
Data3.plot_eigvct= 0;

[~, sol3, fem3, ~] = MainHMZ(Data3, 200);
x3 = fem3.coord(:,1);
u3 = sol3.uh;

figure('Name','Amplitude — Point 3','NumberTitle','off');
plot(x3, abs(u3), '-', 'Color',c1, 'LineWidth',lw);
xline(0.5,'--k','LineWidth',1.2,'Label','x = 0.5(\mu jump)',...
'LabelVerticalAlignment','bottom','FontSize',fontsize,'HandleVisibility','off');
xlabel('x [m]','FontSize',fontsize); ylabel('|u(x)|','FontSize',fontsize);
title('Amplitude of the FEM Solution — Impedance Jump at x = 0.5',...
'FontSize',fontsizeT,'FontWeight','bold');
grid on; box on; ax=gca; ax.FontSize=fontsize; xlim([0 1]);

figure('Name','Phase — Point 3','NumberTitle','off');
plot(x3, unwrap(angle(u3)*180/pi), '-', 'Color',c2, 'LineWidth',lw);
xline(0.5,'--k','LineWidth',1.2,'Label','x = 0.5(\mu jump)',...
'LabelVerticalAlignment','top','FontSize',fontsize,'HandleVisibility','off');
xlabel('x [m]','FontSize',fontsize); ylabel('Phase [deg]','FontSize',fontsize);
title('Phase of the FEM Solution — Impedance Jump at x = 0.5',...
'FontSize',fontsizeT,'FontWeight','bold');
grid on; box on; ax=gca; ax.FontSize=fontsize; xlim([0 1]);

%% =========================================================
%POINT 4 — Numerical Dispersion Analysis
%=========================================================

Data4 = DataTest('HW1_P4');
Data4.visual_graph = 0;
Data4.calc_errors= 0;
Data4.eig_analysis = 0;
Data4.plot_eigvct= 0;

nElVec4 = [10, 20, 40, 80, 160];
hOverLambda_all = [];
khOverK_all = [];

for i = 1:length(nElVec4)
    N = nElVec4(i);
    Region= CreateMesh(Data4, N);
    Femregion = CreateFemregion(Data4, Region);
    h = Femregion.h;
    [A_nbc, M_nbc] = Matrix1D(Data4, Femregion);
    ndof= Femregion.ndof;
    idx = 1:ndof-1;
    A_red = A_nbc(idx,idx);
    M_red = M_nbc(idx,idx);
    [~, D] = eig(full(A_red), full(M_red));
    k2_h = diag(D); k2_h = real(k2_h);
    k2_h = sort(k2_h(k2_h > 0));
    k_h= sqrt(k2_h);
    n_vec= (1:length(k_h))';
    k_ex = (2*n_vec - 1)*pi/2;
    lambda = 2*pi ./ k_ex;
    hOverLambda = h ./ lambda;
    khOverK = k_h ./ k_ex;
    mask = hOverLambda < 0.5;
    hOverLambda_all = [hOverLambda_all; hOverLambda(mask)];
    khOverK_all = [khOverK_all; khOverK(mask)];
end

figure('Name','Dispersion Curve — Point 4','NumberTitle','off');
hLam_ref = linspace(0, 0.5, 200);
plot(hLam_ref, ones(size(hLam_ref)), 'k--', 'LineWidth',1.5,...
 'DisplayName','Exactk_h/k = 1');
hold on;
scatter(hOverLambda_all, khOverK_all, 18, 'b', 'filled',...
'DisplayName','FEM P1 — k_h/k');
xlabel('h/\lambda(mesh size / wavelength)','FontSize',fontsize);
ylabel('k_h / k','FontSize',fontsize);
title('Numerical Dispersion Curve — FEM P1','FontSize',fontsizeT,'FontWeight','bold');
legend('Location','best','FontSize',fontsize);
xlim([0 0.5]); ylim([0.8 1.05]);
grid on; box on; ax=gca; ax.FontSize=fontsize;

%% =========================================================
%POINT 5 — Perfectly Matched Layer (PML)
%=========================================================

Data5 = DataTest('HW1_P5');
Data5.visual_graph = 0;
Data5.calc_errors= 0;
Data5.eig_analysis = 0;
Data5.plot_eigvct= 0;

[~, sol5, fem5, ~] = MainHMZ(Data5, 240);
x5 = fem5.coord(:,1);
u5 = sol5.uh;
mask_phys = x5 <= 1 + 1e-10;
x5_phys= x5(mask_phys);
u5_phys= u5(mask_phys);

figure('Name','Amplitude — P3 vs P5','NumberTitle','off');
plot(x3,abs(u3),'-','Color',c1, 'LineWidth',lw,...
 'DisplayName','Impedance BC (Point 3)');
hold on;
plot(x5_phys, abs(u5_phys), '--', 'Color',c2, 'LineWidth',lw,...
 'DisplayName','PML (Point 5)');
xline(0.5,'--k','LineWidth',1.2,'Label','x = 0.5(\mu jump)',...
'LabelVerticalAlignment','bottom','FontSize',fontsize,'HandleVisibility','off');
xlabel('x [m]','FontSize',fontsize); ylabel('|u(x)|','FontSize',fontsize);
title('Amplitude Comparison — Impedance BC vs PML',...
'FontSize',fontsizeT,'FontWeight','bold');
legend('Location','best','FontSize',fontsize);
grid on; box on; ax=gca; ax.FontSize=fontsize; xlim([0 1]);

figure('Name','Phase — P3 vs P5','NumberTitle','off');
plot(x3,unwrap(angle(u3)*180/pi), '-','Color',c1, 'LineWidth',lw,...
 'DisplayName','Impedance BC (Point 3)');
hold on;
plot(x5_phys, unwrap(angle(u5_phys)*180/pi),'--', 'Color',c2, 'LineWidth',lw,...
 'DisplayName','PML (Point 5)');
xline(0.5,'--k','LineWidth',1.2,'Label','x = 0.5(\mu jump)',...
'LabelVerticalAlignment','top','FontSize',fontsize,'HandleVisibility','off');
xlabel('x [m]','FontSize',fontsize); ylabel('Phase [deg]','FontSize',fontsize);
title('Phase Comparison — Impedance BC vs PML',...
'FontSize',fontsizeT,'FontWeight','bold');
legend('Location','best','FontSize',fontsize);
grid on; box on; ax=gca; ax.FontSize=fontsize; xlim([0 1]);