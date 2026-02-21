% Fractal-Fractional derivative of Pehlivan system
% Caputo-derivative sense
clc; clear; close all;

%% ================= NUMERICAL PARAMETERS =================
h = 0.01;
t(1) = 0.1;
tfinal = 300;
t = t(1):h:tfinal;
N = ceil((tfinal - t(1))/h);

x = zeros(1, N+1);
y = zeros(1, N+1);
z = zeros(1, N+1);

x(1) = 0.4;
y(1) = 0.5;
z(1) = 0.5;

% ===== Choose (alpha, tau) =====
alpha = 0.70;  tau = 0.70;      % D1
% alpha=0.8; tau=0.8;           % D2
% alpha=0.85;  tau=0.995;       % D3
% alpha=1; tau=0.95;            % D4
% alpha=1; tau=1;               % D5

%% ================= SYSTEM PARAMETERS =================
a = 0.5;
b = 0.5;

f1 = @(t,x,y,z) (y - x);
f2 = @(t,x,y,z) (-x.*z + a.*y);
f3 = @(t,x,y,z) (-b + x.*y);

%% ================= SOLVER =================
tic
for n = 1:N
    j = 2:n;

    x(n+1)=x(1)+(tau*h^alpha/gamma(alpha+2))*sum( ...
        ((n+1-j).^alpha.*(n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2*alpha)).* ...
        f1(t(j),x(j),y(j),z(j)).*t(j).^(tau-1) ...
      -((n+1-j).^(alpha+1)-(n-j).^alpha.*(n-j+1+alpha)).* ...
        f1(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));

    y(n+1)=y(1)+(tau*h^alpha/gamma(alpha+2))*sum( ...
        ((n+1-j).^alpha.*(n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2*alpha)).* ...
        f2(t(j),x(j),y(j),z(j)).*t(j).^(tau-1) ...
      -((n+1-j).^(alpha+1)-(n-j).^alpha.*(n-j+1+alpha)).* ...
        f2(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));

    z(n+1)=z(1)+(tau*h^alpha/gamma(alpha+2))*sum( ...
        ((n+1-j).^alpha.*(n-j+2+alpha)-(n-j).^alpha.*(n-j+2+2*alpha)).* ...
        f3(t(j),x(j),y(j),z(j)).*t(j).^(tau-1) ...
      -((n+1-j).^(alpha+1)-(n-j).^alpha.*(n-j+1+alpha)).* ...
        f3(t(j-1),x(j-1),y(j-1),z(j-1)).*t(j-1).^(tau-1));

    t(n+1) = t(n) + h;
end
cpuTime = toc;

%% ================= PLOT SETTINGS =================
markerStyle = 'p';
markerSize  = 1;

%% ================= SAVE SETTINGS =================
outputDir = 'Figures_Pehlivan_FF';
if ~exist(outputDir,'dir')
    mkdir(outputDir);
end

dpi = 600;
tag = sprintf('alpha_%0.2f_tau_%0.2f_h_%g', alpha, tau, h);

%% ================= PHASE PORTRAITS =================

figure('Color','w');
plot3(x,y,z,'r','Marker',markerStyle,'MarkerSize',markerSize);
xlabel('x'); ylabel('y'); zlabel('z');
box on;
print(gcf, fullfile(outputDir,['Pehlivan_3D_',tag,'.png']),'-dpng',['-r',num2str(dpi)]);

figure('Color','w');
plot(x,y,'b','Marker',markerStyle,'MarkerSize',markerSize);
xlabel('x'); ylabel('y');
box on;
print(gcf, fullfile(outputDir,['Pehlivan_xy_',tag,'.png']),'-dpng',['-r',num2str(dpi)]);

figure('Color','w');
plot(x,z,'g','Marker',markerStyle,'MarkerSize',markerSize);
xlabel('x'); ylabel('z');
box on;
print(gcf, fullfile(outputDir,['Pehlivan_xz_',tag,'.png']),'-dpng',['-r',num2str(dpi)]);

figure('Color','w');
plot(y,z,'m','Marker',markerStyle,'MarkerSize',markerSize);
xlabel('y'); ylabel('z');
box on;
print(gcf, fullfile(outputDir,['Pehlivan_yz_',tag,'.png']),'-dpng',['-r',num2str(dpi)]);

%% ================= TIME-SERIES PLOT =================

figure('Color','w');
plot(t, x, 'r', 'LineWidth', 1.5); hold on;
plot(t, y, 'b', 'LineWidth', 1.5);
plot(t, z, 'g', 'LineWidth', 1.5);
xlabel('Time');
ylabel('Value');
legend('x(t)','y(t)','z(t)','Location','best');
grid on;
box on;

print(gcf, fullfile(outputDir,['Pehlivan_TimeSeries_',tag,'.png']), ...
      '-dpng',['-r',num2str(dpi)]);

fprintf('Simulation completed in %.2f seconds\n', cpuTime);
