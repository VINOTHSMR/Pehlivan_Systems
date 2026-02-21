clc; clear; close all;

%% Fractional order
alpha = 0.9;

%% Time interval
t0 = 0;
tfinal = 20;

%% Step sizes
hlist = [1e-2, 5e-3, 1e-3];

%% Initial conditions
y0 = [0.1; 0.2; 0.3];

%% Parameters (two only)
param.a = 0.5;
param.b = 0.8;

%% Preallocate CPU times
cpu_fde12 = zeros(length(hlist),1);
cpu_abm   = zeros(length(hlist),1);

for k = 1:length(hlist)
    h = hlist(k);

    % ---- FDE12 ----
    tic
    fde12(alpha,@pehlivan_cf,t0,tfinal,y0,h,param,1);
    cpu_fde12(k) = toc;

    % ---- Classical ABM ----
    tic
    pehlivan_abm_direct(alpha,@pehlivan_cf,t0,tfinal,y0,h,param);
    cpu_abm(k) = toc;
end

%% Display table
CPU_Table = table(hlist',cpu_fde12,cpu_abm, ...
    'VariableNames',{'StepSize_h','CPU_FDE12_sec','CPU_ABM_sec'});

disp(CPU_Table);
