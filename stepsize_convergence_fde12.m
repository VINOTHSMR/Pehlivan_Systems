clc; clear; close all;

%% Fractional order
alpha = 0.9;

%% Time interval
t0 = 0;
tfinal = 20;

%% Initial conditions
y0 = [0.1; 0.2; 0.3];

%% System parameters (two only)
param.a = 0.5;
param.b = 0.8;

%% Step sizes
h_list  = [1e-1, 5e-2, 1e-2, 5e-3];
h_ref   = 1e-4;     % reference step size

%% Corrector
mu = 1;

%% Reference solution
[t_ref, y_ref] = fde12(alpha,@pehlivan_cf,t0,tfinal,y0,h_ref,param,mu);

yT_ref = y_ref(:,end);

%% Preallocate
Err = zeros(length(h_list),1);
CPU = zeros(length(h_list),1);

%% Loop over step sizes
for k = 1:length(h_list)
    h = h_list(k);

    tic
    [t,y] = fde12(alpha,@pehlivan_cf,t0,tfinal,y0,h,param,mu);
    CPU(k) = toc;

    yT = y(:,end);
    Err(k) = norm(yT - yT_ref, inf);
end

%% Experimental order of convergence (EOC)
EOC = zeros(length(h_list)-1,1);
for k = 1:length(EOC)
    EOC(k) = log(Err(k)/Err(k+1))/log(h_list(k)/h_list(k+1));
end

%% Display table
Convergence_Table = table( ...
    h_list', Err, [EOC; NaN], CPU, ...
    'VariableNames',{'StepSize_h','Error_inf','EOC','CPU_time_sec'} );

disp(Convergence_Table);
