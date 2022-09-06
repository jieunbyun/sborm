clear;
filename = 'cantilever_exact';
rng(1)

bpf_target = 1e-3; 

cov_target = 0.05; 
nSample = ceil( (1-bpf_target)/bpf_target/cov_target^2 );

[pars.X.A, pars.X.b, pars.X.Aeq, pars.X.beq, pars.X.lb, pars.X.ub, info, samples_u] = example.dataFeasSet_cantilever(nSample, bpf_target);

% Problem
L = 5;

X_mean = 150;
X_std = 30;

M_mean_bound = [500 1500];
M_std = 300;

T_mean_bound = [50 150];
T_std = 20;

costFunc = @(M_mean,T_mean) 2*M_mean + T_mean;

bpf_target = 1e-3;
cov_target = 0.1;
nSample = ceil( (1-bpf_target)/bpf_target/cov_target^2 );

samples_X = X_mean + X_std * samples_u(:, info.X_ind);
samples_normal_M = samples_u(:, info.M_ind);
samples_normal_T = samples_u(:, info.T_ind);

% Grid for evaluation
nGrid = 1e2;
M_mean_grid = linspace( M_mean_bound(1), M_mean_bound(2), nGrid );
T_mean_grid = linspace( T_mean_bound(1), T_mean_bound(2), nGrid );
[M_mean_grid_2d,T_mean_grid_2d] = meshgrid(M_mean_grid,T_mean_grid);

% System evaluation
nGfun = 6;

bpf_mean = zeros( nGrid, nGrid );
for iRow = 1:nGrid
    for jCol = 1:nGrid
        ijM_mean = M_mean_grid_2d( iRow, jCol );
        ijT_mean = T_mean_grid_2d( iRow, jCol );

        ijSamples_M = ijM_mean + M_std * samples_normal_M;
        ijSamples_T = ijT_mean + T_std * samples_normal_T;

        ijG1 = -( ijSamples_T - 5*samples_X / 16 );
        ijG2 = -( ijSamples_M - L*samples_X );
        ijG3 = -( ijSamples_M - 3*L*samples_X/8 );
        ijG4 = -( ijSamples_M - L * samples_X / 3 );
        ijG5 = -( ijSamples_M + 2*L*ijSamples_T - L*samples_X );

        ijCutset1 = min( [ijG1 ijG2], [], 2 );
        ijCutset2 = min( [ijG3 ijG4], [], 2 );
        ijCutset3 = min( [ijG3 ijG5], [], 2 );

        ijSys = max( [ijCutset1 ijCutset2 ijCutset3], [], 2 );

        ijBpf = fun.evalBpfFromSamples(ijSys);

        bpf_mean(iRow, jCol) = ijBpf;

    end

    if ~rem(iRow, 10)
        disp(['Row ' num2str(iRow) ' (/' num2str(nGrid) ') done.'])
    end
end

% Figure
figure;
imagesc(M_mean_grid,T_mean_grid,bpf_mean)
colormap jet; axis xy;
colorbar;

% Comparison with Song and Der Kiureghian (2003)
M_mean = 1000;
T_mean = 110;

M_sample = M_mean + M_std * samples_normal_M;
T_sample = T_mean + T_std * samples_normal_T;

G1 = -( T_sample - 5*samples_X / 16 );
G2 = -( M_sample - L*samples_X );
G3 = -( M_sample - 3*L*samples_X/8 );
G4 = -( M_sample - L * samples_X / 3 );
G5 = -( M_sample + 2*L*T_sample - L*samples_X );

cutset1 = min( [G1 G2], [], 2 );
cutset2 = min( [G3 G4], [], 2 );
cutset3 = min( [G3 G5], [], 2 );

sys = max( [cutset1 cutset2 cutset3], [], 2 );
pf = mean( sys>0 );
disp(['Failure prob: ' num2str(pf)])

% Bpf
target_bpf = 0.001;
sample_grid2d_feasibleSol = (bpf_mean < target_bpf);
cost_grid2d = costFunc(M_mean_grid_2d, T_mean_grid_2d);
cost_grid2d(~sample_grid2d_feasibleSol) = inf;
[minCost, minCost_ind] = min(cost_grid2d(:));
[minCost_ind_row, minCost_ind_col] = find(cost_grid2d == minCost);

M_mean_optim = M_mean_grid_2d(sub2ind(size(M_mean_grid_2d), minCost_ind_row, minCost_ind_col));
T_mean_optim = T_mean_grid_2d(sub2ind(size(T_mean_grid_2d), minCost_ind_row, minCost_ind_col));
disp( ['Optimal solution (M mean, T mean): ' num2str(M_mean_optim(:)') ', ' num2str(T_mean_optim(:)')] )
disp( ['Bpf with optimal design: ' num2str(bpf_mean(sub2ind(size(M_mean_grid_2d), minCost_ind_row, minCost_ind_col)))] )

save(filename)