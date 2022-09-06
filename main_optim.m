clear;
scriptVersionName = '_optim';
rng(1)

bpf_target = 1e-3; 

cov_target = 0.05; 
nSample = ceil( (1-bpf_target)/bpf_target/cov_target^2 );

% [pars.X.A, pars.X.b, pars.X.Aeq, pars.X.beq, pars.X.lb, pars.X.ub, info, samples_u] = example.dataFeasSet_cantilever(nSample, bpf_target);
% [pars.X.A, pars.X.b, pars.X.Aeq, pars.X.beq, pars.X.lb, pars.X.ub, info, samples_u] = example.dataFeasSet_truss(nSample, bpf_target);
[pars.X.A, pars.X.b, pars.X.Aeq, pars.X.beq, pars.X.lb, pars.X.ub, info, samples_u] = example.dataFeasSet_power(nSample, bpf_target);

lb_rat = 0.5;
x_init = lb_rat * pars.X.lb + (1-lb_rat) * pars.X.ub;
samples = feval(info.evalSample_name, samples_u, x_init, info );
Gvals = fun.evalAllG( info.evalG_name, x_init, samples, info );
[Gsys, G_ind] = fun.evalGFromCutsets( Gvals, info.cutset );
[bpf, gamma] = fun.evalBpfFromSamples( Gsys );
if bpf == 1
    gamma = min(Gsys)*1.1;
elseif bpf == 0
    gamma = max(Gsys)*1.1;
end

pars.x = [x_init;gamma]; 
pars.n = info.nDv+1;
pars.X.lb   = [pars.X.lb(:);-inf];
pars.X.ub   = [pars.X.ub(:);inf];
pars.QPdual = false;    
pars.matlab = false;
pars.nBun   = min(5*pars.n,300);%10*pars.n+100;
pars.nBun2  = 1;
pars.MaxIt  = 3000;
pars.tol    = 1e-2;
pars.kappa  = 0.01;
pars.tmin   = 1e-5;
pars.t      = 1;
pars.tmax   = 1e2;
pars.flow    =-inf;  

% for outer loop
pars.theta0 = 1; 
pars.theta_max = 1e5;
pars.lambda0 = 0.01; 
pars.tol_out = 1e-2;
pars.kappa_out = 0.01;
pars.nActiveRatio = 2;
nActiveSample = ceil( nSample * bpf_target * pars.nActiveRatio );

sol_diff = 1;
pars.maxOutloop = 1000;
pars.nOutloop = 1;
pars.nGradEvalRound = 0;
serious_step = 1;

cost = feval( info.evalCost_name, pars.x, info );

pars.theta = pars.theta0;
pars.lambda = pars.lambda0;
pars.xnu = pars.x;
Gvals_nu = Gvals;
Gsys_nu = Gsys;
cost_nu = cost;
samples_nu = samples;
tic;
while 1
    if pars.nOutloop > pars.maxOutloop
        break
    end
    disp(['[Outer loop ' num2str(pars.nOutloop) '] x: ' num2str( pars.xnu(1:end-1).' ) '; lambda: ' num2str(pars.lambda) '; theta: ' num2str(pars.theta) ' ..'])

    if serious_step
        pars.nGradEvalRound = pars.nGradEvalRound+1;

        [~, Gsys_sortInd] = sort( Gsys, 'descend' );
        activeSampleInd = Gsys_sortInd(1:nActiveSample);
        samples_nu_active = samples_nu( activeSampleInd,: );

        alpha = fun.evalAllGgrad( info.evalGgrad_name, pars.xnu, samples_nu_active, info );
        beta = zeros(nActiveSample, info.nG);
        for iGInd = 1:info.nG
            beta(:,iGInd) = Gvals_nu(activeSampleInd,iGInd) - cellfun( @(dg) dg(:).'* pars.xnu(1:end-1), alpha(:,iGInd) );
        end
    end

    out = fun.PBMDC(pars, nSample, info, bpf_target, alpha, beta);

    Fnu = cost_nu + pars.theta * max( [0, pars.xnu(end) + 1/bpf_target/nSample*sum(max( [zeros(nSample,1), Gsys_nu - pars.xnu(end)], [], 2 ))] );

    cost = feval( info.evalCost_name, out.sol, info );
    Gvals_approx = fun.evalAllGFromAlphaBeta( alpha, beta, out.sol );
    Gsys_approx = fun.evalGFromCutsets( Gvals_approx, info.cutset );
    F_approx = cost + pars.theta * max( [0, out.sol(end) + 1/bpf_target/nSample*sum(max( [zeros(nActiveSample,1), Gsys_approx - out.sol(end)], [], 2 ))] );

    del = Fnu - F_approx - 0.5 * pars.lambda * norm(out.sol - pars.xnu)^2;

    samples = feval(info.evalSample_name, samples_u, out.sol(1:end-1), info );
    Gvals = fun.evalAllG( info.evalG_name, out.sol(1:end-1), samples, info );
    Gsys = fun.evalGFromCutsets( Gvals, info.cutset );
    F = cost + pars.theta * max( [0, out.sol(end) + 1/bpf_target/nSample*sum(max( [zeros(nSample,1), Gsys - out.sol(end)], [], 2 ))] );

    sol_diff = norm(out.sol(1:end-1) - pars.xnu(1:end-1))^2 + ( out.sol(end) - pars.xnu(end) )^2; % Q to Welington: right?
    if sol_diff < pars.tol_out
        bpf = fun.evalBpfFromSamples( Gsys );

        if bpf < bpf_target 
            pars.x = out.sol;
            break;
        end
    end

    if F < Fnu - pars.kappa_out * del % serious step

        serious_step = 1;

        pars.xnu = out.sol;
        pars.x = pars.xnu;    
        Gvals_nu = Gvals;
        Gsys_nu = Gsys;
        cost_nu = cost;
        samples_nu = samples;

        pars.lambda = pars.lambda0;

    else % null step

        serious_step = 0;
        pars.lambda = 2*pars.lambda; 
        
    end
    
    pars.theta = min([1.5*pars.theta, pars.theta_max]);
    pars.nOutloop = pars.nOutloop + 1;
    disp( ['Seious step: ' num2str(serious_step), '; del: ' num2str(del)] )

end
computeSec = toc;

Gvals = fun.evalAllG( info.evalG_name, pars.x, samples, info );
[Gsys, G_ind] = fun.evalGFromCutsets( Gvals, info.cutset );
[bpf, gamma] = fun.evalBpfFromSamples( Gsys );
cost = feval( info.evalCost_name, pars.x, info );
pf = mean(Gsys>0);
disp(['Final bpf: ' num2str(bpf) '; target bpf: ' num2str(bpf_target) '; cost: ' num2str(cost)])

% Save output
params_save = table( bpf_target, cov_target, nSample, {x_init}, pars.theta0, pars.lambda0, pars.nActiveRatio, nActiveSample );
params_save.Properties.VariableNames = {'Target_bpf', 'Target_cov', 'No_sample', 'Init_sol', 'theta0', 'lambda0', 'Rat_active', 'No_active'};
result_save = table(computeSec, pars.nOutloop, pars.nGradEvalRound, {pars.x(:)'}, cost, bpf, pf);
result_save.Properties.VariableNames = {'Time_sec', 'No_outLoop', 'No_gradEvalRound', 'Optim_sol', 'Optim_cost', 'Optim_bpf', 'Optim_pf'};
save(strcat(info.filename, scriptVersionName), 'params_save', 'result_save')