function [A, b, Aeq, beq, lb, ub, info, samples_u] = dataFeasSet_cantilever(nSample, bpf_target)

info.bpf_target = bpf_target;

A = [];
b = [];
Aeq = [];
beq = [];
lb = [50; 500];
ub = [150; 1500];

info.L = 5;
info.X_mean = 150;
info.X_std = 30;
info.T_std = 20;
info.M_std = 300;

info.nRv = 3; % number of random variables
info.nDv = 2; % number of decision variables
info.nG = 5; % number of limit-state functions

info.T_ind = 1;
info.M_ind = 2;
info.X_ind = 3;

info.cutset = {[1 2], [3 4], [3 5]};
info.evalG_name = 'example.evalG_cantilever';
info.evalGgrad_name = 'example.evalGgrad_cantilever';
info.evalSample_name = 'example.evalSample_cantilever';
info.evalCost_name = 'example.evalCost_cantilever';
info.filename = 'cantilever';

samples_u = randn(nSample, info.nRv);