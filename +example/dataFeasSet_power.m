function [A, b, Aeq, beq, lb, ub, info, samples_u] = dataFeasSet_power(nSample, bpf_target)

info.bpf_target = bpf_target;

A = [];
b = [];
Aeq = [];
beq = [];
lb = [1 1 1 1 1 1].'; % days
ub = [10 10 10 10 10 10].'; % days

info.nRv = 12; % number of random variables
info.nDv = 6; % number of decision variables
info.nG = 12; % number of limit-state functions

info.type_ds = 1; info.type_cb = 2; info.type_pt = 3; info.type_db = 4; info.type_tb = 5; info.type_fb = 6;
info.comp2x = [info.type_ds info.type_ds info.type_ds info.type_cb info.type_cb info.type_pt ...
    info.type_pt info.type_db info.type_db info.type_tb info.type_fb info.type_fb];

info.TargetTime = 365; % days
info.alpha = 9 * ones(info.nDv, 1);
info.beta = 2 * ones(info.nDv, 1);

info.cutset = {[1 2], [4 5], [4 7], [4 9], [5 6], [6 7], [6 9], [5 8], [7 8], [8 9], ...
    [11 12], [1 3 5], [1 3 7], [1 3 9], [2 3 4], [2 3 6], [2 3 8], [4 10 12], [6 10 12], ...
    [8 10 12], [5 10 11], [7 10 11], [9 10 11], [1 3 10 12], [2 3 10 11]};
info.evalG_name = 'example.evalG_power';
info.evalGgrad_name = 'example.evalGgrad_power';
info.evalSample_name = 'example.evalSample_power';
info.evalCost_name = 'example.evalCost_power';
info.filename = 'power';

samples_u = rand(nSample, info.nRv);