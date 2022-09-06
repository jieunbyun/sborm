function [A, b, Aeq, beq, lb, ub, info, samples_u] = dataFeasSet_truss(nSample, bpf_target)

info.bpf_target = bpf_target;

nDv = 4; % number of decision variables

A = [];
b = [];
Aeq = [];
beq = [];
lb = 1 * ones( nDv,1 ); % 10-3 * m2
ub = 2 * ones( nDv,1 ); % 10-3 * m2

info.nDv = nDv; % number of decision variables
info.x2member = {[1 2 9 10] [3 8] [4 7] [5 6]};
info.nMember = max( cellfun(@max, info.x2member) );
info.member2x = zeros(info.nMember,1);
for iMemberInd = 1:info.nMember
    info.member2x(iMemberInd) = find( cellfun( @(x) ismember(iMemberInd, x), info.x2member ) );
end

info.nMember_x = cellfun( @length, info.x2member );
length_h = 2; length_v = 1.6; length_diag = norm([length_h, length_v]);
info.memberLengths = [length_diag length_h length_v length_h length_diag length_diag length_h length_v length_diag length_h]';
info.xCostWeight = zeros(nDv,1);
for iDvInd = 1:nDv
    info.xCostWeight( iDvInd ) = sum( info.memberLengths(info.x2member{iDvInd}) );
end
info.nRv = 1 + info.nMember; % number of random variables


% Get failure modes
load( 'truss_structuralAnalysis.mat', 'F0', 'F1', 'failMode1', 'sysFail1' )

info.cutset = {};
info.G_memberForce = [];
info.G_memberIndex = [];

iGInd = 0; iCutInd = 0;

% % 1 member fails
for iModeInd = 1:size(F1,1)
    iMemberForce = F1(iModeInd,:);
    if all(~iMemberForce)
        iGInd = iGInd + 1;
        info.G_memberForce = [info.G_memberForce; abs( F0(iModeInd) )];
        info.G_memberIndex = [info.G_memberIndex; failMode1(iModeInd)];

        iCutInd = iCutInd + 1;
        info.cutset = [info.cutset; {[iGInd]}];
    end
end

% % 2 members fail
for iMemberInd = 1:info.nMember
    if ~sysFail1( iMemberInd )
        for jMemberInd = setdiff( 1:info.nMember, iMemberInd )
            ijForcePair = [abs( F0(iMemberInd) ) abs( F1( iMemberInd, jMemberInd ) )];
            if all(ijForcePair > 0)
                iGInd = iGInd + 1;
                info.G_memberIndex = [info.G_memberIndex; iMemberInd];
                info.G_memberForce = [info.G_memberForce; ijForcePair(1)];

                iGInd = iGInd + 1;
                info.G_memberIndex = [info.G_memberIndex; jMemberInd];
                info.G_memberForce = [info.G_memberForce; ijForcePair(2)];

                iCutInd = iCutInd + 1;
                info.cutset = [info.cutset; {[iGInd-1 iGInd]}];
            end
        end
    end
end

info.nG = length(info.G_memberIndex);

% Random variables
info.L_ind = 1; % r.v. representing load
info.Y_ind = 1+(1:info.nMember); % r.v. representing yield strength
info.L_mean = 190;
% info.L_mean = 160;
info.L_std = 0.1 * info.L_mean;
info.Y_mean = 276; % MPa
info.Y_std = 0.05 * info.Y_mean;

% File names
info.evalG_name = 'example.evalG_truss';
info.evalGgrad_name = 'example.evalGgrad_truss';
info.evalSample_name = 'example.evalSample_truss';
info.evalCost_name = 'example.evalCost_truss';
info.filename = 'truss';

samples_u = randn(nSample, info.nRv);