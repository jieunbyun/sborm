function [Gamma, GGamma, Lambda, GLambda] = bbGammaLambda( x, nSample, info, bpf_target, alpha, beta )

%{
Input:
- x: (decision variable #) x 1 column vector (last element is gamma)
- nSample: the total number of samples -- can be different from the lengths
of alpha and beta when an active stategy is employed
- info: a structure that contains information required to evaluate
limit-state functions
- bpf_target: a scalar

Output: 
- Gamma: scalar
- GGamma: |X|+1 column vector
- Lambda: scalar
- GLambda: |X|+1 column vector
%}

x = x(:);
alpha_x = cellfun( @(a) a(:).' * x(1:info.nDv), alpha );
nActiveSample = size(alpha,1);

cutsets = info.cutset;
nCutset = length(cutsets);
p_val = zeros(nActiveSample, nCutset);
p_ind = zeros(nActiveSample, nCutset);
for iCutInd = 1:nCutset
    iCutset = cutsets{iCutInd};

    [iPval, iPval_ind] = min(beta(:,iCutset) + alpha_x(:,iCutset), [], 2 );
    iPval = -( iPval - x(end) );

    p_val(:,iCutInd) = iPval;
    p_ind(:,iCutInd) = iCutset(iPval_ind);
end

p_grad = cell(nActiveSample, nCutset);
for iCutInd = 1:nCutset
    iAlpha_ind = sub2ind( size(alpha), (1:nActiveSample).', p_ind(:,iCutInd) );
    p_grad(:,iCutInd) = alpha(iAlpha_ind);
end
p_grad = cellfun( @(x) [-x(:);1], p_grad, 'UniformOutput', false );

phi = sum(p_val, 2);
[p_min, p_min_ind] = min( p_val, [], 2 );
psi = phi - p_min;

phi_grad = cell(nActiveSample,1);
psi_grad = cell(nActiveSample,1);
for iSampleInd = 1:nActiveSample
    iP_grad = cell2mat(p_grad(iSampleInd,:));
    phi_grad{iSampleInd} = sum(iP_grad, 2);

    psi_grad{iSampleInd} = sum(iP_grad, 2) - iP_grad(:, p_min_ind(iSampleInd));
end

psi_loc = 1; phi_loc = 2;
[psi_phi_max, psi_phi_max_ind] = max( [psi phi], [], 2 );
Gamma = x(end) + 1/bpf_target/nSample * sum(psi_phi_max);
Lambda = 1/bpf_target/nSample * sum(phi);

GGamma = zeros( size(x) );
GLambda = zeros( size(x) );
for iSampleInd = 1:nActiveSample
    switch psi_phi_max_ind(iSampleInd)
        case psi_loc
            GGamma = GGamma + 1/bpf_target/nSample * psi_grad{iSampleInd};
        case phi_loc
            GGamma = GGamma + 1/bpf_target/nSample * phi_grad{iSampleInd};
    end

    GLambda = GLambda + 1/bpf_target/nSample * phi_grad{iSampleInd};
end
GGamma(end) = GGamma(end)+1;