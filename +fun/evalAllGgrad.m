function G_grads = evalAllGgrad( evalGgrad_name, x, samples, info )
%{
Input:
- evalG_name: a string 
- x: decision variable values to compute gradient upon
- samples: (sample #) x (random variable #) array
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- G_grads: (sample #) x (Limit-state function #) cell array (in each cell,
(decision variable #)-dim vector)
%}

nSample = size(samples,1);
G_grads = cell(nSample, info.nG);

for iGInd = 1:info.nG
    iG_grad = feval( evalGgrad_name, iGInd, x, samples, info);

    G_grads(:,iGInd) = mat2cell( iG_grad, ones(nSample,1) );
end