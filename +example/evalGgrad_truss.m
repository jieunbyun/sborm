function g_grad = evalGgrad_truss( G_ind, x, samples, info )
%{
Input:
- G_ind: a scalar
- x: decision variable values to compute gradient upon
- samples: (sample #) x (random variable #) array
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- g_grad: (sample #) x (decision variable #) array
%}

memberInd = info.G_memberIndex(G_ind);
Y = samples(:, info.Y_ind(memberInd));
xInd = info.member2x(memberInd);

nSample = size(samples,1);
g_grad = zeros( nSample, info.nDv );
g_grad(:, xInd) = - Y;