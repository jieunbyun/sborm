function g_grad = evalGgrad_cantilever( G_ind, x, samples, info )
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

L = info.L;

nSample = size(samples, 1);

switch G_ind
    case 1
        g_grad = repmat([-1 0], nSample, 1);
    case 2
        g_grad = repmat([0 -1], nSample, 1);
    case 3
        g_grad = repmat([0 -1], nSample, 1);
    case 4
        g_grad = repmat([0 -1], nSample, 1);
    case 5
        g_grad = repmat([-2*L -1], nSample, 1);
end

