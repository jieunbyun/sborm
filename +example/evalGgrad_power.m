function g_grad = evalGgrad_power( G_ind, x, samples, info )
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

samples1 = samples(:,G_ind);

x1_ind = info.comp2x(G_ind);
alpha1 = info.alpha(x1_ind);
beta1 = info.beta(x1_ind);
x1 = x( x1_ind );

faultIntensity = alpha1 * beta1 * exp( -beta1 * x1 );

nSample = size(samples,1);
g_grad  = zeros(nSample, info.nDv);

g_grad(:,x1_ind) = (log(samples1)) / faultIntensity * beta1;

