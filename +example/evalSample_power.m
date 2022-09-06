function samples = evalSample_power( samples_u, x, info )
%{
Input:
- samples_u: (Sample #) x (random variable #) array -- uniform random variables
x: (decision variable #) x 1 column vector
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- samples: (Sample #) x (random variable #) array -- actual values
%}

samples = samples_u;