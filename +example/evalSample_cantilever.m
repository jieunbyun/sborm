function samples = evalSample_cantilever( samples_u, x, info )
%{
Input:
- samples_u: (Sample #) x (random variable #) array -- standard normal
random variables
x: (decision variable #) x 1 column vector
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- samples: (Sample #) x (random variable #) array -- actual values
%}

samples = zeros( size(samples_u) );

samples(:, info.T_ind) = x(info.T_ind) + info.T_std * samples_u(:, info.T_ind);
samples(:, info.M_ind) = x(info.M_ind) + info.M_std * samples_u(:, info.M_ind);
samples(:, info.X_ind) = info.X_mean + info.X_std * samples_u(:, info.X_ind);