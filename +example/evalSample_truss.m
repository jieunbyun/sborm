function samples = evalSample_truss( samples_u, x, info )
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

samples(:, info.L_ind) = info.L_mean + info.L_std * samples_u( :, info.L_ind );
samples(:, info.Y_ind) = info.Y_mean + info.Y_std * samples_u(:, info.Y_ind);