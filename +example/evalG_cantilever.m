function g_val = evalG_cantilever( G_ind, x, samples, info )
%{
Input:
- G_ind: a scalar
- x: decision variable values to compute upon
- samples: (sample #) x (random variable #) array
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- g_val: (sample #) x 1 array
- g_grad: (sample #) x (decision variable #) array
%}

L = info.L;
T = samples(:, info.T_ind);
M = samples(:, info.M_ind);
X = samples(:, info.X_ind);

switch G_ind
    case 1
        g_val = -( T - 5*X/16 );
    case 2
        g_val = -( M - L*X );
    case 3
        g_val = -( M - 3*L*X/8 );
    case 4
        g_val = -( M - L*X/3 );
    case 5
        g_val = -( M + 2*L*T - L*X );
end
