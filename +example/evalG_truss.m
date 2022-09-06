function g_val = evalG_truss( G_ind, x, samples, info )
%{
Input:
- G_ind: a scalar
- x: decision variable values to compute upon
- samples: (sample #) x (random variable #) array
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- g_val: (sample #) x 1 array
%}

L = samples(:, info.L_ind);

memberInd = info.G_memberIndex(G_ind);
Y = samples(:, info.Y_ind(memberInd));
xInd = info.member2x(memberInd);

memberForce = info.G_memberForce(G_ind);
g_val = L * memberForce - Y * x(xInd);