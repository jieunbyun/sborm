function [cost, cost_grad] = evalCost_cantilever( x, info )
%{
Input:
- x: (decision variable #) x (solution #) array (whether
it includes gamma does not matter)
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- cost: 1 x (solution #) vector
- cost_grad: (decision variable #) x (solution #) vector
%}

cost = 2*x(info.M_ind, :) + x(info.T_ind, :);
if nargout > 1
    cost_grad = repmat( [1; 2], 1, size(x,2) );
end