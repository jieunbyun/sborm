function [cost, cost_grad] = evalCost_power( x, info )
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

cost = sum( x(1:info.nDv, :), 1 );
if nargout > 1
    cost_grad = ones( info.nDv, size(x,2) );
end