function [cost, cost_grad] = evalCost_truss( x, info )
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

cost = x(1:info.nDv, :) .* repmat( info.xCostWeight(:), 1, size(x,2) );
cost = sum(cost,1);
if nargout > 1
    cost_grad = repmat( info.xCostWeight(:), 1, size(x,2) );
end