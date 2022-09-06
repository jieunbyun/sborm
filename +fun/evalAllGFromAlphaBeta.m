function Gvals = evalAllGFromAlphaBeta( alpha, beta, x )
%{
Input:
- alpha: a cell array of size (sample #) x (limit-state function #), each cell having size (decision variable #) x 1
- beta: a double array of size (sample #) x (limit-state function #)
- x: decision variable values to compute the values upon (does not matter whether it includes gamma)

Output:
- Gvals: (sample #) x (Limit-state function #) array
%}

nG = size(alpha,2);
nX = length(alpha{1});
nSample = size(alpha,1);
Gvals = zeros(nSample, nG);

for iGInd = 1:nG
    Gvals(:,iGInd) = cellfun( @(al) al(:).' * x(1:nX), alpha(:,iGInd) );
end
Gvals = Gvals + beta;