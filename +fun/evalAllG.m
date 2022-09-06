function Gvals = evalAllG( evalG_name, x, samples, info )
%{
Input:
- evalG_name: a string 
- x: decision variable values to compute gradient upon
- samples: (sample #) x (random variable #) array
- info: a structure that contains information required to evaluate
limit-state functions

Output:
- Gvals: (sample #) x (Limit-state function #) array
%}

nSample = size(samples,1);
Gvals = zeros(nSample, info.nG);

for iGInd = 1:info.nG
    Gvals(:,iGInd) = feval( evalG_name, iGInd, x, samples, info);
end