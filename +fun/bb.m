function [f1,g1,f2,g2,ind] = bb(nSample, info, bpf_target, pars, alpha, beta)
[Gamma, GGamma, Lambda, GLambda] = fun.bbGammaLambda( pars.x, nSample, info, bpf_target, alpha, beta );
[f1, g1]  = feval( info.evalCost_name, pars.x, info );
g1  = [g1;0];
lambda = pars.lambda;
xnu    = pars.xnu;
x      = pars.x;
if Gamma>=Lambda
    f1 = f1 + pars.theta*Gamma     + 0.5*lambda*norm(x-xnu)^2;
    g1 = g1 + pars.theta*GGamma(:) + lambda*(x-xnu);
else
    f1 = f1 + pars.theta*Lambda    + 0.5*lambda*norm(x-xnu)^2;
    g1 = g1 + pars.theta*GLambda(:)+ lambda*(x-xnu);
end
%
f2 = pars.theta*Lambda;
g2 = pars.theta*GLambda(:);
%
ind= 1;
return