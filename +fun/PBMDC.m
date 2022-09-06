function out = PBMDC(pars, nSample, info, bpf_target, alpha, beta)
%--------------------------------------------------------------------------
% Proximal bundle method for DC programming using a DC model
%
% by Welington de Oliveira, May of 2017
%--------------------------------------------------------------------------
%   
%     min  f1(x) - f_2(x)    s.t. x \in X
%
%  where fi: R^n \to R  are convex but nonsmooth functions and X is a polyhedron 
%
%==========================================================================
tic;                                                                        % Marca o tempo de processamento
ns = 0;ss=0;out.ind   = -1;
tol         = pars.tol;%sqrt(pars.n)*pars.tol                                     
out.msg     = 'Ocorreu um erro no oraculo';                                 % Mensagem de erro
[f1,g1,f2,g2,ind] = fun.bb(nSample, info, bpf_target, pars, alpha, beta); % Primeira chamada ao oraculo
if ind < 0;  fprintf(1,'%s',out.msg); return;end                            % Confere se tudo foi bem no oraculo
pars.xc     = pars.x;                                                       % candidato a solucao
fxc         = f1 - f2;
%--------------------------------------------------------------------------
% Planos cortantes: f1(x) + g1'(y-x) \geq r    ===>  g1'y - r \geq g'x - f1(x)
modeloPC(1).G   = g1';
modeloPC(1).rhs = g1'*pars.x - f1;
% Planos cortantes: f1(x) + g1'(y-x) \geq r    ===>  g1'y - r \geq g'x - f1(x)
modeloPC(2).G   = g2';
modeloPC(2).rhs = g2'*pars.x - f2;
%--------------------------------------------------------------------------
% Ciclo
k = 0;nbun=1;nbun2=1;in=0;ds=0;
while (k <= pars.MaxIt)
    k = k + 1;
    if pars.QPdual
       [pars,mu,fcheck] = SubProbProximalDualDC(pars,modeloPC);
    else
       [pars,mu,fcheck] = SubProbProximalDC(pars,modeloPC);
    end
    %----------------------------------------------------------------------       
    % Teste de parada    
    normdif = norm(pars.x(1:end-1)-pars.xc(1:end-1));% / (1+norm(pars.x(1:end-1)));
    normdif_ratio = normdif / norm(pars.xc(1:end-1));
    if normdif <=tol
         fprintf(1,'k: %3d, |xkk -xck|/|xck|: %4.1d  fup: %4.1d \n',k,normdif_ratio,fxc);                  
        break
    end
    pars.mode = 1;
    [f1,g1,f2,g2,ind] = fun.bb(nSample, info, bpf_target, pars, alpha, beta);                                      
    if ind < 0; fprintf(1,'%s',out.msg); return;end 
    aux = fxc - 0.5*(pars.kappa/pars.tmax)*normdif^2;    
    ps=0;
    if f1-f2 <= aux
        ds = ds+1;
        pars.xc = pars.x;
        fxc    = f1 - f2;
        if ind < 0; fprintf(1,'%s',out.msg); return;end 
        ss = ss+1;
        if ss>3, pars.t = min(pars.tmax,pars.t*(1.25));ss=0;end
        ps=1;
        ns = 0;
        in = 0;
%
       if (nbun2<pars.nBun2)
            modeloPC(2).G   = [modeloPC(2).G; g2'];
            modeloPC(2).rhs = [modeloPC(2).rhs; g2'*pars.x - f2];       
            nbun2           = nbun2 + 1;
       else
            modeloPC(2).G   = [modeloPC(2).G(2:end,:); g2'];
            modeloPC(2).rhs = [modeloPC(2).rhs(2:end); g2'*pars.x - f2];       
            nbun2 = pars.nBun2;
       end
    else
        ns = ns+1;
        ss = 0;
         if ns>2
             pars.t=max(pars.tmin,pars.t/(in+2));
             in = in+1;
             ns=0;
         end
    end
    %----------------------------------------------------------------------
%    fprintf(pars.fid,'k: %3d, |xkk -xck|: %4.1d  fup: %4.1d,  \n',...% Imprime dados no arquivo
%               k,       normdif,   fxc );
    fprintf(1,'k: %3d, |xkk -xck|/|xck|: %4.1d  fup: %8.4f, t %4.2d ,  nb %3d,   seriou step %3d  \n',...% Imprime dados na tela
                k, normdif_ratio,   fxc, pars.t, nbun, ps );    
    disp( pars.x(:).' )
   %----------------------------------------------------------------------
   J  = (abs(mu) > 1e-7);
   nAc=sum(J);nNAc = nbun-nAc;
   if nbun>10*pars.nBun
       ps=1;
   end
   if (ps==0)||(nbun<pars.nBun)
      modeloPC(1).G   = [modeloPC(1).G; g1'];
      modeloPC(1).rhs = [modeloPC(1).rhs; g1'*pars.x - f1];       
   elseif (nAc>=pars.nBun)
       i = nAc-pars.nBun;
       aux = [1:nbun];J=aux(J);
       if i>3,J = J(i+3:end);else J=J(3:end);end
       d = modeloPC(1).G'*mu;
       modeloPC(1).G   = [modeloPC(1).G(J,:);d'; g1'];
       modeloPC(1).rhs = [modeloPC(1).rhs(J);d'*pars.x - fcheck; g1'*pars.x - f1];
   else  
       i = pars.nBun-nAc -1;
       aux = 1:nbun;I=aux(~J);I=I(nNAc-i+1:nNAc);
       modeloPC(1).G   = [modeloPC(1).G(I,:);modeloPC(1).G(J,:);g1'];
       modeloPC(1).rhs = [modeloPC(1).rhs(I);modeloPC(1).rhs(J); g1'*pars.x - f1];
   end
   nbun = length(modeloPC(1).rhs);
end
%==========================================================================
% Se voce chegou aqui, signifca que tudo ocorreu bem, ou o numero de
% iteracoes foi atingido:
if k > pars.MaxIt, out.msg = 'Numero maximo de iteracoes atingido';
else out.ind=1; out.msg = 'Teste de parada satisfeito'; end
%
out.k   = k;
out.sol = pars.xc;
out.val = fxc;
out.cpu = toc;
out.ds = ds;
    
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [pars,Dualvar,fcheck] = SubProbProximalDualDC(pars,modeloPC)
%
ncut1 = length(modeloPC(1).rhs);
ncut2 = length(modeloPC(2).rhs);
Aeq   = ones(1,ncut1);
beq   = 1;
lb    = zeros(ncut1,1);
model = inf;
Dualvar = ones(ncut1,1)/ncut1;
for i=1:ncut2    
   K     = (repmat(modeloPC(2).G(i,:),ncut1,1) - modeloPC(1).G  )';
   Q     = pars.t*(K'*K);
   condK = 1./condest(Q);
   tr    = 0; 
   sigma = 1.d-12;
   while condK < 1e-15 && tr<= 10
     tr    = tr + 1;
     Q     = Q + sigma*eye(ncut1);
     sigma = 10*sigma;
     condK = 1./condest(Q);
   end
   q   = modeloPC(1).rhs + K'*pars.xc;
   %-------------------------------------------------------------------------
   % resolve o QP para obter o novo iterado
   if pars.matlab
     [mu,~,ind]  = quadprog(Q,q,[],[],Aeq,beq,lb);
   else
     prob.obj   = q;
     prob.Q     = 0.5*sparse(Q);
     prob.A     = sparse(Aeq);
     prob.rhs   = beq;
     prob.sense = '=';
     prob.lb    = lb;
     par.OutputFlag =0; % 1 to print output
     res = gurobi(prob,par);
     if strcmp(res.status, 'OPTIMAL')|| strcmp(res.status, 'SUBOPTIMAL')
        mu   = res.x;
        ind = 1;
     else
       ind=-1;
     end
   end
   if ind>=0
     x       = pars.xc + pars.t*K*mu;
     fcheck1 = max(modeloPC(1).G*x - modeloPC(1).rhs);
     fcheck2 = max(modeloPC(2).G*x - modeloPC(2).rhs);
     %
     if fcheck1-fcheck2< model
        model  = fcheck1-fcheck2;
        fcheck = fcheck1;
        pars.x = x;
        Dualvar= mu;
     end
   end
end
return

function [pars,Dualvar,fcheck] = SubProbProximalDC(pars,modeloPC)
n    = pars.n;                     % dimensao do problema
ncut1 = length(modeloPC(1).rhs);       % numero de cortes f1
ncut2 = length(modeloPC(2).rhs);       % numero de cortes f2
model = inf;
Dualvar = ones(ncut1,1)/ncut1;
%-------------------------------------------------------------------------
% restricoes de igualdade 
Aeq = pars.X.Aeq;
beq = pars.X.beq;
if ~isempty(Aeq)
    m   = size(Aeq,1);
    Aeq = [Aeq, zeros(m,1)];
end
%-------------------------------------------------------------------------
% restricoes de desigualdade 
A = pars.X.A;
b = pars.X.b;
nineq = length(b);
if ~isempty(A)
    m   = size(A,1);
    A   = [A,        zeros(m,1);
           modeloPC(1).G, -ones(ncut1,1)];
    b   = [b;modeloPC(1).rhs];
else
    A   = [modeloPC(1).G, -ones(ncut1,1)];
    b   = [modeloPC(1).rhs];
end
%-------------------------------------------------------------------------
% limite inferior
lb = pars.X.lb;
if ~isempty(lb)
    lb = [lb;pars.flow];
else
    lb      = -inf(n+1,1);
    lb(n+1) = pars.flow;
end
%-------------------------------------------------------------------------
% limite superior
ub = pars.X.ub;
if ~isempty(ub)
    ub = [ub;inf];
else
    ub      = inf(n+1,1);
end
%-------------------------------------------------------------------------
% Funcao objetivo do QP
aux = 1/pars.t;
Q = diag(ones(n+1,1)*aux);
Q(n+1,n+1)=0;
m = length(b);
prob.Q   = 0.5*sparse(Q);
prob.A   = sparse([A;Aeq]);
prob.rhs = [b;beq];
m1 = length(prob.rhs);
sense=[];for i=1:m;sense=strcat(sense,'<');end
for j=m+1:m1; sense=strcat(sense,'=');end
prob.sense = sense;
prob.lb = lb;
prob.ub = ub;
par.OutputFlag =0; % 1 to print output
opts = optimset('display','off');
for i=1:ncut2
    %-------------------------------------------------------------------------
    % resolve o PL para obter o novo iterado
    q =  [-modeloPC(2).G(i,:)'-aux*pars.xc;1];
    if pars.matlab
        [x,~,ind,~,lambdap]  = quadprog(Q,q,A,b,Aeq,beq,lb,ub,[],opts);
        mu                         = lambdap.ineqlin;
        mu                         = mu(nineq+1:end);
    else
        clear lb ub Q A Aeq b beq
        prob.obj = q;
        res = gurobi(prob,par);
        if strcmp(res.status, 'OPTIMAL')|| strcmp(res.status, 'SUBOPTIMAL')
            x=res.x;
            ind = 1;
            mu  = -res.pi;
            mu  = mu(nineq+1:end);
        else
            ind=-1;
        end
%        prob.start = x;
    end
   if ind>=0
     fcheck1  = x(end);
     x       = x(1:end-1);
     fcheck2 = max(modeloPC(2).G*x - modeloPC(2).rhs);
          %
     if fcheck1-fcheck2< model
        model  = fcheck1-fcheck2;
        fcheck = fcheck1;
        pars.x = x;
        Dualvar= mu;
     end
   end
end
return