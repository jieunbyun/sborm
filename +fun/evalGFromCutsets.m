function [Gsys, G_ind] = evalGFromCutsets( Gvals, cutset )
%{
Input:
- Gvals: (sample #) x (limit-state function #) array
- cutset: (cutset #)-dim cell array

Output:
- Gsys: (sample #) x 1 array
- G_ind: (sample #) x 1 array
%}

nSample = size(Gvals,1);
nCut = length(cutset);
G_cuts = zeros(nSample, nCut);
G_cuts_ind = zeros(nSample, nCut);

for iCutInd = 1:nCut
    iCutset = cutset{iCutInd};
    [G_cuts(:, iCutInd), iGcut_ind] = min( Gvals(:,iCutset), [] ,2 );

    if nargout > 1
        G_cuts_ind(:,iCutInd) = iCutset(iGcut_ind);
    end
end

[Gsys, Gsys_ind_max] = max( G_cuts, [], 2 );

if nargout > 1
    Gsys_ind_mat = sub2ind( size(G_cuts_ind), (1:nSample).', Gsys_ind_max );
    G_ind = G_cuts_ind(Gsys_ind_mat);
end
