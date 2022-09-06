function [bpf, gamma] = evalBpfFromSamples( gSampleVals )

gValSorted = sort( gSampleVals, 'descend' );
gValSortedCumsum = cumsum( gValSorted );
if gValSortedCumsum(1) < 0 
    bpf = 0;
    gamma = inf;
elseif gValSortedCumsum(end) > 0
    bpf = 1;
    gamma = -inf;
else
    gammaId = find( gValSortedCumsum > 0 , 1, 'last' );
    nSample = length( gSampleVals );
    bpf = gammaId / nSample;
    
    gamma = gValSorted( gammaId+1 );
end