function [Pmat,Pnorm,stateVals] = trans_probs(symbSeq,stateSz)

if nargin < 2
    stateSz = 'uniq' ;
end

nT = length(symbSeq) ;

switch stateSz
    case 'max'
        stateSz = max(symbSeq) ;
        stateVals = 1:stateSz ; 
        rectsymbseq = symbSeq ; 
    case 'uniq'
        stateVals = unique(symbSeq) ; 
        stateSz = length(stateVals) ;

        rectsymbseq = nan(nT,1) ;
        for idx = 1:stateSz
            tmp = stateVals(idx) ;
            rectsymbseq(symbSeq==tmp) = idx ;
        end
    otherwise % the user has specified a number
        stateVals = 1:stateSz ;
        rectsymbseq = symbSeq ; 
end

Pmat = zeros(stateSz) ;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% populate matrix
s1 = rectsymbseq(1) ;
for idx = 1:(nT-1)
    s2 = rectsymbseq(idx+1) ;
    Pmat(s1,s2) = Pmat(s1,s2) + 1;
    s1 = s2 ;
end

% normalize rows to add to one!
Pnorm = normalize(Pmat,2,"norm",1) ;
