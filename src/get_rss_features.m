function [binSeg, binSim, binMats,rss] = get_rss_features(ts,nbins,bintype) 

if nargin < 2
    nbins = 10 ;
end

if nargin < 3
    bintype = 'quantiles' ;
end

% ets = get_ets(ts) ;
% rss = sqrt(sum(ets.^2,2)) ;
[ets,rss] = get_etsrss(ts) ; 

%Set up bins by distance 
switch bintype
    case 'sturges'
        [~,binedges] = histcounts(rss,'BinMethod','sturges') ;
    case 'quantiles'
        %bin by quantiles (bins contain equal numbers of edges)
        if nbins>2
            binedges=[min(rss) quantile(rss,nbins-1) max(rss)+eps]; 
        elseif nbins==2
            binedges=[min(rss) median(rss) max(rss)+eps];
        else
            binedges=[min(rss) max(rss)+eps];
        end
    case 'equalwidth'
        %bin into fixed-width bins 
        binedges=linspace(min(rss),max(rss)+eps,nbins+1);
    otherwise
        error('invalid bintype: %s',bintype)
end

%% segment the timepoints into bins
binSeg = discretize(rss,binedges) ;

%% get similarity to whole fc for each bin

binSim = nan(nbins,1) ;
binMats = cell(nbins,1) ;

origfc = corr(ts) ;
trium = logical(triu(ones(size(origfc)),1)) ;

for idx = 1:nbins

    dat = ets(binSeg==idx,:) ;
    mdat = mean(dat,1) ;
    
    binSim(idx) = corr(mdat(:),origfc(trium),'type','s') ;
   
    tmp = zeros(size(origfc)) ;
    tmp(trium) = mdat ; 
    binMats{idx} = tmp + tmp' ;

end
