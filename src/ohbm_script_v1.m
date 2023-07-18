%% OHBM 2023 script
% Josh Faskowitz, 2023, MIT License for this script and all functions
% contained here
% 
% need to download and install WSBM to run that stuff: 
% https://aaronclauset.github.io/wsbm/

%% fresh mind
clc 
clearvars

%% path stuff yo

addpath(genpath(pwd))
addpath('/SOMEWHERE/WSBM_v1.3/analysis tools/')

%% load up the data

load('../data/ts.mat') ; %loads a variable called ts
load('schaefer_2_yeo.mat') ;
yeosys = schaefer_2_yeo17.dictNoUkn('200') ;

%% loadup some singular data

[ntp,nnodes] = size(ts) ;
L = 1000 ;
TR = 0.72 ;
hrfdur = 2 ; 
onstate = [0.1 0.9] ; 

%% can we match the markov model to the empirical?

onD = logspace(log10(0.99),log10(0.5),20) ; 
offD = logspace(log10(0.01),log10(0.5),20) ; 
nreps = 100 ;
distMat = nan(length(onD),length(offD)) ; 

% first the empirical power
[~,tmp] = arrayfun(@(x_) quick_pspec(ts(:,x_),1/TR) , 1:size(ts,2) , 'UniformOutput' , false) ;
emp_spec = normalize(mean(cell2mat(tmp),2),'norm',1) ; 
L = size(ts,1) ; 

for idx = 1:length(onD) ; disp(idx)
    for jdx = 1:length(offD)

        on_d = onD(idx) ; 
        off_d = offD(jdx) ; 

        targMarkov = [on_d off_d ; onstate] ; 

        tmp_d = nan(nreps,1) ; 
        for kdx = 1:nreps

            tmp_ts = sim_bold_markov([L,1],TR,[0.008 0.08],hrfdur,targMarkov) ;
            
            [~,p1] = quick_pspec(tmp_ts,1/TR) ;
            tmp_d(kdx) = pdist([ emp_spec normalize(p1,'norm',1) ]') ; 

        end

        distMat(idx,jdx) = mean(tmp_d) ; 

    end
end

%% just take a looksie at it

[~,ii] = min(distMat,[],'all') ;
[a,b] = ind2sub([20 20],ii) ;

on_d = onD(a) ; 
off_d = offD(b) ; 

targMarkov = [on_d off_d ; onstate] ; 
rng(42)
tmp_ts = sim_bold_markov([L,1],TR,[0.008 0.08],hrfdur,targMarkov) ;
plot(tmp_ts)

%% first, why are there even spikes?
% how can two slow signals make faster spikes? we show here that it happens
% after the multiplication step! 

rng(42)

nreps = 200 ;
ncoups = 100 ; 
coupparams = 0:(1/ncoups):(1-(1/ncoups)) ;
powerspec_dist = nan(ncoups,100) ;

ts_p = nan((L/2)+1,ncoups) ;
ets_p = nan((L/2)+1,ncoups) ;
ts_env = nan(L,ncoups) ; 
ets_env = nan(L,ncoups) ; 

for idx = 1:ncoups

    disp(idx)

    ts_p_tmp = nan((L/2)+1,nreps) ;
    ets_p_tmp = nan((L/2)+1,nreps) ;

    ts_s_tmp = nan(L,nreps) ;
    ets_s_tmp = nan(L,nreps) ;

    for jdx = 1:nreps

        s1 = sim_bold_markov([L,2],TR,[0.008 0.08],hrfdur) ;
        sig = [ 1 coupparams(idx) ; coupparams(idx) 1 ] ;
        s2 = s1 * chol(sig) ;
    
        pp = prod(zscore(s2),2) ; 
    
        [~,p1] = quick_pspec(s1(:,1),1/TR) ;
        [~,p2] = quick_pspec(s1(:,2),1/TR) ;
  
%         [p1] = pspectrum(s1(:,1),1/TR) ;
%         [p2] = pspectrum(s1(:,2),1/TR) ;

        % mean power
        powTS = mean([p1 p2],2) ;
    
        [~,powETS] = quick_pspec(pp,TR) ;
%         [powETS,~] = pspectrum(pp,TR) ;

        [~,~,powerspec_dist(idx,jdx)] = kstest2(powETS,powTS) ; % KS statistic

        ts_p_tmp(:,jdx) = powTS ; 
        ets_p_tmp(:,jdx) = powETS ; 

%         ts_s_tmp(:,jdx) = s1  ;
        ets_s_tmp(:,jdx) = pp ;
        ts_s_tmp(:,jdx) = max([s1 s2],[],2) ;
    end

    ts_p(:,idx) = mean(ts_p_tmp,2) ;
    ets_p(:,idx) = mean(ets_p_tmp,2) ;
    
    ts_env(:,idx) = max(ts_s_tmp,[],2) ;
    ets_env(:,idx) = max(ets_s_tmp,[],2) ; 
    

end

%% and plot that

colmap = summer(ncoups) ;
colmap2 = winter(ncoups) ; 
timeVec = (1:L).*TR ;

tiledlayout(3,2)

rng(41)
nexttile([1 2])
% make a sample timeseries to plot
[s1,m1] = sim_bold_markov([L,1],TR,[0.008 0.08],hrfdur,targMarkov) ;
[f,p] = quick_pspec(s1,1/TR) ; 
% [p,f] = pspectrum(s1,1/TR) ; 

plot(timeVec,m1,'LineWidth',2)
hold on 
plot(timeVec,s1,'LineWidth',3,'Color',[colmap(1,:) 0.7])
hold off
xlim([0 720])
xlabel('Time (s)')
ylabel('Arbitary units')
title('Time series simulation')

nexttile

for idx = (1:ncoups)
%     plot(f,abs(ts_p(:,idx).^2)/L,'Color',[colmap(idx,:) 0.1],'LineWidth',2)
    plot(f,ts_p(:,idx),'Color',[colmap(idx,:) 0.1],'LineWidth',2)
    hold on
%     plot(f,abs(ets_p(:,idx).^2)/L,'Color',[ colmap(idx,:) 0.1],'LineWidth',2)
    plot(f,ets_p(:,idx),'Color',[ colmap(idx,:) 0.1],'LineWidth',2)

    xlim([0 0.15])
    %ylim([-0.00001 0.0004])
end
hold off
% and plot bandpass
line([ 0.008 0.008], [0.7 1],'Color','red','LineStyle','--')
line([ 0.08 0.08], [0.7 1],'Color','red','LineStyle','--')

cb = colorbar() ; colormap(colmap)
cb.Label.String = 'Coupling (r)' ; 
xlabel('Frequency (Hz)')
ylabel('Single-sided amp. of DFT')
title('Power spectrum TS & ETS')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

plot_smokey(coupparams,mean(powerspec_dist,2),std(powerspec_dist,[],2),...
    colmap(1,:))
xlabel('Time series coupling')
ylabel("Pow spect dist (KS statistic)")

title('Power spect. distance')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

% overall ylim
maxPlotYLim = [0 20] ;

for idx = (1:ncoups)
    plot(timeVec,ts_env(:,idx),'Color',[colmap(idx,:) 0.5],'LineWidth',2)
    hold on
end
hold off
cb = colorbar() ; colormap(colmap2)
cb.Label.String = 'Coupling (r)' ; 
xlabel('Time (s)')
ylabel('Signal max (a.u.)')
title('Time series max (both channels)')
ylim(maxPlotYLim)
xlim([0 720])

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nexttile

for idx = (1:ncoups)
    plot(timeVec,ets_env(:,idx),'Color',[colmap(idx,:) 0.5],'LineWidth',2)
    hold on
end
hold off
cb = colorbar() ; colormap(colmap)
cb.Label.String = 'Coupling (r)' ; 
xlabel('Time (s)')
ylabel('Signal max (a.u.)')
title('Edge time series max')
ylim(maxPlotYLim)
xlim([0 720])

%% now save
outfolder = './viz/ohbm_pdf/' ; 
mkdir(outfolder)

outfile = 'sim_n_coupling.pdf' ; 
orient(gcf,'landscape')
print(gcf,[outfolder outfile],'-dpdf','-vector','-bestfit')


%% NOW THE BLOCKMODEL!?!?
% you'll need it downloaded, installed, and in your path

addpath('/somewhere/WSBM_v1.3/')
addpath('/somewhere/WSBM_v1.3/analysis tools/')

wsbmCommCoup = 0.2:0.4:2  ;
wsbmReps = 100 ;
plantedComms = repmat(1:5,4,1) ; 
plantedComms = plantedComms(:);
% Generate Data Parameters
group_sizes = [5;5;5;5];
nnodes = sum(group_sizes) ;
R = [1,2,3,4; 
     2,1,4,3; 
     3,4,1,2;
     4,3,2,1]; % community structure
theta_e = [1; 1; 1; 1]; % make edge existence dense

reps = 200 ;
mod_vals = nan(reps,length(wsbmCommCoup)) ; 
rss_vals = nan(L,reps,length(wsbmCommCoup)) ; 

exampleMats = cell(5,1) ; 

for idx = 1:length(wsbmCommCoup)

    coupVal = wsbmCommCoup(idx) ;

    disp(idx)

    for jdx = 1:reps

        theta_w = [     coupVal,0.20      ;
                        0.05,0.10      ; 
                        0.0,0.05      ; 
                        0.0,0.01     ];
        [E,True_Model] = generateEdges('Normal','Bernoulli',R,theta_w,theta_e,group_sizes);
        synthMat = (Edg2Adj(E) + Edg2Adj(E)') ./ 2 ;
        synthMat(~~eye(size(synthMat,1))) = 1 ;
        synthTargCov = (synthMat*synthMat')^0.5 ;
        [V,D] = eig(synthTargCov) ;
        
        s1 = sim_bold_markov([L,20],TR,[0.008 0.08],5) ;
        synthTs = ((V*sqrt(D))*s1')' ;
        
        [~,mod_vals(jdx,idx)] = fcn_get_q(synthTargCov,plantedComms,'potts',1) ; 
        [~,rss_vals(:,jdx,idx)] = get_etsrss(synthTs) ; 

        if jdx==1
            exampleMats{idx} = synthTargCov ; 
        end

    end
end

%% plot the WSBM stuff 

colmap = winter(5) ;

set(gcf,'Position', [0 0 1200 600]);
TL = tiledlayout(5,4) ; 
mylayout = reshape(1:20,4,5)' ;

nexttile(2,[5 1])

rng(42)
reps = 200 ;
tmpCommC = repmat(wsbmCommCoup,reps,1) ; 
fcn_boxpts(mod_vals(:),tmpCommC(:),winter(5),0,cellstr(num2str(wsbmCommCoup')))
hold off
ylabel('Modularity (Potts, \gamma 1)')
xlabel('On-diagonal blockmodel weighting')
title('Modularity across synth. networks')

% tmpPlc = [ 10 8 6 4 2 ] ;
tmpPlc = flipud(mylayout(:,3)) ;
for idx = 1:length(wsbmCommCoup)
    nexttile(tmpPlc(idx),[1 2] )
    % plot(timeVec,rss_vals(:,:,idx),'Color',[ colmap(idx,:) 0.1],'LineWidth',1)
    % function [ h , pbin_spans , cmap ] = plot_manylines_aspatch(inlines,patchcolor,pbins,varargin)
    plot_manylines_aspatch({timeVec rss_vals(:,:,idx)},[ colmap(idx,:) ])
    hold on
    plot(timeVec,max(rss_vals(:,:,idx),[],2),'Color',[0.5 0.5 0.5])
    ylim([0 75])
    xlim([0 720])
    % also put the 95% line
    pp = prctile(rss_vals(:,:,idx),100,2) ;
    yline(mean(pp),':','Color',[1 0 0 0.25],'LineWidth',1.5)

    text(1.05,0.5,num2str(wsbmCommCoup(idx)),...
        'HorizontalAlignment','left','Rotation',0,'Units','normalized')
end
nexttile(tmpPlc(end))
title('Edge time series RSS')
nexttile(tmpPlc(1))
xlabel('Time (s)')
nexttile(tmpPlc(ceil(length(tmpPlc)/2)))
ylabel('RSS signal (a.u.)')

tmpPlc = flipud(mylayout(:,1)) ;
for idx = 1:5
    nexttile(tmpPlc(idx))

    imagesc(exampleMats{idx})
    clim([0 2])
    yticks([])
    xticks([])
    axis square
    text(-1,0.5,['On-diag val: ' num2str(wsbmCommCoup(idx))],...
        'HorizontalAlignment','left','Rotation',0,'Units','normalized')
end
nexttile(tmpPlc(end))
title('Exemplar networks')

%% 
outfolder = './viz/ohbm_pdf/' ; 
mkdir(outfolder)

% helps to fit paper
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
ppos(3:4) = pos(3:4);
% pos(1:2) = [1 1];
set(gcf,'paperposition',ppos);
set(gcf,'units',unis);

outfile = 'wsbm_sim.pdf' ; 
orient(gcf,'landscape')
print(gcf,[outfolder outfile],'-dpdf','-vector')

%% The decomposition work, to compare to previous papers
% based on decompositions

ctarg = cov(ts) ;
ptarg = mean((abs(fft(ts',[],2)).^2),1);

% get the ts
rng(42)
[~,simts_nocov] = simulate_BOLD_timecourse_func_v3(1100,0.72,0.72,ctarg,ptarg) ;

[V,D] = eig(cov(ts)) ;
[~,sortEigs] = sort(diag(D),'descend') ;
nbins = 10 ; 
nnodes = size(ts,2) ;

simRes = nan(nbins,nnodes) ; 
simRSS = nan(ntp,nnodes) ;
simFC = cell(nnodes) ;

for idx = 1:nnodes
    disp(idx)

    % get the largest eigvals
    sInds = sortEigs(1:(idx-1)) ;

    d = diag(D.*1) ;
    d(sInds) = 0 ; % nullify the eigenvectors
    
    redu_ts1 = (V*sqrt(diag(d))*simts_nocov')' ; % reduced timeseries

    [~,simRes(:,idx),~,rss] = get_rss_features(redu_ts1,nbins) ;

    simFC{idx} = corr(redu_ts1) ; 
    simRSS(:,idx) = rss ;
end

%% image it

% timeVec = (1:1100)' .* 0.72 ; 

set(gcf,'Position', [0 0 800 1200]);
tiledlayout(4,3)
mylayout = reshape(1:12,3,4)' ;

nexttile(1,[4 1])

h = imagesc(simRes(:,1:(nnodes/2))) ;
xlabel('# components reduced')
ylabel('Percentile bin')
labs = strtrim([cellstr(num2str((0:10:90)')) ...
    [ cellstr(repmat('--',9,1)) ; {'-'} ] ...
    cellstr(num2str((10:10:100)'))] ) ;
% yticks(1:10) ;
yticklabels(arrayfun(@(x_)strcat(labs{x_,:}),1:10,'UniformOutput',false))
cb = colorbar ;
cb.Label.String = "Rank corr. to data across bins" ; 
clim([0.5 1])

title('Amp. bin similarity to time-averaged corr.')

nexttile(mylayout(1,2),[2 2])

normRSS = normalize(simRSS,'norm',2) ;
imagesc(normRSS)
clim([0 0.13])
cb = colorbar() ;
cb.Title.String = 'Normalized RSS' ;
xlabel('# components reduced')
ylabel('Time (TR)')
title('ETS RSS after dim. reduction')

nexttile(mylayout(3,2),[1 2])

plot(timeVec,normRSS(:,1),'Color',[ 0.1 0.1 0.1 0.2],'LineWidth',1)
% hold on
% plot(timeVec,normRSS(:,1),'Color',[ 0.5 0.5 0.5 0.2],'LineWidth',5)
% hold off
ylabel('Normalized RSS')
title('ETS RSS')

nexttile(mylayout(4,2),[1 2])

colmap = flipud(parula(nnodes)) ;
for idx = fliplr(1:nnodes)
    plot(timeVec,normRSS(:,idx),'Color',[ colmap(idx,:) 0.2],'LineWidth',2)
    hold on 
end
xlim([0 timeVec(end)])
xlabel('Time (s)')
ylabel('Normalized RSS')
cb = colorbar ;
cb.Label.String = '# components redu.' ; 
cb.Ticks = 0:0.1:1 ;
cb.TickLabels = flipud(cellstr(num2str((0:20:200)'))) ;

title('ETS RSS after dim. redu. stacked')

%% 
outfolder = './viz/ohbm_pdf/' ; 
mkdir(outfolder)

% helps to fit paper
unis = get(gcf,'units');
ppos = get(gcf,'paperposition');
set(gcf,'units',get(gcf,'paperunits'));
pos = get(gcf,'position');
ppos(3:4) = pos(3:4);
% pos(1:2) = [1 1];
set(gcf,'paperposition',ppos);
set(gcf,'units',unis);

outfile = 'dim_redu.pdf' ; 
orient(gcf,'landscape')
print(gcf,[outfolder outfile],'-dpdf','-vector')


%% make system transition matrix

randlevels = [0 floor(logspace(log10(1),log10(200),7))] ;
nLevels = length(randlevels) ; 
reps = 100 ; 

transMats = cell(nLevels,1) ; 
rssVecs = cell(nLevels,1) ; 

for idx = 1:nLevels
    rr = randlevels(idx) ;
    
    disp([ num2str(idx) ' - ' num2str(rr)])

    if rr % randomization applied

        tmpPNorm = nan(17,17,reps) ;
        tmpRSS = nan(ntp,reps) ;

        for jdx = 1:reps
    
            disp(jdx)

            surrts = gen_phaserand_partial(ts,rr) ;
            newts = ts_2_sysload(surrts,yeosys) ;
            
            % make the average transprob
            
            [~,tmp] = arrayfun(@(x_) trans_probs(newts(:,x_),17),1:nnodes,'UniformOutput',false) ;
            % nanless = cellfun(@(c) fillmissing(c,'constant',0), tmp,'UniformOutput',false) ;
            avgpnorm = mean(carray_stack3(tmp),3,'omitnan') ;
    
            tmpPNorm(:,:,jdx) = avgpnorm ; 
    
            [~,~,~,tmpRSS(:,jdx)] = get_rss_features(surrts) ;

        end

        transMats{idx} = tmpPNorm ;
        rssVecs{idx} = tmpRSS ; 

    else % no randomization 

        newts = ts_2_sysload(ts,yeosys) ;

        [~,tmp] = arrayfun(@(x_) trans_probs(newts(:,x_),17),1:nnodes,'UniformOutput',false) ;
        avgpnorm = mean(carray_stack3(tmp),3,'omitnan') ;

        transMats{idx} = avgpnorm ; 
        [~,~,~,rssVecs{idx}] = get_rss_features(ts) ; 

    end

end

[~,actualRSS] = get_etsrss(ts) ;

%% 

% imsc_grid_comm(avgpnorm,1:17,2,[0.7 0.7 0.7],1,schaefer_2_yeo17.yeo17names)
% cb = colorbar() ;
% cb.Label.String = 'log10(probability)' ;
% % clim([-5 0])

%% 


%%

set(gcf,'Position', [0 0 1000 800]);
tiledlayout(3,5)

colmap = autumn(length(randlevels)*2) ; 
% timeVec = (1:1100)' .* 0.72 ; 

nexttile
plot(timeVec,actualRSS,'Color', [0.3 0.3 0.3 ])
title('Unperturbed')

for idx = 2:2:8

    nexttile

    % plot(timeVec,rssVecs{idx},'Color',[colmap(idx,:) 0.1])
    plot_manylines_aspatch({timeVec rssVecs{idx}},colmap(idx,:))
    hold on
    plot(timeVec,actualRSS,'Color', [0.1 0.1 0.1])
    xlim([0 792])
    title(randlevels(idx))

%     plot(timeVec,min(rssVecs{idx},[],2),'Color', [0.9 0.9 0.9])
%     plot(timeVec,max(rssVecs{idx},[],2),'Color', [0.9 0.9 0.9])

    hold off

    if idx == 4
        xlabel('Time (s)')
    end

end
    
nexttile
imsc_grid_comm(log10(avgpnorm),1:17,2,[0.7 0.7 0.7],[0.7 0.7 0.7 0.3],schaefer_2_yeo17.yeo17names)
clim([-3 0])
axis square

for idx = 2:2:8

    nexttile

    imsc_grid_comm(log10(mean(transMats{idx},3)),1:17,2,[0.7 0.7 0.7],[0.7 0.7 0.7 0.3])
    clim([-3 0])
    yticklabels([])   ;
    axis square

    colormap(hot)

    if idx == 4
        title('Transition matrices')
        xlabel('Canonical systems')
    end

end

cb = colorbar();
cb.Label.String = 'Transition probability (log10 scale)' ;

for idx = [1 2:2:8]
    nexttile

    if idx == 1
        tmp = avgpnorm ;
    else
        tmp = mean(transMats{idx},3) ;
    end
    s = scatter(avgpnorm(:),tmp(:),...
        50,repmat(colmap(idx,:),length(tmp(:)),1),...
        'filled') ;
    s.AlphaData = ones(length(tmp(:)),1) .* 0.2 ; 
    s.MarkerFaceAlpha = 'flat' ; 
    hold on
    % and highlight which are the on-diag
    s = scatter(diag(avgpnorm),diag(tmp)) ; 
    s.MarkerFaceAlpha = 1 ; 
    s.MarkerEdgeColor = [ 0.2 0.2 0.2 ] ; 
    s.LineWidth = 1 ;
    hold off

    line([0 1],[0 1])

    sqrt(sum((avgpnorm(:)-tmp(:)).^2))

    if idx == 1
        ylabel('Peturbed prob.')
    end

    if idx == 4
        title('Observed vs. Perturbed')
        xlabel('Observed prob.')
    end
end

%% 
outfolder = './viz/ohbm_pdf/' ; 
mkdir(outfolder)

outfile = 'trans_probs.pdf' ; 
print(gcf,[outfolder outfile],'-dpdf','-bestfit','-vector')

