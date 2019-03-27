function varargout = sh1_imaging(what, varargin);
%%  function varargout = sh1_imaging(what, varargin);
% Do imaging analyses of Yokoi&Diedricsen (2019)
%   - run pattern component modeling (pcm) on each cortical surface patch
%   - run representational clustering
%   - display results on flat surface
% 
% Note:
%   - for visualization (e.g., flattned map) smoothing is applied for
%       purely visual aesthetic purpose.
%   - no smoothing is applied for statistical analyses
%   - as a default setting, it only does visualization from pre-calculated results. 
% 
% To run analyses:
%   - when the data 'pcmResult.mat' is missing in the /data
%     folder, it runs pattern component modeling (PCM) using the pre-whitened
%     activity patterns for each surface patch ('SurfPatch_pwhBeta_s%02d.mat')
%   - when the data 'Cluster_10corr.mat' is missing in the /data folder, it
%     runs representational clustering using the pre-calculated RDMs for
%     each cortical patch ('SurfPatch_RDM.mat').
%   - when the data for cluster-wise noise-ceiling is missing (either
%     'noiseceiling_10cluster.mat' or 'noiseceiling_10cluster_relabel.mat')
% 	  it also runs PCM to estimate it from ('ClusterCorr10_pwhBeta_s%02d.mat').
% 
% Usage
%   - sh1_imaging(action, options)
%   - action can be following string inputs;
%       - 'Figure5'
%       - 'Figure6'
%       - etc.
%   - option can be 'bgcolor' and value, such as
%     sh1_imaging('Figure4c','bgcolor','k');
% 
% 
% ayokoi (2019) at.yokoi.work@gmail.com


%% Path etc.
% baseDir         = '/Volumes/G_Thunderbolt/Yokoi_Research/data/SequenceLearning/sh1/gittoshare/data'; 
% path to data/result (for release, uncomment following lines)
thisfile = mfilename('fullpath');
baseDir = fileparts(thisfile);
p={fullfile(baseDir,'code'), fullfile(baseDir,'code','caret')};
addpath(p); % path to necessary functions

% File structure
F.baseDir = baseDir;
F.pwhBeta_individ = fullfile(baseDir,'SurfPatch_pwhBeta_s%02d.mat'); % pre-whitened beta for each cortical surface patch
% As Github has 100MB size limitation for single file, concatenated pre-whitened activity data
% for all participants is not saved.
F.RDM = fullfile(baseDir, 'SurfPatch_RDM.mat'); % pre-calculated RDM for all patches, for all subjects
F.modeldesign = fullfile(baseDir,'pcmDesign.mat'); % representational models and experimental design info (necessary for pcm)
F.pcmResult = fullfile(baseDir,'pcmResult.mat'); % result file for pcm
F.coord = {fullfile(baseDir,'lh.FLAT.coord'), fullfile(baseDir,'rh.FLAT.coord')}; % flattened surface (fsaverage_sym)
F.topo = {fullfile(baseDir,'lh.CUT.topo'), fullfile(baseDir,'rh.CUT.topo')}; % flattened surface (fsaverage_sym)
F.shape = {fullfile(baseDir,'lh.surface_shape'),fullfile(baseDir,'rh.surface_shape')}; % flattened surface (fsaverage_sym)
F.border = {fullfile(baseDir,'lh.All.borderproj'),fullfile(baseDir,'rh.All.borderproj')}; % flattened surface (fsaverage_sym)
F.p2n = fullfile(baseDir,'patch2node.mat'); % node assignment for each patch
F.meandist = fullfile(baseDir, {'lh.meanDist.metric','rh.meanDist.metric'}); % mean crossnobis distance computed on continuous searchlight
F.pwhBeta_cluster_individ = fullfile(baseDir, 'ClusterCorr10_pwhBeta_s%02d.mat'); % prewhitened beta for cluster (10 solutions)
F.cluster = fullfile(baseDir, 'Cluster_10corr.mat');
F.clustercolor = fullfile(baseDir,'10ClusterColors.mat');
F.ceilingCluster = fullfile(baseDir,'noiseceiling_10cluster.mat');
F.ceilingCluster_relabel = fullfile(baseDir,'noiseceiling_10cluster_relabel.mat');
F.design_relabel = fullfile(baseDir,'pcmDesign_visualcue.mat');

% set path for functions

% other parameters for figure
plotrange = {{[-79.8252,104.7103],[-64.1677,81.5173-10]},...
    {[-104.7103,79.8252],[-64.1677,81.5173-10]}};
goldorange  = [243,152,0]/255;
magenta     = [228,0,127]/255;
pink = [255,192,203]/255;
cyan        = [0,104,183]/255;
maxclust = 10;
clabel10 = [6,5,2,7,9,1,3,8,4,10]; % labeling for cluster (10 solution)

bgcolor = 'w'; % background color
vararginoptions(varargin,{'bgcolor'});
switch lower(bgcolor)
    case {'w',[1,1,1]};
        bgcolor=[1,1,1];
        axcolor=[0,0,0];
    case {'k',[0,0,0]};
        bgcolor=[0,0,0];
        axcolor=[1,1,1];
    otherwise
        bgcolor=[1,1,1];
        axcolor=[0,0,0];
end

%% Main
switch (what)
    case 'Figure4c'             % Plot mean distance map
        % load result
        for h=1:2
            M = caret_load(F.meandist{h}); % mean crossnobis distance on continuous searchlight
            crossnobis = M.data(:,1);
            crossnobis(isnan(crossnobis)) = 0; % pad
            nodedata(:,h) = crossnobis;            
        end
        % map mean distance
        plotrange = {{[-126, 140],[-124,82]}, {[-140, 126],[-124,82]}}; % all surface
        map_fsaverage(nodedata, F.coord, F.topo, F.shape, F.border,...
            'threshold', [0.03,0.10], 'MAP', 'parula','plotrange', plotrange,...
            'label','mean crossnobis dist.','bgcolor',bgcolor);
        % save?
        varargout = {};
    case 'Figure5c'             % Plot noise-ceiling map
        % load result
        R = sh1_imaging('run_pcm');
        
        load(F.p2n); % patch-node mapping
        
        % condence to plot group mean
        R = tapply(R, {'hemis','patch','nodeID'}, {'noiseceiling','nanmean(x,1)', 'name', 'noiseceiling'},...
            {'pxp_ceiling','nanmean(x,1)','name','pxp'});
        
        % map pxp and get mask
        patchdata = R.pxp;
        nodedatapxp = assign_patch2node(patchdata, R.hemis, R.patch, SLnodes,'smooth',1,'F',F);
        nodedatapxp(nodedatapxp<=0.75)=NaN;
        map_fsaverage(nodedatapxp, F.coord, F.topo, F.shape, F.border,...
            'threshold', [0.75,1.0], 'MAP', 'autumn','plotrange', plotrange,...
            'label','PXP(noise-ceiling)','distfile',F.meandist,'bgcolor',bgcolor);
        
        % map logBF (noise-ceiling)
        patchdata = R.noiseceiling;
        nodedatabf = assign_patch2node(patchdata, R.hemis, R.patch, SLnodes,'smooth',1,'F',F);
        nodedatabf(nodedatapxp<=0.75)=NaN; % mask with pxp>0.75
        map_fsaverage(nodedatabf, F.coord, F.topo, F.shape, F.border,...
            'threshold', [1,10], 'MAP', 'jet','plotrange', plotrange,...
            'label','logBF(noise-ceiling)','distfile',F.meandist,'bgcolor',bgcolor);
        
        % print figure
        
        varargout = {};
    case 'Figure6abc'           % Plot logBFc map
        % load result
        R = sh1_imaging('run_pcm');
        
        load(F.p2n); % patch-node mapping
        
        % condence to plot group mean
        R = tapply(R, {'hemis','patch','nodeID'}, {'logBFc','nanmean(x,1)', 'name', 'logBFc'},...
            {'pxp','nanmean(x,1)','name','pxp'});
        
        % do flatmap
        models = {'First-finger','Chunk','Sequence'};
        for i=1:3
            % map pxp and get mask
            patchdata = R.pxp(:,i);
            nodedatapxp = assign_patch2node(patchdata, R.hemis, R.patch, SLnodes,'F',F,'smooth',1);
            map_fsaverage(nodedatapxp, F.coord, F.topo, F.shape, F.border,...
                'threshold', [0.75,1.0], 'MAP', 'autumn','plotrange', plotrange,...
                'label',sprintf('PXP(%s)',models{i}),'bgcolor',bgcolor);
            
            % map logBF (noise-ceiling)
            patchdata = R.logBFc(:,i);
            %patchdata(R.pxp(:,i)<0.75)=-inf;
            nodedatabfc = assign_patch2node(patchdata, R.hemis, R.patch, SLnodes,'F',F,'smooth',1);
            nodedatabfc(nodedatapxp<=0.75) = NaN;
            map_fsaverage(nodedatabfc, F.coord, F.topo, F.shape, F.border,...
                'threshold', [1,3], 'MAP', 'parula','plotrange', plotrange,...
                'label',sprintf('logBFc(%s)',models{i}),'distfile',F.meandist,'bgcolor',bgcolor);
        end
    case 'Figure6d'             % Plot merged map
        % load
        if exist(F.pcmResult, 'file');
            R = load(F.pcmResult);
        else
            R = sh1_imaging('run_pcm');
        end
        load(F.p2n); % patch-node mapping
        
        % condence to plot group mean
        R = tapply(R, {'hemis','patch','nodeID'}, {'logBFc','nanmean(x,1)', 'name', 'logBFc'},...
            {'pxp','nanmean(x,1)','name','pxp'});
        
        % do merged flatmap
        models = {'First-finger','Chunk','Sequence'};
        for i=1:3
            nodedatapxp = assign_patch2node(R.pxp(:,i), R.hemis, R.patch, SLnodes,'smooth',1,'F',F);
            % map logBF (noise-ceiling)
            patchdata = R.logBFc(:,i); %.*double(R.pxp(:,i)>0.75);
            nodedatabfc = assign_patch2node(patchdata, R.hemis, R.patch, SLnodes,'smooth',1,'F',F);
            nodedatabfc(nodedatapxp<=0.75)=NaN; % mask with pxp>0.75
            nodedata(:,:,i) = nodedatabfc;
        end
        nodedata = permute(nodedata, [1,3,2]);
        MAPs = {repmat(cyan,100,1), repmat(magenta,100,1),repmat(goldorange,100,1)};
        map_fsaverage_multi(nodedata, F.coord, F.topo, F.shape, F.border,...
            'threshold', [1,3], 'MAP', MAPs,'plotrange', plotrange,'label',sprintf('logBFc(merged)'),...
            'distfile',F.meandist,'bgcolor',bgcolor);
    case 'Figure6fg'           % Plot scatter plot (overlap analyses)
        % load result
        R = sh1_imaging('run_pcm');
        
        % condence to plot group mean
        R = tapply(R, {'hemis','patch','nodeID'}, {'noiseceiling','nanmean(x,1)', 'name', 'noiseceiling'},...
            {'logBFc','nanmean(x,1)','name','bfc'},{'pxp','nanmean(x,1)','name','pxp'});
        
        % cases
        pxpok(:,1) = R.pxp(:,1)>0.75; % first finger
        pxpok(:,2) = R.pxp(:,2)>0.75; % chunk
        pxpok(:,3) = R.pxp(:,3)>0.75; % sequence
        bfcok(:,1) = R.bfc(:,1)>1; % first finger
        bfcok(:,2) = R.bfc(:,2)>1; % chunk
        bfcok(:,3) = R.bfc(:,3)>1; % sequence
        present = pxpok.*bfcok;
        
        category=present(:,1)+present(:,2)*2+present(:,3)*4; % Binary-based code 0-7
        % 0: None
        % 1: F
        % 2: C
        % 3: F+C
        % 4: S
        % 5: F+S
        % 6: C+S
        % 7: F+C+S
        
        % test overlap
        pivottable(category,[],category,'length','subset',ismember(R.hemis,[1,2]) & R.noiseceiling>1);
        [M,C]=pivottable(category,[],category,'length','subset',ismember(R.hemis,[1,2]) & R.noiseceiling>1);
        % first-finger and others
        [G_f, p_f, exp_f, obs_f] = Gtest([M(C==0), M(C==1), sum(M(C==2|C==4|C==6)), sum(M(C==3|C==5|C==7))]);
        % chunk and sequence
        [G_cs, p_cs, exp_cs, obs_cs] = Gtest([M(C==0), M(C==2), M(C==4), M(C==6)]);
        stats.name = {'first-finger vs. others','chunk vs. sequence'}';
        stats.expected_overlap = [exp_f(end);exp_cs(end)];
        stats.observed_overlap = [obs_f(end);obs_cs(end)];
        stats.Gvalue = [G_f; G_cs];
        stats.Pvalue = [p_f; p_cs];
        struct2table(stats)
        
        % figure;
        figsize=[25*3/2,15/sqrt(2)];
        handle = figure('unit','centimeters',...
            'papersize',figsize    ,...
            'position',[1,1,figsize(1),figsize(2)],...
            'paperposition',[0,0,figsize(1),figsize(2)],...
            'paperpositionmode','manual',...
            'Color',bgcolor,'InvertHardcopy','off');
        
        subplot(1,3,1); % venn diagram
        axis square; 
        title({'Venn diagram', '(not handled for this implementation)'},'color',axcolor);
        axis off;
        % You can draw it manually
        
        % Scatter plot
        subplot(1,3,2); % ff vs seq
        set(gca,'xcolor',axcolor);
        colors = {cyan, magenta, pink, goldorange};
        catname = {'F', 'C', 'F+C', 'C+S', 'S'};
        catname = {'F', 'C', 'F+C', 'S', 'S+F', 'C+S', 'F+C+S'};
        colors = {cyan, magenta, 'k', goldorange, 'k', pink, 'w'};
        subset = R.noiseceiling>1;
        R=getrow(R, category>0); % remove 'N' patches
        subset=subset(category>0);
        category=category(category>0);
        scatterplot(R.bfc(subset,1),R.bfc(subset,3),...
            'split', category(subset),...
            'markertype','o',...
            'markersize', {10,10,10,10},...
            'markerfill',colors(C(C>0)),...
            'markercolor',{axcolor},...
            'leg',catname(C(C>0)),'leglocation','best'); box off; hold on;
        axis square; %title('title');
        xlabel('logBFc (first-finger)','color',axcolor); 
        ylabel('logBFc (sequence)','color',axcolor);
        set(gca, 'tickdir', 'out', 'ticklength', [0.03, 0.01]);
        set(gca,'xcolor',axcolor,'ycolor',axcolor,'color','none');
        
        subplot(1,3,3); % ch vs seq
        scatterplot(R.bfc(:,2),R.bfc(:,3),...
            'split', category.*subset,...
            'markertype','o',...
            'markersize', {10,10,10,10},...
            'markerfill',colors(C(C>0)),...
            'markercolor',{axcolor}); box off; hold on;
        axis square; %title('title');
        xlabel('logBFc (chunk)','color',axcolor); 
        ylabel('logBFc (sequence)','color',axcolor);
        set(gca, 'tickdir', 'out', 'ticklength', [0.03, 0.01]);
        set(gca,'xcolor',axcolor,'ycolor',axcolor,'color','none');
        
        varargout = {stats};
    case 'Figure7b'             % Plot clustering result map
        % load clustering result
        [clusters,Lap,LapEigVec,Idx] = sh1_imaging('run_clustering');
        
        % load color (for more or less clusters than 10, colors should be additionally defined)
        load(F.clustercolor);
        
        % map
        load(F.p2n); % patch-node mapping
        for i=1:max(clusters)
            nodedata(:,:,i) = assign_patch2node(double(clusters==i), Idx.hemis, Idx.patch, SLnodes,'smooth',1,'F',F);
            MAPs{i} = colors(i,:);
        end
        nodedata = permute(nodedata, [1,3,2]);
        map_fsaverage_multi(nodedata, F.coord, F.topo, F.shape, F.border,...
            'threshold', [0.7,1], 'MAP', MAPs,'plotrange', plotrange,'label',sprintf('Clusters'),...
            'distfile',F.meandist,'alpha',0.9);
    case 'Figure7c'             % Plot Marimekko chart (mosaic plot)
        % load clustering result
        [clusters,Lap,LapEigVec,Idx] = sh1_imaging('run_clustering');
        Idx.cluster=clabel10(clusters)';
        Idx=tapply(Idx,{'hemis','patch','nodeID'}, {'cluster','nanmean','name','cluster'});
        
        % load pcm result
        if exist(F.pcmResult, 'file');
            R = load(F.pcmResult);
        else
            R = sh1_imaging('run_pcm');
        end
        % count cases
        R = tapply(R, {'hemis','patch','nodeID'}, {'noiseceiling','nanmean(x,1)', 'name', 'noiseceiling'},...
            {'logBFc','nanmean(x,1)','name','bfc'},{'pxp','nanmean(x,1)','name','pxp'});
        pxpok(:,1) = R.pxp(:,1)>0.75; % first finger
        pxpok(:,2) = R.pxp(:,2)>0.75; % chunk
        pxpok(:,3) = R.pxp(:,3)>0.75; % sequence
        bfcok(:,1) = R.bfc(:,1)>1; % first finger
        bfcok(:,2) = R.bfc(:,2)>1; % chunk
        bfcok(:,3) = R.bfc(:,3)>1; % sequence
        present = pxpok.*bfcok;
        
        category=present(:,1)+present(:,2)*2+present(:,3)*4; % Binary-based code 0-7
        % 0: None
        % 1: F
        % 2: C
        % 3: F+C
        % 4: S
        % 5: F+S
        % 6: C+S
        % 7: F+C+S
        category(category==6)=3; % just for reordering
        
        % plot
        edgecolor = 'none';
        linewidth = 1;
        row=category;
        col = Idx.cluster;
        ilc = R.noiseceiling>1;
        facecolor = {cyan, magenta, pink, goldorange, [0.8,0.8,0.8]};
        rowleg = {'F','C','C+S','S'};
        for c=1:max(clusters);
            colleg{c} = sprintf('%d', c);
        end
        colleg{7} = ''; % for preventing visual cluttering
        colleg{9} = '';
        figure('color',bgcolor,'inverthardcopy','off');
        set(gca,'xcolor',axcolor,'color','none');
        [M, Rh,Ch,xcenter,xrange]=mosaicplot(row, col, ones(size(row)),...
            'facecolor', facecolor,...
            'rowlabel', [],...
            'columnlabel', [],...
            'columnlabelrotation', 0,...
            'edgecolor', edgecolor,...
            'linewidth',linewidth,...
            'textcolor','w',...
            'subset', row>0&ilc,...
            'leg', rowleg,...
            'leglocation','southoutside');
        title({'Cortical "storage" of sequence representations',''},'color',axcolor);
        axis on;
        xlabel('Clusters','color',axcolor);
        set(gca,'ycolor','none','xtick',xcenter,'xticklabel',colleg,'xlim',xrange,...
            'tickdir','out','ylim',[-0.01,1.05]);
    case 'Figure7d'             % Plot relabeling result
        % load relabeling result
        S=sh1_imaging('est_noiseceiling_cluster_relabel');
        S.cluster=clabel10(S.cluster)';
        
        % load original result
        US=sh1_imaging('est_noiseceiling_cluster');
        US.cluster=clabel10(US.cluster)';
        
        % load color (for more or less clusters than 10, colors should be additionally defined)
        load(F.clustercolor);
        for c=1:maxclust
            cols{clabel10(c)} = colors(c,:);
            label{c} = sprintf('%d',c);
        end
        % plot
        figure('color',bgcolor,'inverthardcopy','off');
        CAT.facecolor = cols;
        CAT.edgecolor = axcolor;
        CAT.linecolor = axcolor;
        CAT.linewidth = 1;
        CAT.errorcolor=axcolor;
        Y=S.noiseceiling./US.noiseceiling;
        
        [xpos,ypos,e]=barplot(S.cluster, Y, 'subset', isfinite(Y),...
            'CAT', CAT,'gapwidth', [0.5,0.5],'split',S.cluster,'capwidth', 0.001);
        drawline(1,'dir','horiz','color',axcolor,'linewidth',1);
        ylabel({'logBF-ratio'}, 'color',axcolor);
        xlabel({'Clusters'}, 'color',axcolor);
        title('Relabeling impact on noise-ceiling','color',axcolor);
        set(gca,'ytick',[0:0.5:2],'ylim',[0,2],...
            'xtick',[xpos],'xticklabel',label,'xcolor',axcolor,'ycolor',axcolor,'color','none');
        set(gca,'tickdir','out','ticklength',[0.02,0.01]);
    
    case 'following parts could have been local functions instead'
    case 'concat_pwhbeta'   % Concatenate pre-whitened beta for all subjects into single file
        T = [];
        for s=1:12
            D = load(sprintf(F.pwhBeta_individ, s));
            T = addstruct(T,D);
        end;
        varargout = {T};
    case 'run_pcm'              % Run pcm on each surface patch
        if exist(F.pcmResult, 'file');
            R = load(F.pcmResult);
        else
            % load data
            Data = sh1_imaging('concat_pwhbeta');%load(F.pwhBeta); % concatenated file
            
            % load design
            load(F.modeldesign); % FCS family
            Ncondition = length(partitionVec{1});
            
            % run group-cross-validation at each patch
            patches = unique(Data.patch)'; % this is not surface node ID
            hemis = [1:2];
            sn = unique(Data.sn)';
            R = [];
            for h=hemis;
                for p=patches;
                    idx=Data.patch==p&Data.hemis==h;
                    
                    pvec=[]; cvec=[]; Yprewh=[];c=0;
                    if sum(idx)>0
                        D=getrow(Data, idx);
                        for s=sn
                            Y=D.pwhBeta{s};
                            if ~isempty(Y)
                                c=c+1;
                                Yprewh{c,1}=Y(1:Ncondition,:); % remove run intercept
                                pvec{c,1}=partitionVec{s};
                                cvec{c,1}=conditionVec{s};
                            else
                            end
                        end
                    end
                    if ~isempty(Yprewh)&&numel(Yprewh)>6;
                        [T,Tcv,theta_cv] = runPCM(Yprewh, cvec, pvec, MF);
                        
                        % calculate normal logBF
                        logBF = bsxfun(@minus, Tcv.likelihood, Tcv.likelihood(:,1));
                        noiseceiling = logBF(:,end);
                        logBF = logBF(:,1:end-1);
                        
                        % calculate component logBF
                        [PP, logBFc] = pcm_componentPosterior(logBF, CompIdx);
                        
                        % get PXP (spm_BMS) for logBFc
                        for c=1:3
                            lme = [logBFc(:,c), zeros(size(logBFc(:,2)))];
                            [~,~,~,pxp_] = spm_BMS (lme);
                            pxp(:,c) = repmat(pxp_(1), size(lme,1),1);
                        end
                        
                        % get PXP for noise-ceiling logBF (vs null)
                        lme = [noiseceiling, zeros(size(noiseceiling))];
                        [~,~,~,pxp_] = spm_BMS (lme);
                        pxp_ceiling(:,1) = repmat(pxp_(1), size(lme,1),1);
                        
                        % concatenate results
                        r.logBF = logBF;
                        r.noiseceiling = noiseceiling;
                        r.pxp_ceiling = pxp_ceiling;
                        r.logBFc = logBFc;
                        r.pxp = pxp;
                        r.sn = repmat(s, length(pxp),1);
                        r.hemis = repmat(h, length(pxp),1);
                        r.patch = repmat(p, length(pxp),1);
                        r.nodeID = repmat(D.nodeID(1), length(pxp),1);
                        R = addstruct(R,r);
                    end
                end;
            end;
            % save result
            save(F.pcmResult, '-struct', 'R');
        end
        varargout = {R};
    case 'run_clustering'       % Run spectral clustering
        distfun='correlation'; % distance function defining adjacency matrices across patches
        W=[];
        if exist(F.cluster,'file')
            load(F.cluster);
        else % do clustering
            % load RDM on each patch
            T = load(F.RDM);
            
            for s=1:12
                D = getrow(T,T.sn==s);
                D = tapply(D, {'nodeID', 'hemis', 'patch'},...
                    {'distance','nanmean(x,1)','name','distance'});
                
                % derive individual w
                Dist = D.distance; % D.data is now Nroi x 28xNsubj matrix
                
                % calculate adjacency matrix
                tmp = squareform(pdist(Dist,distfun));
                M=tmp; % now in the form of distance
                % transform into similarity
                sigma = quantile(squareform(tmp),0.05);
                M = exp(-M.^2 ./ (2*sigma^2)); % Gaussian similarity transformation
                W=cat(3,W,M);
                fprintf('...done.\n');
            end
            
            % take group-average first rather than individual clustering
            W_avrg = nanmean(W,3);
            
            % do clustering
            [clusters, Lap, LapEigVec] = clusterRDM_spectralClustering(W_avrg, maxclust);
            
            % save result
            Patch.hemis = D.hemis;
            Patch.patch = D.patch;
            Patch.nodeID = D.nodeID;
            save(F.cluster, 'clusters','Lap','LapEigVec','Patch');
        end
        varargout = {clusters,Lap,LapEigVec,Patch};
    case 'est_noiseceiling_cluster'       % Run pcm (estimate noise-ceiling)
        if exist(F.ceilingCluster, 'file');
            S=load(F.ceilingCluster);
        else % do it from scratch
            % load prewhitened beta for cluster and concat across subjects
            Data = [];
            for s=1:12
                D=load(sprintf(F.pwhBeta_cluster_individ,s));
                Data=addstruct(Data,D);
            end
            
            % load design and model
            load(F.modeldesign);
            MF=MF(1); % leave only null model
            
            % load relabelled design
            %load(F.design_relabel); % this adds 'conditionVec' and 'partitionVec'
            
            % do PCM to estimate noise-ceiling
            S = [];
            for c=1:maxclust
                data=getrow(Data, Data.cluster==c);
                Yprewh=data.beta;
                [~,Tcv] = runPCM(Yprewh, conditionVec, partitionVec, MF);
                
                % calc ceiling
                like = bsxfun(@minus, Tcv.likelihood, Tcv.likelihood(:,1));
                s.noiseceiling = like(:,end);
                s.SN = Tcv.SN;
                s.cluster = repmat(c,size(s.SN));
                S = addstruct(S,s);
            end
            S.sort = ones(size(S.SN));
            
            % save relabelling result
            save(F.ceilingCluster, '-struct','S','-v7.3');
        end
        varargout = {S};
    case 'est_noiseceiling_cluster_relabel'       % Run pcm (estimate noise-ceiling) with relabeling
        if exist(F.ceilingCluster_relabel, 'file');
            S=load(F.ceilingCluster_relabel);
        else % do it from scratch
            % load prewhitened beta for cluster and concat across subjects
            Data = [];
            for s=1:12
                D=load(sprintf(F.pwhBeta_cluster_individ,s));
                Data=addstruct(Data,D);
            end
            
            % load design and model
            load(F.modeldesign);
            MF=MF(1); % leave only null model
            
            % load relabelled design
            load(F.design_relabel); % this adds 'conditionVec' and 'partitionVec'
            
            % do PCM to estimate noise-ceiling
            S = [];
            for c=1:maxclust
                data=getrow(Data, Data.cluster==c);
                Yprewh=data.beta;
                [~,Tcv] = runPCM(Yprewh, conditionVec, partitionVec, MF);
                
                % calc ceiling
                like = bsxfun(@minus, Tcv.likelihood, Tcv.likelihood(:,1));
                s.noiseceiling = like(:,end);
                s.SN = Tcv.SN;
                s.cluster = repmat(c,size(s.SN));
                S = addstruct(S,s);
            end
            S.sort = ones(size(S.SN));
            
            % save relabelling result
            save(F.ceilingCluster_relabel, '-struct','S','-v7.3');
        end
        varargout = {S};        
    otherwise
        warning('No such case.');
end
end

% local functions
function varargout = runPCM(Yprewh, conditionVec, partitionVec, MF, varargin) % run PCM on prewhitened beta for each patch
runEffect   = 'fixed'; % run effect
verbose     = 1;
fitAlgorithm= 'NR';
numIter = 1e3;
pcm_vararginoptions(varargin,{'runEffect', 'verbose','fitAlgorithm', 'prior','numIter'});

% append noise ceiling model
MF{end+1}.type = 'freedirect';
MF{end}.numGparams = 0;
MF{end}.theta0 = [];
MF{end}.name = 'noise_ceiling';

% Run group-fit
[T, theta] = pcm_fitModelGroup(Yprewh, MF, partitionVec, conditionVec,...
    'runEffect', runEffect, ...
    'fitAlgorithm', 'NR', ...
    'verbose', verbose, 'MaxIteration', numIter);

% Do the crossvalidated group-fit
[Tcv, theta_cv] = pcm_fitModelGroupCrossval(Yprewh, MF, partitionVec, conditionVec,...
    'runEffect', runEffect, ...
    'fitAlgorithm', fitAlgorithm, ...
    'groupFit', theta, 'MaxIteration', numIter);

varargout = {T, Tcv, theta_cv};
end
function varargout = map_fsaverage(data, coords, topos, shapes, border, varargin) % map result on fsaverage_sym template 
% map data onto fsaverage flat surface
threshold = [0,0.2];
alpha=0.7;
bordersize=10;
fontsize = 15;
bgcolor = [0 0 0]; % must be vector
MAP = 'parula';
label = '';
distfile = [];
plotrange = {{[-126, 140],[-124,82]}, {[-140, 126],[-124,82]}};
%isSave = 0;
vararginoptions(varargin(1:end), {'threshold', 'bgcolor','template','MAP','label','isSave','plotrange','distfile'});

borderalighment={'left','right'};

% loop over hemisphere
[axespositions, hmap, hscale] = setupFlatMap(plotrange,1,1,bgcolor); % make panels
figure(hmap);
for h=1:2
    % show result on flatmap
    coord = coords{h};
    topo = topos{h};
    shape = shapes{h};
    
    S=caret_load(shape);
    depth = S.data(:,2); % cortical depth information
    
    F=caret_load(coord);
    xlims = [min(F.data(:,1)), max(F.data(:,1))]; % plot range (x)
    ylims = [min(F.data(:,2)), max(F.data(:,2))]; % plot range (y)
    
    %subplot(1,2,h);
    ax(h) = axes('position',axespositions{h});    
    if ~isempty(distfile)
        M=caret_load(distfile{h});
        caret_plotflatmap_rgb('coord',coord,'topo',topo,...
            'underlay', depth, 'data', M.data(:,1), 'dscale', [0.03,10],'threshold',0.03,...
            'xlims', xlims, 'ylims', ylims,'dmap',hot(100),'alpha',alpha);
        idx = ~isnan(data(:,h))&data(:,h)>threshold(1);
        [~,p]=caret_plotflatmap_rgb('coord',coord,'topo',topo,...
            'underlay', depth, 'data', data(:,h), 'dscale', threshold,'threshold',threshold(1),...
            'xlims', xlims, 'ylims', ylims,'dmap',eval(MAP),'alpha',alpha,'idx',idx);
        %set(p,'facealpha',alpha);
    else
        idx = ~isnan(data(:,h));
        caret_plotflatmap_rgb('coord',coord,'topo',topo,...
            'underlay', depth, 'data', data(:,h), 'dscale', threshold,'threshold',threshold(1),...
            'xlims', xlims, 'ylims', ylims,'dmap',eval(MAP),'alpha',alpha,'idx',idx);
    end
    
    axis equal; axis tight; axis off; hold on;
    set(gca, 'xlim',plotrange{h}{1},'ylim',plotrange{h}{2});
    
    % draw sulci onto flatmap
    B = caret_load(border{h});
    for b=1:length(B.Border)
        %borderidx = cat(1,B.Border.vertex);
        borderidx = cat(1,B.Border(b).vertex);
        borderX = F.data(borderidx(:,1),1);
        borderY = F.data(borderidx(:,1),2);
        plot(borderX,borderY,'w.', 'markersize', bordersize);
        text(borderX(end),borderY(end), B.Border(b).name,...
            'horizontalalignment',borderalighment{h},...
            'verticalalignment','top',...
            'fontsize',fontsize,'color','w'); % border name
    end
end
% draw colorbar
figure(hscale); colormap(MAP);
ax(end+1)=axes('position', axespositions{end});
colorscale = linspace(threshold(1),threshold(2),100);
imagesc(colorscale); caxis(threshold);
set(gca,'xcolor',[1,1,1]-bgcolor,'ycolor',[1,1,1]-bgcolor,'fontsize',fontsize,'tickdir','out');
set(gca,'ytick',[],'xtick',[1,100],'xticklabel',{num2str(threshold(1)), num2str(threshold(2))});
title(label,'color',[1,1,1]-bgcolor,'fontsize',fontsize);

% save figures
% fname_map = sprintf('Map.%s.%s.png',template,label);
% fname_scale = sprintf('Colorbar.%s.%s.png',template,label);
% mySaveFig(hmap, fullfile(figDir, fname_map), isSave, '-dpng','-r400');
% mySaveFig(hscale, fullfile(figDir, fname_scale), isSave, '-dpng','-r400');

varargout = {};
end
function varargout = map_fsaverage_multi(data, coords, topos, shapes, border, varargin) % map multiple results on fsaverage_sym template
% map multiple data onto fsaverage flat surface to make merged image
threshold = [0,0.2];
alpha=0.7;
bordersize=10;
fontsize = 15;
bgcolor = [0 0 0]; % must be vector
MAP = parula(100);
label = '';
plotrange = {{[-126, 140],[-124,82]}, {[-140, 126],[-124,82]}};
%isSave = 0;
distfile = [];
vararginoptions(varargin(1:end), {'threshold', 'bgcolor','template','MAP','label','isSave','plotrange','distfile','alpha'});

borderalighment={'left','right'};

% loop over hemisphere
[axespositions, hmap, hscale] = setupFlatMap(plotrange,1,1,bgcolor); % make panels
close(hscale);
figure(hmap);
for h=1:2
    % show result on flatmap
    coord = coords{h};
    topo = topos{h};
    shape = shapes{h};
    
    S=caret_load(shape);
    depth = S.data(:,2); % cortical depth information
    
    F=caret_load(coord);
    xlims = [min(F.data(:,1)), max(F.data(:,1))]; % plot range (x)
    ylims = [min(F.data(:,2)), max(F.data(:,2))]; % plot range (y)
    if ~iscell(MAP)
        MAPs = {MAP,MAP,MAP};
    else
        MAPs=MAP;
    end
    
    %subplot(1,2,h);
    ax(h) = axes('position',axespositions{h});
    if ~isempty(distfile)
        M=caret_load(distfile{h});
        % underlay
        caret_plotflatmap_rgb('coord',coord,'topo',topo,...
            'underlay', depth, 'data', M.data(:,1), 'dscale', [0.03,10],'threshold',0.03,...
            'xlims', xlims, 'ylims', ylims,'dmap',hot,'alpha',alpha); hold on;
        
        % overlay
        for col=1:size(data,2)
            idx = ~isnan(data(:,col,h))&data(:,col,h)>threshold(1);
            [~,p(col)]=caret_plotflatmap_rgb('coord',coord,'topo',topo,...
                'underlay', depth, 'data', data(:,col,h), 'dscale', threshold,'threshold',threshold(1),...
                'xlims', xlims, 'ylims', ylims,'dmap',MAPs{col},'alpha',alpha,'idx',idx);
            set(p(col),'facealpha',alpha);
        end
    else
        for col=1:size(data,2)
            idx = ~isnan(data(:,col,h))&data(:,col,h)>threshold(1);
            [~,p(col)]=caret_plotflatmap_rgb('coord',coord,'topo',topo,...
                'underlay', depth, 'data', data(:,col,h), 'dscale', threshold,'threshold',threshold(1),...
                'xlims', xlims, 'ylims', ylims,'dmap',MAPs{col},'alpha',alpha,'idx',idx);
            set(p(col),'facealpha',alpha);
        end
    end
    axis equal; axis tight; axis off; hold on;
    set(gca, 'xlim',plotrange{h}{1},'ylim',plotrange{h}{2});
    
    % draw sulci onto flatmap
    B = caret_load(border{h});
    for b=1:length(B.Border)
        %borderidx = cat(1,B.Border.vertex);
        borderidx = cat(1,B.Border(b).vertex);
        borderX = F.data(borderidx(:,1),1);
        borderY = F.data(borderidx(:,1),2);
        plot(borderX,borderY,'w.', 'markersize', bordersize);
        text(borderX(end),borderY(end), B.Border(b).name,...
            'horizontalalignment',borderalighment{h},...
            'verticalalignment','top',...
            'fontsize',fontsize,'color','w'); % border name
    end
end
% draw colorbar
% figure(hscale); colormap(MAP);
% ax(end+1)=axes('position', axespositions{end});
% colorscale = linspace(threshold(1),threshold(2),100);
% imagesc(colorscale); caxis(threshold);
% set(gca,'xcolor',[1,1,1]-bgcolor,'ycolor',[1,1,1]-bgcolor,'fontsize',fontsize,'tickdir','out');
% set(gca,'ytick',[],'xtick',[1,100],'xticklabel',{num2str(threshold(1)), num2str(threshold(2))});
% title(label,'color',[1,1,1]-bgcolor,'fontsize',fontsize);

% save figures
% fname_map = sprintf('Map.%s.%s.png',template,label);
% fname_scale = sprintf('Colorbar.%s.%s.png',template,label);
% mySaveFig(hmap, fullfile(figDir, fname_map), isSave, '-dpng','-r400');
% mySaveFig(hscale, fullfile(figDir, fname_scale), isSave, '-dpng','-r400');

varargout = {};
end
function varargout = setupFlatMap(plotrange,makeFig,printScale,bgcolor) % setup for flat cortical surface
% Make figure
if nargin==0
    plotrange = {{[-79.8252,104.7103],[-64.1677,81.5173-10]},...
        {[-104.7103,79.8252],[-64.1677,81.5173-10]}};
    makeFig=1;
    printScale=1;
    bgcolor='k';
end

aspectratio = diff(plotrange{1}{2}) / diff(plotrange{1}{1}) /2;
width = 40;
heightMap= width*aspectratio*1.05; % 20
heightScale=20-heightMap;

if makeFig==1
    figmap = figure('units','centimeters',...
        'papersize',[width heightMap]/2,...
        'position',[5 5 width heightMap],...
        'paperposition',[0 0 width heightMap]/2,...
        'paperpositionmode','manual',...
        'renderer','opengl'); % [width heightMap]
    figmap.Color = bgcolor;
    figmap.InvertHardcopy = 'off';
    if printScale==1
        figscale = figure('units','centimeters',...
            'papersize',[width heightScale]/2,...
            'position',[5 5 width heightScale],...
            'paperposition',[0 0 width heightScale]/2,...
            'paperpositionmode','manual',...
            'renderer','opengl'); % [width heightMap]
        figscale.Color = bgcolor;
        figscale.InvertHardcopy = 'off';
    end
end

topmargin = 0.01; %0.2
colorbarsize = [0.2 0.2]; % width height
centermargin = 0.005;
axespositions = {[0.0 topmargin 0.5-centermargin 1.0-topmargin],...
    [0.5+centermargin topmargin 0.5-centermargin 1.0-topmargin], ...
    [0.5-colorbarsize(1)/2 0.5-colorbarsize(2)/2 colorbarsize(1), colorbarsize(2)]};
% third one is for color scale

if ~exist('figmap','var');
    figmap = [];
end;
if ~exist('figscale','var');
    figscale = [];
end;
varargout = {axespositions,figmap,figscale};
end
function varargout = assign_patch2node(patchdata, hemis, patch, p2n, varargin) % assign surface nodes to patches
nodedata = NaN(163842,2);
smooth = 0;
F=[];
deleteMetric=1;
smoothoptions = {'algorithm','AN','iterations',20,'strength',0.5}; % 15
vararginoptions(varargin, {'smooth','F','smoothoptions'});

for h=1:2
    pa=patch(hemis==h);
    for p=pa'
        currNode=p2n{h}{p};
        nodedata(currNode,h) = patchdata(hemis==h&patch==p);
    end
    if (smooth==1)
        M = caret_load(F.meandist{h});
        M.data = nodedata(:,h);
        M.data(isnan(M.data)) = 0; % force nans to zero
        M.num_cols = 1;
        M.column_name = {'tmp'};
        metricName = fullfile(F.baseDir, 'tmp.metric');
        caret_savemetric(metricName, M);
        S = caret_smooth(metricName,...
            'coord', F.coord{h},...
            'topo', F.topo{h},...
            smoothoptions{:}); % smooth .metric file
        sM = caret_load(S{1}); % load smoothed data
        nodedata(:,h) = sM.data; % swap node data with smoothed data
        % delete .metric file
        if deleteMetric
            delete(metricName); delete(S{1});
        end
        %nodedata=ceil(nodedata);
    end
end
varargout  = {nodedata};

end
function varargout = Gtest(counts) % perform likelihood-ratio G-test
%counts = varargin{1}; % [00 10 01 11] ;  Counts for no, A present, B present, both present
N=sum(counts);
pA = (counts(2)+counts(4)) / N;  % Overall probability of A present
pB = (counts(3)+counts(4)) / N;  % Overall probability of B present
ExpCounts = [N*(1-pA)*(1-pB) N*pA*(1-pB) N*(1-pA)*pB N*pA*pB];
Chi2 =sum((counts -ExpCounts).^2./ExpCounts);
G    = 2*nansum(counts.*log(counts./ExpCounts));
criticalValue = chi2inv(0.95,1);
pval = 1-chi2cdf(G,1);
fprintf('Obs:');
fprintf('%2.2f   ',counts);fprintf('\nExp:');
fprintf('%2.2f   ',ExpCounts);fprintf('\n');
fprintf('Chi2:%2.2f\n',Chi2);
fprintf('G:   %2.2f\n',G);
fprintf('Crit:%2.2f\n',criticalValue);
fprintf('p-val:%f\n',pval);
varargout = {G, pval, ExpCounts, counts};
end
