function varargout=sh1_behavior(what,varargin)
%% function varargout=sh1_behavior(what,varargin)
% Do all behavioral data analyses and reproduce figures 
% for Yokoi&Diedrichsen(2019) paper.
%
%
%
%
%
%

%% Settings
% path to data/result
baseDir         = '/Volumes/G_Thunderbolt/Yokoi_Research/data/SequenceLearning/sh1/gittoshare'; % external HDD 2
analyzeDir 	= fullfile(baseDir, 'data');
figDir      = fullfile(baseDir, 'figure');

% subject info
subject_name = {{'sh1_t01','sh1_t02','sh1_s03','sh1_s04','sh1_s05','sh1_s07','sh1_s22'},... 
    {'sh1_s09','sh1_s10','sh1_s14','sh1_s15','sh1_s17','sh1_s19','sh1_s23','sh1_s25'}};
subject_excluded = {'sh1_s05','sh1_s09','sh1_s14'}; 
% s05, s09, s14 were due to poor performance in the scanner
for s=1:numel(subject_excluded);
    subject_excluded_post{s} = [subject_excluded{s},'_post'];
end
for s=1:numel(subject_excluded);
    subject_excluded_scan{s} = [subject_excluded{s},'_scan'];
end

%% Main operations
switch (what)
    case 'Figure_all' % draw all panels and save stat table
        ipitype = 'ipi_maxtime'; % 'ipi_onset'
        bgcolor = 'w'; % 'k'
        vararginoptions(varargin,{'ipitype','bgcolor'});
        
        % Create Figures
        % paper size: 8.5 x 11 inch (1 inch=2.54 cm)
        % 1 column, 85 mm; 1.5 column, 114 mm; and 2 column, 174 mm (the full width of the page)
        fig1 = figure('units','centimeters',...
            'position',[1,1,11.4,11.4/sqrt(2)*1.2],...
            'paperposition',[1,1,11.4,11.4/sqrt(2)*1.2],...
            'papersize',[8.5,11]*2.54,'name','Figure1'); % 1.5 column size        
        fig2 = figure('units','centimeters',...
            'position',[2,2,17.4,17.4/sqrt(2)*2],...
            'paperposition',[2,2,17.4,17.4/sqrt(2)*2],...
            'papersize',[8.5,11]*2.54,'name','Figure2'); % 2 column size
        
        % training data (within- vs. between-chunk IPIs for training sessions)
        % Fig. 1D
        w1 = 11.4*2/3; % cm
        h1 = w1/sqrt(2);
        x1 = 11.4/5;
        y1 = h1/2;
        figure(fig1);
        [stats1, ax1]=sh1_behavior('panel_training', 'ipitype', ipitype, 'fig', [x1,y1,w1,h1],'bgcolor',bgcolor);
        
        % followup data 1 (Inter-subject correlation of IPI patterns for physically identical sequences)
        % Fig. 2B,C
        figure(fig2);
        w2 = 17.4*0.65*0.9; % cm
        h2 = w2/sqrt(2)/2;
        x2 = 17.4/7;
        y2 = fig2.Position(4) - 1.5*h2;
        stats2=sh1_behavior('panel_followup_1', 'ipitype', ipitype,...
            'fig', {[x2,y2,w2,h2],[2*x2+w2,y2,w2/6,h2]},'bgcolor',bgcolor);
        
        % followup data 2 (within- vs. non-within-chunk IPIs for followup session)
        % Fig. 2E
        figure(fig2);
        %subplot(2,15,[12:15])
        w3 = w2*0.5; % cm
        h3 = w3/sqrt(2)*1.3;
        x3 = 17.4/7;
        y3 = y2-1.5*h3;
        
        stats3=sh1_behavior('panel_followup_2', 'ipitype', ipitype, 'fig', [x3,y3,w3,h3],'bgcolor',bgcolor);
        
        % followup data 3 (model-averaged regression coefficients with protected exceedance probability)
        % Fig. 2F
        figure(fig2);
        w4 = w3*0.5; % cm
        h4 = w4/sqrt(2)*2.5;
        x4 = x3+w3*1.7;
        y4 = y3;
        stats4=sh1_behavior('panel_followup_3', 'ipitype', ipitype,...
            'fig', {[x4,y4,w4,h4],[x4+1.7*w4,y4,w4/6,h4]},'bgcolor',bgcolor);
        
        % show stats
        Stats=addstruct(stats1,stats2);
        Stats=addstruct(Stats,stats3);
        struct2table(Stats)
        struct2table(stats4)
        
        % save figure and stats
        print(fig1,'-dpsc2','-r400',fullfile(figDir,sprintf('sh1_Fig1_%s.eps',ipitype)));
        print(fig2,'-dpsc2','-r400',fullfile(figDir,sprintf('sh1_Fig2_%s.eps',ipitype)));
        dsave(fullfile(figDir,sprintf('sh1_behav_stats_%s.txt',ipitype)),Stats);
        dsave(fullfile(figDir,sprintf('sh1_behav_bms_stats_%s.txt',ipitype)),stats4);
        
        varargout = {Stats,stats4};
    case 'panel_training' % plot within- and between-chunk IPIs
        ipitype = 'ipi_maxtime';  % 'ipi_onset'
        markersize  = 4;
        linewidth   = 1.25;
        errorwidth  = 1.25;
        fig = 1;
        bgcolor = 'w';
        markercolors = {'b','r',[0.8,0.8,0.8]};
        vararginoptions(varargin, {'ipitype','markersize','linewidth','errorwidth','fig','bgcolor','markercolors'});
                
        % ===================================== %
        % load data for training and calc summary
        % ===================================== %
        fname_training = fullfile(analyzeDir, 'All_training.mat');
        T = load(fname_training);
        T = getrow(T, ~ismember(T.subjname, subject_excluded));
        % some stuff to do
        T.ipi = T.(ipitype);
        T.ipi_norm = T.([ipitype, '_norm']);
        T.error = nansum(T.errorVec,2);
        for t=1:length(T.SN) % .ipi_within and .ipi_between (un-normalized)
            T.ipi_within(t,1) = nanmean(T.ipi(t, T.withinchunk(t,:)==1));
            if ~isempty(T.ipi(t, T.withinchunk(t,:)~=1));
                T.ipi_between(t,1) = nanmean(T.ipi(t, T.withinchunk(t,:)~=1));
            else
                T.ipi_between(t,1) = NaN;
            end
        end;
        % tapply
        T.errorRate=T.error>0;
        subset = T.error==0;
        Tplot = tapply(T, {'BN','Day','SN','group','subjname'},... % ,'seqType'
            {T.ipi_within, 'nanmean', 'subset', subset,'name', 'ipi_within'},...
            {T.ipi_between, 'nanmean', 'subset', subset, 'name' 'ipi_between'},...
            {T.error, 'nanmean', 'name', 'seqError'},...
            {T.errorRate, 'nanmean', 'name', 'errorRate'});
        
        % ===================================== %
        % load data for scanning session and calc summary
        % ===================================== %
        fname_scan = fullfile(analyzeDir, 'All_scan.mat');
        Tscan = load(fname_scan);
        Tscan = getrow(Tscan, ~ismember(Tscan.subjname, subject_excluded_scan));
        % some stuff to do
        Tscan.ipi = Tscan.(ipitype);
        Tscan.ipi_norm = Tscan.([ipitype, '_norm']);
        Tscan.withinchunk = logical(Tscan.withinchunk);
        Tscan.error = nansum(Tscan.errorVec,2);
        
        for t=1:length(Tscan.SN) % .ipi_within and .ipi_between (un-normalized)
            Tscan.ipi_within(t,1) = nanmean(Tscan.ipi(t, Tscan.withinchunk(t,:)));
            Tscan.ipi_between(t,1) = nanmean(Tscan.ipi(t, ~Tscan.withinchunk(t,:)));
        end
        Tscan.errorRate = Tscan.error>0;
        subset = Tscan.error==0;
        Tplot_scan = tapply(Tscan,{'SN','Day','group','subjname'},...
            {'ipi_within','nanmean','subset', subset, 'name','ipi_within'},...
            {'ipi_between','nanmean','subset', subset, 'name','ipi_between'},...
            {'error','nanmean','name','seqError'},...
            {'errorRate', 'nanmean', 'name', 'errorRate'});
        Tplot_scan.BN = repmat(55, size(Tplot_scan.SN));
        
        % --------------------------------------------------------------- %
        % Plot block-mean of between- and within-chunk intervals
        % --------------------------------------------------------------- %
        figName = 'Between- and within-chunk intervals (w.o. cue, block-averaged)';
        if fig==1
            % plot in single axis
            figure('name',sprintf('%s (N=%d)',figName, numel(unique(T.SN))),...
                'color','w',...
                'Units','centimeters','Position',2*[5 5 8.9*1.2 8.9/sqrt(2)*1.2]);
            ax{1} = subplot(2,15,[1:7]); % ipi
        elseif length(fig)>1
            ax{1} = axes('units','centimeters', 'position', fig);
            set(ax{1},'ActivePositionProperty','position');
        end
        set(ax{1}, 'yaxislocation', 'left', 'color', 'none', 'hittest', 'off');
        ax{2} = axes('handlevisibility',get(ax{1},'handlevisibility'),...
            'units',get(ax{1},'units'),...
            'position',get(ax{1},'position'),...
            'parent',get(ax{1},'parent')); % error
        set(ax{2}, 'yaxislocation', 'right', 'color', 'none', 'hittest', 'off');
        
        % deal with background/axis colors
        set(gcf,'Color',bgcolor,'InvertHardcopy','off');            
        switch bgcolor
            case 'w'
                axcolor = [0,0,0];
            case 'k'
                axcolor = [1,1,1];
            otherwise
                axcolor = [0,0,0];
        end
        set(ax{1},'xcolor',axcolor,'ycolor',axcolor);
        set(ax{2},'xcolor',axcolor,'ycolor',axcolor);
        
        T = addstruct(Tplot,Tplot_scan);
        
        % (t01,t02 data is also excluded from visualization due to slightly different block assignment)
        subject_excluded{end+1} = 'sh1_t01';
        subject_excluded{end+1} = 'sh1_t02';
        subject_excluded{end+1} = 'sh1_t01_scan';
        subject_excluded{end+1} = 'sh1_t02_scan';
        subset = ~ismember(T.subjname, subject_excluded);
        % get block number
        figure(99);
        [x,blocks] = lineplot([T.Day, T.BN], T.BN,'subset',subset,'gap',[1.5 0.5]);
        close(99);
        
        % 1. plot ipi
        axes(ax{1});
        colors = markercolors(1:2); %{'b','r'};
        leg = {'Within','Between'};
        [xtick,y,e] = lineplot([T.Day T.BN], [T.ipi_within, T.ipi_between],...
            'subset', subset,...
            'plotfcn','nanmean','leg', leg, 'leglocation','northeast',...
            'style_thickline','markertype',{'o','o'},'gap',[1.5 0.5],...
            'markersize',markersize,'linewidth',linewidth,'markerfill',colors,...
            'markercolor',colors,'linecolor',colors,...
            'errorwidth',errorwidth,'errorcap',0.0001,'errorcolor',colors); hold on
        ylabel('Inter press interval (msec)','color',axcolor);
        xlabel('Block no.','color',axcolor);
        
        % 2. plot seqError
        axes(ax{2});
        colors = markercolors(3);
        leg = {'Error'};
        [x,y,e]=lineplot([T.Day T.BN], [T.seqError],...
            'subset', subset,'split',true(size(subset)),...
            'plotfcn','nanmean',...
            'style_thickline','markertype',{'s'},'gap',[1.5 0.5],...
            'markersize',markersize,'linewidth',linewidth,'markerfill',colors,...
            'markercolor',colors,'linecolor',colors,...
            'errorwidth',errorwidth,'errorcap',0.0001,'errorcolor',colors,...
            'leg', leg, 'leglocation','east'); hold on
        ylabel('#Incorrect presses','color',axcolor);
        xlabel('');
        
        % Axis properties
        ticklocs=[1,10,18,22,30,32,33];
        Xticks = xtick(ticklocs);
        for i=1:length(ticklocs)
            Xticklabels{i} = sprintf('%d', blocks(ticklocs(i)));
        end
        Xticklabels{end}='S';
        %Xticklabels = {'2','15','17','28','30','35','37','47','49','50',''};
        Xlim = [Xticks(1)-1,Xticks(end)+1];
        % 1
        set(ax{1},'XTick',Xticks,'Xlim',Xlim,'XTicklabel',Xticklabels,'ylim',[0 2000],...
            'ytick',[0:400:2000],'tickdir','out','ticklength',[0.02,0.02],'linewidth',1);
        % 2
        set(ax{2},'XTick',[],'Xlim',Xlim,'ylim',[0, 5],'ytick',[0:5],...
            'tickdir','out','ticklength',[0.02,0.02],'linewidth',1);
        axes(ax{1}); % make ax1 top
        % align legend
        h=get(get(gca,'parent'),'children');
        h(3).Position(1)=h(1).Position(1);
        h(3).Position(2)=h(1).Position(2)-h(1).Position(4)/2;
        
        % ============================== %
        % Stats: compare ipi on day 5
        % ============================== %
        D = tapply(T,{'SN'},{'ipi_within','nanmean','subset',T.Day==5,'name','within'},...
            {'ipi_between','nanmean','subset',T.Day==5,'name','between'});
        fprintf('---- ttest within vs between on day 5 ----\n');
        [t,p]=ttest_mc(D.within,D.between,2,'paired');
        
        Stats.comparison{1} = 'within- vs. between-chunk ipi on day 5';
        Stats.mean1 = nanmean(D.within);
        Stats.SD1 = nanstd(D.within);
        Stats.mean2 = nanmean(D.between);
        Stats.SD2 = nanstd(D.between);
        Stats.df = numel(D.SN)-1;
        Stats.tval = t;
        Stats.pval = p;
        struct2table(Stats)
        
        varargout = {Stats, ax};
    case 'panel_followup_1' % plot across- vs within-group correlation of ipi (normalized)
        ipitype = 'ipi_maxtime';  % 'ipi_onset'
        titlename = {'Different chunkings on common sequences','(followup session)'};
        corrRange = [0 0.7];
        fig=1;
        bgcolor = 'w'; % 'k'
        vararginoptions(varargin,{'ipitype','titlename','corrRange','fig','bgcolor'});
        
        % ===================================== %
        % load data for followup session and calc summary
        % ===================================== %
        fname_followup = fullfile(analyzeDir, 'All_posttest.mat');
        T = load(fname_followup);
        T = getrow(T, ~ismember(T.subjname, subject_excluded_post));
        % some stuff to do
        T.ipi = T.(ipitype);
        T.ipi_norm = T.([ipitype, '_norm']);
        T.error = nansum(T.errorVec,2);
        
        % Take only trained&identical sequences and go
        seq_identity = [1,2,4,6,8]; % these are phisically the same sequences across set A and set B sbj
        D = getrow(T, T.seqCat==1&T.error==0&ismember(T.seqType, seq_identity));
        D.ipi_norm = bsxfun(@rdivide,D.ipi,sum(D.ipi,2)); % normalized ipi
        D = tapply(D, {'group','SN','seqType'},{'ipi_norm','nanmean(x,1)','name','ipi'});
        
        % ===================================== %
        % Calc within- and across-group correlation of ipi pattern
        % ===================================== %
        % some preparation
        Nsubj   = length(unique(D.SN));
        NsetA   = length(unique(D.SN(D.group==1)));
        NsetB   = Nsubj-NsetA;
        Idxw    = logical(blockdiag(true(NsetA),true(NsetB)));
        count = 0;
        for seq = seq_identity
            count=count+1;
            Ds = getrow(D,D.seqType==seq);
            
            % get correlation and z-transform
            tCorrs(:,:) = fisherz(corr(Ds.ipi','type','Kendall'));
            pCorrs(:,:) = fisherz(corr(Ds.ipi','type','Pearson'));
            
            % kill diagonal part
            tCorrs(logical(eye(Nsubj))) = NaN;
            pCorrs(logical(eye(Nsubj))) = NaN;
            
            % split into within and between parts
            tmpc = tCorrs;
            tmpc(~Idxw) = NaN;
            tWcorr(:,:,count) = tmpc;
            tmpc = tCorrs;
            tmpc(Idxw) = NaN;
            tBcorr(:,:,count) = tmpc;
            tmpc = pCorrs;
            tmpc(~Idxw) = NaN;
            pWcorr(:,:,count) = tmpc;
            tmpc = pCorrs;
            tmpc(Idxw) = NaN;
            pBcorr(:,:,count) = tmpc;
        end
        tWcorr = nanmean(tWcorr,3); % average over sequences
        tBcorr = nanmean(tBcorr,3);
        pWcorr = nanmean(pWcorr,3);
        pBcorr = nanmean(pBcorr,3);
        
        % ===================================== %
        % Figure (ipi patterns and within- vs -between correlation)
        % ===================================== %
        if length(fig)>1
            % create axes
            ax{1} = axes('units','centimeters','position',fig{1});
            ax{2} = axes('units','centimeters','position',fig{2});
        elseif fig==1
            figure('name','setA vs setB','color','w',...
                'Units','centimeters','Position',2*[5 5 8.9*1.5 8.9/sqrt(2)*1.0]);
            ax{1} = subplot(2,10,[1:8]); ax{1}.Position(2)=ax{1}.Position(2)/2;
            ax{2} = subplot(2,10,9.8); ax{2}.Position(2)=ax{2}.Position(2)/2;
        end
        % deal with background/axis colors
        set(gcf,'Color',bgcolor,'InvertHardcopy','off');
        switch bgcolor
            case 'w'
                axcolor = [0,0,0];
            case 'k'
                axcolor = [1,1,1];
            otherwise
                axcolor = [0,0,0];
        end
        set(ax{1},'xcolor',axcolor,'ycolor',axcolor,'color','none');
        set(ax{2},'xcolor',axcolor,'ycolor',axcolor,'color','none');
        
        % average IPI pattern for two groups
        data = D.ipi;
        x = repmat(D.group,1,10);
        data = data';
        x = x';
        data = data(:);
        seq = repmat(D.seqType,1,10);
        seq = seq';
        seq = seq(:);
        press = repmat([1:10]',1,length(D.SN));
        press = press(:);
        
        axes(ax{1}); % ipi patterns
        xcat = lineplot([seq,press],data,'split',x(:),'linecolor',{'r','b'},...
            'markertype','non','errorcolor',{'r','b'},'errorcap',0.0001,...
            'linewidth',1.5,'errorwidth',1.5,'leg',{'G1','G2'},...
            'leglocation','northeast');
        
        ylabel('Normalized IPI','fontsize',15,'color',axcolor);
        xlabel('Interval no.','fontsize',15,'color',axcolor);
        title(titlename,'fontsize',15,'color',axcolor);
        set(gca,'fontsize',15,'tickdir','out','ticklength',[0.02 0.02],...
            'xtick',xcat([1 10 11 20 21 30 31 40 41 50]),'linewidth',1,...
            'xticklabel',{'1','10','','','','','','','',''},...
            'ylim',[0.05 0.2],'ytick',[0.05:0.05:0.2]);
        
        % within- vs between-group correlation
        r.Kendall   = [nanmean(tWcorr,2); nanmean(tBcorr,2)];
        r.Pearson   = [nanmean(pWcorr,2); nanmean(pBcorr,2)];
        r.type      = [repmat(1,length(tWcorr),1);repmat(2,length(tBcorr),1)];
        
        axes(ax{2}); % within vs. between correlation
        barplot(r.type,fisherinv(r.Pearson),...
            'facecolor', {[0.9 0.9 0.9],[0.1 0.1 0.1]},...
            'edgecolor', {axcolor}, 'linewidth', 1,...
            'errorwidth', 1, 'errorcolor', axcolor,...
            'capwidth',0.0001);
        
        drawline(0,'dir','horz','linewidth',1,'color',axcolor);
        ylabel('Correlation coeff.','fontsize',15,'color',axcolor);
        %title({'Within vs Between set', 'correlation of ipi (posttest)'},'fontsize',15);
        set(gca,'fontsize',15,'tickdir','out','ticklength',[0.02 0.02]*4,...
            'linewidth',1,'xticklabel',{'W','B'},...
            'ylim',corrRange,'ytick',[0 0.2 0.4 0.6 0.8]);
        set(ax{2},'xcolor',axcolor,'ycolor',axcolor,'color','none');
        
        % ===================================== %
        % Stats within-set and between-set across-subject correlation
        % ===================================== %
        % Kendall's tau
        data1 = nanmean(tWcorr,2);
        data2 = nanmean(tBcorr,2);
        [t,p]=ttest_mc(data1, data2, 2,'paired');
        
        Stats.comparison{1,1} = 'Paired t-test for Kendall''s \tau (z-transformed)';
        Stats.mean1(1,1) = nanmean(data1);
        Stats.SD1(1,1) = nanstd(data1);
        Stats.mean2(1,1) = nanmean(data2);
        Stats.SD2(1,1) = nanstd(data2);
        Stats.df(1,1) = length(data1)-1;
        Stats.tval(1,1) = t;
        Stats.pval(1,1) = p;
        
        % Pearson's r
        data1 = nanmean(pWcorr,2);
        data2 = nanmean(pBcorr,2);
        [t,p]=ttest_mc(data1, data2, 2,'paired');
        
        Stats.comparison{2,1} = 'Paired t-test for Pearson''s r (z-transformed)';
        Stats.mean1(2,1) = nanmean(data1);
        Stats.SD1(2,1) = nanstd(data1);
        Stats.mean2(2,1) = nanmean(data2);
        Stats.SD2(2,1) = nanstd(data2);
        Stats.df(2,1) = length(data1)-1;
        Stats.tval(2,1) = t;
        Stats.pval(2,1) = p;
        
        varargout = {Stats,ax};
    case 'panel_followup_2' % plot within-chunk vs non-within-chunk ipis for followup session
        ipitype = 'ipi_maxtime';  % 'ipi_onset'
        fig=1;
        bgcolor = 'w'; 
        facecolor = {'b','r'};
        vararginoptions(varargin, {'ipitype','fig', 'bgcolor','facecolor'});
        
        % ===================================== %
        % load data for followup session and calc summary
        % ===================================== %
        fname_followup = fullfile(analyzeDir, 'All_posttest.mat');
        T = load(fname_followup);
        T = getrow(T, ~ismember(T.subjname, subject_excluded_post));
        % some stuff to do
        T.ipi = T.(ipitype);
        T.ipi_norm = T.([ipitype, '_norm']);
        T.error = nansum(T.errorVec,2);
        
        for t=1:length(T.SN) % .ipi_within and .ipi_between (un-normalized)
            within = T.withinchunk(t,:)==1;
            nowithin = T.withinchunk(t,:)~=1;
            if ~isempty(within)
                T.ipi_within(t,1) = nanmean(T.ipi(t, within));
            else
                T.ipi_within(t,1) = NaN;
            end
            if isempty(nowithin)
                T.ipi_between(t,1) = NaN;
            else
                T.ipi_between(t,1) = nanmean(T.ipi(t, nowithin));
            end
        end
        T.errorRate = T.error>0;
        
        % Merge Chunk and Super-chunk categories
        T.seqCat(T.seqCat==3) = 2;
        subset = T.error==0;
        T = tapply(T,{'group','SN','seqCat'},...
            {'ipi_within','nanmean','subset',subset,'name','ipi_wChunk'},...
            {'ipi_between','nanmean','subset',subset,'name','ipi_oChunk'},...
            {'error','nanmean','name','seqError'});
        
        % =============================== %
        % Figure
        % =============================== %
        if length(fig)>1
            ax{1} = axes('units','centimeters','position',fig);
        elseif fig==1
            figure('name','Within- vs. non-within-chunk IPIs (followup session)',...
                'color','w',...
                'Units','centimeters','Position',2*[5 5 8.9/2 8.9/sqrt(2)*1.3]);            
            ax{1}=subplot(2,3,[2,3]);
        end
        % deal with background/axis colors
        set(gcf,'Color',bgcolor,'InvertHardcopy','off');
        switch bgcolor
            case 'w'
                axcolor = [0,0,0];
            case 'k'
                axcolor = [1,1,1];
            otherwise
                axcolor = [0,0,0];
        end
        set(ax{1},'xcolor',axcolor,'ycolor',axcolor,'color','none');
        
        seqCat = T.seqCat;
        seqCat(seqCat==2)=6; % re-ordering
        seqCat(seqCat==3)=2;
        seqCat(seqCat==6)=3;
        
        xtick = barplot(seqCat,[T.ipi_wChunk,T.ipi_oChunk],...
            'barwidth',0.8,'linewidth',1.0,...
            'facecolor',facecolor,...
            'edgecolor',axcolor,'gapwidth',[1 0.3],...
            'errorcolor', axcolor,...
            'subset',T.seqError>=0,'capwidth',0.0001); hold on;
        
        set(gca,'xlim',[xtick(1)-2,xtick(end)+2],'ylim',[150 500]);
        drawline(ax{1}.YLim(1),'dir','horz','color',axcolor,'linewidth',ax{1}.LineWidth);
        set(gca,'tickdir','out','ticklength',[0.03 0.03],'linewidth',1);
        set(gca,'xtick',(xtick(1:2:end)+xtick(2:2:end))/2,'ytick',[200:100:500]);
        set(gca,'xticklabel',{'Trained','Chunk','Chunk+Random','Random'},...
            'xticklabelrotation',45);
        ylabel('Inter pres interval (msec)','color',axcolor)%,'fontsize',13);
        set(ax{1},'color','none','xcolor',axcolor,'ycolor',axcolor);
        
        % =============================== %
        % Stats
        % =============================== %
        comparisons = {'within vs. non-within (trained)',...
            'within vs. non-within (chunk reordered)',...
            'within vs. non-within (chunk+new)',...
            'non-within (trained) vs. non-within (others)'};
        % 1. t-test: within vs outside chunk intervals for each
        % sequence types
        fprintf('----- t-test: Within-chunk interval vs outside-of-chunk interval -----\n\n');
        c=0;
        for seq=[1,2,4]
            c=c+1;
            subset = T.seqCat==seq;
            [t,p]=ttest_mc(T.ipi_wChunk(subset),T.ipi_oChunk(subset),2,'paired');
            
            Stats.comparison{c,1} = comparisons{c};
            Stats.mean1(c,1) = nanmean(T.ipi_wChunk(subset));
            Stats.SD1(c,1) = nanstd(T.ipi_wChunk(subset));
            Stats.mean2(c,1) = nanmean(T.ipi_oChunk(subset));
            Stats.SD2(c,1) = nanstd(T.ipi_oChunk(subset));
            Stats.df(c,1) = sum(subset)-1;
            Stats.tval(c,1) = t;
            Stats.pval(c,1) = p;
        end
        
        % 2. t-test: non-chunk interval of trained vs that of other
        D = tapply(T,{'SN'},...
            {'ipi_oChunk','nanmean','subset',T.seqCat==1,'name','trained'},...
            {'ipi_oChunk','nanmean','subset',ismember(T.seqCat,[2,3,4]),'name','novel'}...
            );
        fprintf('---- ttest: nonchunk(trained) vs nonchunk(novel) ----\n');
        [t,p]=ttest_mc(D.trained,D.novel,2,'paired');
        
        Stats.comparison{end+1} = comparisons{end};
        Stats.mean1(end+1) = nanmean(D.trained);
        Stats.SD1(end+1) = nanstd(D.trained);
        Stats.mean2(end+1) = nanmean(D.novel);
        Stats.SD2(end+1) = nanstd(D.novel);
        Stats.df(end+1) = length(D.novel)-1;
        Stats.tval(end+1) = t;
        Stats.pval(end+1) = p;
        if nargout==0;
            struct2table(Stats)
        end
        varargout = {Stats, ax};
    case 'panel_followup_3' % plot results for ipi regression with k-fold cross-validation
        ipitype = 'ipi_maxtime';  % 'ipi_onset'; % 
        CVfolds = {[1],[2],[3],[4],[5],[6],[7],[8],[9],[10],[11],[12],[13],[14],[15],[16]};%{[1,3,5],[2,4]};
        fig = 1;
        models = [1:3];
        bgcolor = 'w';
        vararginoptions(varargin, {'ipitype','CVfolds', 'fig','bgcolor'});
        
        % ===================================== %
        % load data for followup session and calc summary
        % ===================================== %
        fname_followup = fullfile(analyzeDir, 'All_posttest.mat');
        T = load(fname_followup);
        T = getrow(T, ~ismember(T.subjname, subject_excluded_post));
        % some stuff to do
        T.ipi = T.(ipitype);
        T.ipi_norm = T.([ipitype, '_norm']);
        
        % normalize bio (choose onset or peakf)
        Biomech = load(fullfile(analyzeDir, 'BiomechanicalEffect.mat'));
        switch ipitype
            case 'ipi_maxtime'
                biomech = Biomech.peakf;
            case 'ipi_onset'
                biomech = Biomech.onset;
        end
        
        % transition frequency during training
        Tran=linregress_get_tranfreq(analyzeDir, subject_excluded);
        
        % loop over individuals
        F = [];
        subjects = unique(T.SN);
        for s=subjects';
            fprintf('s=%d\n', s);
            D = getrow(T, T.SN==s);
            tran = getrow(Tran,Tran.SN==s);
            
            % get X
            [X, Xname] = linregress_getX(D, biomech, tran.freq);
            freq = X(:,3);
            meanfreq(s) = mean(freq(freq~=0)); % this is for later visualization
            
            % get Y
            [Y, isOK, execution] = linregress_getY(D);
            
            % define model family and do regression with CV
            [Result, compidx] = fitModelFamilyCrossval(X,Xname,Y,isOK,execution,CVfolds);
            Result.SN = repmat(s, size(Result.model,1),1);
            Result.group = repmat(D.group(1), size(Result.model,1),1);
            
            % do model averaging
            r=getrow(Result, Result.SN==s);
            
            % averaged parameter
            dlike = r.loglike-nanmax(r.loglike);
            w = exp(dlike)/nansum(exp(dlike));
            avrgbeta = bsxfun(@times, r.beta, w);
            fac.avrgbeta = nansum(avrgbeta,1);
            
            % logBFc
            [PP,logBF]=pcm_componentPosterior(r.loglike', r.compidx);
            fac.PP = PP;
            fac.logBF=logBF;
            fac.SN=s;
            fac.group=D.group(1);
            F=addstruct(F,fac);
        end
        % ===================================== %
        % Stats (bms on logBFc and ttest on avrgbeta)
        % ===================================== %
        for com=1:size(compidx,2)-1;
            beta = F.avrgbeta(:,com);
            [t,p] = ttest_mc(beta,0,1,'onesample');
            Stats.modelnum(com,1) = com;
            Stats.model{com,1} = Xname{com};
            Stats.Beta(com,1) = nanmean(beta);
            Stats.SD_beta(com,1) = nanstd(beta);
            Stats.df(com,1) = numel(beta)-1;
            Stats.tval(com,1) = t;
            Stats.pval(com,1) = p;
            
            y = F.logBF(:,com); % individual factor BF for a factor
            y = y(isfinite(y));
            [alpha(:,com),exp_r(:,com,1),xp(:,com,1),pxp,bor(:,com,1)] = spm_BMS ([y, zeros(size(y))]);
            Stats.logBFc(com,1) = nanmean(y);
            Stats.SD_logBFc(com,1) = nanstd(y);
            Stats.pxp(com,1) = pxp(1);
        end
        varargout = {Stats};
        Stats.Beta(3) = -Stats.Beta(3)*mean(meanfreq(meanfreq~=0)); % scale transition weight up to its mean effect
        Stats.SD_beta(3) = -Stats.SD_beta(3)*mean(meanfreq(meanfreq~=0));
        
        % ===================================== %
        % Figure
        % ===================================== %
        % set color for pxp
        minpxp=0.4;
        pxpcolor = autumn(100);
        pxpedge = linspace(minpxp,1,100);
        pxpbins = discretize(Stats.pxp,pxpedge);
        facecolor = {pxpcolor(pxpbins(1),:), pxpcolor(pxpbins(2),:), pxpcolor(pxpbins(3),:), pxpcolor(pxpbins(4),:)};
        if length(fig)>1
            ax{1} = axes('units','centimeters','position',fig{1});
            ax{2} = axes('units','centimeters','position',fig{2});
        elseif fig==1 % create window or axes
            %myFigure([52/2.5 17/2.2]);
            figure('name','Linear regression on followup IPI data (cross-validated)',...
                'color','w',...
                'Units','centimeters','Position',[5, 5, 52/2.5, 17]);
            ax{1} = subplot(2,18,[1:4]); % average regression weights with pxp
            ax{2} = subplot(2,18,7); % colorbar
        end
        % deal with background/axis colors
        set(gcf,'Color',bgcolor,'InvertHardcopy','off');
        switch bgcolor
            case 'w'
                axcolor = [0,0,0];
            case 'k'
                axcolor = [1,1,1];
            otherwise
                axcolor = [0,0,0];
        end
        
        axes(ax{1}); % main graph
        xpos=barplot([],Stats.Beta(models),'split',Stats.modelnum(models),...
            'errorfcn',Stats.SD_beta(models)./(sqrt(length(subjects))),...
            'capwidth',0.001,'errorwidth',1,'linewidth',1,...
            'facecolor', facecolor, 'gapwidth', [1.5,0.3,0.1],...
            'edgecolor',axcolor,'errorcolor',axcolor);
        
        ylabel({'Avrg. regression coeff.', '(msec)'},'color',axcolor);
        set(gca,'xtick', [xpos],'xcolor',axcolor,'ycolor',axcolor,'color','none');
        set(gca,'xticklabel',Stats.model(models),'ylim', [-10,160],'ytick',[0, 50, 100, 150]);
        set(gca, 'tickdir', 'out', 'ticklength', [0.03 0.01],'xticklabelrotation',45);
        title({''},'color',axcolor);
        %set(gca,'view',[90,90],'yaxislocation','left');
        drawline([0],'dir','horz','color',axcolor); hold on;
        set(ax{1},'xcolor',axcolor,'ycolor',axcolor,'color','none','linewidth',1);        
        
        axes(ax{2}); % color bar
        myColorbar('MAP','autumn','scaling',[minpxp,1],'dir','vert','ylabelname','PXP',...
            'bgcolor',bgcolor);
        set(ax{2},'xcolor',axcolor,'ycolor',axcolor,'color','none','linewidth',1);
        
        if nargout==0
            struct2table(Stats)
        end
        varargout = {Stats,ax};
    otherwise
        warning('Case not found.');
end

end
%% Local functions
function [Y, isOK, execution] = linregress_getY(D, maxfold) % prepare IPIs for regression
if nargin<2
    maxfold=0;
end
% =========================== %
% Organize ipi data into single vector (Y)
% =========================== %
Y = reshape(D.ipi',10*size(D.BN,1),1);
isOK = ~isnan(Y);

% Define maximally possible number of cross-validation folds
% by treating repetition x successive execusion as single cv fold
% (this yields similar result)
if (maxfold==1)
    uniseq=[D.seqType(1)];
    nrep=zeros(length(unique(D.seqType)));
    maxexe = max(D.exe);
    for t=2:length(D.exe)
        preseq= D.seqType(t-1);
        currseq= D.seqType(t);
        if preseq==currseq
        else
            if ismember(currseq, uniseq)
                %fprintf('trial=%d, seqtype %d appeared before.\n', t, currseq);
                nrep(currseq) = nrep(currseq)+1;
                D.exe(t:t+maxexe-1) = D.exe(t:t+maxexe-1)+maxexe*(nrep(currseq));
                %keyboard();
            else
                uniseq = unique([uniseq, D.seqType(t)]);
            end
        end
    end
end
execution = reshape(repmat(D.exe,1,size(D.ipi,2))', [],1);
end
function [X, Xname] =  linregress_getX(D,biomech,tranfreq) % prepare design matrix for regression
% =========================== %
% Make design matrix X
% =========================== %
% 1. Within-chunk
X_chunk = reshape(D.withinchunk',[],1);

% 2. Higher-order sequence (known chunk-transition)
% collect trained chunk-transitions
ctran=[];
for seq=1:8
    chunk=D.chunk(D.seqType==seq,:);
    idx=find(chunk(1,:)==0);
    for i=1:numel(idx)
        tmp(1,1) = chunk(1,idx(i)-1);
        tmp(1,2) = chunk(1,idx(i)+1);
        ctran=[ctran; tmp];
    end
end
% handle trained sequences
CT = zeros(size(D.withinchunk));
idx=bsxfun(@and, D.seqCat==1, D.withinchunk==0);
CT(idx)=1;
% handle new sequqnces which may contain trained chunk transitions
trials=find(D.seqCat==3|D.seqCat==2);
for trial=trials';
    % get chunk transition
    chunk = D.chunk(trial,:);
    idx = find(chunk==0);
    for i=1:numel(idx)
        tmp(1,1) = chunk(1,idx(i)-1);
        tmp(1,2) = chunk(1,idx(i)+1);
        tmp2=zeros(size(ctran,1),1);
        for c=1:size(ctran,1);
            tmp2(c)=isequal(tmp, ctran(c,:));
        end
        if any(tmp2)
            CT(trial, idx(i))=1;
        end
    end
end
X_chunk_transition = reshape(CT', [],1);

% 3. Transigion frequency (only frequency matters)
% (e.g., Stadler 1992).
tranfreq=tranfreq/max(tranfreq(tranfreq>0));
AllTransitions = reshape(D.transition',[],1);
X_transition_freq = tranfreq(AllTransitions)';

% 4. Biomechanical effect
biomech = reshape(biomech',[],1);
biomech = biomech-nanmean(biomech(:));
biomech = biomech/nanstd(biomech(:));
X_biomech = biomech(AllTransitions);

% 5. Post-error slowing
% (e.g., Botvinick 2001)
X_error = reshape(D.errorVec(:,1:10)',[],1);

% 6. Explicit knowledge
% (e.g., Verwey 2010; Wong 2015)
EK=zeros(size(D.withinchunk));
EK(D.seqCat==1,:) = 1;
X_explicit = reshape(EK',[],1);

% 7. Intercept
X_intercept = ones(size(X_chunk));

% Adjust direction of effect
X_chunk=-X_chunk;
X_chunk_transition=-X_chunk_transition;
X_transition_freq=-X_transition_freq;
X_explicit=-X_explicit;

% Concatenate Xs
X = [X_chunk, X_chunk_transition, X_transition_freq, X_explicit, X_biomech, X_error, X_intercept];
Xname = {'Chunk','Sequence','FingerTransitionFreq','ExplicitKnowledge','Biomech.Difficulty','PostErrorSlowing','Intercept'};
%X = [X_chunk, X_chunk_transition, X_transition_freq, X_biomech, X_error, X_intercept];
%Xname = {'Chunk','Sequence','FingerTransitionFreq','Biomech.Difficulty','PostErrorSlowing','Intercept'};

end
function [Tran] = linregress_get_tranfreq(analyzeDir,subject_excluded) % calculate transition frequency from training phase data
% ===================================== %
% load data for training and calc summary
% ===================================== %
fname_training = fullfile(analyzeDir, 'All_training.mat');
T = load(fname_training);
T = getrow(T, ~ismember(T.subjname, subject_excluded));
% Individually calculate transition frequency based on goal press
edges=[1:26];
Tran=[];
for s=unique(T.SN)'
    D=getrow(T, T.SN==s);
    transitions = D.transition(D.transition>0&isfinite(D.transition));
    [freq, edge] = histcounts(transitions, edges);
    
    tran.SN=s;
    tran.freq=freq;
    Tran=addstruct(Tran, tran);
end
end
function [Result, compidx] = fitModelFamilyCrossval(X,Xname,Y,isOK,execution,folds) % Fit all combinations (intercept are always in) with cross-validation
[Nrow,Ncol] = size(X);
compidx = [];
for n=1:Ncol-1
    components = zeros(1,Ncol-1);
    components(1:n) = 1;
    compidx = cat(1, compidx, unique(perms(components), 'rows'));
end
% add intercept to all and create null model (only intercept)
compidx = cat(1,compidx, zeros(1,size(compidx,2)));
compidx = cat(2,compidx, ones(size(compidx,1),1));
Result = [];
for mf=1:size(compidx,1);
    %fprintf('%d',mf);
    % Setup model families
    modelidx = find(compidx(mf,:)==1);
    modelnames = '';
    for i=1:numel(modelidx)
        if i>1
            if modelidx(i)==numel(modelidx) % intercept
            else
                modelnames=[modelnames, '+', Xname{modelidx(i)}];
            end
        else
            if modelidx(i)==numel(modelidx) % intercept
                modelnames = 'Intercept';
            else
                modelnames = [modelnames, Xname{modelidx(i)}];
            end
        end
    end
    Xfam(mf).X = X(:,modelidx);
    Xfam(mf).name = modelnames;
    
    try
        % Do regression
        y = Y(isOK,1);
        x = Xfam(mf).X(isOK,:);
        k = execution(isOK,1);
        [loglike, beta, sigma, res_cv] = crossvalKfold(x,y,k,folds);
        
        % concat result
        result.ndata = size(x,1);
        result.model = mf;
        result.compidx = [compidx(mf,:)];
        result.name = {modelnames};
        result.beta = zeros(1,Ncol);
        result.beta(find(result.compidx==1)) = nanmean(beta,2);
        result.sig = sigma;
        result.predR2 = 1-(res_cv'*res_cv)/(y'*y);
        result.loglike = loglike;
        Result = addstruct(Result,result);
    catch
    end
end
%fprintf('\n');
end
function [loglike, beta, sigma, res_cv] = crossvalKfold(x,y,k,folds) % Do k-fold cross-validation
res_cv = zeros(size(y));
beta = [];%zeros(size(x,2), numel(folds));
sigma = [];%zeros(numel(folds,2),1);
kk=0;
for f = 1:numel(folds)
    test = ismember(k,folds{f});
    train = ~test;
    if sum(test)>0
        kk=kk+1;
        % make test and training data
        ytrain = y(train);
        ytest = y(test);
        xtrain = x(train,:);
        xtest = x(test,:);
        
        % fitting
        pX = (xtrain'*xtrain)\xtrain';
        beta(:, kk) = pX*ytrain;
        
        % residual
        res_uncv = ytrain-xtrain*beta(:, kk);
        res_cv(test,1) = ytest-xtest*beta(:, kk);
        
        % est error variance (uncv)
        [ndata, nmodel] = size(xtrain);
        sigma(kk,1) = (res_uncv'*res_uncv)/(ndata-nmodel);
    end
end

% calc crossvalidated log-likelihood
[n, nk] = size(x);
%sig = mean(sigma);
sig = (res_cv'*res_cv)/(n-nk);
loglike = -0.5*n*log(2*pi) - 0.5*n*log(sig) - 0.5*(1/sig) * (res_cv'*res_cv);
end
function ax=myColorbar(varargin)
%% Draw color bar in current axes

bgcolor = 'w';
Nlevel = 100;
scaling = [0,1];
MAP = [];
ylabelname = [];
barname = [];
dir = 'horz';
vararginoptions(varargin, {'bgcolor','Nlevel','scaling','MAP','ylabelname','barname','dir'});

% handle color
if isempty(MAP)
    MAP = parula(Nlevel);
elseif ischar(MAP)
    MAP = eval(sprintf('%s(Nlevel)', MAP));
end

% handle scaling
xtickls=[scaling(1), scaling(2)];
xticks=[1, Nlevel];
xtickls=[scaling(1), scaling(1)+diff(scaling)/2, scaling(2)];
xticks=[1,ceil(Nlevel/2),Nlevel];
xticklabel = cell(length(xticks),1);
for i=1:length(xticks)
    xticklabel{i} = sprintf('%1.2f',xtickls(i));
end

switch bgcolor
    case 'w'
        axcolor = [0 0 0];%[0.3,0.3,0.3];
    case 'k'
        axcolor = [1 1 1];
    otherwise
        axcolor = [0.5 0.5 0.5];
end

% directly draw color bar
for i=1:Nlevel
    patch([i-1,i-1,i,i],[0 1 1 0],MAP(i,:),'edgecolor','non');hold on;
end;
% draw edge
rectangle('position', [1,0,Nlevel-1,1],'edgecolor','k','facecolor','non');

ylabel(ylabelname); %xlabel(valname);
title(barname,'color',axcolor);
set(gca,'xlim',[1,Nlevel],'xtick',xticks,...
    'xticklabel',xticklabel,...
    'ylim',[0,1],'ytick',[],'tickdir','out','ticklength',[0.03,0.1]);
set(gca,'xcolor',axcolor,'ycolor',axcolor);
switch dir
    case 'vert'
        set(gca,'view',[90,-90]);
    otherwise
end
%axis off;
set(gcf,'color',bgcolor);
ax=gca;
end