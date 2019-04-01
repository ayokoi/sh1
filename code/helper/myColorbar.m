function varargout=myColorbar(varargin)
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

if ischar(bgcolor)
    switch bgcolor
        case 'w'
            axcolor = [0 0 0];%[0.3,0.3,0.3];
        case 'k'
            axcolor = [1 1 1];
        otherwise
            axcolor = [0.5 0.5 0.5];
    end
else
    axcolor = [1,1,1]-bgcolor;
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

end