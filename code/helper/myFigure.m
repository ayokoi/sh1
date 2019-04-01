function handle = myFigure(figsize,varargin) 
%% function handle = myFigure(figsize,varargin) 
% 
% Creates figure window ready for save in .eps with minimal input.
% 
% Inputs: 
%   figsize:    1x2 (or 2x1) vector of figure size [width, height] in cm
%   varargin:   additional options passed to figure.m
% 
% Output:
%   handle:     figure handle
%   
% Example:
%   x = 0:360;
%   y = sind(x);
%   h = myFigure([10,10],'name','hoge'); % create figure window
%   plot(x,y);
%   mySaveFig(h,'hoge.eps'); % save figure in .eps format
% 
% see also: mySaveFig.m
% 
% a-yokoi (2017)

switch isempty(varargin)
    case 1
        handle = figure('units','centimeters',...
            'paperunits','centimeters',...
            'papersize',figsize    ,...
            'position',[1,1,figsize(1),figsize(2)],...
            'paperposition',[0,0,figsize(1),figsize(2)],...
            'paperpositionmode','manual');
    otherwise
        handle = figure('units','centimeters',...
            'paperunits','centimeters',...
            'papersize',figsize    ,...
            'position',[1,1,figsize(1),figsize(2)],...
            'paperposition',[0,0,figsize(1),figsize(2)],...
            'paperpositionmode','manual',...
            varargin{:});
end
% other possible options
% 'Color',bgcolor,'InvertHardcopy','off',...
% 'PaperPositionMode','auto'
end
