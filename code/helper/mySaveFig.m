function mySaveFig(handle,figname,flag,format,resolution) 
%% function mySaveFig(handle,figname,flag,format) 
% 
% Saves figure in .eps or other format.
% 
% Input:
%   handle:     figure handle to be saved
%   figname:    figure name to be saved, could be either full-path or just name
%   flag:       0/1 flag, if set to 0, program does not save figure. default is 1
%               frees you from writing 'if' every time
%   format:     format that the figure is saved as
%                   '-dpsc2': .eps
%                   '-djpeg': .jpeg
%                   '-dpng': .png
% 
%   resolution: resolution
%                   '-r300' gives 300 dpi resolution
% 
% Example:
%   x = 0:360;
%   y = sind(x);
%   h = myFigure([10,10],'name','hoge'); % create figure window
%   plot(x,y);
%   mySaveFig(h,'hoge.eps',1,'dpsc2'); % save figure in .eps format
% 
% see also myFigure.m, print.m
% 
% a-yokoi (2017)

[fdir,fname,fext] = fileparts(figname);
switch (nargin)
    case 2 % default is .eps
        flag        = 1;
        format      = '-dpsc2';
        fext        = '.eps';
        resolution  = '';
    case 3
        format      = '-dpsc2';
        fext        = '.eps';
        resolution  = '';
    case 4        
        switch (format)
            case '-dpsc2'
                fext = '.eps';
            case '-dpsc'
                fext = '.eps';
            case '-dpdf'
                fext = '.pdf';
            case '-dtiff'
                fext = '.tiff';
            otherwise
                fext = [];
        end
        resolution = '';
    case 5
        resolution = resolution;
    otherwise
        error('too many or too less input.');
end
if flag==1;
    set(handle,'Renderer','painters');
    figname = fullfile(fdir,[fname,fext]);
    print(handle,format,resolution,figname);
    set(handle,'Renderer','opengl'); % for faster rendering of patch objects (e.g., brain map)
end;

end