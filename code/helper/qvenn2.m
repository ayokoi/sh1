function varargout = qvenn2(sets, varargin)
%% varargout = venn3(cases, varargin)
% Draw a simple quantitative venn diagram of 2 sets.
% 
% 
% 
% 
% 
% 
% 
% 
% 
tolf = 0.0001;
linewidth = 2;
linecolor = {[1,0,0],[0,1,0]};
orientation = 90;
fontname = 'Helvetica';
fontsize = 15;
leg = {'Set1','Set2'};
printnum = [1,1];
vararginoptions(varargin, {'tolf','orientation','linewidth','linecolor','fontname','fontsizse','leg','printnum'});

% input check
[Nobs, Nset] = size(sets);
if Nobs==1|Nset<2
    error;
end

% check set size and intersection
for s=1:Nset
    setsize(s) = sum(sets(:,s)==1);
    radius(s) = sqrt(setsize(s)/pi);
end
uniquen(1) = sum(sets(:,1)==1&~sets(:,2)==1);
uniquen(2) = sum(sets(:,2)==1&~sets(:,1)==1);
intersec12 = sum(sets(:,1)==1&sets(:,2)==1);

setvec = [uniquen(1)/setsize(1), uniquen(2)/setsize(2)];

% find optimal distance between two circles that matches the observation
dist12 = sum(radius);
e = 1; iter = 0; preve = 1e5; hold off;
if intersec12>0
    while e>tolf % simple search
        iter = iter+1;
        [n1, n2, n11, n22, n12] = countoverlap(radius(1), radius(2), dist12);
        currvec = [n11/n1, n22/n2];
        e = (setvec-currvec);
        e = e*e';
        de = preve-e;                
        dist12 = 0.99*dist12;
        preve = e;
        if (0)
            figure(99); 
            subplot(2,1,1);
            plot(iter,e,'ko'); hold on
            ylabel('Squared error');
            subplot(2,1,2);
            plot(iter, de, 'ko'); hold on
            xlabel('Iterations');
            ylabel('e(n)/e(n-1)');
        end
        if mod(iter,100)==0
            %keyboard();
            disp('100 iterations passed. Reset.');
            dist12=sum(radius);
        end
    end
    
else
    dist12 = 1.3*sum(radius);
end

% draw circles
centre = {[0,0],[dist12*cosd(orientation),dist12*sind(orientation)]};
theta = [0:360];
hold on; axis off; axis equal; 
for c=1:Nset    
    x = centre{c}(1)+cosd(theta)*radius(c);
    y = centre{c}(2)+sind(theta)*radius(c);
    plot(x,y,'-', 'linewidth', linewidth,'color',linecolor{c}); hold on    
end

% text
unitvec = {-centre{2}/dist12, centre{2}/dist12};
for i=1:2
    % labels
    p = centre{i} + unitvec{i}*radius(i)*1.3;
    text(p(1), p(2), leg{i},...
    'horizontalalignment', 'center',...
    'verticalalignment', 'middle',...
    'fontname', fontname,...
    'fontsize', fontsize,...
    'color', linecolor{i});
end
if ~isempty(printnum)
    for i=1:2
        if printnum(i)==1
            % only current set is true
            p = centre{i} + unitvec{i}*radius(i)*0.7;
            text(p(1), p(2), sprintf('%d', uniquen(i)),...
                'horizontalalignment', 'center',...
                'verticalalignment', 'middle',...
                'fontname', fontname,...
                'fontsize', fontsize,...
                'color',get(gca,'xcolor'));
        end
    end
    
    % intersection
    if intersec12>1
        y=dist12/(sum(radius))*radius(1);
        text(0, y, sprintf('%d', intersec12),...
            'horizontalalignment', 'center',...
            'verticalalignment', 'middle',...
            'fontname', fontname,...
            'fontsize', fontsize,...
            'color',get(gca,'xcolor'));
    end
end


varargout = {};

end
%% local
function [n1, n2, n11, n22, n12] = countoverlap(r1, r2, d0, varargin)
Ndots = 5000;

% define range
Xrange = [-r1, r2+d0];
Yrange = [-max([r1,r2]), max([r1,r2])];

% generate particles
p(:,1) = unifrnd(Xrange(1), Xrange(2), Ndots,1);
p(:,2) = unifrnd(Yrange(1), Yrange(2), Ndots,1);

% count overlaps (center1:[0,0], center2:[d0,0])
p1 = p; 
p2 = [p(:,1)-d0, p(:,2)];
idx1 = sum(p1.^2,2) <= r1^2;
idx2 = sum(p2.^2,2) <= r2^2;

n1 =  sum(idx1);
n2 = sum(idx2);
n11 = sum(idx1==1&~idx2==1);
n22 = sum(idx2==1&~idx1==1);
n12 = sum(idx1==1&idx2==1);
end