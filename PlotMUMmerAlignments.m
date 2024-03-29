function PlotMUMmerAlignments(COOR, QRY, REF, varargin)
% This function plots sequence alignment between 2 genomes from MUMmer
% nucmer alignments using coordinates file generated by show-coords
% option.
% % %
% % Example command line -basic use-
% %
% PlotMUMmerAlignments('nucmerAlignment.coor', 'Query.fasta', 'Reference.fasta')
%
% % %
% % Example command line -advanced use-
% %
% PlotMUMmerAlignments('nucmerAlignment.coor','Query.fasta','Reference.fasta',...
%   'sortfileQRY','QRY_seq_sort.txt','sortfileREF','REF_seq_sort.txt', ...
%   'HMThreshold',90,'allticks',0,'LastXTick',20,'LastYTick',30,'HeaderLines',5)
%
% % %
% % optional arguments [type] : description (deafult)
% sortfileQRY    [string] : single column text file with fasta headers QRY 
%                           (QRY order)
% sortfileREF    [string] : single column text file with fasta headers of REF 
%                           (REF order)
% HMThreshold    [double] : sequence identitiy % down limit for colorbar (75)
% allTicks      [boolean] : show sequence borders as ticks (true)
% grid          [boolean] : show grids (true)
% LastXTick     [integer] : remove ticks on x-axis after nth sequence, 0 for 
%                           show all ticks (0)
% LastYTicks    [integer] : remove ticks on y-axis after nth sequence, 0 for 
%                           show all ticks (0)
% HeaderLines   [integer] : nth line is the last line of the headers in the 
%                           coordinates file (5)
% figVisibility [boolean] : matlab figure visiblility (1)
%
% % %
% version 0.2.2
% Dilay Ayhan, 10/2021

%% input parsing
defaultsortfileQRY = '';
defaultsortfileREF = '';
defaulthmThreshold = 75;
defaultallTicks = 1;
defaultGRID = 1;
defaultLastXTick = 0;
defaultLastYTick = 0;
defaultHeaderLines = 5;
defaultfigVisibility = true;

p = inputParser;

validScalarPosNum = @(x) isnumeric(x) && isscalar(x) && (x >= 0);
validDubNum = @(x) isnumeric(x) && (x <=100) && (x >= 0);
validLogical = @(x) islogical(x) || (x==1) || (x==0);
validStringChar = @(x) isstring(x) || ischar(x);

addRequired(p,'COOR',validStringChar);
addRequired(p,'QRY',validStringChar);
addRequired(p,'REF',validStringChar);
addParameter(p,'sortfileQRY',defaultsortfileQRY, validStringChar);
addParameter(p,'sortfileREF',defaultsortfileREF, validStringChar);
addParameter(p,'HMThreshold',defaulthmThreshold, validDubNum);
addParameter(p,'allTicks',defaultallTicks, validLogical);
addParameter(p,'grid',defaultGRID, validLogical)
addParameter(p,'LastXTick',defaultLastXTick,validScalarPosNum);
addParameter(p,'LastYTick',defaultLastYTick,validScalarPosNum);
addParameter(p,'HeaderLines',defaultHeaderLines,validScalarPosNum);
addParameter(p,'figVisibility', defaultfigVisibility, validLogical);

parse(p,COOR,QRY,REF,varargin{:});

allTicks = p.Results.allTicks;
GRID = p.Results.grid;
HeaderLines = p.Results.HeaderLines;
hmThreshold = p.Results.HMThreshold;
LastXTick = p.Results.LastXTick;
LastYTick = p.Results.LastYTick;
sortfileQRY = p.Results.sortfileQRY;
sortfileREF = p.Results.sortfileREF;
figVisibility = p.Results.figVisibility;

%% create sorted fai like structure
for k = 1:2 % for REF and QRY
    if k == 1
        f = fastaread(REF);
        sortfile = sortfileREF;
    else
        f = fastaread(QRY);
        sortfile = sortfileQRY;
    end
    fH = {f.Header}';
    for i = 1: length(fH)
        tem = strsplit(fH{i},' ');  % remove string after space in the header
        fH{i} = tem{1};
   end
    if ~isempty(sortfile) % sort fasta
        SF = importdata(sortfile);
        for i = 1: length(SF)
            tem = strsplit(SF{i},' ');
            SF{i} = tem{1};
        end
        [lia, loc] = ismember(SF,fH);
        g = f(loc(lia));
        fH = SF;
    else
        g = f;
    end
    
    for i = 1:length(g)  % cumulative sum for sorted sequences
        F(k).l(i,1) = length(g(i).Sequence);
    end
    F(k).cl = cumsum(F(k).l);
    F(k).H = fH;
end

clearvars('f','i','k', 'g', 'lia', 'loc', 'SF', 'sortfile')

%% read coordinates file generated by mummer/show-coor
fileID = fopen(COOR);
C = textscan(fileID,'%d %d | %d %d | %d %d | %f | %s\t%s','HeaderLines',HeaderLines);
fclose(fileID);
pos = [C{1} C{2} C{3} C{4} C{5} C{6}];
id = C{7};
seq = [C{8} C{9}];

clearvars('fileID', 'C')

%% remove sequences that are not in the sort files
removetheserows = ~(ismember(seq(:,1),F(1).H) & ismember(seq(:,2),F(2).H));
seq(removetheserows,:) = [];
id(removetheserows,:) = [];
pos(removetheserows,:) = [];

%% change orientation of QRY elments if necessary
pos_new = pos;
for i = 1: length(F(2).H)
    q = F(2).H{i};
    [Lia, ~] = ismember(seq(:,2),q);
    if any(Lia)
        r = seq(Lia,1);
        [uniqueC,~,~] = unique(r);
        uniqueCcount=[];
        for j = 1: length(uniqueC)
            uniqueCcount(j,1)=sum(pos(Lia & ismember(seq(:,1),uniqueC{j}),5));
        end
        r = uniqueC{find(uniqueCcount ==  max(uniqueCcount),1)};
        qr = Lia & ismember(seq(:,1),r);
        if sum(pos(qr,4) - pos(qr,3)) < 0
            % the orientation of the QRY is wrong I need to flip it in the
            % coordinates file for all of the instances
            len = F(2).l(i);
            
            pos_new(Lia,3) = len - pos(Lia,3);
            pos_new(Lia,4) = len - pos(Lia,4);
        end
    end
end
pos = pos_new;

%% plot
if figVisibility
    Visible = 'on';
else
    Visible = 'off';
end

f = figure('PaperSize', [10 10],'Visible', Visible);

ax = axes();
fhot = flipud(hot);
lh = length(fhot);

for i = 1:size(id,1)
    P = pos(i,:);
    I = id(i,1);
    if I < hmThreshold
        I = hmThreshold;
    end
    S = seq(i,:);
    
    % position
    for j = 1:2
        ind = find(ismember(F(j).H, S(j)));
        if ind < 2
            add = 0;
        else
            add = F(j).cl(ind-1);
        end
        P(2*j-1:2*j) = P(2*j-1:2*j) + add;
    end
    
    % color
    k = floor(((lh-1)./(100-hmThreshold)).*I + (100 - lh.*hmThreshold)./(100-hmThreshold));
    
    line(ax, P(1:2),P(3:4), ...
        'Marker','.','Color', fhot(k,:), 'LineStyle','-', 'LineWidth',2)
end
colormap(fhot)
caxis([hmThreshold 100])

if allTicks && LastXTick == 0
    LastXTick = length(F(1).cl);
end
if allTicks && LastYTick == 0
    LastYTick = length(F(2).cl);
end
if GRID
    ax.XGrid = 'on';
    ax.YGrid = 'on';
else
    ax.XGrid = 'off';
    ax.YGrid = 'off';
end
ax.XTick = F(1).cl(1:LastXTick);
ax.YTick = F(2).cl(1:LastYTick);
ax.XTickLabel = F(1).H(1:LastXTick);
ax.XTickLabel = F(1).H(1:LastXTick);
ax.YTickLabel = F(2).H(1:LastYTick);

ax.XLim = [1, F(1).cl(end)];
ax.YLim = [1, F(2).cl(end)];

ax.DataAspectRatio = [1 1 1];
ax.Box = 'on';
ax.XTickLabelRotation = 90;
colorbar
ax.TickLabelInterpreter='none';

print(f,[COOR '.pdf'],'-dpdf', '-bestfit','-painters')
end
