function scaffolding(FASTA,FAI,BNDvcf,insertSize,OUTFASTA)
% scaffolding
% Breakend vcf file should be sorted by coordinates

% FASTA = '26365_pilon.fasta';
% OUTFASTA = '26365_scaffolding.fasta';
% BNDvcf = '26365_v3.SV.vcf';
% FAI = '26365_pilon.fasta.fai';
% insertSize = 350;
BNDpairTag = '(?:PARID|MATEID)';
ScaffoldName = 'Scaffold_';

% load fasta file 
fa = fastaread(FASTA);
newfa = fa;

% load index file
formatSpec = '%s%n%n%n%[^\n]';
fileID = fopen(FAI,'r');
dataArray = textscan(fileID, formatSpec, 'Delimiter', '\t', ...
    'TextType', 'string', 'ReturnOnError', false, 'EndOfLine', '\n');
fclose(fileID);
chr = dataArray{1,1};
len = dataArray{1,2};

% load SV.VCF file
formatSpec = '%s%n%q%s%s%n%q%q%q%q%[^\n]';
fileID = fopen(BNDvcf,'r');
dataArray = textscan(fileID,formatSpec,'Delimiter', '\t', ...
    'TextType', 'string','CommentStyle','#',...
    'ReturnOnError', false, 'EndOfLine', '\n');
fclose(fileID);
SV.text = [dataArray{1,1}, dataArray{1,3}, dataArray{1,4}, dataArray{1,5}, dataArray{1,7}, dataArray{1,8}];
% SV.num colums: {1)chr, 2)pos, 3)chr_len, 4)qual, 5)mate_chr, 6)mate_pos, 7)mate_chr_len}
SV.num = [nan(length(dataArray{1,2}),1), dataArray{1,2}, nan(length(dataArray{1,2}),1), dataArray{1,6}];
% clear('ans','dataArray','fileID','formatSpec')

% find chr index
[Lia, Locb] = ismember(SV.text(:,1),chr);
if ~all(Lia)
    disp('ERROR: At least one sequence name does not match with index file')
    return
end
SV.num(:,3) = len(Locb);
SV.num(:,1) = Locb;
for i=1:length(chr)
    SV.num(contains(SV.text(:,4),chr(i)), 5) = i;
end

% remove one of the pairs in the vcf record list
SV.F_num = SV.num;
SV.F_text = SV.text;
i = 1;
check = 1;
while check
    t = SV.F_text{i,6};
    mate = regexp(t,[BNDpairTag '=(.*?);'],'tokens');
    if isempty(mate)
        disp(['WARNING: Mate info tag for ' SV.F_text{i,2} ' is not found. Exiting pair removal process.'])
        break
    end
    [Li, in] = ismember(mate{1},SV.F_text(:,2));
    if Li
        SV.F_text(in,:) = [];
        SV.F_num(in,:) = [];
    else
        disp(['WARING: Mate record for ' SV.F_text{i,2} ' was not found.']);
    end
    i = i + 1;
    if i >= length(SV.F_num)
        check = 0;
    end
end

% add link orientation and intermediate seqs.
% SV.F_num colums: {1)chr, 2)pos, 3)chr_len, 4)qual, 5)mate_chr, 6)mate_pos,
%        7)mate_chr_len, 8)link_orientation}
for i = 1:size(SV.F_num,1)
    LO = DetermineLinkOrientation(SV.F_text{i,4});
    SV.F_num(i,8) = LO{1};
    SV.F_text{i,7} = LO{2};
end

% Find mate positions
for i = 1:size(SV.F_text,1)
    alt = SV.F_text(i,4);
    t = char(regexp(alt,':\d*','match'));
    SV.F_num(i,6) = str2num(t(2:end));
end
SV.F_num(:,7) = len(SV.F_num(:,5));
% clear('Lia','Locb','check', 'b', 't', 'mate', 'Li', 'in','LO','alt')


% filters
%%% remove intercintig links
filter1 = SV.F_num(:,1) ~= SV.F_num(:,5);
%%% remove links to the middle of mate
filter2 = (SV.F_num(:,6) < insertSize) | ((SV.F_num(:,7) - SV.F_num(:,6)) < insertSize);
%%% remove links based on case and pos (if beginning>case 3 or 4; if end>case
%%%  1 or 2; remove links from middle)
filter3 = ((SV.F_num(:,2) < insertSize) & ((SV.F_num(:,8) == 3) | (SV.F_num(:,8) == 4))) | ...
    (((SV.F_num(:,3) - SV.F_num(:,2)) < insertSize) & ((SV.F_num(:,8) == 1) | (SV.F_num(:,8) == 2)));
filter0 = filter1 & filter2 & filter3;
if ~any(filter0)
    disp('ERROR: No links were found.')
    return;
end
fSV.num = SV.F_num(filter0,:);
fSV.text = SV.F_text(filter0,:);

% select the best record in multiple linkage records for same contigs pair
ucontigs = unique(fSV.num(:,1));
removeDuplicates = [];
for i = ucontigs'
    f = find(fSV.num(:,1) == i);
    records = fSV.num(f,5:7);
    u = unique(records(:,1));
    for k = 1:length(u)
        ind = records(:,1) == u(k);
        if sum(ind) > 1
        r = PosValueToEndRelative(records(ind,2),records(ind,3));
        if all(r < 0) || all(r > 0)
            [~,maxin] = nanmax(fSV.num(f(ind),4));
            kk = find(ind);
            ind(kk(maxin)) = 0;
            removeDuplicates = [removeDuplicates; f(ind)];
        else 
            removeDuplicates = [removeDuplicates; f(ind)];                
        end
        end
    end
end
fSV.num(removeDuplicates,:) = [];
fSV.text(removeDuplicates,:) = [];
% clear('ucontigs','f','records','ind','r','maxin','kk','removeDuplicates')

% reformat
m1 = PosValueToEndRelative(fSV.num(:,2),fSV.num(:,3));
m2 = PosValueToEndRelative(fSV.num(:,6),fSV.num(:,7));
d = [fSV.num(:,1), m1, fSV.num(:,5), m2, fSV.num(:,8)];

% remove multiple linkages
kk = [d(:,1:2);d(:,3:4)];
[~,I]=sort(kk(:,1));
B=kk(I,:);
uB = unique(B(:,1));
removeMulti = [];
for i = uB'
    f = find(B(:,1) == i);
    n = B(f,2) < 0;
    p = ~n;
    if sum(n) > 1
        removeMulti = [removeMulti; f(n)];
    end
    if sum(p) > 1
        removeMulti = [removeMulti; f(p)];
    end
end
IremoveMulti = I(removeMulti);
t= IremoveMulti-size(d,1);
IremoveMulti( t > 0) = t(t>0);
uIremoveMulti = unique(IremoveMulti);
d(uIremoveMulti,:) = [];
fSV.num(uIremoveMulti,:) = [];
fSV.text(uIremoveMulti,:) = [];
% clearvars('kk','I','B','uB','f','n','p','removeMulti','t','uRemoveMulti')

% Cluster the contigs for scaffolding
Clusters = [];
for i = 1:size(d,1)
    m1 = d(i,1);
    m2 = d(i,3);
    s = size(Clusters);
    if min(s) ~= 0
        [m1l, m1i] = ismember(Clusters,m1);
        [m2l, m2i] = ismember(Clusters,m2);
        
        if any(m1l(:)) && ~any(m2l(:))
            [a, b] = ind2sub(s,find(m1i));
            z = find(Clusters(a,:) == 0);
            if isempty(z)
                Clusters(a,end+1) = m2;
            else
                Clusters(a,z(1)) = m2;
            end
        elseif ~any(m1l(:)) && any(m2l(:))
            [a, b] = ind2sub(s,find(m2i));
            z = find(Clusters(a,:) == 0);
            if isempty(z)
                Clusters(a,end+1) = m1;
            else
                Clusters(a,z(1)) = m1;
            end
        elseif ~any(m1l(:)) && ~any(m2l(:))
            Clusters(end+1,1:2) = [m1,m2];
        else
            [b, ~] = ind2sub(s,find(m1i));
            [a, ~] = ind2sub(s,find(m2i));
            z = find(Clusters(a,:) == 0);
            w = find(Clusters(b,:));
            if isempty(z)
                Clusters(a,end+1:end+length(w)) = Clusters(b,w);
            else
                Clusters(a,z(1):z(1)+length(w)-1) = Clusters(b,w);
            end
            Clusters(b,:) = [];
        end
    else
        Clusters = [m1,m2];
    end
end
disp(['Total possible scaffolds: ' num2str(size(Clusters,1)) '.'])
% clearvars('m1','m2','s','m1l','m1i','m2l','m2i','a','b','z','w')

% combine the contigs
dcontigs = d(:,[1,3]);
RecordKeep = cell(size(Clusters,1),1);
for i = 1:size(Clusters,1)
    t = Clusters(i,Clusters(i,:)~=0);
    l = any(ismember(dcontigs, t),2);
    D = d(l,:);
    Dseq = fSV.text(l,7);
    x = D(:,[1,3]);
    x = x(:);
    y = sum(x==x');
    % z = find(y == 1);
    s = x(y == 1); % find the ends, assign beginning and end contigs.
    contigsorted = nan(size(D,1)+1,3);
    contigsorted(1,1) = s(1);
    tsorted = strings(size(contigsorted,1)-1,1);
    Dsorted = [];
    for j = 1:size(contigsorted)-1
        ia = any(ismember(D(:,[1 3]),contigsorted(j,1)),2);
        con_link = D(ia,:);
        if sum(ia,1) > 1
            [con_link, ib] = setdiff(con_link,old_con_link,'rows');
            tt= find(ia);
            ia(ia == 1) = 0;
            ia(tt(ib)) = 1;
        end
        tsorted(j) = Dseq(ia,1);
        Dsorted(j,:) = D(ia,:);
        if con_link(1) == contigsorted(j,1)
            contigsorted(j,3) = con_link(2);
            contigsorted(j+1,1) = con_link(3);
            contigsorted(j+1,2) = con_link(4);
        else
            contigsorted(j,3) = con_link(4);
            contigsorted(j+1,1) = con_link(1);
            contigsorted(j+1,2) = con_link(2);
        end
        old_con_link = con_link;
    end
    contigsorted_i = ~((contigsorted(:,2) > 0) | (contigsorted(:,3) < 0)); % reverse ones
    % replace NaNs
    if contigsorted(1,3) < 0
        contigsorted(1,2) = 1;
    else
        contigsorted(1,2) = -1;
    end
    if contigsorted(end,2) < 0
        contigsorted(end,3) = 1;
    else 
        contigsorted(end,3) = -1;
    end
    
    %combine sequences
    newfa(end+1).Header = [ScaffoldName num2str(i)];
    newseq = '';
    rk = '';
    for j = 1: size(contigsorted,1)
        sPos = max(contigsorted(j,[2,3]));
        ePos = min(contigsorted(j,[2,3]));
        cn = fa(contigsorted(j,1)).Header;
        s = fa(contigsorted(j,1)).Sequence(sPos:end+ePos+1);
        if j ~= size(contigsorted,1)
        extention = tsorted{j};
        if (Dsorted(j,5) == 1) || (Dsorted(j,5) == 2)
            s = [s, extention];
        else
            s = [extention, s];
        end
        end
        rk = [rk,cn,':',num2str(sPos),'-',num2str(len(contigsorted(j,1))+ePos+1),'...'];
        if contigsorted_i(j)
            s = seqrcomplement(s);
        end
        
        newseq = [newseq s];
    end
    newfa(end).Sequence = newseq;
    RecordKeep{i,1} = rk(1:end-3);
end
disp([num2str(size(Clusters,1)) ' clusters from ' num2str(numel(dcontigs)) ...
    ' contigs were generated.'])
% remove the contigs that are present in the scaffolds.
newfa(dcontigs(:)) = [];
disp(['Writing fasta file with ' num2str(length(newfa)) ' sequences.'])
myfastawrite(OUTFASTA, newfa)
fileID = fopen([OUTFASTA '.links'],'w');
formatSpec = '%s\n';
for i = 1:size(RecordKeep,1)
    fprintf(fileID,formatSpec,['Scaffold_' num2str(i) ' ' RecordKeep{i}]);
end
fclose(fileID);
disp('Done.')
end