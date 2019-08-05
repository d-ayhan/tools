function assemblystats(FASTA,REF)
% prints basic stats and plots histograms for a given fasta file (FASTA),
% and/or prints stats of comparison the assembly to a reference (REF)
% % example usage:
% assemblystats('/path/to/Genome.fasta')
% Dilay Ayhan, 2019

if nargin > 1
    ref = fastaread(REF);
else
    ref = [];
end
%%
a = fastaread(FASTA);
n=length(a);
len=nan(n,1);
Ns=nan(n,1);
Gaps=nan(n,1);
ATs=nan(n,1);
GCs=nan(n,1);

for i=1:n
    seq=a(i).Sequence;
    len(i,1)=length(seq);
    Ns(i,1)=sum(seq=='N');
    Gaps(i,1)=floor(length(SplitVec(seq=='N'))/2);
    ATs(i,1)=sum(seq=='A'|seq=='T');
    GCs(i,1)=sum(seq=='G'|seq=='C');
end
GCcontent=GCs./(GCs+ATs);

g=sum(len);
sortedlengts=sort(len,'descend');
L50=find(g/2 < cumsum(sortedlengts), 1, 'first');
N50=sortedlengts(L50);
N90=sortedlengts(find(g.*.9 < cumsum(sortedlengts), 1, 'first'));
smallest=min(len);
largest=max(len);
largerThan1000 = sum(len>1000);
%% esimated genome size from alignments to reference genome. %ignores high coverage regions
if ~isempty(ref)
    refgenomeSize=length([ref(:).Sequence])-sum([ref(:).Sequence]=='N'); %nonNs
    half=0.128;
    cov1x=1-half;
    estimated_g=((refgenomeSize.*half)./2)+refgenomeSize.*cov1x;
    covered_genome_percent=g/estimated_g.*100;
    NG50=sortedlengts(find(estimated_g/2 < cumsum(sortedlengts), 1, 'first'));
end

%%
figure
subplot(2,1,1)
histogram(GCcontent)
title('GC% Distribution of Contigs')
xlabel('GC%')
ylabel('Number of Contigs')
xlim([0 1])
subplot(2,1,2)
plot(len,GCcontent,'ok')
xlabel('Length')
ylabel('GC%')

figure
histogram(len)
title('Size Distribution of Contigs')
xlabel('Contig Size')
ylabel('Number of Contigs')

%%
display(['Assembly size = ' num2str(g)])
if ~isempty(ref)
    display(['Reference genome size = ' num2str(refgenomeSize)])
    display(['Estimated genome size = ' num2str(estimated_g)])
    display(['% of covered genome = ' num2str(covered_genome_percent)])
end

display(['Number of contigs & supercontigs = ' num2str(length(len))])
display(['Size of largest supercontig = ' num2str(largest)])
display(['Size of smallest contig = ' num2str(smallest)])
display(['Number of contigs with size larger than 1kb = ' num2str(largerThan1000)])
display(['N50 = ' num2str(N50)])
display(['N90 = ' num2str(N90)])

if ~isempty(ref)
    display(['NG50 = ' num2str(NG50)])
end

display(['L50 = ' num2str(L50)])
display(['%GC = ' num2str(100.*sum(GCs)./(sum(GCs)+sum(ATs)))])
end
