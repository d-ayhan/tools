function FASTAeditWithVCF(VCF,FASTA,OUTFASTA)
% This script edits given FASTA file with the variants in the VCF file.
% VCF has to be sorted by coordinates.
% if OUTFASTA exist, the script will expand the file. Delete the previously
% exising file before running.
% % Example command line
% FASTAeditWithVCF('Variants.vcf','in.fasta','out.fasta')
% Dilay Ayhan, 2018

formatSpec = '%s%n%q%s%s%n%q%q%q%q%[^\n]';
fileID = fopen(VCF,'r');
dataArray = textscan(fileID,formatSpec,'Delimiter', '\t', ...
    'TextType', 'string','CommentStyle','#',...
    'ReturnOnError', false, 'EndOfLine', '\n');
fclose(fileID);

FA=fastaread(FASTA);
newFA=FA;
FAheaders={FA.Header}';

variants = {flipud(dataArray{1}), flipud(dataArray{2}), flipud(dataArray{4}),...
    flipud(dataArray{5})}; % reverse order for fixing without causing problems

[l,p] = ismember(variants{1},FAheaders); % find entries
if ~min(l) % if l has 0
    disp('VCF has following CHROM entries that are not present in the FASTA');
    disp(unique(dataArray{1}(~l)))
    return
end

for i=1:length(variants{1,1})
    newFA(p(i)).Sequence = replaceBetween(newFA(p(i)).Sequence,...
        variants{1,2}(i),variants{1,2}(i)+strlength(variants{1,3}(i))-1,...
        variants{1,4}(i));
end

fastawrite(OUTFASTA,newFA);

end
