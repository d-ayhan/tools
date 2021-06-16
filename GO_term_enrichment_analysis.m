function GOresults = GO_term_enrichment_analysis(upgenes, GO)

% upgenes: the subset of the genes to analyze for GOterm enrichement
%          vector of gene ids 
% GO:      all annotations in the chipset 
%          table with following headers: ID, gotermId, goName, gotermType, goAcc
% upgenes needs to match GO.ID

t = {'bio','mol','cel'};
go.bio = GO(GO.gotermType == "biological_process",:);
go.mol = GO(GO.gotermType == "molecular_function",:);
go.cel = GO(GO.gotermType == "cellular_component",:);

for i = 1:3
    [goterms,ia,~] = unique(go.(t{i}).gotermId);  % all go terms
    gotermsdef = go.(t{i})(ia,:);
    m = size(goterms,1);
    
    chipgenesCount.(t{i}) = zeros(m,1);  % total gene number counts for each goterm
    upgenesCount.(t{i})  = zeros(m,1); % de gene count for each goterm
    
    for goid = 1:length(goterms)
        subsetgenes = go.(t{i}).ID(go.(t{i}).gotermId == goterms(goid)); % genes in go term
        chipgenesCount.(t{i})(goid,1) = size(subsetgenes,1);   % number of genes in go term
        upgenesCount.(t{i})(goid,1) = sum(ismember(upgenes,subsetgenes)); % number of de genes in go term
    end
    x = upgenesCount.(t{i}); % number of DE genes in specific GO term
    k = chipgenesCount.(t{i}); % number of genes in specific GO term
    N = sum(upgenesCount.(t{i})); % number of DE genes with GO term
    M = sum(chipgenesCount.(t{i})); % number of all genes with GO terms
    gopvalues = nan(size(x));
    for j = 1:length(x)
        gopvalues(j,1) = hygecdf(x(j,1), M, k(j,1), N, 'upper');
    end
    [pval, idx] = sort(gopvalues);
    
    counts_subset = upgenesCount.(t{i})(idx);
    counts_chipset = chipgenesCount.(t{i})(idx);

    desc = gotermsdef.goName(idx);
    goAcc = gotermsdef.goAcc(idx);
    GOresults.(t{i}) = table(desc, goAcc, pval, counts_subset,counts_chipset);
    GOresults.(t{i})(GOresults.(t{i}).counts_subset == 0,:) = [];
end
end

