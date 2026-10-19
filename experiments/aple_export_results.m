function aple_export_results(folder)
% Derive lightweight tables; checkpoints remain the authoritative resume state.
families={'snr','grid','aperture','range','rf','iterations','coverage'}; states={};
for k=1:numel(families)
    file=fullfile(folder,[families{k} '.mat']);
    if isfile(file) || isfile([file '.previous.mat']) || isfile([file '.pending.mat'])
        states{end+1}=aple_checkpoint('load',file,''); %#ok<AGROW>
    end
end
[summary,cdf_table,meta]=aple_tables(states);
if ~isempty(summary), writetable(summary,fullfile(folder,'summary.csv'),'Encoding','UTF-8'); end
writetable(cdf_table,fullfile(folder,'cdf.csv'),'Encoding','UTF-8');
file=fullfile(folder,'metadata.json'); fid=fopen(file,'w','n','UTF-8'); assert(fid>=0);
clean=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(meta,'PrettyPrint',true));
end
