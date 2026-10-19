function output=plot_results(source,varargin)
%PLOT_RESULTS Draw M1-M8 in journal order from reference CSVs or a completed run.
% plot_results('reference','Language','zh') redraws the accepted thesis results.
% Semantic figure names (grid, snr, cdf, aperture, rf, range, coverage,
% iterations) are also accepted; see aple_figure_catalog for numeric IDs.
root=setup_aple;
if nargin<1 || strcmp(source,'reference'), source=fullfile(root,'reference'); end
source=char(source); p=inputParser;
addParameter(p,'OutputDir','',@(x)ischar(x)||isstring(x));
addParameter(p,'Language','en',@(x)any(strcmp(x,{'en','zh'})));
catalog=aple_figure_catalog;
addParameter(p,'Figures',{catalog.code});
addParameter(p,'AllowPartial',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'SaveFig',false,@(x)islogical(x)&&isscalar(x));
parse(p,varargin{:}); o=p.Results;
has_checkpoints=~isempty(dir(fullfile(source,'*.mat')));
if has_checkpoints || isfolder(fullfile(source,'.aple-lock'))
    guard=aple_lock(source); %#ok<NASGU>
    if has_checkpoints, aple_export_results(source); end
end
meta=jsondecode(fileread(fullfile(source,'metadata.json')));
assert(o.AllowPartial || (meta.complete && meta.failed_trials==0),'APLE:IncompleteResults', ...
    'Run is incomplete or contains failed estimates. Resume it, or explicitly use AllowPartial=true for inspection.');
t=readtable(fullfile(source,'summary.csv'),'TextType','string');
cdf=readtable(fullfile(source,'cdf.csv'),'TextType','string');
clear guard; % The plot now owns a consistent in-memory snapshot.
output=char(o.OutputDir);
if isempty(output)
    if strcmp(source,fullfile(root,'reference')), output=fullfile(root,'runs','reference_plots');
    else, output=fullfile(source,'figures'); end
end
if ~isfolder(output), mkdir(output); end
figures=upper(cellstr(string(o.Figures))); audit=cell(1,numel(figures));
for k=1:numel(figures)
    audit{k}=aple_plot_main(figures{k},t,cdf,output,o);
end
fid=fopen(fullfile(output,'plot_manifest.json'),'w','n','UTF-8'); assert(fid>=0);
clean=onCleanup(@()fclose(fid)); fprintf(fid,'%s',jsonencode(struct('dataset',meta,'figures',{audit}),'PrettyPrint',true));
fprintf('Plots: %s\n',output);
end
