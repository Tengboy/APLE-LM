function output=run_experiments(figures,profile,varargin)
%RUN_EXPERIMENTS Reproduce M1-M8 without historical data files.
% run_experiments('all','smoke') checks selected cases with one seed.
% run_experiments('all','full') runs every case with seeds 1:1000.
% Options: OutputDir, Trials, Seeds, Cases, BatchSize, MaxBatches,
%          UseParallel, Workers. See docs/REPRODUCIBILITY.md.
root=setup_aple;
if nargin<1, figures='all'; end
if nargin<2, profile='smoke'; end
profile=validatestring(profile,{'smoke','full'});
p=inputParser;
addParameter(p,'OutputDir',fullfile(root,'runs',profile),@(x)ischar(x)||isstring(x));
addParameter(p,'Trials',[],@(x)isempty(x)||(isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x));
addParameter(p,'Seeds',[],@(x)isnumeric(x)&&isvector(x)&&all(isfinite(x)&x>=0&x<=2^32-1&fix(x)==x)||isempty(x));
addParameter(p,'Cases',[],@(x)isnumeric(x)&&(isempty(x)||(isvector(x)&&all(isfinite(x)&x>=1&fix(x)==x))));
addParameter(p,'BatchSize',25,@(x)isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x);
addParameter(p,'MaxBatches',Inf,@(x)isscalar(x)&&x>=0&&(isinf(x)||fix(x)==x));
addParameter(p,'UseParallel',false,@(x)islogical(x)&&isscalar(x));
addParameter(p,'Workers',2,@(x)isscalar(x)&&isfinite(x)&&x>=1&&fix(x)==x);
parse(p,varargin{:}); o=p.Results;
if isempty(o.Trials), o.Trials=1000; if strcmp(profile,'smoke'), o.Trials=1; end, end
seeds=o.Seeds; if isempty(seeds), seeds=1:o.Trials; end
seeds=reshape(seeds,1,[]); assert(numel(unique(seeds))==numel(seeds),'APLE:DuplicateSeeds','Seeds must be unique.');
families=select_families(figures); output=char(o.OutputDir);
if ~isfolder(output), mkdir(output); end
output=char(java.io.File(output).getCanonicalPath()); guard=aple_lock(output); %#ok<NASGU>
if o.UseParallel
    assert(license('test','Distrib_Computing_Toolbox'),'APLE:ParallelUnavailable','Parallel Computing Toolbox is required for UseParallel=true.');
    pool=gcp('nocreate');
    if isempty(pool), pool=parpool('local',o.Workers); end
    assert(pool.NumWorkers==o.Workers,'APLE:PoolSize','Existing pool size differs from Workers; use its actual size or close the pool.');
end
batches=0;
for f=1:numel(families)
    d=aple_cases(families{f}); cid=o.Cases;
    if isempty(cid), cid=1:numel(d.cases); if strcmp(profile,'smoke'), cid=d.smoke_cases; end, end
    assert(all(cid<=numel(d.cases)),'APLE:CaseIndex','Case index exceeds %s case count.',d.family);
    cid=unique(cid,'stable'); identity=aple_signature(d,seeds);
    file=fullfile(output,[d.family '.mat']); state=aple_checkpoint('load',file,identity);
    if isempty(state), state=new_state(d,seeds,identity,profile); end
    bundle=[]; active_key='';
    for k=cid
        if batches>=o.MaxBatches, break; end
        pending=find(~state.done(:,k));
        if isempty(pending), continue; end
        c=d.cases(k); key=jsonencode(rmfield(c,{'theta','radius','key'}));
        if ~strcmp(key,active_key)
            clear bundle;
            fprintf('Building %s / %s in memory...\n',d.family,c.key);
            bundle=aple_system(d,k); active_key=key;
        end
        for offset=1:o.BatchSize:numel(pending)
            if batches>=o.MaxBatches, break; end
            index=pending(offset:min(offset+o.BatchSize-1,numel(pending))); records=cell(numel(index),1);
            if o.UseParallel
                parfor j=1:numel(index), records{j}=safe_trial(seeds(index(j)),d,k,bundle); end
            else
                for j=1:numel(index), records{j}=safe_trial(seeds(index(j)),d,k,bundle); end
            end
            for j=1:numel(index)
                row=index(j); r=records{j};
                state.mse(row,:,:,k)=reshape(r.mse,1,numel(d.methods),numel(d.points));
                state.bound(row,:,k)=r.bound;
                state.iterations(row,:,:,k)=reshape(r.iterations,1,numel(d.methods),numel(d.points));
                if strcmp(d.family,'snr')
                    state.nmse(row,:,:,k)=reshape(r.nmse,1,4,6); state.channel_bound(row,:,k)=r.channel_bound;
                end
                state.seconds(row,k)=single(r.seconds); state.done(row,k)=true;
                bad=any(~isfinite(r.mse),'all') || (strcmp(d.family,'snr')&&any(~isfinite(r.nmse),'all'));
                state.failed(row,k)=bad;
                if bad && numel(state.failure_examples)<10
                    state.failure_examples{end+1}=struct('case',k,'seed',seeds(row),'message',r.error); %#ok<AGROW>
                end
            end
            state.generation=state.generation+1;
            aple_checkpoint('save',file,state); batches=batches+1;
            fprintf('%s case %d/%d: %d/%d trials committed, %d failed.\n', ...
                d.family,k,numel(d.cases),nnz(state.done(:,k)),numel(seeds),nnz(state.failed(:,k)));
        end
        if batches>=o.MaxBatches, break; end
    end
    if batches>=o.MaxBatches, break; end
end
aple_export_results(output);
fprintf('Results: %s\n',output);
end
function state=new_state(d,seeds,identity,profile)
n=numel(seeds); m=numel(d.methods); p=numel(d.points); c=numel(d.cases);
state=struct('schema',1,'signature',identity,'profile',profile,'definition',d,'seeds',seeds, ...
    'matlab_release',version('-release'),'generation',0,'done',false(n,c),'failed',false(n,c), ...
    'mse',nan(n,m,p,c),'bound',nan(n,p,c),'nmse',[],'channel_bound',[], ...
    'iterations',zeros(n,m,p,c,'uint16'),'seconds',zeros(n,c,'single'),'failure_examples',{{}});
if strcmp(d.family,'snr'), state.nmse=nan(n,m,p,c); state.channel_bound=nan(n,p,c); end
end
function r=safe_trial(seed,d,cid,bundle)
timer=tic; r=[]; caught=[];
% Printed warnings are transient; failed jobs retain a bounded short message.
evalc('try; r=aple_trial(seed,d,cid,bundle); catch err; caught=err; end');
if ~isempty(caught)
    r=struct('mse',nan(numel(d.methods),numel(d.points)),'bound',nan(1,numel(d.points)), ...
        'nmse',nan(numel(d.methods),numel(d.points)),'channel_bound',nan(1,numel(d.points)), ...
        'iterations',zeros(numel(d.methods),numel(d.points),'uint16'),'error',[caught.identifier ': ' caught.message]);
    r.error=r.error(1:min(300,numel(r.error)));
end
r.seconds=toc(timer);
end
function families=select_families(figures)
names=upper(cellstr(string(figures))); families={};
catalog=aple_figure_catalog;
for k=1:numel(names)
    index=find(strcmpi(names{k},{catalog.code})|strcmpi(names{k},{catalog.key}));
    if strcmp(names{k},'ALL'), add={catalog.family};
    elseif isscalar(index), add={catalog(index).family};
    else, error('APLE:UnknownFigure','Unknown figure/family: %s',names{k});
    end
    families=[families,add]; %#ok<AGROW>
end
families=unique(families,'stable');
end
