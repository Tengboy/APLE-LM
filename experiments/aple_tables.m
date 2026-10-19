function [summary,cdf_table,meta]=aple_tables(states)
% Single statistical definition for both bundled references and fresh runs.
rows={}; cdf_rows={}; total=0; complete=0; failures=0; families={}; reference_sized=true;
for f=1:numel(states)
    st=states{f}; d=st.definition; families{end+1}=d.family; %#ok<AGROW>
    n=numel(st.seeds); reference_sized=reference_sized && isequal(st.seeds,1:1000);
    total=total+numel(st.done); complete=complete+nnz(st.done); failures=failures+nnz(st.failed);
    for c=1:numel(d.cases)
        mask=st.done(:,c); if ~any(mask), continue; end
        v=d.cases(c);
        for p=1:numel(d.points)
            b=st.bound(mask,p,c); bound=sqrt(mean(b));
            for m=1:numel(d.methods)
                values=st.mse(mask,m,p,c); s=aple_metric(values,3);
                nmse=NaN; cb=NaN;
                if ~isempty(st.nmse)
                    nmse=10*log10(mean(st.nmse(mask,m,p,c))); cb=10*log10(mean(st.channel_bound(mask,p,c)));
                end
                rows{end+1}=struct('Experiment',d.family,'Condition',v.key,'Method',d.methods{m}, ...
                    'Point',d.points(p),'M',v.M,'Ns',v.Ns,'RF',v.RF,'Grid',v.grid, ...
                    'RangeMin',v.range(1),'RangeMax',v.range(2),'Theta',v.theta,'Radius',v.radius, ...
                    'Trials',n,'Attempted',s.attempted,'Valid',s.valid,'Included',s.included,'Excluded',s.exceeded, ...
                    'RMSE_m',s.rmse,'Conditional_RMSE_m',s.conditional_rmse,'Exceedance_probability',s.outage, ...
                    'Wilson95_low',s.wilson_low,'Wilson95_high',s.wilson_high,'BCRB_m',bound, ...
                    'Channel_NMSE_dB',nmse,'Channel_BCRB_dB',cb); %#ok<AGROW>
                if strcmp(d.family,'snr') && ismember(d.points(p),[5,20])
                    e=sort(sqrt(values(isfinite(values) & values>=0))); count=numel(e);
                    if count>0
                        cdf_rows{end+1}=table(repmat(d.points(p),count+1,1), ...
                            repmat(string(d.methods{m}),count+1,1),[e(1);e(:)],(0:count).'/s.attempted, ...
                            'VariableNames',{'SNR_dB','Method','Error_m','Probability'}); %#ok<AGROW>
                    end
                end
            end
        end
    end
end
if isempty(rows), summary=table; else, summary=struct2table([rows{:}]); end
if isempty(cdf_rows)
    cdf_table=table(zeros(0,1),strings(0,1),zeros(0,1),zeros(0,1),'VariableNames',{'SNR_dB','Method','Error_m','Probability'});
else, cdf_table=vertcat(cdf_rows{:}); end
meta=struct('schema',1,'families',{families},'expected_trials',total,'completed_trials',complete, ...
    'failed_trials',failures,'complete',total>0&&complete==total,'reference_sized',~isempty(states)&&reference_sized, ...
    'threshold_m',3,'iteration_display_limit',10,'raw_iterations',20, ...
    'BCRB_definition','Label follows the paper; full-sample mean of per-trial bounds, no explicit prior FIM', ...
    'position_matching','one-to-one minimum total squared distance', ...
    'user_separation','historical 2/Ns per axis, OR acceptance rule; not fixed across all groups');
end
