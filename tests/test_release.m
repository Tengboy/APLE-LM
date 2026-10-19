function tests=test_release
tests=functiontests(localfunctions);
end
function setupOnce(t)
root=fileparts(fileparts(mfilename('fullpath')));
addpath(root,fullfile(root,'src'),fullfile(root,'experiments'),fullfile(root,'plotting'));
t.TestData.root=root;
end
function testPublicAPIExists(t)
verifyEqual(t,exist('run_experiments','file'),2,'The independent experiment entry point is missing.');
verifyEqual(t,exist('plot_results','file'),2,'The independent plotting entry point is missing.');
end
function testCanonicalCases(t)
c=aple_cases('rf');
verifyEqual(t,numel(c.cases),15);
verifyEqual(t,[c.cases(6).M,c.cases(6).Ns,c.cases(6).RF,c.cases(6).W_seed],[5,9,60,1]);
verifyEqual(t,c.cases(6).grid,36);
a=aple_cases('aperture');
verifyEqual(t,[[a.cases.M].',[a.cases.Ns].'],[3,15;4,15;5,15;3,20;3,25]);
i=aple_cases('iterations'); verifyEqual(t,i.points,1:20);
h=aple_cases('coverage'); verifyEqual(t,numel(h.cases),143);
end
function testCheckpointIdentityAndRecovery(t)
folder=tempname; mkdir(folder); clean=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
file=fullfile(folder,'state.mat');
s=struct('signature','identity-a','generation',1,'done',[true;false],'mse',[.25;NaN]);
aple_checkpoint('save',file,s);
verifyEqual(t,aple_checkpoint('load',file,'identity-a'),s);
verifyError(t,@()aple_checkpoint('load',file,'identity-b'),'APLE:CheckpointIdentity');
s.generation=2; s.done(:)=true; s.mse(2)=.5; aple_checkpoint('save',file,s);
fid=fopen(file,'w'); fwrite(fid,'interrupted'); fclose(fid);
warning('off','APLE:CheckpointRecovered'); cw=onCleanup(@()warning('on','APLE:CheckpointRecovered')); %#ok<NASGU>
old=aple_checkpoint('load',file,'identity-a');
verifyEqual(t,old.generation,1); verifyEqual(t,old.done,[true;false]);
% Saving after recovery must not replace the valid backup with corrupt bytes.
old.generation=3; aple_checkpoint('save',file,old);
fid=fopen(file,'w'); fwrite(fid,'interrupted again'); fclose(fid);
recovered=aple_checkpoint('load',file,'identity-a');
verifyEqual(t,recovered.generation,1);
verifyEqual(t,aple_checkpoint('load',file,''),recovered);
verifyLessThanOrEqual(t,numel(dir(fullfile(folder,'*.mat'))),3);
end
function testConditionalStatisticsKeepFailures(t)
s=aple_metric([1,4,16,NaN],3);
verifyEqual(t,s.attempted,4); verifyEqual(t,s.valid,3);
verifyEqual(t,s.included,2); verifyEqual(t,s.exceeded,2);
verifyEqual(t,s.conditional_rmse,sqrt(2.5),'AbsTol',1e-12);
verifyEqual(t,s.outage,.5); verifyTrue(t,isnan(s.rmse));
verifyLessThan(t,s.wilson_low,.5); verifyGreaterThan(t,s.wilson_high,.5);
end
function testExclusiveRunLock(t)
folder=tempname; mkdir(folder); clean=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
guard=aple_lock(folder);
verifyError(t,@()aple_lock(folder),'APLE:RunLocked');
clear guard;
guard=aple_lock(folder); clear guard;
end
function testThreadPolicyReproducesStressCase(t)
% Compare the same low-RF problem under two caller thread settings. Kernels
% must run with the historical single-thread policy, then restore the caller.
old=maxNumCompThreads; clean=onCleanup(@()maxNumCompThreads(old)); %#ok<NASGU>
d=aple_cases('rf');
maxNumCompThreads(1); b1=aple_system(d,6); evalc('r1=aple_trial(1,d,6,b1);');
verifyEqual(t,maxNumCompThreads,1);
maxNumCompThreads(2); b2=aple_system(d,6); evalc('r2=aple_trial(1,d,6,b2);');
verifyEqual(t,maxNumCompThreads,2);
verifyEqual(t,b2.accurate.noise_precision,b1.accurate.noise_precision);
verifyEqual(t,r2.mse,r1.mse);
verifyTrue(t,all(isfinite(r1.mse),'all'));
end
function testPartialIterationPlot(t)
folder=iteration_fixture; clean=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
out=plot_results(folder,'Figures',{'M8'},'AllowPartial',true,'SaveFig',true);
verifyTrue(t,isfile(fullfile(out,'M8_iteration_behavior.pdf')));
fig=openfig(fullfile(out,'M8_iteration_behavior.fig'),'invisible'); cf=onCleanup(@()close(fig)); %#ok<NASGU>
axes_handles=findall(fig,'Type','axes');
ax=axes_handles(arrayfun(@(x)x.Layout.Tile==2,axes_handles));
labels=findall(ax,'Type','text'); messages={labels.String};
verifyTrue(t,any(cellfun(@(x)isequal(x,'No valid trials'),messages)));
end
function testPlotRespectsRunLock(t)
folder=iteration_fixture; clean=onCleanup(@()rmdir(folder,'s')); %#ok<NASGU>
guard=aple_lock(folder); %#ok<NASGU>
verifyError(t,@()plot_results(folder,'Figures',{'M8'},'AllowPartial',true),'APLE:RunLocked');
end
function folder=iteration_fixture
folder=tempname; mkdir(folder); d=aple_cases('iterations');
s=struct('signature','test-fixture','definition',d,'seeds',1,'done',[true,false], ...
    'failed',[false,false],'mse',nan(1,1,20,2),'bound',ones(1,20,2),'nmse',[],'channel_bound',[]);
s.mse(1,1,:,1)=linspace(.04,.02,20);
aple_checkpoint('save',fullfile(folder,'iterations.mat'),s);
end
