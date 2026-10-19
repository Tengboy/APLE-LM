function result=aple_checkpoint(action,file,value)
% One small current checkpoint and one previous checkpoint per experiment.
% A failed write never replaces the last complete checkpoint with partial data.
previous=[file '.previous.mat']; pending=[file '.pending.mat']; result=[];
switch action
    case 'save'
        % Never overwrite a valid recovery copy with a corrupt current file.
        [current,valid]=read_state(file);
        if valid
            assert(strcmp(current.signature,value.signature),'APLE:CheckpointIdentity', ...
                'Refusing to overwrite a checkpoint with a different identity.');
        end
        state=value; %#ok<NASGU>
        save(pending,'state','-v7');
        check=load(pending,'state');
        assert(isequaln(check.state,value),'APLE:CheckpointWrite','Checkpoint verification failed.');
        if valid, copyfile(file,previous,'f'); end
        [ok,message]=movefile(pending,file,'f'); assert(ok,'APLE:CheckpointWrite','%s',message);
    case 'load'
        candidates={file,previous}; recovered=false;
        for k=1:2
            if ~isfile(candidates{k}), continue; end
            [saved,valid]=read_state(candidates{k});
            if ~valid, recovered=true; continue; end
            assert(isempty(value) || strcmp(saved.signature,value),'APLE:CheckpointIdentity', ...
                'Configuration, seeds, MATLAB release, or numerical source changed. Use a different output directory.');
            result=saved;
            if recovered || k==2, warning('APLE:CheckpointRecovered','Using previous complete checkpoint: %s',candidates{k}); end
            return;
        end
        assert(~isfile(file) && ~isfile(previous) && ~isfile(pending),'APLE:CheckpointCorrupt', ...
            'No readable complete checkpoint. Preserve this folder and use a new output directory.');
    otherwise
        error('APLE:CheckpointAction','Unknown checkpoint action.');
end
end
function [state,valid]=read_state(file)
state=[]; valid=false;
if ~isfile(file), return; end
try
    saved=load(file,'state');
    valid=isfield(saved,'state') && isstruct(saved.state) && isfield(saved.state,'signature');
    if valid, state=saved.state; end
catch
    % A previous completed checkpoint may still be usable.
end
end
