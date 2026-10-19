function guard=aple_lock(folder)
% Atomic directory lock: separate MATLAB processes must not write one run.
file=fullfile(folder,'.aple-lock'); lock=java.io.File(file);
assert(lock.mkdir(),'APLE:RunLocked', ...
    'This run is locked. If its previous MATLAB process has stopped, remove only %s and resume.',file);
guard=onCleanup(@()release(file));
end
function release(file)
if isfolder(file), rmdir(file); end
end
