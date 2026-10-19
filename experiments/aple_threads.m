function guard=aple_threads
% Match the frozen experiment's BLAS thread policy without changing the caller.
previous=maxNumCompThreads(1);
guard=onCleanup(@()maxNumCompThreads(previous));
end
