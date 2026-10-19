function root=setup_aple
% Add only this release's directories, without recursively adding run data.
% Copyright (c) 2026 Boyu Teng. MIT License; see LICENSE.
root=fileparts(mfilename('fullpath'));
addpath(root,fullfile(root,'src'),fullfile(root,'experiments'),fullfile(root,'plotting'),'-begin');
assert(startsWith(which('pos_ini_est'),[root filesep]),'APLE:PathConflict', ...
    'A different APLE-LM copy shadows this release. Start a clean MATLAB session.');
end
