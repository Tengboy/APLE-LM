function value=aple_signature(d,seeds)
% Portable identity: no absolute paths; plotting changes do not invalidate runs.
root=fileparts(fileparts(mfilename('fullpath')));
files=[dir(fullfile(root,'src','*.m'));dir(fullfile(root,'experiments','*.m'));dir(fullfile(root,'run_experiments.m'))];
md=java.security.MessageDigest.getInstance('SHA-256');
md.update(uint8(unicode2native(jsonencode(struct('protocol',d,'seeds',seeds,'release',version('-release'))),'UTF-8')));
[~,order]=sort({files.name}); files=files(order);
for k=1:numel(files)
    fid=fopen(fullfile(files(k).folder,files(k).name),'r'); assert(fid>=0);
    bytes=fread(fid,Inf,'*uint8'); fclose(fid); md.update(uint8(files(k).name)); md.update(bytes);
end
value=lower(reshape(dec2hex(typecast(md.digest(),'uint8'),2).',1,[]));
end
