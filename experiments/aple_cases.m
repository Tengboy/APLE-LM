function d=aple_cases(family)
% Canonical protocol, including the historical separation rule in gene_pos.
family=lower(char(family));
c=struct('M',3,'Ns',15,'RF',160,'W_seed',0,'grid',60,'range',[5,10], ...
    'users',3,'theta',NaN,'radius',NaN,'key','');
d=struct('family',family,'points',5:3:20,'methods',{{'APLE-LM'}},'cases',c,'smoke_cases',1);
switch family
    case 'snr'
        d.methods={'APLE-LM','APLE-LM-ACM','ES-GA','P-SIGW'};
        c.key='A45_M3_S15'; d.cases=c;
    case 'grid'
        d.points=15; d.methods={'APLE-LM','ES-GA','P-SIGW'}; d.cases=repmat(c,1,12);
        for g=1:2
            for j=1:6
                k=(g-1)*6+j; v=c; v.M=2*g+1; v.grid=15*j;
                v.key=sprintf('A%d_M%d_S15_grid%d',15*v.M,v.M,v.grid); d.cases(k)=v;
            end
        end
        d.smoke_cases=[1,6,7,12];
    case 'aperture'
        configs=[3,15;4,15;5,15;3,20;3,25]; d.cases=repmat(c,1,5);
        for k=1:5
            v=c; v.M=configs(k,1); v.Ns=configs(k,2); v.grid=4*v.Ns;
            v.key=sprintf('A%d_M%d_S%d',v.M*v.Ns,v.M,v.Ns); d.cases(k)=v;
        end
        d.smoke_cases=1:5;
    case {'range','iterations'}
        count=3; if strcmp(family,'iterations'), count=2; d.points=1:20; end
        d.cases=repmat(c,1,count);
        for k=1:count, v=c; v.range=[5*k,5*k+5]; v.key=sprintf('range_%d_%d',v.range); d.cases(k)=v; end
        d.smoke_cases=1:count;
    case 'rf'
        configs=[3,15;5,9;4,15]; d.points=15; d.cases=repmat(c,1,15);
        for g=1:3
            for j=1:5
                v=c; v.M=configs(g,1); v.Ns=configs(g,2); v.RF=30+30*j; v.W_seed=j; v.grid=4*v.Ns;
                v.key=sprintf('A%d_M%d_S%d_RF%d',v.M*v.Ns,v.M,v.Ns,v.RF); d.cases((g-1)*5+j)=v;
            end
        end
        d.smoke_cases=[1,5,6,10,11,15];
    case 'coverage'
        d.points=15; [angle,radius]=ndgrid(0:5:60,5:15); d.cases=repmat(c,1,numel(angle));
        for k=1:numel(angle)
            v=c; v.users=1; v.range=[5,15]; v.theta=angle(k); v.radius=radius(k);
            v.key=sprintf('theta%03d_r%02d',v.theta,v.radius); d.cases(k)=v;
        end
        d.smoke_cases=[1,13,72,131,143];
    otherwise
        error('APLE:UnknownExperiment','Unknown experiment family: %s',family);
end
end
