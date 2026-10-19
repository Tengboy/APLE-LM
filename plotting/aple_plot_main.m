function audit=aple_plot_main(code,t,cdf,out,o)
% Figure styling and statistics are separated: this file never clips samples.
catalog=aple_figure_catalog;
index=find(strcmpi(code,{catalog.code})|strcmpi(code,{catalog.key}));
assert(isscalar(index),'APLE:UnknownFigure','Unknown figure: %s',code);
item=catalog(index); code=item.code; topic=item.key;
C=[0,.36,.70;.18,.56,.34;.84,.32,.10;.53,.25,.67];
styles={'-o','--d','-.s',':^'}; names={'APLE-LM','APLE-LM-ACM','ES-GA','P-SIGW'};
font='Arial'; zh=strcmp(o.Language,'zh');
if zh
    fonts=listfonts;
    preferred={'SimSun','宋体','Noto Sans CJK SC','Microsoft YaHei','微软雅黑','PingFang SC','WenQuanYi Micro Hei'};
    found=find(ismember(preferred,fonts),1);
    assert(~isempty(found),'APLE:ChineseFont','Install a Chinese font or use Language=en.');
    font=preferred{found};
end
snr=5:3:20; rf=60:30:180; cn=C([1,3,2],:);
switch topic
    case 'snr'
        require(t,{'snr'}); [fig,tiles]=canvas(1,2,19,8.7);
        ax=nexttile(tiles); hold(ax,'on'); h=gobjects(1,5);
        for m=1:4, r=getrows(t,'snr',names{m}); h(m)=lineplot(ax,r.Point,r.RMSE_m,C(m,:),styles{m}); end
        r=getrows(t,'snr','APLE-LM'); h(5)=plot(ax,r.Point,r.BCRB_m,'k--');
        set(ax,'YScale','log','XTick',snr); labels(ax,'接收信噪比 / dB','Received SNR / dB','位置 RMSE / m','Position RMSE / m',1);
        ax=nexttile(tiles); hold(ax,'on');
        for m=1:4, r=getrows(t,'snr',names{m}); h(m)=lineplot(ax,r.Point,r.Channel_NMSE_dB,C(m,:),styles{m}); end
        r=getrows(t,'snr','APLE-LM'); h(5)=plot(ax,r.Point,r.Channel_BCRB_dB,'k--');
        set(ax,'XTick',snr); labels(ax,'接收信噪比 / dB','Received SNR / dB','等效信道 NMSE / dB','Effective-channel NMSE / dB',2);
        legend(h,[names,{'BCRB'}],'Location','southwest');
    case 'cdf'
        require(t,{'snr'}); [fig,tiles]=canvas(1,2,19,8.7);
        finite=cdf.Error_m(isfinite(cdf.Error_m)&cdf.Error_m>0); lim=[min(.006,min(finite)/1.08),max(20,max(finite)*1.08)];
        for p=1:2
            vals=[5,20]; ax=nexttile(tiles); hold(ax,'on'); h=gobjects(1,4);
            for m=1:4
                r=cdf(cdf.SNR_dB==vals(p)&cdf.Method==names{m},:);
                h(m)=stairs(ax,r.Error_m,r.Probability,'Color',C(m,:),'LineStyle',styles{m}(1:end-1),'LineWidth',1.35);
            end
            set(ax,'XScale','log'); xlim(ax,lim); ylim(ax,[0,1]);
            labels(ax,'单次试验位置 RMS / m','Trial position RMS / m','累积概率','Cumulative probability',p);
            legend(h,names,'Location','southeast');
        end
    case 'grid'
        require(t,{'grid'}); [fig,tiles]=canvas(1,2,19,10.4); ids=[1,3,4];
        % Like journal Fig. 4: aperture identifies color/marker, method identifies line style.
        colors=[.80,.16,.12;0,.36,.70]; marks={'o','v'}; ls={'-','-.','--'}; bs={'--','-.'};
        h=gobjects(1,8); legend_text=cell(1,8);
        ax=nexttile(tiles,1); hold(ax,'on');
        for g=1:2
            side=30*g+15;
            for m=1:3
                r=getrows(t,'grid',names{ids(m)}); r=sortrows(r(r.M==2*g+1,:),'Grid');
                j=(g-1)*4+m;
                h(j)=lineplot(ax,r.Grid,r.Conditional_RMSE_m,colors(g,:),[ls{m},marks{g}]);
                legend_text{j}=sprintf('%s, %d×%d',names{ids(m)},side,side);
            end
            r=getrows(t,'grid','APLE-LM'); r=r(r.M==2*g+1,:); b=r.BCRB_m(isfinite(r.BCRB_m));
            if isempty(b), value=NaN; else
                assert(max(abs(b-b(1)))<1e-10*max(1,abs(b(1))),'APLE:NonconstantGridBound', ...
                    'BCRB varies with the grid; a horizontal reference would misrepresent this data.');
                value=b(1);
            end
            h(g*4)=plot(ax,[15,90],value*[1,1],bs{g},'Color',[.1,.1,.1],'LineWidth',1.15);
            legend_text{g*4}=sprintf('BCRB, %d×%d',side,side);
        end
        set(ax,'YScale','log','XTick',15:15:90); padlog(ax);
        labels(ax,'名义角网格参数 M_x','Nominal angular-grid parameter M_x','条件位置 RMSE / m','Conditional position RMSE / m',1);
        order=[1,5,2,6,3,7,4,8]; lg=legend(ax,h(order),legend_text(order),'NumColumns',4,'Interpreter','none');
        lg.Layout.Tile='south';
        ax=nexttile(tiles,2); hold(ax,'on');
        for g=1:2
            for m=1:3
                r=getrows(t,'grid',names{ids(m)}); r=sortrows(r(r.M==2*g+1,:),'Grid');
                probability(ax,r.Grid,r,colors(g,:),[ls{m},marks{g}]);
            end
        end
        set(ax,'XTick',15:15:90); ylim(ax,[-1,102]);
        labels(ax,'名义角网格参数 M_x','Nominal angular-grid parameter M_x','超限概率 P(E > 3 m) / %','Exceedance probability P(E > 3 m) / %',2);
    case 'aperture'
        require(t,{'aperture'}); [fig,tiles]=canvas(1,1,12.5,10.5);
        ax=nexttile(tiles); hold(ax,'on'); h=gobjects(1,6); a=[3,15;4,15;5,15;3,20;3,25]; ci=[1,2,3,2,3]; st={'-o','-s','-^','--d','--v'};
        for g=1:5
            r=getrows(t,'aperture','APLE-LM'); r=r(r.M==a(g,1)&r.Ns==a(g,2),:);
            h(g)=lineplot(ax,r.Point,r.Conditional_RMSE_m,cn(ci(g),:),st{g});
            if g<=3, plot(ax,r.Point,r.BCRB_m,':','Color',cn(g,:),'HandleVisibility','off'); end
        end
        h(6)=plot(ax,nan,nan,'k:');
        legend(h,{'M=3, N_s=15','M=4, N_s=15','M=5, N_s=15','M=3, N_s=20','M=3, N_s=25',word('BCRB（同色点线）','BCRB (matching colors)')}, ...
            'Interpreter','tex','Location','northeast','NumColumns',2);
        set(ax,'YScale','log','XTick',snr); padlog(ax); lim=ylim(ax); ylim(ax,[min(.026,lim(1)),max(.85,1.5*lim(2))]);
        labels(ax,'接收信噪比 / dB','Received SNR / dB','条件位置 RMSE / m','Conditional position RMSE / m',0);
    case 'range'
        require(t,{'range'}); [fig,tiles]=canvas(1,1,12.5,10.5);
        ax=nexttile(tiles); hold(ax,'on'); h=gobjects(1,4); st={'-o','--s','-.^'};
        for g=1:3
            r=getrows(t,'range','APLE-LM'); r=r(r.RangeMin==5*g,:);
            h(g)=lineplot(ax,r.Point,r.RMSE_m,cn(g,:),st{g}); plot(ax,r.Point,r.BCRB_m,':','Color',cn(g,:));
        end
        h(4)=plot(ax,nan,nan,'k:'); legend(h,{'5-10 m','10-15 m','15-20 m',word('BCRB（同色点线）','BCRB (matching colors)')},'Location','northeast');
        set(ax,'YScale','log','XTick',snr); padlog(ax); lim=ylim(ax); ylim(ax,[lim(1),1.4*lim(2)]);
        labels(ax,'接收信噪比 / dB','Received SNR / dB','位置 RMSE / m','Position RMSE / m',0);
    case 'rf'
        require(t,{'rf'}); [fig,tiles]=canvas(1,2,19,8.7); ax=nexttile(tiles); hold(ax,'on');
        h=gobjects(1,5); a=[3,15;5,9;4,15]; colors=C([1,4,3],:); st={'-o','-d','-s'};
        lg={'APLE-LM, M=3, N_s=15','APLE-LM, M=5, N_s=9','APLE-LM, M=4, N_s=15','BCRB, M=3, N_s=15','BCRB, M=4, N_s=15'};
        for g=1:3
            r=getrows(t,'rf','APLE-LM'); r=sortrows(r(r.M==a(g,1)&r.Ns==a(g,2),:),'RF');
            h(g)=lineplot(ax,r.RF,r.Conditional_RMSE_m,colors(g,:),st{g});
            if g~=2, b=1+(g==3); h(3+b)=plot(ax,r.RF,r.BCRB_m,'--','Color',colors(g,:)); end
        end
        set(ax,'YScale','log','XTick',rf); padlog(ax); legend(h,lg,'Interpreter','tex','Location','northeast');
        labels(ax,'RF 链数量','Number of RF chains','条件位置 RMSE / m','Conditional position RMSE / m',1);
        ax=nexttile(tiles); hold(ax,'on');
        for g=1:3
            r=getrows(t,'rf','APLE-LM'); r=sortrows(r(r.M==a(g,1)&r.Ns==a(g,2),:),'RF');
            h(g)=probability(ax,r.RF,r,colors(g,:),st{g});
        end
        set(ax,'XTick',rf); ylim(ax,[-1,102]); legend(h(1:3),lg(1:3),'Interpreter','tex','Location','northeast');
        labels(ax,'RF 链数量','Number of RF chains','超限概率 P(E > 3 m) / %','Exceedance probability P(E > 3 m) / %',2);
    case 'iterations'
        require(t,{'iterations'}); [fig,tiles]=canvas(2,1,11.5,11.5);
        for g=1:2
            ax=nexttile(tiles); hold(ax,'on'); r=getrows(t,'iterations','APLE-LM'); r=r(r.RangeMin==5*g & r.Point<=10,:);
            h=lineplot(ax,r.Point,r.RMSE_m,C(1,:),'-o'); set(ax,'XLim',[1,10],'XTick',1:10);
            values=r.RMSE_m(isfinite(r.RMSE_m));
            if isempty(values)
                ylim(ax,[0,1]);
                text(ax,.5,.5,word('暂无有效结果','No valid trials'),'Units','normalized', ...
                    'HorizontalAlignment','center','FontName',font,'Interpreter','none');
            else
                padding=max(1e-8,.12*(max(values)-min(values)));
                ylim(ax,[min(values)-padding,max(values)+padding]);
            end
            lg=sprintf(word('APLE-LM，距离 %d-%d m','APLE-LM, range %d-%d m'),5*g,5*g+5);
            if ~isempty(h), legend(h,{lg},'Location','northeast'); end
            labels(ax,'外层迭代次数','Outer iteration','位置 RMSE / m','Position RMSE / m',g);
        end
    case 'coverage'
        require(t,{'coverage'}); [fig,tiles]=canvas(1,2,19,8.7); r=getrows(t,'coverage','APLE-LM');
        maps=nan(11,13,2);
        for k=1:height(r), maps(r.Radius(k)-4,r.Theta(k)/5+1,1)=r.RMSE_m(k); maps(r.Radius(k)-4,r.Theta(k)/5+1,2)=r.BCRB_m(k); end
        v=maps(isfinite(maps)&maps>0); lim=[min(v),max(v)]; if lim(1)==lim(2), lim=lim.*[.9,1.1]; end
        for p=1:2
            ax=nexttile(tiles); imagesc(ax,0:5:60,5:15,maps(:,:,p)); axis(ax,'xy');
            set(ax,'ColorScale','log','XTick',0:15:60,'YTick',5:2:15); caxis(ax,lim); colormap(ax,parula(256));
            cb=colorbar(ax); cb.Label.String=word('误差 / m','Error / m'); cb.Label.FontName=font;
            cb.Ticks=[.04,.08,.16,.32]; cb.TickLabels={'0.04','0.08','0.16','0.32'};
            labels(ax,'扫描角度 / °','Scan angle / deg','用户距离 / m','User distance / m',p);
        end
    otherwise, error('APLE:UnknownFigure','Unknown figure: %s',code);
end
filename=item.filename;
axes_handles=findall(fig,'Type','axes');
for k=1:numel(axes_handles)
    ax=axes_handles(k); ax.Title.String=''; ax.FontName='Times New Roman'; ax.FontSize=10.5;
    ax.XLabel.FontName=font; ax.YLabel.FontName=font;
    ax.XLabel.Interpreter='tex'; ax.YLabel.Interpreter='none'; ax.LineWidth=.65;
    box(ax,'on'); if ~strcmp(topic,'coverage'), grid(ax,'on'); end
    ax.GridAlpha=.16; ax.Layer='top'; ax.TickDir='out';
end
legends=findall(fig,'Type','legend');
for k=1:numel(legends)
    legends(k).FontName=font; legends(k).FontSize=9; legends(k).Box='on';
    legends(k).Color=[1,1,1]; legends(k).EdgeColor=[.2,.2,.2]; legends(k).LineWidth=.65;
    legends(k).ItemTokenSize=[22,10];
end
drawnow;
exportgraphics(fig,fullfile(out,[filename '.pdf']),'ContentType','vector','BackgroundColor','white');
exportgraphics(fig,fullfile(out,[filename '.png']),'Resolution',240,'BackgroundColor','white');
if o.SaveFig, savefig(fig,fullfile(out,[filename '.fig'])); end
audit=struct('code',code,'key',topic,'previous_code',item.previous_code,'filename',filename, ...
    'size_cm',fig.Position(3:4),'axes_count',numel(axes_handles),'legend_count',numel(legends),'axes',{{}},'legends',{{}});
for k=1:numel(legends)
    audit.legends{k}=struct('box',legends(k).Box,'color',legends(k).Color,'labels',{legends(k).String});
end
for k=1:numel(axes_handles)
    ax=axes_handles(k); obj=[findall(ax,'Type','line');findall(ax,'Type','stair');findall(ax,'Type','errorbar')]; data=cell(1,numel(obj));
    for j=1:numel(obj)
        data{j}=struct('type',obj(j).Type,'x',obj(j).XData,'y',obj(j).YData,'color',obj(j).Color,'style',obj(j).LineStyle);
        if strcmp(obj(j).Type,'errorbar'), data{j}.lower=obj(j).YNegativeDelta; data{j}.upper=obj(j).YPositiveDelta; end
    end
    pictures=findall(ax,'Type','image'); pixels=cell(1,numel(pictures));
    for j=1:numel(pictures), pixels{j}=struct('x',pictures(j).XData,'y',pictures(j).YData,'values',pictures(j).CData); end
    audit.axes{k}=struct('tile',ax.Layout.Tile,'xlim',ax.XLim,'ylim',ax.YLim,'title',ax.Title.String,'curves',{data},'images',{pixels});
end
close(fig);
    function s=word(chinese,english)
        if zh, s=chinese; else, s=english; end
    end
    function labels(ax,xc,xe,yc,ye,p)
        if p==0, xlabel(ax,word(xc,xe)); ylabel(ax,word(yc,ye)); return; end
        panel=sprintf('(%c)',96+p);
        switch topic
            case 'snr'
                detail={word('用户定位','User localization'),word('RF 合并后等效信道','RF-combined channel')};
                panel=[panel ' ' detail{p}];
            case 'cdf'
                values=[5,20]; panel=sprintf('%s SNR = %d dB',panel,values(p));
            case 'grid'
                detail={word('门限内定位精度','Accuracy within the threshold'),word('定位可靠性','Localization reliability')};
                panel=[panel ' ' detail{p}];
            case 'rf'
                detail={word('门限内定位精度','Accuracy within the threshold'),word('定位可靠性','Localization reliability')};
                panel=[panel ' ' detail{p}];
            case 'iterations'
                panel=sprintf(word('%s 距离 %d-%d m','%s Range %d-%d m'),panel,5*p,5*p+5);
            case 'coverage'
                detail={word('APLE-LM 位置 RMSE','APLE-LM position RMSE'),'BCRB'}; panel=[panel ' ' detail{p}];
        end
        xlabel(ax,{word(xc,xe),panel}); ylabel(ax,word(yc,ye));
    end
end
function [fig,t]=canvas(nr,nc,w,h)
fig=figure('Visible','off','Color','w','Units','centimeters','Position',[2,2,w,h],'Renderer','painters');
t=tiledlayout(fig,nr,nc,'TileSpacing','compact','Padding','compact');
end
function h=lineplot(ax,x,y,color,style)
h=plot(ax,x,y,style,'Color',color,'MarkerSize',4.6,'MarkerFaceColor','w','LineWidth',1.35);
end
function h=probability(ax,x,r,color,style)
y=100*r.Exceedance_probability;
h=plot(ax,x,y,style,'Color',color,'MarkerSize',4.6, ...
    'MarkerFaceColor','w','LineWidth',1.25);
end
function r=getrows(t,family,method)
r=sortrows(t(t.Experiment==family & t.Method==method,:),'Point');
end
function require(t,families)
assert(all(ismember(string(families),unique(t.Experiment))),'APLE:MissingExperiment', ...
    'Required experiment is absent; select available Figures or run the missing experiment.');
end
function padlog(ax)
objects=findall(ax,'Type','line'); v=[];
for k=1:numel(objects), v=[v,reshape(objects(k).YData,1,[])]; end %#ok<AGROW>
v=v(isfinite(v)&v>0); assert(~isempty(v),'APLE:NoValidValues','No positive values to plot.');
ylim(ax,[.82*min(v),1.22*max(v)]);
end
