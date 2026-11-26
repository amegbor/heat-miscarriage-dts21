suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
root<-Sys.getenv('DTS_ROOT');ptab<-file.path(root,'PIF_QGComp_v3','tables')
out<-file.path(root,'Supplementary_DTS_v1','figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
read_four<-function(pattern){f<-list.files(ptab,pattern=pattern,full.names=TRUE);if(length(f)!=4)stop('Expected four files: ',pattern);rbindlist(lapply(f,fread),fill=TRUE)}
save_selected<-function(p,name,w=11,h=8){
 for(ext in c('png','pdf','tiff')){
  a<-list(filename=file.path(out,paste0(name,'.',ext)),plot=p,width=w,height=h,bg='white',limitsize=FALSE)
  if(ext=='pdf'&&capabilities('cairo'))a$device<-grDevices::cairo_pdf
  if(ext!='pdf')a$dpi<-if(ext=='tiff')600 else 300
  if(ext=='tiff')a$compression<-'lzw'
  do.call(ggsave,a)
 }
}
theme_m<-function()theme_classic(base_size=12)+theme(legend.position='bottom',legend.title=element_blank(),strip.background=element_rect(fill='grey96'))
ml<-c('UTCI','MRT','WBGT','Heat Index')
sc<-c('cap_p90','cap_p75','cap_p50','shift_down_1iqr')
sl<-c(cap_p90='Cap at country 90th percentile',cap_p75='Cap at country 75th percentile',cap_p50='Cap at country median',shift_down_1iqr='Reduce by one country IQR')
col<-c(cap_p90='#0072B2',cap_p75='#009E73',cap_p50='#D55E00',shift_down_1iqr='#CC79A7')
pif<-read_four('_PIF_results[.]csv$');pif[,metric:=factor(metric,levels=ml)]
fwrite(pif,file.path(ptab,'PIF_combined.csv'))
z<-copy(pif[population=='Pooled']);z[,scenario:=factor(scenario,levels=sc)]
z[,label:=sprintf('%.2f (%.2f-%.2f)',PIF_percent,lower_95,upper_95)]
right<-max(z$upper_95);span<-diff(range(c(0,z$lower_95,z$upper_95)));z[,label_x:=right+.06*span]
p<-ggplot(z,aes(PIF_percent,scenario,colour=scenario))+
 geom_vline(xintercept=0,linetype='dashed',colour='grey60')+
 geom_segment(aes(x=lower_95,xend=upper_95,yend=scenario),linewidth=.7)+geom_point(size=2.8)+
 geom_text(aes(x=label_x,label=label),hjust=0,colour='black',size=3)+
 facet_wrap(~metric,ncol=2)+scale_colour_manual(values=col,labels=sl)+
 scale_y_discrete(labels=sl)+scale_x_continuous(limits=c(min(0,z$lower_95),right+.85*span))+
 labs(x='Potential impact fraction (%, 95% CI)',y=NULL)+theme_m()+theme(legend.position='none')
save_selected(p,'PIF_scenario_forest',14,8)
cc<-copy(pif[population!='Pooled']);cc[,country:=country_labels(population)]
cc[,scenario:=factor(scenario,levels=sc)]
cc[,country:=factor(country,levels=rev(sort(unique(country))))]
p<-ggplot(cc,aes(scenario,country,fill=PIF_percent))+geom_tile(colour='white',linewidth=.2)+
 facet_wrap(~metric,ncol=2)+scale_fill_gradient(low='#FFF5EB',high='#A63603',name='Potential impact fraction (%)')+
 scale_x_discrete(labels=c('90th percentile','75th percentile','Median','One-IQR shift'))+
 labs(x=NULL,y=NULL)+theme_m()+theme(legend.title=element_text(size=10),axis.text.x=element_text(angle=30,hjust=1))
save_selected(p,'PIF_country_heatmap',13,12)
curves<-read_four('_PIF_net_risk_curves[.]csv$');curves[,metric:=factor(metric,levels=ml)]
fwrite(curves,file.path(ptab,'PIF_net_risk_curves_combined.csv'))
p<-ggplot(curves,aes(gestational_interval,100*net_risk,colour=scenario))+
 geom_line(linewidth=.7)+geom_point(size=1.8)+facet_wrap(~metric,ncol=2)+
 scale_colour_manual(values=c(natural='#222222',col),labels=c(natural='Observed exposure',sl))+
 labs(x='Gestational interval',y='Standardised cumulative net miscarriage risk (%)')+theme_m()
save_selected(p,'PIF_net_risk_curves',11,8)

# Main Figure 4: package joint heat-air HRs only.
hr<-read_four('^package_Air_(WBGT|UTCI|MRT|HI)_cloglog_results[.]csv$')
hr[,Heat_metric:=ml[match(mixture,c('Air_UTCI','Air_MRT','Air_WBGT','Air_HI'))]]
hr[,Heat_metric:=factor(Heat_metric,levels=ml)]
hr[,label:=sprintf('%.3f (%.3f-%.3f)',HR,lower_95,upper_95)]
x<-max(hr$upper_95)+.025
p<-ggplot(hr,aes(HR,Heat_metric,colour=Heat_metric))+
 geom_vline(xintercept=1,linetype='dashed',colour='grey55')+
 geom_segment(aes(x=lower_95,xend=upper_95,yend=Heat_metric),linewidth=.8)+geom_point(size=3)+
 geom_text(aes(x=x,label=label),hjust=0,colour='black',size=3.7)+
 scale_colour_manual(values=c(UTCI='#D55E00',MRT='#009E73',WBGT='#0072B2','Heat Index'='#CC79A7'))+
 scale_x_continuous(limits=c(min(1,hr$lower_95)-.01,x+.16))+
 labs(x='Joint hazard ratio per simultaneous one-quartile increase (95% CI)',y=NULL)+theme_m()+theme(legend.position='none')
save_selected(p,'DTS_joint_air_heat_HR_forest',11,5)

weights<-read_four('^package_Air_(WBGT|UTCI|MRT|HI)_cloglog_weights[.]csv$')
plots<-list();component_names<-c(CO='Carbon monoxide',DUST='Mineral dust',PM='PM2.5',SO='Sulphur dioxide',SU='Sulphate')
for(j in seq_along(ml)){
 key<-c('UTCI','MRT','WBGT','HI')[j];z<-copy(weights[mixture==paste0('Air_',key)])
 z[,label:=unname(component_names[sub('_T1.*$','',component)])]
 z[is.na(label),label:=ml[j]]
 z[,signed_weight:=ifelse(direction=='Negative',-weight,weight)]
 z[,label:=factor(label,levels=label[order(signed_weight)])]
 z[,Heat_metric:=ml[j]]
 p<-ggplot(z,aes(signed_weight,label,fill=direction))+geom_col(width=.7)+geom_vline(xintercept=0,colour='grey50')+
 facet_wrap(~Heat_metric)+scale_fill_manual(values=c(Positive='#D55E00',Negative='#0072B2',Zero='grey70'))+
 scale_x_continuous(limits=c(-1,1))+labs(x='Directional weight',y=NULL)+theme_m()
 plots[[j]]<-p
}
save_selected(wrap_plots(plots,ncol=2,guides='collect'),'Mixture_directional_weights_package',13,10)

