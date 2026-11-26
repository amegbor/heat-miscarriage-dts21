# Source in OOD RStudio. Reads completed tables; no models are fitted.
suppressPackageStartupMessages({library(data.table);library(ggplot2)})
root<-Sys.getenv('DTS_ROOT',file.path(getwd(),'output'))
base<-file.path(root,'Supplementary_DTS_v1')
tables<-file.path(base,'tables');figures<-file.path(base,'figures')
dir.create(figures,recursive=TRUE,showWarnings=FALSE)
metrics<-c('UTCI','MRT','WBGT','Heat Index')
keys<-c('UTCI','MRT','WBGT','HI')
models<-c('Model 1','Model 2','Model 3')
engines<-c('package','mixed')
model_colours<-c('Model 1'='#D55E00','Model 2'='#0072B2','Model 3'='#009E73')
engine_colours<-c(package='#D55E00',mixed='#0072B2')
engine_labels<-c(package='Quantile g-computation',mixed='Mixed-effects sensitivity model')
theme0<-function() theme_classic(base_size=12)+theme(legend.position='bottom',
 legend.title=element_blank(),strip.background=element_rect(fill='grey96',colour='grey75'),
 strip.text=element_text(face='bold'),plot.margin=margin(8,20,8,8))
save_fig<-function(p,stem,w=11,h=6) {
 for(ext in c('png','pdf','tiff')) {
  args<-list(filename=file.path(figures,paste0(stem,'.',ext)),plot=p,width=w,height=h,bg='white',limitsize=FALSE)
  if(ext=='pdf' && capabilities('cairo')) args$device<-grDevices::cairo_pdf
  if(ext!='pdf') args$dpi<-if(ext=='tiff') 600 else 300
  if(ext=='tiff') args$compression<-'lzw'
  do.call(ggsave,args)
 }
 message('Saved: ',file.path(figures,stem))
}
read_expected<-function(names0) {
 p<-file.path(tables,names0)
 if(any(!file.exists(p))) stop('Missing tables: ',paste(basename(p[!file.exists(p)]),collapse=', '))
 rbindlist(lapply(p,function(f) {z<-fread(f);z[,source_file:=basename(f)];z}),fill=TRUE)
}
check_estimates<-function(z,column) {
 if(any(!is.finite(z[[column]])|!is.finite(z$lower_95)|!is.finite(z$upper_95))) stop('Nonfinite estimates')
 if(any(z$lower_95<=0|z[[column]]<z$lower_95|z[[column]]>z$upper_95)) stop('Invalid estimates/CIs')
}
forest<-function(z,estimate,group,levels0,colours,labels,xlab) {
 z<-copy(z);z[,Heat_metric:=factor(Heat_metric,levels=metrics)]
 z[,(group):=factor(get(group),levels=levels0)]
 offsets<-seq(.22,-.22,length.out=length(levels0))
 z[,y:=5-as.integer(Heat_metric)+offsets[as.integer(get(group))]]
 z[,value:=get(estimate)]
 z[,text:=sprintf('%.3f (%.3f-%.3f)',value,lower_95,upper_95)]
 left<-min(c(z$lower_95,1))*.98;right<-max(c(z$upper_95,1))
 span<-log(right)-log(left)
 text_x<-exp(log(right)+.12*span);edge<-exp(log(right)+.90*span)
 z[,label_x:=text_x]
 ticks<-sort(unique(c(1,pretty(c(left,right),n=4))))
 ticks<-ticks[ticks>0 & ticks>=left & ticks<=right]
 ggplot(z,aes(x=value,y=y,colour=.data[[group]],shape=.data[[group]]))+
  geom_vline(xintercept=1,linetype='dashed',colour='grey55')+
  geom_segment(aes(x=lower_95,xend=upper_95,yend=y),linewidth=.65)+
  geom_point(size=2.7)+geom_text(aes(x=label_x,label=text),hjust=0,size=3.3,show.legend=FALSE)+
  scale_y_continuous(breaks=4:1,labels=metrics,limits=c(.5,4.5))+
  scale_x_log10(breaks=ticks)+coord_cartesian(xlim=c(left,edge),clip='off')+
  scale_colour_manual(values=colours,breaks=levels0,labels=labels)+
  scale_shape_manual(values=setNames(c(16,17,15)[seq_along(levels0)],levels0),breaks=levels0,labels=labels)+
  labs(x=xlab,y=NULL)+theme0()
}

# 1. Confirm pooled-year provenance, then replace the two five-year forests.
five<-read_expected(as.vector(outer(keys,c('IQR','RAW'),function(k,s)paste0('Five_year_reported_',k,'_',s,'_GM1_GM2_GM3.csv'))))
if(nrow(five)!=24 || anyDuplicated(five[,.(Heat_metric,Exposure_scale,Model)])) stop('Expected 24 distinct five-year estimates')
if(!'calendar_year_adjustment' %in% names(five)) stop('Pooling provenance column is missing')
if(any(!grepl('2015-2016 pooled',five[Model %in% c('Model 2','Model 3'),calendar_year_adjustment],fixed=TRUE))) stop('Five-year tables do not all use revised year pooling')
check_estimates(five,'HR')
fwrite(five,file.path(tables,'Five_year_reported_verified_recovery_combined.csv'))
for(s in c('IQR','RAW')) {
 z<-five[grepl(paste0('_',s,'_GM1_GM2_GM3.csv$'),source_file)]
 p<-forest(z,'HR','Model',models,model_colours,models,
  if(s=='IQR') 'Hazard ratio per country-specific IQR (95% CI)' else 'Hazard ratio per 1 degree Celsius (95% CI)')
 save_fig(p,paste0('Five_year_reported_',s,'_forest'))
}

# 2. Conditional joint HR: simultaneous one-quartile-category increase.
hr<-read_expected(as.vector(outer(keys,engines,function(k,e)paste0('Air_',k,'_',e,'_net_RR_conditional_HR.csv'))))
hr[,Heat_metric:=metrics[match(sub('^Air_','',mixture),keys)]]
if(nrow(hr)!=8||anyNA(hr$Heat_metric)||anyDuplicated(hr[,.(Heat_metric,engine)])) stop('Expected eight joint HRs')
check_estimates(hr,'HR')
fwrite(hr,file.path(tables,'Mixture_conditional_HR_comparison_combined.csv'))
save_fig(forest(hr,'HR','engine',engines,engine_colours,engine_labels,
 'Joint hazard ratio per simultaneous one-quartile increase (95% CI)'),
 'Mixture_conditional_HR_package_vs_mixed')

# 3. Standardized cumulative net-risk RRs over intervals 1-6.
rr<-read_expected(as.vector(outer(keys,engines,function(k,e)paste0('Air_',k,'_',e,'_net_RR_results.csv'))))
if(nrow(rr)!=40||anyDuplicated(rr[,.(Heat_metric,engine,scenario)])) stop('Expected 40 mixture risk rows')
if(any(rr$horizon!='1,2,3,4,5,6')) stop('Unexpected standardization horizon')
if(any(abs(rr$RR-rr$standardized_net_risk/rr$reference_net_risk)>1e-9)) stop('Risk-ratio arithmetic mismatch')
check_estimates(rr,'RR')
fwrite(rr,file.path(tables,'Mixture_standardized_net_RR_combined.csv'))
z<-copy(rr[scenario %in% paste0('joint_',0:3)])
z[,quartile:=as.integer(sub('joint_','',scenario))+1L]
z[,Heat_metric:=factor(Heat_metric,levels=metrics)]
z[,engine:=factor(engine,levels=engines)]
p<-ggplot(z,aes(quartile,RR,colour=engine,fill=engine,group=engine))+
 geom_hline(yintercept=1,linetype='dashed',colour='grey55')+
 geom_ribbon(aes(ymin=lower_95,ymax=upper_95),alpha=.12,colour=NA)+
 geom_line(linewidth=.75)+geom_point(aes(shape=engine),size=2.3)+
 facet_wrap(~Heat_metric,ncol=2)+scale_x_continuous(breaks=1:4,labels=paste0('Q',1:4))+
 scale_colour_manual(values=engine_colours,labels=engine_labels)+
 scale_fill_manual(values=engine_colours,labels=engine_labels)+
 scale_shape_manual(values=c(package=16,mixed=17),labels=engine_labels)+
 labs(x='Quartile assigned to every mixture component',y='Cumulative net-risk ratio versus all Q1 (95% CI)')+theme0()
save_fig(p,'Mixture_standardized_net_RR_joint_quartile_curves',11,8)
writeLines(capture.output(sessionInfo()),file.path(base,'qc','Recovery_figure_sessionInfo.txt'))

