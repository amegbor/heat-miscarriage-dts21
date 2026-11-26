suppressPackageStartupMessages({library(data.table);library(ggplot2);library(patchwork)})
source(file.path(Sys.getenv('DTS_BUNDLE'),'supp','00_settings.R'), local=TRUE)
write_tab <- function(x,name,folder='tables')fwrite(x,file.path(OUT,folder,paste0(name,'.csv')))
COUNTRY_NAMES <- c(AO='Angola',BD='Bangladesh',BF='Burkina Faso',KH='Cambodia',CD='Democratic Republic of the Congo',CI="Côte d’Ivoire",ET='Ethiopia',GH='Ghana',JO='Jordan',KE='Kenya',LS='Lesotho',MW='Malawi',ML='Mali',MZ='Mozambique',NP='Nepal',NG='Nigeria',PH='Philippines',SN='Senegal',TJ='Tajikistan',TZ='Tanzania',ZM='Zambia')
country_labels <- function(x) {
 raw <- as.character(x); key <- toupper(sub('^([A-Za-z]{2}).*$', '\\1', raw))
 ans <- unname(COUNTRY_NAMES[key]); ans[is.na(ans)] <- raw[is.na(ans)]
 if(any(is.na(COUNTRY_NAMES[key]))) message('Country labels retained verbatim for: ',paste(unique(raw[is.na(COUNTRY_NAMES[key])]),collapse=', '))
 ans
}
# Run after completed stages; reports missing analyses rather than inventing rows.
suppressPackageStartupMessages({library(ggplot2);library(patchwork)})
tabs <- list.files(file.path(OUT,'tables'),full.names=TRUE,pattern='\\.csv$')
colours <- c('Model 1'='#E41A1C','Model 2'='#377EB8','Model 3'='#4DAF4A')
metric_levels <- c('UTCI','MRT','WBGT','Heat Index')
metric_colours <- c('UTCI'='#D55E00','MRT'='#009E73','WBGT'='#0072B2','Heat Index'='#CC79A7')
theme0 <- function() theme_classic(base_size=12)+theme(legend.position='bottom',plot.title=element_blank(),plot.subtitle=element_blank(),plot.caption=element_blank(),strip.background=element_rect(fill='grey95'),strip.text=element_text(face='bold'))
save_fig <- function(p,name,w=12,h=8) {
 for(ext in c('png','pdf','tiff')) {
  args <- list(filename=file.path(OUT,'figures',paste0(name,'.',ext)),plot=p,width=w,height=h,bg='white',limitsize=FALSE)
  if(ext=='pdf' && capabilities('cairo'))args$device<-grDevices::cairo_pdf
  if(ext!='pdf') args$dpi <- if(ext=='tiff') 600 else 300
  if(ext=='tiff') args$compression <- 'lzw'
  do.call(ggsave,args)
 }
}
forest <- function(z,label,name) {
 z <- copy(z[is.finite(HR)&is.finite(lower_95)&is.finite(upper_95)])
 z[, Heat_metric:=factor(Heat_metric,levels=metric_levels)]
 z[, Model:=factor(Model,levels=names(colours))]
 z[, y:=5-as.integer(Heat_metric)+c(.22,0,-.22)[as.integer(Model)]]
 z[, text_label:=sprintf('%.3f (%.3f–%.3f)',HR,lower_95,upper_95)]
 bounds <- range(c(1,z$lower_95,z$upper_95)); span <- diff(bounds)
 if(span==0) span <- .1
 text_x <- bounds[2]+.12*span
 p <- ggplot(z,aes(HR,y,colour=Model,shape=Model))+
  geom_hline(yintercept=c(1.5,2.5,3.5),colour='grey90')+
  geom_vline(xintercept=1,linetype='dashed',colour='grey55')+
  geom_segment(aes(x=lower_95,xend=upper_95,yend=y),linewidth=.7)+
  geom_point(size=3)+geom_text(aes(x=text_x,label=text_label),colour='black',hjust=0,size=3.4)+
  annotate('text',x=text_x,y=4.6,label='HR (95% CI)',hjust=0,fontface='bold',size=3.5)+
  scale_colour_manual(values=colours)+scale_shape_manual(values=c(16,17,15))+
  scale_y_continuous(breaks=1:4,labels=rev(metric_levels),limits=c(.5,4.7))+
  scale_x_continuous(breaks=pretty(bounds,n=5))+
  coord_cartesian(xlim=c(bounds[1]-.05*span,bounds[2]+.65*span),clip='off')+
  labs(x=label,y=NULL)+theme0()
 save_fig(p,name,12,6)
}
for(prefix in c('Raw_linear')) for(scale in c('RAW','IQR')) {
 files <- tabs[grepl(paste0('^',prefix,'_.*_',scale,'_GM1_GM2_GM3\\.csv$'),basename(tabs))]
 if(length(files)) {
  z <- rbindlist(lapply(files,fread),fill=TRUE);write_tab(z,paste0(prefix,'_',scale,'_combined'))
  if(uniqueN(z$Heat_metric)==4 && nrow(z)==12) forest(z,if(scale=='RAW') 'Hazard ratio per degree Celsius (95% CI)' else 'Hazard ratio per country-specific IQR (95% CI)',paste0(prefix,'_',scale,'_forest'))
  else message('Incomplete forest: ',prefix,' ',scale)
 }
}
for(scale in c('IQR','RAW')) {
 files <- tabs[grepl(paste0('^Nonlinear_(WBGT|UTCI|MRT|HI)_',scale,'_predictions\\.csv$'),basename(tabs))]
 histfiles <- tabs[grepl(paste0('^Nonlinear_(WBGT|UTCI|MRT|HI)_',scale,'_distribution\\.csv$'),basename(tabs))]
 if(length(files)!=4 || length(histfiles)!=4) {message('Need all four completed nonlinear metrics for ',scale);next}
 nl <- rbindlist(lapply(files,fread));hist <- rbindlist(lapply(histfiles,fread))
 if(nrow(unique(nl[,.(Model,Heat_metric)]))!=12) {message('Incomplete GM1-GM3 nonlinear results');next}
 write_tab(nl,paste0('Nonlinear_',scale,'_predictions_combined'));write_tab(hist,paste0('Nonlinear_',scale,'_distribution_combined'))
 nl[,Model:=factor(Model,levels=names(colours))];nl[,Heat_metric:=factor(Heat_metric,levels=metric_levels)]
 hist[,Heat_metric:=factor(Heat_metric,levels=metric_levels)]
 medians <- unique(nl[,.(Heat_metric,reference_exposure)])
 if(anyDuplicated(medians$Heat_metric)) stop('Inconsistent nonlinear reference medians')
 write_tab(medians,paste0('Nonlinear_',scale,'median_lines'))
 xlabel <- if(scale=='RAW') 'Heat exposure (degrees Celsius)' else 'Heat exposure (country-specific IQR units)'
 # Shared x limits from both curve grid and pregnancy exposures.
 anchors <- rbind(nl[,.(Heat_metric,exposure)],hist[,.(Heat_metric,exposure)])[,.(exposure=range(exposure)),by=Heat_metric]
 p <- ggplot(nl,aes(exposure,estimate,colour=Model,fill=Model))+geom_ribbon(aes(ymin=lower_95,ymax=upper_95),alpha=.15,colour=NA)+
  geom_line(linewidth=.8)+geom_vline(data=medians,aes(xintercept=reference_exposure),linetype='dotted',colour='black',linewidth=.6)+geom_hline(yintercept=1,linetype='dashed',colour='grey55')+
  geom_blank(data=anchors,aes(exposure,NA_real_),inherit.aes=FALSE)+facet_grid(Model~Heat_metric,scales='free_x')+
  scale_colour_manual(values=colours)+scale_fill_manual(values=colours)+scale_x_continuous(expand=expansion(mult=.02))+
  labs(x=xlabel,y='Hazard ratio (95% CI)')+theme0()+theme(legend.position='none')
 ph <- ggplot(hist,aes(exposure,fill=Heat_metric))+geom_histogram(bins=40,colour='white',linewidth=.15)+geom_vline(data=medians,aes(xintercept=reference_exposure),linetype='dotted',colour='black',linewidth=.6)+
  geom_blank(data=anchors,aes(exposure,NA_real_),inherit.aes=FALSE)+facet_grid(.~Heat_metric,scales='free_x')+
  scale_fill_manual(values=metric_colours,guide='none')+scale_x_continuous(expand=expansion(mult=.02))+
  labs(x=xlabel,y='Pregnancies')+theme0()
 # Only the version with exposure histograms is exported.
 save_fig(wrap_plots(p+labs(x=NULL),ph,ncol=1,heights=c(3,1)),paste0('Nonlinear_',scale,'_GM1_GM2_GM3_with_histograms'),12,9.8)
}
# Country forest includes failures in its output table; points plot stable estimates.
files <- tabs[grepl('^Country_(WBGT|UTCI|MRT|HI)_GM3_results\\.csv$',basename(tabs))]
if(length(files)) {
 z <- rbindlist(lapply(files,fread),fill=TRUE);write_tab(z,'Country_GM3_combined')
 if(!'HR' %in% names(z)) message('No stable country estimates') else {
  z[,Heat_metric:=factor(Heat_metric,levels=metric_levels)];z[,country:=country_labels(country)];allcc <- sort(unique(z$country))
  z[,country:=factor(country,levels=rev(allcc))]
  p <- ggplot(z,aes(HR,country,colour=Heat_metric))+geom_vline(xintercept=1,linetype='dashed',colour='grey55')+
   geom_errorbar(aes(xmin=lower_95,xmax=upper_95),orientation='y',width=.2,na.rm=TRUE)+geom_point(na.rm=TRUE)+
   geom_text(data=z[status!='passed'],aes(x=1,y=country,label='No stable estimate'),inherit.aes=FALSE,size=2.5)+
   facet_wrap(~Heat_metric,ncol=2,scales='free_x')+scale_y_discrete(drop=FALSE)+scale_colour_manual(values=metric_colours,guide='none')+
   labs(x='Hazard ratio per country-specific IQR (95% CI)',y='Country')+theme0()
  save_fig(p,'Country_GM3_forest',14,max(10,length(allcc)*.6))
 }
}

primary_tables<-file.path(ROOT,'tables','Linear_Nonlinear')
f<-file.path(primary_tables,paste0('DTS_Linear_',c('WBGT','UTCI','MRT','Heat_Index'),'_GM1_GM2_GM3_CountryFE_MotherRE_v2.csv'))
if(any(!file.exists(f)))stop('Primary linear tables incomplete')
forest(rbindlist(lapply(f,fread)), 'Hazard ratio per country-specific IQR (95% CI)', 'DTS_linear_IQR_GM1_GM2_GM3_forest')

source(file.path(Sys.getenv('DTS_BUNDLE'),'Recovery_figures.R'), local=TRUE)
source(file.path(Sys.getenv('DTS_BUNDLE'),'PIF_and_weights_figures.R'), local=TRUE)
