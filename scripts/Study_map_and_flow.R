suppressPackageStartupMessages({library(data.table);library(ggplot2)})
root<-Sys.getenv('DTS_ROOT');out<-file.path(root,'Supplementary_DTS_v1','figures')
dir.create(out,recursive=TRUE,showWarnings=FALSE)
save_plot<-function(p,stem,w,h){
 for(ext in c('png','pdf','tiff')){
  a<-list(filename=file.path(out,paste0(stem,'.',ext)),plot=p,width=w,height=h,bg='white',limitsize=FALSE)
  if(ext=='pdf'&&capabilities('cairo'))a$device<-grDevices::cairo_pdf
  if(ext!='pdf')a$dpi<-if(ext=='tiff')600 else 300
  if(ext=='tiff')a$compression<-'lzw'
  do.call(ggsave,a)
 }
}
f<-fread(file.path(root,'qc','DTS_selection_flow.csv'))
stopifnot(all(c('step','rows','pregnancies')%in%names(f)))
f[,retained:=fifelse(is.na(pregnancies),rows,pregnancies)]
f[,excluded:=shift(retained)-retained]
if(any(f$excluded<0,na.rm=TRUE))stop('Sample flow counts increase; review selection log')
f[,y:=rev(seq_len(.N))]
f[,label:=paste0(vapply(step,function(s)paste(strwrap(s,width=49),collapse='\n'),character(1)),
 '\nN = ',format(retained,big.mark=',',trim=TRUE))]
p<-ggplot()+geom_rect(data=f,aes(xmin=0,xmax=6,ymin=y-.32,ymax=y+.32),fill='#EFF5FA',colour='#234E70')+
 geom_text(data=f,aes(x=3,y=y,label=label),size=3.2,lineheight=.9)+
 geom_segment(data=f[-.N],aes(x=3,xend=3,y=y-.32,yend=y-.68),arrow=grid::arrow(length=grid::unit(.10,'inches')))+
 geom_text(data=f[!is.na(excluded)&excluded>0],aes(x=6.2,y=y+.48,label=paste0('Excluded: ',format(excluded,big.mark=',',trim=TRUE))),hjust=0,size=3)+
 coord_cartesian(xlim=c(0,8.5),ylim=c(.5,nrow(f)+.5),clip='off')+theme_void()+
 labs(caption='Counts follow the reconstructed pregnancy-cohort selection.')+
 theme(plot.caption=element_text(size=9,hjust=0),plot.margin=margin(12,12,12,12))
save_plot(p,'DTS21_sample_selection_flow',10,max(8,nrow(f)*.9))
# Use bundled maps package data only; no boundary download during finalization.
if(requireNamespace('maps',quietly=TRUE)){
 world<-ggplot2::map_data('world')
 selected<-c('Angola','Bangladesh','Burkina Faso','Democratic Republic of the Congo',
 'Ivory Coast','Ethiopia','Ghana','Jordan','Kenya','Cambodia','Lesotho','Mali','Malawi',
 'Mozambique','Nigeria','Nepal','Philippines','Senegal','Tajikistan','Tanzania','Zambia')
 missing<-setdiff(selected,unique(world$region))
 if(length(missing)){
  stop('Unmatched map names: ',paste(missing,collapse=', '))
 }else{
  world$included<-ifelse(world$region%in%selected,'Study country','Other country')
  p<-ggplot(world,aes(long,lat,group=group,fill=included))+geom_polygon(colour='white',linewidth=.12)+
   scale_fill_manual(values=c('Study country'='#196B8C','Other country'='#E7E7E7'),name=NULL)+
   coord_quickmap(xlim=c(-20,130),ylim=c(-36,42),expand=FALSE)+theme_void()+
   theme(legend.position='bottom')
  save_plot(p,'DTS21_study_countries_map',12,7)
 }
}else stop('Install maps to generate the study map')
