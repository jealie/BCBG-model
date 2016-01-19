# must set the working directory of the test file before launching this script
## For example: setwd("/media/jean/ext4/phd/Work/lastcode/code/BCBG/jneurocomp/")


plot_radar = function(start_ind, label)
{
  require('fmsb')
  data_msn = data.frame(matrix(NA,nrow=2+nrow(tst),ncol=8))
  data_msn[1,] = 2
  data_msn[2,] = 0
  start_ind_ref = start_ind-8*5-9
  for (j in 1:nrow(tst)) {
    data_msn[j+2,] = (tst[j,start_ind:(start_ind+7)] / tst[j,start_ind_ref:(start_ind_ref+7)])[c(4, 2, 1, 3, 5, 6, 8, 7)]
  }
  radarchart(data.frame(data_msn), plty = 1, pcol = 'red', title = label, vlabels = '',xpd=T,axistype=0, seg=2)
}

tst=read.csv('model_output',sep=" ",header=F)

n_channels = 8

str_nuclei = c('msn', 'fsi', 'stn', 'gpe', 'gpi')

i = 1
for (nuc_i in seq_along(str_nuclei)) {
  nuc = str_nuclei[nuc_i]
  for (n_i in 1:n_channels) {
    colnames(tst)[i] = paste0(nuc, n_i)
    i = i + 1
  }
}
i = i - 1

colnames(tst)[i+1] = 'AMPA antagonist in GPe' # GPe
colnames(tst)[i+2] = 'AMPA and GABAA antagonists in GPe' # GPe
colnames(tst)[i+3] = 'NMDA antagonist in GPe' # GPe
colnames(tst)[i+4] = 'GABAA antagonist in GPe' # GPe

colnames(tst)[i+5] = 'NMDA antagonist in GPi' # GPi
colnames(tst)[i+6] = 'NMDA and AMPA antagonists in GPi' # GPi
colnames(tst)[i+7] = 'AMPA antagonist in GPi' # GPi
colnames(tst)[i+8] = 'GABAA antagonist in GPi' #GPi
colnames(tst)[i+9] = 'GABAA, NMDA and AMPA antagonists in GPi' # GPi

i = i + 10

for (nuc_i in seq_along(str_nuclei)) {
  nuc = str_nuclei[nuc_i]
  for (n_i in 1:n_channels) {
    colnames(tst)[i] = paste0('selec_',nuc, n_i)
    i = i + 1
  }
}




par(mai=c(0.35,0.5,0.1,0.05))
layout(matrix(1:8,ncol=2))

lab = 'AMPA antagonist in GPe'
plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,120),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('gpe1',lab)],col='darkblue',type='b',pch=0))
axis(side=1,at=c(1,2),labels=c('Control',lab),lwd=0, lwd.ticks = 1)

lab = 'AMPA and GABAA antagonists in GPe'
plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,140),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('AMPA antagonist in GPe',lab)],col='darkblue',type='b',pch=0))
axis(side=1,at=c(1,2),labels=c('AMPA antagonist in GPe',lab),lwd=0, lwd.ticks = 1)

lab = 'NMDA antagonist in GPe'
plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,120),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('gpe1',lab)],col='darkblue',type='b',pch=0))
axis(side=1,at=c(1,2),labels=c('Control',lab),lwd=0, lwd.ticks = 1)

lab = 'GABAA antagonist in GPe'
plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,250),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('gpe1',lab)],col='darkblue',type='b',pch=0))
axis(side=1,at=c(1,2),labels=c('Control',lab),lwd=0, lwd.ticks = 1)


lab = 'NMDA antagonist in GPi'
plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,140),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('gpi1',lab)],col='darkblue',type='b',pch=0))
axis(side=1,at=c(1,2),labels=c('Control',lab),lwd=0, lwd.ticks = 1)

lab = 'NMDA and AMPA antagonists in GPi'
plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,120),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('NMDA antagonist in GPi',lab)],col='darkblue',type='b',pch=0))
axis(side=1,at=c(1,2),labels=c('NMDA antagonist in GPi',lab),lwd=0, lwd.ticks = 1)

# lab = 'AMPA antagonist in GPi'
# plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,250),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
# sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('gpi1',lab)],col='darkblue',type='b',pch=0))
# axis(side=1,at=c(1,2),labels=c('Control',lab),lwd=0, lwd.ticks = 1)

lab = 'GABAA antagonist in GPi'
plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,250),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('gpi1',lab)],col='darkblue',type='b',pch=0))
axis(side=1,at=c(1,2),labels=c('Control',lab),lwd=0, lwd.ticks = 1)

# lab = 'GABAA, NMDA and AMPA antagonists in GPi'
# plot(1,1,type='n',xlim=c(0.66,2.33),ylim=c(0,250),xlab='',xaxt='n',ylab='',las=1,xaxs='i',yaxs='i')
# sapply(1:nrow(tst),function(iii)lines(c(1,2),tst[iii,c('gpi1',lab)],col='darkblue',type='b',pch=0))
# axis(side=1,at=c(1,2),labels=c('Control',lab),lwd=0, lwd.ticks = 1)




## SELECTION TEST

first_msn = which(colnames(tst) == 'selec_msn1')
first_fsi = which(colnames(tst) == 'selec_fsi1')
first_stn = which(colnames(tst) == 'selec_stn1')
first_gpe = which(colnames(tst) == 'selec_gpe1')
first_gpi = which(colnames(tst) == 'selec_gpi1')

par(mai=c(0.35,0.5,0.1,0.05))
layout(matrix(1:6,ncol=2))

plot(1,1,type='n',xaxt='n',yaxt='n',bty='n') # should replace by CTX
plot_radar(first_msn, 'MSN')
plot_radar(first_fsi, 'FSI')
plot_radar(first_stn, 'STN')
plot_radar(first_gpe, 'GPe')
plot_radar(first_gpi, 'GPi')
