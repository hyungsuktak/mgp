# SDSS Stripe 82: triangular PSD and squared-coherence comparison
bands <- c('u','g','r','i','z')
tau_mat <- c(424.6893,353.1032,689.4087,827.1391,825.1932)
sd_mat <- c(0.16563299,0.06047612,0.06484202,0.04908173,0.04936020)
rho_mat <- matrix(c(1,.9915415,.9432275,.8731145,.8703337,.9915415,1,.8994494,.8244922,.8220550,.9432275,.8994494,1,.9739943,.9713308,.8731145,.8244922,.9739943,1,.9998767,.8703337,.8220550,.9713308,.9998767,1),5,byrow=TRUE)
# Updated MCARMA(1,0) BFGS fit, logLik = 487.855845654774
# Stationary marginal SDs, damping timescales, and zero-lag correlations.
tau_sde <- c(441.5315, 351.17098, 698.56327, 1154.2732, 1098.5143)
sd_sde <- c(0.1593024, 0.055834073, 0.067183569, 0.054503761, 0.054191356)
sigma_sde <- sd_sde * sqrt(2 / tau_sde)

# Pairs ordered by combn(1:5, 2): ug, ur, ui, uz, gr, gi, gz, ri, rz, iz.
pairs <- combn(1:5, 2)
rho_values <- c(0.99348184, 0.97425527, 0.89471643, 0.90444063,
                0.94365406, 0.84582196, 0.8568769, 0.96928186,
                0.97491983, 0.99969365)
rho_sde <- diag(5)
for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]; j <- pairs[2, k]
  rho_sde[i, j] <- rho_sde[j, i] <- rho_values[k]
}

# For diagonal drift a_j = 1/tau_j, stationary rho_ij equals
# driver correlation r_ij * 2*sqrt(tau_i*tau_j)/(tau_i+tau_j).
# The printed rho estimates are rounded to 8 digits and imply
# driver correlations within ~1e-8 of 1; clamp rounding overshoot.
r_sde <- diag(5)
for (k in seq_len(ncol(pairs))) {
  i <- pairs[1, k]; j <- pairs[2, k]
  raw_r <- rho_sde[i, j] * (tau_sde[i] + tau_sde[j]) /
           (2 * sqrt(tau_sde[i] * tau_sde[j]))
  if (abs(raw_r) > 1 + 1e-6) stop('Inconsistent MCARMA stationary correlation')
  r_sde[i, j] <- r_sde[j, i] <- max(-1, min(1, raw_r))
}
omega <- 10^seq(-5,-1,length.out=600)
psd <- function(w,s,t) 2*s^2*t/(1+(w*t)^2)
mat_cross <- function(w,i,j) {tc<-sqrt(2/(tau_mat[i]^(-2)+tau_mat[j]^(-2)));2*rho_mat[i,j]*sd_mat[i]*sd_mat[j]*tc/(1+(w*tc)^2)}
sde_cross <- function(w,i,j) sigma_sde[i]*sigma_sde[j]*r_sde[i,j]/((1/tau_sde[i]+1i*w)*(1/tau_sde[j]-1i*w))
mat_psd<-lapply(1:5,function(j)psd(omega,sd_mat[j],tau_mat[j]))
sde_psd<-lapply(1:5,function(j)psd(omega,sd_sde[j],tau_sde[j]))
all_psd<-unlist(c(mat_psd,sde_psd)); ylim_psd<-range(all_psd)
cols<-c('#27689B','#BD4D40')

plot_triangle <- function(file,device=c('pdf','png')) {
 op<-par(mfrow=c(5,5),mar=c(2.4,2.6,1.5,.6),oma=c(3,3,1,1),mgp=c(1.5,.45,0),tcl=-.22,cex.axis=.75)
 for(i in 1:5) for(j in 1:5) {
  if(j>i){plot.new();next}
  if(i==j) {
   plot(omega,mat_psd[[i]],type='l',log='xy',ylim=ylim_psd,xlim=range(omega),col=cols[1],lwd=2,axes=FALSE,xlab='',ylab='')
   lines(omega,sde_psd[[i]],col=cols[2],lwd=2,lty=2)
   axis(1);axis(2,las=1);box();mtext(bands[i],side=3,line=.1,font=2,cex=.95)
   if(i==1)legend('topright',legend=c('Matern','SDE mDRW'),col=cols,lty=c(1,2),lwd=2,bty='n',cex=.63)
  } else {
   m<-Mod(mat_cross(omega,i,j))^2/(mat_psd[[i]]*mat_psd[[j]])
   s<-Mod(sde_cross(omega,i,j))^2/(sde_psd[[i]]*sde_psd[[j]])
   plot(omega,m,type='l',log='x',ylim=c(0,1.03),xlim=range(omega),col=cols[1],lwd=2,axes=FALSE,xlab='',ylab='')
   lines(omega,s,col=cols[2],lwd=2,lty=2)
   axis(1);axis(2,at=c(0,.5,1),las=1);box();mtext(paste0(bands[i],'-',bands[j]),side=3,line=.1,cex=.85)
  }
 }
 mtext(expression(omega~'(rad/day)'),side=1,outer=TRUE,line=1.1,cex=1.1)
 mtext('Marginal PSD (diagonal) / squared coherence (lower triangle)',side=2,outer=TRUE,line=1.1,cex=.9)
}
plot_triangle('/mnt/data/sdss_psd_coherence_triangle.pdf','pdf')
plot_triangle('/mnt/data/sdss_psd_coherence_triangle.png','png')
