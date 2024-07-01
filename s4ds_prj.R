library(stats)
library(poweRlaw)
set.seed(123)

linear_reg = function(x, y, weights=NULL,  main, plot = TRUE){
  
  # prm for log() need to be non-negative
  x = x[y>0]
  y = y[y>0]

  fit = lm(log(y)~log(x), weights = weights)
  alpha = fit[2]
  
  if (plot){
    # plot regression line
    plot(log(x), log(y))
    mtext(bquote(y == .(round(fit$coefficients[1],4)) + .(round(fit$coefficients[2],4)) * x),  adj=1, padj=0, col='#9f3e35', cex = 0.8) 
    abline(fit, col='#b6410f', lw=2)
    # plot residuals
    fitted1 = fitted(fit)
    res = resid(fit)
    plot(fitted1, res, xlab="fitted", ylab="residuals", main = main) 
  }
  
  fit
  
}

linear_bins = function(data, b_min, width){
  
  # R will fail in using seq() when the end numbers are extremely large
  # so we use for loop to create breaks 
  # as if the bins with 0 observation are filtered
  
    breaks = c(b_min )
    prev = b_min
    
    for (r in sort(data)){
      b = floor(r/width)*width
      if (b != prev) {
        breaks = c(breaks, b)
        prev = b
      }
    }
 
  #  breaks[length(breaks)] = Inf
  bins = cut(sort(data), breaks = breaks, right = FALSE)
  return(list(bins=bins, b=breaks[1:(length(breaks)-1)]))
}

log_bins = function(data, c, xmin){

  if (floor(xmin)==0){
    max = ceiling(log(max(data)) / log(c))
    breaks = exp(0:max*log(c))
    breaks = c(xmin, exp(0:max*log(2)))
  } else{
    max = ceiling(log(max(data) / xmin) / log(c))
    breaks = floor(xmin)*exp(0:max*log(c))
  }
  #  breaks[length(breaks)+1] = Inf
  bins = cut(sort(data), breaks = breaks, right = FALSE)

  return(list(bins=bins, b=breaks[1:(length(breaks)-1)]))
}

# create binned data
binned_data = function(data, bin_method, c = 2, xmin = min(data), b_min = min(data), width=10, N=length(data)){
  
  if (bin_method == 'log') bins = log_bins(data, c=c, xmin=xmin)
  if (bin_method == 'linear') bins = linear_bins(data, b_min =b_min, width=width)
  
  # h: # observations in the bins
  # b: lower boundary of the bins
  binned = table(bins$bins) 
  h = as.numeric(binned)
  b = bins$b
  binned = table(bins$bins)/N
  
  return(list(h = h, b = b, binned = binned))
}

N = 10^4
xmin = 10
alpha = 3.5

# continuos rpl
rpl = rplcon(N, xmin, alpha)
linear_binned = binned_data(rpl, 'linear', b_min=10)
log2_binned = binned_data(rpl, 'log', xmin=xmin)

# linear regression not suitable for powerlaw
par(mfrow = c(2, 2))
par(bg = "#f3f1ec")
fit = linear_reg(linear_binned$b, as.numeric(linear_binned$binned), main = 'OLS: linear binned')
fit = linear_reg(log2_binned$b, as.numeric(log2_binned$binned), main = 'OLS: log binned')

# MLE for estimating alpha
alpha_est = function(data, binned, n, method, c=2, width = 10, bmin_pos=1, interval=seq(1.5, 3.5, 0.01)){
  
  b = binned$b[binned$h>0]
  h = binned$h[binned$h>0]
  bmin = b[bmin_pos]
  if (method == 'log') b_1 = c(b[-1], b[length(b)]*c)
  if (method == 'linear') b_1 = c(b[-1], b[length(b)]+width)
  
  LL = function(alpha, data) {
    sum = sum(h*log( b^(1-alpha) - b_1^(1-alpha)))
    return(  n*(alpha-1)*log(bmin) + sum )
  }
#  alpha_est = optimize(LL, interval, data=data,  maximum=T)
#  alpha_est$maximum
  alpha_hat = interval[which.max(sapply(interval, FUN=function(alpha) LL(alpha, data)))]
  alpha_hat
}

# experiment on the accuracy of MLE techniques wrt OLS/WLS
par(mfrow = c(1, 1))
alpha_range = seq(1.5, 3.5, 0.1)
xmin=10
alpha_log = numeric()
alpha_log_reg = numeric()
alpha_log_wreg = numeric()
for (alpha in alpha_range){
  rpl = rplcon(N, xmin, alpha)
  log2_binned = binned_data(rpl, 'log', c=2, xmin=xmin)
  
  alpha_hat_log = alpha_est(rpl, log2_binned, N, 'log', c=2)
  alpha_log = c(alpha_log, alpha_hat_log)
 
  coef = -linear_reg(log2_binned$b, as.numeric(log2_binned$binned), plot = FALSE)$coefficients[2]
  alpha_log_reg = c(alpha_log_reg, coef)
  
  wcoef = -linear_reg(log2_binned$b, as.numeric(log2_binned$binned), plot = FALSE, weights = log2_binned$h[log2_binned$h>0])$coefficients[2]
  alpha_log_wreg = c(alpha_log_wreg, wcoef)
  
  
}

alpha_linear = numeric()
alpha_linear_reg = numeric()
alpha_linear_wreg = numeric()
for (alpha in alpha_range){
  rpl = rplcon(N, xmin, alpha)
  linear_binned = binned_data(rpl, 'linear', b_min = xmin)
  
  alpha_hat_linear = alpha_est(rpl, linear_binned, N, 'linear' )
  alpha_linear = c(alpha_linear, alpha_hat_linear)
  
  coef = -linear_reg(linear_binned$b, as.numeric(linear_binned$binned), plot = FALSE)$coefficients[2]
  alpha_linear_reg = c(alpha_linear_reg, coef)
  
  wcoef = -linear_reg(linear_binned$b, as.numeric(linear_binned$binned), plot = FALSE, weights = linear_binned$h)$coefficients[2]
  alpha_linear_wreg = c(alpha_linear_wreg, wcoef)
  
}

par(mar=c(4,4,1,1))
par(mfrow = c(1, 1))

plot(alpha_range, alpha_log, type = 'b', lwd = 2, col='#0B425E',cex = 2, xlab='true α', ylab = 'estimated α')
points(alpha_range, alpha_linear, pch = 4, type = 'b', lwd = 2, col='#934833',cex = 1.3)
points(alpha_range, alpha_log_reg, pch = 2, type = 'b', lwd = 2, col='#6f8278',cex = 1.6)
points(alpha_range, alpha_linear_reg, pch = 0, type = 'b', lwd = 2, col='#d5ba82',cex = 1.6)
points(alpha_range, alpha_log_wreg, pch = 5, type = 'b', lwd = 2, col='#b6410f',cex = 1.6)
points(alpha_range, alpha_linear_wreg, pch = 6, type = 'b', lwd = 2, col='#522b5b',cex = 1.6)
abline(a=0, b=1, lwd=2, lty = 'dashed')
legend("topleft", 
       legend = c("MLE: Log", "MLE: Linear", "OLS: Log", "OLS: Linear", "WLS: Log", "WLS: Linear"), 
       col = c('#0B425E', '#934833', '#6f8278', '#d5ba82', '#b6410f', '#522b5b'), 
       pch = c(1, 4, 2, 0, 5, 6), 
       lwd = 1.5,
       cex = 0.8) 

# ks distance to choose best bmin
empirical_cdf = function(binned){

  temp = cumsum(rev(binned$h))
  temp_flipped = rev(temp)
  n = sum(binned$h)
  cx = 1 - temp_flipped / n
  
  cx
}
ks_dist = function(binned, bmin, alpha_hat){
  
  if (is.null(bmin) || is.null(alpha_hat) || length(binned$b) == 0) {
    return(Inf)
  }
  
  empcdf = empirical_cdf(binned)
  theorcdf = 1 - (binned$b/bmin)^(1 - alpha_hat)
  max(abs(empcdf - theorcdf))
}

fit_pl = function(orig_data, binned, method='log', c=2, width=10, interval=seq(1.5, 3.5, 0.01) ){
  
  min_ks=Inf
  b = binned$b
  kss = numeric()
  
  # require fit to span at at least 2 bins
  for (bmin in (b[1:(length(b)-1)])) {
    
    data = orig_data[orig_data>=bmin]
    binned = binned_data(data, method, c=c, xmin=bmin)
    alpha_hat = alpha_est(data, binned, length(data), method, c=c, interval = interval)
    ks = ks_dist(binned, bmin, alpha_hat)
    kss = c(kss, ks)
    
    if (ks < min_ks) {
      min_ks = ks
      best_bmin = bmin
      best_alpha = alpha_hat
    }
    
  }
  
  return(list(best_bmin = best_bmin, best_alpha = best_alpha, kss=kss, min_ks=min_ks))
  
}
# semi-parametric bootstrap
generate_synthetic_data <- function(data, binned, bmin, alpha) {
  
  b = binned$b
  N <- length(data)
  if (which(b == bmin) == 1){
    x = rplcon(N, bmin, alpha)
  }else{
  # proportion of power law
    nz <- sum(binned$h[b >= bmin])
    pz <- nz / N
    
    # data below bmin
    ind <- which(b>=bmin)[1]
    by <- b[1:ind] # bins < b_min
    ly <- by[1:(length(by)-1)]
    uy <- by[2:length(by)] 
    y <- binned$h[1:(ind-1)]
    
    n1 <- sum(runif(N) > pz)
    temp <- (ly + uy) / 2
    temp2 <- unlist(mapply(rep, temp, y)) 
    temp2 <- sample(temp2)
    x1 <- temp2[ceiling(length(temp2) * runif(n1))]
    
    n2 <- N - n1
    #x2 <- bmin * (1 - runif(n2))^(-1 / (alpha - 1))
    x2 <- rplcon(n2, bmin, alpha)
    x <- c(x1, x2)
  }
  return(x)
}

# calculate p-value
p_value = function(data, binned, bmin, alpha, ks, method='log', significant_level=0.03,  c=2){
  
  range = 1:ceiling(0.25*significant_level^-2)
  min_ks = numeric()
  for (i in range){
    synthetic_data <- generate_synthetic_data(data, binned, bmin, alpha)
    binned_synthetic <- binned_data(synthetic_data, method, c, xmin=min(synthetic_data))
    ksd = fit_pl(synthetic_data, binned_synthetic, method=method)$min_ks
    min_ks = c(min_ks, ksd)
    
  }
  
  p = sum(min_ks>=ks)/length(min_ks)
  return(list(min_ks = min_ks, p = p))
  
}

# real data application
us_pop = read.csv("~/us_pop.txt", sep="")
attach(us_pop)
binned_pop = binned_data(population, 'log', 2, xmin=min(population))
fit = fit_pl(population, binned_pop)

print(paste("Best bmin:", fit$best_bmin))
print(paste("Best alpha:", fit$best_alpha))

# distribution of ks distance for each possible bmin
plot(fit$kss, type = 'b',  pch = 16, lwd = 2, col='#103d3e', cex = 2, xlab='candidate bmins', ylab = 'KS distance')

# best bmin in paper
bmin = binned_pop$b[17]
data = population[population>=bmin]
binned = binned_data(data, 'log', c=2, xmin=bmin)
alpha_hat = alpha_est(data, binned, length(data), 'log', c=2)
alpha_hat

# pvalue
bmin = fit$best_bmin
alpha = fit$best_alpha
ks = fit$min_ks
p = p_value(population, binned_pop, bmin, alpha, ks)
p$p
# visualize the bootstraps ks distance
par(mar=c(4,4,1,1))
par(mfrow = c(1, 2))
min_ks = p$min_ks
ks = fit$min_ks
plot(min_ks, col = '#0B425E', lwd = 1.5, cex = 1.5, xlab="")
abline(h=ks, col = '#b6410f', lwd = 3)

# using best_bmin in paper
p_b17 = p_value(population, binned_pop, ks=fit$kss[17], bmin=binned_pop$b[17], alpha=2.38)
p_b17$p
min_ks_b17 = p_b17$min_ks
ks_b17 = fit$kss[17]
plot(min_ks_b17, col = '#0B425E', lwd = 1.5, cex = 1.5,xlab="")
abline(h=ks_b17, col = '#b6410f', lwd = 3)

# average p-value
# computational expensive
# just try N=1000
bmin=16
alpha=2.5
sample_data <- rplcon(1000, bmin, alpha)

p_vals = numeric(20) # Number of repetitions for averaging
binned_sample = binned_data(sample_data, 'log', c = 2, xmin = bmin)
fit_sample = fit_pl(sample_data, binned_sample, method = 'log', c = 2)
ks_sample = fit_sample$min_ks
for (j in 1:20) {
  p_val = p_value(sample_data, binned_sample, fit_sample$best_bmin, fit_sample$best_alpha, ks_sample)$p
  p_vals[j] = p_val
}


print(paste('average p-value: ', mean(p_vals)))

          