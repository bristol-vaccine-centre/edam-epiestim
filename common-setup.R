# BRMSM0070 Epidemic Data Analysis and Modelling: 
# Uses and limitations of EpiEstim in practice

## Check dependencies
if (!require("devtools")) install.packages("devtools")
if (!require("tidyverse")) install.packages("tidyverse")
if (!require("roogledocs")) devtools::install_github("terminological/roogledocs")
if (!require("EpiEstim")) install.packages("EpiEstim")

## Serial interval data from meta-analysis paper.

here::i_am("common-setup.R")

# Set general theme
ggplot2::theme_set(ggplot2::theme_bw(base_size = 14)+
   ggplot2::theme(
     axis.text.x.bottom = ggplot2::element_text(angle=30,vjust = 1,hjust=1),
     axis.text.x.top = ggplot2::element_text(angle=30,vjust = 0,hjust=0),
     legend.key.size = ggplot2::unit(0.1,"inches"),
     legend.justification = c(1,1),
     legend.position = c(0.95,0.95),
     legend.background = element_blank()
   )
)

doplot = function(p) {
  loc1 = roogledocs::ggplot_to_png(p, width = dims$width*0.8, height = dims$height)
  loc2 = fs::path_ext_remove(loc1) %>% fs::path_ext_set("clip.png")
  system(sprintf("convert %s -resize %s %s", loc1, sprintf("%1.0f%%",96/300*100), loc2))
  system(sprintf("xclip -selection clipboard -t image/png -i %s",loc2))
  # loc2 %>% viewer()
  knitr::include_graphics(loc1)
}


geom_gamma = function(colour,mean=shape/rate,sd=sqrt(shape/rate^2),shape = mean^2/sd^2,rate = mean/sd^2,position_x=0,position_y=0) {
  
  f = function(x) dgamma(x,shape=shape,rate=rate)
  q = function(p) qgamma(p,shape=shape,rate=rate)
  
  shape = mean^2/sd^2
  rate = mean/sd^2
  
  data = tibble::tibble(
    x = seq(0,mean+sd*3,length.out=1001)
  ) %>% dplyr::mutate(
    ymax=f(x)+position_y,
    ymin=position_y,
    x=x+position_x
  )
  
  label = sprintf("mean: %1.2g, sd: %1.3g;\n\u03B1: %1.2g, \u03B2: %1.2g; \u03BA: %1.3g",mean,sd,shape,rate,sd/mean)
  
  return(list(
    ggplot2::geom_line(data=data,mapping=ggplot2::aes(x=x,y=ymax,colour=colour)),
    ggplot2::geom_ribbon(data=data,mapping=ggplot2::aes(x=x,ymax=ymax,ymin=ymin,fill=colour),alpha=0.05),
    ggplot2::geom_segment(ggplot2::aes(x=position_x+mean, y=position_y, xend = position_x+mean, yend=position_y+f(mean), colour=colour)),
    ggplot2::geom_point(ggplot2::aes(x=position_x+q(0.025), y=position_y+f(q(0.025)), colour=colour)),
    ggplot2::geom_segment(ggplot2::aes(x=position_x+q(0.5), y=position_y, xend = position_x+q(0.5), yend=position_y+f(q(0.5)), colour=colour), linetype="dashed"),
    ggplot2::geom_point(ggplot2::aes(x=position_x+q(0.975), y=position_y+f(q(0.975)), colour=colour)),
    ggplot2::geom_label(mapping=ggplot2::aes(colour=colour),x=mean+position_x, y=f(mean)+position_y,label=label,vjust=-0.05,hjust=-0.05, size=14/ggplot2::.pt),
    scale_color_identity(aesthetics=c("fill","colour"))
  ))
}



geom_discrete = function(colour, ymatrix, position_x=0, position_y=0, reversed = FALSE) {
  
  single=!is.matrix(ymatrix)
  
  d = if (is.matrix(ymatrix)) dim(ymatrix) else c(length(ymatrix),1)
  data = tibble::tibble(
    y = as.vector(ymatrix),
    xmin = (rep(1:d[1], d[2])-1)*(1-2*reversed),
    xmax = rep(1:d[1], d[2])*(1-2*reversed),
    group = unlist(lapply(1:d[2],rep,d[1]))
  )
  
  if (single) {
    return(list(
      if (!reversed) {
        ggplot2::geom_step(data = data %>% dplyr::bind_rows(tibble::tibble(y=0,xmin=max(data$xmax))), mapping=ggplot2::aes(x=xmin+position_x, y = y+position_y, colour=colour), alpha=1)
      } else {
        ggplot2::geom_step(data = data %>% dplyr::bind_rows(tibble::tibble(y=0,xmin=min(data$xmax))), mapping=ggplot2::aes(x=xmin+position_x, y = y+position_y, colour=colour),direction = "vh", alpha=1)
      },
      # ggplot2::geom_segment(data = data, mapping=ggplot2::aes(x=xmin+position_x, xend=xmax+position_x, y = y+position_y, yend = y+position_y), colour=colour, alpha=1),
      ggplot2::geom_rect(data = data, mapping=ggplot2::aes(xmin=xmin+position_x, ymin=position_y, xmax = xmax+position_x, ymax=y+position_y, fill=colour), colour=NA, alpha=0.2),
      scale_color_identity(aesthetics=c("fill","colour"))
    ))
  } else {
    return(list(
      ggplot2::geom_segment(data = data, mapping=ggplot2::aes(x=xmin+position_x, xend=xmax+position_x, y = y+position_y, yend = y+position_y), colour="black", alpha=0.25),
      ggplot2::geom_boxplot(data = data, mapping=ggplot2::aes(x=(xmin+xmax)/2+position_x, y=y+position_y, group=as.factor(xmin), colour=colour),outlier.shape = NA),
      scale_color_identity(aesthetics=c("fill","colour"))
    ))
  }
}


# hist(rtnorm_epiestim(10000,5,3,2,10),breaks = 100)
rtnorm_epiestim = function(n, mean, sd, low, high) {
  x=rep(-Inf,n)
  while (any(x<low | x>high)) {
      x = ifelse(x<low | x>high, rnorm(n,mean,sd), x)
  }
  return(x)
  
  # mean_si_sample <- rep(-1, config$n1)
  # std_si_sample <- rep(-1, config$n1)
  # for (k in seq_len(config$n1)) {
  #   while (mean_si_sample[k] < config$min_mean_si || 
  #          mean_si_sample[k] > config$max_mean_si) {
  #     mean_si_sample[k] <- rnorm(1, mean = config$mean_si, 
  #                                sd = config$std_mean_si)
  #   }
  #   while (std_si_sample[k] < config$min_std_si || 
  #          std_si_sample[k] > config$max_std_si) {
  #     std_si_sample[k] <- rnorm(1, mean = config$std_si, 
  #                               sd = config$std_std_si)
  #   }
  # }
}



# ggplot2::ggplot() + geom_discrete(si_sample_data, colour="red")
# ggplot2::ggplot() + geom_discrete(si_sample_data_mean, colour="blue")
# 
# mean = sum(EpiEstim::discr_si(0:30,5,3)*0:30)
# sd = sqrt(sum((0:30*EpiEstim::discr_si(0:30,5,3) - mean)^2/31))



geom_incidence = function(colour, x, i, dates = c(NA,NA)) {
  data = tibble::tibble(x=x,i=i)
  dates = as.numeric(as.Date(dates))
  min_date = as.Date(max(c(dates[1],min(x)),na.rm = TRUE),origin="1970-01-01")
  max_date = as.Date(min(c(dates[2],max(x)),na.rm = TRUE),origin="1970-01-01")
  data= data %>% filter(x >= min_date & x <= max_date)
  
  return(list(
    ggplot2::geom_bar(data=data, mapping=ggplot2::aes(x=x,y=i,fill=colour), stat="identity", colour="black", alpha =0.2),
    scale_color_identity(aesthetics=c("fill","colour")),
    xlab(NULL),
    ylab("count")
  ))
}


geom_rt = function(colour, epiestim_out, dates = c(NA,NA)) {
  data = epiestim_out$R %>% mutate(date = epiestim_out$dates[t_end])
  dates = as.numeric(as.Date(dates))
  min_date = as.Date(max(c(dates[1],min(data$date)),na.rm = TRUE),origin="1970-01-01")
  max_date = as.Date(min(c(dates[2],max(data$date)),na.rm = TRUE),origin="1970-01-01")
  data= data %>% filter(date >= min_date & date <= max_date)
  return(list(
    ggplot2::geom_line(data=data,mapping=ggplot2::aes(x=date, y=`Mean(R)`, colour=colour)),
    ggplot2::geom_ribbon(data=data,mapping=ggplot2::aes(x=date, ymin=`Quantile.0.25(R)`,ymax=`Quantile.0.75(R)`, fill=colour), colour=NA, alpha=0.2),
    ggplot2::geom_ribbon(data=data,mapping=ggplot2::aes(x=date, ymin=`Quantile.0.025(R)`,ymax=`Quantile.0.975(R)`, fill=colour), colour=NA, alpha=0.1),
    ggplot2::geom_hline(yintercept = 1, colour="grey50"),
    scale_color_identity(aesthetics=c("fill","colour")),
    xlab(NULL),
    ylab(expression(R[t]))
  ))
}

# plot_incidence_rt = function(ts, rt_estimate, dates = c(NA,NA), date_breaks = "1 week",date_labels = "%d %b") {
#   p1 = ggplot()+
#     geom_incidence("black",ts$dates, ts$I, dates)+
#     coord_cartesian(xlim=as.Date(dates))+
#     scale_x_date(date_breaks = "1 week",date_labels = "%d %b")+
#     theme(axis.text.x.bottom = element_blank())
#   
#   p2 = ggplot()+
#     geom_rt("black",rt_estimate, dates)+
#     coord_cartesian(ylim=c(0.5,1.75), xlim=as.Date(dates))+
#     scale_x_date(date_breaks = "1 week",date_labels = "%d %b")
#   
#   p1+p2+patchwork::plot_layout(ncol = 1)
# }

# Renewal equation

latex_labels = function(breaks) {
  tmp = lapply(breaks, function(b) latex2exp::TeX(sprintf("$%s$",b)))
  names(tmp) = breaks
  return(tmp)
}

plotmath_labels = function(breaks) {
  tmp = lapply(breaks, function(b) parse(text=b))
  names(tmp) = breaks
  return(tmp)
}
