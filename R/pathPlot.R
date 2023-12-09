#' @title pathPlot
#' @description
#' Plot the path of \code{beta}, \code{theta} or \code{sigma} over \code{spike_params}.
#'
#' @param fit_obj The fitted object from 'lsum' function.
#' @param path_of Character. \code{path_of} = \code{"beta"}, \code{"theta"} or \code{"sigma"}.
#' Default is \code{"beta"}.
#' @import ggplot2
#' @import reshape2
#' @export

pathPlot = function(fit_obj, path_of = "beta"){
  x_axis = log(fit_obj$spike_params)
  if (path_of == "beta"){
    indices = fit_obj$beta_indices
    beta_path = fit_obj$beta_path
    p = dim(beta_path)[2]
    Nonzero_beta = data.frame(cbind(x_axis, beta_path[,indices]))
    Zero_beta = data.frame(cbind(x_axis, beta_path[,-indices]))
    ymax = max(beta_path)
    ymin = min(beta_path)
    
    melt_Nonzero = melt(Nonzero_beta, id=names(Nonzero_beta)[1])
    melt_Zero = melt(Zero_beta, id=names(Zero_beta)[1])
    
    Plot = ggplot(melt_Nonzero,aes(x = x_axis, y = value,group = variable))+
      geom_line(aes(color = "Nonzero",linetype = "Nonzero"),size = 0.5)+
      geom_point(data = melt_Nonzero,aes(x = x_axis, y = value, group = variable,color = "Nonzero"), size = 2)+
      geom_line(data = melt_Zero, aes(x = x_axis, y = value,group = variable, color = "Zero",linetype = "Zero"),size = 0.3)+
      geom_point(data = melt_Zero,aes(x = x_axis, y = value, group = variable,color = "Zero"), size = 2)+ ylim(ymin,ymax)+
      xlab("log(spike_param)")+ylab("beta")+
      labs(title ="Path of beta")+
      theme(plot.title = element_text(size=13,hjust=0.5))
      
    Plot = Plot + scale_color_manual("",values = c("Nonzero" = "brown","Zero" = "grey"))+ scale_linetype_manual("", values = c("Nonzero" = "solid","Zero" = "dashed"))
    Plot
  }else if(path_of == "theta"){
    theta_path = fit_obj$theta_path
    ggplot()+
      geom_line(aes(x = x_axis, y = theta_path, colour = "theta"),size = 1)+
      geom_point(aes(x = x_axis, y = theta_path, colour = "theta"),size = 3)+
      xlab("log(spike_param)")+ylab("theta")+
      labs(title ="Path of theta")+
      theme(plot.title = element_text(size=13,hjust=0.5))+
      scale_color_manual("",values = c("theta" = "lightblue3"))
  }else if(path_of == "sigma"){
    sigma_path = fit_obj$sigma_path
    ggplot()+
      geom_line(aes(x = x_axis, y = sigma_path, colour = "sigma"),size = 1)+
      geom_point(aes(x = x_axis, y = sigma_path, colour = "sigma"),size = 3)+
      xlab("log(spike_param)")+ylab("sigma")+
      labs(title ="Path of sigma")+
      theme(plot.title = element_text(size=13,hjust=0.5))+
      scale_color_manual("",values = c("sigma" = "lightblue3"))
  }else{
    stop("path_of should be chosen from 'beta','theta' and 'sigma'.")
  }
}