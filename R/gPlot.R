#' @title gPlot
#' @description
#' Plot maximized \eqn{g(\bold{\beta},\theta,\sigma|\bold{\beta^{(t)}},\theta^{(t)},
#' \sigma^{(t)})} over iterations at a specified spike parameter.
#' The function \eqn{g} is a concave lower bound for objective log-likelihood function at
#' \eqn{\bold{\beta^{(t)}}}, \eqn{\theta^{(t)}} and \eqn{\sigma^{(t)}} by Jensen's Inequality.
#' \eqn{\bold{\beta^{(t+1)}}}, \eqn{\theta^{(t+1)}} and \eqn{\sigma^{(t+1)}}
#' are obtained by maximizing \eqn{g(\bold{\beta},\theta,\sigma|\bold{\beta^{(t)}},\theta^{(t)},
#' \sigma^{(t)})}. See Vignette for more details.
#'
#'
#' @param fit_obj Fitted object from function 'lsum'.
#' @param spike_param The value of the spike parameter to be specified, which must be in
#' \code{spike_params} designated in 'lsum' function. Each \code{spike_param} corresponds
#' to a unique gPlot. If default \code{spike_params} is used in function 'lsum', one can
#' check the value of \code{spike_params} by calling `$spike_params` from the returned
#' object from function 'lsum'.
#' @export

gPlot = function(fit_obj, spike_param){
  
  spike_index = which(fit_obj$spike_params == spike_param)
  g_vals = fit_obj$g_List[[spike_index]]
  iter_num = fit_obj$iter_nums[spike_index]
  
  theta_path = fit_obj$theta_path
  ggplot()+
    geom_line(aes(x = 1:iter_num, y = g_vals, colour = "g"),size = 1)+
    geom_point(aes(x = 1:iter_num, y = g_vals, colour = "g"),size = 3)+
    xlab("iteration")+ylab("g")+
    labs(title ="g over iterations")+
    theme(plot.title = element_text(size=13,hjust=0.5))+
    scale_color_manual("",values = c("g" = "lightblue3"))
}