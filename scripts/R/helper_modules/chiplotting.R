library(evd)
library(ggplot2)
library(dplyr)
library(tidync)

plotChi <- function(netcdf_path, save_path, cols, ltype = 'twodash', lsize = 0.7, 
                    confreg_lsize = 0.7, confreg_col = 'blue', confreg_alpha = 0.1, 
                    title = 'chi, estimate + 95% confidence band', fsize_title = 17, fsize_axlabs = 20, fsize_axticks = 15) {

    # Function to plot the chi statistic along with 95% confidence band, as a function of quantile
    #
    # Arguments:
    #    netcdf_path, save_path: strings representing path + filename for the netcdf file and where you want to save the resulting plot
    #	 cols: a vector of indices corresponding to columns of interest in dataframe
    #    ltype: type of line to use for estimate, see ggplot documentation for alternatives
    #    lsize, confreg_lsize: thickness of estimate line, as well as lines bounding confidence region
    #    confreg_col: color of shaded confidence band
    #    confreg_alpha: opacity of shaded confidence band
    #    title: string, simply the title
    #    fsize_title, fsize_axlabs, fsize_axticks: font size of plot title, axis labels, and numbers on axis ticks
   
    # convert netcdf file to a dataframe	
    dat <- tidync(netcdf_path)
    dat <- dat %>% hyper_tibble()
    # extract relevant columns
    dat <- dat[,cols]
    
    # create chiplot, but turn graphics off so we can make our own plot
    chiplot_res <- chiplot(dat, which = 1)
    graphics.off()

    # make a dataframe of output from chiplot
    plt_df <- data.frame(chiplot_res[[2]])
    plt_df <- plt_df %>% mutate(quantile = chiplot_res[[1]])
    
    # make and save plot
    plt <- ggplot(plt_df) + geom_line(aes(x = quantile, y = chi), linetype = ltype, size = lsize) + 
        geom_line(aes(x = quantile, y = chiupp), size = confreg_lsize) + geom_line(aes(x = quantile, y = chilow), size = confreg_lsize) + 
        geom_ribbon(aes(x = quantile, ymin = chilow, ymax = chiupp), fill = confreg_col, alpha = confreg_alpha) +
        coord_cartesian(ylim=c(-1, 1), xlim = c(0, 1)) + ggtitle(title) + theme_light() + 
        theme(axis.title = element_text(size = fsize_axlabs), axis.text = element_text(size = fsize_axticks), text = element_text(size = fsize_title))
    
    ggsave(save_path, plt)
}

checkNetCDFCols <- function(path) {

    # Function that returns the column names for a data frame extracted from a netcdf file.
    # To be used to find indices corresponding to columns of interest, and then passed as an arg to plotChi.
    #
    # path: string, representing path + filename for netcdf file of interest


    dat <- tidync(path)
    dat <- dat %>% hyper_tibble()
    return(colnames(dat))

}
