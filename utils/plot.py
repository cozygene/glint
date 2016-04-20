from numpy import std, log10, linspace, sqrt, ceil
import matplotlib.pyplot as plot
from utils import pca
import logging

def draw_setup(function):
    """
    decorator that runs only on functions that belong to Plot class or it's derivatives
    those functions need to handle the plot drawing and titling

    the decorator opens a new subplot (if there is more than one plot to draw on the same screen)
    then returns to the drawing function (the calling function)
    and finally saves the plots if it draw the last one 
    """
    def wrapper(self, *args, **kwargs):
        if self.plots_number > 1:
            # open a subplot
            plot.subplot(self.nrows, self.ncols, self.current_draw_index)
        
        #run function
        output = function(self, *args, **kwargs)

        self.current_draw_index += 1

        # if finished drawing all plots, save plots
        if self.current_draw_index -1 == self.plots_number: # draw all plots
            if self.save_file:
                plot.savefig(self.save_file)
            # plot.show()

        return output
        
    return wrapper


class Plot(object):
    def __init__(self, save_file=None, plots_number=1):

        assert plots_number > 0
        self.save_file = save_file

        self.plots_number = plots_number
        self.current_draw_index = 1
        if plots_number > 1:
            self.nrows = ceil(sqrt(plots_number))
            self.ncols = ceil(self.plots_number / self.nrows)
            # create windows
            fig, axes = plot.subplots()
            fig.subplots_adjust(hspace=.5)


    def add_title(self, title = None, xtitle = None, ytitle = None):
        if xtitle:
            plot.xlabel(xtitle)
        if ytitle:
            plot.ylabel(ytitle)
        if title:
            plot.title(title)



class QQPlot(Plot):
    def __init__(self, save_file=None, plots_number=1):
        super(QQPlot, self).__init__(save_file, plots_number)

    # Generates a QQ-plot for a given vector of p-values.
    @draw_setup
    def draw(self, y, title = None, xtitle = None, ytitle = None, style = 'b.'):
        #x
        x = -log10(linspace(0,1,len(y)+1))
        x = x[1:]
        x.sort()

        # y
        y = -log10(y)
        y.sort()

        # qqplot
        plot.plot(x, y, style)

        # add y=x trend
        xlim = int(y.max()) + 2
        plot.plot(range(xlim), range(xlim), 'r-') 

        # axis limit
        plot.xlim(0, xlim)
        plot.ylim(0, round(y.max()))

        self.add_title(title, xtitle, ytitle)



class PCAScatterPlot(Plot):
    STD_LINES_SYMMETRIC_LIMIT = 6

    def __init__(self, pca_out, save_file=None, plots_number=1):
        if plots_number  >= pca_out.P.shape[1] - 1:
            logging.warning("there are no %s pcs to show, will show all existing pcs" % plots_number)
            plots_number = pca_out.P.shape[1] - 1

        super(PCAScatterPlot, self).__init__(save_file, plots_number)
        self.pca_out = pca_out

    def draw(self):
        for i in range(self.plots_number):
            self.draw_std_outliers(x = self.pca_out.P[:,i], y = self.pca_out.P[:,i+1], title = "pc%s vs pc%s"%(i+1, i+2), xtitle = "pc%s"%(i+1), ytitle = "pc%s"%(i+2))


    @draw_setup
    def draw_std_outliers(self, x, y, title=None, xtitle=None, ytitle=None):
        plot.scatter(x, y)
        # the top title will be set seperatly
        # this must be called before ax2 initiation
        self.add_title(xtitle = xtitle, ytitle = ytitle)

        # add vertical lines for each std to easly see outliers
        x_std = std(x)
        lines = range(-1* self.STD_LINES_SYMMETRIC_LIMIT, self.STD_LINES_SYMMETRIC_LIMIT + 1)
        lines.remove(0) 
        
        # draw the lines
        [plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]
        
        # set lines titles on second x axis
        # this is a patch that works (I dont understand why)
        
        left,right = plot.xlim()   # xlim must be taken before calling plot.twiny(), since calling to twiny() changes the limit
        ax2 = plot.twiny()         # set another x axis
        ax2.set_xlim(left, right)  # set xlim to the one before duplicatin x-axis: that line is what causes the std-lines title to adjust window size changing
        
        ax2.set_xticks([x_std*i for i in  lines])           # set ticks in the std places (dont show them with grid since this grid is wrong, that is why we use axvlines above)
        ax2.set_xticklabels(["%d std"%i for i in lines])    # set title for each tick
       
        plot.text(0.5, 1.11, title,
             horizontalalignment='center',
             fontsize=14,
             transform = ax2.transAxes)

        
        # it's also possible to draw using "ax2.xaxis.grid()"
            # but it doesn't work (lines are drawn bad) without the "[plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]".
            # moreover it only works if the axvline line is called before the .set_xticks line. otherwise it won't work.
            # I also noticed that if I set the xlim (out.set_xlim) it works with or without "[plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]"
            # so I think that changing the limit before calling xticks is what makes it work (plot.axvline makes the limit wider)
            # Elior can you help and check if the plot is OK?

