from numpy import std, log10, linspace, sqrt, ceil
import matplotlib.pyplot as plot
from utils import pca


def draw_setup(function):
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
            plot.show()

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
        super(PCAScatterPlot, self).__init__(save_file, plots_number)
        self.pca_out = pca_out

    def draw(self):
        for i in range(self.plots_number + 1):
            # plot.subplot(nrows, ncols,i +1)
            self.draw_std_outliers(x = self.pca_out.P[:,i], y = self.pca_out.P[:,i+1], title = "pc%s vs pc%s"%(i+1, i+2), xtitle = "pc%s"%(i+1), ytitle = "pc%s"%(i+2))
            # self.add_title()

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
        # left,right = out.get_xlim()
        ax2 = plot.twiny()
        # ax2.set_xlim(left, right)
        ax2.set_xticks([x_std*i for i in  lines])
        ax2.set_xticklabels(["%d std"%i for i in lines])
        
        plot.text(0.5, 1.11, title,
             horizontalalignment='center',
             fontsize=14,
             transform = ax2.transAxes)

        
        # TODO - it's also possible to draw using "ax2.xaxis.grid()"
            # but it doesn't work (lines are drawn bad) without the "[plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]".
            # moreover it only works if the axvline line is called before the .set_xticks line. otherwise it won't work.
            # I also noticed that if I set the xlim (out.set_xlim) it works with or without "[plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]"
            # so I think that changing the limit before calling xticks is what makes it work (plot.axvline makes the limit wider)
            # Elior can you help and check if the plot is OK?





# # Generates a QQ-plot for a given vector of p-values.
# def draw_qqplot(y, title, xtitle, ytitle, style = 'b.', save_file = None):
#     # x
#     x = -log10(linspace(0,1,len(y)+1))
#     x = x[1:]
#     x.sort()

#     # y
#     y = -log10(y)
#     y.sort()

#     plot.plot(x, y, style)

#     # add y=x trend
#     xlim = int(y.max()) + 2
#     plot.plot(range(xlim), range(xlim), 'r-') 

#     # axis limit
#     plot.xlim(0, xlim)
#     plot.ylim(0, round(y.max()))

#     # titles
#     plot.xlabel(xtitle)
#     plot.ylabel(ytitle)
#     plot.title(title)

#     if save_file:
#         plot.savefig(save_file)

# def draw_multiple_plots(number_pf_plots, plot_function, save_file = None, **kwargs):
#     # calc the plot cells number
#     nrows = ceil(sqrt(number_pf_plots))
#     ncols = ceil(sqrt(number_pf_plots))+1 # cannot use cell no.0 so we need to add more cells 
#     # create windows
#     fig, axes = plot.subplots()
#     fig.subplots_adjust(hspace=.5)
    
#     # draw plots from window index 1
#     for i in range(number_pf_plots + 1):
#         plot.subplot(nrows, ncols,i +1)
#         plot_function(**kwargs)
    
#     if save_file:
#         # mng = plot.get_current_fig_manager()
#         # mng.frame.Maximize(True)
#         plot.savefig(save_file, dpi=1)



# STD_LINES_SYMMETRIC_LIMIT = 6
# def draw_scatter_plot(x, y, title=None,  xtitle=None, ytitle=None, save_file=None):
#     plot.scatter(x,y)
#     if xtitle:
#         plot.xlabel(xtitle)
#     if ytitle:
#         plot.ylabel(ytitle)
#     if title:
#         plot.title(title)

#     if save_file:
#         plot.savefig(save_file)




# def coords(s):
#     try:
#         x, y = map(int, s.split(','))
#         return x, y
#     except:
#         raise argparse.ArgumentTypeError("Coordinates must be x,y")

# # # parser.add_argument('--cord', help="Coordinate", dest="cord", type=coords, nargs='+') # + - expects at least one argument
# args = parser.parse_args(["--maxpcstd 1,2 2,4 3,6"])

# # parser.add_argument('--cord', help="Coordinate", dest="cord", type=coords, action='append') # + - expects at least one argument
# args = parser.parse_args(["--maxpcstd 1,2 --maxpcstd 2,4 --maxpcstd 3,6"])

# parser.add_argument('--cord', help="Coordinate", dest="cord", action='append', type=int, nargs =2) # 2 - expects at least one argument
# args = parser.parse_args(['--maxpcstd 1 2 --maxpcstd 3 6'])#["--cord","1,2","--cord","2,4","--cord","3,6"])

# print args.cord


# #----------------- Rdata to numpy
# import rpy2.robjects as robjects
# import numpy as np

# # load your file
# robjects.r['load'](RData_file_path)

# # retrieve the matrix that was loaded from the file
# matrix = robjects.r[RData_argument_name]

# # turn the R matrix into a numpy array
# a = np.array(matrix)
