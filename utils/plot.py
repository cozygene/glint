from numpy import std, linspace, sqrt, ceil, log10
import matplotlib.pyplot as plot
# import matplotlib
import logging
from pandas import DataFrame
from scipy.stats import uniform
from scipy.stats import randint
from itertools import cycle
from math import ceil
from utils import tools

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

        output = function(self, *args, **kwargs)

        self.current_draw_index += 1

        # if finished drawing all plots, save plots
        if self.current_draw_index -1 == self.plots_number: # draw all plots
            if self.save_file:
                logging.info("Saving plot to {filename}.png and {filename}.eps".format(filename =self.save_file))
                plot.savefig(self.save_file + ".png" , dpi = 300) # low-quality, can be changed with dpi =600
                plot.savefig(self.save_file + ".eps", format='eps') # has no dpi param
                # formats supported eps, pdf, pgf, png, ps, raw, rgba, svg, svgz
                # png has a dpi but eps and svg no (vector)
            plot.close()

        return output
        
    return wrapper

class Plot(object):
    """
    testes on ubuntu, macOS and windows
    """
    def __init__(self, save_file=None, plots_number=1):

        assert plots_number > 0
        self.save_file = save_file


        self.plots_number = plots_number
        self.current_draw_index = 1
        if plots_number > 1:
            self.nrows = ceil(sqrt(plots_number))
            self.ncols = ceil(self.plots_number / self.nrows)
            # create windows
            self.fig, axes = plot.subplots(figsize = (self.ncols*8,self.nrows*4))
            self.fig.subplots_adjust(hspace=.5)

        # else: this is not woring on ubuntu
        #     self.fig, axes = plot.subplots()
        #     pass

  



    def add_title(self, title = None, xtitle = None, ytitle = None):
        if xtitle:
            plot.xlabel(xtitle)
        if ytitle:
            plot.ylabel(ytitle)
        if title:
            plot.title(title, y = 1.08)



class QQPlot(Plot):
    X_LABEL = "Expected -log(p-values)"
    Y_LABEL = "Observed -log(p-values)"

    def __init__(self, save_file=None, plots_number=1):
        super(QQPlot, self).__init__(save_file, plots_number)

    # Generates a QQ-plot for a given vector of p-values.
    @draw_setup
    def draw(self, y, title = None, xtitle = None, ytitle = None, style = 'b.'):
        #x
        x = tools.minusLog10(linspace(0.001,1.001,len(y)+1)) # 0.001 and 1.001 instead of 0 and 1 in order to avoid dividing by zero warning
        x = x[1:]
        x.sort()

        # y
        y = tools.minusLog10(y) #changes y in place
        y.sort()

        # qqplot
        plot.plot(x, y, style)

        # add y=x trend
        xlim = int(ceil(x.max()) + 2)
        ylim = int(ceil(y.max()) + 2)
        xy = max(ylim, xlim)
        plot.plot(range(xy+1) , range(xy+1), 'r-') 

        # axis limit
        plot.xlim(0, xlim)
        plot.ylim(0, ylim)

        if ytitle is None:
            ytitle = self.Y_LABEL
        if xtitle is None:
            xtitle = self.X_LABEL

        ax = plot.axes()
        pos1 = ax.get_position() # get the original position 
        pos2 = [pos1.x0 + 0.02, pos1.y0 + 0.1,  pos1.width/1.1, pos1.height/1.3] 
        ax.set_position(pos2)
        self.add_title(title, xtitle, ytitle)



class PCAScatterPlot(Plot):
    STD_LINES_SYMMETRIC_LIMIT = 6

    def __init__(self, pca_out, save_file=None, plots_number=1):
        if plots_number  >= pca_out.P.shape[1] - 1:
            logging.warning("There are no %s PCs to show, will show all existing pcs." % plots_number)
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
        self.add_title(title, xtitle = xtitle, ytitle = ytitle)

        # add vertical lines for each std to easly see outliers
        x_std = std(x)
        lines = range(-1* self.STD_LINES_SYMMETRIC_LIMIT, self.STD_LINES_SYMMETRIC_LIMIT + 1)
        lines.remove(0) 
        lines.remove(1) 
        lines.remove(-1) 
        
        # draw the lines
        [plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]
        
        # set lines titles on second x axis
        # this is a patch that works (I dont understand why)
        
        left,right = plot.xlim()   # xlim must be taken before calling plot.twiny(), since calling to twiny() changes the limit
        ax2 = plot.twiny()         # set another x axis
        ax2.set_xlim(left, right)  # set xlim to the one before duplicatin x-axis: that line is what causes the std-lines title to adjust window size changing
        
        ax2.set_xticks([x_std*i for i in  lines])           # set ticks in the std places (dont show them with grid since this grid is wrong, that is why we use axvlines above)
        ax2.set_xticklabels(["%dsd"%i for i in lines])    # set title for each tick
       
        
        # it's also possible to draw using "ax2.xaxis.grid()"
            # but it doesn't work (lines are drawn bad) without the "[plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]".
            # moreover it only works if the axvline line is called before the .set_xticks line. otherwise it won't work.
            # I also noticed that if I set the xlim (out.set_xlim) it works with or without "[plot.axvline(x=x_std*i, color='r',linestyle='-') for i in  lines]"
            # so I think that changing the limit before calling xticks is what makes it work (plot.axvline makes the limit wider)
            # Elior can you help and check if the plot is OK?



class ManhattanPlot(Plot):
    Y_LABEL = "-log(p-value)"

    def __init__(self, save_file=None, plots_number=1):
        super(ManhattanPlot, self).__init__(save_file, plots_number)

    @draw_setup
    def draw(self, sites, pvalues, chromosomes, positions, title = None, xtitle = None, ytitle = None, style = 'b.'):

        ax = plot.axes()
        self.manhattan(sites, pvalues, chromosomes, positions, ax)

        if ytitle is None:
            ytitle = self.Y_LABEL
        self.add_title(title, xtitle, ytitle)

        # option #1
        pos1 = ax.get_position() # get the original position 
        pos2 = [pos1.x0 + 0.02, pos1.y0 + 0.1,  pos1.width/1.1, pos1.height/1.6] 
        ax.set_position(pos2) # set a new position

    def manhattan(self, sites, pvalues, chromosomes, positions, ax):
        # thiese two lines make sure that chromosomes will be sorted by int value and not string value
        # (chromosomes are list from the values [1,2,..., X, Y] )
        chromosomes = [int(ch) if ch.isdigit() else  ch for ch in chromosomes]
        all_chromosomes = list(set(chromosomes))
        all_chromosomes.sort() # all the chromosomes sorted
        
        number_of_sites = len(sites)
        number_of_chr = len(all_chromosomes)
        space_for_each_chr = float(number_of_sites/number_of_chr)

        df = DataFrame({'sites' : sites,
                        'minuslog10pvalue' : tools.minusLog10(pvalues),
                        'positions' : positions,
                        'chromosome' : chromosomes})
        
        df.chromosome = df.chromosome.astype('category')
        df.chromosome = df.chromosome.cat.set_categories(all_chromosomes, ordered=True)
        
        df['ind'] = range(len(df))
        df_grouped = df.groupby(('chromosome'))
        
        colors = cycle(['darkblue','skyblue'])
        x_labels = []
        for (name, group) in df_grouped:
            clr = colors.next()
            # the position of this chromosome on the x-asix is according to its index at the sorted chromosomes name list
            x_pos = all_chromosomes.index(name) 
            # x-axis is determined by the position of the chromosome (the positions are normalized to the space
            # that is assigned to the chromosome) - that way, the plot is grouped by chromosome and than by sorted position
            # y-axis is the -log10(p-value)
            x = group['positions']%space_for_each_chr  + x_pos*space_for_each_chr
            ax.scatter(x, group['minuslog10pvalue'], color=clr, edgecolors='none' )
            # # to set the x-axis space of each chromoseme by the amount of sites at this chromosome (more space for a chromosome
            # # with more sites) - do something like this:
            # # that option is good when you want to see the samples more clear at the wider chromosomes)
            # group.plot(kind='scatter', x='ind', y='minuslog10pvalue',color=clr, ax=ax, edgecolors='none')
            # x_labels_pos.append((group['ind'].iloc[-1] - (group['ind'].iloc[-1] - group['ind'].iloc[0])/2))

            x_labels.append('ch%s' % name)

        ax.axhline(y=-log10(0.05/number_of_sites), color='r',linestyle='-')
        
        ax.set_xticks([space_for_each_chr*i+space_for_each_chr/2 for i in range(number_of_chr)])
        ax.set_xticklabels(x_labels,  rotation='vertical')
        ax.set_xlim([0, len(df)])
        ax.set_ylim([0, max(df.minuslog10pvalue) + 1.5])

