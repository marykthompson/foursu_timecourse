#Define colors and plot helper functions
import matplotlib.pyplot as plt
import seaborn as sns
import yaml

plt.style.use('../scripts/style_sheets/Ribopop.mplstyle')

#config file copied from snakemake will be here:
with open('../thisconfig.yml', 'r') as file:
    config = yaml.load(file, Loader=yaml.FullLoader)

results_dir = config['results_dir']
probe_designer_dir = config['probe_designer_dir']
gffutils_db = config['gffutils_db']

carto12_hex = ['#7F3C8D', '#11A579', '#3969AC', '#F2B701', '#E73F74', '#80BA5A', '#E68310',
'#008695', '#CF1C90', '#f97b72', '#4b4b8f', '#A5AA99']

#https://gist.github.com/JoachimGoedhart/b2f91652de2b9e3b393c6c28be843e00
#these are the tol_muted colors
color_dict = {'indigo': '#332288', 'cyan':'#88CCEE', 'teal':'#44AA99', 'green':'#117733',
'olive':'#999933', 'sand':'#DDCC77', 'rose':'#CC6677', 'wine':'#882255', 'purple':'#AA4499',
 'pale_grey':'#DDDDDD', 'grey':'#BBBBBB'}

ordered_colors = [color_dict[i] for i in ['grey', 'purple', 'cyan', 'indigo', 'olive', 'rose', 'sand', 'wine', 'teal', 'green']]
#tol_muted = ['#332288', '#88CCEE', '#44AA99', '#117733', '#999933', '#DDCC77', '#CC6677', '#882255', '#AA4499', '#DDDDDD']

##selected_colors = carto12_hex
selected_colors = ordered_colors

#put the grey color first
#selected_colors.reverse()
sns.set_palette(sns.color_palette(selected_colors))

#Set up the figsize and axis position
x_size = 3.5
y_size = 3.5

#column widths, converting to inches for matplotlib
scol = 84/25.4
dcol = 178/25.4

sfig = scol/2
dfig = dcol/4
#I think generally the subpanel will be a quarter of page width, so
#single_col/2 or double_col/4

figsizes = {'single_col': (sfig, sfig), 'double_col': (dfig, dfig)}

#For some reason ax.text is not respecting font.size in the rc params file, so set it here.
label_fontsize = 8

class Plotter(object):
    def __init__(self, figtype = 'single_col', corners = [0.2, 0.2, 0.75, 0.75], figsize = None):
        if figsize:
            self.figsize = figsize
        else:
            self.figsize = figsize[figtype]
        #keep track of original corner position so that we can calculate differences later
        self.original_corners = corners
        self.corners = corners

    def nudge_corners(self, top = False, bottom = False, left = False, right = False, tb_in = 0.03, lr_in = 0.03):

        x0, y0, width, height = self.corners

        #get the x position of the label. Because it is scaled to the figsize, need to check
        size_x, size_y = self.figsize
        xfrac = lr_in/size_x
        yfrac = tb_in/size_y

        #adjust left/right
        #if right, only need to adjust the width of the axes
        if right:
            width = width - xfrac
        if left:
            x0 = x0 + xfrac
            width = width - xfrac
        if top:
            height = height - yfrac
        if bottom:
            y0 = y0 + yfrac
            height = height - yfrac

        self.corners = [x0, y0, width, height]

    def setup_axis(self, new_fig = True):
        '''
        Set up an axes. Normally this will be as part of a new figure.
        If new_fig == True, use the current figure to add the axis.
        '''
        if new_fig:
            self.fig = plt.figure(figsize = self.figsize)
        self.ax = self.fig.add_axes(self.corners)

    def set_xlabel(self,text, pos = None, nudge = None):
        '''
        Label x-axis with respect to the figure coordinates
        '''
        self.ax.set_xlabel('')
        if pos:
            xpos, ypos = pos
        else:
            axpos = self.ax.get_position()
            xpos = axpos.x0 + axpos.width/2
            ##try putting the x-axis label a little higher
            ypos = 0.07

        if nudge:
            xpos = xpos + nudge[0]
            ypos = ypos + nudge[1]

        lpos = (xpos, ypos)
        self.xlabel_pos = lpos
        self.ax.text(*lpos, text, horizontalalignment='center', verticalalignment='center',
        transform = self.fig.transFigure, fontsize = label_fontsize)


    def set_ylabel(self, text, pos = None, nudge = None):

        '''
        Label y-axis with respect to the figure coordinates
        Nudge provides a way to move the axis label slightly to accommodate the
        subpanel label, generally down, i.e. a negative number to move down slightly.
        '''
        self.ax.set_ylabel('')
        if pos:
            xpos, ypos = pos
        else:
            axpos = self.ax.get_position()
            ypos = axpos.y0 + axpos.height/2

            #assume if not sure, x position needs to be scaled
            size_x, size_y = self.fig.get_size_inches()
            if size_x != size_y:
                sfactor = size_y/size_x
            else:
                sfactor = 1

            xpos = 0.04*sfactor

        if self.original_corners != self.corners:
            x_diff = self.corners[0] - self.original_corners[0]
            xpos = xpos + x_diff

        if nudge:
            ypos = ypos + nudge[1]
            xpos = xpos + nudge[0]
        lpos = (xpos, ypos)
        self.ylabel_pos = lpos
        self.ax.text(*lpos, text, horizontalalignment='center', verticalalignment='center',
        transform = self.fig.transFigure, fontsize = label_fontsize, rotation = 90)

    def add_letter(self, letter, ha = 'center'):
        '''
        Add the subpanel letter at the top left, i.e. A, B, C
        ha (horizontalalignment) is best at center for single line labels,
        but better as right for multi-line labels.
        '''
        if self.original_corners != self.corners:
            y_diff = self.corners[-1] - self.original_corners[-1]
        else:
            y_diff = 0
        pos = self.ylabel_pos[0], 0.95 + y_diff
        #try changing horizontalalignment to right to accommodate multiline labels
        self.ax.text(*pos, letter, horizontalalignment = ha, verticalalignment='center',
        transform = self.fig.transFigure, fontsize = 11, fontweight = 'bold')
