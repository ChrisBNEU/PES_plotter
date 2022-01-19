import cantera as ct
import os
import sys
sys.path.append(f'{os.getcwd()}/tools/PyEnergyDiagram')
print(os.getcwd())
from energydiagram import ED
import logging
from matplotlib.lines import Line2D

############################################
#
#   Plots a potential energy surface 
#   (enthalpy vs rxn coordinate) for a 
#   given cti file mechanism
#   
#   uses https://github.com/giacomomarchioro/PyEnergyDiagrams
#
############################################

class single_species():
    """
    takes a single "species" cantera object. 

    Arguments:
    Energy - default is enthalpy, for y placement on energy diagram
    level - x position on energy diagram. integer 0, 1, 2, etc. 
    """
    def __init__(
            self,
            species,
            position,
            reverse=False,
            ):

        self.energy = (species.thermo.h(temp)/1000**2)/96
        self.name = species.name

        self.position = position # 1
        self.top_text = ''
        self.color = 'k'  
        self.right_text = ''
        self.left_text = ''
        self.uid = hash(self)

    
class combined_species():
    """
    makes a grouping of reactants or products. can be used with a single species. 
    """
    def __init__(
            self,
            species,
            position,
            reverse=False,
            ):

        name = ''
        for i in species: 
            self.energy += (species.thermo.h(temp)/1000**2)/96
            self.name = name + ' + ' + species.name

        self.species = species
        self.position = position
        self.top_text = ''
        self.color = 'k'  
        self.right_text = ''
        self.left_text = ''
        self.uid = hash(self)

class reaction():
    """
    makes a reaction object. 
    """
    def __init__(
            self,
            reaction,
            position,
            reverse=False,
            ):
        
                name = ''
        for i in species: 
            self.energy += (species.thermo.h(temp)/1000**2)/96
        equation = species.name

        self.species = species

        self.energy = h_eV  # 0
        self.position = position # 1
        self.bottom_text = name  # 2
        self.top_text = ''
        self.color = 'k'  
        self.right_text = ''
        self.left_text = ''
        self.uid = hash(self)


class ED:
    def __init__(self, aspect='equal'):
        # plot parameters
        self.ratio = 1.6181
        self.dimension = 'auto'
        self.space = 'auto'
        self.offset = 'auto'
        self.offset_ratio = 0.02
        self.color_bottom_text = 'blue'
        self.aspect = aspect
        self.data = []
        self.pos_number = 0
        self.energies = []
        self.positions = []
        self.colors = []
        self.top_texts = []
        self.bottom_texts = []
        self.left_texts = []
        self.right_texts = []
        self.links = []
        self.arrows = []
        self.electons_boxes = []
        # matplotlib fiugre handlers
        self.fig = None
        self.ax = None



    def create_data(self):
        '''
        Method of ED class
        instantiates "data" object, so we do not have to plot each time we need to change that structure

        Returns
        -------
        fig (plt.figure) and ax (fig.add_subplot())

        '''


        self.data = list(zip(self.energies,  # 0
                        self.positions,  # 1
                        self.bottom_texts,  # 2
                        self.top_texts,  # 3
                        self.colors,  # 4
                        self.right_texts,  # 5
                        self.left_texts,))  # 6

        # for idx, link in enumerate(self.links):
        #     # here we connect the levels with the links
        #     # x1, x2   y1, y2
        #     for i in link:
        #         # i is a tuple: (end_level_id,ls,linewidth,color)
        #         start = self.positions[idx]*(self.dimension+self.space)
        #         x1 = start + self.dimension
        #         x2 = self.positions[i[0]]*(self.dimension+self.space)
        #         y1 = self.energies[idx]
        #         y2 = self.energies[i[0]]
        #         line = Line2D([x1, x2], [y1, y2],
        #                       ls=i[1],
        #                       linewidth=i[2],
        #                       color=i[3])




    def plot(
        self, 
        show_IDs=False, 
        ylabel="Energy / $kcal$ $mol^{-1}$", 
        width = 10,
        height = 10,
        ax: plt.Axes = None
        ):
        '''
        Method of ED class
        Plot the energy diagram. Use show_IDs=True for showing the IDs of the
        energy levels and allowing an easy linking.
        E|          4__
        n|   2__    /  \
        e|1__/  \__/5   \
        r|  3\__/       6\__
        g|
        y|

        Parameters
        ----------
        show_IDs : bool
            show the IDs of the energy levels
        ylabel : str
            The label to use on the left-side axis. "Energy / $kcal$
            $mol^{-1}$" by default.
        ax : plt.Axes
            The axes to plot onto. If not specified, a Figure and Axes will be
            created for you.

        Returns
        -------
        fig (plt.figure) and ax (fig.add_subplot())

        '''
        # Create a figure and axis if the user didn't specify them.
        if not ax:
            fig = plt.figure()
            ax = fig.add_subplot(111, aspect=self.aspect)
        
        # Otherwise register the axes and figure the user passed.
        else:
            self.ax = ax
            self.fig = ax.figure

            # Constrain the target axis to have the proper aspect ratio
            self.ax.set_aspect(self.aspect)

        ax.set_ylabel(ylabel)
        ax.axes.get_xaxis().set_visible(False)
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)
        ax.spines['bottom'].set_visible(False)

        self.__auto_adjust()

        fig.set_figwidth(width)
        fig.set_figheight(height)
        self.data = list(zip(self.energies,  # 0
                        self.positions,  # 1
                        self.bottom_texts,  # 2
                        self.top_texts,  # 3
                        self.colors,  # 4
                        self.right_texts,  # 5
                        self.left_texts,))  # 6

        for level in self.data:
            start = level[1]*(self.dimension+self.space)
            ax.hlines(level[0], start, start + self.dimension, color=level[4])
            ax.text(start+self.dimension/2.,  # X
                    level[0]+self.offset,  # Y
                    level[3],  # self.top_texts
                    horizontalalignment='center',
                    verticalalignment='bottom')

            ax.text(start + self.dimension,  # X
                    level[0],  # Y
                    level[5],  # self.bottom_text
                    horizontalalignment='left',
                    verticalalignment='center',
                    color=self.color_bottom_text)

            ax.text(start,  # X
                    level[0],  # Y
                    level[6],  # self.bottom_text
                    horizontalalignment='right',
                    verticalalignment='center',
                    color=self.color_bottom_text)

            ax.text(start + self.dimension/2.,  # X
                    level[0] - self.offset*2,  # Y
                    level[2],  # self.bottom_text
                    horizontalalignment='center',
                    verticalalignment='top',
                    color=self.color_bottom_text)
        if show_IDs:
            # for showing the ID allowing the user to identify the level
            for ind, level in enumerate(self.data):
                start = level[1]*(self.dimension+self.space)
                ax.text(start, level[0]+self.offset, str(ind),
                        horizontalalignment='right', color='red')

        for idx, arrow in enumerate(self.arrows):
            # by Kalyan Jyoti Kalita: put arrows between to levels
            # x1, x2   y1, y2
            for i in arrow:
                start = self.positions[idx]*(self.dimension+self.space)
                x1 = start + 0.5*self.dimension
                x2 = start + 0.5*self.dimension
                y1 = self.energies[idx]
                y2 = self.energies[i]
                gap = y1-y2
                gapnew = '{0:.2f}'.format(gap)
                middle = y1-0.5*gap  # warning: this way works for negative HOMO/LUMO energies
                ax.annotate("", xy=(x1, y1), xytext=(x2, middle), arrowprops=dict(
                    color='green', width=2.5, headwidth=5))
                ax.annotate(s=gapnew, xy=(x2, y2), xytext=(x1, middle), color='green', arrowprops=dict(width=2.5, headwidth=5, color='green'),
                            bbox=dict(boxstyle='round', fc='white'),
                            ha='center', va='center')

        for idx, link in enumerate(self.links):
            # here we connect the levels with the links
            # x1, x2   y1, y2
            for i in link:
                # i is a tuple: (end_level_id,ls,linewidth,color)
                start = self.positions[idx]*(self.dimension+self.space)
                x1 = start + self.dimension
                x2 = self.positions[i[0]]*(self.dimension+self.space)
                y1 = self.energies[idx]
                y2 = self.energies[i[0]]
                line = Line2D([x1, x2], [y1, y2],
                              ls=i[1],
                              linewidth=i[2],
                              color=i[3])
                ax.add_line(line)

        for box in self.electons_boxes:
            # here we add the boxes
            # x,y,boxes,electrons,side,spacing_f
            x, y, boxes, electrons, side, spacing_f = box
            plot_orbital_boxes(ax, x, y, boxes, electrons, side, spacing_f)

    def __auto_adjust(self):
        '''
        Method of ED class
        This method use the ratio to set the best dimension and space between
        the levels.

        Affects
        -------
        self.dimension
        self.space
        self.offset

        '''
        # Max range between the energy
        Energy_variation = abs(max(self.energies) - min(self.energies))
        if self.dimension == 'auto' or self.space == 'auto':
            # Unique positions of the levels
            unique_positions = float(len(set(self.positions)))
            space_for_level = Energy_variation*self.ratio/unique_positions
            self.dimension = space_for_level*0.7
            self.space = space_for_level*0.3

        if self.offset == 'auto':
            self.offset = Energy_variation*self.offset_ratio




            

        