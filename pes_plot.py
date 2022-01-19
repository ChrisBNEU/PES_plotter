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
    Energy - enthalpy, for y placement on energy diagram
    level - x position on energy diagram. integer 0, 1, 2, etc. 
    """
    def __init__(
            self,
            species,
            position,
            reverse=False,
            ):

        self.h_eV = (species.thermo.h(temp)/1000**2)/96
        self.name = species.name

        self.energy = h_eV  # 0
        self.position = position # 1
        self.bottom_text = name  # 2
        self.top_text = ''
        self.color = 'k'  
        self.right_text = ''
        self.left_text = ''

class comb_species():
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
            self.h_eV += (species.thermo.h(temp)/1000**2)/96
            name = name + ' + ' + species.name
        
        self.name = name

        self.species = species
        self.energy = h_eV  # 0
        self.position = position # 1
        self.top_text = ''
        self.color = 'k'  
        self.right_text = ''
        self.left_text = ''

class pes_reaction():
    """
    feed in a cantera reaction, get an object out of it that we can use for making a chart

    arguments: 
    reaction - a ct reaction object
    phase_gas - gas phase in mechanism (for looking up species)
    phase_surf - solid phase in mechanism (for looking up species)
    reverse - bool, True if we swap reactants and products, and change Ea to Ea from products

    properties:
    reactants - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    products - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    barrier - float, Ea for reaction as value
    equation - string, cantera reaction equation
    links - list of ids for connecting the diagram
            [reactant id, Ea id, product id]
    positions - int, list of positions for reactants, reactions, and products. 
            [reactant pos, Reaction Ea pos, product pos]
    """
    def __init__(
        self,
        reaction,
        phase_gas,
        phase_surf,
        reverse=False,
        ):

        self.equation = reaction.equation
        self.reactants  = {}
        self.products  = {}
        self.links = [-1,-1,-1]
        self.positions = [-1, -1, -1]


        if reverse:
            reactants_obj = reaction.products
            products_obj = reaction.reactants

            # flip reaction equation, e.g. A <=> B is now B <=> A
            str_orig = reaction.equation
            split_list = str_orig.split("<=>")
            str1 = split_list[1] + " <=> " + split_list[0]
            str1 = str1.strip() 
            self.equation = str1
            print("flipped equation: ", str_orig, self.equation)

        else:
            reactants_obj = reaction.reactants
            products_obj = reaction.products
            self.equation = reaction.equation

        # lookup each reactant, put in dictionary as 
        # {species_name : enthalpy at phase temp (in Ev) * stoich coeff}
        total_reac_enth = 0
        reac_str = ""
        for i in reactants_obj: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            reac_str += f"{i} "
            total_reac_enth += reactants_obj[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96

        reac_str = reac_str.strip()
        self.reactants[reac_str] = total_reac_enth 
        self.reactants[reac_str] = round(self.reactants[reac_str],3)

        total_prod_enth = 0
        prod_str = ""
        for i in products_obj: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            prod_str += f"{i} "
            total_prod_enth += products_obj[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
        
        prod_str = prod_str.strip()
        self.products[prod_str] = total_prod_enth 
        self.products[prod_str] = round(self.products[prod_str],3)

        # reaction activation energy. need to add to 
        # reactant enthalpy to get barrier
        # if reversed, need to get barrier from the products
        if reverse: 
            self.barrier = (reaction.rate.activation_energy/1000**2)/96 + total_prod_enth
        else: 
            self.barrier = (reaction.rate.activation_energy/1000**2)/96 + total_reac_enth

        self.barrier = round(self.barrier, 3)

        self.energy = self.barrier
        self.name = self.equation
        self.position = position
        self.bottom_text =
        self.top_text = ''
        self.color = 'k'  
        self.right_text = ''
        self.left_text = ''

class pes_plot():
    """
    Plots a potential energy surface 
    (enthalpy vs rxn coordinate) for a 
    given cti file mechanism
    """
    def __init__(
        self,
        yaml_file,
        temp=528,
        press=75,
        ):
        """
        initialize model
        yaml_file = cti or yaml file for mechanism
        temp = temperature (K)
        press = pressure (atm)
        """

        # set initial temps & pressures
        self.temp = temp # kelvin
        self.pressure =  press * ct.one_atm  # Pascals

        # create thermo phases
        self.yaml_file = yaml_file
        self.gas = ct.Solution(yaml_file, "gas")
        self.surf = ct.Interface(yaml_file,"surface1", [self.gas])

        # initialize T and P for each phase
        self.gas.TP = self.temp, self.pressure
        self.surf.TP = self.temp, self.pressure

        # create reaction diagram object
        self.diagram = ED()

        # initialize the reaction object dictionary
        self.pes_rxn_dict = {}

    def find_reactions(self, species, temp):
        """
        find all reactions that involve a certain species or set of species.
        species is a species object
        pes_rxn_dict is a dictionary, reaction equation is the key, PES reaction object is the value
        """
        pes_rxn_dict = {}
        species_names = [i.name for i in species]
        print(species_names)
        # get combined H for species as the starting point for Ea
        
        for index,rxn in enumerate(self.gas.reactions()):
            if all(x in rxn.reactants.keys() for x in species_names):
                pes_obj = pes_reaction(rxn, self.gas, self.surf)
                pes_rxn_dict[pes_obj.equation] = pes_obj

            # if we want to show the reverse reaction, specify that 
            # in call to pes_reaction
            if all(x in j.products.keys() for x in species_names):
                pes_obj = pes_reaction(j, self.gas, self.surf, reverse=True)
                pes_rxn_dict[pes_obj.equation] = pes_obj
                
        for i,j in enumerate(self.surf.reactions()):
            if all(x in j.reactants.keys() for x in species_names):
                pes_obj = pes_reaction(j, self.gas, self.surf)
                pes_rxn_dict[pes_obj.equation] = pes_obj

            # if we want to show the reverse reaction, specify that 
            # in call to pes_reaction
            if all(x in j.products.keys() for x in species_names):
                pes_obj = pes_reaction(j, self.gas, self.surf, reverse=True)
                pes_rxn_dict[pes_obj.equation] = pes_obj
        
        # if no reactions are found, throw an error
        if len(pes_rxn_dict)==0:
            raise Exception(f"no reactions found with reactants {species}")
                
        return pes_rxn_dict

    def update_diag_data():
        


        self.diagram.data = list(zip(
            self.energies,  # 0
            self.positions,  # 1
            self.bottom_texts,  # 2
            self.top_texts,  # 3
            self.colors,  # 4
            self.right_texts,  # 5
            self.left_texts,))  # 6
    

    def update_diagram_links():
        


    def plot_pes_diagram(
        self, 
        species,
        products=None, 
        width=20, 
        height=40, 
        offset=None,
        dimension=None,
        space=None,
        combined=True,
        ):
        """
        plots a potential energy surface given an input for species.
        the "species" are the starting point for the mechanism. each successive 
        run will take an input of the species and get all reactions for that pair. 

        inputs:
        species - (str or [strs]) matching starting reactant species name in cantera mechanism.
        products - (str or [strs]) if specified, required products for the output reactions
        width - (float) matplotlib plot width in inches
        height - (float) matplotlib plot height in inches
        offset - (float) vertical distance that energy level and upper/lower labels are spaced
        dimension - (float) width of platform used for energy level 
        space - (float) distance between bars for energy levels
        combined - (bool) if true combine all reactants to a single energy level. do the same for products.
        """
        # self.data = list(zip(self.energies,  # 0
        #                 self.positions,  # 1
        #                 self.bottom_texts,  # 2
        #                 self.top_texts,  # 3
        #                 self.colors,  # 4
        #                 self.right_texts,  # 5
        #                 self.left_texts,))  # 6


        # get a list of all reactions containing the two species identified
        species_obj = []
        for i in species:
            if i in self.gas.species_names:
                species_obj.append(self.gas.species(i))
            elif i in self.surf.species_names:
                species_obj.append(self.surf.species(i))
            else:
                print(f'species {i} not found!')

        rxns = self.find_reactions(species_obj, self.temp)
        self.pes_rxn_dict.update(rxns)

        link_num = 0
        for i,j in self.pes_rxn_dict.items():
            # generate a pes plot for each pes_reaction reactant
            for k,l in j.reactants.items():
                reac = k
                H_r = l

                # make a new energy level
                self.diagram.add_level(H_r, k, 0)
                j.positions[0] = 0
                    
            j.links[0] = link_num
            link_num+=1

        for i,j in self.pes_rxn_dict.items():
            # plot rxn Ea. for it to show up between species, should be here
            rxn_eq = j.equation
            rxn_Ea = j.barrier

            # make a new energy level
            self.diagram.add_level(rxn_Ea, rxn_eq, 1)
            j.positions[1] = 1
            
            j.links[1] = link_num
            link_num+=1

        for i,j in self.pes_rxn_dict.items():
            # generate a pes plot for each pes_reaction product

            for k,l in j.products.items():
                prod = k
                H_p = l

                # make a new energy level
                self.diagram.add_level(H_p, prod,2)
                j.positions[2] = 2

            # add link id for line drawing
            j.links[2] = link_num
            link_num+=1

        
        for i in self.pes_rxn_dict.values():
            # get connections between each reac - Ea and each Ea - product
            self.diagram.add_link(i.links[0],i.links[1])
            self.diagram.add_link(i.links[1],i.links[2])

        # optional arguments 
        if space: 
            self.diagram.space = space
        if offset:
            self.diagram.offset = offset
        if dimension:
            self.diagram.dimension = dimension

        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=width, height=height)

    def add_next_reaction(
        self, 
        species, 
        width, 
        height, 
        offset=None,
        dimension=None,
        space=None,
        ):
        """
        adds the next reaction to the plot. 
        
        species is the product species selected for the next step
        """

        # get a list of all reactions containing the species identified
        species_obj = []
        for i in species:
            if i in self.gas.species_names:
                species_obj.append(self.gas.species(i))
            elif i in self.surf.species_names:
                species_obj.append(self.surf.species(i))

        rxns = self.find_reactions(species_obj, self.temp)
        self.pes_rxn_dict.update(rxns)

        # ED.position can be assigned to an integer (1,2,3,4, etc) 
        # so we do not need to use "l"
        # get starting position (whatever the last energy level position was)
        starting_pos = max(self.diagram.positions)

        link_num = len(self.diagram.data)

        for i,j in rxns.items():
            # generate a pes plot for each pes_reaction
            for k,l in self.pes_rxn_dict[i].reactants.items():
                reac = k
                H_r = l

                # make a new energy level
                self.diagram.add_level(H_r, k, starting_pos)
                self.pes_rxn_dict[i].positions[0] = starting_pos
                    
            self.pes_rxn_dict[i].links[0] = link_num
            link_num+=1

        for i,j in rxns.items():   
            # plot rxn Ea. for it to show up between species, should be here
            rxn_eq = self.pes_rxn_dict[i].equation
            rxn_Ea = self.pes_rxn_dict[i].barrier

            # make a new energy level
            self.diagram.add_level(rxn_Ea, rxn_eq, starting_pos + 1)
            self.pes_rxn_dict[i].positions[1] = starting_pos + 1
            
            self.pes_rxn_dict[i].links[1] = link_num
            link_num+=1

        for i,j in rxns.items():  
            # generate a pes plot for each pes_reaction 
            for k,l in self.pes_rxn_dict[i].products.items():
                prod = k
                H_p = l

                # make a new energy level
                self.diagram.add_level(H_p, prod,starting_pos + 2)
                self.pes_rxn_dict[i].positions[2] = starting_pos + 2

            # add link id for line drawing
            self.pes_rxn_dict[i].links[2] = link_num
            link_num+=1

        
        for i,j in rxns.items():
            # get connections between each reac - Ea and each Ea-product
            self.diagram.add_link(self.pes_rxn_dict[i].links[0],self.pes_rxn_dict[i].links[1])
            self.diagram.add_link(self.pes_rxn_dict[i].links[1],self.pes_rxn_dict[i].links[2])

        # optional arguements 
        if space: 
            self.diagram.space = space
        if offset:
            self.diagram.offset = offset
        if dimension:
            self.diagram.dimension = dimension
        
        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=width, height=height)



    def _redraw(
        self,
        width=20, 
        height=40,
        offset=None,
        dimension=None,
        space=None,
        ):

        """ redraw after a trim"""
        
        # create new reaction diagram object 
        # (is there a more efficient way to erase the old one?)
        self.diagram = ED()

        link_num = 0
        for i,j in self.pes_rxn_dict.items():

            # generate a pes plot for each pes_reaction reactant
            for k,l in j.reactants.items():
                reac = k
                H_r = l

                # make a new energy level
                self.diagram.add_level(H_r, k, j.positions[0])
                
            rxn_eq = j.equation
            rxn_Ea = j.barrier

            # make a new energy level
            self.diagram.add_level(rxn_Ea, rxn_eq, j.positions[1])

            # generate a pes plot for each pes_reaction product
            for k,l in j.products.items():
                prod = k
                H_p = l

                # make a new energy level
                self.diagram.add_level(H_p, prod, j.positions[2])

        self.diagram.create_data()

        # go through reaction dictionary and match reactant and product entries
        # match only adjacent ones
        for i,j in self.pes_rxn_dict.items():
            for m,n in enumerate(self.diagram.data):

                # there has to be a better way to do this ".keys())[0]"" nonsense
                if str(list(j.reactants.keys())[0]) == n[2] and j.positions[0] == n[1]:
                    j.links[0] = m
                    position = n[1]
                
                if j.equation == n[2] and j.positions[1] == n[1]:
                    j.links[1] = m

                if str(list(j.products.keys())[0]) == n[2] and j.positions[2] == n[1]:
                    j.links[2] = m


        # need to figure out how to add links
        for i in self.pes_rxn_dict.values():
            # get connections between each reac - Ea and each Ea - product
            self.diagram.add_link(i.links[0],i.links[1])
            self.diagram.add_link(i.links[1],i.links[2])

        # optional arguements 
        if space: 
            self.diagram.space = space
        if offset:
            self.diagram.offset = offset
        if dimension:
            self.diagram.dimension = dimension

        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=width, height=height)

    def trim(
        self, 
        reac, 
        width=20, 
        height=40,
        offset=None,
        dimension=None,
        space=None,
        ):
        """
        trims the specified reactions and their species from plot. 

        updates the pes_rxn_object to remove reactions we do not want. 
        runs through diagram.data and pes_rxn_object to update all of the links
        
        reac is the reaction to be trimmed (will remove reactants + products)
        """
        
        for key in list(self.pes_rxn_dict.keys()):
            if reac == key: 
                del self.pes_rxn_dict[key]
        
        self._redraw()


            

        