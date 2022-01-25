import cantera as ct
import os
import sys
sys.path.append(f'{os.getcwd()}/PyEnergyDiagram')
from energydiagram import ED
import logging
from matplotlib.lines import Line2D
import itertools
import copy

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
    takes a single "species" cantera object and makes an object for 
    diagramming that species on a potential energy plot.

    species - a string with the species name as it is written in the Cantera 
            file. 
    position - an integer determining what the x position the species will
            occupy on the pes plot
    temp - the temperature of the system being studied. 
    """
    def __init__(
            self,
            species,
            position,
            temp,
            ):

        self.h_eV = (species.thermo.h(temp)/1000**2)/96
        self.energy = round(self.h_eV, 3) 
        self.name = species.name
        self.position = position
        self.index = -1 #index in diagram. used for connections


class comb_species():
    """
    makes a grouping of reactants or products. Can be used with a single species. 
    species - can be either of the following: 
            a list of tuples with the arguments for the reactant name and 
            stoichiometry in the form [(str, int),], or a list of strings 
            a list of strings that are simply the species names.

    """
    def __init__(
            self,
            species,
            position,
            temp,
            ):

        self.h_eV = 0
        self.name_list = []
        self.stoich_name_list = []
        self.reactants = {}
        if isinstance(species[0], ct.Species):
            for spec in species: 
                self.h_eV += (spec.thermo.h(temp)/1000**2)/96
                self.name_list.append(spec.name)
                self.reactants.update({spec.name:1})
                self.stoich_name_list.append(spec.name)

        elif isinstance(species[0], tuple):
            for spec,stoich in species: 
                self.h_eV += stoich*(spec.thermo.h(temp)/1000**2)/96
                self.name_list.append(spec.name)
                if stoich > 1:
                    self.stoich_name_list.append(str(stoich) + spec.name)
                    self.reactants.update({spec.name:stoich})
                else: 
                    self.stoich_name_list.append(spec.name)
                    self.reactants.update({spec.name:1})

        

        self.energy = round(self.h_eV, 3)
        self.name = "+".join(self.stoich_name_list)
        self.species = species
        self.position = position
        self.index = -1 #index in diagram. used for connections

    def link_rxn(self,reaction):
        """
        checks if combined species object is in a given reaction
        and if it is adjacent to it (position = reaction position +/-1)
         0 if it is not appropriate to link, 
         1 if it creates a forward link,
        -1 if it is a reverse link
        """
        link = 0
        # if all the species in the combined reactants are in the reaction
        if self.reactants == reaction.reac_stoich and reaction.position == self.position + 1:
            link = 1
        elif self.reactants == reaction.prod_stoich and reaction.position == self.position - 1:
            link = -1

        return link           

class pes_reaction():
    """
    feed in a cantera reaction, get an object out of it that we can use for 
    making a chart

    arguments: 
    reaction - a ct reaction object
    phase_gas - gas phase in mechanism (for looking up species)
    phase_surf - solid phase in mechanism (for looking up species)
    reverse - bool, True if we swap reactants and products, and change Ea to Ea 
            from products

    properties:
    reactants - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    products - dict, reactant string as key (e.g. "CO2+H2O"), combined Hf as value
    barrier - float, Ea for reaction as value
    equation - string, cantera reaction equation
    """
    def __init__(
        self,
        reaction,
        position,
        phase_gas,
        phase_surf,
        reverse=False,
        ):

        self.orig_equation = reaction.equation
        self.equation = reaction.equation
        self.reactants  = {}  # combined reactants
        self.products  = {}   # combined products

        if reverse:
            self.reac_stoich = reaction.products   # dict, {spec: stoich}
            self.prod_stoich  = reaction.reactants # dict, {spec: stoich}

            # flip reaction equation, e.g. A <=> B is now B <=> A
            str_orig = reaction.equation
            split_list = str_orig.split("<=>")
            str1 = split_list[1] + " <=> " + split_list[0]
            str1 = str1.strip() 
            self.equation = str1
            print("flipped equation: ", str_orig, self.equation)

        else:
            self.reac_stoich = reaction.reactants
            self.prod_stoich  = reaction.products
            self.equation = reaction.equation

        # lookup each reactant, put in dictionary as 
        # {species_name : enthalpy at phase temp (in Ev) * stoich coeff}
        total_reac_enth = 0
        reac_str = "+".join(self.reac_stoich)

        for i in self.reac_stoich: 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            total_reac_enth += self.reac_stoich[i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96

        self.reactants[reac_str] = total_reac_enth 
        self.reactants[reac_str] = round(self.reactants[reac_str],3)

        total_prod_enth = 0
        prod_str = "+".join(self.prod_stoich )

        for i in self.prod_stoich : 
            if i in phase_gas.species_names:
                phase = phase_gas
            elif i in phase_surf.species_names:
                phase = phase_surf
            else: 
                logging.error(f"species {i} cannot be found")

            total_prod_enth += self.prod_stoich [i] * (phase.species(i).thermo.h(phase.T)/1000**2)/96
        
        self.products[prod_str] = total_prod_enth 
        self.products[prod_str] = round(self.products[prod_str],3)

        # reaction activation energy. need to add to 
        # reactant enthalpy to get barrier
        # if reversed, need to get barrier from the products
        if reverse: 
            self.barrier = (reaction.rate.activation_energy/1000**2)/96 + total_prod_enth
        else: 
            self.barrier = (reaction.rate.activation_energy/1000**2)/96 + total_reac_enth

        self.energy = round(self.barrier, 3)
        self.position = position
        self.name = self.equation.strip()
        self.index = -1 #index in diagram. used for connections

class pes_plot():
    """
    Plots a potential energy surface 
    (enthalpy vs rxn coordinate) for a 
    given cti file mechanism
    """
    def __init__(
        self,
        yaml_file,
        temp=300,
        press=1,
        ):
        """
        initialize cantera model for PES plot
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
        self.species_list = []
        self.comb_species_list = []
        self.reaction_list = []

    def purge(self, position):
        """
        purge all components at a given position
        make copies for iterating
        """
        self.species_list = list(filter(lambda num: num.position != position ,self.species_list))
        self.comb_species_list = list(filter(lambda num: num.position != position ,self.comb_species_list))
        self.reaction_list = list(filter(lambda num: num.position != position ,self.reaction_list))

    def update_diag_data(self):
        """
        updates the PyEnergyDiagram object's data. called each time "plot" is called
        """
        self.diagram.reset()
        mechanism_list = self.species_list+self.comb_species_list+self.reaction_list
        for index,item in enumerate(mechanism_list):
            self.diagram.add_level(
                item.energy,
                bottom_text=item.name, 
                position=item.position, 
                color='k',
                top_text='Energy', 
                right_text='', 
                left_text=''
                )
            item.index = index
        
        # draw links as appropriate.
        # if species have a "combined_species" to the left or right that they're in, add link
        # if combined species have a reaction to the left that they are in, add link
        for spec in self.species_list:
            for spec_c in self.comb_species_list:
                if spec.name in spec_c.name_list:
                    if spec_c.position == spec.position + 1:
                        self.diagram.add_link(
                            spec.index, 
                            spec_c.index,
                            color='k',
                            ls='--',
                            linewidth=1,
                            )
                    elif spec_c.position == spec.position - 1:
                        self.diagram.add_link(
                            spec_c.index, 
                            spec.index,
                            color='k',
                            ls='--',
                            linewidth=1,
                            )
        # for combined species objects, automatically draw a link between 
        # reactions where it is either a reactant or a product. 
        for spec_c in self.comb_species_list:
            for rxn in self.reaction_list:
                if spec_c.link_rxn(rxn) == 1:
                    self.diagram.add_link(
                        spec_c.index, 
                        rxn.index,
                        color='k',
                        ls='--',
                        linewidth=1,
                        )
                elif spec_c.link_rxn(rxn) == -1:
                    self.diagram.add_link(
                        rxn.index, 
                        spec_c.index,
                        color='k',
                        ls='--',
                        linewidth=1,
                        )


    
    def plot(self):
        """
        plot the potential energy surface
        """
        self.update_diag_data()
        self.diagram.offset = 0.01
        self.diagram.plot(show_IDs=True, ylabel="Energy / $eV$", width=20, height=20)
        
    
    def add_specie(self, species, position):
        """
        adds a single species to your PES plot.
        
        species - a string for a species 
        position - an integer for the position of the species
        """

        if any((species==spec.name and position==spec.position) for spec in self.species_list):
                print(f'species {species} already in this diagram position')

        elif species in self.gas.species_names:
            new_spec = single_species(
                self.gas.species(species),
                position,
                self.gas.T
                )
            self.species_list.append(new_spec)

        elif species in self.surf.species_names:
            new_spec = single_species(
                self.surf.species(species),
                position,
                self.gas.T
                )
            self.species_list.append(new_spec)
        else: 
            print(f'species {species} could not be found, so it was not added')

    def add_combined_specie(self, species, position):
        """
        adds a single group of species (i.e. reactants or products) to your 
        potential energy plot. will be displayed as 
        species - can be either of the following: 
                    a list of strings for the species to be added (will assume 
                    that the stoichiometric coefficient is = 1)
                    a list of tuples, in the following format with the species_name 
                    and the stoichiometric coefficient, e.g.
                    [("CH4, 1), (H2,2)]
        position - an integer for the position of the species
        """

        if len(species) > 3:
            logging.error("4 reactants and beyond not supported")

        # check for combined species 
        # species1 + species2 
        # and 
        # species2 + species 1
        # it would be good to figure out how to standardize the arguments so 
        # that is can accept (species, stoich) instead of only (species, stoich)
        reac_perms = []  

        # if each entry is just the species name, permutate w/o stoichiometry
        if isinstance(species[0], str):
            perms = itertools.permutations(species)
            for s in list(perms):
                reac_str = "+".join(s)
                reac_perms.append(reac_str)
        
        # add stoichiometry if specified. 
        elif isinstance(species[0], tuple):
            stoich_spec = []
            for spec, stoich in species: 
                if stoich > 1:
                    stoich_spec.append(str(stoich) + spec)
                else:
                    stoich_spec.append(spec)
            if len(species) > 1:
                perms = itertools.permutations(stoich_spec)
                for s in list(perms):
                    reac_str = "+".join(s)
                    reac_perms.append(reac_str)
            else: 
                perms = stoich_spec
                reac_str = perms
                reac_perms = reac_str                  

        if any((spec_c.name in reac_perms and spec_c.position == position) for spec_c in self.comb_species_list):
            print(f'species {reac_perms[0]} already in this diagram position')

        else:
            spec_list = []
            if isinstance(species[0], str):
                for spec in species:
                    if spec in self.gas.species_names:
                        spec_list.append(self.gas.species(spec))
                    elif spec in self.surf.species_names:
                        spec_list.append(self.surf.species(spec))
                    else: 
                        print(f'species {species} could not be found, so it was not added')
            if isinstance(species[0], tuple): 
                for spec, stoich in species:                             
                    if spec in self.gas.species_names:
                        spec_list.append((self.gas.species(spec), stoich))
                    elif spec in self.surf.species_names:
                        spec_list.append((self.surf.species(spec),stoich))
                    else: 
                        print(f'species {species} could not be found, so it was not added')

            if len(spec_list) == len(species):
                new_comb_species = comb_species(
                    spec_list,
                    position, 
                    self.gas.T
                    )
                self.comb_species_list.append(new_comb_species)
            else: 
                print(f'could not make combined species object. ensure all species names are correct')


    def add_reaction(self, reaction, position, reverse=False):
        """
        adds a reaction to the Potential energy surface plot. 
        accepts 
        reaction - a string for the reaction
        position - an integer for the position of the reaction
        reverse - whether or not the reaction should be "flipped" so it 
                displays in reverse.
        """

        if any((reaction==rxn.orig_equation and position==rxn.position) for rxn in self.reaction_list):
                print(f'reaction {reaction} already in this diagram position')

        elif reaction in self.gas.reaction_equations():
            rxn_index = self.gas.reaction_equations().index(reaction)
            new_rxn = pes_reaction(
                self.gas.reaction(rxn_index),
                position,
                self.gas,
                self.surf,
                reverse,
                )
            self.reaction_list.append(new_rxn)

        elif reaction in self.surf.reaction_equations():
            rxn_index = self.surf.reaction_equations().index(reaction)
            new_rxn = pes_reaction(
                self.surf.reaction(rxn_index),
                position,
                self.gas,
                self.surf,
                reverse,
                )
            self.reaction_list.append(new_rxn)
        else: 
            print(f'species {species} could not be found, so it was not added')


    def add_species(self, species, position):
        """
        add multiple species
        """
        for spec in species:
            self.add_specie(spec, position)

    def add_combined_species(self, species, position):
        """
        add multiple combined species objects. accepts list of lists, 
        elements in sublist can be tuples or strings, e.g.
        [
            [(CH4, 2), (H2,2)],
            [(CO2,2), CO,1)],
        ]
        """
        for spec_c in species:
            self.add_combined_specie(spec_c, position)

    def add_reactions(self, reactions, position, reverse=False):

        """
        add multiple reactions, accepts a list of strings.
        """
        for rxn in reactions:
            self.add_reaction(rxn, position, reverse)
