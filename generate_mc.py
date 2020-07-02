import mbuild as mb
import numpy as np
import copy

# Set the defaults for all the required arguments
# It doesn't seem like defaults_dict covers ALL of the req arguments
# The default values are pass into and pulled from argparse entries
'''
defaults_dict = {'stoichiometry': {'Mo': 1, 'V': 0.15, 'Nb': 0.13, 'Te': 0.12},
                 'dimensions': [1, 1, 1], # x, y, z repeats
                 'template': 'M1UnitCell.pdb',  # Path to file to load into mbuild. 
                 'reactant_composition': {'C2H6': 1},
                 'crystal_separation': 25.0,
                 'z_reactor_size': 20.0,
                 'reactant_num_mol': None,
                 'reactant_density': None,
                 'forcefield': None,
                 'integrate_crystal': False}
'''
class UnitCell(mb.Compound):
    # This class will contain the unit cell for manipulation and replication
    def __init__(self, template, stoichiometry_dict):
        # Call the mb.Compound initialization
        super().__init__()
        # Load the unit cell
        mb.load(template, compound=self)
        # Replacable atoms in the matrix are assigned as type `X'
        atom_types, atom_probs, atom_mass_dict = calculate_probabilities(stoichiometry_dict)
        for particle in self.particles():  # Replace all X atoms with Mo, V, Nb, Te atoms.
            if particle.name == 'X':
                # `Randomly' select atom type based on stoichiometry_dict
                particle.name = np.random.choice(atom_types, p=atom_probs)

class Surface(mb.Compound):
    # This class will describe the surface and consist of several
    # UnitCell instances in a specified dimension
    # Default stoichiometry found in: Nanostructured Catalysts: Selective
    # Oxidations (Hess and Schl\"ogl, 2011, RSC)
    def __init__(self, surface_dimensions, template, stoichiometry_dict,
                 crystal_bonds, x_extent, y_extent, z_extent):
        # Call the mb.Compound initialization
        super().__init__()
        # OUTER LOOP: Create multiple layers based on the input dimensions
        for z_repeat in range(surface_dimensions[2]):
            # MIDDLE LOOP: Multiply up each x_row to create as many y repeats as specified
            # complete_cell_matrix is required to keep track of new bonds across diagonal elements
            complete_cell_matrix = []
            previous_row = None
            for y_repeat in range(surface_dimensions[1]):
                current_row = [] # List of UnitCell instances
                # INNER LOOP: First, create as many x repeats as specified
                # Note: Each cell has 159 atoms in it
                previous_cell = None
                for x_repeat in range(surface_dimensions[0]):
                    print("\rAdding " + repr([x_repeat, y_repeat, z_repeat])
                          + " to system...", end=" ")
                    current_cell = UnitCell(template, # Creates new UnitCell instance
                                            stoichiometry_dict)
                    current_row.append(current_cell)
                    current_cell.translate([x_repeat * x_extent,
                                            y_repeat * y_extent,
                                            z_repeat * z_extent])
                    self.add(current_cell)
                    if crystal_bonds and (previous_cell is not None):
                        self.add_x_connecting_bonds(previous_cell,
                                                    current_cell)
                    previous_cell = current_cell
                if crystal_bonds and (previous_row is not None):
                    for cell_ID, current_cell_y in enumerate(current_row):
                        previous_cell_y = previous_row[cell_ID]
                        self.add_y_connecting_bonds(previous_cell_y,
                                                    current_cell_y)
                complete_cell_matrix.append(current_row)
                previous_row = current_row
            # Now that the cell_matrix is complete for this layer, there might
            # be a few more bonds to add in.
            # Go across all rows first
            for y_coord in range(surface_dimensions[1]):
                # Now each column
                for x_coord in range(surface_dimensions[0]):
                    # Bonds located across the diagonals (i.e. [0, 0] bonded
                    # to [1, 1]; [0, 1] bonded to [1, 2] etc.)
                    if crystal_bonds and (x_coord + 1 < surface_dimensions[0])\
                       and (y_coord + 1 < surface_dimensions[1]):
                        first_cell = complete_cell_matrix[x_coord][y_coord]
                        second_cell = complete_cell_matrix[
                            x_coord + 1][y_coord + 1]
                        self.add_diagonal_connecting_bonds(first_cell,
                                                           second_cell)
        print()

    def add_x_connecting_bonds(self, cell1, cell2):
        self.add_bond([cell1[60], cell2[21]])
        self.add_bond([cell1[137], cell2[13]])
        self.add_bond([cell1[134], cell2[16]])
        self.add_bond([cell1[122], cell2[16]])
        self.add_bond([cell1[19], cell2[119]])
        self.add_bond([cell1[66], cell2[33]])
        self.add_bond([cell1[18], cell2[65]])

    def add_y_connecting_bonds(self, cell1, cell2):
        self.add_bond([cell1[72], cell2[27]])
        self.add_bond([cell1[1], cell2[58]])
        self.add_bond([cell1[1], cell2[73]])
        self.add_bond([cell1[4], cell2[123]])
        self.add_bond([cell1[4], cell2[141]])
        self.add_bond([cell1[6], cell2[141]])
        self.add_bond([cell1[6], cell2[156]])
        self.add_bond([cell1[114], cell2[12]])
        self.add_bond([cell1[159], cell2[12]])

    def add_diagonal_connecting_bonds(self, cell1, cell2):
        self.add_bond([cell1[61], cell2[21]])


#########################

def crystal_system(bottom_crystal, top_crystal=None, crystal_separation=25.0):
    system = mb.Compound()
    if not top_crystal:
        top_crystal = mb.clone(bottom_crystal)
    top_COM = copy.deepcopy(top_crystal.pos)
    bottom_COM = copy.deepcopy(bottom_crystal.pos)
    top_crystal.translate(-top_crystal.pos)
    bottom_crystal.translate(-bottom_crystal.pos)
    top_crystal.rotate(np.pi, [1, 0, 0])
    bottom_crystal.translate([0, 0, crystal_separation / 20.0])
    top_crystal.translate([0, 0, -crystal_separation / 20.0])
    system.add(bottom_crystal)
    system.add(top_crystal)
    return system

def morphology(args):
    '''
    Creates a surface --> Cover it --> Unblock active sites
    Make a clone of the surface.
    Pass both surfaces into crystal_system() function
    Save a file with the top and bottom surfaces
    '''
    output_file = create_output_file_name(args)
    surface1 = Surface(args.dimensions, args.template,
                               args.stoichiometry, args.crystal_bonds,
                               args.crystal_x, args.crystal_y, args.crystal_z)
    surface2 = mb.clone(surface1)
    system = crystal_system(surface1, surface2, args.crystal_separation)
    system.save(output_file, overwrite=True, box=None)
    return system


def simulation_system(morphology, args): # Corresponds with the create_morphology() function in rhaco
    # Get the crystal IDs because we're going to need them later so that HOOMD
    # knows not to integrate them.
    crystal_IDs = range(system.n_particles)
    # Now we can populate the box with reactant
    reactant_components, reactant_probs, reactant_masses = calculate_probabilities(
        args.reactant_composition, ratio_type='number')
    
    # Define the regions that the hydrocarbons can go in, so we don't end
    # up with them between layers
    box_top = mb.Box(mins=[-(args.crystal_x * args.dimensions[0]) / 2.0,
                           -(args.crystal_y * args.dimensions[1]) / 2.0,
                           args.crystal_separation / 20.0
                           + (args.crystal_z * args.dimensions[2])],
                     maxs=[(args.crystal_x * args.dimensions[0]) / 2.0,
                           (args.crystal_y * args.dimensions[1]) / 2.0,
                           args.z_reactor_size / 2.0])
    box_bottom = mb.Box(mins=[-(args.crystal_x * args.dimensions[0]) / 2.0,
                              -(args.crystal_y * args.dimensions[1]) / 2.0,
                              -args.z_reactor_size / 2.0],
                        maxs=[(args.crystal_x * args.dimensions[0]) / 2.0,
                              (args.crystal_y * args.dimensions[1]) / 2.0,
                              -args.crystal_separation / 20.0
                              - (args.crystal_z * args.dimensions[2])])
    box_top_vol = np.prod(box_top.maxs - box_top.mins)
    box_bottom_vol = np.prod(box_bottom.maxs - box_bottom.mins)
    reactor_vol = box_top_vol + box_bottom_vol

    # No reactants given
    if (args.reactant_density is None):
        number_of_reactant_mols = 0
    elif (args.reactant_density is not None):
        # Work backwards to come up with how many reactant molecules are needed
        # to get the specified density.
        # Get the average mass for each molecule based on reactant probabilities
        mass_per_n = np.sum([reactant_masses[key] * reactant_probs[index] for index,
                             key in enumerate(reactant_components)])
        # Given the reactor volume and the specified reactant density, calculate
        # the total number of reactant molecules needed.
        # to convert from CGS (g/cm^{3} -> AMU/nm^{3})
        reactant_density_conv = args.reactant_density * G_TO_AMU / (CM_TO_NM**3)
        number_of_reactant_mols = int(reactant_density_conv * reactor_vol / mass_per_n)

        reactant_compounds = []
        n_compounds = []
        for compound_index, reactant_molecule in enumerate(reactant_components):
            reactant_compounds.append(mbuild_template(reactant_molecule))
            n_compounds.append(int(np.round(np.round(
                reactant_probs[compound_index] * number_of_reactant_mols) / 2.0)))
        reactant_top = mb.packing.fill_box(reactant_compounds, n_compounds, box_top,
                                      seed=np.random.randint(0, 2**31 - 1)) 
        reactant_bottom = mb.packing.fill_box(reactant_compounds, n_compounds,
                                         box_bottom,
                                         seed=np.random.randint(0, 2**31 - 1)) 
        system.add(reactant_top)
        system.add(reactant_bottom)


    if "M1UnitCell.pdb" in args.template:
        # Check the separation of crystal and reactant that we will use later is
        # correct. Get the set of atom types that produce the crystal (and
        # don't include the base atom type, which we asusme to be oxygen).
        names = [particle.name for particle_ID, particle in
                 enumerate(system.particles()) if (particle_ID in crystal_IDs)
                 and (particle.name != 'O')]
        # Ensure that this is the same as the stoichiometry dictionary keys
        assert(np.array_equal(args.stoichiometry.keys(), set(names)))

    # Generate the morphology box based on the input parameters
    system_box = mb.Box(mins=[-(args.crystal_x * args.dimensions[0]) / 2.0,
                              -(args.crystal_y * args.dimensions[1]) / 2.0,
                              -args.z_reactor_size / 2.0],
                        maxs=[(args.crystal_x * args.dimensions[0]) / 2.0,
                              (args.crystal_y * args.dimensions[1]) / 2.0,
                              args.z_reactor_size / 2.0])
    print("Morphology generated.")
    

    if (args.forcefield is None) or ((len(args.forcefield[0]) == 0) and (len(args.forcefield[1]) == 0)):
        print("Saving morphology...")
        system.save(output_file, overwrite=True, box=system_box)
        # Fix the images because mbuild doesn't set them correctly
        morphology = fix_images(output_file)
    else:
        print("Applying forcefield...")





### UTILITY / HELPER FUNCTIONS:

def get_masses(reactant_names):
    number_ratio = {}
    mass_dict = {}
    # First split out the key names into atoms
    for reactant_name in reactant_names:
        # Consult the mass lookup table
        total_mass = mbuild_template(reactant_name).mass
        mass_dict[reactant_name] = total_mass
    # Return dictionary of number ratios
    return mass_dict


def zero_out_images(morphology):
    morphology['image_text'] = [['0', '0', '0']]\
        * len(morphology['position_text'])
    morphology['image_attrib'] = {'num': morphology['position_attrib']['num']}
    return morphology


def get_bond_dict(morphology):
    bond_dict = {atom_id: [] for atom_id, atom_type in
                 enumerate(morphology['type_text'])}
    for bond in morphology['bond_text']:
        bond_dict[int(bond[1])].append(int(bond[2]))
        bond_dict[int(bond[2])].append(int(bond[1]))
    return bond_dict


def load_morphology_xml(xml_file_name):
    morphology_dictionary = OrderedDict()
    with open(xml_file_name, 'r') as xml_file:
        xml_tree = ET.parse(xml_file)
    root = xml_tree.getroot()
    morphology_dictionary['root_tag'] = root.tag
    morphology_dictionary['root_attrib'] = root.attrib
    morphology_dictionary['root_text'] = root.text
    for config in root:
        morphology_dictionary['config_tag'] = config.tag
        morphology_dictionary['config_attrib'] = config.attrib
        morphology_dictionary['config_text'] = config.text
        for child in config:
            if len(child.attrib) > 0:
                morphology_dictionary[child.tag + '_attrib'] = {
                    key.lower(): val for key, val in child.attrib.items()}
            else:
                morphology_dictionary[child.tag + '_attrib'] = {}
            if child.text is not None:
                morphology_dictionary[child.tag + '_text'] = [
                    x.split() for x in child.text.split('\n') if len(x) > 0]
            else:
                morphology_dictionary[child.tag + '_text'] = []
    return morphology_dictionary



def write_morphology_xml(morphology_dictionary, output_file_name):
    # morphology_dictionary is a bunch of keys with the tagnames given for
    # both attributes and text: tag + '_attrib', tag + '_text'
    # The only boilerplate bits are the 'root_tag', 'root_attrib', and
    # 'root_text', which is (obviously) the outmost layer of the xml.
    # Immediately inside are the 'config_tag', 'config_attrib', and
    # 'config_text'. Everything else is a child of config.
    morphology_dictionary = check_wrapped_positions(morphology_dictionary)
    # Build the xml tree.
    root = ET.Element(morphology_dictionary['root_tag'],
                      **morphology_dictionary['root_attrib'])
    root.text = morphology_dictionary['root_text']
    config = ET.Element(morphology_dictionary['config_tag'],
                        **morphology_dictionary['config_attrib'])
    config.text = morphology_dictionary['config_text']
    # Find the remaining elements to make (set is easier here, but a
    # disordered structure, so instead we use lists to keep the order
    # consistent with reading in).
    all_child_tags = ['_'.join(key.split('_')[:-1]) for key in
                      morphology_dictionary.keys()
                      if '_'.join(key.split('_')[:-1]) not in
                      ['root', 'config']]
    child_tags = []
    for tag in all_child_tags:
        # The list comprehension makes two blank entries for some reason and
        # I can't work out why. This will just skip those two entries, as well
        # as make the set.
        if (tag not in child_tags) and (len(tag) > 0):
            child_tags.append(tag)
    for child_tag in child_tags:
        child = ET.Element(child_tag,
                           **morphology_dictionary[child_tag + '_attrib'])
        if child_tag != "external_forcefields":
            data_to_write = '\n'.join(['\t'.join(el) for el in
                                       morphology_dictionary[
                                           child_tag + '_text']])
        else:
            data_to_write = '\n'.join([el for el in morphology_dictionary[child_tag + "_text"]])
        if len(data_to_write) > 0:
            child.text = '\n' + data_to_write + '\n'
        child.tail = '\n'
        config.append(child)
    root.insert(0, config)
    tree = ET.ElementTree(root)
    tree.write(output_file_name, xml_declaration=True, encoding='UTF-8')
    print("XML file written to", str(output_file_name) + "!")


def calculate_probabilities(input_dictionary, ratio_type='stoic'):
    '''
    This function takes an input dictionary corresponding to the relative
    ratios of some parameter, then returns normalized probabilities for each
    option that can be used to make choices with appropriate bias.
    '''
    choices = list(input_dictionary.keys())
    number_ratios = np.array(list(input_dictionary.values()))
    mass_dict = None
    if ratio_type == 'number':
        mass_dict = get_masses(input_dictionary.keys())
    probabilities = list(number_ratios / np.sum(number_ratios))
    return choices, probabilities, mass_dict


def create_output_file_name(args, file_type='hoomdxml'):
    if args.signac is True:
        return 'output.hoomdxml'
    else:
        output_file = "out"
        for (arg_name, arg_val) in sorted(args._get_kwargs()):
            try:
                if (arg_val == defaults_dict[arg_name]):
                    continue
            except KeyError:
                continue
            output_file += "-"
            if arg_name == "stoichiometry":
                output_file += "S_"
                for key, val in arg_val.items():
                    output_file += str(key) + ":" + str(val) + "_"
                output_file = output_file[:-1]
            elif arg_name == "reactant_composition":
                output_file += "RC_"
                for key, val in arg_val.items():
                    output_file += str(key) + ":" + str(val) + "_"
                output_file = output_file[:-1]
            elif arg_name == "dimensions":
                output_file += "D_" + "x".join(list(map(str, arg_val)))
            elif arg_name == "template":
                output_file += "T_" + args.template.split("/")[-1].split(
                    ".")[0]
            elif arg_name == "forcefield":
                if len(args.forcefield[0]) > 0:
                    output_file += "F1"
                    for FF in args.forcefield[0]:
                        output_file += "_" + os.path.split(FF)[1]
                if len(args.forcefield[1]) > 0:
                    if len(args.forcefield[0]) > 0:
                        output_file += "-"
                    output_file += "F2"
                    for FF in args.forcefield[1]:
                        output_file += "_" + os.path.split(FF)[1]
            elif arg_val is False:
                output_file += arg_name[0].upper() + "_Off"
            elif arg_val is True:
                output_file += arg_name[0].upper() + "_On"
            else:
                output_file += arg_name[0].upper() + "_" + str(arg_val)
        return output_file + '.' + file_type


