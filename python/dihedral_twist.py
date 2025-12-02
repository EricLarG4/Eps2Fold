import pymol
from pymol import cmd
import statistics


def twist_calc(pdb, chain, residue_1, residue_2, state=None):
    """
    Calculate the helical twist between two consecutive stacked bases using PyMOL.

    Parameters:
    pdb (str): PDB accession code.
    chain (str): Chain identifier.
    residue_1 (int): First residue number.
    residue_2 (int): Second residue number.
    state (int, optional): Specific state (model) to calculate angle for. If None, will use current state.

    Returns:
    float: Helical twist angle in degrees (rounded to 2 decimal places).
    """
    base_atoms = {"purine": ["C1'", "N9"], "pyrimidine": ["C1'", "N1"]}

    def get_base_atoms(residue):
        # Note: state is used in count_atoms
        is_purine = (
            cmd.count_atoms(
                f"{pdb} and chain {chain} and resi {residue} and name N9", state=state
            )
            > 0
        )
        return base_atoms["purine"] if is_purine else base_atoms["pyrimidine"]

    atoms_res1, atoms_res2 = get_base_atoms(residue_1), get_base_atoms(residue_2)

    # Note: order is important to get the correct angle
    ordered_atoms = [
        f"{pdb} and chain {chain} and resi {residue_1} and name {atoms_res1[0]}",
        f"{pdb} and chain {chain} and resi {residue_1} and name {atoms_res1[1]}",
        f"{pdb} and chain {chain} and resi {residue_2} and name {atoms_res2[1]}",
        f"{pdb} and chain {chain} and resi {residue_2} and name {atoms_res2[0]}",
    ]

    # Pass the state to get_dihedral
    angle = round(cmd.get_dihedral(*ordered_atoms, state=state), 2)
    return angle


def get_glycosidic_bond_angle(pdb, chain, residue, state=None):
    """
    Calculate the glycosidic bond angle for a given residue.
    The angle is used to determine whether the base is in a syn (s) or anti (a) conformation.

    Parameters:
    pdb (str): PDB accession code.
    chain (str): Chain identifier.
    residue (int): Residue number.
    state (int, optional): Specific state (model) to calculate angle for. If None, will use current state.

    Returns:
    float: GBA dihedral angle in degrees (rounded to 2 decimal places).
    """
    base_atoms = {"purine": "N9", "pyrimidine": "N1"}

    # Use state parameter with count_atoms
    is_purine = (
        cmd.count_atoms(
            f"{pdb} and chain {chain} and resi {residue} and name N9", state=state
        )
        > 0
    )
    key_atom = base_atoms["purine"] if is_purine else base_atoms["pyrimidine"]

    # Note: order is important to get the correct angle
    ordered_atoms = [
        f"{pdb} and chain {chain} and resi {residue} and name O4'",
        f"{pdb} and chain {chain} and resi {residue} and name C1'",
        f"{pdb} and chain {chain} and resi {residue} and name {key_atom}",
        f"{pdb} and chain {chain} and resi {residue} and name C2",
    ]

    # Pass the state to get_dihedral
    gba_angle = round(cmd.get_dihedral(*ordered_atoms, state=state), 2)
    return gba_angle


import os

def twist_calc_fetch(pdb, chain, residue_1, residue_2):
    """
    Calculate the average helical twist and glycosidic bond angles between a
    single stack of two bases across all models of a given PDB.

    Parameters:
    pdb (str): PDB accession code or file path.
    chain (str): Chain identifier.
    residue_1 (int): First residue number.
    residue_2 (int): Second residue number.

    Returns:
    tuple: (float) The average helical twist angle in degrees with corresponding average gba values in
    degrees (rounded to 2 decimal places).
    """
    cmd.delete("all")

    # Check if the input is a file path
    if os.path.isfile(pdb):
        # Load the local file
        cmd.load(pdb)
        # Extract the base name of the file without extension to use as the object name
        pdb_object_name = os.path.splitext(os.path.basename(pdb))[0]
    else:
        # Fetch the PDB structure
        cmd.fetch(pdb)
        pdb_object_name = pdb

    # Get the number of states (models)
    state_count = cmd.count_states(pdb_object_name)
    print(f"Found {state_count} states (models) in {pdb_object_name}")

    twist_angles = []
    gba1_angles = []
    gba2_angles = []

    # Calculate angles for each state (model)
    for state in range(1, state_count + 1):
        try:
            twist = twist_calc(pdb_object_name, chain, residue_1, residue_2, state=state)
            gba1 = get_glycosidic_bond_angle(pdb_object_name, chain, residue_1, state=state)
            gba2 = get_glycosidic_bond_angle(pdb_object_name, chain, residue_2, state=state)

            twist_angles.append(twist)
            gba1_angles.append(gba1)
            gba2_angles.append(gba2)
            print(f"  State {state}: twist={twist}, gba1={gba1}, gba2={gba2}")
        except Exception as e:
            print(f"  Error calculating angles for state {state} of {pdb_object_name}: {e}")
            continue

    # Calculate averages
    if twist_angles:
        avg_twist = round(statistics.mean(twist_angles), 2)
        avg_gba1 = round(statistics.mean(gba1_angles), 2)
        avg_gba2 = round(statistics.mean(gba2_angles), 2)
        print(
            f"  Calculated averages from {len(twist_angles)} states: twist={avg_twist}, gba1={avg_gba1}, gba2={avg_gba2}"
        )
    else:
        # Fall back to the first state if no valid states were found
        print(
            f"  Warning: No valid states found for {pdb_object_name} {chain} {residue_1}-{residue_2}. Trying again with state=1 only."
        )
        try:
            avg_twist = twist_calc(pdb_object_name, chain, residue_1, residue_2, state=1)
            avg_gba1 = get_glycosidic_bond_angle(pdb_object_name, chain, residue_1, state=1)
            avg_gba2 = get_glycosidic_bond_angle(pdb_object_name, chain, residue_2, state=1)
        except Exception as e:
            print(f"  Error with fallback calculation: {e}")
            avg_twist = 0.0
            avg_gba1 = 0.0
            avg_gba2 = 0.0

    return avg_twist, avg_gba1, avg_gba2



def twist_calc_iterate(pdb, chain, stack_pairs):
    """
    Compute the average helical twist angles and glycosidic bond orientation for
    all given stacked bases across all models for a given PDB.

    Parameters:
    pdb (str): PDB accession code.
    chain (str): Chain identifier.
    stack_pairs (list of tuples): List of (residue_1, residue_2) pairs to compute angles for.

    Returns:
    list: List of dictionaries containing the PDB name, stack pair, average twist angle, and glycosidic bond annotation.
    """
    results = []
    for residue_1, residue_2 in stack_pairs:
        print(f"Processing {pdb} {chain} {residue_1}-{residue_2}")
        avg_angle, avg_gba1, avg_gba2 = twist_calc_fetch(
            pdb, chain, residue_1, residue_2
        )

        # Determine base conformation from average GBA angles
        gba_1 = "s" if 10 <= avg_gba1 <= 120 else "a"
        gba_2 = "s" if 10 <= avg_gba2 <= 120 else "a"
        gba_combined = f"{gba_1}{gba_2}"

        results.append(
            {
                "pdb": pdb,
                "stack_name": f"{residue_1}-{residue_2}",
                "twist_angle": avg_angle,
                "gba_angle_1": avg_gba1,
                "gba_angle_2": avg_gba2,
                "gba": gba_combined,
            }
        )
    return results
