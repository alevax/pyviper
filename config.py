config = {
    "regulators_species":"human",
    "regulators_filepaths":{
        "human":{
            "tfs": None,
            "cotfs": None,
            "sig": None,
            "surf": None
        },
        "mouse":{
            "tfs": None,
            "cotfs": None,
            "sig": None,
            "surf": None
        }
    }
}

def set_regulators_filepath(group, species, new_filepath):
    """\
    Allows the user to use a custom list of regulatory proteins instead of the
    default ones within pyVIPER's data folder.

    Parameters
    ----------
    group
        A group of regulatory proteins of either: "tfs", "cotfs", "sig" or "surf".
    species
        The species to which the group of proteins belongs to: "human" or "mouse".
    new_filepath
        The new filepath that should be used to retrieve these sets of proteins.

    Returns
    -------
    None
    """
    if not species in ["human", "mouse"]:
        raise ValueError("Unsupported species: " + str(species))
    if not group in ["tfs", "cotfs", "sig", "surf"]:
        raise ValueError("Unsupported species: " + str(group))
    config['regulators_filepaths'][species][group] = new_filepath

def set_regulators_species_to_use(species):
    """\
    Allows the user to specify which species they are currently studying, so the
    correct sets of regulatory proteins will be used during analysis.

    Parameters
    ----------
    species
        The species to which the group of proteins belongs to: "human" or "mouse".

    Returns
    -------
    None
    """
    if not species in ["human", "mouse"]:
        raise ValueError("Unsupported species: " + str(species))
    config['regulators_species'] = species
