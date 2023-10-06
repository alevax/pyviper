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
    if not species in ["human", "mouse"]:
        raise ValueError("Unsupported species: " + str(species))
    if not group in ["tfs", "cotfs", "sig", "surf"]:
        raise ValueError("Unsupported species: " + str(group))
    config['regulators_filepaths'][species][group] = new_filepath

def set_regulators_species_to_use(species):
    if not species in ["human", "mouse"]:
        raise ValueError("Unsupported species: " + str(species))
    config['regulators_species'] = species
