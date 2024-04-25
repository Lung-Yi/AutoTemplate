import os
import random

seed = 42
all_rxn_class = ["AdamsDecarboxylation",
                "Baylis-HillmanReaction",
                "Buchwald-HartwigCross-Coupling",
                "Chan_LamCoupling",
                "DielsAlder",
                "FischerIndoleSynthesis",
                "Friedel-CraftsAcylation",
                "Friedel-CraftsAlkylation",
                "GrignardReaction",
                "HiyamaCoupling",
                "HuisgenCycloaddition",
                "Hydrogenation",
                "Kabachnik-FieldsReaction",
                "KumadaCoupling",
                "MannichReaction",
                "NegishiCoupling",
                "PausonKhandReaction",
                "ReductiveAmination",
                "SuzukiCoupling",
                "WittigReaction",
                ]

data_processed_path = "../data"
fractions = [10, 25, 50, 75]

for rxn in all_rxn_class:
    input_rxn_file = os.path.join(os.path.join(data_processed_path, rxn), f"MappingResult_{rxn}.txt")
    with open(input_rxn_file, "r") as f:
        input_data = f.readlines()
    for fraction in fractions:
        data_output_path = f"../data_{fraction}perc/{rxn}"
        random.seed(seed)
        sample_size = len(input_data) * fraction // 100
        output_data = [input_data[i] for i in sorted(random.sample(range(len(input_data)), sample_size))]
        os.makedirs(data_output_path, exist_ok=True)
        with open(os.path.join(data_output_path, f"MappingResult_{rxn}.txt"), "w") as g:
            g.writelines(output_data)