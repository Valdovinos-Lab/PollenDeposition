Code and data accompanying Valdovinos 2025 manuscript.
Main Matlab function is RunValdovinos2013_800M_PollServ.m, which runs Valdovinos et al (2013) model for
800 networks of varying connnectance, richness, and nestedness levels.
The code outputs 3 tables with key variables for plant and animal species and for the plant-animal interactions.
RunValdovinos2013_800M_PollServ.m calls functions:
1) J_zero_pattern.m,
2) IntegrateValdovinos2013_dataset.m, which calls create_metadata.m, uniform_rand.m, Valdovinos2013_rhs.m, unpack.m, etc. 
3) calValMechs_matrices,
4) effective_degree.m

Csv fils include:
1) NetProperties_800m.csv: Input for the Matlb code with the properties of the 800 networks.
2) PollinationServicesReceivedByP_muAP3_800m.csv: Output of the Matlab code, with the plant variables.
3) PolServByA_muAP3_800m.csv: Output of the Matlab code, with the pollinator variables.
4) Missing file is EffortAtoP_muAP3_800m, which is an output for every interaction in every network, which GitHub could not handle because it is too large. It can be obtained by running the Matlab code.
