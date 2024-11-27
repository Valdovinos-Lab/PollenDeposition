
% Code accompanying Valdovinos 2025 manuscript.
% Runs Valdovinos et al (2013) model for 800 networks of varying connnectance, richness, 
% and nestedness levels. The code outputs 3 tables with key variables for plant and animal 
% species and for the plant-animal interactions.
% Last Modification 10/25/2024, Davis.
% Calls functions: 1) J_zero_pattern.m, 2) IntegrateValdovinos2013_dataset.m,
% 3) calValMechs_matrices, 4) effective_degree.m
%

muAP=2; % setting up mortality level, which in the ms is case 3, which is low mortality for 
               % both animal and plant species.
sem=0; % seeting up the seed (semilla in Spanish) so random draws are reproducible.

load('800m.mat') % Loads the 800 netoworks, which are sorted by ascending degree for both
                          % plant (rows) and animal (columns) species.
             
NetProperties=readtable('NetProperties_800m.csv'); % loads properties of all 88 networks as 
                                                                             % table.
NetPropertiesData = table2array(NetProperties); % Converts table to array.

nSpecA=zeros(800,1); % Pre-allocating 1 column with 800 rows to enter the number of
                                  % specialist pollinators (animals, with only 1 link) for each network

PolServP_All=[]; % Starting the matrix that will contain all plant variables as empty.
PolServA_All=[]; % Starting the matrix that will contain all animal variables as empty.
EffortAP_All=[]; % Starting the matrix that will contain all foraging effort that pollinators assign
                        % to plants as empty.

for i=1:800
    
   rand('seed',sem+i); % choosing the seed in a predictable way
    
    In=m800{i}; % selecting the i matrix, representing the i network
    [nP, nA]=size(In);% number of plant (nP) and pollinator or animal (nA) species
    J_pattern = J_zero_pattern(In) ;% Calling the function J_zero_pattern, which outputs a
                                                   % matrix for the Jacobian so the integrator knows where
                                                   % the non-zero elements are in this model, to speed up
                                                   % the numerical integration.

    abundP=zeros(nP,2); % preallocating matrix for recording abundance of each plant species
                                    % per network, without (first column) and with (second column) AF.
    abundA=zeros(nA,2);% preallocating matrix for recording abundance of each animal species
                                    % per network, without (first column) and with (second column) AF.
    sPolDepP=zeros(nP,2);% preallocating matrix for recording total per-capita pollen deposition
                                    % of each plant species per network, without (first column) and with
                                    % (second column) AF.
    sPolDepA=zeros(nA,2);% preallocating matrix for recording total per-capita pollen deposition
                                    % performed by each animal species per network, without (first column) 
                                    % and with (second column) AF.
    sVisits_perP=zeros(nP,2);% preallocating matrix for recording total per-capita visits received
                                    % by each plant species per network, without (first column) and with
                                    % (second column) AF.
    sVisits_perPbyA=zeros(nA,2);% preallocating matrix for recording total per-capita visits
                                    % performed by each pollinator species to all its plant species per
                                    % network, without (first column) and with (second column) AF.
    EffDegA=zeros(nA,2);% preallocating matrix for recording effective degree of each pollinator
                                    % species per network, without (first) and with (second column) AF.
                                    % See function "effective_degree" for more details.
    
    degreeP=sum(In,2); % Calculates degree of each plant species
    degreeA=sum(In)'; % Calculates degree of each pollinator species
    nSpecA(i)=sum(degreeA==1); % Calculates number of specialist pollinator species (1 link)

    [fP, fA]=find(In); % Finds each pair of plant and animal species that interacts.
    degAnP=[degreeA(fA) degreeP(fP)];% degree of each interacting animal and plant species pair
    L=sum(degreeP); % number of links in the network
    EffortAP=zeros(L, 2); % Preallocates the effort assign by each pollinator species to each
                                    % plant species in its diet. In the model, it is the alpha variable,
                                    % which is one per interaction or link.

    for frG=[0,1] % Indicates without frG=0, or with frG=1 adaptive foraging (AF)
        
        % Vector indicating which pollinators have AF, in this case all eaither do no or do have AF.
        vectG=frG*ones(1,nA);

        % Numerical integration of Valdovinos et al (2013) model:
        [p, N, a, Alpha, network_metadata]=IntegrateValdovinos2013_dataset(vectG, In, muAP, J_pattern);
       % Model outputs are:
        % p = vector with density of each plant species in the network.
        % N = vector with the density of floral rewards for each plant species.
        % a = vector with density of each animal species in the network.
        % Alpha = matrix with foraging effort that each animal species assign to each plant species.
        % network_metadata = all parameter values stored in a Matlab structure.
        
        % Using model outputs to calculate key reponse variables of this study:
        [M_Visits_perP, M_PolDep]= calValMechs_matrices(Alpha, p, a, network_metadata);
        
        abundP(:, frG+1)=p;
        abundA(:, frG+1)=a;
        
        sPolDepA(:, frG+1)=sum(M_PolDep)';   % Total per-plant (i.e., per-capita) pollen deposition rate.
        sPolDepP(:, frG+1)=sum(M_PolDep,2); % Total per-plant pollen deposition rate by pollinators.
        
        sVisits_perP(:, frG+1)=sum(M_Visits_perP,2); % Total per-plant visits each plant species receives.
        sVisits_perPbyA(:, frG+1)=sum(M_Visits_perP)'; % Total per-plant visits by each pollinator species.

        EffDegA(:, frG+1)=effective_degree(Alpha)'; % Interesting metric but finally not used in this study.

        EffortAP(:, frG+1)=Alpha(In==1);                  
    
    end
    
    % NetPropertiesData_i=[NetPropertiesData(i,:) nSpecA(i)];
    NetPropertiesData_i=NetPropertiesData(i,:);
    MNetPropertiesData_iP=repmat(NetPropertiesData_i, nP, 1);
    MNetPropertiesData_iA=repmat(NetPropertiesData_i, nA, 1);
    MNetPropertiesData_iL=repmat(NetPropertiesData_i, L, 1);
    
    PolServP=[MNetPropertiesData_iP degreeP abundP sVisits_perP sPolDepP];
    PolServA=[MNetPropertiesData_iA degreeA EffDegA abundA sVisits_perPbyA sPolDepA];
    EffortAPnVars=[MNetPropertiesData_iL degAnP EffortAP abundA(fA,:) abundP(fP,:) sVisits_perPbyA(fA,:) sVisits_perP(fP,:) sPolDepA(fA,:) sPolDepP(fP,:)];

    PolServP_All=[PolServP_All; PolServP];
    PolServA_All=[PolServA_All; PolServA];
    EffortAP_All=[EffortAP_All; EffortAPnVars];
    
end

ColumnNames_NetProperties = NetProperties.Properties.VariableNames;

% Write table of variables for plant species
extraColNames_PolServP = {'degreeP', 'abundP0', 'abundP1', 'sVisits_perP0', 'sVisits_perP1', 'sPollenDepP0', 'sPollenDepP1'};
ColNames_PolServP = [ColumnNames_NetProperties, extraColNames_PolServP];
PolServP = array2table(PolServP_All);
PolServP.Properties.VariableNames = ColNames_PolServP;
writetable(PolServP, sprintf('PollinationServicesReceivedByP_muAP%d_800m.csv',muAP))

% Write table of variables for animals species
extraColNames_PolServA = {'degreeA', 'effective_degree0', 'effective_degree1', 'abundA0', 'abundA1','sVisits_perPbyA0',...
    'sVisits_perPbyA1','sPollenDep0', 'sPollenDep1'};
ColNames_PolServA = [ColumnNames_NetProperties, extraColNames_PolServA];
PolServA = array2table(PolServA_All);
PolServA.Properties.VariableNames = ColNames_PolServA;
writetable(PolServA, sprintf('PolServByA_muAP%d_800m.csv',muAP))

% Write table of variable for plant-animal interactions
EffortAP = array2table(EffortAP_All);
extraColNames_EffortAP = {'degreeA', 'degreeP', 'EffortAP0', 'EffortAP1', 'abundA0', 'abundA1', 'abundP0', 'abundP1', 'sVisits_perPbyA0', 'sVisits_perPbyA1',...
     'sVisits_perP0', 'sVisits_perP1', 'sPollenDepByA0', 'sPollenDepByA1', 'sPollenDepP0', 'sPollenDepP1'};
ColNames_EffortAP = [ColumnNames_NetProperties, extraColNames_EffortAP];
EffortAP.Properties.VariableNames = ColNames_EffortAP;
writetable(EffortAP, sprintf('EffortAtoP_muAP%d_800m.csv',muAP))

