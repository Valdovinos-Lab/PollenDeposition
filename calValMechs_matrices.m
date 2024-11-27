function [Visits_perP, pol_event]= calValMechs_matrices(Alpha, p, a, network_metadata)

tau     = network_metadata.tau ;
epsilon = network_metadata.epsilon ;

Visits_perP = Alpha * diag(a .*tau) ;

sigma = diag(p.*epsilon) * Alpha ;
sigma = sigma * diag(1./(sum(sigma)+realmin)) ;

pol_event = sigma .* Visits_perP ; % Matrix of pollination events per visit (includes animal abundance)


end