function [s, species, n_species] = species_map()

% specify the names of all the species and use those as keys that uses a
% map 's' to link them to integers. these integers will be their position
% in the reaction matrices, etc.
species = {'O', 'C', 'N+', 'N', 'N-', 'S+', 'S', 'S-', 'H', 'H2O', 'N2', 'CO2', 'HCO3', 'null', 'zero'};
n_species = length(species);
s = containers.Map(species, 1: n_species);

end
