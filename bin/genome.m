function [div_mat, ma_op_rxns] = test_genome(Cmax, Gmax, Pmax, carbon_precipitation)

%% Simulation parameters
% These are constants that affect the simulation but do not make
% assertions about the actual system.
n_x = 17;   % number of compartments
n_time_slices = 100;

%% Species map
% import the species list using the separate function file:
% "s" is a hash from a string that names the species to its index in the
% concentration matrix
[s, ~, n_species] = species_map();

%% Precipitation constants
% assert the precipitation constants: a constant 0 means equal diffusion up
% and down; +0.1 means a molecule is 10% as likely to go down as to go up;
% -0.1 means 10% more likley to go up than down
precipitation_constant_input = [
    s('C'), carbon_precipitation
];

%% Reaction constants
% there are mass action reactions and the primary oxidations
% reactions from Reed et al PNAS (2014) 111(5):1879-1884
% SI table S2
cox_sp_growth_rate=0.28;
cox_half_sat_C=0.7;
cox_half_sat_O=0.121;
nar_sp_growth_rate=0.151;
nar_half_sat_C=0.7;
nar_half_sat_N=0.3;
nir_sp_growth_rate=0.247;
nir_half_sat_C=0.7;
nir_half_sat_N=0.3;
nrf_sp_growth_rate=0.162;
nrf_half_sat_C=0.7;
nrf_half_sat_N=0.3;
dsr_sp_growth_rate=0.0636;
dsr_half_sat_C=0.7;
dsr_half_sat_S=3;
amo_sp_growth_rate=0.432;
amo_half_sat_N=107;
amo_half_sat_O=18.75;
hzo_sp_growth_rate=0.864;
hzo_half_sat_Nm=5;
hzo_half_sat_N=5;
nor_sp_growth_rate=0.432;
nor_half_sat_N=64.3;
nor_half_sat_O=16.9;
sox_sp_growth_rate=0.864;
sox_half_sat_S=0.121;
sox_half_sat_O=0.121;
nap_sp_growth_rate=0.864;
nap_half_sat_N=0.121;
nap_half_sat_S=0.121;
%The deltaG0 are maade up for now
%in kJ/mol
%G from 298K
%Calculated as G from CHNOSZ with
cox_deltaG0=-477.2734; %subcrt(c(1646, 65, 1576, 1), c(-1/6, -1, 1, 1), T=seq(298, 299, 300))
nar_deltaG0=-313.2163; %subcrt(c(1646, 16, 1576, 17, 1), c(-1/6, -2, 1, 2, 1), T=seq(298, 299, 300))
nir_deltaG0=-573.6341; %subcrt(c(1646, 17, 3, 1576, 74, 1), c(-1/6, -4/3, -4/3, 1, 2/3, 5/3), T=seq(298, 299, 300))
nrf_deltaG0=-343.9649; %subcrt(c(1646, 17, 3, 1576, 18, 1), c(-1/6, -2/3, -4/3, 1, 2/3, 1/3), T=seq(298, 299, 300))
dsr_deltaG0=-76.10632; %subcrt(c(1646, 24, 13, 67), c(-1/6, -1/2, 1, 1/2), T=seq(298, 299, 300))
nap_deltaG0=-100.4539; %subcrt(c(67, 65, 13, 1576, 24, 1), c(-1, -2, -2, 2, 1, 1), T=seq(298, 299, 300))
sox_deltaG0=-822.0793; %subcrt(c(18, 65, 17, 1, 3), c(-1, -3/2, 1, 1, 2), T=seq(298, 299, 300));
amo_deltaG0=-214.7716; %subcrt(c(18, 17, 74, 1), c(-1, -1, 1, 2), T=seq(298, 299, 300))
hzo_deltaG0=-344.5038; %subcrt(c(18, 17, 74, 1), c(-1, -1, 1, 2), T=seq(298, 299, 300))
nor_deltaG0=-173.9297; %subcrt(c(65, 17, 16), c(-1, -2, 2), T=seq(298, 299, 300))

% mass action, one product reactions (ma_op_rxns)
ma_op_rxns = [
    % all reactions
    % From Reed et al PNAS (2014) 111(5):1879–1884
    % aA + bB + cC -> dD + eE + fF 
    % (columns: a, A, b, B, c, C, d, D, e, E, f, F, specific_growth_rate, half_sat_donor, half_sat_acceptor, deltaG0)
    % specific growth rates (days-1)
    % half_sat (uM or nM as in SI of Reed et al)
    % deltaG0 (kJ/mol)
    % separated by three spaces
    0.167   s('C')   1   s('O')   1   s('null')   1   s('CO2')   1   s('H2O')   1   s('null')   cox_sp_growth_rate   cox_half_sat_C   cox_half_sat_O   cox_deltaG0
    0.167   s('C')   2   s('N+')   1   s('null')   1   s('CO2')   1   s('N')   1   s('H2O')   nar_sp_growth_rate   nar_half_sat_C   nar_half_sat_N   nar_deltaG0
    0.167   s('C')   1.333   s('N')   1.333   s('H')   1   s('CO2')   0.67   s('N2')   1.67   s('H2O')   nir_sp_growth_rate   nir_half_sat_C   nir_half_sat_N   nir_deltaG0
    0.167   s('C')   0.67   s('N-')   1.333   s('H')   1   s('CO2')   0.67   s('N+')   0.333   s('H2O')   nrf_sp_growth_rate   nrf_half_sat_C   nrf_half_sat_N   nrf_deltaG0
    0.167   s('C')   0.5   s('S+')   1   s('null')   1   s('HCO3')   0.5   s('S-')   1   s('null')   dsr_sp_growth_rate   dsr_half_sat_C   dsr_half_sat_S   dsr_deltaG0
    0.25   s('S-')   1   s('N+')   1   s('null')   1   s('N')   0.25   s('S+')   0.5   s('H')   nap_sp_growth_rate   nap_half_sat_S   nap_half_sat_N   nap_deltaG0
    1   s('S-')   2   s('O')   2   s('HCO3')   2   s('CO2')   1   s('S+')   1   s('H2O')   sox_sp_growth_rate   sox_half_sat_S   sox_half_sat_O   sox_deltaG0
    1   s('N-')   1.5   s('O')   1   s('null')   1   s('N')   1   s('H2O')   2   s('H')   amo_sp_growth_rate   amo_half_sat_N   amo_half_sat_O   amo_deltaG0
    1   s('N-')   1   s('N')   1   s('null')   1   s('N2')   2   s('H2O')   1   s('null')   hzo_sp_growth_rate   hzo_half_sat_Nm   hzo_half_sat_N   hzo_deltaG0
    2    s('N')   1   s('O')   1   s('null')   2   s('N+')   1   s('null')   1   s('null')   nor_sp_growth_rate   nor_half_sat_N   nor_half_sat_O   nor_deltaG0
];
%[n_ma_op_rxns, ~] = size(ma_op_rxns);


%% Internal parameters
% These are parameters that simplify the code but make no assertions about
% the mechanics of the simulation or the actual system

% make the precipitation constants
%precipitation_constants = zeros([1 n_species]);
%for i = 1: size(precipitation_constant_input, 1)
%    idx = precipitation_constant_input(i, 1);
%    val = precipitation_constant_input(i, 2);
%
%%    precipitation_constants(idx) = val;
%end

% make new precipitation values
%D = diffusion_constant;
%D_plus = (1.0 + precipitation_constants) * D;
%D_minus = (1.0 - precipitation_constants) * D;

% grab the unchanging columns from the reaction matrix
ma_op_reac1_c = ma_op_rxns(:, 1)';
ma_op_reac1_i = ma_op_rxns(:, 2);
ma_op_reac2_c = ma_op_rxns(:, 3)';
ma_op_reac2_i = ma_op_rxns(:, 4);
ma_op_reac3_c = ma_op_rxns(:, 5)’;
ma_op_reac3_i = ma_op_rxns(:, 6);
ma_op_prod1_c = ma_op_rxns(:, 7)’;
ma_op_prod1_i = ma_op_rxns(:, 8);
ma_op_prod2_c = ma_op_rxns(:, 9)’;
ma_op_prod2_i = ma_op_rxns(:, 10);
ma_op_prod3_c = ma_op_rxns(:, 11)’;
ma_op_prod3_i = ma_op_rxns(:, 12);
%% are these supposed to be transformed?
% the rc is transposed in the original
ma_op_sp_growth_rate = ma_op_rxns(:, 13)’;
ma_op_half_sat_1 = ma_op_rxns(:, 14)’;
ma_op_half_sat_2 = ma_op_rxns(:, 15)’;
ma_op_deltaG0 = ma_op_rxns(:, 16)’;

% Initiate community structure
% The community structure will be stored
% as a matrix with the following information
% Genes to mediate any of the 10 ma_op reactions? 0, no or 1, yes
% Also, the total number of genes in the genome
% Initiate all with 1 copy
% matrix structre [genome_copies, genome_length, cox, nar, nir, nrf, dsr nap, sox, amo, hzo, nor]
% Foreach row from 1-Cmax
% genome_copies = 1 to start

div_mat=zeros(Cmax, 11);

for x = 1: Cmax
   %set the copies to 1;
   div_mat(x, 1)=1;
   genome_length = randi(Gmax);
   div_mat(x, 2)=genome_length;
   genome_compo = randi(Pmax,genome_length,1);
   % set up the rows of div_mat according to genome composition
   %glut1 is gene 13
   if ismember(13,genome_compo)
      %cox is gene 1
      if ismember(1,genome_compo)
         div_mat(x,3)=1;
      end
      %nar is gene 2
      if ismember(2,genome_compo)
         div_mat(x,4)=1;
      end
      %nir is gene 3
      if ismember(3, genome_compo)
         div_mat(x,5)=1;
      end
      %nrf is gene 4
      if ismember(4, genome_compo)
         div_mat(x, 6)=1;
      end
      %dsr is gene 5
      if ismember(5, genome_compo)
         div_mat(x, 7)=1;
      end
   end
   %rbcl is gene 6
   if ismember(6, genome_compo)
      %nap is gene 7 
      if ismember(7, genome_compo)
         div_mat(x, 8)=1;
      end
      %sox is gene 8
      if ismember(8, genome_compo)
         div_mat(x, 9)=1;
      end
      %amoA is gene 9
      if ismember(9, genome_compo)
         div_mat(x, 10)=1;
      end
      %hzo is gene 10
      if ismember(10, genome_compo)
         div_mat(x,11)=1;
      end
      %nor is gene 11
      if ismember(11, genome_compo)
         div_mat(x, 12)=1;
      end
   end
end

% choose a rand set of numbers from 1-Pmax (max size of gene pool)

%% Define the flux functions

% -- rates --
% This function takes a row from the concentration matrix (i.e., a
% horizontal slice from the lake) and computes the rates of the mass action
% and primary oxidation reactions

n_total = n_x * n_species;

function [ma_op_rates, po_carbon_rate] = rates(concs_row)
    % compute the mass action rates
    % for both chemicals and genomes
    % equations 1 and 3 in Reed et al
    % get the concentration of chemical species
    % for both reactants and products
    ma_op_reac1 = concs_row(ma_op_reac1_i);
    ma_op_reac2 = concs_row(ma_op_reac2_i);
    ma_op_reac3 = concs_row(ma_op_reac3_i);
    ma_op_prod1 = concs_row(ma_op_prod1_i);
    ma_op_prod2 = concs_row(ma_op_prod2_i);
    ma_op_prod3 = concs_row(ma_op_prod3_i);

    % Calculate deltaG
    % R is gas constant kJ K-1 mole -1
    % room temp in Kelvin
    %%Q%% can I multiply and add the vectors by these numbers?

    ma_op_deltaG = ma_op_deltaG0 + 8.3144598 *273*ln((ma_op_prod1*ma_op_prod2*ma_op_prod3)/(ma_op_reac1*ma_op_reac2*ma_op_reac3));
    
    %%Q%% can I multiple the vector (ma_op_deltaG) by a number (i.e. 2.08)?
    Y = 2.08 - 0.0211*ma_op_deltaG;

    % Ft calculated from each sample
    %%Q%% Check if this is doing what I Think it is
    Ft = 1/(exp((deltaG + 96.485*0.120)/(8.3144598*273))+1);

    % foreach genome, figure out which reaction is most favorable given the deltaG calc above and their genome content
    % something like go through this and figure out for each how many

    % Then calculate gamma (i.e. the number of specific genones)
    % Gamma = number of genomes cycling each reaction

    % Now use gamma to calculate the rate of the reaction
    ma_op_rates = Gamma*Ft*ma_op_sp_growth*(ma_op_reac1/(ma_op_reac1 + ma_op_half_sat_1))*(ma_op_reac2/(ma_op_reac2 + ma_op_half_sat_2));


end

% -- flux --
% This function is taken as an argument by the ODE solver, which feeds this
% function the concentration matrix (flattened into a concentration
% vector). The function computes the rates of reactions, diffusions, and
% precipitations and outputs the fluxes for each metabolite at each depth.

% twiddle because this is time-independent
function [conc_fluxes] = flux(~, concs_vector)
    % Extract both concs and genomes from concs_vector
    concs = reshape(concs_vector, [n_x, n_species]);
    %something like this, where concs in line above replaced by concs_extract
    % concs=concs_extract(1,n_chemicals);
    % genomes=concs_extraction(n_chemicals+1,end);

    conc_fluxes = zeros(n_x, n_species);
    % genome_fluxes = zeros(n_x, n_genomes)

    % apply the oxygen bubbles
    conc_fluxes(:, s('O')) = conc_fluxes(:, s('O')) + oxygen_bubble_rate;

    % apply the fixed source terms
    conc_fluxes(1, s('O')) = conc_fluxes(1, s('O')) + oxygen_source;
    conc_fluxes(1, s('C')) = conc_fluxes(1, s('C')) + carbon_source;
    conc_fluxes(1, s('N+')) = conc_fluxes(1, s('N+')) + nitrogen_source;
    conc_fluxes(end, s('CH4')) = conc_fluxes(end, s('CH4')) + methane_source;
    
    for x = 1: n_x
        [ma_op_rates, po_carbon_rate, tea_rates] = rates(concs(x, :));

        % apply the mass action rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac1_i, ma_op_reac1_c .* ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac2_i, ma_op_reac2_c .* ma_op_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(ma_op_prod_i, ma_op_prod_c .* ma_op_rates, [n_species, 1])';

        % apply the primary oxidation rates
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(po_tea_i, tea_rates, [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(po_tea_prod_i, tea_rates, [n_species, 1])';
        
        conc_fluxes(x, s('C')) = conc_fluxes(x, s('C')) - po_carbon_rate;
        conc_fluxes(x, s('N-')) = conc_fluxes(x, s('N-')) + nitrogen_ratio * po_carbon_rate;
        
        % diffusion      
        if x > 1
            conc_fluxes(x, :) = conc_fluxes(x, :) + D_plus .* concs(x - 1, :) - D_minus .* concs(x, :);
        end

        if x < n_x
            conc_fluxes(x, :) = conc_fluxes(x, :) - D_plus .* concs(x, :) + D_minus .* concs(x + 1, :);
        end

    end % for x
    
    conc_fluxes(:, s('null')) = 0.0;
    conc_fluxes = reshape(conc_fluxes, [n_total, 1]);
end


%% ODE solver
% This section feeds the flux function to the ODE solver.

% all concentrations are constrained to be nonnegative
options = odeset('NonNegative', 1: n_total);

% initially flatten the concentration matrix
concs0_vector = reshape(concs0, [n_total, 1]);

% run the ODE solver (ode15s)
% t is the times at which the ODE solver gives output. They are not evenly
% spaced! y is a matrix whose rows are the flattened concentration matrices
% at each time step
[time_slices, y] = ode15s(@flux, linspace(0.0, t_max, n_time_slices), concs0_vector, options);

% unfold the result y, putting it into a 3D space whose dimensions
% correspond to time, depth, and metabolite
[n_time_slices, ~] = size(y);
concs_history = reshape(y, n_time_slices, n_x, n_species);

%% get all the reaction rates for all timepoints
% ma_op then teas
rates_history = zeros(n_time_slices, n_x, n_ma_op_rxns + n_po_teas);
for time = 1: n_time_slices
    for x = 1: n_x
        [rates_history(time, x, 1:n_ma_op_rxns), ~, rates_history(time, x, n_ma_op_rxns + 1:end)] = rates(squeeze(concs_history(time, x, :))');
    end
end

%end
