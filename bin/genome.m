function [time_slices, y] = genome(t_max, concs0)

%% Simulation parameters
% These are constants that affect the simulation but do not make
% assertions about the actual system.
n_x = 17;   % number of compartments
n_time_slices = 100;
%t_max=100.0;
Cmax=11;
Pmax=500;
Gmax=250;
carbon_precipitation = 0.1;
diffusion_constant=1;
oxygen_bubble_rate=0;
oxygen_source=1;
carbon_source=10;
nitrogen_source=0.0;
%Currently there isn't any release of ammonia upon degredation
%Maybe add this somehow with this from the amount of C removed 
%or the amount of primary oxidaion with rates 1-5
nitrogen_ratio=0.1;
lambda=0.001; %from Reed et al table S2
D_cell_plus=0.001;
D_cell_minus=0.001;

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
cox_hs_inhib=1; %This is arbitrary non-zero number to make ma_op_inhib_eq 1
nar_sp_growth_rate=0.151;
nar_half_sat_C=0.7;
nar_half_sat_N=0.3;
nar_hs_inhib=15; %%based on dsr inhib half sat constant
nir_sp_growth_rate=0.247;
nir_half_sat_C=0.7;
nir_half_sat_N=0.3;
nir_hs_inhib=15; %%based on dsr inhib half sat constant
nrf_sp_growth_rate=0.162;
nrf_half_sat_C=0.7;
nrf_half_sat_N=0.3;
nrf_hs_inhib=15; %based on dsr inhib half sat constant
dsr_sp_growth_rate=0.0636;
dsr_half_sat_C=0.7;
dsr_half_sat_S=3;
dsr_hs_inhib=15;
amo_sp_growth_rate=0.432;
amo_half_sat_N=107;
amo_half_sat_O=18.75;
amo_hs_inhib=1; %arbitrary non-zero number to make ma_op_inhib_eq 1
hzo_sp_growth_rate=0.864;
hzo_half_sat_Nm=5;
hzo_half_sat_N=5;
hzo_hs_inhib=0.2; %From Reed et al
nor_sp_growth_rate=0.432;
nor_half_sat_N=64.3;
nor_half_sat_O=16.9;
nor_hs_inhib=1; %arbitrary non-zero number to make ma_op_inhib_eq 1
sox_sp_growth_rate=0.4; %changed to 0.4 because it's dominating the community
sox_half_sat_S=0.121;
sox_half_sat_O=0.121;
sox_hs_inhib=1; %arbitrary non-zero number to make ma_op_inhib_eq 1
nap_sp_growth_rate=0.864;
nap_half_sat_N=0.121;
nap_half_sat_S=0.121;
nap_hs_inhib=0.2; %based on hzo inhibition hs constant
%Test with sulfur disproport
sdp_sp_growth_rate=0.864;
sdp_half_sat_H=0.121;
sdp_half_sat_S=0.121;
sdp_hs_inhib=0.2; %based on hzo inhibition hs constant
%The deltaG0 are
%in kJ/mol
%G from 298K
%Calculated as G from CHNOSZ in R with
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
sdp_deltaG0=120.51; %From Aquatic Geomicrobiology vol 48 Canfield
% mass action, one product reactions (ma_op_rxns)
ma_op_rxns = [
    % all reactions
    % From Reed et al PNAS (2014) 111(5):1879–1884
    % aA + bB + cC -> dD + eE + fF 
    % (columns: a, A, b, B, c, C, d, D, e, E, f, F, specific_growth_rate, half_sat_donor, half_sat_acceptor, half_sat_inhib, inhib_index, deltaG0)
    % specific growth rates (days-1)
    % half_sat (uM or nM as in SI of Reed et al)
    % deltaG0 (kJ/mol)
    % separated by three spaces
    0.167   s('C')   1   s('O')   1   s('null')   1   s('CO2')   1   s('H2O')   1   s('null')   cox_sp_growth_rate   cox_half_sat_C   cox_half_sat_O   cox_hs_inhib   s('zero')   cox_deltaG0
    0.167   s('C')   2   s('N+')   1   s('null')   1   s('CO2')   1   s('N')   1   s('H2O')   nar_sp_growth_rate   nar_half_sat_C   nar_half_sat_N   nar_hs_inhib   s('O')   nar_deltaG0
    0.167   s('C')   1.333   s('N')   1.333   s('H')   1   s('CO2')   0.67   s('N2')   1.67   s('H2O')   nir_sp_growth_rate   nir_half_sat_C   nir_half_sat_N   nir_hs_inhib   s('O')   nir_deltaG0
    0.167   s('C')   0.67   s('N')   1.333   s('H')   1   s('CO2')   0.67   s('N-')   0.333   s('H2O')   nrf_sp_growth_rate   nrf_half_sat_C   nrf_half_sat_N   nrf_hs_inhib   s('O')   nrf_deltaG0
    0.167   s('C')   0.5   s('S+')   1   s('null')   1   s('HCO3')   0.5   s('S-')   1   s('null')   dsr_sp_growth_rate   dsr_half_sat_C   dsr_half_sat_S   dsr_hs_inhib   s('O')   dsr_deltaG0
    0.25   s('S-')   1   s('N+')   1   s('null')   1   s('N')   0.25   s('S+')   0.5   s('H')   nap_sp_growth_rate   nap_half_sat_S   nap_half_sat_N   nap_hs_inhib   s('O')   nap_deltaG0
    1   s('S-')   2   s('O')   2   s('HCO3')   2   s('CO2')   1   s('S+')   1   s('H2O')   sox_sp_growth_rate   sox_half_sat_S   sox_half_sat_O   sox_hs_inhib   s('zero') sox_deltaG0
    1   s('N-')   1.5   s('O')   1   s('null')   1   s('N')   1   s('H2O')   2   s('H')   amo_sp_growth_rate   amo_half_sat_N   amo_half_sat_O   amo_hs_inhib   s('zero')   amo_deltaG0
    1   s('N-')   1   s('N')   1   s('null')   1   s('N2')   2   s('H2O')   1   s('null')   hzo_sp_growth_rate   hzo_half_sat_Nm   hzo_half_sat_N   hzo_hs_inhib   s('O')   hzo_deltaG0
    2   s('N')   1   s('O')   1   s('null')   2   s('N+')   1   s('null')   1   s('null')   nor_sp_growth_rate   nor_half_sat_N   nor_half_sat_O   nor_hs_inhib   s('zero')   nor_deltaG0
    4   s('H2O')   4   s('S')   1   s('null')   3   s('S-')   1   s('S+')   2   s('H')   sdp_sp_growth_rate   sdp_half_sat_H   sdp_half_sat_S   sdp_hs_inhib   s('O')   sdp_deltaG0   
];
[n_ma_op_rxns, ~] = size(ma_op_rxns);


%% Internal parameters
% These are parameters that simplify the code but make no assertions about
% the mechanics of the simulation or the actual system

% make the precipitation constants
precipitation_constants = zeros([1 n_species]);
for i = 1: size(precipitation_constant_input, 1)
    idx = precipitation_constant_input(i, 1);
    val = precipitation_constant_input(i, 2);

    precipitation_constants(idx) = val;
end

% make new precipitation values
D = diffusion_constant;
D_plus = (1.0 + precipitation_constants) * D;
D_minus = (1.0 - precipitation_constants) * D;

%Change to per second, rather than per year
%D_plus = D_plus/3.1536E7;
%D_minus = D_minus/3.1536E7;

% grab the unchanging columns from the reaction matrix
ma_op_reac1_c = ma_op_rxns(:, 1)';
ma_op_reac1_i = ma_op_rxns(:, 2);
ma_op_reac2_c = ma_op_rxns(:, 3)';
ma_op_reac2_i = ma_op_rxns(:, 4);
ma_op_reac3_c = ma_op_rxns(:, 5)';
ma_op_reac3_i = ma_op_rxns(:, 6);
ma_op_prod1_c = ma_op_rxns(:, 7)';
ma_op_prod1_i = ma_op_rxns(:, 8);
ma_op_prod2_c = ma_op_rxns(:, 9)';
ma_op_prod2_i = ma_op_rxns(:, 10);
ma_op_prod3_c = ma_op_rxns(:, 11)';
ma_op_prod3_i = ma_op_rxns(:, 12);
ma_op_sp_growth_rate = ma_op_rxns(:, 13)';
ma_op_half_sat_1 = ma_op_rxns(:, 14)';
ma_op_half_sat_2 = ma_op_rxns(:, 15)';
ma_op_hs_inhib = ma_op_rxns(:, 16)';
ma_op_hs_inhib_i = ma_op_rxns(:, 17)';
ma_op_deltaG0 = ma_op_rxns(:, 18)';

% Initiate community structure
% The community structure will be stored
% as a matrix with the following information
% Genes to mediate any of the 10 ma_op reactions? 0, no or 1, yes
% Also, the total number of genes in the genome
% Initiate all with 1 copy
% matrix structre [genome_copies, genome_length, cox, nar, nir, nrf, dsr nap, sox, amo, hzo, nor]
% Foreach row from 1-Cmax
% genome_copies = 1 to start
div_mat=zeros(11, 13);
%Set up temp structure for now with each function represented onces
div_mat(1,3)=1;
div_mat(2,4)=1;
div_mat(3,5)=1;
div_mat(4,6)=1;
div_mat(5,7)=1;
div_mat(6,8)=1;
div_mat(7,9)=1;
div_mat(8,10)=1;
div_mat(9,11)=1;
div_mat(10,12)=1;
%turn off sdp by removing genes from the community
div_mat(11,13)=0;
%This will be for random structure
%for x = 1: Cmax
%   %set the copies to 1;
%   div_mat(x, 1)=1;
%   genome_length = randi(Gmax);
%   div_mat(x, 2)=genome_length;
%   genome_compo = randi(Pmax,genome_length,1);
%   % set up the rows of div_mat according to genome composition
%   %glut1 is gene 13
%   if ismember(13,genome_compo)
%      %cox is gene 1
%      if ismember(1,genome_compo)
%         div_mat(x,3)=1;
%      end
%      %nar is gene 2
%      if ismember(2,genome_compo)
%         div_mat(x,4)=1;
%      end
%      %nir is gene 3
%      if ismember(3, genome_compo)
%         div_mat(x,5)=1;
%      end
%      %nrf is gene 4
%     if ismember(4, genome_compo)
%         div_mat(x, 6)=1;
%      end
%      %dsr is gene 5
%      if ismember(5, genome_compo)
%         div_mat(x, 7)=1;
%      end
%   end
%   %rbcl is gene 6
%   if ismember(6, genome_compo)
%      %nap is gene 7 
%      if ismember(7, genome_compo)
%         div_mat(x, 8)=1;
%      end
%      %sox is gene 8
%      if ismember(8, genome_compo)
%         div_mat(x, 9)=1;
%      end
%      %amoA is gene 9
%      if ismember(9, genome_compo)
%         div_mat(x, 10)=1;
%      end
%      %hzo is gene 10
%      if ismember(10, genome_compo)
%         div_mat(x,11)=1;
%      end
%      %nor is gene 11
%      if ismember(11, genome_compo)
%         div_mat(x, 12)=1;
%      end
%   end
%end

n_total_chem = n_x * n_species;
%for now diversity remains constant so div_mat doesn't change, but the number of organisms of each does
n_total_div = n_x * Cmax;
div0=ones(n_x, Cmax);
%% Define the flux functions

% -- rates --
% This function takes a row from the concentration matrix (i.e., a
% horizontal slice from the lake) and computes the rates of the mass action
% and primary oxidation reactions
function [ma_op_rates, ma_op_deltaG, Y] = rates(concs_row, Gamma)
    %to change this to dynamic community, we need to pass comunity too
    %Community will change like concs will change
    %it will be just as dynamic 
    % compute the mass action rates
    % for chemicals
    %assume unchaning community for now
    % use equations 1 and 3 in Reed et al
    % get the concentration of chemical species
    % for both reactants and products
    %Here null should be 1 so it doesn't affect Q calc
%    concs_row(s('null'))=ones(n_x);
%        concs_row=concs(x, :);
        %change the null and H2O to 1M so it doesn't affect Q calc.
	%1E6/1E6=1 see below it is divided by 1E6 to make into moles
        concs_row(s('null'))=1E6;
	concs_row(s('H2O'))=1E6;	
	%Change the H so it is buffered at pH 7.5
	concs_row(s('H'))=0.0316;
	%N2 and CO2, HCO3 should remain steady
	%I'm not sure what to keep it at for now, but these should remain steady
	%Approximations for now, but get refs and confirm in future
	%From sat of 29.14 mg/L N2
	concs_row(s('N2'))=1E3;
	%From 200 mg/L alk and 1.2192 mg/L of HCO3
	concs_row(s('CO2'))=4E3;
	concs_row(s('HCO3'))=4E3;

        ma_op_reac1 = concs_row(ma_op_reac1_i);
        ma_op_reac2 = concs_row(ma_op_reac2_i);
        ma_op_reac3 = concs_row(ma_op_reac3_i);
        ma_op_prod1 = concs_row(ma_op_prod1_i);
        ma_op_prod2 = concs_row(ma_op_prod2_i);
        ma_op_prod3 = concs_row(ma_op_prod3_i);
	ma_op_inhib = concs_row(ma_op_hs_inhib_i);
        % Calculate deltaG
        % R is gas constant kJ K-1 mole -1
        % room temp in Kelvin
        % Equation SI5 in Reed et al
	numerator=times(times(power(ma_op_prod3/1E6,ma_op_prod3_c), power(ma_op_prod2/1E6,ma_op_prod2_c)), power(ma_op_prod1/1E6,ma_op_prod1_c));
	denominator=times(times(power(ma_op_reac3/1E6,ma_op_reac3_c), power(ma_op_reac2/1E6,ma_op_reac2_c)), power(ma_op_reac1/1E6, ma_op_reac1_c));
	Q=rdivide(numerator,denominator);
	g = 0.0083144598 *298*log(Q);
	ma_op_deltaG = plus(ma_op_deltaG0, g);
    
        %calculate the biomass yield for each reaction related to free energy
        %This is not per mole, but per uM
	%Check if that is correct
	%Equation SI4 of Reed et al 2016
        Y = 2.08 - 0.0211*rdivide(ma_op_deltaG,ma_op_reac1);

        % Ft calculated from each sample
	%Thermodynamic potential factor (unitless)
	%Equation S1 of Reed et al
	Ft=rdivide(1,exp(rdivide((ma_op_deltaG + 96.485*0.120),0.083144598*298))+1);
        % foreach genome, figure out which reaction is most favorable given the deltaG calc above and their genome content
        % something like go through this and figure out for each how many

        % Then calculate gamma (i.e. the number of specific genones)
        % Gamma = number of individuals cycling each reaction
        % For now assume all individuals with capability contribute 
        % Now use gamma to calculate the rate of the reaction
        %Gamma is a vector of the sum of genes for each reaction
	%Based on equation 1 of Reed et al
	half_sat_eq_1 = rdivide(ma_op_reac1,plus(ma_op_reac1,ma_op_half_sat_1));
	half_sat_eq_2 = rdivide(ma_op_reac2,plus(ma_op_reac2,ma_op_half_sat_2));
	hs_inhib_eq = rdivide(ma_op_hs_inhib,plus(ma_op_inhib,ma_op_hs_inhib));
        ma_op_rates=times(Gamma,times(Ft,times(ma_op_sp_growth_rate, times(half_sat_eq_1, times(half_sat_eq_2, hs_inhib_eq)))));
%        ma_op_rates_mat(x, :)=ma_op_rates;

end

% -- flux --
% This function is taken as an argument by the ODE solver, which feeds this
% function the concentration matrix (flattened into a concentration
% vector). The function computes the rates of reactions, diffusions, and
% precipitations and outputs the fluxes for each metabolite at each depth.

% twiddle because this is time-independent
function [merged_fluxes] = flux(~, merged_vector)
    % Extractconcs from concs_vector
    %Divide vector into concs_vector and div_vector
    %
    concs_vector = merged_vector(1:n_total_chem);
    div_vector = merged_vector(n_total_chem+1:end);
    concs = reshape(concs_vector, [n_x, n_species]);
    div = reshape(div_vector, [n_x, Cmax]);

    conc_fluxes = zeros(n_x, n_species);
    div_fluxes = zeros(n_x, Cmax);
    % apply the oxygen bubbles
    conc_fluxes(:, s('O')) = conc_fluxes(:, s('O')) + oxygen_bubble_rate;

    % apply the fixed source terms
    conc_fluxes(1, s('O')) = conc_fluxes(1, s('O')) + oxygen_source;
    conc_fluxes(1, s('C')) = conc_fluxes(1, s('C')) + carbon_source;
    conc_fluxes(1, s('N+')) = conc_fluxes(1, s('N+')) + nitrogen_source;
    %conc_fluxes(end, s('CH4')) = conc_fluxes(end, s('CH4')) + methane_source;
    
    for x = 1: n_x
%        %set 'null' to 1 for Q calculation, so it doesn't influence deltaG
	%Each species has unchaning a vecotr of genes
	%If the community diversity is dynamic too, it needs to be treated differently
        Gamma=sum(div(x, :)'.*div_mat(:,3:end));
        [ma_op_rates, ma_op_deltaG, Y] = rates(concs(x, :), Gamma);

        % apply the mass action rates
	%Based on Equation 3 of Reed et al ?
	%removing reactants
	%adding products
	%MISSING the stoichiometry for these reactions
	%React 1 is the electron donor, so it doesn't need to be adjusted by stochiometry
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac1_i, rdivide(ma_op_rates,Y), [n_species, 1])';
	%Others need to be adjusted by the relationship between the stoichiometry of 1 and itself
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac2_i, times(rdivide(ma_op_reac2_c,ma_op_reac1_c),rdivide(ma_op_rates,Y)), [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) - accumarray(ma_op_reac3_i, times(rdivide(ma_op_reac3_c,ma_op_reac1_c),rdivide(ma_op_rates,Y)), [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(ma_op_prod1_i, times(rdivide(ma_op_prod1_c,ma_op_reac1_c),rdivide(ma_op_rates,Y)), [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(ma_op_prod2_i, times(rdivide(ma_op_prod2_c,ma_op_reac1_c),rdivide(ma_op_rates,Y)), [n_species, 1])';
        conc_fluxes(x, :) = conc_fluxes(x, :) + accumarray(ma_op_prod3_i, times(rdivide(ma_op_prod3_c,ma_op_reac1_c),rdivide(ma_op_rates,Y)), [n_species, 1])';

	%aaply the growth terms
	%This is based on equation 2 of Reed et al?
	%Right now, metabolic plasticity isn't included 
	%This should change the gene rate based on genome structure
	%For now, for any process that is available in the genome
	%These genes are arrayed in the same way as the genome
	%So for anything that has a gene, it should change based on that rate
	%cycle through all of the organisms and increase their numbers according to gene content and rates
	for gx = 1 : Cmax
		div_fluxes(x, gx) = div_fluxes(x, gx) + sum(times(ma_op_rates,div_mat(gx,3:end))) - lambda*div(x,gx);
	end
        
        % diffusion      
        if x > 1
            conc_fluxes(x, :) = conc_fluxes(x, :) + D_plus .* concs(x - 1, :) - D_minus .* concs(x, :);
            for gx = 1 : Cmax
               div_fluxes(x,gx) = div_fluxes(x, gx) + D_cell_plus .* div(x-1,gx) - D_cell_minus .* div(x, gx);
            end
        end

        if x < n_x
            conc_fluxes(x, :) = conc_fluxes(x, :) - D_plus .* concs(x, :) + D_minus .* concs(x + 1, :);
	    for gx = 1: Cmax
                div_fluxes(x,gx)=div_fluxes(x, gx) - D_cell_plus .* div(x, gx) + D_cell_minus .* div(x+1, gx);
	    end
        end
        %For any carbon that is used for energy, I want to release N
	%This should be for any carbon oxidizing processes
	temp_flux(x, :) = accumarray(ma_op_reac1_i, rdivide(ma_op_rates,Y), [n_species, 1])';
	%16 times amount of oxidized carbon
	conc_fluxes(x, s('N-'))=conc_fluxes(x, s('N-')) + 16*(temp_flux(x, s('C')));
    end % for x

    %Fix buffered values so they don't change
    conc_fluxes(:, s('null'))=0.0;
    conc_fluxes(:, s('H2O'))=0.0;
    conc_fluxes(:, s('H'))=0.0;
    conc_fluxes(:, s('N2'))=0.0;
    conc_fluxes(:, s('CO2'))=0.0;
    conc_fluxes(:, s('HCO3'))=0.0;
    conc_fluxes(:, s('zero'))=0.0;
    % reshape into long vector
    conc_fluxes = reshape(conc_fluxes, [n_total_chem, 1]);
    div_fluxes = reshape(div_fluxes, [n_total_div, 1]);
    merged_fluxes_2 = [conc_fluxes' div_fluxes'];
    merged_fluxes = merged_fluxes_2';
end


%% ODE solver
% This section feeds the flux function to the ODE solver.

% all concentrations are constrained to be nonnegative
n_total=n_total_chem+n_total_div;
options = odeset('NonNegative', 1: n_total);

% initially flatten the concentration matrix
concs0_vector = reshape(concs0, [n_total_chem, 1]);
%check to see if Cmax is the right thing here
div0_vector = reshape(div0, [n_total_div, 1]);
init_vector = [concs0_vector' div0_vector'];
init_vector_2=init_vector';
% run the ODE solver (ode15s)
% t is the times at which the ODE solver gives output. They are not evenly
% spaced! y is a matrix whose rows are the flattened concentration matrices
% at each time step
% add div_mat to concs0_vector
[time_slices, y] = ode15s(@flux, linspace(0.0, t_max, n_time_slices), init_vector_2, options);

% unfold the result y, putting it into a 3D space whose dimensions
% correspond to time, depth, and metabolite
[n_time_slices, ~] = size(y);
%concs_history = reshape(y, n_time_slices, n_x, n_species);

%% get all the reaction rates for all timepoints
% ma_op then teas
%rates_history = zeros(n_time_slices, n_x, n_ma_op_rxns);
%for time = 1: n_time_slices
%    for x = 1: n_x
%        [rates_history(time, x, 1:n_ma_op_rxns), ~, rates_history(time, x, n_ma_op_rxns + 1:end)] = rates(squeeze(concs_history(time, x, :))');
%    end
%end

end
