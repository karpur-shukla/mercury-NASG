dodecane_data = table2array(readtable('./../Dodecaneproperties/Dodecane_NIST_CSV_Vukalovich_order_ASCII.csv'));
%save('./Hgproperties/Vukalovich_Saturation_Line_by_Temperature_ASCII.mat',"water_data_all");

%%
% Defining the functions used for the liquid phase parameters
data_size = 201; 

syms f1_all(p_inf_liq) f2_all(p_inf_liq) f3_all(p_inf_liq);
syms A_term_all(p_inf_liq) B_term_all(p_inf_liq);

% parameters from ideal gas approx (Pg. 9 of LeMétayer-Saurel, below Eq. 45 and 46)
b_gas_all = 0;
p_inf_gas_all = 0;

% Defining the average temperature, pressure, average specific enthalpy (gas and liquid) for both the full and the ranges
T_avg_all = mean(dodecane_data(:,1));
p_avg_all = mean(dodecane_data(:,2));
h_gas_avg_all = mean(dodecane_data(:,8));
h_liq_avg_all = mean(dodecane_data(:,7));
v_liq_avg_all = mean(dodecane_data(:,4));

%c_3 constants
T_avg_norm_factor_all = sum(dodecane_data(:,1) .* (dodecane_data(:,1) - T_avg_all));

% Using Eq. 50, 51, 55 in LeMétayer-Saurel to calculate c_p,g; q_g; c_v,g; gamma_g for both full data set and range. c_p,g in J/(g K); c_v,g in m^3/kg.
c_p_gas_all = sum(dodecane_data(:,1) .* (dodecane_data(:,8) - h_gas_avg_all))/T_avg_norm_factor_all;
q_gas_all = h_gas_avg_all - c_p_gas_all .* T_avg_all;
c_v_gas_all = c_p_gas_all - sum(dodecane_data(:,1) .* dodecane_data(:,5) ./ dodecane_data(:,2))/sum((dodecane_data(:,1) ./ dodecane_data(:,2)) .^ 2);
gamma_gas_all = c_p_gas_all/c_v_gas_all;

% reference state: 20 C point on saturation curve (using Vukalovich data). Density taken from Wikipedia.
reference_p = 1.128;
reference_dens = 589.73;
reference_sound_speed_sq = (620.4).^2;
reference_dens_sound_speed_sq = reference_dens .* reference_sound_speed_sq;
ideal_gas_const = 8.31446261815324./100000;

A_term_relative_coeff_all = 1 - (reference_dens .* v_liq_avg_all);

% Setting up the root-finding problem using LeMétayer-Saurel Eq. 60, 61, 64, 65, 68.
% c_1 coefficients
c_p_liq_numer_1st_term_all = sum(dodecane_data(:,1) .* (dodecane_data(:,7) - h_liq_avg_all)); %c_1 (all)

% c_2 coefficients
c_p_liq_numer_2nd_term_coeff_all = sum(dodecane_data(:,1) .* (dodecane_data(:,2) - p_avg_all)); %c_2 (all)

% Constant term in the expanded version of Eq 68 in LeMS
LeMS_Eq_68_const_term_all = reference_p - reference_dens_sound_speed_sq + (((reference_dens).^2) .* reference_sound_speed_sq .* v_liq_avg_all);

% 
f1_all = 0;
f2_all = 0;
f3_all = 0;

for i=1:data_size
  f1_all = f1_all + (vpa(dodecane_data(i,1)) * (vpa(dodecane_data(i,4)) - v_liq_avg_all)/(vpa(dodecane_data(i,2)) + p_inf_liq));
  f2_all = f2_all + (1/data_size) * (vpa(dodecane_data(i,1))/(vpa(dodecane_data(i,2)) + p_inf_liq));
  f3_all = f3_all + vpa(dodecane_data(i,1))/(vpa(dodecane_data(i,2)) + p_inf_liq) * (vpa(dodecane_data(i,1))/(vpa(dodecane_data(i,2)) + p_inf_liq) - (1/data_size) * (vpa(dodecane_data(i,1))/(vpa(dodecane_data(i,2)) + p_inf_liq)));
end

clear i;

A_term_all = (T_avg_norm_factor_all .* f1_all)/(c_p_liq_numer_1st_term_all - (c_p_liq_numer_2nd_term_coeff_all .* v_liq_avg_all .* f3_all) + (c_p_liq_numer_2nd_term_coeff_all .* f1_all .* f2_all));
B_term_all = (f1_all .* f2_all)/(f3_all);


%
v_gas_analytic = zeros(data_size,1);
h_gas_analytic = zeros(data_size,1);
compress_factor_gas_analytic = zeros(data_size,1);
dens_gas_analytic = zeros(data_size,1);
dens_gas_vukalovich = zeros(data_size,1);

for i=1:data_size
  v_gas_analytic(i) = ((gamma_gas_all - 1) * c_v_gas_all * dodecane_data(i,1))/(dodecane_data(i,2));
  compress_factor_gas_analytic(i) = (dodecane_data(i,2) .* v_gas_analytic(i)) ./ (ideal_gas_const .* dodecane_data(i,1) .* 4.9856360761414);
  h_gas_analytic(i) = gamma_gas_all * c_v_gas_all * dodecane_data(i,1) + q_gas_all;
  dens_gas_analytic(i) = 1 ./ v_gas_analytic(i);
  dens_gas_vukalovich(i) = 1 ./ dodecane_data(i,5);
end

clear i;


%%
figure(1);
  semilogy(dodecane_data(:,1),dodecane_data(:,5),'.','Color',[0 0.4470 0.7410],'MarkerSize',16);
  hold on;
  semilogy(dodecane_data(:,1),v_gas_analytic,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
  fontsize(15,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${v}_{g} \left( \mathrm{{m}^{3} / kg} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(NIST)','NASG gas phase' + string(newline) + 'theoretical prediction','Location','northeast','Interpreter','latex','FontSize',14);
  set(gca,'TickLabelInterpreter','latex');

figure(2);
  semilogy(dodecane_data(:,1),dens_gas_vukalovich,'.','Color',[0.3010 0.7450 0.9330],'MarkerSize',16);
  hold on;
  semilogy(dodecane_data(:,1),dens_gas_analytic,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
  fontsize(15,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${\rho}_{g} \left( \mathrm{kg / {m}^{3}} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(NIST)','NASG gas phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',14);
  set(gca,'TickLabelInterpreter','latex');

figure(3);
  plot(dodecane_data(:,1),dodecane_data(:,8),'.','Color',[0.4940 0.1840 0.5560],'MarkerSize',16);
  hold on;
  plot(dodecane_data(:,1),h_gas_analytic,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(15,"points");
  %title('\bf{Specific Enthalpy vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${h}_{g} \left( \mathrm{kJ/kg} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(NIST)','NASG gas phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',14);
  set(gca,'TickLabelInterpreter','latex');

figure(5);
  plot(dodecane_data(:,1),dodecane_data(:,6),'.','Color',[0.7 0 1],'MarkerSize',16);
  hold on;
  plot(dodecane_data(:,1),compress_factor_gas_analytic,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(15,"points");
  %title('Compressibility Factor vs. Temperature, Mercury Gas, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${z}_{g}$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(NIST)','NASG gas phase' + string(newline) + 'theoretical prediction','Location','southwest','Interpreter','latex','FontSize',14);
  set(gca,'TickLabelInterpreter','latex');


%% setting up Eq 68 in LeM-S to find p_inf_liq
syms LeMS_Eq_68_all(p_inf_liq) LeMS_Eq_68_range_small(p_inf_liq) LeMS_Eq_68_range_medium(p_inf_liq);

LeMS_Eq_68_all = LeMS_Eq_68_const_term_all - reference_dens_sound_speed_sq .* B_term_all + reference_dens_sound_speed_sq .* A_term_relative_coeff_all .* A_term_all + reference_dens_sound_speed_sq .* A_term_all .* B_term_all + p_inf_liq;

%% finding p_inf_liq
p_inf_sol_all_sym = vpasolve(LeMS_Eq_68_all == 0, p_inf_liq);

%% convert p_inf_liq from symbol to number
p_inf_sol_all_num = double(p_inf_sol_all_sym);

%% take all of the values of p_inf_sol that are positive and real. There's always only one! LeM Eq 64: c_p,l - c_v,l
p_inf_liq_result_all = 24.521478323038426; % taken from p_inf_sol_all_num. 6.284993546892424 for trunc dec
LeM_Eq_64_T_P_frac_all = dodecane_data(:,1) ./ (dodecane_data(:,2) + p_inf_liq_result_all);
LeM_Eq_64_num_all = sum(LeM_Eq_64_T_P_frac_all .* (dodecane_data(:,4) - v_liq_avg_all));
LeM_Eq_64_denom_all = sum(LeM_Eq_64_T_P_frac_all .* (LeM_Eq_64_T_P_frac_all - mean(LeM_Eq_64_T_P_frac_all)));
c_p_liq_minus_c_v_liq_all = LeM_Eq_64_num_all ./ LeM_Eq_64_denom_all;
b_liq_all = v_liq_avg_all - (c_p_liq_minus_c_v_liq_all .* mean(LeM_Eq_64_T_P_frac_all));

%%
v_liq_analytic_all = zeros(data_size_full,1);
dens_liq_analytic_all = zeros(data_size_full,1);
dens_liq_vukalovich = zeros(data_size_full,1);

for i=1:data_size_full
  %v_liq_analytic_all(i) = ((gamma_liq_all - 1) .* c_p_liq_minus_c_v_liq_all .* dodecane_data(i,1)) ./ (dodecane_data(i,2) + p_inf_liq_result_all) + b_liq_all;
  %v_liq_analytic_all(i) = (((c_p_liq_all - c_v_liq_all) .* dodecane_data(i,1))./(dodecane_data(i,2) + p_inf_liq_result_all)) + b_liq_all;
  v_liq_analytic_all(i) = (((c_p_liq_minus_c_v_liq_all) .* dodecane_data(i,1))./(dodecane_data(i,2) + p_inf_liq_result_all)) + b_liq_all;
  dens_liq_analytic_all(i) = 1./v_liq_analytic_all(i);
  dens_liq_vukalovich(i) = 1./dodecane_data(i,4);
end

clear i;

%%
figure(6);
  plot(dodecane_data(:,1),dodecane_data(:,4),'.','Color',[0.4660 0.6740 0.1880],'MarkerSize',16);
  hold on;
  plot(dodecane_data(:,1),v_liq_analytic_all,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Volume vs. Temperature, Mercury Liquid, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${v}_{l} \left( \mathrm{{m}^{3} / kg} \right) $','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Grigull)','NASG liquid phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(7);
  semilogy(dodecane_data(:,1),dens_liq_vukalovich,'.','Color',[0.3010 0.7450 0.9330],'MarkerSize',16);
  hold on;
  semilogy(dodecane_data(:,1),dens_liq_analytic_all,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
  fontsize(30,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${\rho}_{l} \left( \mathrm{kg / {m}^{3}} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Grigull)','NASG liquid phase' + string(newline) + 'theoretical prediction','Location','northeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');
%%
c_p_liq_all = (sum(dodecane_data(:,1) .* (dodecane_data(:,7) - h_liq_avg_all)) - (b_liq_all .* sum(dodecane_data(:,1) .* (dodecane_data(:,2) - p_avg_all))))/sum(dodecane_data(:,1) .* (dodecane_data(:,1) - T_avg_all));
c_v_liq_all = c_p_liq_all - c_p_liq_minus_c_v_liq_all;
q_liq_all = h_liq_avg_all - c_p_liq_all .* T_avg_all - b_liq_all .* p_avg_all;
h_liq_analytic_all = zeros(data_size,1);

for i=1:data_size_full
  h_liq_analytic_all(i) = (c_p_liq_all .* dodecane_data(i,1)) + (b_liq_all .* dodecane_data(i,2)) + q_liq_all;
end

clear i;

%%
gamma_liq_all = c_p_liq_all/c_v_liq_all;

%%
figure(8);
  plot(dodecane_data(:,1),dodecane_data(:,7),'.','Color',[0.6350 0.0780 0.1840],'MarkerSize',16);
  hold on;
  plot(dodecane_data(:,1),h_liq_analytic_all,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Enthalpy vs. Temperature, Mercury Liquid, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${h}_{l} \left( \mathrm{kJ/kg} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Grigull)','NASG liquid phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

%% LeM Eq 41 coefficients, to determine specific entropies
LeM_Eq_41_B_coeff_all = (q_liq_all - q_gas_all)./(c_p_gas_all - c_v_gas_all);
LeM_Eq_41_C_coeff_all = (c_p_gas_all - c_p_liq_all)./(c_p_gas_all - c_v_gas_all);
LeM_Eq_41_D_coeff_all = c_p_liq_minus_c_v_liq_all./(c_p_gas_all - c_v_gas_all);
LeM_Eq_41_E_coeff_all = (b_liq_all - b_gas_all)./(c_p_gas_all - c_v_gas_all);
q_prime_gas_all = (c_p_gas_all - c_p_liq_all) + (c_p_gas_all - c_v_gas_all) .* mean(log(dodecane_data(:,2)) - (LeM_Eq_41_B_coeff_all + (LeM_Eq_41_E_coeff_all .* water_data_full(:,2)))./water_data_full(:,1) - LeM_Eq_41_C_coeff_all .* log(water_data_full(:,2)) - LeM_Eq_41_D_coeff_all .* log((water_data_full(:,2) + p_inf_liq_result_all)));
q_prime_liq_all = 0;

s_liq_analytic_all = zeros(data_size,1);
s_gas_analytic_all = zeros(data_size,1);
cs_sq_liq_analytic_all = zeros(data_size,1);

for i=1:data_size
  s_liq_analytic_all = c_v_liq_all .* log(((dodecane_data(:,1) .^ gamma_liq_all) ./ ((dodecane_data(:,2) + p_inf_liq_result_all) .^(gamma_liq_all - 1)))) + q_prime_liq_all;
  s_gas_analytic_all = c_v_gas_all .* log(((dodecane_data(:,1) .^ gamma_gas_all) ./ ((dodecane_data(:,2) + p_inf_gas_all) .^(gamma_gas_all - 1)))) + q_prime_gas_all;
  cs_sq_liq_analytic_all = (gamma_liq_all .* (dodecane_data(:,4) .^ 2) .* (dodecane_data(:,2) - p_inf_liq_result_all)) ./ (dodecane_data(:,4) - b_liq_all);
end

clear i;

%%
figure(4);
  plot(dodecane_data(:,1),dodecane_data(:,11),'.','Color',[0.4390 0.1650 0.5290],'MarkerSize',16);
  hold on;
  plot(dodecane_data(:,1),s_gas_analytic_all,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Entropy vs. Temperature, Mercury Gas, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${s}_{g} \left( \mathrm{kJ/ \left(kg \cdot K \right)} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Grigull)','NASG gas phase' + string(newline) + 'theoretical prediction','Location','northeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(9);
  plot(dodecane_data(:,1),dodecane_data(:,10),'.','Color',[0.4390 0.1650 0.5290],'MarkerSize',16);
  hold on;
  plot(dodecane_data(:,1),s_liq_analytic_all,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Entropy vs. Temperature, Mercury Gas, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${s}_{l} \left( \mathrm{kJ/ \left(kg \cdot K \right)} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Grigull)','NASG liquid phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

%%
v_liq_analytic_first_way = zeros(data_size,1);
v_liq_analytic_second_way = zeros(data_size,1);
v_liq_analytic_third_way = zeros(data_size,1);

for i=1:data_size_full
  v_liq_analytic_first_way(i) = ((gamma_liq_all - 1) .* c_p_liq_minus_c_v_liq_all .* dodecane_data(i,1)) ./ (dodecane_data(i,2) + p_inf_liq_result_all) + b_liq_all;
  v_liq_analytic_second_way(i) = (((c_p_liq_all - c_v_liq_all) .* dodecane_data(i,1))./(dodecane_data(i,2) + p_inf_liq_result_all)) + b_liq_all;
  v_liq_analytic_third_way(i) = (((c_p_liq_minus_c_v_liq_all) .* dodecane_data(i,1))./(dodecane_data(i,2) + p_inf_liq_result_all)) + b_liq_all;
end

%%
figure(10);
  semilogy(dodecane_data(:,1),dodecane_data(:,4),'.','Color',[0 0.4470 0.7410],'MarkerSize',16);
  hold on;
  semilogy(dodecane_data(:,1),v_liq_analytic_first_way,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
  hold on;
  semilogy(dodecane_data(:,1),v_liq_analytic_second_way,'-','Color',[0.4390 0.1650 0.5290],'LineWidth',1.5);
  hold on;
  semilogy(dodecane_data(:,1),v_liq_analytic_third_way,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
  fontsize(15,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${v}_{g} \left( \mathrm{{m}^{3} / kg} \right)$','Interpreter','latex');
  legend('Experimental values (Grigull)','NASG liquid version 1','NASG liquid version 2','NASG liquid version 3','Location','northeast','Interpreter','latex','FontSize',14);
  set(gca,'TickLabelInterpreter','latex');