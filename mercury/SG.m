vukalovich_data = table2array(readtable('./Vukalovich_Saturation_Line_by_Temperature_ASCII.csv'));
%save('./Hgproperties/Vukalovich_Saturation_Line_by_Temperature_ASCII.mat',"vukalovich_data_all");

%%
% Defining the functions used for the liquid phase parameters
syms LeMS_Eq_64_num(p_inf_liq) LeMS_Eq_64_denom(p_inf_liq) LeMS_Eq_64_T_over_P_avg(p_inf_liq) LeMS_Eq_68_all(p_inf_liq);

% parameters from ideal gas approx (Pg. 9 of LeMétayer-Saurel, below Eq. 45 and 46)
b_gas_all = 0;
b_liq_all = 0;
p_inf_gas_all = 0;

% Defining the average temperature, pressure, average specific enthalpy (gas and liquid) for both the full and the ranges
T_avg_all = mean(vukalovich_data(:,1));
p_avg_all = mean(vukalovich_data(:,2));
h_gas_avg_all = mean(vukalovich_data(:,8));
h_liq_avg_all = mean(vukalovich_data(:,7));
v_liq_avg_all = mean(vukalovich_data(:,4));

%c_3 constants
T_avg_norm_factor_all = sum(vukalovich_data(:,1) .* (vukalovich_data(:,1) - T_avg_all));

% Using Eq. 50, 51, 55 in LeMétayer-Saurel to calculate c_p,g; q_g; c_v,g; gamma_g for both full data set and range. c_p,g in J/(g K); c_v,g in m^3/kg.
c_p_gas_all = sum(vukalovich_data(:,1) .* (vukalovich_data(:,8) - h_gas_avg_all))/T_avg_norm_factor_all;
q_gas_all = h_gas_avg_all - c_p_gas_all .* T_avg_all;
c_v_gas_all = c_p_gas_all - sum(vukalovich_data(:,1) .* vukalovich_data(:,5) ./ vukalovich_data(:,2))/sum((vukalovich_data(:,1) ./ vukalovich_data(:,2)) .^ 2);
gamma_gas_all = c_p_gas_all/c_v_gas_all;

% reference state: 20 C point on saturation curve (using Vukalovich data). Density taken from Wikipedia.
reference_p = 1.729e-6;
reference_dens = 13545.83;
reference_sound_speed_sq = (1451.4).^2;
reference_dens_sound_speed_sq = reference_dens .* reference_sound_speed_sq;
ideal_gas_const = 8.31446261815324./100000;

%
v_gas_analytic = zeros(88,1);
h_gas_analytic = zeros(88,1);
compress_factor_gas_analytic = zeros(88,1);
dens_gas_analytic = zeros(88,1);
dens_gas_vukalovich = zeros(88,1);

for i=1:88
  v_gas_analytic(i) = ((gamma_gas_all - 1) * c_v_gas_all * vukalovich_data(i,1))/(vukalovich_data(i,2));
  compress_factor_gas_analytic(i) = (vukalovich_data(i,2) .* v_gas_analytic(i)) ./ (ideal_gas_const .* vukalovich_data(i,1) .* 4.9856360761414);
  h_gas_analytic(i) = gamma_gas_all * c_v_gas_all * vukalovich_data(i,1) + q_gas_all;
  dens_gas_analytic(i) = 1 ./ v_gas_analytic(i);
  dens_gas_vukalovich(i) = 1 ./ vukalovich_data(i,5);
end

clear i;


%%
figure(1);
  semilogy(vukalovich_data(:,1),vukalovich_data(:,5),'.','Color',[0 0.4470 0.7410],'MarkerSize',16);
  hold on;
  semilogy(vukalovich_data(:,1),v_gas_analytic,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
  fontsize(30,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${v}_{g} \left( \mathrm{{m}^{3} / kg} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG gas phase' + string(newline) + 'theoretical prediction','Location','northeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(2);
  semilogy(vukalovich_data(:,1),dens_gas_vukalovich,'.','Color',[0.3010 0.7450 0.9330],'MarkerSize',16);
  hold on;
  semilogy(vukalovich_data(:,1),dens_gas_analytic,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
  fontsize(30,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${\rho}_{g} \left( \mathrm{kg / {m}^{3}} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG gas phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(3);
  plot(vukalovich_data(:,1),vukalovich_data(:,8),'.','Color',[0.4940 0.1840 0.5560],'MarkerSize',16);
  hold on;
  plot(vukalovich_data(:,1),h_gas_analytic,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(30,"points");
  %title('\bf{Specific Enthalpy vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${h}_{g} \left( \mathrm{kJ/kg} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG gas phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(5);
  plot(vukalovich_data(:,1),vukalovich_data(:,6),'.','Color',[0.7 0 1],'MarkerSize',16);
  hold on;
  plot(vukalovich_data(:,1),compress_factor_gas_analytic,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Compressibility Factor vs. Temperature, Mercury Gas, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${z}_{g}$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG gas phase' + string(newline) + 'theoretical prediction','Location','southwest','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

%%
c_p_liq_all = sum(vukalovich_data(:,1) .* (vukalovich_data(:,7) - h_liq_avg_all)) ./ T_avg_norm_factor_all;
q_liq_all = h_liq_avg_all - c_p_liq_all .* T_avg_all;

%% setting up Eq 68 in LeM-S to find p_inf_liq
LeMS_Eq_68_const_term = vpa(reference_p - reference_dens_sound_speed_sq);
LeMS_Eq_64_num = sum(vpa(vukalovich_data(:,1) .* (vukalovich_data(:,4) - v_liq_avg_all))./vpa(vukalovich_data(:,2) - p_inf_liq));
LeMS_Eq_64_T_over_P_avg = (1/88) .* sum(vpa((vukalovich_data(:,1) ./ (vukalovich_data(:,2) + p_inf_liq))));
LeMS_Eq_64_denom = sum(vpa((vukalovich_data(:,1)./(vukalovich_data(:,2) - p_inf_liq))) .* vpa((vukalovich_data(:,1)./(vukalovich_data(:,2) - p_inf_liq)) - LeMS_Eq_64_T_over_P_avg));
LeMS_Eq_68_all = LeMS_Eq_68_const_term + p_inf_liq + (LeMS_Eq_64_num ./ (c_p_liq_all .* LeMS_Eq_64_denom));

%% finding p_inf_liq
p_inf_sol_all_sym = vpasolve(LeMS_Eq_68_all == 0, p_inf_liq);

%% convert p_inf_liq from symbol to number
p_inf_sol_all_num = double(p_inf_sol_all_sym);

%%
for i=1:263
  if imag(p_inf_sol_all_num(i)) == 0 && real(p_inf_sol_all_num(i)) > 0
    fprintf('%d\n',p_inf_sol_all_num(i));
  end
end

clear i;

%%
p_inf_liq_result_all = 2.853513e+10;
c_p_liq_minus_c_v_liq_all_num = sum((vukalovich_data(:,1) .* (vukalovich_data(:,4) - v_liq_avg_all)) ./ (vukalovich_data(:,2) - p_inf_liq_result_all));
T_over_P_avg = (1/88) .* sum((vukalovich_data(:,1) ./ (vukalovich_data(:,2) + p_inf_liq_result_all)));
c_p_liq_minus_c_v_liq_all_denom = sum((vukalovich_data(:,1) ./ (vukalovich_data(:,2) + p_inf_liq_result_all)) .* ((vukalovich_data(:,1) ./ (vukalovich_data(:,2) + p_inf_liq_result_all)) - T_over_P_avg));
c_p_liq_minus_c_v_liq_all = c_p_liq_minus_c_v_liq_all_num ./ c_p_liq_minus_c_v_liq_all_denom;

c_v_liq_all = c_p_liq_all - c_p_liq_minus_c_v_liq_all;

gamma_liq_all = c_p_liq_all/c_v_liq_all;

%%
v_liq_analytic_all = zeros(1,88);
dens_liq_analytic_all = zeros(1,88);
dens_liq_vukalovich = zeros(1,88);
h_liq_analytic_all = zeros(1,88);
s_liq_analytic_all = zeros(88,1);
s_gas_analytic_all = zeros(88,1);
cs_sq_liq_analytic_all = zeros(88,1);

for i=1:88
  v_liq_analytic_all(i) = (c_p_liq_minus_c_v_liq_all .* vukalovich_data(i,1)) ./ (vukalovich_data(i,2) + p_inf_liq_result_all);
  %v_liq_analytic_all(i) = ((gamma_liq_all - 1) .* c_v_liq_all .* vukalovich_data(i,1)) ./ (vukalovich_data(i,2) + p_inf_liq_result_all);
  h_liq_analytic_all(i) = c_p_liq_all .* vukalovich_data(i,1) + q_liq_all;
  dens_liq_analytic_all(i) = 1 ./ v_liq_analytic_all(i);
  dens_liq_vukalovich(i) = 1./vukalovich_data(i,4);
end

clear i;

%%
figure(6);
  plot(vukalovich_data(:,1),vukalovich_data(:,4),'.','Color',[0.4660 0.6740 0.1880],'MarkerSize',16);
  hold on;
  plot(vukalovich_data(:,1),v_liq_analytic_all,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Volume vs. Temperature, Mercury Liquid, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${v}_{l} \left( \mathrm{{m}^{3} / kg} \right) $','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG liquid phase' + string(newline) + 'theoretical prediction','Location','northwest','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(7);
  plot(vukalovich_data(:,1),dens_liq_vukalovich,'.','Color',[0.3010 0.7450 0.9330],'MarkerSize',16);
  hold on;
  plot(vukalovich_data(:,1),dens_liq_analytic_all,'-','Color',[0.8500 0.3250 0.0980],'LineWidth',1.5);
  fontsize(30,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${\rho}_{l} \left( \mathrm{kg / {m}^{3}} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG liquid phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(8);
  plot(vukalovich_data(:,1),vukalovich_data(:,7),'.','Color',[0.6350 0.0780 0.1840],'MarkerSize',16);
  hold on;
  plot(vukalovich_data(:,1),h_liq_analytic_all,'-','Color',[0.9290 0.6940 0.1250],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Enthalpy vs. Temperature, Mercury Liquid, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${h}_{l} \left( \mathrm{kJ/kg} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG liquid phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

%%
LeM_Eq_41_B_coeff_all = (q_liq_all - q_gas_all)./(c_p_gas_all - c_v_gas_all);
LeM_Eq_41_C_coeff_all = (c_p_gas_all - c_p_liq_all)./(c_p_gas_all - c_v_gas_all);
LeM_Eq_41_D_coeff_all = c_p_liq_minus_c_v_liq_all./(c_p_gas_all - c_v_gas_all);
q_prime_gas_all = (c_p_gas_all - c_p_liq_all) + (c_p_gas_all - c_v_gas_all) .* mean(log(vukalovich_data(:,2)) - (LeM_Eq_41_B_coeff_all)./vukalovich_data(:,1) - LeM_Eq_41_C_coeff_all .* log(vukalovich_data(:,2)) - LeM_Eq_41_D_coeff_all .* log((vukalovich_data(:,2) + p_inf_liq_result_all)));
q_prime_liq_all = 0;

for i=1:88
  s_liq_analytic_all = c_v_liq_all .* log(((vukalovich_data(:,1) .^ gamma_liq_all) ./ ((vukalovich_data(:,2) + p_inf_liq_result_all) .^(gamma_liq_all - 1)))) + q_prime_liq_all;
  s_gas_analytic_all = c_v_gas_all .* log(((vukalovich_data(:,1) .^ gamma_gas_all) ./ ((vukalovich_data(:,2) + p_inf_gas_all) .^(gamma_gas_all - 1)))) + q_prime_gas_all;
  cs_sq_liq_analytic_all = (gamma_liq_all .* (vukalovich_data(:,4) .^ 2) .* (vukalovich_data(:,2) - p_inf_liq_result_all)) ./ (vukalovich_data(:,4) - b_liq_all);
end

clear i;

%%
figure(4);
  semilogy(vukalovich_data(:,1),vukalovich_data(:,11),'.','Color',[0.4390 0.1650 0.5290],'MarkerSize',16);
  hold on;
  semilogy(vukalovich_data(:,1),s_gas_analytic_all,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Entropy vs. Temperature, Mercury Gas, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${s}_{g} \left( \mathrm{kJ/ \left(kg \cdot K \right)} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG gas phase' + string(newline) + 'theoretical prediction','Location','northeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');

figure(9);
  semilogy(vukalovich_data(:,1),vukalovich_data(:,10),'.','Color',[0.4390 0.1650 0.5290],'MarkerSize',16);
  hold on;
  semilogy(vukalovich_data(:,1),s_liq_analytic_all,'-','Color',[0.6350 0.0780 0.1840],'LineWidth',1.5);
  fontsize(30,"points");
  %title('Specific Entropy vs. Temperature, Mercury Gas, NA');
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${s}_{l} \left( \mathrm{kJ/ \left(kg \cdot K \right)} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','SG liquid phase' + string(newline) + 'theoretical prediction','Location','southeast','Interpreter','latex','FontSize',28);
  set(gca,'TickLabelInterpreter','latex');