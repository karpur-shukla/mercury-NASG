%original
Vukalovich_data = table2array(readtable('./Vukalovich_Saturation_Line_by_Temperature_ASCII.csv'));
Fokin_data_original = table2array(readtable('./Fokin_Saturation_Line_by_Temperature_ASCII.csv'));
NIST_primary_data_original = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/_NIST_primary_all.csv'));
NIST_primary_KS_additional_data_original = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/_NIST_primary_KS_additional_all.csv'));

%primary data, individual (if necessary)
% NIST_primary_Ambrose_Sprake_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/Ambrose-Sprake_(1972).csv'));
% NIST_primary_Beattie_Blaisdell_Kaminsky_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/Beattie-Blaisdell-Kaminsky_(1937).csv'));
% NIST_primary_Ernsberger_Pitman_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/Ernsberger-Pitman_(1955).csv'));
% NIST_primary_Menzies_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/Menzies_(1920-1927).csv'));
% NIST_primary_Schonherr_Hensel_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/Schönherr-Hensel_(1981).csv'));
% NIST_primary_Shpilrain_Nikanorov_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/Shpilrain-Nikanorov_(1971).csv'));
% NIST_primary_Spedding_Dye_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets/Spedding-Dye_(1955).csv'));

%secondary main data, individual (if necessary)
% NIST_primary_Douglas_Ball_Ginnings_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Douglas-Ball-Ginnings_(1951).csv'));
% NIST_primary_Millar_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Millar_(1927).csv'));
% NIST_primary_Murgulescu_Topor_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Murgulescu-Topor_(1966).csv'));
% NIST_primary_Pedder_Barratt_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Pedder-Barratt_(1933).csv'));
% NIST_primary_Rodebush_Dixon_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Rodebush-Dixon_(1925).csv'));
% NIST_primary_Roeder_Morawietz_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Roeder-Morawietz_(1956).csv'));
% NIST_primary_Schmahl_Barthel_Kaloff_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Schmahl-Barthel-Kaloff_(1965).csv'));
% NIST_primary_Scott_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Scott_(1924).csv'));
% NIST_primary_Sugawara_Sato_Minamiyama_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Sugawara-Sato-Minamiyama_(1962).csv'));
% NIST_primary_Volmer_Kirchhoff_data = table2array(readtable('./DoE_refs/NIST_T_vs_p_CSVs/primary_data_sets_KS_additional/Volmer-Kirchhoff_(1925).csv'));
%%
figure(1);
  semilogy(vukalovich_data(:,1),vukalovich_data(:,2),'.','Color',[0 0.4470 0.7410],'MarkerSize',12);
  hold on;
  semilogy(vukalovich_data(:,1),v_gas_analytic,'.','Color',[0.8500 0.3250 0.0980],'LineWidth',12);
  fontsize(20,"points");
  %title('\bf{Specific Volume vs. Temperature, Mercury Gas, NA}','Interpreter','latex','FontSize', 40);
  xlabel('$T \left( \mathrm{K} \right)$','Interpreter','latex');
  ylabel('${v}_{g} \left( \mathrm{{m}^{3} / kg} \right)$','Interpreter','latex');
  legend('Experimental values' + string(newline) + '(Vukalovich)','NA gas phase' + string(newline) + 'theoretical prediction','Location','northeast','Interpreter','latex','FontSize',20);
  set(gca,'TickLabelInterpreter','latex');