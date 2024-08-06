vukalovich_data = table2array(readtable('./Vukalovich_Saturation_Line_by_Temperature_ASCII.csv'));

% reference state: 20 C point on saturation curve (using Vukalovich data). Density taken from Wikipedia.
reference_p = 1.729e-6;
reference_dens = 13545.83;
reference_sound_speed_sq = (1451.4).^2;
reference_dens_sound_speed_sq = reference_dens .* reference_sound_speed_sq;
R = 8.31446261815324./100000;

T_data = vukalovich_data(:,1);
p_data = vukalovich_data(:,2);
v_data_liq = vukalovich_data(:,4);
v_data_gas = vukalovich_data(:,5);
%%
liq_Tv_data = [T_data v_data_liq];
gas_Tv_data = [T_data v_data_gas];

mercury_liq_RK_consts = [0.000217617; 7.24549e-6];
mercury_gas_RK_consts = [0; 0];

RK_press_liq = RK_eq(mercury_liq_RK_consts,liq_Tv_data);
RK_press_gas = RK_eq(mercury_gas_RK_consts,gas_Tv_data);

%%
RK_liq_param_result_Wiki_params = nlinfit(liq_Tv_data,p_data,@RK_eq,mercury_liq_RK_consts)

%%
RK_gas_param_result_Wiki_params = nlinfit(gas_Tv_data,p_data,@RK_eq,mercury_gas_RK_consts)
