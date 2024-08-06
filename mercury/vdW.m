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

mercury_liq_vdW_consts = [8.2; 0.01696];
mercury_gas_vdW_consts = [0; 0];

vdW_press_liq = vdW_eq(mercury_liq_vdW_consts,liq_Tv_data);
vdW_press_gas = vdW_eq(mercury_gas_vdW_consts,gas_Tv_data);

%%
vdW_liq_param_result_Wiki_params = nlinfit(liq_Tv_data,p_data,@vdW_eq,mercury_liq_vdW_consts)
vdW_liq_param_result_zero_params = nlinfit(liq_Tv_data,p_data,@vdW_eq,[0; 0])
vdW_liq_param_result_params_1 = nlinfit(liq_Tv_data,p_data,@vdW_eq,[0.001; 0.001])
vdW_liq_param_result_params_2 = nlinfit(liq_Tv_data,p_data,@vdW_eq,[0.00001; 1])
vdW_liq_param_result_params_3 = nlinfit(liq_Tv_data,p_data,@vdW_eq,[0.000000001; 1])
vdW_liq_param_result_params_4 = nlinfit(liq_Tv_data,p_data,@vdW_eq,[0; 1])
vdW_liq_param_result_params_5 = nlinfit(liq_Tv_data,p_data,@vdW_eq,[0.000001; 0.0000005])
%%
vdW_gas_param_result_Wiki_params = nlinfit(gas_Tv_data,p_data,@vdW_eq,mercury_gas_vdW_consts)
vdW_gas_param_result_zero_params = nlinfit(gas_Tv_data,p_data,@vdW_eq,[0; 0])
vdW_gas_param_result_params_1 = nlinfit(gas_Tv_data,p_data,@vdW_eq,[0.001; 0.001])
vdW_gas_param_result_params_2 = nlinfit(gas_Tv_data,p_data,@vdW_eq,[0.00001; 1])
vdW_gas_param_result_params_3 = nlinfit(gas_Tv_data,p_data,@vdW_eq,[0.000000001; 1])
vdW_gas_param_result_params_4 = nlinfit(gas_Tv_data,p_data,@vdW_eq,[0; 1])
vdW_gas_param_result_params_5 = nlinfit(gas_Tv_data,p_data,@vdW_eq,[0.000001; 0.0000005])
vdW_gas_param_result_params_6 = nlinfit(gas_Tv_data,p_data,@vdW_eq,[1; 1])
