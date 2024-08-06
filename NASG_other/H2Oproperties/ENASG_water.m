water_data_trunc = table2array(readtable('./../H2Oproperties/Grigull_CSV_Vukalovich_order_ASCII_truncated_300_500K.csv'));
water_data_trunc_dec = table2array(readtable('./../H2Oproperties/Grigull_CSV_Vukalovich_order_ASCII_decimated_10_truncated_300_500K.csv'));
water_data_full = table2array(readtable('./../H2Oproperties/Grigull_CSV_Vukalovich_order_ASCII'));
%save('./Hgproperties/Vukalovich_Saturation_Line_by_Temperature_ASCII.mat',"water_data_all");

%%
