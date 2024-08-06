function vdW_press = vdW_eq(vdW_param_in,Tv_data_in)
  R = 8.31446261815324./100000;

  a_in = vdW_param_in(1);
  b_in = vdW_param_in(2);

  T_in = Tv_data_in(:,1);
  v_in = Tv_data_in(:,2);

  vdW_press = (R .* T_in) ./ (v_in - b_in) - (a_in ./ (v_in.^2));
end