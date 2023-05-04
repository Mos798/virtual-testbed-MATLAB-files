t = out.simlog.Battery.Input_Current.i.series.time;
current = out.simlog.Battery.Input_Current.i.series.values('A');
T1 = out.simlog.Battery.Pack_1.T.series.values('degC');
T2 = out.simlog.Battery.Pack_2.T.series.values('degC');
T3 = out.simlog.Battery.Pack_3.T.series.values('degC');
T4 = out.simlog.Battery.Pack_4.T.series.values('degC');
T_ave = (T1+T2+T3+T4)/4;
charging_current = [t,current];
charging_temp = [t,T_ave];
T_ave_tot = trapz(t,T_ave)/t(end);
save('charge_profile','charging_temp',"charging_current")