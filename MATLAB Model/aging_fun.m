clc;clear all;close all

% Combined aging function
% cycle aging model from 'Wang, John, et al. "Cycle-life model for graphite-LiFePO4 cells." Journal of power sources 196.8 (2011): 3942-3948.'
% calendar aging model from 'Naumann, Maik, et al. "Analysis and modeling of calendar aging of a commercial LiFePO4/graphite cell." Journal of Energy Storage 17 (2018): 153-169.'

% either load the charge, discharge file into workspace or uncomment the
% next two lines to load the charge and discharge profiles each inluding current
% and temperature with time
%load('ch_c2_T45.mat')
%load('d_T=45.mat')


% forming the current and temperature profile for one charging + one discharge cycles
I = [charging_current;discharge_current(:,1)+charging_current(end,1),discharge_current(:,2)];
T = [charging_temp(:,1),charging_temp(:,2);discharge_temp(:,1)+charging_temp(end,1),discharge_temp(:,2)];

[loss_ratio,n,q_plt,q_cyc_n,q_cal_n] = cap_loss(I,T);

function [loss_ratio,n,q_plt,q_cyc_n,q_cal_n] = cap_loss(I,T)
% q_cyc: capacity loss due to cycling(as percentage of initial capacity)
% q_cal: capacity loss due to storage(as percentage of initial capacity)
% n: cycle number
% q_plt: capacity loss at the end of each cycle


% Coefficients for Nuaman Capacity loss
k_ref = 0.0012571;c=2.8575;d=0.60225;
soc = 0.5;
k_soc = c*(soc-0.5)^3+d;
E=17126;
R=8.314;
T_ref = 298.15;
T_cal =25;
z_c =0.5;
% time at rest for each cycle
time_cal = 365/88*24*3600-I(end,1);% seconds at rest for each one cycle
% cell assumed to be at rest any time other than charge discharge

% capacity of the cell has been derated to 2 Ah at Wang et al.
c_rate = abs(I(:,2))/2;

q_cyc = 0;
%q_cal_30 = 0;
%q_cal_65 = 0;
q_cal=0;
loss_ratio = 0;
n = 0; 
while (q_cyc+q_cal) < 20
    n = n + 1;
    for j=11:10:size(I,1)
        
        % interpolation function for preexponenetial
        p1=-47.84;p2=1215;p3=-9419;p4=3.604e4;
        B = p1*c_rate(j)^3+p2*c_rate(j)^2+p3*c_rate(j)+p4;
        z=0.55;
        dt = I(j,1)-I(j-10,1);
        if (c_rate(j)+c_rate(j-10))/2 > 0.5
            I_ave = (abs(I(j,2))+abs(I(j-10,2)))/2; % avergae of initial and final I in each delta t
            T_ave = (T(j,2)+T(j-10,2))/2+273.15; % avergae of initial and final T in each delta t
            disp(T_ave);
            X = I_ave*(B*exp((-31700+370.3*(c_rate(j)+c_rate(j-10))/2)/(8.314*(T_ave))))^(1/z);
           
            q_cyc = (X*dt/3600+q_cyc^(1/z))^z;
            
        end
        
    end
    % calendar loss for SoC = 0.30
    X_p = (k_ref*k_soc * exp(-E/R*(1/(T_cal+273.15)-1/T_ref)))^(1/z_c);
    q_cal = (X_p*time_cal+q_cal^(1/z_c))^z_c;
    
    q_plt(n)=q_cyc+q_cal;
    q_cyc_n(n) = q_cyc;
    q_cal_n(n) = q_cal;
    loss_ratio = q_cyc/20;
    
end
end


