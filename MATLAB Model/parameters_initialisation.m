clear;
clc;

%Add Subfolders to Path
addpath(genpath(cd))

%Parameters Initialisation
% load the autonomie drive cycle data (we only use speed and power, speed
% goes to condenser and radiator, power is used to calculate current from
% voltage of the battery pack at that particular temperature
load autonomie_drive_cycles_data.mat;

% load the the 2 RC cell model from 
% Lin, Xinfan, et al. "A lumped-parameter electro-thermal model for cylindrical batteries." Journal of Power Sources 257 (2014): 1-11.
load 2RC_battery_parameters.mat;

% simulation time is set to 10000 sec but the simulation stops at SOC=0.8
% during charge and SOC=0.2 during discharge
sim_t_end=10000;


%% Battery 
n_parallel=84; % number of cells in parralel
Ncells=120/4;  % dividing by 4 as there are four packs in simulink model
cell_capacity=2.5; % Ah
cell_mass=0.076; %kg
cell_Cp_heat=1130; % H. Maleki et al, “Thermal Properties of Lithium-Ion Battery and Components”, Journal of The Electrochemical Society, 146 (3) 947-954 (1999)
c_rate = 2;
const_current=c_rate*cell_capacity*n_parallel;
charging='Y'; %Y for charging, N for Discharging
if charging=='Y'
    Qe_init=.8*cell_capacity; 
    C12_SOC=C12_SOC_c;
    R12_SOC=R12_SOC_c;
    C22_SOC=C22_SOC_c;
    R22_SOC=R22_SOC_c;
    Rs2_SOC=Rs2_SOC_c;
    SOC_Temp=SOC_Temp_c;
else
    Qe_init=.2*cell_capacity;
    C12_SOC=C12_SOC_d;
    R12_SOC=R12_SOC_d;
    C22_SOC=C22_SOC_d;
    R22_SOC=R22_SOC_d;
    Rs2_SOC=Rs2_SOC_d;
    SOC_Temp=SOC_Temp_d;
end
clear R12_SOC_c R12_SOC_d C12_SOC_c C12_SOC_d R22_SOC_c R22_SOC_d C22_SOC_c C22_SOC_d Rs2_SOC_c Rs2_SOC_d SOC_Temp_c SOC_Temp_d

%% Operating and Critical Temperatures
T_battery_nom=30+273;  %K  set point teperature of battery

%% Env conditions
summer_T=298; %K
extreme_summer_T=318; %K
winter_T=273; %K
extreme_winter_T= 253; %K

ambient_temp=273+45;% change this temperature to ambient temperature
ambient_pressure=0.101325; %MPa
ambient_rh=.65;  %annual average
ambient_co2=412*10^-6; 

T_env=ambient_temp;
%T_init=T_env-0;
%T_init_bat=303.15; % the initial battery temperature, equal to precooled temperature
% if precooled, otherwise equal to ambient temperature
T_init_bat = T_env;

%% Coolant and Tank Properties
coolant_initial_T=ambient_temp;
coolant_initial_P=ambient_pressure;
coolant_ethylene_glycol_vol_frac = 0.5; 
coolant_tank_vol=2; %litres
coolant_tank_initial_liq_vol=0.3*coolant_tank_vol; 
coolant_tank_cross_sec_area=.02; % m^2
coolant_pipe_cross_sec_area= 0.0000785; % m^2

%% Centrifugal Pump 
%https://www.vovyopump.com/product/12-volt-electric-water-pump-automotive/#Curve
pump_nom_capacity=25; %lpm
pump_max_capacity=50; %lpm
pump_nom_head=7.5; %m
pump_max_head=10; %m
pump_nom_brake_power=60; %W
pump_ref_speed=2000; %RPM

%% HX Coolant Pipe Channel 
coolant_channel_segments = 1;
coolant_channel_length = .5; %metre
coolant_channel_cs_area = coolant_pipe_cross_sec_area;
coolant_channel_hydraulic_dia = 10^-2; %m
%Haaland correlation parameters
coolant_channel_eqv_resistance_length=2*coolant_channel_length; %m
coolant_channel_int_surf_roughness=1.5^10^-5; %m
coolant_channel_lam_friction_const=64;
coolant_channel_lam_flow_Re=100;
coolant_channel_turb_flow_Re=500;
%Dittus-Boelter correlation-Nusselt=a*Re^b*Pr*c
coolant_channel_ht_coeff_a=0.3;
coolant_channel_ht_coeff_b=0.8;
coolant_channel_ht_coeff_c=0.33;
coolant_channel_ht_Nusselt=3.66;

%% Radiator 
al_conductivity = 239;
radiator_total_width=0.5; % m
radiator_tube_length=0.9; % m
radiator_tube_width=.022; % m
radiator_tube_height=.0015; % m
radiator_tube_pitch= 0.009; %.01107-.003048; % m
radiator_n_tubes=round(radiator_total_width/radiator_tube_pitch); %width/tube pitch
radiator_tube_eqv_resistance_length=0.2*(radiator_n_tubes*radiator_tube_length); % m
radiator_air_cross_sec_area=radiator_total_width*radiator_tube_length; % m^2
radiator_wall_thermal_conductivity= 239; % W/m/K  Aluminium thermal conductivity
radiator_wall_thickness=0.0001; %m

radiator_air_min_flow_area=radiator_air_cross_sec_area-radiator_tube_height*radiator_tube_length*radiator_n_tubes-.00015*(radiator_tube_pitch-radiator_tube_height)*582*radiator_tube_length*radiator_n_tubes;  
    % radiator c/s area - tubes c/s area - fins c/s area  m^2
radiator_air_ht_area=radiator_tube_length*(radiator_tube_width+radiator_tube_height)*2*(radiator_n_tubes-1);
radiator_ht_surf_area_without_fins= radiator_tube_width*radiator_tube_length*2*radiator_n_tubes;
radiator_fin_surf_area=(radiator_tube_pitch-radiator_tube_height)*1.1*radiator_tube_width*2*582*radiator_tube_length*radiator_n_tubes;  
    % fin height*fin height factor for extra length*fin width*two surfaces*fin pitch*radiator length* number of tubes m^2
radiator_air_fin_eff=0.7;
radiator_wall_thermal_resistance=radiator_wall_thickness/radiator_wall_thermal_conductivity/radiator_air_ht_area;
%Correlation input
rad_air_side_coef=[0.27, .67, 0.36]; % Coefficients [a, b, c] for a*Re^b*Pr^c

%% R134a Referigernat and Tank Properties
refrigerant_initial_P=1.5; %MPa
refrigerant_initial_P_cond = 1.5; %MPa
refrigerant_initial_P_chiller = 0.3; %MPa
refrigerant_pipe_cross_sec_area= 0.0000785; % m^2
refrigerant_tank_vol=3; %litres make it around 2
refrigerant_initial_void_fraction = 0.6; %vapor void fraction

%refrigerant_initial_T=ambient_temp;
%refrigerant_tank_cross_sec_area=.02; % m^2


%% Condenser 
condenser_height=0.4466; % m
condenser_width=0.611; % m
condenser_depth=2.2*.0185; % m

condenser_tube_height=.0015; % m
condenser_tube_pitch= 0.009; %.01107-.003048; % m

condenser_coolant_pressure_loss_coef=1750;  
condenser_air_pressure_loss_coef=5;  

%Calcuated parameters
condenser_tube_length=condenser_width; % m
condenser_tube_width=condenser_depth; % m
condenser_n_tubes=round(condenser_height/condenser_tube_pitch); %width/tube pitch
condenser_tube_eqv_resistance_length=0.18*(condenser_n_tubes*condenser_tube_length); % m
condenser_air_cross_sec_area=condenser_height*condenser_width; % m^2
condenser_wall_thermal_conductivity= al_conductivity; % W/m/K  Aluminium thermal conductivity
condenser_wall_thickness=0.0001; %m
%rad_coolant_side_coef=[0.023, .8, 0.3]; % Coefficients [a, b, c] for a*Re^b*Pr^c

%Air Side Generic Parameterization
condenser_air_min_flow_area=condenser_air_cross_sec_area-condenser_tube_height*condenser_tube_length*condenser_n_tubes-.00015*(condenser_tube_pitch-condenser_tube_height)*582*condenser_tube_length*condenser_n_tubes;  
    % condenser c/s area - tubes c/s area - fins c/s area  m^2
condenser_air_ht_area=condenser_tube_length*(condenser_tube_width+condenser_tube_height)*2*(condenser_n_tubes-1);
condenser_ht_surf_area_without_fins= condenser_tube_width*condenser_tube_length*2*condenser_n_tubes;
condenser_fin_surf_area=(condenser_tube_pitch-condenser_tube_height)*1.2*condenser_tube_width*2*(582)*condenser_tube_length*condenser_n_tubes;  %582
    % fin height*fin height factor for extra length*fin width*two surfaces*fin pitch*condenser length* number of tubes m^2
condenser_air_fin_eff=1; %0.7;
condenser_wall_thermal_resistance=0;%condenser_wall_thickness/condenser_wall_thermal_conductivity/condenser_air_ht_area;

%Correlation input
condenser_air_side_coef=[0.27, .67, 0.36];%[4*0.27, .67, 0.36]; % Coefficients [a, b, c] for a*Re^b*Pr^c, Can be looked over from Kays and London % Fig 10.30

%% Evaporator 
evaporator_n_tubes=20;
evaporator_tube_length=0.23; % m
evaporator_total_width=0.18; % m
evaporator_tube_width=.1; % m
evaporator_tube_height=.0015; % m
evaporator_tube_pitch=evaporator_total_width/evaporator_n_tubes; % m
evaporator_tube_eqv_resistance_length=.2*(evaporator_tube_length*evaporator_n_tubes); % m
evaporator_air_cross_sec_area=evaporator_total_width*evaporator_tube_length; % m^2
evaporator_wall_thermal_conductivity= 239; % W/m/K; Aluminium thermal conductivity
evaporator_wall_thickness=0.0001; %m

evaporator_air_min_flow_area=evaporator_air_cross_sec_area-evaporator_tube_height*evaporator_tube_length*evaporator_n_tubes-.00015*(evaporator_tube_pitch-evaporator_tube_height)*582*evaporator_tube_length*evaporator_n_tubes;  % m^2
evaporator_ht_surf_area_without_fins= evaporator_tube_width*evaporator_tube_length*2*evaporator_n_tubes;
evaporator_fin_surf_area=(evaporator_tube_pitch-evaporator_tube_height)*1.1*evaporator_tube_width*2*582*evaporator_tube_length*evaporator_n_tubes;  
    %fin height*fin height factor for extra length*fin width*two surfaces*fin pitch*evaporator length* number of tubes m^2
evaporator_air_fin_eff=0.7;
evaporator_wall_thermal_resistance=evaporator_wall_thickness/evaporator_wall_thermal_conductivity/evaporator_ht_surf_area_without_fins;
%evaporator_air_flow_area=.7;  % m^2

%% Chiller  (Tube and Shell HX)
chiller_tube_length=.5;  %length along tube flow (Refrigerant) *4
chiller_width=.1;   %length transverse to shell fluid flow (TL)
chiller_height=chiller_width;  %length along shell fluid flow (TL)
chiller_tubes_n_rows= 10; %along shell flow direction (TL)
chiller_tubes_n_columns= chiller_tubes_n_rows; %transverse to shell flow direction (TL)
chiller_tubes_n_total=chiller_tubes_n_rows*chiller_tubes_n_columns;
chiller_baffles_n=4;
chiller_long_pitch=chiller_height/chiller_tubes_n_rows;
chiller_transverse_pitch=chiller_width/chiller_tubes_n_columns;
chiller_baffles_area= chiller_width*(chiller_height*0.8)*2*chiller_baffles_n;
chiller_bafffles_fin_eff=0.7;
chiller_tube_dia=0.004; %m
chiller_tube_eqv_resistance_length=0;%chiller_tube_length*.25; 

chiller_tube_width=chiller_tube_dia;
chiller_wall_thermal_conductivity= al_conductivity;
chiller_wall_thickness=0.0001; %m
chiller_ht_surf_area_without_fins= chiller_tube_width*chiller_tube_length*2*chiller_tubes_n_total;
chiller_wall_thermal_resistance=chiller_wall_thickness/chiller_wall_thermal_conductivity/chiller_ht_surf_area_without_fins;

%% Expansion valve
txv_nom_capacity=5000; %W
txv_max_capacity=7500; %W
evap_nom_temp=273; %K
evap_min_superheat_delta=5; % delta K
evap_nom_superheat_delta=10; %delta K
cond_nom_temp=313; %K
cond_nom_subcool_delta=5; %K

%% Compressor 
%load('compressor_parameters.mat')  %taken from Matlab tutorial 
%load('matlab_model_data.mat')
compressor_fixed_m_dot=0.018; %kg/s fixed refrigerant mass flow rate
eta_mech_elec = 0.85; % mechanical and electrical efficiency of compressor
eta_isen = 0.7; % isentropic efficiencyog compressor
%% Vavles
three_w_valve_max_opening=.005; %m
three_w_valve_long_len=.05; %m
three_w_valve_max_open_area_factor=0.98;
%% condenser and radiator fan volumetric flow rate
fan_vdot = 3000*0.00047194745; % [m3/s]
%% compressor mass flow rate PID control
%PID_prop = 1/10; % set the proportional part of PID to 1/10 in case of fast-charging 
PID_prop = 1/30; % set to 1/30 for discharge and precooling
