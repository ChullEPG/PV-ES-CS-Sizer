import numpy as np
########### NPV ########### 


def calculate_npv(initial_investment, cash_flows, discount_rate):
    present_values = [cf / (1 + discount_rate) ** idx for idx, cf in enumerate(cash_flows)]
    npv = sum(present_values) - initial_investment
    return npv



def calculate_pv_capital_cost(pv_capacity, a):
    inverter_cost = pv_capacity * a['inverter_cost_per_kw'] 
        
    panel_cost = a['pv_cost_per_kw'] * pv_capacity 

    # Components cost
    component_cost = panel_cost + inverter_cost # $/kW

    # Other costs
    peripherals_cost = component_cost * 0.20 # $/kW
    installation_cost = (component_cost + peripherals_cost) * 0.10 # $/kW
    markup_cost = (component_cost + peripherals_cost + installation_cost) * 0.33 # $/kW
    
    total_pv_capital_cost = component_cost + peripherals_cost + installation_cost + markup_cost
    
    return total_pv_capital_cost

def calculate_crf(a):
    # Calculate annual payments   
    return (a['interest rate'] * (1 + a['interest rate'])**a['Rproj']) / ((1 + a['interest rate'])**a['Rproj'] - 1)


def get_energy_savings(cost_of_trickle_charging, load_profile, net_load_profile, year, a):
    
    #load_profile = a['load_profile']
    #net_load_profile = load_profile - pv_with_battery_output_profile

    peak_hours = a['time_periods']['peak_hours']
    standard_hours = a['time_periods']['standard_hours']
    off_peak_hours = a['time_periods']['off_peak_hours']
    
    # Initialize total cost variables
    total_cost_no_pv = np.zeros(len(load_profile))
    total_cost_with_pv = np.zeros(len(net_load_profile))
    
    # Calculate total cost of energy with and without PV
    for i in range(len(total_cost_no_pv)):      
        curr_hour_of_week = i % 168 
        curr_hour_of_day = i % 24
        
        if (i > a['high_period_start']) & (i <= a['high_period_end']): # high period (all peak)
            peak_cost = a['time_of_use_tariffs_high']['peak'] * (1 + a['inflation rate'])**(year - 1)
            standard_cost = a['time_of_use_tariffs_high']['standard'] * (1 + a['inflation rate'])**(year - 1)
            off_peak_cost = a['time_of_use_tariffs_high']['off_peak'] * (1 + a['inflation rate'])**(year - 1) 
        else:
            peak_cost = a['time_of_use_tariffs_low']['peak'] * (1 + a['inflation rate'])**(year - 1)
            standard_cost = a['time_of_use_tariffs_low']['standard'] * (1 + a['inflation rate'])**(year - 1)
            off_peak_cost = a['time_of_use_tariffs_low']['off_peak'] * (1 + a['inflation rate'])**(year - 1)       
             
        if curr_hour_of_week > 120: # weekend is all off-peak
            total_cost_no_pv[i] = load_profile[i] * off_peak_cost 
            if net_load_profile[i] < 0: # If PV output is greater than EV load
                total_cost_with_pv[i] = 0 # No feed-in tariff, so excess supply is curtailed 
            else:
                total_cost_with_pv[i] = net_load_profile[i] * off_peak_cost        
        elif curr_hour_of_day in peak_hours: #peak 
            total_cost_no_pv[i] = load_profile[i] * peak_cost 
            if net_load_profile[i] < 0: 
                total_cost_with_pv[i] = 0 
            else:
                total_cost_with_pv[i] = net_load_profile[i] * peak_cost         
        elif curr_hour_of_day in standard_hours: #standard 
            total_cost_no_pv[i] = load_profile[i] * standard_cost      
            if net_load_profile[i] < 0:
                total_cost_with_pv[i] = 0
            else:
                total_cost_with_pv[i] = net_load_profile[i] * standard_cost     
        else: # off peak       
            total_cost_no_pv[i] = load_profile[i] * off_peak_cost      
            if net_load_profile[i] < 0:
                total_cost_with_pv[i] = 0
            else:
                total_cost_with_pv[i] = net_load_profile[i] * off_peak_cost      
                
    return total_cost_no_pv.sum() - total_cost_with_pv.sum() - cost_of_trickle_charging





def get_value_of_resiliency(year: int, kWh_affected_by_loadshedding: list, a: dict):

    # Obtain energy costs for each time period of the day
    peak_hours = a['time_periods']['peak_hours']
    standard_hours = a['time_periods']['standard_hours']
    off_peak_hours = a['time_periods']['off_peak_hours']


    # Initialize total cost variables
    val_kwh_missed = np.zeros(len(kWh_affected_by_loadshedding))
    
    
    for hour, kWh in enumerate(kWh_affected_by_loadshedding):
        # Find cost of doing with ICE
        kwh_L = a['kwh_km'] * (1/a['L_km'])  # kWh/L
        cost_of_ICE_operation_per_kwh = (a['cost_diesel'] / kwh_L) * (1 + a['inflation rate'])**(year - 1) # $/kWh
        
        curr_hour_of_week = hour % 168 
        curr_hour_of_day = hour % 24
        
        if (hour > a['high_period_start']) & (hour <= a['high_period_end']): # high period (all peak)
            peak_cost = a['time_of_use_tariffs_high']['peak'] * (1 + a['inflation rate'])**(year - 1)
            standard_cost = a['time_of_use_tariffs_high']['standard'] * (1 + a['inflation rate'])**(year - 1)
            off_peak_cost = a['time_of_use_tariffs_high']['off_peak'] * (1 + a['inflation rate'])**(year - 1) 
        else:
            peak_cost = a['time_of_use_tariffs_low']['peak'] * (1 + a['inflation rate'])**(year - 1)
            standard_cost = a['time_of_use_tariffs_low']['standard'] * (1 + a['inflation rate'])**(year - 1)
            off_peak_cost = a['time_of_use_tariffs_low']['off_peak'] * (1 + a['inflation rate'])**(year - 1)       
            
        if curr_hour_of_week > 120: # weekend 
            val_kwh_missed[hour] = kWh * (cost_of_ICE_operation_per_kwh - off_peak_cost)
        
        elif curr_hour_of_day in peak_hours:
            val_kwh_missed[hour] = kWh * (cost_of_ICE_operation_per_kwh - peak_cost)

        elif curr_hour_of_day in standard_hours:
            val_kwh_missed[hour] = kWh * (cost_of_ICE_operation_per_kwh - standard_cost)
                
        else:
            val_kwh_missed[hour] = kWh * (cost_of_ICE_operation_per_kwh - off_peak_cost)
                         
    return val_kwh_missed
