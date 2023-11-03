import numpy as np
import calculations
import generate_data
########### NPV ########### 

def calculate_npv(initial_investment, cash_flows, discount_rate):
    present_values = [cf / (1 + discount_rate) ** idx for idx, cf in enumerate(cash_flows)]
    npv = sum(present_values) - initial_investment
    return npv

def calculate_pv_capital_cost(pv_capacity, a):
    
    # find the inverter cost using a['inverter_cost_schedule'] and the pv_capacity 
    # to do so, must find the closest pv capacity in the schedule
    
    # if linearize_inverter_cost:
    inverter_cost = pv_capacity * a['inverter_cost_per_kw'] 
        
    # else:  # Use cost schedule
    #     pv_capacities = a['inverter_cost_schedule'].keys()
    #     closest_pv_capacity = min(pv_capacities, key=lambda x:abs(x-pv_capacity))
    #     inverter_cost = a['inverter_cost_schedule'][closest_pv_capacity]
        
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



def get_energy_served_by_system(pv_capacity, battery_capacity, load_profile, a):
    '''
    Energy from PV directly to EV load, not including energy that goes to the battery'''
    total_energy_served_by_pv = 0 
    total_energy_served_by_batt =0
    
    load_profile = np.array(load_profile) # for usage in np.where function
    
    battery_exists = True
    this_yr_battery_throughput=0
    repurchase_battery = True
    
    year_of_battery_lifetime = 0
    total_energy_served_by_system = 0
    
    for year in range(a['Rproj']):
 # Update PV and battery capacity after degradation   
        usable_pv_capacity = calculations.get_usable_pv_capacity(pv_capacity, year, a) 
        
        if battery_exists:
            annual_deg = (100 * (1 - a['battery_end_of_life_perc']))  / a['battery_lifetime_years'] / 100   # annual degradation rate (%)
            usable_battery_capacity =  battery_capacity * (1 - (annual_deg * year_of_battery_lifetime))
            year_of_battery_lifetime += 1 # update year of battery lifetime
            
        # Generate PV Output profile 
        pv_output_profile = generate_data.get_pv_output(a['annual_capacity_factor'], usable_pv_capacity) 
        
        # Initialize battery repurchae cost, residual value, and cost of trickle charging (these are all zero at beginning and change throughout)
        battery_repurchase_cost = 0
        battery_residual_value = 0
        cost_of_trickle_charging = 0
            
         ###########################################################################################
          # Battery
          ###########################################################################################  
        
        # Check if the battery is still alive 
        if year < a['battery_lifetime_years']:

            pv_with_battery_output_profile, this_yr_battery_throughput = generate_data.simulate_battery_storage_and_get_battery_throughput(pv_output_profile, usable_battery_capacity, a)
            
            
        else:
            if repurchase_battery: 
                year_of_battery_lifetime = 0
                pv_with_battery_output_profile, this_yr_battery_throughput = generate_data.simulate_battery_storage_and_get_battery_throughput(pv_output_profile, battery_capacity, a)
                
                if a['limit_battery_repurchases']:
                    repurchase_battery = False
            else:
                pv_with_battery_output_profile = pv_output_profile 
                battery_exists = False
                usable_battery_capacity = 0
                
        energy_served_by_pv = np.where((pv_with_battery_output_profile > 0) & (load_profile > 0), np.minimum(load_profile, pv_with_battery_output_profile), 0)
        energy_served_by_system = np.where((pv_with_battery_output_profile > 0) & (load_profile > 0), np.minimum(load_profile, pv_with_battery_output_profile), 0)
        total_energy_served_by_system += energy_served_by_system.sum()
      #  total_energy_served_by_batt += this_yr_battery_throughput
       # total_energy_served_by_pv += energy_served_by_pv.sum() - this_yr_battery_throughput
        
        
    return total_energy_served_by_system

def value_of_resiliency(pv_capacity, battery_capacity, a):
    total_val_kwh_ls=0
    repurchase_battery = a['repurchase_battery']
    battery_exists = True 
    year_of_battery_lifetime = 0
    
    for year in range(a['Rproj']):
        usable_pv_capacity = calculations.get_usable_pv_capacity(pv_capacity, year, a)   
        if battery_exists:
            annual_deg = (100 * (1 - a['battery_end_of_life_perc']))  / a['battery_lifetime_years'] / 100   # annual degradation rate (%)
            usable_battery_capacity =  battery_capacity * (1 - (annual_deg * year_of_battery_lifetime))
            year_of_battery_lifetime += 1 # update year of battery lifetime
                      
        # Generate PV Output profile 
        pv_output_profile = generate_data.get_pv_output(a['annual_capacity_factor'], usable_pv_capacity)
        
        # Check if the battery is still alive 
        if year < a['battery_lifetime_years']:
            pv_with_battery_output_profile, cost_of_trickle_charging = generate_data.simulate_battery_storage_v5(pv_output_profile, usable_battery_capacity, a)

        else:
            if repurchase_battery: 
                year_of_battery_lifetime = 0
                pv_with_battery_output_profile, cost_of_trickle_charging = generate_data.simulate_battery_storage_v5(pv_output_profile, battery_capacity, a)
                if a['limit_battery_repurchases']:
                    repurchase_battery = False
            else:
                pv_with_battery_output_profile = pv_output_profile 
                battery_exists = False
                usable_battery_capacity = 0
                
                
        loadshedding_schedule = a['load_shedding_schedule']
        net_load_profile = a['load_profile'] - pv_with_battery_output_profile
        gross_load_lost_to_loadshedding = np.array([a['load_profile'][i] if is_shedding else 0 for i, is_shedding in enumerate(loadshedding_schedule)])
        
        # Profile of kWh that would have been lost to loadshedding but are saved by the solar + battery generation [these are beneficial, and not to be charged $$ for]
        saved_free_kWh = [min(pv_with_battery_output_profile[i], gross_load_lost_to_loadshedding[i]) if is_shedding else 0 for i, is_shedding in enumerate(loadshedding_schedule)]
        
        # Profile of kWh that would be lost to load shedding WITH solar and battery
        net_load_lost_to_loadshedding = np.array([net_load_profile[i] if is_shedding and net_load_profile[i] > 0 else 0 for i, is_shedding in enumerate(loadshedding_schedule)])
        gross_load_minus_loadshedding = a['load_profile'] - gross_load_lost_to_loadshedding
        net_load_minus_loadshedding = net_load_profile - net_load_lost_to_loadshedding 
        value_of_charging_saved_by_pv_from_loadshedding = get_cost_of_missed_passengers_from_loadshedding_v2(year, saved_free_kWh, a)
        val_kwh_ls = value_of_charging_saved_by_pv_from_loadshedding.sum()
        
        total_val_kwh_ls += val_kwh_ls
    return total_val_kwh_ls


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






