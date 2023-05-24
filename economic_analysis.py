import numpy as np


########### NPV ########### 

def calculate_npv(initial_investment, cash_flows, discount_rate):
    values = []
    for idx, cash_flow in enumerate(cash_flows):
        this_year_value = cash_flow /(1 + discount_rate)**idx
        values.append(this_year_value)
    total_benefits = sum(values)
    total_costs = initial_investment
    npv = total_benefits - total_costs
    return npv



########### Cost of charging ########### 

def get_cost_of_charging(load_profile: np.ndarray, net_load_profile: np.ndarray,
                         time_of_use_tariffs: dict, time_periods: dict,
                         feed_in_tariff: int,
                         feed_in_tariff_bool: bool):
    
    # Obtain energy costs for each time period of the day
    morning_cost = time_of_use_tariffs['morning']
    afternoon_cost = time_of_use_tariffs['afternoon']
    evening_cost = time_of_use_tariffs['evening']
    night_cost = time_of_use_tariffs['night']
    
    # Obtain time periods for each time period of the day
    morning_start = time_periods['morning_start']
    afternoon_start = time_periods['afternoon_start']
    evening_start = time_periods['evening_start']
    night_start = time_periods['night_start']
    
    # Initialize total cost variables
    total_cost_no_pv = np.zeros(len(load_profile))
    total_cost_with_pv = np.zeros(len(net_load_profile))
    
    # Calculate total cost of energy with and without PV

    for i in range(len(total_cost_no_pv)):
        curr_hour_of_day = i % 24
        
        if morning_start <= curr_hour_of_day < afternoon_start:
            # Without PV
            total_cost_no_pv[i] = load_profile[i] * morning_cost
            # With PV
            if net_load_profile[i] < 0: # PV output is greater than EV load
                
                if feed_in_tariff_bool: # if there is a feed-in-tariff 
                    total_cost_with_pv[i] = net_load_profile[i] * feed_in_tariff # apply feed-in tariff
                else:
                    total_cost_with_pv[i] = 0 # otherwise, energy is simply curtailed and there is no cost
            else:
                total_cost_with_pv[i] = net_load_profile[i] * morning_cost
            
        elif afternoon_start <= curr_hour_of_day < evening_start:
            
            total_cost_no_pv[i] = load_profile[i] * afternoon_cost
            
            if net_load_profile[i] < 0:
                if feed_in_tariff_bool:
                    total_cost_with_pv[i] = net_load_profile[i] * feed_in_tariff
                else:
                    total_cost_with_pv[i] = 0
            else:
                total_cost_with_pv[i] = net_load_profile[i] * afternoon_cost
                
        elif evening_start <= curr_hour_of_day < night_start:
            
            total_cost_no_pv[i] = load_profile[i] * evening_cost
            
            if net_load_profile[i] < 0:
                if feed_in_tariff_bool:
                    total_cost_with_pv[i] = net_load_profile[i] * feed_in_tariff
                else:
                    total_cost_with_pv[i] = 0
            else:
                total_cost_with_pv[i] = net_load_profile[i] * evening_cost
                
        else:
            
            total_cost_no_pv[i] = load_profile[i] * night_cost
            
            if net_load_profile[i] < 0:
                if feed_in_tariff_bool:
                    total_cost_with_pv[i] = net_load_profile[i] * feed_in_tariff
                else: 
                    total_cost_with_pv[i] = 0
                
            else:
                total_cost_with_pv[i] = net_load_profile[i] * night_cost
                
        
    return total_cost_no_pv, total_cost_with_pv 

########### Cost of loadshedding ########### 

def get_cost_of_missed_passengers_from_loadshedding(kWh_affected_by_loadshedding: list,
                                                     cost_per_passenger: float,
                                                     time_passenger_per_kWh: float, 
                                                     time_periods: dict):

    # Obtain energy costs for each time period of the day
    morning_passenger_per_kWh = time_passenger_per_kWh['morning'] # (float) number of passengers per kWh in the morning
    afternoon_passenger_per_kWh= time_passenger_per_kWh['afternoon']
    evening_passenger_per_kWh = time_passenger_per_kWh['evening']
    night_passenger_per_kWh = time_passenger_per_kWh['night']
    
    # Obtain time periods for each time period of the day
    morning_start = time_periods['morning_start']
    afternoon_start = time_periods['afternoon_start']
    evening_start = time_periods['evening_start']
    night_start = time_periods['night_start']
    
    # Initialize total cost variables
    passengers_missed = np.zeros(len(kWh_affected_by_loadshedding))
    
    # Calculate total cost of energy with and without PV

    for hour, kWh in enumerate(kWh_affected_by_loadshedding):
        
        curr_hour_of_day = hour % 24
        
        if morning_start <= curr_hour_of_day < afternoon_start:
            passengers_missed[hour] = kWh * morning_passenger_per_kWh

        elif afternoon_start <= curr_hour_of_day < evening_start:
            passengers_missed[hour] = kWh * afternoon_passenger_per_kWh
                
        elif evening_start <= curr_hour_of_day < night_start:
            passengers_missed[hour] = kWh * evening_passenger_per_kWh
            
        else:
            passengers_missed[hour] = kWh * night_passenger_per_kWh
            
    return passengers_missed * cost_per_passenger




########### Carbon offsets ########### 

def get_value_of_carbon_offsets(load_profile, net_load_profile, grid_carbon_intensity, carbon_price):
    carbon_cost_gross_load = load_profile.sum() * grid_carbon_intensity *  carbon_price  # kWh * kgCO2/kWh * $/kgCO2 = $
    carbon_cost_net_load = net_load_profile.sum() * grid_carbon_intensity *  carbon_price  # kWh * kgCO2/kWh * $/kgCO2 = $ 
    return carbon_cost_gross_load - carbon_cost_net_load