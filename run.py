import time
import os
import multiprocessing as mp
import optimization
import economic_analysis
from constraints import cons  
import input_params as input 
import time 
import numpy as np
from scipy.optimize import minimize 
from scipy.optimize import differential_evolution
from tqdm import tqdm

###### OPTIMIZE SYSTEM varying load profile and solar and storage ##########
def optimize_system(solar_cost, battery_cost, load_profile, grid_load_shedding_schedule, scenario_name, bounds, a):
    this_run_start_time = time.time()
    a['load_profile'] = load_profile
    a['pv_cost_per_kw'] = solar_cost
    a['battery_cost_per_kWh'] = battery_cost
    a['load_shedding_schedule'] = grid_load_shedding_schedule
  
    if sum(grid_load_shedding_schedule) > 0:
        load_shedding_scenario = "Grid LS on"
    else:
        load_shedding_scenario = "Grid LS off"
    print(load_shedding_scenario)
    
    initial_guess = [500, 500]
    
    if '100%' in scenario_name:
        a['full_ev_fleet'] = True

    # Run 
    result = minimize(optimization.objective_function, x0 = initial_guess, args = (a,), bounds=bounds, constraints = cons , method= 'SLSQP')
     

    # Extract results         
    optimal_pv_capacity = result.x[0]
    optimal_battery_capacity = result.x[1]

    max_npv = result.fun

    this_run_end_time = time.time()
    this_run_time = (this_run_end_time - this_run_start_time) / 60

    return optimal_pv_capacity, optimal_battery_capacity, -max_npv, this_run_time, scenario_name, load_shedding_scenario



###############  STARMAP #################### 
if __name__ == "__main__":
    
    start_time = time.time()

    # Loading in the load profiles.
    annual_25_perc_ev = np.loadtxt(f"processed_ev_schedule_data/annual_25_perc_ev.txt")
    annual_50_perc_ev = np.loadtxt(f"processed_ev_schedule_data/annual_50_perc_ev.txt")
    annual_75_perc_ev = np.loadtxt(f"processed_ev_schedule_data/annual_75_perc_ev.txt")
    annual_100_perc_ev = np.loadtxt(f"processed_ev_schedule_data/annual_100_perc_ev.txt")
    
    annual_25_perc_ls_1 = np.loadtxt(f"processed_ev_schedule_data/25_perc/annual_ls_1.txt") 
    annual_50_perc_ls_1 = np.loadtxt(f"processed_ev_schedule_data/50_perc/annual_ls_1.txt") 
    annual_75_perc_ls_1 = np.loadtxt(f"processed_ev_schedule_data/75_perc/annual_ls_1.txt") 
    annual_100_perc_ls_1 = np.loadtxt(f"processed_ev_schedule_data/100_perc/annual_ls_1.txt") 
    
    annual_25_perc_ls_2 = np.loadtxt(f"processed_ev_schedule_data/25_perc/annual_ls_2.txt") 
    annual_50_perc_ls_2 = np.loadtxt(f"processed_ev_schedule_data/50_perc/annual_ls_2.txt") 
    annual_75_perc_ls_2 = np.loadtxt(f"processed_ev_schedule_data/75_perc/annual_ls_2.txt") 
    annual_100_perc_ls_2 = np.loadtxt(f"processed_ev_schedule_data/100_perc/annual_ls_2.txt") 
    
    annual_25_perc_ls_3 = np.loadtxt(f"processed_ev_schedule_data/25_perc/annual_ls_3.txt") 
    annual_50_perc_ls_3 = np.loadtxt(f"processed_ev_schedule_data/50_perc/annual_ls_3.txt") 
    annual_75_perc_ls_3 = np.loadtxt(f"processed_ev_schedule_data/75_perc/annual_ls_3.txt") 
    annual_100_perc_ls_3 = np.loadtxt(f"processed_ev_schedule_data/100_perc/annual_ls_3.txt") 
    
    
    
    
     # Possible input for load_profiles
    load_profile_list_ls_1 = [annual_25_perc_ls_1,annual_50_perc_ls_1, annual_75_perc_ls_1, annual_100_perc_ls_1]
    load_profile_list_ls_2 = [annual_25_perc_ls_2, annual_50_perc_ls_2, annual_75_perc_ls_2, annual_100_perc_ls_2]
    load_profile_list_ls_3 = [annual_25_perc_ls_3, annual_50_perc_ls_3, annual_75_perc_ls_3, annual_100_perc_ls_3]
    
    # Possible inputs for grid_load_shedding_schedules = [input.ls_annual_empty, input.ls_annual_1, input.ls_annual_2, input.ls_annual_3]
    grid_load_shedding_schedules = [input.ls_annual_1, input.ls_annual_2, input.ls_annual_3]
    
    # Ensure that scenario names line up with load profile + grid load shedding schedules used.
    scenario_names = ["100% LS1", "100% LS2", "100% LS3"]
   
    load_profile_list = [annual_100_perc_ev]

    # Solar panel costs per kW to consider
    solar_costs = [500]
    # Battery costs per kWh to consider
    battery_costs = [300]
   
    # Limit on [(solar),(battery)] sizes
    bounds = [(0,1000), (0,1000)]
    a = input.a 
    pool = mp.Pool(processes=mp.cpu_count())

    total_combinations = len(load_profile_list) * len(solar_costs) * len(battery_costs) #len(grid_load_shedding_schedules)
    combinations = []

    # Running the combos 
    for idx, load_profile in enumerate(load_profile_list):
        for jdx, grid_load_shedding_schedule in enumerate(grid_load_shedding_schedules): 
            for solar_cost in solar_costs:
                for battery_cost in battery_costs:
                    combinations.append((solar_cost, battery_cost, load_profile, grid_load_shedding_schedule, scenario_names[jdx], bounds, a))

    # Use starmap to parallelize the optimization process
    results = list(tqdm(pool.starmap(optimize_system, combinations), total=total_combinations, desc="Progress"))

    pool.close()
    pool.join()
    
        
    for idx, result in enumerate(results):
        solar_cost, battery_cost, load_profile, grid_load_shedding_schedule, scenario_name, _, _, = combinations[idx]
        optimal_pv_capacity, optimal_battery_capacity, npv_value, execution_time, scenario_name, grid_load_shedding_scenario = result
        
        
        #Make directory scenario_name if doesn't exist
        if not os.path.exists(f"../results/{grid_load_shedding_scenario}/{scenario_name}"):   
            os.makedirs(f"../results/{grid_load_shedding_scenario}/{scenario_name}")
            
        # Create the file name
        file_name = f"../results/{grid_load_shedding_scenario}/{scenario_name}/solar={solar_cost}_battery={battery_cost}.txt"
        
        # Write the results to the file
        with open(file_name, "w") as f2:
            f2.write(f"Optimal PV Capacity: {optimal_pv_capacity} kW\n")
            f2.write(f"Optimal Battery Capacity: {optimal_battery_capacity} kWh\n")
            f2.write(f"NPV: {npv_value}\n")
            # Calculate investment cost, ROI, and LCOE
            investment_cost = economic_analysis.calculate_pv_capital_cost(optimal_pv_capacity, a) + (a['battery_cost_per_kWh'] * optimal_battery_capacity)
            f2.write(f"Investment cost: ${investment_cost}\n") 
            roi = npv_value / investment_cost
            f2.write(f"ROI: {roi} %\n")
            lcoe_pv = economic_analysis.calculate_lcoe_pv(optimal_pv_capacity, optimal_battery_capacity, load_profile, a)
            lcoe_batt = economic_analysis.calculate_lcoe_batt(optimal_pv_capacity, optimal_battery_capacity, a)
            f2.write(f"LCOE_PV: {lcoe_pv}$/kWh \n")
            f2.write(f"LCOE_Batt: {lcoe_batt}$/kWh \n")
            p_grid = economic_analysis.compute_p_grid(load_profile, a)
            f2.write(f"P_grid: {p_grid}$/kWh \n")
            lps = load_profile.sum()
            f2.write(f"Load profile sum {lps} \n")
            kwh_ls, op_savings = economic_analysis.get_kwh_ls_and_op_savings(optimal_pv_capacity, optimal_battery_capacity, a)
            f2.write(f"kwh_ls: {kwh_ls} \n")
            f2.write(f"op_savings: {op_savings} \n")
            # demand served by system
            d_sys = economic_analysis.get_demand_served_by_system(optimal_pv_capacity, optimal_battery_capacity, a)
            f2.write(f"d_sys: {d_sys} \n")
            carbon_savings = d_sys * 0.95 # SA grid intensity = 0.95kgCO2/kWh
            f2.write(f"carbon_savings: {carbon_savings} \n")
            
            coe_w_pv, coe_no_pv = economic_analysis.get_cost_of_energy(optimal_pv_capacity, optimal_battery_capacity, a)
            f2.write(f"coe_w_pv: {coe_w_pv} \n")
            f2.write(f"coe_no_pv: {coe_no_pv} \n")
            
        
        # Print the results
        print(f"Combination {idx + 1}/{total_combinations}:")
        print("Solar Cost:", solar_cost)
        print("Battery Cost:", battery_cost)
        print("Scenario:", scenario_name)  # You might want to print the name or other information about the load profile
        print("Optimal PV Capacity:", optimal_pv_capacity, 'kW')
        print("Optimal Battery Capacity:", optimal_battery_capacity, 'kWh')
        print("NPV:", npv_value)
        print("Investment cost: $", investment_cost) 
        print("LCOE PV", lcoe_pv)
        print("LCOE Batt", lcoe_batt)
        print("Avg grid price", p_grid)
        print("Execution Time (minutes):", execution_time)
        print("kwh_ls:", kwh_ls)
        print('op_savings:', op_savings)
        print("d_sys:", d_sys)
        print("carbon_savings:", carbon_savings)
        print()

            
    end_time = time.time()
    
    total_execution_time = (end_time - start_time)/60
    
    print("Total program execution time:", total_execution_time)
    # print number of true observations in a['load_shedding_schdule'] array
    

    