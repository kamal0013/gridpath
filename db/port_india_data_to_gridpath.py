#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

"""
Port India data to GridPath
"""

from collections import OrderedDict
import os
import sqlite3
import pandas as pd
from calendar import monthrange
import numpy as np

# Data-import modules

from db.utilities import temporal, geography, project_list, project_zones, \
    project_operational_chars, project_availability, fuels, \
    project_portfolios, project_existing_params, project_new_costs, \
    project_new_potentials, project_prm, transmission_portfolios, \
    transmission_zones, transmission_capacities, simultaneous_flow_groups, \
    simultaneous_flows, transmission_hurdle_rates, carbon_cap, system_load, \
    system_reserves, system_prm, rps, scenario

## MODULES FOR PORTING INDIA DATA TO SQL DATABASE
from db.india_data_functions import port_data_to_gridpath_demand_modules, port_data_to_gridpath_project_modules



## INPUTS
# SQL database
sql_database = 'india.db'

# Input csv path
inputPath = os.path.join(os.getcwd(), "db", "india_data")
dbPath = os.path.join(os.getcwd(), "db")

# csv with periods data
periods_csv = os.path.join(inputPath, 'periods_4periods.csv')

# csv with horizons data
horizons_csv = os.path.join(inputPath, 'horizons_12_1dayPerMonth.csv')

# csv with existing generators data
existing_projects_csv = os.path.join(inputPath, 'gen_all_input_cc_ccgt_diesel.csv')

# csv with new variable generation projects data
new_projects_wind_csv = os.path.join(inputPath, 'variable_gen_projects/wind_candidate_projects.csv')
new_projects_solar_csv = os.path.join(inputPath, 'variable_gen_projects/solar_candidate_projects.csv')
new_projects_wind_timeseries_csv = os.path.join(inputPath, 'variable_gen_projects/wind_candidate_project_timeseries.csv')
new_projects_solar_timeseries_csv = os.path.join(inputPath, 'variable_gen_projects/solar_candidate_project_timeseries.csv')

# csv with hourly demand forecasts
load_csv = os.path.join(inputPath, 'demand/demand_2014_2032_hourly.csv')

# csv with year, month, days selected - sample days
year_month_days_selected_csv = os.path.join(inputPath, 'demand/year_month_days_selected.csv')
year_month_days_selected = pd.read_csv(year_month_days_selected_csv)

# csv with hydro data
hydro_max_energy_csv = os.path.join(inputPath, 'hydro_max_energy.csv')
hydro_min_gen_csv = os.path.join(inputPath, 'hydro_min_gen.csv')
mustrun_gen_csv = os.path.join(inputPath, 'mustrun_gen.csv')
# This is the curated hydro parameters output file
projects_hydro_horizon_chars_csv = os.path.join(inputPath, 'projects_hydro_horizon_chars.csv')
# This is the curated variable non-curtailable time series
projects_variable_noncurtailable_timeseries_csv = os.path.join(inputPath, 'variable_gen_projects/variable_noncurtailable_timeseries.csv')

# data with heat rate curves and variable om costs
projects_hr_var_costs_excel = os.path.join(inputPath, 'projects_hr_var_costs.xlsx')

# data with fuel costs
fuels_costs_csv = os.path.join(inputPath, 'Sys_Fuel_Costs.csv')

# data with new projects costs and operational characteristics
new_projects_excel = os.path.join(inputPath, 'new_projects.xlsx')

# main_scenarios to run
main_scenarios_csv = os.path.join(inputPath, 'main_scenarios.csv')

# RPS targets by zones and scenarios
rps_targets_excel = os.path.join(inputPath, 'rps_targets.xlsx')

discount_rate = 0.07
base_year = 2018
horizons_per_year = 12
number_of_hours_in_timepoint = 1  # 0.25 for 15 minutes
inr_to_usd = 70

# TIME INDEXES
day_array_365 = np.arange(1,366)
day_array_8760 = np.repeat(day_array_365, 24)
hour_array_8760 = np.arange(1,8761)
hour_of_day_array = np.tile(np.arange(1,25), 365)

# Connect to database
io = sqlite3.connect(
    os.path.join(dbPath, sql_database)
)

c2 = io.cursor()


#### Periods and Horizons and Timepoints ####

def create_periods():
    """
    Create periods dataframe
    :return periods_df:
    """
    periods_df = pd.DataFrame()
    periods_df = pd.read_csv(periods_csv, sep=',')

    # Calculate the discount factors
    periods_df['discount_factor'] = 1 / (
                (1 + discount_rate) ** (periods_df['period'] - base_year)) # CHECK THIS FORMULA. Removed +1 after base_year
    return periods_df

def create_horizons():
    """
    Create horizons array and horizon_weights_and_months_df
    :return horizons array:
    :return horizon_weights_and_months_df:
    """
    # global horizon_array  # So other modules can use this variable
    # global month_array
    # global horizon_weights_and_months_df
    horizon_weights_and_months_df = pd.DataFrame()
    if horizons_csv != '':
        horizon_weights_and_months_df = pd.read_csv(horizons_csv, sep=',')
        horizon_array = horizon_weights_and_months_df['horizon'].values
        #month_array = horizon_weights_and_months_df['month'].values
    else:
        horizon_array = np.arange(1, horizons_per_year + 1, 1)
        weights_array = np.repeat(1, horizons_per_year)
        # weights_array = np.repeat(1, (366 if calendar.isleap(base_year) else 365))
        month_array = []
        # This month array is for 1 day horizon with 365 days. For less number of days, use a csv.
        for m in range(1, 12 + 1, 1):
            month_array = np.append(month_array, np.repeat(m, monthrange(base_year, m)[1])).astype(int)
        horizon_weights_and_months_df = pd.DataFrame(
            {'horizon': horizon_array, 'weight': weights_array, 'month': month_array})

    return (horizon_array, horizon_weights_and_months_df)


def create_period_horizon_timepoints():
    """
    Create period, horizon, timepoint dataframe
    :return period_horizon_timepoints:
    """

    period_horizon_timepoints = pd.DataFrame()

    for p in periods_df['period'].to_list():
        print(p)
        period_horizon_timepoints_p = pd.DataFrame()
        period_horizon_timepoints_p['period'] = np.repeat(p, horizons_per_year * 24 * 1 / number_of_hours_in_timepoint)
        period_horizon_timepoints_p['horizonstep'] = np.repeat(range(1, horizons_per_year + 1, 1),
                                                               24 * 1 / number_of_hours_in_timepoint)
        # Multiply by 3 if horizons are 3 digit e.g. for production cost models; else 2
        period_horizon_timepoints_p['horizon'] = period_horizon_timepoints_p['period'] * 10 ** 2 + \
                                                 period_horizon_timepoints_p['horizonstep']
        period_horizon_timepoints_p['timestep'] = np.tile(range(1, int(24 * 1 / number_of_hours_in_timepoint + 1), 1), horizons_per_year)
        # Multiply by 5 if horizons are 3 digit e.g. for production cost models; else 4
        period_horizon_timepoints_p['tmp'] = period_horizon_timepoints_p['period'] * 10 ** 4 + \
                                             period_horizon_timepoints_p['horizonstep'] * 10 ** 2 + \
                                             period_horizon_timepoints_p['timestep']
        period_horizon_timepoints = period_horizon_timepoints.append(period_horizon_timepoints_p, ignore_index=True)
        return period_horizon_timepoints

def create_period_horizon():
    """
    Create period, horizon dataframe
    :return period_horizon:
    """

    period_horizon = pd.DataFrame()

    for p in periods_df['period'].to_list():
        print(p)
        period_horizon_p = pd.DataFrame()
        period_horizon_p['period'] = np.repeat(p, horizons_per_year)
        period_horizon_p['horizonstep'] = range(1, horizons_per_year + 1, 1)
        period_horizon_p['horizon'] = period_horizon_p['period'] * 10 ** 2 + \
                                      period_horizon_p['horizonstep']
        period_horizon = period_horizon.append(period_horizon_p, ignore_index=True)

    return period_horizon


# Using the same variable names inside and outside function. Need to change.
periods_df = create_periods()
(horizon_array, horizon_weights_and_months_df) = create_horizons()
horizon_weights_and_months_df = horizon_weights_and_months_df.set_index("horizon")
period_horizon_timepoints = create_period_horizon_timepoints()
period_horizon = create_period_horizon()

def load_temporal_data():
    """
    Add timepoints/days into database
    Investment periods are 2018, 2022, 2026, and 2030
    Discount factors are calculated using an input discount rate and a base year
    Number of years represented by a period is 4 but arbitrarily set to a
    higher number for the last period
    Horizon boundary is circular
    :return:
    """

    # PERIODS
    periods = dict(dict())
    for p in periods_df['period']:
        periods[p] = {}
        periods[p]['discount_factor'] = periods_df.loc[
            periods_df["period"] == p]['discount_factor'].iloc[0]
        periods[p]['number_years_represented'] = float(periods_df.loc[
            periods_df["period"] == p]['number_years_represented'].iloc[0])

    # TODO: CAN WE USE THIS ALTERNATE CODE TO POPULATE DICTIONARIES?
    periods_alt = dict(dict())
    periods_df_alt = periods_df.copy()
    periods_df_alt = periods_df_alt.set_index("period")
    periods_alt = periods_df_alt.to_dict(orient='index')

    # HORIZONS
    horizon_weights_and_months = horizon_weights_and_months_df.to_dict(orient='index')

    # Subproblems and stages (one subproblem, one stage)
    subproblems = [1]
    subproblem_stages = {sid: [(1, "single stage")] for sid in subproblems}

    # Timepoints
    subproblem_stage_timepoints = dict()
    for subproblem_id in subproblem_stages.keys():
        print(subproblem_id)
        subproblem_stage_timepoints[subproblem_id] = dict()
        for stage in subproblem_stages[subproblem_id]:
            stage_id = stage[0]
            subproblem_stage_timepoints[subproblem_id][stage_id] = dict()
            for _period in periods.keys():
                for _day in horizon_weights_and_months.keys():
                    for hour in range(1, 25):
                        timepoint = _period * 10**4 + _day * 10**2 + hour
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint] = dict()
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint]["period"] = _period
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint][
                            "number_of_hours_in_timepoint"] = 1
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint]["timepoint_weight"] = \
                            horizon_weights_and_months[_day]["weight"]
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint][
                            "previous_stage_timepoint_map"] = None
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint][
                            "spinup_or_lookahead"] = None
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint]["month"] = \
                            int(horizon_weights_and_months[_day]["month"])
                        subproblem_stage_timepoints[subproblem_id][
                            stage_id][timepoint]["hour_of_day"] = hour

    # Horizons
    subproblem_horizons = dict()
    for subproblem_id in subproblem_stages.keys():
        subproblem_horizons[subproblem_id] = dict()
        for period in periods.keys():
            for day in horizon_weights_and_months.keys():
                horizon = period * 10**2 + day
                subproblem_horizons[subproblem_id][horizon] = dict()
                subproblem_horizons[subproblem_id][horizon]["period"] = period
                subproblem_horizons[subproblem_id][horizon]["boundary"] = \
                    "circular"
                subproblem_horizons[subproblem_id][horizon][
                    "balancing_type_horizon"] = "day"

    # Timepoint horizons
    subproblem_stage_timepoint_horizons = dict()
    for subproblem_id in subproblem_stage_timepoints.keys():
        subproblem_stage_timepoint_horizons[subproblem_id] = dict()
        for stage_id in subproblem_stage_timepoints[subproblem_id].keys():
            subproblem_stage_timepoint_horizons[subproblem_id][stage_id] = \
                dict()
            for timepoint in subproblem_stage_timepoints[subproblem_id][
                    stage_id].keys():
                subproblem_stage_timepoint_horizons[subproblem_id][
                    stage_id][timepoint] = [(int(timepoint/10**2), 'day')]

    # Load data into GridPath database
    temporal.temporal(
            io=io, c=c2,
            temporal_scenario_id=1,
            scenario_name="default_4_periods_12_days_24_hours",
            scenario_description="2018, 2022, 2026, 2030; 12 average days, "
                                 "24 hours each",
            periods=periods,
            subproblems=[1],
            subproblem_stages={1: [(1, "single stage")]},
            subproblem_stage_timepoints=subproblem_stage_timepoints,
            subproblem_horizons=subproblem_horizons,
            subproblem_stage_timepoint_horizons=subproblem_stage_timepoint_horizons
    )


def load_geography_load_zones():
    """
    Only single zone
    :return:
    """

    # Load data into GridPath database
    geography.geography_load_zones(
        io=io, c=c2,
        load_zone_scenario_id=1,
        scenario_name='default_1_zone',
        scenario_description = '1 zone: India',
        zones=['India'],
        zone_overgen_penalties={
            'India': 50000
        },
        zone_unserved_energy_penalties={
            'India': 50000
        }
    )


def load_geography_rps_zones():
    """
    single rps zone
    :return:
    """
    geography.geography_rps_zones(
        io=io, c=c2,
        rps_zone_scenario_id=1,
        scenario_name='india_rps',
        scenario_description='INDIA RPS',
        zones=['India']
    )


#### PROCESS RAW PROJECT DATA ####

def create_projects_data():
    """
    Create the projects data table with existing and new generators and their properties
    :return new_projects_var_gen
    :return projects_all:
    :return existing_projects_capacities:
    """

    # Existing project list
    existing_project_units = pd.read_csv(existing_projects_csv, sep=',')

    # Reclassify projects
    existing_project_units.loc[(existing_project_units['gen_capacity'] <= 500) & (existing_project_units['type'] == 'coal'), 'type'] = 'Existing_Subcritical_Coal'
    existing_project_units.loc[(existing_project_units['gen_capacity'] > 500) & (
                existing_project_units['type'] == 'coal'), 'type'] = 'Existing_Supercritical_Coal'
    existing_project_units.loc[(existing_project_units['type'] == 'gas_ccgt'), 'type'] = 'Existing_CCGT'
    existing_project_units.loc[(existing_project_units['type'] == 'gas_ct'), 'type'] = 'Existing_CT'
    existing_project_units.loc[(existing_project_units['generator'] == 'HYDRO-ROR'), 'type'] = 'Existing_Hydro_ROR'
    existing_project_units.loc[(existing_project_units['generator'] == 'HYDRO-PONDAGE'), 'type'] = 'Existing_Hydro_Pondage'
    existing_project_units.loc[(existing_project_units['generator'] == 'HYDRO-STORAGE'), 'type'] = 'Existing_Hydro_Storage'
    existing_project_units.loc[(existing_project_units['generator'] == 'NUCLEAR'), 'type'] = 'Existing_Nuclear'
    existing_project_units.loc[(existing_project_units['type'] == 'diesel'), 'type'] = 'Existing_Diesel'
    existing_project_units.loc[(existing_project_units['type'] == 'other'), 'type'] = 'Existing_Other'

    # Group existing generator units and estimate mean variable costs, mean unit size, and sum of capacity # THIS CODE IS FOR CAPACITY EXPANSION AGGREGATION
    mean = lambda x: x.mean()
    existing_project_units.loc[:, 'variable_om_cost_per_mwh'] = existing_project_units.groupby(['type']).var_cost.transform(mean) / inr_to_usd
    existing_project_units.loc[:, 'unit_size_mw_mean'] = existing_project_units.groupby(
        ['type']).gen_capacity.transform(mean)
    existing_project_units.loc[:, 'gen_capacity_mw_sum'] = existing_project_units.groupby(
        ['type']).gen_capacity.transform(sum)
    # Drop generator unit level data and then drop duplicates to create the aggregate dataframe for existing projects
    # Use copy so the new dataframe is not a view of the original df. Also avoids "SettingWithCopyWarning".
    existing_projects = existing_project_units[['type', 'variable_om_cost_per_mwh', 'unit_size_mw_mean']].copy()
    existing_projects.rename(index = str, columns = {'type': 'project'}, inplace = True)
    # Drop duplicates leaving only aggregated projects
    existing_projects = existing_projects.drop_duplicates()

    # Set unit size as the nearest multiple of 10.
    existing_projects['unit_size_mw'] = np.ceil(existing_projects['unit_size_mw_mean']/10)*10

    # TODO: This is arbitrarily setting hydro and nuclear unit sizes to 500 MW. Not sure how this will be treated in GP. The total capacity is not a multiple of 500.
    existing_projects.loc[existing_projects['project'].isin(['Existing_Hydro_Storage', 'Existing_Hydro_ROR', 'Existing_Hydro_Pondage', 'Existing_Nuclear']), 'unit_size_mw'] = 500

    # Add technology and fuel to existing generators

    #existing_projects = existing_projects.reindex(columns = ['project', 'variable_om_cost_per_mwh', 'capacity_type', 'operational_type', 'technology', 'fuel', 'unit_size_mw'])
    existing_projects = existing_projects.reindex(
        columns=['project', 'variable_om_cost_per_mwh', 'capacity_type', 'operational_type', 'balancing_type_project', 'technology', 'fuel',
                       'heat_rate_curves_scenario_id', 'min_stable_level', 'unit_size_mw',
                       'startup_cost_per_mw', 'shutdown_cost_per_mw', 'startup_plus_ramp_up_rate', 'shutdown_plus_ramp_down_rate',
                       'ramp_up_when_on_rate', 'ramp_down_when_on_rate', 'min_up_time_hours', 'min_down_time_hours',
                       'charging_efficiency', 'discharging_efficiency', 'minimum_duration_hours', 'regulation_up_ba',
                       'regulation_up_derate', 'regulation_down_ba', 'regulation_down_derate', 'spinning_reserves_ba',
                       'spinning_reserves_derate', 'regulation_up_ramp_rate', 'regulation_down_ramp_rate',
                       'spinning_reserves_ramp_rate', 'carbon_cap_zone'])

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_Supercritical_Coal'), ['capacity_type', 'operational_type',
                                                                                            'balancing_type_project',
                                                                                            'technology', 'fuel',
                                                                                            'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'dispatchable_capacity_commit', 'day', 'Coal', 'Coal_existing', 640]

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_Subcritical_Coal'), ['capacity_type', 'operational_type',
                                                                                          'balancing_type_project',
                                                                                          'technology', 'fuel',
                                                                                          'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'dispatchable_capacity_commit', 'day', 'Coal', 'Coal_existing', 200]

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_CCGT'), ['capacity_type', 'operational_type',
                                                            'balancing_type_project', 'technology', 'fuel',
                                                            'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'dispatchable_capacity_commit', 'day', 'CCGT', 'Gas_domestic', 120]

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_CT'), ['capacity_type', 'operational_type',
                                                          'balancing_type_project', 'technology', 'fuel',
                                                          'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'dispatchable_capacity_commit', 'day', 'Peaker', 'Gas_domestic', 30]

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_Diesel'), ['capacity_type', 'operational_type',
                                                              'balancing_type_project', 'technology', 'fuel',
                                                              'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'dispatchable_capacity_commit', 'day', 'Peaker', 'Diesel', 30]

    ## TODO: Double check this one
    existing_projects.loc[(existing_projects['project'] == 'Existing_Hydro_ROR') | (
            existing_projects['project'] == 'Existing_Hydro_Pondage'), ['capacity_type', 'operational_type',
                                                                        'balancing_type_project', 'technology',
                                                                        'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'variable_no_curtailment', 'day', 'Hydro', 500]

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_Hydro_Storage'), ['capacity_type', 'operational_type',
                                                                     'balancing_type_project', 'technology',
                                                                     'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'hydro_noncurtailable', 'day', 'Hydro', 500]

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_Nuclear'), ['capacity_type', 'operational_type',
                                                               'balancing_type_project', 'technology', 'fuel',
                                                               'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'must_run', 'day', 'Nuclear', 'Uranium', 500]

    existing_projects.loc[
        (existing_projects['project'] == 'Existing_Other'), ['capacity_type', 'operational_type',
                                                             'balancing_type_project', 'technology',
                                                             'unit_size_mw']] = [
        'existing_gen_no_economic_retirement', 'must_run', 'day', 'Biomass', 200]


    ## Add heat rates, start up shutdown costs to existing projects

    existing_projects.loc[(existing_projects['project'] == 'Existing_Supercritical_Coal'), ['min_stable_level',
                                                                                            'min_up_time_hours',
                                                                                            'min_down_time_hours',
                                                                                          'startup_cost_per_mw',
                                                                                          'shutdown_cost_per_mw']] = [
        0.6, 12, 12, 69.3, 69.3]

    existing_projects.loc[(existing_projects['project'] == 'Existing_Subcritical_Coal'), ['min_stable_level',
                                                                                          'min_up_time_hours',
                                                                                          'min_down_time_hours',
                                                                                          'startup_cost_per_mw',
                                                                                          'shutdown_cost_per_mw']] = [
        0.6, 12, 12, 69.3, 69.3]

    existing_projects.loc[(existing_projects['project'] == 'Existing_CCGT'), ['min_stable_level',
                                                                                            'min_up_time_hours',
                                                                                            'min_down_time_hours',
                                                                                          'startup_cost_per_mw',
                                                                                          'shutdown_cost_per_mw']] = [
        0.5, 6, 6, 100, 100]

    existing_projects.loc[(existing_projects['project'] == 'Existing_CT'), ['min_stable_level',
                                                                                            'min_up_time_hours',
                                                                                            'min_down_time_hours',
                                                                                          'startup_cost_per_mw',
                                                                                          'shutdown_cost_per_mw']] = [
        0.4, 1, 1, 100, 100]

    existing_projects.loc[(existing_projects['project'] == 'Existing_Diesel'), ['min_stable_level',
                                                                                            'min_up_time_hours',
                                                                                            'min_down_time_hours',
                                                                                          'startup_cost_per_mw',
                                                                                          'shutdown_cost_per_mw']] = [
        0.4, 1, 1, 100, 100]

    ## Add heat rate curve id for conventional generators
    existing_projects.loc[existing_projects['project'].isin(['Existing_CCGT', 'Existing_CT', 'Existing_Diesel',
                                                             'Existing_Subcritical_Coal', 'Existing_Supercritical_Coal',
                                                             'Existing_Nuclear']),
                          ['heat_rate_curves_scenario_id']] = 1

    ## Add reserve zones
    existing_projects.loc[existing_projects['project'].isin(['Existing_CCGT', 'Existing_CT', 'Existing_Diesel',
                                                             'Existing_Subcritical_Coal', 'Existing_Supercritical_Coal',
                                                             'Existing_Hydro_Storage']),
                          ['regulation_up_ba', 'regulation_down_ba', 'spinning_reserves_ba']] =['India', 'India', 'India']

    ## Add carbon cap zone
    existing_projects.loc[existing_projects['project'].isin(['Existing_CCGT', 'Existing_CT', 'Existing_Diesel',
                                                             'Existing_Subcritical_Coal', 'Existing_Supercritical_Coal']),
                          'carbon_cap_zone'] ='India'

    ## Add load zone [for single load zone
    existing_projects['load_zone'] = 'India'



    ## New project operational chars data
    new_projects = pd.read_excel(new_projects_excel, sheet_name='new_projects_op_chars')

    ## New variable generation project list
    # global new_projects_var_gen
    new_projects_wind = pd.read_csv(new_projects_wind_csv)
    new_projects_wind['technology'] = 'Wind'
    new_projects_solar = pd.read_csv(new_projects_solar_csv)
    new_projects_solar['technology'] = 'Solar'
    new_projects_var_gen = pd.concat([new_projects_solar, new_projects_wind], axis=0)

    new_projects_var_gen['capacity_type'] = 'new_build_generator'
    new_projects_var_gen['operational_type'] = 'variable'
    new_projects_var_gen['variable_om_cost_per_mwh'] = 0
    new_projects_var_gen['balancing_type_project'] = 'day'
    # new_projects_var_gen['unit_size_mw'] = 100


    ## Concat existing and new projects
    # global projects_all
    projects_all = pd.DataFrame()
    projects_all = pd.concat([existing_projects, new_projects, new_projects_var_gen], sort = False)
    projects_all = projects_all.reset_index(drop=True)


    ### UPDATE VARIABLE OM COSTS ### NOTE: THIS IS ONLY IF YOU HAVE FUEL COSTS.
    projects_var_om_costs = pd.read_excel(projects_hr_var_costs_excel, sheet_name='project_var_om_costs')

    for prj in projects_var_om_costs['project'].to_list():
        projects_all.loc[projects_all['project'] == prj, 'variable_om_cost_per_mwh'] = projects_var_om_costs.loc[projects_var_om_costs['project'] == prj, 'variable_om_cost_per_mwh'].iloc[0]

    ## Rename the variable_om_cost_per_mwh to variable_cost_per_mwh, because the sql database column is named as that.
    # # TODO: Suggest adding the 'om' in that name to avoid confusion
    projects_all.rename(index=str, columns={'variable_om_cost_per_mwh': 'variable_cost_per_mwh'},
                                                     inplace=True)

    ### CREATE EXISTING PROJECT CAPACITIES ###
    # Drop generator unit level data and then drop duplicates to create the aggregate dataframe for existing projects capacities
    # Use copy so the new dataframe is not a view of the original df. Also avoids "SettingWithCopyWarning".
    # global existing_projects_capacities
    existing_projects_capacities = existing_project_units[['type', 'gen_capacity_mw_sum']].copy()
    existing_projects_capacities.rename(index = str, columns = {'type': 'project', 'gen_capacity_mw_sum': 'existing_capacity_mw'}, inplace = True)
    # Drop duplicates leaving only aggregated projects
    existing_projects_capacities = existing_projects_capacities.drop_duplicates()
    existing_projects_capacities = existing_projects_capacities.reset_index()


    ### CREATE MUST RUN GENERATORS CAPACITY FACTORS TIME SERIES ###

    variable_noncurtailable_cap_factors_daily = pd.read_csv(mustrun_gen_csv)
    # Limit data to 365 days i.e. no lookahead days
    variable_noncurtailable_cap_factors_daily = variable_noncurtailable_cap_factors_daily[: 365]
    # rename projects
    variable_noncurtailable_cap_factors_daily.rename(index = str, columns = {'HYDRO-ROR': 'Existing_Hydro_ROR', 'HYDRO-PONDAGE': 'Existing_Hydro_Pondage'}, inplace = True)

    projects_variable_noncurtailable_timeseries = pd.DataFrame()

    #for prj in projects_all.loc[projects_all['operational_type'] == 'variable_noncurtailable', 'project'].to_list():
    for prj in projects_all.loc[projects_all['project'].str.contains('ROR|Pondage'), 'project'].to_list():
        print(prj)
        projects_variable_noncurtailable_timeseries_prj = pd.DataFrame()
        variable_noncurtailable_cap_factors_daily_prj = variable_noncurtailable_cap_factors_daily[prj]
        variable_noncurtailable_cap_factors_hourly_prj = np.tile(variable_noncurtailable_cap_factors_daily_prj.values, 24)
        projects_variable_noncurtailable_timeseries_prj['cap_factor'] = variable_noncurtailable_cap_factors_hourly_prj
        projects_variable_noncurtailable_timeseries_prj['project'] = prj
        projects_variable_noncurtailable_timeseries_prj['day_of_year'] = day_array_8760
        projects_variable_noncurtailable_timeseries_prj['hour_of_year'] = hour_array_8760
        projects_variable_noncurtailable_timeseries_prj['hour'] = hour_of_day_array

        projects_variable_noncurtailable_timeseries = pd.concat([projects_variable_noncurtailable_timeseries, projects_variable_noncurtailable_timeseries_prj], axis = 0)

    # Write csv
    projects_variable_noncurtailable_timeseries.reset_index(drop = True)
    projects_variable_noncurtailable_timeseries.to_csv(projects_variable_noncurtailable_timeseries_csv, index=False)

    return (new_projects_var_gen, projects_all, existing_projects_capacities)

(new_projects_var_gen, projects_all, existing_projects_capacities) = create_projects_data()


def create_hydro_opchar():
    hydro_max_energy = pd.read_csv(hydro_max_energy_csv)
    hydro_min_gen = pd.read_csv(hydro_min_gen_csv)
    mustrun_gen = pd.read_csv(mustrun_gen_csv)

    projects_hydro_horizon_chars = pd.DataFrame()

    # This month array is for 1 day horizon with 365 days. For less number of days, use a csv.
    month_array = []
    for m in range(1, 12 + 1, 1):
        month_array = np.append(month_array, np.repeat(m, monthrange(base_year, m)[1])).astype(int)

    ## Hydro storage max energy and average power generation
    # Take only the first 365 days (avoid the lookahead days)
    hydro_max_energy = hydro_max_energy[:365]
    hydro_max_energy.rename(columns = {'HYDRO-STORAGE': 'daily_energy_mwh'}, inplace = True)
    hydro_max_energy['month'] = pd.DataFrame(month_array)
    # Take the monthly mean so that all average hourly generation for a month/horizon is the same
    mean = lambda x: x.mean()
    hydro_max_energy['average_daily_energy_mwh'] = hydro_max_energy.groupby(['month']).daily_energy_mwh.transform(mean)
    # Average power generation by dividing average daily energy by 24 hours
    hydro_max_energy['mwa'] = hydro_max_energy['average_daily_energy_mwh'] / 24
    hydro_max_energy['average_power_fraction'] = hydro_max_energy['mwa'] / existing_projects_capacities.loc[
        existing_projects_capacities['project'] == 'Existing_Hydro_Storage', 'existing_capacity_mw'].iloc[0]
    # drop duplicates and create average power generation by horizon dataframe
    projects_hydro_horizon_chars = hydro_max_energy[['month', 'average_power_fraction']].drop_duplicates().reset_index(drop=True)


    ## Hydro storage min generation
    # Take only the first 365 days (avoid the lookahead days)
    hydro_min_gen = hydro_min_gen[:365]
    hydro_min_gen.rename(columns = {'HYDRO-STORAGE': 'min_generation_perc'}, inplace = True)
    hydro_min_gen['month'] = pd.DataFrame(month_array)
    # Take the monthly mean so that all minimum generation for a month/horizon is the same
    hydro_min_gen['average_min_generation_perc'] = hydro_min_gen.groupby(['month']).min_generation_perc.transform(mean)
    # Minimum power generation by multiplying minimum generation capacity factor by rated capacity of hydro storage
    hydro_min_gen['min_power_mw'] = hydro_min_gen['average_min_generation_perc'] * existing_projects_capacities.loc[existing_projects_capacities['project'] == 'Existing_Hydro_Storage', 'existing_capacity_mw'].iloc[0]
    # drop duplicates and create average minimum power generation fraction by horizon dataframe
    projects_hydro_horizon_chars['min_power_fraction'] = hydro_min_gen[['average_min_generation_perc']].drop_duplicates().reset_index(drop=True)

    ## Hydro storage max generation as rated capacity
    # projects_hydro_horizon_chars['max_mw'] = existing_projects_capacities.loc[existing_projects_capacities['project'] == 'Existing_Hydro_Storage', 'existing_capacity_mw'].iloc[0]
    projects_hydro_horizon_chars['max_power_fraction'] = 1
    projects_hydro_horizon_chars['project'] = 'Existing_Hydro_Storage'

    ## Write out the csv
    projects_hydro_horizon_chars.to_csv(projects_hydro_horizon_chars_csv, index = False)


def calculate_rps(
        rps_scenario,
        load_for_rps,
        start_year,
        end_year
):
    """

    :param rps_scenario:
    :param load_for_rps:
    :param start_year:
    :param end_year:
    :return:
    """

    rps_targets_all = pd.read_excel(rps_targets_excel, sheet_name='rps_targets')

    rps_targets = OrderedDict()

    for rps_zone in rps_targets_all.columns.difference(
            ['year', 'rps_scenario']).to_list():  # Nice way of getting specific columns in a list
        print(rps_zone)
        rps_targets[rps_zone] = OrderedDict()
        for y in range(start_year, end_year + 1):
            rps_targets[rps_zone][y] = load_for_rps.loc[load_for_rps['year'] == y, rps_zone].iloc[0] * rps_targets_all.loc[(rps_targets_all['year'] == y) & (rps_targets_all['rps_scenario'] == rps_scenario), rps_zone].iloc[0]

    return rps_targets

def load_rps_targets():
    """

    :return:
    """

    # Specify scenarios to include, start and end year, annual load for each load zone that is also an rps zone
    rps_targets = pd.read_excel(rps_targets_excel, sheet_name = 'rps_targets')

    rps_scenarios_to_include = {1: '2030_rps_10_percent', 2: '2030_rps_30_percent', 3: '2030_rps_50_percent', 4: '2030_rps_70_percent'}
    st_year = 2018
    ed_year = 2030

    for rps_sc_id in rps_scenarios_to_include.keys():
        rps_scenario_name = rps_scenarios_to_include[rps_sc_id]
        rps_targets_for_scenario = calculate_rps(rps_scenario_name, load_annual, st_year, ed_year)

        # Insert into DB
        rps.insert_rps_targets(
            io=io, c=c2,
            rps_target_scenario_id=rps_sc_id,
            scenario_name=rps_scenario_name,
            scenario_description=rps_scenario_name.replace('_', ' '),
            zone_period_targets=rps_targets_for_scenario,
        )


def tuning():
    """

    :return:
    """
    print("tuning")
    # subscenarios = [
    #     (0, 'no_tuning', 'No tuning (tuning params = 0)', 0, 0, 0),
    #     (1, 'tune_carbon_imports', 'Tune carbon imports only', 10e-9, 0, 0),
    #     (2, 'tune_ramps', 'Tune ramps only', 0, 10e-9, 0),
    #     (3, 'tune_carbon_imports_and_ramps',
    #      'Tune carbon imports and ramps', 10e-9, 10e-9, 0),
    #     (4, 'tune_elcc_surface', 'Tune dynamic ELCC only', 0, 0, 10e-3),
    #     (5, 'tune_carbon_imports_and_elcc_surface',
    #      'Tune carbon imports and ELCC surface', 10e-9, 0, 10e-3),
    #     (6, 'tune_carbon_imports_ramps_and_elcc_surface',
    #      'Tune carbon imports, ramps, and ELCC surface', 10e-9, 10e-9, 10e-3)
    # ]

    subscenarios = [
        (0, 'no_tuning', 'No tuning (tuning params = 0)', 0, 0, 0),
    ]

    # Subscenarios
    for subscenario in subscenarios:
        c2.execute(
            """INSERT INTO subscenarios_tuning (tuning_scenario_id, name,
            description)
            VALUES ({}, '{}', '{}');""".format(
                subscenario[0], subscenario[1], subscenario[2]
            )
        )
    io.commit()

    # Data
    for subscenario in subscenarios:
        c2.execute(
            """INSERT INTO inputs_tuning
            (tuning_scenario_id, import_carbon_tuning_cost, 
            ramp_tuning_cost, dynamic_elcc_tuning_cost)
            VALUES ({}, {}, {}, {});""".format(
                subscenario[0], subscenario[3], subscenario[4], subscenario[5]
            )
        )
    io.commit()


# Default subscenario IDs for
defaults = {
    "of_fuels": 1,
    "of_multi_stage": 0,
    "of_transmission": 0,
    "of_transmission_hurdle_rates": 0,
    "of_simultaneous_flow_limits": 0,
    "of_lf_reserves_up": 0,
    "of_lf_reserves_down": 0,
    "of_regulation_up": 0,
    "of_regulation_down": 0,
    "of_frequency_response": 0,
    "of_spinning_reserves": 0,
    "of_rps": 1,
    "of_carbon_cap": 0,
    "of_track_carbon_imports": 0,
    "of_prm": 0,
    "of_local_capacity": 0,
    "of_elcc_surface": 0,
    "temporal_scenario_id": 1,
    "load_zone_scenario_id": 1,
    "lf_reserves_up_ba_scenario_id": 'NULL',
    "lf_reserves_down_ba_scenario_id": 'NULL',
    "regulation_up_ba_scenario_id": 'NULL',
    "regulation_down_ba_scenario_id": 'NULL',
    "frequency_response_ba_scenario_id": 'NULL',
    "spinning_reserves_ba_scenario_id": 'NULL',
    "rps_zone_scenario_id": 1,
    "carbon_cap_zone_scenario_id": 'NULL',
    "prm_zone_scenario_id": 'NULL',
    "local_capacity_zone_scenario_id": 'NULL',
    "project_portfolio_scenario_id": 1,
    "project_operational_chars_scenario_id": 1,
    "project_availability_scenario_id": 'NULL',
    "fuel_scenario_id": 1,
    "project_load_zone_scenario_id": 1,
    "project_lf_reserves_up_ba_scenario_id": 'NULL',
    "project_lf_reserves_down_ba_scenario_id": 'NULL',  # Allow renewables
    "project_regulation_up_ba_scenario_id": 'NULL',
    "project_regulation_down_ba_scenario_id": 'NULL',
    "project_frequency_response_ba_scenario_id": 'NULL',
    "project_spinning_reserves_ba_scenario_id": 'NULL',
    "project_rps_zone_scenario_id": 1,
    "project_carbon_cap_zone_scenario_id": 'NULL',
    "project_prm_zone_scenario_id": 'NULL',
    "project_elcc_chars_scenario_id": 'NULL',
    "prm_energy_only_scenario_id": 'NULL',
    "project_local_capacity_zone_scenario_id": 'NULL',
    "project_local_capacity_chars_scenario_id": 'NULL',
    "project_existing_capacity_scenario_id": 1,
    "project_existing_fixed_cost_scenario_id": 1,
    "fuel_price_scenario_id": 1,
    "project_new_cost_scenario_id": 1,
    "project_new_potential_scenario_id": 1,
    "transmission_portfolio_scenario_id": 'NULL',
    "transmission_load_zone_scenario_id": 'NULL',
    "transmission_existing_capacity_scenario_id": 'NULL',
    "transmission_operational_chars_scenario_id": 'NULL',
    "transmission_hurdle_rate_scenario_id": 'NULL',
    "transmission_carbon_cap_zone_scenario_id": 'NULL',
    "transmission_simultaneous_flow_limit_scenario_id": 'NULL',
    "transmission_simultaneous_flow_limit_line_group_scenario_id": 'NULL',
    "load_scenario_id": 1,
    "lf_reserves_up_scenario_id": 'NULL',
    "lf_reserves_down_scenario_id": 'NULL',
    "regulation_up_scenario_id": 'NULL',
    "regulation_down_scenario_id": 'NULL',
    "frequency_response_scenario_id": 'NULL',
    "spinning_reserves_scenario_id": 'NULL',
    "rps_target_scenario_id": 1,
    "carbon_cap_target_scenario_id": 'NULL',
    "prm_requirement_scenario_id": 'NULL',
    "elcc_surface_scenario_id": 'NULL',
    "local_capacity_requirement_scenario_id": 'NULL',
    "tuning_scenario_id": 0  # No tuning
}


def create_scenarios():
    """

    :return:
    """

    # Create 'Base_42MMT' scenario
    create_scenario.create_scenario(
        io=io, c=c2,
        scenario_name=main_sc,
        of_fuels=defaults["of_fuels"],
        of_multi_stage=defaults["of_multi_stage"],
        of_transmission=defaults["of_transmission"],
        of_transmission_hurdle_rates=defaults["of_transmission_hurdle_rates"],
        of_simultaneous_flow_limits=defaults["of_simultaneous_flow_limits"],
        of_lf_reserves_up=defaults["of_lf_reserves_up"],
        of_lf_reserves_down=defaults["of_lf_reserves_down"],
        of_regulation_up=defaults["of_regulation_up"],
        of_regulation_down=defaults["of_regulation_down"],
        of_frequency_response=defaults["of_frequency_response"],
        of_spinning_reserves=defaults["of_spinning_reserves"],
        of_rps=defaults["of_rps"],
        of_carbon_cap=defaults["of_carbon_cap"],
        of_track_carbon_imports=defaults["of_track_carbon_imports"],
        of_prm=defaults["of_prm"],
        of_local_capacity=defaults["of_local_capacity"],
        of_elcc_surface=defaults["of_elcc_surface"],
        temporal_scenario_id=defaults["temporal_scenario_id"],
        load_zone_scenario_id=defaults["load_zone_scenario_id"],
        lf_reserves_up_ba_scenario_id=
        defaults["lf_reserves_up_ba_scenario_id"],
        lf_reserves_down_ba_scenario_id=
        defaults["lf_reserves_down_ba_scenario_id"],
        regulation_up_ba_scenario_id=defaults["regulation_up_ba_scenario_id"],
        regulation_down_ba_scenario_id=
        defaults["regulation_down_ba_scenario_id"],
        frequency_response_ba_scenario_id=
        defaults["frequency_response_ba_scenario_id"],
        spinning_reserves_ba_scenario_id=
        defaults["spinning_reserves_ba_scenario_id"],
        rps_zone_scenario_id=defaults["rps_zone_scenario_id"],
        carbon_cap_zone_scenario_id=defaults["carbon_cap_zone_scenario_id"],
        prm_zone_scenario_id=defaults["prm_zone_scenario_id"],
        local_capacity_zone_scenario_id=defaults[
            "local_capacity_zone_scenario_id"],
        project_portfolio_scenario_id=
        defaults["project_portfolio_scenario_id"],
        project_operational_chars_scenario_id=
        defaults["project_operational_chars_scenario_id"],
        project_availability_scenario_id=
        defaults["project_availability_scenario_id"],
        fuel_scenario_id=defaults["fuel_scenario_id"],
        project_load_zone_scenario_id=
        defaults["project_load_zone_scenario_id"],
        project_lf_reserves_up_ba_scenario_id=
        defaults["project_lf_reserves_up_ba_scenario_id"],
        project_lf_reserves_down_ba_scenario_id=
        defaults["project_lf_reserves_down_ba_scenario_id"],
        project_regulation_up_ba_scenario_id=
        defaults["project_regulation_up_ba_scenario_id"],
        project_regulation_down_ba_scenario_id=
        defaults["project_regulation_down_ba_scenario_id"],
        project_frequency_response_ba_scenario_id=
        defaults["project_frequency_response_ba_scenario_id"],
        project_spinning_reserves_ba_scenario_id=
        defaults["project_spinning_reserves_ba_scenario_id"],
        project_rps_zone_scenario_id=defaults["project_rps_zone_scenario_id"],
        project_carbon_cap_zone_scenario_id=
        defaults["project_carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=defaults[
            "project_local_capacity_zone_scenario_id"],
        project_local_capacity_chars_scenario_id=
        defaults["project_local_capacity_chars_scenario_id"],
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=project_new_cost_scenario_id, # defaults["project_new_cost_scenario_id"]
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        transmission_portfolio_scenario_id=
        defaults["transmission_portfolio_scenario_id"],
        transmission_load_zone_scenario_id=
        defaults["transmission_load_zone_scenario_id"],
        transmission_existing_capacity_scenario_id=
        defaults["transmission_existing_capacity_scenario_id"],
        transmission_operational_chars_scenario_id=
        defaults["transmission_operational_chars_scenario_id"],
        transmission_hurdle_rate_scenario_id=
        defaults["transmission_hurdle_rate_scenario_id"],
        transmission_carbon_cap_zone_scenario_id=
        defaults["transmission_carbon_cap_zone_scenario_id"],
        transmission_simultaneous_flow_limit_scenario_id=
        defaults["transmission_simultaneous_flow_limit_scenario_id"],
        transmission_simultaneous_flow_limit_line_group_scenario_id=
        defaults[
            "transmission_simultaneous_flow_limit_line_group_scenario_id"],
        load_scenario_id=defaults["load_scenario_id"],
        lf_reserves_up_scenario_id=defaults["lf_reserves_up_scenario_id"],
        lf_reserves_down_scenario_id=defaults["lf_reserves_down_scenario_id"],
        regulation_up_scenario_id=defaults["regulation_up_scenario_id"],
        regulation_down_scenario_id=defaults["regulation_down_scenario_id"],
        frequency_response_scenario_id=
        defaults["frequency_response_scenario_id"],
        spinning_reserves_scenario_id=
        defaults["spinning_reserves_scenario_id"],
        rps_target_scenario_id=rps_target_scenario_id, #defaults["rps_target_scenario_id"]
        carbon_cap_target_scenario_id=defaults["carbon_cap_target_scenario_id"],
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"]
    )





load_temporal_data()
load_geography_load_zones()
load_geography_rps_zones()
create_projects_data()
create_hydro_opchar()
create_loads()
load_loads()
load_projects()
load_project_load_zones()
load_project_operational_chars()
load_project_hr_curves()
load_project_variable_profiles()
load_project_hydro_opchar()
load_project_portfolios()
load_project_capacities()
load_project_new_potentials()
load_project_fixed_costs()
load_project_new_costs()
load_fuels()
load_fuel_prices()
load_project_rps_zones()
load_rps_targets()

tuning()

### scenarios ###
# RPS0_VRElow_SThigh_CONVhigh
# RPS0_VREhigh_SThigh_CONVhigh
# RPS0_VRElow_STlow_CONVhigh
# RPS0_VREhigh_STlow_CONVhigh
# RPS30_VRElow_SThigh_CONVhigh
# RPS30_VREhigh_SThigh_CONVhigh
# RPS30_VRElow_STlow_CONVhigh
# RPS30_VREhigh_STlow_CONVhigh
# RPS50_VRElow_SThigh_CONVhigh
# RPS50_VREhigh_SThigh_CONVhigh
# RPS50_VRElow_STlow_CONVhigh
# RPS50_VREhigh_STlow_CONVhigh
# RPS70_VRElow_SThigh_CONVhigh
# RPS70_VREhigh_SThigh_CONVhigh
# RPS70_VRElow_STlow_CONVhigh
# RPS70_VREhigh_STlow_CONVhigh

main_scenarios = pd.read_csv(main_scenarios_csv)

for main_sc in main_scenarios['main_scenario_name'].to_list():
    if main_scenarios.loc[main_scenarios['main_scenario_name'] == main_sc, 'include'].iloc[0] == 1:
        project_new_cost_scenario_id = main_scenarios.loc[main_scenarios['main_scenario_name'] == main_sc, 'project_new_cost_scenario_id'].iloc[0]
        rps_target_scenario_id = main_scenarios.loc[main_scenarios['main_scenario_name'] == main_sc, 'rps_target_scenario_id'].iloc[0]
        print(main_sc)
        create_scenarios()


