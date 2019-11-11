
"""
Modules to import project data
"""

from collections import OrderedDict
import csv
import pandas as pd
import numpy as np
from db.utilities import project_list, project_zones, \
    project_operational_chars, fuels, \
    project_portfolios, project_existing_params, project_new_costs, \
    project_new_potentials


def load_projects():
    """
    Get a list of all projects
    Get projects from the projects_all dataframe.
    :return:
    """

    projects = list()

    projects = projects_all['project'].tolist()

    project_list.project_list(
        io=io, c=c2,
        projects=projects
    )


def load_project_load_zones():
    """
    Pull the load zone from the projects_all dataframe.
    :return:
    """

    # Load zones
    project_load_zones = dict()

    # Get projects and load_zone from dataframe and convert to dictionary with
    # projects as key
    project_load_zones = projects_all.loc[:, ['project', 'load_zone']].set_index(
        'project')['load_zone'].to_dict()

    project_zones.project_load_zones(
        io=io, c=c2,
        load_zone_scenario_id=1,  # Single load zone
        project_load_zone_scenario_id=1,
        scenario_name=
        'single load zone',
        scenario_description=
        "Default single load zone",
        project_load_zones=project_load_zones
    )


def load_project_operational_chars():
    """
    Operational characteristics of projects
    :return:
    """

    # Make subscenario and insert all projects into operational
    # characteristics table; we'll then update that table with the
    # operational characteristics each project needs
    project_operational_chars.make_scenario_and_insert_all_projects(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        scenario_name='default_op_char',
        scenario_description='Default operational characteristics.'
    )

    # ### Operational types ### #
    # Pull operational_type from projects_all

    project_op_types = dict()
    project_op_types = projects_all.loc[:, ['project', 'operational_type']].set_index(
        'project')['operational_type'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="operational_type",
        project_char=project_op_types
    )

    project_balancing_types = dict()
    project_balancing_types = projects_all.loc[:, ['project', 'balancing_type_project']].set_index(
        'project')['balancing_type_project'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="balancing_type_project",
        project_char=project_balancing_types
    )

    # ### Technologies ### #
    # Pull technology from projects_all
    project_tech = dict()
    project_tech = projects_all.loc[:, ['project', 'technology']].set_index(
        'project')['technology'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="technology",
        project_char=project_tech
    )

    # ### Variable O&M costs ### #
    # Eventually break this charge into fuel and VO&M charge assuming some heat rates.
    project_vom = dict()
    project_vom = projects_all.loc[:, ['project', 'variable_cost_per_mwh']].dropna().set_index(
        'project')['variable_cost_per_mwh'].to_dict()

    # Update 'conventional' projects
    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="variable_cost_per_mwh",
        project_char=project_vom
    )

    ###########################################################################
    # Example of SQL code to modify columns directly in the sql database.
    # c2.execute(
    #     """UPDATE inputs_project_operational_chars
    #     SET variable_cost_per_mwh = 600.0
    #     WHERE project LIKE 'CAISO_Shed_DR%'
    #     AND project_operational_chars_scenario_id = 1;"""
    # )
    # io.commit()
    #
    #
    # # Set variable cost for all remaining projects to 0
    # c2.execute(
    #     """UPDATE inputs_project_operational_chars
    #     SET variable_cost_per_mwh = 0.0
    #     WHERE variable_cost_per_mwh IS NULL
    #     AND project_operational_chars_scenario_id = 1;"""
    # )
    # io.commit()
    #
    # c2.execute(
    #     """UPDATE inputs_project_operational_chars
    #     SET fuel = 'CA_Natural_Gas'
    #     WHERE project LIKE 'CAISO_Peaker%_Planned'
    #     OR project = 'CAISO_CCGT1_Planned'
    #     OR project = 'CAISO_Reciprocating_Engine_Candidate'
    #     AND project_operational_chars_scenario_id = 1;"""
    # )
    # io.commit()
    ###########################################################################

    # ### Fuel ### #
    project_fuel = dict()
    project_fuel = projects_all.loc[:, ['project', 'fuel']].dropna().set_index(
        'project')['fuel'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="fuel",
        project_char=project_fuel
    )

    # ### Min Stable Level ### #
    project_min_stable_level = dict()
    project_min_stable_level = projects_all.loc[:, ['project', 'min_stable_level']].dropna().set_index(
        'project')['min_stable_level'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="min_stable_level",
        project_char=project_min_stable_level
    )

    # ### Heat rates ### #
    project_hr_curve_id = dict()
    project_hr_curve_id_df = projects_all.loc[:, ['project', 'heat_rate_curves_scenario_id']].dropna()
    project_hr_curve_id_df['heat_rate_curves_scenario_id'] = project_hr_curve_id_df['heat_rate_curves_scenario_id'].astype('Int64').astype('str') # Otherwise they are floats
    project_hr_curve_id = project_hr_curve_id_df.set_index(
        'project')['heat_rate_curves_scenario_id'].to_dict()

    project_hr_curve_id = projects_all.loc[:, ['project', 'heat_rate_curves_scenario_id']].dropna().set_index(
        'project')['heat_rate_curves_scenario_id'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="heat_rate_curves_scenario_id",
        project_char=project_hr_curve_id
    )

    # project_inc_hr = dict()
    # project_inc_hr = projects_all.loc[:, ['project', 'inc_heat_rate_mmbtu_per_mwh']].dropna().set_index(
    #     'project')['inc_heat_rate_mmbtu_per_mwh'].to_dict()
    #
    # project_operational_chars.update_project_opchar_column(
    #     io=io, c=c2,
    #     project_operational_chars_scenario_id=1,
    #     column="inc_heat_rate_mmbtu_per_mwh",
    #     project_char=project_inc_hr
    # )
    #
    # project_min_input = dict()
    # project_min_input = projects_all.loc[:, ['project', 'minimum_input_mmbtu_per_hr']].dropna().set_index(
    #     'project')['minimum_input_mmbtu_per_hr'].to_dict()
    #
    # project_operational_chars.update_project_opchar_column(
    #     io=io, c=c2,
    #     project_operational_chars_scenario_id=1,
    #     column="minimum_input_mmbtu_per_hr",
    #     project_char=project_min_input
    # )

    # ### Unit size ### #
    project_unit_size = projects_all.loc[:, ['project', 'unit_size_mw']].dropna().set_index(
        'project')['unit_size_mw'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="unit_size_mw",
        project_char=project_unit_size
    )

    # ### Startup and shutdown costs ### #
    project_startup_cost = dict()
    project_shutdown_cost = dict()
    project_startup_cost = projects_all.loc[:, ['project', 'startup_cost_per_mw']].dropna().set_index(
        'project')['startup_cost_per_mw'].to_dict()
    project_shutdown_cost = projects_all.loc[:, ['project', 'shutdown_cost_per_mw']].dropna().set_index(
        'project')['shutdown_cost_per_mw'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="startup_cost_per_mw",
        project_char=project_startup_cost
    )

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="shutdown_cost_per_mw",
        project_char=project_shutdown_cost
    )

    # ### Ramp rates ### #
    # project_off_on_plus_ramp_up = dict()
    # project_off_on_plus_ramp_down = dict()
    # project_when_on_ramp_up = dict()
    # project_when_on_ramp_down = dict()
    #
    # project_off_on_plus_ramp_up = projects_all.loc[:, ['project', 'startup_plus_ramp_up_rate']].dropna().set_index(
    #     'project')['startup_plus_ramp_up_rate'].to_dict()
    #
    # project_off_on_plus_ramp_down = projects_all.loc[:, ['project', 'shutdown_plus_ramp_down_rate']].dropna().set_index(
    #     'project')['shutdown_plus_ramp_down_rate'].to_dict()
    #
    # project_when_on_ramp_up = projects_all.loc[:, ['project', 'ramp_up_when_on_rate']].dropna().set_index(
    #     'project')['ramp_up_when_on_rate'].to_dict()
    #
    # project_when_on_ramp_down = projects_all.loc[:, ['project', 'ramp_down_when_on_rate']].dropna().set_index(
    #     'project')['ramp_down_when_on_rate'].to_dict()
    #
    # project_operational_chars.update_project_opchar_column(
    #     io=io, c=c2,
    #     project_operational_chars_scenario_id=1,
    #     column="startup_plus_ramp_up_rate",
    #     project_char=project_off_on_plus_ramp_up
    # )
    #
    # project_operational_chars.update_project_opchar_column(
    #     io=io, c=c2,
    #     project_operational_chars_scenario_id=1,
    #     column="shutdown_plus_ramp_down_rate",
    #     project_char=project_off_on_plus_ramp_down
    # )
    #
    # project_operational_chars.update_project_opchar_column(
    #     io=io, c=c2,
    #     project_operational_chars_scenario_id=1,
    #     column="ramp_up_when_on_rate",
    #     project_char=project_when_on_ramp_up
    # )
    # project_operational_chars.update_project_opchar_column(
    #     io=io, c=c2,
    #     project_operational_chars_scenario_id=1,
    #     column="ramp_down_when_on_rate",
    #     project_char=project_when_on_ramp_down
    # )

    # Not used
    # project_operational_chars.update_project_opchar_column(
    #     io=io, c=c2,
    #     project_operational_chars_scenario_id=1,
    #     column="lf_reserves_down_derate",
    #     project_char=var_proj_lf_down_derate
    #     )

    # ### Min up and down time ### #
    project_min_up_time = dict()
    project_min_down_time = dict()
    project_min_up_time = projects_all.loc[:, ['project', 'min_up_time_hours']].dropna().set_index(
        'project')['min_up_time_hours'].to_dict()
    project_min_down_time = projects_all.loc[:, ['project', 'min_down_time_hours']].dropna().set_index(
        'project')['min_down_time_hours'].to_dict()

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="min_up_time_hours",
        project_char=project_min_up_time
    )

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="min_down_time_hours",
        project_char=project_min_down_time
    )

    ### Charging and discharging efficiency ### #
    proj_eff_charging = dict()
    proj_eff_discharging = dict()
    proj_eff_charging = projects_all.loc[:, ['project', 'charging_efficiency']].dropna().set_index(
        'project')['charging_efficiency'].to_dict()
    proj_eff_discharging = projects_all.loc[:, ['project', 'discharging_efficiency']].dropna().set_index(
        'project')['discharging_efficiency'].to_dict()

    # proj_eff_charging = {
    #     'Existing_Pumped_Storage': 0.81**(1/2.0),
    #     'Storage_Generic': 0.85**(1/2.0),
    #     'New_Pumped_Storage': 0.81**(1/2.0),
    #     'New_Flow_Battery': 0.70**(1/2.0),
    #     'New_Li_Battery': 0.85**(1/2.0)
    # }
    # proj_eff_discharging = proj_eff_charging

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="charging_efficiency",
        project_char=proj_eff_charging
    )

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="discharging_efficiency",
        project_char=proj_eff_discharging
    )

    # ### Minimum duration ### #
    # LCR-eligible batteries must have at least 4 hours of duration
    proj_min_dur = dict()
    proj_min_dur = projects_all.loc[:, ['project', 'minimum_duration_hours']].dropna().set_index(
        'project')['minimum_duration_hours'].to_dict()

    # proj_min_dur = {
    #     'Existing_Pumped_Storage': 12,
    #     'Storage_Generic': 1,
    #     'New_Pumped_Storage': 12,
    #     'New_Flow_Battery': 1,
    #     'New_Li_Battery': 1
    # }

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="minimum_duration_hours",
        project_char=proj_min_dur
    )

    # ### Variable generator profiles ### #
    project_operational_chars. \
        update_project_opchar_variable_gen_profile_scenario_id(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        variable_generator_profile_scenario_id=1
    )

    # ### Hydro operational characteristics ### #
    project_operational_chars.update_project_opchar_hydro_opchar_scenario_id(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        hydro_operational_chars_scenario_id=1
    )


def load_project_hr_curves():
    """
    Heat rate curves of projects
    :return:
    """
    projects_heatrate_curves = pd.read_excel(projects_hr_var_costs_excel, sheet_name='project_heat_rate_curves')
    heat_rate_curves_scenarios = pd.read_excel(projects_hr_var_costs_excel, sheet_name='heat_rate_curves_scenarios')

    for prj in projects_heatrate_curves['project'].to_list():
        projects_heatrate_curves.loc[projects_heatrate_curves['project'] == prj, 'load_point_mw'] = \
            projects_all.loc[projects_all['project'] == prj, 'unit_size_mw'].iloc[0] * \
                                                    projects_heatrate_curves.loc[projects_heatrate_curves['project'] == prj, 'load_point_percentage']

    proj_heat_rate_curves_dict = OrderedDict()
    proj_heat_rate_curves_names_dict = OrderedDict()

    for prj in projects_heatrate_curves['project'].to_list():
        proj_heat_rate_curves_dict[prj] = OrderedDict()
        proj_heat_rate_curves_names_dict[prj] = OrderedDict()

        for scenario in projects_heatrate_curves['heat_rate_curves_scenario_id'].unique():
            proj_heat_rate_curves_dict[prj][scenario] = OrderedDict()

            hr_curve_points_project = projects_heatrate_curves.loc[(projects_heatrate_curves['project'] == prj) & (projects_heatrate_curves['heat_rate_curves_scenario_id'] == scenario), 'heat_rate_curve_point'].to_list()

            for hr_curve_point in hr_curve_points_project:
                proj_heat_rate_curves_dict[prj][scenario][hr_curve_point] = OrderedDict()


                proj_heat_rate_curves_dict[prj][scenario][hr_curve_point] = (projects_heatrate_curves.loc[
                                                                (projects_heatrate_curves['project'] == prj) &
                                                                                         (projects_heatrate_curves['heat_rate_curves_scenario_id'] == scenario), 'load_point_mw'].iloc[0],
                                                            projects_heatrate_curves.loc[
                                                                (projects_heatrate_curves['project'] == prj) &
                                                                (projects_heatrate_curves[
                                                                     'heat_rate_curves_scenario_id'] == scenario), 'average_heat_rate_mmbtu_per_mwh'].iloc[0]
                                                            )

            proj_heat_rate_curves_names_dict[prj][scenario] = (heat_rate_curves_scenarios.loc[heat_rate_curves_scenarios['heat_rate_curves_scenario_id'] == scenario, 'heat_rate_curves_scenario'].iloc[0],
                                                               heat_rate_curves_scenarios.loc[
                                                                   heat_rate_curves_scenarios[
                                                                       'heat_rate_curves_scenario_id'] == scenario, 'heat_rate_curves_scenario_description'].iloc[0])

    project_operational_chars.update_project_hr_curves(
        io=io, c=c2,
        proj_opchar_names=proj_heat_rate_curves_names_dict,
        proj_hr_chars=proj_heat_rate_curves_dict
    )


def load_project_variable_profiles():
    """
    Profiles for 'variable' generators
    :return:
    """

    new_projects_wind_timeseries = pd.read_csv(new_projects_wind_timeseries_csv)
    new_projects_solar_timeseries = pd.read_csv(new_projects_solar_timeseries_csv)
    projects_variable_noncurtailable_timeseries = pd.read_csv(projects_variable_noncurtailable_timeseries_csv)
    new_projects_var_gen_timeseries = pd.concat([new_projects_wind_timeseries, new_projects_solar_timeseries, projects_variable_noncurtailable_timeseries], axis = 0, ignore_index = True)


    #projects_w_tmp = projects_all.loc[projects_all['technology'].str.contains('Wind|Solar')]
    projects_w_tmp = projects_all.loc[projects_all['operational_type'].str.contains('variable')]
    proj_tmp_profiles_all_dict = OrderedDict()
    proj_profile_names_dict = OrderedDict()

    for prj in projects_w_tmp['project'].to_list():
        print(prj)
        proj_tmp_profile = pd.DataFrame()
        # prj_orig = prj.replace("_", " ")
        #proj_tmp_profile = new_projects_var_gen_timeseries.loc[new_projects_var_gen_timeseries['project'] == prj]

        for y in year_month_days_selected['year'].unique():
            selected_days = year_month_days_selected.loc[year_month_days_selected['year'] == y, 'day_of_year']
            proj_tmp_profile_y = new_projects_var_gen_timeseries.loc[((new_projects_var_gen_timeseries['project'] == prj) &
                                                                                             (new_projects_var_gen_timeseries['day_of_year'].isin(selected_days)))]
            proj_tmp_profile_y['year'] = y
            proj_tmp_profile = pd.concat([proj_tmp_profile, proj_tmp_profile_y], axis = 0, ignore_index=True)




        # proj_tmp_profile = proj_tmp_profile.loc[pd.to_numeric(proj_tmp_profile['YEAR']) == base_year,] # Select data based on dates
        # proj_tmp_profile['cap_factor'] = proj_tmp_profile['VALUE'] # Either rename the VALUE column or create a new one.
        # Add timepoint and zone
        proj_tmp_profile = pd.concat([proj_tmp_profile.reset_index(), period_horizon_timepoints['tmp']], axis=1)
        proj_tmp_profile = proj_tmp_profile.set_index(['tmp'])
        # Delete all other columns
        proj_tmp_profile = proj_tmp_profile.filter(['cap_factor'])
        # Convert cap_factor to dictionary
        proj_tmp_profile_dict = OrderedDict()
        proj_tmp_profile_dict = proj_tmp_profile['cap_factor'].to_dict()
        ## TODO: Modify this code to accept multiple stages and scenarios from inputs. Currently hardcoded.
        # Add to stage dictionary
        proj_tmp_profile_stage_dict = OrderedDict()
        proj_tmp_profile_stage_dict[1] = proj_tmp_profile_dict
        # Add to scenario id dictionary
        proj_tmp_profile_stage_scenario_dict = OrderedDict()
        proj_tmp_profile_stage_scenario_dict[1] = proj_tmp_profile_stage_dict

        # Add vre profile to project as key in a nested dictionary of all vre profiles
        proj_tmp_profiles_all_dict[prj] = proj_tmp_profile_stage_scenario_dict

        # Add project to project names list along with scenario id and description
        proj_profile_names_scenario_dict = OrderedDict()
        proj_profile_names_scenario_dict[1] = ('default_cap_factors', 'Default profiles for variable generators')
        proj_profile_names_dict[prj] = proj_profile_names_scenario_dict

    project_operational_chars.update_project_variable_profiles(
        io=io, c=c2,
        proj_profile_names=proj_profile_names_dict,
        proj_tmp_profiles=proj_tmp_profiles_all_dict
    )


def load_project_hydro_opchar():
    """
    Energy budget, min, and max by horizon for hydro projects
    :return:
    """
    projects_hydro_horizon_chars = pd.read_csv(projects_hydro_horizon_chars_csv)

    projects_hydro_horizon_chars_dict = OrderedDict()

    projects_hydro_horizon_chars_prj = pd.DataFrame()

    proj_opchar_names_dict = OrderedDict()

    for prj in projects_hydro_horizon_chars['project'].unique():
        print(prj)
        for p in period_horizon['period'].unique():
            projects_hydro_horizon_chars_prj_p = pd.concat([projects_hydro_horizon_chars.loc[projects_hydro_horizon_chars['project'] == prj].set_index('month'), period_horizon.loc[period_horizon['period'] == p].set_index('horizonstep')], axis=1)
            projects_hydro_horizon_chars_prj = projects_hydro_horizon_chars_prj.append(projects_hydro_horizon_chars_prj_p, ignore_index=True)

        projects_hydro_horizon_chars_prj = projects_hydro_horizon_chars_prj[['horizon', 'period', 'average_power_fraction', 'max_power_fraction', 'min_power_fraction']].set_index(['horizon'])
        projects_hydro_horizon_chars_prj_dict = OrderedDict()
        # TODO: Creating only a single scenario here. Modify code to accept multiple scenarios from input file
        hydro_operational_chars_scenario_id = 1
        projects_hydro_horizon_chars_dict[prj] = OrderedDict()
        # Add scenario
        projects_hydro_horizon_chars_dict[prj][hydro_operational_chars_scenario_id] = OrderedDict()
        # Add the "day" as the balancing horizon
        projects_hydro_horizon_chars_dict[prj][hydro_operational_chars_scenario_id]["day"] = OrderedDict()
        # Add data
        projects_hydro_horizon_chars_dict[prj][hydro_operational_chars_scenario_id]["day"] = projects_hydro_horizon_chars_prj.to_dict(orient='index')

        # Add project to project names list along with scenario id and description
        proj_opchar_names_scenario_dict = OrderedDict()
        proj_opchar_names_scenario_dict[hydro_operational_chars_scenario_id] = ('default_hydro_budget_min_max', 'Default hydro horizon energy budgets, minimum, maximum')
        proj_opchar_names_dict[prj] = proj_opchar_names_scenario_dict

    # Insert data
    project_operational_chars.update_project_hydro_opchar(
        io=io, c=c2,
        proj_opchar_names=proj_opchar_names_dict,
        proj_horizon_chars=projects_hydro_horizon_chars_dict
    )


def load_project_portfolios():
    """
    Project portfolios
    :return:
    """

    portfolios = OrderedDict()
    scenario_list = ['base_portfolio']

    for sc in scenario_list:
        portfolios[sc] = OrderedDict()
        # for proj in projects_all['project'].to_list():
        #     portfolios[sc][proj] = projects_all.loc[projects_all['project'] == proj, 'capacity_type']
        #
        portfolios[sc] = projects_all.loc[:, ['project', 'capacity_type']].dropna().set_index(
            'project')['capacity_type'].to_dict()

    scenario_id = 1
    for portfolio in portfolios.keys():
        project_portfolios.update_project_portfolios(
            io=io, c=c2,
            project_portfolio_scenario_id=scenario_id,
            scenario_name=portfolio,
            scenario_description=portfolio.replace("_", " "),
            project_cap_types=portfolios[portfolio]
        )
        scenario_id += 1


def load_project_capacities():
    """

    :return:
    """

    proj_capacities = OrderedDict()
    scenario_list = ['existing_capacity_no_retirements']

    for sc in scenario_list:
        proj_capacities[sc] = OrderedDict()
        for proj in existing_projects_capacities['project'].to_list():
            proj_capacities[sc][proj] = OrderedDict()
            # TODO: Modify this code to import from csv and include retirements and additional scenarios
            for p in periods_df['period'].to_list():
                proj_capacities[sc][proj][p] = (existing_projects_capacities.loc[existing_projects_capacities['project'] == proj, 'existing_capacity_mw'].iloc[0], None)


    # Insert data
    scenario_id = 1
    for scenario in proj_capacities.keys():
        project_existing_params.update_project_capacities(
            io=io, c=c2,
            project_existing_capacity_scenario_id=scenario_id,
            scenario_name=scenario,
            scenario_description=scenario,
            project_capacities=proj_capacities[scenario]
        )
        scenario_id += 1


def load_project_new_potentials():
    """
    Max potentials for candidate projects
    :return:
    """

    proj_new_potentials = OrderedDict()
    scenario_list = ['base_potentials']

    for sc in scenario_list:
        proj_new_potentials[sc] = OrderedDict()
        for proj in projects_all.loc[projects_all['technology'].str.contains('Wind|Solar'), 'project'].to_list():
            proj_new_potentials[sc][proj] = OrderedDict()
            for p in periods_df['period'].to_list():
                proj_new_potentials[sc][proj][p] = (None, None, projects_all.loc[projects_all['project'] == proj, 'maximum_cumulative_new_build_mw'].iloc[0], None)


    # Insert data
    scenario_id = 1
    for scenario in proj_new_potentials.keys():
        project_new_potentials.update_project_potentials(
            io=io, c=c2,
            project_new_potential_scenario_id=scenario_id,
            scenario_name=scenario,
            scenario_description=scenario.replace("_", ""),
            project_period_potentials=proj_new_potentials[scenario]
        )

        scenario_id += 1

def load_project_fixed_costs():
    """

    :return:
    """
    project_fixed_costs = OrderedDict()

    for prj in existing_projects_capacities['project'].to_list():
        project_fixed_costs[prj] = OrderedDict()
        for p in np.arange(2018,2031):
            project_fixed_costs[prj][p] = (0, None)

    project_existing_params.update_project_fixed_costs(
        io=io, c=c2,
        project_existing_fixed_cost_scenario_id=1,
        scenario_name="default_fixed_costs",
        scenario_description='Default fixed costs',
        project_fixed_costs=project_fixed_costs
    )

def load_project_new_costs():
    """
    Costs for candidate projects
    :return:
    """

    new_battery_costs_all = pd.read_excel(new_projects_excel, sheet_name='new_battery_costs')
    new_vre_costs_all = pd.read_excel(new_projects_excel, sheet_name='new_vre_costs')
    new_conventional_costs_all = pd.read_excel(new_projects_excel, sheet_name='new_conventional_costs')
    new_costs_scenarios = pd.read_excel(new_projects_excel, sheet_name='new_costs_scenarios')


    # sc = 1

    for sc_id in new_costs_scenarios['new_costs_scenario_id'].to_list():
        new_battery_costs_scenario = new_costs_scenarios.loc[
            new_costs_scenarios['new_costs_scenario_id'] == sc_id, 'new_battery_costs_scenario'].iloc[0]
        new_vre_costs_scenario = new_costs_scenarios.loc[
            new_costs_scenarios['new_costs_scenario_id'] == sc_id, 'new_vre_costs_scenario'].iloc[0]
        new_conventional_costs_scenario = new_costs_scenarios.loc[
            new_costs_scenarios['new_costs_scenario_id'] == sc_id, 'new_conventional_costs_scenario'].iloc[0]
        new_costs_scenario = new_costs_scenarios.loc[
            new_costs_scenarios['new_costs_scenario_id'] == sc_id, 'new_costs_scenario'].iloc[0]
        new_costs_scenario_description = new_costs_scenarios.loc[
            new_costs_scenarios['new_costs_scenario_id'] == sc_id, 'new_costs_scenario_description'].iloc[0]
        print(new_costs_scenario)

        new_battery_costs = new_battery_costs_all.loc[(new_battery_costs_all['scenario'] == new_battery_costs_scenario), :]
        new_conventional_costs = new_conventional_costs_all.loc[
                            (new_conventional_costs_all['scenario'] == new_conventional_costs_scenario), :]


        # Assign new technology costs to specific projects
        new_vre_project_costs = pd.DataFrame()

        # for scenario_vre in new_vre_costs['scenario'].unique():

        for proj in projects_all.loc[projects_all['project'].str.contains('wind'), 'project'].to_list():
            new_vre_project_costs_prj = new_vre_costs_all.loc[(new_vre_costs_all['project'] == 'wind') & (new_vre_costs_all['scenario'] == new_vre_costs_scenario), :]
            new_vre_project_costs_prj.loc[:, 'project'] = proj
            new_vre_project_costs = new_vre_project_costs.append(new_vre_project_costs_prj, ignore_index=True)

        for proj in projects_all.loc[projects_all['project'].str.contains('solarPV'), 'project'].to_list():
            new_vre_project_costs_prj = new_vre_costs_all.loc[(new_vre_costs_all['project'] == 'solarPV') & (new_vre_costs_all['scenario'] == new_vre_costs_scenario), :]
            new_vre_project_costs_prj.loc[:, 'project'] = proj
            new_vre_project_costs = new_vre_project_costs.append(new_vre_project_costs_prj, ignore_index=True)


        new_project_costs = pd.concat([new_conventional_costs, new_battery_costs, new_vre_project_costs], ignore_index=True)
        new_project_costs = new_project_costs.loc[:, ['project', 'period', 'lifetime_yrs', 'annualized_real_cost_per_kw_yr', 'annualized_real_cost_per_kwh_yr']]

        project_period_lifetimes_costs = OrderedDict()

        # TODO: Creating only a single scenario here. Modify code to accept multiple scenarios from input file
        project_new_cost_scenario_id = 1
        for prj in new_project_costs['project'].unique():
            #project_period_lifetimes_costs_prj = OrderedDict()
            new_project_costs_prj = new_project_costs.loc[
                new_project_costs['project'] == prj, ['period', 'lifetime_yrs', 'annualized_real_cost_per_kw_yr',
                                                      'annualized_real_cost_per_kwh_yr']]
            project_period_lifetimes_costs[prj] = OrderedDict()
            for p in new_project_costs_prj['period'].to_list():
                project_period_lifetimes_costs[prj][p] = (new_project_costs_prj.loc[new_project_costs_prj['period'] == p, 'lifetime_yrs'].iloc[0],
                                                          new_project_costs_prj.loc[
                                                              new_project_costs_prj['period'] == p, 'annualized_real_cost_per_kw_yr'].iloc[0],
                                                          new_project_costs_prj.loc[
                                                              new_project_costs_prj['period'] == p, 'annualized_real_cost_per_kwh_yr'].iloc[0])

        project_new_costs.update_project_new_costs(
            io=io, c=c2,
            project_new_cost_scenario_id=sc_id,
            scenario_name=new_costs_scenario,
            scenario_description=new_costs_scenario_description,
            project_period_lifetimes_costs=project_period_lifetimes_costs
        )



def load_fuels():
    """
    Fuels and CO2 intensity
    :return:
    """
    fuel_chars = dict()

    with open(fuels_costs_csv, "r") as f:
        rows_list = list(csv.reader(f))
        # CO2 intensity
        for row in range(2-1, 7):  # skip Conventional_DR, not fuel in GP
            fuel_chars[rows_list[row][2-1]] = rows_list[row][3-1]

    # Insert data
    fuels.update_fuels(
        io=io, c=c2,
        fuel_scenario_id=1,
        scenario_name="default_fuels_and_co2_intensities",
        scenario_description="Default fuels and their CO2 intensities.",
        fuel_chars=fuel_chars
    )


def load_fuel_prices():
    """
    Fuel prices
    1 base price scenario1 X 1 carbon cost scenarios
    :return:
    """

    # modify this code if more fuel and carbon price scenarios need to be added. Change the csv input.
    scenario_fuel_month_prices = OrderedDict()

    with open(fuels_costs_csv, "r") as f:
        rows_list = list(csv.reader(f))

        # Get the monthly price shape
        fuel_month_price_shape = OrderedDict()
        for row in range(2 - 1, 7): # Rows with data on all fuels
            fuel_month_price_shape[rows_list[row][6 - 1]] = OrderedDict() # Column with fuel
            for column in range(9 - 1, 20):
                fuel_month_price_shape[
                    rows_list[row][6 - 1] # Column with fuel
                ][
                    int(float(rows_list[1 - 1][column])) # Row with month number
                ] = \
                    float(rows_list[row][column])

        # Get the base fuel prices by scenario
        scenario_fuel_year_base_price = OrderedDict()
        for base_scenario in ["Zero"]: # Add all the fuel price scenarios
            scenario_fuel_year_base_price[base_scenario] = OrderedDict()
        # Zero
        for row in range(12 - 1, 17): # Rows with fuel price data for all fuels
            scenario_fuel_year_base_price["Zero"][rows_list[row][3 - 1]] = OrderedDict() # Column with fuel names
            for column in range(8 - 1, 23): # Columns for all the years
                scenario_fuel_year_base_price[
                    "Zero"
                ][
                    rows_list[row][3 - 1]
                ][
                    int(float(rows_list[9 - 1][column])) # Row with all the years
                ] = \
                    float(rows_list[row][column])

        # Get the carbon adder by carbon price scenario
        scenario_year_carbon_price = OrderedDict()
        for row in range(21 - 1, 21):
            scenario_year_carbon_price[rows_list[row][3 - 1]] = OrderedDict()
            for column in range(8 - 1, 23):
                scenario_year_carbon_price[
                    rows_list[row][3 - 1]
                ][
                    int(float(rows_list[19 - 1][column])) # Row with years
                ] = \
                    float(rows_list[row][column])

        # Get the fuel CO2 intensity for prices
        fuel_co2_intensities = dict()
        for row in range(2-1, 7):  # Rows with fuels and their CO2 intensities. Can skip a fuel if you
            fuel_co2_intensities[rows_list[row][2-1]] = \
                float(rows_list[row][4-1]) # First Column with fuel names and then Column with fuel CO2 intensities

        # Make fuel price scenarios dictionary
        for base_price_scenario in scenario_fuel_year_base_price.keys():
            for carbon_price_scenario in \
                    scenario_year_carbon_price.keys():
                scenario_name = \
                    base_price_scenario + "_" + carbon_price_scenario
                scenario_fuel_month_prices[scenario_name] = \
                    OrderedDict()
                for fuel in scenario_fuel_year_base_price[base_price_scenario]:
                    scenario_fuel_month_prices[scenario_name][fuel] = \
                        OrderedDict()
                    fuel_co2_intensity = fuel_co2_intensities[fuel]
                    for year in scenario_fuel_year_base_price[
                        base_price_scenario
                    ][
                        fuel
                    ].keys():
                        scenario_fuel_month_prices[scenario_name][fuel][year]\
                            = OrderedDict()
                        carbon_price = scenario_year_carbon_price[
                            carbon_price_scenario][year]
                        for month in fuel_month_price_shape[fuel].keys():
                            base_monthly_price = scenario_fuel_year_base_price[
                                base_price_scenario
                            ][
                                fuel
                            ][
                                year
                            ] * float(fuel_month_price_shape[fuel][month])
                            carbon_adjusted_price = \
                                base_monthly_price + \
                                carbon_price * fuel_co2_intensity

                            scenario_fuel_month_prices[scenario_name][fuel][
                                year][month] = carbon_adjusted_price

    # Insert scenarios
    scenario_id = 1  # scenario ID to start with
    for base_price_scenario in ["Zero"]:
        for carbon_price_scenario in ["Zero"]:
            dict_scenario_name = base_price_scenario + "_" \
                                 + carbon_price_scenario
            fuels.update_fuel_prices(
                io=io, c=c2,
                fuel_price_scenario_id=scenario_id,
                scenario_name=base_price_scenario + "_Base_Fuel_Price_w_" +
                              carbon_price_scenario + "_Carbon_Adder",
                scenario_description=
                base_price_scenario + "_Base_Fuel_Price_w_" +
                carbon_price_scenario + "_Carbon_Adder",
                fuel_month_prices=
                scenario_fuel_month_prices[dict_scenario_name]
            )

            scenario_id += 1


def load_project_rps_zones():
    """
    Project RPS zones
    :return:
    """
    project_rps_zones = dict()

    for prj in projects_all.loc[projects_all['project'].str.contains('solarPV|wind'), 'project'].to_list():
        project_rps_zones[prj] = 'India'


    # Insert data
    project_zones.project_policy_zones(
        io=io, c=c2,
        policy_zone_scenario_id=1,
        project_policy_zone_scenario_id=1,
        scenario_name="wind_solar_only",
        scenario_description="wind and solar only",
        project_zones=project_rps_zones,
        policy_type="rps"
    )


if __name__ == "__main__":
    pass

# load_projects(projects_all = projects_all)
# load_project_load_zones(projects_all)
# load_project_operational_chars(projects_all)
# load_project_hr_curves()
# load_project_variable_profiles()
# load_project_hydro_opchar()
# load_project_portfolios()
# load_project_capacities()
# load_project_new_potentials()
# load_project_fixed_costs()
# load_project_new_costs()
# load_fuels()
# load_fuel_prices()
# load_project_rps_zones()