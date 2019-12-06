#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

"""
Port the IRP RESOLVE data to GridPath
"""

from collections import OrderedDict
import csv
import os.path
import sqlite3
import warnings

# Data-import modules

from db.utilities import temporal, geography, project_list, project_zones, \
    project_operational_chars, project_availability, fuels, \
    project_portfolios, project_existing_params, project_new_costs, \
    project_new_potentials, project_prm, transmission_portfolios, \
    transmission_zones, transmission_capacities, simultaneous_flow_groups, \
    simultaneous_flows, transmission_hurdle_rates, carbon_cap, system_load, \
    system_reserves, system_prm, rps, scenario


resolve = sqlite3.connect(
    os.path.join(os.getcwd(), 'resolve.db')
)

io = sqlite3.connect(
    os.path.join(os.getcwd(), 'io.db')
)

c1 = resolve.cursor()
c2 = io.cursor()


def load_temporal_data():
    """
    Add RESOLVE timepoints/days into database
    Investment periods are 2018, 2022, 2026, and 2030
    Discount factors and number of years represented are the same as in RESOLVE
    (the GridPath discount_factor is the RESOLVE weighted divided by the
    years in period)
    Horizon boundary is circular
    :return:
    """

    # Make a dictionary with the discount factor and number of years
    # represented for each investment period
    periods = dict(dict())
    with open(os.path.join("cpuc_irp_data", "csvs", "Lists.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(4 - 1, 7):
            for column in range(32 - 1, 34):
                periods[
                    int(float(rows_list[row][32 - 1]))
                ] = {}
                periods[
                    int(float(rows_list[row][32 - 1]))
                ]["discount_factor"] = float(rows_list[row][33 - 1]) / \
                          float(rows_list[row][34 - 1])
                periods[
                    int(float(rows_list[row][32 - 1]))
                ]["number_years_represented"] = float(rows_list[row][34 - 1])

    # Make a dictionary of horizons with their weights and months
    # We'll also keep these in the RESOLVE database

    c1.execute(
        """DROP TABLE IF EXISTS resolve_days;"""
    )
    c1.execute(
        """CREATE TABLE resolve_days(
        day_id INTEGER PRIMARY KEY,
        day_weight FLOAT,
        month_of_year INTEGER,
        hydro_year INTEGER
        );"""
    )

    # Make dictionary for GridPath database
    horizon_weights_and_months = dict()

    for row in range(4 - 1, 40):
        # Save in RESOLVE database
        c1.execute(
            """INSERT INTO resolve_days (day_id, day_weight,
            month_of_year, hydro_year) VALUES ({}, {}, {}, {});""".format(
                int(float(rows_list[row][15 - 1])), rows_list[row][16 - 1],
                int(float(rows_list[row][17 - 1])),
                int(float(rows_list[row][18 - 1]))
            )
        )

        # Populate dictionary for GridPath database
        for column in range(15 - 1, 17):
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ] = {}
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ]["weight"] = float(rows_list[row][16 - 1])
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ]["month"] = float(rows_list[row][17 - 1])

    resolve.commit()

    # Subproblems and stages (one subproblem, one stage)
    subproblems = [1]
    subproblem_stages = {sid: [(1, "single stage")] for sid in subproblems}

    # Timepoints
    subproblem_stage_timepoints = dict()
    for subproblem_id in subproblem_stages.keys():
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
            scenario_name="default_4_periods_37_days_24_hours",
            scenario_description="2018, 2022, 2026, 2030; 37 RESOLVE days, "
                                 "24 hours each",
            periods=periods,
            subproblems=[1],
            subproblem_stages={1: [(1, "single stage")]},
            subproblem_stage_timepoints=subproblem_stage_timepoints,
            subproblem_horizons=subproblem_horizons,
            subproblem_stage_timepoint_horizons=subproblem_stage_timepoint_horizons
    )


def load_temporal_data_2030_only():
    """
    Add RESOLVE timepoints/days into database
    Investment periods are 2018, 2022, 2026, and 2030
    Discount factors and number of years represented are the same as in RESOLVE
    (the GridPath discount_factor is the RESOLVE weighted divided by the
    years in period)
    Horizon boundary is circular
    :return:
    """

    # Make a dictionary with the discount factor and number of years
    # represented for each investment period
    periods = dict(dict())
    with open(os.path.join("cpuc_irp_data", "csvs", "Lists.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(4 - 1, 7):
            for column in range(32 - 1, 34):
                periods[
                    int(float(rows_list[row][32 - 1]))
                ] = {}
                periods[
                    int(float(rows_list[row][32 - 1]))
                ]["discount_factor"] = float(rows_list[row][33 - 1]) / \
                          float(rows_list[row][34 - 1])
                periods[
                    int(float(rows_list[row][32 - 1]))
                ]["number_years_represented"] = float(rows_list[row][34 - 1])

    # Make a dictionary of horizons with their weights and months
    # We'll also keep these in the RESOLVE database

    c1.execute(
        """DROP TABLE IF EXISTS resolve_days;"""
    )
    c1.execute(
        """CREATE TABLE resolve_days(
        day_id INTEGER PRIMARY KEY,
        day_weight FLOAT,
        month_of_year INTEGER,
        hydro_year INTEGER
        );"""
    )

    # Make dictionary for GridPath database
    horizon_weights_and_months = dict()

    for row in range(4 - 1, 40):
        # Save in RESOLVE database
        c1.execute(
            """INSERT INTO resolve_days (day_id, day_weight,
            month_of_year, hydro_year) VALUES ({}, {}, {}, {});""".format(
                int(float(rows_list[row][15 - 1])), rows_list[row][16 - 1],
                int(float(rows_list[row][17 - 1])),
                int(float(rows_list[row][18 - 1]))
            )
        )

        # Populate dictionary for GridPath database
        for column in range(15 - 1, 17):
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ] = {}
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ]["weight"] = float(rows_list[row][16 - 1])
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ]["month"] = float(rows_list[row][17 - 1])

    resolve.commit()

    # Subproblems and stages (one subproblem, one stage)
    subproblems = [1]
    subproblem_stages = {sid: [(1, "single stage")] for sid in subproblems}

    # Timepoints
    subproblem_stage_timepoints = dict()
    for subproblem_id in subproblem_stages.keys():
        subproblem_stage_timepoints[subproblem_id] = dict()
        for stage in subproblem_stages[subproblem_id]:
            stage_id = stage[0]
            subproblem_stage_timepoints[subproblem_id][stage_id] = dict()
            for _period in [2030]:
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
        for period in [2030]:
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
            temporal_scenario_id=2,
            scenario_name="2030 only",
            scenario_description="2030 only; 37 RESOLVE days, "
                                 "24 hours each",
            periods={2030: periods[2030]},
            subproblems=[1],
            subproblem_stages={1: [(1, "single stage")]},
            subproblem_stage_timepoints=subproblem_stage_timepoints,
            subproblem_horizons=subproblem_horizons,
            subproblem_stage_timepoint_horizons=subproblem_stage_timepoint_horizons
    )


def load_temporal_data_2030_only_1horizon():
    """
    Add RESOLVE timepoints/days into database
    Investment periods are 2030
    Period is represented by just one horizon (day)
    Discount factors and number of years represented are the same as in RESOLVE
    (the GridPath discount_factor is the RESOLVE weighted divided by the
    years in period)
    Horizon boundary is circular

    The point of this temporal scenario is to create a small toy model that
    can be run very quickly.
    :return:
    """

    # Make a dictionary with the discount factor and number of years
    # represented for each investment period
    periods = dict(dict())
    with open(os.path.join("cpuc_irp_data", "csvs", "Lists.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(4 - 1, 7):
            for column in range(32 - 1, 34):
                periods[
                    int(float(rows_list[row][32 - 1]))
                ] = {}
                periods[
                    int(float(rows_list[row][32 - 1]))
                ]["discount_factor"] = float(rows_list[row][33 - 1]) / \
                          float(rows_list[row][34 - 1])
                periods[
                    int(float(rows_list[row][32 - 1]))
                ]["number_years_represented"] = float(rows_list[row][34 - 1])

    # Make a dictionary of horizons with their weights and months
    # We'll also keep these in the RESOLVE database

    c1.execute(
        """DROP TABLE IF EXISTS resolve_days;"""
    )
    c1.execute(
        """CREATE TABLE resolve_days(
        day_id INTEGER PRIMARY KEY,
        day_weight FLOAT,
        month_of_year INTEGER,
        hydro_year INTEGER
        );"""
    )

    # Make dictionary for GridPath database
    horizon_weights_and_months = dict()

    for row in range(4 - 1, 40):
        # Save in RESOLVE database
        c1.execute(
            """INSERT INTO resolve_days (day_id, day_weight,
            month_of_year, hydro_year) VALUES ({}, {}, {}, {});""".format(
                int(float(rows_list[row][15 - 1])), rows_list[row][16 - 1],
                int(float(rows_list[row][17 - 1])),
                int(float(rows_list[row][18 - 1]))
            )
        )

        # Populate dictionary for GridPath database
        for column in range(15 - 1, 17):
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ] = {}
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ]["weight"] = float(rows_list[row][16 - 1])
            horizon_weights_and_months[
                int(float(rows_list[row][15 - 1]))
            ]["month"] = float(rows_list[row][17 - 1])

    resolve.commit()

    # Subproblems and stages (one subproblem, one stage)
    subproblems = [1]
    subproblem_stages = {sid: [(1, "single stage")] for sid in subproblems}

    # Timepoints
    subproblem_stage_timepoints = dict()
    for subproblem_id in subproblem_stages.keys():
        subproblem_stage_timepoints[subproblem_id] = dict()
        for stage in subproblem_stages[subproblem_id]:
            stage_id = stage[0]
            subproblem_stage_timepoints[subproblem_id][stage_id] = dict()
            for _period in [2030]:
                _day = list(horizon_weights_and_months.keys())[0]
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
                        stage_id][timepoint]["timepoint_weight"] = 365
                    subproblem_stage_timepoints[subproblem_id][
                        stage_id][timepoint][
                        "previous_stage_timepoint_map"] = 'NULL'
                    subproblem_stage_timepoints[subproblem_id][
                        stage_id][timepoint][
                        "spinup_or_lookahead"] = 'NULL'
                    subproblem_stage_timepoints[subproblem_id][
                        stage_id][timepoint]["month"] = \
                        int(horizon_weights_and_months[_day]["month"])
                    subproblem_stage_timepoints[subproblem_id][
                        stage_id][timepoint]["hour_of_day"] = hour

    # Horizons
    subproblem_horizons = dict()
    for subproblem_id in subproblem_stages.keys():
        subproblem_horizons[subproblem_id] = dict()
        for period in [2030]:
            day = list(horizon_weights_and_months.keys())[0]  # pick first day
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
            temporal_scenario_id=3,
            scenario_name="2030 only, 1 horizon",
            scenario_description="2030 only; 1 RESOLVE day, "
                                 "24 hours each",
            periods={2030: periods[2030]},
            subproblems=[1],
            subproblem_stages={1: [(1, "single stage")]},
            subproblem_stage_timepoints=subproblem_stage_timepoints,
            subproblem_horizons=subproblem_horizons,
            subproblem_stage_timepoint_horizons=subproblem_stage_timepoint_horizons
    )


def load_geography_load_zones():
    """
    Six load zones from the IRP
    The penalties are taken from the Inputs2Write tab (system_params)
    :return:
    """

    # Load data into GridPath database
    geography.geography_load_zones(
        io=io, c=c2,
        load_zone_scenario_id=1,
        scenario_name='default_6_zones',
        scenario_description = '6 zones: CAISO, NW, SW, LDWP, BANC, IID',
        zones=['CAISO', 'NW', 'SW', 'LDWP', 'BANC', 'IID'],
        zone_overgen_penalties={
            'CAISO': (1, 50000), 'NW': (1, 50000), 'SW': (1, 50000),
            'LDWP': (1, 50000), 'BANC': (1, 50000), 'IID': (1, 50000)
        },
        zone_unserved_energy_penalties={
            'CAISO': (1, 50000), 'NW': (1, 50000), 'SW': (1, 50000),
            'LDWP': (1, 50000), 'BANC': (1, 50000), 'IID': (1, 50000)
        }
    )


def load_geography_lf_reserves_up_bas():
    """
    CAISO-only LF up reserves from IRP
    Violation penalties are taken from Inputs2Write tab (system_params)
    :return:
    """
    # Load data into GridPath database
    geography.geography_lf_reserves_up_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        scenario_name='default_caiso_only_pt2_reserve_to_energy_adj',
        scenario_description='CAISO only, no reserve-to-energy adjustment',
        bas=['CAISO'],
        ba_penalties={'CAISO': (1, 10000)},
        reserve_to_energy_adjustments={'CAISO': 0.2}
    )


def load_geography_lf_reserves_down_bas():
    """
    CAISO-only LF down reserves from IRP
    Violation penalties are taken from Inputs2Write tab (system_params)
    :return:
    """
    geography.geography_lf_reserves_down_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        scenario_name='default_caiso_only_pt2_reserve_to_energy_adj',
        scenario_description='CAISO only, no reserve-to-energy adjustment',
        bas=['CAISO'],
        ba_penalties={'CAISO': (1, 10000)},
        reserve_to_energy_adjustments={'CAISO': 0.2}
    )


def load_geography_regulation_up_bas():
    """
    CAISO-only reg up reserves from IRP
    Violation penalties are taken from Inputs2Write tab (system_params)
    :return:
    """
    geography.geography_regulation_up_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        scenario_name='default_caiso_only_pt2_reserve_to_energy_adj',
        scenario_description='CAISO only, no reserve-to-energy adjustment',
        bas=['CAISO'],
        ba_penalties={'CAISO': (1, 10000)},
        reserve_to_energy_adjustments={'CAISO': 0.2}
    )


def load_geography_regulation_down_bas():
    """
    CAISO-only reg down reserves from IRP
    Violation penalties are taken from Inputs2Write tab (system_params)
    :return:
    """
    geography.geography_regulation_down_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        scenario_name='default_caiso_only_pt2_reserve_to_energy_adj',
        scenario_description='CAISO only, no reserve-to-energy adjustment',
        bas=['CAISO'],
        ba_penalties={'CAISO': (1, 10000)},
        reserve_to_energy_adjustments={'CAISO': 0.2}
    )


def load_geography_spinning_reserves_bas():
    """
    CAISO-only spinning reserves from IRP
    Violation penalties are taken from Inputs2Write tab (system_params)
    :return:
    """
    geography.geography_spinning_reserves_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        scenario_name='default_caiso_only_pt2_reserve_to_energy_adj',
        scenario_description='CAISO only, no reserve-to-energy adjustment',
        bas=['CAISO'],
        ba_penalties={'CAISO': (1, 10000)},
        reserve_to_energy_adjustments={'CAISO': 0.2}
    )


def load_geography_frequency_response_bas():
    """
    CAISO-only frequency response reserves from IRP
    For frequency response, assume penalty is 10000 (no data from RESOLVE)
    :return:
    """
    geography.geography_frequency_response_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        scenario_name='default_caiso_only_pt2_reserve_to_energy_adj',
        scenario_description='CAISO only, no reserve-to-energy adjustment',
        bas=['CAISO'],
        ba_penalties={'CAISO': (1, 10000)},
        reserve_to_energy_adjustments={'CAISO': 0.2}
    )


def load_geography_rps_zones():
    """
    CAISO rps zone
    :return:
    """
    geography.geography_rps_zones(
        io=io, c=c2,
        rps_zone_scenario_id=1,
        scenario_name='caiso_rps',
        scenario_description='CAISO RPS only',
        zones=['CAISO'],
        zone_penalties={'CAISO': (0, 0)}
    )


def load_geography_carbon_cap_zones():
    """
    CAISO carbon cap zone
    :return:
    """
    geography.geography_carbon_cap_zones(
        io=io, c=c2,
        carbon_cap_zone_scenario_id=defaults["carbon_cap_zone_scenario_id"],
        scenario_name='caiso_carbon_cap',
        scenario_description='CAISO carbon cap only',
        zones=['CAISO'],
        zone_penalties={'CAISO': (0, 0)}
    )


def load_geography_prm_zones():
    """
    CAISO PRM zone
    :return:
    """
    geography.geography_prm_zones(
        io=io, c=c2,
        prm_zone_scenario_id=1,
        scenario_name='caiso_only_prm',
        scenario_description='CAISO PRM only',
        zones=['CAISO'],
        zone_penalties={'CAISO': (0, 0)}
    )


def load_geography_local_capacity_zones():
    """
    CAISO PRM zone
    :return:
    """
    geography.geography_local_capacity_zones(
        io=io, c=c2,
        local_capacity_zone_scenario_id=1,
        scenario_name='caiso_local_zones',
        scenario_description='CAISO local zones',
        zones=[
            'Bay_Area',
            'Big_Creek_Ventura',
            'Fresno',
            'Humboldt',
            'Kern',
            'LA_Basin',
            'NCNB',
            'San_Diego',
            'Sierra',
            'Stockton'
        ],
        zone_penalties = {
            'Bay_Area': (1, 1540000),
            'Big_Creek_Ventura': (1, 1540000),
            'Fresno': (1, 1540000),
            'Humboldt': (1, 1540000),
            'Kern': (1, 1540000),
            'LA_Basin': (1, 1540000),
            'NCNB': (1, 1540000),
            'San_Diego': (1, 1540000),
            'Sierra': (1, 1540000),
            'Stockton': (1, 1540000)
        }
    )
    io.commit()


def load_projects():
    """
    Get a list of all projects
    Get the projects from the Inputs2Write tab in the RESOLVE sceanario tool
    (hardcoded in that tab), with the exception of the candidate renewable
    resources, which we'll get from the REN_Candidate tab
    :return:
    """

    projects = list()

    # All but candidate renewables
    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        i2w_rows_list = list(csv.reader(f))

        for row in range(3 - 1, 94):
            # GridPath needs two reciprocating engine projects (existing and
            # candidate)
            # We'll also split the CAISO_Peaker1, CAISO_Peaker1,
            # and CAISO_CCGT1 fleets into existing (same name) and
            # planned for the retirement cases (only existing allowed to
            # retire, not what's added during the study period)
            if i2w_rows_list[row][91 - 1] == 'CAISO_Reciprocating_Engine':
                for tail in ["", "_Candidate"]:
                    projects.append(i2w_rows_list[row][91 - 1] + tail)
            elif i2w_rows_list[row][91 - 1] in [
                'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
            ]:
                for tail in ["", "_Planned"]:
                    projects.append(i2w_rows_list[row][91 - 1] + tail)
            else:
                projects.append(i2w_rows_list[row][91 - 1])

    # Candidate renewables
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        # This will make a list of lists where each (sub)list is a row and
        # each element of the (sub)list is a column
        rencand_rows_list = list(csv.reader(f))

        for row in range(11 - 1, 52):
            projects.append(rencand_rows_list[row][3 - 1])

    # Add disaggregated projects for peaker analysis
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            projects.append(row[0])

    # # Add projects for CCGT analysis
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         print(row[0])
    #         projects.append(row[0])

    project_list.project_list(
        io=io, c=c2,
        projects=projects
    )


def load_project_load_zones():
    """
    Assign 'modeled' load zone to all projects
    It appears the IRP is only modeling out-of-state renewables as delivered
    to the California border, not just locally -- that's the scenario we're
    creating here
    :return:
    """

    # Load zones
    project_load_zones = dict()

    # Get hard-coded load_zones directly from spreadsheet for all but
    # candidate renewables
    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        i2w_rows_list = list(csv.reader(f))

        for row in range(3 - 1, 94):
            # GridPath needs two reciprocating engine projects (existing and
            # candidate)
            # We'll also split the CAISO_Peaker1, CAISO_Peaker2,
            # and CAISO_CCGT1 fleets into existing and
            # planned for the retirement cases (only existing allowed to
            # retire, not what's added during the study period)
            if i2w_rows_list[row][91 - 1] == 'CAISO_Reciprocating_Engine':
                for tail in ["", "_Candidate"]:
                    project_load_zones[i2w_rows_list[row][91 - 1] + tail] = \
                        i2w_rows_list[row][93 - 1]
            elif i2w_rows_list[row][91 - 1] in [
                'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
            ]:
                for tail in ["", "_Planned"]:
                    project_load_zones[i2w_rows_list[row][91 - 1] + tail] = \
                        i2w_rows_list[row][93 - 1]
            else:
                project_load_zones[i2w_rows_list[row][91 - 1]] = \
                i2w_rows_list[row][93 - 1]

    # Assign CAISO as load_zone for candidate renewables (all modeled as
    # delivered to CAISO)
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        # This will make a list of lists where each (sub)list is a row and
        # each element of the (sub)list is a column
        rencand_rows_list = list(csv.reader(f))

        for row in range(11 - 1, 52):
            project_load_zones[rencand_rows_list[row][3 - 1]] = 'CAISO'

    # Add disaggregated gas projects for peaker analysis
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            project_load_zones[row[0]] = 'CAISO'

    # # Add disaggregated gas projects for CCGT analysis
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_load_zones[row[0]] = 'CAISO'

    # Add LCR-eligible batteries
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "storage_lcr.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            project_load_zones[row[0]] = 'CAISO'

    project_zones.project_load_zones(
        io=io, c=c2,
        load_zone_scenario_id=1, # Default 6 zones
        project_load_zone_scenario_id=1,
        scenario_name=
        'oos_ren_delivered_to_ca_border',
        scenario_description=
        "Default six zones from IRP, OOS renewables are delivered to the CA "
        "border (modeled as if in CAISO).",
        project_load_zones=project_load_zones
    )


def load_project_lf_reserves_bas():
    """
    Assign BA to projects that can provide the reserve type
    In RESOLVE, the same projects can provide up and down LF
    :return:
    """
    project_reserve_bas_no_ren = dict()

    # Get hard-coded booleans (can/cannot provide) directly from
    # spreadsheet; assign 'CAISO' as BA if tech/project can provide
    # In RESOLVE, the flag is by technology, so we have to assign the
    # reserve zone to all projects from that technology

    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(3 - 1, 39):
            tech = rows_list[row][106 - 1]
            can_provide_reserves = int(float(rows_list[row][112 - 1]))

            # The pumped hydro projecs actually have "Pumped_Storage" not
            # "Pumped_Hydro" in the name, so update the tech
            if tech == 'Pumped_Hydro':
                tech = 'Pumped'
            else:
                pass

            # If tech can provide, find projects of that tech with 'CAISO'
            # in the name (all projects except for candidate renewables
            # have their physical load zone in their name in the data) and
            # assign "CAISO" as reserve BA
            # Need to exclude 'Small_Hydro' manually, as it gets picked up
            # with the 'Hydro' technology
            if can_provide_reserves:
                projects = [
                    p[0] for p in c2.execute(
                        """SELECT project
                        FROM inputs_project_all
                        WHERE project LIKE '%{}%'
                        AND project LIKE '%CAISO%'
                        AND project NOT LIKE '%Small_Hydro%'""".format(
                            tech
                        )
                    ).fetchall()
                ]

                for project in projects:
                    project_reserve_bas_no_ren[project] = "CAISO"

    project_reserve_bas_no_ren_lf_up = project_reserve_bas_no_ren
    project_reserve_bas_no_ren_lf_down = project_reserve_bas_no_ren

    # Add disaggregated gas projects not flagged as 'baseload'
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[13] != 'CHP':
                project_reserve_bas_no_ren[row[0]] = 'CAISO'

    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         if int(float(row[28 - 1])):
    #             project_reserve_bas_no_ren_lf_up[row[0]] = 'CAISO'
    #
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         if int(float(row[29 - 1])):
    #             project_reserve_bas_no_ren_lf_down[row[0]] = 'CAISO'

    # # Add LCR-eligible batteries
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "storage_lcr.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_reserve_bas_no_ren[row[0]] = 'CAISO'


    # LF up
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=1,
        scenario_name="thermal_hydro_storage_only",
        scenario_description="Only thermal, hydro, and storage can provide "
                             "reserves.",
        project_bas=project_reserve_bas_no_ren_lf_up,
        reserve_type="lf_reserves_up"
    )

    # LF down
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=1,
        scenario_name="thermal_hydro_storage_only",
        scenario_description="Only thermal, hydro, and storage can provide "
                             "reserves.",
        project_bas=project_reserve_bas_no_ren_lf_down,
        reserve_type="lf_reserves_down"
    )

    # Add scenarios where variable renewables can provide reserves
    # Start with thermal, hydro, and storage
    project_reserve_bas_w_ren_lf_up = project_reserve_bas_no_ren_lf_up
    project_reserve_bas_w_ren_lf_down = project_reserve_bas_no_ren_lf_down

    # Exclude all non-CAISO projects
    # 'CAISO_Wind_for_Other' and 'CAISO_Solar_for_Other' are included
    caiso_var_ren_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE (project LIKE '%Solar%' OR project LIKE '%Wind%')
            AND project NOT LIKE 'BANC%'
            AND project NOT LIKE 'IID%'
            AND project NOT LIKE 'LDWP%'
            AND project NOT LIKE 'NW%'
            AND project NOT LIKE 'SW%';"""
        ).fetchall()
    }

    # Manually add back the SW and NW existing Tx wind projects
    existing_tx_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE project LIKE '%Ext_Tx%'"""
        ).fetchall()
    }

    # Add exsiting Tx projects to the variable renewables dictionary
    caiso_var_ren_projects.update(existing_tx_projects)
    # Add the variable renewables to the conventional projects
    project_reserve_bas_no_ren_lf_up.update(caiso_var_ren_projects)
    project_reserve_bas_no_ren_lf_down.update(caiso_var_ren_projects)

    # Load data for variable renewables providing reserves scenario
    # LF up
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=2,
        scenario_name="var_renewables_thermal_hydro_storage",
        scenario_description="Variable renewables can provide reserves in "
                             "addition to thermal, hydro, and storage.",
        project_bas=project_reserve_bas_no_ren_lf_up,
        reserve_type="lf_reserves_up"
    )

    # LF down
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=2,
        scenario_name="var_renewables_thermal_hydro_storage",
        scenario_description="Variable renewables can provide reserves in "
                             "addition to thermal, hydro, and storage.",
        project_bas=project_reserve_bas_no_ren_lf_down,
        reserve_type="lf_reserves_down"
    )


def load_project_regulation_bas():
    """
    Assign BA to projects that can provide the reserve type
    In RESOLVE, the same projects can provide up and down regulation
    :return:
    """
    project_reserve_bas_no_ren = dict()

    # Get hard-coded booleans (can/cannot provide) directly from
    # spreadsheet; assign 'CAISO' as BA if tech/project can provide
    # In RESOLVE, the flag is by technology, so we have to assign the
    # reserve zone to all projects from that technology

    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(3 - 1, 39):
            tech = rows_list[row][106 - 1]
            can_provide_reserves = int(float(rows_list[row][111 - 1]))

            # The pumped hydro projecs actually have "Pumped_Storage" not
            # "Pumped_Hydro" in the name, so update the tech
            if tech == 'Pumped_Hydro':
                tech = 'Pumped'
            else:
                pass

            # If tech can provide, find projects of that tech with 'CAISO'
            # in the name (all projects except for candidate renewables
            # have their physical load zone in their name in the data) and
            # assign "CAISO" as reserve BA
            # Need to exclude 'Small_Hydro' manually, as it gets picked up
            # with the 'Hydro' technology
            if can_provide_reserves:
                projects = [
                    p[0] for p in c2.execute(
                        """SELECT project
                        FROM inputs_project_all
                        WHERE project LIKE '%{}%'
                        AND project LIKE '%CAISO%'
                        AND project NOT LIKE '%Small_Hydro%'""".format(
                            tech
                        )
                    ).fetchall()
                ]

                for project in projects:
                    project_reserve_bas_no_ren[project] = "CAISO"

    project_reserve_bas_no_ren_reg_up = project_reserve_bas_no_ren
    project_reserve_bas_no_ren_reg_down = project_reserve_bas_no_ren

    # Add disaggregated gas projects
    # Add disaggregated gas projects not flagged as 'baseload'
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[13] != 'CHP':
                project_reserve_bas_no_ren[row[0]] = 'CAISO'
    #
    # # Add LCR-eligible batteries
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "storage_lcr.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_reserve_bas_no_ren[row[0]] = 'CAISO'

    # Regulation up
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=1,
        scenario_name="thermal_hydro_storage_only",
        scenario_description="Only thermal, hydro, and storage can provide "
                             "reserves.",
        project_bas=project_reserve_bas_no_ren_reg_up,
        reserve_type="regulation_up"
    )

    # Regulation down
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=1,
        scenario_name="thermal_hydro_storage_only",
        scenario_description="Only thermal, hydro, and storage can provide "
                             "reserves.",
        project_bas=project_reserve_bas_no_ren_reg_down,
        reserve_type="regulation_down"
    )

    # Add scenarios where variable renewables can provide reserves
    # Start with thermal, hydro, and storage
    project_reserve_bas_w_ren_reg_up = project_reserve_bas_no_ren_reg_up
    project_reserve_bas_w_ren_reg_down = project_reserve_bas_no_ren_reg_down

    # Exclude all non-CAISO projects
    caiso_var_ren_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE (project LIKE '%Solar%' OR project LIKE '%Wind%')
            AND project NOT LIKE 'BANC%'
            AND project NOT LIKE 'IID%'
            AND project NOT LIKE 'LDWP%'
            AND project NOT LIKE 'NW%'
            AND project NOT LIKE 'SW%';"""
        ).fetchall()
    }

    # Manually add back the SW and NW existing Tx wind projects
    existing_tx_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE project LIKE '%Ext_Tx%'"""
        ).fetchall()
    }

    # Add exsiting Tx projects to the variable renewables dictionary
    caiso_var_ren_projects.update(existing_tx_projects)
    # Add the variable renewables to the conventional projects
    project_reserve_bas_w_ren_reg_up.update(caiso_var_ren_projects)
    project_reserve_bas_w_ren_reg_down.update(caiso_var_ren_projects)

    # Load data for variable renewables providing reserves scenario
    # Regulation up
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=2,
        scenario_name="var_renewables_thermal_hydro_storage",
        scenario_description="Variable renewables can provide reserves in "
                             "addition to thermal, hydro, and storage.",
        project_bas=project_reserve_bas_w_ren_reg_up,
        reserve_type="regulation_up"
    )

    # Regulation down
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=2,
        scenario_name="var_renewables_thermal_hydro_storage",
        scenario_description="Variable renewables can provide reserves in "
                             "addition to thermal, hydro, and storage.",
        project_bas=project_reserve_bas_w_ren_reg_down,
        reserve_type="regulation_down"
    )


def load_project_spinning_reserves_bas():
    """
    Assign BA to projects that can provide the reserve type
    :return:
    """
    project_reserve_bas_no_ren = dict()

    # Get hard-coded booleans (can/cannot provide) directly from
    # spreadsheet; assign 'CAISO' as BA if tech/project can provide
    # In RESOLVE, the flag is by technology, so we have to assign the
    # reserve zone to all projects from that technology

    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(3 - 1, 39):
            tech = rows_list[row][106 - 1]
            can_provide_reserves = int(float(rows_list[row][110 - 1]))

            # The pumped hydro projecs actually have "Pumped_Storage" not
            # "Pumped_Hydro" in the name, so update the tech
            if tech == 'Pumped_Hydro':
                tech = 'Pumped'
            else:
                pass

            # If tech can provide, find projects of that tech with 'CAISO'
            # in the name (all projects except for candidate renewables
            # have their physical load zone in their name in the data) and
            # assign "CAISO" as reserve BA
            # Need to exclude 'Small_Hydro' manually, as it gets picked up
            # with the 'Hydro' technology
            if can_provide_reserves:
                projects = [
                    p[0] for p in c2.execute(
                        """SELECT project
                        FROM inputs_project_all
                        WHERE project LIKE '%{}%'
                        AND project LIKE '%CAISO%'
                        AND project NOT LIKE '%Small_Hydro%'""".format(
                            tech
                        )
                    ).fetchall()
                ]

                for project in projects:
                    project_reserve_bas_no_ren[project] = "CAISO"

    project_reserve_bas_no_ren_spin = project_reserve_bas_no_ren

    # Add disaggregated gas projects not flagged as 'baseload'
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[13] != 'CHP':
                project_reserve_bas_no_ren[row[0]] = 'CAISO'
    #
    # # Add LCR-eligible batteries
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "storage_lcr.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_reserve_bas_no_ren[row[0]] = 'CAISO'

    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=1,
        scenario_name="thermal_hydro_storage_only",
        scenario_description="Only thermal, hydro, and storage can provide "
                             "reserves.",
        project_bas=project_reserve_bas_no_ren_spin,
        reserve_type="spinning_reserves"
    )

    # Add scenarios where variable renewables can provide reserves
    # Start with thermal, hydro, and storage
    project_reserve_bas_w_ren = project_reserve_bas_no_ren_spin

    # Exclude all non-CAISO projects
    caiso_var_ren_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE (project LIKE '%Solar%' OR project LIKE '%Wind%')
            AND project NOT LIKE 'BANC%'
            AND project NOT LIKE 'IID%'
            AND project NOT LIKE 'LDWP%'
            AND project NOT LIKE 'NW%'
            AND project NOT LIKE 'SW%';"""
        ).fetchall()
    }

    # Manually add back the SW and NW existing Tx wind projects
    existing_tx_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE project LIKE '%Ext_Tx%'"""
        ).fetchall()
    }

    # Add exsiting Tx projects to the variable renewables dictionary
    caiso_var_ren_projects.update(existing_tx_projects)
    # Add the variable renewables to the conventional projects
    project_reserve_bas_no_ren_spin.update(caiso_var_ren_projects)

    # Load data for variable renewables providing reserves scenario
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=2,
        scenario_name="var_renewables_thermal_hydro_storage",
        scenario_description="Variable renewables can provide reserves in "
                             "addition to thermal, hydro, and storage.",
        project_bas=project_reserve_bas_no_ren_spin,
        reserve_type="spinning_reserves"
    )


def load_project_frequency_response_bas():
    """
    Assign BA to projects that can provide the reserve type
    :return:
    """
    project_reserve_bas_no_ren = dict()
    project_reserve_bas_no_ren_freq_partial = dict()

    # Get hard-coded booleans (can/cannot provide) directly from
    # spreadsheet; assign 'CAISO' as BA if tech/project can provide
    # In RESOLVE, the flag is by technology, so we have to assign the
    # reserve zone to all projects from that technology

    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(3 - 1, 39):
            tech = rows_list[row][106 - 1]
            can_provide_reserves = int(float(rows_list[row][114 - 1]))
            partial = int(float(rows_list[row][115 - 1]))

            # The pumped hydro projecs actually have "Pumped_Storage" not
            # "Pumped_Hydro" in the name, so update the tech
            if tech == 'Pumped_Hydro':
                tech = 'Pumped'
            else:
                pass

            # If tech can provide, find projects of that tech with 'CAISO'
            # in the name (all projects except for candidate renewables
            # have their physical load zone in their name in the data) and
            # assign "CAISO" as reserve BA
            # Need to exclude 'Small_Hydro' manually, as it gets picked up
            # with the 'Hydro' technology
            if can_provide_reserves:
                projects = [
                    p[0] for p in c2.execute(
                        """SELECT project
                        FROM inputs_project_all
                        WHERE project LIKE '%{}%'
                        AND project LIKE '%CAISO%'
                        AND project NOT LIKE '%Small_Hydro%'""".format(
                            tech
                        )
                    ).fetchall()
                ]

                for project in projects:
                    project_reserve_bas_no_ren[project] = "CAISO"
                    project_reserve_bas_no_ren_freq_partial[project] = partial

    project_reserve_bas_no_ren_freq = project_reserve_bas_no_ren

    # Add disaggregated gas projects not flagged as 'baseload'
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[13] != 'CHP':
                project_reserve_bas_no_ren[row[0]] = 'CAISO'
                project_reserve_bas_no_ren_freq_partial[row[0]] = 1

    #
    # # Add LCR-eligible batteries
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "storage_lcr.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_reserve_bas_no_ren[row[0]] = 'CAISO'
    #         project_reserve_bas_partial_no_ren[row[0]] = 1


    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=1,
        scenario_name="thermal_hydro_storage_only",
        scenario_description="Only thermal, hydro, and storage can provide "
                             "reserves.",
        project_bas=project_reserve_bas_no_ren_freq,
        reserve_type="frequency_response"
    )

    # Update partial frequency response flag
    for prj in project_reserve_bas_no_ren_freq_partial.keys():
        c2.execute(
            """UPDATE inputs_project_frequency_response_bas
            SET contribute_to_partial = {}
            WHERE project = '{}'
            AND frequency_response_ba_scenario_id = {}
            AND project_frequency_response_ba_scenario_id = {};""".format(
                project_reserve_bas_no_ren_freq_partial[prj], prj, 1, 1
            )
        )
    io.commit()

    # Add scenarios where variable renewables can provide reserves
    # Start with thermal, hydro, and storage
    project_reserve_bas_w_ren = project_reserve_bas_no_ren

    # Exclude all non-CAISO projects
    caiso_var_ren_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE (project LIKE '%Solar%' OR project LIKE '%Wind%')
            AND project NOT LIKE 'BANC%'
            AND project NOT LIKE 'IID%'
            AND project NOT LIKE 'LDWP%'
            AND project NOT LIKE 'NW%'
            AND project NOT LIKE 'SW%';"""
        ).fetchall()
    }

    # Manually add back the SW and NW existing Tx wind projects
    existing_tx_projects = {
        p[0]: 'CAISO' for p in c2.execute(
            """SELECT project
            FROM inputs_project_all
            WHERE project LIKE '%Ext_Tx%'"""
        ).fetchall()
    }

    # Add exsiting Tx projects to the variable renewables dictionary
    caiso_var_ren_projects.update(existing_tx_projects)
    # Add the variable renewables to the conventional projects
    project_reserve_bas_no_ren_freq.update(caiso_var_ren_projects)

    # Load data for variable renewables providing reserves scenario
    project_zones.project_reserve_bas(
        io=io, c=c2,
        reserve_ba_scenario_id=1,
        project_reserve_scenario_id=2,
        scenario_name="var_renewables_thermal_hydro_storage",
        scenario_description="Variable renewables can provide reserves in "
                             "addition to thermal, hydro, and storage.",
        project_bas=project_reserve_bas_w_ren,
        reserve_type="frequency_response"
    )

    # Update partial frequency response flag
    for prj in project_reserve_bas_no_ren_freq_partial.keys():
        c2.execute(
            """UPDATE inputs_project_frequency_response_bas
            SET contribute_to_partial = {}
            WHERE project = '{}'
            AND frequency_response_ba_scenario_id = {}
            AND project_frequency_response_ba_scenario_id = {};""".format(
                project_reserve_bas_no_ren_freq_partial[prj], prj, 1, 2
            )
        )

    # Don't let the variable renewables provide partial frequency response
    for prj in caiso_var_ren_projects.keys():
        c2.execute(
            """UPDATE inputs_project_frequency_response_bas
            SET contribute_to_partial = {}
            WHERE project = '{}'
            AND frequency_response_ba_scenario_id = {}
            AND project_frequency_response_ba_scenario_id = {};""".format(
                0, prj, 1, 2
            )
        )
    io.commit()


def load_project_rps_zones():
    """
    Project RPS zones
    :return:
    """
    project_rps_zones = dict()

    # Get hard-coded 'contract' directly from spreadsheet for all but
    # candidate renewables, which we'll get from REN_Candidate
    # Existing/planned renewables are in rows 50-93 (Customer_PV is excluded)
    # We'll only add CAISO renewables to the dictionary, as the rest don't
    # matter since we're not enforcing RPS outside of CAISO
    # We also need to model CAISO storage as counting toward the RPS to
    # prevent storage from acting as load
    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        i2w_rows_list = list(csv.reader(f))

        # Renewables
        for row in range(50 - 1, 93):
            if i2w_rows_list[row][94 - 1] == 'CAISO':
                project_rps_zones[i2w_rows_list[row][91 - 1]] = \
                    i2w_rows_list[row][94 - 1]

        # Storage
        for row in range(38 - 1, 42):
            if i2w_rows_list[row][94 - 1] == 'CAISO':  # might as well check
                project_rps_zones[i2w_rows_list[row][91 - 1]] = \
                    i2w_rows_list[row][94 - 1]

    # LCR-eligible batteries
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "storage_lcr.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            project_rps_zones[row[0]] = 'CAISO'

    # Assign CAISO as rps_zone for candidate renewables (all modeled as
    # delivered to CAISO)
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        # This will make a list of lists where each (sub)list is a row and
        # each element of the (sub)list is a column
        rencand_rows_list = list(csv.reader(f))

        for row in range(11 - 1, 52):
            project_rps_zones[rencand_rows_list[row][3 - 1]] = 'CAISO'

    # Insert data
    project_zones.project_policy_zones(
        io=io, c=c2,
        policy_zone_scenario_id=1,
        project_policy_zone_scenario_id=1,
        scenario_name="renewables_and_storage_losses",
        scenario_description="CAISO RPS, storage losses included in the RPS.",
        project_zones=project_rps_zones,
        policy_type="rps"
    )


def load_project_prm_zones():
    """
    Project PRM zones
    :return:
    """
    project_prm_zones = dict()

    # All projects listed in SYS_Planning_Reserve tab will contribute via
    # simple ELCC fraction
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Planning_Reserve.csv"), "r") as f:
        # This will make a list of lists where each (sub)list is a row and
        # each element of the (sub)list is a column
        rows_list = list(csv.reader(f))

        for row in range(22 - 1, 50):
            if rows_list[row][3 - 1] == 'CAISO_Reciprocating_Engine':
                for tail in ["", "_Candidate"]:
                    project_prm_zones[rows_list[row][3 - 1] + tail] = \
                        'CAISO'
            elif rows_list[row][3 - 1] in [
                'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
            ]:
                for tail in ["", "_Planned"]:
                    project_prm_zones[rows_list[row][3 - 1] + tail] = \
                        'CAISO'
            project_prm_zones[rows_list[row][3 - 1]] = 'CAISO'

    # CAISO storage will also contribute (with a duration derate)
    # CAISO variable renewables located in CAISO (physical or via
    # transmission) will contribute via an ELCC surface or simple fraction
    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        i2w_rows_list = list(csv.reader(f))

        # Storage
        for row in range(38 - 1, 42):
            project_prm_zones[i2w_rows_list[row][91 - 1]] = 'CAISO'

        # Existing renewables (located in CAISO, including those contracted
        # to other zones)
        for row in range(49 - 1, 93):
            if i2w_rows_list[row][93 - 1] == 'CAISO':
                project_prm_zones[i2w_rows_list[row][91 - 1]] = \
                    i2w_rows_list[row][93 - 1]

    # Candidate renewables are all delivered to CAISO, so can contribute to PRM
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        # This will make a list of lists where each (sub)list is a row and
        # each element of the (sub)list is a column
        rencand_rows_list = list(csv.reader(f))

        for row in range(11 - 1, 52):
            project_prm_zones[rencand_rows_list[row][3 - 1]] = 'CAISO'

    # Customer_PV contributes via the surface
    project_prm_zones['Customer_PV'] = 'CAISO'

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            project_prm_zones[row[0]] = 'CAISO'

    # # Add disaggregated gas projects
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_prm_zones[row[0]] = 'CAISO'

    #
    # # LCR-eligible batteries
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "storage_lcr.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_prm_zones[row[0]] = 'CAISO'

    # Insert data
    project_zones.project_policy_zones(
        io=io, c=c2,
        policy_zone_scenario_id=1,
        project_policy_zone_scenario_id=1,
        scenario_name="caiso_located_projects_incl_via_tx_contribute",
        scenario_description="CAISO PRM, projects located in CAISO ("
                             "including via transmission) contribute, "
                             "so CAISO_X_for_Other projects are included.",
        project_zones=project_prm_zones,
        policy_type="prm"
    )


def load_project_local_capacity_zones():
    """
    Local capacity zones of projects modeled as contributing to the local
    capacity requirements
    :return:
    """

    c2.execute(
        """INSERT INTO subscenarios_project_local_capacity_zones 
        (local_capacity_zone_scenario_id, 
        project_local_capacity_zone_scenario_id, name, description) 
        VALUES (1, 1, 'default', 'default');"""
    )
    io.commit()

    # project_lc_zones = OrderedDict()
    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "gas_disagg_plants.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         project_lc_zones[row[0]] = row[17]

    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[17] == 'CAISO_System':
                pass
            else:
                c2.execute(
                    """INSERT INTO inputs_project_local_capacity_zones
                    (local_capacity_zone_scenario_id, 
                    project_local_capacity_zone_scenario_id, project, 
                    local_capacity_zone)
                    VALUES (1, 1, '{}', '{}');""".format(
                        row[0], row[17]
                    )
                )

    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         if row[17] == 'CAISO_System':
    #             pass
    #         else:
    #             c2.execute(
    #                 """INSERT INTO inputs_project_local_capacity_zones
    #                 (local_capacity_zone_scenario_id,
    #                 project_local_capacity_zone_scenario_id, project,
    #                 local_capacity_zone)
    #                 VALUES (1, 1, '{}', '{}');""".format(
    #                     row[0], row[17]
    #                 )
    #             )

    io.commit()


def load_project_carbon_cap_zones():
    """
    Project carbon cap zones
    :return:
    """
    project_carbon_cap_zones = dict()

    # All thermal fossil projects in CAISO
    # No need to add zero carbon nuclear
    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        i2w_rows_list = list(csv.reader(f))

        for row in range(3 - 1, 12):
            if i2w_rows_list[row][91 - 1] == 'CAISO_Nucleaer':
                pass
            elif i2w_rows_list[row][91 - 1] == 'CAISO_Reciprocating_Engine':
                for tail in ["", "_Candidate"]:
                    project_carbon_cap_zones[
                        i2w_rows_list[row][91 - 1] + tail
                    ] = 'CAISO'
            elif i2w_rows_list[row][91 - 1] in [
                'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
            ]:
                for tail in ["", "_Planned"]:
                    project_carbon_cap_zones[i2w_rows_list[
                        row][91 - 1] + tail
                    ] = 'CAISO'
            else:
                project_carbon_cap_zones[i2w_rows_list[row][91 - 1]] = 'CAISO'

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                project_carbon_cap_zones[row[0]] = 'CAISO'

    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         if row[0].startswith('Battery'):
    #             pass
    #         else:
    #             project_carbon_cap_zones[row[0]] = 'CAISO'

    # Add data
    project_zones.project_policy_zones(
        io=io, c=c2,
        policy_zone_scenario_id=1,
        project_policy_zone_scenario_id=1,
        scenario_name="in_caiso_thermal_fossil_generation",
        scenario_description="Thermal fossil generation in CAISO.",
        project_zones=project_carbon_cap_zones,
        policy_type="carbon_cap"
    )


def load_project_availability():
    """

    :return:
    """
    project_types_and_char_ids = dict()
    all_projects = c2.execute(
        "SELECT project FROM inputs_project_all;"
    ).fetchall()
    for prj in all_projects:
        project_types_and_char_ids[prj[0]] = {}
        project_types_and_char_ids[prj[0]]["type"] = "exogenous"
        project_types_and_char_ids[prj[0]]["exogenous_availability_id"] = None
        project_types_and_char_ids[prj[0]]["endogenous_availability_id"] = None

    project_availability.make_scenario_and_insert_types_and_ids(
        io=io, c=c2,
        project_availability_scenario_id=1,
        scenario_name="default",
        scenario_description="Default availabilities.",
        project_types_and_char_ids=project_types_and_char_ids
    )

    # Update the projects that are actually derated
    project_avail_scenarios = {}

    timepoints = c2.execute(
        """SELECT timepoint, month FROM inputs_temporal_timepoints
        WHERE temporal_scenario_id = {}
        AND subproblem_id = {}
        AND stage_id = {};""".format(1, 1, 1)
    ).fetchall()

    av_by_prj_month = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"),
              "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(38 - 1, 44):
            prj = rows_list[row][3 - 1]
            c2.execute("""
            UPDATE inputs_project_availability_types
            SET exogenous_availability_scenario_id = 1
            WHERE project = '{}';
            """.format(prj))
            av_by_prj_month[prj] = OrderedDict()
            for column in range(6 - 1, 17):
                month = int(float(rows_list[37 - 1][column]))
                av_by_prj_month[prj][month] = float(rows_list[row][column])

    avail_by_prj_tmp = OrderedDict()
    for prj in av_by_prj_month.keys():
        # exogenous_availability_scenario_id
        project_avail_scenarios[prj] = {}
        project_avail_scenarios[prj][1] = ("default", "default derate")

        # Derates
        avail_by_prj_tmp[prj] = OrderedDict()
        avail_by_prj_tmp[prj][
            1] = OrderedDict()  # exogenous_availability_scenario_id
        avail_by_prj_tmp[prj][1][1] = OrderedDict()  # stage
        for tmp_row in timepoints:
            tmp = tmp_row[0]
            month = tmp_row[1]
            avail_by_prj_tmp[prj][1][1][tmp] = av_by_prj_month[prj][month]

    project_availability.insert_project_availability_exogenous(
        io=io, c=c2,
        project_avail_scenarios=project_avail_scenarios,
        project_avail=avail_by_prj_tmp
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
    # 'Must run' gas disaggregated projects manually applied (should match
    # 'baseload' flag)
    project_op_types = dict()
    project_balancing_types = dict()
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "project_operational_types_gas_disagg.csv"),
            "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            project_op_types[row[0]] = row[1]
            project_balancing_types[row[0]] = "day"

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="operational_type",
        project_char=project_op_types
    )

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="balancing_type_project",
        project_char=project_balancing_types
    )

    # ### Technologies ### #
    # Tech for disaggregated gas projects is based on the E3 type
    project_tech = dict()
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "project_technologies_gas_disagg.csv"),
            "r") as f:
        reader = csv.reader(f)
        next(reader)  # skip header
        for row in reader:
            project_tech[row[0]] = row[1]

    # with open(os.path.join(
    #         "gridpath_specific", "plants_to_add_for_ccgt_analysis.csv"),
    #         "r") as f:
    #     reader = csv.reader(f)
    #     next(reader)  # skip header
    #     for row in reader:
    #         project_tech[row[0]] = row[13]

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="technology",
        project_char=project_tech
    )

    # ### Variable O&M costs ### #
    project_vom = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            project_vom[rows_list[row][3 - 1]] = float(rows_list[row][25 - 1])

    # Update 'conventional' projects
    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="variable_cost_per_mwh",
        project_char=project_vom
    )

    # The list from CONV_OpChar is actually by RESOLVE 'technology,' which
    # is not the same as the project in the case of DR
    # We have to also divide CAISO_Peaker1 and CAISO_Reciprocating_Engine into
    # two projects each for GridPath, so these also need to be updated
    # separately
    # Update DR, CAISO_Peaker1_Planned, CAISO_Peaker2_Planned,
    # CAISO_CCGT1_Planned, and CAISO_Reciprocating_Engine_Candidate separately
    c2.execute(
        """UPDATE inputs_project_operational_chars
        SET variable_cost_per_mwh = 600.0
        WHERE project LIKE 'CAISO_Shed_DR%'
        AND project_operational_chars_scenario_id = 1;"""
    )
    io.commit()

    c2.execute(
        """UPDATE inputs_project_operational_chars
        SET variable_cost_per_mwh = 5.0
        WHERE project LIKE 'CAISO_Peaker%_Planned'
        AND project_operational_chars_scenario_id = 1;"""
    )
    io.commit()

    c2.execute(
        """UPDATE inputs_project_operational_chars
        SET variable_cost_per_mwh = 5.7
        WHERE project = 'CAISO_CCGT1_Planned'
        AND project_operational_chars_scenario_id = 1;"""
    )
    io.commit()

    c2.execute(
        """UPDATE inputs_project_operational_chars
        SET variable_cost_per_mwh = 6.0
        WHERE project = 'CAISO_Reciprocating_Engine_Candidate'
        AND project_operational_chars_scenario_id = 1;"""
    )
    io.commit()

    # Of the other projects, only biomass projects have non-zero VOM cost
    # Set biomass project VOM to 4
    c2.execute(
        """UPDATE inputs_project_operational_chars
        SET variable_cost_per_mwh = 4.0
        WHERE project LIKE '%Biomass%'
        AND project_operational_chars_scenario_id = 1;"""
    )
    io.commit()

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            c2.execute(
                """UPDATE inputs_project_operational_chars
                SET variable_cost_per_mwh = {}
                WHERE project = '{}'
                AND project_operational_chars_scenario_id = 1;""".format(
                    row[11], row[0]
                )
            )

    # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
    #                        "plants_to_add_for_ccgt_analysis.csv"),
    #           "r") as f:
    #     reader = csv.reader(f, delimiter=",")
    #     next(reader)
    #     for row in reader:
    #         c2.execute(
    #             """UPDATE inputs_project_operational_chars
    #             SET variable_cost_per_mwh = {}
    #             WHERE project = '{}'
    #             AND project_operational_chars_scenario_id = 1;""".format(
    #                 row[11], row[0]
    #             )
    #         )

    io.commit()
    
    # Set variable cost for all remaining projects to 0
    c2.execute(
        """UPDATE inputs_project_operational_chars
        SET variable_cost_per_mwh = 0.0
        WHERE variable_cost_per_mwh IS NULL
        AND project_operational_chars_scenario_id = 1;"""
    )
    io.commit()

    # ### Fuel ### #
    project_fuel = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            project_fuel[rows_list[row][3 - 1]] = rows_list[row][5 - 1]

    # Update 'conventional' projects
    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="fuel",
        project_char=project_fuel
    )

    # The list from CONV_OpChar is actually by RESOLVE 'technology,' which
    # is not the same as the project in the case of DR
    # However, DR does not need a fuel in GridPath, so we'll skip updating it
    # We have to also divide CAISO_Peaker1 and CAISO_Reciprocating_Engine into
    # two projects each for GridPath, so these also need to be updated
    # separately
    # Update CAISO_Peaker1_Planned, CAISO_Peaker2_Planned, CAISO_CCGT1_Planned,
    # and CAISO_Reciprocating_Engine_Candidate separately
    c2.execute(
        """UPDATE inputs_project_operational_chars
        SET fuel = 'CA_Natural_Gas'
        WHERE project LIKE 'CAISO_Peaker%_Planned'
        OR project = 'CAISO_CCGT1_Planned'
        OR project = 'CAISO_Reciprocating_Engine_Candidate'
        AND project_operational_chars_scenario_id = 1;"""
    )
    io.commit()
    
    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith("Battery"):
                pass
            else:
                c2.execute(
                    """UPDATE inputs_project_operational_chars
                    SET fuel = 'CA_Natural_Gas'
                    WHERE project = '{}'
                    AND project_operational_chars_scenario_id = 1;""".format(
                        row[0]
                    )
                )

    io.commit()

    # ### Heat rates ### #
    # Assign a heat_curve_scenario_id to all fuel projects
    project_hr_curve_id = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            if rows_list[row][3 - 1] == 'CAISO_Conventional_DR':
                pass
            # If Pmax=Pmin, assign 0 to min input and heat rate at Pmax
            else:
                project_hr_curve_id[rows_list[row][3 - 1]] = 1

    project_hr_curve_id["CAISO_Peaker1_Planned"] = 1
    project_hr_curve_id["CAISO_Peaker2_Planned"] = 1
    project_hr_curve_id["CAISO_CCGT1_Planned"] = 1
    project_hr_curve_id["CAISO_Reciprocating_Engine_Candidate"] = 1
    # TODO: seems like this is a duplicate?
    project_hr_curve_id["CAISO_Peaker1_Planned"] = 1
    project_hr_curve_id["CAISO_Peaker2_Planned"] = 1
    project_hr_curve_id["CAISO_CCGT1_Planned"] = 1
    project_hr_curve_id["CAISO_Reciprocating_Engine_Candidate"] = 1

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                project_hr_curve_id[row[0]] = 1

    project_operational_chars.update_project_opchar_column(
        io, c2, 1, 'heat_rate_curves_scenario_id', project_hr_curve_id
    )


    # NOTE: RESOLVE data table header says Btu/MWh for the heat rates at
    # Pmin and Pmax, but this looks like it's MMBtu/MWh
    # Add manually

    hr_names = {
        'CAISO_CHP': {1: ('default', 'default')},
        'CAISO_Nuclear': {1: ('default', 'default')},
        'CAISO_CCGT1': {1: ('default', 'default')},
        'CAISO_CCGT2': {1: ('default', 'default')},
        'CAISO_Peaker1': {1: ('default', 'default')},
        'CAISO_Peaker2': {1: ('default', 'default')},
        'CAISO_Advanced_CCGT': {1: ('default', 'default')},
        'CAISO_Aero_CT': {1: ('default', 'default')},
        'CAISO_Reciprocating_Engine': {1: ('default', 'default')},
        'CAISO_ST': {1: ('default', 'default')},
        'NW_Nuclear': {1: ('default', 'default')},
        'NW_Coal': {1: ('default', 'default')},
        'NW_CCGT': {1: ('default', 'default')},
        'NW_Peaker': {1: ('default', 'default')},
        'SW_Nuclear': {1: ('default', 'default')},
        'SW_Coal': {1: ('default', 'default')},
        'SW_CCGT': {1: ('default', 'default')},
        'SW_Peaker': {1: ('default', 'default')},
        'LDWP_Nuclear': {1: ('default', 'default')},
        'LDWP_Coal': {1: ('default', 'default')},
        'LDWP_CCGT': {1: ('default', 'default')},
        'LDWP_Peaker': {1: ('default', 'default')},
        'IID_CCGT': {1: ('default', 'default')},
        'IID_Peaker': {1: ('default', 'default')},
        'BANC_CCGT': {1: ('default', 'default')},
        'BANC_Peaker': {1: ('default', 'default')},
        'CAISO_Peaker1_Planned': {1: ('default', 'default')},
        'CAISO_Peaker2_Planned': {1: ('default', 'default')},
        'CAISO_CCGT1_Planned': {1: ('default', 'default')},
        'CAISO_Reciprocating_Engine_Candidate': {1: ('default', 'default')}
    }

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                hr_names[row[0]] = {1: ('default', 'default')}

    # The heat rate curves; hard-coded for the aggregated projects
    hr_dict = dict()
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "heat_rate_curves.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery') or row[0] == \
                    'CAISO_Conventional_DR':
                pass
            else:
                hr_curve = {1: (row[1], row[2])}
                if hr_curve[1] == (row[3], row[4]):  # single bp for must-run
                    pass
                else:
                    hr_curve[2] = (row[3], row[4])

                hr_dict[row[0]] = {1: hr_curve}


    # hr_dict = {
    #     'CAISO_CHP': {1: {1: (19, 7.606)}},
    #     'CAISO_Nuclear': {1: {1: (423, 13.008), 2: (584, 12.554)}},
    #     'CAISO_CCGT1': {1: {1: (291, 7.280), 2: (484, 6.865)}},
    #     'CAISO_CCGT2': {1: {1: (129, 7.996), 2: (248, 7.381)}},
    #     'CAISO_Peaker1': {1: {1: (29, 12.904), 2: (62, 9.308)}},
    #     'CAISO_Peaker2': {1: {1: (22, 14.988), 2: (45, 11.955)}},
    #     'CAISO_Advanced_CCGT': {1: {1: (120, 10.17), 2: (600, 6.833)}},
    #     'CAISO_Aero_CT': {1: {1: (30, 17.63), 2: (100, 9.663)}},
    #     'CAISO_Reciprocating_Engine': {1: {1: (1, 10.893), 2: (5, 9.151)}},
    #     'CAISO_ST': {1: {1: (27, 17.117), 2: (337, 9.151)}},
    #     'NW_Nuclear': {1: {1: (1170, 10.907)}},
    #     'NW_Coal': {1: {1: (129, 11.259), 2: (305, 10.609)}},
    #     'NW_CCGT': {1: {1: (178, 7.721), 2: (337, 7.141)}},
    #     'NW_Peaker': {1: {1: (10, 12.500), 2: (28, 10.591)}},
    #     'SW_Nuclear': {1: {1: (1403, 10.544)}},
    #     'SW_Coal': {1: {1: (174, 11.211), 2: (414, 10.374)}},
    #     'SW_CCGT': {1: {1: (205, 7.661), 2: (372, 7.143)}},
    #     'SW_Peaker': {1: {1: (26, 14.269), 2: (71, 10.554)}},
    #     'LDWP_Nuclear': {1: {1: (152, 10.544)}},
    #     'LDWP_Coal': {1: {1: (328, 10.289), 2: (900, 9.608)}},
    #     'LDWP_CCGT': {1: {1: (123, 7.095), 2: (215, 6.995)}},
    #     'LDWP_Peaker': {1: {1: (36, 10.532), 2: (74, 9.042)}},
    #     'IID_CCGT': {1: {1: (61, 9.209), 2: (128, 7.905)}},
    #     'IID_Peaker': {1: {1: (14, 16.208), 2: (41, 12.140)}},
    #     'BANC_CCGT': {1: {1: (124, 8.037), 2: (234, 7.677)}},
    #     'BANC_Peaker': {1: {1: (16, 12.121), 2: (40, 10.392)}},
    #     'CAISO_Peaker1_Planned': {1: {1: (29, 12.904), 2: (62, 9.308)}},
    #     'CAISO_Peaker2_Planned': {1: {1: (22, 14.988), 2: (45, 11.955)}},
    #     'CAISO_CCGT1_Planned': {1: {1: (291, 7.280), 2: (484, 6.865)}},
    #     'CAISO_Reciprocating_Engine_Candidate': {1: {1: (1, 10.893), 2: (5, 9.151)}}
    # }


    project_operational_chars.update_project_hr_curves(
        io, c2, hr_names, hr_dict
    )
    io.commit()

    # ### Min stable level ### #
    project_min_stable_level = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            if rows_list[row][3 - 1] == 'CAISO_Conventional_DR':
                pass
            # If Pmax=Pmin, we're modeling as must-run, so no Pmin
            elif float(rows_list[row][17 - 1]) == 1:
                pass
            else:
                project_min_stable_level[rows_list[row][3 - 1]] = \
                    float(rows_list[row][17 - 1])

    # Add GridPath split projects
    project_min_stable_level["CAISO_Peaker1_Planned"] = \
        project_min_stable_level["CAISO_Peaker1"]
    project_min_stable_level["CAISO_Peaker2_Planned"] = \
        project_min_stable_level["CAISO_Peaker2"]
    project_min_stable_level["CAISO_CCGT1_Planned"] = \
        project_min_stable_level["CAISO_CCGT1"]
    project_min_stable_level["CAISO_Reciprocating_Engine_Candidate"] = \
        project_min_stable_level["CAISO_Reciprocating_Engine"]

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="min_stable_level",
        project_char=project_min_stable_level
    )

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                # If this is a 'must run' project, min stable level remains NULL
                if int(float(row[2])):
                    pass
                else:
                    c2.execute(
                        """UPDATE inputs_project_operational_chars
                        SET min_stable_level = {}
                        WHERE project = '{}'
                        AND project_operational_chars_scenario_id = 1;""".format(
                            row[5], row[0]
                        )
                    )

    io.commit()

    # ### Unit size ### #
    project_unit_size = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            if rows_list[row][3 - 1] == 'CAISO_Conventional_DR':
                pass
            # If Pmax=Pmin, we're modeling as must-run, so don't need unit size
            elif float(rows_list[row][17 - 1]) == 1:
                pass
            else:
                project_unit_size[rows_list[row][3 - 1]] = \
                    float(rows_list[row][15 - 1])

    # Add GridPath split projects
    project_unit_size["CAISO_Peaker1_Planned"] = \
        project_unit_size["CAISO_Peaker1"]
    project_unit_size["CAISO_Peaker2_Planned"] = \
        project_unit_size["CAISO_Peaker2"]
    project_unit_size["CAISO_CCGT1_Planned"] = \
        project_unit_size["CAISO_CCGT1"]
    project_unit_size["CAISO_Reciprocating_Engine_Candidate"] = \
        project_unit_size["CAISO_Reciprocating_Engine"]

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="unit_size_mw",
        project_char=project_unit_size
    )

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                c2.execute(
                    """UPDATE inputs_project_operational_chars
                    SET unit_size_mw = {}
                    WHERE project = '{}'
                    AND project_operational_chars_scenario_id = 1;""".format(
                        row[6], row[0]
                    )
                )
    io.commit()

    # ### Min up and down time ### #
    project_min_up_down_time = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"),
              "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            if rows_list[row][3 - 1] == 'CAISO_Conventional_DR':
                pass
            # If Pmax=Pmin, we're modeling as must-run
            elif float(rows_list[row][17 - 1]) == 1:
                pass
            # CAISO_Nuclear modeled as always-on
            elif rows_list[row][3 - 1] == 'CAISO_Nuclear':
                pass
            else:
                project_min_up_down_time[rows_list[row][3 - 1]] = \
                    float(rows_list[row][21 - 1])

        # Add GridPath split projects
        project_min_up_down_time["CAISO_Peaker1_Planned"] = \
            project_min_up_down_time["CAISO_Peaker1"]
        project_min_up_down_time["CAISO_Peaker2_Planned"] = \
            project_min_up_down_time["CAISO_Peaker2"]
        project_min_up_down_time["CAISO_CCGT1_Planned"] = \
            project_min_up_down_time["CAISO_CCGT1"]
        project_min_up_down_time["CAISO_Reciprocating_Engine_Candidate"] = \
            project_min_up_down_time["CAISO_Reciprocating_Engine"]

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="min_up_time_hours",
        project_char=project_min_up_down_time
    )

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="min_down_time_hours",
        project_char=project_min_up_down_time
    )

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                type = c2.execute(
                    """SELECT operational_type FROM
                    inputs_project_operational_chars
                    WHERE project_operational_chars_scenario_id = 1
                    AND project = '{}'""".format(row[0])
                ).fetchone()[0]
                if type == 'must_run':
                    pass
                else:
                    tech = c2.execute(
                        """SELECT technology FROM inputs_project_operational_chars
                        WHERE project_operational_chars_scenario_id = 1
                        AND project = '{}'""".format(row[0])
                    ).fetchone()[0]

                    c2.execute(
                        """UPDATE inputs_project_operational_chars
                        SET min_up_time_hours = {}, min_down_time_hours = {}
                        WHERE project = '{}'
                        AND project_operational_chars_scenario_id = 1;""".format(
                            6.0 if tech == 'CCGT' else 1.0,
                            6.0 if tech == 'CCGT' else 1.0,
                            row[0]
                        )
                    )

    io.commit()

    # ### Shutdown costs ### #
    # Costs in CONV_OpChar are divided by 2 in Inputs2Write
    project_startup_shutdown_cost = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            if rows_list[row][3 - 1] == 'CAISO_Conventional_DR':
                pass
            # If Pmax=Pmin, we're modeling as must-run, so no startup/shutdown
            elif float(rows_list[row][17 - 1]) == 1:
                pass
            # CAISO_Nuclear modeled as always-on, so no startup/shutdown
            elif rows_list[row][3 - 1] == 'CAISO_Nuclear':
                pass
            else:
                project_startup_shutdown_cost[rows_list[row][3 - 1]] = \
                    float(rows_list[row][22 - 1]) / 2.0

    # Add GridPath split projects
    project_startup_shutdown_cost["CAISO_Peaker1_Planned"] = \
        project_startup_shutdown_cost["CAISO_Peaker1"]
    project_startup_shutdown_cost["CAISO_Peaker2_Planned"] = \
        project_startup_shutdown_cost["CAISO_Peaker2"]
    project_startup_shutdown_cost["CAISO_CCGT1_Planned"] = \
        project_startup_shutdown_cost["CAISO_CCGT1"]
    project_startup_shutdown_cost["CAISO_Reciprocating_Engine_Candidate"] = \
        project_startup_shutdown_cost["CAISO_Reciprocating_Engine"]

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="shutdown_cost_per_mw",
        project_char=project_startup_shutdown_cost
    )

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                # If this is a 'must run' baseload project, no startup/shutdown
                # costs
                if int(float(row[2])):
                    pass
                elif row[18] == 'always_on':
                    pass
                else:
                    c2.execute(
                        """UPDATE inputs_project_operational_chars
                        SET shutdown_cost_per_mw = {}
                        WHERE project = '{}'
                        AND project_operational_chars_scenario_id = 1;""".format(
                            row[8], row[0]
                        )
                    )
    io.commit()

    # ### Startup Chars ### #
    # 1. Assign a startup_chars_scenario_id to all startup projects
    startup_chars_scenario_id = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            project = rows_list[row][3 - 1]
            if project.startswith('Battery') or project == \
                    'CAISO_Conventional_DR':
                pass
            elif int(float(rows_list[row][7-1])):  # must-run flag
                pass
            # If Pmax=Pmin, we're modeling as must-run, so no startups
            elif float(rows_list[row][17-1]) == 1:
                pass
            else:
                startup_chars_scenario_id[project] = 1

    # Gridpath split projects
    startup_chars_scenario_id["CAISO_Peaker1_Planned"] = 1
    startup_chars_scenario_id["CAISO_Peaker2_Planned"] = 1
    startup_chars_scenario_id["CAISO_CCGT1_Planned"] = 1
    startup_chars_scenario_id["CAISO_Reciprocating_Engine_Candidate"] = 1

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery') or row[0] == \
                    'CAISO_Conventional_DR':
                pass
            elif row[18] == "must_run" or row[18] == "always_on":
                pass
            elif int(float(row[2])):
                pass
            else:
                startup_chars_scenario_id[row[0]] = 1

    project_operational_chars.update_project_opchar_column(
        io, c2, 1, 'startup_chars_scenario_id', startup_chars_scenario_id
    )

    # 2. Create startup_scenarios for each project in startup_names variable
    # which is used in update_project_startup_chars
    startup_names = {
        'CAISO_CHP': {1: ('default', 'default')},
        'CAISO_Nuclear': {1: ('default', 'default')},
        'CAISO_CCGT1': {1: ('default', 'default')},
        'CAISO_CCGT2': {1: ('default', 'default')},
        'CAISO_Peaker1': {1: ('default', 'default')},
        'CAISO_Peaker2': {1: ('default', 'default')},
        'CAISO_Advanced_CCGT': {1: ('default', 'default')},
        'CAISO_Aero_CT': {1: ('default', 'default')},
        'CAISO_Reciprocating_Engine': {1: ('default', 'default')},
        'CAISO_ST': {1: ('default', 'default')},
        'NW_Nuclear': {1: ('default', 'default')},
        'NW_Coal': {1: ('default', 'default')},
        'NW_CCGT': {1: ('default', 'default')},
        'NW_Peaker': {1: ('default', 'default')},
        'SW_Nuclear': {1: ('default', 'default')},
        'SW_Coal': {1: ('default', 'default')},
        'SW_CCGT': {1: ('default', 'default')},
        'SW_Peaker': {1: ('default', 'default')},
        'LDWP_Nuclear': {1: ('default', 'default')},
        'LDWP_Coal': {1: ('default', 'default')},
        'LDWP_CCGT': {1: ('default', 'default')},
        'LDWP_Peaker': {1: ('default', 'default')},
        'IID_CCGT': {1: ('default', 'default')},
        'IID_Peaker': {1: ('default', 'default')},
        'BANC_CCGT': {1: ('default', 'default')},
        'BANC_Peaker': {1: ('default', 'default')},
        'CAISO_Peaker1_Planned': {1: ('default', 'default')},
        'CAISO_Peaker2_Planned': {1: ('default', 'default')},
        'CAISO_CCGT1_Planned': {1: ('default', 'default')},
        'CAISO_Reciprocating_Engine_Candidate': {1: ('default', 'default')}
    }

    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                startup_names[row[0]] = {1: ('default', 'default')}

    # 3. The startup char dicts, based on 2 different files:
    #  one is the conv_opchar.csv, and one is gas_disagg_plants.csv
    startup_char_dict = dict()
    with open(os.path.join("cpuc_irp_data", "csvs",
                           "CONV_OpChar.csv"),
              "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            project = rows_list[row][3 - 1]
            if project.startswith('Battery') or project == \
                    'CAISO_Conventional_DR':
                pass
            elif int(float(rows_list[row][7-1])):  # must-run flag
                pass
            # If Pmax=Pmin, we're modeling as must-run, so no startups
            elif float(rows_list[row][17-1]) == 1:
                pass
            else:
                # Note: the key is the startup type id
                # Startup costs are divided by 2 in Inputs2Write
                startup_chars = {1: (float(rows_list[row][21-1]),
                                     min(float(rows_list[row][18 - 1]), 1.0),
                                     float(rows_list[row][22-1])/2.,
                                     None)}
                startup_char_dict[project] = {1: startup_chars}  # 1 is scen id

    # Add GridPath split projects
        startup_char_dict["CAISO_Peaker1_Planned"] = \
            startup_char_dict["CAISO_Peaker1"]
        startup_char_dict["CAISO_Peaker2_Planned"] = \
            startup_char_dict["CAISO_Peaker2"]
        startup_char_dict["CAISO_CCGT1_Planned"] = \
            startup_char_dict["CAISO_CCGT1"]
        startup_char_dict["CAISO_Reciprocating_Engine_Candidate"] = \
            startup_char_dict["CAISO_Reciprocating_Engine"]

    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery') or row[0] == \
                    'CAISO_Conventional_DR':
                pass
            elif row[18] == "must_run" or row[18] == "always_on":
                pass
            elif int(float(row[2])):
                pass
            else:
                # Note: make sure down_times are set already!
                # If there are none, set to zero or 1
                down_time = c2.execute(
                    """SELECT min_down_time_hours FROM 
                    inputs_project_operational_chars
                    WHERE project_operational_chars_scenario_id = 1
                    AND project = '{}'""".format(row[0])
                ).fetchone()[0]
                startup_chars = {1: (down_time,
                                     min(float(row[9]), 1),
                                     float(row[7]),
                                     None)}
                startup_char_dict[row[0]] = {1: startup_chars}

    # 4. Update the table
    project_operational_chars.update_project_startup_chars(
        io, c2, startup_names, startup_char_dict
    )
    io.commit()

    # ### Ramp rates ### #
    project_off_on_plus_ramp = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(8 - 1, 34):
            if rows_list[row][3 - 1] == 'CAISO_Conventional_DR':
                pass
            # If Pmax=Pmin, we're modeling as must-run, so no ramping
            elif float(rows_list[row][17 - 1]) == 1:
                pass
            # CAISO_Nuclear modeled as always-on, so we'll apply ramp when on
            elif rows_list[row][3 - 1] == 'CAISO_Nuclear':
                pass
            else:
                # Ramp rates greater than 1 will be set to 1, as GridPath
                # requires this param to be PercentFraction
                if float(rows_list[row][18 - 1]) > 1:
                    project_off_on_plus_ramp[rows_list[row][3 - 1]] = 1.0
                else:
                    project_off_on_plus_ramp[rows_list[row][3 - 1]] = \
                        float(rows_list[row][18 - 1])

    # Add GridPath split projects
    project_off_on_plus_ramp["CAISO_Peaker1_Planned"] = \
        project_off_on_plus_ramp["CAISO_Peaker1"]
    project_off_on_plus_ramp["CAISO_Peaker2_Planned"] = \
        project_off_on_plus_ramp["CAISO_Peaker2"]
    project_off_on_plus_ramp["CAISO_CCGT1_Planned"] = \
        project_off_on_plus_ramp["CAISO_CCGT1"]
    project_off_on_plus_ramp["CAISO_Reciprocating_Engine_Candidate"] = \
        project_off_on_plus_ramp["CAISO_Reciprocating_Engine"]

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="shutdown_plus_ramp_down_rate",
        project_char=project_off_on_plus_ramp
    )
    
    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                # If this is a 'must run' project, no ramp rates
                if int(float(row[2])):
                    pass
                elif row[18] == 'always_on':
                    pass
                else:
                    c2.execute(
                        """UPDATE inputs_project_operational_chars
                        SET shutdown_plus_ramp_down_rate = {},
                        ramp_up_when_on_rate = {},
                        ramp_down_when_on_rate = {}
                        WHERE project = '{}'
                        AND project_operational_chars_scenario_id = 1;""".format(
                            1.0 if float(row[10]) > 1.0 else row[10],
                            1.0 if float(row[22]) > 1.0 else row[22],
                            1.0 if float(row[23]) > 1.0 else row[23],
                            row[0]
                        )
                    )
    io.commit()

    # Apply 'when on' ramp rate to CAISO_Nuclear, modeled as always-on,
    # and 'CAISO_Hydro'
    project_when_on_ramp = {
        'CAISO_Nuclear': 0.18446,
        'CAISO_Hydro': 967.24/7844.85
    }
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            # If this is a 'must run' project, no ramp rates
            if row[18] == 'always_on':
                project_when_on_ramp[row[0]] = \
                    1 if float(row[9]) > 1 else float(row[9])
            else:
                pass

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="ramp_up_when_on_rate",
        project_char=project_when_on_ramp
    )
    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="ramp_down_when_on_rate",
        project_char=project_when_on_ramp
    )

    # Ramp rate limits on reserves
    # 10-min ramp rate assumed for LF up, LF down, reg up, reg down, spin
    proj_reserve_ramp_rate_limits = dict()
    for row in range(8 - 1, 34):
        if rows_list[row][3 - 1] == 'CAISO_Conventional_DR':
            pass
        # If Pmax=Pmin, we're modeling as must-run, so no reserves
        elif float(rows_list[row][17 - 1]) == 1:
            pass
        # CAISO_Nuclear modeled does not provide reserves
        elif rows_list[row][3 - 1] == 'CAISO_Nuclear':
            pass
        else:
            proj_reserve_ramp_rate_limits[rows_list[row][3 - 1]] = \
                float(rows_list[row][18 - 1]) * 1/6.0 \
                if float(rows_list[row][18 - 1]) * 1/6.0 < 1 \
                else 1

    for reserve_product in [
        "lf_reserves_up", "lf_reserves_down",
        "regulation_up", "regulation_down",
        "spinning_reserves"
    ]:
        project_operational_chars.update_project_opchar_column(
            io=io, c=c2,
            project_operational_chars_scenario_id=1,
            column="{}_ramp_rate".format(reserve_product),
            project_char=proj_reserve_ramp_rate_limits
        )

    # 8% of committed capacity assumed for frequency response
    proj_freq_resp_ramp_rate_limits = dict()
    for proj in project_off_on_plus_ramp.keys():
        proj_freq_resp_ramp_rate_limits[proj] = 0.08
    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="frequency_response_ramp_rate",
        project_char=proj_freq_resp_ramp_rate_limits
    )

    # TODO: this should be based on ramp_up_when_on rate
    # Add disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery'):
                pass
            else:
                # If this is a 'must run' project, no ramp rates
                if int(float(row[2])):
                    pass
                else:
                    c2.execute(
                        """UPDATE inputs_project_operational_chars
                        SET lf_reserves_up_ramp_rate = {},
                        lf_reserves_down_ramp_rate = {},
                        regulation_up_ramp_rate = {},
                        regulation_down_ramp_rate = {},
                        spinning_reserves_ramp_rate = {},
                        frequency_response_ramp_rate = {}
                        WHERE project = '{}'
                        AND project_operational_chars_scenario_id = 1;""".format(
                            1.0 if float(row[9])/6.0 > 1.0 else float(row[9])/6.0,
                            1.0 if float(row[9]) / 6.0 > 1.0 else float(
                                row[9]) / 6.0,
                            1.0 if float(row[9]) / 6.0 > 1.0 else float(
                                row[9]) / 6.0,
                            1.0 if float(row[9]) / 6.0 > 1.0 else float(
                                row[9]) / 6.0,
                            1.0 if float(row[9]) / 6.0 > 1.0 else float(
                                row[9]) / 6.0,
                            0.08,
                            row[0]
                        )
                    )
    io.commit()

    # ### Derate for renewables providing LF down reserves ### #
    lf_reserves_down_projects = [p[0] for p in c2.execute(
        """SELECT project
        FROM inputs_project_lf_reserves_down_bas
        WHERE lf_reserves_down_ba_scenario_id = 1
        AND project_lf_reserves_down_ba_scenario_id = 2;"""
    ).fetchall()]
    var_proj_lf_down_derate = dict()
    for project in lf_reserves_down_projects:
        if project_op_types[project] == 'variable':
            var_proj_lf_down_derate[project] = 0.5

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="lf_reserves_down_derate",
        project_char=var_proj_lf_down_derate
        )

    # ### Charging and discharging efficiency ### #
    proj_eff = {
        'CAISO_Existing_Pumped_Storage': 0.81**(1/2.0),
        'CAISO_Storage_Mandate': 0.85**(1/2.0),
        'CAISO_New_Pumped_Storage': 0.81**(1/2.0),
        'CAISO_New_Flow_Battery': 0.70**(1/2.0),
        'CAISO_New_Li_Battery': 0.85**(1/2.0),
        'Battery_LCR_Bay_Area': 0.85**(1/2.0),
        'Battery_LCR_Big_Creek_Ventura': 0.85**(1/2.0),
        'Battery_LCR_Fresno': 0.85**(1/2.0),
        'Battery_LCR_Humboldt': 0.85**(1/2.0),
        'Battery_LCR_Kern': 0.85**(1/2.0),
        'Battery_LCR_LA_Basin': 0.85**(1/2.0),
        'Battery_LCR_NCNB': 0.85**(1/2.0),
        'Battery_LCR_San_Diego': 0.85**(1/2.0),
        'Battery_LCR_Sierra': 0.85**(1/2.0),
        'Battery_LCR_Stockton': 0.85**(1/2.0)
    }

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="charging_efficiency",
        project_char=proj_eff
    )

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="discharging_efficiency",
        project_char=proj_eff
    )

    # ### Minimum duration ### #
    # LCR-eligible batteries must have at least 4 hours of duration
    proj_min_dur = {
        'CAISO_Existing_Pumped_Storage': 12,
        'CAISO_Storage_Mandate': 1,
        'CAISO_New_Pumped_Storage': 12,
        'CAISO_New_Flow_Battery': 1,
        'CAISO_New_Li_Battery': 1,
        'Battery_LCR_Bay_Area': 4,
        'Battery_LCR_Big_Creek_Ventura': 4,
        'Battery_LCR_Fresno': 4,
        'Battery_LCR_Humboldt': 4,
        'Battery_LCR_Kern': 4,
        'Battery_LCR_LA_Basin': 4,
        'Battery_LCR_NCNB': 4,
        'Battery_LCR_San_Diego': 4,
        'Battery_LCR_Sierra': 4,
        'Battery_LCR_Stockton': 4
    }

    project_operational_chars.update_project_opchar_column(
        io=io, c=c2,
        project_operational_chars_scenario_id=1,
        column="minimum_duration_hours",
        project_char=proj_min_dur
    )

    # ### Variable generator profiles ### #
    project_operational_chars.\
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

    # Create operational portfolios with 6- and 8-hour LCR battery min
    # duration requirements
    project_operational_chars.make_scenario_and_insert_all_projects(
        io=io, c=c2,
        project_operational_chars_scenario_id=2,
        scenario_name='default_lcr6h',
        scenario_description='LCR batteries required 6 hours of duration; '
                             'default for rest.'
    )

    project_operational_chars.make_scenario_and_insert_all_projects(
        io=io, c=c2,
        project_operational_chars_scenario_id=3,
        scenario_name='default_lcr8h',
        scenario_description='LCR batteries required 8 hours of duration; '
                             'default for rest.'
    )

    opchar_scenario_4hmin = c2.execute(
        "SELECT * FROM inputs_project_operational_chars WHERE project_operational_chars_scenario_id = 1;"
    ).fetchall()

    # TODO: REMOVE STARTUP HERE
    for duration_min in [6, 8]:
        for row in opchar_scenario_4hmin:
            c2.execute(
                """UPDATE inputs_project_operational_chars
                SET technology = '{}',
                operational_type = '{}',
                balancing_type_project = '{}',
                variable_cost_per_mwh = {},
                fuel = {},
                heat_rate_curves_scenario_id = {},
                min_stable_level = {},
                unit_size_mw = {},
                startup_chars_scenario_id = {},
                shutdown_cost_per_mw = {},
                shutdown_plus_ramp_down_rate = {},
                ramp_up_when_on_rate = {},
                ramp_down_when_on_rate = {},
                min_up_time_hours = {},
                min_down_time_hours = {},
                charging_efficiency = {},
                discharging_efficiency = {},
                minimum_duration_hours = {},
                variable_generator_profile_scenario_id = {},
                hydro_operational_chars_scenario_id = {},
                lf_reserves_up_derate = {},
                lf_reserves_down_derate = {},
                regulation_up_derate = {},
                regulation_down_derate = {},
                frequency_response_derate = {},
                spinning_reserves_derate = {},
                lf_reserves_up_ramp_rate = {},
                lf_reserves_down_ramp_rate = {},
                regulation_up_ramp_rate = {},
                regulation_down_ramp_rate = {},
                frequency_response_ramp_rate = {},
                spinning_reserves_ramp_rate = {}
                WHERE project_operational_chars_scenario_id = {}
                AND project = '{}';""".format(
                    'NULL' if row[2] is None else row[2],
                    'NULL' if row[3] is None else row[3],
                    'NULL' if row[4] is None else row[4],
                    'NULL' if row[5] is None else row[5],
                    'NULL' if row[6] is None else "'" + row[6] + "'",
                    'NULL',
                    'NULL' if row[8] is None else row[8],
                    'NULL' if row[9] is None else row[9],
                    'NULL',
                    'NULL' if row[11] is None else row[11],
                    'NULL' if row[12] is None else row[12],
                    'NULL' if row[13] is None else row[13],
                    'NULL' if row[14] is None else row[14],
                    'NULL' if row[15] is None else row[15],
                    'NULL' if row[16] is None else row[16],
                    'NULL' if row[17] is None else row[17],
                    'NULL' if row[18] is None else row[18],
                    'NULL' if row[19] is None else row[19],
                    duration_min if row[1].startswith('Battery_LCR')
                    else ('NULL' if row[20] is None else row[20]),
                    'NULL' if row[21] is None else row[21],
                    'NULL' if row[22] is None else row[22],
                    'NULL' if row[23] is None else row[23],
                    'NULL' if row[24] is None else row[24],
                    'NULL' if row[25] is None else row[25],
                    'NULL' if row[26] is None else row[26],
                    'NULL' if row[27] is None else row[27],
                    'NULL' if row[28] is None else row[28],
                    'NULL' if row[29] is None else row[29],
                    'NULL' if row[30] is None else row[30],
                    'NULL' if row[31] is None else row[31],
                    'NULL' if row[32] is None else row[32],
                    'NULL' if row[33] is None else row[33],
                    'NULL' if row[34] is None else row[34],
                    2 if duration_min == 6 else 3,
                    row[1]
                )
            )
    io.commit()


def load_project_variable_profiles():
    """
    Profiles for 'variable' generators
    :return:
    """
    td_losses = 0.0733  # for scaling customer PV profile

    proj_profile_names = dict()
    proj_tmp_profiles = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Profiles.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for column in range(5 - 1, 61):
            proj_profile_names[rows_list[26 - 1][column]] = dict()
            proj_profile_names[rows_list[26 - 1][column]][1] = \
                ("default_cap_factors",
                 "Default profiles for variable generators.")
            proj_tmp_profiles[rows_list[26 - 1][column]] = OrderedDict()
            # Add scenario and stage
            proj_tmp_profiles[rows_list[26 - 1][column]][1] = OrderedDict()
            proj_tmp_profiles[rows_list[26 - 1][column]][1][1] = OrderedDict()
            for period in [2018, 2022, 2026, 2030]:
                for row in range(27 - 1, 914):
                    day = int(float(rows_list[row][2 - 1]))
                    hour = \
                        int(float(rows_list[row][3 - 1])) + 1  # start at 1
                    tmp = period * 10**4 + day * 10**2 + hour
                    if rows_list[26 - 1][column] == 'Customer_PV':
                        proj_tmp_profiles[rows_list[26 - 1][column]][1][1][tmp] = \
                            float(rows_list[row][column]) / (1 - td_losses)
                    else:
                        proj_tmp_profiles[rows_list[26 - 1][column]][1][1][tmp] = \
                            float(rows_list[row][column])

    project_operational_chars.update_project_variable_profiles(
        io=io, c=c2,
        proj_profile_names=proj_profile_names,
        proj_tmp_profiles=proj_tmp_profiles
    )


def load_project_hydro_opchar():
    """
    Energy budget, min, and max by horizon for hydro projects
    :return:
    """
    proj_capacities = dict()
    proj_opchar_names = dict()
    proj_horizon_chars = OrderedDict()

    with open(os.path.join("cpuc_irp_data", "csvs", "HYD_OpChar.csv"), "r") as f:
        # This will make a list of lists where each (sub)list is a row and
        # each element of the (sub)list is a column
        rows_list = list(csv.reader(f))

        # Capacities
        for column in range(5 - 1, 10):
            proj_capacities[rows_list[4 - 1][column]] = \
                float(rows_list[5 -1][column])

        # Hydro budgets
        for column in range(5 - 1, 10):
            proj = rows_list[14 - 1][column]
            proj_opchar_names[proj] = dict()
            proj_opchar_names[proj][1] = (
                "default_hydro_budget_min_max",
                "Default hydro horizon energy budgets, minimum, and maximum."
            )
            proj_horizon_chars[proj] = OrderedDict()
            # Add scenario
            proj_horizon_chars[proj][1] = OrderedDict()
            # Add the "day" as the balancing horizon
            proj_horizon_chars[proj][1]["day"] = OrderedDict()
            for period in [2018, 2022, 2026, 2030]:
                for row in range(15 - 1, 51):
                    day = int(float(rows_list[row][2 - 1]))
                    horizon = period * 10 ** 2 + day
                    proj_horizon_chars[proj][1]["day"][horizon] = \
                        OrderedDict()
                    mwa = float(rows_list[row][column]) / 24  # MWh to MWa
                    proj_horizon_chars[proj][1]["day"][horizon]["period"] = \
                        period
                    proj_horizon_chars[proj][1]["day"][horizon]["avg"] = \
                        mwa / proj_capacities[proj]


        # Min
        for column in range(11 - 1, 16):
            proj = rows_list[14 - 1][column]
            for period in [2018, 2022, 2026, 2030]:
                for row in range(15 - 1, 51):
                    day = int(float(rows_list[row][2 - 1]))
                    horizon = period * 10 ** 2 + day
                    min_mw = float(rows_list[row][column])
                    proj_horizon_chars[proj][1]["day"][horizon]["min"] = \
                        min_mw / proj_capacities[proj]
        # Max
        for column in range(17 - 1, 22):
            proj = rows_list[14 - 1][column]
            for period in [2018, 2022, 2026, 2030]:
                for row in range(15 - 1, 51):
                    day = int(float(rows_list[row][2 - 1]))
                    horizon = period * 10 ** 2 + day
                    max_mw = float(rows_list[row][column])
                    proj_horizon_chars[proj][1]["day"][horizon]["max"] = \
                        max_mw / proj_capacities[proj]

    # Insert data
    project_operational_chars.update_project_hydro_opchar(
            io=io, c=c2,
            proj_opchar_names=proj_opchar_names,
            proj_horizon_chars=proj_horizon_chars
    )


def load_fuels():
    """
    Fuels and CO2 intensity
    :return:
    """
    fuel_chars = dict()

    with open(os.path.join("cpuc_irp_data", "csvs", "Sys_Fuel_Costs.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        # CO2 intensity
        for row in range(8-1, 13):  # skip Conventional_DR, not fuel in GP
            fuel_chars[rows_list[row][3-1]] = rows_list[row][4-1]

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
    3 base price scenarios X 4 carbon cost scenarios
    :return:
    """
    scenario_fuel_month_prices = OrderedDict()

    with open(os.path.join("cpuc_irp_data", "csvs", "Sys_Fuel_Costs.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # Get the monthly price shape
        fuel_month_price_shape = OrderedDict()
        for row in range(8 - 1, 13):
            fuel_month_price_shape[rows_list[row][8 - 1]] = OrderedDict()
            for column in range(12 - 1, 23):
                fuel_month_price_shape[
                    rows_list[row][8 - 1]
                ][
                    int(float(rows_list[7 - 1][column]))
                ] = \
                    float(rows_list[row][column])

        # Get the base fuel prices by scenario
        scenario_fuel_year_base_price = OrderedDict()
        for base_scenario in ["Low", "Mid", "High"]:
            scenario_fuel_year_base_price[base_scenario] = OrderedDict()
        # Low
        for row in range(19 - 1, 24):
            scenario_fuel_year_base_price["Low"][rows_list[row][3 - 1]] = \
                OrderedDict()
            for column in range(8 - 1, 43):
                scenario_fuel_year_base_price[
                    "Low"
                ][
                    rows_list[row][3 - 1]
                ][
                    int(float(rows_list[16 - 1][column]))
                ] = \
                    float(rows_list[row][column])
        # Mid
        for row in range(27 - 1, 32):
            scenario_fuel_year_base_price["Mid"][
                rows_list[row][3 - 1]] = \
                OrderedDict()
            for column in range(8 - 1, 43):
                scenario_fuel_year_base_price[
                    "Mid"
                ][
                    rows_list[row][3 - 1]
                ][
                    int(float(rows_list[16 - 1][column]))
                ] = \
                    float(rows_list[row][column])
        # High
        for row in range(35 - 1, 40):
            scenario_fuel_year_base_price["High"][
                rows_list[row][3 - 1]] = \
                OrderedDict()
            for column in range(8 - 1, 43):
                scenario_fuel_year_base_price[
                    "High"
                ][
                    rows_list[row][3 - 1]
                ][
                    int(float(rows_list[16 - 1][column]))
                ] = \
                    float(rows_list[row][column])

        # Get the carbon adder by carbon price scenario
        scenario_year_carbon_price = OrderedDict()
        for row in range(53 - 1, 56):
            scenario_year_carbon_price[rows_list[row][3 - 1]] = OrderedDict()
            for column in range(8 - 1, 43):
                scenario_year_carbon_price[
                    rows_list[row][3 - 1]
                ][
                    int(float(rows_list[51 - 1][column]))
                ] = \
                    float(rows_list[row][column])

        # Get the fuel CO2 intensity for prices
        fuel_co2_intensities = dict()
        for row in range(8-1, 13):  # skip Conventional_DR, not fuel in GP
            fuel_co2_intensities[rows_list[row][3-1]] = \
                float(rows_list[row][5-1])

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
    for base_price_scenario in ["Low", "Mid", "High"]:
        for carbon_price_scenario in ["Low", "Mid", "High", "Zero"]:
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


# def load_project_portfolios():
#     """
#     Project portfolios
#     :return:
#     """
#
#     # ### Scenarios with aggregated gas ### #
#     # Conventional projects first
#     with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
#         i2w_rows_list = list(csv.reader(f))
#
#         # No retirement scenarios conventional projects
#         no_ret_conv_prj_cap_types = OrderedDict()
#         for row in range(3 - 1, 94):
#             # GridPath needs two reciprocating engine projects (existing and
#             # candidate)
#             # We'll also split the CAISO_Peaker1, CAISO_Peaker2,
#             # and CAISO_CCGT1 fleets into existing and planned for the
#             # retirement cases (only existing allowed to retire,
#             # not planned, i.e. capacity added during the study period)
#             if i2w_rows_list[row][91 - 1] == 'CAISO_Reciprocating_Engine':
#                 no_ret_conv_prj_cap_types[
#                     i2w_rows_list[row][91 - 1]
#                  ] = "existing_gen_no_economic_retirement"
#                 no_ret_conv_prj_cap_types[
#                     i2w_rows_list[row][91 - 1] + "_Candidate"
#                     ] = "new_build_generator"
#             elif i2w_rows_list[row][91 - 1] in [
#                 'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
#             ]:
#                 for tail in ["", "_Planned"]:
#                     no_ret_conv_prj_cap_types[
#                         i2w_rows_list[row][91 - 1] + tail
#                         ] = "existing_gen_no_economic_retirement"
#             # Other candidate resources (Reciprocating_Engine_Candidate treated
#             # above)
#             elif i2w_rows_list[row][91 - 1] in [
#                 "CAISO_Advanced_CCGT",
#                 "CAISO_Aero_CT",
#                 "CAISO_Shed_DR_Tranche1",
#                 "CAISO_Shed_DR_Tranche2",
#                 "CAISO_Shed_DR_Tranche3",
#                 "CAISO_Shed_DR_Tranche4",
#                 "CAISO_Shed_DR_Tranche5",
#                 "CAISO_Shed_DR_Tranche6",
#                 "CAISO_Shed_DR_Tranche7",
#                 "CAISO_Shed_DR_Tranche8"
#             ]:
#                 no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "new_build_generator"
#             # Storage projects get a different capacity_type
#             # Existing storage
#             elif i2w_rows_list[row][91 - 1] in [
#                 "CAISO_Existing_Pumped_Storage",
#                 "CAISO_Storage_Mandate"
#             ]:
#                 no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "storage_specified_no_economic_retirement"
#             # New storage
#             elif i2w_rows_list[row][91 - 1] in [
#                 "CAISO_New_Pumped_Storage",
#                 "CAISO_New_Flow_Battery",
#                 "CAISO_New_Li_Battery"
#             ]:
#                 no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "new_build_storage"
#             else:
#                 no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "existing_gen_no_economic_retirement"
#
#         # Retirement scenarios existing projects
#         # Allowed retirements: CAISO_CCGT1, CAISO_CCGT1,
#         # CAISO_Peaker1_Existing, CAISO_Peaker2,
#         # CAISO_Reciprocating_Engine_Existing, CAISO_ST
#         ret_conv_prj_cap_types = OrderedDict()
#         for row in range(3 - 1, 94):
#             # GridPath needs two reciprocating engine projects (existing and
#             # candidate)
#             # We'll also split the CAISO_Peaker1, CAISO_Peaker2,
#             # and CAISO_CCGT1 fleets into existing and
#             # planned for the retirement cases (only existing allowed to
#             # retire, not what's added during the study period)
#             if i2w_rows_list[row][91 - 1] == 'CAISO_Reciprocating_Engine':
#                 ret_conv_prj_cap_types[
#                     i2w_rows_list[row][91 - 1]
#                  ] = "existing_gen_linear_economic_retirement"
#                 ret_conv_prj_cap_types[
#                     i2w_rows_list[row][91 - 1] + "_Candidate"
#                     ] = "new_build_generator"
#             elif i2w_rows_list[row][91 - 1] in [
#                 'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
#             ]:
#                 ret_conv_prj_cap_types[
#                     i2w_rows_list[row][91 - 1]
#                  ] = "existing_gen_linear_economic_retirement"
#                 ret_conv_prj_cap_types[
#                     i2w_rows_list[row][91 - 1] + "_Planned"
#                     ] = "existing_gen_no_economic_retirement"
#             elif i2w_rows_list[row][91 - 1] in [
#                 "CAISO_CCGT2", "CAISO_ST", "CAISO_CHP"
#             ]:
#                 ret_conv_prj_cap_types[
#                     i2w_rows_list[row][91 - 1]
#                     ] = "existing_gen_linear_economic_retirement"
#             # Other candidate resources (Reciprocating_Engine_New treated
#             # above)
#             elif i2w_rows_list[row][91 - 1] in [
#                 "CAISO_Advanced_CCGT",
#                 "CAISO_Aero_CT",
#                 "CAISO_Shed_DR_Tranche1",
#                 "CAISO_Shed_DR_Tranche2",
#                 "CAISO_Shed_DR_Tranche3",
#                 "CAISO_Shed_DR_Tranche4",
#                 "CAISO_Shed_DR_Tranche5",
#                 "CAISO_Shed_DR_Tranche6",
#                 "CAISO_Shed_DR_Tranche7",
#                 "CAISO_Shed_DR_Tranche8"
#             ]:
#                 ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "new_build_generator"
#             # Storage projects get a different capacity_type
#             # Existing storage
#             elif i2w_rows_list[row][91 - 1] in [
#                 "CAISO_Existing_Pumped_Storage",
#                 "CAISO_Storage_Mandate"
#             ]:
#                 ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "storage_specified_no_economic_retirement"
#             # New storage
#             elif i2w_rows_list[row][91 - 1] in [
#                 "CAISO_New_Pumped_Storage",
#                 "CAISO_New_Flow_Battery",
#                 "CAISO_New_Li_Battery"
#             ]:
#                 ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "new_build_storage"
#             else:
#                 ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
#                     "existing_gen_no_economic_retirement"
#
#     # Peaker-only retirement scenarios
#     ret_peaker_cap_types = OrderedDict()
#     for proj in no_ret_conv_prj_cap_types.keys():
#         ret_peaker_cap_types[proj] = \
#             no_ret_conv_prj_cap_types[proj] \
#             if proj not in ['CAISO_Peaker1', 'CAISO_Peaker2'] \
#             else 'existing_gen_linear_economic_retirement'
#
#     # Peaker and CCGT retirement scenarios
#     ret_peaker_ccgt_cap_types = OrderedDict()
#     for proj in no_ret_conv_prj_cap_types.keys():
#         ret_peaker_ccgt_cap_types[proj] = \
#             no_ret_conv_prj_cap_types[proj] \
#             if proj not in ['CAISO_Peaker1', 'CAISO_Peaker2',
#                             'CAISO_CCGT1', 'CAISO_CCGT2'] \
#             else 'existing_gen_linear_economic_retirement'
#
#     # Candidate renewable projects
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
#         cand_ren_rows_list = list(csv.reader(f))
#
#         cand_ren_oos_none = OrderedDict()
#         cand_ren_oos_ex_tx_only = OrderedDict()
#         cand_ren_oos_ex_and_new_tx = OrderedDict()
#
#         for row in range(11 - 1, 52):
#             if int(float(cand_ren_rows_list[row][12 - 1])):
#                 cand_ren_oos_none[cand_ren_rows_list[row][3 - 1]] = \
#                     "new_build_generator"
#             else:
#                 pass
#
#             if int(float(cand_ren_rows_list[row][13 - 1])):
#                 cand_ren_oos_ex_tx_only[cand_ren_rows_list[row][3 - 1]] = \
#                     "new_build_generator"
#             else:
#                 pass
#
#             if int(float(cand_ren_rows_list[row][14 - 1])):
#                 cand_ren_oos_ex_and_new_tx[cand_ren_rows_list[row][3 - 1]] = \
#                     "new_build_generator"
#             else:
#                 pass
#
#     # Final portfolios
#     no_ret_oos_none = no_ret_conv_prj_cap_types.copy()
#     no_ret_oos_none.update(cand_ren_oos_none)
#
#     no_ret_oos_ex_tx = no_ret_conv_prj_cap_types.copy()
#     no_ret_oos_ex_tx.update(cand_ren_oos_ex_tx_only)
#
#     no_ret_oos_ex_and_new_tx = no_ret_conv_prj_cap_types.copy()
#     no_ret_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)
#
#     ret_oos_none = ret_conv_prj_cap_types.copy()
#     ret_oos_none.update(cand_ren_oos_none)
#
#     ret_oos_ex_tx = ret_conv_prj_cap_types.copy()
#     ret_oos_ex_tx.update(cand_ren_oos_ex_tx_only)
#
#     ret_oos_ex_and_new_tx = ret_conv_prj_cap_types.copy()
#     ret_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)
#
#     ret_peakers_oos_none = ret_peaker_cap_types.copy()
#     ret_peakers_oos_none.update(cand_ren_oos_none)
#
#     ret_peakers_oos_ex_tx = ret_peaker_cap_types.copy()
#     ret_peakers_oos_ex_tx.update(cand_ren_oos_ex_tx_only)
#
#     ret_peakers_oos_ex_and_new_tx = ret_peaker_cap_types.copy()
#     ret_peakers_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)
#
#     ret_peakers_ccgt_oos_none = ret_peaker_ccgt_cap_types.copy()
#     ret_peakers_ccgt_oos_none.update(cand_ren_oos_none)
#
#     ret_peakers_ccgt_oos_ex_tx = ret_peaker_ccgt_cap_types.copy()
#     ret_peakers_ccgt_oos_ex_tx.update(cand_ren_oos_ex_tx_only)
#
#     ret_peakers_ccgt_oos_ex_and_new_tx = ret_peaker_ccgt_cap_types.copy()
#     ret_peakers_ccgt_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)
#
#     portfolios = OrderedDict()
#     portfolios["No_Retirements_No_OOS_Renewables"] = no_ret_oos_none
#     portfolios["No_Retirements_Existing_Tx_Only_Renewables"] = \
#         no_ret_oos_ex_tx
#     portfolios["No_Retirements_Existing_and_New_Tx_Renewables"] = \
#         no_ret_oos_ex_and_new_tx
#     # portfolios["Retirements_Allowed_No_OOS_Renewables"] = ret_oos_none
#     # portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables"] = \
#     #     ret_oos_ex_tx
#     # portfolios["Retirements_Allowed_Existing_and_New_Tx_Renewables"] = \
#     #     ret_oos_ex_and_new_tx
#     portfolios["Retirements_Allowed_No_OOS_Renewables"] = \
#         ret_peakers_ccgt_oos_none
#     portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables"]\
#         = ret_peakers_ccgt_oos_ex_tx
#     portfolios[
#         "Retirements_Allowed_Existing_and_New_Tx_Renewables"] = \
#         ret_peakers_ccgt_oos_ex_and_new_tx
#     portfolios["Peaker_Retirements_Allowed_No_OOS_Renewables"] = \
#         ret_peakers_oos_none
#     portfolios["Peaker_Retirements_Allowed_Existing_Tx_Only_Renewables"] = \
#         ret_peakers_oos_ex_tx
#     portfolios["Peaker_Retirements_Allowed_Existing_and_New_Tx_Renewables"] = \
#         ret_peakers_oos_ex_and_new_tx
#
#     # # Add shiftable load project
#     # no_ret_oos_none_dr = no_ret_oos_none.copy()
#     # no_ret_oos_none_dr.update(
#     #     {"Shift_DR":"new_shiftable_load_supply_curve"}
#     # )
#     # no_ret_oos_ex_tx_dr = no_ret_oos_ex_tx.copy()
#     # no_ret_oos_ex_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # no_ret_oos_ex_and_new_tx_dr = no_ret_oos_ex_and_new_tx.copy()
#     # no_ret_oos_ex_and_new_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_oos_none_dr = ret_oos_none.copy()
#     # ret_oos_none_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_oos_ex_tx_dr = ret_oos_ex_tx.copy()
#     # ret_oos_ex_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_oos_ex_and_new_tx_dr = ret_oos_ex_and_new_tx.copy()
#     # ret_oos_ex_and_new_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_peakers_oos_none_dr = ret_peakers_oos_none.copy()
#     # ret_peakers_oos_none_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_peakers_oos_ex_tx_dr = ret_peakers_oos_ex_tx.copy()
#     # ret_peakers_oos_ex_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_peakers_oos_ex_and_new_tx_dr = ret_peakers_oos_ex_and_new_tx.copy()
#     # ret_peakers_oos_ex_and_new_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_peakers_ccgt_oos_none_dr = ret_peakers_ccgt_oos_none.copy()
#     # ret_peakers_ccgt_oos_none_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_peakers_ccgt_oos_ex_tx_dr = ret_peakers_ccgt_oos_ex_tx.copy()
#     # ret_peakers_ccgt_oos_ex_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     # ret_peakers_ccgt_oos_ex_and_new_tx_dr = ret_peakers_ccgt_oos_ex_and_new_tx.copy()
#     # ret_peakers_ccgt_oos_ex_and_new_tx_dr.update(
#     #     {"Shift_DR": "new_shiftable_load_supply_curve"}
#     # )
#     #
#     # portfolios["No_Retirements_No_OOS_Renewables_Shift_DR"] = \
#     #     no_ret_oos_none_dr
#     # portfolios["No_Retirements_Existing_Tx_Only_Renewables_Shift_DR"] = \
#     #     no_ret_oos_ex_tx_dr
#     # portfolios["No_Retirements_Existing_and_New_Tx_Renewables_Shift_DR"] = \
#     #     no_ret_oos_ex_and_new_tx_dr
#     # portfolios["Retirements_Allowed_No_OOS_Renewables_Shift_DR"] = \
#     #     ret_oos_none_dr
#     # portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables_Shift_DR"] = \
#     #     ret_oos_ex_tx_dr
#     # portfolios["Retirements_Allowed_Existing_and_New_Tx_Renewables_Shift_DR"] \
#     #     = ret_oos_ex_and_new_tx_dr
#     # portfolios["Peaker_Retirements_Allowed_No_OOS_Renewables_Shift_DR"] = \
#     #     ret_peakers_oos_none_dr
#     # portfolios["Peaker_Retirements_Allowed_Existing_Tx_Only_Renewables_Shift_DR"] = \
#     #     ret_peakers_oos_ex_tx_dr
#     # portfolios["Peaker_Retirements_Allowed_Existing_and_New_Tx_Renewables_Shift_DR"] = \
#     #     ret_peakers_oos_ex_and_new_tx_dr
#     # portfolios["Retirements_Allowed_No_OOS_Renewables_Shift_DR"] = \
#     #     ret_peakers_ccgt_oos_none_dr
#     # portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables_Shift_DR"]\
#     #     = ret_peakers_ccgt_oos_ex_tx_dr
#     # portfolios[
#     #     "Retirements_Allowed_Existing_and_New_Tx_Renewables_Shift_DR"] = \
#     #     ret_peakers_ccgt_oos_ex_and_new_tx_dr
#
#     # ### Portfolios with disaggregated gas projects ### #
#     # Start with aggregated gas scenarios
#     no_hyb_oos_none_gas_disagg = no_ret_oos_none.copy()
#     no_hyb_oos_ex_tx_gas_disagg = no_ret_oos_ex_tx.copy()
#     no_hyb_oos_ex_and_new_tx_gas_disagg = no_ret_oos_ex_and_new_tx.copy()
#     hyb_oos_none_gas_disagg = ret_peakers_ccgt_oos_none.copy()
#     hyb_peakers_ccgt_oos_ex_tx_gas_disagg = ret_peakers_ccgt_oos_ex_tx.copy()
#     hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg = \
#         ret_peakers_ccgt_oos_ex_and_new_tx.copy()
#
#     # No retirement gas disaggregation scenarios
#     # First will need to remove the aggregated projects from our dictionaries
#     for agg_proj in ['CAISO_CCGT1', 'CAISO_CCGT2', 'CAISO_Peaker1',
#                      'CAISO_Peaker2', 'CAISO_CCGT1_Planned',
#                      'CAISO_Peaker1_Planned', 'CAISO_Peaker2_Planned']:
#         no_hyb_oos_none_gas_disagg.pop(agg_proj)
#         no_hyb_oos_ex_tx_gas_disagg.pop(agg_proj)
#         no_hyb_oos_ex_and_new_tx_gas_disagg.pop(agg_proj)
#         hyb_oos_none_gas_disagg.pop(agg_proj)
#         hyb_peakers_ccgt_oos_ex_tx_gas_disagg.pop(agg_proj)
#         hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.pop(agg_proj)
#
#     # Then add the disaggregated plants
#     with open(os.path.join("cpuc_irp_data", "gridpath_specific",
#                            "gas_disagg_capacity_types.csv"), "r") as f:
#         reader = csv.reader(f, delimiter=",")
#         next(reader)
#
#         for row in reader:
#             # Skip hybridized resources for no-hybridization scenarios
#             if str(row[0][-3:]) == "Hyb":
#                 pass
#             else:
#                 no_hyb_oos_none_gas_disagg[row[0]] = \
#                     'existing_gen_no_economic_retirement'
#                 no_hyb_oos_ex_tx_gas_disagg[row[0]] = \
#                     'existing_gen_no_economic_retirement'
#                 no_hyb_oos_ex_and_new_tx_gas_disagg[row[0]] = \
#                     'existing_gen_no_economic_retirement'
#             # Add all resources for hybridization scenarios
#             hyb_oos_none_gas_disagg[row[0]] = row[1]
#             hyb_peakers_ccgt_oos_ex_tx_gas_disagg[row[0]] = row[1]
#             hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg[row[0]] = row[1]
#
#     portfolios["No_Hyb_No_OOS_Renewables_Gas_Disagg"] = \
#         no_hyb_oos_none_gas_disagg
#     portfolios["No_Hyb_Existing_Tx_Only_Renewables_Gas_Disagg"] = \
#         no_hyb_oos_ex_tx_gas_disagg
#     portfolios["No_Hyb_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
#         no_hyb_oos_ex_and_new_tx_gas_disagg
#     portfolios["Hyb_Allowed_No_OOS_Renewables_Gas_Disagg"] = \
#         hyb_oos_none_gas_disagg
#     portfolios["Hyb_Allowed_Existing_Tx_Only_Renewables_Gas_Disagg"]\
#         = hyb_peakers_ccgt_oos_ex_tx_gas_disagg
#     portfolios[
#         "Hyb_Allowed_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
#         hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg
#
#     # CCGT analysis portfolios
#     ccgt_analysis_no_hyb_oos_none_gas_disagg = \
#         no_hyb_oos_none_gas_disagg.copy()
#     ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg = \
#         no_hyb_oos_ex_tx_gas_disagg.copy()
#     ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg = \
#         no_hyb_oos_ex_and_new_tx_gas_disagg.copy()
#     ccgt_analysis_hyb_oos_none_gas_disagg = \
#         hyb_oos_none_gas_disagg.copy()
#     ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg = \
#         hyb_peakers_ccgt_oos_ex_tx_gas_disagg.copy()
#     ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg = \
#         hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.copy()
#
#     # Remove the two aggregated CCGT fleets that include the hybrid CCGT
#     # candidates and add the new aggregated fleets
#     for agg_proj in ['CCGT_AverageInflex', 'CCGT_AverageFlex']:
#         ccgt_analysis_no_hyb_oos_none_gas_disagg.pop(agg_proj)
#         ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg.pop(agg_proj)
#         ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg.pop(agg_proj)
#         ccgt_analysis_hyb_oos_none_gas_disagg.pop(agg_proj)
#         ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg.pop(agg_proj)
#         ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.pop(
#             agg_proj)
#     # Then add the new plants
#     with open(os.path.join("cpuc_irp_data", "gridpath_specific",
#                            "plants_to_add_for_ccgt_analysis.csv"), "r") as f:
#         reader = csv.reader(f, delimiter=",")
#         next(reader)
#
#         for row in reader:
#             # Skip hybridized resources for no-hybridization scenarios
#             if str(row[0][-3:]) == "Hyb":
#                 pass
#             else:
#                 ccgt_analysis_no_hyb_oos_none_gas_disagg[row[0]] = \
#                     'existing_gen_no_economic_retirement'
#                 ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg[row[0]] = \
#                     'existing_gen_no_economic_retirement'
#                 ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg[row[0]] = \
#                     'existing_gen_no_economic_retirement'
#             # Add all resources for hybridization scenarios
#             ccgt_analysis_hyb_oos_none_gas_disagg[row[0]] = row[19]
#             ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg[row[0]] = \
#                 row[19]
#             ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg[row[
#                 0]] = row[19]
#
#     portfolios["CCGT_Analysis_No_Hyb_No_OOS_Renewables_Gas_Disagg"] = \
#         ccgt_analysis_no_hyb_oos_none_gas_disagg
#     portfolios["CCGT_Analysis_No_Hyb_Existing_Tx_Only_Renewables_Gas_Disagg"] = \
#         ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg
#     portfolios["CCGT_Analysis_No_Hyb_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
#         ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg
#     portfolios["CCGT_Analysis_Hyb_Allowed_No_OOS_Renewables_Gas_Disagg"] = \
#         ccgt_analysis_hyb_oos_none_gas_disagg
#     portfolios["CCGT_Analysis_Hyb_Allowed_Existing_Tx_Only_Renewables_Gas_Disagg"]\
#         = ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg
#     portfolios[
#         "CCGT_Analysis_Hyb_Allowed_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
#         ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg
#
#     # Retirements for CCGT retirements
#     ret_ccgt_analysis_no_hyb_oos_none_gas_disagg = \
#         ccgt_analysis_no_hyb_oos_none_gas_disagg.copy()
#     ret_ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg = \
#         ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg.copy()
#     ret_ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg = \
#         ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg.copy()
#     ret_ccgt_analysis_hyb_oos_none_gas_disagg = \
#         ccgt_analysis_hyb_oos_none_gas_disagg.copy()
#     ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg = \
#         ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg.copy()
#     ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg = \
#         ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.copy()
#
#     # Remove the retired projects (will need to adjust capacities for
#     # retirements part of the aggregated fleets
#     for agg_proj in ['Mandalay', 'GilroyPeaker1', 'GilroyPeaker2',
#                      'GilroyPeaker3', 'Henrietta1', 'Henrietta2',
#                      'GooseHaven', 'Lambie']:
#         ret_ccgt_analysis_no_hyb_oos_none_gas_disagg.pop(agg_proj)
#         ret_ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg.pop(agg_proj)
#         ret_ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg.pop(agg_proj)
#         ret_ccgt_analysis_hyb_oos_none_gas_disagg.pop(agg_proj)
#         ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg.pop(agg_proj)
#         ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.pop(
#             agg_proj)
#
#     for agg_proj in ['GilroyPeaker1_Hyb', 'GilroyPeaker2_Hyb',
#                      'GilroyPeaker3_Hyb', 'Henrietta1_Hyb', 'Henrietta2_Hyb',
#                      'GooseHaven_Hyb', 'Lambie_Hyb']:
#         ret_ccgt_analysis_hyb_oos_none_gas_disagg.pop(agg_proj)
#         ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg.pop(agg_proj)
#         ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.pop(
#             agg_proj)
#
#     portfolios["Ret_CCGT_Analysis_No_Hyb_No_OOS_Renewables_Gas_Disagg"] = \
#         ret_ccgt_analysis_no_hyb_oos_none_gas_disagg
#     portfolios["Ret_CCGT_Analysis_No_Hyb_Existing_Tx_Only_Renewables_Gas_Disagg"] = \
#         ret_ccgt_analysis_no_hyb_oos_ex_tx_gas_disagg
#     portfolios["Ret_CCGT_Analysis_No_Hyb_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
#         ret_ccgt_analysis_no_hyb_oos_ex_and_new_tx_gas_disagg
#     portfolios["Ret_CCGT_Analysis_Hyb_Allowed_No_OOS_Renewables_Gas_Disagg"] = \
#         ret_ccgt_analysis_hyb_oos_none_gas_disagg
#     portfolios["Ret_CCGT_Analysis_Hyb_Allowed_Existing_Tx_Only_Renewables_Gas_Disagg"]\
#         = ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_tx_gas_disagg
#     portfolios[
#         "Ret_CCGT_Analysis_Hyb_Allowed_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
#         ret_ccgt_analysis_hyb_peakers_ccgt_oos_ex_and_new_tx_gas_disagg
#
#     scenario_id = 1
#     for portfolio in portfolios.keys():
#         project_portfolios.update_project_portfolios(
#             io=io, c=c2,
#             project_portfolio_scenario_id=scenario_id,
#             scenario_name=portfolio,
#             scenario_description=portfolio.replace("_", " "),
#             project_cap_types=portfolios[portfolio]
#         )
#         scenario_id += 1

def load_project_portfolios():
    """
    Project portfolios
    :return:
    """

    # ### Scenarios with aggregated gas ### #
    # Conventional projects first
    with open(os.path.join("cpuc_irp_data", "csvs", "Inputs2Write.csv"), "r") as f:
        i2w_rows_list = list(csv.reader(f))

        # No retirement scenarios conventional projects
        no_ret_conv_prj_cap_types = OrderedDict()
        for row in range(3 - 1, 94):
            # GridPath needs two reciprocating engine projects (existing and
            # candidate)
            # We'll also split the CAISO_Peaker1, CAISO_Peaker2,
            # and CAISO_CCGT1 fleets into existing and planned for the
            # retirement cases (only existing allowed to retire,
            # not planned, i.e. capacity added during the study period)
            if i2w_rows_list[row][91 - 1] == 'CAISO_Reciprocating_Engine':
                no_ret_conv_prj_cap_types[
                    i2w_rows_list[row][91 - 1]
                ] = "existing_gen_no_economic_retirement"
                no_ret_conv_prj_cap_types[
                    i2w_rows_list[row][91 - 1] + "_Candidate"
                    ] = "new_build_generator"
            elif i2w_rows_list[row][91 - 1] in [
                'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
            ]:
                for tail in ["", "_Planned"]:
                    no_ret_conv_prj_cap_types[
                        i2w_rows_list[row][91 - 1] + tail
                        ] = "existing_gen_no_economic_retirement"
            # Other candidate resources (Reciprocating_Engine_Candidate treated
            # above)
            elif i2w_rows_list[row][91 - 1] in [
                "CAISO_Advanced_CCGT",
                "CAISO_Aero_CT",
                "CAISO_Shed_DR_Tranche1",
                "CAISO_Shed_DR_Tranche2",
                "CAISO_Shed_DR_Tranche3",
                "CAISO_Shed_DR_Tranche4",
                "CAISO_Shed_DR_Tranche5",
                "CAISO_Shed_DR_Tranche6",
                "CAISO_Shed_DR_Tranche7",
                "CAISO_Shed_DR_Tranche8"
            ]:
                no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "new_build_generator"
            # Storage projects get a different capacity_type
            # Existing storage
            elif i2w_rows_list[row][91 - 1] in [
                "CAISO_Existing_Pumped_Storage",
                "CAISO_Storage_Mandate"
            ]:
                no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "storage_specified_no_economic_retirement"
            # New storage
            elif i2w_rows_list[row][91 - 1] in [
                "CAISO_New_Pumped_Storage",
                "CAISO_New_Flow_Battery",
                "CAISO_New_Li_Battery"
            ]:
                no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "new_build_storage"
            else:
                no_ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "existing_gen_no_economic_retirement"

        # Retirement scenarios existing projects
        # Allowed retirements: CAISO_CCGT1, CAISO_CCGT1,
        # CAISO_Peaker1_Existing, CAISO_Peaker2,
        # CAISO_Reciprocating_Engine_Existing, CAISO_ST
        ret_conv_prj_cap_types = OrderedDict()
        for row in range(3 - 1, 94):
            # GridPath needs two reciprocating engine projects (existing and
            # candidate)
            # We'll also split the CAISO_Peaker1, CAISO_Peaker2,
            # and CAISO_CCGT1 fleets into existing and
            # planned for the retirement cases (only existing allowed to
            # retire, not what's added during the study period)
            if i2w_rows_list[row][91 - 1] == 'CAISO_Reciprocating_Engine':
                ret_conv_prj_cap_types[
                    i2w_rows_list[row][91 - 1]
                ] = "existing_gen_linear_economic_retirement"
                ret_conv_prj_cap_types[
                    i2w_rows_list[row][91 - 1] + "_Candidate"
                    ] = "new_build_generator"
            elif i2w_rows_list[row][91 - 1] in [
                'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
            ]:
                ret_conv_prj_cap_types[
                    i2w_rows_list[row][91 - 1]
                ] = "existing_gen_linear_economic_retirement"
                ret_conv_prj_cap_types[
                    i2w_rows_list[row][91 - 1] + "_Planned"
                    ] = "existing_gen_no_economic_retirement"
            elif i2w_rows_list[row][91 - 1] in [
                "CAISO_CCGT2", "CAISO_ST", "CAISO_CHP"
            ]:
                ret_conv_prj_cap_types[
                    i2w_rows_list[row][91 - 1]
                ] = "existing_gen_linear_economic_retirement"
            # Other candidate resources (Reciprocating_Engine_New treated
            # above)
            elif i2w_rows_list[row][91 - 1] in [
                "CAISO_Advanced_CCGT",
                "CAISO_Aero_CT",
                "CAISO_Shed_DR_Tranche1",
                "CAISO_Shed_DR_Tranche2",
                "CAISO_Shed_DR_Tranche3",
                "CAISO_Shed_DR_Tranche4",
                "CAISO_Shed_DR_Tranche5",
                "CAISO_Shed_DR_Tranche6",
                "CAISO_Shed_DR_Tranche7",
                "CAISO_Shed_DR_Tranche8"
            ]:
                ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "new_build_generator"
            # Storage projects get a different capacity_type
            # Existing storage
            elif i2w_rows_list[row][91 - 1] in [
                "CAISO_Existing_Pumped_Storage",
                "CAISO_Storage_Mandate"
            ]:
                ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "storage_specified_no_economic_retirement"
            # New storage
            elif i2w_rows_list[row][91 - 1] in [
                "CAISO_New_Pumped_Storage",
                "CAISO_New_Flow_Battery",
                "CAISO_New_Li_Battery"
            ]:
                ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "new_build_storage"
            else:
                ret_conv_prj_cap_types[i2w_rows_list[row][91 - 1]] = \
                    "existing_gen_no_economic_retirement"

    # Peaker-only retirement scenarios
    ret_peaker_cap_types = OrderedDict()
    for proj in no_ret_conv_prj_cap_types.keys():
        ret_peaker_cap_types[proj] = \
            no_ret_conv_prj_cap_types[proj] \
                if proj not in ['CAISO_Peaker1', 'CAISO_Peaker2'] \
                else 'existing_gen_linear_economic_retirement'

    # Peaker and CCGT retirement scenarios
    ret_peaker_ccgt_cap_types = OrderedDict()
    for proj in no_ret_conv_prj_cap_types.keys():
        ret_peaker_ccgt_cap_types[proj] = \
            no_ret_conv_prj_cap_types[proj] \
                if proj not in ['CAISO_Peaker1', 'CAISO_Peaker2',
                                'CAISO_CCGT1', 'CAISO_CCGT2'] \
                else 'existing_gen_linear_economic_retirement'

    # Candidate renewable projects
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        cand_ren_rows_list = list(csv.reader(f))

        cand_ren_oos_none = OrderedDict()
        cand_ren_oos_ex_tx_only = OrderedDict()
        cand_ren_oos_ex_and_new_tx = OrderedDict()

        for row in range(11 - 1, 52):
            if int(float(cand_ren_rows_list[row][12 - 1])):
                cand_ren_oos_none[cand_ren_rows_list[row][3 - 1]] = \
                    "new_build_generator"
            else:
                pass

            if int(float(cand_ren_rows_list[row][13 - 1])):
                cand_ren_oos_ex_tx_only[cand_ren_rows_list[row][3 - 1]] = \
                    "new_build_generator"
            else:
                pass

            if int(float(cand_ren_rows_list[row][14 - 1])):
                cand_ren_oos_ex_and_new_tx[cand_ren_rows_list[row][3 - 1]] = \
                    "new_build_generator"
            else:
                pass

    # Final portfolios
    no_ret_oos_none = no_ret_conv_prj_cap_types.copy()
    no_ret_oos_none.update(cand_ren_oos_none)

    no_ret_oos_ex_tx = no_ret_conv_prj_cap_types.copy()
    no_ret_oos_ex_tx.update(cand_ren_oos_ex_tx_only)

    no_ret_oos_ex_and_new_tx = no_ret_conv_prj_cap_types.copy()
    no_ret_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)

    ret_oos_none = ret_conv_prj_cap_types.copy()
    ret_oos_none.update(cand_ren_oos_none)

    ret_oos_ex_tx = ret_conv_prj_cap_types.copy()
    ret_oos_ex_tx.update(cand_ren_oos_ex_tx_only)

    ret_oos_ex_and_new_tx = ret_conv_prj_cap_types.copy()
    ret_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)

    ret_peakers_oos_none = ret_peaker_cap_types.copy()
    ret_peakers_oos_none.update(cand_ren_oos_none)

    ret_peakers_oos_ex_tx = ret_peaker_cap_types.copy()
    ret_peakers_oos_ex_tx.update(cand_ren_oos_ex_tx_only)

    ret_peakers_oos_ex_and_new_tx = ret_peaker_cap_types.copy()
    ret_peakers_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)

    ret_peakers_ccgt_oos_none = ret_peaker_ccgt_cap_types.copy()
    ret_peakers_ccgt_oos_none.update(cand_ren_oos_none)

    ret_peakers_ccgt_oos_ex_tx = ret_peaker_ccgt_cap_types.copy()
    ret_peakers_ccgt_oos_ex_tx.update(cand_ren_oos_ex_tx_only)

    ret_peakers_ccgt_oos_ex_and_new_tx = ret_peaker_ccgt_cap_types.copy()
    ret_peakers_ccgt_oos_ex_and_new_tx.update(cand_ren_oos_ex_and_new_tx)

    portfolios = OrderedDict()
    portfolios["No_Retirements_No_OOS_Renewables"] = no_ret_oos_none
    portfolios["No_Retirements_Existing_Tx_Only_Renewables"] = \
        no_ret_oos_ex_tx
    portfolios["No_Retirements_Existing_and_New_Tx_Renewables"] = \
        no_ret_oos_ex_and_new_tx
    # portfolios["Retirements_Allowed_No_OOS_Renewables"] = ret_oos_none
    # portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables"] = \
    #     ret_oos_ex_tx
    # portfolios["Retirements_Allowed_Existing_and_New_Tx_Renewables"] = \
    #     ret_oos_ex_and_new_tx
    portfolios["Retirements_Allowed_No_OOS_Renewables"] = \
        ret_peakers_ccgt_oos_none
    portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables"] \
        = ret_peakers_ccgt_oos_ex_tx
    portfolios[
        "Retirements_Allowed_Existing_and_New_Tx_Renewables"] = \
        ret_peakers_ccgt_oos_ex_and_new_tx
    portfolios["Peaker_Retirements_Allowed_No_OOS_Renewables"] = \
        ret_peakers_oos_none
    portfolios["Peaker_Retirements_Allowed_Existing_Tx_Only_Renewables"] = \
        ret_peakers_oos_ex_tx
    portfolios["Peaker_Retirements_Allowed_Existing_and_New_Tx_Renewables"] = \
        ret_peakers_oos_ex_and_new_tx

    # # Add shiftable load project
    # no_ret_oos_none_dr = no_ret_oos_none.copy()
    # no_ret_oos_none_dr.update(
    #     {"Shift_DR":"new_shiftable_load_supply_curve"}
    # )
    # no_ret_oos_ex_tx_dr = no_ret_oos_ex_tx.copy()
    # no_ret_oos_ex_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # no_ret_oos_ex_and_new_tx_dr = no_ret_oos_ex_and_new_tx.copy()
    # no_ret_oos_ex_and_new_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_oos_none_dr = ret_oos_none.copy()
    # ret_oos_none_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_oos_ex_tx_dr = ret_oos_ex_tx.copy()
    # ret_oos_ex_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_oos_ex_and_new_tx_dr = ret_oos_ex_and_new_tx.copy()
    # ret_oos_ex_and_new_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_peakers_oos_none_dr = ret_peakers_oos_none.copy()
    # ret_peakers_oos_none_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_peakers_oos_ex_tx_dr = ret_peakers_oos_ex_tx.copy()
    # ret_peakers_oos_ex_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_peakers_oos_ex_and_new_tx_dr = ret_peakers_oos_ex_and_new_tx.copy()
    # ret_peakers_oos_ex_and_new_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_peakers_ccgt_oos_none_dr = ret_peakers_ccgt_oos_none.copy()
    # ret_peakers_ccgt_oos_none_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_peakers_ccgt_oos_ex_tx_dr = ret_peakers_ccgt_oos_ex_tx.copy()
    # ret_peakers_ccgt_oos_ex_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    # ret_peakers_ccgt_oos_ex_and_new_tx_dr = ret_peakers_ccgt_oos_ex_and_new_tx.copy()
    # ret_peakers_ccgt_oos_ex_and_new_tx_dr.update(
    #     {"Shift_DR": "new_shiftable_load_supply_curve"}
    # )
    #
    # portfolios["No_Retirements_No_OOS_Renewables_Shift_DR"] = \
    #     no_ret_oos_none_dr
    # portfolios["No_Retirements_Existing_Tx_Only_Renewables_Shift_DR"] = \
    #     no_ret_oos_ex_tx_dr
    # portfolios["No_Retirements_Existing_and_New_Tx_Renewables_Shift_DR"] = \
    #     no_ret_oos_ex_and_new_tx_dr
    # portfolios["Retirements_Allowed_No_OOS_Renewables_Shift_DR"] = \
    #     ret_oos_none_dr
    # portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables_Shift_DR"] = \
    #     ret_oos_ex_tx_dr
    # portfolios["Retirements_Allowed_Existing_and_New_Tx_Renewables_Shift_DR"] \
    #     = ret_oos_ex_and_new_tx_dr
    # portfolios["Peaker_Retirements_Allowed_No_OOS_Renewables_Shift_DR"] = \
    #     ret_peakers_oos_none_dr
    # portfolios["Peaker_Retirements_Allowed_Existing_Tx_Only_Renewables_Shift_DR"] = \
    #     ret_peakers_oos_ex_tx_dr
    # portfolios["Peaker_Retirements_Allowed_Existing_and_New_Tx_Renewables_Shift_DR"] = \
    #     ret_peakers_oos_ex_and_new_tx_dr
    # portfolios["Retirements_Allowed_No_OOS_Renewables_Shift_DR"] = \
    #     ret_peakers_ccgt_oos_none_dr
    # portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables_Shift_DR"]\
    #     = ret_peakers_ccgt_oos_ex_tx_dr
    # portfolios[
    #     "Retirements_Allowed_Existing_and_New_Tx_Renewables_Shift_DR"] = \
    #     ret_peakers_ccgt_oos_ex_and_new_tx_dr

    # ### Portfolios with disaggregated gas projects ### #
    # Start with aggregated gas scenarios
    no_ret_oos_none_gas_disagg = no_ret_oos_none.copy()
    no_ret_oos_ex_tx_gas_disagg = no_ret_oos_ex_tx.copy()
    no_ret_oos_ex_and_new_tx_gas_disagg = no_ret_oos_ex_and_new_tx.copy()
    ret_peakers_ccgt_oos_none_gas_disagg = ret_peakers_ccgt_oos_none.copy()
    ret_peakers_ccgt_oos_ex_tx_gas_disagg = ret_peakers_ccgt_oos_ex_tx.copy()
    ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg = \
        ret_peakers_ccgt_oos_ex_and_new_tx.copy()

    # No retirement gas disaggregation scenarios
    # First will need to remove the aggregated projects from our dictionaries
    for agg_proj in ['CAISO_CCGT1', 'CAISO_CCGT2', 'CAISO_Peaker1',
                     'CAISO_Peaker2', 'CAISO_CCGT1_Planned',
                     'CAISO_Peaker1_Planned', 'CAISO_Peaker2_Planned']:
        no_ret_oos_none_gas_disagg.pop(agg_proj)
        no_ret_oos_ex_tx_gas_disagg.pop(agg_proj)
        no_ret_oos_ex_and_new_tx_gas_disagg.pop(agg_proj)
        ret_peakers_ccgt_oos_none_gas_disagg.pop(agg_proj)
        ret_peakers_ccgt_oos_ex_tx_gas_disagg.pop(agg_proj)
        ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.pop(agg_proj)

    # Then add the disaggregated plants
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_capacity_types.csv"), "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)

        for row in reader:
            # don't add LCR resources yet
            if row[0].startswith("Battery") or row[0].startswith("Advanced") or \
                    row[0].startswith("Aero"):
                pass
            else:
                no_ret_oos_none_gas_disagg[row[0]] = \
                    'existing_gen_no_economic_retirement'
                no_ret_oos_ex_tx_gas_disagg[row[0]] = \
                    'existing_gen_no_economic_retirement'
                no_ret_oos_ex_and_new_tx_gas_disagg[row[0]] = \
                    'existing_gen_no_economic_retirement'
                ret_peakers_ccgt_oos_none_gas_disagg[row[0]] = row[1]
                ret_peakers_ccgt_oos_ex_tx_gas_disagg[row[0]] = row[1]
                ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg[row[0]] = row[1]

    portfolios["No_Retirements_No_OOS_Renewables_Gas_Disagg"] = \
        no_ret_oos_none_gas_disagg
    portfolios["No_Retirements_Existing_Tx_Only_Renewables_Gas_Disagg"] = \
        no_ret_oos_ex_tx_gas_disagg
    portfolios["No_Retirements_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
        no_ret_oos_ex_and_new_tx_gas_disagg
    portfolios["Retirements_Allowed_No_OOS_Renewables_Gas_Disagg"] = \
        ret_peakers_ccgt_oos_none_gas_disagg
    portfolios["Retirements_Allowed_Existing_Tx_Only_Renewables_Gas_Disagg"] \
        = ret_peakers_ccgt_oos_ex_tx_gas_disagg
    portfolios[
        "Retirements_Allowed_Existing_and_New_Tx_Renewables_Gas_Disagg"] = \
        ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg

    # Portfolios with LCR-eligible resources
    no_ret_oos_none_gas_disagg_lcr_resources = \
        no_ret_oos_none_gas_disagg.copy()
    no_ret_oos_ex_tx_gas_disagg_lcr_resources = \
        no_ret_oos_ex_tx_gas_disagg.copy()
    no_ret_oos_ex_and_new_tx_gas_disagg_lcr_resources = \
        no_ret_oos_ex_and_new_tx_gas_disagg.copy()
    ret_peakers_ccgt_oos_none_gas_disagg_lcr_resources = \
        ret_peakers_ccgt_oos_none_gas_disagg.copy()
    ret_peakers_ccgt_oos_ex_tx_gas_disagg_lcr_resources = \
        ret_peakers_ccgt_oos_ex_tx_gas_disagg.copy()
    ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg_lcr_resources = \
        ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg.copy()

    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_capacity_types.csv"), "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)

        for row in reader:
            # Add LCR resources yet
            if row[0].startswith("Battery") or row[0].startswith("Advanced") or \
                    row[0].startswith("Aero"):
                no_ret_oos_none_gas_disagg_lcr_resources[row[0]] = row[1]
                no_ret_oos_ex_tx_gas_disagg_lcr_resources[row[0]] = row[1]
                no_ret_oos_ex_and_new_tx_gas_disagg_lcr_resources[row[0]] = \
                    row[1]
                ret_peakers_ccgt_oos_none_gas_disagg_lcr_resources[row[0]] = \
                    row[1]
                ret_peakers_ccgt_oos_ex_tx_gas_disagg_lcr_resources[row[0]] = \
                    row[1]
                ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg_lcr_resources[
                    row[0]] = row[1]
            else:
                pass
    portfolios["No_Retirements_No_OOS_Renewables_Gas_Disagg_LCR_Resources"] = \
        no_ret_oos_none_gas_disagg_lcr_resources
    portfolios[
        "No_Retirements_Existing_Tx_Only_Renewables_Gas_Disagg_LCR_Resources"] = \
        no_ret_oos_ex_tx_gas_disagg_lcr_resources
    portfolios[
        "No_Retirements_Existing_and_New_Tx_Renewables_Gas_Disagg_LCR_Resources"] = \
        no_ret_oos_ex_and_new_tx_gas_disagg_lcr_resources
    portfolios[
        "Retirements_Allowed_No_OOS_Renewables_Gas_Disagg_LCR_Resources"] = \
        ret_peakers_ccgt_oos_none_gas_disagg_lcr_resources
    portfolios[
        "Retirements_Allowed_Existing_Tx_Only_Renewables_Gas_Disagg_LCR_Resources"] \
        = ret_peakers_ccgt_oos_ex_tx_gas_disagg_lcr_resources
    portfolios[
        "Retirements_Allowed_Existing_and_New_Tx_Renewables_Gas_Disagg_LCR_Resources"] = \
        ret_peakers_ccgt_oos_ex_and_new_tx_gas_disagg_lcr_resources

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


# def load_project_capacities():
#     """
# 
#     :return:
#     """
#     # Load the data into resolve.db, as it will be easier to split the
#     # existing and planned projects
#     c1.execute(
#         """DROP TABLE IF EXISTS CONV_planned_capacity;"""
#     )
#     c1.execute(
#         """CREATE TABLE CONV_planned_capacity(
#             project VARCHAR(32),
#             proj_zone VARCHAR(16),
#             forecast_year INTEGER,
#             capacity_mw FLOAT,
#             PRIMARY KEY (project, forecast_year)
#             );"""
#     )
#     resolve.commit()
# 
#     with open(os.path.join("cpuc_irp_data", "csvs", "CONV_Baseline_no_early_retirements.csv"),
#               "r") as f:
#         rows_list = list(csv.reader(f))
# 
#         for row in range(55 - 1, 81):
#             if all(x == "" for x in rows_list[row]):
#                 pass  # skip if no data in row
#             else:
#                 for column in range(6-1, 41):
#                     c1.execute(
#                         """INSERT INTO CONV_planned_capacity 
#                         (project, proj_zone, forecast_year, capacity_mw)
#                         VALUES ('{}', '{}', {}, {});""".format(
#                             rows_list[row][3 - 1],
#                             rows_list[row][4 - 1],
#                             int(float(rows_list[6-1][column])),
#                             rows_list[row][column]
#                         )
#                     )
#         resolve.commit()
# 
#     # Split any projects that have capacity added in the future
#     # Make key table
#     c1.execute("DROP TABLE IF EXISTS PROCESSED_split_project_list;")
#     c1.execute(
#         """CREATE TABLE PROCESSED_split_project_list (
#         split_project VARCHAR(64),
#         parent_project VARCHAR(64),
#         PRIMARY KEY (split_project, parent_project)
#         );"""
#     )
#     resolve.commit()
# 
#     # Table with new capacities
#     c1.execute("DROP TABLE IF EXISTS PROCESSED_CONV_planned_capacity;")
#     c1.execute(
#         """CREATE TABLE PROCESSED_CONV_planned_capacity (
#         project VARCHAR(32),
#         proj_zone VARCHAR(16),
#         forecast_year INTEGER,
#         capacity_mw FLOAT,
#         PRIMARY KEY (project, forecast_year)
#         );"""
#     )
#     resolve.commit()
# 
#     # Insert non-CAISO project capacities
#     c1.execute(
#         """INSERT INTO PROCESSED_CONV_planned_capacity
#         (project, proj_zone, forecast_year, capacity_mw)
#         SELECT project, proj_zone, forecast_year, capacity_mw
#         FROM CONV_planned_capacity
#         WHERE project NOT LIKE 'CAISO%';"""
#     )
#     resolve.commit()
# 
#     # CAISO projects
#     # Process data -- split planned from existing apacity
#     parent_projects = c1.execute(
#         """SELECT DISTINCT project
#         FROM CONV_planned_capacity
#         WHERE project LIKE 'CAISO%';"""
#     ).fetchall()
# 
#     # Start with year 2015; insert data from CONV_planned_capacity table
#     c1.execute(
#         """INSERT INTO PROCESSED_CONV_planned_capacity
#         (project, proj_zone, forecast_year, capacity_mw)
#         SELECT project, proj_zone, forecast_year, capacity_mw
#         FROM CONV_planned_capacity
#         WHERE forecast_year = 2015
#         AND project LIKE 'CAISO%';"""
#     )
#     resolve.commit()
# 
#     for p in parent_projects:
#         for yr in range(2016, 2051):
#             # Get the project zone
#             zone = c1.execute(
#                 """SELECT proj_zone 
#                 FROM CONV_planned_capacity 
#                 WHERE project = '{}'
#                 AND forecast_year = {}""".format(p[0], yr)
#             ).fetchone()
# 
#             split_project_name = p[0] + '_Planned'
#             current_capacity = c1.execute(
#                 """SELECT capacity_mw
#                 FROM CONV_planned_capacity
#                 WHERE project = '{}'
#                 AND forecast_year = {}""".format(p[0], yr)
#             ).fetchone()
#             previous_capacity = c1.execute(
#                 """SELECT capacity_mw
#                 FROM CONV_planned_capacity
#                 WHERE project = '{}'
#                 AND forecast_year = {}""".format(p[0], yr-1)
#             ).fetchone()
# 
#             # Increase or decrease in capacity from previous year
#             delta = current_capacity[0] - previous_capacity[0]
# 
#             # If we see increase in capacity, we'll need to split
#             # the project
#             if delta > 0:
#                 # Insert split project into our key table if not already
#                 # there
#                 c1.execute(
#                     """INSERT OR IGNORE INTO PROCESSED_split_project_list
#                     (split_project, parent_project)
#                     VALUES ('{}', '{}')""".format(split_project_name, p[0])
#                 )
#                 # Insert capacities for this year into capacities table
#                 # First, get last year's capacities
#                 parent_project_previous_capacity = c1.execute(
#                     """SELECT capacity_mw
#                     FROM PROCESSED_CONV_planned_capacity
#                     WHERE project = '{}'
#                     AND forecast_year = {};""".format(p[0], yr-1)
#                 ).fetchone()
# 
#                 split_project_previous_capacity = c1.execute(
#                     """SELECT capacity_mw 
#                     FROM PROCESSED_CONV_planned_capacity
#                     WHERE project = '{}'
#                     AND forecast_year = {};""".format(
#                         split_project_name, yr-1
#                     )
#                 ).fetchone()
# 
#                 # If there was no capacity for the split project in the
#                 # previous year, it was 0
#                 if split_project_previous_capacity is None:
#                     split_project_previous_capacity = 0
#                     c1.execute(
#                         """INSERT INTO PROCESSED_CONV_planned_capacity
#                         (project, proj_zone, forecast_year, capacity_mw)
#                         VALUES ('{}', '{}', {}, {});""".format(
#                             split_project_name, zone[0], yr - 1,
#                             split_project_previous_capacity
#                         )
#                     )
#                 # Insert previous capacity + the delta for the split prj
#                 c1.execute(
#                     """INSERT INTO PROCESSED_CONV_planned_capacity
#                     (project, proj_zone, forecast_year, capacity_mw)
#                     VALUES ('{}', '{}', {}, {});""".format(
#                         split_project_name, zone[0], yr,
#                         (0 if split_project_previous_capacity == 0
#                          else split_project_previous_capacity[0]) + delta
#                     )
#                 )
#                 # Insert the previous capacity for the parent project
#                 c1.execute(
#                     """INSERT INTO PROCESSED_CONV_planned_capacity
#                     (project, proj_zone, forecast_year, capacity_mw)
#                     VALUES ('{}', '{}', {}, {});""".format(
#                         p[0], zone[0], yr,
#                         parent_project_previous_capacity[0]
#                     )
#                 )
#             # If we have the same as the previous year,
#             # simply insert the calculated previous year capacity
#             # + delta for both the parent and the previous year capacity
#             # for the split project
#             elif delta <= 0:
#                 # Insert the previous capacity for the parent project
#                 c1.execute(
#                     """INSERT INTO PROCESSED_CONV_planned_capacity
#                     (project, proj_zone, forecast_year, capacity_mw)
#                     SELECT project, proj_zone, forecast_year + 1, 
#                     capacity_mw + {}
#                     FROM PROCESSED_CONV_planned_capacity
#                     WHERE project = '{}'
#                     AND forecast_year = {};""".format(
#                         delta, p[0], yr - 1
#                     )
#                 )
# 
#                 # If no previous year capacity, we won't insert anything
#                 c1.execute(
#                     """INSERT INTO PROCESSED_CONV_planned_capacity
#                     (project, proj_zone, forecast_year, capacity_mw)
#                     SELECT project, proj_zone, forecast_year + 1, 
#                     capacity_mw
#                     FROM PROCESSED_CONV_planned_capacity
#                     WHERE project = '{}'
#                     AND forecast_year = {};""".format(
#                         split_project_name, yr - 1
#                     )
#                 )
# 
#     resolve.commit()
# 
#     # Make final dictionary
#     proj_capacities = OrderedDict()
# 
#     # Add conventional project capacities
#     distinct_conv_projects = c1.execute(
#         """SELECT DISTINCT project
#         FROM PROCESSED_CONV_planned_capacity
#         ORDER BY project;"""
#     ).fetchall()
# 
#     conv_projects = c1.execute(
#         """SELECT project, forecast_year, capacity_mw
#         FROM PROCESSED_CONV_planned_capacity
#         ORDER BY project, forecast_year;"""
#     ).fetchall()
# 
#     for project in distinct_conv_projects:
#         p = project[0]
#         proj_capacities[p] = OrderedDict()
# 
#     for proj_yr_cap in conv_projects:
#         project = proj_yr_cap[0]
#         year = proj_yr_cap[1]
#         capacity = proj_yr_cap[2]
#         proj_capacities[project][year] = (capacity, None)
# 
#     # Add CAISO_Shed_DR_Existing capacities
#     with open(os.path.join("cpuc_irp_data", "csvs", "DR_Conventional.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#         proj_capacities['CAISO_Shed_DR_Existing'] = OrderedDict()
#         for clmn in range(10 - 1, 45):
#             dr_yr = int(float(rows_list[31 - 1][clmn]))
#             proj_capacities['CAISO_Shed_DR_Existing'][dr_yr] = \
#                 (float(rows_list[35 - 1][clmn]), None)
# 
#     # Add hydro project capacities
#     with open(os.path.join("cpuc_irp_data", "csvs", "HYD_OpChar.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#         for column in range(5 - 1, 10):
#             proj_capacities[rows_list[4 - 1][column]] = OrderedDict()
#             for yr in range(2015, 2050 + 1):
#                 proj_capacities[rows_list[4 - 1][column]][yr] = \
#                     (float(rows_list[5 - 1][column]), None)
# 
#     # Add existing/planned renewable capacities
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Baseline.csv"), "r") as f:
#         # This will make a list of lists where each (sub)list is a row and
#         # each element of the (sub)list is a column
#         rows_list = list(csv.reader(f))
# 
#         # CAISO resources
#         for row in range(59 - 1, 73):
#             project = rows_list[row][3 - 1]
#             proj_capacities[project] = OrderedDict()
#             for column in range(10 - 1, 45):
#                 year = int(float(rows_list[3 - 1][column]))
#                 proj_capacities[project][year] \
#                     = (float(rows_list[row][column]), None)
# 
#          # Non-CAISO resources
#         for row in range(339 - 1, 368):
#             project = rows_list[row][3 - 1]
#             proj_capacities[project] = OrderedDict()
#             for column in range(10 - 1, 45):
#                 year = int(float(rows_list[3 - 1][column]))
#                 proj_capacities[project][year] \
#                     = (float(rows_list[row][column]), None)
# 
#     # Disaggregated gas projects
#     with open(os.path.join("cpuc_irp_data", "gridpath_specific", "gas_disagg_plants.csv"),
#               "r") as foo:
#         reader = csv.reader(foo, delimiter=",")
#         next(reader)
#         for row in reader:
#             if str(row[0][-3]) == "Hyb":
#                 pass
#             else:
#                 proj_capacities[row[0]] = OrderedDict()
#                 if row[20] == '':
#                     start_year = 2015
#                 else:
#                     start_year = int(float(row[20]))
#                 if row[21] == '':
#                     end_year = 2050 + 1
#                 else:
#                     end_year = int(float(row[21])) + 1
#                 for yr in range(2015, 2050 + 1):
#                     if start_year <= yr <= end_year:
#                         proj_capacities[row[0]][yr] = (float(row[16]), None)
#                     else:
#                         proj_capacities[row[0]][yr] = (0, None)
# 
#     # with open(os.path.join("cpuc_irp_data", "gridpath_specific",
#     #                        "plants_to_add_for_ccgt_analysis.csv"),
#     #           "r") as foo:
#     #     reader = csv.reader(foo, delimiter=",")
#     #     next(reader)
#     #     for row in reader:
#     #         if str(row[0][-3]) == "Hyb":
#     #             pass
#     #         else:
#     #             proj_capacities[row[0]] = OrderedDict()
#     #             if row[20] == '':
#     #                 start_year = 2015
#     #             else:
#     #                 start_year = int(float(row[20]))
#     #             if row[21] == '':
#     #                 end_year = 2050 + 1
#     #             else:
#     #                 end_year = int(float(row[21])) + 1
#     #             for yr in range(2015, 2050 + 1):
#     #                 if start_year <= yr <= end_year:
#     #                     proj_capacities[row[0]][yr] = (float(row[16]), None)
#     #                 else:
#     #                     proj_capacities[row[0]][yr] = (0, None)
# 
#     # Add storage capacities
#     with open(os.path.join("cpuc_irp_data", "csvs", "STOR_Inputs.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
# 
#         # Existing pumped storage
#         ps_project = rows_list[28 - 1][3 - 1]
#         proj_capacities[ps_project] = OrderedDict()
#         for column in range(6-1, 41):
#             year = int(float(rows_list[27 - 1][column]))
#             capacity = rows_list[28 - 1][column]
#             energy = rows_list[36 - 1][column]
#             proj_capacities[ps_project][year] = (capacity, energy)
# 
#         # 3 storage mandate scenarios
#         # Storage mandate assumed to have 4 hours duration
#         scenario_storage_mandate_capacities = OrderedDict()
#         sm_project = rows_list[29 - 1][3 - 1]
#         for scenario_row in range(21 - 1, 23):
#             scenario = rows_list[scenario_row][3 - 1]
#             scenario_storage_mandate_capacities[scenario] = OrderedDict()
#             scenario_storage_mandate_capacities[scenario][sm_project] = \
#                 OrderedDict()
#             for column in range(6-1, 41):
#                 sm_year = int(float(rows_list[20 - 1][column]))
#                 sm_capacity = float(rows_list[scenario_row][column])
#                 sm_energy = 4 * sm_capacity
#                 scenario_storage_mandate_capacities[
#                     scenario
#                 ][
#                     sm_project
#                 ][
#                     sm_year] = (sm_capacity, sm_energy)
# 
#     # Add Customer_PV capacities
#     # 4 scenarios
#     # We'll need the shape capacity factor to convert to MW from GWh
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Profiles.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#         customer_pv_cap_factor = float(rows_list[23 - 1][42 - 1])
# 
#     with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_Forecast.csv"), "r") as f:
#         scenario_customer_pv_capacities = OrderedDict()
#         rows_list = list(csv.reader(f))
#         for scenario_row in range(29 - 1, 32):
#             custpv_scenario = rows_list[scenario_row][3 - 1]
#             scenario_customer_pv_capacities[custpv_scenario] = \
#                 OrderedDict()
#             scenario_customer_pv_capacities[custpv_scenario]['Customer_PV'
#             ] = OrderedDict()
#             for column in range(10 - 1, 45):
#                 custpv_year = int(float(rows_list[28 - 1][column]))
#                 custpv_capacity = \
#                     float(rows_list[scenario_row][column]) * 1000 / 8760 / \
#                     customer_pv_cap_factor
#                 scenario_customer_pv_capacities[
#                     custpv_scenario
#                 ][
#                     "Customer_PV"
#                 ][
#                     custpv_year] = \
#                     (custpv_capacity, None)    
#     
#     # Make final 12 scenarios (base X 3 storage mandate cases X 4 customer
#     # PV cases)
#     scenario_projects = OrderedDict()
# 
#     # 1: no storage mandate, no additional customer PV
#     no_storage_mandate_no_addtl_customer_pv = \
#         proj_capacities.copy()
#     no_storage_mandate_no_addtl_customer_pv.update(
#         scenario_storage_mandate_capacities['None']
#     )
#     no_storage_mandate_no_addtl_customer_pv.update(
#         scenario_customer_pv_capacities['No Addtl. BTM']
#     )
#     scenario_projects["No_Storage_Mandate_No_Additional_BTM_PV"] = \
#         no_storage_mandate_no_addtl_customer_pv
# 
#     # 2: no storage mandate, low additional customer PV
#     no_storage_mandate_low_customer_pv = \
#         proj_capacities.copy()
#     no_storage_mandate_low_customer_pv.update(
#         scenario_storage_mandate_capacities['None']
#     )
#     no_storage_mandate_low_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - Low PV']
#     )
#     scenario_projects["No_Storage_Mandate_Low_Additional_BTM_PV"] = \
#         no_storage_mandate_low_customer_pv
# 
#     # 3: no storage mandate, mid additional customer PV
#     no_storage_mandate_mid_customer_pv = \
#         proj_capacities.copy()
#     no_storage_mandate_mid_customer_pv.update(
#         scenario_storage_mandate_capacities['None']
#     )
#     no_storage_mandate_mid_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - Mid PV']
#     )
#     scenario_projects["No_Storage_Mandate_Mid_Additional_BTM_PV"] = \
#         no_storage_mandate_mid_customer_pv
# 
#     # 4: no storage mandate, high additional customer PV
#     no_storage_mandate_high_customer_pv = \
#         proj_capacities.copy()
#     no_storage_mandate_high_customer_pv.update(
#         scenario_storage_mandate_capacities['None']
#     )
#     no_storage_mandate_high_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - High PV']
#     )
#     scenario_projects["No_Storage_Mandate_High_Additional_BTM_PV"] = \
#         no_storage_mandate_high_customer_pv
# 
#     # 5: 1325 MW by 2024 storage mandate, no additional customer PV
#     storage_mandate_1325_by_2024_no_addtl_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1325_by_2024_no_addtl_customer_pv.update(
#         scenario_storage_mandate_capacities['1325 MW by 2024']
#     )
#     storage_mandate_1325_by_2024_no_addtl_customer_pv.update(
#         scenario_customer_pv_capacities['No Addtl. BTM']
#     )
#     scenario_projects["Storage_Mandate_1325MW_by_2024_No_Additional_BTM_PV"] \
#         = \
#         storage_mandate_1325_by_2024_no_addtl_customer_pv
# 
#     # 6: 1325 MW by 2024 storage mandate, low additional customer PV
#     storage_mandate_1325_by_2024_low_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1325_by_2024_low_customer_pv.update(
#         scenario_storage_mandate_capacities['1325 MW by 2024']
#     )
#     storage_mandate_1325_by_2024_low_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - Low PV']
#     )
#     scenario_projects["Storage_Mandate_1325MW_by_2024_Low_Additional_BTM_PV"] \
#         = \
#         storage_mandate_1325_by_2024_low_customer_pv
# 
#     # 7: 1325 MW by 2024 storage mandate, mid additional customer PV
#     storage_mandate_1325_by_2024_mid_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1325_by_2024_mid_customer_pv.update(
#         scenario_storage_mandate_capacities['1325 MW by 2024']
#     )
#     storage_mandate_1325_by_2024_mid_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - Mid PV']
#     )
#     scenario_projects["Storage_Mandate_1325MW_by_2024_Mid_Additional_BTM_PV"] \
#         = \
#         storage_mandate_1325_by_2024_mid_customer_pv
# 
#     # 8: 1325 MW by 2024 storage mandate, high additional customer PV
#     storage_mandate_1325_by_2024_high_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1325_by_2024_high_customer_pv.update(
#         scenario_storage_mandate_capacities['1325 MW by 2024']
#     )
#     storage_mandate_1325_by_2024_high_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - High PV']
#     )
#     scenario_projects["Storage_Mandate_1325MW_by_2024_High_Additional_BTM_PV"
#         ] = \
#         storage_mandate_1325_by_2024_high_customer_pv
# 
#     # 9: 1825 MW by 2024 storage mandate, no additional customer PV
#     storage_mandate_1825_by_2024_no_addtl_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1825_by_2024_no_addtl_customer_pv.update(
#         scenario_storage_mandate_capacities['1825 MW by 2024']
#     )
#     storage_mandate_1825_by_2024_no_addtl_customer_pv.update(
#         scenario_customer_pv_capacities['No Addtl. BTM']
#     )
#     scenario_projects["Storage_Mandate_1825MW_by_2024_No_Additional_BTM_PV"
#         ] = \
#         storage_mandate_1825_by_2024_no_addtl_customer_pv
# 
#     # 10: 1825 MW by 2024 storage mandate, low additional customer PV
#     storage_mandate_1825_by_2024_low_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1825_by_2024_low_customer_pv.update(
#         scenario_storage_mandate_capacities['1825 MW by 2024']
#     )
#     storage_mandate_1825_by_2024_low_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - Low PV']
#     )
#     scenario_projects["Storage_Mandate_1825MW_by_2024_Low_Additional_BTM_PV"
#      ] = \
#         storage_mandate_1825_by_2024_low_customer_pv
# 
#     # 11: 1825 MW by 2024 storage mandate, mid additional customer PV
#     storage_mandate_1825_by_2024_mid_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1825_by_2024_mid_customer_pv.update(
#         scenario_storage_mandate_capacities['1825 MW by 2024']
#     )
#     storage_mandate_1825_by_2024_mid_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - Mid PV']
#     )
#     scenario_projects["Storage_Mandate_1825MW_by_2024_Mid_Additional_BTM_PV"
#         ] = \
#         storage_mandate_1825_by_2024_mid_customer_pv
# 
#     # 12: 1825 MW by 2024 storage mandate, high additional customer PV
#     storage_mandate_1825_by_2024_high_customer_pv = \
#         proj_capacities.copy()
#     storage_mandate_1825_by_2024_high_customer_pv.update(
#         scenario_storage_mandate_capacities['1825 MW by 2024']
#     )
#     storage_mandate_1825_by_2024_high_customer_pv.update(
#         scenario_customer_pv_capacities['CEC 2016 IEPR - High PV']
#     )
#     scenario_projects["Storage_Mandate_1825MW_by_2024_High_Additional_BTM_PV"
#         ] = \
#         storage_mandate_1825_by_2024_high_customer_pv
# 
#     # Insert data
#     scenario_id = 1
#     for scenario in scenario_projects.keys():
#         project_existing_params.update_project_capacities(
#             io=io, c=c2,
#             project_existing_capacity_scenario_id=scenario_id,
#             scenario_name=scenario,
#             scenario_description=scenario,
#             project_capacities=scenario_projects[scenario]
#         )
#         scenario_id += 1
# 
#     # Hybridized gas counts toward storage mandate
#     c2.execute(
#         """INSERT INTO subscenarios_project_existing_capacity
#         (project_existing_capacity_scenario_id, name, description)
#         VALUES (13, 'Storage_Mandate_1103MW_by_2024_Mid_Additional_BTM_PV', 
#         '222 MW of hybridized gas counts toward storage mandate');"""
#     )
# 
#     c2.execute(
#         """
#         INSERT INTO inputs_project_existing_capacity
#         (project_existing_capacity_scenario_id, project, period, 
#         existing_capacity_mw, existing_capacity_mwh)
#         SELECT 13, project, period, existing_capacity_mw, existing_capacity_mwh
#         FROM inputs_project_existing_capacity
#         WHERE project_existing_capacity_scenario_id = 7
#         AND project != 'CAISO_Storage_Mandate';
#         """
#     )
# 
#     c2.execute(
#         """
#         INSERT INTO inputs_project_existing_capacity
#         (project_existing_capacity_scenario_id, project, period, 
#         existing_capacity_mw, existing_capacity_mwh)
#         SELECT 13, project, period, existing_capacity_mw - 222, 
#         (existing_capacity_mw - 222) * 4
#         FROM inputs_project_existing_capacity
#         WHERE project_existing_capacity_scenario_id = 7
#         AND project = 'CAISO_Storage_Mandate';
#         """
#     )
# 
#     io.commit()
# 
#     # CCGT analysis retirements scenario fleet capacities
#     capacity_to_subtract_by_fleet = {
#         'CCGT_SuperInefficient': 751.09,
#         'Peaker_Average2': 991.79,
#         'Peaker_SuperInefficient': 1140.88,
#         'CCGT_Inefficient': 1909.63,
#         'CCGT_AverageFlex': 263.00,
#         'Peaker_Average1': 65.00
#     }
#     scenario_projects['ccgt_ret_analysis'] = \
#         scenario_projects["Storage_Mandate_1825MW_by_2024_Mid_Additional_BTM_PV"]
# 
#     for fleet in capacity_to_subtract_by_fleet.keys():
#         for year in range(2015, 2051):
#             scenario_projects['ccgt_ret_analysis'][fleet][year] = \
#                 (scenario_projects['ccgt_ret_analysis'][fleet][year][0] -
#                  capacity_to_subtract_by_fleet[fleet],
#                  scenario_projects['ccgt_ret_analysis'][fleet][year][1])
# 
#     project_existing_params.update_project_capacities(
#         io=io, c=c2,
#         project_existing_capacity_scenario_id=14,
#         scenario_name=
#         'CCGT_Analysis_Ret_Storage_Mandate_1825MW_by_2024_Mid_Additional_BTM_PV',
#         scenario_description=
#         'CCGT_Analysis_Ret_Storage_Mandate_1825MW_by_2024_Mid_Additional_BTM_PV',
#         project_capacities=scenario_projects['ccgt_ret_analysis']
#     )
def load_project_capacities():
    """

    :return:
    """
    # Load the data into resolve.db, as it will be easier to split the
    # existing and planned projects
    c1.execute(
        """DROP TABLE IF EXISTS CONV_planned_capacity;"""
    )
    c1.execute(
        """CREATE TABLE CONV_planned_capacity(
            project VARCHAR(32),
            proj_zone VARCHAR(16),
            forecast_year INTEGER,
            capacity_mw FLOAT,
            PRIMARY KEY (project, forecast_year)
            );"""
    )
    resolve.commit()

    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_Baseline_no_early_retirements.csv"),
              "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(55 - 1, 81):
            if all(x == "" for x in rows_list[row]):
                pass  # skip if no data in row
            else:
                for column in range(6 - 1, 41):
                    c1.execute(
                        """INSERT INTO CONV_planned_capacity 
                        (project, proj_zone, forecast_year, capacity_mw)
                        VALUES ('{}', '{}', {}, {});""".format(
                            rows_list[row][3 - 1],
                            rows_list[row][4 - 1],
                            int(float(rows_list[6 - 1][column])),
                            rows_list[row][column]
                        )
                    )
        resolve.commit()

    # Split any projects that have capacity added in the future
    # Make key table
    c1.execute("DROP TABLE IF EXISTS PROCESSED_split_project_list;")
    c1.execute(
        """CREATE TABLE PROCESSED_split_project_list (
        split_project VARCHAR(64),
        parent_project VARCHAR(64),
        PRIMARY KEY (split_project, parent_project)
        );"""
    )
    resolve.commit()

    # Table with new capacities
    c1.execute("DROP TABLE IF EXISTS PROCESSED_CONV_planned_capacity;")
    c1.execute(
        """CREATE TABLE PROCESSED_CONV_planned_capacity (
        project VARCHAR(32),
        proj_zone VARCHAR(16),
        forecast_year INTEGER,
        capacity_mw FLOAT,
        PRIMARY KEY (project, forecast_year)
        );"""
    )
    resolve.commit()

    # Insert non-CAISO project capacities
    c1.execute(
        """INSERT INTO PROCESSED_CONV_planned_capacity
        (project, proj_zone, forecast_year, capacity_mw)
        SELECT project, proj_zone, forecast_year, capacity_mw
        FROM CONV_planned_capacity
        WHERE project NOT LIKE 'CAISO%';"""
    )
    resolve.commit()

    # CAISO projects
    # Process data -- split planned from existing apacity
    parent_projects = c1.execute(
        """SELECT DISTINCT project
        FROM CONV_planned_capacity
        WHERE project LIKE 'CAISO%';"""
    ).fetchall()

    # Start with year 2015; insert data from CONV_planned_capacity table
    c1.execute(
        """INSERT INTO PROCESSED_CONV_planned_capacity
        (project, proj_zone, forecast_year, capacity_mw)
        SELECT project, proj_zone, forecast_year, capacity_mw
        FROM CONV_planned_capacity
        WHERE forecast_year = 2015
        AND project LIKE 'CAISO%';"""
    )
    resolve.commit()

    for p in parent_projects:
        for yr in range(2016, 2051):
            # Get the project zone
            zone = c1.execute(
                """SELECT proj_zone 
                FROM CONV_planned_capacity 
                WHERE project = '{}'
                AND forecast_year = {}""".format(p[0], yr)
            ).fetchone()

            split_project_name = p[0] + '_Planned'
            current_capacity = c1.execute(
                """SELECT capacity_mw
                FROM CONV_planned_capacity
                WHERE project = '{}'
                AND forecast_year = {}""".format(p[0], yr)
            ).fetchone()
            previous_capacity = c1.execute(
                """SELECT capacity_mw
                FROM CONV_planned_capacity
                WHERE project = '{}'
                AND forecast_year = {}""".format(p[0], yr - 1)
            ).fetchone()

            # Increase or decrease in capacity from previous year
            delta = current_capacity[0] - previous_capacity[0]

            # If we see increase in capacity, we'll need to split
            # the project
            if delta > 0:
                # Insert split project into our key table if not already
                # there
                c1.execute(
                    """INSERT OR IGNORE INTO PROCESSED_split_project_list
                    (split_project, parent_project)
                    VALUES ('{}', '{}')""".format(split_project_name, p[0])
                )
                # Insert capacities for this year into capacities table
                # First, get last year's capacities
                parent_project_previous_capacity = c1.execute(
                    """SELECT capacity_mw
                    FROM PROCESSED_CONV_planned_capacity
                    WHERE project = '{}'
                    AND forecast_year = {};""".format(p[0], yr - 1)
                ).fetchone()

                split_project_previous_capacity = c1.execute(
                    """SELECT capacity_mw 
                    FROM PROCESSED_CONV_planned_capacity
                    WHERE project = '{}'
                    AND forecast_year = {};""".format(
                        split_project_name, yr - 1
                    )
                ).fetchone()

                # If there was no capacity for the split project in the
                # previous year, it was 0
                if split_project_previous_capacity is None:
                    split_project_previous_capacity = 0
                    c1.execute(
                        """INSERT INTO PROCESSED_CONV_planned_capacity
                        (project, proj_zone, forecast_year, capacity_mw)
                        VALUES ('{}', '{}', {}, {});""".format(
                            split_project_name, zone[0], yr - 1,
                            split_project_previous_capacity
                        )
                    )
                # Insert previous capacity + the delta for the split prj
                c1.execute(
                    """INSERT INTO PROCESSED_CONV_planned_capacity
                    (project, proj_zone, forecast_year, capacity_mw)
                    VALUES ('{}', '{}', {}, {});""".format(
                        split_project_name, zone[0], yr,
                        (0 if split_project_previous_capacity == 0
                         else split_project_previous_capacity[0]) + delta
                    )
                )
                # Insert the previous capacity for the parent project
                c1.execute(
                    """INSERT INTO PROCESSED_CONV_planned_capacity
                    (project, proj_zone, forecast_year, capacity_mw)
                    VALUES ('{}', '{}', {}, {});""".format(
                        p[0], zone[0], yr,
                        parent_project_previous_capacity[0]
                    )
                )
            # If we have the same as the previous year,
            # simply insert the calculated previous year capacity
            # + delta for both the parent and the previous year capacity
            # for the split project
            elif delta <= 0:
                # Insert the previous capacity for the parent project
                c1.execute(
                    """INSERT INTO PROCESSED_CONV_planned_capacity
                    (project, proj_zone, forecast_year, capacity_mw)
                    SELECT project, proj_zone, forecast_year + 1, 
                    capacity_mw + {}
                    FROM PROCESSED_CONV_planned_capacity
                    WHERE project = '{}'
                    AND forecast_year = {};""".format(
                        delta, p[0], yr - 1
                    )
                )

                # If no previous year capacity, we won't insert anything
                c1.execute(
                    """INSERT INTO PROCESSED_CONV_planned_capacity
                    (project, proj_zone, forecast_year, capacity_mw)
                    SELECT project, proj_zone, forecast_year + 1, 
                    capacity_mw
                    FROM PROCESSED_CONV_planned_capacity
                    WHERE project = '{}'
                    AND forecast_year = {};""".format(
                        split_project_name, yr - 1
                    )
                )

    resolve.commit()

    # Make final dictionary
    proj_capacities = OrderedDict()

    # Add conventional project capacities
    distinct_conv_projects = c1.execute(
        """SELECT DISTINCT project
        FROM PROCESSED_CONV_planned_capacity
        ORDER BY project;"""
    ).fetchall()

    conv_projects = c1.execute(
        """SELECT project, forecast_year, capacity_mw
        FROM PROCESSED_CONV_planned_capacity
        ORDER BY project, forecast_year;"""
    ).fetchall()

    for project in distinct_conv_projects:
        p = project[0]
        proj_capacities[p] = OrderedDict()

    for proj_yr_cap in conv_projects:
        project = proj_yr_cap[0]
        year = proj_yr_cap[1]
        capacity = proj_yr_cap[2]
        proj_capacities[project][year] = (capacity, None)

    # Add CAISO_Shed_DR_Existing capacities
    with open(os.path.join("cpuc_irp_data", "csvs", "DR_Conventional.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        proj_capacities['CAISO_Shed_DR_Existing'] = OrderedDict()
        for clmn in range(10 - 1, 45):
            dr_yr = int(float(rows_list[31 - 1][clmn]))
            proj_capacities['CAISO_Shed_DR_Existing'][dr_yr] = \
                (float(rows_list[35 - 1][clmn]), None)

    # Add hydro project capacities
    with open(os.path.join("cpuc_irp_data", "csvs", "HYD_OpChar.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for column in range(5 - 1, 10):
            proj_capacities[rows_list[4 - 1][column]] = OrderedDict()
            for yr in range(2015, 2050 + 1):
                proj_capacities[rows_list[4 - 1][column]][yr] = \
                    (float(rows_list[5 - 1][column]), None)

    # Add existing/planned renewable capacities
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Baseline.csv"), "r") as f:
        # This will make a list of lists where each (sub)list is a row and
        # each element of the (sub)list is a column
        rows_list = list(csv.reader(f))

        # CAISO resources
        for row in range(59 - 1, 73):
            project = rows_list[row][3 - 1]
            proj_capacities[project] = OrderedDict()
            for column in range(10 - 1, 45):
                year = int(float(rows_list[3 - 1][column]))
                proj_capacities[project][year] \
                    = (float(rows_list[row][column]), None)

        # Non-CAISO resources
        for row in range(339 - 1, 368):
            project = rows_list[row][3 - 1]
            proj_capacities[project] = OrderedDict()
            for column in range(10 - 1, 45):
                year = int(float(rows_list[3 - 1][column]))
                proj_capacities[project][year] \
                    = (float(rows_list[row][column]), None)

    # Disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific", "gas_disagg_plants.csv"),
              "r") as foo:
        reader = csv.reader(foo, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery') \
                    or row[0].startswith('Advanced_CCGT') \
                    or row[0].startswith('Aero_CT'):
                pass
            else:
                proj_capacities[row[0]] = OrderedDict()
                if row[20] == '':
                    start_year = 2015
                else:
                    start_year = int(float(row[20]))
                if row[21] == '':
                    end_year = 2050 + 1
                else:
                    end_year = int(float(row[21])) + 1
                for yr in range(2015, 2050 + 1):
                    if start_year <= yr <= end_year:
                        proj_capacities[row[0]][yr] = (float(row[16]), None)
                    else:
                        proj_capacities[row[0]][yr] = (0, None)

    # Add storage capacities
    with open(os.path.join("cpuc_irp_data", "csvs", "STOR_Inputs.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # Existing pumped storage
        ps_project = rows_list[28 - 1][3 - 1]
        proj_capacities[ps_project] = OrderedDict()
        for column in range(6 - 1, 41):
            year = int(float(rows_list[27 - 1][column]))
            capacity = rows_list[28 - 1][column]
            energy = rows_list[36 - 1][column]
            proj_capacities[ps_project][year] = (capacity, energy)

        # 3 storage mandate scenarios
        # Storage mandate assumed to have 4 hours duration
        scenario_storage_mandate_capacities = OrderedDict()
        sm_project = rows_list[29 - 1][3 - 1]
        for scenario_row in range(21 - 1, 23):
            scenario = rows_list[scenario_row][3 - 1]
            scenario_storage_mandate_capacities[scenario] = OrderedDict()
            scenario_storage_mandate_capacities[scenario][sm_project] = \
                OrderedDict()
            for column in range(6 - 1, 41):
                sm_year = int(float(rows_list[20 - 1][column]))
                sm_capacity = float(rows_list[scenario_row][column])
                sm_energy = 4 * sm_capacity
                scenario_storage_mandate_capacities[
                    scenario
                ][
                    sm_project
                ][
                    sm_year] = (sm_capacity, sm_energy)

    # Add Customer_PV capacities
    # 4 scenarios
    # We'll need the shape capacity factor to convert to MW from GWh
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Profiles.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        customer_pv_cap_factor = float(rows_list[23 - 1][42 - 1])

    with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_Forecast.csv"), "r") as f:
        scenario_customer_pv_capacities = OrderedDict()
        rows_list = list(csv.reader(f))
        for scenario_row in range(29 - 1, 32):
            custpv_scenario = rows_list[scenario_row][3 - 1]
            scenario_customer_pv_capacities[custpv_scenario] = \
                OrderedDict()
            scenario_customer_pv_capacities[custpv_scenario]['Customer_PV'
            ] = OrderedDict()
            for column in range(10 - 1, 45):
                custpv_year = int(float(rows_list[28 - 1][column]))
                custpv_capacity = \
                    float(rows_list[scenario_row][column]) * 1000 / 8760 / \
                    customer_pv_cap_factor
                scenario_customer_pv_capacities[
                    custpv_scenario
                ][
                    "Customer_PV"
                ][
                    custpv_year] = \
                    (custpv_capacity, None)

                # Make final 12 scenarios (base X 3 storage mandate cases X 4 customer
    # PV cases)
    scenario_projects = OrderedDict()

    # 1: no storage mandate, no additional customer PV
    no_storage_mandate_no_addtl_customer_pv = \
        proj_capacities.copy()
    no_storage_mandate_no_addtl_customer_pv.update(
        scenario_storage_mandate_capacities['None']
    )
    no_storage_mandate_no_addtl_customer_pv.update(
        scenario_customer_pv_capacities['No Addtl. BTM']
    )
    scenario_projects["No_Storage_Mandate_No_Additional_BTM_PV"] = \
        no_storage_mandate_no_addtl_customer_pv

    # 2: no storage mandate, low additional customer PV
    no_storage_mandate_low_customer_pv = \
        proj_capacities.copy()
    no_storage_mandate_low_customer_pv.update(
        scenario_storage_mandate_capacities['None']
    )
    no_storage_mandate_low_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - Low PV']
    )
    scenario_projects["No_Storage_Mandate_Low_Additional_BTM_PV"] = \
        no_storage_mandate_low_customer_pv

    # 3: no storage mandate, mid additional customer PV
    no_storage_mandate_mid_customer_pv = \
        proj_capacities.copy()
    no_storage_mandate_mid_customer_pv.update(
        scenario_storage_mandate_capacities['None']
    )
    no_storage_mandate_mid_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - Mid PV']
    )
    scenario_projects["No_Storage_Mandate_Mid_Additional_BTM_PV"] = \
        no_storage_mandate_mid_customer_pv

    # 4: no storage mandate, high additional customer PV
    no_storage_mandate_high_customer_pv = \
        proj_capacities.copy()
    no_storage_mandate_high_customer_pv.update(
        scenario_storage_mandate_capacities['None']
    )
    no_storage_mandate_high_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - High PV']
    )
    scenario_projects["No_Storage_Mandate_High_Additional_BTM_PV"] = \
        no_storage_mandate_high_customer_pv

    # 5: 1325 MW by 2024 storage mandate, no additional customer PV
    storage_mandate_1325_by_2024_no_addtl_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1325_by_2024_no_addtl_customer_pv.update(
        scenario_storage_mandate_capacities['1325 MW by 2024']
    )
    storage_mandate_1325_by_2024_no_addtl_customer_pv.update(
        scenario_customer_pv_capacities['No Addtl. BTM']
    )
    scenario_projects["Storage_Mandate_1325MW_by_2024_No_Additional_BTM_PV"] \
        = \
        storage_mandate_1325_by_2024_no_addtl_customer_pv

    # 6: 1325 MW by 2024 storage mandate, low additional customer PV
    storage_mandate_1325_by_2024_low_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1325_by_2024_low_customer_pv.update(
        scenario_storage_mandate_capacities['1325 MW by 2024']
    )
    storage_mandate_1325_by_2024_low_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - Low PV']
    )
    scenario_projects["Storage_Mandate_1325MW_by_2024_Low_Additional_BTM_PV"] \
        = \
        storage_mandate_1325_by_2024_low_customer_pv

    # 7: 1325 MW by 2024 storage mandate, mid additional customer PV
    storage_mandate_1325_by_2024_mid_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1325_by_2024_mid_customer_pv.update(
        scenario_storage_mandate_capacities['1325 MW by 2024']
    )
    storage_mandate_1325_by_2024_mid_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - Mid PV']
    )
    scenario_projects["Storage_Mandate_1325MW_by_2024_Mid_Additional_BTM_PV"] \
        = \
        storage_mandate_1325_by_2024_mid_customer_pv

    # 8: 1325 MW by 2024 storage mandate, high additional customer PV
    storage_mandate_1325_by_2024_high_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1325_by_2024_high_customer_pv.update(
        scenario_storage_mandate_capacities['1325 MW by 2024']
    )
    storage_mandate_1325_by_2024_high_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - High PV']
    )
    scenario_projects["Storage_Mandate_1325MW_by_2024_High_Additional_BTM_PV"
    ] = \
        storage_mandate_1325_by_2024_high_customer_pv

    # 9: 1825 MW by 2024 storage mandate, no additional customer PV
    storage_mandate_1825_by_2024_no_addtl_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1825_by_2024_no_addtl_customer_pv.update(
        scenario_storage_mandate_capacities['1825 MW by 2024']
    )
    storage_mandate_1825_by_2024_no_addtl_customer_pv.update(
        scenario_customer_pv_capacities['No Addtl. BTM']
    )
    scenario_projects["Storage_Mandate_1825MW_by_2024_No_Additional_BTM_PV"
    ] = \
        storage_mandate_1825_by_2024_no_addtl_customer_pv

    # 10: 1825 MW by 2024 storage mandate, low additional customer PV
    storage_mandate_1825_by_2024_low_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1825_by_2024_low_customer_pv.update(
        scenario_storage_mandate_capacities['1825 MW by 2024']
    )
    storage_mandate_1825_by_2024_low_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - Low PV']
    )
    scenario_projects["Storage_Mandate_1825MW_by_2024_Low_Additional_BTM_PV"
    ] = \
        storage_mandate_1825_by_2024_low_customer_pv

    # 11: 1825 MW by 2024 storage mandate, mid additional customer PV
    storage_mandate_1825_by_2024_mid_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1825_by_2024_mid_customer_pv.update(
        scenario_storage_mandate_capacities['1825 MW by 2024']
    )
    storage_mandate_1825_by_2024_mid_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - Mid PV']
    )
    scenario_projects["Storage_Mandate_1825MW_by_2024_Mid_Additional_BTM_PV"
    ] = \
        storage_mandate_1825_by_2024_mid_customer_pv

    # 12: 1825 MW by 2024 storage mandate, high additional customer PV
    storage_mandate_1825_by_2024_high_customer_pv = \
        proj_capacities.copy()
    storage_mandate_1825_by_2024_high_customer_pv.update(
        scenario_storage_mandate_capacities['1825 MW by 2024']
    )
    storage_mandate_1825_by_2024_high_customer_pv.update(
        scenario_customer_pv_capacities['CEC 2016 IEPR - High PV']
    )
    scenario_projects["Storage_Mandate_1825MW_by_2024_High_Additional_BTM_PV"
    ] = \
        storage_mandate_1825_by_2024_high_customer_pv

    # Insert data
    scenario_id = 1
    for scenario in scenario_projects.keys():
        project_existing_params.update_project_capacities(
            io=io, c=c2,
            project_existing_capacity_scenario_id=scenario_id,
            scenario_name=scenario,
            scenario_description=scenario,
            project_capacities=scenario_projects[scenario]
        )
        scenario_id += 1

def load_project_fixed_costs():
    """

    :return:
    """
    fixed_costs_by_tech = dict()

    with open(os.path.join("cpuc_irp_data", "csvs", "COSTS_Resource_Char.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for column in range(9 - 1, 38):
            tech = rows_list[13 - 1][column]
            fixed_cost = float(rows_list[28 - 1][column])
            fixed_costs_by_tech[tech] = fixed_cost

    fc_tech_by_project = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "project_tech_for_fixed_costs.csv"),
              "r") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            fc_tech_by_project[row[0]] = (row[1], row[2])

    # We'll apply fixed costs only to projects that will be allowed to
    # retire; the rest will get 0s
    all_projects = [p[0] for p in c2.execute(
        """SELECT DISTINCT project
        FROM inputs_project_existing_capacity;"""
    ).fetchall()]
    relevant_projects = [p[0] for p in c2.execute(
        """SELECT DISTINCT project
        FROM inputs_project_existing_capacity
        WHERE project IN (
        'CAISO_CCGT1', 
        'CAISO_CCGT2',
        'CAISO_CHP',
        'CAISO_Peaker1',
        'CAISO_Peaker2',
        'CAISO_Reciprocating_Engine',
        'CAISO_ST'
        );"""
    ).fetchall()]

    project_fixed_costs = dict()
    # Find project and assign fixed cost based on tech
    for project in all_projects:
        project_fixed_costs[project] = OrderedDict()
        for period in range(2015, 2050+1):  # no escalation
            if project in relevant_projects:
                if fc_tech_by_project[project][0] == '':
                    pass  # This will skip CAISO_CHP and CAISO_ST
                else:
                    project_fixed_costs[project][period] = \
                        (fixed_costs_by_tech[fc_tech_by_project[project][0]],
                         None)
            else:
                # If storage, fixed costs per kW-yr and kWh-yr
                if project in [
                    "CAISO_Storage_Mandate",
                    "CAISO_Existing_Pumped_Storage"
                ]:
                    project_fixed_costs[project][period] = \
                        (0, 0)
                else:
                    project_fixed_costs[project][period] = (0, None)

    # Separately update CHP and STs, for which there is no data in the
    # scenario tool
    with open(os.path.join("cpuc_irp_data", "csvs", "chp_and_st_fixed_costs.csv"),
              "r") as f:
        reader = csv.reader(f)
        next(reader)
        for row in reader:
            project_fixed_costs[row[0]] = OrderedDict()
            for period in range(2015, 2050+1):
                project_fixed_costs[row[0]][period] = (float(row[1]), None)

    # Disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific", "gas_disagg_plants.csv"),
              "r") as foo:
        reader = csv.reader(foo, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith('Battery') \
                    or row[0].startswith('Advanced_CCGT') \
                    or row[0].startswith('Aero_CT'):
                pass
            else:
                project_fixed_costs[row[0]] = OrderedDict()
                for yr in range(2015, 2050+1):
                    project_fixed_costs[row[0]][yr] = (float(row[15]), None)

    # Insert data
    project_existing_params.update_project_fixed_costs(
        io=io, c=c2,
        project_existing_fixed_cost_scenario_id=1,
        scenario_name="default_fixed_costs_for_caiso_gas",
        scenario_description='Default fixed costs applied to CAISO existing '
                             'gas generators',
        project_fixed_costs=project_fixed_costs
)


def load_project_new_costs():
    """
    Costs for candidate projects
    :return:
    """

    # GridPath needs project lifetimes (financing lifetime assumed the same
    # as project lifetime)
    # These are financing lifetimes from COSTS_Resource_Char (manual)
    # 'Biomass': 20
    # 'Geothermal': 20
    # 'Solar': 25
    # 'Wind': 25
    # 'Battery': 20
    # 'Pumped Storage': 25
    # 'CCGT': 20
    # 'CT': 20
    # 'Reciprocating_Engine': 20

    project_period_lifetimes_costs = OrderedDict()

    # Conventional candidate projects
    with open(os.path.join("cpuc_irp_data", "csvs", "CONV_Candidate.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(11 - 1, 13):
            if rows_list[row][3 - 1] == 'CAISO_Reciprocating_Engine':
                conv_proj = 'CAISO_Reciprocating_Engine_Candidate'
            else:
                conv_proj = rows_list[row][3 - 1]
            project_period_lifetimes_costs[conv_proj] = OrderedDict()
            for column in range(6-1, 41):
                conv_yr = int(float(rows_list[10 - 1][column]))
                project_period_lifetimes_costs[conv_proj][conv_yr] = \
                    (20, float(rows_list[row][column]), None, None)

    # Renewable candidate projects
    # Take the 'Total Fixed Technology Costs for Candidate Resources ($/kW-yr)
    # - sum of Tx and Resource Cost' table
    # These are the 'mid' costs for solar
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(240 - 1, 281):
            ren_proj = rows_list[row][3 - 1]
            project_period_lifetimes_costs[ren_proj] = OrderedDict()
            for column in range(10-1, 45):
                ren_yr = int(float(rows_list[239 - 1][column]))
                # Lifetime is 25 years if solar or wind, 20 years if biomass
                #  or geothermal
                if 'Wind' or 'Solar' in ren_proj:
                    ren_lifetime = 25
                elif 'Biomass' or 'Geothermal' in ren_proj:
                    ren_lifetime = 20
                else:  # this shouldn't happen
                    ren_lifetime = None

                project_period_lifetimes_costs[ren_proj][ren_yr] = \
                    (ren_lifetime, float(rows_list[row][column]), None, None)
    # Shed DR
    with open(os.path.join("cpuc_irp_data", "csvs", "DR_Conventional.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(40 - 1, 47):
            dr_proj = rows_list[row][3 - 1]
            project_period_lifetimes_costs[dr_proj] = OrderedDict()
            for column in range(10 - 1, 45):
                dr_yr = int(float(rows_list[39 - 1][column]))
                # Give DR a lifetime of 20 years, so that it's present the
                # entire study
                project_period_lifetimes_costs[dr_proj][dr_yr] = \
                    (20, float(rows_list[row][column]), None, None)

    # Storage candidate projects
    # These are the 'mid' costs for storage
    # Low and high available in Table 29, p. 42 of the IRP Staff Proposal
    with open(os.path.join("cpuc_irp_data", "csvs", "STOR_Inputs.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(74 - 1, 76):
            stor_proj = rows_list[row][3 - 1]
            project_period_lifetimes_costs[stor_proj] = OrderedDict()
            for column in range(6 - 1, 41):
                stor_yr = int(float(rows_list[73 - 1][column]))
                # Lifetime is 25 yrs for pumped storage, 20 yrs for batteries
                if 'Pumped' in stor_proj:
                    stor_lifetime = 25
                elif 'Battery' in stor_proj:
                    stor_lifetime = 20
                else:
                    stor_lifetime = None  # this shouldn't happen

                project_period_lifetimes_costs[stor_proj][stor_yr] = \
                    (stor_lifetime, float(rows_list[row][column]),)

        for row in range(80 - 1, 82):
            for column in range(6 - 1, 41):
                project_period_lifetimes_costs[
                    rows_list[row][3 - 1]
                ][
                    int(float(rows_list[79 - 1][column]))
                ] = \
                    project_period_lifetimes_costs[
                        rows_list[row][3 - 1]
                    ][
                        int(float(rows_list[79 - 1][column]))
                    ] + (float(rows_list[row][column]), None)

    # Shift DR

    # Manually update storage costs to the newer version of the data
    # from Sep 2017
    # Newer storage costs
    project_period_lifetimes_costs['CAISO_New_Flow_Battery'][2018] = \
        (20, 273.769086905514, 30.8468596528025, None)
    project_period_lifetimes_costs['CAISO_New_Flow_Battery'][2022] = \
        (20, 222.155995290217, 25.0313681694496, None)
    project_period_lifetimes_costs['CAISO_New_Flow_Battery'][2026] = \
        (20, 196.454628548388, 22.1354734512705, None)
    project_period_lifetimes_costs['CAISO_New_Flow_Battery'][2030] = \
        (20, 190.036799973108, 21.4123463094331, None)
    project_period_lifetimes_costs['CAISO_New_Li_Battery'][2018] = \
        (20, 50.4869109291559, 69.1935130107353, None)
    project_period_lifetimes_costs['CAISO_New_Li_Battery'][2022] = \
        (20, 35.6553801527522, 48.7557195080645, None)
    project_period_lifetimes_costs['CAISO_New_Li_Battery'][2026] = \
        (20, 29.0961953560795, 39.7865885444321, None)
    project_period_lifetimes_costs['CAISO_New_Li_Battery'][2030] = \
        (20, 27.550798250155, 37.6733885868876, None)

    # LCR-eligible batteries
    for lcr_battery in [
        'Battery_LCR_Bay_Area',
        'Battery_LCR_Big_Creek_Ventura',
        'Battery_LCR_Fresno',
        'Battery_LCR_Humboldt',
        'Battery_LCR_Kern',
        'Battery_LCR_LA_Basin',
        'Battery_LCR_NCNB',
        'Battery_LCR_San_Diego',
        'Battery_LCR_Sierra',
        'Battery_LCR_Stockton'
    ]:
        project_period_lifetimes_costs[lcr_battery] = OrderedDict()
        project_period_lifetimes_costs[lcr_battery][2018] = \
            (20, 50.4869109291559, 69.1935130107353, None)
        project_period_lifetimes_costs[lcr_battery][2022] = \
            (20, 35.6553801527522, 48.7557195080645, None)
        project_period_lifetimes_costs[lcr_battery][2026] = \
            (20, 29.0961953560795, 39.7865885444321, None)
        project_period_lifetimes_costs[lcr_battery][2030] = \
            (20, 27.550798250155, 37.6733885868876, None)

    # LCR-eligible CCGTs
    for lcr_ccgt in [
        'Advanced_CCGT_LCR_Bay_Area',
        'Advanced_CCGT_LCR_Big_Creek_Ventura',
        'Advanced_CCGT_LCR_Fresno',
        'Advanced_CCGT_LCR_Humboldt',
        'Advanced_CCGT_LCR_Kern',
        'Advanced_CCGT_LCR_LA_Basin',
        'Advanced_CCGT_LCR_NCNB',
        'Advanced_CCGT_LCR_San_Diego',
        'Advanced_CCGT_LCR_Sierra',
        'Advanced_CCGT_LCR_Stockton'
    ]:
        project_period_lifetimes_costs[lcr_ccgt] = OrderedDict()
        for yr in [2018, 2022, 2026, 2030]:
            project_period_lifetimes_costs[lcr_ccgt][yr] = \
                (20, 201.848094258, None, None)
    
    # LCR-eligible CTs
    for lcr_ct in [
        'Aero_CT_LCR_Bay_Area',
        'Aero_CT_LCR_Big_Creek_Ventura',
        'Aero_CT_LCR_Fresno',
        'Aero_CT_LCR_Humboldt',
        'Aero_CT_LCR_Kern',
        'Aero_CT_LCR_LA_Basin',
        'Aero_CT_LCR_NCNB',
        'Aero_CT_LCR_San_Diego',
        'Aero_CT_LCR_Sierra',
        'Aero_CT_LCR_Stockton'
    ]:
        project_period_lifetimes_costs[lcr_ct] = OrderedDict()
        for yr in [2018, 2022, 2026, 2030]:
            project_period_lifetimes_costs[lcr_ct][yr] = \
                (20, 196.97067234, None, None)

    project_new_costs.update_project_new_costs(
        io=io, c=c2,
        project_new_cost_scenario_id=1,
        scenario_name='default_costs_mid_solar_mid_storage',
        scenario_description='Default costs (mid solar PV costs and mid '
                             'storage costs)',
        project_period_lifetimes_costs=project_period_lifetimes_costs
    )

    # Make subscenarios with low and high storage costs
    # Got the costs manually from the RESOLVE cases
    # li_ion_low = {
    #     2018: (23.9260355369971, 59.3201444241859),
    #     2022: (19.6749510093235, 49.904248483839),
    #     2026: (17.7709234922922, 44.2547109196308),
    #     2030: (17.1362476532818, 43.3131213255961)
    # }
    #
    # li_ion_high = {
    #     2018: (34.4197542789376, 122.240311404952),
    #     2022: (26.5260420128839, 95.1959947224402),
    #     2026: (22.8325424921026, 82.2147227148347),
    #     2030: (21.8252244409804, 78.9694047129334)
    # }
    # These are the costs from Sep 2017 cases
    li_ion_low = {
        2018: (35.9306356878573, 33.9643604761634),
        2022: (23.386462997225, 22.0876370062929),
        2026: (18.2138947010904, 17.2023402887501),
        2030: (17.0358748885025, 16.0897447667274)
    }
    flow_low = {
        2018: (206.738244254675, 22.6506044502134),
        2022: (158.160600766124, 17.3283526735794),
        2026: (135.156802945179, 14.8080162589339),
        2030: (129.552788585568, 14.1940306219248)
    }

    li_ion_high = {
        2018: (66.2280716152593, 120.81764409051),
        2022: (50.6663468678193, 92.1932374734895),
        2026: (43.2971386451174, 78.7841167125319),
        2030: (41.5019068742408, 75.5174863118014)
    }
    flow_high = {
        2018: (344.76696800075, 39.5322903749105),
        2022: (296.497138843087, 33.9974883790188),
        2026: (271.242515618114, 31.1016973337571),
        2030: (264.780329768902, 30.3607185534259)
    }

    low_storage_costs = dict()
    high_storage_costs = dict()
    for prj in project_period_lifetimes_costs.keys():
        low_storage_costs[prj] = OrderedDict()
        high_storage_costs[prj] = OrderedDict()
        for yr in [2018, 2022, 2026, 2030]:
            if prj in [
                'CAISO_New_Li_Battery',
                'Battery_LCR_Bay_Area',
                'Battery_LCR_Big_Creek_Ventura',
                'Battery_LCR_Fresno',
                'Battery_LCR_Humboldt',
                'Battery_LCR_Kern',
                'Battery_LCR_LA_Basin',
                'Battery_LCR_NCNB',
                'Battery_LCR_San_Diego',
                'Battery_LCR_Sierra',
                'Battery_LCR_Stockton'
            ]:
                low_storage_costs[prj][yr] = (
                    project_period_lifetimes_costs[prj][yr][0],
                    li_ion_low[yr][0], li_ion_low[yr][1],
                    None
                )
                high_storage_costs[prj][yr] = (
                    project_period_lifetimes_costs[prj][yr][0],
                    li_ion_high[yr][0], li_ion_high[yr][1],
                    None
                )
            elif prj == 'CAISO_New_Flow_Battery':
                low_storage_costs[prj][yr] = (
                    project_period_lifetimes_costs[prj][yr][0],
                    flow_low[yr][0], flow_low[yr][1],
                    None
                )
                high_storage_costs[prj][yr] = (
                    project_period_lifetimes_costs[prj][yr][0],
                    flow_high[yr][0], flow_high[yr][1],
                    None
                )
            else:
                low_storage_costs[prj][yr] \
                    = project_period_lifetimes_costs[prj][yr]
                high_storage_costs[prj][yr] \
                    = project_period_lifetimes_costs[prj][yr]

    project_new_costs.update_project_new_costs(
        io=io, c=c2,
        project_new_cost_scenario_id=2,
        scenario_name='low_li_io_costs',
        scenario_description='Like default but low Li-ion costs',
        project_period_lifetimes_costs=low_storage_costs
    )

    project_new_costs.update_project_new_costs(
        io=io, c=c2,
        project_new_cost_scenario_id=3,
        scenario_name='high_li_io_costs',
        scenario_description='Like default but high Li-ion costs',
        project_period_lifetimes_costs=high_storage_costs
    )


# TODO: some discrepancy between spreadsheet and IRP Staff proposal -- make
#  sure it's not a big issue
def load_project_new_potentials():
    """
    Max potentials for candidate projects
    :return:
    """
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        project_first_year = dict()
        for row in range(11 - 1, 52):
            project_first_year[rows_list[row][3 - 1]] = \
                int(float(rows_list[row][29- 1]))

    with open(os.path.join("cpuc_irp_data", "csvs",
                           "environmental_screens_from_irp_staff_proposal.csv"
                           ), "r") as f:
        rows_list = list(csv.reader(f))
        scenario_proj_potential = OrderedDict()

        for column in range(3 - 1, 8):
            scenario = rows_list[1 - 1][column]
            scenario_proj_potential[scenario] = OrderedDict()
            for row in range(2 - 1, 43):
                project = rows_list[row][9 - 1]
                scenario_proj_potential[scenario][project] = OrderedDict()
                for yr in range(2015, 2050 + 1):
                    scenario_proj_potential[scenario][project][yr] = \
                        (None, None, float(rows_list[row][column]), None) \
                        if yr >= project_first_year[project] \
                        else (None, None, 0, None)

    # Add pumped storage
    for scenario in scenario_proj_potential.keys():
        scenario_proj_potential[scenario]['CAISO_New_Pumped_Storage'] = \
            OrderedDict()
        for yr in range(2015, 2050 + 1):
            if yr >= 2025:
                scenario_proj_potential[scenario][
                    'CAISO_New_Pumped_Storage'][yr] = (None, None, 4000, None)
            elif yr < 2020:
                scenario_proj_potential[scenario][
                    'CAISO_New_Pumped_Storage'][yr] = (None, None, 0, None)
            else:
                scenario_proj_potential[scenario][
                    'CAISO_New_Pumped_Storage'][yr] = (None, None, 2000, None)

    # Add DR
    with open(os.path.join("cpuc_irp_data", "csvs", "DR_Conventional.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for scenario in scenario_proj_potential.keys():
            for row in range(60 - 1, 67):
                project = rows_list[row][3 - 1]
                scenario_proj_potential[scenario][project] = OrderedDict()
                for column in range(10 - 1, 45):
                    yr = int(float(rows_list[59 - 1][column]))
                    scenario_proj_potential[scenario][project][yr] = \
                        (None, None, float(rows_list[row][column]), None)

    # Update 'default' with numbers directly from RESOLVE Sep 2017 cases
    for prj in scenario_proj_potential['DRECP_SJV']:
        for yr in scenario_proj_potential['DRECP_SJV'][prj]:
            scenario_proj_potential['DRECP_SJV'][prj][yr] = \
                (None, None, None, None)
    scenario_proj_potential['DRECP_SJV']['CAISO_New_Li_Battery'] = \
        OrderedDict()
    scenario_proj_potential['DRECP_SJV']['CAISO_New_Flow_Battery'] = \
        OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "sep2017_drecp_sjv.csv"), "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            scenario_proj_potential['DRECP_SJV'][row[0]][int(row[1])] = \
                (None, None, float(row[2]), None)

    # Greater Imperial Wind is not a resource in Sep 2017 RESOLVE data,
    # so make potential 0 here
    for yr in range(2015, 2050 + 1):
        scenario_proj_potential['DRECP_SJV']['Greater_Imperial_Wind'][yr] = \
            (None, None, 0, None)
        scenario_proj_potential['DRECP_SJV']['Southern_California_Desert_Wind'
        ][yr] = \
            (None, None, 0, None)
        scenario_proj_potential['DRECP_SJV']['Kramer_Inyokern_Wind'][yr] = \
            (None, None, 0, None)

    # Create scenario with DRECP_SJV and additional limits preventing LCR
    # resources from being built before 2026
    lcr_limits_scenario = OrderedDict()
    for prj in scenario_proj_potential['DRECP_SJV']:
        lcr_limits_scenario[prj] = OrderedDict()
        for yr in [2018, 2022, 2026, 2030]:
            lcr_limits_scenario[prj][yr] = scenario_proj_potential[
                    'DRECP_SJV'][prj][yr]

    for prj in [
        'Battery_LCR_Bay_Area',
        'Battery_LCR_Big_Creek_Ventura',
        'Battery_LCR_Fresno',
        'Battery_LCR_Humboldt',
        'Battery_LCR_Kern',
        'Battery_LCR_LA_Basin',
        'Battery_LCR_NCNB',
        'Battery_LCR_San_Diego',
        'Battery_LCR_Sierra',
        'Battery_LCR_Stockton',
        'Advanced_CCGT_LCR_Bay_Area',
        'Advanced_CCGT_LCR_Big_Creek_Ventura',
        'Advanced_CCGT_LCR_Fresno',
        'Advanced_CCGT_LCR_Humboldt',
        'Advanced_CCGT_LCR_Kern',
        'Advanced_CCGT_LCR_LA_Basin',
        'Advanced_CCGT_LCR_NCNB',
        'Advanced_CCGT_LCR_San_Diego',
        'Advanced_CCGT_LCR_Sierra',
        'Advanced_CCGT_LCR_Stockton',
        'Aero_CT_LCR_Bay_Area',
        'Aero_CT_LCR_Big_Creek_Ventura',
        'Aero_CT_LCR_Fresno',
        'Aero_CT_LCR_Humboldt',
        'Aero_CT_LCR_Kern',
        'Aero_CT_LCR_LA_Basin',
        'Aero_CT_LCR_NCNB',
        'Aero_CT_LCR_San_Diego',
        'Aero_CT_LCR_Sierra',
        'Aero_CT_LCR_Stockton'
    ]:
        lcr_limits_scenario[prj] = OrderedDict()
        for yr in [2018, 2022]:
            lcr_limits_scenario[prj][yr] = (None, None, 0, None)

    scenario_proj_potential['DRECP_SJV_w_LCR_project_limits'] = \
        lcr_limits_scenario

    # Insert data
    scenario_id = 1
    for scenario in scenario_proj_potential.keys():
        project_new_potentials.update_project_potentials(
            io=io, c=c2,
            project_new_potential_scenario_id=scenario_id,
            scenario_name=scenario,
            scenario_description=scenario,
            project_period_potentials=scenario_proj_potential[scenario]
        )

        scenario_id += 1

    # Update New_Mexico and Wyoming wind to be 1500 MW of potential starting
        # in 2026
    c2.execute(
        """UPDATE inputs_project_new_potential
        SET maximum_cumulative_new_build_mw = 1500
        WHERE (project = 'Wyoming_Wind' or project = 'New_Mexico_Wind')
        AND period > 2025;"""
    )
    io.commit()


# def load_project_elcc_chars():
#     """
#     ELCC charecteristics of projects
#     :return:
#     """
#     proj_prm_type = dict()
#     proj_elcc_simple_fraction = dict()
#     proj_elcc_surface = dict()
#     proj_min_duration_for_full = dict()
#     proj_deliv_group = dict()
#
#     # Start with the simple fraction
#     with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Planning_Reserve.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#
#         for row in range(22 - 1, 50):
#             if rows_list[row][3 - 1] == 'CAISO_Reciprocating_Engine':
#                 for tail in ["", "_Candidate"]:
#                     proj_elcc_simple_fraction[
#                         rows_list[row][3 - 1] + tail
#                     ] = float(rows_list[row][4 - 1])
#             elif rows_list[row][3 - 1] in [
#                 'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
#             ]:
#                 for tail in ["", "_Planned"]:
#                     proj_elcc_simple_fraction[
#                         rows_list[row][3 - 1] + tail
#                     ] = float(rows_list[row][4 - 1])
#             proj_elcc_simple_fraction[
#                 rows_list[row][3 - 1]
#             ] = float(rows_list[row][4 - 1])
#
#     # Disaggregated gas projects
#     with open(os.path.join("cpuc_irp_data", "gridpath_specific", "gas_disagg_plants.csv"),
#               "r") as foo:
#         reader = csv.reader(foo, delimiter=",")
#         next(reader)
#         for row in reader:
#             if row[0].startswith("Battery"):
#                 pass  # will add the batteries separately with storage below
#             else:
#                 for yr in range(2015, 2050+1):
#                     proj_elcc_simple_fraction[row[0]] = \
#                         float(row[14]) if float(row[14]) < 1 else 1.0
#
#     with open(os.path.join("cpuc_irp_data", "gridpath_specific",
#                            "plants_to_add_for_ccgt_analysis.csv"),
#               "r") as foo:
#         reader = csv.reader(foo, delimiter=",")
#         next(reader)
#         for row in reader:
#             for yr in range(2015, 2050+1):
#                 proj_elcc_simple_fraction[row[0]] = \
#                     float(row[14]) if float(row[14]) < 1 else 1.0
#
#     # PRM TYPE
#     # All simple-fraction projects so far are 'fully_deliverable'
#     for proj in proj_elcc_simple_fraction.keys():
#         proj_prm_type[proj] = 'fully_deliverable'
#
#     # Storage projects will be of the type
#     # 'fully_deliverable_energy_limited," have a simple ELCC fraction of 1,
#     # and get a min duration for full credit of 4
#     storage_projects = [
#         'CAISO_Existing_Pumped_Storage',
#         'CAISO_Storage_Mandate',
#         'CAISO_New_Pumped_Storage',
#         'CAISO_New_Flow_Battery',
#         'CAISO_New_Li_Battery',
#         'Battery_LCR_Bay_Area',
#         'Battery_LCR_Big_Creek_Ventura',
#         'Battery_LCR_Fresno',
#         'Battery_LCR_Humboldt',
#         'Battery_LCR_Kern',
#         'Battery_LCR_LA_Basin',
#         'Battery_LCR_NCNB',
#         'Battery_LCR_San_Diego',
#         'Battery_LCR_Sierra',
#         'Battery_LCR_Stockton'
#     ]
#
#     for stor_proj in storage_projects:
#         proj_prm_type[stor_proj] = 'fully_deliverable_energy_limited'
#         proj_elcc_simple_fraction[stor_proj] = 1
#         proj_min_duration_for_full[stor_proj] = 4
#
#     # Existing renewable resources are 'fully_deliverable' and have been
#     # included in the SYS_Planning_Reserve tab, imported above, except for
#     # CAISO_Solar_for_CAISO, CAISO_Wind_for_CAISO, CAISO_Solar_for_Other,
#     # and CAISO_Wind_for_Other
#     # Add the latter two here
#     proj_prm_type['CAISO_Solar_for_CAISO'] = 'fully_deliverable'
#     proj_prm_type['CAISO_Wind_for_CAISO'] = 'fully_deliverable'
#     proj_prm_type['CAISO_Solar_for_Other'] = 'fully_deliverable'
#     proj_prm_type['CAISO_Wind_for_Other'] = 'fully_deliverable'
#
#     # New renewable resources will be 'energy_only_allowed'
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Tx_Costs.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#         for row in range(28-1, 69):
#             proj_prm_type[rows_list[row][3 - 1]] = 'energy_only_allowed'
#
#     # Customer PV will be fully deliverable
#     proj_prm_type['Customer_PV'] = 'fully_deliverable'
#
#     # ELCC SURFACE
#     # Solar and wind resources will contribute to the ELCC surface
#     # Existing CAISO resources
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Baseline.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#         for row in range(59 - 1, 73):
#             proj = rows_list[row][3 - 1]
#             zone = rows_list[row][5 - 1]
#             tech = rows_list[row][7 - 1]
#             if zone == 'CAISO' and tech in ["Wind", "Solar"]:
#                 proj_elcc_surface[proj] = 1
#             else:
#                 pass
#
#     # Also add CAISO_Solar_for_Other and CAISO_Wind_for_Other
#     proj_elcc_surface['CAISO_Solar_for_Other'] = 1
#     proj_elcc_surface['CAISO_Wind_for_Other'] = 1
#
#     # New resources
#     # We'll add non-OOS resources and existing Tx OOS resources to the surface
#     # OOS resources that require new transmission will get a simple fraction
#     # for the import capacity instead
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#         for row in range(11 - 1, 52):
#             proj = rows_list[row][3 - 1]
#             tech = rows_list[row][9 - 1]
#             existing_tx = int(float(rows_list[row][13 - 1]))
#             if tech in ["Wind", "Solar"] and existing_tx:
#                 proj_elcc_surface[proj] = 1
#             else:
#                 pass
#     # It also looks like Baja resources are treated via the surface, not via
#     #  the new transmission simple fraction
#     proj_elcc_surface['Baja_California_Solar'] = 1
#     proj_elcc_surface['Baja_California_Wind'] = 1
#
#     # Northern California wind is missing because it's not allowed in any
#     # portfolios, but add it here nonetheless
#     proj_elcc_surface['Northern_California_Wind'] = 1
#
#     # Customer_PV
#     proj_elcc_surface['Customer_PV'] = 1
#
#     # OOS projects requiring new transmission contribute via a simple
#     # fraction (credit for the additional import capacity)
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Tx_Costs.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#
#         for row in range(28-1, 69):
#             proj = rows_list[row][3 - 1]
#             apply_credit = int(float(rows_list[row][7 - 1]))
#             deliverability_group = rows_list[row][5 - 1]
#
#             # Additional credit for transmission
#             if apply_credit:
#                 # check that we haven't already included this project
#                 if proj in proj_elcc_simple_fraction.keys():
#                     warnings.warn(
#                         """Project {} already has an ELCC simple fraction. This
#                         will overwrite it."""
#                     ).format(proj)
#                 proj_elcc_simple_fraction[proj] = 0.6
#             else:
#                 pass
#
#             # Deliverability groups
#             proj_deliv_group[proj] = deliverability_group
#
#     # CLEAN-UP
#     # Go over all projects and assign 0 for ELCC surface contribution if not
#     # contributing via surface and 0 for simple contribution if not
#     # contributing via simple fraction
#     for proj in proj_prm_type.keys():
#         # Surface
#         if proj in proj_elcc_surface.keys():
#             pass
#         else:
#             proj_elcc_surface[proj] = 0
#         # Simple
#         if proj in proj_elcc_simple_fraction.keys():
#             pass
#         else:
#             proj_elcc_simple_fraction[proj] = 0
#
#     project_prm.project_elcc_chars(
#         io=io, c=c2,
#         project_elcc_chars_scenario_id=1,
#         scenario_name="default",
#         scenario_description="""Conventional projects and baseload
#         renewables contribute via simple fraction; wind and solar
#         contribute via surface; additional contribution for imports for OOS
#         wind and solar; new renewables can be energy only; existing
#         resources are fully deliverable; storage min duration for full
#         credit is 4 hours.""",
#         proj_prm_type=proj_prm_type,
#         proj_elcc_simple_fraction=proj_elcc_simple_fraction,
#         proj_elcc_surface=proj_elcc_surface,
#         proj_min_duration_for_full=proj_min_duration_for_full,
#         proj_deliv_group=proj_deliv_group
#     )
#
#
# def load_local_capacity_chars():
#     """
#
#     :return:
#     """
#     print("local capaciy chars")
#     c2.execute(
#         """INSERT INTO subscenarios_project_local_capacity_chars
#         (project_local_capacity_chars_scenario_id, name, description)
#         VALUES (1, 'default', 'default');"""
#     )
#
#     with open(os.path.join("cpuc_irp_data", "gridpath_specific",
#                            "gas_disagg_plants.csv"),
#               "r") as f:
#         reader = csv.reader(f, delimiter=",")
#         next(reader)
#         for row in reader:
#             c2.execute(
#                 """INSERT INTO inputs_project_local_capacity_chars
#                 (project_local_capacity_chars_scenario_id, project,
#                 local_capacity_fraction)
#                 VALUES ({}, '{}', {});""".format(
#                     1, row[0], float(row[14]) if float(row[14]) < 1 else 1.0
#                 )
#             )
#     io.commit()
#
#
# def load_deliverability_group_params():
#     """
#
#     :return:
#     """
#     deliv_group_params = dict()
#     with open(os.path.join("cpuc_irp_data", "csvs", "REN_Tx_Costs.csv"), "r") as f:
#         rows_list = list(csv.reader(f))
#
#         for row in range(77 - 1, 88):
#             group = rows_list[row][3 - 1]
#             no_cost_cap = float(rows_list[row][6 - 1])
#             cost = float(rows_list[row][7 - 1]) * 1000  # conver to per MW-yr
#             eo_cap_lim = float(rows_list[row][8 - 1])
#
#             deliv_group_params[group] = (no_cost_cap, cost, eo_cap_lim)
#
#     # Insert data
#     project_prm.deliverability_groups(
#     io=io, c=c2,
#     prm_energy_only_scenario_id=1,
#     scenario_name='default',
#     scenario_description='Default deliverability groups and parameters.',
#     deliv_group_params=deliv_group_params
#     )
def load_project_elcc_chars():
    """
    ELCC charecteristics of projects
    :return:
    """
    proj_prm_type = dict()
    proj_elcc_simple_fraction = dict()
    proj_elcc_surface = dict()
    proj_min_duration_for_full = dict()
    proj_deliv_group = dict()

    # Start with the simple fraction
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Planning_Reserve.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(22 - 1, 50):
            if rows_list[row][3 - 1] == 'CAISO_Reciprocating_Engine':
                for tail in ["", "_Candidate"]:
                    proj_elcc_simple_fraction[
                        rows_list[row][3 - 1] + tail
                        ] = float(rows_list[row][4 - 1])
            elif rows_list[row][3 - 1] in [
                'CAISO_Peaker1', 'CAISO_Peaker2', 'CAISO_CCGT1'
            ]:
                for tail in ["", "_Planned"]:
                    proj_elcc_simple_fraction[
                        rows_list[row][3 - 1] + tail
                        ] = float(rows_list[row][4 - 1])
            proj_elcc_simple_fraction[
                rows_list[row][3 - 1]
            ] = float(rows_list[row][4 - 1])

    # Disaggregated gas projects
    with open(os.path.join("cpuc_irp_data", "gridpath_specific", "gas_disagg_plants.csv"),
              "r") as foo:
        reader = csv.reader(foo, delimiter=",")
        next(reader)
        for row in reader:
            if row[0].startswith("Battery"):
                pass  # will add the batteries separately with storage below
            else:
                for yr in range(2015, 2050 + 1):
                    proj_elcc_simple_fraction[row[0]] = \
                        float(row[14]) if float(row[14]) < 1 else 1.0

    # PRM TYPE
    # All simple-fraction projects so far are 'fully_deliverable'
    for proj in proj_elcc_simple_fraction.keys():
        proj_prm_type[proj] = 'fully_deliverable'

    # Storage projects will be of the type
    # 'fully_deliverable_energy_limited," have a simple ELCC fraction of 1,
    # and get a min duration for full credit of 4
    storage_projects = [
        'CAISO_Existing_Pumped_Storage',
        'CAISO_Storage_Mandate',
        'CAISO_New_Pumped_Storage',
        'CAISO_New_Flow_Battery',
        'CAISO_New_Li_Battery',
        'Battery_LCR_Bay_Area',
        'Battery_LCR_Big_Creek_Ventura',
        'Battery_LCR_Fresno',
        'Battery_LCR_Humboldt',
        'Battery_LCR_Kern',
        'Battery_LCR_LA_Basin',
        'Battery_LCR_NCNB',
        'Battery_LCR_San_Diego',
        'Battery_LCR_Sierra',
        'Battery_LCR_Stockton'
    ]

    for stor_proj in storage_projects:
        proj_prm_type[stor_proj] = 'fully_deliverable_energy_limited'
        proj_elcc_simple_fraction[stor_proj] = 1
        proj_min_duration_for_full[stor_proj] = 4

    # Existing renewable resources are 'fully_deliverable' and have been
    # included in the SYS_Planning_Reserve tab, imported above, except for
    # CAISO_Solar_for_CAISO, CAISO_Wind_for_CAISO, CAISO_Solar_for_Other,
    # and CAISO_Wind_for_Other
    # Add the latter two here
    proj_prm_type['CAISO_Solar_for_CAISO'] = 'fully_deliverable'
    proj_prm_type['CAISO_Wind_for_CAISO'] = 'fully_deliverable'
    proj_prm_type['CAISO_Solar_for_Other'] = 'fully_deliverable'
    proj_prm_type['CAISO_Wind_for_Other'] = 'fully_deliverable'

    # New renewable resources will be 'energy_only_allowed'
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Tx_Costs.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(28 - 1, 69):
            proj_prm_type[rows_list[row][3 - 1]] = 'energy_only_allowed'

    # Customer PV will be fully deliverable
    proj_prm_type['Customer_PV'] = 'fully_deliverable'

    # ELCC SURFACE
    # Solar and wind resources will contribute to the ELCC surface
    # Existing CAISO resources
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Baseline.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(59 - 1, 73):
            proj = rows_list[row][3 - 1]
            zone = rows_list[row][5 - 1]
            tech = rows_list[row][7 - 1]
            if zone == 'CAISO' and tech in ["Wind", "Solar"]:
                proj_elcc_surface[proj] = 1
            else:
                pass

    # Also add CAISO_Solar_for_Other and CAISO_Wind_for_Other
    proj_elcc_surface['CAISO_Solar_for_Other'] = 1
    proj_elcc_surface['CAISO_Wind_for_Other'] = 1

    # New resources
    # We'll add non-OOS resources and existing Tx OOS resources to the surface
    # OOS resources that require new transmission will get a simple fraction
    # for the import capacity instead
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Candidate.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(11 - 1, 52):
            proj = rows_list[row][3 - 1]
            tech = rows_list[row][9 - 1]
            existing_tx = int(float(rows_list[row][13 - 1]))
            if tech in ["Wind", "Solar"] and existing_tx:
                proj_elcc_surface[proj] = 1
            else:
                pass
    # It also looks like Baja resources are treated via the surface, not via
    #  the new transmission simple fraction
    proj_elcc_surface['Baja_California_Solar'] = 1
    proj_elcc_surface['Baja_California_Wind'] = 1

    # Northern California wind is missing because it's not allowed in any
    # portfolios, but add it here nonetheless
    proj_elcc_surface['Northern_California_Wind'] = 1

    # Customer_PV
    proj_elcc_surface['Customer_PV'] = 1

    # OOS projects requiring new transmission contribute via a simple
    # fraction (credit for the additional import capacity)
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Tx_Costs.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(28 - 1, 69):
            proj = rows_list[row][3 - 1]
            apply_credit = int(float(rows_list[row][7 - 1]))
            deliverability_group = rows_list[row][5 - 1]

            # Additional credit for transmission
            if apply_credit:
                # check that we haven't already included this project
                if proj in proj_elcc_simple_fraction.keys():
                    warnings.warn(
                        """Project {} already has an ELCC simple fraction. This 
                        will overwrite it."""
                    ).format(proj)
                proj_elcc_simple_fraction[proj] = 0.6
            else:
                pass

            # Deliverability groups
            proj_deliv_group[proj] = deliverability_group

    # CLEAN-UP
    # Go over all projects and assign 0 for ELCC surface contribution if not
    # contributing via surface and 0 for simple contribution if not
    # contributing via simple fraction
    for proj in proj_prm_type.keys():
        # Surface
        if proj in proj_elcc_surface.keys():
            pass
        else:
            proj_elcc_surface[proj] = 0
        # Simple
        if proj in proj_elcc_simple_fraction.keys():
            pass
        else:
            proj_elcc_simple_fraction[proj] = 0

    project_prm.project_elcc_chars(
        io=io, c=c2,
        project_elcc_chars_scenario_id=1,
        scenario_name="default",
        scenario_description="""Conventional projects and baseload 
        renewables contribute via simple fraction; wind and solar 
        contribute via surface; additional contribution for imports for OOS 
        wind and solar; new renewables can be energy only; existing 
        resources are fully deliverable; storage min duration for full 
        credit is 4 hours.""",
        proj_prm_type=proj_prm_type,
        proj_elcc_simple_fraction=proj_elcc_simple_fraction,
        proj_elcc_surface=proj_elcc_surface,
        proj_min_duration_for_full=proj_min_duration_for_full,
        proj_deliv_group=proj_deliv_group
    )


def load_local_capacity_chars():
    """

    :return:
    """
    print("local capaciy chars")
    c2.execute(
        """INSERT INTO subscenarios_project_local_capacity_chars
        (project_local_capacity_chars_scenario_id, name, description)
        VALUES (1, 'default', 'default');"""
    )

    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "gas_disagg_plants.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            c2.execute(
                """INSERT INTO inputs_project_local_capacity_chars
                (project_local_capacity_chars_scenario_id, project, 
                local_capacity_fraction)
                VALUES ({}, '{}', {});""".format(
                    1, row[0], float(row[14]) if float(row[14]) < 1 else 1.0
                )
            )
    io.commit()


def load_deliverability_group_params():
    """

    :return:
    """
    deliv_group_params = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Tx_Costs.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(77 - 1, 88):
            group = rows_list[row][3 - 1]
            no_cost_cap = float(rows_list[row][6 - 1])
            cost = float(rows_list[row][7 - 1]) * 1000  # conver to per MW-yr
            eo_cap_lim = float(rows_list[row][8 - 1])

            deliv_group_params[group] = (no_cost_cap, cost, eo_cap_lim)

    # Insert data
    project_prm.deliverability_groups(
        io=io, c=c2,
        prm_energy_only_scenario_id=1,
        scenario_name='default',
        scenario_description='Default deliverability groups and parameters.',
        deliv_group_params=deliv_group_params
    )


def load_elcc_surface():
    """
    GridPath multiplies the coefficients by MW capacity
    In RESOLVE the coefficients are multiplied by:
    main zone peak load x solar fraction of annual load =
    main zone peak load x MW x annual cap factor x 8760 / annual load
    We'll therefore adjust the coefficients for GridPath as follows:
    GridPath coefficient =
    RESOLVE coefficient x peak load x cap factor x 8760 / annual load

    The intercept in GridPath has units of MW and in RESOLVE it's fraction
    of peak, so:
    GridPath intercept = RESOLVE intercept x peak load
    :return:
    """
    # Get the raw surface as input to RESOLVE
    resolve_surface = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Planning_Reserve.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(59 - 1, 74):
            resolve_surface[int(float(rows_list[row][3 - 1]))] = {
                "solar_coefficient": float(rows_list[row][4 - 1]),
                "wind_coefficient": float(rows_list[row][5 - 1]),
                "intercept": float(rows_list[row][6 - 1])
            }

    # Manually get the peak load requirement at the default load assumptions
    peak_load = {
        2018: 46404.0,
        2022: 45815.0,
        2026: 45568.0,
        2030: 45624.0
    }

    # Manually get the annual load at the default load assumptions
    annual_load = {
        2018: 235158857.0,
        2022: 235493690.0,
        2026: 238486926.0,
        2030: 242473023.0

    }

    # Get the capacity factors for variable projects
    proj_cf = dict()
    with open(os.path.join("cpuc_irp_data", "csvs", "REN_Profiles.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for column in range(5 - 1, 61):
            proj_cf[rows_list[21 - 1][column]] = \
                float(rows_list[23 - 1][column])

    # Make the GridPath intercepts
    gridpath_intercepts = OrderedDict()
    for yr in peak_load.keys():
        gridpath_intercepts[yr] = OrderedDict()
        for facet in resolve_surface.keys():
            gridpath_intercepts[yr][facet] = \
                peak_load[yr] * resolve_surface[facet]["intercept"]

    zone_period_facet_intercepts = {'CAISO': gridpath_intercepts}

    # Assign coefficients to projects
    proj_period_facet_coeff = OrderedDict()
    solar_projects = c2.execute(
        """SELECT project
        FROM inputs_project_elcc_chars
        WHERE contributes_to_elcc_surface = 1
        AND project LIKE '%Solar%';"""
    ).fetchall()
    wind_projects = c2.execute(
        """SELECT project
        FROM inputs_project_elcc_chars
        WHERE contributes_to_elcc_surface = 1
        AND project LIKE '%Wind%';"""
    ).fetchall()

    for p in solar_projects:
        project = p[0]
        proj_period_facet_coeff[project] = OrderedDict()
        for period in [2018, 2022, 2026, 2030]:
            proj_period_facet_coeff[project][period] = OrderedDict()
            for facet in resolve_surface.keys():
                proj_period_facet_coeff[project][period][facet] = \
                    resolve_surface[facet]["solar_coefficient"] \
                    * peak_load[period] \
                    * proj_cf[project] * 8760.0 \
                    / annual_load[period]

    # Customer_PV
    proj_period_facet_coeff['Customer_PV'] = OrderedDict()
    for period in [2018, 2022, 2026, 2030]:
        proj_period_facet_coeff['Customer_PV'][period] = OrderedDict()
        for facet in resolve_surface.keys():
            proj_period_facet_coeff['Customer_PV'][period][facet] = \
                resolve_surface[facet]["solar_coefficient"] \
                    * peak_load[period] \
                    * proj_cf['Customer_PV'] * 8760.0 \
                    / annual_load[period]

    for p in wind_projects:
        project = p[0]
        proj_period_facet_coeff[project] = OrderedDict()
        for period in [2018, 2022, 2026, 2030]:
            proj_period_facet_coeff[project][period] = OrderedDict()
            for facet in resolve_surface.keys():
                proj_period_facet_coeff[project][period][facet] = \
                    resolve_surface[facet]["wind_coefficient"] \
                    * peak_load[period] \
                    * proj_cf[project] * 8760.0 \
                    / annual_load[period]

    project_prm.elcc_surface(
        io=io, c=c2,
        prm_zone_scenario_id=1,
        project_prm_zone_scenario_id=1,
        elcc_surface_scenario_id=1,
        scenario_name='default',
        scenario_description='Default surface contributions.',
        zone_period_facet_intercepts=zone_period_facet_intercepts,
        proj_period_facet_coeff=proj_period_facet_coeff
    )


def load_transmission_portfolios():
    """
    Transmission portfolios
    :return:
    """
    tx_line_cap_types = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Regional_Settings.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(11 - 1, 23):
            tx_line_cap_types[rows_list[row][3 - 1]] = 'specified_transmission'

    transmission_portfolios.insert_transmission_portfolio(
        io=io, c=c2,
        transmission_portfolio_scenario_id=1,
        scenario_name='default_tx_lines',
        scenario_description='Default transmission lines',
        tx_line_cap_types=tx_line_cap_types
    )


def load_transmission_load_zones():
    """

    :return:
    """
    tx_line_load_zones = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Regional_Settings.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(11 - 1, 23):
            tx_line_load_zones[rows_list[row][3 - 1]] = \
                (rows_list[row][4 - 1], rows_list[row][5 - 1])

    transmission_zones.insert_transmission_load_zones(
        io=io, c=c2,
        load_zone_scenario_id=1,
        transmission_load_zone_scenario_id=1,
        scenario_name='default_tx_load_zone_geography',
        scenario_description='Default transmission line geography.',
        tx_line_load_zones=tx_line_load_zones
    )


def load_transmission_capacities():
    """

    :return:
    """
    tx_line_period_capacities = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Regional_Settings.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(11 - 1, 23):
            tx_line_period_capacities[rows_list[row][3 - 1]] = OrderedDict()
            for yr in range(2015, 2050+1):
                tx_line_period_capacities[rows_list[row][3 - 1]][yr] = \
                    (rows_list[row][6 - 1], rows_list[row][7 - 1])

    # Insert data
    transmission_capacities.insert_transmission_capacities(
        io=io, c=c2,
        transmission_existing_capacity_scenario_id=1,
        scenario_name='default_tx_capacities',
        scenario_description='Default transmission line capacities.',
        tx_line_period_capacities=tx_line_period_capacities
    )


def load_simultaneous_flow_groups():
    """
    Groups of transmission lines over which simultaneous flow constraints
    will be applied
    :return:
    """
    group_lines = OrderedDict()

    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Regional_Settings.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        resolve_to_gridpath_group_names = {
            'NW_to_CA': [('NW_to_CA_Imports', 1), ('CA_to_NW_Exports', -1)],
            'CAISO_Simultaneous_Import': [
                ("CAISO_Imports", 1), ('CAISO_Exports', -1)]
        }

        for resolve_name in resolve_to_gridpath_group_names.keys():
            for gridpath_name in resolve_to_gridpath_group_names[resolve_name]:
                group_lines[gridpath_name[0]] = list()
                for row in range(95 - 1, 102):
                    if rows_list[row][3 - 1] == resolve_name:
                        group_lines[gridpath_name[0]].append(
                            (rows_list[row][4 - 1],
                             int(float(rows_list[row][5 - 1]))
                             * gridpath_name[1])
                        )

    # Insert data
    simultaneous_flow_groups.insert_transmission_simultaneous_flow_groups(
        io=io, c=c2,
        transmission_simultaneous_flow_limit_line_group_scenario_id=1,
        scenario_name="CA_to_NW_Imports_Exports_and_CAISO_Imports_Exports",
        scenario_description="Lines between CA and NW, lines in and out of "
                             "CAISO",
        group_lines=group_lines
    )


def load_simultaneous_flow_limits():
    """
    Simultanous flow limits on groups of transmission lines
    :return:
    """
    scenario_group_period_limits = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Regional_Settings.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # 'Low' scenario
        scenario_group_period_limits["Low_Limits"] = OrderedDict()

        # CA_to_NW_Exports
        scenario_group_period_limits["Low_Limits"]["CA_to_NW_Exports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Low_Limits"]["CA_to_NW_Exports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                - float(rows_list[108 - 1][column])

        # NW_to_CA_Imports
        scenario_group_period_limits["Low_Limits"]["NW_to_CA_Imports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Low_Limits"]["NW_to_CA_Imports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                float(rows_list[109 - 1][column])

        # CAISO_Exports
        scenario_group_period_limits["Low_Limits"]["CAISO_Exports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Low_Limits"]["CAISO_Exports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                - float(rows_list[110 - 1][column])

        # CAISO_Imports
        scenario_group_period_limits["Low_Limits"]["CAISO_Imports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Low_Limits"]["CAISO_Imports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                float(rows_list[111 - 1][column])

        # 'Mid' scenario
        scenario_group_period_limits["Mid_Limits"] = OrderedDict()

        # CA_to_NW_Exports
        scenario_group_period_limits["Mid_Limits"]["CA_to_NW_Exports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Mid_Limits"]["CA_to_NW_Exports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                - float(rows_list[114 - 1][column])

        # NW_to_CA_Imports
        scenario_group_period_limits["Mid_Limits"]["NW_to_CA_Imports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Mid_Limits"]["NW_to_CA_Imports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                float(rows_list[115 - 1][column])

        # CAISO_Exports
        scenario_group_period_limits["Mid_Limits"]["CAISO_Exports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Mid_Limits"]["CAISO_Exports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                - float(rows_list[116 - 1][column])

        # CAISO_Imports
        scenario_group_period_limits["Mid_Limits"]["CAISO_Imports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["Mid_Limits"]["CAISO_Imports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                float(rows_list[117 - 1][column])

        # 'High' scenario
        scenario_group_period_limits["High_Limits"] = OrderedDict()

        # CA_to_NW_Exports
        scenario_group_period_limits["High_Limits"]["CA_to_NW_Exports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["High_Limits"]["CA_to_NW_Exports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                - float(rows_list[120 - 1][column])

        # NW_to_CA_Imports
        scenario_group_period_limits["High_Limits"]["NW_to_CA_Imports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["High_Limits"]["NW_to_CA_Imports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                float(rows_list[121 - 1][column])

        # CAISO_Exports
        scenario_group_period_limits["High_Limits"]["CAISO_Exports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["High_Limits"]["CAISO_Exports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                - float(rows_list[122 - 1][column])

        # CAISO_Imports
        scenario_group_period_limits["High_Limits"]["CAISO_Imports"] = \
            OrderedDict()
        for column in range(6 - 1, 41):
            scenario_group_period_limits["High_Limits"]["CAISO_Imports"][
                int(float(rows_list[105 - 1][column]))
            ] = \
                float(rows_list[123 - 1][column])

    # Insert data
    scenario_id = 1
    for scenario in scenario_group_period_limits.keys():
        simultaneous_flows.insert_transmission_simultaneous_flow_limits(
            io=io, c=c2,
            transmission_simultaneous_flow_limit_scenario_id=scenario_id,
            scenario_name=scenario,
            scenario_description=scenario,
            group_period_limits=scenario_group_period_limits[scenario]
        )
        scenario_id += 1


def load_transmission_carbon_cap_zones():
    """
    Transmission carbon cap zones, direction, and intensity
    Positive direction into CAISO for all lines
    The emissions intensity does not vary by year
    :return:
    """
    tx_carbon_cap_zones = OrderedDict()

    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Regional_Settings.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        for row in range(85 - 1, 89):
            tx_carbon_cap_zones[rows_list[row][3-1]] = (
                'CAISO', 'positive', 0.42716568
            )

    # Insert data
    transmission_zones.insert_transmission_carbon_cap_zones(
        io=io, c=c2,
        carbon_cap_zone_scenario_id=defaults["carbon_cap_zone_scenario_id"],
        transmission_carbon_cap_zone_scenario_id=defaults["carbon_cap_zone_scenario_id"],
        scenario_name="default",
        scenario_description="All lines into CAISO, intensity = 0.427",
        tx_line_carbon_cap_zones=tx_carbon_cap_zones
    )


def load_transmission_hurdle_rates():
    """
    Hurdle rates (including carbon adder)
    :return:
    """

    # Base hurdle rate
    # Values will be tuples with the first element of the tuple the positive
    # direction base hurdle rate and the second element of the tuple the
    # negative direction base hurdle rate
    tx_line_base_hurdle_rates = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Regional_Settings.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # Base hurdle rates
        for row in range(45 - 1, 57):
            tx_line_base_hurdle_rates[rows_list[row][3 - 1]] = \
                (float(rows_list[row][6 - 1]),
                 float(rows_list[row][7 - 1]))

    # Calculate carbon adder
    # Get carbon prices
    with open(os.path.join("cpuc_irp_data", "csvs", "Sys_Fuel_Costs.csv"), "r") as f:
        carbon_rows_list = list(csv.reader(f))
    # Get the carbon adder by carbon price scenario
    scenario_year_carbon_price = OrderedDict()
    for row in range(53 - 1, 56):
        scenario_year_carbon_price[carbon_rows_list[row][3 - 1]] = OrderedDict()
        for column in range(8 - 1, 43):
            scenario_year_carbon_price[
                carbon_rows_list[row][3 - 1]
            ][
                int(float(carbon_rows_list[51 - 1][column]))
            ] = \
                float(carbon_rows_list[row][column])
    # Get carbon intensities
    # Lines into California will get the adder (but not lines within CA)
    # All lines have the same adder in all years for now
    carbon_adder_lines_intensities = {
        'SW_to_CAISO': 0.42716568,
        'NW_to_CAISO': 0.42716568,
        'NW_to_LDWP': 0.42716568,
        'SW_to_LDWP': 0.42716568,
        'SW_to_IID': 0.42716568,
        'NW_to_BANC': 0.42716568,
        'SW_to_BANC': 0.42716568
    }

    # Three hurdle rate scenarios
    # Base X 3 carbon price scenarios

    # Low carbon price
    low_carbon_price_hurdle_rates = OrderedDict()
    for tx_line in tx_line_base_hurdle_rates.keys():
        low_carbon_price_hurdle_rates[tx_line] = OrderedDict()
        for period in scenario_year_carbon_price["Low"].keys():
            # Adder applied in positive direction only, as all lines are
            # defined with the positive direction into California,
            # so no need to adjust any negative directions
            low_carbon_price_hurdle_rates[tx_line][period] = (
                (tx_line_base_hurdle_rates[tx_line][0] + (
                 scenario_year_carbon_price["Low"][period] *
                 carbon_adder_lines_intensities[tx_line]
                 if tx_line in carbon_adder_lines_intensities.keys() else 0)),
                tx_line_base_hurdle_rates[tx_line][1]
            )

    # Mid carbon price
    mid_carbon_price_hurdle_rates = OrderedDict()
    for tx_line in tx_line_base_hurdle_rates.keys():
        mid_carbon_price_hurdle_rates[tx_line] = OrderedDict()
        for period in scenario_year_carbon_price["Mid"].keys():
            # Adder applied in positive direction only, as all lines are
            # defined with the positive direction into California,
            # so no need to adjust any negative directions
            mid_carbon_price_hurdle_rates[tx_line][period] = (
                (tx_line_base_hurdle_rates[tx_line][0] + (
                 scenario_year_carbon_price["Mid"][period] *
                 carbon_adder_lines_intensities[tx_line]
                 if tx_line in carbon_adder_lines_intensities.keys() else 0)),
                tx_line_base_hurdle_rates[tx_line][1]
            )

    # High carbon price
    high_carbon_price_hurdle_rates = OrderedDict()
    for tx_line in tx_line_base_hurdle_rates.keys():
        high_carbon_price_hurdle_rates[tx_line] = OrderedDict()
        for period in scenario_year_carbon_price["High"].keys():
            # Adder applied in positive direction only, as all lines are
            # defined with the positive direction into California,
            # so no need to adjust any negative directions
            high_carbon_price_hurdle_rates[tx_line][period] = (
                (tx_line_base_hurdle_rates[tx_line][0] + (
                 scenario_year_carbon_price["High"][period] *
                 carbon_adder_lines_intensities[tx_line]
                 if tx_line in carbon_adder_lines_intensities.keys() else 0)),
                tx_line_base_hurdle_rates[tx_line][1]
            )

    # No carbon price
    no_carbon_price_hurdle_rates = OrderedDict()
    for tx_line in tx_line_base_hurdle_rates.keys():
        no_carbon_price_hurdle_rates[tx_line] = OrderedDict()
        for period in scenario_year_carbon_price["High"].keys():
            no_carbon_price_hurdle_rates[tx_line][period] = (
                tx_line_base_hurdle_rates[tx_line][0],
                tx_line_base_hurdle_rates[tx_line][1]
            )

    # Insert data
    # Low
    transmission_hurdle_rates.insert_transmission_hurdle_rates(
        io=io, c=c2,
        transmission_hurdle_rate_scenario_id=1,
        scenario_name="low_carbon_prices",
        scenario_description="Base hurdle rates plus low carbon price carbon "
                             "adders into California",
        tx_line_period_hurdle_rates=low_carbon_price_hurdle_rates
    )

    # Mid
    transmission_hurdle_rates.insert_transmission_hurdle_rates(
        io=io, c=c2,
        transmission_hurdle_rate_scenario_id=2,
        scenario_name="mid_carbon_prices",
        scenario_description="Base hurdle rates plus mid carbon price carbon "
                             "adders into California",
        tx_line_period_hurdle_rates=mid_carbon_price_hurdle_rates
    )

    # High
    transmission_hurdle_rates.insert_transmission_hurdle_rates(
        io=io, c=c2,
        transmission_hurdle_rate_scenario_id=3,
        scenario_name="high_carbon_prices",
        scenario_description="Base hurdle rates plus high carbon price carbon "
                             "adders into California",
        tx_line_period_hurdle_rates=high_carbon_price_hurdle_rates
    )

    # No carbon adder
    transmission_hurdle_rates.insert_transmission_hurdle_rates(
        io=io, c=c2,
        transmission_hurdle_rate_scenario_id=4,
        scenario_name="no_carbon_prices",
        scenario_description="Base hurdle rates only",
        tx_line_period_hurdle_rates=no_carbon_price_hurdle_rates
    )


def load_carbon_cap_targets():
    """
    Carbon cap targets by period
    :return:
    """
    scenario_zone_period_targets = OrderedDict()

    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_RPS_GHG_Targets_JULY_UPDATE.csv"),
              "r") as f:
        rows_list = list(csv.reader(f))

        for row in list(range(58 - 1, 61)) + [66 - 1]:
            if rows_list[row][3 - 1] == 'Very Large':
                scenario = 'Very_Large'
            elif rows_list[row][3 - 1] == 'None':
                scenario = 'Non_Binding'
            else:
                scenario = rows_list[row][3 - 1]

            scenario_zone_period_targets[scenario] = OrderedDict()
            scenario_zone_period_targets[scenario]['CAISO'] = OrderedDict()

            for column in range(5 - 1, 40):
                period = int(float(rows_list[57 - 1][column]))
                scenario_zone_period_targets[scenario]['CAISO'][period] = dict()
                for subproblem in [1]:
                    scenario_zone_period_targets[scenario]['CAISO'][
                        period][subproblem] = dict()
                    for stage in [1]:
                        if scenario == 'Non_Binding':
                            scenario_zone_period_targets[scenario]['CAISO'][
                                period][subproblem][stage] \
                                = 99
                        else:
                            # Apply specified imports adder and convert to MMT
                            target = (float(rows_list[row][column])
                                      + float(rows_list[69 - 1][column])) \
                                     * 10**-6

                            scenario_zone_period_targets[scenario]['CAISO'][
                                period][subproblem][stage] = \
                                target

    # From carbon_caps.xlsx
    additional_targets = {
        "40_MMT": {
            "CAISO": {
                2015: 64.2,
                2016: 62.3,
                2017: 60.3,
                2018: 58.4,
                2019: 56.5,
                2020: 54.5,
                2021: 52.6,
                2022: 50.7,
                2023: 48.7,
                2024: 46.8,
                2025: 44.9,
                2026: 42.9,
                2027: 41.0,
                2028: 39.1,
                2029: 37.1,
                2030: 35.2
            }
        },
        "38_MMT": {
            "CAISO": {
                2015: 64.2,
                2016: 62.2,
                2017: 60.1,
                2018: 58.1,
                2019: 56.0,
                2020: 54.0,
                2021: 52.0,
                2022: 49.9,
                2023: 47.9,
                2024: 45.8,
                2025: 43.8,
                2026: 41.8,
                2027: 39.7,
                2028: 37.7,
                2029: 35.6,
                2030: 33.6
            }
        },
        "36_MMT": {
            "CAISO": {
                2015: 64.2,
                2016: 62.1,
                2017: 59.9,
                2018: 57.8,
                2019: 55.6,
                2020: 53.5,
                2021: 51.3,
                2022: 49.2,
                2023: 47.0,
                2024: 44.9,
                2025: 42.7,
                2026: 40.6,
                2027: 38.4,
                2028: 36.3,
                2029: 34.1,
                2030: 32.0
            }
        },
        "34_MMT": {
            "CAISO": {
                2015: 64.2,
                2016: 62.0,
                2017: 59.7,
                2018: 57.4,
                2019: 55.2,
                2020: 52.9,
                2021: 50.7,
                2022: 48.4,
                2023: 46.1,
                2024: 43.9,
                2025: 41.6,
                2026: 39.4,
                2027: 37.1,
                2028: 34.9,
                2029: 32.6,
                2030: 30.3
            }
        },
        "32_MMT": {
            "CAISO": {
                2015: 64.2,
                2016: 61.8,
                2017: 59.5,
                2018: 57.1,
                2019: 54.7,
                2020: 52.4,
                2021: 50.0,
                2022: 47.7,
                2023: 45.3,
                2024: 42.9,
                2025: 40.6,
                2026: 38.2,
                2027: 35.8,
                2028: 33.5,
                2029: 31.1,
                2030: 28.7
            }
        },
    }

    for scenario in additional_targets.keys():
        scenario_zone_period_targets[scenario] = OrderedDict()
        for zone in additional_targets[scenario].keys():
            scenario_zone_period_targets[scenario][zone] = OrderedDict()
            for period in additional_targets[scenario][zone].keys():
                scenario_zone_period_targets[scenario][zone][period] = dict()
                for subproblem in [1]:
                    scenario_zone_period_targets[scenario][zone][period][
                        subproblem] = dict()
                    for stage in [1]:
                        scenario_zone_period_targets[scenario][zone][
                            period][subproblem][stage] = \
                            additional_targets[scenario][zone][period]


    # Insert data
    scenario_id = 1
    for scenario in scenario_zone_period_targets.keys():
        carbon_cap.insert_carbon_cap_targets(
            io=io, c=c2,
            carbon_cap_target_scenario_id=scenario_id,
            scenario_name=scenario,
            scenario_description=scenario,
            zone_period_targets=scenario_zone_period_targets[scenario]
        )
        scenario_id += 1


def load_loads():
    """

    :return:
    """

    # ### CAISO ### #
    # Variables for totals to which profiles will be applied:
    # 1. Baseline load: equals the CEC IEPR total 'baseline' forecast minus
    # non-PV self gen
    # 2. EVs
    # 3. Building electrification
    # 4. Energy efficiency
    # 5. TOU adjustment: this is actually a percentage, which equals 0 before
    # 2018 and the ratio of:
    # (current year total baseline - current year Mid AAEE) to
    # (2025 total baseline - 2025 total Mid AAEE)

    # Variables for profiles
    # 1. Baseline load
    # 2. EVs: two profiles EVHome and EVHomeAndWork -- must get relative
    # penetration of EVHomeAndWork
    # 3. Building electrification
    # 4. Energy efficiency: two profiles, 2015 and 2030, weighted depending
    # on current year
    # 5. TOU: TOU adjustment for the year from above times profile for MW
    # impact

    # Forecasts
    with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_Forecast.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # Baseline forecast -- single scenario
        # Baseline - Non-PV Self Generation (CEC 2016 IEPR-Mid Demand)
        baseline_forecast = OrderedDict()
        for column in range(10 - 1, 45):
            yr = int(float(rows_list[12 - 1][column]))
            baseline_forecast[yr] = \
                (float(rows_list[13 - 1][column])
                 - float(rows_list[37 - 1][column]))

        # EV forecast -- three scenarios
        scenario_ev_forecast = OrderedDict()
        for row in range(17 - 1, 19):
            ev_scenario = rows_list[row][3 - 1]
            scenario_ev_forecast[ev_scenario] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[16 - 1][column]))
                scenario_ev_forecast[ev_scenario][yr] = \
                    float(rows_list[row][column])

        # Building electrification -- three scenarios
        scenario_building_electrification_forecast = OrderedDict()
        for row in range(23 - 1, 25):
            building_electrification_scenario = rows_list[row][3 - 1]
            scenario_building_electrification_forecast[
                building_electrification_scenario
            ] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[22 - 1][column]))
                scenario_building_electrification_forecast[
                    building_electrification_scenario
                ][yr] = \
                    float(rows_list[row][column])

        # Energy efficiency -- four scenarios
        scenario_ee_forecast = OrderedDict()
        for row in range(41 - 1, 45):
            ee_scenario = rows_list[row][3 - 1]
            scenario_ee_forecast[ee_scenario] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[40 - 1][column]))
                scenario_ee_forecast[ee_scenario][yr] = \
                    float(rows_list[row][column])

        # TOU adjustment percentage -- single scenario
        tou_adjustment_forecast = OrderedDict()
        for column in range(10 - 1, 45):
            yr = int(float(rows_list[48 - 1][column]))
            tou_adjustment_forecast[yr] = \
                float(rows_list[49 - 1][column])

    # Profiles
    with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_Profiles.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # Baseline profile
        baseline_profile = OrderedDict()
        for row in range(28 - 1, 915):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in baseline_profile.keys():
                baseline_profile[day][hour] = float(rows_list[row][7 - 1])
            else:
                baseline_profile[day] = OrderedDict()
                baseline_profile[day][hour] = float(rows_list[row][7 - 1])

        # EV home-only profile
        ev_home_only_profile = OrderedDict()
        for row in range(28 - 1, 915):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in ev_home_only_profile.keys():
                ev_home_only_profile[day][hour] = float(
                    rows_list[row][8 - 1])
            else:
                ev_home_only_profile[day] = OrderedDict()
                ev_home_only_profile[day][hour] = float(
                    rows_list[row][8 - 1])

        # EV home-only profile
        ev_home_and_work_profile = OrderedDict()
        for row in range(28 - 1, 915):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in ev_home_and_work_profile.keys():
                ev_home_and_work_profile[day][hour] = float(
                    rows_list[row][9 - 1])
            else:
                ev_home_and_work_profile[day] = OrderedDict()
                ev_home_and_work_profile[day][hour] = float(
                    rows_list[row][9 - 1])

        # Building electrification profile
        building_electrification_profile = OrderedDict()
        for row in range(28 - 1, 915):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in building_electrification_profile.keys():
                building_electrification_profile[day][hour] = float(
                    rows_list[row][10 - 1])
            else:
                building_electrification_profile[day] = OrderedDict()
                building_electrification_profile[day][hour] = float(
                    rows_list[row][10 - 1])

        # EE2015 profile
        energy_efficiency_2015_profile = OrderedDict()
        for row in range(28 - 1, 915):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in energy_efficiency_2015_profile.keys():
                energy_efficiency_2015_profile[day][hour] = float(
                    rows_list[row][11 - 1])
            else:
                energy_efficiency_2015_profile[day] = OrderedDict()
                energy_efficiency_2015_profile[day][hour] = float(
                    rows_list[row][11 - 1])

        # EE2030 profile
        energy_efficiency_2030_profile = OrderedDict()
        for row in range(28 - 1, 915):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in energy_efficiency_2030_profile.keys():
                energy_efficiency_2030_profile[day][hour] = float(
                    rows_list[row][12 - 1])
            else:
                energy_efficiency_2030_profile[day] = OrderedDict()
                energy_efficiency_2030_profile[day][hour] = float(
                    rows_list[row][12 - 1])

        # TOU MW impact profile -- 4 scenarios
        scenario_tou_mw_impact_profile = OrderedDict()
        for column in range(23 - 1, 26):
            tou_mw_impact_scenario = rows_list[27 - 1][column]
            scenario_tou_mw_impact_profile[tou_mw_impact_scenario] = \
                OrderedDict()
            for row in range(28 - 1, 915):
                day = int(float(rows_list[row][20 - 1]))
                hour = int(float(rows_list[row][21 - 1]))
                if day in scenario_tou_mw_impact_profile[
                    tou_mw_impact_scenario
                ].keys():
                    scenario_tou_mw_impact_profile[
                        tou_mw_impact_scenario
                    ][day][hour] = float(rows_list[row][column])
                else:
                    scenario_tou_mw_impact_profile[
                        tou_mw_impact_scenario
                    ][day] = OrderedDict()
                    scenario_tou_mw_impact_profile[
                        tou_mw_impact_scenario
                    ][day][hour] = float(rows_list[row][column])

    # EV home-and-work percent forecast
    with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_EV_Char.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        scenario_ev_home_and_work_charging_percent_forecast = OrderedDict()

        for row in range(18 - 1, 20):
            ev_home_and_work_charging_percent_forecast_scenario = \
                rows_list[row][3 - 1]
            scenario_ev_home_and_work_charging_percent_forecast[
                ev_home_and_work_charging_percent_forecast_scenario
            ] = OrderedDict()
            for column in range(4 - 1, 39):
                yr = int(float(rows_list[17 - 1][column]))
                scenario_ev_home_and_work_charging_percent_forecast[
                    ev_home_and_work_charging_percent_forecast_scenario
                ][yr] = float(rows_list[row][column])

    # Create the 'Default' CAISO load
    caiso_default_load = create_caiso_load_shape(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
            'CEC 2016 IEPR - Mid Demand'
        ],
        energy_efficiency_forecast=scenario_ee_forecast[
            'CEC 2016 IEPR - Mid AAEE + AB802'
        ],
        tou_adjustment_forecast=tou_adjustment_forecast,
        baseline_profile=baseline_profile,
        ev_home_and_work_charging_percent_forecast=
        scenario_ev_home_and_work_charging_percent_forecast['Mid'],
        ev_home_only_profile=ev_home_only_profile,
        ev_home_and_work_profile=ev_home_and_work_profile,
        building_electrification_profile=building_electrification_profile,
        energy_efficiency_2015_profile=energy_efficiency_2015_profile,
        energy_efficiency_2030_profile=energy_efficiency_2030_profile,
        tou_mw_impact_profile=scenario_tou_mw_impact_profile[
            'High (MRW S4 x1.5)'
        ],
        td_losses=0.0733
    )

    caiso_2x_aaee_load = create_caiso_load_shape(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
            'CEC 2016 IEPR - Mid Demand'
        ],
        energy_efficiency_forecast=scenario_ee_forecast[
            'CEC 2016 IEPR - Mid AAEE (x2)'
        ],
        tou_adjustment_forecast=tou_adjustment_forecast,
        baseline_profile=baseline_profile,
        ev_home_and_work_charging_percent_forecast=
        scenario_ev_home_and_work_charging_percent_forecast['Mid'],
        ev_home_only_profile=ev_home_only_profile,
        ev_home_and_work_profile=ev_home_and_work_profile,
        building_electrification_profile=building_electrification_profile,
        energy_efficiency_2015_profile=energy_efficiency_2015_profile,
        energy_efficiency_2030_profile=energy_efficiency_2030_profile,
        tou_mw_impact_profile=scenario_tou_mw_impact_profile[
            'High (MRW S4 x1.5)'
        ],
        td_losses=0.0733
    )

    caiso_mid_aaee_load = create_caiso_load_shape(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
            'CEC 2016 IEPR - Mid Demand'
        ],
        energy_efficiency_forecast=scenario_ee_forecast[
            'CEC 2016 IEPR - Mid AAEE'
        ],
        tou_adjustment_forecast=tou_adjustment_forecast,
        baseline_profile=baseline_profile,
        ev_home_and_work_charging_percent_forecast=
        scenario_ev_home_and_work_charging_percent_forecast['Mid'],
        ev_home_only_profile=ev_home_only_profile,
        ev_home_and_work_profile=ev_home_and_work_profile,
        building_electrification_profile=building_electrification_profile,
        energy_efficiency_2015_profile=energy_efficiency_2015_profile,
        energy_efficiency_2030_profile=energy_efficiency_2030_profile,
        tou_mw_impact_profile=scenario_tou_mw_impact_profile[
            'High (MRW S4 x1.5)'
        ],
        td_losses=0.0733
    )

    # Other regions -- single scenario
    with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_Forecast.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        region_year_forecast = OrderedDict()
        for row in range(148 - 1, 152):
            region = rows_list[row][3 - 1]
            region_year_forecast[region] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[146 -1][column]))
                region_year_forecast[region][yr] = \
                    float(rows_list[row][column])

    with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_Profiles.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # Region profile
        region_profile = OrderedDict()
        for column in range(13 - 1, 17):
            region = rows_list[27 - 1][column]
            region_profile[region] = OrderedDict()
            for row in range(28 - 1, 915):
                day = int(float(rows_list[row][3 - 1]))
                hour = int(float(rows_list[row][4 - 1]))
                if day in region_profile[region].keys():
                    region_profile[region][day][hour] = \
                        float(rows_list[row][column])
                else:
                    region_profile[region][day] = OrderedDict()
                    region_profile[region][day][hour] = \
                        float(rows_list[row][column])

    region_loads = OrderedDict()

    for region in region_year_forecast.keys():
        region_loads[region] = OrderedDict()
        region_loads[region][1] = OrderedDict()  # add stage
        for year in range(2015, 2050+1):
            for day in range(1, 37 + 1):
                for hour in range(1, 24 + 1):
                    tmp = year * 10**4 + day * 10**2 + hour
                    # Calculate the load
                    region_loads[region][1][tmp] = \
                        region_year_forecast[region][year] \
                        * region_profile[region][day][hour] * 10**3

    # Create final dictionary with the loads for all regions
    all_regions_loads = OrderedDict()
    all_regions_loads['CAISO'] = OrderedDict()
    all_regions_loads['CAISO'][1] = caiso_default_load
    for region in region_loads.keys():
        all_regions_loads[region] = OrderedDict()
        all_regions_loads[region][1] = OrderedDict()  # add stage
        all_regions_loads[region][1] = region_loads[region][1]

    # Insert data
    system_load.insert_system_static_loads(
        io=io, c=c2,
        load_scenario_id=1,
        scenario_name='IRP_Default_CAISO_Load_Mid_AAEE_plus_AB802',
        scenario_description='Default CAISO load and other regions loads '
                             'from the IRP',
        zone_stage_timepoint_static_loads=all_regions_loads
    )

    # 2x AAEE
    all_regions_loads['CAISO'][1] = caiso_2x_aaee_load
    # Insert data
    system_load.insert_system_static_loads(
        io=io, c=c2,
        load_scenario_id=2,
        scenario_name='IRP_CAISO_2x_AAEE_Load',
        scenario_description='Default CAISO load and other regions loads '
                             'from the IRP',
        zone_stage_timepoint_static_loads=all_regions_loads
    )

    # Mid AAEE only
    all_regions_loads['CAISO'][1] = caiso_mid_aaee_load
    # Insert data
    system_load.insert_system_static_loads(
        io=io, c=c2,
        load_scenario_id=3,
        scenario_name='IRP_CAISO_Mid_AAEE_Load',
        scenario_description='Default CAISO load and other regions loads '
                             'from the IRP',
        zone_stage_timepoint_static_loads=all_regions_loads
    )


def create_caiso_load_shape(
        baseline_forecast,
        ev_forecast,
        building_electrification_forecast,
        energy_efficiency_forecast,
        tou_adjustment_forecast,
        baseline_profile,
        ev_home_and_work_charging_percent_forecast,
        ev_home_only_profile,
        ev_home_and_work_profile,
        building_electrification_profile,
        energy_efficiency_2015_profile,
        energy_efficiency_2030_profile,
        tou_mw_impact_profile,
        td_losses
):
    """
    Forecasts are before td_losses and in GWh (except TOU adjustment)
    :param baseline_forecast:
    Baseline load forecast by year
    :param ev_forecast:
    EV load forecast by year
    :param building_electrification_forecast:
    Building electrification load forecast by year
    :param energy_efficiency_forecast:
    EE load forecast by year
    :param tou_adjustment_forecast:
    Percent adjustment to MW TOU impact by year
    :param baseline_profile:
    37 days x 24 hours normalized profile for baseline load
    :param ev_home_and_work_charging_percent_forecast:
    Percent of total EV load with access to workplace charging by year
    :param ev_home_only_profile:
    37 days x 24 hours normalized profile for EV load with home-only charging
    :param ev_home_and_work_profile:
    37 days x 24 hours normalized profile for EV load with home-and-work charging
    :param building_electrification_profile:
    37 days x 24 hours normalized profile for EV building electrification load
    :param energy_efficiency_2015_profile:
    37 days x 24 hours normalized profile for 2015 EE load
    :param energy_efficiency_2030_profile:
    37 days x 24 hours normalized profile for 2030 EE load
    :param tou_mw_impact_profile:
    37 days x 24 hours MW profile for TOU impact
    :return:
    """

    caiso_load = OrderedDict()

    for year in range(2015, 2050+1):
        for day in range(1, 37 + 1):
            for hour in range(1, 24 + 1):
                tmp = year * 10**4 + day * 10**2 + hour

                # Weights for two EE profiles
                weight_2030 = \
                    ((float(year) - 2015.0) / (2030.0 - 2015.0)) \
                    if year < 2030 else 1.0
                weight_2015 = \
                    1.0 - weight_2030 if year < 2030 else 0

                # Calculate the load
                caiso_load[tmp] = (
                    # Baseline
                    baseline_forecast[year] * baseline_profile[day][hour]
                    / (1 - td_losses)
                    # EV Home Only
                    + ev_forecast[year]
                    * (1 - ev_home_and_work_charging_percent_forecast[year])
                    * ev_home_only_profile[day][hour] / (1 - td_losses)
                    # EV Home and Work
                    + ev_forecast[year]
                    * ev_home_and_work_charging_percent_forecast[year]
                    * ev_home_and_work_profile[day][hour] / (1 - td_losses)
                    # EE 2015
                    - energy_efficiency_forecast[year]
                    * weight_2015
                    * energy_efficiency_2015_profile[day][hour]
                    / (1 - td_losses)
                    # EE 2030
                    - energy_efficiency_forecast[year]
                    * weight_2030
                    * energy_efficiency_2030_profile[day][hour]
                    / (1 - td_losses)
                    # Building electrification
                    + building_electrification_forecast[year]
                    * building_electrification_profile[day][hour]
                    / (1 - td_losses)
                ) * 10**3 \
                    + tou_adjustment_forecast[year] \
                    * tou_mw_impact_profile[day][hour]

    return caiso_load


def load_rps_targets():
    """

    :return:
    """
    # Forecasts
    with open(os.path.join("cpuc_irp_data", "csvs", "LOADS_Forecast.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # Baseline forecast -- single scenario
        # Note this is different from the baseline forecast in the load
        # generation method where non-PV self gen is subtracted
        baseline_forecast = OrderedDict()
        for column in range(10 - 1, 45):
            yr = int(float(rows_list[12 - 1][column]))
            baseline_forecast[yr] = \
                float(rows_list[13 - 1][column]) * 10**3

        # EV forecast -- three scenarios
        scenario_ev_forecast = OrderedDict()
        for row in range(17 - 1, 19):
            ev_scenario = rows_list[row][3 - 1]
            scenario_ev_forecast[ev_scenario] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[16 - 1][column]))
                scenario_ev_forecast[ev_scenario][yr] = \
                    float(rows_list[row][column]) * 10**3

        # Building electrification -- three scenarios
        scenario_building_electrification_forecast = OrderedDict()
        for row in range(23 - 1, 25):
            building_electrification_scenario = rows_list[row][3 - 1]
            scenario_building_electrification_forecast[
                building_electrification_scenario
            ] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[22 - 1][column]))
                scenario_building_electrification_forecast[
                    building_electrification_scenario
                ][yr] = \
                    float(rows_list[row][column]) * 10**3

        # Energy efficiency -- four scenarios
        scenario_ee_forecast = OrderedDict()
        for row in range(41 - 1, 45):
            ee_scenario = rows_list[row][3 - 1]
            scenario_ee_forecast[ee_scenario] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[40 - 1][column]))
                scenario_ee_forecast[ee_scenario][yr] = \
                    float(rows_list[row][column]) * 10**3

        # TOU effects -- three scenarios
        scenario_tou_effects_forecast = OrderedDict()
        for row in range(52 - 1, 55):
            tou_effects_scenario = rows_list[row][3 - 1]
            scenario_tou_effects_forecast[tou_effects_scenario] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[51 - 1][column]))
                scenario_tou_effects_forecast[tou_effects_scenario][yr] = \
                    float(rows_list[row][column]) * 10**3

        # Customer PV -- four scenarios
        scenario_customer_pv_forecast = OrderedDict()
        for row in range(29 - 1, 32):
            customer_pv_scenario = rows_list[row][3 - 1]
            scenario_customer_pv_forecast[customer_pv_scenario] = OrderedDict()
            for column in range(10 - 1, 45):
                yr = int(float(rows_list[28 - 1][column]))
                scenario_customer_pv_forecast[customer_pv_scenario][yr] = \
                    float(rows_list[row][column]) * 10**3

        # Pumping loads -- single scenario
        pumping_loads_forecast = OrderedDict()
        for column in range(10 - 1, 45):
            yr = int(float(rows_list[70 - 1][column]))
            pumping_loads_forecast[yr] = \
                float(rows_list[71 - 1][column]) * 10**3

        # Non-PV self generation -- single scenario
        self_gen_forecast = OrderedDict()
        for column in range(10 - 1, 45):
            yr = int(float(rows_list[36 - 1][column]))
            self_gen_forecast[yr] = \
                float(rows_list[37 - 1][column]) * 10**3

    # RPS targets
    scenario_rps_targets = OrderedDict()
    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_RPS_GHG_Targets.csv"), "r") as f:
        rows_list = list(csv.reader(f))
        for row in range(11 - 1, 16):
            rps_target_scenario = rows_list[row][3 - 1]
            scenario_rps_targets[rps_target_scenario] = OrderedDict()
            for column in range(5 - 1, 40):
                yr = int(float(rows_list[10 - 1][column]))
                scenario_rps_targets[rps_target_scenario][yr] = \
                    float(rows_list[row][column])

    # Adjustments
    # Non-modeled resources
    nonmodeled_resources = {
        yr: 2655000 for yr in range(2015, 2050 + 1)
    }

    # Banked resources
    # Only know this for 2026 and 2030; assign 0 to other years
    rps_bank = dict()
    for yr in range(2015, 2050 + 1):
        if yr == 2026:
            rps_bank[yr] = 145768.458533837
        elif yr == 2030:
            rps_bank[yr] = 8441051.90351601
        else:
            rps_bank[yr] = 0

    # Create RPS targets with default assumptions
    # Non-binding
    # rps_0 = {'CAISO': create_caiso_rps_target(
    #     baseline_forecast=baseline_forecast,
    #     ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
    #     building_electrification_forecast
    #     =scenario_building_electrification_forecast[
    #             'CEC 2016 IEPR - Mid Demand'
    #         ],
    #     energy_efficiency_forecast=scenario_ee_forecast[
    #             'CEC 2016 IEPR - Mid AAEE + AB802'
    #         ],
    #     customer_pv_forecast=scenario_customer_pv_forecast[
    #         'CEC 2016 IEPR - Mid PV'
    #     ],
    #     other_self_gen_forecast=self_gen_forecast,
    #     tou_adjustment_forecast=scenario_tou_effects_forecast[
    #             'High (MRW S4 x1.5)'],
    #     pumping_loads_forecast=pumping_loads_forecast,
    #     rps_percent_target=scenario_rps_targets['no RPS']
    # )}
    # Make it non-zero to avoid division by 0 issues
    rps_0 = {
        'CAISO': {yr: {1: {1: 100.0}} for yr in range(2015, 2050 + 1)}
    }

    rps_50 = {'CAISO': create_caiso_rps_target(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
                'CEC 2016 IEPR - Mid Demand'
            ],
        energy_efficiency_forecast=scenario_ee_forecast[
                'CEC 2016 IEPR - Mid AAEE + AB802'
            ],
        customer_pv_forecast=scenario_customer_pv_forecast[
            'CEC 2016 IEPR - Mid PV'
        ],
        other_self_gen_forecast=self_gen_forecast,
        tou_adjustment_forecast=scenario_tou_effects_forecast[
                'High (MRW S4 x1.5)'],
        pumping_loads_forecast=pumping_loads_forecast,
        rps_percent_target=scenario_rps_targets['50% by 2030'],
        nonmodeled_resources=nonmodeled_resources,
        rps_bank=rps_bank
    )}

    rps_55 = {'CAISO': create_caiso_rps_target(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
            'CEC 2016 IEPR - Mid Demand'
        ],
        energy_efficiency_forecast=scenario_ee_forecast[
            'CEC 2016 IEPR - Mid AAEE + AB802'
        ],
        customer_pv_forecast=scenario_customer_pv_forecast[
            'CEC 2016 IEPR - Mid PV'
        ],
        other_self_gen_forecast=self_gen_forecast,
        tou_adjustment_forecast=scenario_tou_effects_forecast[
            'High (MRW S4 x1.5)'],
        pumping_loads_forecast=pumping_loads_forecast,
        rps_percent_target=scenario_rps_targets['55% by 2030'],
        nonmodeled_resources=nonmodeled_resources,
        rps_bank=rps_bank
    )}

    rps_60 = {'CAISO': create_caiso_rps_target(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
            'CEC 2016 IEPR - Mid Demand'
        ],
        energy_efficiency_forecast=scenario_ee_forecast[
            'CEC 2016 IEPR - Mid AAEE + AB802'
        ],
        customer_pv_forecast=scenario_customer_pv_forecast[
            'CEC 2016 IEPR - Mid PV'
        ],
        other_self_gen_forecast=self_gen_forecast,
        tou_adjustment_forecast=scenario_tou_effects_forecast[
            'High (MRW S4 x1.5)'],
        pumping_loads_forecast=pumping_loads_forecast,
        rps_percent_target=scenario_rps_targets['60% by 2030'],
        nonmodeled_resources=nonmodeled_resources,
        rps_bank=rps_bank
    )}

    rps_80 = {'CAISO': create_caiso_rps_target(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
            'CEC 2016 IEPR - Mid Demand'
        ],
        energy_efficiency_forecast=scenario_ee_forecast[
            'CEC 2016 IEPR - Mid AAEE + AB802'
        ],
        customer_pv_forecast=scenario_customer_pv_forecast[
            'CEC 2016 IEPR - Mid PV'
        ],
        other_self_gen_forecast=self_gen_forecast,
        tou_adjustment_forecast=scenario_tou_effects_forecast[
            'High (MRW S4 x1.5)'],
        pumping_loads_forecast=pumping_loads_forecast,
        rps_percent_target=scenario_rps_targets['80% by 2030'],
        nonmodeled_resources=nonmodeled_resources,
        rps_bank=rps_bank
    )}

    rps_100 = {'CAISO': create_caiso_rps_target(
        baseline_forecast=baseline_forecast,
        ev_forecast=scenario_ev_forecast['CARB Scoping Plan - SP'],
        building_electrification_forecast
        =scenario_building_electrification_forecast[
            'CEC 2016 IEPR - Mid Demand'
        ],
        energy_efficiency_forecast=scenario_ee_forecast[
            'CEC 2016 IEPR - Mid AAEE + AB802'
        ],
        customer_pv_forecast=scenario_customer_pv_forecast[
            'CEC 2016 IEPR - Mid PV'
        ],
        other_self_gen_forecast=self_gen_forecast,
        tou_adjustment_forecast=scenario_tou_effects_forecast[
            'High (MRW S4 x1.5)'],
        pumping_loads_forecast=pumping_loads_forecast,
        rps_percent_target=scenario_rps_targets['100% by 2030'],
        nonmodeled_resources=nonmodeled_resources,
        rps_bank=rps_bank
    )}

    # Insert into DB
    rps.insert_rps_targets(
        io=io, c=c2,
        rps_target_scenario_id=1,
        scenario_name='default_non_binding',
        scenario_description='2030 RPS target of 0% at default load '
                             'assumptions.',
        zone_period_targets=rps_0
    )

    rps.insert_rps_targets(
        io=io, c=c2,
        rps_target_scenario_id=2,
        scenario_name='default_50_percent',
        scenario_description='2030 RPS target of 50% at default load '
                             'assumptions.',
        zone_period_targets=rps_50
    )

    rps.insert_rps_targets(
        io=io, c=c2,
        rps_target_scenario_id=3,
        scenario_name='default_55_percent',
        scenario_description='2030 RPS target of 55% at default load '
                             'assumptions.',
        zone_period_targets=rps_55
    )

    rps.insert_rps_targets(
        io=io, c=c2,
        rps_target_scenario_id=4,
        scenario_name='default_60_percent',
        scenario_description='2030 RPS target of 60% at default load '
                             'assumptions.',
        zone_period_targets=rps_60
    )

    rps.insert_rps_targets(
        io=io, c=c2,
        rps_target_scenario_id=5,
        scenario_name='default_80_percent',
        scenario_description='2030 RPS target of 80% at default load '
                             'assumptions.',
        zone_period_targets=rps_80
    )

    rps.insert_rps_targets(
        io=io, c=c2,
        rps_target_scenario_id=6,
        scenario_name='default_100_percent',
        scenario_description='2030 RPS target of 100% at default load '
                             'assumptions.',
        zone_period_targets=rps_100
    )


def create_caiso_rps_target(
        baseline_forecast,
        ev_forecast,
        building_electrification_forecast,
        energy_efficiency_forecast,
        customer_pv_forecast,
        other_self_gen_forecast,
        tou_adjustment_forecast,
        pumping_loads_forecast,
        rps_percent_target,
        nonmodeled_resources,
        rps_bank

):
    """
    MWh
    :param baseline_forecast:
    :param ev_forecast:
    :param building_electrification_forecast:
    :param energy_efficiency_forecast:
    :param customer_pv_forecast:
    :param other_self_gen_forecast:
    :param tou_adjustment_forecast:
    :param pumping_loads_forecast:
    :param rps_percent_target:
    :param nonmodeled_resources:
    :param rps_bank:
    :return:
    """

    rps_target = OrderedDict()
    for year in range(2015, 2050 + 1):
        rps_target[year] = OrderedDict()
        for subproblem in [1]:
            rps_target[year][subproblem] = OrderedDict()
            for stage in [1]:
                rps_target[year][subproblem][stage] = (
                    baseline_forecast[year]
                    + ev_forecast[year]
                    + building_electrification_forecast[year]
                    - energy_efficiency_forecast[year]
                    - customer_pv_forecast[year]
                    - other_self_gen_forecast[year]
                    + tou_adjustment_forecast[year]
                    - pumping_loads_forecast[year]
                ) * rps_percent_target[year] \
                    - nonmodeled_resources[year] \
                    - rps_bank[year]

    return rps_target


def load_reg_requirement():
    """
    1% of load
    :return:
    """
    percentage = 0.01
    load_scenario_id = 1
    reg_req_query = c2.execute(
        """SELECT timepoint, {} * load_mw
        FROM inputs_system_load
        WHERE load_zone = 'CAISO'
        AND stage_id = 1
        AND load_scenario_id = {};""".format(
            percentage, load_scenario_id
        )
    ).fetchall()

    reg_req = OrderedDict()
    reg_req['CAISO'] = OrderedDict()
    reg_req['CAISO'][1] = OrderedDict()

    for i in reg_req_query:
        tmp = i[0]
        req = i[1]
        reg_req['CAISO'][1][tmp] = req

    for direction in ["up", "down"]:
        system_reserves.insert_system_reserves(
            io=io, c=c2,
            reserve_scenario_id=load_scenario_id,
            scenario_name='1_pcnt_of_default_load_scenario_id_{}'.format(
                load_scenario_id
            ),
            scenario_description='One percent of the default IRP CAISO load',
            ba_stage_timepoint_reserve_req=reg_req,
            reserve_type="regulation_{}".format(direction)
        )


def load_spin_requirement():
    """
    3% of load
    :return:
    """
    percentage = 0.03
    load_scenario_id = 1
    spin_req_query = c2.execute(
        """SELECT timepoint, {} * load_mw
        FROM inputs_system_load
        WHERE load_zone = 'CAISO'
        AND stage_id = 1
        AND load_scenario_id = {};""".format(
            percentage, load_scenario_id
        )
    ).fetchall()

    spin_req = OrderedDict()
    spin_req['CAISO'] = OrderedDict()
    spin_req['CAISO'][1] = OrderedDict()

    for i in spin_req_query:
        tmp = i[0]
        req = i[1]
        spin_req['CAISO'][1][tmp] = req

    system_reserves.insert_system_reserves(
        io=io, c=c2,
        reserve_scenario_id=load_scenario_id,
        scenario_name='3_pcnt_of_default_load_scenario_id_{}'.format(
            load_scenario_id
        ),
        scenario_description='Three percent of the default IRP CAISO load',
        ba_stage_timepoint_reserve_req=spin_req,
        reserve_type="spinning_reserves"
    )


def load_lf_requirement():
    """
    Linearly interpolate between 2015 w/o Flex and the three 2030 scenarios
    for each year through 2050
    This is the approach taken in the IRP, although it may make more sense
    to interpolate between 2015 and 33% in 2020, then from 33% in 2020 to 2030
    NOTE: this probably overestimates reserve requirements after 2030
    :return:
    """

    with open(os.path.join("cpuc_irp_data", "csvs", "SYS_Reserves.csv"), "r") as f:
        rows_list = list(csv.reader(f))

        # LF up
        lf_up_2015_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_up_2015_req.keys():
                pass
            else:
                lf_up_2015_req[day] = OrderedDict()

            lf_up_2015_req[day][hour] = float(rows_list[row][5 - 1])

        # 33 percent
        lf_up_2030_33pcnt_req = OrderedDict()
        slope_intercept_33pcnt_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_up_2030_33pcnt_req.keys():
                pass
            else:
                lf_up_2030_33pcnt_req[day] = OrderedDict()
                slope_intercept_33pcnt_req[day] = OrderedDict()

            lf_up_2030_33pcnt_req[day][hour] = \
                float(rows_list[row][6 - 1])
            slope = \
                (lf_up_2030_33pcnt_req[day][hour]
                 - lf_up_2015_req[day][hour]) / (2030 - 2015)
            intercept = \
                lf_up_2015_req[day][hour] - 2015 * slope
            slope_intercept_33pcnt_req[day][hour] = (slope, intercept)

        # High solar
        lf_up_2030_high_solar_req = OrderedDict()
        slope_intercept_high_solar_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_up_2030_high_solar_req.keys():
                pass
            else:
                lf_up_2030_high_solar_req[day] = OrderedDict()
                slope_intercept_high_solar_req[day] = OrderedDict()

            lf_up_2030_high_solar_req[day][hour] = \
                float(rows_list[row][7 - 1])
            slope = \
                (lf_up_2030_high_solar_req[day][hour]
                 - lf_up_2015_req[day][hour]) / (2030 - 2015)
            intercept = \
                lf_up_2015_req[day][hour] - 2015 * slope
            slope_intercept_high_solar_req[day][hour] = (slope, intercept)

        # Diverse
        lf_up_2030_diverse_req = OrderedDict()
        slope_intercept_diverse_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_up_2030_diverse_req.keys():
                pass
            else:
                lf_up_2030_diverse_req[day] = OrderedDict()
                slope_intercept_diverse_req[day] = OrderedDict()

            lf_up_2030_diverse_req[day][hour] = \
                float(rows_list[row][8 - 1])
            slope = \
                (lf_up_2030_diverse_req[day][hour]
                 - lf_up_2015_req[day][hour]) / (2030 - 2015)
            intercept = \
                lf_up_2015_req[day][hour] - 2015 * slope
            slope_intercept_diverse_req[day][hour] = (slope, intercept)


        # Calculate and insert final requirements
        lf_reserves_up_33pcnt_req = OrderedDict()
        lf_reserves_up_33pcnt_req['CAISO'] = OrderedDict()
        lf_reserves_up_33pcnt_req['CAISO'][1] = OrderedDict()
        for yr in range(2015, 2050 + 1):
            for day in slope_intercept_33pcnt_req.keys():
                for hour in slope_intercept_33pcnt_req[day].keys():
                    tmp = yr * 10**4 + day * 10**2 + hour
                    lf_reserves_up_33pcnt_req['CAISO'][1][tmp] = \
                        slope_intercept_33pcnt_req[day][hour][0] * yr \
                        + slope_intercept_33pcnt_req[day][hour][1]

        lf_reserves_up_high_solar_req = OrderedDict()
        lf_reserves_up_high_solar_req['CAISO'] = OrderedDict()
        lf_reserves_up_high_solar_req['CAISO'][1] = OrderedDict()
        for yr in range(2015, 2050 + 1):
            for day in slope_intercept_high_solar_req.keys():
                for hour in slope_intercept_high_solar_req[day].keys():
                    tmp = yr * 10**4 + day * 10**2 + hour
                    lf_reserves_up_high_solar_req['CAISO'][1][tmp] = \
                        slope_intercept_high_solar_req[day][hour][0] * yr \
                        + slope_intercept_high_solar_req[day][hour][1]

        lf_reserves_up_diverse_req = OrderedDict()
        lf_reserves_up_diverse_req['CAISO'] = OrderedDict()
        lf_reserves_up_diverse_req['CAISO'][1] = OrderedDict()
        for yr in range(2015, 2050 + 1):
            for day in slope_intercept_diverse_req.keys():
                for hour in slope_intercept_diverse_req[day].keys():
                    tmp = yr * 10**4 + day * 10**2 + hour
                    lf_reserves_up_diverse_req['CAISO'][1][tmp] = \
                        slope_intercept_diverse_req[day][hour][0] * yr \
                        + slope_intercept_diverse_req[day][hour][1]

        # Insert data
        system_reserves.insert_system_reserves(
            io=io, c=c2,
            reserve_scenario_id=1,
            scenario_name='high_solar_lf_up',
            scenario_description='High solar portfolio LF reserves up',
            ba_stage_timepoint_reserve_req=lf_reserves_up_high_solar_req,
            reserve_type="lf_reserves_up"
        )

        system_reserves.insert_system_reserves(
            io=io, c=c2,
            reserve_scenario_id=2,
            scenario_name='diverse_lf_up',
            scenario_description='Diverse portfolio LF reserves up',
            ba_stage_timepoint_reserve_req=lf_reserves_up_diverse_req,
            reserve_type="lf_reserves_up"
        )

        system_reserves.insert_system_reserves(
            io=io, c=c2,
            reserve_scenario_id=3,
            scenario_name='33pcnt_lf_up',
            scenario_description='33 percent portfolio LF reserves up',
            ba_stage_timepoint_reserve_req=lf_reserves_up_33pcnt_req,
            reserve_type="lf_reserves_up"
        )

        # LF down
        lf_down_2015_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_down_2015_req.keys():
                pass
            else:
                lf_down_2015_req[day] = OrderedDict()

            lf_down_2015_req[day][hour] = float(rows_list[row][14 - 1])

        # 33 percent
        lf_down_2030_33pcnt_req = OrderedDict()
        slope_intercept_33pcnt_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_down_2030_33pcnt_req.keys():
                pass
            else:
                lf_down_2030_33pcnt_req[day] = OrderedDict()
                slope_intercept_33pcnt_req[day] = OrderedDict()

            lf_down_2030_33pcnt_req[day][hour] = \
                float(rows_list[row][15 - 1])
            slope = \
                (lf_down_2030_33pcnt_req[day][hour]
                 - lf_down_2015_req[day][hour]) / (2030 - 2015)
            intercept = \
                lf_down_2015_req[day][hour] - 2015 * slope
            slope_intercept_33pcnt_req[day][hour] = (slope, intercept)

        # High solar
        lf_down_2030_high_solar_req = OrderedDict()
        slope_intercept_high_solar_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_down_2030_high_solar_req.keys():
                pass
            else:
                lf_down_2030_high_solar_req[day] = OrderedDict()
                slope_intercept_high_solar_req[day] = OrderedDict()

            lf_down_2030_high_solar_req[day][hour] = \
                float(rows_list[row][16 - 1])
            slope = \
                (lf_down_2030_high_solar_req[day][hour]
                 - lf_down_2015_req[day][hour]) / (2030 - 2015)
            intercept = \
                lf_down_2015_req[day][hour] - 2015 * slope
            slope_intercept_high_solar_req[day][hour] = (slope, intercept)

        # Diverse
        lf_down_2030_diverse_req = OrderedDict()
        slope_intercept_diverse_req = OrderedDict()
        for row in range(37 - 1, 924):
            day = int(float(rows_list[row][3 - 1]))
            hour = int(float(rows_list[row][4 - 1]))
            if day in lf_down_2030_diverse_req.keys():
                pass
            else:
                lf_down_2030_diverse_req[day] = OrderedDict()
                slope_intercept_diverse_req[day] = OrderedDict()

            lf_down_2030_diverse_req[day][hour] = \
                float(rows_list[row][17 - 1])
            slope = \
                (lf_down_2030_diverse_req[day][hour]
                 - lf_down_2015_req[day][hour]) / (2030 - 2015)
            intercept = \
                lf_down_2015_req[day][hour] - 2015 * slope
            slope_intercept_diverse_req[day][hour] = (slope, intercept)

        # Calculate and insert final requirements
        lf_reserves_down_33pcnt_req = OrderedDict()
        lf_reserves_down_33pcnt_req['CAISO'] = OrderedDict()
        lf_reserves_down_33pcnt_req['CAISO'][1] = OrderedDict()
        for yr in range(2015, 2050 + 1):
            for day in slope_intercept_33pcnt_req.keys():
                for hour in slope_intercept_33pcnt_req[day].keys():
                    tmp = yr * 10 ** 4 + day * 10 ** 2 + hour
                    lf_reserves_down_33pcnt_req['CAISO'][1][tmp] = \
                        slope_intercept_33pcnt_req[day][hour][0] * yr \
                        + slope_intercept_33pcnt_req[day][hour][1]

        lf_reserves_down_high_solar_req = OrderedDict()
        lf_reserves_down_high_solar_req['CAISO'] = OrderedDict()
        lf_reserves_down_high_solar_req['CAISO'][1] = OrderedDict()
        for yr in range(2015, 2050 + 1):
            for day in slope_intercept_high_solar_req.keys():
                for hour in slope_intercept_high_solar_req[day].keys():
                    tmp = yr * 10 ** 4 + day * 10 ** 2 + hour
                    lf_reserves_down_high_solar_req['CAISO'][1][tmp] = \
                        slope_intercept_high_solar_req[day][hour][0] * yr \
                        + slope_intercept_high_solar_req[day][hour][1]

        lf_reserves_down_diverse_req = OrderedDict()
        lf_reserves_down_diverse_req['CAISO'] = OrderedDict()
        lf_reserves_down_diverse_req['CAISO'][1] = OrderedDict()
        for yr in range(2015, 2050 + 1):
            for day in slope_intercept_diverse_req.keys():
                for hour in slope_intercept_diverse_req[day].keys():
                    tmp = yr * 10 ** 4 + day * 10 ** 2 + hour
                    lf_reserves_down_diverse_req['CAISO'][1][tmp] = \
                        slope_intercept_diverse_req[day][hour][0] * yr \
                        + slope_intercept_diverse_req[day][hour][1]

        # Insert data
        system_reserves.insert_system_reserves(
            io=io, c=c2,
            reserve_scenario_id=1,
            scenario_name='high_solar_lf_down',
            scenario_description='High solar portfolio LF reserves down',
            ba_stage_timepoint_reserve_req=lf_reserves_down_high_solar_req,
            reserve_type="lf_reserves_down"
        )

        system_reserves.insert_system_reserves(
            io=io, c=c2,
            reserve_scenario_id=2,
            scenario_name='diverse_lf_down',
            scenario_description='Diverse portfolio LF reserves down',
            ba_stage_timepoint_reserve_req=lf_reserves_down_diverse_req,
            reserve_type="lf_reserves_down"
        )

        system_reserves.insert_system_reserves(
            io=io, c=c2,
            reserve_scenario_id=3,
            scenario_name='33pcnt_lf_down',
            scenario_description='33 percent portfolio LF reserves down',
            ba_stage_timepoint_reserve_req=lf_reserves_down_33pcnt_req,
            reserve_type="lf_reserves_down"
        )


def load_frequency_response_requirement():
    """

    :return:
    """
    print ("frequency response")

    freq_resp_scenario = 1
    freq_resp_scenario_name = 'default_caiso'
    description = 'Default from IRP draft data, total 770 MW, partial 385 MW ' \
                  'at all times.'
    # Subscenarios
    c2.execute(
        """INSERT INTO subscenarios_system_frequency_response
        (frequency_response_scenario_id,
        name, description)
        VALUES ({}, '{}', '{}');""".format(
            freq_resp_scenario,
            freq_resp_scenario_name,
            description
        )
    )
    io.commit()

    for yr in range(2015, 2050 + 1):
        for day in range(1, 37 + 1):
            for hour in range(1, 24 + 1):
                tmp = yr * 10 ** 4 + day * 10 ** 2 + hour

                c2.execute(
                    """INSERT INTO inputs_system_frequency_response
                    (frequency_response_scenario_id, frequency_response_ba, 
                    stage_id,
                    timepoint, frequency_response_mw, 
                    frequency_response_partial_mw)
                    VALUES ({}, '{}', {}, {}, {}, {});""".format(
                        freq_resp_scenario, 'CAISO', 1, tmp, 770, 385
                    )
                )
    io.commit()


def load_prm_requirement():
    """

    :return:
    """
    # Manually get the peak load requirement at the default load assumptions
    peak_load = {
        2018: 46404,
        2022: 45815,
        2026: 45568,
        2030: 45582
    }

    effective_import_capacity = 9891

    prm = 0.15

    zone_period_prm_requirement = {'CAISO': OrderedDict()}

    for period in peak_load.keys():
        zone_period_prm_requirement['CAISO'][period] = \
            (1 + prm) * peak_load[period] - effective_import_capacity

    system_prm.prm_requirement(
        io=io, c=c2,
        prm_requirement_scenario_id=1,
        scenario_name='default_load_15pcnt_margin',
        scenario_description='Default load level, 15 percent margin, adjust '
                             'for imports minus OOS gen capacity modeled in '
                             'CAISO',
        zone_period_requirement=zone_period_prm_requirement
)


def load_local_capacity_requirement():
    """

    :return:
    """
    print("local capacity requirement")

    # Subscenarios
    c2.execute(
        """INSERT INTO subscenarios_system_local_capacity_requirement
        (local_capacity_requirement_scenario_id, name, description)
        VALUES ({}, '{}', '{}');""".format(
            1, "default", "default"
        )
    )
    io.commit()

    zone_period_requirement = dict()
    with open(os.path.join("cpuc_irp_data", "gridpath_specific",
                           "lcr.csv"),
              "r") as f:
        reader = csv.reader(f, delimiter=",")
        next(reader)
        for row in reader:
            zone_period_requirement[row[0]] = dict()
            zone_period_requirement[row[0]][2018] = \
                0 if float(row[1]) < 0 else float(row[1])
            zone_period_requirement[row[0]][2022] = \
                0 if float(row[2]) < 0 else float(row[2])
            zone_period_requirement[row[0]][2026] = \
                0 if float(row[3]) < 0 else float(row[3])
            zone_period_requirement[row[0]][2030] = \
                0 if float(row[4]) < 0 else float(row[4])

    # Insert data
    for zone in zone_period_requirement.keys():
        for period in zone_period_requirement[zone].keys():
            c2.execute(
                """INSERT INTO inputs_system_local_capacity_requirement
                (local_capacity_requirement_scenario_id, 
                local_capacity_zone, period, local_capacity_requirement_mw)
                VALUES ({}, '{}', {}, {});""".format(
                    1, zone, period,
                    zone_period_requirement[zone][period]
                )
            )
    io.commit()


def tuning():
    """

    :return:
    """
    print("tuning")
    subscenarios = [
        (0, 'no_tuning', 'No tuning (tuning params = 0)', 0, 0, 0),
        (1, 'tune_carbon_imports', 'Tune carbon imports only', 10e-9, 0, 0),
        (2, 'tune_ramps', 'Tune ramps only', 0, 10e-9, 0),
        (3, 'tune_carbon_imports_and_ramps',
         'Tune carbon imports and ramps', 10e-9, 10e-9, 0),
        (4, 'tune_elcc_surface', 'Tune dynamic ELCC only', 0, 0, 10e-3),
        (5, 'tune_carbon_imports_and_elcc_surface',
         'Tune carbon imports and ELCC surface', 10e-9, 0, 10e-3),
        (6, 'tune_carbon_imports_ramps_and_elcc_surface',
         'Tune carbon imports, ramps, and ELCC surface', 10e-9, 10e-9, 10e-3)
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
            (tuning_scenario_id, import_carbon_tuning_cost_per_ton, 
            ramp_tuning_cost_per_mw, dynamic_elcc_tuning_cost_per_mw)
            VALUES ({}, {}, {}, {});""".format(
                subscenario[0], subscenario[3], subscenario[4], subscenario[5]
            )
        )
    io.commit()


def transmission_operational_chars():
    """
    Currently does nothing
    Single transmission portfolio from IRP  data
    No new transmission
    :return:
    """
    print("transmission operational chars")

    transmission_operational_chars_scenario_id = 1
    name = 'default_lines_does_nothing'
    description = 'Default transmission lines, but no operating chars'

    # Subscenarios
    c2.execute(
        """INSERT INTO subscenarios_transmission_operational_chars
        (transmission_operational_chars_scenario_id, name, description)
        VALUES ({}, '{}', '{}');""".format(
            transmission_operational_chars_scenario_id, name, description
        )
    )
    io.commit()

    tx_lines = c2.execute(
        """SELECT DISTINCT transmission_line
        FROM inputs_transmission_portfolios
        WHERE transmission_portfolio_scenario_id = 1;"""
    ).fetchall()

    for tx_line in tx_lines:
        c2.execute(
            """INSERT INTO inputs_transmission_operational_chars
               (transmission_operational_chars_scenario_id,
               transmission_line, operational_type)
               VALUES ({}, '{}', 'simple_transmission');""".format(
                transmission_operational_chars_scenario_id,
                tx_line[0]
            )
        )
    io.commit()


def options_solver():
    c2.execute(
        """INSERT INTO options_solver_descriptions
        (solver_options_id, solver, description)
        VALUES (1, 'cplex', 'cplex, barrier, 4 threads');"""
    )

    c2.execute(
        """INSERT INTO options_solver_values 
        (solver_options_id, solver, solver_option_name, solver_option_value) 
        VALUES (1, 'cplex', 'lpmethod', 4), (1, 'cplex', 'threads', '4')"""
    )
    io.commit()

# def hybridization_costs():
#     """
#
#     :return:
#     """
#     print("hybridization costs")
#     c2.execute("""DROP TABLE IF EXISTS inputs_project_hybridization_costs;""")
#     c2.execute(
#         """CREATE TABLE inputs_project_hybridization_costs (
#         project VARCHAR(64),
#         vintage INTEGER,
#         hybridization_annualized_real_cost_per_mw_yr FLOAT,
#         PRIMARY KEY (project, vintage)
#         );""")
#     io.commit()
#
#     with open(os.path.join("cpuc_irp_data", "gridpath_specific",
#                            "hybridization_costs.csv"),
#               "r") as f:
#         reader = csv.reader(f, delimiter=",")
#         next(reader)
#         for row in reader:
#             proj, yr, cost = row[0], row[1], row[2]
#             c2.execute(
#                 """INSERT INTO inputs_project_hybridization_costs
#                 (project, vintage,
#                 hybridization_annualized_real_cost_per_mw_yr)
#                 VALUES ('{}', {}, {});""".format(proj, yr, cost)
#             )
#     io.commit()


# Default subscenario IDs for
defaults = {
    "of_fuels": 1,
    "of_multi_stage": 0,
    "of_transmission": 1,
    "of_transmission_hurdle_rates": 1,
    "of_simultaneous_flow_limits": 1,
    "of_lf_reserves_up": 1,
    "of_lf_reserves_down": 1,
    "of_regulation_up": 1,
    "of_regulation_down": 1,
    "of_frequency_response": 1,
    "of_spinning_reserves": 1,
    "of_rps": 1,
    "of_carbon_cap": 1,
    "of_track_carbon_imports": 1,
    "of_prm": 1,
    "of_local_capacity": 1,
    "of_elcc_surface": 1,
    "of_tuning": 0,
    "temporal_scenario_id": 1,
    "load_zone_scenario_id": 1,
    "lf_reserves_up_ba_scenario_id": 1,
    "lf_reserves_down_ba_scenario_id": 1,
    "regulation_up_ba_scenario_id": 1,
    "regulation_down_ba_scenario_id": 1,
    "frequency_response_ba_scenario_id": 1,
    "spinning_reserves_ba_scenario_id": 1,
    "rps_zone_scenario_id": 1,
    "carbon_cap_zone_scenario_id": 1,
    "prm_zone_scenario_id": 1,
    "local_capacity_zone_scenario_id": 1,
    "project_portfolio_scenario_id": 11,  # No retirements, Existing Tx only
    "project_operational_chars_scenario_id": 1,
    "project_availability_scenario_id": 1,
    "fuel_scenario_id": 1,
    "project_load_zone_scenario_id": 1,
    "project_lf_reserves_up_ba_scenario_id": 1,
    "project_lf_reserves_down_ba_scenario_id": 2,  # Allow renewables
    "project_regulation_up_ba_scenario_id": 1,
    "project_regulation_down_ba_scenario_id": 1,
    "project_frequency_response_ba_scenario_id": 1,
    "project_spinning_reserves_ba_scenario_id": 1,
    "project_rps_zone_scenario_id": 1,
    "project_carbon_cap_zone_scenario_id": 1,
    "project_prm_zone_scenario_id": 1,
    "project_elcc_chars_scenario_id": 1,
    "prm_energy_only_scenario_id": 1,
    "project_local_capacity_zone_scenario_id": 1,
    "project_local_capacity_chars_scenario_id": 1,
    "project_existing_capacity_scenario_id": 7,  # Storage mandate, mid PV BTM
    "project_existing_fixed_cost_scenario_id": 1,
    "fuel_price_scenario_id": 5,  # 'Mid' fuel prices, 'Low' carbon adder
    "project_new_cost_scenario_id": 1,
    "project_new_potential_scenario_id": 7,  # DRECP/SJV w LCR proj limits
    "project_new_binary_build_size_scenario_id": None,
    "transmission_portfolio_scenario_id": 1,
    "transmission_load_zone_scenario_id": 1,
    "transmission_existing_capacity_scenario_id": 1,
    "transmission_operational_chars_scenario_id": 1,
    "transmission_hurdle_rate_scenario_id": 1,  # Base hurdle + 'Low' CO2
    "transmission_carbon_cap_zone_scenario_id": 1,
    "transmission_simultaneous_flow_limit_scenario_id": 2,  # "Mid"
    "transmission_simultaneous_flow_limit_line_group_scenario_id": 1,
    "load_scenario_id": 1,
    "lf_reserves_up_scenario_id": 1,
    "lf_reserves_down_scenario_id": 1,
    "regulation_up_scenario_id": 1,
    "regulation_down_scenario_id": 1,
    "frequency_response_scenario_id": 1,
    "spinning_reserves_scenario_id": 1,
    "rps_target_scenario_id": 2,  # 50% by 2030
    "carbon_cap_target_scenario_id": 5,
    "prm_requirement_scenario_id": 1,
    "elcc_surface_scenario_id": 1,
    "local_capacity_requirement_scenario_id": 1,
    "tuning_scenario_id": None,
    "solver_options_id": 1
}


def create_base_scenarios():
    """

    :return:
    """

    # Create 'Base_42MMT' scenario
    # Implement carbon cap
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT",
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
        of_tuning=defaults["of_tuning"],
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
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT' scenario
    # Implement carbon cap
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT",
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
        of_tuning=defaults["of_tuning"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        transmission_carbon_cap_zone_scenario_id=defaults[
            "carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Peakers and CCGTs allowed to retire
    # Base_42MMT_Ret
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_Ret",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=14,
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
        defaults["carbon_cap_zone_scenario_id"],
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
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT_Ret' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_Ret",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=14,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=
        defaults["project_local_capacity_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=defaults[
            "project_local_capacity_zone_scenario_id"],
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Base_42MMT_Ret_NoLCR
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_Ret_NoLCR",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=14,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_local_capacity_chars_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT_Ret_NoLCR' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_Ret_NoLCR",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=14,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=None,
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Base_42MMT_AggFleet
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_AggFleet",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=2,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_local_capacity_chars_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT_AggFleet' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_AggFleet",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=2,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=None,
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Base_42MMT_AggFleet_Ret
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_AggFleet_Ret",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=5,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_local_capacity_chars_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT_AggFleet_Ret' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_AggFleet_Ret",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=5,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=None,
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Peakers and CCGTs allowed to retire, LCR-eligible batteries & gas
    # included
    # Base_42MMT_Ret_LCR_Resources_4h_Batteries
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_Ret_LCR_Resources_4h_Batteries",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=20,
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
        defaults["carbon_cap_zone_scenario_id"],
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
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT_Ret_LCR_Resources_4h_Batteries' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_Ret_LCR_Resources_4h_Batteries",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=20,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=
        defaults["project_local_capacity_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=defaults[
            "project_local_capacity_zone_scenario_id"],
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # 6h batteries
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_Ret_LCR_Resources_6h_Batteries",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=20,
        project_operational_chars_scenario_id=2,
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
        defaults["carbon_cap_zone_scenario_id"],
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
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT_Ret_LCR_Resources_6h_Batteries' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_Ret_LCR_Resources_6h_Batteries",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=20,
        project_operational_chars_scenario_id=2,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=
        defaults["project_local_capacity_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=defaults[
            "project_local_capacity_zone_scenario_id"],
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # 8h batteries
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_Ret_LCR_Resources_8h_Batteries",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=20,
        project_operational_chars_scenario_id=3,
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
        defaults["carbon_cap_zone_scenario_id"],
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
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )

    # Create 'Base_30MMT_Ret_LCR_Resources_8h_Batteries' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_Ret_LCR_Resources_8h_Batteries",
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
        of_tuning=defaults["of_tuning"],
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
        project_portfolio_scenario_id=20,
        project_operational_chars_scenario_id=3,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=
        defaults["project_local_capacity_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=defaults[
            "project_local_capacity_zone_scenario_id"],
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=defaults[
            "local_capacity_requirement_scenario_id"],
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=defaults["solver_options_id"]
    )


def create_2030_scenario():
    # Base_42MMT_AggFleet_2030
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_42MMT_AggFleet_2030",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
        temporal_scenario_id=2,
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=2,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_local_capacity_chars_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=2,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=None
    )

    # Create 'Base_30MMT_AggFleet_2030' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_AggFleet_2030",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
        temporal_scenario_id=2,
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=2,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=None,
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=None
    )

    # Create 'Base_30MMT_AggFleet_2030_1horizon' scenario
    scenario.create_scenario_all_args(
        io=io, c=c2,
        scenario_name="Base_30MMT_AggFleet_2030_1horizon",
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
        of_local_capacity=0,
        of_elcc_surface=defaults["of_elcc_surface"],
        of_tuning=defaults["of_tuning"],
        temporal_scenario_id=3,
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
        local_capacity_zone_scenario_id=None,
        project_portfolio_scenario_id=2,
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
        defaults["carbon_cap_zone_scenario_id"],
        project_prm_zone_scenario_id=defaults["project_prm_zone_scenario_id"],
        project_elcc_chars_scenario_id=
        defaults["project_elcc_chars_scenario_id"],
        project_local_capacity_chars_scenario_id=None,
        prm_energy_only_scenario_id=
        defaults["prm_energy_only_scenario_id"],
        project_local_capacity_zone_scenario_id=None,
        project_existing_capacity_scenario_id=
        defaults["project_existing_capacity_scenario_id"],
        project_existing_fixed_cost_scenario_id=
        defaults["project_existing_fixed_cost_scenario_id"],
        fuel_price_scenario_id=defaults["fuel_price_scenario_id"],
        project_new_cost_scenario_id=defaults["project_new_cost_scenario_id"],
        project_new_potential_scenario_id=
        defaults["project_new_potential_scenario_id"],
        project_new_binary_build_size_scenario_id=
        defaults["project_new_binary_build_size_scenario_id"],
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
        defaults["carbon_cap_zone_scenario_id"],
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
        rps_target_scenario_id=defaults["rps_target_scenario_id"],
        carbon_cap_target_scenario_id=1,
        prm_requirement_scenario_id=defaults["prm_requirement_scenario_id"],
        elcc_surface_scenario_id=defaults["elcc_surface_scenario_id"],
        local_capacity_requirement_scenario_id=None,
        tuning_scenario_id=defaults["tuning_scenario_id"],
        solver_options_id=None
    )


if __name__ == "__main__":
    # ### Temporal ### #
    load_temporal_data()
    load_temporal_data_2030_only()
    load_temporal_data_2030_only_1horizon()

    # ### Geography ### #
    load_geography_load_zones()
    load_geography_lf_reserves_up_bas()
    load_geography_lf_reserves_down_bas()
    load_geography_regulation_up_bas()
    load_geography_regulation_down_bas()
    load_geography_spinning_reserves_bas()
    load_geography_frequency_response_bas()
    load_geography_rps_zones()
    load_geography_carbon_cap_zones()
    load_geography_prm_zones()
    load_geography_local_capacity_zones()

    # ### Projects ### #
    load_projects()  # project list

    load_project_load_zones()  # project load zones
    load_project_lf_reserves_bas()  # project lf reserves up and down bas
    load_project_regulation_bas()  # project regulation up and down bas
    load_project_spinning_reserves_bas()  # project spinning reserves bas
    load_project_frequency_response_bas()  # project frequency response bas
    load_project_rps_zones()  # project rps zones
    load_project_prm_zones()  # project prm zones
    load_project_local_capacity_zones()  # project local capacity zones
    load_project_carbon_cap_zones()  # project carbon cap zones

    load_project_operational_chars()  # operational characteristics
    load_project_variable_profiles()  # variable generator cap factors
    load_project_hydro_opchar()  # hydro operational characteristics
    load_project_availability()  # project availability (i.e. maintenance)

    load_fuels()  # fuels and fuel characteristics
    load_fuel_prices()  # fuel prices

    load_project_portfolios()  # project portfolios
    load_project_capacities()  # existing/planned project capacities
    load_project_fixed_costs()  # existing/planned project fixed costs

    load_project_new_costs()  # candidate project costs
    load_project_new_potentials()  # candidate project potentials

    load_project_elcc_chars()  # elcc characteristics
    load_deliverability_group_params()  # deliverability group parameters
    load_elcc_surface()  # elcc surface

    load_local_capacity_chars()  # local capacity contribution chars

    # ### Transmission ### #
    load_transmission_portfolios()  # transmission portfolios
    load_transmission_load_zones()  # transmission load zones
    load_transmission_capacities()  # transmission capacities
    transmission_operational_chars()

    load_simultaneous_flow_groups()  # simultaneous flow group lines
    load_simultaneous_flow_limits()

    load_transmission_carbon_cap_zones()
    load_transmission_hurdle_rates()

    # ### Policy ### #
    load_carbon_cap_targets()
    load_rps_targets()

    # ### Loads ### #
    load_loads()

    # # ### Reserves ### #
    load_reg_requirement()
    load_spin_requirement()
    load_lf_requirement()
    load_frequency_response_requirement()

    # ### PRM requirement ### #
    load_prm_requirement()

    # ### Local capacity requirement ### #
    load_local_capacity_requirement()

    # ### Tuning ### #
    tuning()

    # # ### Solver options ### # #
    options_solver()

    # ### Hybridization costs ### #
    # hybridization_costs()

    # ### Create scenarios ### #
    create_base_scenarios()
    create_2030_scenario()
    # create_2x_aaee_scenarios()
    # create_mid_aaee_scenarios()
    # create_high_export_scenarios()
    # create_low_export_scenarios()
    # create_oos_scenarios()
    # create_no_oos_scenarios()
    # create_low_battery_cost_scenarios()
    # create_high_battery_cost_scenarios()

    # Vacuum
    io.execute("VACUUM")
    io.close()
