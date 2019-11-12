#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

"""
Modules to import load data
"""
from collections import OrderedDict
import pandas as pd
from db.utilities import system_load


def create_loads():
    """

    :return:
    """
    # Import demand timeseries from csv
    global load_all
    load_all = pd.read_csv(load_csv)

    # Select demand for selected years and days of month (for last year 2030)
    year_month_days_selected = pd.read_csv(year_month_days_selected_csv)
    year_month_days_selected = year_month_days_selected.set_index(['year', 'month', 'day'])
    load_all = load_all.set_index(['year', 'month', 'day'])
    global load
    load = load_all[(load_all.index.isin(year_month_days_selected.index))].reset_index()

    # Write to csv
    # load.to_csv(inputFolder + 'demand_hourly_selected.csv')

    # Reset index for year month day
    year_month_days_selected = year_month_days_selected.reset_index()
    load_all = load_all.reset_index()

    # Estimate annual demand
    global load_annual
    load_annual = load_all.groupby(['year']).sum()
    load_annual = load_annual.drop(['month', 'day', 'hour'], axis=1).reset_index()


def load_loads():
    """
    :return:
    """

    load_dict = OrderedDict()

    for z in load.columns.to_list()[4:]:
        print(z)
        load_zone = load[z]
        load_zone = pd.concat([load_zone, period_horizon_timepoints['tmp']], axis=1)
        load_zone[['tmp']] = load_zone[['tmp']].astype(int)
        load_zone = load_zone.set_index(['tmp'])
        load_zone_dict = load_zone[z].to_dict()
        # Add stage id TODO: modify to be flexible in adding multiple stages
        load_zone_stage_dict = OrderedDict()
        load_zone_stage_dict[1] = load_zone_dict
        # load_zone_stage_dict = {1: load_zone_dict}
        load_dict[z] = load_zone_stage_dict

    # Insert data
    system_load.insert_system_static_loads(
        io=io, c=c2,
        load_scenario_id=1,
        scenario_name='base_profile_2014_forecast',
        scenario_description='Base profile forecast using 2014 base load',
        zone_stage_timepoint_static_loads=load_dict
    )


if __name__ == "__main__":
    pass
