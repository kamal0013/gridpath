#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

"""
Simple results dashboard

Plan:
Different tabs, each one has something (cost, capacity, etc.)
Also have selection menu to pick a zone (or multiple?) of interest
Also have selection menu to pick a stage
Also have selection to pick period of interest (or multiple)

TODO:
 for testing can simply mock up some dataframes so we can test the callbacks
 and test the layout stuff. Can then later worry about how to get data
 in exactly the right format

Seems like we can use CDSView and JSfilter to do simple filters

Next steps:
 - mock up something with different tabs and a chart in each tab
 mock up a cpaacity expansion example with multiple zones
 (mock up a prod cost example with multiple zones)
 (cost tab, capacity tab, etc.)
 - after that, try and add some sort of customjs filter to a tab so we can
 figure out which zone we want to see for e.g. capacity results
 -
"""

from argparse import ArgumentParser
from bokeh.embed import json_item

import pandas as pd
import sys

# GridPath modules
from db.common_functions import connect_to_database
from gridpath.auxiliary.auxiliary import get_scenario_id_and_name
from viz.common_functions import create_stacked_bar_plot, show_plot, \
    get_parent_parser, get_unit


def create_parser():
    """

    :return:
    """
    parser = ArgumentParser(add_help=True, parents=[get_parent_parser()])
    parser.add_argument("--scenario_id", help="The scenario ID. Required if "
                                              "no --scenario is specified.")
    parser.add_argument("--scenario", help="The scenario name. Required if "
                                           "no --scenario_id is specified.")

    return parser


def parse_arguments(arguments):
    """

    :return:
    """
    parser = create_parser()
    parsed_arguments = parser.parse_args(args=arguments)

    return parsed_arguments


def get_cost_data(conn, scenario_id, **kwargs):
    """
    Get costs results by period and component for a given scenario

    **kwargs needed, so that an error isn't thrown when calling this
    function with extra arguments from the UI.

    :param conn:
    :param scenario_id:
    :return:
    """
    # Note: fuel cost and variable O&M cost are actually cost *rates* in $/hr
    #  and should be multiplied by the timepoint duration to get the actual
    #  cost.

    # System costs by scenario and period -- by source and total
    # Spinup/lookahead timepoints are ignored by adding the resp. column tag
    # through inner joins and adding a conditional to ignore those timepoints

    # TODO: transmission hurdle rates are now assigned to the load_zone_to
    #  zone but should treat this smarter
    sql = """SELECT period, stage_id, load_zone,
        capacity_cost/1000000 as Capacity,
        fuel_cost/1000000 as Fuel,
        variable_om_cost/1000000 as Variable_OM,
        startup_cost/1000000 as Startups,
        shutdown_cost/1000000 as Shutdowns,
        hurdle_cost/1000000 as Hurdle_Rates

        FROM
        
        (SELECT scenario_id, period, stage_id, load_zone, 
        sum(capacity_cost) AS capacity_cost
        FROM results_project_costs_capacity
        WHERE scenario_id = ?
        GROUP BY scenario_id, period, stage_id, load_zone) AS cap_costs

        LEFT JOIN

        (SELECT scenario_id, period, stage_id, load_zone,
        sum(fuel_cost * timepoint_weight * number_of_hours_in_timepoint) 
        AS fuel_cost,
        sum(variable_om_cost * timepoint_weight * number_of_hours_in_timepoint) 
        AS variable_om_cost,
        sum(startup_cost * timepoint_weight) AS startup_cost,
        sum(shutdown_cost * timepoint_weight) AS shutdown_cost
        FROM results_project_costs_operations
        
        -- add temporal scenario id so we can join timepoints table
        INNER JOIN
        
        (SELECT temporal_scenario_id, scenario_id FROM scenarios)
        USING (scenario_id)
        
        -- filter out spinup_or_lookahead timepoints
        INNER JOIN
        
        (SELECT temporal_scenario_id, stage_id, subproblem_id, timepoint, 
        spinup_or_lookahead
        FROM inputs_temporal)
        USING (temporal_scenario_id, stage_id, subproblem_id, timepoint)
        
        WHERE scenario_id = ?
        AND spinup_or_lookahead is NULL
        
        GROUP BY scenario_id, period, stage_id, load_zone) AS operational_costs
        USING (scenario_id, period, stage_id, load_zone)

        LEFT JOIN

        (SELECT scenario_id, period, stage_id, load_zone_to AS load_zone,
        sum((hurdle_cost_positive_direction + hurdle_cost_negative_direction) * 
        timepoint_weight * number_of_hours_in_timepoint) AS hurdle_cost
        FROM
        results_transmission_hurdle_costs
        
        -- add temporal scenario id so we can join timepoints table
        INNER JOIN
        
        (SELECT temporal_scenario_id, scenario_id FROM scenarios)
        USING (scenario_id)
        
        -- filter out spinup_or_lookahead timepoints
        INNER JOIN
        
        (SELECT temporal_scenario_id, stage_id, subproblem_id, timepoint, 
        spinup_or_lookahead
        FROM inputs_temporal)
        USING (temporal_scenario_id, stage_id, subproblem_id, timepoint)
        
        WHERE scenario_id = ?
        AND spinup_or_lookahead is NULL
        
        GROUP BY scenario_id, period, stage_id, load_zone) AS hurdle_costs
        USING (scenario_id, period, stage_id, load_zone)
        ;"""

    df = pd.read_sql(
        sql,
        con=conn,
        params=(scenario_id, scenario_id, scenario_id)
    )

    return df


def main(args=None):
    """
    Parse the arguments, get the data in a df, and create the plot

    :return: if requested, return the plot as JSON object
    """
    if args is None:
        args = sys.argv[1:]
    parsed_args = parse_arguments(arguments=args)

    conn = connect_to_database(db_path=parsed_args.database)
    c = conn.cursor()

    scenario_id, scenario = get_scenario_id_and_name(
        scenario_id_arg=parsed_args.scenario_id,
        scenario_name_arg=parsed_args.scenario,
        c=c,
        script="cost_plot"
    )

    cost_unit = "million " + get_unit(c, "cost")

    plot_title = "test"
    plot_name = "Dashboard"


    cost_df = pd.DataFrame(
        data={'period': [2020, 2020, 2030, 2030, 2040, 2040],
              'stage': [1, 1, 1, 1, 1, 1],
              'load_zone': ['zone1', 'zone2', 'zone1', 'zone2', 'zone1', 'zone2'],
              'capacity_cost': [10, 0, 20, 0, 30, 20],
              'fuel_cost': [10, 0, 20, 0, 30, 20],
              'vom_cost': [10, 0, 20, 0, 30, 20],
              'startup_cost': [10, 0, 20, 0, 30, 20],
              'shutdown_cost': [10, 0, 20, 0, 30, 20]
              }
    )
    # cost_df = get_cost_data(
    #     conn=conn,
    #     scenario_id=scenario_id
    # )

    print(cost_df.head())

    cap_df = pd.DataFrame(
        columns=['period', 'stage', 'load_zone', 'technology',
                 'new_cap', 'cumulative_new_cap', 'retired_cap', 'cumulative_retired_cap', 'total_cap'],
        data=[[2020, 1, 'Zone1', 'Solar', 10, 10, 0, 0, 10],
              [2030, 1, 'Zone1', 'Solar', 10, 20, 0, 0, 20],
              [2040, 1, 'Zone1', 'Solar', 10, 30, 0, 0, 30],
              [2020, 1, 'Zone1', 'Wind', 10, 10, 0, 0, 10],
              [2030, 1, 'Zone1', 'Wind', 10, 20, 0, 0, 20],
              [2040, 1, 'Zone1', 'Wind', 10, 30, 0, 0, 30],
              [2020, 1, 'Zone2', 'Solar', 0, 0, 0, 0, 0],
              [2030, 1, 'Zone2', 'Solar', 10, 10, 0, 0, 20],
              [2040, 1, 'Zone2', 'Solar', 10, 20, 0, 0, 30],
              [2020, 1, 'Zone2', 'Wind', 20, 20, 0, 0, 20],
              [2030, 1, 'Zone2', 'Wind', 10, 30, 0, 0, 30],
              [2040, 1, 'Zone2', 'Wind', 10, 40, 0, 0, 40]
              ]
    )
    print(cap_df.head())

    tabs = Tabs(tabs=[cost, cap])

    # plot = create_stacked_bar_plot(
    #     df=df,
    #     title=plot_title,
    #     y_axis_column="Cost (MW)",
    #     x_axis_column="period",
    #     group_column="Cost Component",
    #     column_mapper={"Cost (MW)": "Cost ({})".format(cost_unit),
    #                    "period": "Period"},
    #     ylimit=parsed_args.ylimit
    # )
    #
    # # Show plot in HTML browser file if requested
    # if parsed_args.show:
    #     show_plot(plot=plot,
    #               plot_name=plot_name,
    #               plot_write_directory=parsed_args.plot_write_directory,
    #               scenario=scenario)
    #
    # # Return plot in json format if requested
    # if parsed_args.return_json:
    #     return json_item(plot, plot_name)


if __name__ == "__main__":
    main()
