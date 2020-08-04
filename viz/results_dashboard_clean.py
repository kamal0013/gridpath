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
 this requires adding a columndatasource and some type of filtering!
 -
"""

from argparse import ArgumentParser
from bokeh.embed import json_item
from bokeh.models import Tabs, Panel, PreText
from bokeh.io import show
from bokeh.layouts import column, row

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


# TODO: step proces
# 1. this function gives you certain zone and you put it in a tab with the costs
# 2. allow for different zones to be picked dyncamically
#   architecture: like in example where you have an update function that updates the CDS?
#   --> problem is that we currently don't run function off a CDS but off a dataframe
#   in update function, simply


# TODO: add ability to pick certain load zone (or period?)
#  and add an "update" function that pulls this data again for different info
#  It seems like using CDSView and CustomJSfilter is kind of a nuisance, see here:
#  https://stackoverflow.com/questions/58868542/update-cdsview-filter-with-customjs-for-line-glyph
#  -> plot all load_zones (or periods or ...), and jscallback decides which one to show
def get_cost_data(conn, scenario_id, **kwargs):
    df = pd.DataFrame(
        data={'period': [2020, 2020, 2030, 2030, 2040, 2040],
              'stage': [1, 1, 1, 1, 1, 1],
              'load_zone': ['Zone1', 'Zone2', 'Zone1', 'Zone2', 'Zone1', 'Zone2'],
              'capacity_cost': [10, 0, 20, 0, 30, 20],
              'fuel_cost': [10, 0, 20, 0, 30, 20],
              'vom_cost': [10, 0, 20, 0, 30, 20],
              'startup_cost': [10, 0, 20, 0, 30, 20],
              'shutdown_cost': [10, 0, 20, 0, 30, 20]
              }
    )
    df = df[df['load_zone'] == 'Zone1']

    df = pd.melt(
        df,
        id_vars=['period', 'stage'],
        value_vars=['capacity_cost', 'fuel_cost', 'vom_cost',
                    'startup_cost', 'shutdown_cost'],
        var_name='Cost Component',
        value_name='Cost ()'
    )

    print(df.head())

    return df


def get_cap_data(conn, scenario_id, **kwargs):
    df = pd.DataFrame(
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

    df = df[df['load_zone'] == 'Zone1']
    df = df[['period', 'technology', 'total_cap']]

    print(df.head())
    return df


def get_energy_data(conn, scenario_id, **kwargs):
    df = pd.DataFrame(
        columns=['period', 'stage', 'load_zone', 'technology', 'energy'],
        data=[[2020, 1, 'Zone1', 'Solar', 10],
              [2030, 1, 'Zone1', 'Solar', 10],
              [2040, 1, 'Zone1', 'Solar', 10],
              [2020, 1, 'Zone1', 'Wind', 10],
              [2030, 1, 'Zone1', 'Wind', 10],
              [2040, 1, 'Zone1', 'Wind', 10],
              [2020, 1, 'Zone2', 'Solar', 0],
              [2030, 1, 'Zone2', 'Solar', 10],
              [2040, 1, 'Zone2', 'Solar', 10],
              [2020, 1, 'Zone2', 'Wind', 20],
              [2030, 1, 'Zone2', 'Wind', 10],
              [2040, 1, 'Zone2', 'Wind', 10]
              ]
    )

    df = df[df['load_zone'] == 'Zone1']
    df = df[['period', 'technology', 'energy']]

    print(df.head())
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
        script="dashboard"
    )

    cost_unit = "million " + get_unit(c, "cost")
    power_unit = get_unit(c, "power")

    plot_name = "Dashboard"

    # TODO: depending on selection, slice out load zone of interest

    cost_df = get_cost_data(conn, scenario_id)
    cap_df = get_cap_data(conn, scenario_id)
    energy_df = get_energy_data(conn, scenario_id)

    cost = create_stacked_bar_plot(
        df=cost_df,
        title="Annual Cost by Component",
        y_axis_column="Cost ()",
        x_axis_column="period",
        group_column="Cost Component",
        column_mapper={"Cost ()": "Cost ({})".format(cost_unit),
                       "period": "Period"},
        ylimit=parsed_args.ylimit
    )

    cap = create_stacked_bar_plot(
        df=cap_df,
        title="Capacity by Technology",
        y_axis_column="total_cap",  # TODO: allows this to change with select button
        x_axis_column="period",
        group_column="technology",
        column_mapper={"total_cap": "Capacity ({})".format(power_unit),
                       "period": "Period"},
        ylimit=parsed_args.ylimit
    )

    energy = create_stacked_bar_plot(
        df=energy_df,
        title="Annual Energy by Technology",
        y_axis_column="energy",  # TODO: allow to have toggles for including curtailment and whatnot?
        x_axis_column="period",
        group_column="technology",
        column_mapper={"total_cap": "Capacity ({})".format(power_unit),
                       "period": "Period"},
        ylimit=parsed_args.ylimit
    )

    # TODO: link to an actual table with summary stats!
    stats = PreText(text=cap_df.to_string())

    # cost = Panel(child=cost, title='cap')
    # cap = Panel(child=cap, title='cost')
    # tabs = Tabs(tabs=[cost, cap])

    # set up layout
    cap_energy = row(cap, energy)
    layout = column(stats, cost, cap_energy)

    show(layout)



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
