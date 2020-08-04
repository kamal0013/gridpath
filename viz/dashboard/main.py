#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

"""
"""

from argparse import ArgumentParser
from bokeh.embed import json_item
from bokeh.models import Tabs, Panel, PreText, Select, ColumnDataSource, Legend
from bokeh.io import curdoc, show
from bokeh.layouts import column, row
from bokeh.plotting import figure
from bokeh.palettes import cividis


import pandas as pd
import sys

# GridPath modules
from db.common_functions import connect_to_database
from gridpath.auxiliary.auxiliary import get_scenario_id_and_name
from viz.common_functions import create_stacked_bar_plot, show_plot, \
    get_parent_parser, get_unit, show_hide_legend

# TODO: base on actual data
ZONE_OPTIONS = ['Zone1', 'Zone2']  # todo: Add "all" option
CAPACITY_OPTIONS = ["new_cap", "cumulative_new_cap", "total_cap"]
COST_COMPONENTS = ["capacity_cost", "fuel_cost", "vom_cost", "startup_cost",
                   "shutdown_cost"]


# Data gathering functions
def get_summary(conn, scenario_id):
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

    return df


def get_all_cost_data(conn, scenario_id):
    # get data for all zones and periods
    # TODO: use COST COMPONENTS var?
    idx_cols = ['period', 'stage', 'load_zone']
    value_cols = ['capacity_cost', 'fuel_cost', 'vom_cost',
                  'startup_cost', 'shutdown_cost']
    df = pd.DataFrame(
        columns=['period', 'stage', 'load_zone',
                 'capacity_cost', 'fuel_cost', 'vom_cost',
                 'startup_cost', 'shutdown_cost'],
        data=[[2020, 1, 'Zone1', 10, 10, 10, 10, 10],
              [2020, 1, 'Zone2', 0, 0, 0, 0, 0],
              [2030, 1, 'Zone1', 20, 20, 20, 20, 20],
              [2030, 1, 'Zone2', 0, 0, 0, 0, 0],
              [2040, 1, 'Zone1', 30, 30, 30, 30, 30],
              [2040, 1, 'Zone2', 20, 20, 20, 20, 20]
              ]
    )

    df['period'] = df['period'].astype(str)  # for categorical axis in Bokeh

    return df, idx_cols, value_cols


def get_all_capacity_data(conn, scenario_id):
    # get data for all zones and periods
    df = pd.DataFrame(
        columns=['period', 'stage', 'load_zone', 'technology',
                 'new_cap', 'cumulative_new_cap', 'retired_cap',
                 'cumulative_retired_cap', 'total_cap'],
        data=[[2020, 1, 'Zone1', 'Solar', 10, 10, 0, 0, 50],
              [2030, 1, 'Zone1', 'Solar', 10, 20, 0, 0, 60],
              [2040, 1, 'Zone1', 'Solar', 10, 30, 0, 0, 70],
              [2020, 1, 'Zone1', 'Wind', 10, 10, 0, 0, 20],
              [2030, 1, 'Zone1', 'Wind', 10, 20, 0, 0, 30],
              [2040, 1, 'Zone1', 'Wind', 10, 30, 0, 0, 40],
              [2020, 1, 'Zone2', 'Solar', 0, 0, 0, 0, 20],
              [2030, 1, 'Zone2', 'Solar', 10, 10, 0, 0, 30],
              [2040, 1, 'Zone2', 'Solar', 10, 20, 0, 0, 40],
              [2020, 1, 'Zone2', 'Wind', 20, 20, 0, 0, 30],
              [2030, 1, 'Zone2', 'Wind', 10, 30, 0, 0, 40],
              [2040, 1, 'Zone2', 'Wind', 10, 40, 0, 0, 50]
              ]
    )

    # 'Unpivot' from wide to long format
    df = pd.melt(
        df,
        id_vars=['period', 'stage', 'load_zone', 'technology'],
        var_name='capacity_metric',
        value_name='capacity'
    )

    # Pivot technologies to wide format (for stack chart)
    # Note: df.pivot does not work with multi-index as of pandas 1.0.5
    idx_cols = ['period', 'stage', 'load_zone', 'capacity_metric']
    value_cols = list(df['technology'].unique())
    df = pd.pivot_table(
        df,
        index=idx_cols,
        columns='technology',
        values='capacity'
    ).fillna(0).reset_index()

    df['period'] = df['period'].astype(str)  # for categorical axis in Bokeh

    return df, idx_cols, value_cols


def get_all_energy_data(conn, scenario_id):
    # get data for all zones and periods
    df = pd.DataFrame(
        columns=['period', 'stage', 'load_zone', 'technology', 'energy'],
        data=[[2020, 1, 'Zone1', 'Solar', 10],
              [2030, 1, 'Zone1', 'Solar', 10],
              [2040, 1, 'Zone1', 'Solar', 10],
              [2020, 1, 'Zone1', 'Wind', 20],
              [2030, 1, 'Zone1', 'Wind', 20],
              [2040, 1, 'Zone1', 'Wind', 20],
              [2020, 1, 'Zone2', 'Solar', 0],
              [2030, 1, 'Zone2', 'Solar', 10],
              [2040, 1, 'Zone2', 'Solar', 20],
              [2020, 1, 'Zone2', 'Wind', 20],
              [2030, 1, 'Zone2', 'Wind', 30],
              [2040, 1, 'Zone2', 'Wind', 40]
              ]
    )

    # Pivot technologies to wide format (for stack chart)
    # Note: df.pivot does not work with multi-index as of pandas 1.0.5
    idx_cols = ['period', 'stage', 'load_zone']
    value_cols = list(df['technology'].unique())
    df = pd.pivot_table(
        df,
        index=idx_cols,
        columns='technology',
        values='energy'
    ).fillna(0).reset_index()

    df['period'] = df['period'].astype(str)  # for categorical axis in Bokeh

    return df, idx_cols, value_cols


# Plot creation functions
def get_stacked_plot_by_period(source, stackers, title):
    print(source.data)
    print(stackers)
    print(source.data['period'])

    # TODO: configure layout better (less tall)
    # TODO: make sure x is categoricla (now source shows numbers)
    # TODO: fix legend (shows numbers?)
    # TODO: add axes
    plot = figure(title=title, x_range=source.data['period'])
    area_renderers = plot.vbar_stack(
        stackers=stackers,
        x='period',
        source=source,
        color=cividis(len(stackers)),
        width=0.5
    )
    # Note: cannot use legend=stackers because there is an issue when the legend
    # and column names are the same, see here:
    # https://github.com/bokeh/bokeh/issues/5365

    # Add Legend
    legend_items = [(y, [area_renderers[i]]) for i, y in enumerate(stackers)
                    if source.data[y].mean() > 0]
    legend = Legend(items=legend_items)
    plot.add_layout(legend, 'right')
    plot.legend.title = 'Legend Title'
    plot.legend[0].items.reverse()  # Reverse legend to match stacked order
    plot.legend.click_policy = 'hide'  # Add interactivity to the legend
    show_hide_legend(plot=plot)  # Hide legend on double click

    return plot


# Define callbacks
def zone_change(attr, old, new):
    """
    When the selected load zone changes, update the cost, energy and capacity.
    """
    update_cost_data(zone=new)
    update_energy_data(zone=new)
    update_capacity_data(zone=new, capacity_metric=capacity_select.value)


def capacity_change(attr, old, new):
    """
    When the selected capacity metric changes, update the capacity.
    """
    update_capacity_data(zone=zone_select.value, capacity_metric=new)


# TODO: deal with default zone more generally (dynamic link or "All" or first
#  one)
def update_cost_data(zone="Zone1"):
    """
    Update the ColumnDataSource object 'cost_source' with the appropriate
    load_zone slice of the data.
    :param zone:
    :return:
    """
    slice = cost[cost["load_zone"] == zone]
    cost_source.data = slice


def update_energy_data(zone="Zone1"):
    """
    Update the ColumnDataSource object 'energy_source' with the appropriate
    load_zone slice of the data.
    :param zone:
    :return:
    """
    slice = energy[energy["load_zone"] == zone]
    energy_source.data = slice


def update_capacity_data(zone="Zone1", capacity_metric="new_cap"):
    """
    Update the ColumnDataSource object 'cost_source' with the appropriate
    load_zone slice of the data, and with the appropriate capacity metric.
    :param zone:
    :param capacity_metric:
    :return:
    """
    slice = capacity[(capacity["load_zone"] == zone)
                     & (capacity["capacity_metric"] == capacity_metric)]
    capacity_source.data = slice


# TODO: update summary too?
# TODO: need to dynamically link initial zone value (or set to "all"?)
# Set up widgets
summary_table = PreText(text='empty PreText', width=500)
zone_select = Select(title="Load Zone:", value='Zone1', options=ZONE_OPTIONS)
capacity_select = Select(title="Capacity Metric:", value='new_cap',
                         options=CAPACITY_OPTIONS)

# Get the data
# TODO: get actual connection and data
conn = "dbtest"
scenario_id = 1
summary = get_summary(conn, scenario_id)
cost, cost_idx_cols, cost_value_cols = get_all_cost_data(conn, scenario_id)
capacity, capacity_idx_cols, capacity_value_cols = \
    get_all_capacity_data(conn, scenario_id)
energy, energy_idx_cols, energy_value_cols = \
    get_all_energy_data(conn, scenario_id)

# Set up CDS and update with default values
# todo: is there better way to set up initial CDS?
summary_source = ColumnDataSource(data=summary)
cost_source = ColumnDataSource()
energy_source = ColumnDataSource()
capacity_source = ColumnDataSource()

# TODO: maybe deal with defaults here? (e.g. update with first zone in list)
update_cost_data()
update_capacity_data()
update_energy_data()
# TODO: update summary table here too!

# TRY AND PASS AROUND DF INSTEAD OF SOURCE AND SEE IF IT WORKS?
# Set up plots
cost_plot = get_stacked_plot_by_period(
    cost_source,
    stackers=cost_value_cols,
    title="Cost by Component"
)
energy_plot = get_stacked_plot_by_period(
    energy_source,
    stackers=energy_value_cols,
    title="Energy by Technology"
)
cap_plot = get_stacked_plot_by_period(
    capacity_source,
    stackers=capacity_value_cols,
    title="Capacity by Technology"
)

# Set up callback behavior
zone_select.on_change('value', zone_change)
capacity_select.on_change('value', capacity_change)

# Set up layout
widgets = column(zone_select, capacity_select)
top_row = row(summary_table, widgets)
middle_row = row(cap_plot, energy_plot)
bottom_row = row(cost_plot)
layout = column(top_row, middle_row, bottom_row)

curdoc().add_root(layout)
curdoc().title = "Dashboard"