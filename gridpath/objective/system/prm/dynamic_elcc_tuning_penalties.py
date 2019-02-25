#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

from builtins import next
import csv
import os.path
from pyomo.environ import Param, Expression, NonNegativeReals

from gridpath.auxiliary.dynamic_components import total_cost_components


def add_model_components(m, d):
    """

    :param m:
    :param d:
    :return:
    """

    m.dynamic_elcc_tuning_cost = Param(default=0)

    def total_elcc_tuning_cost_rule(mod):
        """
        Set dynamic elcc to max available by subtracting a small amount from 
        the objective function when Dynamic_ELCC is higher
        :param mod:
        :return:
        """
        if mod.dynamic_elcc_tuning_cost == 0:
            return 0
        else:
            return - sum(
                mod.Dynamic_ELCC_MW[z, p]
                * mod.dynamic_elcc_tuning_cost
                * mod.number_years_represented[p]
                * mod.discount_factor[p]
                for (z, p)
                in mod.PRM_ZONE_PERIODS_WITH_REQUIREMENT
            )

    m.Total_Dynamic_ELCC_Tuning_Cost = Expression(
        rule=total_elcc_tuning_cost_rule
    )
    getattr(d, total_cost_components).append("Total_Dynamic_ELCC_Tuning_Cost")
    
    
def load_model_data(m, d, data_portal, scenario_directory, horizon, stage):
    """
    Get tuning param value from file if file exists
    :param m:
    :param d:
    :param data_portal:
    :param scenario_directory:
    :param horizon:
    :param stage:
    :return:
    """
    tuning_param_file = os.path.join(
        scenario_directory, horizon, stage, "inputs", "tuning_params.tab"
    )

    if os.path.exists(tuning_param_file):
        data_portal.load(filename=tuning_param_file,
                         select=("dynamic_elcc_tuning_cost",),
                         param=m.dynamic_elcc_tuning_cost
                         )
    else:
        pass


def get_inputs_from_database(subscenarios, c, inputs_directory):
    """

    :param subscenarios
    :param c:
    :param inputs_directory:
    :return:
    """

    dynamic_elcc_tuning_cost = c.execute(
        """SELECT dynamic_elcc_tuning_cost
        FROM inputs_tuning
        WHERE tuning_scenario_id = {}""".format(
            subscenarios.TUNING_SCENARIO_ID
        )
    ).fetchone()[0]

    # If tuning params file exists, add column to file, else create file and
    #  writer header and tuning param value
    if os.path.isfile(os.path.join(inputs_directory, "tuning_params.tab")):
        with open(os.path.join(inputs_directory, "tuning_params.tab"), "r"
                  ) as tuning_params_file_in:
            reader = csv.reader(tuning_params_file_in, delimiter="\t")

            new_rows = list()

            # Append column header
            header = next(reader)
            header.append("dynamic_elcc_tuning_cost")
            new_rows.append(header)

            # Append tuning param value
            param_value = next(reader)
            param_value.append(dynamic_elcc_tuning_cost)
            new_rows.append(param_value)

        with open(os.path.join(inputs_directory, "tuning_params.tab"),
                  "w") as \
                tuning_params_file_out:
            writer = csv.writer(tuning_params_file_out, delimiter="\t")
            writer.writerows(new_rows)

    else:
        with open(os.path.join(inputs_directory, "tuning_params.tab"),
                  "w") as \
                tuning_params_file_out:
            writer = csv.writer(tuning_params_file_out, delimiter="\t")
            writer.writerows(["dynamic_elcc_tuning_cost"])
            writer.writerows([dynamic_elcc_tuning_cost])