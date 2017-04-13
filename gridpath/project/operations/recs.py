#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

"""
Get RECs for each project
"""
import csv
import os.path
from pyomo.environ import Param, Set, Expression, value

from gridpath.auxiliary.dynamic_components import required_operational_modules
from gridpath.auxiliary.auxiliary import load_operational_type_modules


def add_model_components(m, d):
    """

    :param m:
    :param d:
    :return:
    """
    # First figure out which projects are RPS-eligible
    m.RPS_PROJECTS = Set(within=m.PROJECTS)
    m.rps_zone = Param(m.RPS_PROJECTS, within=m.RPS_ZONES)

    def determine_rps_generators_by_rps_zone(mod, rps_z):
        return [p for p in mod.RPS_PROJECTS if mod.rps_zone[p] == rps_z]

    m.RPS_PROJECTS_BY_RPS_ZONE = \
        Set(m.RPS_ZONES, within=m.RPS_PROJECTS,
            initialize=determine_rps_generators_by_rps_zone)

    # Get operational RPS projects - timepoints combinations
    m.RPS_PROJECT_OPERATIONAL_TIMEPOINTS = Set(
        within=m.PROJECT_OPERATIONAL_TIMEPOINTS,
        rule=lambda mod: [(p, tmp) for (p, tmp) in
                          mod.PROJECT_OPERATIONAL_TIMEPOINTS
                          if p in mod.RPS_PROJECTS]
    )
    # Import needed operational modules
    imported_operational_modules = \
        load_operational_type_modules(getattr(d, required_operational_modules))

    def scheduled_recs_rule(mod, g, tmp):
        """
        This how many RECs are scheduled to be delivered at the timepoint
        (hourly) schedule
        :param mod:
        :param g:
        :param tmp:
        :return:
        """
        gen_op_type = mod.operational_type[g]
        return imported_operational_modules[gen_op_type]. \
            rec_provision_rule(mod, g, tmp)

    m.Scheduled_RPS_Energy_MW = Expression(
        m.RPS_PROJECT_OPERATIONAL_TIMEPOINTS, 
        rule=scheduled_recs_rule
    )

    # Keep track of curtailment
    def scheduled_curtailment_rule(mod, g, tmp):
        """
        Keep track of curtailment to make it easier to calculate total
        curtailed RPS energy for example -- this is the scheduled
        curtailment component
        :param mod:
        :param g:
        :param tmp:
        :return:
        """
        gen_op_type = mod.operational_type[g]
        return imported_operational_modules[gen_op_type]. \
            scheduled_curtailment_rule(mod, g, tmp)

    m.Scheduled_Curtailment_MW = Expression(
        m.RPS_PROJECT_OPERATIONAL_TIMEPOINTS, rule=scheduled_curtailment_rule
    )

    def subhourly_curtailment_rule(mod, g, tmp):
        """
        Keep track of curtailment to make it easier to calculate total
        curtailed RPS energy for example -- this is the subhourly
        curtailment component
        :param mod:
        :param g:
        :param tmp:
        :return:
        """
        gen_op_type = mod.operational_type[g]
        return imported_operational_modules[gen_op_type]. \
            subhourly_curtailment_rule(mod, g, tmp)

    m.Subhourly_Curtailment_MW = Expression(
        m.RPS_PROJECT_OPERATIONAL_TIMEPOINTS, rule=subhourly_curtailment_rule
    )

    def subhourly_recs_delivered_rule(mod, g, tmp):
        """
        Keep track of curtailment to make it easier to calculate total
        curtailed RPS energy for example -- this is the subhourly
        curtailment component
        :param mod:
        :param g:
        :param tmp:
        :return:
        """
        gen_op_type = mod.operational_type[g]
        return imported_operational_modules[gen_op_type]. \
            subhourly_energy_delivered_rule(mod, g, tmp)

    m.Subhourly_RPS_Energy_Delivered_MW = Expression(
        m.RPS_PROJECT_OPERATIONAL_TIMEPOINTS,
        rule=subhourly_recs_delivered_rule
    )


def load_model_data(m, d, data_portal, scenario_directory, horizon, stage):
    """

    :param m:
    :param d:
    :param data_portal:
    :param scenario_directory:
    :param horizon:
    :param stage:
    :return:
    """
    data_portal.load(filename=os.path.join(scenario_directory,
                                           "inputs", "projects.tab"),
                     select=("project", "rps_zone"),
                     param=(m.rps_zone,)
                     )

    data_portal.data()['RPS_PROJECTS'] = {
        None: data_portal.data()['rps_zone'].keys()
    }


def export_results(scenario_directory, horizon, stage, m, d):
    """

    :param scenario_directory:
    :param horizon:
    :param stage:
    :param m:
    :param d:
    :return:
    """
    with open(os.path.join(scenario_directory, horizon, stage, "results",
                           "rps_by_project.csv"), "wb") as rps_results_file:
        writer = csv.writer(rps_results_file)
        writer.writerow(["project", "load_zone", "rps_zone",
                         "timepoint", "period", "horizon", "horizon_weight",
                         "scheduled_rps_energy_mw",
                         "scheduled_curtailment_mw",
                         "subhourly_rps_energy_delivered_mw",
                         "subhourly_curtailment_mw"])
        for (p, tmp) in m.RPS_PROJECT_OPERATIONAL_TIMEPOINTS:
            writer.writerow([
                p,
                m.load_zone[p],
                m.rps_zone[p],
                tmp,
                m.period[tmp],
                m.horizon[tmp],
                m.horizon_weight[m.horizon[tmp]],
                value(m.Scheduled_RPS_Energy_MW[p, tmp]),
                value(m.Scheduled_Curtailment_MW[p, tmp]),
                value(m.Subhourly_RPS_Energy_Delivered_MW[p, tmp]),
                value(m.Subhourly_Curtailment_MW[p, tmp])
            ])

    # Export list of RPS projects and their zones for later use
    with open(os.path.join(scenario_directory, horizon, stage, "results",
                           "rps_project_zones.csv"), "wb") as \
            rps_project_zones_file:
        writer = csv.writer(rps_project_zones_file)
        writer.writerow(["project", "rps_zone"])
        for p in m.RPS_PROJECTS:
            writer.writerow([p, m.rps_zone[p]])


def get_inputs_from_database(subscenarios, c, inputs_directory):
    """

    :param subscenarios
    :param c:
    :param inputs_directory:
    :return:
    """

    project_zones = c.execute(
        """SELECT project, rps_zone
        FROM inputs_project_rps_zones
            WHERE rps_zone_scenario_id = {}
            AND project_rps_zone_scenario_id = {}""".format(
            subscenarios.RPS_ZONE_SCENARIO_ID,
            subscenarios.PROJECT_RPS_ZONE_SCENARIO_ID
        )
    ).fetchall()

    # Make a dict for easy access
    prj_zone_dict = dict()
    for (prj, zone) in project_zones:
        prj_zone_dict[str(prj)] = "." if zone is None else str(zone)

    with open(os.path.join(inputs_directory, "projects.tab"), "r"
              ) as projects_file_in:
        reader = csv.reader(projects_file_in, delimiter="\t")

        new_rows = list()

        # Append column header
        header = reader.next()
        header.append("rps_zone")
        new_rows.append(header)

        # Append correct values
        for row in reader:
            # If project specified, check if BA specified or not
            if row[0] in prj_zone_dict.keys():
                row.append(prj_zone_dict[row[0]])
                new_rows.append(row)
            # If project not specified, specify no BA
            else:
                row.append(".")
                new_rows.append(row)

    with open(os.path.join(inputs_directory, "projects.tab"), "w") as \
            projects_file_out:
        writer = csv.writer(projects_file_out, delimiter="\t")
        writer.writerows(new_rows)


def import_results_into_database(scenario_id, c, db, results_directory):
    """
    Assign rps_zone to appropriate results tables
    :param scenario_id:
    :param c:
    :param db:
    :param results_directory:
    :return:
    """
    print("update rps zones")
    # Figure out RPS zone for each project
    prj_rps_zones = dict()
    with open(os.path.join(results_directory, "rps_project_zones.csv"),
              "r") as prj_rps_zones_file:
        reader = csv.reader(prj_rps_zones_file)
        reader.next()
        for row in reader:
            prj_rps_zones[row[0]] = row[1]

    # Update tables with RPS zone
    tables_to_update = [
        "results_project_capacity_all",
        "results_project_dispatch_all",
        "results_project_dispatch_variable",
        "results_project_dispatch_capacity_commit",
        "results_project_costs_capacity",
        "results_project_costs_operations_variable_om",
        "results_project_costs_operations_fuel",
        "results_project_costs_operations_startup",
        "results_project_costs_operations_shutdown"
    ]
    for prj in prj_rps_zones.keys():
        for tbl in tables_to_update:
            c.execute(
                """UPDATE {}
                SET rps_zone = '{}'
                WHERE scenario_id = {}
                AND project = '{}';""".format(
                    tbl,
                    prj_rps_zones[prj],
                    scenario_id,
                    prj
                )
            )
    db.commit()

    # Set rps_zone to 'no_rps' for all other projects
    # This helps for later joins (can't join on NULL values)
    for tbl in tables_to_update:
        c.execute(
            """UPDATE {}
                            SET rps_zone = 'no_rps'
                            WHERE scenario_id = {}
                            AND rps_zone IS NULL;""".format(
                tbl,
                scenario_id
            )
        )
    db.commit()