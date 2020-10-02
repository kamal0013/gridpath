#!/usr/bin/env python
# Copyright 2017 Blue Marble Analytics LLC. All rights reserved.

"""
This script iterates over all modules required for a GridPath scenario and
calls their *validate_inputs()* method, which performs various validations
of the input data and scenario setup.
"""

import pandas as pd
import sqlite3
import sys
from argparse import ArgumentParser

from db.common_functions import connect_to_database, spin_on_database_lock
from gridpath.auxiliary.auxiliary import get_scenario_id_and_name
from gridpath.auxiliary.validations import write_validation_to_database, \
    validate_cols_equal
from gridpath.common_functions import get_db_parser
from gridpath.auxiliary.module_list import determine_modules, load_modules
from gridpath.auxiliary.scenario_chars import OptionalFeatures, SubScenarios, \
    SubProblems


def validate_inputs(subproblems, loaded_modules, scenario_id, subscenarios, conn):
    """"
    For each module, load the inputs from the database and validate them

    :param subproblems: SubProblems object with info on the subproblem/stage
        structure
    :param loaded_modules: list of imported modules (Python <class 'module'>
        objects)
    :param subscenarios: SubScenarios object with all subscenario info
    :param conn: database connection
    :return:
    """

    # TODO: check if we even need database cursor (and subscenarios?)
    #   since presumably we already have our data in the inputs.
    #   need to go through each module's input validation to check this

    # TODO: see if we can do some sort of automatic dtype validation for
    #  each table in the database? Problem is that you don't necessarily want
    #  to check the full table but only the appropriate subscenario

    subproblems_list = subproblems.SUBPROBLEMS
    for subproblem in subproblems_list:
        stages = subproblems.SUBPROBLEM_STAGE_DICT[subproblem]
        for stage in stages:
            # 1. input validation within each module
            for m in loaded_modules:
                if hasattr(m, "validate_inputs"):
                    m.validate_inputs(
                        scenario_id=scenario_id,
                        subscenarios=subscenarios,
                        subproblem=subproblem,
                        stage=stage,
                        conn=conn
                    )
                else:
                    pass

            # 2. input validation across modules
            #    make sure geography and projects are in line
            #    ... (see Evernote validation list)
            #    create separate function for each validation that you call here

    # Validation across subproblems and stages:
    validate_hours_in_subproblem_period(scenario_id, subscenarios, conn)


def validate_hours_in_subproblem_period(scenario_id, subscenarios, conn):
    """

    :param subscenarios:
    :param conn:
    :return:
    """

    sql = """
        SELECT stage_id, period, n_hours, hours_in_full_period
        FROM
            (SELECT temporal_scenario_id, stage_id, period,
            sum(number_of_hours_in_timepoint * timepoint_weight) as n_hours
            FROM inputs_temporal
            WHERE spinup_or_lookahead IS NULL
            AND temporal_scenario_id = {}
            GROUP BY temporal_scenario_id, stage_id, period
            ) AS tmp_table
    
        INNER JOIN
        
        (SELECT temporal_scenario_id, period, hours_in_full_period
        FROM inputs_temporal_periods) AS period_tbl
    
        USING (temporal_scenario_id, period)
    
    """.format(subscenarios.TEMPORAL_SCENARIO_ID)

    df = pd.read_sql(sql, con=conn)

    msg = """Total number of hours in timepoints adjusted for timepoint weight
          and duration and excluding spinup and lookahead timepoints should 
          match hours_in_full_period."""
    write_validation_to_database(
        conn=conn,
        scenario_id=scenario_id,
        subproblem_id="N/A",
        stage_id="N/A",
        gridpath_module="N/A",
        db_table="inputs_temporal",
        severity="Mid",
        errors=validate_cols_equal(
            df=df,
            col1="n_hours",
            col2="hours_in_full_period",
            idx_col=["stage_id", "period"],
            msg=msg
        )
    )

    # TODO: check that subproblems don't straddle periods


def validate_subscenario_ids(scenario_id, subscenarios, optional_features, conn):
    """
    Check whether subscenarios_ids are consistent with:
     - core required subscenario_ids
     - data dependent subscenario_ids (e.g. new build)
     - optional features
    data.

    :param subscenarios:
    :param optional_features:
    :param conn:
    :return: Boolean, True is all subscenario IDs are valid.
    """

    valid_core = validate_required_subscenario_ids(scenario_id, subscenarios, conn)

    if valid_core:
        valid_data_dependent = validate_data_dependent_subscenario_ids(
            scenario_id, subscenarios, conn)
    else:
        valid_data_dependent = False

    valid_feature = validate_feature_subscenario_ids(
        scenario_id, subscenarios, optional_features, conn)

    return valid_core and valid_data_dependent and valid_feature


def validate_feature_subscenario_ids(scenario_id, subscenarios, optional_features, conn):
    """

    :param subscenarios:
    :param optional_features:
    :param conn:
    :return:
    """

    subscenario_ids_by_feature = determine_subscenarios_by_feature(conn)
    feature_list = optional_features.get_active_features()

    errors = {"High": [], "Low": []}  # errors by severity
    for feature, subscenario_ids in subscenario_ids_by_feature.items():
        if feature not in ["core", "optional", "data_dependent"]:
            for sc_id in subscenario_ids:
                # If the feature is requested, and the associated subscenarios
                # are not specified, raise a validation error
                if feature in feature_list and \
                        getattr(subscenarios, sc_id) == "NULL":
                    errors["High"].append(
                        "Requested feature '{}' requires an input for '{}'"
                        .format(feature, sc_id)
                    )

                # If the feature is not requested, and the associated
                # subscenarios are specified, raise a validation error
                elif feature not in feature_list and \
                        getattr(subscenarios, sc_id) != "NULL":
                    errors["Low"].append(
                        "Detected inputs for '{}' while related feature '{}' "
                         "is not requested".format(sc_id, feature)
                    )
    for severity, error_list in errors.items():
        write_validation_to_database(
            conn=conn,
            scenario_id=scenario_id,
            subproblem_id="N/A",
            stage_id="N/A",
            gridpath_module="N/A",
            db_table="scenarios",
            severity=severity,
            errors=error_list
        )

    # Return True if all required subscenario_ids are valid (list is empty)
    return not bool(sum(errors.values(), []))


def validate_required_subscenario_ids(scenario_id, subscenarios, conn):
    """
    Check whether the required subscenario_ids are specified in the db
    :param subscenarios:
    :param conn:
    :return: boolean, True if all required subscenario_ids are specified
    """

    required_subscenario_ids = determine_subscenarios_by_feature(conn)["core"]

    errors = []
    for sc_id in required_subscenario_ids:
        if getattr(subscenarios, sc_id) is None:
            errors.append(
                "'{}' is a required input in the 'scenarios' table"
                .format(sc_id)
            )

    write_validation_to_database(
        conn=conn,
        scenario_id=scenario_id,
        subproblem_id="N/A",
        stage_id="N/A",
        gridpath_module="N/A",
        db_table="scenarios",
        severity="High",
        errors=errors
    )

    # Return True if all required subscenario_ids are valid (list is empty)
    return not bool(errors)


def validate_data_dependent_subscenario_ids(scenario_id, subscenarios, conn):
    """

    :param subscenarios:
    :param conn:
    :return:
    """

    assert subscenarios.PROJECT_PORTFOLIO_SCENARIO_ID is not None

    req_cap_types = set(
        get_required_capacity_type_modules(conn, scenario_id)
    )

    new_build_types = {
        "gen_new_lin,",
        "gen_new_bin",
        "stor_new_lin",
        "stor_new_bin",
        "dr_new"
    }
    existing_build_types = {
        "gen_spec",
        "gen_ret_bin",
        "gen_ret_lin",
        "stor_spec",
    }
    dr_types = {
        "dr_new"
    }

    # Determine required subscenario_ids
    sc_id_type = []
    if bool(req_cap_types & new_build_types):
        sc_id_type.append(("PROJECT_NEW_COST_SCENARIO_ID",
                           "New Build"))
    if bool(req_cap_types & existing_build_types):
        sc_id_type.append(("PROJECT_SPECIFIED_CAPACITY_SCENARIO_ID",
                           "Existing"))
        sc_id_type.append(("PROJECT_SPECIFIED_FIXED_COST_SCENARIO_ID",
                           "Existing"))
    if bool(req_cap_types & dr_types):
        sc_id_type.append(("PROJECT_NEW_POTENTIAL_SCENARIO_ID",
                           "Demand Response (DR)"))

    # Check whether required subscenario_ids are present
    errors = []
    for sc_id, build_type in sc_id_type:
        if getattr(subscenarios, sc_id) is None:
            errors.append(
                 "'{}' is a required input in the 'scenarios' table if there "
                 "are '{}' resources in the portfolio".format(sc_id, build_type)
            )

    write_validation_to_database(
        conn=conn,
        scenario_id=scenario_id,
        subproblem_id="N/A",
        stage_id="N/A",
        gridpath_module="N/A",
        db_table="scenarios",
        severity="High",
        errors=errors
    )

    # Return True if all required subscenario_ids are valid (list is empty)
    return not bool(errors)


def reset_input_validation(conn, scenario_id):
    """
    Reset input validation: delete old input validation outputs and reset the
    input validation status.
    :param conn: database connection
    :param scenario_id: scenario_id
    :return: 
    """
    c = conn.cursor()

    sql = """
        DELETE FROM status_validation
        WHERE scenario_id = ?;
        """
    spin_on_database_lock(conn=conn, cursor=c, sql=sql, data=(scenario_id,),
                          many=False)

    sql = """
        UPDATE scenarios
        SET validation_status_id = 0
        WHERE scenario_id = ?;
        """
    spin_on_database_lock(conn=conn, cursor=c, sql=sql,  data=(scenario_id,),
                          many=False)


def update_validation_status(conn, scenario_id):
    """

    :param conn:
    :param scenario_id:
    :return:
    """
    c = conn.cursor()
    validations = c.execute(
        """SELECT scenario_id 
        FROM status_validation
        WHERE scenario_id = {}""".format(str(scenario_id))
    ).fetchall()

    if validations:
        status = 2
    else:
        status = 1

    sql = """
        UPDATE scenarios
        SET validation_status_id = ?
        WHERE scenario_id = ?;
        """
    spin_on_database_lock(conn=conn, cursor=c, sql=sql,
                          data=(status, scenario_id), many=False)


def parse_arguments(args):
    """
    :param args: the script arguments specified by the user
    :return: the parsed known argument values (<class 'argparse.Namespace'>
    Python object)

    Parse the known arguments.
    """
    parser = ArgumentParser(
        add_help=True,
        parents=[get_db_parser()]
    )

    # Add quiet flag which can suppress run output
    parser.add_argument("--quiet", default=False, action="store_true",
                        help="Don't print run output.")

    parsed_arguments = parser.parse_known_args(args=args)[0]

    return parsed_arguments


def determine_subscenarios_by_feature(conn):
    """

    :param conn:
    :return:
    """
    c = conn.cursor()

    feature_sc = c.execute(
        """SELECT feature, subscenario_id
        FROM mod_feature_subscenarios"""
    ).fetchall()
    feature_sc_dict = {}
    for f, sc in feature_sc:
        if f in feature_sc_dict:
            feature_sc_dict[f].append(sc.upper())
        else:
            feature_sc_dict[f] = [sc.upper()]
    return feature_sc_dict


# TODO: refactor this in capacity_types/__init__? (similar functions are
#   used in prm_types/operational_types etc.
def get_required_capacity_type_modules(conn, scenario_id):
    """
    :param c: database cursor
    :return: List of the required capacity type submodules

    Get the required capacity type submodules based on the database inputs
    for the specified scenario_id. Required modules are the unique set of
    generator capacity types in the scenario's portfolio. Get the list based
    on the project_operational_chars_scenario_id of the scenario_id.

    This list will be used to know for which capacity type submodules we
    should validate inputs, get inputs from database , or save results to
    database. It is also used to figure out which suscenario_ids are required
    inputs (e.g. cost inputs are required when there are new build resources)

    Note: once we have determined the dynamic components, this information
    will also be stored in the DynamicComponents class object.

    """
    c = conn.cursor()

    project_portfolio_scenario_id = c.execute(
        """SELECT project_portfolio_scenario_id 
        FROM scenarios 
        WHERE scenario_id = {}""".format(scenario_id)
    ).fetchone()[0]

    required_capacity_type_modules = [
        p[0] for p in c.execute(
            """SELECT DISTINCT capacity_type 
            FROM inputs_project_portfolios
            WHERE project_portfolio_scenario_id = {}""".format(
                project_portfolio_scenario_id
            )
        ).fetchall()
    ]

    return required_capacity_type_modules


def main(args=None):
    """

    :return:
    """

    # Retrieve scenario_id and/or name from args + "quiet" flag
    if args is None:
        args = sys.argv[1:]
    parsed_arguments = parse_arguments(args=args)

    if not parsed_arguments.quiet:
        print("Validating inputs...")

    db_path = parsed_arguments.database
    scenario_id_arg = parsed_arguments.scenario_id
    scenario_name_arg = parsed_arguments.scenario

    conn = connect_to_database(db_path=db_path,
                               detect_types=sqlite3.PARSE_DECLTYPES)
    c = conn.cursor()

    scenario_id, scenario_name = get_scenario_id_and_name(
        scenario_id_arg=scenario_id_arg,
        scenario_name_arg=scenario_name_arg,
        c=c,
        script="validate_inputs"
    )

    # Reset input validation status and results
    reset_input_validation(conn, scenario_id)

    # TODO: this is very similar to what's in get_scenario_inputs.py,
    #  so we should consolidate
    # Get scenario characteristics (features, scenario_id, subscenarios, subproblems)
    optional_features = OptionalFeatures(conn=conn, scenario_id=scenario_id)
    subscenarios = SubScenarios(conn=conn, scenario_id=scenario_id)
    subproblems = SubProblems(conn=conn, scenario_id=scenario_id)

    # Check whether subscenario_ids are valid
    is_valid = validate_subscenario_ids(scenario_id, subscenarios, optional_features, conn)

    # Only do the detailed input validation if all required subscenario_ids
    # are specified (otherwise will get errors when loading data)
    if is_valid:
        # Load modules for all requested features
        feature_list = optional_features.get_active_features()
        # If any subproblem's stage list is non-empty, we have stages, so set
        # the stages_flag to True to pass to determine_modules below
        # This tells the determine_modules function to include the
        # stages-related modules
        stages_flag = any([
            len(subproblems.SUBPROBLEM_STAGE_DICT[subp]) > 1 for subp in
            subproblems.SUBPROBLEM_STAGE_DICT.keys()
        ])
        modules_to_use = determine_modules(
            features=feature_list, multi_stage=stages_flag
        )
        loaded_modules = load_modules(modules_to_use=modules_to_use)

        # Read in inputs from db and validate inputs for loaded modules
        validate_inputs(subproblems, loaded_modules, scenario_id, subscenarios, conn)
    else:
        if not parsed_arguments.quiet:
            print("Invalid subscenario ID(s). Skipped detailed input validation.")

    # Update validation status:
    update_validation_status(conn, scenario_id)

    # Close the database connection explicitly
    conn.close()


if __name__ == "__main__":
    main()
