#!/usr/bin/env python
# Copyright 2020 Blue Marble Analytics LLC. All rights reserved.

import os
import pandas as pd


def get_availabilities_as_df():
    all_profiles_df = pd.read_csv(
        "./project/resolve_maintenance.csv"
    )

    return all_profiles_df


def create_gp_profile(resolve_df, project):
    print(project)

    # Make dataframe for the GridPath profile with headers
    gp_profile_df = pd.DataFrame(
        columns=["stage_id", "timepoint", "availability_derate"]
    )

    row = 0
    # Make profiles for every year until 2050, for each of the 37 RESOLVE days
    for year in range(2020, 2051):
        print("... {}".format(year))

        for resolve_day in range(1, 38):
            derate = resolve_df[
                (resolve_df["day"] == resolve_day) &
                (resolve_df["resource"] == project)
                ].iloc[0]["value"]
            for hour in range(1, 25):
                tmp = year * 10**4 + resolve_day * 10**2 + hour
                gp_profile_df.loc[row] = [1, tmp, derate]
                row += 1

    return gp_profile_df


def save_project_av_to_csv(prj, profile):
    profile_filename = os.path.join(
        "./project/availability/exogenous/",
        "{}-1-resolve_availability.csv".format(prj)
    )
    profile.to_csv(profile_filename, index=False)


if __name__ == "__main__":
    availabilities_df = get_availabilities_as_df()
    projects = pd.unique(availabilities_df["resource"])

    for p in projects:
        av_df = create_gp_profile(resolve_df=availabilities_df, project=p)
        save_project_av_to_csv(prj=p, profile=av_df)


