import os
import pandas as pd
from mylogger import get_logger

logger = get_logger(__name__)


def get_experiment_info():
    # Sets path to ~/MATseq
    path = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

    conditions_pd = pd.read_excel(
        os.path.join(path, "experiment/exp_settings.xlsx"), index_col=("field")
    ).fillna("")

    try:
        experiment_information = conditions_pd.iloc[
            : conditions_pd.index.get_loc("sample_id")
        ].to_dict()["values"]
        if experiment_information["project_number"] == "":
            raise ValueError("Can not find project number")
    except ValueError as e:
        logger.error("Can not find project number, provided by a sequencing company\n ")
        print(e)
        exit()

    sample_dict = conditions_pd.iloc[
        conditions_pd.index.get_loc("sample_id") + 1 :
    ].to_dict()["values"]

    logger.info(
        "Project "
        + str(experiment_information["project_number"])
        + " information loaded"
    )

    return (experiment_information, sample_dict)
