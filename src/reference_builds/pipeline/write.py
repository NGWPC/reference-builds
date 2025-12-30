"""Contains all code for downloading hydrofabric data"""

import logging
from typing import Any, cast

from reference_builds.configs import PRVI
from reference_builds.task_instance import TaskInstance

logger = logging.getLogger(__name__)


def write_nhd_data(**context: dict[str, Any]) -> dict[str, Any]:
    """Opens local / downloads for the reference-build process

    Parameters
    ----------
    **context : dict
        Airflow-compatible context containing:
        - ti : TaskInstance for XCom operations
        - config : HFConfig with pipeline configuration
        - task_id : str identifier for this task
        - run_id : str identifier for this pipeline run
        - ds : str execution date
        - execution_date : datetime object

    Returns
    -------
    dict[str, Any]
        The reference flowpath and divides references in memory
    """
    cfg = cast(PRVI, context["config"])
    ti = cast(TaskInstance, context["ti"])
    cfg.output_reference_flowpaths_path.unlink(missing_ok=True)  # deletes files that exist with the same name
    cfg.output_reference_divides_path.unlink(missing_ok=True)  # deletes files that exist with the same name

    final_flowpaths = ti.xcom_pull(task_id="build_reference", key="reference_flowpaths")
    final_divides = ti.xcom_pull(task_id="build_reference", key="reference_divides")

    if cfg.write_gpkg:
        cfg.output_reference_gpkg_path.unlink(missing_ok=True)
        final_flowpaths.to_file(cfg.output_reference_gpkg_path, layer="reference_flowpaths", driver="GPKG")
        final_divides.to_file(cfg.output_reference_gpkg_path, layer="reference_divides", driver="GPKG")
        logger.info(f"write_nhd_data task: wrote geopackage layers to {cfg.output_reference_gpkg_path}")

    final_flowpaths.to_parquet(cfg.output_reference_flowpaths_path)
    final_flowpaths.to_parquet(cfg.output_reference_divides_path)
    return {}
