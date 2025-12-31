"""An end-to-end build file that will take the NHD ReferenceConfig and turn it into a reference fabric"""

import argparse
import logging

from pydantic import ValidationError

from reference_builds.configs import ReferenceConfig
from reference_builds.local_runner import LocalRunner
from reference_builds.pipeline import build_graphs, build_reference, download_nhd_data, write_reference

logger = logging.getLogger(__name__)


def main() -> int:
    """Main entry point for the hydrofabric-build pipeline CLI.

    Returns
    -------
    int
        Exit code: 0 for success, 1 for failure.
    """
    parser = argparse.ArgumentParser(description="A local runner for hydrofabric data processing")
    parser.add_argument("--config", required=False, help="Config file")
    args = parser.parse_args()

    try:
        config = ReferenceConfig.from_yaml(args.config)
    except ValidationError as e:
        print("Configuration validation failed:")
        for error in e.errors():
            print(f"  {error['loc']}: {error['msg']}")
        return 1
    except FileNotFoundError:
        logger.error(f"Config file not found: {args.config}")
        return 1
    except TypeError as e:
        logger.error("Config file not specified.")
        raise TypeError("Config file not specified.") from e

    with LocalRunner(config) as runner:
        runner.run_task(task_id="download", python_callable=download_nhd_data, op_kwargs={})
        runner.run_task(task_id="build_graphs", python_callable=build_graphs, op_kwargs={})
        runner.run_task(task_id="build_reference", python_callable=build_reference, op_kwargs={})
        runner.run_task(task_id="write_reference", python_callable=write_reference, op_kwargs={})

        print("Pipeline completed")
        print("=" * 60)
        for task_id, info in runner.results.items():
            status = "✓" if info["status"] == "success" else "✗"
            print(f"  {status} {task_id}: {info['status']}")
        print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
