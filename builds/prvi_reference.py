"""An end-to-end build file that will take the NHD PRVI and turn it into a reference fabric"""

import argparse
import logging

from pydantic import ValidationError

from reference_builds.configs import PRVI
from reference_builds.local_runner import LocalRunner
from reference_builds.pipeline import download_nhd_data

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
        config = PRVI.from_yaml(args.config)
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

        print("Pipeline completed")
        print("=" * 60)
        for task_id, info in runner.results.items():
            status = "✓" if info["status"] == "success" else "✗"
            print(f"  {status} {task_id}: {info['status']}")
        print("=" * 60)

    return 0


if __name__ == "__main__":
    exit(main())
