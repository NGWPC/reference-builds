
"""Local runner for building the NGWPC hydrofabric"""

from collections.abc import Callable
from datetime import datetime
from typing import Any, Self

from reference_builds.configs import PRVI
from reference_builds.task_instance import TaskInstance
from reference_builds.logs import setup_logging

class LocalRunner:
    """Execute pipeline tasks locally with Airflow-like interface.

    Parameters
    ----------
    config : HFConfig
        Pipeline configuration containing build settings and parameters.
    run_id : str or None, default=None
        Unique identifier for this pipeline run. If None, generated from
        current timestamp in format 'YYYYMMDD_HHMMSS'.

    Attributes
    ----------
    config : HFConfig
        The pipeline configuration.
    run_id : str
        Unique identifier for this run.
    ti : TaskInstance
        TaskInstance for XCom operations.
    results : dict[str, dict[str, Any]]
        Execution results for each task, keyed by task_id.
    """

    def __init__(
        self,
        config: PRVI,
        run_id: str | None = None,
    ) -> None:
        """Initialize the LocalRunner.

        Parameters
        ----------
        config : HFConfig
            Pipeline configuration.
        run_id : str or None, default=None
            Optional run identifier. Auto-generated if not provided.
        """
        self.config: PRVI = config
        self.run_id: str = run_id or datetime.now().strftime("%Y%m%d_%H%M%S")
        self.ti: TaskInstance = TaskInstance()
        self.results: dict[str, dict[str, Any]] = {}
        self.logger = setup_logging()

    def cleanup(self) -> None:
        """Clean up resources"""
        self.logger.info("runner: Closing processes")

    def __enter__(self: Self) -> Self:
        """Context manager entry."""
        return self

    def __exit__(self: Self, *args: str, **kwargs: str) -> None:
        """Context manager exit - ensures cleanup."""
        self.cleanup()

    def run_task(
        self,
        task_id: str,
        python_callable: Callable[..., Any],
        op_kwargs: dict[str, Any] | None = None,
    ) -> Any:
        """Execute a single task.

        Parameters
        ----------
        task_id : str
            Unique identifier for this task. Used in XCom keys and result tracking.
        python_callable : Callable[..., Any]
            The function to execute. Must accept **kwargs to receive context.
        op_kwargs : dict[str, Any] or None, default=None
            Additional keyword arguments to pass to the callable.

        Returns
        -------
        Any
            The return value from the callable.
        """
        self.logger.info(f"Running task: {task_id}")

        context: dict[str, Any] = {
            "ti": self.ti,
            "task_id": task_id,
            "run_id": self.run_id,
            "ds": datetime.now().strftime("%Y-%m-%d"),
            "execution_date": datetime.now(),
            "config": self.config,
        }

        kwargs = {**(op_kwargs or {}), **context}

        result = python_callable(**kwargs)

        for k, v in result.items():
            self.ti.xcom_push(f"{task_id}.{k}", v)
        self.results[task_id] = {"status": "success", "result": result}

        self.logger.info(f"✓ Task {task_id} completed")
        return result

    def get_result(self, task_id: str) -> dict[str, Any]:
        """Retrieve execution results for a specific task.

        Parameters
        ----------
        task_id : str
            The identifier of the task to get results for.

        Returns
        -------
        dict[str, Any] or None
            Dictionary containing 'status' and either 'result' (on success)
            or 'error' (on failure). Returns None if task_id not found.
        """
        result = self.results.get(task_id)
        if result is None:
            raise ValueError("Cannot find result from task")
        return result