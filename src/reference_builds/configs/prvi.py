"""A file to host all Hydrofabric Schemas"""

from pathlib import Path
from typing import Self

import yaml
from pydantic import BaseModel, Field
from pyprojroot import here

from reference_builds import __version__


class PRVI(BaseModel):
    """Configs for building the PRVI reference"""

    output_dir: Path = Field(
        default=here() / "data/",
        description="The directory for output files to be saved from Hydrofabric builds",
    )

    input_file_regex: Path = Field(
        default_factory=lambda data: data["output_dir"] / "PRVI/NHDPLUS_H_2101_HU4_GPKG.gpkg",
        description="input file to be converted into a reference product",
    )

    crs: str = Field(
        default="EPSG:6566",
        description="Coordinate Reference System for the PRVI reference builds. Defaults to https://epsg.io/6566",
    )

    output_reference_divides_path: Path = Field(
        default_factory=lambda data: data["output_dir"] / f"prvi_{__version__}_reference_divides.parquet",
        description="Save directory for the PRVI reference divides",
    )

    output_reference_flowpaths_path: Path = Field(
        default_factory=lambda data: data["output_dir"] / f"prvi_{__version__}_reference_flowpaths.parquet",
        description="Save directory for the PRVI reference flowpaths",
    )

    @classmethod
    def from_yaml(cls, path: str | Path) -> Self:
        """An internal method to read a config from a YAML file

        Parameters
        ----------
        path : str | Path
            The path to the provided YAML file

        Returns
        -------
        HFConfig
            A configuration object validated
        """
        with open(path) as f:
            data = yaml.safe_load(f)

        return cls(**data)
