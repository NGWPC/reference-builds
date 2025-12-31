"""A file to host all Hydrofabric Schemas"""

from pathlib import Path
from typing import Self

import yaml
from pydantic import BaseModel, Field
from pyprojroot import here

from reference_builds import __version__


class ReferenceConfig(BaseModel):
    """Configs for building the ReferenceConfig reference"""

    output_dir: Path = Field(
        default=here() / "data/",
        description="The directory for output files to be saved from Hydrofabric builds",
    )

    domain: str = Field(description="The domain used for the building your reference")

    input_file_regex: str = Field(
        description="input file to be converted into a reference product",
    )

    crs: str = Field(
        default="EPSG:4326",
        description="Coordinate Reference System for the domain reference builds. Defaults to https://epsg.io/4326",
    )

    vpu_id: str = Field(description="The VPUID for the domain")

    write_gpkg: bool = Field(
        default=True, description="Writes a geopackage in addition to parquet files for output"
    )

    permitted_fcodes: list[str] = Field(
        default_factory=lambda: [
            "Stream/River: Hydrographic Category = Intermittent",
            "Artificial Path",
            "Connector",
            "Stream/River: Hydrographic Category = Perennial",
            "Canal/Ditch",
            "Canal Ditch: Canal Ditch Type = Stormwater",
        ],
        description="The permitted fcode descriptions for the reference",
    )

    output_reference_gpkg_path: Path = Field(
        default_factory=lambda data: data["output_dir"] / f"{data['domain']}_{__version__}_reference.gpkg",
        description="Save directory for the domain's reference (in .gpkg form)",
    )

    output_reference_divides_path: Path = Field(
        default_factory=lambda data: data["output_dir"]
        / f"{data['domain']}_{__version__}_reference_divides.parquet",
        description="Save directory for the domain's reference divides",
    )

    output_reference_flowpaths_path: Path = Field(
        default_factory=lambda data: data["output_dir"]
        / f"{data['domain']}_{__version__}_reference_flowpaths.parquet",
        description="Save directory for the domain's reference flowpaths",
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
