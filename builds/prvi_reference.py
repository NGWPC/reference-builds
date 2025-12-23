"""An end-to-end build file that will take the v2.2 Hydrofabric for PRVI and turn it into the reference fabric"""

from reference_builds.logs import setup_logging

logger = setup_logging()

if __name__ == "__main__":
    logger.info("Here")
