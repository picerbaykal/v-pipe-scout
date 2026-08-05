"""Centralized configuration loading utilities."""

import yaml
import pathlib
import logging
from typing import Tuple, Dict, Any

logger = logging.getLogger(__name__)


def load_config() -> Dict[str, Any]:
    """
    Load configuration from config.yaml.
    
    Returns:
        Dictionary containing the configuration
    """
    config_path = pathlib.Path(__file__).parent.parent / "config.yaml"
    try:
        with open(config_path, 'r') as file:
            config = yaml.safe_load(file)
            logger.info(f"Loaded configuration from {config_path}")
            return config
    except FileNotFoundError:
        logger.error(f"Configuration file not found at {config_path}")
        raise
    except yaml.YAMLError as e:
        logger.error(f"Error parsing YAML configuration: {e}")
        raise


def get_api_urls() -> Tuple[str, str]:
    """
    Get API URLs from configuration.
    
    Returns:
        Tuple of (wiseloculus_url, covspectrum_url)
    """
    config = load_config()
    
    wise_url = config.get('server', {}).get('lapis_address', 'http://default_ip:8000')
    covspectrum_url = config.get('server', {}).get('cov_spectrum_api', 'https://lapis.cov-spectrum.org')
    
    logger.info(f"API URLs - WiseLoculus: {wise_url}, CovSpectrum: {covspectrum_url}")
    
    return wise_url, covspectrum_url


def get_wiseloculus_url() -> str:
    """
    Get WiseLoculus API URL from configuration.
    
    Returns:
        WiseLoculus API URL
    """
    wise_url, _ = get_api_urls()
    return wise_url


def get_covspectrum_url() -> str:
    """
    Get CovSpectrum API URL from configuration.
    
    Returns:
        CovSpectrum API URL
    """
    _, covspectrum_url = get_api_urls()
    return covspectrum_url

# ── Cooc-specific configuration ──────────────────────────────────────────

_COOC_CONFIG_CACHE: Dict[str, Any] | None = None


def load_cooc_config() -> Dict[str, Any]:
    """
    Load co-occurrence analysis configuration from config/cooc_config.yaml.

    Cached after first read.

    Returns:
        Dictionary containing the cooc config.
    """
    global _COOC_CONFIG_CACHE
    if _COOC_CONFIG_CACHE is not None:
        return _COOC_CONFIG_CACHE

    config_path = pathlib.Path(__file__).parent.parent / "config" / "cooc_config.yaml"
    try:
        with open(config_path, "r") as file:
            _COOC_CONFIG_CACHE = yaml.safe_load(file)
            logger.info(f"Loaded cooc config from {config_path}")
            return _COOC_CONFIG_CACHE
    except FileNotFoundError:
        logger.error(f"Cooc config file not found at {config_path}")
        raise
    except yaml.YAMLError as e:
        logger.error(f"Error parsing cooc YAML config: {e}")
        raise


def get_cooc_setting(key_path: str, default: Any = None) -> Any:
    """
    Look up a nested key in cooc config, e.g. "query.max_position_distance_bp".

    Args:
        key_path: Dot-separated path to the setting.
        default: Value to return if the key is missing.

    Returns:
        The setting value, or default if not found.
    """
    config = load_cooc_config()
    val: Any = config
    for part in key_path.split("."):
        if not isinstance(val, dict) or part not in val:
            return default
        val = val[part]
    return val
