"""Implements functions to load pre-defined variant signature mutations from GitHub."""

import requests
import yaml
import os
import re
from typing import Dict, List, Optional, Any, Tuple
from pathlib import Path
# pydantic 2 (2026-10-08; was pydantic 1: validator / parse_obj). The worker
# imports this file through /app_shared and needs pydantic 2 for covvfit.
from pydantic import BaseModel, field_validator
import logging

# Configure logging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

# Load configuration from config.yaml
def load_config() -> Dict[str, Any]:
    """Load configuration from config.yaml."""
    # Try multiple potential config locations
    config_paths = [
        Path("app/config.yaml"),      # When running from repository root
        Path("config.yaml"),          # When running from inside app directory
        Path("../config.yaml")        # When running from inside app/api directory
    ]
    
    for config_path in config_paths:
        if config_path.exists():
            logger.info(f"Loading config from {config_path.absolute()}")
            with open(config_path, 'r') as f:
                return yaml.safe_load(f)
    
    logger.warning("Config file not found. Using default values.")
    return {}

config = load_config()

# Get GitHub URL from config
variant_config = config.get("curated_variant_definitions", {})

# Default GitHub URL
DEFAULT_GITHUB_URL = "https://github.com/cbg-ethz/cowwid/tree/master/voc"
# Get the GitHub URL from config with fallback to default
GITHUB_URL = variant_config.get("github_url", DEFAULT_GITHUB_URL)



LOCAL_CACHE_DIR = Path("data/known_variants")

class Mutation(BaseModel): 
    """A mutation in nucleotide format."""
    position: int
    ref: str
    alt: str

    @field_validator('position')
    @classmethod
    def validate_position(cls, v):
        if v <= 0:
            raise ValueError("Position must be a positive integer")
        return v

    @field_validator('ref')
    @classmethod
    def validate_ref(cls, v):
        if not (len(v) == 0 or len(v) == 1):
            raise ValueError("Reference must be empty or a single nucleotide")
        if v and v not in "ACGTN":
            raise ValueError("Reference must be one of A, C, G, T, or N")
        return v

    @field_validator('alt')
    @classmethod
    def validate_alt(cls, v):
        if len(v) != 1:
            raise ValueError("Alternative must be a single character")
        if v not in "ACGTN-":
            raise ValueError("Alternative must be one of A, C, G, T, N, or - (deletion)")
        return v
    
    @classmethod
    def validate_mutation_string(cls, mutation_str: str) -> tuple[bool, str, dict]:
        """
        Validate a mutation string and return a tuple of:
        (is_valid, error_message, mutation_data)
        
        Valid formats:
        - REF+POS+ALT (e.g., C123T)
        - POS+ALT (e.g., 123T)
        - POS+- (e.g., 123-, for deletion)
        
        Args:
            mutation_str: Mutation string to validate
            
        Returns:
            Tuple of (is_valid, error_message, mutation_data)
            If valid, error_message will be empty and mutation_data will contain the parsed data
            If invalid, mutation_data will be empty
        """
        # Regex to capture: optional REF, mandatory POS, mandatory ALT
        # REF: A, C, G, T, N (optional)
        # POS: digits
        # ALT: A, C, G, T, N, -
        match = re.match(r"^([ACGTN]?)(\d+)([ACGTN-])$", mutation_str.upper())
        
        if not match:
            return False, f"Invalid format for '{mutation_str}'. Expected format is like 'C123T', '123-', or '123A'.", {}
        
        ref_char = match.group(1)
        pos_str = match.group(2)
        alt_char = match.group(3)
        
        try:
            position = int(pos_str)
            mutation_data = {
                "position": position,
                "ref": ref_char,
                "alt": alt_char
            }
            
            # Create an instance to validate
            cls(**mutation_data)
            
            return True, "", mutation_data
            
        except ValueError as e:
            if "position" in str(e):
                return False, f"Invalid position in '{mutation_str}': {str(e)}", {}
            return False, f"Invalid value in '{mutation_str}': {str(e)}", {}
            
        except Exception as e:
            return False, f"Error validating '{mutation_str}': {str(e)}", {}


# Pydantic models that match the YAML structure
class VariantInfo(BaseModel):
    """ A variant information type to load from YAML."""
    short: str
    pangolin: str
    nextstrain: str = ""

class VariantDefinition(BaseModel):
    """A variant definition type to load from YAML.
        VariantInfo: Contains variant information.
        mut: Dictionary of mutations with position as key and change as value.
        The mutation format is expected to be in the form of REF>ALT (e.g., C>T).
        The mutation position is 1-based.
    """
    variant: VariantInfo
    mut: Dict[int, str]

    @field_validator('mut', mode='before')
    @classmethod
    def convert_string_keys_to_int(cls, v):
        """Convert dictionary string keys to integers for mutation positions."""
        return {int(k): v for k, v in v.items()}

    def format_mutation(self, position: int, change: str) -> List[str]:
        """
        Format a mutation from position and change (e.g., 241: 'C>T') to the new format (e.g., 'C241T').
        
        If multiple mutations are concatenated:
          28881: GGG>AAC
          29734: '--'
        They are split into separate mutations:
            G28881A
            G28882A
            G28883C

        Deletions are ignored.
            
        Args:
            position: Position of the mutation
            change: Change in format 'REF>ALT' or 'REF>-' for deletions
            
        Returns:
            List of formatted mutation strings
        """
        mutations = []

        # handle deletion
        if all(c == '-' for c in change):
            ref = ""
            alt = "-"

            for i in range(len(change)):
                mutations.append(f"{ref}{position}{alt}")
            
            return mutations
                

        # Check if it's in the REF>ALT format
        if '>' in change:
            parts = change.split('>')
            if len(parts) != 2:
                logger.warning(f"Invalid mutation format at position {position}: {change}. Expected REF>ALT format. In variant {self.variant.pangolin}.")
                return mutations
                
            ref, alt = parts
                
            # Handle multiple mutations (e.g., GGG>AAC)
            if len(ref) > 1 and len(alt) > 1 and len(ref) == len(alt):
                for i in range(len(ref)):
                    if i < len(alt):  #  ensure index is valid
                        mutations.append(f"{ref[i]}{position + i}{alt[i]}")
            # Handle single ref multiple alt (e.g., G>AAC) or multiple ref single alt (e.g., GGG>A)
            elif len(ref) != len(alt):
                logger.warning(f"Complex mutation at position {position}: {change}. Skipping.")
            # Handle single mutation (e.g., C>T)
            elif ref and alt:
                mutations.append(f"{ref}{position}{alt}")
        else:
            logger.warning(f"Unexpected mutation format at position {position}: {change}. Expected REF>ALT format. In variant {self.variant.pangolin}.")
            
        return mutations
    

class Variant(BaseModel):
    name: str
    short_name: str
    nextstrain_name: str = ""
    signature_mutations: List[str]
    
    @classmethod
    def from_variant_definition(cls, variant_def: VariantDefinition) -> "Variant":
        """Create a Variant from a VariantDefinition."""
        # Convert mutations from {position: change} format to list of strings in new format
        mutations = []
        for pos, change in variant_def.mut.items():
            mutations.extend(variant_def.format_mutation(pos, change))
            
        return cls(
            name=variant_def.variant.pangolin,
            short_name=variant_def.variant.short,
            nextstrain_name=variant_def.variant.nextstrain,
            signature_mutations=mutations
        )

class VariantList(BaseModel):
    variants: List[Variant] = []
    
    def add_variant(self, variant: Variant) -> None:
        """Add a variant to the list."""
        self.variants.append(variant)
    
    def remove_variant(self, variant: Variant) -> None:
        """Remove a variant from the list."""
        self.variants.remove(variant)
    
    def get_variant_by_name(self, name: str) -> Optional[Variant]:
        """Get a variant by its name."""
        for variant in self.variants:
            if variant.name == name:
                return variant
        return None

def ensure_cache_dir() -> None:
    """Ensure that the local cache directory exists."""
    os.makedirs(LOCAL_CACHE_DIR, exist_ok=True)

def _parse_github_url(github_url: str) -> dict:
    """Parse GitHub URL into components needed for API and raw content access."""
    # Default values
    result = {
        "owner": "cbg-ethz",
        "repo": "cowwid",
        "branch": "master",
        "path": "voc"
    }
    
    try:
        # Extract owner/repo from github URL
        if "github.com" in github_url:
            # Extract owner and repo from URL
            url_parts = github_url.split("github.com/")[-1].split("/")
            
            if len(url_parts) >= 2:
                result["owner"] = url_parts[0]
                result["repo"] = url_parts[1]
            
            # Extract branch and path if present
            if len(url_parts) >= 4 and url_parts[2] in ["tree", "blob"]:
                result["branch"] = url_parts[3]
                if len(url_parts) >= 5:
                    result["path"] = "/".join(url_parts[4:])
    except Exception as e:
        logger.warning(f"Error parsing GitHub URL {github_url}: {e}. Using default values.")
        
    return result

def list_github_files() -> List[Dict[str, Any]]:
    """List all YAML files in the GitHub repository's VOC folder."""
    try:
        # Parse GitHub URL components
        github_components = _parse_github_url(GITHUB_URL)
        owner = github_components["owner"]
        repo = github_components["repo"]
        path = github_components["path"]
        
        # Construct API URL
        api_url = f"https://api.github.com/repos/{owner}/{repo}/contents/{path}"
        
        logger.info(f"Fetching files from GitHub API: {api_url}")
        response = requests.get(api_url)
        response.raise_for_status()
        
        # Filter for YAML files only
        files = [file for file in response.json() if file['name'].endswith('.yaml')]
        logger.info(f"Found {len(files)} YAML files in GitHub repository")
        return files
    except requests.RequestException as e:
        logger.error(f"Error fetching files from GitHub: {e}")
        return []

def download_yaml_file(file_name: str) -> Optional[Dict[str, Any]]:
    """Download and parse a YAML file from GitHub."""
    # Check if file exists in local cache
    cache_path = LOCAL_CACHE_DIR / file_name
    
    if cache_path.exists():
        logger.info(f"Loading {file_name} from local cache")
        try:
            with open(cache_path, 'r') as f:
                return yaml.safe_load(f)
        except Exception as e:
            logger.error(f"Error loading cached file {file_name}: {e}")
    
    # If not in cache or error occurred, download from GitHub
    try:
        # Parse GitHub URL components
        github_components = _parse_github_url(GITHUB_URL)
        owner = github_components["owner"]
        repo = github_components["repo"]
        branch = github_components["branch"]
        path = github_components["path"]
        
        # Construct raw content URL
        raw_url = f"https://raw.githubusercontent.com/{owner}/{repo}/{branch}/{path}/{file_name}"
        
        logger.info(f"Downloading {file_name} from {raw_url}")
        response = requests.get(raw_url)
        response.raise_for_status()
        
        # Parse the YAML content
        content = yaml.safe_load(response.text)
        
        # Save to cache
        ensure_cache_dir()
        with open(cache_path, 'w') as f:
            yaml.dump(content, f)
        
        return content
    except (requests.RequestException, yaml.YAMLError) as e:
        logger.error(f"Error downloading or parsing {file_name}: {e}")
        return None

def load_variant_definition(yaml_data: Dict[str, Any]) -> Optional[VariantDefinition]:
    """Convert YAML data to VariantDefinition."""
    try:
        # NOTE: Old manual Variant Definitions reature `shared` and `mut` keys,
        # which need to be merged correctly as they are part of the signature mutations.
        # Merge 'shared' mutations into 'mut' if present
        mut = yaml_data.get('mut', {})
        shared = yaml_data.get('shared', {})
        # Avoid overwriting any keys in mut with shared
        merged_mut = {**shared, **mut}
        yaml_data['mut'] = merged_mut
        return VariantDefinition.model_validate(yaml_data)
    except Exception as e:
        logger.error(f"Error parsing variant definition: {e}")
        return None

def get_all_variant_definitions() -> List[VariantDefinition]:
    """Get all variant definitions from GitHub, or load
       from local cache if GitHub is unavailable."""
    files = list_github_files()
    variant_defs = []
  
    if files is None or len(files) == 0:
            logger.error("No YAML files found in the GitHub repository.")
            # load all locally cached files
            logger.info("Loading all locally cached variant definitions")
            ensure_cache_dir()
            for file in LOCAL_CACHE_DIR.glob("*.yaml"):
                try:
                    with open(file, 'r') as f:
                        yaml_data = yaml.safe_load(f)
                        variant_def = load_variant_definition(yaml_data)
                        if variant_def:
                            variant_defs.append(variant_def)
                except Exception as e:
                    logger.error(f"Error loading cached file {file.name}: {e}") 

            return variant_defs
    else:
        logger.info(f"Found {len(files)} YAML files to process")
        for file in files:
            yaml_data = download_yaml_file(file['name'])
            if yaml_data:
                variant_def = load_variant_definition(yaml_data)
                if variant_def:
                    variant_defs.append(variant_def)
    
    logger.info(f"Successfully loaded {len(variant_defs)} variant definitions")
    return variant_defs

def get_variant_list() -> VariantList:
    """Get a VariantList containing all variants from GitHub."""
    variant_defs = get_all_variant_definitions()
    variant_list = VariantList()
    
    for variant_def in variant_defs:
        variant = Variant.from_variant_definition(variant_def)
        variant_list.add_variant(variant)
    
    return variant_list

def get_variant_by_name(name: str) -> Optional[Variant]:
    """Get a variant by its name (pangolin designation)."""
    variant_list = get_variant_list()
    return variant_list.get_variant_by_name(name)

def get_variant_names() -> List[str]:
    """Get a list of all available variant names."""
    variant_list = get_variant_list()
    return [variant.name for variant in variant_list.variants]

def validate_mutation_strings(mutations_str_list: List[str]) -> Tuple[bool, List[str], List[str]]:
    """
    Validate a list of mutation strings.
    
    Args:
        mutations_str_list: List of mutation strings to validate
        
    Returns:
        Tuple of (all_valid, valid_mutations, error_messages)
    """
    valid_mutations = []
    error_messages = []
    all_valid = True
    
    for mutation_str in mutations_str_list:
        is_valid, error_message, _ = Mutation.validate_mutation_string(mutation_str)
        if is_valid:
            valid_mutations.append(mutation_str)
        else:
            error_messages.append(error_message)
            all_valid = False
            
    return all_valid, valid_mutations, error_messages


if __name__ == "__main__":

    # Test the functions
    logger.info("Testing variant signature loading...")
    variant_list = get_variant_list()
    logger.info(f"Loaded {len(variant_list.variants)} variants")
    
    # Test the new mutation format
    logger.info("Examples of reformatted mutations:")
    for variant in variant_list.variants[:10]:  # Show first 3 variants as examples
        logger.info(f"Variant: {variant.name}, Signature mutations (first 5): {variant.signature_mutations[:5]}")

