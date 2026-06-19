# V-Pipe Scout: Rapid Interactive Viral Variant Detection 

![POC](https://img.shields.io/badge/status-POC-yellow)
![Python](https://img.shields.io/badge/python-3.13-blue)
![Streamlit](https://img.shields.io/badge/streamlit-1.51.0-brightgreen)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)


## Overview

Recognizing and quantifying viral variants from wastewater requires expert human judgment in the final steps.
V-Pipe Scout allows for rapid exploration of wastewater viral sequences down to the single read level. 

Its aim: Discover novel viral threats a few weeks earlier than traditional methods.

This Proof-of-Concept is set up for SARS-CoV-2, yet is built to be virus-agnostic and will be expanded to RSV and Influenza soon.

This is an effort of the V-Pipe team.
For more information about V-Pipe, visit the [V-Pipe website](https://cbg-ethz.github.io/V-pipe/).

<div align="center">
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="./app/images/index/POC_DeployForInternal_inverted.png">
    <img src="./app/images/index/POC_DeployForInternal.png" alt="Fast Query Visualization" width="800"/>
  </picture>
  <p><em>Technical architecture for real-time visualization of viral sequencing data & rapid on-demand analysis</em></p>
</div>

Specifically, V-Pipe Scout enables:
- **Exploration of mutations at the read level**  
    - For known resistance mutations  
    - Guided by smart filters and variant signatures
- **Composition of variant signatures for abundance estimates**  
    - Leveraging clinical sequence databases (e.g., [CovSpectrum](https://cov-spectrum.org/))  
    - Using curated variant signatures
- **Variant fitness inference**  
    - Estimates relative fitness advantages with [covvfit](https://github.com/cbg-ethz/covvfit)  
    - Pooling data across selected locations, with forecasts of future variant dynamics
- **URL-based session sharing**  
    - Share analysis configurations via URLs
    - Collaborate by sharing specific page setups
    - Bookmark and resume analysis sessions

Further, we will implement:
- On-demand variant abundance estimates by [Lollipop](https://github.com/cbg-ethz/LolliPop)

V-Pipe Scout brings together:
- [V-pipe](https://github.com/cbg-ethz/V-pipe) - our prime Wastewater Viral Analysis Pipeline, see [publication](https://www.biorxiv.org/content/10.1101/2023.10.16.562462v1.full).
- [covvfit](https://github.com/cbg-ethz/covvfit) - fitness estimates of SARS-CoV-2 variants from variant abundance data, see [publication](https://doi.org/10.1016/j.watres.2026.126018)
- [GenSpectrum](https://genspectrum.org/) - in particular the novel fast database for genomic sequences [LAPIS-SILO](https://github.com/GenSpectrum/LAPIS-SILO), see [publication](https://bmcbioinformatics.biomedcentral.com/articles/10.1186/s12859-023-05364-3)


This application relies on two other repos as connecting infrastructure:
- [WisePulse](https://github.com/cbg-ethz/WisePulse) - to pre-process and run the SILO database, powering read-level queries
- [sr2silo](https://github.com/cbg-ethz/sr2silo) - large scale data-wrangler of nucleotide alignments, to amino-acids and SILO input format


## Deployment

The current deployment of this project can be accessed at [dev.vpipe.ethz.ch](http://dev.vpipe.ethz.ch).
_Only accessible within ETH Zürich Networks._

### Installation

1. Clone the repository:
    ```sh
    git clone https://github.com/cbg-ethz/v-pipe-scout.git
    cd v-pipe-scout
    ```

2. **Setup environment:**
    ```sh
    ./setup.sh  # Creates .env with secure Redis password (single source of truth)
    ```

3. **Configure LAPIS connection** in `app/config.yaml`:
    ```yaml
    server:
      lapis_address: "http://host.docker.internal:8083"  # For local LAPIS
    ```

4. **Run the application:**
    ```sh
    docker compose up --build
    ```

### Automatic Deployment

For production deployments on VMs or servers, you can set up automatic deployment to eliminate the need for manual updates. See [DEPLOYMENT.md](DEPLOYMENT.md) for detailed instructions on:

- Setting up automatic deployment with cron jobs
- Monitoring and logging deployment activities  
- Configuring rollback mechanisms
- Troubleshooting deployment issues

### Architecture

- **Streamlit Frontend**: Interactive web interface
- **Celery Worker**: Background task processing  
- **Redis**: Message broker (password-protected, internal only)


## Project Origin

This project was initiated as part of a hackathon project at the [BioHackathon Europe 2024](https://biohackathon-europe.org/).

## CI/CD

**Testing**: Automated tests for frontend, worker, and full Docker Compose stack run on every push/PR.

**Deployment**:
- **Development**: [`auto-deploy.sh`](scripts/auto-deploy.sh) runs every 5 minutes via cron → [dev.vpipe.ethz.ch](http://dev.vpipe.ethz.ch) (ETH network only)
- **Production**: [`deploy.yml`](.github/workflows/deploy.yml) triggers on GitHub releases for production deployment


## Contributing

Contributions are welcome! Please fork the repository and submit a pull request with your changes. For major changes, please open an issue first to discuss what you would like to change.

## License

This project is licensed under the MIT License. See the [LICENSE](LICENSE) file for details.
