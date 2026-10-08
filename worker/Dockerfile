FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy requirements file
COPY environment.yml .

# Create the conda environment from environment.yml
RUN conda env create -f environment.yml && \
    conda run -n v-pipe-scout-worker python -c "import aiohttp, requests" || \
    (echo "Sanity check failed" && exit 1)

# covvfit (growth advantages, after the deconvolution), from GitHub at a fixed
# commit: the PyPI release 0.3.1 predates covvfit's move to numpy 2 and pins
# numpy < 2, which has no build for Python 3.13. covvfit declares Python < 3.12
# but runs on 3.13: tested 2026-10-08 with this commit (jax 0.4.38, numpy 2.5),
# its 99 tests pass and `covvfit infer` runs end to end. --ignore-requires-python
# skips pip's Python check, hence the fixed commit; `covvfit check` fails the
# build if the install is broken. git comes from environment.yml.
ARG COVVFIT_COMMIT=80b12755116a59af49029393444c25eea3a46cd4
RUN conda run -n v-pipe-scout-worker pip install --no-cache-dir \
        --ignore-requires-python \
        "covvfit @ git+https://github.com/cbg-ethz/covvfit@${COVVFIT_COMMIT}" && \
    conda run -n v-pipe-scout-worker covvfit check

# Copy all worker files (scripts, configs, and data)
COPY . .

# Set environment variables
ENV PYTHONUNBUFFERED=1

# Run the Celery worker using the conda environment
CMD ["/bin/bash", "/app/start.sh"]