FROM continuumio/miniconda3:latest

WORKDIR /app

# Copy requirements file
COPY environment.yml .

# Create the conda environment from environment.yml
RUN conda env create -f environment.yml && \
    conda run -n v-pipe-scout-worker python -c "import aiohttp, requests" || \
    (echo "Sanity check failed" && exit 1)

# Copy all worker files (scripts, configs, and data)
COPY . .

# Set environment variables
ENV PYTHONUNBUFFERED=1

# Run the Celery worker using the conda environment
CMD ["conda", "run", "--no-capture-output", "-n", "v-pipe-scout-worker", \
     "celery", "-A", "tasks", "worker", "--concurrency", "4", "--loglevel=info"]