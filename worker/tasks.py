
from celery import Celery
import os
import json
import time
import redis
import pickle
import base64
import pandas as pd
from deconvolve import devconvolve
from covvfit import run_covvfit_inference

# Initialize Celery
app = Celery(
    'tasks',
    broker=os.environ.get('CELERY_BROKER_URL', 'redis://redis:6379/0'),
    backend=os.environ.get('CELERY_RESULT_BACKEND', 'redis://redis:6379/0')
)

# Initialize Redis client for storing progress updates
redis_client = redis.Redis(
    host=os.environ.get('REDIS_HOST', 'redis'),
    port=int(os.environ.get('REDIS_PORT', 6379)),
    password=os.environ.get('REDIS_PASSWORD', 'defaultpassword123'),
    db=0
)

# EXAMPLE TASK: that progressively updates its status in Redis
@app.task(bind=True)
def long_running_task(self, n_iterations, sleep_time):
    """
    A simple long-running task that simulates work by sleeping.
    Updates progress in Redis to allow the frontend to track status.
    """
    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"
    
    # Initialize result container
    result = {
        "iterations_completed": 0,
        "total_iterations": n_iterations,
        "results": []
    }
    
    # Process each iteration
    for i in range(n_iterations):
        # Simulate work
        time.sleep(sleep_time)
        
        # Calculate some dummy result
        iteration_result = {
            "iteration": i + 1,
            "timestamp": time.time(),
            "value": (i + 1) * sleep_time
        }
        
        # Add to results
        result["results"].append(iteration_result)
        result["iterations_completed"] = i + 1
        
        # Update progress in Redis
        progress_data = {
            "current": i + 1,
            "total": n_iterations,
            "status": f"Processing iteration {i + 1}/{n_iterations}",
            "partial_results": result["results"]
        }
        
        redis_client.set(
            progress_key,
            json.dumps(progress_data),
            ex=3600  # Expire after 1 hour
        )
    
    # Task completed
    progress_data = {
        "current": n_iterations,
        "total": n_iterations,
        "status": "Completed",
        "partial_results": result["results"]
    }
    
    redis_client.set(
        progress_key,
        json.dumps(progress_data),
        ex=3600  # Expire after 1 hour
    )
    
    return result

@app.task(bind=True)
def run_deconvolve(self, mutation_counts_df, mutation_variant_matrix_df, 
                   bootstraps=None, bandwidth=None, regressor=None, 
                   regressor_params=None, deconv_params=None, locationName=None):
    """
    A task that runs the deconvolve function with progress tracking.
    
    Args:
        mutation_counts_df (pd.DataFrame): DataFrame containing mutation counts data (required)
        mutation_variant_matrix_df (pd.DataFrame): DataFrame containing mutation variant matrix data (required)
        bootstraps (int, optional): Number of bootstrap iterations
        bandwidth (int, optional): Bandwidth parameter for kernel
        regressor (str, optional): Regressor type
        regressor_params (dict, optional): Parameters for the regressor
        deconv_params (dict, optional): Parameters for deconvolution
        locationName (str, optional): Name of the location for proper result structuring
    """
    
    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"
    
    # Initialize progress tracking
    progress_data = {
        "current": 0,
        "total": 5,  # We'll track progress in 5 stages
        "status": "Preparing input data",
        "partial_results": None
    }
    
    redis_client.set(
        progress_key,
        json.dumps(progress_data),
        ex=3600  # Expire after 1 hour
    )
    
    try:
        # Update progress
        progress_data["current"] = 1
        progress_data["status"] = f"Preparing deconvolution (bootstraps={bootstraps if bootstraps is not None else 'default'})"
        redis_client.set(progress_key, json.dumps(progress_data), ex=3600)
        
        # Function to update progress
        def update_progress(stage, message):
            progress_data["current"] = stage
            progress_data["status"] = message
            redis_client.set(progress_key, json.dumps(progress_data), ex=3600)
        
        # Convert serialized DataFrames back to pandas DataFrames if needed
        try:
            # Add debug info about the input types
            update_progress(1.5, f"Input types: mutation_counts_df: {type(mutation_counts_df)}, mutation_variant_matrix_df: {type(mutation_variant_matrix_df)}")
            
            # Check if inputs are already DataFrames or need to be deserialized
            if isinstance(mutation_counts_df, pd.DataFrame) and isinstance(mutation_variant_matrix_df, pd.DataFrame):
                update_progress(2, "Inputs are already DataFrames, no parsing needed")
            else:
                # Try to deserialize if needed
                try:
                    # If they're base64 encoded pickle strings
                    if isinstance(mutation_counts_df, str):
                        try:
                            mutation_counts_df = pickle.loads(base64.b64decode(mutation_counts_df))
                            update_progress(2, f"Successfully unpickled counts DataFrame, shape: {mutation_counts_df.shape}")
                        except:
                            update_progress(2, "Failed to unpickle counts DataFrame as base64")
                    
                    if isinstance(mutation_variant_matrix_df, str):
                        try:
                            mutation_variant_matrix_df = pickle.loads(base64.b64decode(mutation_variant_matrix_df))
                            update_progress(2, f"Successfully unpickled matrix DataFrame, shape: {mutation_variant_matrix_df.shape}")
                        except:
                            update_progress(2, "Failed to unpickled matrix DataFrame as base64")
                    
                except Exception as e:
                    update_progress(2, f"Error parsing DataFrames: {str(e)}")
                    raise ValueError(f"Failed to deserialize DataFrames: {str(e)}")
        except Exception as e:
            update_progress(2, f"Error processing DataFrames: {str(e)}")
            raise ValueError(f"Failed to process DataFrames: {str(e)}")
        
        # Create kwargs dict with required parameters and optional parameters if provided
        kwargs = {
            'mutation_counts_df': mutation_counts_df,
            'mutation_variant_matrix_df': mutation_variant_matrix_df
        }
        
        # Add optional parameters only if they're not None (exclude locationName)
        if bootstraps is not None:
            kwargs['bootstraps'] = bootstraps
        if bandwidth is not None:
            kwargs['bandwidth'] = bandwidth
        if regressor is not None:
            kwargs['regressor'] = regressor
        if regressor_params is not None:
            kwargs['regressor_params'] = regressor_params
        if deconv_params is not None:
            kwargs['deconv_params'] = deconv_params
            
        # Update progress before running deconvolution
        update_progress(3, "Running deconvolution algorithm")
        
        # Run the deconvolution with only the provided parameters
        deconvolved_data = devconvolve(**kwargs)
        
        # Update progress after deconvolution is complete
        update_progress(4, "Processing results")
        
        # If locationName is provided, restructure the result to use the actual location name
        if locationName and isinstance(deconvolved_data, dict) and "location" in deconvolved_data:
            # Replace the generic "location" key with the actual location name
            location_data = deconvolved_data.pop("location")  # Remove and get the data
            deconvolved_data[locationName] = location_data   # Add with proper name
            update_progress(4.5, f"Restructured result for location: {locationName}")
        
        # Stage 5: Finalize results
        progress_data["current"] = 5
        progress_data["total"] = 5
        progress_data["status"] = "Completed"
        progress_data["partial_results"] = {"summary": "Deconvolution completed successfully"}
        redis_client.set(progress_key, json.dumps(progress_data), ex=3600)
        
        return deconvolved_data
    except Exception as e:
        # If there's an error, report it
        error_message = str(e)
        progress_data["status"] = f"Error: {error_message}"
        redis_client.set(progress_key, json.dumps(progress_data), ex=3600)
        raise

@app.task(bind=True)
def run_covvfit(self, location_data, matrix_df, max_days=180, horizon=90):
    """
   A task that runs covvfit fitness inference with progress tracking.

   Runs no-smoothing deconvolution per location, pools the results, and runs
   covvfit inference, returning the generated figures.

   Args:
       location_data (dict): {location_name: counts_df} — each counts_df is a
           base64-encoded pickle string of the mutation counts DataFrame.
       matrix_df (str): base64-encoded pickle string of the mutation-variant matrix.
       max_days (int): number of past days to restrict the analysis to.
       horizon (int): number of future days to forecast.
   """

    task_id = self.request.id
    progress_key = f"task_progress:{task_id}"

    def update_progress(current, total, message):
        redis_client.set(
            progress_key,
            json.dumps({"current": current, "total": total, "status": message}),
            ex=3600
        )

    # deserialize inputs and run
    try:
        update_progress(1, 3, "Deserializing inputs")

        # matrix_df arrives as a base64-encoded pickle string (same as run_deconvolve)
        if isinstance(matrix_df, str):
            matrix_df = pickle.loads(base64.b64decode(matrix_df))

        # location_data is {location: counts_df}, each counts_df a pickle string
        deserialized_location_data = {}
        for location, counts_df in location_data.items():
            if isinstance(counts_df, str):
                counts_df = pickle.loads(base64.b64decode(counts_df))
            deserialized_location_data[location] = counts_df

        update_progress(2, 3, "Running no-smoothing deconvolution and covvfit inference")

        result = run_covvfit_inference(
            location_data=deserialized_location_data,
            matrix_df=matrix_df,
            max_days=max_days,
            horizon=horizon,
        )

        update_progress(3, 3, "Completed")

        return result

    except Exception as e:
        update_progress(0, 3, f"Error: {str(e)}")
        raise




