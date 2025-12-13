/Users/kevintu/Library/Caches/org.R-project.R/R/reticulate/uv/cache/archive-v0/VH8Np4G5fHN_6YqbSBEn9/bin/python3 -m pip install requests anndata pandas numpy
import requests
import anndata as ad
import os 
import pandas as pd 
import numpy as np 
from scipy.sparse import issparse # To check if matrix is sparse

def load_h5ad_data(local_filename):
    """
    Loads an H5AD file into an AnnData object.
    """
    if not os.path.exists(local_filename):
        print(f"Error: File not found at {local_filename}. Cannot load.")
        return None
    
    print(f"Attempting to load H5AD file: {local_filename}")
    try:
        adata = ad.read_h5ad(local_filename)
        print("H5AD file loaded successfully into an AnnData object.")
        return adata
    except Exception as e:
        print(f"Error loading H5AD file: {e}")
        return None

def save_dataframe_to_text(df, path, component_name):
    """Saves a pandas DataFrame to a tab-separated text file."""
    try:
        print(f"Exporting {component_name} to {path} ({df.shape[0]} rows, {df.shape[1]} cols)...")
        df.to_csv(path, sep='\t', index=True) # Keep index for obs, var, etc.
        print(f"{component_name} saved successfully.")
    except Exception as e:
        print(f"Error exporting {component_name} to {path}: {e}")

def save_matrix_to_text(matrix_data, index_names, column_names, path, component_name):
    """
    Converts matrix data (potentially sparse) to a DataFrame and saves it to a text file.
    Matrix is expected to be genes (rows) x cells (columns).
    """
    print(f"Preparing {component_name} for export...")
    if issparse(matrix_data):
        print(f"{component_name} is sparse with shape: {matrix_data.shape}")
        # Convert to DataFrame, handling sparse data efficiently
        df = pd.DataFrame.sparse.from_spmatrix(
            data=matrix_data, 
            index=index_names, 
            columns=column_names
        )
    else: # If data is already dense (e.g., numpy array)
        print(f"{component_name} is dense with shape: {matrix_data.shape}")
        df = pd.DataFrame(
            data=matrix_data, 
            index=index_names, 
            columns=column_names
        )
    
    # Reset index to make the current index (gene names) a column named 'V1'
    df.reset_index(inplace=True)
    df.rename(columns={'index': 'V1'}, inplace=True)
    save_dataframe_to_text(df, path, component_name)


# --- Configuration ---
target_directory = os.path.expanduser("~/Library/CloudStorage/OneDrive-Personal/Research/Duke Amgen/Glioblastoma TME")
local_filename_only = "GBMap_data.h5ad"
local_h5ad_full_path = os.path.join(target_directory, local_filename_only)

# --- Main execution logic ---
print(f"Attempting to import data from: {local_h5ad_full_path}")
adata_object = load_h5ad_data(local_h5ad_full_path)

if adata_object:
    print("\n--- AnnData Object Summary ---")
    print(adata_object)
    
    # Ensure output directory exists
    export_dir = os.path.join(target_directory, "adata_export")
    os.makedirs(export_dir, exist_ok=True)
    print(f"\nExporting components to directory: {export_dir}")

    # 1. Export adata_object.obs (cell metadata)
    if not adata_object.obs.empty:
        obs_path = os.path.join(export_dir, "obs_data.txt")
        save_dataframe_to_text(adata_object.obs, obs_path, "obs (cell_metadata)")
    else:
        print("adata_object.obs is empty. Skipping export.")

    # 2. Export adata_object.var (gene metadata)
    if not adata_object.var.empty:
        var_path = os.path.join(export_dir, "var_data.txt")
        save_dataframe_to_text(adata_object.var, var_path, "var (gene_metadata)")
    else:
        print("adata_object.var is empty. Skipping export.")

    # 3. Export adata_object.X (primary expression matrix)
    # Format: genes x cells, with gene names as 'V1' column
    if adata_object.X is not None:
        x_path = os.path.join(export_dir, "X_data.txt")
        # adata.X is obs x vars (cells x genes)
        save_matrix_to_text(adata_object.X.T, adata_object.var_names.astype(str), adata_object.obs_names.astype(str), x_path, "X (primary_matrix)")
    else:
        print("adata_object.X is None. Skipping export.")

    # 4. Export adata_object.layers
    if adata_object.layers:
        print("\n--- Exporting Layers ---")
        for layer_name, layer_data in adata_object.layers.items():
            print(f"Processing layer: {layer_name}")
            if layer_data is not None:
                layer_path = os.path.join(export_dir, f"layer_{layer_name}_data.txt")
                # Layers are typically obs x vars (cells x genes)
                save_matrix_to_text(layer_data.T, adata_object.var_names.astype(str), adata_object.obs_names.astype(str), layer_path, f"layer_{layer_name}")
            else:
                print(f"Layer '{layer_name}' is None. Skipping.")
    else:
        print("No layers found in adata_object.layers.")

    # 5. Export adata_object.obsm (multi-dimensional observation annotations)
    if adata_object.obsm:
        print("\n--- Exporting Obsm ---")
        for key, data_array in adata_object.obsm.items():
            print(f"Processing obsm key: {key}")
            if data_array is not None:
                obsm_path = os.path.join(export_dir, f"obsm_{key}_data.txt")
                # obsm arrays are typically n_obs x n_features
                try:
                    df = pd.DataFrame(data_array, index=adata_object.obs_names.astype(str))
                    # Add column names if it's a typical embedding like UMAP (usually 2 or more components)
                    df.columns = [f"{key}_{i}" for i in range(data_array.shape[1])]
                    save_dataframe_to_text(df, obsm_path, f"obsm_{key}")
                except Exception as e:
                    print(f"Could not convert obsm key '{key}' to DataFrame and save: {e}. Shape: {data_array.shape if hasattr(data_array, 'shape') else 'N/A'}")
            else:
                print(f"Obsm key '{key}' is None. Skipping.")

    else:
        print("No data found in adata_object.obsm.")
        
    # 6. Export adata_object.varm (multi-dimensional variable annotations)
    if adata_object.varm:
        print("\n--- Exporting Varm ---")
        for key, data_array in adata_object.varm.items():
            print(f"Processing varm key: {key}")
            if data_array is not None:
                varm_path = os.path.join(export_dir, f"varm_{key}_data.txt")
                # varm arrays are typically n_vars x n_features
                try:
                    df = pd.DataFrame(data_array, index=adata_object.var_names.astype(str))
                    df.columns = [f"{key}_{i}" for i in range(data_array.shape[1])]
                    save_dataframe_to_text(df, varm_path, f"varm_{key}")
                except Exception as e:
                     print(f"Could not convert varm key '{key}' to DataFrame and save: {e}. Shape: {data_array.shape if hasattr(data_array, 'shape') else 'N/A'}")
            else:
                print(f"Varm key '{key}' is None. Skipping.")
    else:
        print("No data found in adata_object.varm.")

    # 7. Export adata_object.uns (unstructured annotations) - simple cases
    if adata_object.uns:
        print("\n--- Exporting Uns (simple cases) ---")
        for key, data_item in adata_object.uns.items():
            print(f"Processing uns key: {key}")
            uns_path = os.path.join(export_dir, f"uns_{key}_data.txt")
            try:
                if isinstance(data_item, pd.DataFrame):
                    save_dataframe_to_text(data_item, uns_path, f"uns_{key} (DataFrame)")
                elif isinstance(data_item, pd.Series):
                     save_dataframe_to_text(data_item.to_frame(), uns_path, f"uns_{key} (Series)")
                elif isinstance(data_item, dict):
                    # Save dictionary as key-value pairs in a text file
                    with open(uns_path, 'w') as f:
                        for k_dict, v_dict in data_item.items():
                            f.write(f"{str(k_dict)}\t{str(v_dict)}\n")
                    print(f"uns_{key} (dict) saved to {uns_path}")
                elif isinstance(data_item, (list, np.ndarray)):
                    # Save list or numpy array, each element on a new line
                     np.savetxt(uns_path, np.array(data_item, dtype=object), fmt='%s', delimiter='\t')
                     print(f"uns_{key} (list/array) saved to {uns_path}")
                elif isinstance(data_item, (str, int, float, bool)):
                    with open(uns_path, 'w') as f:
                        f.write(str(data_item))
                    print(f"uns_{key} (scalar) saved to {uns_path}")
                else:
                    print(f"Skipping uns key '{key}': data type {type(data_item)} not handled for simple text export.")
            except Exception as e:
                print(f"Error exporting uns key '{key}': {e}")
    else:
        print("No data found in adata_object.uns.")

    print("\n--- Full object export process complete ---")

else:
    print("\nFailed to load the data into an AnnData object.")
    print("Please ensure the file exists at the specified path and is a valid H5AD file.")
    print("Path checked:", local_h5ad_full_path)
