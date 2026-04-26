import os
import subprocess

# Define your paths and parameters
input_pdb = "/Users/shantanu.yadav19/Downloads/MTP2Offline/Data/long_chain_pdb_new.pdb"
base_out_dir = "/Users/shantanu.yadav19/Downloads/MTP2Offline/Data/clamp_tensors_121_A30_34"

# Loop through the desired binder lengths (115 to 125 inclusive)
for binder_len in range(121, 122):
    
    # Create a specific output directory for this length to keep things organized
    out_dir = os.path.join(base_out_dir, f"len_{binder_len}")
    os.makedirs(out_dir, exist_ok=True)
    
    # Construct the command as a list of arguments
    command = [
        "python", "/Users/shantanu.yadav19/Downloads/MTP2Offline/Code/interface_tensors/make_interface_tensor.py",
        "--input_pdb", input_pdb,
        "--out_dir", f"{out_dir}/",
        "--binderlen", str(binder_len),
        "--target_adj", "A30-34",   # colon-separated targets
        "--binder_ss", "E",              # comma-separated, one per target
        "--binder_ss_len", "4",          # comma-separated, one per target
    ]
    
    print(f"Generating tensors for binder length {binder_len}...")
    
    try:
        # Execute the command
        subprocess.run(command, check=True)
    except subprocess.CalledProcessError as e:
        print(f"An error occurred while generating tensors for length {binder_len}.")
        print(e)
        break # Stops the loop if an error occurs

print("Tensor generation complete!")