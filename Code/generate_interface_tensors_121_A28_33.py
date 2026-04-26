import os
import subprocess

# Define your paths and parameters
input_pdb = "/Users/shantanu.yadav19/Downloads/MTP2Offline/Data/long_chain_pdb_new.pdb"
base_out_dir = "/Users/shantanu.yadav19/Downloads/MTP2Offline/Data/clamp_tensors_121_A28_33"

# We only need one length (121 aa)
binder_len = 121

# Create output directory
out_dir = os.path.join(base_out_dir, f"len_{binder_len}")
os.makedirs(out_dir, exist_ok=True)

# Corrected command
command = [
    "python",
    "/Users/shantanu.yadav19/Downloads/MTP2Offline/Code/interface_tensors/make_interface_tensor.py",
    "--input_pdb", input_pdb,
    "--out_dir", f"{out_dir}/",
    "--binderlen", str(binder_len),
    "--target_adj", "A28-33",      # 6-residue interface
    "--binder_ss", "E",
    "--binder_ss_len", "6",        # ← fixed: must match target length
]

print(f"Generating tensors for binder length {binder_len} with 6-residue interface (A28-33)...")

try:
    subprocess.run(command, check=True)
    print("Tensor generation completed successfully!")
except subprocess.CalledProcessError as e:
    print(f"Error during tensor generation: {e}")