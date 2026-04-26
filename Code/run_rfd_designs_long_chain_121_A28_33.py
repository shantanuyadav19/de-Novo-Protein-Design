import os
import subprocess

INPUT_PDB        = "/home/sbl_lab1/protein_design/long_chain_design/long_chain_pdb_new.pdb"
BASE_TENSOR_DIR  = "/home/sbl_lab1/protein_design/long_chain_design/clamp_tensors_121_A28_33"
BASE_OUTPUT_DIR  = "/home/sbl_lab1/protein_design/long_chain_design/RFDiffusion_designs_121_A28_33"
RFDIFF_SCRIPT    = "/home/sbl_lab1/protein_design/RFdiffusion/scripts/run_inference.py"

BINDER_LENGTHS      = [121]
DESIGNS_PER_LENGTH  = 100
HOTSPOTS            = "['A28','A29','A30','A31','A32','A33']"   # ← corrected (target-first)

os.makedirs(BASE_OUTPUT_DIR, exist_ok=True)

for blen in BINDER_LENGTHS:
    tensor_dir = os.path.join(BASE_TENSOR_DIR, f"len_{blen}")
    if not os.path.isdir(tensor_dir):
        print(f"[skip] {tensor_dir} not found")
        continue

    stems = sorted({
        f[:-len("_adj.pt")] for f in os.listdir(tensor_dir) if f.endswith("_adj.pt")
    })
    if not stems:
        print(f"[skip] no *_adj.pt tensors in {tensor_dir}")
        continue

    scaffold_list_path = os.path.join(BASE_TENSOR_DIR, f"len_{blen}_scaffolds.txt")
    with open(scaffold_list_path, "w") as fh:
        fh.write("\n".join(stems) + "\n")

    out_prefix = os.path.join(BASE_OUTPUT_DIR, f"rfdesigns_len_{blen}")

    env = os.environ.copy()
    env["CUDA_VISIBLE_DEVICES"]    = "GPU-2df4463e-120d-947a-2053-d23ef13eb073"
    env["PYTORCH_CUDA_ALLOC_CONF"] = "expandable_segments:True"

    command = [
        "python", RFDIFF_SCRIPT,
        f"inference.input_pdb={INPUT_PDB}",
        f"inference.output_prefix={out_prefix}",
        f"contigmap.contigs=[A1-71/0 {blen}]",           # ← target-first (cleaner numbering)
        "inference.model_runner=ScaffoldedSampler",
        f"ppi.hotspot_res={HOTSPOTS}",
        "scaffoldguided.scaffoldguided=True",
        f"scaffoldguided.scaffold_dir={tensor_dir}",
        f"scaffoldguided.scaffold_list={scaffold_list_path}",
        f"inference.num_designs={DESIGNS_PER_LENGTH}",
        "denoiser.noise_scale_ca=0.5",
        "denoiser.noise_scale_frame=0.5",
        "scaffoldguided.mask_loops=False",
    ]

    print(f"\n=== Binder length {blen}: {DESIGNS_PER_LENGTH} designs from {len(stems)} scaffolds ===")
    try:
        subprocess.run(command, env=env, check=True)
    except subprocess.CalledProcessError as e:
        print(f"Error at length {blen}: {e}")
        break

print("\nAll RFdiffusion runs complete.")