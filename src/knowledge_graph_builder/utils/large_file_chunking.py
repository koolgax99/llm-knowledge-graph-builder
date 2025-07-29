import os

# ==== Configuration ====
input_file = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data/output/data_output_uv_unep/txt/unep_UV_radiation_report.txt"
output_dir = "/home/exouser/masters-thesis/ai-knowledge-graph-main/data/output/data_output_uv_unep/txt/split_parts"

# ==== Prepare Output Folder ====
os.makedirs(output_dir, exist_ok=True)

# ==== Read and Split ====
lines_per_file = 500_000 // 24

with open(input_file, 'r') as f:
    lines = f.readlines()

for i in range(24):
    start = i * lines_per_file
    end = (i + 1) * lines_per_file if i < 23 else None  # last chunk gets the rest
    chunk = lines[start:end]

    output_path = os.path.join(output_dir, f"unep_UV_radiation_report_part_{i + 1}.txt")
    with open(output_path, 'w') as f_out:
        f_out.writelines(chunk)

print(f"File successfully split into 6 parts in: {output_dir}")