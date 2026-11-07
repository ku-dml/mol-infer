import sys
import pandas as pd
from pathlib import Path

# Command-line arguments
csv_file = sys.argv[1]
n_samples = int(sys.argv[2])

# Optional random seed
random_seed = 75

# Load dataset
df = pd.read_csv(csv_file)

if n_samples > len(df):
    print(f"Error: Dataset mein sirf {len(df)} rows hain.")
    sys.exit(1)

# Random sampling
sampled_df = df.sample(
    n=n_samples,
    random_state=random_seed
)

# Output filename
input_path = Path(csv_file)

output_file = input_path.with_name(
    f"{input_path.stem}_seed_{random_seed}_{n_samples}.csv"
)

# Save output
sampled_df.to_csv(output_file, index=False)

print(f"Input dataset: {csv_file}")
print(f"Total rows: {len(df)}")
print(f"Selected samples: {n_samples}")
print(f"Random seed: {random_seed}")
print(f"Output file: {output_file}")