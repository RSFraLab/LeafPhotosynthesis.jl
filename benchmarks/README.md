# Independent leaf throughput

This benchmark evaluates C3 photosynthesis over independent leaf environments,
using the same scalar function on CPU and GPU. It excludes stomatal dynamics,
oxygen diagnostics, canopy aggregation, radiative transfer, and ESM coupling.
Settings are in [workload.toml](workload.toml). The proposed workload is
180 × 360 × 0.3 × 3000 = 58,320,000 leaves; the land fraction is a sizing estimate,
not a geographic mask.

Run from the package root, selecting an idle GPU explicitly:

```sh
CUDA_VISIBLE_DEVICES=1 julia --threads=8 --project=test/gpu \
  benchmarks/leaf_batches.jl /tmp/leaf-throughput.csv
```

## Measured on wurst

November 5, 2026; physical GPU 1, NVIDIA L40S; Julia 1.12.6; CUDA.jl 6.5.1.
Host: two Intel Xeon Platinum 8462Y+ sockets, 128 logical CPUs. The threaded
comparison uses eight Julia threads, not the entire host.

| Precision | Leaves | Input + output payload | Serial CPU | 8-thread CPU | GPU evaluation | Upload + evaluation + download |
| --- | ---: | ---: | ---: | ---: | ---: | ---: |
| FP32 | 58,320,000 | 2.33 GB | 8.68 s | 1.02 s | 3.46 ms | 0.216 s |
| FP64 | 58,320,000 | 4.67 GB | 8.25 s | 1.00 s | 69.49 ms | 0.476 s |

GPU-resident evaluation was about 295× and 14.4× faster than the eight-thread
CPU comparison. Including both transfers reduces those ratios to approximately
4.7× and 2.1×. Keep this distinction when considering an ESM implementation: data
residency and the amount of retained diagnostics can dominate total runtime.

[Raw results](results/wurst_l40s_2026-11-05.csv) also include smaller batches.
These are exploratory measurements, with one warmed CPU sample and the minimum
of three warmed GPU event timings per case. They are not robust benchmark statistics
or an end-to-end model speedup. GPU timings exclude allocation, compilation, and
transfers; the final column sums separately measured upload, evaluation, and download.
Host preparation and allocations are excluded from all reported evaluation times.

The input array repeats the heterogeneous environments from all three Bonan 11.1
curves. Every leaf's temperature responses are evaluated; there is no layer-level
rate reuse or canopy reduction. The output retains six scalar diagnostics per leaf.
The benchmark checks sampled CPU/GPU outputs after downloading; the separate CUDA
test suite checks every reference-environment result, both precisions, dark/light
ordering, and the new optional components. Larger coupled kernels need their own
timing and memory budgets before predicting global-model performance.
