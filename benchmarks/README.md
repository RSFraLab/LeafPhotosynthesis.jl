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

## Coupled A–Ci equilibrium

`equilibrium_batches.jl` separately measures coupled C3–CO₂ diffusion–Ball–Berry
solving with safeguarded Newton and bisection. The 1,000,000-leaf input repeats
108 combinations of three temperatures, four light levels, three ambient CO₂
values, and three VPD values specified in `workload.toml`. Ball–Berry uses the
prescribed RH, so VPD variation does not change that particular law. This is a
synthetic batch, not a realistic canopy population. One quarter of the leaves
are dark and take the exact zero-iteration path. No warm-start guess is supplied.

```sh
CUDA_VISIBLE_DEVICES=1 julia --threads=8 --project=test/gpu \
  benchmarks/equilibrium_batches.jl /tmp/equilibrium-throughput.csv
```

Measured on the same host/GPU/software listed above, November 5, 2026:

| Precision | Method | 8-thread CPU | GPU evaluation | Mean CPU iterations | Maximum CPU iterations |
| --- | --- | ---: | ---: | ---: | ---: |
| FP32 | Safeguarded Newton | 29.70 ms | 0.203 ms | 2.92 | 8 |
| FP32 | Bisection | 71.40 ms | 0.233 ms | 14.06 | 24 |
| FP64 | Safeguarded Newton | 35.26 ms | 5.349 ms | 3.94 | 28 |
| FP64 | Bisection | 120.19 ms | 6.949 ms | 24.42 | 37 |

Newton reduced CPU time by about 2.4×/3.4× and GPU time by 1.15×/1.30× in FP32/FP64.
Fewer iterations do not translate proportionally into GPU speed: per-step work,
branch divergence, output traffic, and the device's FP64 throughput also matter.
The Newton GPU evaluation was about 146×/6.6× faster than eight CPU threads.
Upload + evaluation + download was 6.12 ms/15.92 ms for Newton, reducing that
comparison to approximately 4.9×/2.2×. Transfers were warmed before timing.

[Raw results](results/wurst-l40s-equilibrium-2026-11-05.csv) include separate transfers
and GPU iteration statistics. All CPU/GPU leaves converged; sampled fluxes were
cross-checked. FP32 bisection has slightly different iteration counts across devices
because rounding can change which iteration meets the same residual tolerances.
Timings are the minimum of three warmed CPU and five GPU event samples; allocations,
compilation, input construction, and checking are excluded. Transfers are single
warmed samples. This shared host was not isolated for statistical benchmarking.
These measurements exclude prognostic integration, sensitivities, O₂, leaf energy
balance, canopy transport, and ESM coupling. No full-model speedup is implied.
