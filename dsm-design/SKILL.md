---
name: dsm-design
description: "Delta-Sigma modulator design using the Delta Sigma Toolbox methodology. Based on the Delta Sigma Toolbox by Richard Schreier from 'Understanding Delta-Sigma Data Converters' (Appendix B). For Octave/MATLAB modeling, use Python/Octave hybrid workflow: Python generates coefficients, Octave simulates."
homepage: https://www.mathworks.com/matlabcentral/fileexchange/
metadata:
  {
    "openclaw":
      {
        "emoji": "📊",
        "requires": { "bins": ["python3", "octave"] },
        "optional": { "bins": ["matlab"] },
        "install": [],
      },
  }
---

# Delta-Sigma Modulator Design Skill

This skill guides you through designing Delta-Sigma modulators using the methodology from the Delta Sigma Toolbox by Richard Schreier, as documented in Appendix B of "Understanding Delta-Sigma Data Converters" (2nd Edition) by Pavan, Schreier, and Temes.

## Important: Recommended Workflow

**For Octave/MATLAB modeling, use the Python/Octave hybrid approach:**

```
Python (synthesizeNTF) → Coefficients → Octave (stuffABCD + simulation)
```

This is necessary because:
1. `synthesizeNTF` in Delta Sigma Toolbox has API incompatibilities with modern Octave
2. Python deltasigma package works reliably (with minor patches for Python 3.10+)
3. Octave can still use `stuffABCD` and simulation functions from toolbox

## Quick Start: Python/Octave Hybrid

### Step 1: Python Generates Coefficients

```python
#!/usr/bin/env python3
import sys
sys.path.insert(0, 'python-deltasigma-master')
from deltasigma import synthesizeNTF, realizeNTF

# Your specifications
order = 3
OSR = 16
H_inf = 2.0
opt = 1  # optimized zeros

# Synthesize NTF
ntf = synthesizeNTF(order, OSR, opt, H_inf, 0)
z, p, k = ntf

# Realize coefficients
a, g, b, c = realizeNTF(ntf, 'CIFF')

# Save for Octave
with open('coefficients.txt', 'w') as f:
    f.write(f"a = [{' '.join([f'{x:.10f}' for x in a])}];\n")
    f.write(f"g = [{g:.10f}];\n")
    f.write(f"b = [{' '.join([f'{x:.10f}' for x in b])}];\n")
    f.write(f"c = [{' '.join([f'{x:.10f}' for x in c])}];\n")
```

### Step 2: Octave Simulates

```matlab
% Load coefficients from Python
run('coefficients.txt');

% Build ABCD using toolbox stuffABCD
addpath('python-deltasigma-master/delsig/');
ABCD = stuffABCD(a, g, b, c, 'CIFF');

% Run simulation (see full example below)
% ... simulation code ...
```

## Full Working Example

```matlab
%% DSM Design - Python/Octave Hybrid Workflow
%% 3rd-Order CIFF, fs=2MHz, OSR=16, 5-bit, H_inf=2.0

clear all; close all; clc;

%% Step 1: Load coefficients from Python
fprintf('Loading coefficients from Python deltasigma...\n');
run('coefficients_from_python.txt');

fprintf('Coefficients:\n');
fprintf('  a = ['); fprintf('%.6f ', a); fprintf(']\n');
fprintf('  g = %.6f\n', g);
fprintf('  b = ['); fprintf('%.4f ', b); fprintf(']\n');
fprintf('  c = ['); fprintf('%.4f ', c); fprintf(']\n\n');

%% Step 2: Build ABCD using toolbox stuffABCD
addpath('python-deltasigma-master/delsig/');
fprintf('[Building ABCD with stuffABCD]\n');
ABCD = stuffABCD(a, g, b, c, 'CIFF');

%% Step 3: Simulation parameters
order = 3;
fs = 2e6;
OSR = 16;
fB = fs / (2*OSR);
A_in = 0.2;

n_bits = 5;
n_levels = 32;
V_fs = 1.0;
LSB = 2*V_fs / (n_levels - 1);
q_levels = linspace(-V_fs, V_fs, n_levels);

%% Step 4: Run simulation
N = 8192;
f_bin = round(sqrt(1/7) * N / (2*OSR));
u = A_in * sin(2*pi*f_bin*(0:N-1)/N);

n_states = size(ABCD, 1) - 1;
x = zeros(n_states, 1);
y = zeros(1, N);
v = zeros(1, N);

A_mat = ABCD(1:n_states, 1:n_states);
B_mat = ABCD(1:n_states, n_states+1:end);
C_mat = ABCD(n_states+1, 1:n_states);
D_mat = ABCD(n_states+1, n_states+1:end);

for i = 1:N
    if i == 1, v_prev = 0; else, v_prev = v(i-1); end
    y(i) = C_mat * x + D_mat(1)*u(i) + D_mat(2)*v_prev;
    
    % Quantize
    y_clip = max(-V_fs, min(V_fs, y(i)));
    idx = round((y_clip - (-V_fs)) / LSB) + 1;
    idx = max(1, min(n_levels, idx));
    v(i) = q_levels(idx);
    
    % Update state
    x = A_mat*x + B_mat(:,1)*u(i) + B_mat(:,2)*v(i);
end

%% Step 5: Calculate SNR
w = 0.5 * (1 - cos(2*pi*(0:N-1)/N));
V_out = fft(v .* w) / (N/4);
V_out_mag = abs(V_out);

[~, sig_idx] = max(V_out_mag(2:N/2));
sig_bin = sig_idx + 1;
fB_bins = ceil(N / (2*OSR));
sig_bins = sig_bin-1:sig_bin+1;
signal_power = sum(V_out_mag(sig_bins).^2);
noise_bins = setdiff(2:fB_bins, sig_bins);
noise_power = sum(V_out_mag(noise_bins).^2);
SNR = 10*log10(signal_power / noise_power);

fprintf('SNR: %.2f dB\n', SNR);
```

## Core Concepts

### Key Parameters

| Parameter | Description | Typical Values |
|-----------|-------------|----------------|
| `order` | Modulator order | 1-5 for stable designs |
| `OSR` | Oversampling Ratio | 32, 64, 128 |
| `H_inf` | Max out-of-band NTF gain | 1.5-2.0 (Lee's rule: < 2) |
| `f0` | Center frequency | 0 for LP, 0.25 for BP at fs/4 |
| `n_bits` | Quantizer bits | 1-8 |
| `n_levels` | Quantizer levels | 2^n_bits |
| `V_fs` | Full scale voltage | +/- 1 (normalized) |

### Design Flow (Hybrid)

```
┌─────────────────────────────────────────────────────────────┐
│  PYTHON: NTF Synthesis                                      │
│  ─────────────────────                                      │
│  1. Call synthesizeNTF(order, OSR, opt, H_inf, f0)         │
│  2. Returns: zeros, poles, gain                            │
│  3. Call realizeNTF(ntf, 'CIFF')                           │
│  4. Returns: a, g, b, c coefficients                      │
│  5. Save to file (coefficients.txt)                         │
└─────────────────────┬───────────────────────────────────────┘
                      │
                      ▼
┌─────────────────────────────────────────────────────────────┐
│  OCTAVE: Simulation                                          │
│  ─────────────────                                           │
│  6. Load coefficients from file                              │
│  7. Call stuffABCD(a, g, b, c, 'CIFF')                     │
│  8. Returns: ABCD matrix                                    │
│  9. Run time-domain simulation                              │
│  10. Calculate SNR from output spectrum                     │
│  11. Plot results                                            │
└─────────────────────────────────────────────────────────────┘
```

## Python Environment Setup

### Install dependencies:
```bash
# Download python-deltasigma
curl -L "https://github.com/ggventurini/python-deltasigma/archive/refs/heads/master.zip" -o deltasigma.zip
unzip deltasigma.zip

# Apply Python 3.10+ compatibility patches:
# 1. In _utils.py: change 'gcd = fractions.gcd' to 'gcd = math.gcd'
# 2. Add 'import math' at top of _utils.py
# 3. Change 'collections.Iterable' to 'collections.abc.Iterable'
```

### Python script template:
```python
#!/usr/bin/env python3
import sys
sys.path.insert(0, 'python-deltasigma-master')
from deltasigma import synthesizeNTF, realizeNTF

# Your specs here
order = 3
OSR = 16
H_inf = 2.0
opt = 1
form = 'CIFF'

# Generate NTF
ntf = synthesizeNTF(order, OSR, opt, H_inf, 0)

# Get coefficients
a, g, b, c = realizeNTF(ntf, form)

# Print for Octave
print(f"a = [{' '.join([f'{x:.10f}' for x in a])}];")
print(f"g = [{g:.10f}];")
print(f"b = [{' '.join([f'{x:.10f}' for x in b])}];")
print(f"c = [{' '.join([f'{x:.10f}' for x in c])}];")
```

## Critical Points

### Why Not Pure Octave?

| Function | Status | Issue |
|----------|--------|-------|
| `synthesizeNTF` | **BROKEN** | Uses `tf.z`, `ss.z` incompatible with octave-control v4.x+ |
| `realizeNTF` | ✓ Works | No control package dependency |
| `stuffABCD` | ✓ Works | Pure MATLAB code, no dependencies |
| `simulateDSM` | ✓ Works | Pure MATLAB code |

**Note:** Even control v3.2.0 has compatibility issues. The `tf` object API changed.

### SNR Calculation (Critical!)

**CORRECT:** Find signal in **OUTPUT spectrum**
```matlab
V_out = fft(v .* w) / (N/4);
[~, sig_idx] = max(V_out_mag(2:N/2));
sig_power = sum(V_out_mag(sig_bins).^2);
noise_power = sum(V_out_mag(noise_bins).^2);
SNR = 10*log10(sig_power / noise_power);
```

**WRONG:** Finding signal in **noise spectrum**
```matlab
noise = v - u;
V_noise = fft(noise .* w);  % DON'T find signal here!
```

## Example Specifications

### Example 1: 2nd-Order CIFF
```python
order = 2
OSR = 64
H_inf = 1.5
opt = 1
fs = 1e6
n_bits = 5
A_in = 0.5

# Expected output:
# a ≈ [0.9996, 0.9984]
# g ≈ 0.0032
# SNR ≈ 70 dB
```

### Example 2: 3rd-Order CIFF
```python
order = 3
OSR = 16
H_inf = 2.0
opt = 1
fs = 2e6
n_bits = 5
A_in = 0.2

# Expected output:
# a ≈ [1.3406, 0.7291, 0.1340]
# g ≈ 0.0230
# SNR ≈ 74 dB
```

## Troubleshooting

### Python deltasigma import fails
```
AttributeError: module 'fractions' has no attribute 'gcd'
```
**Fix:** Edit `_utils.py`:
- Add `import math`
- Change `gcd = fractions.gcd` to `gcd = math.gcd`
- Change `collections.Iterable` to `collections.abc.Iterable`

### Octave stuffABCD C matrix looks wrong
If C = [a(1), a(2), ...] instead of proper feedforward:
- Use hardcoded values from successful runs
- Or use manual ABCD construction

### Simulation unstable (states diverge)
- Reduce input amplitude (try 0.1V instead of 0.5V)
- Check H_inf not too high for order
- Verify coefficients from Python

## References

1. Pavan, Schreier, Temes - "Understanding Delta-Sigma Data Converters" (2nd Ed.), Appendix B
2. Schreier, R. - "The Delta-Sigma Toolbox" (MATLAB Central)
3. Python deltasigma: https://github.com/ggventurini/python-deltasigma
4. Delta Sigma Toolbox: https://github.com/ggventurini/python-deltasigma

## Reference Designs

Located in `references/` directory:

### 4th-Order CIFF Design (2MSPS, OSR=16, 5-bit, H_inf=4.0)

**Files:**
- `generate_4th_order.py` - Python coefficient generator
- `simulate_4th_order.m` - Octave simulation script
- `coefficients_4th_order.txt` - Generated coefficients
- `dsm_4th_order.png` - Result plots

**Specifications:**
| Parameter | Value |
|-----------|-------|
| Order | 4 |
| fs | 2 MHz |
| OSR | 16 (BW = 62.5 kHz) |
| H_inf | 4.0 |
| Topology | CIFF |
| Quantizer | 5-bit (32 levels) |

**Coefficients (from Python deltasigma):**
```matlab
a = [2.566151 2.736364 1.365746 0.218951];  % feedback coefficients
g = [0.004450 0.028318];                      % resonator coefficients
b = [1 0 0 0 1];                              % input coefficients
c = [1 1 1 1];                                % inter-stage coefficients
```

**Results:**
| Input | SNR | ENOB | Status |
|-------|-----|------|--------|
| 0.1V | 89.82 dB | 14.63 bits | ✓ Stable |
| 0.5V | 103.90 dB | 16.97 bits | ✓ Stable |

**ABCD Matrix:**
```
[  1.000000  -0.004450   0.000000   0.000000   1.000000  -1.000000 ]
[  1.000000   1.000000   0.000000   0.000000   0.000000   0.000000 ]
[  0.000000   1.000000   1.000000  -0.028318   0.000000   0.000000 ]
[  0.000000   0.000000   1.000000   1.000000   0.000000   0.000000 ]
[  2.566151   2.736364   1.365746   0.218951   1.000000   0.000000 ]
```

**Usage:**
```bash
# Generate coefficients
cd /home/ubuntu/DSM_test
python3 references/generate_4th_order.py

# Simulate in Octave (includes NTF/STF plots)
octave --no-gui references/simulate_4th_order.m
```

**Output Files:**
- `dsm_4th_order.png` - Time domain, output spectrum, noise spectrum, histogram, status, ABCD matrix
- `dsm_4th_order_ntf_stf.png` - NTF and STF frequency response plots

**NTF/STF Features:**
The merged simulation script includes:
1. **NTF Calculation** - Uses exact zeros/poles from Python `synthesizeNTF`
   - Deep notches at optimized zero locations
   - H_inf marked (~12 dB)
   - Bandwidth edge marked
2. **STF Calculation** - CIFF has unity STF (0 dB)
   - Flat response in signal band
   - DC gain marked

---
*Reference designs demonstrate proven working configurations*

