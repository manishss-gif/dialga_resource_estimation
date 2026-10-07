# Quantum Resource Estimation of the Dialga Cipher Family

This repository contains standalone, reproducible Clifford+$T$ quantum implementations, Grover phase-oracle constructions, and complete key-search resource estimators for the **Dialga** tweakable block-cipher family.

---

## 1. Overview of Cipher Variants & Implementation Profiles

The Dialga family comprises four configurations:
- **Dialga-128**: 128-bit block, 256-bit key, 128-bit tweak (20 rounds).
- **Dialga-128† (128d)**: 128-bit block, 256-bit key, 128-bit tweak (16 rounds).
- **Dialga-256**: 128-bit block, 256-bit key, 256-bit tweak (20 rounds).
- **Dialga-256† (256d)**: 128-bit block, 256-bit key, 256-bit tweak (16 rounds).

For each variant, two complementary reversible circuit organizations are provided:
- **Low-Width (LW)**: Minimizes qubit width through sequential ancilla reuse across S-box evaluations.
- **Low-Depth (LD)**: Minimizes $T$-depth and circuit depth by evaluating all 32 S-boxes concurrently with dedicated ancillae.

Both profiles implement the exact classical Dialga transformation and produce identical additive gate counts.

---

## 2. Requirements & Installation

- **Python**: 3.14.4
- **Dependencies**: ProjectQ 0.8.0, mpmath 1.4.1

Install dependencies:
```bash
pip install projectq==0.8.0 mpmath==1.4.1
```

**Note:** ProjectQ 0.8.0 builds a C++ extension, and recent `setuptools` releases break that build (`AttributeError: 'Compiler' object has no attribute 'dry_run'`). If the installation fails, install with an older `setuptools` as a build constraint:
```bash
echo "setuptools<75" > constraints.txt
pip install --build-constraint constraints.txt -r requirements.txt   # pip >= 25.3
# with older pip:  PIP_CONSTRAINT=constraints.txt pip install -r requirements.txt
# or, with uv:
uv pip install --build-constraint constraints.txt -r requirements.txt
```

---

## 3. Quick Start / One-Click Execution

- **Reproduce All Full Resource Estimation Reports (Tables 2, 3, 4, 5, 6, 7)**:
  ```bash
  ./run_all_resources.sh
  ```
  *Outputs and logs for all 8 configurations are automatically displayed and saved in `results/`.*

- **Run Shared Correctness & Regression Tests**:
  ```bash
  ./run_all_tests.sh
  ```

*Expected run time on a laptop: about 10 minutes for `run_all_resources.sh` and about 15 minutes for `run_all_tests.sh`.*

---

## 4. How to Reproduce Paper Tables (Tables 2, 3, 4, 5, 6, 7)

### Table 2 — Three-Tuple Grover Phase Oracle

Reports resources for the phase oracle: 3 parallel encryptions ($\text{ENC}$), 384-bit ciphertext comparison, phase inversion ($Z_{\text{mc}}$), 384-bit comparator uncomputation, and 3 parallel uncomputations ($\text{ENC}^\dagger$).

*These are reported in section **[2]** of each `*_full_grover.py` report:*

| Cipher Variant & Profile | Command to Run | Expected Resources |
| :--- | :--- | :--- |
| **Dialga-128 (LW)** | `python3 dialga128_ct_lw_full_grover.py --resource` *([2])* | 361,206 CNOT, 96,764 T, T-depth 7,936, Width 2,117 |
| **Dialga-128 (LD)** | `python3 dialga128_ct_ld_full_grover.py --resource` *([2])* | 361,206 CNOT, 96,764 T, T-depth 177, Width 2,495 |
| **Dialga-128† (LW)** | `python3 dialga128d_ct_lw_full_grover.py --resource` *([2])* | 292,086 CNOT, 78,332 T, T-depth 6,400, Width 2,117 |
| **Dialga-128† (LD)** | `python3 dialga128d_ct_ld_full_grover.py --resource` *([2])* | 292,086 CNOT, 78,332 T, T-depth 145, Width 2,495 |
| **Dialga-256 (LW)** | `python3 dialga256_ct_lw_full_grover.py --resource` *([2])* | 456,438 CNOT, 121,340 T, T-depth 9,984, Width 2,501 |
| **Dialga-256 (LD)** | `python3 dialga256_ct_ld_full_grover.py --resource` *([2])* | 456,438 CNOT, 121,340 T, T-depth 177, Width 2,879 |
| **Dialga-256† (LW)** | `python3 dialga256d_ct_lw_full_grover.py --resource` *([2])* | 364,278 CNOT, 96,764 T, T-depth 7,936, Width 2,501 |
| **Dialga-256† (LD)** | `python3 dialga256d_ct_ld_full_grover.py --resource` *([2])* | 364,278 CNOT, 96,764 T, T-depth 145, Width 2,879 |

---

### Table 3 — 256-bit Grover Diffusion Resources ($D_{256}$)

Reports logical Clifford+$T$ resources for the standalone 256-bit Grover diffusion operator $D_{256}$ acting on the key register:

| Operator | Verification Command | Expected Resources |
| :--- | :--- | :--- |
| **Diffusion ($D_{256}$)** | `python3 dialga128_ct_lw_full_grover.py --resource` *([3])* | 2,550 CNOT, 3,065 1qClifford, 1,020 T, 255 Measure, T-depth 8, Full depth 63, Width 639 |

*Note: In the complete Grover iteration (Table 5), the diffusion operator reuses the 384-bit comparator workspace, leaving the peak circuit width unchanged.*

---

### Table 4 — Reversible Encryption Resources (Single Primitive)

Reports `#CNOT`, `#1qClifford`, `#T`, `#Measure`, `T-depth`, `Full depth`, and `Width` for a single encryption:

| Cipher Variant & Profile | Command to Run | Expected Resources |
| :--- | :--- | :--- |
| **Dialga-128 (LW)** | `python3 qdialga128_ct_lw.py --resource` | 59,392 CNOT, 15,872 T, T-depth 3,968, Width 514 |
| **Dialga-128 (LD)** | `python3 qdialga128_ct_ld.py --resource` | 59,392 CNOT, 15,872 T, T-depth 88, Width 640 |
| **Dialga-128† (LW)** | `python3 qdialga128d_ct_lw.py --resource` | 47,872 CNOT, 12,800 T, T-depth 3,200, Width 514 |
| **Dialga-128† (LD)** | `python3 qdialga128d_ct_ld.py --resource` | 47,872 CNOT, 12,800 T, T-depth 72, Width 640 |
| **Dialga-256 (LW)** | `python3 qdialga256_ct_lw.py --resource` | 75,264 CNOT, 19,968 T, T-depth 4,992, Width 642 |
| **Dialga-256 (LD)** | `python3 qdialga256_ct_ld.py --resource` | 75,264 CNOT, 19,968 T, T-depth 88, Width 768 |
| **Dialga-256† (LW)** | `python3 qdialga256d_ct_lw.py --resource` | 59,904 CNOT, 15,872 T, T-depth 3,968, Width 642 |
| **Dialga-256† (LD)** | `python3 qdialga256d_ct_ld.py --resource` | 59,904 CNOT, 15,872 T, T-depth 72, Width 768 |

---

### Table 5 — One Complete Grover Iteration

Reports resources for One Grover Iteration = Phase Oracle ($U_f$) + 256-bit Diffusion Operator ($D_{256}$):

*These are reported in section **[4]** of each `*_full_grover.py` report:*

| Cipher Variant & Profile | Command to Run | Expected Resources |
| :--- | :--- | :--- |
| **Dialga-128 (LW)** | `python3 dialga128_ct_lw_full_grover.py --resource` *([4])* | 363,756 CNOT, 97,784 T, T-depth 7,944, Width 2,117 |
| **Dialga-128 (LD)** | `python3 dialga128_ct_ld_full_grover.py --resource` *([4])* | 363,756 CNOT, 97,784 T, T-depth 185, Width 2,495 |
| **Dialga-128† (LW)** | `python3 dialga128d_ct_lw_full_grover.py --resource` *([4])* | 294,636 CNOT, 79,352 T, T-depth 6,408, Width 2,117 |
| **Dialga-128† (LD)** | `python3 dialga128d_ct_ld_full_grover.py --resource` *([4])* | 294,636 CNOT, 79,352 T, T-depth 153, Width 2,495 |
| **Dialga-256 (LW)** | `python3 dialga256_ct_lw_full_grover.py --resource` *([4])* | 458,988 CNOT, 122,360 T, T-depth 9,992, Width 2,501 |
| **Dialga-256 (LD)** | `python3 dialga256_ct_ld_full_grover.py --resource` *([4])* | 458,988 CNOT, 122,360 T, T-depth 185, Width 2,879 |
| **Dialga-256† (LW)** | `python3 dialga256d_ct_lw_full_grover.py --resource` *([4])* | 366,828 CNOT, 97,784 T, T-depth 7,944, Width 2,501 |
| **Dialga-256† (LD)** | `python3 dialga256d_ct_ld_full_grover.py --resource` *([4])* | 366,828 CNOT, 97,784 T, T-depth 153, Width 2,879 |

---

### Table 6 — Depth–Width Cost and MAXDEPTH Figures

Table 6 and the NIST MAXDEPTH figures in the text follow by direct arithmetic from Table 5 and $N_G$:
$\mathrm{DW} = \text{Full depth} \times \text{Width}$ of one Grover iteration, and for a depth bound $D_{\max}$ the number of parallel Grover instances is $p \approx (D_{\text{search}}/D_{\max})^2$, where $D_{\text{search}}$ is the full depth reported in section **[7]**.

---

### Table 7 — Complete Grover Key Search ($N_G = \lfloor \frac{\pi}{4} 2^{128} \rfloor$ Iterations)

Reports total resources including initial state preparation, all $N_G = 267,257,146,016,241,686,964,920,093,290,467,695,825$ Grover iterations, and final key readout:

*These are reported in section **[7]** of each `*_full_grover.py` report:*

| Cipher Variant & Profile | Command to Run | Expected Full Attack Resources |
| :--- | :--- | :--- |
| **Dialga-128 (LW)** | `python3 dialga128_ct_lw_full_grover.py --resource` | $1.0898 \times 2^{146}$ CNOT, $1.1719 \times 2^{144}$ T, Width 2,117 |
| **Dialga-128 (LD)** | `python3 dialga128_ct_ld_full_grover.py --resource` | $1.0898 \times 2^{146}$ CNOT, $1.1719 \times 2^{144}$ T, Width 2,495 |
| **Dialga-128† (LW)** | `python3 dialga128d_ct_lw_full_grover.py --resource` | $1.7655 \times 2^{145}$ CNOT, $1.9019 \times 2^{143}$ T, Width 2,117 |
| **Dialga-128† (LD)** | `python3 dialga128d_ct_ld_full_grover.py --resource` | $1.7655 \times 2^{145}$ CNOT, $1.9019 \times 2^{143}$ T, Width 2,495 |
| **Dialga-256 (LW)** | `python3 dialga256_ct_lw_full_grover.py --resource` | $1.3752 \times 2^{146}$ CNOT, $1.4664 \times 2^{144}$ T, Width 2,501 |
| **Dialga-256 (LD)** | `python3 dialga256_ct_ld_full_grover.py --resource` | $1.3752 \times 2^{146}$ CNOT, $1.4664 \times 2^{144}$ T, Width 2,879 |
| **Dialga-256† (LW)** | `python3 dialga256d_ct_lw_full_grover.py --resource` | $1.0990 \times 2^{146}$ CNOT, $1.1719 \times 2^{144}$ T, Width 2,501 |
| **Dialga-256† (LD)** | `python3 dialga256d_ct_ld_full_grover.py --resource` | $1.0990 \times 2^{146}$ CNOT, $1.1719 \times 2^{144}$ T, Width 2,879 |

---

## 5. File Manifest

| File | Description |
| :--- | :--- |
| `corefunctions.py` | Shared Clifford+$T$ primitives, S-box ($\mathsf{Sb}_0$), circulant linear layer ($M$), and permutations. |
| `dialga_grover_common.py` | Shared 384-bit comparator, 256-bit diffusion operator ($D_{256}$, Table 3), $N_G$ calculation, and the test and report driver used by all `*_full_grover.py` scripts. |
| `test_advanced_counter.py` | `AdvancedResourceCounter` DAG-based logical depth and gate counter. |
| `classical_dialga_reference.py` | Classical reference for Dialga-128 / Dialga-128†. |
| `classical_dialga256_reference.py` | Classical reference for Dialga-256 / Dialga-256†. |
| `qdialga*_ct_lw.py` | Low-Width reversible encryption circuits (Table 4). |
| `qdialga*_ct_ld.py` | Low-Depth reversible encryption circuits (Table 4). |
| `dialga*_ct_*_grover.py` | 3-tuple Grover phase oracle and uncomputation circuits (Table 2). |
| `dialga*_ct_*_full_grover.py` | Per-variant definitions (test tuples, circuits, reference values) run by the shared driver: phase oracle (Table 2), iteration (Table 5), full key search (Table 7). |
| `test_ct_primitives.py` | Unit tests for Clifford+$T$ primitives (AND gate, uncomputation, SubCell, 8-bit S-box). |
| `test_ct_profiles.py` | Cross-verification suite of all 8 quantum profiles against classical reference. |
| `test_all_ct_profiles.py` | Runs `test_ct_primitives.py` and `test_ct_profiles.py` together. |
| `test_dialga_grover_common.py` | Unit tests for comparator and diffusion operators (Table 3). |
| `run_all_resources.sh` | One-step automated execution of all resource estimators. |
| `run_all_tests.sh` | One-step automated execution of all test suites. |
| `LICENSE` | Anonymous MIT License. |
