# 💻 Source Code and Experiments for “Towards Optimal Differential Attacks on FLY and PIPO”

This repository contains the full implementation and experimental scripts for the paper:

> **Insung Kim, et al.**,  
> *“Towards Optimal Differential Attacks on FLY and PIPO”*,  
> submitted to **Designs, Codes and Cryptography**, 2025.

It provides the codebase for enumerating high-probability differential trails, performing clustering effect analysis, and estimating the overall key recovery complexity for the lightweight block ciphers **FLY** and **PIPO**.

---

## ⚙️ Environment Compatibility

### ✅ Windows
- Visual Studio 2022  
- MSVC v143 (cl version: 14.33.31629)  
- Python 3.10  

### ✅ Linux
- Ubuntu 22.04  
- GCC 9.4, glibc 2.31  
- Python 3.10  

> 🔧 The project builds shared libraries:  
> - Windows: `.dll` via MSVC  
> - Linux: `.so` via GCC and glibc  
> Please ensure your environment supports these binary interface requirements.

---

## 📁 Directory Structure

Each cipher has its own folder with the following layout:

### 🔹 `PIPO/` and `FLY/`

```
├── Handle_LINUX/      # Shared library interface (Linux)
├── Handle_Windows/    # DLL interface and test harness (Windows)
├── Lib_Source/        # Common C source code: trail search, clustering, recovery logic
├── PIPO_LINUX/        # Platform-specific build scripts (or FLY_LINUX/)
├── PIPO_Windows/      # Visual Studio project (or FLY_Windows/)
├── Python/            # Python ctypes wrappers and analysis scripts
├── x64/               # Windows build artifacts
└── PIPO.sln           # Visual Studio solution (or FLY.sln)
```

---

## 🔬 Analysis Workflow

Our analysis proceeds as follows:

1. **Differential Trail Enumeration**  
   - Enumerate all trails with expected differential probability (EDP) ≥ 2⁻⁶³.

2. **Distinguisher Filtering**  
   - Retain trails whose complexity (under known bounds) is low enough to mount a key-recovery attack using a single right pair.

3. **Clustering Effect Estimation**  
   - Compute groups of trails sharing the same I/O difference, each with EDP ≥ (threshold + gap).

4. **Key Recovery Complexity Calculation**  
   - Under the key ranking paradigm, calculate the number of right pairs and advantage needed to achieve ≥90% success, and reflect this in the final key recovery complexity.

---

## 🐍 Python Wrappers

Python bindings (via `ctypes`) are provided to simplify analysis and integration.

### Key Functions

- **`get_BDP`**  
  Returns log₂ of the best differential probability up to the given round using Matsui’s branch-and-bound algorithm.

- **`get_all_BDP`**  
  Returns all input-output difference pairs with EDP ≥ threshold.  
  Duplicate entries imply multiple paths with the same I/O difference.

- **`get_DT_IO`**  
  Returns all trails matching the input-output pair with EDP ≥ `p + gap`.  
  Used to quantify clustering effects.

---

## 🧪 Example Script

A working example (based on **FLY**) is provided in the following Jupyter notebook:

```
Python/Differential-attack.ipynb
```


---

## ✉️ Contact

For any questions regarding the implementation or analysis, please contact the authors at:
> **cmcom35@korea.ac.kr**
