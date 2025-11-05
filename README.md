# V10_iScaling — Multi-Dataset 1D NMR Plotter with Per-Trace Scaling 📈
Written by Uno Bolin Hartl (with chatGPT 5.0), BSc candidate Molecular Biology (Stockholm University, MBW) and Research Assistant
in Katja Petzold Lab (RNA Dynamics by NMR - Uppsala University)

This script loads **TopSpin-processed 1D spectra** (`pdata/<procno>`) and generates **overlay** or **stacked** comparison plots.  
It adds **per-trace intensity scaling**, allowing visual normalization across datasets with different experimental conditions.

The script is useful for:
- plot classic 1D spectra
- Comparing titration points
- Monitoring spectral shifts
- Visualizing signal intensity trends
- Generating clean *publication-ready* spectral overlays

---
##  Requirements 🚧
- Python 3.8+
- Install dependencies: pandas, matplotlib, and numpy.
- Alternatively use commands below for discription.
  ```bash
  pip install -r requirements.txt
---

## ✨ Key Features

| Feature | Description |
|--------|-------------|
| Works directly on **Bruker pdata** directories | No conversion required |
| Overlay or stacked plotting modes | `MULTI_MODE = "overlay"` or `"stack"` |
| **Per-trace scale control** | `SCALE_FACTORS` list or global `SCALE_FACTOR` |
| Adjustable ppm window | `PPM_MIN`, `PPM_MAX` |
| Optional colormap scaling | Intensity-based color grading |
| Stack direction and spacing configurable | Vertical offsets managed cleanly |
| Automatic or fixed Y-axis limits | Consistent figure comparison |
| Safe filename templating | Prevents overwriting unless allowed |

---

## 🔧 Example Configuration (edit inside script)

```python
CONFIG = {
    "INPUT_PATHS": [
        r"/path/to/sampleA/pdata/1",
        r"/path/to/sampleB/pdata/1",
        r"/path/to/sampleC/pdata/1",
    ],
    "MULTI_MODE": "stack",     # "overlay" or "stack"
    "SCALE_FACTORS": [1, 1, 4],
    "SCALE_FACTORS_MODE": "multiply",
    "PPM_MIN": 14,
    "PPM_MAX": 9.5,
    "HIDE_Y_AXIS": True,
    "USE_COLORMAP": True,
}
