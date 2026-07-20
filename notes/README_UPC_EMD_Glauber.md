# UPC EMD Glauber Probability Macro

`calc_UPC_EMD_Glauber_unified.C` generates impact-parameter probability tables for
UPC calculations. The output ROOT file is used by the rho-prime grid generator to
apply hadronic survival and EMD tag probabilities.

## What It Computes

The macro computes:

- hadronic survival probability, `P_surv(b)`
- single-side EMD probabilities, `P_1n(b)` and `P_Xn(b)`
- mutual EMD probabilities, `P_1n1n(b)` and `P_XnXn(b)`
- survival-weighted SCD/MCD products for compatibility with older workflows

PbPb and AuAu parameters are built into the macro, including Woods-Saxon
parameters and the photonuclear cross-section tables used for EMD.

## Requirements

Run with ROOT:

```bash
root -l -b -q 'calc_UPC_EMD_Glauber_unified.C("PbPb",5360,68,"UPC_Probabilities_PbPb_5360GeV_Glauber_Final.root")'
```

No external data files are required.

## Function Signature

```cpp
void calc_UPC_EMD_Glauber_unified(
    const char* system = "PbPb",
    double sqrt_s_NN_GeV = 5360.0,
    double sigma_NN_mb = 68.0,
    const char* outputFileName = ""
)
```

Parameters:

- `system`: `"PbPb"`, `"Pb"`, `"AuAu"`, or `"Au"`
- `sqrt_s_NN_GeV`: nucleon-nucleon center-of-mass energy in GeV
- `sigma_NN_mb`: inelastic nucleon-nucleon cross section in mb
- `outputFileName`: output ROOT file path. If empty, a default name is used.

For AuAu, if `sigma_NN_mb` is left at the default `68.0`, the macro internally
uses `42.0 mb`.

## Common Commands

Generate the PbPb 5.36 TeV table:

```bash
root -l -b -q 'calc_UPC_EMD_Glauber_unified.C("PbPb",5360,68,"UPC_Probabilities_PbPb_5360GeV_Glauber_Final.root")'
```

Generate the AuAu 200 GeV table:

```bash
root -l -b -q 'calc_UPC_EMD_Glauber_unified.C("AuAu",200,42,"UPC_Probabilities_AuAu_200GeV_Glauber_Final.root")'
```

## Output Histograms

The ROOT file contains:

```text
hProb_Surv
hProb_EMD_1n
hProb_EMD_Xn
hProb_EMD_1n1n
hProb_EMD_XnXn
hProb_SCD_1n
hProb_SCD_Xn
hProb_MCD_1n1n
hProb_MCD_XnXn
```

The rho-prime grid generator reads:

- `hProb_Surv`
- `hProb_EMD_<tag>` for tagged runs, for example `hProb_EMD_XnXn`

For `NoTag`, only `hProb_Surv` is used.

## Notes

The macro prints integrated MCD and SCD cross sections in barns. The histograms
store probabilities as functions of impact parameter `b` in fm.
