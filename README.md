# Processing_Maps

# Dynamic Material Model (DMM) Processing Maps for Nickel-Based Superalloys

This Jupyter Notebook computes and visualizes **Dynamic Material Model (DMM) Processing Maps** for four Nickel-based superalloys:

- **Inconel 600 (IN 600)**
- **Inconel 625 (IN 625)**
- **Inconel 718 (IN 718)**
- **Inconel 100 (IN 100)**

Processing Maps integrate the **Power Dissipation Efficiency (η)** and **Flow Instability Parameter (ξ)** to identify safe and unsafe hot-working domains across temperature and strain-rate space.

The methodology follows the classical DMM framework and uses experimental flow stress data at multiple temperatures and strain rates.

---

## Project Goals

1. **Calculate DMM Parameters**  
   - Strain-rate sensitivity index (**m**)  
   - Power dissipation efficiency (**η**)  
   - Flow instability parameter (**ξ**)  

2. **Generate Contour Plots**  
   - η-maps (efficiency)  
   - ξ-maps (instability)  

3. **Superimpose Processing Maps**  
   - Combine η and ξ contours to visualize stable (ξ ≥ 0) and unstable (ξ < 0) deformation domains.

---

## Methodology (Dynamic Material Model)

### 1. Strain-Rate Sensitivity (m)

\[
m = \frac{\partial \log \sigma}{\partial \log \dot{\epsilon}}
\]

Computed by fitting a **third-degree polynomial** to  
`log(stress)` vs. `log(strain rate)` at constant temperature.

---

### 2. Power Dissipation Efficiency (η)

\[
\eta = \frac{2m}{m+1} \times 100\%
\]

Represents the material’s ability to dissipate energy via microstructural mechanisms  
(e.g., dynamic recovery, dynamic recrystallization).  
Higher η values (> 30%) typically indicate favorable hot-working conditions.

---

### 3. Flow Instability Parameter (ξ)

\[
\xi = m + \frac{\partial \log\left(\frac{m}{m+1}\right)}{\partial \log \dot{\epsilon}}
\]

A negative ξ value (**ξ < 0**) indicates flow instability zones associated with:

- Flow localization  
- Shear banding  
- Intergranular cracking  

Stable domains are regions where **ξ ≥ 0**.

---

## Notebook Structure

The notebook is organized alloy-wise. Each section includes:

- Polynomial fitting for m-calculation  
- Efficiency (η) contour map  
- Instability (ξ) contour map  
- Final Processing Map using superposition  

---

## I. INCONEL 600 (IN 600)

| Cell | Description |
|------|-------------|
| **In[65]** | Polynomial fits: log(σ) vs. log(ε̇) for all temperatures using cubic polynomials. |
| **In[53]** | Instability parameter (ξ): computed and visualized using contour plots. |
| **In[54]** | Efficiency (η) map generated using `SmoothBivariateSpline`. |
| **In[55]** | Final Processing Map by superimposing η and ξ maps using `skimage`. |

---

## II. INCONEL 625 (IN 625)

| Cell | Description |
|------|-------------|
| **In[57]** | Instability parameter (ξ) for IN 625. |
| **In[58]** | Efficiency (η) contour map. |
| **In[43]** | Superimposed Processing Map for IN 625. |

---

## III. INCONEL 718 (IN 718)

(Flow stress at engineering strain ε = 0.4)

| Cell | Description |
|------|-------------|
| **In[59]** | Efficiency (η) contour map. |
| **In[60]** | Instability (ξ) contour map. |
| **In[69]** | Final Processing Map. |

---

## IV. INCONEL 100 (IN 100)

(Flow stress at engineering strain ε = 0.3)

| Cell | Description |
|------|-------------|
| **In[61]** | Efficiency (η) map. |
| **In[62]** | Instability (ξ) contour map. |
| **In[51]** | Superimposed Processing Map. |

---

## Dependencies

```bash
numpy
pandas
matplotlib
scipy
skimage
