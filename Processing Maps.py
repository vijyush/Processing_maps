#!/usr/bin/env python
# coding: utf-8

# ### IN 600

# In[65]:


import numpy as np
import matplotlib.pyplot as plt
from numpy import polyfit, polyval

# Data
strain_rate = np.array([0.001, 0.01, 0.1, 1, 10, 100])
temperatures = np.array([850, 900, 950, 1000, 1050, 1100, 1150, 1200])
stress = np.array([
    [181.0, 143.0, 106.0, 83.4, 66.6, 49.5, 36.6, 28.5],
    [275.0, 225.0, 175.0, 128.0, 99.4, 72.5, 59.8, 26.1],
    [428.0, 295.0, 241.0, 204.0, 150.0, 120.0, 94.2, 47.6],
    [526.0, 433.0, 351.0, 303.0, 232.0, 196.0, 156.0, 124.0],
    [841.0, 580.0, 470.0, 411.0, 336.0, 303.0, 217.0, 181.0],
    [799.0, 623.0, 529.0, 482.0, 430.0, 379.0, 316.0, 246.0]
])

# Log-transformed strain rate
strainlog = np.log10(strain_rate)
strainlog_fine = np.linspace(np.min(strainlog), np.max(strainlog), 200)

# Plot setup
plt.figure(figsize=(10, 6))

# Loop over all temperatures
for k in range(len(temperatures)):
    stress_at_temp = stress[:, k]
    
    # Fit 3rd-degree polynomial in log-log space
    coeffs = polyfit(strainlog, np.log10(stress_at_temp), 3)
    stresslog_fit = polyval(coeffs, strainlog_fine)
    
    # Plot fit
    plt.plot(strainlog_fine, stresslog_fit, label=f'{temperatures[k]}°C')

# Scatter actual data points for reference (optional)
for k in range(len(temperatures)):
    plt.scatter(strainlog, np.log10(stress[:, k]), s=20)

# Final plot settings
plt.title('3rd-Degree Polynomial Fits: log(Stress) vs log(Strain Rate)')
plt.xlabel('log10(Strain Rate)')
plt.ylabel('log10(Stress)')
plt.legend(title='Temperature')
plt.grid(True)
plt.tight_layout()
plt.show()


# In[53]:


import numpy as np
import matplotlib.pyplot as plt
from numpy import polyfit, polyder, polyval


strain_rate = np.array([0.001, 0.01, 0.1, 1, 10, 100])  # Strain rates
temperatures = np.array([850, 900, 950, 1000, 1050, 1100, 1150, 1200])  # Temperatures (°C)
stress = np.array([
    [181.0, 143.0, 106.0, 83.4, 66.6, 49.5, 36.6, 28.5],
    [275.0, 225.0, 175.0, 128.0, 99.4, 72.5, 59.8, 26.1],
    [428.0, 295.0, 241.0, 204.0, 150.0, 120.0, 94.2, 47.6],
    [526.0, 433.0, 351.0, 303.0, 232.0, 196.0, 156.0, 124.0],
    [841.0, 580.0, 470.0, 411.0, 336.0, 303.0, 217.0, 181.0],
    [799.0, 623.0, 529.0, 482.0, 430.0, 379.0, 316.0, 246.0]
])

# Logarithmic transformation of strain rate
strainlog = np.log10(strain_rate)

# Create an empty array to store the strain rate sensitivity values
m = np.zeros((len(stress), len(stress[0])))

# Loop to calculate strain rate sensitivity for each temperature
for k in range(len(stress[0])):  # Loop over columns (temperatures)
    for p in range(len(strain_rate)):  # Loop over strain rate points
        # Fit a third-degree polynomial to the log of strain rate vs log of stress
        c1 = polyfit(np.log10(strain_rate), np.log10(stress[:, k]), 3)
        
        # Compute the derivative of the polynomial
        deriv = polyder(c1)
        
        # Evaluate the derivative at the logarithmic strain rate point
        m[p, k] = polyval(deriv, strainlog[p])

# Calculate 'e' and instability zone 'iz'
e = (2 * m) / (m + 1) * 100  # e calculation
e1 = e / 200  # Scale e for instability zone calculation

# Ensure no negative or zero values before applying log10
e1 = np.maximum(e1, 1e-10)  # Replace non-positive values with a small positive value

# Initialize slope matrix
slope = np.zeros_like(m)

# Loop to calculate slope for instability zone
for k in range(len(stress[0])):  # Loop over columns (temperatures)
    for p in range(len(strain_rate)):  # Loop over strain rate points
        # Fit a third-degree polynomial to the log of strain rate vs e1
        c2 = polyfit(np.log10(strain_rate), np.log10(e1[:, k]), 3)
        
        # Compute the derivative of the polynomial
        deriv = polyder(c2)
        
        # Evaluate the derivative at the logarithmic strain rate point
        slope[p, k] = polyval(deriv, strainlog[p])

# Calculate instability zone (iz)
iz = slope + m

# Create the contour plot with larger figure size (matching the 'e' plot size)
plt.figure(figsize=(12, 8))  # Increased figsize

# Correct meshgrid for matching iz shape
X, Y = np.meshgrid(temperatures, strainlog)  # Y should correspond to strain_rate

# Filled contour plot for Instability Zone (IZ) with shading for unstable regions (below 0)
contour_filled = plt.contourf(X, Y, iz, levels=np.linspace(iz.min(), 0, 15), cmap='coolwarm')

# Add a colorbar for contour levels
plt.colorbar(contour_filled, label="Instability Zone (IZ)")

# Customize plot
plt.xlabel("Temperature (°C)")
plt.ylabel("Log(Strain Rate) (s^-1)")
plt.title(" Instability Zone for IN  600 at strain=0.5 ")

# Show the plot
plt.show()


# In[54]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import SmoothBivariateSpline
from numpy import polyfit, polyder, polyval

# --- Step 1: Define your original data ---
temperatures = np.array([850, 900, 950, 1000, 1050, 1100, 1150, 1200])  # Temperature (°C)
strain_rate = np.array([0.001, 0.01, 0.1, 1, 10, 100])  # Strain Rate (1/s)
strainlog = np.log10(strain_rate)

# Stress data
stress = np.array([
    [181.0, 143.0, 106.0, 83.4, 66.6, 49.5, 36.6, 28.5],
    [275.0, 225.0, 175.0, 128.0, 99.4, 72.5, 59.8, 26.1],
    [428.0, 295.0, 241.0, 204.0, 150.0, 120.0, 94.2, 47.6],
    [526.0, 433.0, 351.0, 303.0, 232.0, 196.0, 156.0, 124.0],
    [841.0, 580.0, 470.0, 411.0, 336.0, 303.0, 217.0, 181.0],
    [799.0, 623.0, 529.0, 482.0, 430.0, 379.0, 316.0, 246.0]
])

# --- Step 2: Calculate m (strain rate sensitivity) and efficiency (e) ---
m = np.zeros((len(strain_rate), len(temperatures)))

for k in range(len(temperatures)):  # Loop over temperatures
    c1 = polyfit(np.log10(strain_rate), np.log10(stress[:, k]), 3)
    deriv = polyder(c1)
    for p in range(len(strain_rate)):
        m[p, k] = polyval(deriv, strainlog[p])

e = (2 * m) / (1 + m) * 100

# --- Step 3: Smooth and interpolate data for finer grid ---
xData, yData = np.meshgrid(temperatures, strainlog)
xData = xData.flatten()
yData = yData.flatten()
zData = e.flatten()

# Smooth the efficiency data
spline_e = SmoothBivariateSpline(xData, yData, zData, s=0.5)

# Create a fine meshgrid
temp_dense = np.linspace(temperatures.min(), temperatures.max(), 100)
strainlog_dense = np.linspace(strainlog.min(), strainlog.max(), 100)
X_dense, Y_dense = np.meshgrid(temp_dense, strainlog_dense)

# Evaluate smoothed efficiency
Z_dense = spline_e.ev(X_dense.ravel(), Y_dense.ravel()).reshape(X_dense.shape)

# --- Step 4: Plot the Efficiency Map ---
fig, ax = plt.subplots(figsize=(12, 8))  # Increased figsize to match IZ plot

# Efficiency contour lines
contour = ax.contour(X_dense, Y_dense, Z_dense, levels=np.linspace(0, 50, 11), cmap='jet')
ax.clabel(contour, inline=True, fontsize=10, fmt='%2.1f')

# Axes labels and title
ax.set_xlabel('Temperature (°C)', fontsize=12)
ax.set_ylabel('Log(Strain Rate) (s⁻¹)', fontsize=12)
ax.set_title(' Efficiency Contours for IN 600 at strain=0.5', fontsize=14)
ax.grid(True)

# Colorbar
fig.colorbar(contour, ax=ax, label='Efficiency (e)', shrink=0.75)

plt.show()


# In[55]:


import numpy as np
import matplotlib.pyplot as plt
import cv2
from skimage import io, img_as_float
from skimage.color import rgb2gray
from skimage.transform import resize
from skimage.exposure import equalize_hist

# Read images
fig1 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 112054.png")  # Efficiency contour image
fig2 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 112109.png")  # Instability zone image

# Fix: Remove alpha channel if present
if fig1.shape[-1] == 4:
    fig1 = fig1[..., :3]
if fig2.shape[-1] == 4:
    fig2 = fig2[..., :3]

# Resize images
fig1_resized = resize(fig1, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)
fig2_resized = resize(fig2, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)

# Normalize images
fig1_norm = img_as_float(fig1_resized)
fig2_norm = img_as_float(fig2_resized)

# Superimpose with false color effect
fusion = np.zeros_like(fig1_norm)
fusion[..., 0] = fig2_norm[..., 0]  # Red channel: instability
fusion[..., 1] = fig1_norm[..., 1]  # Green channel: efficiency
fusion[..., 2] = (fig1_norm[..., 2] + fig2_norm[..., 2]) / 2  # Blue: average

# Convert to grayscale
fusion_gray = rgb2gray(fusion)

# Step 1: Apply histogram equalization to improve overall contrast
fusion_eq = equalize_hist(fusion_gray)  # Enhance overall contrast

# Step 2: Further brighten the instability region and darken the contours if necessary
fusion_eq = np.clip(fusion_eq, 0, 0.95)  # Avoid clipping the highlights too much

# Final resize if needed
final_image = resize(fusion_eq, (1968, 2622), preserve_range=True, anti_aliasing=True)

# --- Step 1: Define your original data ---
temperatures = np.array([850, 900, 950, 1000, 1050, 1100, 1150, 1200])  # Temperature (°C)
strain_rate = np.array([0.001, 0.01, 0.1, 1, 10, 100])  # Strain Rate (1/s)

# Logarithmic strain rates
strainlog = np.log10(strain_rate)

# Round off the temperatures and logarithmic strain rate to two decimal places
temperatures_rounded = np.round(temperatures, 2)
strainlog_rounded = np.round(strainlog, 2)

# Plotting
plt.figure(figsize=(10, 12))
plt.imshow(final_image, cmap='gray', vmin=0, vmax=1)  # Ensure proper brightness range

# Set x and y axis labels
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Log(Strain Rate) (s⁻¹)', fontsize=14)

# Set x and y ticks (for proper scaling)
plt.xticks(np.linspace(0, final_image.shape[1], len(temperatures_rounded)), temperatures_rounded)
plt.yticks(np.linspace(0, final_image.shape[0], len(strainlog_rounded)), strainlog_rounded)

# Add title
plt.title('Processing Map of IN-600 at a Strain of 0.5', fontsize=16)

# Hide or show axis grid
plt.axis('on')  # or 'off' to hide axes

# Display the plot
plt.show()


# In[ ]:





# ### IN 625

# In[57]:


import numpy as np
import matplotlib.pyplot as plt
from numpy import polyfit, polyder, polyval

# Define the data
strain_rate = np.array([0.001, 0.01, 0.1, 1, 10, 100])  # Strain rates (s^-1)
temperatures = np.array([900, 950, 1000, 1050, 1100, 1150, 1200])  # Temperatures (°C)
stress = np.array([
    [233.0, 153.9, 103.1, 72.4, 48.7, 31.8, 22.9],
    [278.7, 242.8, 165.5, 117.1, 90.4, 71.4, 47.6],
    [525.6, 382.1, 256.1, 182.0, 139.4, 118.2, 86.2],
    [748.9, 568.4, 416.2, 264.3, 205.2, 155.5, 128.5],
    [1214.4, 865.4, 633.3, 429.0, 302.9, 229.0, 200.7],
    [1040.3, 891.1, 772.7, 633.1, 492.7, 413.2, 272.4]
])

# Logarithmic transformation of strain rate
strainlog = np.log10(strain_rate)

# Create an empty array to store the strain rate sensitivity values
m = np.zeros((len(stress), len(stress[0])))

# Loop to calculate strain rate sensitivity for each temperature
for k in range(len(stress[0])):  # Loop over columns (temperatures)
    for p in range(len(strain_rate)):  # Loop over strain rate points
        # Fit a third-degree polynomial to the log of strain rate vs log of stress
        c1 = polyfit(np.log10(strain_rate), np.log10(stress[:, k]), 3)
        
        # Compute the derivative of the polynomial
        deriv = polyder(c1)
        
        # Evaluate the derivative at the logarithmic strain rate point
        m[p, k] = polyval(deriv, strainlog[p])

# Calculate 'e' and instability zone 'iz'
e = (2 * m) / (m + 1) * 100  # e calculation
e1 = e / 200  # Scale e for instability zone calculation

# Ensure no negative or zero values before applying log10
e1 = np.maximum(e1, 1e-10)  # Replace non-positive values with a small positive value

# Initialize slope matrix
slope = np.zeros_like(m)

# Loop to calculate slope for instability zone
for k in range(len(stress[0])):  # Loop over columns (temperatures)
    for p in range(len(strain_rate)):  # Loop over strain rate points
        # Fit a third-degree polynomial to the log of strain rate vs e1
        c2 = polyfit(np.log10(strain_rate), np.log10(e1[:, k]), 3)
        
        # Compute the derivative of the polynomial
        deriv = polyder(c2)
        
        # Evaluate the derivative at the logarithmic strain rate point
        slope[p, k] = polyval(deriv, strainlog[p])

# Calculate instability zone (iz)
iz = slope + m

# Create the contour plot with larger figure size (matching the 'e' plot size)
plt.figure(figsize=(12, 8))  # Increased figsize

# Correct meshgrid for matching iz shape
X, Y = np.meshgrid(temperatures, strainlog)  # Y should correspond to strain_rate

# Filled contour plot for Instability Zone (IZ) with shading for unstable regions (below 0)
contour_filled = plt.contourf(X, Y, iz, levels=np.linspace(iz.min(), 0, 15), cmap='coolwarm')

# Add a colorbar for contour levels
plt.colorbar(contour_filled, label="Instability Zone (IZ)")

# Customize plot
plt.xlabel("Temperature (°C)")
plt.ylabel("Log(Strain Rate) (s^-1)")
plt.title(" Instability Zone (IZ) for IN 625 at Strain=0.5 ")

# Show the plot
plt.show()


# In[58]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import SmoothBivariateSpline
from numpy import polyfit, polyder, polyval

# --- Step 1: Define your new data ---
temperatures = np.array([900, 950, 1000, 1050, 1100, 1150, 1200])  # Temperature (°C)
strain_rate = np.array([0.001, 0.01, 0.1, 1, 10, 100])  # Strain Rate (1/s)
strainlog = np.log10(strain_rate)

# Stress data for Nickel-Base Superalloy IN 625 (PM)
stress = np.array([
    [233.0, 153.9, 103.1, 72.4, 48.7, 31.8, 22.9],
    [278.7, 242.8, 165.5, 117.1, 90.4, 71.4, 47.6],
    [525.6, 382.1, 256.1, 182.0, 139.4, 118.2, 86.2],
    [748.9, 568.4, 416.2, 264.3, 205.2, 155.5, 128.5],
    [1214.4, 865.4, 633.3, 429.0, 302.9, 229.0, 200.7],
    [1040.3, 891.1, 772.7, 633.1, 492.7, 413.2, 272.4]
])

# --- Step 2: Calculate m (strain rate sensitivity) and efficiency (e) ---
m = np.zeros((len(strain_rate), len(temperatures)))

for k in range(len(temperatures)):  # Loop over temperatures
    c1 = polyfit(np.log10(strain_rate), np.log10(stress[:, k]), 3)
    deriv = polyder(c1)
    for p in range(len(strain_rate)):
        m[p, k] = polyval(deriv, strainlog[p])

e = (2 * m) / (1 + m) * 100  # Calculate efficiency map

# --- Step 3: Smooth and interpolate data for finer grid ---
xData, yData = np.meshgrid(temperatures, strainlog)
xData = xData.flatten()
yData = yData.flatten()
zData = e.flatten()

# Smooth the efficiency data using a Bivariate Spline
spline_e = SmoothBivariateSpline(xData, yData, zData, s=0.5)

# Create a fine meshgrid for interpolation
temp_dense = np.linspace(temperatures.min(), temperatures.max(), 100)
strainlog_dense = np.linspace(strainlog.min(), strainlog.max(), 100)
X_dense, Y_dense = np.meshgrid(temp_dense, strainlog_dense)

# Evaluate smoothed efficiency on the fine meshgrid
Z_dense = spline_e.ev(X_dense.ravel(), Y_dense.ravel()).reshape(X_dense.shape)

# --- Step 4: Plot the Efficiency Map ---
fig, ax = plt.subplots(figsize=(12, 8))  # Increased figsize for better visualization

# Efficiency contour lines
contour = ax.contour(X_dense, Y_dense, Z_dense, levels=np.linspace(0, 60, 13), cmap='jet')  # Adjusted levels
ax.clabel(contour, inline=True, fontsize=10, fmt='%2.1f')

# Axes labels and title
ax.set_xlabel('Temperature (°C)', fontsize=12)
ax.set_ylabel('Log(Strain Rate) (s⁻¹)', fontsize=12)
ax.set_title(' Efficiency Contours for IN 625 at strain=0.5', fontsize=14)
ax.grid(True)

# Colorbar to show efficiency values
fig.colorbar(contour, ax=ax, label='Efficiency (e)', shrink=0.75)

plt.show()


# In[43]:


import numpy as np
import matplotlib.pyplot as plt
import cv2
from skimage import io, img_as_float
from skimage.color import rgb2gray
from skimage.transform import resize
from skimage.exposure import equalize_hist

# Read images
fig1 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 111431.png")  # Efficiency contour image
fig2 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 111454.png")  # Instability zone image

# Fix: Remove alpha channel if present
if fig1.shape[-1] == 4:
    fig1 = fig1[..., :3]
if fig2.shape[-1] == 4:
    fig2 = fig2[..., :3]

# Resize images
fig1_resized = resize(fig1, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)
fig2_resized = resize(fig2, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)

# Normalize images
fig1_norm = img_as_float(fig1_resized)
fig2_norm = img_as_float(fig2_resized)

# Superimpose with false color effect
fusion = np.zeros_like(fig1_norm)
fusion[..., 0] = fig2_norm[..., 0]  # Red channel: instability
fusion[..., 1] = fig1_norm[..., 1]  # Green channel: efficiency
fusion[..., 2] = (fig1_norm[..., 2] + fig2_norm[..., 2]) / 2  # Blue: average

# Convert to grayscale
fusion_gray = rgb2gray(fusion)

# Step 1: Apply histogram equalization to improve overall contrast
fusion_eq = equalize_hist(fusion_gray)  # Enhance overall contrast

# Step 2: Further brighten the instability region and darken the contours if necessary
fusion_eq = np.clip(fusion_eq, 0, 0.95)  # Avoid clipping the highlights too much

# Final resize if needed
final_image = resize(fusion_eq, (1968, 2622), preserve_range=True, anti_aliasing=True)

# --- Step 1: Define your new data ---
temperatures = np.array([900, 950, 1000, 1050, 1100, 1150, 1200])  # Temperature (°C)
strain_rate = np.array([0.001, 0.01, 0.1, 1, 10, 100])  # Strain Rate (1/s)

# Logarithmic strain rates
strainlog = np.log10(strain_rate)

# Round off the temperatures and logarithmic strain rate to two decimal places
temperatures_rounded = np.round(temperatures, 2)
strainlog_rounded = np.round(strainlog, 2)

# Plotting
plt.figure(figsize=(10, 12))
plt.imshow(final_image, cmap='gray', vmin=0, vmax=1)  # Ensure proper brightness range

# Set x and y axis labels
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Log(Strain Rate) (s⁻¹)', fontsize=14)

# Set x and y ticks (for proper scaling)
plt.xticks(np.linspace(0, final_image.shape[1], len(temperatures_rounded)), temperatures_rounded)
plt.yticks(np.linspace(0, final_image.shape[0], len(strainlog_rounded)), strainlog_rounded)

# Add title
plt.title('Processing Map of IN-625 at a Strain of 0.5', fontsize=16)

# Hide or show axis grid
plt.axis('on')  # or 'off' to hide axes

# Display the plot
plt.show()


# In[ ]:





# ### IN 718

# In[59]:


# Forged  IN718 (25 μm)
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import SmoothBivariateSpline
from numpy import polyfit, polyder, polyval

# --- Step 1: Define the new data for IN 718 at strain = 0.4 ---
temperatures = np.array([900, 950, 975, 1000, 1025, 1050, 1100])  # Temperature (°C)
strain_rate = np.array([0.001, 0.01, 0.1, 1.0])  # Strain Rate (1/s)
strainlog = np.log10(strain_rate)

# Flow Stress values at 0.4 strain
stress = np.array([
    [306.8, 159.3, 119.8, 100.5,  81.9,  60.5,  44.8],   # 0.001 1/s
    [408.4, 258.5, 193.0, 153.8, 133.3, 118.2, 102.9],   # 0.01 1/s
    [555.9, 348.0, 286.9, 252.1, 220.4, 189.3, 153.4],   # 0.1 1/s
    [686.7, 486.3, 412.5, 342.8, 346.1, 282.5, 215.7],   # 1.0 1/s
])

# --- Step 2: Calculate m (strain rate sensitivity) and efficiency (e) ---
m = np.zeros((len(strain_rate), len(temperatures)))

for k in range(len(temperatures)):  # Loop over temperatures
    c1 = polyfit(np.log10(strain_rate), np.log10(stress[:, k]), 3)  # Fit 3rd degree polynomial
    deriv = polyder(c1)  # Differentiate the polynomial
    for p in range(len(strain_rate)):
        m[p, k] = polyval(deriv, strainlog[p])  # Evaluate the derivative

e = (2 * m) / (1 + m) * 100  # Calculate efficiency map

# --- Step 3: Smooth and interpolate data for finer grid ---
xData, yData = np.meshgrid(temperatures, strainlog)
xData = xData.flatten()
yData = yData.flatten()
zData = e.flatten()

# Smooth the efficiency data using a Bivariate Spline
spline_e = SmoothBivariateSpline(xData, yData, zData, s=0.5)

# Create a fine meshgrid for interpolation
temp_dense = np.linspace(temperatures.min(), temperatures.max(), 100)
strainlog_dense = np.linspace(strainlog.min(), strainlog.max(), 100)
X_dense, Y_dense = np.meshgrid(temp_dense, strainlog_dense)

# Evaluate smoothed efficiency on the fine meshgrid
Z_dense = spline_e.ev(X_dense.ravel(), Y_dense.ravel()).reshape(X_dense.shape)

# --- Step 4: Plot the Efficiency Map ---
fig, ax = plt.subplots(figsize=(12, 8))  # Increased figsize for better visualization

# Efficiency contour lines
contour = ax.contour(X_dense, Y_dense, Z_dense, levels=np.linspace(0, 60, 13), cmap='jet')  # Adjusted levels
ax.clabel(contour, inline=True, fontsize=10, fmt='%2.1f')

# Axes labels and title
ax.set_xlabel('Temperature (°C)', fontsize=12)
ax.set_ylabel('Log(Strain Rate) (s⁻¹)', fontsize=12)
ax.set_title(' Efficiency Contours for IN 718 (Strain = 0.5)', fontsize=14)
ax.grid(True)

# Colorbar to show efficiency values
fig.colorbar(contour, ax=ax, label='Efficiency (e)', shrink=0.75)

plt.show()


# In[60]:


import numpy as np
import matplotlib.pyplot as plt
from numpy import polyfit, polyder, polyval

# --- Step 1: Define the data ---
temperatures = np.array([900, 950, 975, 1000, 1025, 1050, 1100])  # Temperature (°C)
strain_rate = np.array([0.001, 0.01, 0.1, 1.0])  # Strain rate (1/s)
strainlog = np.log10(strain_rate)

# Stress data at 0.4 strain
stress = np.array([
    [306.8, 159.3, 119.8, 100.5,  81.9,  60.5,  44.8],   # 0.001 1/s
    [408.4, 258.5, 193.0, 153.8, 133.3, 118.2, 102.9],   # 0.01 1/s
    [555.9, 348.0, 286.9, 252.1, 220.4, 189.3, 153.4],   # 0.1 1/s
    [686.7, 486.3, 412.5, 342.8, 346.1, 282.5, 215.7],   # 1.0 1/s
])


# --- Step 2: Calculate strain rate sensitivity 'm' ---
m = np.zeros_like(stress)

for k in range(stress.shape[1]):  # Loop over temperatures
    c1 = polyfit(strainlog, np.log10(stress[:, k]), 3)
    deriv = polyder(c1)
    for p in range(len(strain_rate)):
        m[p, k] = polyval(deriv, strainlog[p])

# --- Step 3: Calculate efficiency 'e' and normalized efficiency 'e1' ---
e = (2 * m) / (1 + m) * 100
e1 = e / 200
e1 = np.maximum(e1, 1e-10)  # Avoid log10 issues

# --- Step 4: Calculate slope of log(e1) w.r.t log(strain_rate) ---
slope = np.zeros_like(m)
for k in range(stress.shape[1]):
    c2 = polyfit(strainlog, np.log10(e1[:, k]), 3)
    deriv = polyder(c2)
    for p in range(len(strain_rate)):
        slope[p, k] = polyval(deriv, strainlog[p])

# --- Step 5: Calculate instability zone ---
iz = slope + m

# --- Step 6: Plot the instability zone contour ---
X, Y = np.meshgrid(temperatures, strainlog)

plt.figure(figsize=(12, 8))
contour_filled = plt.contourf(X, Y, iz, levels=np.linspace(iz.min(), 0, 15), cmap='coolwarm')
plt.colorbar(contour_filled, label="Instability Zone (IZ)")

# Customize plot
plt.xlabel("Temperature (°C)")
plt.ylabel("Log(Strain Rate) (s⁻¹)")
plt.title("Instability Zone (IZ) for IN 718 at 0.5 Strain ")
plt.grid(True)
plt.show()


# In[69]:


import numpy as np
import matplotlib.pyplot as plt
import cv2
from skimage import io, img_as_float
from skimage.color import rgb2gray
from skimage.transform import resize
from skimage.exposure import equalize_hist

# Read images
fig1 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 122802.png")  # Efficiency contour image
fig2 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 122819.png")  # Instability zone image



# Fix: Remove alpha channel if present
if fig1.shape[-1] == 4:
    fig1 = fig1[..., :3]
if fig2.shape[-1] == 4:
    fig2 = fig2[..., :3]

# Resize images
fig1_resized = resize(fig1, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)
fig2_resized = resize(fig2, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)

# Normalize images
fig1_norm = img_as_float(fig1_resized)
fig2_norm = img_as_float(fig2_resized)

# Superimpose with false color effect
fusion = np.zeros_like(fig1_norm)
fusion[..., 0] = fig2_norm[..., 0]  # Red channel: instability
fusion[..., 1] = fig1_norm[..., 1]  # Green channel: efficiency
fusion[..., 2] = (fig1_norm[..., 2] + fig2_norm[..., 2]) / 2  # Blue: average

# Convert to grayscale
fusion_gray = rgb2gray(fusion)

# Step 1: Apply histogram equalization to improve overall contrast
fusion_eq = equalize_hist(fusion_gray)  # Enhance overall contrast

# Step 2: Further brighten the instability region and darken the contours if necessary
fusion_eq = np.clip(fusion_eq, 0, 0.95)  # Avoid clipping the highlights too much

# Final resize if needed
final_image = resize(fusion_eq, (1968, 2622), preserve_range=True, anti_aliasing=True)

# Define your new data
temperatures = np.array([900, 950, 975, 1000, 1025, 1050, 1100])  # Temperature (°C)
strain_rate = np.array([0.001, 0.01, 0.1, 1.0])  # Strain rate (1/s)

# Compute logarithm of strain rate for the y-axis
log_strain_rate = np.log10(strain_rate)

# Round the temperatures and logarithmic strain rate to two decimal places
temperatures_rounded = np.round(temperatures, 2)
log_strain_rate_rounded = np.round(log_strain_rate, 2)

# Plotting
plt.figure(figsize=(10, 12))
plt.imshow(final_image, cmap='gray', vmin=0, vmax=1)  # Ensure proper brightness range

# Set x and y axis labels
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Log(Strain Rate) (s⁻¹)', fontsize=14)

# Set x and y ticks (for proper scaling)
plt.xticks(np.linspace(0, final_image.shape[1], len(temperatures_rounded)), temperatures_rounded)
plt.yticks(np.linspace(0, final_image.shape[0], len(log_strain_rate_rounded)), log_strain_rate_rounded)

# Add title
plt.title('Processing Map of Forged  IN718 (25 μm)  at a Strain of 0.5', fontsize=16)

# Hide or show axis grid
plt.axis('on')  # or 'off' to hide axes

# Display the plot
plt.show()


# ### IN 100

# In[61]:


import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import SmoothBivariateSpline
from numpy import polyfit, polyder, polyval

# --- Step 1: Define the data for IN 718 at strain = 0.3 ---
temperatures = np.array([1000, 1050, 1100, 1150, 1200])  # Temperature (°C)
strain_rate = np.array([0.0003, 0.001, 0.01, 0.1, 1.0, 10.0])  # Strain rate (1/s)
strainlog = np.log10(strain_rate)

stress = np.array([
    [154,  82,  60,  40,  20],   # 0.0003 1/s
    [173,  82,  74,  52,  22],   # 0.001 1/s
    [286, 175, 122,  72,  36],   # 0.01 1/s
    [337, 226, 172,  82,  53],   # 0.1 1/s
    [438, 353, 258, 163,  83],   # 1.0 1/s
    [455, 402, 362, 278, 131],   # 10.0 1/s
])


# --- Step 2: Calculate m (strain rate sensitivity) and efficiency (e) ---
m = np.zeros((len(strain_rate), len(temperatures)))

for k in range(len(temperatures)):  # Loop over temperatures
    c1 = polyfit(strainlog, np.log10(stress[:, k]), 3)  # Fit 3rd degree polynomial
    deriv = polyder(c1)  # Differentiate the polynomial
    for p in range(len(strain_rate)):
        m[p, k] = polyval(deriv, strainlog[p])  # Evaluate the derivative

e = (2 * m) / (1 + m) * 100  # Calculate efficiency map

# --- Step 3: Smooth and interpolate data for finer grid ---
xData, yData = np.meshgrid(temperatures, strainlog)
xData = xData.flatten()
yData = yData.flatten()
zData = e.flatten()

# Smooth the efficiency data using a Bivariate Spline
spline_e = SmoothBivariateSpline(xData, yData, zData, s=0.5)

# Create a fine meshgrid for interpolation
temp_dense = np.linspace(temperatures.min(), temperatures.max(), 100)
strainlog_dense = np.linspace(strainlog.min(), strainlog.max(), 100)
X_dense, Y_dense = np.meshgrid(temp_dense, strainlog_dense)

# Evaluate smoothed efficiency on the fine meshgrid
Z_dense = spline_e.ev(X_dense.ravel(), Y_dense.ravel()).reshape(X_dense.shape)

# --- Step 4: Plot the Efficiency Map ---
fig, ax = plt.subplots(figsize=(12, 8))  # Increased figsize for better visualization

# Efficiency contour lines
contour = ax.contour(X_dense, Y_dense, Z_dense, levels=np.linspace(0, 60, 13), cmap='jet')  # Adjusted levels
ax.clabel(contour, inline=True, fontsize=10, fmt='%2.1f')

# Axes labels and title
ax.set_xlabel('Temperature (°C)', fontsize=12)
ax.set_ylabel('Log(Strain Rate) (s⁻¹)', fontsize=12)
ax.set_title(' Efficiency Contours for IN 100 (Strain = 0.5)', fontsize=14)
ax.grid(True)

# Colorbar to show efficiency values
fig.colorbar(contour, ax=ax, label='Efficiency (e)', shrink=0.75)

plt.show()


# In[62]:


import numpy as np
import matplotlib.pyplot as plt
from numpy import polyfit, polyder, polyval
from matplotlib.colors import LinearSegmentedColormap

# IN-100 data (True strain = 0.3)
temperatures = np.array([1000, 1050, 1100, 1150, 1200])  # °C
strain_rate = np.array([0.0003, 0.001, 0.01, 0.1, 1.0, 10.0])  # 1/s

# Flow stress at 0.3 strain for each strain rate and temperature
stress = np.array([
    [154,  82,  60,  40,  20],   # 0.0003 1/s
    [173,  82,  74,  52,  22],   # 0.001 1/s
    [286, 175, 122,  72,  36],   # 0.01 1/s
    [337, 226, 172,  82,  53],   # 0.1 1/s
    [438, 353, 258, 163,  83],   # 1.0 1/s
    [455, 402, 362, 278, 131],   # 10.0 1/s
])


# Logarithmic strain rates
strainlog = np.log10(strain_rate)

# Strain rate sensitivity (m) matrix
m = np.zeros_like(stress, dtype=float)

# Fit log-log stress vs strain rate with cubic polynomial and compute m
for j in range(stress.shape[1]):  # Over temperatures
    log_stress = np.log10(stress[:, j])
    coeffs = polyfit(strainlog, log_stress, 3)
    deriv = polyder(coeffs)
    m[:, j] = polyval(deriv, strainlog)

# Power dissipation efficiency
e = (2 * m) / (m + 1) * 100
e1 = np.maximum(e / 200, 1e-10)  # Normalize and clip to avoid log(0)

# Instability slope matrix
slope = np.zeros_like(m)
for j in range(stress.shape[1]):
    log_e = np.log10(e1[:, j])
    coeffs2 = polyfit(strainlog, log_e, 3)
    deriv2 = polyder(coeffs2)
    slope[:, j] = polyval(deriv2, strainlog)

# Instability criterion
iz = slope + m  # Negative values indicate instability

# Plotting
plt.figure(figsize=(12, 8))
X, Y = np.meshgrid(temperatures, strainlog)

# Create custom colormap: Red for instability (ξ < 0), White for stability (ξ > 0)
cmap = LinearSegmentedColormap.from_list("instability_map", [(1, 0, 0), (1, 1, 1)])  # Red to White

# Define levels: We want to highlight instability region in red (ξ < 0) and stable region in white (ξ > 0)
levels = np.linspace(np.min(iz), 0, 100)  # From the minimum instability value (ξ < 0) to 0

# Plot the instability region (ξ < 0) shaded red and stability region (ξ > 0) white
contour = plt.contourf(X, Y, iz, levels=levels, cmap=cmap, extend='both')

# Add colorbar
plt.colorbar(contour, label="Instability Parameter (ξ)")

# Labels and plot formatting
plt.xlabel("Temperature (°C)")
plt.ylabel("Log(Strain Rate) (s⁻¹)")
plt.title("Instability Map for IN-100 Superalloy (True Strain = 0.5)")
plt.grid(True)
plt.tight_layout()
plt.show()


# In[51]:


import numpy as np
import matplotlib.pyplot as plt
import cv2
from skimage import io, img_as_float
from skimage.color import rgb2gray
from skimage.transform import resize
from skimage.exposure import equalize_hist

# Read images
fig1 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 123106.png")  # Efficiency contour image
fig2 = io.imread(r"C:\Users\burra\OneDrive\Pictures\Screenshots\Screenshot 2025-04-29 123121.png")  # Instability zone image

# Fix: Remove alpha channel if present
if fig1.shape[-1] == 4:
    fig1 = fig1[..., :3]
if fig2.shape[-1] == 4:
    fig2 = fig2[..., :3]

# Resize images
fig1_resized = resize(fig1, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)
fig2_resized = resize(fig2, (1968, 2622), preserve_range=True, anti_aliasing=True).astype(np.uint8)

# Normalize images
fig1_norm = img_as_float(fig1_resized)
fig2_norm = img_as_float(fig2_resized)

# Superimpose with false color effect
fusion = np.zeros_like(fig1_norm)
fusion[..., 0] = fig2_norm[..., 0]  # Red channel: instability
fusion[..., 1] = fig1_norm[..., 1]  # Green channel: efficiency
fusion[..., 2] = (fig1_norm[..., 2] + fig2_norm[..., 2]) / 2  # Blue: average

# Convert to grayscale
fusion_gray = rgb2gray(fusion)

# Step 1: Apply histogram equalization to improve overall contrast
fusion_eq = equalize_hist(fusion_gray)  # Enhance overall contrast

# Step 2: Further brighten the instability region and darken the contours if necessary
fusion_eq = np.clip(fusion_eq, 0, 0.95)  # Avoid clipping the highlights too much

# Final resize if needed
final_image = resize(fusion_eq, (1968, 2622), preserve_range=True, anti_aliasing=True)

# Define your data
temperatures = np.array([1000, 1050, 1100, 1150, 1200])  # °C
strain_rate = np.array([0.0003, 0.001, 0.01, 0.1, 1.0, 10.0])  # 1/s

# Compute logarithm of strain rate for the y-axis
log_strain_rate = np.log10(strain_rate)

# Round the temperatures and logarithmic strain rate to two decimal places
temperatures_rounded = np.round(temperatures, 2)
log_strain_rate_rounded = np.round(log_strain_rate, 2)

# Plotting
plt.figure(figsize=(10, 12))
plt.imshow(final_image, cmap='gray', vmin=0, vmax=1)  # Ensure proper brightness range

# Set x and y axis labels
plt.xlabel('Temperature (°C)', fontsize=14)
plt.ylabel('Log(Strain Rate) (s⁻¹)', fontsize=14)

# Set x and y ticks (for proper scaling)
plt.xticks(np.linspace(0, final_image.shape[1], len(temperatures_rounded)), temperatures_rounded)
plt.yticks(np.linspace(0, final_image.shape[0], len(log_strain_rate_rounded)), log_strain_rate_rounded)

# Add title
plt.title('Processing Map of IN-100 at a Strain of 0.5', fontsize=16)

# Hide or show axis grid
plt.axis('on')  # or 'off' to hide axes

# Display the plot
plt.show()


# In[ ]:




