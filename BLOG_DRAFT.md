# Can satellites predict soil moisture beyond seasonality?

## Lessons from Sentinel-1, Sentinel-2 and SMAP in Bagré, Burkina Faso

Ibrahim Kabore | November 2026 | Independent exploratory study

In a seasonal agricultural landscape, a high soil-moisture prediction score can be misleading. A model that knows the month of the year can predict much of the wet-dry cycle without seeing an image. Can radar and optical imagery provide additional information, especially spatial variability within a month?

### Study design

I used monthly observations from January 2024 through December 2025, Sentinel-1 ascending relative orbit 74 (VV/VH), Sentinel-2 NDVI, ESA WorldCover 2021, and quality-screened SMAP L3 enhanced 9 km AM surface soil moisture. Radar power was averaged in linear units. Fine predictors were aggregated at 250 m and then to the actual SMAP grid. Of 2,880 potential core cell-month rows, 1,750 passed the common-support and reference criteria. August 2024 had only one eligible row. The training, validation, spatial, temporal and joint sets had 463, 107, 263, 629 and 288 rows respectively. These are correlated observations, not independent field samples.

### What did we find?

The seasonal harmonic ridge outperformed all satellite-only models. Joint-test RMSE was 0.0342 m³/m³ (R² 0.8142); SAR ridge yielded 0.0628 and optical ridge 0.0631. A calendar-only prediction is constant within each month, so a high R² cannot establish spatial skill.

The exploratory extension fitted satellite residual models to the same training set. The validation-selected NDVI residual model improved joint RMSE to 0.0308 (R² 0.8489; 18.7% squared-error reduction). Fusion residual yielded 0.0290 (R² 0.8666; 28.2% reduction) but was not selected by validation. These are post-hoc results on already inspected holdouts.

### The spatial test changes the story

Removing each evaluated month's mean reveals whether predictions discriminate locations. The selected NDVI residual model had joint within-month skill +7.5%, temporal +5.6%, spatial −10.9%. SAR and fusion residual models had negative centered skill across all three test splits. The anomaly scatterplots show shrunken predictions, especially for large deviations.

### Maps and uncertainty

January and June 2025 were selected by predictor coverage. The QGIS maps juxtapose the uniform seasonal estimate, paired SMAP observations, differences and a dimensionless, scene-standardized 250 m linear score restricted to crop-dominated pixels. The latter is not a fine-resolution moisture retrieval. White areas are masked or unsupported. SMAP is not field ground truth and its effective radiometer footprint exceeds the posted grid spacing. The original seasonal model's validation residual envelope ±0.0535 m³/m³ covered 76.8%, 86.0%, and 87.5% on spatial, temporal and joint tests; it is descriptive, not a guaranteed interval.

### Figure guide

* `01_dielectric.png`: conceptual dielectric reflection, not radar backscatter.
* `01_toy_canopy.png`: illustrative attenuation and canopy scattering.
* `02_smoke_test.png`: first-month sensor coverage, grid contrast and land-cover context.
* `03_training_scatter.png`: seasonal confounding in VV, VH and NDVI.
* `04_heldout_scatter.png`: calendar-only predictions form horizontal bands.
* `05_maps_2025-01.png`, `05_maps_2025-06.png`: dry/wet maps with relative score.
* `06_monthly_rmse.png`: rainy-season error peaks.
* `07_anomaly_skill.png` / PDF: overall versus within-month skill.
* `07_anomaly_scatter.png` / PDF: residual prediction shrinkage.

### Conclusion

Seasonality must be a strong baseline for agricultural remote sensing. Satellite predictors provided incremental overall skill, but fine spatial discrimination remains unproven. Independent future-year or field observations, roughness and rainfall information, and better wet-season optical availability are priorities for follow-up. No root-zone or irrigation prescription claims.
