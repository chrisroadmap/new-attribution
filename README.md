# new-attribution
Attribution of 2023 and 2024 warming

## emissions counterfactual
new dataset prepared for ScenarioMIP CMIP7, with constant emissions of all GHGs and SLCFs from 2019 to 2024.

## solar maximum
counterfactual: IGCC with zero solar forcing from the point where the ERF crosses zero, circa 2022
experiment: IGCC

## Hunga Tonga
counterfactual: IGCC with HTHH removed in 2023, 2024 and 2025 timebounds
- stratospheric water vapour easy to back out as is present in IGCC time series and nothing else really affected stratospheric water vapour
- sulphate more difficult, but I have a method to isolate it from the total SAOD time series
  - look at monthly SAOD prepared for IGCC, from GloSSAC, and assume that HTHH was only significant eruption between Dec 2021 and Apr 2022
  - we see that SAOD increased from 0.008 to 0.015 (about 0.007 units) peaking in Apr 2022
  - assume that SAOD anomaly decays away with an e-folding lifetime of 24 months starting from t=0 in Apr 2022. This is conservative, since large eruptions (is HTHH "large"?) decay away with an e-folding lifetime more like 16 months, see https://acp.copernicus.org/articles/23/921/2023/
  - calculate ERF contribution based on -20 x SAOD + WV (-20 from AR6, should be sampled, as should the other uncertainties)

A back of the envelope sense check for the 0.007 AOD comes from:
- total SO4 burden increase 0.66 Tg (https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL105076)
- AOD per TgSO4, in a SRM context (couldn't find a paper on volcanoes, quickly) is 0.004-0.02 AOD/Tg depending on latitude and injection rate (fig 2d, https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2023GL104417). Taking 0.01 as a mid-range guess and multiply 0.066 x 0.01 = 0.0066 AOD units (rounds to 0.007). Fig. 2c in this paper also gives another line of evidence for a 16 month lifetime.

experiment: IGCC (includes HTHH)