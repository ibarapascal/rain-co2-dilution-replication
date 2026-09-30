# Third-party data and code (attribution)

The derived data in this repository are released under CC BY 4.0 and the code under MIT, but these licenses do not replace the terms of
the third-party sources. Anyone reusing the derived tables should also cite the sources below. Terms were checked on 2026-09-30 from
the sources named in the "Terms checked from" column; "not checked" means no primary statement was found and the entry is an
inference. Identifiers are DOIs where the provider issues them; otherwise the landing page or archive named.

"Redistributed here" says whether values of the source appear in this repository: "no" = only used to compute derived quantities;
"derived values" = per-event or aggregated quantities computed from it; "extract" = source values themselves (small subsets).

## Observations

| Dataset (version) | Identifier | Provider | Terms (summary) | Terms checked from | Acknowledgment / citation requested | Redistributed here |
|---|---|---|---|---|---|---|
| MAPCO2 mooring pCO₂ and salinity, 13 moorings (ERDDAP `pmel_co2_moorings_*`; NCEI accessions 0100071-0162473) | NCEI accession DOIs of the 13 moorings (Sutton et al., 2012-2019), listed on the NCEI OCADS accession pages | NOAA PMEL Carbon Program | Data freely available; PMEL asks to be informed at the outset of work intended for publication, notes that co-authorship may be appropriate when the data are essential, and asks that manuscripts be sent to PMEL for review before submission | `license` attribute of the ERDDAP data sets (NOAA text) | Credit PMEL; cite the NCEI accessions | derived values (event ΔS at the MAPCO2 channel, pCO₂ changes) |
| GTMBA (TAO/TRITON, RAMA) rain, salinity, wind, OceanSITES delayed/mixed mode | OceanSITES GDAC (NDBC) | NOAA PMEL GTMBA Project Office / NDBC | File attribute: collected and made freely available by NDBC; no license attribute in the file checked | `citation` attribute of an OceanSITES TAO file (T0N165E) | GTMBA Project Office acknowledgment wording not checked | derived values |
| KEO and Papa (OceanSITES) | OceanSITES GDAC | NOAA PMEL Ocean Climate Stations | Freely available without restriction | `license` attribute (OS_KEO_*, OS_PAPA_* files) | Acknowledge the OCS Project Office of NOAA/PMEL; preprints requested | derived values |
| WHOTS (OceanSITES) | OceanSITES GDAC | WHOI / University of Hawaii | CLIVAR data policy: free of charge; user must display the citation in any publication; contact PI before commercial use | `license` attribute (OS_WHOTS_* files) | Citation text of the `citation` attribute (WHOTS Ocean Reference Station, Plueddemann and Weller, WHOI; NOAA funding) | derived values |
| Stratus (OceanSITES) | OceanSITES GDAC | WHOI | As WHOTS (CLIVAR policy) | `license` attribute (OS_Stratus_* file) | Citation text of the `citation` attribute (Stratus Ocean Reference Station, Weller, WHOI; NOAA OGP support) | no qualifying events (station metadata only) |
| SOFS (OceanSITES) | OceanSITES GDAC | IMOS (Australia) | CC BY 4.0 | `license` attribute (OS_SOTS_SOFS05 file) | Required statement (from the `acknowledgement` attribute): "Data was sourced from the Integrated Marine Observing System (IMOS) - IMOS is supported by the Australian Government through the National Collaborative Research Infrastructure Strategy and the Super Science Initiative." | derived values |
| SPURS-2 salinity snake, R/V Revelle meteorology, central mooring, Lady Amber (V1.0) | 10.5067/SPUR2-SNAKE, 10.5067/SPUR2-MET00, 10.5067/SPUR2-MOOR1, 10.5067/SPUR2-LAMBR | NASA PO.DAAC | NASA Earth science data: generally unrestricted (CC0 at the program level); cite the data; no implied NASA endorsement | NASA Earthdata "Data Use and Citation Guidance" page; DataCite records carry no rights field | Dataset citations | derived values (per-event ΔS) |
| SPURS-1 WHOI central mooring (V1.0) | 10.5067/SPUR1-MOOR1 | NASA PO.DAAC | as above | as above | Dataset citation | derived values |
| PISTON 2018 R/V Thompson, 2019 R/V Sally Ride | 10.5067/SUBORBITAL/PISTON2018-ONR-NOAA/RVTHOMPSON/DATA001, .../PISTON2019-ONR-NOAA/RVSALLYRIDE/DATA001 | NASA ASDC | as above | as above | Dataset citations | derived values; extract: hourly ship positions (`piston_track_hourly.csv`) |
| NCEI OCADS per-deployment metadata (XML), OceanSITES metadata (.das/.dds) | NCEI accessions above; OceanSITES GDAC | NOAA NCEI; OceanSITES | Public metadata; not checked | not checked | — | no (SHA-256 list only) |

## Satellite and gridded products

| Dataset (version) | Identifier | Provider | Terms (summary) | Terms checked from | Acknowledgment / citation requested | Redistributed here |
|---|---|---|---|---|---|---|
| CMORPH CDR, 8 km / 30 min, V1.0 (bias-corrected) | 10.25921/w9va-q159 | NOAA NCEI | No access or use restrictions stated; no warranty; cite the dataset | NCEI landing page (access/use constraints) and CDR product page | Dataset citation | extract: `data/products/p2/cmorph_pixels.csv` (rain rates at 13 mooring pixels, 14.5 MB); derived values |
| ERA5 (NSF NCAR mirror, 0.25°) | 10.5065/BH6N-5N20 | ECMWF / Copernicus C3S, via NSF NCAR | CC BY 4.0 | DataCite record of the DOI | Attribution to Copernicus C3S/ECMWF (exact statement not checked) | derived aggregates only (fluxes by band/basin/month) |
| HYCOM GOFS 3.1 reanalysis GLBv0.08 expt_53.X | hycom.org | HYCOM Consortium | Approved for public release, distribution unlimited; provided as is | hycom.org data server page | No statement found | derived aggregates only |
| OISST v2.1 | 10.25921/RE9P-PT57 | NOAA NCEI | not checked (no rights field in DataCite; NOAA data, inferred unrestricted) | — | Dataset citation | derived aggregates only |
| Watson et al. fCO₂ product (RECCAP2-ocean data collection) | 10.5281/zenodo.7990823 | Zenodo (Müller, 2023) | CC BY 4.0 | DataCite record | Cite Müller (2023) and Watson et al. (2020) | derived aggregates only |
| GLODAPv2.2016b mapped climatology | glodap.info (Lauvset et al., 2016) | GLODAP | CC BY 4.0 stated for the website contents; applies to the product files by inference | glodap.info mapped-data page | Cite Lauvset et al. (2016) | extract: surface TAlk, TCO₂, salinity, temperature at the 13 mooring cells (`p2_glodap.json`); derived aggregates |
| WOA09 basin mask (`basin.msk`) | NCEI WOA09 | NOAA NCEI | not checked (NOAA data, inferred unrestricted) | — | Cite NODC (2010) / Locarnini et al. (2010) | derived: basin codes in `p3_regions.csv` |
| SMAP RSS L2C SSS V6.0 | 10.5067/SMP60-2SOCS | Remote Sensing Systems, via NASA PO.DAAC | RSS makes a specific citation a condition of use and asks for an acknowledgment that the data are produced by RSS and sponsored by the NASA Ocean Salinity Science Team (www.remss.com); NASA program-level terms as above | RSS SMAP salinity web page | Citation: Meissner, Wentz, Manaster, Lindsley, Brewer, Densberger (2024), RSS SMAP Ocean Surface Salinities, Version 6.0 validated release; acknowledgment statement as on the RSS page | derived values (per-overpass satellite ΔSSS, IMERG rain supplied in the RSS files) |
| SMAP JPL L2B CAP SSS V5.0 | 10.5067/SMP50-2TOCS | NASA/JPL, via PO.DAAC | NASA program-level terms as above | NASA Earthdata guidance | Dataset citation | derived values |
| SOCATv2025 | 10.25921/648f-fv35 | NOAA NCEI | CC BY 4.0 | DataCite record | Cite Bakker et al. (2025) | no (reference for the constant c used in P3; not downloaded) |

## Code

| Code | Identifier | Terms | Checked from | Redistributed here |
|---|---|---|---|---|
| Witte, Zappa and McGillis (2026a), "Code For: On the Importance of Rain-Induced Dilution to the Ocean Carbon Sink", Code Ocean capsule 6552187 (files `code/CO2_Rain_Flux_Toolbox.py`, `code/main.py`) | 10.24433/CO.9378898.v1 | MIT License (DataCite rights list; the record also lists CC0, presumably for data). The license file and copyright line in the capsule were not retrieved (the capsule API refused the request on 2026-09-30) | DataCite record | no; `code/tools/fetch_witte_capsule.py` downloads the two files and checks their SHA-256. `code/pipeline/p3_global.py` re-implements the capsule's formulas (equations and the d₀ table are transcribed with source line references in the docstrings of `p3_global.py` and `p2_rim_test.py`) |

## Obligations when reusing the data

- PMEL MAPCO2: the terms ask users to inform PMEL at the outset of work intended for publication and to send manuscripts to PMEL for
  review before submission; credit PMEL and cite the NCEI accessions.
- RSS SMAP: cite Meissner et al. (2024) and include the RSS acknowledgment statement.
- IMOS (SOFS): the acknowledgment statement above is required, including for re-packaged data.
- WHOTS and Stratus: display the citation text given in the file attributes.
- KEO/Papa (PMEL OCS): acknowledge the OCS Project Office of NOAA/PMEL. GTMBA: acknowledgment wording not checked.
- ERA5: attribute Copernicus C3S/ECMWF (exact statement not checked).
- CMORPH, GLODAP, PISTON: `data/products/p2/cmorph_pixels.csv`, `data/products/p2/p2_glodap.json` and
  `data/products/piston-track/piston_track_hourly.csv` contain source values; cite those datasets when reusing these files.
