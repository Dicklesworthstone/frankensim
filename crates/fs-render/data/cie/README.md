# CIE 1931 observer reference data

The numerical data in this directory are provided by the **International
Commission on Illumination (CIE)** under the **Creative Commons
Attribution-ShareAlike 4.0 International** license:
<https://creativecommons.org/licenses/by-sa/4.0/>.

This license applies to the original CSV and its metadata and to the adapted
numeric table `xyz_1931_2deg_5nm_midpoints.in`. These data retain their own
license separately from FrankenSim's original implementation code.

## Attribution and original sources

CIE (2019), *Colour-matching functions of CIE 1931 standard colorimetric
observer*, International Commission on Illumination (CIE), Vienna, AT.
DOI: <https://doi.org/10.25039/CIE.DS.xvudnb9b>.

- Dataset description:
  <https://cie.co.at/datatable/cie-1931-colour-matching-functions-2-degree-observer>
- Original CSV:
  <https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv>
- Original metadata:
  <https://files.cie.co.at/Publications-datasets/CIE_xyz_1931_2deg.csv_metadata.json>

The CSV and metadata were retrieved on 2026-10-30 and are retained unchanged,
including their original CRLF line endings. The original CSV contains the
471 one-nanometre rows from 360 through 830 nm, taken from CIE 018:2019,
Table 6. Its SHA-256 matches the checksum published in the metadata:

```text
fa663e3535a7e0763a745993a1f0a192eb0275ac46ad2d1befd7626841e713c1
```

Metadata file SHA-256, recorded at retrieval:

```text
03abcaecf4e63d77045ef57c4514b52c8bb1a46dd18e1f93a50044a0f4f481c8
```

## Adaptation for the existing photometric representation

`xyz_1931_2deg_5nm_midpoints.in` contains `[x_bar, y_bar, z_bar]` at
`382.5 + 5*b` nm for `b = 0..79`. Each entry is the arithmetic mean of the
two adjacent original one-nanometre rows. This follows the metadata's linear
interpolation rule. The transformation selects and interpolates the data; it
does not fit a Gaussian, normalize XYZ to a display white, or expand the
renderer beyond its existing 380--780 nm band. The adapted table is distributed
under the same CC BY-SA 4.0 license with the attribution above. CIE does not
endorse this adaptation or FrankenSim.

The original data remain available to independent tests, which parse the CSV
and interpolate each of the 80 reference wavelengths separately from the
production lookup. The implementation loads a fixed numeric table without
runtime parsing or allocation. The legacy renderer's analytic observer is a
separate numerical convention and is not replaced by these data.

## Physical scope

The metadata describes these data as **nominal** and supplies no uncertainty
enclosure. Agreement with this table removes the analytic-observer fit error
at the represented wavelengths; it does not bound spectral discretization,
unresolved narrow emission or absorption lines, geometry, material data,
finite path depth, or fixture calibration. Wavelength-convention differences
between tabulated observer data and physical material laws are not corrected
here. Neither this dataset nor its inclusion establishes CIE 171 conformance
or agreement with measurements of real lighting fixtures.
