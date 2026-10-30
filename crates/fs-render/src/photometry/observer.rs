//! Tabulated CIE 1931 observer for calibrated photometry.
//!
//! The three matching functions come from CIE's published 1 nm table,
//! linearly interpolated at the existing 80 spectral-bin midpoints. This
//! removes the analytic observer fit from calibrated source normalization and
//! XYZ measurements without changing the renderer's represented 380--780 nm
//! spectra. It does not resolve sub-bin spectral structure, bound spectral
//! quadrature error, or establish fixture accuracy or CIE 171 conformance.
//!
//! The original CSV, its metadata and the separately licensed derived table
//! are retained in `data/cie`. The data are nominal reference values, not
//! uncertainty enclosures. No air/vacuum wavelength conversion is inferred.

use crate::spectral::{LAMBDA_MAX, LAMBDA_MIN, QUAD_BINS};

/// Published observer and the exact sampling convention used here.
///
/// Data citation: CIE (2019), *Colour-matching functions of CIE 1931 standard
/// colorimetric observer*, <https://doi.org/10.25039/CIE.DS.xvudnb9b>.
pub const OBSERVER_ID: &str = "CIE.DS.xvudnb9b:1931-2deg:linear-1nm:380-780nm:5nm-midpoints:v1";

// Numeric data are an attributed CC BY-SA 4.0 adaptation, isolated from this
// implementation. See data/cie/README.md for the original reference and hash.
const XYZ_MIDPOINTS: [[f64; 3]; QUAD_BINS] =
    include!("../../data/cie/xyz_1931_2deg_5nm_midpoints.in");

/// Dimensionless CIE `[x_bar, y_bar, z_bar]` at spectral bin `bin`.
///
/// Bin `b` represents the midpoint `382.5 + 5*b` nm. Multiply these functions
/// by spectral radiance, bin width in nm and 683 lm/W to form absolute XYZ,
/// with luminance in cd/m² in Y. There is no display-white normalization.
///
/// # Panics
///
/// Panics if `bin` is outside the existing 80-bin spectral domain.
#[must_use]
pub fn xyz_at_bin(bin: usize) -> [f64; 3] {
    XYZ_MIDPOINTS[bin]
}

/// Observer for the calibrated tracer's exact midpoint wavelength packets.
///
/// Like the represented source spectrum, this lookup is constant inside each
/// 5 nm bin; it is not an arbitrary-wavelength interpolation API. The upper
/// endpoint belongs to the final bin. Nonfinite or out-of-band wavelengths
/// have zero response. Physical transport supplies the bin midpoints.
#[must_use]
pub(crate) fn xyz_at_wavelength(lambda_nm: f64) -> [f64; 3] {
    if !lambda_nm.is_finite() || !(LAMBDA_MIN..=LAMBDA_MAX).contains(&lambda_nm) {
        return [0.0; 3];
    }
    let width_nm = (LAMBDA_MAX - LAMBDA_MIN) / QUAD_BINS as f64;
    let bin = ((lambda_nm - LAMBDA_MIN) / width_nm) as usize;
    xyz_at_bin(bin.min(QUAD_BINS - 1))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn original_rows() -> Vec<[f64; 4]> {
        include_str!("../../data/cie/CIE_xyz_1931_2deg.csv")
            .lines()
            .map(|line| {
                let values: Vec<f64> = line
                    .split(',')
                    .map(|value| value.parse().expect("published decimal data"))
                    .collect();
                values
                    .try_into()
                    .expect("wavelength and three CIE functions")
            })
            .collect()
    }

    #[test]
    fn g1_all_midpoints_match_independent_interpolation_of_original_cie_csv() {
        let rows = original_rows();
        assert_eq!(rows.len(), 471);
        assert_eq!(rows.first().unwrap()[0], 360.0);
        assert_eq!(rows.last().unwrap()[0], 830.0);
        for (index, row) in rows.iter().enumerate() {
            assert_eq!(row[0], 360.0 + index as f64);
        }
        for bin in 0..80 {
            // Deliberately find the bracketing original rows independently of
            // the production lookup and source wavelength helper.
            let wavelength = 382.5 + 5.0 * bin as f64;
            let upper = rows.partition_point(|row| row[0] < wavelength);
            let lower = upper - 1;
            let t = (wavelength - rows[lower][0]) / (rows[upper][0] - rows[lower][0]);
            let actual = xyz_at_bin(bin);
            for channel in 0..3 {
                let expected = (1.0 - t) * rows[lower][channel + 1] + t * rows[upper][channel + 1];
                let residual = (actual[channel] - expected).abs();
                // Decimal parsing/interpolation roundoff only. This is not a
                // physical uncertainty bound for the underlying CIE data.
                assert!(
                    residual <= 4.0 * f64::EPSILON * expected.abs(),
                    "bin={bin}, channel={channel}, actual={}, reference={expected}",
                    actual[channel]
                );
            }
        }
    }

    #[test]
    fn g1_published_blue_green_red_anchors_keep_absolute_observer_scale() {
        // Adjacent 1 nm rows averaged at 447.5, 552.5 and 652.5 nm.
        // The red z_bar is exactly zero in the published table.
        for (bin, expected) in [
            (13, [0.343_035_6, 0.033_702_4, 1.780_449]),
            (34, [0.472_2, 0.998_605_15, 0.007_083_5]),
            (54, [0.249_654, 0.093_655_8, 0.0]),
        ] {
            for (actual, expected) in xyz_at_bin(bin).into_iter().zip(expected) {
                assert!((actual - expected).abs() <= 4.0 * f64::EPSILON * expected.abs());
            }
        }
        // This unnormalized sum is in nanometres. A display-normalized
        // observer would incorrectly turn it into one.
        let integral: f64 = (0..80).map(|bin| 5.0 * xyz_at_bin(bin)[1]).sum();
        assert!((integral - 106.856_341_122_3).abs() < 2.0e-12);
    }

    #[test]
    fn g0_midpoint_packets_and_visible_band_boundaries_use_the_same_observer() {
        for bin in 0..80 {
            assert_eq!(xyz_at_wavelength(382.5 + 5.0 * bin as f64), xyz_at_bin(bin));
        }
        assert_eq!(xyz_at_wavelength(380.0), xyz_at_bin(0));
        assert_eq!(xyz_at_wavelength(780.0), xyz_at_bin(79));
        for wavelength in [379.9, 780.1, f64::NAN, f64::NEG_INFINITY, f64::INFINITY] {
            assert_eq!(xyz_at_wavelength(wavelength), [0.0; 3]);
        }
    }
}
