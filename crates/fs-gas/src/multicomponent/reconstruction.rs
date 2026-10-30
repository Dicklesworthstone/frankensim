//! Compatible mass-fraction slopes for conservative multicomponent flow.
//!
//! Density and composition are independent reconstruction coordinates. Limiting
//! partial densities separately can destroy a uniform transported mass fraction
//! when density changes sharply. Limiting each fraction separately can also
//! destroy linear composition constraints, including elemental fractions. One
//! common limiter scales the complete centered composition slope instead.

use super::{Error, MAX_SPECIES, positive, sum};

/// Reconstruct all canonical mass-fraction slopes on a uniform three-cell
/// stencil. Inputs are species mass densities, including exact zero species.
///
/// Each unlimited slope is `(Y_right-Y_left)/2`. The smallest admissible MC
/// factor over all species scales the entire slope vector. Consequently a
/// linear composition constraint shared by the stencil, including `sum(Y)=1`
/// and a uniform inert or elemental fraction, is retained to floating-point
/// roundoff. No species state, fraction or conserved inventory is clipped,
/// seeded, renormalized or replaced by a remainder species.
///
/// Centered differences no larger than
/// `32*EPSILON*max(abs(Y_left),abs(Y_center),abs(Y_right))`
/// do not constrain the common factor. This relative-only roundoff guard keeps
/// a constant fraction's division noise from flattening every other species;
/// it does not discard small but resolved trace-species gradients. Their raw
/// centered slopes are still scaled by the same factor. Output is scratch until
/// success. All slices must have the same nonzero length, at most `MAX_SPECIES`.
pub fn mass_fraction_slopes(
    left: &[f64],
    center: &[f64],
    right: &[f64],
    out: &mut [f64],
) -> Result<(), Error> {
    slopes_with_geometry(left, center, right, [1.0, 0.5, 1.0], out)
}

// Geometry factors multiply left, centered and right fraction differences.
// For cell width h and neighboring centre distances dl/dr, these factors are
// h/dl, h/(dl+dr), h/dr. Uniform cells give [1, 1/2, 1].
pub(super) fn slopes_with_geometry(
    left: &[f64],
    center: &[f64],
    right: &[f64],
    factors: [f64; 3],
    out: &mut [f64],
) -> Result<(), Error> {
    let n = center.len();
    if n == 0
        || n > MAX_SPECIES
        || left.len() != n
        || right.len() != n
        || out.len() != n
        || !factors.into_iter().all(positive)
    {
        return Err(Error::Invalid(
            "mass-fraction reconstruction dimensions or geometry",
        ));
    }
    let density = |values: &[f64]| {
        if !values.iter().all(|v| v.is_finite() && *v >= 0.0) {
            return Err(Error::Invalid(
                "mass-fraction reconstruction species density",
            ));
        }
        let rho = sum(values.iter().copied());
        if !positive(rho) {
            return Err(Error::Invalid(
                "mass-fraction reconstruction mixture density",
            ));
        }
        Ok(rho)
    };
    let rho = [density(left)?, density(center)?, density(right)?];
    let mut theta: f64 = 1.0;
    for k in 0..n {
        let yl = left[k] / rho[0];
        let yc = center[k] / rho[1];
        let yr = right[k] / rho[2];
        let dc = yr - yl;
        let slope = factors[1] * dc;
        out[k] = slope;
        let noise = 32.0 * f64::EPSILON * yl.abs().max(yc.abs()).max(yr.abs());
        if dc.abs() <= noise {
            continue;
        }
        let dl = yc - yl;
        let dr = yr - yc;
        // Compare signs without products: a resolved 1e-200 trace gradient
        // must not become an artificial extremum through product underflow.
        if dl == 0.0
            || dr == 0.0
            || dl.is_sign_positive() != dc.is_sign_positive()
            || dr.is_sign_positive() != dc.is_sign_positive()
        {
            theta = 0.0;
        } else {
            theta = theta
                .min(2.0 * factors[0] * (dl / slope))
                .min(2.0 * factors[2] * (dr / slope));
        }
    }
    for value in out {
        *value *= theta;
        if !value.is_finite() {
            return Err(Error::Invalid(
                "mass-fraction reconstruction slope overflow",
            ));
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn common_limiter_preserves_inert_and_linear_composition_constraints() {
        let fractions = [
            [0.1, 0.425, 0.05, 0.125, 0.3],
            [0.2, 0.3, 0.1, 0.1, 0.3],
            [0.28, 0.2525, 0.105, 0.0625, 0.3],
        ];
        let states: [[f64; 5]; 3] =
            std::array::from_fn(|i| fractions[i].map(|y| y * [1.1, 0.3, 5.7][i]));
        let mut out = [0.0; 5];
        mass_fraction_slopes(&states[0], &states[1], &states[2], &mut out).unwrap();
        assert!(out[0].abs() > 0.01);
        assert!(out[4].abs() < 3e-16);
        assert!(sum(out).abs() < 3e-16);
        assert!((out[0] + 2.0 * out[1] + 3.0 * out[2]).abs() < 5e-16);
        for k in 0..5 {
            let lower = fractions.iter().map(|q| q[k]).fold(f64::INFINITY, f64::min);
            let upper = fractions
                .iter()
                .map(|q| q[k])
                .fold(f64::NEG_INFINITY, f64::max);
            for sign in [-0.5, 0.5] {
                let trace = fractions[1][k] + sign * out[k];
                assert!(trace >= lower - 3e-16 && trace <= upper + 3e-16);
            }
        }
        let mut nonuniform = [0.0; 5];
        slopes_with_geometry(
            &states[0],
            &states[1],
            &states[2],
            [0.8, 0.6, 1.3],
            &mut nonuniform,
        )
        .unwrap();
        assert!(sum(nonuniform).abs() < 3e-16);
        assert!((nonuniform[0] + 2.0 * nonuniform[1] + 3.0 * nonuniform[2]).abs() < 5e-16);
    }

    #[test]
    fn relative_roundoff_guard_retains_resolved_trace_species_and_structural_zeros() {
        let mut out = [0.0; 3];
        mass_fraction_slopes(
            &[1.0, 1e-200, 0.0],
            &[1.0, 2e-200, 0.0],
            &[1.0, 4e-200, 0.0],
            &mut out,
        )
        .unwrap();
        assert_eq!(out[0], 0.0);
        assert!(out[1] > 1e-200 && out[1] < 2e-200);
        assert_eq!(out[2], 0.0);
        assert!(mass_fraction_slopes(&[0.0], &[0.0], &[0.0], &mut [0.0]).is_err());
        assert!(mass_fraction_slopes(&[1.0], &[1.0], &[1.0], &mut []).is_err());
    }
}
