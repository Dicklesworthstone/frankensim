//! Concentric annular Poiseuille resistance, prepared for the original gas owner.
//!
//! Newtonian specialization of Huilgol and Georgiou (2020), equations (18),
//! (22), (23), <https://doi.org/10.1016/j.jnnfm.2020.104401>.

use super::{Length, Mesh1d, WallExchangeError, positive, same_geometry};
use fs_math::det;

/// Quasi-steady no-slip resistance between concentric circular cylinders.
///
/// Both surfaces share the original contact's axial velocity and temperature.
/// Their complete wetted area is retained. A distinct inner/outer temperature,
/// differential cylinder motion, eccentricity, rotating flow or wall slip needs
/// a different model; none is inferred from the hydraulic diameter.
///
/// The explicitly authored Nusselt number uses `Dh = 2 (outer - inner)` and
/// must apply to the declared radius ratio and thermal boundary conditions.
/// Current viscosity and density still come from the original gas source.
/// This developed-flow approximation does not resolve entrance development,
/// unsteady transverse profiles, rarefaction or turbulent annular flow.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AnnularLaminarWall {
    inner_radius: Length,
    outer_radius: Length,
    area_m2: f64,
    hydraulic_diameter: Length,
    viscous_length_m: f64,
    viscous_factor: f64,
    nusselt_number: f64,
    reynolds_number: [f64; 2],
    maximum_mach_number: f64,
}

impl AnnularLaminarWall {
    /// Prepare exact Newtonian annular geometry and declared applicability.
    ///
    /// Radii must satisfy finite `0 <= inner < outer`. Zero inner radius is
    /// the circular-pipe limit. Nu is finite and nonnegative; zero disables
    /// radial heat exchange. Reynolds support lies within `[0, 2300]`, with
    /// positive upper limit, and the wall-relative Mach ceiling is in `(0,0.3]`.
    /// These are screening bounds, not proof of developed laminar flow.
    pub fn try_new(
        inner_radius: Length,
        outer_radius: Length,
        nusselt_number: f64,
        reynolds_number: [f64; 2],
        maximum_mach_number: f64,
    ) -> Result<Self, WallExchangeError> {
        let (inner, outer) = (inner_radius.value(), outer_radius.value());
        let [minimum, maximum] = reynolds_number;
        if !inner.is_finite()
            || inner < 0.0
            || !positive(outer)
            || inner >= outer
            || !nusselt_number.is_finite()
            || nusselt_number < 0.0
            || !minimum.is_finite()
            || minimum < 0.0
            || !positive(maximum)
            || minimum > maximum
            || maximum > 2300.0
            || !positive(maximum_mach_number)
            || maximum_mach_number > 0.3
        {
            return Err(WallExchangeError::Invalid(
                "annular wall geometry or applicability",
            ));
        }
        let gap = outer - inner;
        let diameter = 2.0 * gap;
        let area_m2 = if inner == 0.0 {
            // Match the existing circular mesh and rate arithmetic exactly.
            core::f64::consts::FRAC_PI_4 * diameter * diameter
        } else {
            (core::f64::consts::PI * gap) * (outer + inner)
        };
        if !positive(area_m2) || !positive(diameter) {
            return Err(WallExchangeError::Invalid(
                "unrepresentable annular wall geometry",
            ));
        }
        let ratio = inner / outer;
        let (viscous_length_m, viscous_factor) = if inner == 0.0 {
            (diameter, 32.0)
        } else if ratio >= 0.5 {
            // Mean velocity is gap² C dp/(12 mu L). A positive series retains
            // C near one without cancellation of almost equal fourth powers.
            let mean = 0.5 * outer + 0.5 * inner;
            (gap, 12.0 / gap_correction(0.5 * (gap / mean)))
        } else {
            let q2 = ratio * ratio;
            let logarithm = if ratio > 0.0 {
                -det::ln(ratio)
            } else {
                // A tiny nonzero core retains logarithmic drag even if r/R
                // underflows. It must not silently become an empty pipe.
                det::ln(outer) - det::ln(inner)
            };
            let difference = 1.0 - q2;
            let factor = 1.0 - q2 * q2 - difference * difference / logarithm;
            (outer, 8.0 * difference / factor)
        };
        if !positive(viscous_factor) {
            return Err(WallExchangeError::Invalid(
                "unrepresentable annular shape factor",
            ));
        }
        Ok(Self {
            inner_radius,
            outer_radius,
            area_m2,
            hydraulic_diameter: Length::new(diameter),
            viscous_length_m,
            viscous_factor,
            nusselt_number,
            reynolds_number,
            maximum_mach_number,
        })
    }

    /// Original inner radius [m]; exactly zero selects an empty circular pipe.
    #[must_use]
    pub const fn inner_radius(self) -> Length {
        self.inner_radius
    }

    /// Original outer radius [m].
    #[must_use]
    pub const fn outer_radius(self) -> Length {
        self.outer_radius
    }

    /// Actual cross-sectional area [m²], evaluated without subtracting squares.
    #[must_use]
    pub const fn flow_area_m2(self) -> f64 {
        self.area_m2
    }

    /// Twice the physical gap [m], used by Reynolds, Nu and wetted area.
    #[must_use]
    pub const fn hydraulic_diameter(self) -> Length {
        self.hydraulic_diameter
    }

    /// Authored `h Dh/k` for this annulus and its common-temperature surfaces.
    #[must_use]
    pub const fn nusselt_number(self) -> f64 {
        self.nusselt_number
    }

    /// Inclusive Reynolds support based on hydraulic diameter and relative speed.
    #[must_use]
    pub const fn reynolds_number(self) -> [f64; 2] {
        self.reynolds_number
    }

    /// Original wall-relative Mach ceiling.
    #[must_use]
    pub const fn maximum_mach_number(self) -> f64 {
        self.maximum_mach_number
    }

    pub(super) fn admit_mesh(self, mesh: &Mesh1d) -> Result<(), WallExchangeError> {
        let Some(area) = mesh.uniform_cross_section() else {
            return Err(WallExchangeError::Invalid(
                "developed annular wall requires constant area",
            ));
        };
        if !same_geometry(self.area_m2, area) {
            return Err(WallExchangeError::Invalid(
                "annular radii differ from original branch area",
            ));
        }
        for (faces, &volume) in mesh.positions().windows(2).zip(mesh.volumes()) {
            if !same_geometry(volume, area * (faces[1] - faces[0])) {
                return Err(WallExchangeError::Invalid(
                    "annular wall volume differs from original geometry",
                ));
            }
        }
        Ok(())
    }

    // U = Q/[pi(R²-r²)] and Q = pi*dp*D/(8*mu*L), where
    // D = R^4-r^4-(R²-r²)²/ln(R/r). No series or logarithm in the hot path.
    pub(super) fn linear_rate(
        self,
        viscosity: f64,
        density: f64,
    ) -> Result<f64, WallExchangeError> {
        let rate = self.viscous_factor * (viscosity / density)
            / self.viscous_length_m
            / self.viscous_length_m;
        if !positive(rate) {
            return Err(WallExchangeError::Invalid(
                "unrepresentable annular laminar wall rate",
            ));
        }
        Ok(rate)
    }
}

fn gap_correction(half_gap_over_mean: f64) -> f64 {
    // e=(R-r)/(R+r)<=1/3, A=atanh(e)/e. The exact normalized conductance
    // C = 3/4 * sum[(1/(2k+1)+1/(2k+3))*e^(2k)] / sum[e^(2k)/(2k+1)]
    // is positive and approaches one. The omitted denominator tail after
    // k=16 is <2e-18; Horner evaluation needs fixed work only at preparation.
    let z = half_gap_over_mean * half_gap_over_mean;
    let (mut numerator, mut denominator) = (0.0, 0.0);
    for k in (0..=16).rev() {
        let n = f64::from(2 * k + 1);
        denominator = 1.0 / n + z * denominator;
        numerator = (1.0 / n + 1.0 / (n + 2.0)) + z * numerator;
    }
    0.75 * numerator / denominator
}

#[cfg(test)]
mod tests {
    use super::*;

    fn wall(inner: f64, outer: f64) -> AnnularLaminarWall {
        AnnularLaminarWall::try_new(
            Length::new(inner),
            Length::new(outer),
            0.0,
            [0.0, 2300.0],
            0.3,
        )
        .unwrap()
    }

    fn close(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 3e-13 * expected.abs(),
            "{actual:.17e} != {expected:.17e}"
        );
    }

    #[test]
    fn g1_annular_rates_match_high_precision_integrated_velocity() {
        // Independent 100-decimal-digit evaluation of 8*(1-r²)/
        // [1-r^4-(1-r²)²/ln(1/r)] for the exact binary64 r, outer radius 1.
        for (inner, expected) in [
            (0.001, 9.354_139_711_611_522),
            (0.1, 13.791_951_037_652_433),
            (0.25, 20.712_687_345_948_765),
            (0.5, 47.625_080_318_225_52),
            (0.9, 1199.778_112_258_746_7),
            (0.99, 119999.797_982_567_12),
            (0.99999, 120000000000.892_24),
            (0.99999999, 1.199_999_987_940_577_8e17),
            (0.999999999999, 1.200_053_093_890_077_3e25),
        ] {
            close(wall(inner, 1.0).linear_rate(1.0, 1.0).unwrap(), expected);
        }
    }

    #[test]
    fn g3_zero_core_matches_original_circular_area_and_rate_bits() {
        for radius in [1e-3, 1.0, 1e3] {
            let model = wall(0.0, radius);
            let diameter = 2.0 * radius;
            let rate = 32.0 * (0.004 / 1.2) / diameter / diameter;
            assert_eq!(
                model.linear_rate(0.004, 1.2).unwrap().to_bits(),
                rate.to_bits()
            );
            assert_eq!(
                model.flow_area_m2(),
                core::f64::consts::FRAC_PI_4 * diameter * diameter
            );
            assert_eq!(model.hydraulic_diameter().value(), diameter);
        }
    }

    #[test]
    fn g1_narrow_gap_retains_parallel_plate_limit_without_subtraction_loss() {
        for inner in [0.99, 0.9999, 0.99999999, 0.999999999999] {
            let model = wall(inner, 1.0);
            let gap = 1.0 - inner;
            let normalized = model.linear_rate(1.0, 1.0).unwrap() * gap * gap;
            assert!((normalized / 12.0 - 1.0).abs() < 2e-6);
            assert!(model.flow_area_m2() > 0.0);
        }
    }

    #[test]
    fn g3_vanishing_core_keeps_logarithmic_drag_even_when_radius_ratio_underflows() {
        for inner in [1e-16_f64, 1e-64, 1e-200] {
            let rate = wall(inner, 1.0).linear_rate(1.0, 1.0).unwrap();
            close(rate, 8.0 / (1.0 + 1.0 / inner.ln()));
            assert!(rate > 8.0);
        }
        let model = wall(1e-250, 1e90);
        close(
            model.linear_rate(1.0, 1.0).unwrap(),
            8e-180 / (1.0 - 1.0 / (340.0 * 10.0_f64.ln())),
        );
    }

    #[test]
    fn g3_material_and_length_scaling_preserve_prepared_resistance() {
        let reference = wall(0.25, 1.0).linear_rate(0.004, 1.2).unwrap();
        for scale in [1e-3, 3.0, 1e3] {
            let model = wall(0.25 * scale, scale);
            close(
                model.linear_rate(0.004, 1.2).unwrap() * scale * scale,
                reference,
            );
            close(
                model.linear_rate(0.012, 2.4).unwrap() * scale * scale,
                1.5 * reference,
            );
        }
    }

    #[test]
    fn g0_invalid_annular_geometry_and_domains_refuse() {
        for (inner, outer) in [
            (-1.0, 1.0),
            (1.0, 1.0),
            (2.0, 1.0),
            (0.0, 0.0),
            (f64::NAN, 1.0),
            (0.0, f64::INFINITY),
            (0.0, 1e-200),
        ] {
            assert!(
                AnnularLaminarWall::try_new(
                    Length::new(inner),
                    Length::new(outer),
                    0.0,
                    [0.0, 2300.0],
                    0.3
                )
                .is_err()
            );
        }
        assert!(
            AnnularLaminarWall::try_new(
                Length::new(0.5),
                Length::new(1.0),
                0.0,
                [0.0, 2301.0],
                0.3
            )
            .is_err()
        );
        assert!(
            AnnularLaminarWall::try_new(
                Length::new(0.5),
                Length::new(1.0),
                -1.0,
                [0.0, 2300.0],
                0.3
            )
            .is_err()
        );
    }
}
