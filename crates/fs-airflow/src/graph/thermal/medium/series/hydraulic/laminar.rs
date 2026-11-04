//! Fully developed no-slip duct resistance from the actual cross section.
//!
//! Rectangular flow follows the Fourier solution of cross-sectional Poisson
//! flow, not circular Poiseuille with a substituted hydraulic diameter.
//! Reference: Mortensen, Okkels and Bruus, Phys. Rev. E 71, 057301 (2005),
//! equations (3), (10), (11), <https://doi.org/10.1103/PhysRevE.71.057301>.

use core::f64::consts::PI;

use fs_math::det;
use fs_qty::{Area, Length};

use super::{MaterialFlowResistance, ResistanceError, positive};

impl MaterialFlowResistance {
    /// Derive Hagen-Poiseuille resistance for a straight circular duct.
    ///
    /// `dp = 128 mu length Q / (pi diameter^4)`. The diameter also owns flow
    /// area and Reynolds length, so pressure loss and convection handoffs use
    /// consistent geometry. Dynamic viscosity remains an actual material input
    /// on every evaluation; it is never frozen into the geometry coefficient.
    ///
    /// This assumes a Newtonian fluid, stationary no-slip walls and a fully
    /// developed profile. The caller selects a Reynolds ceiling in `(0,2300]`;
    /// that ceiling is a screening limit, not proof of laminarity or developed
    /// flow. Entrance losses, slip, rarefaction and unsteady-profile effects
    /// require their own models. Gas use retains the original frozen-density
    /// approximation and its caller-owned pressure/Mach admission.
    ///
    /// # Errors
    /// Nonpositive/nonfinite geometry, unsupported ceiling, or unrepresentable
    /// area or resistance. Signed and zero flow remain admitted at evaluation.
    pub fn circular_laminar(
        length: Length,
        diameter: Length,
        maximum_reynolds: f64,
    ) -> Result<Self, ResistanceError> {
        admit(&[length.value(), diameter.value()], maximum_reynolds)?;
        let d = diameter.value();
        let area = Area::new((PI / 4.0) * d * d);
        let viscous = (128.0 / PI) * (length.value() / d) / d / d / d;
        Self::try_new(viscous, 0.0, area, diameter, (0.0, maximum_reynolds))
    }

    /// Derive fully developed resistance between concentric stationary cylinders.
    ///
    /// For inner/outer radii r/R, the no-slip solution integrates to
    /// `Q = pi dp [R^4-r^4-(R^2-r^2)^2/ln(R/r)] / (8 mu length)`.
    /// This is the Newtonian limit of Huilgol and Georgiou (2020), equations
    /// (18), (22), (23), <https://doi.org/10.1016/j.jnnfm.2020.104401>.
    ///
    /// The actual annular area supplies speed; `Dh = 2 (R-r)` supplies Reynolds
    /// number. A narrow gap uses a convergent positive series instead of
    /// subtracting nearly equal fourth-power terms. Zero inner radius selects
    /// the circular-pipe limit exactly; a nonzero core retains its logarithmic
    /// drag even when r/R underflows. Geometry is compiled once.
    ///
    /// The same Newtonian, developed-flow and applicability assumptions as
    /// [`Self::circular_laminar`] apply. Eccentricity, moving cylinders, slip
    /// and entrance effects are not modeled. An annular heat-transfer law
    /// still needs its own wall conditions; circular Nusselt numbers do not
    /// become applicable merely by substituting this hydraulic diameter.
    ///
    /// # Errors
    /// Requires finite `0 <= r < R`, positive length and a Reynolds ceiling
    /// in `(0,2300]`. Unrepresentable area, diameter or resistance refuses.
    pub fn annular_laminar(
        length: Length,
        inner_radius: Length,
        outer_radius: Length,
        maximum_reynolds: f64,
    ) -> Result<Self, ResistanceError> {
        admit(&[length.value(), outer_radius.value()], maximum_reynolds)?;
        let inner = inner_radius.value();
        let outer = outer_radius.value();
        if !inner.is_finite() || inner < 0.0 || inner >= outer {
            return Err(ResistanceError::Invalid(
                "annular radii must satisfy finite 0 <= inner < outer",
            ));
        }
        if inner == 0.0 {
            return Self::circular_laminar(length, Length::new(2.0 * outer), maximum_reynolds);
        }
        let gap = outer - inner;
        let area = Area::new((PI * gap) * (outer + inner));
        let diameter = Length::new(2.0 * gap);
        let ratio = inner / outer;
        let viscous = if ratio >= 0.5 {
            let mean = 0.5 * outer + 0.5 * inner;
            let correction = annular_gap_correction(0.5 * (gap / mean));
            (6.0 / (PI * correction)) * (length.value() / mean) / gap / gap / gap
        } else {
            let ratio_squared = ratio * ratio;
            let logarithm = if ratio > 0.0 {
                -det::ln(ratio)
            } else {
                det::ln(outer) - det::ln(inner)
            };
            let difference = 1.0 - ratio_squared;
            let factor =
                1.0 - ratio_squared * ratio_squared - difference * difference / logarithm;
            (8.0 / (PI * factor)) * (length.value() / outer) / outer / outer / outer
        };
        Self::try_new(viscous, 0.0, area, diameter, (0.0, maximum_reynolds))
    }

    /// Derive fully developed resistance for a straight rectangular duct.
    ///
    /// With smaller side `h`, larger side `w`, and `a=h/w`,
    /// `Q = w h^3 dp C(a) / (12 mu length)`, where
    /// `C(a) = 1 - 192 a/pi^5 sum_odd(tanh(n pi/(2a))/n^5)`.
    /// The entire side-wall correction is retained, including square ducts;
    /// the parallel-plate approximation is only the `a -> 0` limit.
    ///
    /// Flow area is `w h` and Reynolds length is the hydraulic diameter
    /// `2 w h/(w+h)`. Swapping width and height leaves the law identical.
    /// The same material, wall, Reynolds and flow-profile assumptions as
    /// [`Self::circular_laminar`] apply. Geometry is compiled once with fixed
    /// work; repeated flow/thermal trials use the original constant-cost law.
    pub fn rectangular_laminar(
        length: Length,
        width: Length,
        height: Length,
        maximum_reynolds: f64,
    ) -> Result<Self, ResistanceError> {
        admit(
            &[length.value(), width.value(), height.value()],
            maximum_reynolds,
        )?;
        let h = width.value().min(height.value());
        let w = width.value().max(height.value());
        let aspect = h / w;
        let correction = rectangular_conductance_factor(aspect);
        let area = Area::new(h * w);
        let diameter = Length::new(h * (2.0 / (1.0 + aspect)));
        let viscous = (12.0 / correction) * (length.value() / w) / h / h / h;
        Self::try_new(viscous, 0.0, area, diameter, (0.0, maximum_reynolds))
    }
}

// With e=(R-r)/(R+r), the exact conductance is pi*mean*gap^3/(6*mu*L)
// times C(e). Write A=atanh(e)/e and avoid subtracting 1/A from 1+e^2:
// C = 3/4 * sum[(1/(2k+1)+1/(2k+3))*e^(2k)] / sum[e^(2k)/(2k+1)].
// Here e <= 1/3; the denominator's omitted tail after k=16 is <2e-18.
// No logarithm, small-gap subtraction or per-evaluation series is needed.
fn annular_gap_correction(half_gap_over_mean: f64) -> f64 {
    let z = half_gap_over_mean * half_gap_over_mean;
    let mut denominator = 0.0;
    let mut numerator = 0.0;
    for k in (0..=16).rev() {
        let n = f64::from(2 * k + 1);
        denominator = 1.0 / n + z * denominator;
        numerator = (1.0 / n + 1.0 / (n + 2.0)) + z * numerator;
    }
    0.75 * numerator / denominator
}

fn admit(lengths: &[f64], maximum_reynolds: f64) -> Result<(), ResistanceError> {
    if !lengths.iter().copied().all(positive)
        || !positive(maximum_reynolds)
        || maximum_reynolds > 2300.0
    {
        return Err(ResistanceError::Invalid(
            "positive finite duct dimensions and Reynolds ceiling in (0,2300] required",
        ));
    }
    Ok(())
}

fn rectangular_conductance_factor(aspect: f64) -> f64 {
    // Sum odd n^-5 = (31/32) zeta(5). Subtract only the exponentially small
    // tanh corrections instead of truncating the algebraically decaying sum.
    // For 0 < aspect <= 1 the omitted correction starting at n=13 is <1e-23,
    // well below f64 roundoff. exp(-2x) avoids cancellation in 1-tanh(x) and
    // overflow for very slender ducts. aspect=0 is the rounded plate limit.
    const ODD_ZETA_FIVE: f64 = 1.004_523_762_795_139_6;
    let mut difference = 0.0;
    for n in [1.0_f64, 3.0, 5.0, 7.0, 9.0, 11.0] {
        let e = det::exp(-PI * n / aspect);
        difference += (2.0 * e / (1.0 + e)) / n.powi(5);
    }
    1.0 - (192.0 / PI.powi(5)) * aspect * (ODD_ZETA_FIVE - difference)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fluid_transport::FluidTransportState;
    use fs_qty::MassFlowRate;

    fn fluid(rho: f64, mu: f64) -> FluidTransportState {
        FluidTransportState::try_new(rho, 4000.0, 0.5, mu).unwrap()
    }

    fn near(actual: f64, expected: f64) {
        assert!(
            (actual - expected).abs() <= 2e-12 * expected.abs().max(1e-30),
            "{actual:e} != {expected:e}"
        );
    }

    #[test]
    fn circular_poiseuille_uses_dynamic_viscosity_signed_flow_and_actual_geometry() {
        let law =
            MaterialFlowResistance::circular_laminar(Length::new(0.8), Length::new(0.004), 1500.0)
                .unwrap();
        near(law.flow_area().value(), PI * 0.002_f64.powi(2));
        assert_eq!(law.reynolds_length(), Length::new(0.004));
        for (rho, mu) in [(1000.0, 0.001), (800.0, 0.004)] {
            for q in [-1e-7, 0.0, 1e-7] {
                let state = law
                    .evaluate(&fluid(rho, mu), MassFlowRate::new(rho * q))
                    .unwrap();
                let pressure = 8.0 * mu * 0.8 * q / (PI * 0.002_f64.powi(4));
                near(state.pressure_drop_pa, pressure);
                near(
                    state.dpressure_dmass_flow,
                    8.0 * mu * 0.8 / (rho * PI * 0.002_f64.powi(4)),
                );
                near(state.dissipated_power_w, pressure * q);
                near(state.reynolds, 4.0 * rho * q.abs() / (PI * mu * 0.004));
            }
        }
    }

    #[test]
    fn rectangular_series_retains_square_sidewalls_and_parallel_plate_limit() {
        for aspect in [1.0, 0.8, 0.5, 0.1, 0.001] {
            // Independent, unaccelerated reference series. Its omitted odd
            // tail is below 1/(8*8191^4), not a fit to the implementation.
            let series: f64 = (0..4096)
                .map(|k| {
                    let n = f64::from(2 * k + 1);
                    (n * PI / (2.0 * aspect)).tanh() / n.powi(5)
                })
                .sum();
            let expected = 1.0 - 192.0 * aspect * series / PI.powi(5);
            near(rectangular_conductance_factor(aspect), expected);
        }
        // Independent square-duct Darcy Poiseuille number, not circular 64.
        let law = MaterialFlowResistance::rectangular_laminar(
            Length::new(1.0),
            Length::new(0.01),
            Length::new(0.01),
            2300.0,
        )
        .unwrap();
        let state = law
            .evaluate(&fluid(1000.0, 0.001), MassFlowRate::new(0.01))
            .unwrap();
        let darcy_re =
            state.pressure_drop_pa * 2.0 * 0.01 / (1000.0 * 0.1_f64.powi(2)) * state.reynolds;
        near(darcy_re, 56.908_307_539_124_56);
        near(rectangular_conductance_factor(1e-14), 1.0);
    }

    #[test]
    fn rectangular_rotation_and_si_rescaling_preserve_the_physical_law() {
        let make = |scale, w, h| {
            MaterialFlowResistance::rectangular_laminar(
                Length::new(0.7 * scale),
                Length::new(w * scale),
                Length::new(h * scale),
                1000.0,
            )
            .unwrap()
        };
        let base = make(1.0, 0.02, 0.004);
        assert_eq!(base, make(1.0, 0.004, 0.02));
        near(
            base.reynolds_length().value(),
            4.0 * 0.02 * 0.004 / (2.0 * (0.02 + 0.004)),
        );
        let material = fluid(1000.0, 0.001);
        let reference = base.evaluate(&material, MassFlowRate::new(0.001)).unwrap();
        for scale in [1e-3, 3.0, 1e3] {
            // Q scales with length to hold Re fixed: dp scales as 1/length².
            let scaled = make(scale, 0.02, 0.004)
                .evaluate(&material, MassFlowRate::new(0.001 * scale))
                .unwrap();
            near(scaled.reynolds, reference.reynolds);
            near(
                scaled.pressure_drop_pa * scale * scale,
                reference.pressure_drop_pa,
            );
            near(
                scaled.dissipated_power_w * scale,
                reference.dissipated_power_w,
            );
        }
    }

    #[test]
    fn annular_pressure_matches_integrated_no_slip_velocity() {
        for ratio in [0.1_f64, 0.25, 0.5, 0.75, 0.9] {
            let law = MaterialFlowResistance::annular_laminar(
                Length::new(1.0), Length::new(ratio), Length::new(1.0), 100.0,
            ).unwrap();
            // Independently integrate the radial Stokes velocity, with both
            // stationary walls enforcing u=0. No conductance-series reuse.
            let count = 8192;
            let dr = (1.0 - ratio) / f64::from(count);
            let mut integral = 0.0;
            for j in 0..=count {
                let r = ratio + f64::from(j) * dr;
                let velocity = 0.25 * (1.0 - r * r
                    + (1.0 - ratio * ratio) * r.ln() / (1.0 / ratio).ln());
                let weight = if j == 0 || j == count { 1.0 }
                    else if j % 2 == 0 { 2.0 } else { 4.0 };
                integral += weight * r * velocity;
            }
            let flow_per_pressure = 2.0 * PI * dr * integral / 3.0;
            let state = law.evaluate(&fluid(1.0, 1.0), MassFlowRate::new(1.0)).unwrap();
            near(state.pressure_drop_pa, 1.0 / flow_per_pressure);
            near(state.dpressure_dmass_flow, state.pressure_drop_pa);
            near(state.dissipated_power_w, state.pressure_drop_pa);
            near(law.flow_area().value(), PI * (1.0 - ratio * ratio));
            near(law.reynolds_length().value(), 2.0 * (1.0 - ratio));
        }
    }

    #[test]
    fn annular_thin_gap_retains_the_cubic_parallel_plate_limit() {
        for relative_gap in [1e-6, 1e-9, 1e-12, f64::EPSILON] {
            let inner = 1.0 - relative_gap;
            let gap = 1.0 - inner;
            let mean = 0.5 * (1.0 + inner);
            let law = MaterialFlowResistance::annular_laminar(
                Length::new(2.0), Length::new(inner), Length::new(1.0), 100.0,
            ).unwrap();
            let q = 1e-6 * gap * gap;
            let state = law.evaluate(&fluid(1000.0, 0.001), MassFlowRate::new(1000.0 * q)).unwrap();
            // Plane Poiseuille between walls separated by gap, unrolled over
            // circumference 2*pi*mean. Curvature changes this by O(gap^2).
            let pressure = 6.0 * 0.001 * 2.0 * q / (PI * mean * gap.powi(3));
            near(state.pressure_drop_pa, pressure);
            near(state.reynolds, 1000.0 * q / (PI * mean * 0.001));
            assert!(state.pressure_drop_pa.is_finite() && state.pressure_drop_pa > 0.0);
        }
    }

    #[test]
    fn annular_vanishing_core_retains_logarithmic_drag_and_exact_empty_core_limit() {
        let circular = MaterialFlowResistance::circular_laminar(
            Length::new(1.0), Length::new(2.0), 100.0,
        ).unwrap();
        assert_eq!(circular, MaterialFlowResistance::annular_laminar(
            Length::new(1.0), Length::new(0.0), Length::new(1.0), 100.0,
        ).unwrap());
        let reference = circular.evaluate(&fluid(1.0, 1.0), MassFlowRate::new(1.0)).unwrap();
        let mut previous = f64::INFINITY;
        for inner in [1e-16_f64, 1e-64, 1e-200] {
            let law = MaterialFlowResistance::annular_laminar(
                Length::new(1.0), Length::new(inner), Length::new(1.0), 100.0,
            ).unwrap();
            let state = law.evaluate(&fluid(1.0, 1.0), MassFlowRate::new(1.0)).unwrap();
            near(state.pressure_drop_pa, reference.pressure_drop_pa / (1.0 + 1.0 / inner.ln()));
            assert!(state.pressure_drop_pa > reference.pressure_drop_pa);
            assert!(state.pressure_drop_pa < previous);
            previous = state.pressure_drop_pa;
        }
        // A representable dimensional law must not lose the inner wall merely
        // because the radius ratio underflows or its reciprocal overflows.
        let law = MaterialFlowResistance::annular_laminar(
            Length::new(1e100), Length::new(1e-250), Length::new(1e90), 100.0,
        ).unwrap();
        let state = law.evaluate(&fluid(1.0, 1.0), MassFlowRate::new(1e90)).unwrap();
        let pressure = (8.0 / PI) * 1e-170 / (1.0 - 1.0 / (340.0 * 10.0_f64.ln()));
        assert!((state.pressure_drop_pa / pressure - 1.0).abs() < 2e-12);
    }

    #[test]
    fn annular_signed_flow_and_si_rescaling_share_the_original_material_law() {
        let make = |scale| MaterialFlowResistance::annular_laminar(
            Length::new(0.4 * scale), Length::new(0.003 * scale),
            Length::new(0.005 * scale), 1000.0,
        ).unwrap();
        let material = fluid(1000.0, 0.001);
        let reference = make(1.0).evaluate(&material, MassFlowRate::new(0.001)).unwrap();
        for scale in [1e-3, 3.0, 1e3] {
            for sign in [-1.0_f64, 0.0, 1.0] {
                let value = make(scale).evaluate(&material, MassFlowRate::new(sign * 0.001 * scale)).unwrap();
                near(value.reynolds, sign.abs() * reference.reynolds);
                near(value.pressure_drop_pa * scale * scale, sign * reference.pressure_drop_pa);
                near(value.dissipated_power_w * scale, sign.abs() * reference.dissipated_power_w);
                near(value.dpressure_dmass_flow * scale.powi(3), reference.dpressure_dmass_flow);
            }
        }
        for (inner, outer) in [(-1.0, 1.0), (1.0, 1.0), (2.0, 1.0),
            (0.0, 0.0), (f64::NAN, 1.0), (0.0, f64::INFINITY)] {
            assert!(MaterialFlowResistance::annular_laminar(
                Length::new(1.0), Length::new(inner), Length::new(outer), 100.0,
            ).is_err());
        }
        assert!(make(1.0).evaluate(&material, MassFlowRate::new(10.0)).is_err());
    }

    #[test]
    fn invalid_geometry_and_out_of_regime_flow_refuse_without_regularization() {
        for value in [0.0, -1.0, f64::NAN, f64::INFINITY] {
            assert!(
                MaterialFlowResistance::circular_laminar(
                    Length::new(1.0),
                    Length::new(value),
                    100.0
                )
                .is_err()
            );
            assert!(
                MaterialFlowResistance::rectangular_laminar(
                    Length::new(value),
                    Length::new(1.0),
                    Length::new(1.0),
                    100.0
                )
                .is_err()
            );
            assert!(
                MaterialFlowResistance::circular_laminar(Length::new(1.0), Length::new(1.0), value)
                    .is_err()
            );
        }
        assert!(
            MaterialFlowResistance::circular_laminar(Length::new(1.0), Length::new(1.0), 2300.1)
                .is_err()
        );
        assert!(
            MaterialFlowResistance::rectangular_laminar(
                Length::new(1.0),
                Length::new(1e-200),
                Length::new(1e-200),
                100.0
            )
            .is_err()
        );
        let law =
            MaterialFlowResistance::circular_laminar(Length::new(1.0), Length::new(0.01), 100.0)
                .unwrap();
        assert!(matches!(
            law.evaluate(&fluid(1000.0, 0.001), MassFlowRate::new(1.0)),
            Err(ResistanceError::Reynolds { .. })
        ));
    }
}
