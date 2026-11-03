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
