//! Explicit continuous pipe friction, including its finite Poiseuille limit.
//!
//! Churchill (1977), "Friction-factor equation spans all fluid-flow regimes",
//! equations (5), (18): the paper's friction factor is Darcy f / 8.
//! <https://files.engineering.com/files/85c0f3a6-a102-4a22-9d35-f15858c0dd2b/CEM_-_Friction-factor_equation_%281977%29.pdf>

use fs_material::fluid_transport::FluidTransportState;
use fs_math::det;
use fs_qty::{Length, MassFlowRate};

use super::{MaterialFlowResistance, PipeFriction, ResistanceError, positive};

impl MaterialFlowResistance {
    /// Explicit Churchill (1977) circular-pipe correlation across laminar,
    /// transitional and turbulent flow, including the exact creeping-flow limit.
    ///
    /// This is an empirical steady Newtonian pipe model. Transition depends on
    /// disturbances and upstream conditions; this smooth engineering fit does
    /// not resolve intermittency or predict its onset. Fully developed flow,
    /// negligible compressibility and an equivalent sand-grain roughness are
    /// caller assumptions; entrance and unsteady effects are not included.
    ///
    /// The authored Reynolds interval must lie within [0, 1e8]; roughness / D
    /// must lie within [0, 0.05]. A zero lower endpoint explicitly admits rest.
    /// Density and dynamic viscosity come from the original material snapshot.
    /// Pressure, its exact fixed-material tangent and nonnegative pressure work
    /// use the same flow-dependent correlation. Existing laminar and Colebrook
    /// constructors retain their own domains; no regime is silently switched.
    ///
    /// Existing series and native steady graph owners retain this nonlinear
    /// law. Constant-coefficient exports reject it instead of freezing friction.
    ///
    /// # Errors
    /// Invalid/unrepresentable geometry or support outside the declared domain.
    pub fn circular_churchill(
        length: Length,
        diameter: Length,
        roughness: Length,
        reynolds_range: (f64, f64),
    ) -> Result<Self, ResistanceError> {
        let d = diameter.value();
        let roughness = roughness.value();
        if !positive(d) || !roughness.is_finite() || roughness < 0.0
            || roughness > 0.05 * d
            || !reynolds_range.0.is_finite() || reynolds_range.0 < 0.0
            || !reynolds_range.1.is_finite() || reynolds_range.0 >= reynolds_range.1
            || reynolds_range.1 > 1e8
        {
            return Err(ResistanceError::Invalid("Churchill pipe geometry or Reynolds support"));
        }
        // Retain the original exact Poiseuille geometry. Only its explicitly
        // selected friction model and Reynolds support change.
        let mut law = Self::circular_laminar(length, diameter, 2300.0)?;
        let relative_roughness = roughness / d;
        if roughness > 0.0 && relative_roughness == 0.0 {
            return Err(ResistanceError::Invalid("Churchill roughness is unrepresentable"));
        }
        law.reynolds_range = reynolds_range;
        law.pipe_friction = Some(PipeFriction::Churchill(relative_roughness));
        Ok(law)
    }
}

/// Return g = f_D Re / 64 and d ln(g) / d ln(Re). Scaling the published
/// equation by its Poiseuille term avoids evaluating a singular f_D at rest.
pub(super) fn laminar_multiplier(reynolds: f64, roughness: f64) -> (f64, f64) {
    // The relative correction is below 2e-122 throughout this interval.
    // Returning its represented limit also avoids enormous inverse powers.
    if reynolds <= 1.0 {
        return (1.0, 0.0);
    }
    let z = det::pow(7.0 / reynolds, 0.9);
    let argument = z + 0.27 * roughness;
    let a = -2.457 * det::ln(argument);
    let da = (2.457 * 0.9) * (z / argument);
    let a2 = a * a;
    let a4 = a2 * a2;
    let a8 = a4 * a4;
    let a15 = ((a8 * a4) * a2) * a;
    let b = 37530.0 / reynolds;
    let b2 = b * b;
    let b4 = b2 * b2;
    let b8 = b4 * b4;
    let b16 = b8 * b8;
    let sum = a8 * a8 + b16;
    let root8 = det::sqrt(det::sqrt(det::sqrt(sum)));
    let ratio = reynolds / (8.0 * root8);
    let ratio2 = ratio * ratio;
    let ratio4 = ratio2 * ratio2;
    let ratio12 = (ratio4 * ratio4) * ratio4;
    let one_plus = 1.0 + ratio12;
    // Correct the addition rounding in ln(1+x), retaining weak viscous flow.
    let log1p = if one_plus == 1.0 { ratio12 }
        else { det::ln(one_plus) * (ratio12 / (one_plus - 1.0)) };
    let multiplier = 1.0 + det::expm1(log1p / 12.0);
    let slope = (ratio12 / one_plus)
        * (1.0 - 2.0 * ((a15 * da - b16) / sum));
    (multiplier, slope)
}

/// Fixed-work safeguarded inverse inside the native graph's original admitted
/// Reynolds endpoints. Every iterate uses the same pressure law and tangent.
pub(super) fn invert(
    law: &MaterialFlowResistance,
    transport: FluidTransportState,
    rho: f64,
    pressure: f64,
    volume_bounds: (f64, f64),
) -> Result<f64, ResistanceError> {
    let (mut lower, mut upper) = volume_bounds;
    if !positive(rho) || !positive(pressure) || !lower.is_finite() || lower < 0.0
        || !positive(upper) || lower >= upper
    {
        return Err(ResistanceError::Invalid("Churchill inverse requires positive pressure and admitted bounds"));
    }
    let mut q = (pressure / (law.viscous_geometry * transport.dynamic_viscosity_pa_s()))
        .clamp(lower, upper);
    for _ in 0..64 {
        let value = law.evaluate_snapshot(transport, MassFlowRate::new(q * rho), rho)?;
        let residual = value.pressure_drop_pa - pressure;
        if residual.abs() <= (32.0 * f64::EPSILON) * pressure {
            return Ok(q);
        }
        if residual > 0.0 { upper = q; } else { lower = q; }
        let newton = q - (residual / value.dpressure_dmass_flow) / rho;
        q = if newton > lower && newton < upper { newton }
            else { lower + 0.5 * (upper - lower) };
    }
    Err(ResistanceError::Invalid("Churchill inverse did not close within its fixed work bound"))
}

#[cfg(test)]
mod tests {
    use super::*;

    fn law(roughness: f64) -> MaterialFlowResistance {
        MaterialFlowResistance::circular_churchill(
            Length::new(2.0), Length::new(0.01), Length::new(roughness * 0.01),
            (0.0, 1e8),
        ).unwrap()
    }
    fn fluid(mu: f64) -> FluidTransportState {
        FluidTransportState::try_new(1000.0, 4000.0, 0.6, mu).unwrap()
    }
    // Independent direct transcription of Churchill equations (5), (18),
    // converted to the Darcy convention; no production scaling or helpers.
    fn published_darcy(re: f64, roughness: f64) -> f64 {
        let a = (-2.457 * ((7.0 / re).powf(0.9) + 0.27 * roughness).ln()).powi(16);
        let b = (37530.0 / re).powi(16);
        8.0 * ((8.0 / re).powi(12) + (a + b).powf(-1.5)).powf(1.0 / 12.0)
    }

    #[test]
    fn churchill_original_equation_and_exact_tangent_through_transition() {
        for roughness in [0.0, 1e-4, 0.01, 0.05] {
            let pipe = law(roughness);
            let area = core::f64::consts::PI * 0.01_f64.powi(2) / 4.0;
            for re in [1.0, 100.0, 1000.0, 2100.0, 2300.0, 3000.0, 4000.0, 1e4, 1e6, 9.99e7] {
                let mass = re * 0.001 * area / 0.01;
                let speed = mass / (1000.0 * area);
                let expected = published_darcy(re, roughness) * (2.0 / 0.01)
                    * 1000.0 * speed * speed / 2.0;
                let value = pipe.evaluate(&fluid(0.001), MassFlowRate::new(mass)).unwrap();
                assert!((value.pressure_drop_pa / expected - 1.0).abs() < 3e-13);
                let h = mass * 1e-5;
                let plus = pipe.evaluate(&fluid(0.001), MassFlowRate::new(mass + h)).unwrap();
                let minus = pipe.evaluate(&fluid(0.001), MassFlowRate::new(mass - h)).unwrap();
                let derivative = (plus.pressure_drop_pa - minus.pressure_drop_pa) / (2.0 * h);
                assert!((derivative / value.dpressure_dmass_flow - 1.0).abs() < 2e-8);
                assert!(value.dpressure_dmass_flow > 0.0);
                let reverse = pipe.evaluate(&fluid(0.001), MassFlowRate::new(-mass)).unwrap();
                assert_eq!(reverse.pressure_drop_pa, -value.pressure_drop_pa);
                assert_eq!(reverse.dissipated_power_w, value.dissipated_power_w);
                assert_eq!(reverse.dpressure_dmass_flow, value.dpressure_dmass_flow);
            }
        }
        let (g, _) = laminar_multiplier(1e5, 1e-4);
        assert!((64.0 * g / 1e5 - 0.018_462_624_566_280_075).abs() < 1e-15);
    }

    #[test]
    fn churchill_rest_and_weak_flow_recover_original_poiseuille_bits() {
        let pipe = law(0.01);
        let laminar = MaterialFlowResistance::circular_laminar(
            Length::new(2.0), Length::new(0.01), 2300.0,
        ).unwrap();
        for mass in [-1e-6, -1e-12, 0.0, 1e-12, 1e-6] {
            let actual = pipe.evaluate(&fluid(0.001), MassFlowRate::new(mass)).unwrap();
            let expected = laminar.evaluate(&fluid(0.001), MassFlowRate::new(mass)).unwrap();
            assert_eq!(actual.pressure_drop_pa, expected.pressure_drop_pa);
            assert_eq!(actual.dpressure_dmass_flow, expected.dpressure_dmass_flow);
            assert_eq!(actual.dissipated_power_w, expected.dissipated_power_w);
        }
        assert!(pipe.volume_coefficients(fluid(0.001), 1000.0).is_err());
    }

    #[test]
    fn churchill_viscosity_changes_regime_without_replacing_the_original_law() {
        let pipe = law(0.0);
        let area = core::f64::consts::PI * 0.01_f64.powi(2) / 4.0;
        let mass = 3000.0 * 0.001 * area / 0.01;
        for (mu, re) in [(0.003, 1000.0), (0.001, 3000.0), (0.0003, 10000.0)] {
            let state = fluid(mu);
            let value = pipe.evaluate(&state, MassFlowRate::new(mass)).unwrap();
            assert!((value.reynolds / re - 1.0).abs() < 1e-14);
            let speed = mass / (1000.0 * area);
            let pressure = published_darcy(re, 0.0) * 200.0 * 500.0 * speed * speed;
            assert!((value.pressure_drop_pa / pressure - 1.0).abs() < 3e-13);
            assert_eq!(value.transport, state);
            assert_eq!(value.law, pipe);
        }
        for band in [(-1.0, 1e8), (0.0, 0.0), (0.0, 1.01e8), (1e4, 100.0)] {
            assert!(MaterialFlowResistance::circular_churchill(
                Length::new(1.0), Length::new(1.0), Length::new(0.0), band,
            ).is_err());
        }
        for roughness in [-1.0, 0.051, f64::NAN, f64::INFINITY] {
            assert!(MaterialFlowResistance::circular_churchill(
                Length::new(1.0), Length::new(1.0), Length::new(roughness), (0.0, 1e8),
            ).is_err());
        }
        let restricted = MaterialFlowResistance::circular_churchill(
            Length::new(2.0), Length::new(0.01), Length::new(0.0), (2100.0, 4000.0),
        ).unwrap();
        assert!(matches!(restricted.evaluate(&fluid(0.003), MassFlowRate::new(mass)),
            Err(ResistanceError::Reynolds { .. })));
    }
}
