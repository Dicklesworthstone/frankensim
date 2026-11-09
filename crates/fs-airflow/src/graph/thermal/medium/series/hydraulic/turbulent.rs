//! Explicit developed rough-pipe resistance on the original fluid snapshot.

use fs_math::det;
use super::{Area, Length, MaterialFlowResistance, PipeFriction, ResistanceError, positive};

impl MaterialFlowResistance {
    /// Fully developed Newtonian circular pipe, using Darcy-Weisbach loss and
    /// the implicit Colebrook-White friction factor at every material/flow trial.
    /// `roughness` is equivalent sand-grain height, not an RMS surface height.
    /// The authored Reynolds interval must lie within [4000, 1e8]; roughness / D
    /// must lie within [0, 0.05]. Flow must already be developed; pipe length
    /// alone does not establish that assumption. Transition, entrance effects,
    /// non-Newtonian stress and compressibility require separate models.
    ///
    /// `dp = f(Re, roughness/D) (L/D) rho v |v| / 2`. Density and viscosity are
    /// sampled together, and the mass-flow tangent differentiates that same
    /// implicit friction equation. Signed loss and `dp Q` remain passive.
    /// Existing liquid pressure, pump and viscous-heating paths accept this law.
    /// Native steady graphs retain this nonlinear law; fixed-coefficient exports
    /// reject it instead of freezing friction.
    ///
    /// Reference: Colebrook (1939), DOI 10.1680/ijoti.1939.13150;
    /// NASA KSC-STD-Z-0017 Rev A, section 7.4.1, equation 4:
    /// <https://standards.nasa.gov/sites/default/files/standards/KSC/A/0/KSC-STD-Z-0017_RevA_Final.pdf>.
    ///
    /// # Errors
    /// Invalid/unrepresentable geometry or support outside the declared domain.
    pub fn circular_turbulent(
        length: Length,
        diameter: Length,
        roughness: Length,
        reynolds_range: (f64, f64),
    ) -> Result<Self, ResistanceError> {
        let (length, diameter, roughness) = (length.value(), diameter.value(), roughness.value());
        if !positive(length) || !positive(diameter) || !roughness.is_finite()
            || roughness < 0.0 || roughness > 0.05 * diameter
            || !reynolds_range.0.is_finite() || !reynolds_range.1.is_finite()
            || reynolds_range.0 < 4000.0 || reynolds_range.0 >= reynolds_range.1
            || reynolds_range.1 > 1e8
        {
            return Err(ResistanceError::Invalid("developed turbulent circular pipe geometry or Reynolds support"));
        }
        let area = core::f64::consts::FRAC_PI_4 * diameter * diameter;
        let geometry = (length / diameter) / area / (2.0 * area);
        let relative_roughness = roughness / diameter;
        if !positive(geometry) || (roughness > 0.0 && relative_roughness == 0.0) {
            return Err(ResistanceError::Invalid("turbulent pipe geometry is unrepresentable"));
        }
        let mut law = Self::try_new(0.0, geometry, Area::new(area), Length::new(diameter), reynolds_range)?;
        law.pipe_friction = Some(PipeFriction::Colebrook(relative_roughness));
        Ok(law)
    }


    /// Numerical pressure-to-flow inverse of the original pipe law. Graph
    /// brackets may use it only inside their precomputed Reynolds endpoints;
    /// accepted flows must still pass `evaluate_snapshot`. Colebrook uses its
    /// algebraic inverse; Churchill uses bounded Newton/bisection with the same
    /// exact derivative as the forward law and the admitted volume-flow bounds.
    pub(crate) fn turbulent_trial_volume_flow(
        &self, transport: fs_material::fluid_transport::FluidTransportState,
        rho: f64, pressure_magnitude: f64,
        volume_bounds: (f64, f64),
    ) -> Result<f64, ResistanceError> {
        let roughness = match self.pipe_friction {
            Some(PipeFriction::Colebrook(roughness)) => roughness,
            Some(PipeFriction::Churchill(_)) => {
                return super::churchill::invert(self, transport, rho, pressure_magnitude, volume_bounds);
            }
            None => return Err(ResistanceError::Invalid("pipe friction required for nonlinear inverse")),
        };
        if !positive(rho) || !positive(pressure_magnitude) {
            return Err(ResistanceError::Invalid("positive density and pressure required for pipe inverse"));
        }
        // dp = (b rho) f Q^2 and Re = rho Q D / (mu A).
        // With t = Q sqrt(f), Colebrook gives Q directly:
        // Q = -2 t log10(eps/(3.7 D) + 2.51 mu A/(rho D t)).
        // Separate roots avoid forming a dimensional b*rho product.
        let t = (det::sqrt(pressure_magnitude) / det::sqrt(self.inertial_geometry))
            / det::sqrt(rho);
        let argument = roughness / 3.7
            + ((2.51 * transport.dynamic_viscosity_pa_s() / rho)
                * (self.flow_area.value() / self.reynolds_length.value())) / t;
        let q = (-2.0 * core::f64::consts::LOG10_E * det::ln(argument)) * t;
        if !positive(t) || !positive(argument) || !positive(q) {
            return Err(ResistanceError::Invalid("Colebrook pressure inverse is unrepresentable"));
        }
        Ok(q)
    }

    /// Equivalent sand roughness / diameter for the selected nonlinear pipe law.
    /// `None` denotes an authored constant linear/quadratic resistance.
    #[must_use]
    pub const fn relative_sand_roughness(&self) -> Option<f64> {
        match self.pipe_friction {
            Some(PipeFriction::Colebrook(roughness) | PipeFriction::Churchill(roughness)) => Some(roughness),
            None => None,
        }
    }
}

// Return Darcy f and d ln(f) / d ln(Re). The inverse-root residual is monotone:
// F(x) = x + 2 log10(eps/(3.7 D) + 2.51 x/Re), x = 1/sqrt(f).
// Differentiating F=0 gives the exact logarithmic slope -2c/(1+c), with
// c = (2/ln(10)) (2.51/Re) / (eps/(3.7D) + 2.51x/Re).
// This bounded scalar solve is allocation-free. No finite differences, lookup
// table or unrelated explicit friction approximation enters the tangent.
pub(super) fn colebrook(reynolds: f64, roughness: f64) -> Result<(f64, f64), ResistanceError> {
    let a = roughness / 3.7;
    let b = 2.51 / reynolds;
    let scale = 2.0 * core::f64::consts::LOG10_E;
    // F(1)<0 and F(32)>0 throughout the admitted roughness/Reynolds domain.
    let (mut lower, mut upper) = (1.0, 32.0);
    let mut x = 7.0;
    for _ in 0..32 {
        let argument = a + b * x;
        let residual = x + scale * det::ln(argument);
        let c = scale * b / argument;
        if residual.abs() <= 8.0 * f64::EPSILON * x {
            return Ok((1.0 / (x * x), -2.0 * c / (1.0 + c)));
        }
        if residual < 0.0 { lower = x; } else { upper = x; }
        let candidate = x - residual / (1.0 + c);
        x = if candidate > lower && candidate < upper { candidate }
            else { lower + 0.5 * (upper - lower) };
    }
    Err(ResistanceError::Invalid("Colebrook friction did not close within its bounded iteration budget"))
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fluid_transport::FluidTransportState;
    use fs_qty::MassFlowRate;

    fn law(roughness: f64) -> MaterialFlowResistance {
        MaterialFlowResistance::circular_turbulent(Length::new(2.0), Length::new(0.01),
            Length::new(roughness * 0.01), (4000.0, 1e8)).unwrap()
    }
    fn fluid(rho: f64, mu: f64) -> FluidTransportState {
        FluidTransportState::try_new(rho, 4000.0, 0.6, mu).unwrap()
    }

    #[test]
    fn original_colebrook_equation_and_implicit_tangent_across_roughness_and_reynolds() {
        for roughness in [0.0, 1e-8, 1e-4, 0.01, 0.05] {
            let pipe = law(roughness);
            let medium = fluid(1000.0, 0.001);
            for re in [4001.0, 1e4, 1e5, 1e6, 9.99e7] {
                let mass = re * 0.001 * pipe.flow_area().value() / 0.01;
                let actual = pipe.evaluate(&medium, MassFlowRate::new(mass)).unwrap();
                let velocity = mass / (1000.0 * pipe.flow_area().value());
                let friction = actual.pressure_drop_pa / ((2.0 / 0.01) * 1000.0 * velocity * velocity / 2.0);
                // Independent original f-variable equation, using platform log.
                let inverse = 1.0 / friction.sqrt();
                let residual = inverse + 2.0 * (roughness / 3.7 + 2.51 * inverse / re).log10();
                assert!(residual.abs() < 2e-12, "{roughness} {re}: {residual}");
                let step = mass * 1e-5;
                let plus = pipe.evaluate(&medium, MassFlowRate::new(mass + step)).unwrap();
                let minus = pipe.evaluate(&medium, MassFlowRate::new(mass - step)).unwrap();
                let tangent = (plus.pressure_drop_pa - minus.pressure_drop_pa) / (2.0 * step);
                assert!((actual.dpressure_dmass_flow / tangent - 1.0).abs() < 1e-8);
                // Constant-f differentiation would overestimate this tangent.
                assert!(actual.dpressure_dmass_flow < 2.0 * actual.pressure_drop_pa / mass);
                let reverse = pipe.evaluate(&medium, MassFlowRate::new(-mass)).unwrap();
                assert_eq!(reverse.pressure_drop_pa, -actual.pressure_drop_pa);
                assert_eq!(reverse.dissipated_power_w, actual.dissipated_power_w);
                assert_eq!(reverse.dpressure_dmass_flow, actual.dpressure_dmass_flow);
                assert!(actual.dissipated_power_w > 0.0);
            }
        }
        let (smooth, _) = colebrook(1e5, 0.0).unwrap();
        assert!((smooth - 0.017_989_773_084_273_84).abs() < 2e-14);
        let (rough, slope) = colebrook(1e8, 0.05).unwrap();
        let limit = 1.0 / (-2.0 * (0.05_f64 / 3.7).log10()).powi(2);
        assert!((rough / limit - 1.0).abs() < 5e-6);
        assert!(slope.abs() < 5e-6);
    }

    #[test]
    fn current_material_changes_friction_and_power_without_refitting_geometry() {
        let pipe = law(1e-4);
        let mass = MassFlowRate::new(0.2);
        let base = pipe.evaluate(&fluid(1000.0, 0.001), mass).unwrap();
        let viscous = pipe.evaluate(&fluid(1000.0, 0.002), mass).unwrap();
        let lighter = pipe.evaluate(&fluid(500.0, 0.001), mass).unwrap();
        assert!(viscous.pressure_drop_pa > base.pressure_drop_pa);
        assert!(viscous.reynolds < base.reynolds);
        assert_eq!(base.law, viscous.law);
        assert_eq!(lighter.reynolds, base.reynolds);
        assert_eq!(lighter.pressure_drop_pa, 2.0 * base.pressure_drop_pa);
        assert_eq!(lighter.dissipated_power_w, 4.0 * base.dissipated_power_w);
        assert!(pipe.volume_coefficients(base.transport, 1000.0).is_err());
    }

    #[test]
    fn refuses_transition_extrapolation_and_roughness_substitution() {
        let pipe = law(0.001);
        let medium = fluid(1000.0, 0.001);
        for re in [0.0, 2300.0, 3999.0, 1.01e8] {
            let mass = re * 0.001 * pipe.flow_area().value() / 0.01;
            assert!(matches!(pipe.evaluate(&medium, MassFlowRate::new(mass)),
                Err(ResistanceError::Reynolds { .. })));
        }
        for roughness in [-0.001, 0.051, f64::NAN, f64::INFINITY] {
            assert!(MaterialFlowResistance::circular_turbulent(Length::new(2.0), Length::new(1.0),
                Length::new(roughness), (4000.0, 1e8)).is_err());
        }
        for band in [(2300.0, 1e8), (4000.0, 1.01e8), (1e5, 1e4)] {
            assert!(MaterialFlowResistance::circular_turbulent(Length::new(2.0), Length::new(1.0),
                Length::new(0.001), band).is_err());
        }
    }
}
