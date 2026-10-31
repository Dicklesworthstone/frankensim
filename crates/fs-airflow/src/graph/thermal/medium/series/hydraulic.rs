//! Material-scaled linear/quadratic pressure loss, independent of a heat solver.
//!
//! `dp = a mu Q + b rho Q |Q|`. Geometry/calibration supplies `a` [1/m^3]
//! and `b` [1/m^4]; the shared material snapshot supplies rho and mu. This is
//! the lumped Darcy-Forchheimer form, not a second fluid database or a universal
//! pipe-friction correlation. The coefficients and Reynolds interval are
//! explicit authored assumptions. No roughness, entrance, phase, compressible
//! flow, or Reynolds-dependent turbulent friction law is inferred.
//!
//! Source form: COMSOL, Forchheimer Flow model documentation,
//! <https://doc.comsol.com/6.4/doc/com.comsol.help.models.ssf.forchheimer_flow/forchheimer_flow.html>.
//! Fixed-state resistance is dissipative; this alone is not a coupled-energy
//! certificate. The thermal driver must account for or explicitly neglect the
//! reported dissipation. A gas snapshot admits only a frozen-density estimate.

use fs_material::MaterialError;
use fs_material::fluid_transport::{FluidTransport, FluidTransportState};
use fs_qty::{Area, Length, MassFlowRate};

/// Source-preserving inversion of complete fixed-state series pressure losses.
pub mod pressure;

/// Immutable geometry/calibration, not fluid-specific pressure-loss constants.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct MaterialFlowResistance {
    viscous_geometry: f64,
    inertial_geometry: f64,
    flow_area: Area,
    reynolds_length: Length,
    reynolds_range: (f64, f64),
}

/// Original material refusal or an unsupported pressure-loss evaluation.
#[derive(Clone, Debug, PartialEq)]
pub enum ResistanceError {
    /// Malformed geometry, flow or unrepresentable arithmetic.
    Invalid(&'static str),
    /// The shared material snapshot could not be admitted.
    Material(MaterialError),
    /// Outside the explicitly authored Reynolds applicability interval.
    Reynolds { value: f64, admitted: (f64, f64) },
}
impl core::fmt::Display for ResistanceError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material flow resistance: {self:?}")
    }
}
impl std::error::Error for ResistanceError {}
impl From<MaterialError> for ResistanceError {
    fn from(value: MaterialError) -> Self { Self::Material(value) }
}

/// One numerical snapshot, retaining every input used by momentum and power.
/// Original source receipts remain on the input medium, not on this raw view.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ResistanceEvaluation {
    /// Original geometry and applicability declaration.
    pub law: MaterialFlowResistance,
    /// Exactly one checked material snapshot.
    pub transport: FluidTransportState,
    /// Density used for continuity and inertia. Normally the source density;
    /// an explicitly bounded Boussinesq solve uses its conserved reference rho.
    pub momentum_density_kg_m3: f64,
    /// Signed mass flow [kg/s].
    pub mass_flow_kg_per_s: f64,
    /// Signed local volume flow [m^3/s], using `momentum_density_kg_m3`.
    pub volume_flow_m3_s: f64,
    /// Nonnegative Reynolds number based on speed magnitude.
    pub reynolds: f64,
    /// Signed pressure drop in the declared forward direction [Pa].
    pub pressure_drop_pa: f64,
    /// Fixed-material derivative of pressure drop with respect to mass flow.
    /// This is NOT the derivative of a temperature-coupled operating point.
    pub dpressure_dmass_flow: f64,
    /// Nonnegative `dp * Q` [W]. Not silently added to a thermal source.
    pub dissipated_power_w: f64,
}

impl MaterialFlowResistance {
    /// Declare shape/calibration factors and their finite Reynolds support.
    /// A capillary owner may supply its viscous factor; an orifice/porous owner
    /// may supply its inertial factor. Neither factor includes rho or mu.
    /// Zero in one factor disables that mechanism explicitly, not by fallback.
    ///
    /// # Errors
    /// Both factors zero, negative/nonfinite factors, invalid area/length or band.
    pub fn try_new(
        viscous_geometry_per_m3: f64,
        inertial_geometry_per_m4: f64,
        flow_area: Area,
        reynolds_length: Length,
        reynolds_range: (f64, f64),
    ) -> Result<Self, ResistanceError> {
        let nonnegative = |v: f64| v.is_finite() && v >= 0.0;
        if !nonnegative(viscous_geometry_per_m3) || !nonnegative(inertial_geometry_per_m4)
            || (viscous_geometry_per_m3 == 0.0 && inertial_geometry_per_m4 == 0.0)
            || !positive(flow_area.value()) || !positive(reynolds_length.value())
            || !nonnegative(reynolds_range.0) || !positive(reynolds_range.1)
            || reynolds_range.0 > reynolds_range.1
        {
            return Err(ResistanceError::Invalid("invalid resistance factors, geometry or Reynolds support"));
        }
        Ok(Self { viscous_geometry: viscous_geometry_per_m3, inertial_geometry: inertial_geometry_per_m4,
            flow_area, reynolds_length, reynolds_range })
    }

    /// Flow-normal area used by velocity and Reynolds admission.
    #[must_use]
    pub const fn flow_area(&self) -> Area { self.flow_area }

    /// Evaluate using one snapshot shared with the thermal/flow material owners.
    /// Signed and zero flow are supported when their Reynolds number is admitted.
    ///
    /// # Errors
    /// Invalid material or flow, unsupported Reynolds number, or unrepresentable
    /// active loss, differential resistance or dissipated power. No partial output.
    pub fn evaluate<M: FluidTransport + ?Sized>(
        &self, medium: &M, mass_flow: MassFlowRate,
    ) -> Result<ResistanceEvaluation, ResistanceError> {
        if !mass_flow.value().is_finite() {
            return Err(ResistanceError::Invalid("nonfinite mass flow"));
        }
        let transport = medium.transport_state()?;
        self.evaluate_snapshot(transport, mass_flow, transport.density_kg_m3())
    }

    // One loss law for ordinary and Boussinesq hydraulics. The original source
    // snapshot stays visible; reference rho is only a declared momentum model.
    fn evaluate_snapshot(
        &self, transport: FluidTransportState, mass_flow: MassFlowRate, rho: f64,
    ) -> Result<ResistanceEvaluation, ResistanceError> {
        let mass = mass_flow.value();
        if !mass.is_finite() || !positive(rho) {
            return Err(ResistanceError::Invalid("finite mass flow and positive momentum density required"));
        }
        let mu = transport.dynamic_viscosity_pa_s();
        let q = mass / rho;
        let speed = q.abs() / self.flow_area.value();
        if !q.is_finite() || !speed.is_finite() || (mass != 0.0 && (q == 0.0 || speed == 0.0)) {
            return Err(ResistanceError::Invalid("volume flow or velocity is unrepresentable"));
        }
        let momentum = if rho == transport.density_kg_m3() { transport } else {
            FluidTransportState::try_new(rho, transport.specific_heat_j_kg_k(),
                transport.thermal_conductivity_w_m_k(), mu)?
        };
        let reynolds = momentum.reynolds(speed, self.reynolds_length.value())?;
        if reynolds < self.reynolds_range.0 || reynolds > self.reynolds_range.1 {
            return Err(ResistanceError::Reynolds { value: reynolds, admitted: self.reynolds_range });
        }
        let viscous = (self.viscous_geometry * mu) * q;
        let inertial = (self.inertial_geometry * mass.abs()) * q;
        let dp = viscous + inertial;
        let derivative = (self.viscous_geometry * mu) / rho
            + (2.0 * self.inertial_geometry * q.abs());
        let power = dp * q;
        if [viscous, inertial, dp, derivative, power].iter().any(|v| !v.is_finite())
            || derivative < 0.0 || power < 0.0
            || (self.viscous_geometry > 0.0 && derivative == 0.0)
            || (mass != 0.0 && (dp == 0.0 || derivative == 0.0 || power == 0.0
                || (self.viscous_geometry > 0.0 && viscous == 0.0)
                || (self.inertial_geometry > 0.0 && inertial == 0.0)))
        {
            return Err(ResistanceError::Invalid("pressure loss, tangent or power is unrepresentable"));
        }
        Ok(ResistanceEvaluation { law: *self, transport, momentum_density_kg_m3: rho,
            mass_flow_kg_per_s: mass,
            volume_flow_m3_s: q, reynolds, pressure_drop_pa: dp,
            dpressure_dmass_flow: derivative, dissipated_power_w: power })
    }
}
fn positive(v: f64) -> bool { v.is_finite() && v > 0.0 }

#[cfg(test)]
mod tests {
    use super::*;
    fn law(a: f64, b: f64) -> MaterialFlowResistance {
        MaterialFlowResistance::try_new(a, b, Area::new(1e-4), Length::new(0.01), (0.0, 1e6)).unwrap()
    }
    fn fluid(rho: f64, mu: f64) -> FluidTransportState {
        FluidTransportState::try_new(rho, 2000.0, 0.5, mu).unwrap()
    }
    fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-10*(1.0+a.abs().max(b.abs())), "{a} != {b}"); }

    #[test]
    fn data_changes_both_viscous_and_inertial_pressure_at_conserved_mass() {
        let r = law(1e10, 2e8);
        let a = r.evaluate(&fluid(1000.0, 0.002), MassFlowRate::new(0.01)).unwrap();
        let b = r.evaluate(&fluid(500.0, 0.001), MassFlowRate::new(0.01)).unwrap();
        near(a.pressure_drop_pa, 220.0);
        near(b.pressure_drop_pa, 240.0);
        near(b.volume_flow_m3_s, 2.0*a.volume_flow_m3_s);
        near(a.dissipated_power_w, 0.0022);
        near(b.dissipated_power_w, 0.0048);
    }
    #[test]
    fn signed_loss_is_passive_and_fixed_material_tangent_matches_perturbation() {
        let r = law(1e10, 2e8);
        let f = fluid(1000.0, 0.002);
        for m in [-0.01, 0.0, 0.01] {
            let v = r.evaluate(&f, MassFlowRate::new(m)).unwrap();
            let reverse = r.evaluate(&f, MassFlowRate::new(-m)).unwrap();
            near(v.pressure_drop_pa, -reverse.pressure_drop_pa);
            near(v.dissipated_power_w, reverse.dissipated_power_w);
            let eps = 1e-8;
            let plus = r.evaluate(&f, MassFlowRate::new(m+eps)).unwrap();
            let minus = r.evaluate(&f, MassFlowRate::new(m-eps)).unwrap();
            assert!(((plus.pressure_drop_pa-minus.pressure_drop_pa)/(2.0*eps)
                /v.dpressure_dmass_flow-1.0).abs() < 1e-6);
        }
        near(law(0.0, 2e8).evaluate(&f, MassFlowRate::new(0.0)).unwrap().dpressure_dmass_flow, 0.0);
    }
    #[test]
    fn samples_material_once_and_preserves_the_shared_snapshot() {
        struct Counted(std::cell::Cell<usize>);
        impl FluidTransport for Counted {
            fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
                self.0.set(self.0.get()+1); Ok(fluid(1000.0, 0.002))
            }
        }
        let f = Counted(std::cell::Cell::new(0));
        let value = law(1e10, 0.0).evaluate(&f, MassFlowRate::new(0.01)).unwrap();
        assert_eq!(f.0.get(), 1);
        assert_eq!(value.transport, fluid(1000.0, 0.002));
        near(value.reynolds, 500.0);
    }
    #[test]
    fn unsupported_states_and_invalid_coefficients_refuse() {
        for (a,b) in [(0.0,0.0),(-1.0,0.0),(0.0,f64::NAN),(f64::INFINITY,0.0)] {
            assert!(MaterialFlowResistance::try_new(a,b,Area::new(1e-4),Length::new(0.01),(0.0,1000.0)).is_err());
        }
        let r = MaterialFlowResistance::try_new(1e10,0.0,Area::new(1e-4),Length::new(0.01),(1.0,400.0)).unwrap();
        assert!(matches!(r.evaluate(&fluid(1000.0,0.002),MassFlowRate::new(0.01)),Err(ResistanceError::Reynolds { .. })));
        assert!(r.evaluate(&fluid(1000.0,0.002),MassFlowRate::new(0.0)).is_err());
        assert!(r.evaluate(&fluid(1000.0,0.002),MassFlowRate::new(f64::NAN)).is_err());
        assert!(law(f64::MAX,0.0).evaluate(&fluid(1000.0,2.0),MassFlowRate::new(0.01)).is_err());
        assert!(law(1e-200,0.0).evaluate(&fluid(1e-20,1e-200),MassFlowRate::new(0.0)).is_err());
    }
}
