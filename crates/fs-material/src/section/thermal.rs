//! Per-layer thermal eigenstrains for the existing parallel-strain section.
//!
//! Each expansion path is integrated by `state_point`, not approximated from
//! an endpoint alpha. The current E/rho and each expansion endpoint must name
//! the same card, coordinates and coordinate conventions. No metal is selected
//! here. For concentric, circumferentially uniform layers there is no thermal
//! bending moment in this one-dimensional model.
//!
//! This is small-strain uniaxial thermoelasticity with a common stress-free
//! reference temperature and frozen section geometry. Radial/Poisson constraint,
//! interlayer slip, yielding, temperature evolution and thermodynamic heat/work
//! closure are not modeled. A strain limit is an applicability assumption, not
//! an experimentally established yield criterion.

/// Fixed, force-controlled and elastic axial support equations.
pub mod constraint;
pub use constraint::AxialConstraint;

use crate::state_point::{IntegratedIsotropicThermalExpansion, MaterialStatePointError};
use super::{ConcentricElasticSection, invalid};

/// A source-bound thermal path for every layer of one elastic section.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricThermoelasticSection {
    section: ConcentricElasticSection,
    expansions: Vec<IntegratedIsotropicThermalExpansion>,
    thermal_force_n: f64,
    zero_resultant_strain: f64,
}

impl ConcentricThermoelasticSection {
    /// Bind paths in spatial layer order, without averaging their material data.
    /// Current temperatures may differ by layer; the common reference is the
    /// initially stress-free bonded configuration. Geometry and elastic data are
    /// supplied at the evaluated state, not transported from that reference.
    ///
    /// # Errors
    /// Refuses missing/extra paths, foreign cards, stale current coordinates,
    /// different coordinate conventions or reference temperatures, and arithmetic
    /// overflow/underflow. No partial source-bound section is published.
    pub fn try_new(
        section: ConcentricElasticSection,
        expansions: Vec<IntegratedIsotropicThermalExpansion>,
    ) -> Result<Self, MaterialStatePointError> {
        if expansions.len() != section.layers().len() {
            return Err(invalid("thermal section needs exactly one expansion path per layer"));
        }
        let reference = expansions[0].reference().resolved();
        let reference_temperature = temperature(reference.query_point())?;
        let mut thermal_force_n = 0.0;
        for ((layer, integral), path) in section.layers().iter()
            .zip(section.layer_integrals()).zip(&expansions)
        {
            let current = path.current().resolved();
            let path_reference = path.reference().resolved();
            if current.card_identity() != layer.material.card_identity()
                || current.axis_quantities() != layer.material.axis_quantities()
                || !same_point(current.query_point(), layer.material.query_point())
            {
                return Err(invalid("thermal section path must match its layer card and current state exactly"));
            }
            if temperature(path_reference.query_point())?.to_bits() != reference_temperature.to_bits()
                || path_reference.axis_quantities().get("T") != reference.axis_quantities().get("T")
            {
                return Err(invalid("thermal section requires a common stress-free reference temperature convention"));
            }
            let strain = path.free_linear_strain();
            let force = checked_product(integral.axial_stiffness_n(), strain)?;
            thermal_force_n = finite(thermal_force_n + force)?;
        }
        let zero_resultant_strain = finite(thermal_force_n / section.axial_stiffness_n())?;
        if thermal_force_n != 0.0 && zero_resultant_strain == 0.0 {
            return Err(invalid("thermal section free extension is unrepresentable"));
        }
        Ok(Self { section, expansions, thermal_force_n, zero_resultant_strain })
    }

    /// Original elastic section, including each layer's source bundle.
    #[must_use]
    pub const fn section(&self) -> &ConcentricElasticSection { &self.section }

    /// Original integrated expansion paths, in spatial layer order.
    #[must_use]
    pub fn expansion_paths(&self) -> &[IntegratedIsotropicThermalExpansion] { &self.expansions }

    /// Thermal term `sum(E_i A_i epsilon_th_i)` in `N = EA epsilon - N_th` [N].
    #[must_use]
    pub const fn thermal_force_n(&self) -> f64 { self.thermal_force_n }

    /// Stiffness-weighted free extension of the assembled member.
    /// Zero resultant is NOT zero layer stress: differential expansion can
    /// leave self-equilibrated tension and compression inside the section.
    #[must_use]
    pub const fn zero_resultant_axial_strain(&self) -> f64 { self.zero_resultant_strain }

    /// Evaluate every layer at one imposed total axial engineering strain.
    /// Compression is valid here; a taut-string consumer must separately refuse
    /// a nonpositive resultant. Every layer's thermal AND elastic strain is
    /// bounded, so a small net force cannot conceal excessive opposing stresses.
    ///
    /// # Errors
    /// Refuses nonfinite strain, a limit outside `(0, 1)`, any total/thermal/
    /// elastic strain outside that limit, and unrepresentable stresses or energy.
    pub fn response_at_axial_strain(
        &self,
        total_strain: f64,
        linear_strain_limit: f64,
    ) -> Result<ConcentricAxialResponse, MaterialStatePointError> {
        let mut layers = Vec::with_capacity(self.expansions.len());
        let (force_n, energy_per_length) = self.evaluate_axial(
            total_strain, linear_strain_limit, |layer| layers.push(layer),
        )?;
        Ok(ConcentricAxialResponse { total_strain, layers, force_n, energy_per_length })
    }

    /// Axial force [N] and stored energy per length [J/m], without allocating
    /// per-layer output. This uses the exact same layer arithmetic and ALL the
    /// strain/stress/representability checks as `response_at_axial_strain`.
    /// Useful for constitutive evaluations inside repeated residual/tangent
    /// probes; it does not replace layer stresses with an averaged material.
    pub fn axial_force_and_energy_at_strain(
        &self, total_strain: f64, linear_strain_limit: f64,
    ) -> Result<(f64, f64), MaterialStatePointError> {
        self.evaluate_axial(total_strain, linear_strain_limit, |_| {})
    }

    fn evaluate_axial(
        &self, total_strain: f64, linear_strain_limit: f64,
        mut emit: impl FnMut(ConcentricLayerResponse),
    ) -> Result<(f64, f64), MaterialStatePointError> {
        if !linear_strain_limit.is_finite() || linear_strain_limit <= 0.0
            || linear_strain_limit >= 1.0 || !total_strain.is_finite()
            || total_strain.abs() > linear_strain_limit
        {
            return Err(invalid("thermal section total strain or linear strain limit is inadmissible"));
        }
        let mut force_n = 0.0;
        let mut energy_per_length = 0.0;
        for (integral, path) in self.section.layer_integrals().iter().zip(&self.expansions) {
            let thermal_strain = path.free_linear_strain();
            let elastic_strain = total_strain - thermal_strain;
            if !elastic_strain.is_finite() || thermal_strain.abs() > linear_strain_limit
                || elastic_strain.abs() > linear_strain_limit
            {
                return Err(invalid("thermal section layer strain exceeds its declared linear regime"));
            }
            let force = checked_product(integral.axial_stiffness_n(), elastic_strain)?;
            let stress = finite(force / integral.area_m2())?;
            if force != 0.0 && stress == 0.0 {
                return Err(invalid("thermal section layer stress is unrepresentable"));
            }
            let energy = checked_product(0.5, checked_product(force, elastic_strain)?)?;
            force_n = finite(force_n + force)?;
            energy_per_length = finite(energy_per_length + energy)?;
            emit(ConcentricLayerResponse { thermal_strain, elastic_strain, force_n: force, stress_pa: stress });
        }
        Ok((force_n, energy_per_length))
    }
}

/// Uniaxial result in one layer; no transverse or interface stress is inferred.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ConcentricLayerResponse {
    thermal_strain: f64,
    elastic_strain: f64,
    force_n: f64,
    stress_pa: f64,
}
impl ConcentricLayerResponse {
    /// Integrated free thermal strain from the retained source path.
    #[must_use]
    pub const fn thermal_strain(&self) -> f64 { self.thermal_strain }
    /// Mechanical strain `epsilon_total - epsilon_th`.
    #[must_use]
    pub const fn elastic_strain(&self) -> f64 { self.elastic_strain }
    /// Signed layer force [N], tension positive.
    #[must_use]
    pub const fn force_n(&self) -> f64 { self.force_n }
    /// Signed uniaxial stress [Pa], tension positive.
    #[must_use]
    pub const fn stress_pa(&self) -> f64 { self.stress_pa }
}

/// All layer responses and their resultant at one imposed axial strain.
#[derive(Clone, Debug, PartialEq)]
pub struct ConcentricAxialResponse {
    total_strain: f64,
    layers: Vec<ConcentricLayerResponse>,
    force_n: f64,
    energy_per_length: f64,
}
impl ConcentricAxialResponse {
    /// Imposed engineering strain in the common reference configuration.
    #[must_use]
    pub const fn total_strain(&self) -> f64 { self.total_strain }
    /// Results in the original spatial layer order.
    #[must_use]
    pub fn layers(&self) -> &[ConcentricLayerResponse] { &self.layers }
    /// Sum of signed layer forces [N], not a positivity-clamped tension.
    #[must_use]
    pub const fn force_n(&self) -> f64 { self.force_n }
    /// Uniaxial stored elastic energy per model length [J/m]. This is not a
    /// coupled thermodynamic balance or dissipated heat.
    #[must_use]
    pub const fn strain_energy_per_length_j_m(&self) -> f64 { self.energy_per_length }
}

fn same_point(a: &[(String, f64)], b: &[(String, f64)]) -> bool {
    a.len() == b.len() && a.iter().zip(b)
        .all(|((ka, va), (kb, vb))| ka == kb && va.to_bits() == vb.to_bits())
}
fn temperature(point: &[(String, f64)]) -> Result<f64, MaterialStatePointError> {
    point.iter().find(|(key, _)| key == "T").map(|(_, value)| *value)
        .ok_or_else(|| invalid("thermal section path is missing its reference T coordinate"))
}
fn finite(value: f64) -> Result<f64, MaterialStatePointError> {
    if value.is_finite() { Ok(value) } else { Err(invalid("thermal section arithmetic is unrepresentable")) }
}
fn checked_product(a: f64, b: f64) -> Result<f64, MaterialStatePointError> {
    let product = finite(a * b)?;
    if a != 0.0 && b != 0.0 && product == 0.0 {
        return Err(invalid("thermal section nonzero product underflowed"));
    }
    Ok(product)
}
