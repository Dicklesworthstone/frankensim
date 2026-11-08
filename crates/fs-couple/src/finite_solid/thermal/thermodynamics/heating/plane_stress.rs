//! Actual heat transactions with fixed in-plane deformation and free transverse
//! tractions. The native exact-volume owner determines thickness, shear and
//! pressure at every temperature on the original source paths.
//!
//! With the in-plane deformation held, physical P:Fdot is zero on this branch:
//! the held supports do no work and transverse tractions vanish. Heat therefore
//! changes the constrained internal energy. Its entropy and capacity include
//! pressure work from the changing natural volume. This uniform material
//! inventory has no inferred spatial conduction, inertia or caloric phase change.

use super::{
    CaloricReference, Deformation, FiniteSolidMaterial, FiniteThermoelasticHeatSolve,
    FiniteThermoelasticSourceError, ReferenceDeformation, ResolvedHeatCapacityCurve,
    ResolvedNeoHookeanThermalCoefficients, Result, Source, ThermalReference, combined_capacity,
    expansion_jet_at_temperature, finite, invalid, kernel, poll,
};
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_material::hyper::thermoelastic::incompressible_plane_stress::IncompressibleThermoelasticPlaneStressState;
use fs_thermal::enthalpy::ThermalStorageBasis;
use std::sync::Arc;

/// Accepted temperature and native exact-volume plane-stress response.
/// Every descendant retains the original source identity, in-plane constraint,
/// reference mass and reference volume. Current density follows actual det(F).
#[derive(Clone, Debug)]
pub struct IncompressiblePlaneCaloricState {
    source: Arc<Source>,
    potential: IncompressibleThermoelasticPlaneStressState,
    base_specific_capacity_j_kg_k: f64,
    heat_capacity_j_per_k: f64,
    jacobian: f64,
}

/// Accepted signed heat, final deformation and independently checked U/S changes.
#[derive(Clone, Debug)]
pub struct IncompressiblePlaneHeatBalance {
    result: kernel::Balance<IncompressiblePlaneCaloricState>,
}

impl FiniteSolidMaterial {
    /// Prepare heat-induced exact-volume plane stress with fixed in-plane F.
    ///
    /// The original engineering expansion, coefficient and caloric paths must
    /// describe the same material and actual temperature. Original mass stays
    /// fixed; current source density is not substituted for reference inventory.
    /// The full original capacity is combined with the signed constrained
    /// mechanical contribution, and their sum must be strictly positive.
    ///
    /// The original two-shear solve and source/domain settings apply at every
    /// temperature. This is a material inventory under explicit mechanical
    /// restraints, not a solved spatial temperature or displacement field.
    #[allow(clippy::too_many_arguments)]
    pub fn incompressible_plane_stress_caloric_state(
        &self,
        in_plane: [f64; 4],
        before_thermal: ReferenceDeformation,
        reference_mass_kg: f64,
        coefficients: &ResolvedNeoHookeanThermalCoefficients,
        caloric: &ResolvedHeatCapacityCurve,
        caloric_reference: CaloricReference,
        settings: PlaneStressSettings,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<IncompressiblePlaneCaloricState> {
        let source = self.prepare_heat_source(
            Deformation::ExactPlane { in_plane, settings },
            before_thermal,
            reference_mass_kg,
            coefficients,
            caloric,
            caloric_reference,
            &mut cancelled,
        )?;
        let state = evaluate(&source, coefficients.temperature_k())?;
        poll(&mut cancelled)?;
        Ok(state)
    }
}

fn evaluate(source: &Arc<Source>, temperature: f64) -> Result<IncompressiblePlaneCaloricState> {
    let Some(ThermalReference::Expansion(expansion)) = source.material.thermal_reference() else {
        unreachable!("prepared thermal source remains immutable")
    };
    let Deformation::ExactPlane { in_plane, settings } = source.deformation else {
        return Err(invalid(
            "plane-stress heat requires its original in-plane inventory",
        ));
    };
    let jet = expansion_jet_at_temperature(expansion, source.before_thermal, temperature, false)?;
    let potential = source
        .coefficients
        .incompressible_plane_stress_state_at_temperature(
            temperature,
            &in_plane,
            jet,
            None,
            settings,
        )
        .map_err(FiniteThermoelasticSourceError::Coefficients)?;
    let (base_specific_capacity_j_kg_k, heat_capacity_j_per_k) = combined_capacity(
        source,
        temperature,
        potential.mechanical_heat_capacity_j_m3_k(),
    )?;
    let f = potential.response().deformation_gradient;
    let jacobian = finite((f[0] * f[4] - f[1] * f[3]) * f[8])?;
    if jacobian <= 0.0
        || finite(source.volume_m3 * jacobian)? <= 0.0
        || finite(source.material.reference_density_kg_m3() / jacobian)? <= 0.0
    {
        return Err(invalid(
            "plane-stress current volume and density must be representable",
        ));
    }
    Ok(IncompressiblePlaneCaloricState {
        source: Arc::clone(source),
        potential,
        base_specific_capacity_j_kg_k,
        heat_capacity_j_per_k,
        jacobian,
    })
}

impl IncompressiblePlaneCaloricState {
    pub const fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    /// Original anchor and density. Evolved physical stress is in thermodynamics().
    pub fn material(&self) -> &FiniteSolidMaterial {
        &self.source.material
    }
    pub fn coefficients(&self) -> &ResolvedNeoHookeanThermalCoefficients {
        &self.source.coefficients
    }
    pub fn caloric_source(&self) -> &ResolvedHeatCapacityCurve {
        &self.source.caloric
    }
    pub fn caloric_reference(&self) -> CaloricReference {
        self.source.caloric_reference
    }
    pub fn reference_mass_kg(&self) -> f64 {
        self.source.mass_kg
    }
    pub fn reference_volume_m3(&self) -> f64 {
        self.source.volume_m3
    }
    pub fn current_volume_m3(&self) -> f64 {
        self.source.volume_m3 * self.jacobian
    }
    pub fn current_density_kg_m3(&self) -> f64 {
        self.source.material.reference_density_kg_m3() / self.jacobian
    }
    pub fn temperature_k(&self) -> f64 {
        self.potential.temperature_k()
    }
    pub const fn thermodynamics(&self) -> &IncompressibleThermoelasticPlaneStressState {
        &self.potential
    }
    pub const fn heat_capacity_j_per_k(&self) -> f64 {
        self.heat_capacity_j_per_k
    }
    pub const fn base_specific_capacity_j_kg_k(&self) -> f64 {
        self.base_specific_capacity_j_kg_k
    }

    /// Evaluate an original source endpoint, without claiming heat or elapsed time.
    pub fn at_temperature(
        &self,
        temperature: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let state = evaluate(&self.source, temperature)?;
        poll(&mut cancelled)?;
        Ok(state)
    }

    fn same_source(&self, other: &Self) -> Result<()> {
        if Arc::ptr_eq(&self.source, &other.source) {
            Ok(())
        } else {
            Err(invalid(
                "plane heat endpoints require the same original constrained inventory",
            ))
        }
    }

    /// Stable native constrained energy increment plus the original Cp integral.
    /// No absolute enthalpy datum enters this closed-inventory heat balance.
    pub fn internal_energy_change_to(&self, other: &Self) -> Result<f64> {
        self.same_source(other)?;
        let base = self
            .source
            .caloric
            .specific_enthalpy_change(self.temperature_k(), other.temperature_k())?;
        let mechanical = self
            .potential
            .internal_energy_increment_to(&other.potential)
            .map_err(FiniteThermoelasticSourceError::Material)?;
        finite(
            self.source
                .mass_kg
                .mul_add(base, self.source.volume_m3 * mechanical),
        )
    }

    /// Original Cp/T integral and the pressure-corrected mechanical entropy.
    pub fn entropy_change_to(&self, other: &Self) -> Result<f64> {
        self.same_source(other)?;
        let base = self
            .source
            .caloric
            .specific_entropy_change(self.temperature_k(), other.temperature_k())?;
        finite(self.source.mass_kg.mul_add(
            base,
            self.source.volume_m3
                * (other.potential.entropy_j_m3_k() - self.potential.entropy_j_m3_k()),
        ))
    }

    /// Solve signed heat while allowing native transverse thermal deformation.
    ///
    /// Reuses the full-F owner's bounded continuation, cached source states,
    /// safeguarded inverse and energy/entropy closure. A local positive-capacity
    /// check is not a global stability certificate. Zero heat preserves the exact
    /// state; source, material, budget or cancellation refusal publishes nothing.
    pub fn add_heat(
        &self,
        heat_j: f64,
        controls: FiniteThermoelasticHeatSolve,
        cancelled: impl FnMut() -> bool,
    ) -> Result<IncompressiblePlaneHeatBalance> {
        let result = kernel::solve(self, controls, |_| finite(heat_j), 0.0, cancelled)?;
        Ok(IncompressiblePlaneHeatBalance { result })
    }
}

impl kernel::Endpoint for IncompressiblePlaneCaloricState {
    fn temperature_k(&self) -> f64 {
        self.temperature_k()
    }
    fn temperature_range_k(&self) -> [f64; 2] {
        self.source.coefficients.temperature_range_k()
    }
    fn heat_capacity_j_per_k(&self) -> f64 {
        self.heat_capacity_j_per_k()
    }
    fn at_sample_temperature(&self, temperature: f64) -> Result<Self> {
        evaluate(&self.source, temperature)
    }
    fn internal_energy_change_to(&self, other: &Self) -> Result<f64> {
        self.internal_energy_change_to(other)
    }
    fn entropy_change_to(&self, other: &Self) -> Result<f64> {
        self.entropy_change_to(other)
    }
    fn mechanical_entropy_scale(&self, other: &Self) -> f64 {
        self.source.volume_m3
            * (self.potential.entropy_j_m3_k().abs() + other.potential.entropy_j_m3_k().abs())
    }
}

impl IncompressiblePlaneHeatBalance {
    pub const fn final_state(&self) -> &IncompressiblePlaneCaloricState {
        &self.result.final_state
    }
    pub const fn heat_j(&self) -> f64 {
        self.result.heat_j
    }
    pub const fn internal_energy_change_j(&self) -> f64 {
        self.result.energy_change_j
    }
    pub const fn entropy_change_j_per_k(&self) -> f64 {
        self.result.entropy_change_j_per_k
    }
    pub fn energy_residual_j(&self) -> f64 {
        self.result.energy_change_j - self.result.heat_j
    }
    pub const fn temperature_bracket_k(&self) -> [f64; 2] {
        self.result.temperature_bracket_k
    }
    pub const fn iterations(&self) -> usize {
        self.result.iterations
    }
    pub const fn state_evaluations(&self) -> usize {
        self.result.state_evaluations
    }
}
