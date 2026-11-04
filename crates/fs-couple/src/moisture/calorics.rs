//! Actual isothermal heat duty of native finite-atmosphere moisture transport.
//!
//! The ORIGINAL diffusion and finite gas owners advance first in a private
//! candidate. A pinned differential isosteric heat q_st(M,T), on the same solid
//! card and fixed context, then accounts for every cell's bound-water storage.
//! For an ideal vapor and negligible bound-phase pV,
//! `partial_u_bound = h_vapor - q_st`. Thus the thermostat heat entering the
//! rigid closed gas+solid system is integral `(R_water*T - q_st) dm_water`.
//! Native gas energy and independent bound-water increments audit that duty.
//!
//! Temperature remains prescribed, with one real external thermostat. There is
//! no second diffusion solver, implicit latent-heat default, thermal feedback,
//! or claim that per-cell storage changes are unique local heat deposition.
//! Bound-water energy fluxes and conduction would be needed for that spatial
//! heat equation. Internal redistribution is nevertheless included in the
//! global duty, even when the atmosphere receives no net water.

use super::atmosphere::{
    AtmosphereMoistureCandidate, AtmosphereMoistureError, AtmosphereMoistureSolve,
    DryBasisSorption, FiniteMoistureAtmosphere,
};
use super::binding::MoistureField;
use super::cartesian::{MoistureCartesian, MoistureCartesianReport};
use super::reservoir::ReservoirField;
use super::{
    MAX_MOISTURE_CELLS, MoistureError, MoistureNonlinearConfig, MoistureSlab, MoistureStepConfig,
    MoistureStepReport, poll, positive,
};
use fs_blake3::ContentHash;
use fs_exec::Cx;
use fs_material::state_point::MoistureSorptionEnergy;
use fs_math::dd::Dd;
use fs_qty::Temperature;
use fs_thermochem::{
    Nasa9EvaluationReceiptV1, ThermochemErrorV1, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R,
};

/// Additional source-integration work and independent closed-energy admission.
#[derive(Clone, Copy, Debug)]
pub struct SorptionEnergyControls {
    /// Sum of actual crossed heat-source segments over all physical cells.
    /// Zero admits unchanged states only. At most 4096 segments per cell.
    pub maximum_source_segments: usize,
    /// Absolute allowed `delta_U_gas + delta_U_bound - Q_thermostat` [J].
    /// This cannot be scaled by arbitrary absolute caloric reference energies.
    pub energy_balance_tolerance_j: f64,
}

#[derive(Debug)]
pub enum SorptionEnergyError {
    Atmosphere(AtmosphereMoistureError),
    Moisture(MoistureError),
    Thermochemistry(ThermochemErrorV1),
    Invalid(&'static str),
    SourceBudget { required: usize, maximum: usize },
    EnergyBalance { residual_j: f64, tolerance_j: f64 },
}
impl core::fmt::Display for SorptionEnergyError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "moisture sorption energy: {self:?}")
    }
}
impl std::error::Error for SorptionEnergyError {}
impl From<AtmosphereMoistureError> for SorptionEnergyError {
    fn from(value: AtmosphereMoistureError) -> Self {
        Self::Atmosphere(value)
    }
}
impl From<MoistureError> for SorptionEnergyError {
    fn from(value: MoistureError) -> Self {
        Self::Moisture(value)
    }
}
impl From<fs_material::state_point::MaterialStatePointError> for SorptionEnergyError {
    fn from(value: fs_material::state_point::MaterialStatePointError) -> Self {
        Self::Moisture(MoistureError::Material(value))
    }
}
impl From<ThermochemErrorV1> for SorptionEnergyError {
    fn from(value: ThermochemErrorV1) -> Self {
        Self::Thermochemistry(value)
    }
}
type Result<T> = core::result::Result<T, SorptionEnergyError>;
fn invalid(reason: &'static str) -> SorptionEnergyError {
    SorptionEnergyError::Invalid(reason)
}

/// Measured caloric consequences of one native conservative transport interval.
#[derive(Clone, Debug)]
pub struct SorptionEnergyStep {
    pub interval_s: [f64; 2],
    pub temperature_k: f64,
    pub sorption_energy_source: ContentHash,
    /// Original vapor species, caloric coefficients and reference convention.
    pub vapor_calorics: Nasa9EvaluationReceiptV1,
    /// Local changes of bound-water internal energy on that vapor datum [J].
    /// These are storage changes, NOT unique per-cell heat sources.
    pub cell_bound_internal_energy_change_j: Vec<f64>,
    /// Local bound storage relative to the original vapor u datum [J]. This
    /// preserves loading-dependent redistribution and avoids large datum loss.
    pub cell_binding_energy_change_j: Vec<f64>,
    pub bound_internal_energy_change_j: f64,
    pub gas_internal_energy_change_j: f64,
    /// External thermostat heat entering the closed rigid gas+solid system [J].
    /// Negative means heat must be removed to maintain the prescribed T.
    pub thermostat_heat_j: f64,
    /// Actual mean external heat duty over this interval; not its time waveform.
    pub mean_thermostat_power_w: f64,
    pub energy_balance_residual_j: f64,
    pub source_segments: usize,
}

/// Transport and energy accepted together; inputs remain untouched on refusal.
#[derive(Clone, Debug)]
pub struct SorptionEnergyCandidate<
    Field = MoistureSlab,
    Report = MoistureStepReport,
    const FACES: usize = 2,
> {
    native: AtmosphereMoistureCandidate<Field, Report, FACES>,
    energy: SorptionEnergyStep,
}
impl<Field, Report, const FACES: usize> SorptionEnergyCandidate<Field, Report, FACES> {
    pub const fn native(&self) -> &AtmosphereMoistureCandidate<Field, Report, FACES> {
        &self.native
    }
    pub const fn energy(&self) -> &SorptionEnergyStep {
        &self.energy
    }
    pub fn into_parts(
        self,
    ) -> (
        AtmosphereMoistureCandidate<Field, Report, FACES>,
        SorptionEnergyStep,
    ) {
        (self.native, self.energy)
    }
}
pub type CartesianSorptionEnergyCandidate =
    SorptionEnergyCandidate<MoistureCartesian, MoistureCartesianReport, 6>;

fn admit(
    field: MoistureField<'_>,
    law: &MoistureSorptionEnergy,
    controls: SorptionEnergyControls,
    cx: &Cx<'_>,
) -> Result<()> {
    poll(cx)?;
    if !positive(controls.energy_balance_tolerance_j)
        || controls.maximum_source_segments > MAX_MOISTURE_CELLS * fs_matdb::MAX_PROPERTY_LINE_KNOTS
    {
        return Err(invalid(
            "finite energy tolerance and bounded source-segment work are required",
        ));
    }
    let context = law
        .source_point()
        .clone()
        .with(
            "moisture_content",
            field.reference.query().axes()["moisture_content"],
        )
        .map_err(|_| invalid("invalid original moisture source context"))?;
    if law.card_identity() != field.card.content_hash() || context != field.reference.query() {
        return Err(invalid(
            "sorption heat must retain the same solid card and all fixed transport conditions",
        ));
    }
    for &m in field.moisture {
        poll(cx)?;
        law.integration_segments(m, m)?;
    }
    Ok(())
}

fn finish<Field: ReservoirField<N>, const N: usize>(
    atmosphere: &FiniteMoistureAtmosphere,
    before: &Field,
    native: AtmosphereMoistureCandidate<Field, Field::Report, N>,
    law: &MoistureSorptionEnergy,
    controls: SorptionEnergyControls,
    cx: &Cx<'_>,
) -> Result<SorptionEnergyCandidate<Field, Field::Report, N>> {
    let old = before.inventory();
    let next = native.moisture().inventory();
    if old.identity != next.identity
        || old.moisture.len() != next.moisture.len()
        || old.dry_mass_per_cell_kg.to_bits() != next.dry_mass_per_cell_kg.to_bits()
    {
        return Err(invalid(
            "native transport changed the original dry skeleton inventory",
        ));
    }
    let catalog = atmosphere.source_catalog();
    let water_id = catalog
        .water_species()
        .ok_or_else(|| invalid("original catalog must explicitly identify water vapor"))?;
    let species = catalog
        .thermodynamics()
        .components()
        .iter()
        .find(|species| species.species() == water_id)
        .ok_or_else(|| invalid("original vapor caloric source is absent"))?;
    let temperature_k = atmosphere.gas().state().state().temperature;
    let vapor = species.evaluate(Temperature::new(temperature_k))?;
    let h = vapor.properties().h().value() / species.molar_mass().value();
    let rt = R * temperature_k / species.molar_mass().value();
    if !h.is_finite() || !positive(rt) {
        return Err(invalid("original vapor calorics are unrepresentable"));
    }
    let mut cell_bound = Vec::with_capacity(old.moisture.len());
    let mut cell_binding = Vec::with_capacity(old.moisture.len());
    let mut bound_sum = Dd::ZERO;
    let mut binding_sum = Dd::ZERO;
    let mut segments = 0usize;
    for (&a, &b) in old.moisture.iter().zip(next.moisture) {
        poll(cx)?;
        segments = segments
            .checked_add(law.integration_segments(a, b)?)
            .ok_or_else(|| invalid("source-segment work overflow"))?;
        if segments > controls.maximum_source_segments {
            return Err(SorptionEnergyError::SourceBudget {
                required: segments,
                maximum: controls.maximum_source_segments,
            });
        }
        let mass = Dd::from_f64(old.dry_mass_per_cell_kg);
        let dm = mass * Dd::from_f64(b - a);
        let release = mass * Dd::from_f64(law.integral(a, b)?);
        let bound = Dd::from_f64(h) * dm - release;
        let binding = Dd::from_f64(rt) * dm - release;
        if !bound.to_f64().is_finite() || !binding.to_f64().is_finite() {
            return Err(invalid(
                "local sorption storage increment is unrepresentable",
            ));
        }
        cell_bound.push(bound.to_f64());
        cell_binding.push(binding.to_f64());
        bound_sum = bound_sum + bound;
        binding_sum = binding_sum + binding;
    }
    let gas_change = native.step().gas_internal_energy_change_j;
    let heat = binding_sum.to_f64();
    let residual = (Dd::from_f64(gas_change) + bound_sum - binding_sum).to_f64();
    let power = heat / (next.time_s - old.time_s);
    if !heat.is_finite() || !power.is_finite() || !bound_sum.to_f64().is_finite() {
        return Err(invalid(
            "integrated sorption storage or thermostat power is unrepresentable",
        ));
    }
    if !residual.is_finite() || residual.abs() > controls.energy_balance_tolerance_j {
        return Err(SorptionEnergyError::EnergyBalance {
            residual_j: residual,
            tolerance_j: controls.energy_balance_tolerance_j,
        });
    }
    let energy = SorptionEnergyStep {
        interval_s: [old.time_s, next.time_s],
        temperature_k,
        sorption_energy_source: law.identity(),
        vapor_calorics: vapor.receipt().clone(),
        cell_bound_internal_energy_change_j: cell_bound,
        cell_binding_energy_change_j: cell_binding,
        bound_internal_energy_change_j: bound_sum.to_f64(),
        gas_internal_energy_change_j: gas_change,
        thermostat_heat_j: heat,
        mean_thermostat_power_w: power,
        energy_balance_residual_j: residual,
        source_segments: segments,
    };
    poll(cx)?;
    Ok(SorptionEnergyCandidate { native, energy })
}

impl FiniteMoistureAtmosphere {
    /// Native finite-vapor/slab transport plus actual closed-system thermostat
    /// duty from the solid's original differential isosteric heat. Source,
    /// work, energy or cancellation failure returns neither endpoint.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_sorption_energy(
        &self,
        slab: &MoistureSlab,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 2],
        transport: MoistureStepConfig,
        atmosphere: AtmosphereMoistureSolve,
        law: &MoistureSorptionEnergy,
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<SorptionEnergyCandidate> {
        admit(slab.field(), law, controls, cx)?;
        let native = self.advance(slab, dt_s, faces, transport, atmosphere, cx)?;
        finish(self, slab, native, law, controls, cx)
    }

    /// Same caloric closure over the actual three-dimensional cell inventory
    /// and native six-face endpoint; no surrogate slab or lumped moisture field.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_cartesian_with_sorption_energy(
        &self,
        field: &MoistureCartesian,
        dt_s: f64,
        faces: &[Option<DryBasisSorption>; 6],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        atmosphere: AtmosphereMoistureSolve,
        law: &MoistureSorptionEnergy,
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<CartesianSorptionEnergyCandidate> {
        admit(field.field(), law, controls, cx)?;
        let native =
            self.advance_cartesian(field, dt_s, faces, transport, nonlinear, atmosphere, cx)?;
        finish(self, field, native, law, controls, cx)
    }
}
