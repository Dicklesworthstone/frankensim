//! One finite atmosphere and all specimen binding-energy stores, counted once.

use super::*;
use crate::moisture::atmosphere::{
    CartesianMoistureParticipant, MoistureParticipant, SharedAtmosphereMoistureCandidate,
    admit_participants,
};
use crate::moisture::cartesian::MoistureCartesianSdirk2Report;

/// Bind one original isosteric-heat source to a stable physical specimen ID.
/// Enumeration order is irrelevant; exactly one source per participant is required.
#[derive(Clone, Copy, Debug)]
pub struct SorptionEnergyBinding<'a> {
    pub id: u64,
    pub law: &'a MoistureSorptionEnergy,
}

/// The same native per-material storage account used for layered bodies.
#[derive(Clone, Debug)]
pub struct SorptionSpecimenEnergy {
    pub id: u64,
    pub storage: SorptionLayerEnergy,
}

/// Complete rigid, isothermal gas-plus-specimens energy balance.
#[derive(Clone, Debug)]
pub struct SharedSorptionEnergyStep {
    pub interval_s: [f64; 2],
    pub temperature_k: f64,
    pub vapor_calorics: Nasa9EvaluationReceiptV1,
    /// Stable-ID order, matching the original shared transport candidate.
    pub specimens: Vec<SorptionSpecimenEnergy>,
    pub bound_internal_energy_change_j: f64,
    /// One actual gas increment, irrespective of the number of specimens.
    pub gas_internal_energy_change_j: f64,
    /// External heat entering the complete system to hold its temperature.
    pub thermostat_heat_j: f64,
    pub mean_thermostat_power_w: f64,
    pub energy_balance_residual_j: f64,
    /// One aggregate allowance across every specimen and physical cell.
    pub source_segments: usize,
}

/// Native transport and complete calorics are accepted or refused together.
#[derive(Clone, Debug)]
pub struct SharedSorptionEnergyCandidate<Field = MoistureSlab, Report = MoistureStepReport> {
    native: SharedAtmosphereMoistureCandidate<Field, Report>,
    energy: SharedSorptionEnergyStep,
}
impl<Field, Report> SharedSorptionEnergyCandidate<Field, Report> {
    pub const fn native(&self) -> &SharedAtmosphereMoistureCandidate<Field, Report> {
        &self.native
    }
    pub const fn energy(&self) -> &SharedSorptionEnergyStep {
        &self.energy
    }
    pub fn into_parts(
        self,
    ) -> (
        SharedAtmosphereMoistureCandidate<Field, Report>,
        SharedSorptionEnergyStep,
    ) {
        (self.native, self.energy)
    }
}

pub type SharedCartesianSorptionEnergyCandidate =
    SharedSorptionEnergyCandidate<MoistureCartesian, MoistureCartesianReport>;
pub type SharedCartesianSdirk2SorptionEnergyCandidate =
    SharedSorptionEnergyCandidate<MoistureCartesian, MoistureCartesianSdirk2Report>;
pub type CartesianSdirk2SorptionEnergyCandidate =
    SorptionEnergyCandidate<MoistureCartesian, MoistureCartesianSdirk2Report, 6>;

type EnergyBody<'a> = (u64, MoistureField<'a>, &'a MoistureSorptionEnergy);

fn bind_shared<'a>(
    fields: impl ExactSizeIterator<Item = (u64, MoistureField<'a>)>,
    bindings: &[SorptionEnergyBinding<'a>],
    maximum_total_cells: usize,
    controls: SorptionEnergyControls,
    cx: &Cx<'_>,
) -> Result<Vec<EnergyBody<'a>>> {
    poll(cx)?;
    admit_participants(fields.len(), maximum_total_cells)?;
    if fields.len() != bindings.len() {
        return Err(invalid(
            "one sorption heat source is required per shared specimen",
        ));
    }
    let mut fields: Vec<_> = fields.collect();
    let mut bindings: Vec<_> = bindings.iter().collect();
    fields.sort_unstable_by_key(|(id, _)| *id);
    bindings.sort_unstable_by_key(|binding| binding.id);
    poll(cx)?;
    if fields.windows(2).any(|pair| pair[0].0 == pair[1].0)
        || bindings.windows(2).any(|pair| pair[0].id == pair[1].id)
    {
        return Err(invalid("shared caloric specimen IDs must be unique"));
    }
    let mut total_cells = 0_usize;
    let mut ordered = Vec::with_capacity(fields.len());
    for ((id, field), binding) in fields.into_iter().zip(bindings) {
        poll(cx)?;
        if id != binding.id {
            return Err(invalid(
                "sorption heat source IDs must match the original specimens",
            ));
        }
        total_cells = total_cells
            .checked_add(field.moisture.len())
            .ok_or_else(|| invalid("shared caloric cell budget overflow"))?;
        if total_cells > maximum_total_cells {
            return Err(invalid("shared caloric cell budget exceeded"));
        }
        admit(field, binding.law, controls, cx)?;
        ordered.push((id, field, binding.law));
    }
    Ok(ordered)
}

fn finish_shared<Field: ReservoirField<N>, Report, const N: usize>(
    atmosphere: &FiniteMoistureAtmosphere,
    ordered: &[EnergyBody<'_>],
    native: SharedAtmosphereMoistureCandidate<Field, Report>,
    controls: SorptionEnergyControls,
    cx: &Cx<'_>,
) -> Result<SharedSorptionEnergyCandidate<Field, Report>> {
    poll(cx)?;
    if ordered.is_empty() || native.solids().len() != ordered.len() {
        return Err(invalid("native shared transport changed the specimen set"));
    }
    let vapor = vapor_datum(atmosphere)?;
    let interval_s = native.step().interval_s;
    let mut bound = Dd::ZERO;
    let mut binding = Dd::ZERO;
    let mut segments = 0;
    let mut specimens = Vec::with_capacity(ordered.len());
    for (&(id, old, law), solid) in ordered.iter().zip(native.solids()) {
        poll(cx)?;
        let next = solid.moisture.inventory();
        if id != solid.id
            || old.time_s.to_bits() != interval_s[0].to_bits()
            || next.time_s.to_bits() != interval_s[1].to_bits()
            || old.temperature_kelvin().to_bits() != vapor.temperature_k.to_bits()
            || next.temperature_kelvin().to_bits() != vapor.temperature_k.to_bits()
        {
            return Err(invalid(
                "shared calorics requires original IDs, clocks and temperature",
            ));
        }
        let first_segment = segments;
        let storage = integrate_storage(old, next, law, &vapor, &mut segments, controls, cx)?;
        if !storage.water.to_f64().is_finite()
            || !storage.bound.to_f64().is_finite()
            || !storage.binding.to_f64().is_finite()
        {
            return Err(invalid("shared specimen storage is unrepresentable"));
        }
        bound = bound + storage.bound;
        binding = binding + storage.binding;
        specimens.push(SorptionSpecimenEnergy {
            id,
            storage: SorptionLayerEnergy {
                sorption_energy_source: law.identity(),
                water_change_kg: storage.water.to_f64(),
                cell_bound_internal_energy_change_j: storage.cell_bound,
                cell_binding_energy_change_j: storage.cell_binding,
                bound_internal_energy_change_j: storage.bound.to_f64(),
                binding_energy_change_j: storage.binding.to_f64(),
                source_segments: segments - first_segment,
            },
        });
    }
    let gas_change = native.step().gas_internal_energy_change_j;
    let (heat, power, residual) = close_energy(bound, binding, gas_change, interval_s, controls)?;
    let energy = SharedSorptionEnergyStep {
        interval_s,
        temperature_k: vapor.temperature_k,
        vapor_calorics: vapor.receipt,
        specimens,
        bound_internal_energy_change_j: bound.to_f64(),
        gas_internal_energy_change_j: gas_change,
        thermostat_heat_j: heat,
        mean_thermostat_power_w: power,
        energy_balance_residual_j: residual,
        source_segments: segments,
    };
    poll(cx)?;
    Ok(SharedSorptionEnergyCandidate { native, energy })
}

impl FiniteMoistureAtmosphere {
    /// Advance original slab participants against one finite vapor inventory,
    /// then account every cell's actual binding-energy change on one vapor datum.
    ///
    /// The original native solve owns the common humidity and water balance.
    /// Calorics uses the complete endpoint, including opposing specimen transfers
    /// and sealed internal redistribution. It neither advances transport again nor
    /// infers a local heat source from a cell's storage change. The gas energy
    /// enters once. All fields and reports are discarded on any caloric refusal.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_shared_with_sorption_energy(
        &self,
        participants: &[MoistureParticipant<'_>],
        bindings: &[SorptionEnergyBinding<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        atmosphere: AtmosphereMoistureSolve,
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<SharedSorptionEnergyCandidate> {
        let ordered = bind_shared(
            participants.iter().map(|body| (body.id, body.slab.field())),
            bindings,
            maximum_total_cells,
            controls,
            cx,
        )?;
        let native =
            self.advance_shared(participants, dt_s, maximum_total_cells, atmosphere, cx)?;
        finish_shared::<_, _, 2>(self, &ordered, native, controls, cx)
    }

    /// Six-face, anisotropic Cartesian transport for all original specimens
    /// with one complete gas-plus-solids heat account and source-work budget.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_shared_cartesian_with_sorption_energy(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        bindings: &[SorptionEnergyBinding<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        atmosphere: AtmosphereMoistureSolve,
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<SharedCartesianSorptionEnergyCandidate> {
        let ordered = bind_shared(
            participants
                .iter()
                .map(|body| (body.id, body.moisture.field())),
            bindings,
            maximum_total_cells,
            controls,
            cx,
        )?;
        let native =
            self.advance_shared_cartesian(participants, dt_s, maximum_total_cells, atmosphere, cx)?;
        finish_shared::<_, _, 6>(self, &ordered, native, controls, cx)
    }

    /// Original coupled SDIRK stages with one original-to-final heat account.
    ///
    /// Each stage retains its own finite humidity and material laws. The first
    /// stage is a transport diagnostic, not another water or energy debit.
    /// Isosteric heat integrates each actual full-interval cell increment;
    /// no nonlinear heat is approximated from only the net vapor transfer.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_shared_cartesian_sdirk2_with_sorption_energy(
        &self,
        participants: &[CartesianMoistureParticipant<'_>],
        bindings: &[SorptionEnergyBinding<'_>],
        dt_s: f64,
        maximum_total_cells: usize,
        atmosphere: AtmosphereMoistureSolve,
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<SharedCartesianSdirk2SorptionEnergyCandidate> {
        let ordered = bind_shared(
            participants
                .iter()
                .map(|body| (body.id, body.moisture.field())),
            bindings,
            maximum_total_cells,
            controls,
            cx,
        )?;
        let native = self.advance_shared_cartesian_sdirk2(
            participants,
            dt_s,
            maximum_total_cells,
            atmosphere,
            cx,
        )?;
        finish_shared::<_, _, 6>(self, &ordered, native, controls, cx)
    }

    /// One Cartesian SDIRK field with the same complete endpoint calorics.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_cartesian_sdirk2_with_sorption_energy(
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
    ) -> Result<CartesianSdirk2SorptionEnergyCandidate> {
        admit(field.field(), law, controls, cx)?;
        let native = self
            .advance_cartesian_sdirk2(field, dt_s, faces, transport, nonlinear, atmosphere, cx)?;
        finish(self, field, native, law, controls, cx)
    }
}
