//! One closed gas/solid heat account across original unlike material layers.
//!
//! Each layer retains its original differential isosteric heat and actual cell
//! dry masses. An internal water transfer may change binding energy even when
//! the whole body and atmosphere exchange no net water. The original gas energy
//! enters once, on one shared vapor datum; material interfaces create no extra
//! gas inventory. Temperature remains prescribed by an external thermostat.

use super::*;
use crate::moisture::{
    layered::{MoistureLayeredReport, MoistureLayeredSdirk2Report, MoistureLayeredSlab},
    radial::layered::MoistureLayeredRadial,
};

mod sdirk;
pub use sdirk::{LayeredRadialSdirk2SorptionEnergyCandidate, LayeredSdirk2SorptionEnergyCandidate};

/// Actual storage changes for one original material, in physical layer order.
/// The binding increment is relative to the common original vapor-u datum.
/// It is not a uniquely located heat source or an additional gas energy.
#[derive(Clone, Debug)]
pub struct SorptionLayerEnergy {
    pub sorption_energy_source: ContentHash,
    pub water_change_kg: f64,
    pub cell_bound_internal_energy_change_j: Vec<f64>,
    pub cell_binding_energy_change_j: Vec<f64>,
    pub bound_internal_energy_change_j: f64,
    pub binding_energy_change_j: f64,
    pub source_segments: usize,
}

/// Independent caloric account of the complete heterogeneous body and one gas.
#[derive(Clone, Debug)]
pub struct LayeredSorptionEnergyStep {
    pub interval_s: [f64; 2],
    pub temperature_k: f64,
    pub vapor_calorics: Nasa9EvaluationReceiptV1,
    pub layers: Vec<SorptionLayerEnergy>,
    pub bound_internal_energy_change_j: f64,
    /// Native shared gas increment, counted once across all original layers.
    pub gas_internal_energy_change_j: f64,
    /// External heat entering the rigid closed gas+body to hold temperature.
    pub thermostat_heat_j: f64,
    /// Mean duty over this interval, not a resolved heat-flow waveform.
    pub mean_thermostat_power_w: f64,
    pub energy_balance_residual_j: f64,
    /// Aggregate actual source integration work, not a per-layer allowance.
    pub source_segments: usize,
}

/// Native transport and the full isothermal energy account publish together.
#[derive(Clone, Debug)]
pub struct LayeredSorptionEnergyCandidate<
    Field = MoistureLayeredSlab,
    Report = MoistureLayeredReport,
> {
    native: AtmosphereMoistureCandidate<Field, Report>,
    energy: LayeredSorptionEnergyStep,
}
impl<Field, Report> LayeredSorptionEnergyCandidate<Field, Report> {
    pub const fn native(&self) -> &AtmosphereMoistureCandidate<Field, Report> {
        &self.native
    }
    pub const fn energy(&self) -> &LayeredSorptionEnergyStep {
        &self.energy
    }
    pub fn into_parts(
        self,
    ) -> (
        AtmosphereMoistureCandidate<Field, Report>,
        LayeredSorptionEnergyStep,
    ) {
        (self.native, self.energy)
    }
}

/// Original spherical or cylindrical coatings, finite vapor and caloric duty.
pub type LayeredRadialSorptionEnergyCandidate =
    LayeredSorptionEnergyCandidate<MoistureLayeredRadial>;

fn admit_layers<'a>(
    fields: impl ExactSizeIterator<Item = MoistureField<'a>>,
    laws: &[MoistureSorptionEnergy],
    controls: SorptionEnergyControls,
    cx: &Cx<'_>,
) -> Result<()> {
    poll(cx)?;
    if fields.len() == 0 || fields.len() != laws.len() {
        return Err(invalid(
            "one original sorption heat source is required per physical layer",
        ));
    }
    for (field, law) in fields.zip(laws) {
        admit(field, law, controls, cx)?;
    }
    Ok(())
}

fn account_layers<'a, 'b>(
    atmosphere: &FiniteMoistureAtmosphere,
    before: impl ExactSizeIterator<Item = MoistureField<'a>>,
    after: impl ExactSizeIterator<Item = MoistureField<'b>>,
    laws: &[MoistureSorptionEnergy],
    gas_change: f64,
    controls: SorptionEnergyControls,
    cx: &Cx<'_>,
) -> Result<LayeredSorptionEnergyStep> {
    if before.len() == 0 || before.len() != laws.len() || after.len() != laws.len() {
        return Err(invalid(
            "native transport changed the original layer layout",
        ));
    }
    let vapor = vapor_datum(atmosphere)?;
    let mut layers = Vec::with_capacity(laws.len());
    let mut bound = Dd::ZERO;
    let mut binding = Dd::ZERO;
    let mut segments = 0;
    let mut interval = None;
    for ((old, next), law) in before.zip(after).zip(laws) {
        poll(cx)?;
        let times = [old.time_s, next.time_s];
        if interval.is_some_and(|span: [f64; 2]| {
            span[0].to_bits() != times[0].to_bits() || span[1].to_bits() != times[1].to_bits()
        }) || old.temperature_kelvin().to_bits() != vapor.temperature_k.to_bits()
            || next.temperature_kelvin().to_bits() != vapor.temperature_k.to_bits()
        {
            return Err(invalid(
                "all material heat accounts require common clocks and temperature",
            ));
        }
        interval = Some(times);
        let first_segment = segments;
        let storage = integrate_storage(old, next, law, &vapor, &mut segments, controls, cx)?;
        if !storage.water.to_f64().is_finite()
            || !storage.bound.to_f64().is_finite()
            || !storage.binding.to_f64().is_finite()
        {
            return Err(invalid("integrated layer storage is unrepresentable"));
        }
        bound = bound + storage.bound;
        binding = binding + storage.binding;
        layers.push(SorptionLayerEnergy {
            sorption_energy_source: law.identity(),
            water_change_kg: storage.water.to_f64(),
            cell_bound_internal_energy_change_j: storage.cell_bound,
            cell_binding_energy_change_j: storage.cell_binding,
            bound_internal_energy_change_j: storage.bound.to_f64(),
            binding_energy_change_j: storage.binding.to_f64(),
            source_segments: segments - first_segment,
        });
    }
    let interval_s = interval.expect("admitted nonempty layers");
    let (heat, power, residual) = close_energy(bound, binding, gas_change, interval_s, controls)?;
    let energy = LayeredSorptionEnergyStep {
        interval_s,
        temperature_k: vapor.temperature_k,
        vapor_calorics: vapor.receipt,
        layers,
        bound_internal_energy_change_j: bound.to_f64(),
        gas_internal_energy_change_j: gas_change,
        thermostat_heat_j: heat,
        mean_thermostat_power_w: power,
        energy_balance_residual_j: residual,
        source_segments: segments,
    };
    poll(cx)?;
    Ok(energy)
}

impl FiniteMoistureAtmosphere {
    /// Advance every original slab layer against this finite gas and account
    /// for the heat needed to hold their common temperature.
    ///
    /// Each heat law must bind its corresponding layer's original card and
    /// complete current transport context. The native nonlinear sorption,
    /// diffusivity and interface/surface resistance determine every candidate.
    /// Source integration then uses actual local moisture increments and dry
    /// masses, including internal transfers invisible to net gas uptake.
    ///
    /// One aggregate source-work cap and one closed-system energy tolerance
    /// apply to the whole body. Native gas energy is counted once on its own
    /// vapor caloric datum. Source, work, energy or cancellation refusal returns
    /// neither endpoint. This computes isothermal thermostat duty; it does not
    /// supply a spatial heat equation, sorption-driven thermal feedback or a
    /// unique per-layer heat deposition.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_with_sorption_energy(
        &self,
        body: &MoistureLayeredSlab,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        atmosphere: AtmosphereMoistureSolve,
        laws: &[MoistureSorptionEnergy],
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<LayeredSorptionEnergyCandidate> {
        admit_layers(
            body.layers().iter().map(|layer| layer.slab().field()),
            laws,
            controls,
            cx,
        )?;
        let native =
            self.advance_layered(body, dt_s, exposed, transport, nonlinear, atmosphere, cx)?;
        let energy = account_layers(
            self,
            body.layers().iter().map(|layer| layer.slab().field()),
            native
                .moisture()
                .layers()
                .iter()
                .map(|layer| layer.slab().field()),
            laws,
            native.step().gas_internal_energy_change_j,
            controls,
            cx,
        )?;
        Ok(LayeredSorptionEnergyCandidate { native, energy })
    }

    /// The same native caloric closure for a coated cylinder or sphere.
    /// Each material retains its original source and every shell's actual dry
    /// mass. Curved paths and internal partition stay in the native diffusion
    /// owner; the shared gas energy is not duplicated at material interfaces.
    ///
    /// Admission, aggregate source work, final energy closure and the held
    /// temperature scope are those of [Self::advance_layered_with_sorption_energy].
    #[allow(clippy::too_many_arguments)]
    pub fn advance_layered_radial_with_sorption_energy(
        &self,
        body: &MoistureLayeredRadial,
        dt_s: f64,
        exposed: [bool; 2],
        transport: MoistureStepConfig,
        nonlinear: MoistureNonlinearConfig,
        atmosphere: AtmosphereMoistureSolve,
        laws: &[MoistureSorptionEnergy],
        controls: SorptionEnergyControls,
        cx: &Cx<'_>,
    ) -> Result<LayeredRadialSorptionEnergyCandidate> {
        admit_layers(
            body.layers().iter().map(|layer| layer.radial().field()),
            laws,
            controls,
            cx,
        )?;
        let native =
            self.advance_layered_radial(body, dt_s, exposed, transport, nonlinear, atmosphere, cx)?;
        let energy = account_layers(
            self,
            body.layers().iter().map(|layer| layer.radial().field()),
            native
                .moisture()
                .layers()
                .iter()
                .map(|layer| layer.radial().field()),
            laws,
            native.step().gas_internal_energy_change_j,
            controls,
            cx,
        )?;
        Ok(LayeredSorptionEnergyCandidate { native, energy })
    }
}
