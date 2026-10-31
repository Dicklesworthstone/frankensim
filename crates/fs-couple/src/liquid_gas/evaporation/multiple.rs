//! Joint phase transfer from many native interfaces into one finite gas.
//!
//! Every rate uses the same original gas. The material chamber owns one shared
//! donor account and one gas caloric solve; native VOF owns the represented
//! masses, one material refresh and one subsequent flow step. No per-patch
//! intermediate gas or duplicate liquid inventory is published.

use std::collections::BTreeMap;

use fs_lbm::material::free_surface::enthalpy::{
    FreeSurfaceMassExchangeSample, FreeSurfaceMassProjection,
};
use fs_material::evaporation::gas::chamber::heat::{
    GasLiquidChamberHeatSolve, GasLiquidChamberHeatStep, GasLiquidChamberHeatSurface,
    GasLiquidChamberSurfaceHeatStep, advance_chamber_with_gas_and_heat_mass_selection,
};
use fs_material::evaporation::gas::chamber::{
    GasLiquidChamberSolve, GasLiquidChamberStep, GasLiquidChamberSurface,
    GasLiquidChamberSurfaceStep, advance_chamber_with_gas_mass_selection,
};

use super::{
    CancelGate, ClosedGasParcel, EvaporatingLiquidGasError, EvaporationError,
    FreeSurfaceEnthalpyReport, FreeSurfaceMassExchange, FreeSurfaceMassExchangeControls,
    FreeSurfaceMassExchangeReport, GasEvaporationError, GasFilmHeatSolve,
    LiquidGasEvaporationSurface, LiquidSurfaceState, MaterialFreeSurfaceEnthalpy3, RefCell,
    ThermochemicalGasModel, poll,
};

type Error = EvaporatingLiquidGasError;

/// One native patch with an explicit additional sensible heat conductance.
/// The coefficient excludes latent and vapor-carried sensible energy.
#[derive(Clone, Copy)]
pub struct LiquidGasHeatSurface<'a> {
    pub surface: LiquidGasEvaporationSurface<'a>,
    pub heat_transfer_coefficient_w_m2_k: f64,
    pub heat_controls: GasFilmHeatSolve,
}

/// One joint source budget and one native transaction/flow budget.
#[derive(Clone, Copy, Debug)]
pub struct DistributedLiquidGasSolve {
    pub chamber: GasLiquidChamberSolve,
    pub native: FreeSurfaceMassExchangeControls,
    /// Whole interval, including absolute local defects and native flow [kg].
    pub mass_tolerance_kg: f64,
    /// Whole sum of absolute source, native storage and flow defects [J].
    pub energy_tolerance_j: f64,
}

/// Original source and actual native storage for one unique liquid cell.
#[derive(Clone, Debug)]
pub struct NativeLiquidGasSurfaceStep {
    pub coordinates: [usize; 3],
    /// Phase-only receipt, starting after optional sensible heating.
    pub phase: GasLiquidChamberSurfaceStep,
    /// Original sensible heat receipt, including the exact phase-initial h.
    pub heat: Option<GasLiquidChamberSurfaceHeatStep>,
    pub native: FreeSurfaceMassExchangeSample,
    /// Native minus original source sensible storage, relative to this cell's
    /// original h, before adding any formation-energy reference [J].
    pub native_caloric_defect_j: f64,
}

/// Joint exchange followed by exactly one native physical step.
#[derive(Clone, Debug)]
pub struct DistributedLiquidGasReport {
    pub duration_s: f64,
    /// Phase-only common-gas account after optional joint sensible heating.
    pub phase: GasLiquidChamberStep,
    /// Combined original-to-final storage, temperatures and defects. Its
    /// energy_into_gas_j is vapor-carried energy only; subtract sensible heat
    /// into liquid to obtain the net gas supply. Equals phase without heat.
    pub exchange: GasLiquidChamberStep,
    /// Original common sensible heat account, already included in exchange.
    pub heat: Option<GasLiquidChamberHeatStep>,
    /// Cartesian coordinate order; source IDs are consecutive in this order.
    pub surfaces: Vec<NativeLiquidGasSurfaceStep>,
    pub mass_exchange: FreeSurfaceMassExchangeReport,
    pub flow: FreeSurfaceEnthalpyReport,
    pub native_liquid_enthalpy_change_j: f64,
    pub native_phase_caloric_defect_j: f64,
    pub mass_residual_kg: f64,
    pub volatile_mass_residual_kg: f64,
    /// Local source/native mass defects, original gas total/species defects,
    /// native transfer defects and subsequent native-flow mass drift.
    pub absolute_mass_defect_kg: f64,
    /// Combined original heat/phase defect plus native caloric and flow defects.
    /// The optional separate heat account is not added a second time.
    pub energy_residual_j: f64,
    /// Absolute local defects cannot cancel between opposing phase streams.
    /// Gas representation rounding remains separately mass-budgeted, matching
    /// the original source-relative gas boundary owner.
    pub absolute_energy_defect_j: f64,
}

pub struct DistributedLiquidGasCandidate<'s, 'm, 'catalog> {
    liquid: MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    report: DistributedLiquidGasReport,
}
impl<'s, 'm, 'catalog> DistributedLiquidGasCandidate<'s, 'm, 'catalog> {
    pub const fn liquid(&self) -> &MaterialFreeSurfaceEnthalpy3<'s, 'm> {
        &self.liquid
    }
    pub const fn gas(&self) -> &ClosedGasParcel {
        &self.gas
    }
    pub const fn catalog(&self) -> &'catalog ThermochemicalGasModel {
        self.catalog
    }
    pub const fn report(&self) -> &DistributedLiquidGasReport {
        &self.report
    }
    pub fn into_parts(
        self,
    ) -> (
        MaterialFreeSurfaceEnthalpy3<'s, 'm>,
        ClosedGasParcel,
        &'catalog ThermochemicalGasModel,
        DistributedLiquidGasReport,
    ) {
        (self.liquid, self.gas, self.catalog, self.report)
    }
}

struct NativeInput {
    coordinates: [usize; 3],
    initial_mass_kg: f64,
    initial_h_j_kg: f64,
}

#[derive(Default)]
struct Total {
    sum: f64,
    correction: f64,
}
impl Total {
    fn add(&mut self, value: f64) -> Result<(), Error> {
        let (sum, error) = fs_math::eft::two_sum(self.sum, value);
        let correction = self.correction + error;
        if ![value, sum, error, correction]
            .into_iter()
            .all(f64::is_finite)
        {
            return Err(Error::Invalid(
                "distributed native account is unrepresentable",
            ));
        }
        self.sum = sum;
        self.correction = correction;
        Ok(())
    }
    fn value(&self) -> f64 {
        self.sum + self.correction
    }
}

/// Exchange all measured native patches with ONE common finite gas, then
/// advance native flow once over its retained timestep.
///
/// Every patch must name a different occupied liquid cell and an exposed
/// Cartesian face. The same native pure liquid requires one volatile identity
/// and one constant gas/liquid datum binding throughout. Multiple faces of a
/// single cell must first be reduced by a geometry/flux owner; treating each
/// face as a second liquid inventory refuses here.
///
/// The original chamber evaluates all frozen rates against the same gas and
/// checks gross vapor withdrawals before any opposing evaporation is credited.
/// Each original native mass is selected before the sole gas caloric solve.
/// One batched native mass/enthalpy transaction refreshes the original material
/// fields; one normal native step transports them. Surface enumeration cannot
/// change the physical result. No grid clone or native step occurs per patch.
///
/// This is a first-order simultaneous phase/flow split; refine the original
/// physical timestep to assess accuracy. Gas volume remains a fixed connected
/// control volume. This entry point adds no sensible gas-film heat, gas traction,
/// sealed-vessel displacement, nucleation, dry-out or interface-overflow remap.
/// The original native owner retains its momentum, electrical-core and solved
/// surface-stress/contact-history admission rules. Any late source refusal,
/// balance failure or cancellation discards both private endpoint owners.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_gas<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[LiquidGasEvaporationSurface<'_>],
    controls: DistributedLiquidGasSolve,
    gate: &CancelGate,
) -> Result<DistributedLiquidGasCandidate<'s, 'm, 'catalog>, Error> {
    advance_native(
        liquid,
        gas,
        catalog,
        surfaces.iter().copied().map(|surface| (surface, None)),
        controls,
        None,
        gate,
    )
}

/// Joint reciprocal sensible heat, simultaneous phase exchange and one native
/// flow step over the SAME retained physical interval.
///
/// The original material owner solves all sensible film exchanges against one
/// common gas endpoint, then evaluates phase transfer on those heated states.
/// Native mass selection occurs only in that final phase stage. The complete
/// original-to-final liquid h change is stored once on each actual native mass;
/// no heat is applied again after the one batched native transaction.
///
/// This first-order heat/phase/flow split needs full-timestep refinement. Its
/// geometry, momentum, source and attached-owner limits match the phase-only
/// entry point. Zero sensible coefficients preserve that entry point's native
/// and gas endpoints while retaining explicit zero heat reports.
#[allow(clippy::too_many_arguments)]
pub fn advance_liquid_surfaces_with_gas_and_heat<'s, 'm, 'catalog>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: &[LiquidGasHeatSurface<'_>],
    controls: DistributedLiquidGasSolve,
    heat_controls: GasLiquidChamberHeatSolve,
    gate: &CancelGate,
) -> Result<DistributedLiquidGasCandidate<'s, 'm, 'catalog>, Error> {
    advance_native(
        liquid,
        gas,
        catalog,
        surfaces.iter().map(|surface| {
            (
                surface.surface,
                Some((
                    surface.heat_transfer_coefficient_w_m2_k,
                    surface.heat_controls,
                )),
            )
        }),
        controls,
        Some(heat_controls),
        gate,
    )
}

#[allow(clippy::too_many_arguments)]
fn advance_native<'s, 'm, 'catalog, 'surface>(
    liquid: &MaterialFreeSurfaceEnthalpy3<'s, 'm>,
    gas: &ClosedGasParcel,
    catalog: &'catalog ThermochemicalGasModel,
    surfaces: impl ExactSizeIterator<
        Item = (
            LiquidGasEvaporationSurface<'surface>,
            Option<(f64, GasFilmHeatSolve)>,
        ),
    >,
    controls: DistributedLiquidGasSolve,
    heat_controls: Option<GasLiquidChamberHeatSolve>,
    gate: &CancelGate,
) -> Result<DistributedLiquidGasCandidate<'s, 'm, 'catalog>, Error> {
    poll(gate)?;
    if [controls.mass_tolerance_kg, controls.energy_tolerance_j]
        .into_iter()
        .any(|x| !x.is_finite() || x < 0.0)
        || surfaces.len() == 0
        || controls.chamber.maximum_surfaces == 0
        || surfaces.len() > controls.chamber.maximum_surfaces
    {
        return Err(Error::Invalid(
            "distributed native exchange needs positive bounded work and finite balance budgets",
        ));
    }
    let mut binding = None;
    let mut ordered = BTreeMap::new();
    for (surface, heat) in surfaces {
        poll(gate)?;
        let common = *binding.get_or_insert(surface.binding);
        if surface.binding != common
            || heat.is_some() != heat_controls.is_some()
            || ordered.insert(surface.liquid, (surface, heat)).is_some()
        {
            return Err(Error::Invalid(
                "native patches need unique liquid cells and one volatile/datum binding",
            ));
        }
    }
    let binding = binding.expect("positive surface count admitted");

    let chart_error = RefCell::new(None);
    let temperature = |h| {
        liquid
            .liquid_temperature_at_specific_enthalpy_k(h)
            .map_err(|error| {
                *chart_error.borrow_mut() = Some(error);
                EvaporationError::Invalid("original native liquid caloric branch refused")
            })
    };
    let mut inputs = Vec::with_capacity(ordered.len());
    let mut sources = Vec::with_capacity(ordered.len());
    let mut source_heat = Vec::with_capacity(if heat_controls.is_some() {
        ordered.len()
    } else {
        0
    });
    for (index, (&coordinates, (surface, heat))) in ordered.iter().enumerate() {
        poll(gate)?;
        liquid.admit_interface_gas_patch(coordinates, surface.gas, surface.area_m2)?;
        let initial = liquid.project_interface_mass_exchange(coordinates, 0.0)?;
        let h0 = liquid
            .specific_enthalpy_at(coordinates)
            .ok_or(Error::Invalid(
                "native phase interface has no retained enthalpy",
            ))?;
        let id = u64::try_from(index)
            .map_err(|_| Error::Invalid("native surface identity is unrepresentable"))?;
        inputs.push(NativeInput {
            coordinates,
            initial_mass_kg: initial.initial_mass_kg,
            initial_h_j_kg: h0,
        });
        sources.push(GasLiquidChamberSurface {
            id,
            law: surface.law,
            binding: surface.binding,
            liquid: LiquidSurfaceState {
                mass_kg: initial.initial_mass_kg,
                specific_enthalpy_j_kg: h0,
                area_m2: surface.area_m2,
            },
            temperature_from_enthalpy: &temperature,
            mass_transfer_coefficient_m_s: surface.mass_transfer_coefficient_m_s,
        });
        if let Some((coefficient, controls)) = heat {
            source_heat.push(GasLiquidChamberHeatSurface {
                surface: *sources.last().expect("inserted original source"),
                heat_transfer_coefficient_w_m2_k: *coefficient,
                heat_controls: *controls,
            });
        }
    }
    let duration_s = liquid.flow().scaling().timestep_s();
    let mut selected: Vec<Option<FreeSurfaceMassProjection>> = vec![None; inputs.len()];
    let mut mass_error = None;
    let select_mass = |id: u64, mass: f64| {
        let index = usize::try_from(id)
            .map_err(|_| GasEvaporationError::Invalid("unknown native surface mass selection"))?;
        if index >= inputs.len() || selected[index].is_some() {
            return Err(GasEvaporationError::Invalid(
                "unknown or repeated native surface mass selection",
            ));
        }
        match liquid.project_interface_mass_exchange(inputs[index].coordinates, -mass) {
            Ok(projection) => {
                selected[index] = Some(projection);
                Ok(-projection.mass_change_kg)
            }
            Err(error) => {
                mass_error = Some(error);
                Err(GasEvaporationError::Invalid(
                    "original native interface mass selection refused",
                ))
            }
        }
    };
    let proposal = match heat_controls {
        None => advance_chamber_with_gas_mass_selection(
            &sources,
            gas,
            catalog,
            duration_s,
            controls.chamber,
            select_mass,
            || gate.is_requested(),
        )
        .map(|candidate| {
            let mass_defect = candidate.gas_exchange().step().mass_residual_kg;
            let species_defect = candidate
                .gas_exchange()
                .step()
                .maximum_species_mass_residual_kg;
            let (surfaces, gas, catalog, phase) = candidate.into_parts();
            (
                surfaces,
                gas,
                catalog,
                phase.clone(),
                phase,
                None,
                Vec::new(),
                mass_defect,
                species_defect,
            )
        }),
        Some(heat_controls) => advance_chamber_with_gas_and_heat_mass_selection(
            &source_heat,
            gas,
            catalog,
            duration_s,
            controls.chamber,
            heat_controls,
            select_mass,
            || gate.is_requested(),
        )
        .map(|candidate| {
            let phase = candidate.phase().step().clone();
            let surfaces = candidate.phase().surfaces().to_vec();
            let mass_defect = candidate.phase().gas_exchange().step().mass_residual_kg;
            let species_defect = candidate
                .phase()
                .gas_exchange()
                .step()
                .maximum_species_mass_residual_kg;
            let (_, gas, catalog, exchange, heat, heat_surfaces) = candidate.into_parts();
            (
                surfaces,
                gas,
                catalog,
                phase,
                exchange,
                Some(heat),
                heat_surfaces,
                mass_defect,
                species_defect,
            )
        }),
    };
    if let Some(error) = chart_error.into_inner().or(mass_error) {
        return Err(error.into());
    }
    let (
        phase_surfaces,
        next_gas,
        catalog,
        phase,
        exchange,
        film_heat,
        heat_surfaces,
        gas_mass_defect,
        gas_species_defect,
    ) = proposal?;
    if phase_surfaces.len() != inputs.len()
        || (film_heat.is_some() && heat_surfaces.len() != inputs.len())
    {
        return Err(Error::Invalid("joint source omitted a native surface"));
    }
    let mut transfers = Vec::with_capacity(inputs.len());
    for (index, (input, phase)) in inputs.iter().zip(&phase_surfaces).enumerate() {
        poll(gate)?;
        let projection = selected[index].ok_or(Error::Invalid(
            "phase owner did not select every native mass",
        ))?;
        if phase.id != index as u64 || -phase.transferred_mass_kg != projection.mass_change_kg {
            return Err(Error::Invalid(
                "joint source changed its native mass selection",
            ));
        }
        if film_heat.is_some()
            && (heat_surfaces[index].id != phase.id
                || heat_surfaces[index].phase_initial_liquid.mass_kg != input.initial_mass_kg)
        {
            return Err(Error::Invalid(
                "joint sensible heat changed a native surface inventory",
            ));
        }
        let delta_h = phase.liquid.specific_enthalpy_j_kg - input.initial_h_j_kg;
        let heat = projection.final_mass_kg * delta_h;
        if !heat.is_finite() || (delta_h != 0.0 && heat == 0.0) {
            return Err(Error::Invalid(
                "native phase caloric increment is unrepresentable",
            ));
        }
        transfers.push(FreeSurfaceMassExchange {
            coordinates: input.coordinates,
            mass_kg: projection.mass_change_kg,
            specific_enthalpy_j_kg: input.initial_h_j_kg,
            heat_j: heat,
        });
    }
    let mut next_liquid = liquid.clone();
    let mass_exchange = next_liquid.exchange_interface_mass(&transfers, controls.native, gate)?;
    let mut native_samples = BTreeMap::new();
    for sample in &mass_exchange.samples {
        poll(gate)?;
        native_samples.insert(sample.coordinates, *sample);
    }
    if native_samples.len() != inputs.len() {
        return Err(Error::Invalid(
            "native transaction omitted a source surface",
        ));
    }
    let mut native_caloric = Total::default();
    let mut absolute_caloric = Total::default();
    let mut absolute_mass = Total::default();
    absolute_mass.add(gas_mass_defect.abs())?;
    absolute_mass.add(gas_species_defect.abs())?;
    absolute_mass.add(mass_exchange.absolute_mass_defect_kg)?;
    let mut reports = Vec::with_capacity(inputs.len());
    for (index, (input, phase)) in inputs.iter().zip(phase_surfaces).enumerate() {
        poll(gate)?;
        // Native samples use their sparse-tile order, not Cartesian/source-ID
        // order. Match physical coordinates rather than zipping report arrays.
        let native = native_samples
            .remove(&input.coordinates)
            .ok_or(Error::Invalid(
                "native transaction returned a different surface",
            ))?;
        if native.mass_change_kg != -phase.transferred_mass_kg {
            return Err(Error::Invalid(
                "native phase transfer differs from its selected gas debit",
            ));
        }
        let projection = selected[index].expect("selected before native transaction");
        let heat = film_heat.as_ref().map(|_| heat_surfaces[index]);
        let source_relative_storage = if let Some(heat) = &heat {
            // Keep the exact heated h supplied to phase. Evaporating material
            // carries that heated value, so shifting its reference back to h0
            // includes the signed mass term rather than adding sensible heat
            // or latent energy a second time.
            let heated_h = heat.phase_initial_liquid.specific_enthalpy_j_kg;
            heat.liquid_enthalpy_change_j
                + phase.liquid.mass_kg * (phase.liquid.specific_enthalpy_j_kg - heated_h)
                - phase.transferred_mass_kg * (heated_h - input.initial_h_j_kg)
        } else {
            phase.liquid.mass_kg * (phase.liquid.specific_enthalpy_j_kg - input.initial_h_j_kg)
        };
        let defect = projection.final_mass_kg
            * (native.specific_enthalpy_j_kg[1] - input.initial_h_j_kg)
            - source_relative_storage;
        native_caloric.add(defect)?;
        absolute_caloric.add(defect.abs())?;
        absolute_mass
            .add((native.mass_change_kg - (phase.liquid.mass_kg - input.initial_mass_kg)).abs())?;
        reports.push(NativeLiquidGasSurfaceStep {
            coordinates: input.coordinates,
            phase,
            heat,
            native,
            native_caloric_defect_j: defect,
        });
    }
    poll(gate)?;
    let flow = next_liquid.step(gate)?;
    let native_transport_mass = flow.mass_after_kg - flow.mass_before_kg;
    absolute_mass.add(native_transport_mass.abs())?;
    let native_mass = mass_exchange.mass_change_kg + native_transport_mass;
    let mass_residual_kg = native_mass + (next_gas.mass_kg() - gas.mass_kg());
    let volatile_mass_residual_kg = native_mass
        + (next_gas
            .species_mass_kg(&binding.species)
            .map_err(GasEvaporationError::from)?
            - gas
                .species_mass_kg(&binding.species)
                .map_err(GasEvaporationError::from)?);
    let native_phase_caloric_defect_j = native_caloric.value();
    let energy_residual_j =
        exchange.energy_residual_j + native_phase_caloric_defect_j + flow.balance_defect_j;
    let absolute_energy_defect_j = exchange.sum_absolute_energy_residual_j
        + absolute_caloric.value()
        + flow.balance_defect_j.abs();
    let absolute_mass_defect_kg = absolute_mass.value();
    let native_liquid_enthalpy_change_j =
        mass_exchange.enthalpy_change_j + (flow.enthalpy_after_j - flow.enthalpy_before_j);
    if ![
        mass_residual_kg,
        volatile_mass_residual_kg,
        absolute_mass_defect_kg,
        native_phase_caloric_defect_j,
        energy_residual_j,
        absolute_energy_defect_j,
        native_liquid_enthalpy_change_j,
    ]
    .into_iter()
    .all(f64::is_finite)
    {
        return Err(Error::Invalid(
            "combined distributed native/gas account is unrepresentable",
        ));
    }
    let worst_mass = mass_residual_kg
        .abs()
        .max(volatile_mass_residual_kg.abs())
        .max(absolute_mass_defect_kg);
    if worst_mass > controls.mass_tolerance_kg {
        return Err(Error::MassBalance {
            residual_kg: worst_mass,
            tolerance_kg: controls.mass_tolerance_kg,
        });
    }
    let worst_energy = energy_residual_j.abs().max(absolute_energy_defect_j);
    if worst_energy > controls.energy_tolerance_j {
        return Err(Error::EnergyBalance {
            residual_j: worst_energy,
            tolerance_j: controls.energy_tolerance_j,
        });
    }
    poll(gate)?;
    Ok(DistributedLiquidGasCandidate {
        liquid: next_liquid,
        gas: next_gas,
        catalog,
        report: DistributedLiquidGasReport {
            duration_s,
            phase,
            exchange,
            heat: film_heat,
            surfaces: reports,
            mass_exchange,
            flow,
            native_liquid_enthalpy_change_j,
            native_phase_caloric_defect_j,
            mass_residual_kg,
            volatile_mass_residual_kg,
            absolute_mass_defect_kg,
            energy_residual_j,
            absolute_energy_defect_j,
        },
    })
}
