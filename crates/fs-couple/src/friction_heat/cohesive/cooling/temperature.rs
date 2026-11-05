//! Shared fixed-jump thermal weakening of deformable bonds and their cooling.
//!
//! The original source evaluates actual caloric-cell temperatures. The two
//! native history owners then perform the SAME coefficient-change/relaxation
//! transition, and publish only when points and physical energy accounts agree.
//! Parameter work and new fracture loss stay explicit; neither is silently heat.

use super::{CohesiveCoolingError, EnergySum, require_energy, require_synchronized_owners};
use fs_conduction::transient::enthalpy::deposition::fracture::{
    FractureCooling, temperature as native_temperature,
};
use fs_exec::CancelGate;
use fs_matdb::InterfaceSystemCard;
use fs_material::cohesive::source::ResolvedBilinearCohesive;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsTemperatureStep};
use fs_solid::cohesive::{CohesivePatch, CohesiveSystem};
use fs_thermal::enthalpy::DepositionError;

/// The existing native local and total energy gates apply to both replays and
/// their independent endpoint balance. Both tolerances must be positive.
pub use fs_solid::cohesive::temperature::TemperatureConfig as CohesiveCoolingTemperatureControls;

/// One original card and declared cell-temperature approximation per patch.
/// Supply exactly the ORIGINAL mechanical patch order; IDs guard that ordering.
/// A/B receiving-cell and material correspondence is the original cooling
/// binding's caller obligation, not inferred from this scalar temperature.
#[derive(Clone, Copy)]
pub struct CohesiveCoolingTemperatureBinding<'c> {
    pub patch_id: u64,
    pub card: &'c InterfaceSystemCard,
    /// T_interface = w T_A + (1-w) T_B, w in [0,1]. No resistance/capacity
    /// weighting or resolved through-thickness temperature is inferred.
    pub temperature_weight_a: f64,
}

/// Coefficient preparation only; accepted displacement/history is unchanged.
///
/// Each original patch queries its card at most once, in patch order, with
/// its original selected claims, nonthermal coordinates and model bounds.
/// An unchanged temperature retains the exact original source receipt.
/// Use the returned sources with the existing
/// `CohesiveSystem::with_temperature_sources_controlled`, then pass that native
/// target to `propose_deformable_cohesive_temperature`. Keeping source storage
/// outside the borrowed target avoids a second system or self-referential state.
pub fn prepare_deformable_cohesive_temperature_sources(
    motion: &CohesiveMotion<'_, '_>,
    cooling: &FractureCooling<'_, '_>,
    bindings: &[CohesiveCoolingTemperatureBinding<'_>],
    temperature_axis: &str,
    gate: &CancelGate,
) -> Result<Vec<ResolvedBilinearCohesive>, CohesiveCoolingError> {
    require_synchronized_owners(motion, cooling, gate)?;
    let patches = motion.state().system().patches();
    if bindings.len() != patches.len() {
        return Err(CohesiveCoolingError::Invalid(
            "one thermal binding per original patch required",
        ));
    }
    let mut sources = Vec::with_capacity(patches.len());
    for (patch, binding) in patches.iter().zip(bindings) {
        check(gate)?;
        let original = source(patch, binding)?;
        let temperature_k =
            cooling.interface_temperature_k(patch.id(), binding.temperature_weight_a)?;
        let next = if original.resolved().state_coordinate(temperature_axis) == Some(temperature_k)
        {
            original
                .require_temperature_successor(original, temperature_axis)
                .map_err(|error| CohesiveCoolingError::Source {
                    patch: patch.id(),
                    error,
                })?;
            original.clone()
        } else {
            original
                .at_temperature(binding.card, temperature_axis, temperature_k)
                .map_err(|error| CohesiveCoolingError::Source {
                    patch: patch.id(),
                    error,
                })?
        };
        sources.push(next);
    }
    check(gate)?;
    Ok(sources)
}

/// One physical transition, observed by both native owners. The thermal report
/// is a replay of the same coefficient work/fracture, never a second loss.
pub struct CohesiveCoolingTemperatureCandidate<'s, 'a, 't> {
    pub mechanical: DynamicsTemperatureStep<'s, 'a>,
    pub cooling: FractureCooling<'t, 'a>,
    pub contact: native_temperature::TemperatureReport,
    /// Independently reduced original-area delta(U_bond + D_fracture) minus
    /// signed parameter work. Motion, bulk storage, time and cell heat stay fixed.
    pub combined_energy_residual_j: f64,
}

/// Transfer actual-cell temperature feedback into BOTH original bond owners.
///
/// The caller-owned target is the original native system with temperature
/// sources. Source/card/claim identity, fixed geometry, mass, bulk operator,
/// displacement, velocity and trajectory controls remain with the native motion
/// transfer. The cooling transfer independently verifies target temperatures
/// against current caloric cells and updates its same bonded/failed contact
/// conductance. Both complete point histories and work/loss accounts must agree.
///
/// All mechanical/thermal clocks must already agree: finish a timed segment
/// before changing its coefficients. No time or displacement work occurs here.
/// Thermal weakening may create damage or open heat contact at fixed jump;
/// cooling cannot rebond. A later `propose_deformable_cohesive_cooling` uses the
/// new traction, viscosity and contact topology. Parameter work must be supplied
/// by the caller's declared physical owner; released fracture energy stays
/// nonthermal. No thermal expansion, strong implicit feedback or entropy/heat
/// closure is inferred. Refine the staggered sequence to resolve its dynamics.
#[allow(clippy::too_many_arguments)]
pub fn propose_deformable_cohesive_temperature<'old, 'a, 's, 'b, 't>(
    motion: &CohesiveMotion<'old, 'a>,
    cooling: &FractureCooling<'t, 'a>,
    target: &'s CohesiveSystem<'b>,
    bindings: &[CohesiveCoolingTemperatureBinding<'_>],
    temperature_axis: &str,
    controls: CohesiveCoolingTemperatureControls,
    gate: &CancelGate,
) -> Result<CohesiveCoolingTemperatureCandidate<'s, 'b, 't>, CohesiveCoolingError>
where
    'old: 's,
{
    require_synchronized_owners(motion, cooling, gate)?;
    let patches = motion.state().system().patches();
    if bindings.len() != patches.len()
        || target.patches().len() != patches.len()
        || [
            controls.local_energy_tolerance_j_m2,
            controls.total_energy_tolerance_j,
        ]
        .iter()
        .any(|v| !v.is_finite() || *v <= 0.0)
    {
        return Err(CohesiveCoolingError::Invalid(
            "complete thermal target and positive energy gates required",
        ));
    }
    let mut updates = Vec::with_capacity(patches.len());
    for ((patch, next), binding) in patches.iter().zip(target.patches()).zip(bindings) {
        check(gate)?;
        if next.id() != patch.id() {
            return Err(CohesiveCoolingError::Invalid(
                "thermal target changed original patch order",
            ));
        }
        updates.push(native_temperature::TemperatureUpdate {
            interface: patch.id(),
            previous: source(patch, binding)?,
            successor: next.source().ok_or(CohesiveCoolingError::Invalid(
                "temperature target requires original sourced patches",
            ))?,
            temperature_weight_a: binding.temperature_weight_a,
        });
    }
    let mechanical =
        motion.transfer_temperature_controlled(target, temperature_axis, controls, &mut || {
            gate.is_requested()
        })?;
    let (next, contact) = cooling.transfer_temperature(
        motion.time_s(),
        &updates,
        temperature_axis,
        native_temperature::TemperatureConfig {
            local_energy_tolerance_j_m2: controls.local_energy_tolerance_j_m2,
            total_energy_tolerance_j: controls.total_energy_tolerance_j,
        },
        || {
            if gate.is_requested() {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        },
    )?;
    require_synchronized_owners(&mechanical.motion, &next, gate)?;

    let mut storage = EnergySum::default();
    let mut fracture = EnergySum::default();
    for ((patch, old), point) in patches
        .iter()
        .zip(motion.state().points())
        .zip(mechanical.motion.state().points())
    {
        check(gate)?;
        storage.add(patch.area_m2() * (point.stored_energy_j_m2() - old.stored_energy_j_m2()))?;
        fracture.add(patch.area_m2() * (point.dissipated_j_m2() - old.dissipated_j_m2()))?;
    }
    let storage = storage.value()?;
    let fracture = fracture.value()?;
    let mut balance = EnergySum::default();
    balance.add(storage)?;
    balance.add(fracture)?;
    balance.add(-mechanical.energy.parameter_work_j)?;
    let combined_energy_residual_j = balance.value()?;
    // These are two observations of ONE physical interface, not additive
    // energies. Independently read its represented history/storage endpoints.
    for residual in [
        contact.parameter_work_j - mechanical.energy.parameter_work_j,
        contact.stored_change_j - mechanical.energy.stored_change_j,
        contact.contact.fracture_dissipation_j - mechanical.energy.damage_dissipated_j,
        storage - mechanical.energy.stored_change_j,
        fracture - mechanical.energy.damage_dissipated_j,
        combined_energy_residual_j,
        combined_energy_residual_j + mechanical.energy.defect_j,
        combined_energy_residual_j + contact.defect_j,
    ] {
        require_energy(residual, controls.total_energy_tolerance_j)?;
    }
    check(gate)?;
    Ok(CohesiveCoolingTemperatureCandidate {
        mechanical,
        cooling: next,
        contact,
        combined_energy_residual_j,
    })
}

fn source<'a>(
    patch: &CohesivePatch<'a>,
    binding: &CohesiveCoolingTemperatureBinding<'_>,
) -> Result<&'a ResolvedBilinearCohesive, CohesiveCoolingError> {
    let source = patch.source().ok_or(CohesiveCoolingError::Invalid(
        "thermal feedback requires original sourced cohesive patches",
    ))?;
    if binding.patch_id != patch.id()
        || binding.card.content_hash() != source.resolved().card_identity()
    {
        return Err(CohesiveCoolingError::Invalid(
            "thermal binding changed original patch or card",
        ));
    }
    Ok(source)
}

fn check(gate: &CancelGate) -> Result<(), CohesiveCoolingError> {
    if gate.is_requested() {
        Err(DepositionError::Cancelled.into())
    } else {
        Ok(())
    }
}
