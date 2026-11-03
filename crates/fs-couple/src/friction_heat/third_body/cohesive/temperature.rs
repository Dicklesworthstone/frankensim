//! Original source-temperature transfer and actual finite adhesive heating.
//!
//! Temperature is prescribed, not solved from a cohesive free energy. The
//! caller owns the original pinned successor sources and their native system;
//! the existing material transfer retains all accepted irreversible history.
//! Signed coefficient work is external work, and fixed-jump thermal fracture
//! retains its surface account. Neither becomes caloric heat. The subsequent
//! native motion supplies only actual surviving-bond viscosity to the existing
//! finite adhesive cell and its two physical conduction paths.

use super::{
    CohesiveHeatError, CohesiveLayerHeatBinding, CohesiveLayerHeatControls,
    CohesiveLayerHeatSource, HeatedCohesiveLayerCandidate, propose_deformable_cohesive_layer_heat,
};
use fs_exec::CancelGate;
use fs_solid::cohesive::CohesiveSystem;
use fs_solid::cohesive::dynamics::{CohesiveMotion, DynamicsStep};
use fs_solid::cohesive::temperature::{TemperatureConfig, TemperatureEnergy};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation};

/// Caller-owned temperature successors on the exact original native geometry.
/// Prepare `target` with `CohesiveSystem::with_temperature_sources`; the native
/// transfer checks every source, patch, trace and elastic operator again against
/// this accepted motion. Original claim choices and nonthermal coordinates stay
/// pinned. These tolerances do not replace the motion or heat acceptance gates.
#[derive(Clone, Copy)]
pub struct CohesiveLayerHeatTemperature<'s, 'a> {
    pub target: &'s CohesiveSystem<'a>,
    pub axis: &'s str,
    pub controls: TemperatureConfig,
}

/// Both complete immutable candidates, including the preceding source change.
pub struct TemperatureHeatedCohesiveLayerCandidate<'s, 'a, 't, P> {
    /// Motion and caloric increment AFTER the zero-duration source transfer.
    /// Its local residuals and fracture field refer only to that increment.
    pub heated: HeatedCohesiveLayerCandidate<'s, 'a, 't, P>,
    /// Fixed-jump source work, new damage release and storage change in joules.
    pub temperature: TemperatureEnergy,
    /// Thermal weakening plus subsequent native fracture, never caloric heat.
    pub fracture_work_j: f64,
    /// Original-to-final delta(K + bulk/bond storage + H) + total fracture
    /// - applied work - prescribed-jump work - parameter work - other heat.
    pub combined_energy_residual_j: f64,
    /// Same balance assembled from the original three owners' residuals.
    pub summed_energy_residual_j: f64,
}

/// Change prescribed bond temperature, advance native motion and retain its
/// actual viscous heat in the finite adhesive layer, as one immutable proposal.
///
/// Transfer occurs once at the accepted jump and retained damage, followed by
/// irreversible relaxation onto the new envelope. It preserves displacement,
/// velocity, inertia and clock. Mechanical coefficients then remain frozen for
/// the increment; thermal Newton trials never repeat transfer or damage. The
/// transport callback sees the resulting actual native damage and trial layer
/// temperature, so original resistance/conductivity sources can respond to both.
/// A payload must remain unpublished until this operation succeeds.
///
/// Signed parameter work is an explicit external input, not a debit from body
/// or adhesive enthalpy. This supplies prescribed thermal weakening, not closed
/// thermoelastic feedback, a cohesive entropy model, or an unseen temperature
/// path. The ordinary layer operation's complete source, fixed-reference,
/// closed-path and heat-location restrictions still apply. Every failure leaves
/// the original motion, source history and finite thermal state unchanged.
#[allow(clippy::too_many_arguments)]
pub fn propose_deformable_cohesive_layer_heat_with_temperature<'m, 'a, 's, 'b, 't, P, E>(
    motion: &CohesiveMotion<'m, 'a>,
    mean_force_n: &[f64],
    endpoint_prescribed_jumps_m: &[[f64; 3]],
    dt_s: f64,
    temperature: CohesiveLayerHeatTemperature<'s, 'b>,
    binding: CohesiveLayerHeatBinding<'_>,
    thermal: &DepositionSimulation<'t>,
    controls: CohesiveLayerHeatControls,
    energy_guess_j: &[f64],
    resolve: impl Fn(
        &DepositionSimulation<'t>,
        f64,
        &DynamicsStep<'s, 'b>,
    ) -> Result<CohesiveLayerHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<TemperatureHeatedCohesiveLayerCandidate<'s, 'b, 't, P>, CohesiveHeatError<E>>
where
    'm: 's,
{
    let check = || {
        if gate.is_requested() {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    };
    check()?;
    let ([patch], [before]) = (motion.state().system().patches(), motion.state().points()) else {
        return Err(CohesiveHeatError::Invalid(
            "one complete original adhesive patch required for source temperature and layer heat",
        ));
    };
    let transferred = motion.transfer_temperature_controlled(
        temperature.target,
        temperature.axis,
        temperature.controls,
        &mut || gate.is_requested(),
    )?;
    let heated = propose_deformable_cohesive_layer_heat(
        &transferred.motion,
        mean_force_n,
        endpoint_prescribed_jumps_m,
        dt_s,
        binding,
        thermal,
        controls,
        energy_guess_j,
        resolve,
        gate,
    )?;
    check()?;
    let mechanical = &heated.mechanical;
    let after = &mechanical.motion.state().points()[0];
    let fracture = transferred.energy.damage_dissipated_j + heated.fracture_work_j;
    let history_fracture = patch.area_m2() * (after.dissipated_j_m2() - before.dissipated_j_m2());
    // Source transfer preserves mass, displacement and velocity exactly: the
    // native motion's initial kinetic/bulk energies are the ORIGINAL energies.
    // Read bond storage, accrued fracture and caloric storage at both accepted
    // endpoints independently of the transfer/increment residuals.
    let combined = mechanical.kinetic_after_j - mechanical.kinetic_before_j
        + mechanical.bulk_stored_after_j
        - mechanical.bulk_stored_before_j
        + patch.area_m2() * (after.stored_energy_j_m2() - before.stored_energy_j_m2())
        + heated.thermal.energy_j()
        - thermal.energy_j()
        + history_fracture
        - mechanical.external_work_j
        - mechanical.prescribed_jump_work_j
        - transferred.energy.parameter_work_j
        - heated.heat.external_heat_j;
    let summed = heated.summed_energy_residual_j - transferred.energy.defect_j;
    if [combined, summed, history_fracture - fracture]
        .iter()
        .any(|v| !v.is_finite() || v.abs() > controls.combined_energy_tolerance_j)
        || !fracture.is_finite()
        || fracture < 0.0
    {
        return Err(CohesiveHeatError::Invalid(
            "source temperature, native adhesive and finite layer energy do not balance",
        ));
    }
    check()?;
    Ok(TemperatureHeatedCohesiveLayerCandidate {
        heated,
        temperature: transferred.energy,
        fracture_work_j: fracture,
        combined_energy_residual_j: combined,
        summed_energy_residual_j: summed,
    })
}
