//! Source-bound wear on the original curved friction/finite-heat trajectory.
//! Endpoint k/H and held normal force multiply the actual unsigned path, even
//! when sliding curves or reverses. Wear is a small-removal, fixed-geometry
//! constitutive diagnostic: no surface recession, mass removal, debris energy,
//! or replacement thermal partition is inferred from its volume.

use super::trajectory::{
    DirectionalTrajectoryHeatError, HeatedDirectionalTrajectoryReport,
    apply_source_directional_trajectory_heat,
};
use super::*;
use fs_contact::interface_binding::wear::{SourceArchardWear, SourceWearError};
use fs_tribo::WearAdvance;
use fs_tribo::tangential::elliptic::relaxation::EllipticRelaxationControls;

/// The original mechanics/heat account plus accepted A/B wear increments.
pub struct HeatedDirectionalWearReport<P> {
    pub phase: HeatedDirectionalTrajectoryReport<P>,
    pub wear: [WearAdvance; 2],
}

#[derive(Debug)]
pub enum DirectionalWearHeatError<E> {
    Invalid(&'static str),
    Phase(DirectionalTrajectoryHeatError<E>),
    Wear(SourceWearError),
    Cancelled,
}
impl<E: core::fmt::Debug> core::fmt::Display for DirectionalWearHeatError<E> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "directional friction/heat/wear: {self:?}")
    }
}
impl<E: std::error::Error + 'static> std::error::Error for DirectionalWearHeatError<E> {}

/// Advance actual sliding, finite heat and both original surface wear volumes.
///
/// The original trajectory owner resolves the full-matrix curved slip and its
/// finite arrest. Native Archard uses that unsigned path, never the norm of net
/// displacement or an impulse-derived proxy. k comes from the SAME ordered
/// interface and H from each original bulk card at its own caloric endpoint.
/// Claims remain pinned across calls. Phase-average speed, endpoint k/H and
/// held load are a time discretization; refine physical time for varying data.
///
/// Only the two participating bodies are copied. Their momenta, original
/// friction history, complete finite thermal state and wear publish together
/// after every source, constitutive and cancellation check. A late hardness
/// refusal therefore cannot leave already-heated or partially worn bodies.
/// Other source callbacks retain the original requirement of private staging.
///
/// Wear volume does not determine a physical energy share. The original heat
/// accounting is retained, and no volume-to-energy or geometry law is supplied.
#[allow(clippy::too_many_arguments)]
pub fn apply_source_directional_trajectory_heat_with_wear<'m, P, E>(
    bodies: &mut [RigidBodyState],
    properties: &[MassProperties],
    history: &mut SourceThermalEllipticState,
    law: &SourceThermalEllipticCoulomb,
    contact: RigidDirectionalContact,
    thermal: &mut DepositionSimulation<'m>,
    route: &FrictionHeatRoute,
    material: DirectionalHeatMaterial<'_>,
    wear: &mut SourceArchardWear,
    dt_s: f64,
    controls: DirectionalHeatSolve,
    relaxation: EllipticRelaxationControls,
    energy_guess_j: &[f64],
    interface_temperature_guess_k: f64,
    resolve: impl Fn(
        &DepositionSimulation<'m>,
        [EquilibriumPhaseState; 2],
    ) -> Result<DirectionalHeatSource<P>, E>,
    gate: &CancelGate,
) -> Result<HeatedDirectionalWearReport<P>, DirectionalWearHeatError<E>> {
    if gate.is_requested() {
        return Err(DirectionalWearHeatError::Cancelled);
    }
    wear.admit_friction(law)
        .map_err(DirectionalWearHeatError::Wear)?;
    let binding = contact.binding;
    if bodies.len() != properties.len()
        || bodies.len() > controls.mechanical.maximum_bodies
        || controls.maximum_cells == 0
        || thermal.cells().len() > controls.maximum_cells
        || binding.body_a >= bodies.len()
        || binding.body_b >= bodies.len()
        || binding.body_a == binding.body_b
    {
        return Err(DirectionalWearHeatError::Invalid(
            "bounded original body pair and finite thermal state required",
        ));
    }
    let mut pair = [bodies[binding.body_a], bodies[binding.body_b]];
    let pair_properties = [properties[binding.body_a], properties[binding.body_b]];
    let mut pair_contact = contact;
    pair_contact.binding.body_a = 0;
    pair_contact.binding.body_b = 1;
    let mut next_history = history.clone();
    let mut next_thermal = thermal.clone();
    let phase = apply_source_directional_trajectory_heat(
        &mut pair,
        &pair_properties,
        &mut next_history,
        law,
        pair_contact,
        &mut next_thermal,
        route,
        material,
        dt_s,
        controls,
        relaxation,
        energy_guess_j,
        interface_temperature_guess_k,
        resolve,
        gate,
    )
    .map_err(DirectionalWearHeatError::Phase)?;
    let query = next_history.friction_sources()[0].query();
    let (next_wear, increments) = wear
        .advance_candidate(
            law,
            &query,
            material.point.materials,
            phase.phases,
            contact.normal_load_n,
            phase.actual_distance_m,
        )
        .map_err(DirectionalWearHeatError::Wear)?;
    if gate.is_requested() {
        return Err(DirectionalWearHeatError::Cancelled);
    }
    // No fallible work follows publication of these original owner states.
    bodies[binding.body_a] = pair[0];
    bodies[binding.body_b] = pair[1];
    *history = next_history;
    *thermal = next_thermal;
    *wear = next_wear;
    Ok(HeatedDirectionalWearReport {
        phase,
        wear: increments,
    })
}
