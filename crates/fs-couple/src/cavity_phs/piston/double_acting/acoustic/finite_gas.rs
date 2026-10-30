//! Finite regional gas transport on the original shared acoustic motion.
//!
//! The native piston or electrical owner advances once against the initial
//! acoustic medium. One original material-transport snapshot then supplies the
//! finite regional endpoints. Each region pays the actual native gas-rebind
//! work on its already-marched load; no acoustic or mechanical history is
//! advanced again inside a caloric trial.
//!
//! These inventories belong to distinct, fixed-geometry acoustic regions, not
//! either moving thermodynamic chamber. Acoustic flow does not imply material
//! transport. Mean pressure requires a balanced reference configuration or an
//! authored compensating force; no additional absolute-pressure force or swept
//! thermodynamic volume is inferred. The original frozen-medium, first-order
//! partition and linear perturbation assumptions remain.

use super::{AcousticMotionError, AcousticMotionFace, AcousticMotionSolve, AcousticMotionStep};
use crate::bernoulli_aperture::wall::material::{
    FiniteGasTransportNode, finite_gas_transport_candidate_with_paths,
};
use crate::cavity_phs::piston::PistonError;
use crate::cavity_phs::piston::double_acting::{
    DoubleActingGasPiston, DoublePistonInput, DoublePistonSolve, DoublePistonStep,
};
use crate::dg_port::DgPortLoad;
use crate::dg_port::finite_gas::transport::{
    FiniteGasDgPortTransportStep, FiniteGasDgTransportError, FiniteGasDgTransportInput,
    FiniteGasDgTransportSolve,
};
use crate::dg_port::finite_gas::{
    FiniteGasDgBackground, FiniteGasDgError, FiniteGasDgPortStep, admit_initial, prepare_after_step,
};
use fs_duct::dg::{DgStepOutput, PortInput};
use fs_material::thermochemical_gas::energy::GasEnergyError;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    snapshot::MAX_SNAPSHOT_PARCELS,
    species::network::transport::{GasTransportError, GasTransportPaths},
};
use std::collections::{BTreeMap, BTreeSet};

/// One original acoustic load and its distinct finite regional inventory.
///
/// The parcel key is the original identity used by material transport. Its
/// volume must equal this load's uniform main bore plus main-gas finite terminal
/// passages. Independent branch inventories remain with their original sources.
#[derive(Clone, Copy, Debug)]
pub struct FiniteGasAcousticFace<'a> {
    pub parcel: u64,
    pub load: &'a DgPortLoad,
    pub signed_area_m2: f64,
}

/// Actual native interval exchanges and the final regional gas/load endpoints.
#[derive(Debug)]
pub struct FiniteGasAcousticStep {
    /// Original force, flow and interval work, with refreshed endpoint ledgers.
    pub motion: AcousticMotionStep,
    /// One original material snapshot and complete finite-gas acoustic accounts.
    pub transport: FiniteGasDgPortTransportStep,
    /// Parcel keys in the authored face order used by the motion report.
    pub face_parcels: Vec<u64>,
    /// Native mechanical/electrical storage plus all regional gas and acoustic
    /// storage/loss, less actual external supply. Internal transfer is counted once.
    pub total_energy_residual_j: f64,
}

/// A complete proposal; publish every returned owner together.
///
/// The native owner remains its original type and need not implement Clone.
/// Regional inventories are additional physical stores, not copies of either
/// piston chamber. No accepted owner changes if any proposal is refused.
pub struct FiniteGasAcousticCandidate<N> {
    native: N,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    loads: BTreeMap<u64, DgPortLoad>,
    step: FiniteGasAcousticStep,
}

impl<N> FiniteGasAcousticCandidate<N> {
    pub const fn native(&self) -> &N {
        &self.native
    }

    pub fn parcels(&self) -> &BTreeMap<u64, ClosedGasParcel> {
        &self.parcels
    }

    pub fn loads(&self) -> &BTreeMap<u64, DgPortLoad> {
        &self.loads
    }

    pub const fn step(&self) -> &FiniteGasAcousticStep {
        &self.step
    }

    pub const fn total_energy_residual_j(&self) -> f64 {
        self.step.total_energy_residual_j
    }

    pub fn into_parts(
        self,
    ) -> (
        N,
        BTreeMap<u64, ClosedGasParcel>,
        BTreeMap<u64, DgPortLoad>,
        FiniteGasAcousticStep,
    ) {
        (self.native, self.parcels, self.loads, self.step)
    }
}

/// A failure retains every accepted native, acoustic and gas owner.
#[derive(Debug)]
pub enum FiniteGasAcousticError<M> {
    Invalid(&'static str),
    Motion(AcousticMotionError<M>),
    Transport(FiniteGasDgTransportError),
    Cancelled,
}

impl<M> From<AcousticMotionError<M>> for FiniteGasAcousticError<M> {
    fn from(error: AcousticMotionError<M>) -> Self {
        match error {
            AcousticMotionError::Cancelled => Self::Cancelled,
            error => Self::Motion(error),
        }
    }
}

impl<M> From<FiniteGasDgTransportError> for FiniteGasAcousticError<M> {
    fn from(error: FiniteGasDgTransportError) -> Self {
        match error {
            FiniteGasDgTransportError::Cancelled
            | FiniteGasDgTransportError::Transport(GasTransportError::Cancelled)
            | FiniteGasDgTransportError::Gas {
                error: GasEnergyError::Cancelled,
                ..
            }
            | FiniteGasDgTransportError::Acoustic {
                error: FiniteGasDgError::Cancelled,
                ..
            }
            | FiniteGasDgTransportError::Acoustic {
                error: FiniteGasDgError::Background(GasEnergyError::Cancelled),
                ..
            } => Self::Cancelled,
            error => Self::Transport(error),
        }
    }
}

impl<M: core::fmt::Debug> core::fmt::Display for FiniteGasAcousticError<M> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite-gas shared acoustic motion: {self:?}")
    }
}

impl<M: std::error::Error + 'static> std::error::Error for FiniteGasAcousticError<M> {}

fn poll<M>(cancelled: &mut impl FnMut() -> bool) -> Result<(), FiniteGasAcousticError<M>> {
    if cancelled() {
        Err(FiniteGasAcousticError::Cancelled)
    } else {
        Ok(())
    }
}

/// An already-marched original load, never a second acoustic integrator.
///
/// Cloning the admitted staging map only copies references and an empty endpoint.
/// Native caloric trials clone the marched load inside prepare_after_step, where
/// the original gas-to-wave work solver owns those private proposals.
#[derive(Clone)]
struct StagedPort<'a> {
    original: &'a DgPortLoad,
    marched: &'a DgPortLoad,
    output: &'a DgStepOutput,
    endpoint: Option<DgPortLoad>,
}

impl StagedPort<'_> {
    fn current(&self) -> &DgPortLoad {
        self.endpoint.as_ref().unwrap_or(self.original)
    }
}

impl FiniteGasTransportNode for StagedPort<'_> {
    type Step<S> = FiniteGasDgPortStep<S>;
    type Input = FiniteGasDgTransportInput;

    fn transport_input(
        input: &Self::Input,
        bound: bool,
    ) -> Result<&FiniteGasDgTransportInput, &'static str> {
        <DgPortLoad as FiniteGasTransportNode>::transport_input(input, bound)
    }

    fn dt(&self) -> f64 {
        self.current().propagator().dt()
    }

    fn step_count(&self) -> u64 {
        self.current().state().step_count()
    }

    fn stored_energy_j(&self) -> f64 {
        self.current().stored_energy_j()
    }

    fn admit<E, C: FnMut() -> bool>(
        &self,
        gas: &ClosedGasParcel,
        input: &Self::Input,
        cancelled: &mut C,
    ) -> Result<(), FiniteGasDgError<E>> {
        <DgPortLoad as FiniteGasTransportNode>::admit(self.original, gas, input, cancelled)
    }

    fn advance<S, E, C: FnMut() -> bool>(
        &mut self,
        gas: &mut ClosedGasParcel,
        input: &Self::Input,
        evolve: impl FnMut(&ClosedGasParcel, f64, f64, &mut C) -> Result<FiniteGasDgBackground<S>, E>,
        cancelled: &mut C,
    ) -> Result<Self::Step<S>, FiniteGasDgError<E>> {
        if self.endpoint.is_some() {
            return Err(FiniteGasDgError::Invalid(
                "one shared-motion acoustic endpoint may be completed only once",
            ));
        }
        let controls = input.acoustic.ok_or(FiniteGasDgError::Invalid(
            "a staged shared-motion load requires finite-gas acoustic controls",
        ))?;
        let candidate = prepare_after_step(
            self.original,
            self.marched,
            self.output.clone(),
            gas,
            controls,
            evolve,
            cancelled,
        )?;
        *gas = candidate.gas;
        self.endpoint = Some(candidate.load);
        Ok(candidate.step)
    }

    fn accounts<S>(step: &Self::Step<S>) -> [f64; 3] {
        <DgPortLoad as FiniteGasTransportNode>::accounts(step)
    }

    fn source<S>(step: &Self::Step<S>) -> &S {
        &step.background
    }
}

/// Admit original regional ownership before any native motion proposal.
pub(super) fn prepare_faces<'a, M>(
    faces: &[FiniteGasAcousticFace<'a>],
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    acoustic_controls: AcousticMotionSolve,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<AcousticMotionFace<'a>>, FiniteGasAcousticError<M>> {
    poll(cancelled)?;
    if faces.is_empty()
        || faces.len() > acoustic_controls.maximum_faces
        || faces.len() > parcels.len()
        || parcels.is_empty()
        || parcels.len() > MAX_SNAPSHOT_PARCELS
        || !parcels.keys().eq(inputs.keys())
    {
        return Err(FiniteGasAcousticError::Invalid(
            "bounded nonempty acoustic faces and a complete original parcel/input map required",
        ));
    }
    let mut identities = BTreeSet::new();
    let mut native_faces = Vec::with_capacity(faces.len());
    for (index, face) in faces.iter().enumerate() {
        poll(cancelled)?;
        if !face.signed_area_m2.is_finite()
            || !parcels.contains_key(&face.parcel)
            || !identities.insert(face.parcel)
            || faces[..index]
                .iter()
                .any(|other| std::ptr::eq(other.load, face.load))
        {
            return Err(FiniteGasAcousticError::Invalid(
                "each finite acoustic face requires a distinct existing parcel and original load",
            ));
        }
        native_faces.push(AcousticMotionFace {
            load: face.load,
            signed_area_m2: face.signed_area_m2,
        });
    }
    for (&id, input) in inputs {
        poll(cancelled)?;
        if !input.ports.is_empty()
            || !input.supplied_heat_j.is_finite()
            || input.acoustic.is_some() != identities.contains(&id)
        {
            return Err(FiniteGasAcousticError::Invalid(
                "caller ports must be empty, heat finite and acoustic controls bound exactly to face parcels",
            ));
        }
        input.gas.validate().map_err(|error| {
            FiniteGasAcousticError::from(FiniteGasDgTransportError::Gas { parcel: id, error })
        })?;
    }
    for face in faces {
        poll(cancelled)?;
        let controls = inputs[&face.parcel]
            .acoustic
            .ok_or(FiniteGasAcousticError::Invalid(
                "a finite acoustic face lost its original gas-work controls",
            ))?;
        admit_initial::<GasEnergyError>(face.load, &parcels[&face.parcel], controls, cancelled)
            .map_err(|error| {
                FiniteGasAcousticError::from(FiniteGasDgTransportError::Acoustic {
                    parcel: face.parcel,
                    error,
                })
            })?;
    }
    poll(cancelled)?;
    Ok(native_faces)
}

/// Complete one original native proposal with regional transport and gas work.
#[allow(clippy::too_many_arguments)]
pub(super) fn finish<N, M, C: FnMut() -> bool>(
    faces: &[FiniteGasAcousticFace<'_>],
    parcels: &BTreeMap<u64, ClosedGasParcel>,
    inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
    paths: GasTransportPaths<'_>,
    duration_s: f64,
    transport_controls: FiniteGasDgTransportSolve,
    acoustic_controls: AcousticMotionSolve,
    native: N,
    marched: Vec<DgPortLoad>,
    mut motion: AcousticMotionStep,
    original_total: f64,
    cancelled: &mut C,
) -> Result<FiniteGasAcousticCandidate<N>, FiniteGasAcousticError<M>> {
    poll(cancelled)?;
    if faces.is_empty()
        || faces.len() != marched.len()
        || faces.len() != motion.faces.len()
        || !original_total.is_finite()
        || !acoustic_controls.energy_tolerance_j.is_finite()
        || acoustic_controls.energy_tolerance_j <= 0.0
        || !transport_controls.energy_tolerance_j.is_finite()
        || transport_controls.energy_tolerance_j <= 0.0
    {
        return Err(FiniteGasAcousticError::Invalid(
            "complete original shared-motion endpoints and finite positive energy tolerances required",
        ));
    }
    let original_acoustic_residual = motion.acoustic_energy_residual_j;
    let mut staged = BTreeMap::new();
    let mut native_inputs = inputs.clone();
    let mut face_parcels = Vec::with_capacity(faces.len());
    for ((face, load), report) in faces.iter().zip(&marched).zip(&motion.faces) {
        poll(cancelled)?;
        if report.signed_area_m2.to_bits() != face.signed_area_m2.to_bits()
            || !report.inlet_flow_m3_s.is_finite()
        {
            return Err(FiniteGasAcousticError::Invalid(
                "regional binding disagrees with its original accepted acoustic face",
            ));
        }
        let input = native_inputs
            .get_mut(&face.parcel)
            .ok_or(FiniteGasAcousticError::Invalid(
                "the original acoustic face lost its finite parcel input",
            ))?;
        if !input.ports.is_empty() {
            return Err(FiniteGasAcousticError::Invalid(
                "shared-motion inlet flow is supplied only by the accepted native motion",
            ));
        }
        input.ports.push(PortInput {
            id: face.load.inlet(),
            effort: report.inlet_flow_m3_s,
        });
        if staged
            .insert(
                face.parcel,
                StagedPort {
                    original: face.load,
                    marched: load,
                    output: &report.duct,
                    endpoint: None,
                },
            )
            .is_some()
        {
            return Err(FiniteGasAcousticError::Invalid(
                "one regional inventory cannot own two shared-motion acoustic loads",
            ));
        }
        face_parcels.push(face.parcel);
    }
    let candidate = finite_gas_transport_candidate_with_paths(
        parcels,
        &staged,
        &native_inputs,
        paths,
        duration_s,
        transport_controls,
        cancelled,
    )?;
    let (next_parcels, next_nodes, transport) = candidate.into_parts();
    let mut loads = BTreeMap::new();
    for (id, node) in next_nodes {
        poll(cancelled)?;
        let load = node.endpoint.ok_or(FiniteGasAcousticError::Invalid(
            "native finite-gas transport omitted a shared-motion acoustic endpoint",
        ))?;
        loads.insert(id, load);
    }
    // End all staging borrows before refreshing the actual returned motion rows.
    drop(staged);
    for (report, &id) in motion.faces.iter_mut().zip(&face_parcels) {
        poll(cancelled)?;
        let endpoint = transport
            .acoustics
            .get(&id)
            .ok_or(FiniteGasAcousticError::Invalid(
                "native finite-gas transport omitted an original acoustic report",
            ))?;
        // Preserve measured interval force, flow, pressure, inlet work and poles.
        // The native finite-gas owner has already refreshed only endpoint data.
        report.duct = endpoint.acoustic.clone();
    }
    motion.acoustic_storage_change_j = transport.acoustic_storage_change_j;
    motion.acoustic_dissipated_j = transport.acoustic_dissipation_j;
    motion.parameter_work_j = transport.parameter_work_j;
    let inlet_work = motion
        .faces
        .iter()
        .fold(0.0, |sum, face| sum + face.inlet_work_j);
    motion.acoustic_energy_residual_j = (((motion.acoustic_storage_change_j
        + motion.acoustic_dissipated_j)
        + motion.terminal_work_j)
        - inlet_work)
        - motion.parameter_work_j;

    // The transport account already contains the complete original propagation,
    // all regional gas storage and native parameter work. Replace the old frozen
    // acoustic residual; retain actual native reciprocal mechanical work once.
    let total =
        (original_total - original_acoustic_residual) + transport.combined_energy_residual_j;
    let tolerance = acoustic_controls
        .energy_tolerance_j
        .min(transport_controls.energy_tolerance_j);
    if !total.is_finite()
        || total.abs() > tolerance
        || !motion.acoustic_energy_residual_j.is_finite()
        || motion.acoustic_energy_residual_j.abs() > tolerance
    {
        return Err(FiniteGasAcousticError::Invalid(
            "native shared motion, finite regional gases and complete acoustic energy do not close",
        ));
    }
    poll(cancelled)?;
    Ok(FiniteGasAcousticCandidate {
        native,
        parcels: next_parcels,
        loads,
        step: FiniteGasAcousticStep {
            motion,
            transport,
            face_parcels,
            total_energy_residual_j: total,
        },
    })
}

impl DoubleActingGasPiston {
    /// Advance the original two-gas motion with additional finite acoustic regions.
    ///
    /// Every authored face has one distinct parcel and original complete load.
    /// Other parcels are finite reservoirs. Caller port lists are empty: the
    /// accepted native displacement supplies each acoustic inlet flow. Material
    /// transport uses only the authored native paths and never infers mass flow
    /// from acoustic swept flow or a density equality between distinct regions.
    ///
    /// The original impulse solve and each acoustic load advance once with the
    /// initial gas sources. Native transport and caloric work feedback update the
    /// endpoint medium for the following interval. This first-order partition
    /// requires time refinement; it is not an implicit mean-flow acoustic junction.
    ///
    /// Regions retain fixed geometry, a uniform main-bore gas and the original
    /// finite-load source/volume limits. They are separate from both piston
    /// chambers. Mean-pressure force requires a balanced reference or authored
    /// compensation. Gas enthalpy already includes pressure-flow work; no jet,
    /// acoustic-loss heat or second pressure-work source is introduced here.
    ///
    /// Consume and publish all returned native, gas and acoustic endpoints
    /// together. Refusal or cancellation leaves every supplied owner unchanged.
    #[allow(clippy::too_many_arguments)]
    pub fn transported_acoustic_candidate<C: FnMut() -> bool>(
        &self,
        faces: &[FiniteGasAcousticFace<'_>],
        parcels: &BTreeMap<u64, ClosedGasParcel>,
        inputs: &BTreeMap<u64, FiniteGasDgTransportInput>,
        paths: GasTransportPaths<'_>,
        input: DoublePistonInput,
        piston_controls: DoublePistonSolve,
        acoustic_controls: AcousticMotionSolve,
        transport_controls: FiniteGasDgTransportSolve,
        cancelled: &mut C,
    ) -> Result<
        FiniteGasAcousticCandidate<(DoubleActingGasPiston, DoublePistonStep)>,
        FiniteGasAcousticError<PistonError>,
    > {
        let native_faces = prepare_faces(faces, parcels, inputs, acoustic_controls, cancelled)?;
        let candidate = self.acoustic_candidate(
            input,
            piston_controls,
            &native_faces,
            acoustic_controls,
            &mut *cancelled,
        )?;
        let (piston, loads, step, motion, total) = candidate.into_parts();
        let joint_controls = AcousticMotionSolve {
            energy_tolerance_j: acoustic_controls
                .energy_tolerance_j
                .min(piston_controls.piston.energy_tolerance_j),
            ..acoustic_controls
        };
        finish(
            faces,
            parcels,
            inputs,
            paths,
            input.duration_s,
            transport_controls,
            joint_controls,
            (piston, step),
            loads,
            motion,
            total,
            cancelled,
        )
    }
}

#[cfg(test)]
mod tests;
