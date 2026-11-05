//! Local error control over the ORIGINAL complete diffusion/chemistry/feedback
//! split. A coarse physical candidate is compared with two fine candidates;
//! only the fine original owners and their source receipts can be accepted.
use super::*;
use crate::cavity_phs::piston::{PistonError, reactive::GasReactionError};
use crate::cavity_phs::piston::thermal::ThermalPistonError;
use crate::cavity_phs::piston::thermal::boundary::ports::network::NetworkError;
use crate::circuit_enthalpy::CircuitEnthalpyError;
use fs_material::thermochemical_gas::energy::parcel::transfer::species::network::{SpeciesNetworkError, SpeciesNetworkSolve};
use fs_thermochem::SpeciesId;

mod error;
mod balance;
mod physical;
#[cfg(feature = "molecular-transport")]
pub(crate) use error::scaled_error;
pub use physical::PhysicalSubstep;
pub(crate) use physical::{PhysicalControlError, PhysicalStepError};
#[cfg(feature = "molecular-transport")]
pub(crate) use physical::PhysicalEstimate;
pub use balance::{CoupledDiffusionBalance, CoupledDiffusionChamberBalance,
    CoupledDiffusionSpeciesBalance, CoupledReactionExtent};

/// Positive absolute LOCAL error scales, separately dimensioned for every
/// independent physical coordinate. Thermal energy is per cell, not a network
/// total: an arbitrary caloric datum and opposite cell errors cannot hide it.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledDiffusionScales {
    /// Each chamber/species inventory [kg], including arriving/consumed species.
    pub species_mass_kg: f64,
    pub gas_temperature_k: f64,
    pub gas_pressure_pa: f64,
    /// Gas volume also controls piston position through its unchanged area [m3].
    pub gas_volume_m3: f64,
    pub piston_velocity_m_s: f64,
    /// Each fixed thermal cell's extensive enthalpy increment [J].
    pub cell_energy_j: f64,
    /// Every native capacitor charge [C], including coupled/nonlinear storage.
    pub charge_c: f64,
    /// Every native inductor flux linkage [Wb], including mutual inductance.
    pub flux_wb: f64,
}

/// Complete coupled time accuracy and work, separate from the native physical
/// residual budgets. No material/source/withdrawal tolerance is loosened.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledDiffusionControls {
    pub absolute: CoupledDiffusionScales,
    /// Add rtol*max(|coarse-start|, |fine-start|) to each absolute allowance.
    /// Finite and in [0,1); absolute state/reference magnitudes are never used.
    pub relative_tolerance: f64,
    /// Largest coarse comparison interval [s]; accepted physical steps are halves.
    pub maximum_step_s: f64,
    /// Smallest retry comparison interval [s]; a final clipped interval may be shorter.
    pub minimum_step_s: f64,
    /// ALL whole-assembly attempts, including failed, coarse and rejected fine
    /// probes. Native controls bound inner work per attempt. Also capped by the
    /// original thermal owner's max_substeps over this entire call.
    pub maximum_evaluations: usize,
}

/// Dimensionless coarse-versus-fine discrepancies; every entry must be <= 1.
/// For equal halves the first-order fine local estimate is their difference;
/// unequal representable durations retain their exact truncation factor.
/// There is no Richardson extrapolation or claimed global error bound.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct CoupledDiffusionEstimate {
    pub species_mass: f64,
    pub gas_temperature: f64,
    pub gas_pressure: f64,
    pub gas_volume: f64,
    pub piston_velocity: f64,
    pub cell_energy: f64,
    pub charge: f64,
    pub flux: f64,
}
impl CoupledDiffusionEstimate {
    pub fn maximum(self) -> f64 {
        [self.species_mass, self.gas_temperature, self.gas_pressure, self.gas_volume,
            self.piston_velocity, self.cell_energy, self.charge, self.flux]
            .into_iter().fold(0.0, f64::max)
    }
    pub(crate) fn accumulate(&mut self, other: Self) {
        self.species_mass = self.species_mass.max(other.species_mass);
        self.gas_temperature = self.gas_temperature.max(other.gas_temperature);
        self.gas_pressure = self.gas_pressure.max(other.gas_pressure);
        self.gas_volume = self.gas_volume.max(other.gas_volume);
        self.piston_velocity = self.piston_velocity.max(other.piston_velocity);
        self.cell_energy = self.cell_energy.max(other.cell_energy);
        self.charge = self.charge.max(other.charge);
        self.flux = self.flux.max(other.flux);
    }
}

/// Complete-assembly refusal retaining the local-accuracy payload and the
/// original physical/constitutive error, including a selected heat-law source.
#[derive(Debug)]
pub enum CoupledPhysicalError<E, S = MaterialFlowError> {
    Invalid(&'static str),
    Step { start_s: f64, end_s: f64, error: S },
    Accuracy { start_s: f64, step_s: f64, estimate: E },
    WorkLimit { time_s: f64, evaluations: usize, accepted_substeps: usize },
    Balance(&'static str),
    Circuit(crate::material_circuit::MaterialCircuitError),
    Gas { chamber: u64, error: fs_material::MaterialError },
    Cancelled,
}
impl<E: core::fmt::Debug, S: core::fmt::Debug> core::fmt::Display for CoupledPhysicalError<E, S> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled molecular diffusion: {self:?}")
    }
}
impl<E: core::fmt::Debug, S: core::fmt::Debug> std::error::Error for CoupledPhysicalError<E, S> {}
/// Original concrete closed-assembly error, preserving its variant namespace.
#[derive(Debug)]
pub enum CoupledDiffusionError {
    Invalid(&'static str),
    Step { start_s: f64, end_s: f64, error: MaterialFlowError },
    Accuracy { start_s: f64, step_s: f64, estimate: CoupledDiffusionEstimate },
    WorkLimit { time_s: f64, evaluations: usize, accepted_substeps: usize },
    Balance(&'static str),
    Circuit(crate::material_circuit::MaterialCircuitError),
    Gas { chamber: u64, error: fs_material::MaterialError },
    Cancelled,
}
impl core::fmt::Display for CoupledDiffusionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled molecular diffusion: {self:?}")
    }
}
impl std::error::Error for CoupledDiffusionError {}
impl From<CoupledPhysicalError<CoupledDiffusionEstimate>> for CoupledDiffusionError {
    fn from(error: CoupledPhysicalError<CoupledDiffusionEstimate>) -> Self {
        match error {
            CoupledPhysicalError::Invalid(reason) => Self::Invalid(reason),
            CoupledPhysicalError::Step { start_s, end_s, error } => Self::Step { start_s, end_s, error },
            CoupledPhysicalError::Accuracy { start_s, step_s, estimate } => Self::Accuracy { start_s, step_s, estimate },
            CoupledPhysicalError::WorkLimit { time_s, evaluations, accepted_substeps } =>
                Self::WorkLimit { time_s, evaluations, accepted_substeps },
            CoupledPhysicalError::Balance(reason) => Self::Balance(reason),
            CoupledPhysicalError::Circuit(error) => Self::Circuit(error),
            CoupledPhysicalError::Gas { chamber, error } => Self::Gas { chamber, error },
            CoupledPhysicalError::Cancelled => Self::Cancelled,
        }
    }
}
type Error = CoupledDiffusionError;

/// One actual accepted half step, retaining every original source and balance.
pub type CoupledDiffusionSubstep<F = DiffusionNetworkStep, I = ResolvedMaterialStatePoint> =
    PhysicalSubstep<MaterialFlowAdvance<F, I>>;
pub struct CoupledDiffusionReport<F = DiffusionNetworkStep, I = ResolvedMaterialStatePoint> {
    pub start_s: f64,
    pub end_s: f64,
    pub evaluations: usize,
    pub accepted_intervals: usize,
    pub rejected_intervals: usize,
    pub minimum_accepted_substep_s: f64,
    pub maximum_accepted_substep_s: f64,
    /// Componentwise maximum over ACCEPTED coarse/fine comparisons only.
    pub estimated_error: CoupledDiffusionEstimate,
    /// Independent original-to-final accounts against accepted physical ledgers.
    pub balance: CoupledDiffusionBalance,
    /// Two fine steps per accepted comparison. Rejected receipts are discarded.
    pub steps: Vec<CoupledDiffusionSubstep<F, I>>,
}

// A private ingress for other original transport candidates. They must expose
// the SAME conservative species owner; the time loop and complete assembly
// error/balance logic are never duplicated for a different path discretization.
pub(crate) trait TransportLedger {
    fn species_exchange(&self) -> &fs_material::thermochemical_gas::energy::parcel::transfer::species::network::SpeciesNetworkStep;
    /// Additional typed timestep-sensitive failures owned by this transport.
    fn refinable(_error: &MaterialFlowError) -> bool { false }
}
impl TransportLedger for DiffusionNetworkStep {
    fn species_exchange(&self) -> &fs_material::thermochemical_gas::energy::parcel::transfer::species::network::SpeciesNetworkStep {
        &self.exchange
    }
}

fn check(cx: &Cx<'_>) -> Result<(), PhysicalControlError> {
    cx.checkpoint().map_err(|_| PhysicalControlError::Cancelled)
}
fn finite(value: f64) -> Result<f64, PhysicalControlError> {
    if value.is_finite() { Ok(value) }
    else { Err(PhysicalControlError::Invalid("coupled diffusion arithmetic is unrepresentable")) }
}
impl CoupledDiffusionControls {
    fn validate(self) -> Result<(), PhysicalControlError> {
        let a = self.absolute;
        if [a.species_mass_kg, a.gas_temperature_k, a.gas_pressure_pa,
            a.gas_volume_m3, a.piston_velocity_m_s, a.cell_energy_j, a.charge_c,
            a.flux_wb, self.maximum_step_s, self.minimum_step_s].into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || self.minimum_step_s > self.maximum_step_s || self.maximum_evaluations < 3
        {
            return Err(PhysicalControlError::Invalid("positive dimensioned error scales, bounded steps and at least three evaluations required"));
        }
        Ok(())
    }
}

impl<'m> MaterialExchangeState<'m> {
    /// Adapt the WHOLE diffusion/reaction/heat/electrical/piston split locally.
    /// Each comparison runs one complete coarse step and two complete fine
    /// steps from the same accepted start. Only the two fine states are used;
    /// no extrapolated inventories, heat, charge or momentum are manufactured.
    ///
    /// Every fine step samples current gas volumes, composition, T/p and original
    /// material laws. Arrivals, reaction heat, wall heat and piston compression
    /// therefore feed the next diffusion solve. Connector areas/distances and
    /// thermal topology remain authored/fixed; moving connectors are not inferred.
    /// Both native molecular force/reference frames are explicitly supported.
    ///
    /// Electrical inputs, exterior pressures and forces are held over this call.
    /// Each drive's integrated heat is apportioned at constant power using
    /// cumulative endpoints; guess_j is scaled as a numerical energy guess only.
    /// Pending events at or before target refuse. Split such histories explicitly.
    /// Source/domain/admission failures are fatal; only typed withdrawal, finite
    /// entropy, consumption and nonlinear step-size failures can retry.
    ///
    /// This controls a local estimate of the existing first-order split, including
    /// every physical storage coordinate. It is not Strang splitting, an implicit
    /// diffusion solver, a global error certificate or a resolved flow model.
    #[allow(clippy::too_many_arguments)]
    pub fn coupled_diffusion_candidate(
        &self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve, driving: GasDiffusionDriving,
        controls: NetworkSolve, guess_j: &[f64], time: CoupledDiffusionControls,
    ) -> Result<(Self, CoupledDiffusionReport), Error> {
        self.coupled_transport_window_candidate(cx, target_s, ports, drives, reactions,
            diffusion_controls.exchange, controls, guess_j, time, |state, end, drives, guess| {
                state.diffusion_step_with_driving_candidate(prepared, cx, end, ports, drives,
                    reactions, gas_patches, links, diffusion_controls, driving, controls, guess)
            })
    }

    #[allow(clippy::too_many_arguments)]
    pub(crate) fn coupled_transport_window_candidate<F: TransportLedger>(
        &self, cx: &Cx<'_>, target_s: f64, ports: &[f64], drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>], species_controls: SpeciesNetworkSolve,
        controls: NetworkSolve, guess_j: &[f64], time: CoupledDiffusionControls,
        evolve: impl FnMut(&Self, f64, &[ChamberDrive], &[f64])
            -> Result<(Self, MaterialFlowAdvance<F>), MaterialFlowError>,
    ) -> Result<(Self, CoupledDiffusionReport<F>), Error> {
        self.coupled_transport_window_with_error_candidate(cx, target_s, ports, drives,
            reactions, species_controls, controls, guess_j, time, evolve, F::refinable)
            .map_err(Into::into)
    }

    /// The same closed-window accounts retain the selected heat law's original
    /// source receipt and typed refusal; only the complete trial ingress varies.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn coupled_transport_window_with_error_candidate<F: TransportLedger, I, S: PhysicalStepError>(
        &self, cx: &Cx<'_>, target_s: f64, ports: &[f64], drives: &[ChamberDrive],
        reactions: &[ChamberReaction<'_>], species_controls: SpeciesNetworkSolve,
        controls: NetworkSolve, guess_j: &[f64], time: CoupledDiffusionControls,
        evolve: impl FnMut(&Self, f64, &[ChamberDrive], &[f64])
            -> Result<(Self, MaterialFlowAdvance<F, I>), S>,
        transport_refinable: impl Fn(&S) -> bool,
    ) -> Result<(Self, CoupledDiffusionReport<F, I>), CoupledPhysicalError<CoupledDiffusionEstimate, S>> {
        let (mut next, window) = self.coupled_physical_window_candidate(
            cx, target_s, drives, guess_j, time, evolve,
            |state_error, _, _, _, _| Ok(state_error), transport_refinable,
        )?;
        let balance = balance::summarize(self, &next, cx, ports, drives, reactions,
            species_controls, controls, &window.steps)?;
        balance::macro_ports(self, &mut next, cx, controls, &balance)?;
        check(cx)?;
        Ok((next, CoupledDiffusionReport {
            start_s: window.start_s, end_s: window.end_s, evaluations: window.evaluations,
            accepted_intervals: window.accepted_intervals, rejected_intervals: window.rejected_intervals,
            minimum_accepted_substep_s: window.minimum_accepted_substep_s,
            maximum_accepted_substep_s: window.maximum_accepted_substep_s,
            estimated_error: window.estimated_error, balance, steps: window.steps,
        }))
    }

    /// Publish only after complete-window accuracy, species/energy/impulse
    /// accounts and the final cancellation gate. Failure changes no live owner.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_coupled_diffusion(
        &mut self, prepared: &PreparedMaterialExchangeNetwork<'_, '_>, cx: &Cx<'_>, target_s: f64,
        ports: &[f64], drives: &[ChamberDrive], reactions: &[ChamberReaction<'_>],
        gas_patches: &[MaterialContact], links: &[DiffusionNetworkLink<'_>],
        diffusion_controls: DiffusionNetworkSolve, driving: GasDiffusionDriving,
        controls: NetworkSolve, guess_j: &[f64], time: CoupledDiffusionControls,
    ) -> Result<CoupledDiffusionReport, Error> {
        let (next, report) = self.coupled_diffusion_candidate(prepared, cx, target_s, ports,
            drives, reactions, gas_patches, links, diffusion_controls, driving, controls, guess_j, time)?;
        check(cx)?;
        *self = next;
        Ok(report)
    }
}

// Only typed timestep-sensitive failures retry. Never inspect error strings or
// reinterpret a material/source-domain error as a request for smaller steps.
fn refinable(error: &MaterialFlowError) -> bool {
    use transport::DiffusionNetworkError as D;
    matches!(error,
        MaterialFlowError::Diffusion(D::EntropyDecrease { .. }
            | D::Exchange(SpeciesNetworkError::Withdrawal { .. })))
        || matches!(error, MaterialFlowError::Exchange(error) if feedback_refinable(error))
}

// A caller's source-error type never affects which original physical failures
// are timestep-sensitive. Constitutive/source refusals remain fatal.
fn feedback_refinable<E>(error: &ReactiveExchangeError<E>) -> bool {
    matches!(error,
        ReactiveExchangeError::Reaction {
            error: GasReactionError::ConsumptionLimit { .. } | GasReactionError::EntropyDecrease { .. }
                | GasReactionError::ImplicitOutsideBracket { .. } | GasReactionError::ImplicitNotConverged { .. }, .. }
        | ReactiveExchangeError::Exchange(
            CircuitEnthalpyError::NotConverged(_)
            | CircuitEnthalpyError::Source(NetworkError::Boundary {
                error: ThermalPistonError::OutsideBracket { .. } | ThermalPistonError::NotConverged { .. }
                    | ThermalPistonError::Piston(PistonError::OutsideBracket { .. }
                        | PistonError::NotConverged { .. }
                        | PistonError::TimeResolution(fs_phs::linear_exact::ResolvedBandError::Outside { .. })), .. })))
}

fn fixed(initial: &MaterialExchangeState<'_>, next: &MaterialExchangeState<'_>) -> Result<(), PhysicalControlError> {
    if !initial.chambers.keys().eq(next.chambers.keys())
        || !initial.thermal.cells().keys().eq(next.thermal.cells().keys())
        || initial.thermal.links() != next.thermal.links()
        || initial.circuit.charges_c().len() != next.circuit.charges_c().len()
        || initial.circuit.fluxes_wb().len() != next.circuit.fluxes_wb().len()
        || initial.thermal.cells().iter().any(|(id, a)| {
            let b = &next.thermal.cells()[id];
            a.mass_kg != b.mass_kg || a.material_id != b.material_id || a.geometry != b.geometry
                || a.boundary != b.boundary || a.body_id != b.body_id || a.segment_id != b.segment_id
        })
    {
        return Err(PhysicalControlError::Invalid("coupled diffusion requires fixed cell/circuit/chamber topology and cell masses"));
    }
    Ok(())
}

fn inventory(gas: &ClosedGasParcel) -> Result<BTreeMap<SpeciesId, f64>, PhysicalControlError> {
    let total_mol = finite(gas.mass_kg() / gas.state().thermodynamics().properties().molar_mass().value())?;
    gas.model().thermodynamics().components().iter()
        .zip(gas.model().thermodynamics().mole_composition().fractions())
        .map(|(species, &fraction)| {
            let mass = finite((total_mol * fraction) * species.molar_mass().value())?;
            if mass < 0.0 { return Err(PhysicalControlError::Invalid("negative actual species inventory")); }
            Ok((species.species().clone(), mass))
        }).collect()
}

fn compare(start: &MaterialExchangeState<'_>, coarse: &MaterialExchangeState<'_>,
    fine: &MaterialExchangeState<'_>, policy: CoupledDiffusionControls, fine_error_factor: f64, cx: &Cx<'_>,
) -> Result<CoupledDiffusionEstimate, PhysicalControlError> {
    let mut result = CoupledDiffusionEstimate::default();
    let observe = |target: &mut f64, a: f64, b: f64, absolute: f64| -> Result<(), PhysicalControlError> {
        let ratio = error::scaled_error(a, b, absolute, policy.relative_tolerance)
            .ok_or(PhysicalControlError::Invalid("coupled local error is unrepresentable"))?;
        // Overflow remains an infinite rejecting estimate, never a zero error.
        *target = target.max(ratio * fine_error_factor);
        Ok(())
    };
    let p = policy.absolute;
    for (&id, before) in &start.chambers {
        check(cx)?;
        let a = &coarse.chambers[&id]; let b = &fine.chambers[&id];
        let old = inventory(before.gas())?; let left = inventory(a.gas())?; let right = inventory(b.gas())?;
        let ids: BTreeSet<_> = old.keys().chain(left.keys()).chain(right.keys()).collect();
        for id in ids {
            let initial = old.get(id).copied().unwrap_or(0.0);
            observe(&mut result.species_mass, left.get(id).copied().unwrap_or(0.0) - initial,
                right.get(id).copied().unwrap_or(0.0) - initial, p.species_mass_kg)?;
        }
        let old = before.gas().state().state(); let left = a.gas().state().state(); let right = b.gas().state().state();
        observe(&mut result.gas_temperature, left.temperature - old.temperature, right.temperature - old.temperature, p.gas_temperature_k)?;
        observe(&mut result.gas_pressure, left.pressure - old.pressure, right.pressure - old.pressure, p.gas_pressure_pa)?;
        observe(&mut result.gas_volume, a.gas().volume_m3() - before.gas().volume_m3(), b.gas().volume_m3() - before.gas().volume_m3(), p.gas_volume_m3)?;
        observe(&mut result.piston_velocity, a.velocity_m_s() - before.velocity_m_s(), b.velocity_m_s() - before.velocity_m_s(), p.piston_velocity_m_s)?;
    }
    for (id, before) in start.thermal.cells() {
        check(cx)?;
        observe(&mut result.cell_energy, before.mass_kg * (coarse.thermal.cells()[id].specific_h_j_kg - before.specific_h_j_kg),
            before.mass_kg * (fine.thermal.cells()[id].specific_h_j_kg - before.specific_h_j_kg), p.cell_energy_j)?;
    }
    for ((&old, &a), &b) in start.circuit.charges_c().iter().zip(coarse.circuit.charges_c()).zip(fine.circuit.charges_c()) {
        check(cx)?; observe(&mut result.charge, a - old, b - old, p.charge_c)?;
    }
    for ((&old, &a), &b) in start.circuit.fluxes_wb().iter().zip(coarse.circuit.fluxes_wb()).zip(fine.circuit.fluxes_wb()) {
        check(cx)?; observe(&mut result.flux, a - old, b - old, p.flux_wb)?;
    }
    Ok(result)
}
