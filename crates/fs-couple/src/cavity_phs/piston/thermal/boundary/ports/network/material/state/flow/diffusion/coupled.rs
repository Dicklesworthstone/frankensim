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
    fn accumulate(&mut self, other: Self) {
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
type Error = CoupledDiffusionError;

/// One actual accepted half step, retaining every original source and balance.
pub struct CoupledDiffusionSubstep<F = DiffusionNetworkStep> {
    pub start_s: f64,
    pub end_s: f64,
    pub advance: MaterialFlowAdvance<F>,
}
pub struct CoupledDiffusionReport<F = DiffusionNetworkStep> {
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
    pub steps: Vec<CoupledDiffusionSubstep<F>>,
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

fn check(cx: &Cx<'_>) -> Result<(), Error> {
    cx.checkpoint().map_err(|_| Error::Cancelled)
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() { Ok(value) }
    else { Err(Error::Invalid("coupled diffusion arithmetic is unrepresentable")) }
}
impl CoupledDiffusionControls {
    fn validate(self) -> Result<(), Error> {
        let a = self.absolute;
        if [a.species_mass_kg, a.gas_temperature_k, a.gas_pressure_pa,
            a.gas_volume_m3, a.piston_velocity_m_s, a.cell_energy_j, a.charge_c,
            a.flux_wb, self.maximum_step_s, self.minimum_step_s].into_iter()
            .any(|x| !x.is_finite() || x <= 0.0)
            || !self.relative_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_tolerance)
            || self.minimum_step_s > self.maximum_step_s || self.maximum_evaluations < 3
        {
            return Err(Error::Invalid("positive dimensioned error scales, bounded steps and at least three evaluations required"));
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
        mut evolve: impl FnMut(&Self, f64, &[ChamberDrive], &[f64])
            -> Result<(Self, MaterialFlowAdvance<F>), MaterialFlowError>,
    ) -> Result<(Self, CoupledDiffusionReport<F>), Error> {
        check(cx)?;
        time.validate()?;
        let start_s = self.time_s();
        let duration = finite(target_s - start_s)?;
        if duration <= 0.0 || guess_j.len() != self.thermal.cells().len()
            || guess_j.iter().any(|x| !x.is_finite())
            || drives.iter().any(|d| !d.heat_j.is_finite())
            || self.thermal.pending_events().first().is_some_and(|e| e.time_s <= target_s)
        {
            return Err(Error::Invalid("one forward fixed-assembly window, finite inputs and no pending events through target required"));
        }
        let cap = time.maximum_evaluations.min(self.thermal.work_budget().max_substeps);
        let mut next = Self::from_parts(self.thermal.clone(), self.circuit.clone(), self.chambers.clone());
        let mut steps = Vec::new();
        let (mut evaluations, mut accepted_intervals, mut rejected_intervals) = (0, 0, 0);
        let (mut minimum, mut maximum) = (f64::INFINITY, 0.0_f64);
        let mut estimated_error = CoupledDiffusionEstimate::default();
        let mut step_s = time.maximum_step_s.min(duration);
        while next.time_s() < target_s {
            check(cx)?;
            let before = next.time_s();
            let remaining = finite(target_s - before)?;
            let dt = step_s.min(remaining);
            let end = if dt == remaining { target_s } else { finite(before + dt)? };
            let dt = finite(end - before)?;
            let midpoint = finite(before + 0.5 * (end - before))?;
            if midpoint <= before || midpoint >= end {
                return Err(Error::Invalid("coupled diffusion half-step clock is unrepresentable"));
            }
            let fine_error_factor = error::fine_error_factor(midpoint - before, end - midpoint)
                .ok_or(Error::Invalid("coupled diffusion half-step ratio is unrepresentable"))?;
            let accepted_substeps = steps.len();
            let mut evaluate = |state: &Self, end: f64, evaluations: &mut usize| {
                check(cx)?;
                if *evaluations >= cap {
                    return Err(Error::WorkLimit { time_s: before, evaluations: *evaluations, accepted_substeps });
                }
                let left = (state.time_s() - start_s) / duration;
                let right = if end == target_s { 1.0 } else { (end - start_s) / duration };
                let subdrives: Vec<_> = drives.iter().map(|drive| ChamberDrive {
                    heat_j: drive.heat_j * right - drive.heat_j * left, ..*drive
                }).collect();
                let fraction = (end - state.time_s()) / duration;
                let guess: Vec<_> = guess_j.iter().map(|value| value * fraction).collect();
                *evaluations += 1;
                let (candidate, report) = evolve(state, end, &subdrives, &guess).map_err(|error| match error {
                    MaterialFlowError::Cancelled => Error::Cancelled,
                    error => Error::Step { start_s: state.time_s(), end_s: end, error },
                })?;
                fixed(self, &candidate)?;
                Ok((candidate, report))
            };
            let comparison = (|| {
                let (coarse, _) = evaluate(&next, end, &mut evaluations)?;
                let (half, first) = evaluate(&next, midpoint, &mut evaluations)?;
                let (fine, second) = evaluate(&half, end, &mut evaluations)?;
                let estimate = compare(&next, &coarse, &fine, time, fine_error_factor, cx)?;
                Ok((fine, first, second, estimate))
            })();
            let reduction = match comparison {
                Ok((fine, first, second, estimate)) if estimate.maximum() <= 1.0 => {
                    estimated_error.accumulate(estimate);
                    accepted_intervals += 1;
                    minimum = minimum.min(midpoint - before).min(end - midpoint);
                    maximum = maximum.max(midpoint - before).max(end - midpoint);
                    steps.push(CoupledDiffusionSubstep { start_s: before, end_s: midpoint, advance: first });
                    steps.push(CoupledDiffusionSubstep { start_s: midpoint, end_s: end, advance: second });
                    next = fine;
                    let factor = error::step_factor(estimate.maximum(), true);
                    step_s = error::bounded_step_s(dt, factor, time.minimum_step_s, time.maximum_step_s);
                    continue;
                }
                Ok((_, _, _, estimate)) => {
                    let failure = Error::Accuracy { start_s: before, step_s: dt, estimate };
                    (error::step_factor(estimate.maximum(), false), failure)
                }
                Err(Error::Step { start_s, end_s, error }) if refinable(&error) || F::refinable(&error) =>
                    (0.5, Error::Step { start_s, end_s, error }),
                Err(error) => return Err(error),
            };
            rejected_intervals += 1;
            if dt <= time.minimum_step_s {
                return Err(reduction.1);
            }
            let smaller = error::bounded_step_s(dt, reduction.0, time.minimum_step_s, dt);
            if !smaller.is_finite() || smaller >= dt || before + smaller >= end {
                return Err(reduction.1);
            }
            step_s = smaller;
        }
        let balance = balance::summarize(self, &next, cx, ports, drives, reactions,
            species_controls, controls, &steps)?;
        balance::macro_ports(self, &mut next, cx, controls, &balance)?;
        check(cx)?;
        Ok((next, CoupledDiffusionReport { start_s, end_s: target_s, evaluations,
            accepted_intervals, rejected_intervals, minimum_accepted_substep_s: minimum,
            maximum_accepted_substep_s: maximum, estimated_error, balance, steps }))
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
            | D::Exchange(SpeciesNetworkError::Withdrawal { .. }))
        | MaterialFlowError::Exchange(ReactiveExchangeError::Reaction {
            error: GasReactionError::ConsumptionLimit { .. } | GasReactionError::EntropyDecrease { .. }
                | GasReactionError::ImplicitOutsideBracket { .. } | GasReactionError::ImplicitNotConverged { .. }, .. })
        | MaterialFlowError::Exchange(ReactiveExchangeError::Exchange(
            CircuitEnthalpyError::NotConverged(_)
            | CircuitEnthalpyError::Source(NetworkError::Boundary {
                error: ThermalPistonError::OutsideBracket { .. } | ThermalPistonError::NotConverged { .. }
                    | ThermalPistonError::Piston(PistonError::OutsideBracket { .. }
                        | PistonError::NotConverged { .. }
                        | PistonError::TimeResolution(fs_phs::linear_exact::ResolvedBandError::Outside { .. })), .. }))))
}

fn fixed(initial: &MaterialExchangeState<'_>, next: &MaterialExchangeState<'_>) -> Result<(), Error> {
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
        return Err(Error::Invalid("coupled diffusion requires fixed cell/circuit/chamber topology and cell masses"));
    }
    Ok(())
}

fn inventory(gas: &ClosedGasParcel) -> Result<BTreeMap<SpeciesId, f64>, Error> {
    let total_mol = finite(gas.mass_kg() / gas.state().thermodynamics().properties().molar_mass().value())?;
    gas.model().thermodynamics().components().iter()
        .zip(gas.model().thermodynamics().mole_composition().fractions())
        .map(|(species, &fraction)| {
            let mass = finite((total_mol * fraction) * species.molar_mass().value())?;
            if mass < 0.0 { return Err(Error::Invalid("negative actual species inventory")); }
            Ok((species.species().clone(), mass))
        }).collect()
}

fn compare(start: &MaterialExchangeState<'_>, coarse: &MaterialExchangeState<'_>,
    fine: &MaterialExchangeState<'_>, policy: CoupledDiffusionControls, fine_error_factor: f64, cx: &Cx<'_>,
) -> Result<CoupledDiffusionEstimate, Error> {
    let mut result = CoupledDiffusionEstimate::default();
    let observe = |target: &mut f64, a: f64, b: f64, absolute: f64| -> Result<(), Error> {
        let ratio = error::scaled_error(a, b, absolute, policy.relative_tolerance)
            .ok_or(Error::Invalid("coupled local error is unrepresentable"))?;
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
