//! Scheduled finite-reservoir exchange during continuing pressure playback.
//!
//! The original nozzle/species owner transports partial enthalpies once, and
//! the original moving-source transaction pays acoustic parameter work from
//! that same candidate. Every finite inventory and the source publish together.
//! This is a fixed-volume, well-mixed material network with explicit restrictions,
//! separate from the source's oscillatory acoustic flow and prescribed pressure.

use super::{AperturePerformance, CoupledAperture, GestureCompileError, invalid, receiver_error};
use super::finite_gas::{FiniteGasHeatProgram, GasHeatControl};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::dg::{ApertureDg, DgApertureFrame};
use crate::bernoulli_aperture::dg::finite_gas::transport::{
    FiniteGasApertureDrive, FiniteGasApertureTransportInput, FiniteGasApertureTransportStep,
    FiniteGasDgTransportError, FiniteGasDgTransportInput, FiniteGasDgTransportSolve,
    scalar::{FiniteGasScalarApertureTransportCandidate,
        finite_gas_scalar_aperture_transport_candidate_with_paths},
};
use crate::bernoulli_aperture::tube::TubeDrive;
use crate::dg_port::finite_gas::FiniteGasDgError;
use crate::render::RenderError;
use fs_exec::CancelGate;
use fs_material::thermochemical_gas::energy::{GasEnergyError, GasEnergySolve};
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;
use fs_material::thermochemical_gas::energy::parcel::transfer::{
    nozzle::{GasNozzle, network::snapshot::SnapshotNozzleLink},
    snapshot::{MAX_SNAPSHOT_PARCELS, MAX_SNAPSHOT_TRANSFERS},
    species::network::transport::{GasTransportError, GasTransportPaths},
};
use std::collections::BTreeMap;
use std::sync::Arc;

/// An actual finite inventory, with its original caloric endpoint admission.
/// Reservoirs are adiabatic here; transferred enthalpy changes their temperature.
#[derive(Clone, Debug, PartialEq)]
pub struct FiniteGasReservoir {
    pub gas: ClosedGasParcel,
    pub caloric: GasEnergySolve,
}

/// Change one existing physical restriction before its named source sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct GasNozzleControl {
    pub sample: u64,
    pub id: u64,
    /// Area [m2] and dimensionless discharge coefficient from the native owner.
    /// Zero area closes the path without resetting either finite inventory.
    pub nozzle: GasNozzle,
}

/// One source inventory and finite reservoirs connected by explicit nozzles.
#[derive(Clone, Debug)]
pub struct FiniteGasTransportProgram {
    /// Identity of the parcel covering the actual main bore and its channels.
    pub source_id: u64,
    /// Existing source-clock heat program and native caloric/work admissions.
    pub source: FiniteGasHeatProgram,
    /// Initial finite inventories; none may reuse source_id.
    pub reservoirs: BTreeMap<u64, FiniteGasReservoir>,
    /// Existing source-resolved isentropic/choked restrictions, with unique IDs.
    /// Connections may join any two declared inventories.
    pub nozzles: Vec<SnapshotNozzleLink>,
    /// Held assignments, admitted and executed in (sample, id) order.
    pub controls: Vec<GasNozzleControl>,
    /// Original gross-species, entropy and enclosing-energy admissions.
    pub transport: FiniteGasDgTransportSolve,
}

pub(super) struct FiniteGasTransportPlayback {
    program: FiniteGasTransportProgram,
    parcels: BTreeMap<u64, ClosedGasParcel>,
    inputs: BTreeMap<u64, FiniteGasApertureTransportInput>,
    held: Vec<SnapshotNozzleLink>,
    candidate: Vec<SnapshotNozzleLink>,
    gate: Arc<CancelGate>,
    next: usize,
    candidate_next: usize,
    heat_next: usize,
    held_power_w: f64,
    last: Option<FiniteGasApertureTransportStep>,
    pub(super) assignments: usize,
    pub(super) work: u64,
}

fn runtime(error: FiniteGasDgTransportError) -> RenderError {
    match error {
        FiniteGasDgTransportError::Cancelled
        | FiniteGasDgTransportError::Transport(GasTransportError::Cancelled)
        | FiniteGasDgTransportError::Gas { error: GasEnergyError::Cancelled, .. }
        | FiniteGasDgTransportError::Acoustic {
            error: FiniteGasDgError::Cancelled
                | FiniteGasDgError::Background(GasEnergyError::Cancelled), ..
        } => RenderError::Voice(AcousticRealizeError::Cancelled),
        FiniteGasDgTransportError::Acoustic {
            error: FiniteGasDgError::Acoustic(error), ..
        } => RenderError::Voice(error),
        error => receiver_error(error.to_string()),
    }
}

fn finite(nozzle: GasNozzle) -> bool {
    nozzle.area_m2.is_finite() && nozzle.area_m2 >= 0.0
        && nozzle.discharge_coefficient.is_finite()
        && (0.0..=1.0).contains(&nozzle.discharge_coefficient)
}

impl FiniteGasTransportPlayback {
    fn preview(&mut self, system: &ApertureDg, drive: TubeDrive, force_n: f64)
        -> Result<FiniteGasScalarApertureTransportCandidate, RenderError>
    {
        let sample = system.aperture().accepted_steps();
        self.candidate.copy_from_slice(&self.held);
        self.candidate_next = self.next;
        while let Some(event) = self.program.controls.get(self.candidate_next)
            .filter(|event| event.sample == sample)
        {
            let index = self.candidate.binary_search_by_key(&event.id, |link| link.id)
                .expect("admitted original restriction");
            self.candidate[index].nozzle = event.nozzle;
            self.candidate_next += 1;
        }
        let power = self.program.source.controls.get(self.heat_next)
            .filter(|event| event.sample == sample)
            .map_or(self.held_power_w, |event| event.power_w);
        let input = self.inputs.get_mut(&self.program.source_id)
            .expect("admitted source parcel");
        input.thermal.supplied_heat_j = power * system.aperture().spec().time_step_s;
        input.source = Some(FiniteGasApertureDrive {
            drive, forces_n: [force_n, 0.0], pressures_pa: [0.0; 2],
        });
        let nodes = BTreeMap::from([(self.program.source_id, system.clone())]);
        finite_gas_scalar_aperture_transport_candidate_with_paths(
            &self.parcels, &nodes, &self.inputs,
            GasTransportPaths { nozzles: &self.candidate, ..Default::default() },
            system.aperture().spec().time_step_s, self.program.transport,
            &mut || self.gate.is_requested(),
        ).map_err(runtime)
    }

    pub(super) fn step(&mut self, system: &mut ApertureDg, drive: TubeDrive, force_n: f64)
        -> Result<DgApertureFrame, RenderError>
    {
        let proposal = self.preview(system, drive, force_n)?;
        let (parcels, mut nodes, report) = proposal.into_parts();
        let frame = report.acoustics[&self.program.source_id].coupled.clone();
        // No fallible operation follows joint publication of all inventories and
        // the complete original source/duct/material/terminal histories.
        *system = nodes.remove(&self.program.source_id).expect("native source endpoint");
        self.parcels = parcels;
        self.last = Some(report);
        Ok(frame)
    }

    pub(super) fn accept(&mut self, sample: u64) {
        std::mem::swap(&mut self.held, &mut self.candidate);
        self.next = self.candidate_next;
        if let Some(event) = self.program.source.controls.get(self.heat_next)
            .filter(|event| event.sample == sample)
        {
            self.held_power_w = event.power_w;
            self.heat_next += 1;
        }
    }
}

impl AperturePerformance {
    /// Attach real finite-reservoir mass/species exchange before playback.
    ///
    /// The source uses the same volume, source, loss and observation admissions
    /// as finite-gas heat playback. Every other parcel has an explicit finite
    /// volume and native caloric solver. The original nozzle owner determines
    /// donor, choking, rate and carried enthalpy on each accepted-state snapshot;
    /// the shared transport owner checks gross withdrawal, species conservation,
    /// entropy, final driving forces and the complete energy account.
    ///
    /// Nozzle geometry and heat controls share the existing source clock and
    /// control budgets. Cancellation and any failed physical step preserve all
    /// inventories and physical owners. A failed callback remains poisoned.
    ///
    /// These are explicit, separate material connections. They do not substitute
    /// an ideal nozzle for a resolved jet, turn acoustic volume flow into mass,
    /// or infer mean-pressure prestress. Independent upstream tracts, computed
    /// exterior radiation and independently sourced branch gases retain the
    /// heat-only route's refusals; their open mean-flow feedback is not supplied.
    pub fn with_finite_gas_transport(
        mut self,
        gas: ClosedGasParcel,
        mut program: FiniteGasTransportProgram,
        gate: Arc<CancelGate>,
    ) -> Result<Self, GestureCompileError> {
        if self.finite_transport.is_some() {
            return Err(invalid("finite gas transport must bind once before playback"));
        }
        // Reuse the complete existing source admission without installing a
        // second evolving gas owner. Its cold probe publishes no physical step.
        self = self.with_finite_gas_heat(gas.clone(), program.source.clone(), gate.clone())?;
        self.finite_gas = None;
        let (capacity, work_capacity) = self.remaining_control_budgets()?;
        let events = program.controls.len();
        let initial = program.nozzles.len().checked_add(program.reservoirs.len())
            .ok_or_else(|| invalid("finite gas control count overflow"))?;
        let assignments = initial.checked_add(events)
            .and_then(|count| count.checked_add(program.source.controls.len()))
            .ok_or_else(|| invalid("finite gas control count overflow"))?;
        let levels = usize::BITS - events.saturating_sub(1).leading_zeros();
        let work = program.source.controls.len() as u128 + 2 * initial as u128
            + events as u128 * (program.nozzles.len() as u128 + 4 * (u128::from(levels) + 1));
        if program.reservoirs.len() >= MAX_SNAPSHOT_PARCELS
            || program.reservoirs.contains_key(&program.source_id)
            || program.nozzles.len() > MAX_SNAPSHOT_TRANSFERS
            || assignments > capacity || work > u128::from(work_capacity)
        {
            return Err(invalid("finite gas reservoirs, paths and controls exceed shared storage/work or reuse the source identity"));
        }
        let mut parcels = BTreeMap::from([(program.source_id, gas)]);
        let mut inputs = BTreeMap::from([(program.source_id, FiniteGasApertureTransportInput {
            thermal: FiniteGasDgTransportInput {
                gas: program.source.caloric, supplied_heat_j: 0.0,
                acoustic: Some(program.source.coupling), ports: Vec::new(),
            },
            source: None,
        })]);
        for (&id, reservoir) in &program.reservoirs {
            if gate.is_requested() {
                return Err(GestureCompileError::Render(RenderError::Voice(AcousticRealizeError::Cancelled)));
            }
            reservoir.caloric.validate().map_err(|error|
                GestureCompileError::Render(receiver_error(error.to_string())))?;
            let temperature = reservoir.gas.state().state().temperature;
            if temperature < reservoir.caloric.minimum_temperature_k
                || temperature > reservoir.caloric.maximum_temperature_k
            {
                return Err(invalid("finite reservoir caloric bracket must contain its original temperature"));
            }
            parcels.insert(id, reservoir.gas.clone());
            inputs.insert(id, FiniteGasApertureTransportInput {
                thermal: FiniteGasDgTransportInput {
                    gas: reservoir.caloric, supplied_heat_j: 0.0, acoustic: None, ports: Vec::new(),
                },
                source: None,
            });
        }
        program.nozzles.sort_unstable_by_key(|link| link.id);
        if program.nozzles.windows(2).any(|pair| pair[0].id == pair[1].id) {
            return Err(invalid("finite gas nozzles require unique original identities"));
        }
        for link in &program.nozzles {
            if !finite(link.nozzle) || link.parcels[0] == link.parcels[1]
                || link.parcels.iter().any(|id| !parcels.contains_key(id))
            {
                return Err(invalid("finite gas nozzles require admitted geometry and two distinct existing inventories"));
            }
            link.nozzle_solve.validate().map_err(|error|
                GestureCompileError::Render(receiver_error(error.to_string())))?;
        }
        program.controls.sort_unstable_by_key(|event| (event.sample, event.id));
        if program.controls.windows(2).any(|pair|
            (pair[0].sample, pair[0].id) == (pair[1].sample, pair[1].id))
        {
            return Err(invalid("finite gas controls may assign each (sample, nozzle) only once"));
        }
        for event in &program.controls {
            if gate.is_requested() {
                return Err(GestureCompileError::Render(RenderError::Voice(AcousticRealizeError::Cancelled)));
            }
            if event.sample >= self.config.samples || !finite(event.nozzle)
                || program.nozzles.binary_search_by_key(&event.id, |link| link.id).is_err()
            {
                return Err(invalid("finite gas controls require an existing nozzle, finite geometry and an in-window sample"));
            }
        }
        let mut playback = FiniteGasTransportPlayback {
            held: program.nozzles.clone(), candidate: program.nozzles.clone(),
            program, parcels, inputs, gate, next: 0, candidate_next: 0,
            heat_next: 0, held_power_w: 0.0, last: None,
            assignments, work: work as u64,
        };
        let CoupledAperture::Dg(system) = &self.system else {
            unreachable!("original finite heat admission requires DG");
        };
        let (_, pressure) = super::pressure_candidate(&self.controls, self.next,
            self.held_pressure_pa, 0);
        let force = self.forces.as_ref().map_or(0.0, |forces| forces.candidate(0));
        // Exact first-sample preflight through the native owners, on private
        // snapshots. Construction never consumes inventory or advances a clock.
        let _ = playback.preview(system, TubeDrive {
            upstream_pressure_pa: pressure, body_flow_m3_s: 0.0,
        }, force).map_err(GestureCompileError::Render)?;
        self.finite_transport = Some(playback);
        Ok(self)
    }

    /// Actual accepted source and reservoir inventories, keyed by their IDs.
    #[must_use]
    pub fn finite_gas_parcels(&self) -> Option<&BTreeMap<u64, ClosedGasParcel>> {
        self.finite_transport.as_ref().map(|state| &state.parcels)
    }

    /// Admitted physical connections, heat program and original finite sources.
    #[must_use]
    pub fn finite_gas_transport_program(&self) -> Option<&FiniteGasTransportProgram> {
        self.finite_transport.as_ref().map(|state| &state.program)
    }

    /// Cancellation gate retained by the native transport runtime.
    #[must_use]
    pub fn finite_gas_transport_cancel_gate(&self) -> Option<&Arc<CancelGate>> {
        self.finite_transport.as_ref().map(|state| &state.gate)
    }

    /// Held physical nozzles after the last successfully observed sample.
    #[must_use]
    pub fn held_gas_nozzles(&self) -> &[SnapshotNozzleLink] {
        self.finite_transport.as_ref().map_or(&[], |state| &state.held)
    }

    /// Nozzle assignments published with successful source samples.
    #[must_use]
    pub fn applied_nozzle_controls(&self) -> &[GasNozzleControl] {
        self.finite_transport.as_ref().map_or(&[], |state| &state.program.controls[..state.next])
    }

    /// Heat assignments published by this conservative transport program.
    #[must_use]
    pub fn applied_gas_transport_heat_controls(&self) -> &[GasHeatControl] {
        self.finite_transport.as_ref().map_or(&[], |state|
            &state.program.source.controls[..state.heat_next])
    }

    /// Unmodified native transport, acoustic and finite reservoir accounts.
    #[must_use]
    pub fn last_finite_gas_transport_step(&self) -> Option<&FiniteGasApertureTransportStep> {
        self.finite_transport.as_ref().and_then(|state| state.last.as_ref())
    }
}
