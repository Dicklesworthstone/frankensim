//! Reversible intervals of the original characteristic NSF owner.
//! No source function is cloned or reevaluated on another clock: characteristic
//! topology already excludes manufactured volume sources. Every physical state,
//! boundary force chart, source material and cumulative ledger remains retained.

use super::{EnergyLedger, NsfBoundaryDrive, NsfBoundaryExchange, NsfError, NsfGrid3,
    NsfOpenState, NsfReservoirExchange, NsfReservoirFlow};
use crate::lattice3::Lattice3;

/// A complete unpublished interval. Dropping it changes no accepted field,
/// characteristic, material cache, boundary history or clock.
pub struct NsfCharacteristicTrial<'a, L: Lattice3<27>> {
    owner: &'a mut NsfGrid3<L>,
    candidate: NsfGrid3<L>,
    exchange: [NsfBoundaryExchange; 2],
    patch_exchange: [Option<NsfBoundaryExchange>; 2],
    ledger: EnergyLedger,
    mass_before_kg: f64,
    energy_before_j: f64,
    reservoir_mass_out_kg: f64,
    reservoir_flows: Vec<NsfReservoirFlow>,
    reservoir_exchanges: Vec<NsfReservoirExchange>,
}

impl<L: Lattice3<27>> NsfCharacteristicTrial<'_, L> {
    /// Original candidate grid, including source thermochemistry and face traces.
    pub const fn grid(&self) -> &NsfGrid3<L> { &self.candidate }
    /// Summed mass, momentum, energy and work; pressure/flow are interval means.
    pub const fn boundary_exchange(&self) -> &[NsfBoundaryExchange; 2] { &self.exchange }
    /// Original exchange when exactly one aperture is attached; `None` for zero
    /// or two. Each aperture is already included in its complete face account.
    pub const fn patch_exchange(&self) -> Option<&NsfBoundaryExchange> {
        match &self.patch_exchange {
            [Some(exchange), None] | [None, Some(exchange)] => Some(exchange),
            _ => None,
        }
    }
    /// Original interval exchange for the selected face's aperture.
    pub const fn patch_exchange_on(&self, high: bool) -> Option<&NsfBoundaryExchange> {
        self.patch_exchange[high as usize].as_ref()
    }
    /// Original step accounts summed over this interval, with final storage and
    /// the interval total-energy residual. Viscous heat keeps its budget-remainder
    /// qualification; no independently resolved heat account is invented here.
    pub const fn energy_ledger(&self) -> &EnergyLedger { &self.ledger }
    /// Accepted mass before any incoming-data or step trial [kg].
    pub const fn previous_mass_kg(&self) -> f64 { self.mass_before_kg }
    /// Accepted sensible-chart total energy before the trial [J].
    pub const fn previous_energy_j(&self) -> f64 { self.energy_before_j }
    /// Original streamed mass leaving the localized reservoir aperture over
    /// this interval [kg], independently counted from characteristic exchange.
    /// Zero when no aperture is attached.
    pub const fn reservoir_mass_out_kg(&self) -> f64 { self.reservoir_mass_out_kg }
    /// Original per-tile pressure and continuum flow averaged over every
    /// original step, with pressure work summed over the complete interval.
    /// Tile identity and signs are retained, independently of mass/energy flux.
    pub fn reservoir_flows(&self) -> &[NsfReservoirFlow] { &self.reservoir_flows }
    /// Aggregate original FV pressure work over the complete trial interval.
    /// Means use every substep; covariance is recomputed for this interval,
    /// not summed from shorter projections with different temporal means.
    pub fn reservoir_exchange(&self) -> Option<&NsfReservoirExchange> {
        if self.candidate.reservoir_aperture_count() == 1 { self.reservoir_exchange_at(0) } else { None }
    }
    /// Original FV interval account of the explicitly selected aperture.
    pub fn reservoir_exchange_at(&self, aperture: usize) -> Option<&NsfReservoirExchange> {
        self.reservoir_exchanges.get(aperture)
    }
    /// Original per-owner interval accounts in declared aperture order.
    pub fn reservoir_exchanges(&self) -> &[NsfReservoirExchange] { &self.reservoir_exchanges }
    /// Original interval-mean tile flows of one owner, preserving signed flow.
    pub fn reservoir_flows_at(&self, aperture: usize) -> &[NsfReservoirFlow] {
        let first = self.reservoir_flows.partition_point(|flow| flow.aperture < aperture);
        let end = self.reservoir_flows.partition_point(|flow| flow.aperture <= aperture);
        &self.reservoir_flows[first..end]
    }
    /// Publish the already admitted complete candidate. No fallible work follows.
    pub fn accept(self) { *self.owner = self.candidate; }
}

impl<L: Lattice3<27>> NsfGrid3<L> {
    /// Preview held incoming characteristics for a positive number of original
    /// lattice steps. The caller bounds cell-step work and can cancel before
    /// the snapshot or between steps. A failed interval changes no accepted state.
    /// Every original per-step material, Mach, thermal and boundary refusal applies.
    pub fn preview_characteristic_interval(
        &mut self,
        incoming_m_s: [f64; 2],
        steps: usize,
        maximum_cell_steps: usize,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        self.preview_boundary_interval(
            incoming_m_s.map(NsfBoundaryDrive::IncomingCharacteristic),
            steps, maximum_cell_steps, cancelled,
        )
    }

    /// Preview original steps under typed incoming-wave or pressure controls.
    /// A pressure control is imposed by the original local reconstruction at
    /// every predictor, not by another flow solve or a rebaked source waveform.
    /// Dropping or refusing the trial preserves both accepted controls and all
    /// original populations, force/energy history and material state.
    /// An attached aperture retains its accepted control throughout the trial.
    pub fn preview_boundary_interval(
        &mut self,
        drives: [NsfBoundaryDrive; 2],
        steps: usize,
        maximum_cell_steps: usize,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        self.preview_boundary_interval_inner(drives, [None; 2], &[], steps, maximum_cell_steps, cancelled)
    }

    /// Preview a selected aperture control and its complementary face controls
    /// in one original transaction. The patch shape and source gas are retained.
    /// A rejected or dropped candidate publishes no controls, fields or accounts.
    pub fn preview_boundary_interval_with_patch(
        &mut self, drives: [NsfBoundaryDrive; 2], patch_drive: NsfBoundaryDrive,
        steps: usize, maximum_cell_steps: usize, cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        let patch_drives = self.single_patch_controls(patch_drive)?;
        self.preview_boundary_interval_inner(drives, patch_drives, &[], steps, maximum_cell_steps, cancelled)
    }

    /// Preview independent aperture controls on the low and high faces, together
    /// with their complementary face drives. `None` retains that face's accepted
    /// aperture control. All original populations, work and material histories
    /// are staged once; refusal or a dropped trial publishes no participant.
    pub fn preview_boundary_interval_with_patches(
        &mut self, drives: [NsfBoundaryDrive; 2],
        patch_drives: [Option<NsfBoundaryDrive>; 2],
        steps: usize, maximum_cell_steps: usize, cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        self.preview_boundary_interval_inner(drives, patch_drives, &[],
            steps, maximum_cell_steps, cancelled)
    }

    /// Hold a new state at the original transverse reservoir while solving
    /// the two characteristic controls and optional selected patch together.
    /// The ghost and its force/energy-flux history change only on the candidate;
    /// dropping or refusing the trial also rolls back the reservoir control.
    pub fn preview_boundary_interval_with_reservoir(
        &mut self, drives: [NsfBoundaryDrive; 2], patch_drive: Option<NsfBoundaryDrive>,
        reservoir_state: NsfOpenState, steps: usize, maximum_cell_steps: usize,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        if self.reservoir_aperture_count() != 1 {
            return Err(NsfError::Input("singular reservoir trial requires exactly one aperture"));
        }
        self.preview_boundary_interval_with_reservoirs(drives, patch_drive, &[(0, reservoir_state)],
            steps, maximum_cell_steps, cancelled)
    }

    /// Stage explicit independent transverse controls with the original face
    /// drives and optional singular characteristic patch. Unknown/duplicate
    /// reservoir IDs refuse; unnamed owners retain their accepted state.
    pub fn preview_boundary_interval_with_reservoirs(
        &mut self, drives: [NsfBoundaryDrive; 2], patch_drive: Option<NsfBoundaryDrive>,
        controls: &[(usize, NsfOpenState)], steps: usize, maximum_cell_steps: usize,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        let patch_drives = match patch_drive {
            Some(drive) => self.single_patch_controls(drive)?,
            None => [None; 2],
        };
        self.preview_boundary_interval_with_patches_and_reservoirs(drives, patch_drives, controls,
            steps, maximum_cell_steps, cancelled)
    }

    /// The same two-aperture transaction with a proposed state at the original
    /// transverse reservoir. No ghost, control or work account is published
    /// until the complete original candidate is accepted.
    pub fn preview_boundary_interval_with_patches_and_reservoir(
        &mut self, drives: [NsfBoundaryDrive; 2],
        patch_drives: [Option<NsfBoundaryDrive>; 2],
        reservoir_state: NsfOpenState, steps: usize, maximum_cell_steps: usize,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        if self.reservoir_aperture_count() != 1 {
            return Err(NsfError::Input("singular reservoir trial requires exactly one aperture"));
        }
        self.preview_boundary_interval_with_patches_and_reservoirs(
            drives, patch_drives, &[(0, reservoir_state)], steps, maximum_cell_steps, cancelled)
    }

    /// Stage independently selected characteristic patches on both faces and
    /// transverse reservoir controls in the same original candidate. No field,
    /// control or work account publishes before the complete interval accepts.
    pub fn preview_boundary_interval_with_patches_and_reservoirs(
        &mut self, drives: [NsfBoundaryDrive; 2],
        patch_drives: [Option<NsfBoundaryDrive>; 2],
        controls: &[(usize, NsfOpenState)], steps: usize, maximum_cell_steps: usize,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        self.preview_boundary_interval_inner(drives, patch_drives, controls,
            steps, maximum_cell_steps, cancelled)
    }

    fn preview_boundary_interval_inner(
        &mut self, drives: [NsfBoundaryDrive; 2], patch_drives: [Option<NsfBoundaryDrive>; 2],
        controls: &[(usize, NsfOpenState)], steps: usize, maximum_cell_steps: usize,
        cancelled: impl Fn() -> bool,
    ) -> Result<NsfCharacteristicTrial<'_, L>, NsfError> {
        if cancelled() {
            return Err(NsfError::Input("characteristic interval cancelled"));
        }
        if self.open.is_none() || self.sources.is_some() || steps == 0
            || self.n.checked_mul(steps).is_none_or(|work| work > maximum_cell_steps)
            || self.steps.checked_add(steps).is_none()
        {
            return Err(NsfError::Input("characteristic interval requires bounded positive source-free cell-step work"));
        }
        let mass_before_kg = self.total_mass_kg();
        let energy_before_j = self.total_energy_j();
        let mut candidate = self.characteristic_copy();
        // Batch all transverse ghosts before the original characteristic
        // rebase, preserving the historical one-aperture arithmetic.
        candidate.set_reservoir_aperture_states(controls)?;
        candidate.set_boundary_drives_with_patches(drives, patch_drives)?;
        let mut exchange = [NsfBoundaryExchange::default(); 2];
        let mut patch_exchange = candidate.characteristic_patch
            .map(|patch| patch.map(|_| NsfBoundaryExchange::default()));
        let mut ledger = EnergyLedger::default();
        let mut reservoir_mass_out_kg = 0.0;
        let mut reservoir_exchanges = candidate.empty_reservoir_exchanges();
        let mut reservoir_flows = candidate.last_reservoir_flows.clone();
        for flow in &mut reservoir_flows {
            flow.outward_volume_flow_m3_s = 0.0;
            flow.mean_pressure_pa = 0.0;
            flow.mean_acoustic_pressure_pa = 0.0;
            flow.pressure_work_out_j = 0.0;
            flow.acoustic_pressure_work_out_j = 0.0;
        }
        for _ in 0..steps {
            if cancelled() {
                return Err(NsfError::Input("characteristic interval cancelled"));
            }
            let step = candidate.step()?;
            reservoir_mass_out_kg += candidate.last_reservoir_mass_out_kg().unwrap_or(0.0);
            for total in &mut reservoir_exchanges {
                let step = candidate.last_reservoir_exchange_at(total.aperture).ok_or(NsfError::Input(
                    "original reservoir work account disappeared during the interval",
                ))?;
                if step.aperture != total.aperture || step.area_m2 != total.area_m2
                    || step.reference_pressure_pa != total.reference_pressure_pa
                {
                    return Err(NsfError::Input("reservoir work changed its original area or reference"));
                }
                total.mean_pressure_pa += step.mean_pressure_pa;
                total.mean_acoustic_pressure_pa += step.mean_acoustic_pressure_pa;
                total.normal_volume_flow_m3_s += step.normal_volume_flow_m3_s;
                total.pressure_work_out_j += step.pressure_work_out_j;
                total.acoustic_pressure_work_out_j += step.acoustic_pressure_work_out_j;
                total.mass_out_kg += step.mass_out_kg;
                total.total_energy_out_j += step.total_energy_out_j;
            }
            for (total, flow) in reservoir_flows.iter_mut().zip(candidate.last_reservoir_flows()) {
                total.outward_volume_flow_m3_s += flow.outward_volume_flow_m3_s;
                total.mean_pressure_pa += flow.mean_pressure_pa;
                total.mean_acoustic_pressure_pa += flow.mean_acoustic_pressure_pa;
                total.pressure_work_out_j += flow.pressure_work_out_j;
                total.acoustic_pressure_work_out_j += flow.acoustic_pressure_work_out_j;
            }
            ledger.internal_j = step.internal_j;
            ledger.kinetic_j = step.kinetic_j;
            ledger.viscous_heating_j += step.viscous_heating_j;
            ledger.unresolved_internal_exchange_j += step.unresolved_internal_exchange_j;
            ledger.pressure_work_j += step.pressure_work_j;
            ledger.wall_heat_j += step.wall_heat_j;
            ledger.open_boundary_energy_j += step.open_boundary_energy_j;
            ledger.open_energy_j += step.open_energy_j;
            ledger.source_energy_j += step.source_energy_j;
            for (total, step) in exchange.iter_mut().zip(
                candidate.last_boundary_exchange().expect("admitted characteristic step"))
            {
                total.mass_out_kg += step.mass_out_kg;
                for a in 0..3 { total.momentum_out_kg_m_s[a] += step.momentum_out_kg_m_s[a]; }
                total.total_energy_out_j += step.total_energy_out_j;
                total.mean_pressure_pa += step.mean_pressure_pa;
                total.normal_volume_flow_m3_s += step.normal_volume_flow_m3_s;
                total.pressure_work_out_j += step.pressure_work_out_j;
                total.mean_acoustic_pressure_pa += step.mean_acoustic_pressure_pa;
                total.acoustic_pressure_work_out_j += step.acoustic_pressure_work_out_j;
            }
            for (side, total) in patch_exchange.iter_mut().enumerate() {
                let Some(total) = total else { continue; };
                let step = candidate.last_patch_exchange_on(side == 1)
                    .expect("admitted characteristic patch");
                total.mass_out_kg += step.mass_out_kg;
                for a in 0..3 { total.momentum_out_kg_m_s[a] += step.momentum_out_kg_m_s[a]; }
                total.total_energy_out_j += step.total_energy_out_j;
                total.mean_pressure_pa += step.mean_pressure_pa;
                total.normal_volume_flow_m3_s += step.normal_volume_flow_m3_s;
                total.pressure_work_out_j += step.pressure_work_out_j;
                total.mean_acoustic_pressure_pa += step.mean_acoustic_pressure_pa;
                total.acoustic_pressure_work_out_j += step.acoustic_pressure_work_out_j;
            }
        }
        for face in &mut exchange {
            face.mean_pressure_pa /= steps as f64;
            face.mean_acoustic_pressure_pa /= steps as f64;
            face.normal_volume_flow_m3_s /= steps as f64;
        }
        for flow in &mut reservoir_flows {
            flow.outward_volume_flow_m3_s /= steps as f64;
            flow.mean_pressure_pa /= steps as f64;
            flow.mean_acoustic_pressure_pa /= steps as f64;
        }
        for total in &mut reservoir_exchanges {
            total.mean_pressure_pa /= steps as f64;
            total.mean_acoustic_pressure_pa /= steps as f64;
            total.normal_volume_flow_m3_s /= steps as f64;
            total.projection_covariance_j = total.acoustic_pressure_work_out_j
                - total.mean_acoustic_pressure_pa * total.normal_volume_flow_m3_s
                    * (steps as f64 * candidate.law.dt_s);
        }
        for face in patch_exchange.iter_mut().flatten() {
            face.mean_pressure_pa /= steps as f64;
            face.mean_acoustic_pressure_pa /= steps as f64;
            face.normal_volume_flow_m3_s /= steps as f64;
        }
        let energy = candidate.total_energy_j();
        ledger.residual_rel = (energy - energy_before_j - ledger.wall_heat_j
            - ledger.open_boundary_energy_j - ledger.open_energy_j - ledger.source_energy_j).abs()
            / energy_before_j.abs().max(energy.abs()).max(1e-300);
        if exchange.iter().chain(patch_exchange.iter().flatten()).any(|face| ![
            face.mass_out_kg, face.total_energy_out_j, face.mean_pressure_pa,
            face.normal_volume_flow_m3_s, face.pressure_work_out_j,
            face.mean_acoustic_pressure_pa, face.acoustic_pressure_work_out_j,
        ].into_iter().chain(face.momentum_out_kg_m_s).all(f64::is_finite))
            || reservoir_flows.iter().any(|flow| !flow.is_finite())
            || reservoir_exchanges.iter().any(|exchange| !exchange.is_finite())
            || ![mass_before_kg, energy_before_j, energy, candidate.total_mass_kg(),
                ledger.internal_j, ledger.kinetic_j, ledger.viscous_heating_j,
                ledger.unresolved_internal_exchange_j, ledger.pressure_work_j,
                ledger.wall_heat_j, ledger.open_boundary_energy_j,
                ledger.open_energy_j, ledger.source_energy_j, ledger.residual_rel,
                reservoir_mass_out_kg]
                .iter().all(|x| x.is_finite())
        {
            return Err(NsfError::Input("characteristic interval exchange overflowed"));
        }
        if cancelled() {
            return Err(NsfError::Input("characteristic interval cancelled"));
        }
        Ok(NsfCharacteristicTrial {
            owner: self, candidate, exchange, patch_exchange, ledger, mass_before_kg, energy_before_j,
            reservoir_mass_out_kg, reservoir_flows, reservoir_exchanges,
        })
    }

    // Characteristic topology prohibits SourceFn. Copy every remaining owner
    // field rather than reconstructing populations from macroscopic values.
    fn characteristic_copy(&self) -> Self {
        Self {
            spec: self.spec.clone(), law: self.law.clone(), material: self.material.clone(),
            sgs: self.sgs, sgs_transport: self.sgs_transport.clone(),
            n: self.n, solid: self.solid.clone(), f: self.f.clone(), e: self.e.clone(),
            theta: self.theta.clone(), force: self.force.clone(), div: self.div.clone(),
            wall_out: self.wall_out, open: self.open,
            characteristic_pressure_offset_pa: self.characteristic_pressure_offset_pa,
            characteristic_inflow_temperature_k: self.characteristic_inflow_temperature_k,
            characteristic_patch: self.characteristic_patch,
            open_faces: self.open_faces.clone(),
            reservoir_faces: self.reservoir_faces, open_out: self.open_out,
            reservoir_apertures: self.reservoir_apertures.clone(),
            last_reservoir_mass_out_kg: self.last_reservoir_mass_out_kg,
            last_reservoir_flows: self.last_reservoir_flows.clone(),
            last_reservoir_exchanges: self.last_reservoir_exchanges.clone(),
            last_boundary_exchange: self.last_boundary_exchange, steps: self.steps,
            last_patch_exchange: self.last_patch_exchange,
            cumulative: self.cumulative, edge_velocity: self.edge_velocity,
            moving_walls: self.moving_walls, sources: None, src_now: self.src_now.clone(),
            marker: core::marker::PhantomData,
        }
    }
}
