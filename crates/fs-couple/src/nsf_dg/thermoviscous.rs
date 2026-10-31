//! Resolved nonlinear flow driving a complete thermoviscous acoustic assembly.
//!
//! The NSF zone and FE region are distinct, non-overlapping physical volumes.
//! The original characteristic solve sees the field's exact held-input port
//! response, including every attached finite DG passage. The final NSF trial,
//! complete FE state, all duct trials and an optional finite material source
//! publish together. There is one resolved source; no jet card, independent
//! force waveform or added integrator enters. A finite source retains its
//! original calorics and actual NSF mass/total-energy exchange; prescribed
//! heat and volume motion use that same material owner's midpoint work law.
//!
//! This remains a zero-reference-flow acoustic projection. The linear field
//! does not receive mean mass, species or advected enthalpy. Original NSF mass
//! and energy exchange and signed pressure/flow covariance stay explicit. Zone
//! location, transverse projection, spatial meshes and clocks require their
//! own refinement; this adapter does not qualify a jet or instrument.

use super::{
    NsfDgDrive, NsfDgError, NsfDgLimits, NsfDgReservoirFrame, NsfFiniteReservoirFrame,
    NsfFiniteReservoirInventory, NsfFiniteReservoirLimits, coupled_exchange,
    finite_reservoir::{self, FiniteReservoir}, poll, same_clock,
    reservoir::multiple::ReservoirRadiations, solve_characteristic_interval, upstream,
};
use crate::cavity_modal::interior::thermoviscous::dg::load::{
    ThermoviscousDgLoad, ThermoviscousDgLoadTrial,
};
use fs_duct::dg::DgPortResponse;
use fs_exec::CancelGate;
use fs_lbm::lattice3::Lattice3;
use fs_lbm::nsf3::{
    EnergyLedger, NsfBoundaryDrive, NsfBoundaryExchange, NsfCharacteristicPatch,
    NsfCharacteristicTrial, NsfGrid3, NsfOpenState, NsfReservoirExchange,
    NsfReservoirFlow,
};
use fs_material::gas::GasState;
use fs_material::thermochemical_gas::energy::parcel::ClosedGasParcel;

/// Original pressure/PCM observation of the complete resolved assembly.
pub mod performance;
mod radiation;

/// Accepted interface and original fluid accounts. Complete accepted FE and
/// duct states remain available through `load().system().last_frame()` without
/// cloning every field merely to return a step receipt.
#[derive(Clone, Debug)]
pub struct NsfThermoviscousDgFrame {
    pub step: u64,
    pub evaluations: usize,
    pub incoming_characteristic_m_s: f64,
    pub drive: NsfDgDrive,
    /// Original complete low/high face accounts, including patch remainders.
    pub faces: [NsfBoundaryExchange; 2],
    pub coupled_patch: Option<NsfBoundaryExchange>,
    pub source_patch: Option<NsfBoundaryExchange>,
    /// Conservative exchange through exactly the selected FE connection.
    /// This is already included in the complete face account above.
    pub field_exchange: NsfBoundaryExchange,
    pub reservoir_mass_out_kg: f64,
    pub reservoir_flows: Vec<NsfReservoirFlow>,
    pub reservoir_exchanges: Vec<NsfReservoirExchange>,
    /// Original shared exterior's reciprocal exchange at each attached mouth,
    /// in attachment order. Its work is not another FE/DG storage debit.
    pub reservoir_radiations: Vec<NsfDgReservoirFrame>,
    /// Original material-owned storage and conservative exchange, when bound.
    pub finite_reservoir: Option<NsfFiniteReservoirFrame>,
    pub flow: EnergyLedger,
    pub flow_mass_change_kg: f64,
    pub flow_energy_change_j: f64,
    /// Mean outward NSF flow, held inward on the selected FE port [m³/s].
    pub inlet_flow_m3_s: f64,
    pub inlet_pressure_pa: f64,
    pub inlet_work_j: f64,
    /// Complete acoustic FE plus every attached original finite duct.
    pub acoustic_stored_energy_j: f64,
    pub acoustic_balance_residual_j: f64,
    pub acoustic_balance_tolerance_j: f64,
    pub pressure_residual_pa: f64,
    /// Resolved integral(p' Q) minus mean(p') integral(Q); never heat.
    pub projection_covariance_j: f64,
    /// Resolved acoustic work minus actual complete-field inlet work.
    pub acoustic_work_residual_j: f64,
    pub port_work_residual_j: f64,
    /// Streamed mass minus reference density times continuum swept volume.
    pub mass_projection_difference_kg: f64,
}

/// Original nonlinear fluid and complete linear acoustic histories on one
/// coupled clock. Construction consumes fresh owners, so it cannot silently
/// discard attachments from an existing `NsfDg` assembly.
pub struct NsfThermoviscousDg<'a, 'c, L: Lattice3<27>> {
    zone: NsfGrid3<L>,
    load: ThermoviscousDgLoad<'a, 'c>,
    limits: NsfDgLimits,
    high: bool,
    reference: GasState,
    area_m2: f64,
    steps: u64,
    patch_background: [Option<NsfBoundaryDrive>; 2],
    finite_reservoir: Option<FiniteReservoir>,
    reservoir_radiations: Option<ReservoirRadiations>,
}

impl<'a, 'c, L: Lattice3<27>> NsfThermoviscousDg<'a, 'c, L> {
    /// Bind the complete selected characteristic face to the external FE port.
    pub fn new(zone: NsfGrid3<L>, load: ThermoviscousDgLoad<'a, 'c>, high: bool,
        limits: NsfDgLimits, gate: &CancelGate) -> Result<Self, NsfDgError>
    {
        Self::bind(zone, load, high, limits, gate, [false; 2])
    }

    /// Bind only the authored downstream aperture, retaining its remainder.
    pub fn new_on_patch(zone: NsfGrid3<L>, load: ThermoviscousDgLoad<'a, 'c>, high: bool,
        limits: NsfDgLimits, gate: &CancelGate) -> Result<Self, NsfDgError>
    {
        let mut selected = [false; 2];
        selected[usize::from(high)] = true;
        Self::bind(zone, load, high, limits, gate, selected)
    }

    /// Select a source aperture while the acoustic field receives a whole face.
    pub fn new_with_source_patch(zone: NsfGrid3<L>, load: ThermoviscousDgLoad<'a, 'c>, high: bool,
        limits: NsfDgLimits, gate: &CancelGate) -> Result<Self, NsfDgError>
    {
        let mut selected = [false; 2];
        selected[usize::from(!high)] = true;
        Self::bind(zone, load, high, limits, gate, selected)
    }

    /// Select both original apertures; both surrounding laws remain unchanged.
    pub fn new_on_patches(zone: NsfGrid3<L>, load: ThermoviscousDgLoad<'a, 'c>, high: bool,
        limits: NsfDgLimits, gate: &CancelGate) -> Result<Self, NsfDgError>
    {
        Self::bind(zone, load, high, limits, gate, [true; 2])
    }

    fn bind(zone: NsfGrid3<L>, load: ThermoviscousDgLoad<'a, 'c>, high: bool,
        limits: NsfDgLimits, gate: &CancelGate, selected: [bool; 2]) -> Result<Self, NsfDgError>
    {
        poll(gate)?;
        limits.validate_for_cells(zone.cell_count())?;
        let open = zone.characteristic_axis().ok_or(NsfDgError::Input(
            "resolved field coupling requires the original characteristic axis"))?;
        let reference = zone.characteristic_reference_gas(high).map_err(NsfDgError::Flow)?;
        let drives = zone.boundary_drives().ok_or(NsfDgError::Input(
            "resolved field coupling requires original boundary drives"))?;
        let mut patch_background = [None; 2];
        for side in 0..2 {
            if selected[side] != zone.characteristic_patch_on(side != 0).is_some() {
                return Err(NsfDgError::Input(
                    "select every source and field aperture explicitly; no whole-face role may hide a patch"));
            }
            if selected[side] { patch_background[side] = Some(drives[side]); }
        }
        let traces = if selected[usize::from(high)] { zone.open_patch_traces_on(high) }
            else { zone.open_boundary_traces(high) }
            .ok_or(NsfDgError::Input("resolved field coupling needs original interface quadrature"))?;
        let area = traces.first().ok_or(NsfDgError::Input("resolved field interface is empty"))?
            .area_m2 * traces.len() as f64;
        if zone.steps() != 0 || load.system().accepted_steps() != 0 || load.mouth_count() != 0
            || open.faces.iter().any(|face| face.reference_velocity_m_s != [0.0; 3])
            || load.system().source().gas() != &reference
            || !area.is_finite() || area <= 0.0
            || (area - load.area_m2()).abs() > 32.0 * f64::EPSILON * area.max(load.area_m2())
            || !same_clock(zone.law().dt_s * limits.substeps as f64, load.duration_s())
        {
            return Err(NsfDgError::Input(
                "NSF/complete field requires cold original owners, equal gas/area/clock, zero reference flow and explicitly attached mouths"));
        }
        poll(gate)?;
        Ok(Self { zone, load, limits, high, reference, area_m2: area, steps: 0,
            patch_background, finite_reservoir: None, reservoir_radiations: None })
    }

    /// Attach the original finite material owner at one quiescent transverse
    /// opening. It exclusively supplies its accepted pressure and temperature.
    /// The same named-mixture or NASA9 source reconstructs storage from actual
    /// NSF total-energy and mass exchange; acoustic work is not debited again.
    ///
    /// Optional heat and volume motion use the original simultaneous midpoint
    /// law after the flow exchange. This held-source partition is first order
    /// and requires interval refinement. No chamber momentum is stored.
    pub fn with_finite_reservoir(
        mut self, aperture: usize, volume_m3: f64,
        limits: NsfFiniteReservoirLimits, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        if self.steps != 0 || self.finite_reservoir.is_some()
            || self.reservoir_radiation_at(aperture).is_some()
        {
            return Err(NsfDgError::Input(
                "finite reservoir needs one cold complete-field source aperture distinct from every passive mouth",
            ));
        }
        self.finite_reservoir = Some(FiniteReservoir::bind(
            &self.zone, aperture, volume_m3, limits, gate,
        )?);
        Ok(self)
    }

    /// Accepted NASA9 parcel; named mixtures use finite_reservoir_inventory.
    pub fn finite_reservoir(&self) -> Option<&ClosedGasParcel> {
        self.finite_reservoir.as_ref().and_then(FiniteReservoir::thermochemical)
    }

    /// Accepted source-owned storage for either supported material family.
    pub fn finite_reservoir_inventory(&self) -> Option<NsfFiniteReservoirInventory> {
        self.finite_reservoir.as_ref().map(FiniteReservoir::inventory)
    }

    pub fn finite_reservoir_aperture(&self) -> Option<usize> {
        self.finite_reservoir.as_ref().map(|owner| owner.aperture)
    }

    pub fn finite_reservoir_limits(&self) -> Option<&NsfFiniteReservoirLimits> {
        self.finite_reservoir.as_ref().map(FiniteReservoir::limits)
    }

    pub const fn zone(&self) -> &NsfGrid3<L> { &self.zone }
    pub const fn load(&self) -> &ThermoviscousDgLoad<'a, 'c> { &self.load }
    pub const fn limits(&self) -> &NsfDgLimits { &self.limits }
    pub const fn accepted_steps(&self) -> u64 { self.steps }
    pub const fn reference_gas(&self) -> &GasState { &self.reference }
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    pub const fn coupled_face_is_high(&self) -> bool { self.high }
    pub fn coupled_patch(&self) -> Option<NsfCharacteristicPatch> {
        self.zone.characteristic_patch_on(self.high)
    }
    pub fn source_patch(&self) -> Option<NsfCharacteristicPatch> {
        self.zone.characteristic_patch_on(!self.high)
    }
    pub const fn patch_background_drive_on(&self, high: bool) -> Option<NsfBoundaryDrive> {
        self.patch_background[high as usize]
    }

    /// Admit a source at its physical location without changing any owner.
    /// An attached finite chamber exclusively owns its source controls and
    /// retains its accepted material state. No pressure waveform replaces it.
    pub fn validate_drive(&self, drive: NsfDgDrive) -> Result<(), NsfDgError> {
        self.source_conditions(drive).map(|_| ())
    }

    fn source_conditions(&self, drive: NsfDgDrive)
        -> Result<(NsfBoundaryDrive, Option<(usize, NsfOpenState)>), NsfDgError>
    {
        if !drive.is_finite() {
            return Err(NsfDgError::Input("resolved field source control must be finite"));
        }
        FiniteReservoir::admit_drive(self.finite_reservoir.as_ref(), drive)?;
        let zero = NsfBoundaryDrive::IncomingCharacteristic(0.0);
        let conditions = match drive {
            NsfDgDrive::Boundary(value) => Ok((value, None)),
            NsfDgDrive::Reservoir { pressure_pa, temperature_k } => Ok((zero, Some((0,
                self.zone.reservoir_aperture_state_at(pressure_pa, temperature_k)
                    .map_err(NsfDgError::Flow)?)))),
            NsfDgDrive::ReservoirAperture { aperture, pressure_pa, temperature_k } => Ok((zero, Some((aperture,
                self.zone.reservoir_aperture_state_for(aperture, pressure_pa, temperature_k)
                    .map_err(NsfDgError::Flow)?)))),
            NsfDgDrive::FiniteReservoir { aperture }
            | NsfDgDrive::FiniteReservoirActuation { aperture, .. } => Ok((zero, Some((aperture,
                self.finite_reservoir.as_ref().expect("admitted finite reservoir")
                    .state(&self.zone)?)))),
            NsfDgDrive::UpstreamPressureRise(_) => Err(NsfDgError::Input(
                "complete-field source has no second NSF-side acoustic load for a series pressure rise")),
        }?;
        if conditions.1.is_some_and(|(aperture, _)| self.reservoir_radiation_at(aperture).is_some()) {
            return Err(NsfDgError::Input(
                "prescribed or finite breath and every passive mouth must have distinct original owners"));
        }
        Ok(conditions)
    }

    pub fn step(&mut self, incoming_m_s: f64, gate: &CancelGate)
        -> Result<NsfThermoviscousDgFrame, NsfDgError>
    {
        self.step_with_boundary(NsfBoundaryDrive::IncomingCharacteristic(incoming_m_s), gate)
    }

    pub fn step_with_boundary(&mut self, boundary: NsfBoundaryDrive, gate: &CancelGate)
        -> Result<NsfThermoviscousDgFrame, NsfDgError>
    {
        self.step_with_drive(NsfDgDrive::Boundary(boundary), gate)
    }

    pub fn step_with_reservoir_aperture(&mut self, aperture: usize, pressure_pa: f64,
        temperature_k: f64, gate: &CancelGate) -> Result<NsfThermoviscousDgFrame, NsfDgError>
    {
        self.step_with_drive(NsfDgDrive::ReservoirAperture { aperture, pressure_pa, temperature_k }, gate)
    }

    /// Solve on the original NSF trial trajectory, then stage the complete
    /// field, all finite passages and finite material storage before the shared
    /// work/cancellation gate. Any refusal retains all accepted histories.
    pub fn step_with_drive(&mut self, drive: NsfDgDrive, gate: &CancelGate)
        -> Result<NsfThermoviscousDgFrame, NsfDgError>
    {
        poll(gate)?;
        if self.steps >= self.limits.maximum_steps {
            return Err(NsfDgError::Input("resolved field coupling exhausted its physical horizon"));
        }
        let (external, reservoir) = self.source_conditions(drive)?;
        if self.reservoir_radiations.is_some() {
            return self.step_with_radiations(drive, external, reservoir, gate);
        }
        let reservoirs: Vec<_> = reservoir.into_iter().collect();
        let relation = self.load.port_response(gate).map_err(NsfDgError::Thermoviscous)?;
        let limits = self.limits;
        let side = usize::from(self.high);
        let reference = self.reference;
        let step = self.steps + 1;
        let load = &mut self.load;
        let finite = &mut self.finite_reservoir;
        let frame = solve_characteristic_interval(&mut self.zone, external, side,
            self.patch_background, &reservoirs, relation, limits, gate,
            |fluid, incoming, evaluations| accept(fluid, load, finite, side, relation, reference,
                limits, incoming, drive, evaluations, step, gate),
        )?;
        self.steps += 1;
        Ok(frame)
    }
}

#[allow(clippy::too_many_arguments)]
fn accept<L: Lattice3<27>>(
    fluid: NsfCharacteristicTrial<'_, L>, load: &mut ThermoviscousDgLoad<'_, '_>,
    finite: &mut Option<FiniteReservoir>,
    side: usize, relation: DgPortResponse, reference: GasState, limits: NsfDgLimits,
    incoming: f64, drive: NsfDgDrive, evaluations: usize, step: u64, gate: &CancelGate,
) -> Result<NsfThermoviscousDgFrame, NsfDgError> {
    let face = coupled_exchange(&fluid, side)?;
    let q = face.normal_volume_flow_m3_s;
    let dt = load.duration_s();
    let acoustic = load.preview_step(q, gate).map_err(NsfDgError::Thermoviscous)?;
    accept_trial(fluid, acoustic, finite, side, relation, reference, limits, dt,
        incoming, drive, evaluations, step, Vec::new(), gate)
}

/// One acceptance owner for scalar and reciprocal-mouth source equations.
/// All caller-specific mouth audits precede this final shared publication.
#[allow(clippy::too_many_arguments)]
fn accept_trial<L: Lattice3<27>>(
    fluid: NsfCharacteristicTrial<'_, L>, acoustic: ThermoviscousDgLoadTrial<'_>,
    finite: &mut Option<FiniteReservoir>,
    side: usize, relation: DgPortResponse, reference: GasState, limits: NsfDgLimits,
    dt: f64, incoming: f64, drive: NsfDgDrive, evaluations: usize, step: u64,
    reservoir_radiations: Vec<NsfDgReservoirFrame>, gate: &CancelGate,
) -> Result<NsfThermoviscousDgFrame, NsfDgError> {
    let face = coupled_exchange(&fluid, side)?;
    let q = face.normal_volume_flow_m3_s;
    let pressure = acoustic.mean_pressure_pa();
    let work = acoustic.inlet_work_j();
    let audit = upstream::audit(face, pressure, work, 0.0, relation, reference, dt, limits)?;
    let mass_change = fluid.grid().total_mass_kg() - fluid.previous_mass_kg();
    let energy_change = fluid.grid().total_energy_j() - fluid.previous_energy_j();
    let reservoir_mass = fluid.reservoir_mass_out_kg();
    if ![mass_change, energy_change, reservoir_mass].iter().all(|x| x.is_finite()) {
        return Err(NsfDgError::Input("resolved field source mass/energy account overflowed"));
    }
    let chamber = finite_reservoir::prepare(finite.as_ref(), &fluid, drive, dt, gate)?;
    let mut frame = NsfThermoviscousDgFrame {
        step, evaluations, incoming_characteristic_m_s: incoming, drive,
        faces: *fluid.boundary_exchange(),
        coupled_patch: fluid.patch_exchange_on(side != 0).copied(),
        source_patch: fluid.patch_exchange_on(side == 0).copied(), field_exchange: face,
        reservoir_mass_out_kg: reservoir_mass, reservoir_flows: fluid.reservoir_flows().to_vec(),
        reservoir_exchanges: fluid.reservoir_exchanges().to_vec(), flow: *fluid.energy_ledger(),
        finite_reservoir: None, reservoir_radiations,
        flow_mass_change_kg: mass_change, flow_energy_change_j: energy_change,
        inlet_flow_m3_s: q, inlet_pressure_pa: pressure, inlet_work_j: work,
        acoustic_stored_energy_j: acoustic.frame().stored_energy_j,
        acoustic_balance_residual_j: acoustic.frame().balance_residual_j,
        acoustic_balance_tolerance_j: acoustic.frame().balance_tolerance_j,
        pressure_residual_pa: audit.pressure, projection_covariance_j: audit.covariance,
        acoustic_work_residual_j: audit.work, port_work_residual_j: audit.port,
        mass_projection_difference_kg: audit.mass,
    };
    poll(gate)?;
    // Original publications cannot refuse. No cancellation poll separates the
    // fluid, complete field, finite ducts and retained source inventory.
    fluid.accept();
    acoustic.accept();
    frame.finite_reservoir = finite_reservoir::publish(finite, chamber);
    Ok(frame)
}
