//! Original resolved NSF flow coupled to the original causal DG bore.
//!
//! A complete planar characteristic face or one explicit aperture patch is projected onto the uniform
//! acoustic mode. Its interval mean pressure and volume flow obey the SAME
//! current DG law p = H + Z Q; only the incoming NSF characteristic is solved.
//! The NSF lattice step, thermochemical owner and DG propagator are unchanged.
//! Both candidates remain private until pressure and acoustic-work checks pass.
//!
//! This bounded, zero-reference-flow interface is an acoustic projection, not
//! full mass/entropy/enthalpy transport into a compressible duct. The actual NSF
//! mass and total-energy flux remain separate. The difference between integral
//! p' u and mean(p') integral(u) is retained as signed projection covariance,
//! never relabeled as heat or silently credited to make a balance close.
//! Face location, transverse modes, time resolution and the original first-order
//! NSF boundary still require independent refinement; no jet qualification or
//! reflection guarantee follows from this composition.

use crate::acoustic_realize::AcousticRealizeError;
use crate::dg_port::{DgPortLoad, DgTerminalFrame};
use fs_duct::dg::{DgPortResponse, DgStepOutput};
use fs_exec::CancelGate;
use fs_lbm::lattice3::Lattice3;
use fs_lbm::nsf3::{
    EnergyLedger, NsfBoundaryDrive, NsfBoundaryExchange, NsfCharacteristicPatch,
    NsfCharacteristicTrial, NsfError, NsfGrid3, NsfOpenState,
};
use fs_material::gas::GasState;

/// Finite physical pressure playback through the shared observer/PCM interface.
pub mod performance;
mod upstream;
pub use upstream::NsfDgUpstreamFrame;
mod reservoir;
pub use reservoir::{NsfDgReservoirFrame, NsfDgReservoirRadiation};
mod finite_reservoir;
pub use finite_reservoir::{NsfFiniteReservoirFrame, NsfFiniteReservoirLimits};

/// Physical drive location and units remain explicit in the coupled source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum NsfDgDrive {
    /// Prescribed control on an otherwise unloaded original NSF face.
    Boundary(NsfBoundaryDrive),
    /// Ideal held pressure rise from a finite upstream DG load into the zone
    /// [Pa]. This is a series effort, not prescribed absolute face pressure.
    UpstreamPressureRise(f64),
    /// Absolute pressure [Pa] and temperature [K] at the existing transverse
    /// reservoir. Density comes from the original gas source; its declared
    /// ghost velocity is retained. Both finite acoustic loads, when present,
    /// solve with zero series effort. Otherwise the other face has J-=0.
    Reservoir { pressure_pa: f64, temperature_k: f64 },
    /// Explicit independent transverse source in the grid's authored order.
    /// A distinct aperture may simultaneously receive passive exterior pressure.
    ReservoirAperture { aperture: usize, pressure_pa: f64, temperature_k: f64 },
    /// The selected finite thermochemical chamber supplies its own held p/T.
    /// Original mass and total-energy exchange update its inventory on acceptance.
    FiniteReservoir { aperture: usize },
}
impl NsfDgDrive {
    pub fn is_finite(self) -> bool {
        match self {
            Self::Boundary(value) => value.is_finite(),
            Self::UpstreamPressureRise(value) => value.is_finite(),
            Self::Reservoir { pressure_pa, temperature_k }
            | Self::ReservoirAperture { pressure_pa, temperature_k, .. } =>
                pressure_pa.is_finite() && temperature_k.is_finite(),
            Self::FiniteReservoir { .. } => true,
        }
    }
}

/// Original owner failures and explicit coupling/refinement refusals.
#[derive(Debug)]
pub enum NsfDgError {
    /// Invalid source, interface geometry, clock or numerical allowance.
    Input(&'static str),
    /// Original fluid material, boundary or step refusal.
    Flow(NsfError),
    /// Original acoustic material, port or step refusal.
    Acoustic(AcousticRealizeError),
    /// Original finite chamber source law or caloric inversion refusal.
    FiniteReservoir(fs_material::thermochemical_gas::energy::GasEnergyError),
    /// Cancellation publishes neither candidate.
    Cancelled,
    /// Caller-supplied incoming-characteristic bracket contains no sign change.
    NoBracket {
        /// Pressure residual at the caller's two bracket endpoints [Pa].
        residual_pa: [f64; 2],
    },
    /// Adjacent representable incoming characteristics still straddle a
    /// nonzero pressure residual gap. No trial is accepted at this interface.
    Resolution {
        /// Remaining incoming-characteristic interval [m/s].
        incoming_bracket_m_s: [f64; 2],
        /// Actual pressure residual at each remaining endpoint [Pa].
        residual_pa: [f64; 2],
    },
    /// Adjacent representable reservoir ghost pressures still straddle a
    /// nonzero pressure residual gap. No physical owner accepts this trial.
    PressureResolution {
        /// Remaining ghost-pressure offsets from the source reference [Pa].
        pressure_bracket_pa: [f64; 2],
        /// Actual effort residual at each remaining endpoint [Pa].
        residual_pa: [f64; 2],
    },
    /// The original owners did not meet the pressure budget within bounded work.
    Convergence {
        /// Last evaluated pressure mismatch [Pa].
        residual_pa: f64,
        /// Actual number of original candidate intervals evaluated.
        evaluations: usize,
    },
    /// The uniform/held-flow projection exceeds the caller's work allowance.
    Projection {
        /// Signed resolved work minus original DG inlet work [J].
        residual_j: f64,
        /// Caller-declared allowed magnitude [J].
        allowance_j: f64,
    },
}
impl core::fmt::Display for NsfDgError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "resolved NSF/DG coupling: {self:?}")
    }
}
impl std::error::Error for NsfDgError {}

/// Explicit work, horizon, root and acoustic-projection limits.
#[derive(Clone, Copy, Debug)]
pub struct NsfDgLimits {
    /// Original NSF steps per held DG interval.
    pub substeps: usize,
    /// Maximum actual grid cells, including solid cells in the retained snapshot.
    pub maximum_cells: usize,
    /// Maximum cell-steps across ALL trials of one coupled interval.
    pub maximum_cell_steps: usize,
    /// Root evaluations per interval, in 3..=256.
    pub maximum_evaluations: usize,
    /// Finite positive interval horizon.
    pub maximum_steps: u64,
    /// Absolute incoming-characteristic bracket [m/s], never an amplitude clamp.
    /// Both endpoints must remain inside the original NSF admission domain.
    pub incoming_bracket_m_s: [f64; 2],
    /// Absolute pressure-matching allowance [Pa].
    pub pressure_absolute_tolerance_pa: f64,
    /// Relative pressure allowance, in (0, 1e-3].
    pub pressure_relative_tolerance: f64,
    /// Absolute difference between resolved acoustic work and DG inlet work [J].
    pub projection_absolute_tolerance_j: f64,
    /// Relative work allowance against the sum of both absolute works, in [0, 1].
    /// Choosing a loose bound does not qualify the interface accuracy.
    pub projection_relative_tolerance: f64,
}

impl NsfDgLimits {
    /// Admit the complete trial work and numerical allowances before allocating
    /// a source volume. File and programmatic callers use this same rule; gas,
    /// area, topology and clock equality still require both physical owners.
    pub fn validate_for_cells(&self, cells: usize) -> Result<(), NsfDgError> {
        let work = cells.checked_mul(self.substeps)
            .and_then(|n| n.checked_mul(self.maximum_evaluations));
        if cells == 0 || cells > self.maximum_cells || self.maximum_cells == 0
            || !(1..=65_536).contains(&self.substeps)
            || !(3..=256).contains(&self.maximum_evaluations)
            || work.is_none_or(|n| n > self.maximum_cell_steps)
            || self.maximum_steps == 0
            || self.maximum_steps > (usize::MAX / self.substeps.max(1)) as u64
            || !self.incoming_bracket_m_s.iter().all(|x| x.is_finite())
            || self.incoming_bracket_m_s[0] >= self.incoming_bracket_m_s[1]
            || !self.pressure_absolute_tolerance_pa.is_finite()
            || self.pressure_absolute_tolerance_pa < 0.0
            || !self.pressure_relative_tolerance.is_finite()
            || !(0.0 < self.pressure_relative_tolerance && self.pressure_relative_tolerance <= 1e-3)
            || !self.projection_absolute_tolerance_j.is_finite()
            || self.projection_absolute_tolerance_j < 0.0
            || !self.projection_relative_tolerance.is_finite()
            || !(0.0..=1.0).contains(&self.projection_relative_tolerance)
        {
            return Err(NsfDgError::Input("resolved coupling requires finite complete cell-work, horizon, root and projection budgets"));
        }
        Ok(())
    }
}

/// One accepted resolved/acoustic interval, with no duplicated excitation owner.
#[derive(Debug, Clone)]
pub struct NsfDgFrame {
    /// One-based accepted DG interval.
    pub step: u64,
    /// Number of original NSF candidate intervals evaluated.
    pub evaluations: usize,
    /// Accepted incoming characteristic on the coupled face [m/s].
    pub incoming_characteristic_m_s: f64,
    /// Actual prescribed physical control, including its location and units.
    pub drive: NsfDgDrive,
    /// Original finite upstream load and source-work accounts, when attached.
    pub upstream: Option<NsfDgUpstreamFrame>,
    /// Original NSF face accounts: low, high. Flux/work are summed over all
    /// substeps; pressure and continuum flow are their interval means.
    pub faces: [NsfBoundaryExchange; 2],
    /// Selected aperture exchange when only part of the downstream face is
    /// coupled. The complete face accounts above also include the independently
    /// driven surrounding opening; that escaping flow never enters the bore.
    pub coupled_patch: Option<NsfBoundaryExchange>,
    /// Selected source/upstream aperture exchange, if that face is patched.
    /// Its surrounding opening remains in `faces` but never enters the finite
    /// upstream duct or its series-source work. This is already part of the
    /// complete face account and must not be added to the fluid totals again.
    pub source_patch: Option<NsfBoundaryExchange>,
    /// Original streamed mass leaving all transverse reservoir apertures [kg].
    /// This is separate from both characteristic faces; zero when absent.
    pub reservoir_mass_out_kg: f64,
    /// Original per-tile reservoir continuum volume flow, averaged over all
    /// original substeps. Kept separate from compressible streamed mass.
    pub reservoir_flows: Vec<fs_lbm::nsf3::NsfReservoirFlow>,
    /// Separate original conservative interval accounts in authored owner order.
    pub reservoir_exchanges: Vec<fs_lbm::nsf3::NsfReservoirExchange>,
    /// Live self/mutual exterior reaction at the independent escape opening.
    pub reservoir_radiation: Option<NsfDgReservoirFrame>,
    /// Finite source storage, accepting the original conservative exchange once.
    pub finite_reservoir: Option<NsfFiniteReservoirFrame>,
    /// Original NSF energy accounts, including reservoir energy in `open_energy_j`.
    /// Mixed-boundary unresolved exchange retains its original classification;
    /// it is never relabeled as resolved viscous heat by this adapter.
    pub flow: EnergyLedger,
    /// Actual change of original NSF mass [kg], independent of the face sum.
    pub flow_mass_change_kg: f64,
    /// Actual change of original NSF sensible-chart total energy [J].
    pub flow_energy_change_j: f64,
    /// Original aggregate DG/branch accounts, retaining their inferred-loss class.
    pub duct: DgStepOutput,
    /// Original downstream terminal order, suitable for its bound microphone.
    pub terminals: Vec<DgTerminalFrame>,
    /// Held continuum volume flow into the DG duct [m³/s].
    pub inlet_flow_m3_s: f64,
    /// Actual original DG step-mean inlet pressure [Pa].
    pub inlet_pressure_pa: f64,
    /// Resolved mean acoustic pressure minus actual DG inlet pressure [Pa].
    pub pressure_residual_pa: f64,
    /// Resolved integral(p' u) minus mean(p') mean(Q) dt [J].
    /// Contains transverse and temporal covariance, not physical heat.
    pub projection_covariance_j: f64,
    /// Resolved acoustic work minus original DG integrated inlet work [J].
    pub acoustic_work_residual_j: f64,
    /// Original DG inlet work minus its mean-pressure/held-flow product [J].
    pub port_work_residual_j: f64,
    /// Actual lattice mass flux minus rho_ref times continuum swept volume [kg].
    /// Includes finite-amplitude and discretization differences; not a correction.
    pub mass_projection_difference_kg: f64,
}

/// The only accepted flow and acoustic owners. An initial resolved pulse or
/// interior flow is retained, while both clocks must be cold at binding.
pub struct NsfDg<L: Lattice3<27>> {
    zone: NsfGrid3<L>,
    load: DgPortLoad,
    limits: NsfDgLimits,
    high: bool,
    reference: GasState,
    area_m2: f64,
    steps: u64,
    upstream: Option<upstream::UpstreamLoad>,
    patch_background: [Option<NsfBoundaryDrive>; 2],
    reservoir_radiation: Option<reservoir::ReservoirRadiation>,
    finite_reservoir: Option<finite_reservoir::FiniteReservoir>,
}

fn poll(gate: &CancelGate) -> Result<(), NsfDgError> {
    if gate.is_requested() { Err(NsfDgError::Cancelled) } else { Ok(()) }
}
pub(super) fn same_clock(a: f64, b: f64) -> bool {
    a.is_finite() && b.is_finite() && a > 0.0 && b > 0.0
        && (a-b).abs() <= 8.0 * f64::EPSILON * a.abs().max(b.abs())
}

impl<L: Lattice3<27>> NsfDg<L> {
    /// Bind one complete face to the DG inlet. Areas agree to floating-point
    /// geometry roundoff; clocks differ by at most derived-clock roundoff.
    /// Both characteristic reference velocities must be zero. The complete gas
    /// at the coupled reference must equal the actual DG inlet gas.
    pub fn new(
        zone: NsfGrid3<L>, load: DgPortLoad, high: bool,
        limits: NsfDgLimits, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        Self::bind(zone, load, high, limits, gate, [false; 2])
    }

    /// Couple the grid's explicitly declared rectangular aperture patch to the
    /// bore. The remainder of that face retains its original boundary drive,
    /// independently of the solved aperture pressure and flow. The other whole
    /// face remains the pressure/characteristic source or finite upstream load.
    ///
    /// This permits resolved flow to split between a bore and its surroundings.
    /// Only the selected area enters the DG inlet. Complete fluid mass/energy
    /// flux and selected-patch acoustic work remain separately observable.
    /// The boundary is still the original plane-normal, first-order interface;
    /// aperture geometry and flow resolution require their own convergence.
    pub fn new_on_patch(
        zone: NsfGrid3<L>, load: DgPortLoad, high: bool,
        limits: NsfDgLimits, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        let mut selected = [false; 2];
        selected[usize::from(high)] = true;
        Self::bind(zone, load, high, limits, gate, selected)
    }

    /// Drive the original aperture on the OTHER face while the main DG bore
    /// receives the complete downstream face. Only the selected source cells
    /// receive each prescribed pressure/characteristic or finite upstream load.
    /// The surrounding cells retain their independently declared boundary drive.
    ///
    /// This constructor selects only the source patch and refuses an additional
    /// downstream patch. Use [`Self::new_on_patches`] to select both. A subsequent
    /// [`Self::with_upstream_load`] matches that duct to the selected source area,
    /// not the whole face; all three original candidates are committed together.
    pub fn new_with_source_patch(
        zone: NsfGrid3<L>, load: DgPortLoad, high: bool,
        limits: NsfDgLimits, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        let mut selected = [false; 2];
        selected[usize::from(!high)] = true;
        Self::bind(zone, load, high, limits, gate, selected)
    }

    /// Couple independent finite source and bore openings on opposite faces.
    /// Each face must carry its own original characteristic patch. The main
    /// bore receives only its selected downstream area; a prescribed source
    /// or finite upstream duct receives only the other selected area. Both
    /// remainders retain their independent boundary laws and full-face fluxes.
    /// The same original NSF/DG solve and work tolerances apply to both patches.
    pub fn new_on_patches(
        zone: NsfGrid3<L>, load: DgPortLoad, high: bool,
        limits: NsfDgLimits, gate: &CancelGate,
    ) -> Result<Self, NsfDgError> {
        Self::bind(zone, load, high, limits, gate, [true; 2])
    }

    fn bind(
        zone: NsfGrid3<L>, load: DgPortLoad, high: bool,
        limits: NsfDgLimits, gate: &CancelGate, selected_patches: [bool; 2],
    ) -> Result<Self, NsfDgError> {
        poll(gate)?;
        if limits.maximum_cells == 0 || zone.cell_count() > limits.maximum_cells {
            return Err(NsfDgError::Input("resolved grid exceeds the declared retained-cell budget"));
        }
        let open = zone.characteristic_axis().ok_or(NsfDgError::Input(
            "resolved coupling requires the original characteristic axis"))?;
        let reference = zone.characteristic_reference_gas(high).map_err(NsfDgError::Flow)?;
        let drives = zone.boundary_drives().ok_or(NsfDgError::Input(
            "aperture coupling requires the original surrounding boundary drives"))?;
        let mut patch_background = [None; 2];
        for side in 0..2 {
            if selected_patches[side] != zone.characteristic_patch_on(side != 0).is_some() {
                return Err(NsfDgError::Input(
                    "select every declared source or coupled aperture explicitly; a whole-face role cannot hide a patch"));
            }
            if selected_patches[side] { patch_background[side] = Some(drives[side]); }
        }
        let traces = if selected_patches[usize::from(high)] { zone.open_patch_traces_on(high) }
            else { zone.open_boundary_traces(high) }
            .ok_or(NsfDgError::Input("resolved coupling requires original interface quadrature"))?;
        // The original Cartesian face has one identical area per trace. A
        // count/product avoids resolution-dependent accumulated-sum roundoff.
        let area = traces.first().ok_or(NsfDgError::Input(
            "resolved coupling requires a nonempty characteristic face"))?.area_m2
            * traces.len() as f64;
        let duct_area = load.duct().inlet_area_m2();
        limits.validate_for_cells(zone.cell_count())?;
        if zone.steps() != 0 || load.state().step_count() != 0 || load.has_driven_mouth()
            || open.faces.iter().any(|face| face.reference_velocity_m_s != [0.0; 3])
            || load.duct().inlet_gas() != &reference
            || !area.is_finite() || area <= 0.0 || !duct_area.is_finite()
            || (area-duct_area).abs() > 32.0 * f64::EPSILON * area.max(duct_area)
            || !same_clock(zone.law().dt_s * limits.substeps as f64, load.propagator().dt())
        {
            return Err(NsfDgError::Input("NSF/DG requires cold matching gas, area and clock, zero reference flow and finite complete work/projection budgets"));
        }
        Ok(Self { zone, load, limits, high, reference, area_m2: area, steps: 0,
            upstream: None, patch_background, reservoir_radiation: None,
            finite_reservoir: None })
    }
    /// Original resolved fluid state and thermochemical source.
    pub const fn zone(&self) -> &NsfGrid3<L> { &self.zone }
    /// Original acoustic state, material walls and exterior loads.
    pub const fn load(&self) -> &DgPortLoad { &self.load }
    /// Complete reference material state used by this interface.
    pub const fn reference_gas(&self) -> &GasState { &self.reference }
    /// Actual coupled area [m²]: the selected aperture or complete face.
    pub const fn area_m2(&self) -> f64 { self.area_m2 }
    /// Physical face carrying the downstream DG inlet; the source is opposite.
    pub const fn coupled_face_is_high(&self) -> bool { self.high }
    /// Original downstream aperture bounds, if explicitly selected for the bore.
    pub fn coupled_patch(&self) -> Option<NsfCharacteristicPatch> {
        self.zone.characteristic_patch_on(self.high)
    }
    /// Original source aperture bounds, if explicitly selected on the other face.
    pub fn source_patch(&self) -> Option<NsfCharacteristicPatch> {
        self.zone.characteristic_patch_on(!self.high)
    }
    /// Original remainder drive when exactly one face is patched. Returns
    /// `None` for zero or two patches; use the per-face accessor in that case.
    pub const fn patch_background_drive(&self) -> Option<NsfBoundaryDrive> {
        match self.patch_background {
            [Some(drive), None] | [None, Some(drive)] => Some(drive),
            _ => None,
        }
    }
    /// Independent original boundary law on this physical face's remainder.
    pub const fn patch_background_drive_on(&self, high: bool) -> Option<NsfBoundaryDrive> {
        self.patch_background[high as usize]
    }
    fn patch_boundaries(&self) -> [Option<NsfBoundaryDrive>; 2] {
        self.patch_background
    }
    /// Original explicit numerical limits.
    pub const fn limits(&self) -> &NsfDgLimits { &self.limits }
    /// Accepted coupled intervals.
    pub const fn accepted_steps(&self) -> u64 { self.steps }

    /// Advance the explicitly located source against its retained participants.
    pub fn step_with_drive(&mut self, drive: NsfDgDrive, gate: &CancelGate)
        -> Result<NsfDgFrame, NsfDgError>
    {
        self.admit_finite_reservoir_drive(drive)?;
        match drive {
            NsfDgDrive::Boundary(value) => self.step_with_boundary(value, gate),
            NsfDgDrive::UpstreamPressureRise(value) => self.step_with_upstream_pressure_rise(value, gate),
            NsfDgDrive::Reservoir { pressure_pa, temperature_k } =>
                self.step_with_reservoir(pressure_pa, temperature_k, gate),
            NsfDgDrive::ReservoirAperture { aperture, pressure_pa, temperature_k } =>
                self.step_with_reservoir_aperture(aperture, pressure_pa, temperature_k, gate),
            NsfDgDrive::FiniteReservoir { aperture } =>
                self.step_with_finite_reservoir(aperture, gate),
        }
    }

    /// Blow through the transverse aperture independently of the two axial
    /// acoustic loads. Held absolute pressure and temperature select the SAME
    /// material-owned reservoir state used by populations, force and energy.
    /// No jet speed, source waveform or extra pressure rise is synthesized.
    ///
    /// Both ends retain finite DG feedback when an upstream load is attached;
    /// otherwise the unused end receives its unperturbed incoming acoustic
    /// characteristic. Reservoir conditions publish only with the complete
    /// accepted fluid and acoustic interval, including on failed root solves.
    pub fn step_with_reservoir(
        &mut self, pressure_pa: f64, temperature_k: f64, gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        poll(gate)?;
        let state = self.zone.reservoir_aperture_state_at(pressure_pa, temperature_k)
            .map_err(NsfDgError::Flow)?;
        let drive = NsfDgDrive::Reservoir { pressure_pa, temperature_k };
        if self.upstream.is_some() {
            self.step_with_upstream_conditions(0.0, &[(0, state)], drive, gate)
        } else {
            self.step_with_boundary_conditions(NsfBoundaryDrive::IncomingCharacteristic(0.0),
                &[(0, state)], drive, gate)
        }
    }

    /// Drive one original transverse aperture while the axial loads and any
    /// distinct passive escape opening solve in the same accepted interval.
    pub fn step_with_reservoir_aperture(
        &mut self, aperture: usize, pressure_pa: f64, temperature_k: f64, gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        poll(gate)?;
        let state = self.zone.reservoir_aperture_state_for(aperture, pressure_pa, temperature_k)
            .map_err(NsfDgError::Flow)?;
        let drive = NsfDgDrive::ReservoirAperture { aperture, pressure_pa, temperature_k };
        if self.upstream.is_some() {
            self.step_with_upstream_conditions(0.0, &[(aperture, state)], drive, gate)
        } else {
            self.step_with_boundary_conditions(NsfBoundaryDrive::IncomingCharacteristic(0.0),
                &[(aperture, state)], drive, gate)
        }
    }

    /// Hold the OTHER face's incoming characteristic for one interval [m/s].
    /// The coupled face is solved inside the caller's bracket. No synthetic jet
    /// source, output gain, limiter or lagged reflection is added.
    pub fn step(&mut self, external_incoming_m_s: f64, gate: &CancelGate)
        -> Result<NsfDgFrame, NsfDgError>
    {
        self.step_with_boundary(NsfBoundaryDrive::IncomingCharacteristic(external_incoming_m_s), gate)
    }

    /// Advance against a physical incoming-wave or prescribed-pressure source.
    ///
    /// Pressure is imposed at every original NSF boundary reconstruction on
    /// the declared source patch, or the complete source face when unpatched. The
    /// returning wave determines source-face flow, while the same scalar solve
    /// matches the opposite face to the live DG bore. No independent density,
    /// sound speed, source flow law or pressure-to-wave calibration is added.
    /// The source is an ideal pressure boundary, not finite airway storage.
    pub fn step_with_boundary(&mut self, external_boundary: NsfBoundaryDrive, gate: &CancelGate)
        -> Result<NsfDgFrame, NsfDgError>
    {
        self.step_with_boundary_conditions(external_boundary, &[],
            NsfDgDrive::Boundary(external_boundary), gate)
    }

    fn step_with_boundary_conditions(
        &mut self, external_boundary: NsfBoundaryDrive, reservoir_states: &[(usize, NsfOpenState)],
        control: NsfDgDrive, gate: &CancelGate,
    ) -> Result<NsfDgFrame, NsfDgError> {
        poll(gate)?;
        self.admit_finite_reservoir_drive(control)?;
        if self.reservoir_radiation.is_some() {
            return self.step_with_reservoir_drive(control, gate);
        }
        if self.upstream.is_some() {
            return Err(NsfDgError::Input("a finite upstream load requires an explicit series pressure-rise control"));
        }
        if !external_boundary.is_finite() || self.steps >= self.limits.maximum_steps {
            return Err(NsfDgError::Input("resolved interface input is nonfinite or its horizon is exhausted"));
        }
        let relation = self.load.port_response().map_err(NsfDgError::Acoustic)?;
        let limits = self.limits;
        let side = usize::from(self.high);
        let patch = self.patch_boundaries();
        let mut bracket = limits.incoming_bracket_m_s;
        let mut values = [0.0; 2];
        let mut incoming = bracket[0];
        let mut last_residual = 0.0;
        let mut evaluations = 0;
        for evaluation in 1..=limits.maximum_evaluations {
            poll(gate)?;
            evaluations = evaluation;
            let mut drive = [external_boundary; 2];
            drive[side] = NsfBoundaryDrive::IncomingCharacteristic(incoming);
            let candidate = preview_interval(&mut self.zone, drive, patch,
                reservoir_states, limits, gate)?;
            let face = coupled_exchange(&candidate, side)?;
            let flow = face.normal_volume_flow_m3_s;
            let residual = face.mean_acoustic_pressure_pa
                - (relation.history_pressure_pa + relation.impedance_pa_s_m3 * flow);
            let scale = face.mean_acoustic_pressure_pa.abs() + relation.history_pressure_pa.abs()
                + (relation.impedance_pa_s_m3 * flow).abs();
            if !residual.is_finite() || !scale.is_finite() {
                return Err(NsfDgError::Input("resolved interface pressure equation overflowed"));
            }
            let allowance = limits.pressure_absolute_tolerance_pa
                + limits.pressure_relative_tolerance * scale;
            if residual.abs() <= allowance {
                let frame = accept_candidate(candidate, &mut self.load, &mut self.finite_reservoir, side, relation,
                    self.reference, limits, incoming, control,
                    evaluation, self.steps + 1, gate)?;
                self.steps += 1;
                return Ok(frame);
            }
            last_residual = residual;
            drop(candidate);
            if evaluation == 1 {
                values[0] = residual;
                incoming = bracket[1];
                continue;
            }
            if evaluation == 2 {
                values[1] = residual;
                if values[0].is_sign_positive() == values[1].is_sign_positive() {
                    return Err(NsfDgError::NoBracket { residual_pa: values });
                }
            } else {
                let side = usize::from(residual.is_sign_positive() != values[0].is_sign_positive());
                bracket[side] = incoming;
                values[side] = residual;
            }
            // A safeguarded secant proposal accelerates the nearly linear
            // characteristic law. Only a strict interior point is evaluated;
            // the physical amplitude is never clamped to a bracket endpoint.
            let proposed = bracket[0] - values[0]
                * ((bracket[1]-bracket[0]) / (values[1]-values[0]));
            incoming = if proposed.is_finite()
                && proposed > bracket[0] + 0.1*(bracket[1]-bracket[0])
                && proposed < bracket[1] - 0.1*(bracket[1]-bracket[0])
            { proposed } else { f64::midpoint(bracket[0], bracket[1]) };
            if incoming <= bracket[0] || incoming >= bracket[1] {
                return Err(NsfDgError::Resolution {
                    incoming_bracket_m_s: bracket,
                    residual_pa: values,
                });
            }
        }
        Err(NsfDgError::Convergence { residual_pa: last_residual,
            evaluations })
    }
}

fn preview_interval<'a, L: Lattice3<27>>(
    zone: &'a mut NsfGrid3<L>, mut drives: [NsfBoundaryDrive; 2],
    patches: [Option<NsfBoundaryDrive>; 2], reservoir_states: &[(usize, NsfOpenState)],
    limits: NsfDgLimits, gate: &CancelGate,
) -> Result<NsfCharacteristicTrial<'a, L>, NsfDgError> {
    let mut patch_drives = [None; 2];
    for side in 0..2 {
        if let Some(background) = patches[side] {
            patch_drives[side] = Some(drives[side]);
            drives[side] = background;
        }
    }
    let result = if !reservoir_states.is_empty() {
        zone.preview_boundary_interval_with_patches_and_reservoirs(drives, patch_drives, reservoir_states,
            limits.substeps, limits.maximum_cell_steps, || gate.is_requested())
    } else {
        zone.preview_boundary_interval_with_patches(drives, patch_drives, limits.substeps,
            limits.maximum_cell_steps, || gate.is_requested())
    };
    result.map_err(|error| if gate.is_requested() { NsfDgError::Cancelled }
        else { NsfDgError::Flow(error) })
}

fn coupled_exchange<L: Lattice3<27>>(
    trial: &NsfCharacteristicTrial<'_, L>, side: usize,
) -> Result<NsfBoundaryExchange, NsfDgError> {
    if trial.grid().characteristic_patch_on(side != 0).is_some() {
        trial.patch_exchange_on(side != 0).copied().ok_or(NsfDgError::Input(
            "original aperture exchange disappeared during the coupled trial"))
    } else { Ok(trial.boundary_exchange()[side]) }
}

#[allow(clippy::too_many_arguments)]
fn accept_candidate<L: Lattice3<27>>(
    fluid: NsfCharacteristicTrial<'_, L>, load: &mut DgPortLoad,
    finite_owner: &mut Option<finite_reservoir::FiniteReservoir>, side: usize,
    relation: DgPortResponse, reference: GasState, limits: NsfDgLimits,
    incoming: f64, control: NsfDgDrive,
    evaluations: usize, step: u64, gate: &CancelGate,
) -> Result<NsfDgFrame, NsfDgError> {
    let face = coupled_exchange(&fluid, side)?;
    let q = face.normal_volume_flow_m3_s;
    let dt = load.propagator().dt();
    let inlet = load.inlet();
    let acoustic = load.preview_step(q).map_err(NsfDgError::Acoustic)?;
    let p = acoustic.output().port_outputs.iter().find(|port| port.id == inlet)
        .ok_or(NsfDgError::Input("original DG inlet observation disappeared"))?.conjugate_rate;
    let pressure_residual = face.mean_acoustic_pressure_pa - p;
    let scale = face.mean_acoustic_pressure_pa.abs() + p.abs()
        + relation.history_pressure_pa.abs() + (relation.impedance_pa_s_m3 * q).abs();
    let inlet_work = acoustic.inlet_work_j();
    let port_defect = inlet_work - p * (q * dt);
    let covariance = face.acoustic_pressure_work_out_j
        - face.mean_acoustic_pressure_pa * (q * dt);
    let work_residual = face.acoustic_pressure_work_out_j - inlet_work;
    let work_scale = face.acoustic_pressure_work_out_j.abs() + inlet_work.abs();
    let allowance = limits.projection_absolute_tolerance_j
        + limits.projection_relative_tolerance * work_scale;
    if ![p, pressure_residual, scale, inlet_work, port_defect, covariance,
        work_residual, work_scale, allowance].iter().all(|x| x.is_finite())
        || pressure_residual.abs() > limits.pressure_absolute_tolerance_pa
            + limits.pressure_relative_tolerance * scale
        || port_defect.abs() > 1e-10 * (inlet_work.abs() + scale*q.abs()*dt).max(1e-300)
    {
        return Err(NsfDgError::Input("original mean pressure or integrated port work disagrees at the resolved interface"));
    }
    if work_residual.abs() > allowance {
        return Err(NsfDgError::Projection { residual_j: work_residual, allowance_j: allowance });
    }
    let faces = *fluid.boundary_exchange();
    let coupled_patch = fluid.patch_exchange_on(side != 0).copied();
    let source_patch = fluid.patch_exchange_on(side == 0).copied();
    let reservoir_mass_out_kg = fluid.reservoir_mass_out_kg();
    let reservoir_flows = fluid.reservoir_flows().to_vec();
    let reservoir_exchanges = fluid.reservoir_exchanges().to_vec();
    let flow = *fluid.energy_ledger();
    let mass_change = fluid.grid().total_mass_kg() - fluid.previous_mass_kg();
    let energy_change = fluid.grid().total_energy_j() - fluid.previous_energy_j();
    let mass_projection = face.mass_out_kg - reference.density * (q * dt);
    let terminals = acoustic.terminal_frames().to_vec();
    if ![mass_change, energy_change, mass_projection, reservoir_mass_out_kg].iter().all(|x| x.is_finite()) {
        return Err(NsfDgError::Input("resolved interval mass or energy difference overflowed"));
    }
    let finite_candidate = finite_reservoir::prepare(finite_owner.as_ref(), &fluid, gate)?;
    poll(gate)?;
    // All source, allocation, pressure, work and cancellation checks precede
    // publication of either original owner.
    fluid.accept();
    let duct = acoustic.accept();
    let finite_reservoir = finite_reservoir::publish(finite_owner, finite_candidate);
    Ok(NsfDgFrame {
        step, evaluations, incoming_characteristic_m_s: incoming,
        drive: control, upstream: None, faces, coupled_patch, source_patch,
        reservoir_mass_out_kg, reservoir_flows, reservoir_exchanges, reservoir_radiation: None,
        finite_reservoir, flow,
        flow_mass_change_kg: mass_change, flow_energy_change_j: energy_change, duct, terminals,
        inlet_flow_m3_s: q, inlet_pressure_pa: p, pressure_residual_pa: pressure_residual,
        projection_covariance_j: covariance, acoustic_work_residual_j: work_residual,
        port_work_residual_j: port_defect, mass_projection_difference_kg: mass_projection,
    })
}

#[cfg(test)]
mod tests;
