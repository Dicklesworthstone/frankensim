//! Physical solid-conductor electrodes with retained 3D vector hysteresis.
//!
//! Voltage or total-current control enters the original nonlinear A–phi solve.
//! The actual current distribution, electrode reaction, pinning heat and Joule
//! heat are recovered together with magnetic storage and endpoint numerical
//! loss. Every Newton trial retains the SAME accepted vector memories; changing
//! electrode causality changes only the charge space, never A or its history.
//!
//! Electrodes use the shared exterior-face and connected fixed-A trace admission.
//! Geometry, magnetic coefficients, material frames and occupancy stay fixed.
//! Original conductivity sources may refresh at admitted temperatures while
//! preserving electrode support and source identity. No uniform-current model,
//! fitted terminal impedance, contact resistance or moving memory is introduced.

use super::super::terminals::TerminalLayout;
/// Native RLC circuits driving physical electrodes with retained vector memory.
pub mod circuit;

use super::*;

fn magnetic_entries(
    field: &SpatialHystereticInduction<'_, impl SpatialConductorMaterial>,
) -> Result<usize> {
    field
        .problem()
        .mesh()
        .tets
        .len()
        .checked_mul(36)
        .ok_or_else(|| invalid("magnetic electrode assembly bound overflows").into())
}

fn with_coupling<'m, C: SpatialConductorMaterial>(
    field: SpatialHystereticInduction<'m, C>,
    coupling: Coupling<C>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialHystereticInduction<'m, C>> {
    // Electrode regrouping changes the coupled dimension. Rebuild only the
    // stationary zero operators; all accepted cochains and memories remain.
    let n = coupling.mass.nrows();
    let stationary = Arc::new(Stationary {
        motion: Csr::from_parts(n, n, vec![0; n + 1], vec![], vec![]),
        drive: vec![0.0; n],
        velocity: field.prepared.stationary.velocity.clone(),
    });
    let prepared = Arc::new(Prepared {
        coupling,
        stationary,
    });
    poll(cancelled)?;
    Ok(SpatialHystereticInduction { prepared, ..field })
}

/// Accepted vector history with prescribed physical electrode voltages.
/// The underlying isolated field cannot be stepped through this grounded owner.
#[derive(Clone)]
pub struct SpatialHystereticTerminalInduction<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    field: SpatialHystereticInduction<'m, C>,
    terminals: Arc<TerminalLayout>,
}

/// Accepted vector history with one solved potential per total-current electrode.
#[derive(Clone)]
pub struct SpatialHystereticCurrentInduction<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    field: SpatialHystereticInduction<'m, C>,
    terminals: Arc<TerminalLayout>,
    terminal_unknown: Arc<[Option<usize>]>,
}

pub struct SpatialHystereticTerminalInductionCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    state: SpatialHystereticTerminalInduction<'m, C>,
    step: SpatialHystereticInductionStep,
}
pub struct SpatialHystereticCurrentInductionCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    state: SpatialHystereticCurrentInduction<'m, C>,
    step: SpatialHystereticInductionStep,
}
impl<'m, C: SpatialConductorMaterial> SpatialHystereticTerminalInductionCandidate<'m, C> {
    pub fn state(&self) -> &SpatialHystereticTerminalInduction<'m, C> {
        &self.state
    }
    pub fn step(&self) -> &SpatialHystereticInductionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialHystereticTerminalInduction<'m, C>,
        SpatialHystereticInductionStep,
    ) {
        (self.state, self.step)
    }
}
impl<'m, C: SpatialConductorMaterial> SpatialHystereticCurrentInductionCandidate<'m, C> {
    pub fn state(&self) -> &SpatialHystereticCurrentInduction<'m, C> {
        &self.state
    }
    pub fn step(&self) -> &SpatialHystereticInductionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialHystereticCurrentInduction<'m, C>,
        SpatialHystereticInductionStep,
    ) {
        (self.state, self.step)
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialHystereticTerminalInduction<'m, C> {
    /// Attach disjoint physical electrodes without re-equilibrating or resetting
    /// accepted magnetic memory, impressed source, edge cochains or clock.
    pub fn new(
        field: SpatialHystereticInduction<'m, C>,
        electrode_faces: &[Vec<usize>],
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let (terminals, coupling) = TerminalLayout::prepare(
            field.problem().geometry(),
            &field.prepared.coupling,
            electrode_faces,
            magnetic_entries(&field)?,
            limits,
            cancelled,
        )?;
        Ok(Self {
            field: with_coupling(field, coupling, cancelled)?,
            terminals,
        })
    }

    pub fn problem(&self) -> &SpatialHystereticProblem {
        self.field.problem()
    }
    pub fn potential_wb(&self) -> &[f64] {
        self.field.potential_wb()
    }
    pub fn cell_states(&self) -> &[VectorHysteresisState] {
        self.field.cell_states()
    }
    pub fn time_s(&self) -> f64 {
        self.field.time_s()
    }
    pub fn conductors(&self) -> &[Option<C>] {
        self.field.conductors()
    }
    pub fn electrode_faces(&self) -> &[Vec<usize>] {
        self.terminals.faces()
    }
    pub fn field(&self, cancelled: &mut impl FnMut() -> bool) -> Result<SpatialInductionField> {
        self.field.field(cancelled)
    }

    /// Replace prescribed potentials by one grouped unknown per non-reference
    /// electrode, retaining the original vector and scalar material histories.
    pub fn into_current_driven(
        self,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticCurrentInduction<'m, C>> {
        let (coupling, terminal_unknown) = self.terminals.current_coupling(
            self.field.problem().geometry(),
            &self.field.prepared.coupling,
            magnetic_entries(&self.field)?,
            limits,
            cancelled,
        )?;
        Ok(SpatialHystereticCurrentInduction {
            field: with_coupling(self.field, coupling, cancelled)?,
            terminals: self.terminals,
            terminal_unknown: terminal_unknown.into(),
        })
    }

    /// Backward Euler uses ENDPOINT voltages; native nodal reactions determine
    /// the actual terminal currents and independently reconstructed source work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_candidate(
        &self,
        endpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticTerminalInductionCandidate<'m, C>> {
        self.step_prepared(
            Arc::clone(&self.field.prepared),
            endpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Refresh original conductivity without changing grouped electrode support
    /// or magnetic memory. Fixed magnetic fits retain their admitted domains.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_with_conductors_candidate(
        &self,
        endpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticTerminalInductionCandidate<'m, C>> {
        let prepared = self.field.prepare_conductors(conductors, cancelled)?;
        self.step_prepared(
            prepared,
            endpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &self,
        prepared: Arc<Prepared<C>>,
        voltage_v: &[f64],
        current_a_m2: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticTerminalInductionCandidate<'m, C>> {
        let drive = self.terminals.voltage_drive(voltage_v, cancelled)?;
        let (field, step) = self
            .field
            .step_driven_candidate(
                prepared,
                current_a_m2,
                duration_s,
                controls,
                Some(&drive),
                cancelled,
            )?
            .into_parts();
        Ok(SpatialHystereticTerminalInductionCandidate {
            state: Self {
                field,
                terminals: Arc::clone(&self.terminals),
            },
            step,
        })
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialHystereticCurrentInduction<'m, C> {
    pub fn problem(&self) -> &SpatialHystereticProblem {
        self.field.problem()
    }
    pub fn potential_wb(&self) -> &[f64] {
        self.field.potential_wb()
    }
    pub fn cell_states(&self) -> &[VectorHysteresisState] {
        self.field.cell_states()
    }
    pub fn time_s(&self) -> f64 {
        self.field.time_s()
    }
    pub fn conductors(&self) -> &[Option<C>] {
        self.field.conductors()
    }
    pub fn electrode_faces(&self) -> &[Vec<usize>] {
        self.terminals.faces()
    }
    pub fn field(&self, cancelled: &mut impl FnMut() -> bool) -> Result<SpatialInductionField> {
        self.field.field(cancelled)
    }

    pub fn into_voltage_driven(
        self,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticTerminalInduction<'m, C>> {
        let coupling = self.terminals.voltage_coupling(
            self.field.problem().geometry(),
            &self.field.prepared.coupling,
            magnetic_entries(&self.field)?,
            limits,
            cancelled,
        )?;
        Ok(SpatialHystereticTerminalInduction {
            field: with_coupling(self.field, coupling, cancelled)?,
            terminals: self.terminals,
        })
    }

    /// Prescribe ENDPOINT total currents, positive into the conductor. Each
    /// conductor balances independently; the nonuniform current and electrode
    /// voltages are outputs of the same retained-memory field solve.
    #[allow(clippy::too_many_arguments)]
    pub fn step_current_candidate(
        &self,
        endpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticCurrentInductionCandidate<'m, C>> {
        self.step_prepared(
            Arc::clone(&self.field.prepared),
            endpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_current_with_conductors_candidate(
        &self,
        endpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticCurrentInductionCandidate<'m, C>> {
        let prepared = self.field.prepare_conductors(conductors, cancelled)?;
        self.step_prepared(
            prepared,
            endpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &self,
        prepared: Arc<Prepared<C>>,
        incoming_current_a: &[f64],
        current_a_m2: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialHystereticInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialHystereticCurrentInductionCandidate<'m, C>> {
        let drive = self.terminals.current_drive(
            incoming_current_a,
            &self.terminal_unknown,
            controls.induction.induction.current_tolerance_a,
            cancelled,
        )?;
        let (field, step) = self
            .field
            .step_driven_candidate(
                prepared,
                current_a_m2,
                duration_s,
                controls,
                Some(&drive),
                cancelled,
            )?
            .into_parts();
        Ok(SpatialHystereticCurrentInductionCandidate {
            state: Self {
                field,
                terminals: Arc::clone(&self.terminals),
                terminal_unknown: Arc::clone(&self.terminal_unknown),
            },
            step,
        })
    }
}
