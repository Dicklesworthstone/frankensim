//! Physical solid-conductor electrodes with reversible 3D magnetic saturation.
//!
//! Electrode potentials or total currents drive the original coupled Whitney
//! A–phi equations. The same resolved current supplies Ampere, charge continuity,
//! endpoint electrode work, conductor-frame Joule heat and Lorentz work. No
//! stranded profile, sampled response matrix or fitted impedance replaces that
//! current distribution. Magnetic saturation and diffusion retain accepted A.
//!
//! Electrode geometry and its connected prescribed-A reference trace are fixed.
//! Voltage/current causality can change without projecting field cochains or
//! restarting time. Source refresh retains the original charge groups, material
//! identity and conductivity basis. Temperature-dependent magnetic refresh
//! reports parameter work separately from physical heat and numerical loss.
//! External lead fields, contact impedance and moving interfaces are not added.

/// Native RLC circuit closure of these physical electrode currents.
pub mod circuit;

use super::super::terminals::TerminalLayout;
use super::*;

fn magnetic_entries(
    field: &SpatialNonlinearInduction<'_, impl SpatialConductorMaterial>,
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
    field: SpatialNonlinearInduction<'m, C>,
    coupling: Coupling<C>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<SpatialNonlinearInduction<'m, C>> {
    let prepared = Arc::new(Prepared {
        coupling,
        permeability_bound: field.prepared.permeability_bound.clone(),
    });
    poll(cancelled)?;
    Ok(SpatialNonlinearInduction { prepared, ..field })
}

/// Accepted nonlinear field with prescribed equipotential electrode voltages.
/// The isolated field owner stays private because its charge space is grounded.
#[derive(Clone)]
pub struct SpatialNonlinearTerminalInduction<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    field: SpatialNonlinearInduction<'m, C>,
    terminals: Arc<TerminalLayout>,
}

/// Accepted nonlinear field with one solved potential per current electrode.
#[derive(Clone)]
pub struct SpatialNonlinearCurrentInduction<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    field: SpatialNonlinearInduction<'m, C>,
    terminals: Arc<TerminalLayout>,
    terminal_unknown: Arc<[Option<usize>]>,
}

pub struct SpatialNonlinearTerminalInductionCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    state: SpatialNonlinearTerminalInduction<'m, C>,
    step: SpatialNonlinearInductionStep,
}
pub struct SpatialNonlinearCurrentInductionCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    state: SpatialNonlinearCurrentInduction<'m, C>,
    step: SpatialNonlinearInductionStep,
}

impl<'m, C: SpatialConductorMaterial> SpatialNonlinearTerminalInductionCandidate<'m, C> {
    pub fn state(&self) -> &SpatialNonlinearTerminalInduction<'m, C> {
        &self.state
    }
    pub fn step(&self) -> &SpatialNonlinearInductionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialNonlinearTerminalInduction<'m, C>,
        SpatialNonlinearInductionStep,
    ) {
        (self.state, self.step)
    }
}
impl<'m, C: SpatialConductorMaterial> SpatialNonlinearCurrentInductionCandidate<'m, C> {
    pub fn state(&self) -> &SpatialNonlinearCurrentInduction<'m, C> {
        &self.state
    }
    pub fn step(&self) -> &SpatialNonlinearInductionStep {
        &self.step
    }
    pub fn into_parts(
        self,
    ) -> (
        SpatialNonlinearCurrentInduction<'m, C>,
        SpatialNonlinearInductionStep,
    ) {
        (self.state, self.step)
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialNonlinearTerminalInduction<'m, C> {
    /// Attach disjoint exterior electrode faces to an accepted field. Every
    /// electrode on a conductor shares its connected fixed-A reference trace.
    /// No accepted A, impressed history, source or physical clock is changed.
    pub fn new(
        field: SpatialNonlinearInduction<'m, C>,
        electrode_faces: &[Vec<usize>],
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        let (terminals, coupling) = TerminalLayout::prepare(
            geometry(field.problem()),
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

    pub fn problem(&self) -> &SpatialNonlinearVectorMagnetostatics {
        self.field.problem()
    }
    pub fn potential_wb(&self) -> &[f64] {
        self.field.potential_wb()
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

    /// Regroup charge nodes into equipotential total-current equations once.
    /// The first electrode per conductor fixes only the scalar reference.
    pub fn into_current_driven(
        self,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearCurrentInduction<'m, C>> {
        let (coupling, terminal_unknown) = self.terminals.current_coupling(
            geometry(self.field.problem()),
            &self.field.prepared.coupling,
            magnetic_entries(&self.field)?,
            limits,
            cancelled,
        )?;
        Ok(SpatialNonlinearCurrentInduction {
            field: with_coupling(self.field, coupling, cancelled)?,
            terminals: self.terminals,
            terminal_unknown: terminal_unknown.into(),
        })
    }

    /// Backward Euler uses ENDPOINT electrode voltages and magnetic response.
    /// Reaction currents and incoming work are recovered from the actual field.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_candidate(
        &self,
        endpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalInductionCandidate<'m, C>> {
        self.step_prepared(
            self.field.problem.clone(),
            Arc::clone(&self.field.prepared),
            endpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Refresh original conductivity while preserving electrode charge maps.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_with_conductors_candidate(
        &self,
        endpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalInductionCandidate<'m, C>> {
        let prepared = self.field.prepare_conductors(conductors, cancelled)?;
        self.step_prepared(
            self.field.problem.clone(),
            prepared,
            endpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Source-bound magnetic/transport refresh retains the original field and
    /// terminal support. A thermal consumer must close magnetic parameter work.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_with_materials_candidate(
        &self,
        endpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Vec<SpatialNonlinearMaterial>,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalInductionCandidate<'m, C>> {
        let (problem, prepared) = self
            .field
            .prepare_materials(materials, conductors, cancelled)?;
        self.step_prepared(
            problem,
            prepared,
            endpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &self,
        problem: MagneticOwner<'m>,
        prepared: Arc<Prepared<C>>,
        voltage_v: &[f64],
        current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalInductionCandidate<'m, C>> {
        let drive = self.terminals.voltage_drive(voltage_v, cancelled)?;
        let (field, step) = self
            .field
            .step_driven_candidate(
                problem,
                prepared,
                current_a_m2,
                velocity_m_s,
                duration_s,
                controls,
                Some(&drive),
                cancelled,
            )?
            .into_parts();
        Ok(SpatialNonlinearTerminalInductionCandidate {
            state: Self {
                field,
                terminals: Arc::clone(&self.terminals),
            },
            step,
        })
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialNonlinearCurrentInduction<'m, C> {
    pub fn problem(&self) -> &SpatialNonlinearVectorMagnetostatics {
        self.field.problem()
    }
    pub fn potential_wb(&self) -> &[f64] {
        self.field.potential_wb()
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

    /// Restore prescribed voltage traces without resetting magnetic storage.
    pub fn into_voltage_driven(
        self,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearTerminalInduction<'m, C>> {
        let coupling = self.terminals.voltage_coupling(
            geometry(self.field.problem()),
            &self.field.prepared.coupling,
            magnetic_entries(&self.field)?,
            limits,
            cancelled,
        )?;
        Ok(SpatialNonlinearTerminalInduction {
            field: with_coupling(self.field, coupling, cancelled)?,
            terminals: self.terminals,
        })
    }

    /// Prescribe total ENDPOINT current, positive into each electrode. Only
    /// totals are constrained; the coupled field solves the spatial distribution
    /// and electrode voltages, with balance checked on every physical component.
    #[allow(clippy::too_many_arguments)]
    pub fn step_current_candidate(
        &self,
        endpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearCurrentInductionCandidate<'m, C>> {
        self.step_prepared(
            self.field.problem.clone(),
            Arc::clone(&self.field.prepared),
            endpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
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
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearCurrentInductionCandidate<'m, C>> {
        let prepared = self.field.prepare_conductors(conductors, cancelled)?;
        self.step_prepared(
            self.field.problem.clone(),
            prepared,
            endpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    pub fn step_current_with_materials_candidate(
        &self,
        endpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Vec<SpatialNonlinearMaterial>,
        conductors: Vec<Option<C>>,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearCurrentInductionCandidate<'m, C>> {
        let (problem, prepared) = self
            .field
            .prepare_materials(materials, conductors, cancelled)?;
        self.step_prepared(
            problem,
            prepared,
            endpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_prepared(
        &self,
        problem: MagneticOwner<'m>,
        prepared: Arc<Prepared<C>>,
        incoming_current_a: &[f64],
        current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialNonlinearInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialNonlinearCurrentInductionCandidate<'m, C>> {
        let drive = self.terminals.current_drive(
            incoming_current_a,
            &self.terminal_unknown,
            controls.induction.current_tolerance_a,
            cancelled,
        )?;
        let (field, step) = self
            .field
            .step_driven_candidate(
                problem,
                prepared,
                current_a_m2,
                velocity_m_s,
                duration_s,
                controls,
                Some(&drive),
                cancelled,
            )?
            .into_parts();
        Ok(SpatialNonlinearCurrentInductionCandidate {
            state: Self {
                field,
                terminals: Arc::clone(&self.terminals),
                terminal_unknown: Arc::clone(&self.terminal_unknown),
            },
            step,
        })
    }
}
