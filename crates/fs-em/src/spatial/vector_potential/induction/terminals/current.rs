//! Total-current electrodes in the original coupled Whitney A–phi equations.
//!
//! All charge nodes on one electrode share a single scalar unknown. Its weak
//! row is the SUM of electrode reactions, so only total incoming current is
//! prescribed; the spatial current profile and electrode voltage are solved.
//! The first electrode on each conductor fixes the scalar reference. Summed
//! incoming current must vanish separately on each face-connected conductor.
//!
//! The electrode space and sparse mass operator are prepared once. Each step
//! uses one native coupled solve, with the original magnetic storage, motion,
//! Ohmic heat, Lorentz load and terminal-work checks. No dense electrode
//! response sampling or imposed uniform current-density approximation is used.

use super::*;

/// Accepted native magnetic history with equipotential total-current terminals.
/// Keeping this owner distinct from voltage control prevents a prescribed
/// voltage from accidentally being applied to the grouped current unknowns.
#[derive(Clone)]
pub struct SpatialCurrentInduction<'m, C: SpatialConductorMaterial = ResolvedConductorState> {
    field: SpatialInduction<'m, C>,
    terminals: Arc<TerminalLayout>,
    terminal_unknown: Arc<[Option<usize>]>,
}

/// Unpublished field update, including solved electrode voltages and work.
pub struct SpatialCurrentInductionCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    state: SpatialCurrentInduction<'m, C>,
    step: SpatialInductionStep,
}

impl<'m, C: SpatialConductorMaterial> SpatialCurrentInduction<'m, C> {
    pub(super) fn from_voltage_driven(
        state: SpatialTerminalInduction<'m, C>,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let p = &state.field.prepared;
        let problem = &*p.problem;
        let (coupling, terminal_unknown) = state.terminals.current_coupling(
            geometry(problem),
            &p.coupling,
            problem.stiffness.nnz(),
            limits,
            cancelled,
        )?;
        let magnetic = coupling.restrict(&problem.stiffness, cancelled)?;
        let prepared = Arc::new(Prepared {
            problem: p.problem.clone(),
            coupling,
            magnetic,
        });
        poll(cancelled)?;
        Ok(Self {
            field: SpatialInduction {
                prepared,
                ..state.field
            },
            terminals: state.terminals,
            terminal_unknown: terminal_unknown.into(),
        })
    }

    pub fn problem(&self) -> &SpatialVectorMagnetostatics {
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
        &self.terminals.faces
    }
    pub fn field(&self, cancelled: &mut impl FnMut() -> bool) -> Result<SpatialInductionField> {
        self.field.field(cancelled)
    }

    /// Re-prepare voltage control without changing the accepted A or clock.
    pub fn into_voltage_driven(
        self,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalInduction<'m, C>> {
        SpatialTerminalInduction::new(self.field, &self.terminals.faces, limits, cancelled)
    }

    /// Resolve peak total electrode currents and their physical AC voltages.
    /// Each electrode retains one grouped scalar-potential unknown; only its
    /// total current is imposed, not its spatial distribution. Incoming peak
    /// phasors use `exp(-i omega t)` and must balance on each conductor.
    ///
    /// This reuses one native sparse A–phi solve and the same real/reactive
    /// power, charge, Ampere, skin-depth and motion checks as voltage drive.
    /// Material, geometry and steady Eulerian velocity are held fixed. Neither
    /// this history's magnetic cochain nor its physical clock is advanced.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic_current(
        &self,
        incoming_current_a: &[fs_math::c64::C64],
        impressed_current_a_m2: &[[fs_math::c64::C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: super::super::harmonic::SpatialHarmonicInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<super::super::harmonic::SpatialHarmonicInductionResponse<'_, 'm, C>> {
        poll(cancelled)?;
        if incoming_current_a.len() != self.terminals.faces.len() {
            return Err(invalid(
                "one peak current per original induction electrode required",
            ));
        }
        let real: Vec<_> = incoming_current_a.iter().map(|v| v.re).collect();
        let imaginary: Vec<_> = incoming_current_a.iter().map(|v| v.im).collect();
        let real_drive = self.terminals.current_drive(
            &real,
            &self.terminal_unknown,
            controls.current_tolerance_a,
            cancelled,
        )?;
        let imaginary_drive = self.terminals.current_drive(
            &imaginary,
            &self.terminal_unknown,
            controls.current_tolerance_a,
            cancelled,
        )?;
        self.field.solve_harmonic_driven(
            impressed_current_a_m2,
            velocity_m_s,
            omega_rad_s,
            controls,
            None,
            Some([&real_drive, &imaginary_drive]),
            cancelled,
        )
    }

    /// Prescribe stage total currents positive INTO each original electrode.
    /// Electrode potential is an output; interior and tangential current remain
    /// part of the same resolved induction/charge system. All input currents,
    /// including the reference electrode's, are independently checked on return.
    #[allow(clippy::too_many_arguments)]
    pub fn step_current_candidate(
        &self,
        midpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialCurrentInductionCandidate<'m, C>> {
        self.step_current_prepared_candidate(
            Arc::clone(&self.field.prepared),
            midpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Prescribe total electrode currents while refreshing original
    /// source-backed conductivity and induced-linear permeability. The single
    /// sparse solve retains the original A history and grouped electrode rows;
    /// signed magnetic parameter work is distinct from Ohmic heat.
    #[allow(clippy::too_many_arguments)]
    pub fn step_current_with_materials_candidate(
        &self,
        midpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<C>>,
        assembly: WeightedAssemblyLimits,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialCurrentInductionCandidate<'m, C>> {
        let prepared = self
            .field
            .prepare_materials(materials, conductors, assembly, cancelled)?;
        self.step_current_prepared_candidate(
            prepared,
            midpoint_incoming_current_a,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_current_prepared_candidate(
        &self,
        prepared: Arc<Prepared<'m, C>>,
        midpoint_incoming_current_a: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialCurrentInductionCandidate<'m, C>> {
        let drive = self.terminals.current_drive(
            midpoint_incoming_current_a,
            &self.terminal_unknown,
            controls.current_tolerance_a,
            cancelled,
        )?;
        let candidate = self.field.step_driven_candidate(
            prepared,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            Some(&drive),
            cancelled,
        )?;
        let (field, step) = candidate.into_parts();
        Ok(SpatialCurrentInductionCandidate {
            state: Self {
                field,
                terminals: Arc::clone(&self.terminals),
                terminal_unknown: Arc::clone(&self.terminal_unknown),
            },
            step,
        })
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialCurrentInductionCandidate<'m, C> {
    pub fn state(&self) -> &SpatialCurrentInduction<'m, C> {
        &self.state
    }
    pub fn step(&self) -> &SpatialInductionStep {
        &self.step
    }
    pub fn into_parts(self) -> (SpatialCurrentInduction<'m, C>, SpatialInductionStep) {
        (self.state, self.step)
    }
    /// Mechanical work uses the same original current recovered by this step.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<SpatialLorentzLoads, SpatialLorentzLoadError> {
        let problem = self.state.problem();
        super::super::mechanical::recover(
            problem.positions_m(),
            &problem.mesh().tets,
            &self.step.cells,
            selected_cells,
            limits,
            cancelled,
        )
    }
}
