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
        if p.cells
            .len()
            .checked_mul(100)
            .and_then(|n| n.checked_add(problem.stiffness.nnz()))
            .is_none_or(|n| n > limits.maximum_assembly_terms)
        {
            return Err(invalid("current electrode assembly budget exceeded"));
        }
        let (coupling, terminal_unknown) = p.coupling.with_current_charges(
            geometry(problem),
            &state.terminals.charge_terminal,
            &state.terminals.reference,
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
        poll(cancelled)?;
        if midpoint_incoming_current_a.len() != self.terminals.faces.len()
            || midpoint_incoming_current_a.iter().any(|i| !i.is_finite())
            || !positive(controls.current_tolerance_a)
        {
            return Err(invalid(
                "one finite total incoming current per electrode required",
            ));
        }
        let mut balance = BTreeMap::new();
        for (&component, &incoming) in self
            .terminals
            .component
            .iter()
            .zip(midpoint_incoming_current_a)
        {
            poll(cancelled)?;
            let sum = balance.entry(component).or_insert(0.0);
            *sum = finite(*sum + incoming)?;
        }
        if balance
            .into_values()
            .any(|i| i.abs() > controls.current_tolerance_a)
        {
            return Err(invalid(
                "total electrode current must balance on every conductor",
            ));
        }
        let drive = TerminalDrive {
            layout: &self.terminals,
            potential_v: vec![],
            relative_voltage_v: vec![],
            current: Some(CurrentDrive {
                incoming_a: midpoint_incoming_current_a,
                unknowns: &self.terminal_unknown,
            }),
        };
        let candidate = self.field.step_driven_candidate(
            Arc::clone(&self.field.prepared),
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
