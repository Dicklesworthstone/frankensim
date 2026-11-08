//! Equipotential electrodes on the actual three-dimensional inducing conductors.
//!
//! The same Whitney A–phi solve determines skin/induction current, electrode
//! reaction current, Joule heat and Lorentz work. Electrodes are disjoint unions
//! of exterior faces with fixed tangential A; their tangential laboratory E is
//! therefore zero. Other conductor boundaries remain naturally insulating.
//! Electrodes on one conductor must share a connected prescribed-A trace: this
//! supplies a fixed-induction reference path and preserves the existing gauge.
//! This is the massive-conductor a-v terminal formulation, not a stranded
//! impressed-current surrogate (see GetDP's 3D circuit tutorial:
//! <https://getdp.info/dev/doc/texinfo/getdp.html>).
//!
//! Voltages are prescribed at the midpoint stage. Each face-connected conductor
//! uses its first terminal as reference, avoiding loss of precision and fictitious
//! power from large common-mode voltages. Physical current is positive INTO the
//! conductor: I_t = -sum_{nodes in t} integral grad(N_i).J dV. Consequently the
//! terminal supply is dt*sum_t (V_t-V_ref)*I_t. The full field update independently
//! checks charge conservation, Ampere's law and magnetic/Ohmic/mechanical balance.
//!
//! Geometry and occupancy stay fixed. The explicit source-pinned refresh path
//! permits temperature-dependent induced-linear permeability and conductivity.
//! This does not add electrode contact impedance, external lead fields,
//! moving interfaces, displacement current or a nonlinear material time step.

use super::mechanical::{SpatialLorentzLoadError, SpatialLorentzLoadLimits, SpatialLorentzLoads};
use super::*;
use fs_feec::weighted::WeightedAssemblyLimits;
use std::collections::BTreeSet;

/// Physical electrode closure against the native RLC circuit and sources.
pub mod circuit;
/// Total-current control with solved equipotential electrode voltages.
pub mod current;

/// Native stage reaction of one physical equipotential electrode.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SpatialInductionTerminalTransfer {
    /// Index in the electrode list supplied at preparation.
    pub terminal: usize,
    /// Original face-connected conductor component identifier.
    pub component: usize,
    /// First electrode of this conductor, used only as voltage reference.
    pub reference_terminal: usize,
    /// Supplied stage voltage relative to the component reference [V].
    pub voltage_to_reference_v: f64,
    /// Positive current enters this electrode's original conductor [A].
    pub incoming_current_a: f64,
    /// Relative voltage times incoming current times duration [J].
    pub supplied_work_j: f64,
}

#[derive(Clone)]
pub(super) struct TerminalLayout {
    faces: Vec<Vec<usize>>,
    component: Vec<usize>,
    reference: Vec<usize>,
    charge_terminal: Vec<Option<usize>>,
}

impl TerminalLayout {
    /// Admit physical electrodes against the original geometry and conductor
    /// support, then build the native prescribed-potential charge space once.
    pub fn prepare<C: SpatialConductorMaterial>(
        g: Geometry<'_>,
        coupling: &Coupling<C>,
        electrode_faces: &[Vec<usize>],
        magnetic_entries: usize,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Arc<Self>, Coupling<C>)> {
        poll(cancelled)?;
        Self::admit_assembly(coupling, magnetic_entries, limits)?;
        if electrode_faces.is_empty()
            || electrode_faces.iter().any(Vec::is_empty)
            || electrode_faces
                .iter()
                .try_fold(0usize, |n, f| n.checked_add(f.len()))
                .is_none_or(|n| n > g.geometry.faces.len())
        {
            return Err(invalid(
                "nonempty bounded physical induction electrodes required",
            ));
        }
        let mut faces = Vec::with_capacity(electrode_faces.len());
        let mut component = Vec::with_capacity(electrode_faces.len());
        let mut charge_terminal = vec![None; coupling.charge_vertices.len()];
        let mut used_faces = BTreeSet::new();
        let trace_components =
            super::coupling::fixed_trace_components(g.mesh(), &coupling.fixed_edges, cancelled)?;
        let topology = super::super::VectorTopology {
            geometry: g.geometry,
            cells: g.cells,
        };
        let mut conductor_trace = BTreeMap::new();
        for (terminal, supplied) in electrode_faces.iter().enumerate() {
            let mut selected = supplied.clone();
            selected.sort_unstable();
            let mut region = None;
            for &face in &selected {
                poll(cancelled)?;
                let info = g
                    .geometry
                    .faces
                    .get(face)
                    .ok_or_else(|| invalid("induction electrode face is outside the mesh"))?;
                if info.cells.len() != 1
                    || !used_faces.insert(face)
                    || !topology
                        .face_edges(face)
                        .is_some_and(|edges| edges.iter().all(|&e| coupling.fixed_edges[e]))
                {
                    return Err(invalid(
                        "distinct exterior electrode faces with fixed tangential A required",
                    ));
                }
                let cell = info.cells[0];
                let index = coupling
                    .cells
                    .binary_search_by_key(&cell, |c| c.cell)
                    .map_err(|_| invalid("induction electrode must touch a conducting cell"))?;
                let conductor = coupling.cells[index].component;
                if region.is_some_and(|r| r != conductor) {
                    return Err(invalid(
                        "one electrode cannot join disconnected conductor components",
                    ));
                }
                region = Some(conductor);
                for &v in &g.mesh().faces[face] {
                    let trace = trace_components[v as usize];
                    let reference = *conductor_trace.entry(conductor).or_insert(trace);
                    if trace != reference {
                        return Err(invalid(
                            "conductor electrodes need a connected fixed-A voltage reference path",
                        ));
                    }
                    let charge = coupling
                        .charge_vertices
                        .binary_search(&(conductor, v as usize))
                        .map_err(|_| {
                            invalid("electrode lost its original conductor charge node")
                        })?;
                    if charge_terminal[charge].is_some_and(|t| t != terminal) {
                        return Err(invalid(
                            "distinct electrodes cannot meet at a conductor vertex",
                        ));
                    }
                    charge_terminal[charge] = Some(terminal);
                }
            }
            faces.push(selected);
            component.push(region.expect("nonempty electrode admitted"));
        }
        let mut references = BTreeMap::new();
        let reference = component
            .iter()
            .enumerate()
            .map(|(t, &c)| *references.entry(c).or_insert(t))
            .collect();
        let layout = Arc::new(Self {
            faces,
            component,
            reference,
            charge_terminal,
        });
        let coupling = layout.voltage_coupling(g, coupling, magnetic_entries, limits, cancelled)?;
        poll(cancelled)?;
        Ok((layout, coupling))
    }

    fn admit_assembly<C: SpatialConductorMaterial>(
        coupling: &Coupling<C>,
        magnetic_entries: usize,
        limits: SpatialInductionLimits,
    ) -> Result<()> {
        if coupling
            .cells
            .len()
            .checked_mul(100)
            .and_then(|n| n.checked_add(magnetic_entries))
            .is_none_or(|n| n > limits.maximum_assembly_terms)
        {
            return Err(invalid("physical electrode assembly budget exceeded"));
        }
        Ok(())
    }

    pub fn voltage_coupling<C: SpatialConductorMaterial>(
        &self,
        g: Geometry<'_>,
        coupling: &Coupling<C>,
        magnetic_entries: usize,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Coupling<C>> {
        poll(cancelled)?;
        Self::admit_assembly(coupling, magnetic_entries, limits)?;
        let prescribed = self
            .charge_terminal
            .iter()
            .map(Option::is_some)
            .collect::<Vec<_>>();
        coupling.with_prescribed_charges(g, &prescribed, limits, cancelled)
    }

    pub fn current_coupling<C: SpatialConductorMaterial>(
        &self,
        g: Geometry<'_>,
        coupling: &Coupling<C>,
        magnetic_entries: usize,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Coupling<C>, Vec<Option<usize>>)> {
        poll(cancelled)?;
        Self::admit_assembly(coupling, magnetic_entries, limits)?;
        coupling.with_current_charges(g, &self.charge_terminal, &self.reference, limits, cancelled)
    }

    pub fn faces(&self) -> &[Vec<usize>] {
        &self.faces
    }

    pub fn references(&self) -> &[usize] {
        &self.reference
    }

    pub fn voltage_drive(
        &self,
        voltage_v: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TerminalDrive<'_>> {
        poll(cancelled)?;
        if voltage_v.len() != self.faces.len() || voltage_v.iter().any(|v| !v.is_finite()) {
            return Err(invalid(
                "one finite stage voltage per original induction electrode required",
            ));
        }
        let relative_voltage_v = voltage_v
            .iter()
            .zip(&self.reference)
            .map(|(&v, &reference)| finite(v - voltage_v[reference]))
            .collect::<Result<Vec<_>>>()?;
        let potential_v = self
            .charge_terminal
            .iter()
            .map(|t| t.map_or(0.0, |t| relative_voltage_v[t]))
            .collect();
        Ok(TerminalDrive {
            layout: self,
            potential_v,
            relative_voltage_v,
            current: None,
        })
    }

    pub fn current_drive<'a>(
        &'a self,
        incoming_a: &'a [f64],
        unknowns: &'a [Option<usize>],
        current_tolerance_a: f64,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<TerminalDrive<'a>> {
        poll(cancelled)?;
        if incoming_a.len() != self.faces.len()
            || unknowns.len() != self.faces.len()
            || incoming_a.iter().any(|i| !i.is_finite())
            || !positive(current_tolerance_a)
        {
            return Err(invalid(
                "one finite total incoming current per electrode required",
            ));
        }
        let mut balance = BTreeMap::new();
        for (&component, &incoming) in self.component.iter().zip(incoming_a) {
            poll(cancelled)?;
            let sum = balance.entry(component).or_insert(0.0);
            *sum = finite(*sum + incoming)?;
        }
        if balance.into_values().any(|i| i.abs() > current_tolerance_a) {
            return Err(invalid(
                "total electrode current must balance on every conductor",
            ));
        }
        Ok(TerminalDrive {
            layout: self,
            potential_v: vec![],
            relative_voltage_v: vec![],
            current: Some(CurrentDrive {
                incoming_a,
                unknowns,
            }),
        })
    }
}

/// An accepted magnetic history with fixed physical electrode locations.
/// The native history remains private so a grounded-electrode charge space
/// cannot accidentally be advanced as an electrically isolated conductor.
#[derive(Clone)]
pub struct SpatialTerminalInduction<'m, C: SpatialConductorMaterial = ResolvedConductorState> {
    field: SpatialInduction<'m, C>,
    terminals: Arc<TerminalLayout>,
}

/// Atomic field candidate and all original physical transfers.
pub struct SpatialTerminalInductionCandidate<
    'm,
    C: SpatialConductorMaterial = ResolvedConductorState,
> {
    state: SpatialTerminalInduction<'m, C>,
    step: SpatialInductionStep,
}

impl<'m, C: SpatialConductorMaterial> SpatialTerminalInductionCandidate<'m, C> {
    pub fn state(&self) -> &SpatialTerminalInduction<'m, C> {
        &self.state
    }
    pub fn step(&self) -> &SpatialInductionStep {
        &self.step
    }
    pub fn into_parts(self) -> (SpatialTerminalInduction<'m, C>, SpatialInductionStep) {
        (self.state, self.step)
    }

    /// Reuse the original consistent mechanical transfer of this SAME current.
    pub fn lorentz_loads(
        &self,
        selected_cells: &[usize],
        limits: SpatialLorentzLoadLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> core::result::Result<SpatialLorentzLoads, SpatialLorentzLoadError> {
        let problem = self.state.problem();
        super::mechanical::recover(
            problem.positions_m(),
            &problem.mesh().tets,
            &self.step.cells,
            selected_cells,
            limits,
            cancelled,
        )
    }
}

impl<'m, C: SpatialConductorMaterial> SpatialTerminalInduction<'m, C> {
    /// Attach physical electrodes to an accepted native induction history.
    /// Each electrode must lie on a single conducting component. Electrodes
    /// cannot overlap even at a vertex of that component. Face lists are
    /// canonicalized, while electrode order retains its terminal identity.
    /// No magnetic cochain, physical field or clock is changed by preparation.
    pub fn new(
        field: SpatialInduction<'m, C>,
        electrode_faces: &[Vec<usize>],
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(cancelled)?;
        let p = &field.prepared;
        let problem = &*p.problem;
        let (terminals, coupling) = TerminalLayout::prepare(
            geometry(problem),
            &p.coupling,
            electrode_faces,
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
            field: SpatialInduction { prepared, ..field },
            terminals,
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
    /// Canonically ordered face sets in the original terminal order.
    pub fn electrode_faces(&self) -> &[Vec<usize>] {
        &self.terminals.faces
    }
    pub fn field(&self, cancelled: &mut impl FnMut() -> bool) -> Result<SpatialInductionField> {
        self.field.field(cancelled)
    }

    /// Prepare one scalar unknown per non-reference electrode once, retaining
    /// this exact accepted magnetic history for subsequent current-driven steps.
    pub fn into_current_driven(
        self,
        limits: SpatialInductionLimits,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<current::SpatialCurrentInduction<'m, C>> {
        current::SpatialCurrentInduction::from_voltage_driven(self, limits, cancelled)
    }

    /// Resolve peak AC electrode voltages in the original Whitney A–phi space.
    /// The conductor current distribution, electrode reactions, skin/proximity
    /// response and prescribed-motion work are solved together. Voltages use
    /// `exp(-i omega t)` and are referenced separately on each conductor before
    /// assembly, preserving complex common-mode invariance.
    ///
    /// This is a fixed-material, fixed-geometry settled AC response. It does
    /// not advance the accepted transient history, model a switch-on transient
    /// or infer contact impedance or lead fields. Its source-bound heat and
    /// Lorentz-load APIs apply to this same accepted harmonic field.
    #[allow(clippy::too_many_arguments)]
    pub fn solve_harmonic_voltage(
        &self,
        voltage_v: &[fs_math::c64::C64],
        impressed_current_a_m2: &[[fs_math::c64::C64; 3]],
        velocity_m_s: &[[f64; 3]],
        omega_rad_s: f64,
        controls: super::harmonic::SpatialHarmonicInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<super::harmonic::SpatialHarmonicInductionResponse<'_, 'm, C>> {
        poll(cancelled)?;
        if voltage_v.len() != self.terminals.faces.len() {
            return Err(invalid(
                "one peak voltage per original induction electrode required",
            ));
        }
        let real: Vec<_> = voltage_v.iter().map(|v| v.re).collect();
        let imaginary: Vec<_> = voltage_v.iter().map(|v| v.im).collect();
        let real = self.terminals.voltage_drive(&real, cancelled)?;
        let imaginary = self.terminals.voltage_drive(&imaginary, cancelled)?;
        self.field.solve_harmonic_driven(
            impressed_current_a_m2,
            velocity_m_s,
            omega_rad_s,
            controls,
            None,
            Some([&real, &imaginary]),
            cancelled,
        )
    }

    /// Advance the original midpoint A–phi equations with prescribed stage
    /// electrode voltages. The actual current distribution is solved, never
    /// supplied as a uniform current density. Source and velocity conventions
    /// are identical to SpatialInduction::step_candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_candidate(
        &self,
        midpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalInductionCandidate<'m, C>> {
        self.step_voltage_prepared_candidate(
            Arc::clone(&self.field.prepared),
            midpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    /// Advance the accepted magnetic history with source-pinned endpoint
    /// conductivity and induced-linear permeability. The electrode layout and
    /// original before-state remain unchanged; parameter work is reported
    /// separately from physical Joule heat and electrode supply.
    #[allow(clippy::too_many_arguments)]
    pub fn step_voltage_with_materials_candidate(
        &self,
        midpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        materials: Vec<SpatialMagneticMaterial>,
        conductors: Vec<Option<C>>,
        assembly: WeightedAssemblyLimits,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalInductionCandidate<'m, C>> {
        let prepared = self
            .field
            .prepare_materials(materials, conductors, assembly, cancelled)?;
        self.step_voltage_prepared_candidate(
            prepared,
            midpoint_voltage_v,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            cancelled,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn step_voltage_prepared_candidate(
        &self,
        prepared: Arc<Prepared<'m, C>>,
        midpoint_voltage_v: &[f64],
        endpoint_impressed_current_a_m2: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        duration_s: f64,
        controls: SpatialInductionControls,
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<SpatialTerminalInductionCandidate<'m, C>> {
        let drive = self
            .terminals
            .voltage_drive(midpoint_voltage_v, cancelled)?;
        let native = self.field.step_driven_candidate(
            prepared,
            endpoint_impressed_current_a_m2,
            velocity_m_s,
            duration_s,
            controls,
            Some(&drive),
            cancelled,
        )?;
        let (field, step) = native.into_parts();
        Ok(SpatialTerminalInductionCandidate {
            state: Self {
                field,
                terminals: Arc::clone(&self.terminals),
            },
            step,
        })
    }
}

pub(super) struct TerminalDrive<'a> {
    layout: &'a TerminalLayout,
    pub potential_v: Vec<f64>,
    relative_voltage_v: Vec<f64>,
    current: Option<CurrentDrive<'a>>,
}
struct CurrentDrive<'a> {
    incoming_a: &'a [f64],
    unknowns: &'a [Option<usize>],
}
impl TerminalDrive<'_> {
    pub fn prescribed_potentials(&self) -> Option<&[f64]> {
        self.current
            .is_none()
            .then_some(self.potential_v.as_slice())
    }

    /// Add only physical electrode forcing to the original coupled charge and
    /// edge equations. Both prescribed-potential and grouped-current spaces use
    /// this same RHS path; the conductivity law remains owned by Coupling.
    pub fn add_rhs<C: SpatialConductorMaterial>(
        &self,
        coupling: &Coupling<C>,
        rhs: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        if rhs.len() != coupling.mass.nrows() {
            return Err(invalid(
                "terminal forcing requires the original coupled space",
            ));
        }
        if let Some(potential) = self.prescribed_potentials() {
            for (rhs, known) in rhs
                .iter_mut()
                .zip(coupling.potential_drive(potential, cancelled)?)
            {
                *rhs = finite(*rhs + known)?;
            }
        }
        self.add_current_rhs(rhs, cancelled)
    }

    pub fn add_current_rhs(
        &self,
        rhs: &mut [f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<()> {
        if let Some(current) = &self.current {
            for (&incoming, &unknown) in current.incoming_a.iter().zip(current.unknowns) {
                poll(cancelled)?;
                if let Some(i) = unknown {
                    rhs[i] = finite(rhs[i] + incoming)?;
                }
            }
        }
        Ok(())
    }

    pub fn transfers(
        &self,
        outward_charge_current_a: &[f64],
        duration_s: f64,
        solution: &[f64],
        cancelled: &mut impl FnMut() -> bool,
    ) -> Result<(Vec<SpatialInductionTerminalTransfer>, f64, f64)> {
        if outward_charge_current_a.len() != self.layout.charge_terminal.len() {
            return Err(invalid("terminal reaction lost its original charge map"));
        }
        let mut incoming = vec![0.0; self.layout.faces.len()];
        let mut residual = 0.0_f64;
        for (&current, &terminal) in outward_charge_current_a
            .iter()
            .zip(&self.layout.charge_terminal)
        {
            poll(cancelled)?;
            if let Some(t) = terminal {
                incoming[t] = finite(incoming[t] - current)?;
            } else {
                residual = residual.max(finite(current)?.abs());
            }
        }
        let mut component_balance = BTreeMap::new();
        let mut work = 0.0;
        let mut transfers = Vec::with_capacity(incoming.len());
        for (terminal, incoming_current_a) in incoming.into_iter().enumerate() {
            poll(cancelled)?;
            let component = self.layout.component[terminal];
            let balance = component_balance.entry(component).or_insert(0.0);
            *balance = finite(*balance + incoming_current_a)?;
            let voltage_to_reference_v = if let Some(current) = &self.current {
                residual =
                    residual.max(finite(incoming_current_a - current.incoming_a[terminal])?.abs());
                current.unknowns[terminal].map_or(Ok(0.0), |i| finite(solution[i] / duration_s))?
            } else {
                self.relative_voltage_v[terminal]
            };
            let supplied_work_j = finite(duration_s * voltage_to_reference_v * incoming_current_a)?;
            work = finite(work + supplied_work_j)?;
            transfers.push(SpatialInductionTerminalTransfer {
                terminal,
                component,
                reference_terminal: self.layout.reference[terminal],
                voltage_to_reference_v,
                incoming_current_a,
                supplied_work_j,
            });
        }
        for balance in component_balance.into_values() {
            residual = residual.max(balance.abs());
        }
        Ok((transfers, work, residual))
    }
}
