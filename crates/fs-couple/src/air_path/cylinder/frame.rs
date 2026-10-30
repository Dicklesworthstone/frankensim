//! Source-defined elastic frames in the original causal cylinder fluid.
//!
//! Dry frame modes/normalization come from fs-solid. The same source section's
//! outer radius and the native cubic transverse trace define its fluid load.
//! Four positive Gauss points exactly integrate products of cubic traces for
//! each constant-section member. The reduced inertia is FULL:
//! M = I + sum(w*m_added*b*b^T), not a list of retuned modal frequencies.
//! Each fluid branch retains its own momentum, with resistance
//! w*c*(b^T*v-v_fluid)^2. Original material viscosity keeps all cross terms.
//!
//! This is frozen, small, planar vibration at fixed supports in a stationary
//! Newtonian exterior. Circular prismatic members, local 2-D Stokes flow and
//! explicit modal/spectral truncation apply. No axial/rotary fluid drag,
//! inter-member hydrodynamic interactions, thermal evolution, moving supports,
//! cavitation, nonlinear wake or acoustic radiation is inferred. Hollow members
//! receive EXTERIOR loading only; internal fluid needs its own owner.

use super::{CylinderFlow, CylinderFlowResolution};
use crate::acoustic_realize::AcousticRealizeError;
use crate::frame_vibration::{FrameVibrationBasis, FrameVibrationBudget, FrameVibrationError};
use crate::modal_acoustic_time::ModalAcousticState;
use fs_exec::CancelGate;
use fs_la::factor::{FactorError, cholesky};
use fs_material::fluid_transport::FluidTransport;
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, StepRecord};
use fs_solid::frame2d::Frame2d;

/// Original member index and original exterior material source; never a radius,
/// mass or stiffness override. Different members may have different fluids.
pub struct FrameFluidBinding<F> {
    pub element: usize,
    pub fluid: F,
}

/// Complete preparation limits, in addition to the native dry-frame budget.
#[derive(Clone, Copy, Debug)]
pub struct FrameFluidBudget {
    /// 2*retained_modes + 4*loaded_members*spectral_intervals.
    pub max_states: usize,
    /// Maximum entries in EACH full dense state matrix; not total byte usage.
    /// Native factorizations and time steppers also require workspace.
    pub max_dense_entries: usize,
    /// Six scalar projection terms per retained mode per fluid station or
    /// original material resistance port. Checked before projected allocation.
    pub max_projection_terms: usize,
    /// Infinity norm allowance for M*M_inverse-I, finite and in (0,1).
    /// Arithmetic admission only, not a continuum-error certificate.
    pub mass_inverse_tolerance: f64,
}

#[derive(Debug)]
pub enum FrameFluidError {
    Invalid(&'static str),
    Frame(FrameVibrationError),
    Fluid(AcousticRealizeError),
    Factor(FactorError),
    Native(PhsError),
    Balance { residual_j: f64, tolerance_j: f64 },
    Cancelled,
}
impl core::fmt::Display for FrameFluidError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "material frame/fluid: {self:?}")
    }
}
impl std::error::Error for FrameFluidError {}
impl From<FrameVibrationError> for FrameFluidError {
    fn from(error: FrameVibrationError) -> Self {
        match error { FrameVibrationError::Cancelled => Self::Cancelled, other => Self::Frame(other) }
    }
}
impl From<fs_solid::SolidError> for FrameFluidError {
    fn from(error: fs_solid::SolidError) -> Self { Self::Frame(error.into()) }
}
impl From<PhsError> for FrameFluidError {
    fn from(error: PhsError) -> Self { Self::Native(error) }
}
type Result<T> = core::result::Result<T, FrameFluidError>;
fn invalid(message: &'static str) -> FrameFluidError { FrameFluidError::Invalid(message) }
fn poll(gate: &CancelGate) -> Result<()> {
    if gate.is_requested() { Err(FrameFluidError::Cancelled) } else { Ok(()) }
}
fn finite(x: f64) -> Result<f64> {
    if x.is_finite() { Ok(x) } else { Err(invalid("unrepresentable frame/fluid coefficient or state")) }
}
fn positive(x: f64) -> Result<f64> {
    finite(x)?;
    if x > 0.0 { Ok(x) } else { Err(invalid("positive frame/fluid coefficient required")) }
}
fn dot(a: &[f64], b: &[f64]) -> Result<f64> {
    if a.len() != b.len() { return Err(invalid("frame/fluid projection dimensions differ")); }
    finite(a.iter().zip(b).map(|(a,b)| a*b).sum())
}

/// One unchanged structural member and its source-preserving exterior operator.
pub struct LoadedFrameMember<F> {
    element: usize,
    flow: CylinderFlow<F>,
}
impl<F> LoadedFrameMember<F> {
    pub const fn element(&self) -> usize { self.element }
    pub const fn flow(&self) -> &CylinderFlow<F> { &self.flow }
}

struct Station {
    member: usize,
    weight_m: f64,
    projection: Vec<f64>,
    first_memory: usize,
}
struct MaterialPort {
    projection: Vec<f64>,
    resistance: f64,
}

/// Frozen complete physical component, retaining the original frame and fluids.
/// State order: [dry modal q, LOADED modal momentum p=M*qdot, fluid momenta].
/// Fluid coordinates are member-major, then Gauss-point-major, then the original
/// spectral-branch order. They are never reset when an interval is accepted.
pub struct MaterialFluidFrame<'a, F> {
    frame: &'a Frame2d,
    basis: FrameVibrationBasis,
    members: Vec<LoadedFrameMember<F>>,
    stations: Vec<Station>,
    material_ports: Vec<MaterialPort>,
    modal_mass: Vec<f64>,
    memory_masses: Vec<f64>,
    inverse_residual: f64,
    system: PortHamiltonian,
}

/// Native proposal plus separate physical dissipation accounts. These losses
/// have NOT been deposited into any thermal owner. `native.x` is the complete
/// next state; publication remains the caller's joint-transaction decision.
pub struct FrameFluidStep {
    pub native: StepRecord,
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    /// Same order as basis().material_damping().ports(), or empty for elastic frames.
    pub material_port_dissipation_j: Vec<f64>,
    /// Same original order as members(), including low-tail and memory-slip loss.
    pub fluid_member_dissipation_j: Vec<f64>,
    /// Held nodal forces/moments dotted with actual reconstructed displacement change.
    pub nodal_work_j: f64,
    pub energy_residual_j: f64,
}

const GAUSS: [(f64, f64); 4] = [
    (0.06943184420297371, 0.17392742256872693),
    (0.33000947820757187, 0.32607257743127307),
    (0.6699905217924281, 0.32607257743127307),
    (0.9305681557970262, 0.17392742256872693),
];

impl<'a, F: FluidTransport> MaterialFluidFrame<'a, F> {
    /// Prepare one coupled model from the actual source frame and fluid data.
    ///
    /// Each bound member must have its original circular/concentric material
    /// section; radius is read from that section. Duplicate members refuse.
    /// Dry mass, stiffness, supports, source damping and modal extraction remain
    /// owned by the original frame. No artificial gas or homogeneous solid is
    /// created. All Gauss points and all spectral branches count toward the
    /// COMPLETE budget before their arrays are allocated; none are silently cut.
    ///
    /// Empty bindings are the exact dry Galerkin model, with its existing loss.
    /// Retained-mode and spectral-cutoff refinement remain caller obligations.
    /// Setup/stepping poll around native synchronous solves, not inside them.
    pub fn new(
        frame: &'a Frame2d,
        retained_modes: usize,
        bindings: Vec<FrameFluidBinding<F>>,
        resolution: CylinderFlowResolution,
        frame_budget: FrameVibrationBudget,
        budget: FrameFluidBudget,
        gate: &CancelGate,
    ) -> Result<Self> {
        poll(gate)?;
        resolution.validate().map_err(FrameFluidError::Fluid)?;
        if retained_modes == 0 || !budget.mass_inverse_tolerance.is_finite()
            || budget.mass_inverse_tolerance <= 0.0 || budget.mass_inverse_tolerance >= 1.0 {
            return Err(invalid("retained modes and a finite inverse residual allowance in (0,1) required"));
        }
        let memory_count = bindings.len().checked_mul(4)
            .and_then(|n| n.checked_mul(resolution.intervals))
            .ok_or_else(|| invalid("complete frame fluid-memory count overflow"))?;
        let n = retained_modes.checked_mul(2).and_then(|n| n.checked_add(memory_count))
            .ok_or_else(|| invalid("complete frame/fluid state count overflow"))?;
        let square = n.checked_mul(n).ok_or_else(|| invalid("frame/fluid matrix size overflow"))?;
        if n > budget.max_states || square > budget.max_dense_entries {
            return Err(invalid("complete frame/fluid state or dense matrix budget exceeded"));
        }
        for (i, binding) in bindings.iter().enumerate() {
            poll(gate)?;
            if bindings[..i].iter().any(|b| b.element == binding.element) {
                return Err(invalid("a physical member cannot receive the same exterior loading twice"));
            }
            frame.reference_transverse_trace(binding.element, 0.5)?;
            if frame.material_section(binding.element).is_none() {
                return Err(invalid("fluid loading needs the original member material section and circular geometry"));
            }
        }
        let basis = FrameVibrationBasis::from_frame(frame, retained_modes, frame_budget,
            || gate.is_requested())?;
        let m = retained_modes;
        let port_count = basis.material_damping().map_or(0, |d| d.ports().len());
        if bindings.len().checked_mul(4).and_then(|n| n.checked_add(port_count))
            .and_then(|n| n.checked_mul(m)).and_then(|n| n.checked_mul(6))
            .is_none_or(|n| n > budget.max_projection_terms) {
            return Err(invalid("frame/fluid physical projection budget exceeded"));
        }
        let mut members = Vec::with_capacity(bindings.len());
        let mut stations = Vec::with_capacity(bindings.len()*4);
        let mut memory_masses = Vec::with_capacity(memory_count);
        let mut modal_mass = vec![0.0; m*m];
        for i in 0..m { modal_mass[i*m+i] = 1.0; }
        for binding in bindings {
            poll(gate)?;
            let section = frame.material_section(binding.element)
                .ok_or_else(|| invalid("original circular material section is absent"))?;
            let flow = CylinderFlow::from_fluid(section.outer_radius_m(), binding.fluid, resolution, gate)
                .map_err(|error| match error {
                    AcousticRealizeError::Cancelled => FrameFluidError::Cancelled,
                    other => FrameFluidError::Fluid(other),
                })?;
            let member = members.len();
            for (xi, weight) in GAUSS {
                poll(gate)?;
                let (dofs, weights, length) = frame.reference_transverse_trace(binding.element, xi)?;
                let projection = project(&basis, &dofs, &weights)?;
                let weight_m = positive(weight*length)?;
                add_outer(&mut modal_mass, m, 0, &projection,
                    positive(weight_m*flow.added_mass_per_length_kg_m())?)?;
                let first_memory = memory_masses.len();
                for branch in flow.branches() {
                    let mass = positive(weight_m*branch.mass_per_length_kg_m)?;
                    positive(1.0/mass)?;
                    memory_masses.push(mass);
                }
                stations.push(Station { member, weight_m, projection, first_memory });
            }
            members.push(LoadedFrameMember { element: binding.element, flow });
        }
        let mut material_ports = Vec::new();
        if let Some(damping) = basis.material_damping() {
            for port in damping.ports() {
                poll(gate)?;
                let resistance = finite(port.damping_n_s_m)?;
                if resistance < 0.0 { return Err(invalid("negative original material resistance")); }
                material_ports.push(MaterialPort {
                    projection: project(&basis, &port.dofs, &port.displacement_weights)?, resistance,
                });
            }
        }
        poll(gate)?;
        let (inverse, inverse_residual) = inverse_mass(&modal_mass, m,
            budget.mass_inverse_tolerance, gate)?;
        let mut q = vec![0.0; square];
        let mut j = vec![0.0; square];
        let mut r = vec![0.0; square];
        let mut g = vec![0.0; n*m];
        for i in 0..m {
            q[i*n+i] = basis.modes()[i].eigenvalue_per_s2;
            j[i*n+m+i] = 1.0;
            j[(m+i)*n+i] = -1.0;
            g[(m+i)*m+i] = 1.0;
            for k in 0..m { q[(m+i)*n+m+k] = inverse[i*m+k]; }
        }
        for port in &material_ports {
            poll(gate)?;
            add_outer(&mut r, n, m, &port.projection, port.resistance)?;
        }
        for station in &stations {
            poll(gate)?;
            let flow = &members[station.member].flow;
            add_outer(&mut r, n, m, &station.projection,
                positive(station.weight_m*flow.low_tail_resistance_per_length_n_s_m2())?)?;
            for (k, branch) in flow.branches().iter().enumerate() {
                poll(gate)?;
                let memory = station.first_memory+k;
                let index = 2*m+memory;
                let c = positive(station.weight_m*branch.resistance_per_length_n_s_m2)?;
                q[index*n+index] = 1.0/memory_masses[memory];
                r[index*n+index] = c;
                add_outer(&mut r, n, m, &station.projection, c)?;
                for i in 0..m {
                    let cross = finite(-c*station.projection[i])?;
                    r[(m+i)*n+index] = cross;
                    r[index*n+m+i] = cross;
                }
            }
        }
        let storage = QuadraticStorage::new(q, n)?;
        poll(gate)?;
        let system = PortHamiltonian::new(n, m, j, r, g, Box::new(storage))?;
        poll(gate)?;
        Ok(Self { frame, basis, members, stations, material_ports, modal_mass,
            memory_masses, inverse_residual, system })
    }
}

impl<F> MaterialFluidFrame<'_, F> {
    pub const fn frame(&self) -> &Frame2d { self.frame }
    pub const fn basis(&self) -> &FrameVibrationBasis { &self.basis }
    pub fn members(&self) -> &[LoadedFrameMember<F>] { &self.members }
    /// Full loaded modal inertia in the ORIGINAL dry mass-normalized basis.
    pub fn modal_mass_matrix(&self) -> &[f64] { &self.modal_mass }
    pub fn memory_masses_kg(&self) -> &[f64] { &self.memory_masses }
    pub const fn mass_inverse_residual(&self) -> f64 { self.inverse_residual }
    /// Native force/velocity system. Ports are dry-basis generalized forces;
    /// use the retained basis for physical force/moment projection. Raw native
    /// stepping bypasses step_candidate's source-clock and independent work gates.
    pub const fn system(&self) -> &PortHamiltonian { &self.system }

    /// Initialize original modal displacements/velocities plus ALL physical
    /// quadrature fluid velocities. No fluid-rest assumption or history reset.
    pub fn state_from_motion(&self, modes: &[ModalAcousticState], fluid_velocities: &[f64]) -> Result<Vec<f64>> {
        let m = self.basis.modes().len();
        if modes.len() != m || fluid_velocities.len() != self.memory_masses.len()
            || modes.iter().any(|s| !s.displacement_m_sqrt_kg.is_finite() || !s.velocity_m_sqrt_kg_per_s.is_finite())
            || fluid_velocities.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite modal motion and all fluid velocities required"));
        }
        let velocity: Vec<_> = modes.iter().map(|s| s.velocity_m_sqrt_kg_per_s).collect();
        let mut state = vec![0.0; self.system.state_dim()];
        for i in 0..m {
            state[i] = modes[i].displacement_m_sqrt_kg;
            state[m+i] = dot(&self.modal_mass[i*m..(i+1)*m], &velocity)?;
        }
        for (k, (&mass, &v)) in self.memory_masses.iter().zip(fluid_velocities).enumerate() {
            state[2*m+k] = finite(mass*v)?;
            if v != 0.0 && state[2*m+k] == 0.0 { return Err(invalid("fluid momentum underflow")); }
        }
        self.admit_state(&state)?;
        Ok(state)
    }

    /// Recover original-basis q and qdot=M_inverse*p, not momentum-as-velocity.
    /// Feed these to basis().reconstruct_into for physical nodal observations.
    pub fn modal_states(&self, state: &[f64]) -> Result<Vec<ModalAcousticState>> {
        self.admit_state(state)?;
        let m = self.basis.modes().len();
        let effort = self.system.effort(state);
        Ok((0..m).map(|i| ModalAcousticState { displacement_m_sqrt_kg: state[i],
            velocity_m_sqrt_kg_per_s: effort[m+i] }).collect())
    }

    /// Advance held incremental nodal forces/moments using the ORIGINAL pHS
    /// stepper. Fixed supports and stationary far-field fluid are explicit.
    /// Do not apply the frozen equilibrium's preload twice. Cut force events at
    /// interval boundaries. Refine time, dry modes and fluid spectrum separately.
    ///
    /// The caller's state/source frame is never mutated. Publication follows
    /// original material clock admission, native solve, finite-state admission,
    /// and independent nodal-work / separate solid-and-fluid loss checks.
    pub fn step_candidate(&self, state: &[f64], nodal_force: &[f64], dt: f64,
        energy_tolerance_j: f64, gate: &CancelGate) -> Result<FrameFluidStep> {
        poll(gate)?;
        positive(dt)?;
        positive(energy_tolerance_j)?;
        self.admit_state(state)?;
        if let Some(damping) = self.basis.material_damping() { damping.admit_timestep(dt)?; }
        let m = self.basis.modes().len();
        let mut input = vec![0.0; m];
        self.basis.project_force_into(nodal_force, [0.0; 2], &mut input, || gate.is_requested())?;
        poll(gate)?;
        let native = fs_phs::step(&self.system, state, &input, dt)?;
        poll(gate)?;
        self.admit_state(&native.x)?;
        let midpoint: Vec<_> = state.iter().zip(&native.x).map(|(a,b)| 0.5*a+0.5*b).collect();
        let effort = self.system.effort(&midpoint);
        let velocity = &effort[m..2*m];
        let mut material_port_dissipation_j = Vec::with_capacity(self.material_ports.len());
        for port in &self.material_ports {
            poll(gate)?;
            let speed = dot(&port.projection, velocity)?;
            material_port_dissipation_j.push(finite(dt*port.resistance*speed*speed)?);
        }
        let mut fluid_member_dissipation_j = vec![0.0; self.members.len()];
        for station in &self.stations {
            poll(gate)?;
            let speed = dot(&station.projection, velocity)?;
            let flow = &self.members[station.member].flow;
            let mut power = station.weight_m*flow.low_tail_resistance_per_length_n_s_m2()*speed*speed;
            for (k, branch) in flow.branches().iter().enumerate() {
                let slip = speed-effort[2*m+station.first_memory+k];
                power += station.weight_m*branch.resistance_per_length_n_s_m2*slip*slip;
            }
            fluid_member_dissipation_j[station.member] = finite(
                fluid_member_dissipation_j[station.member]+dt*power)?;
        }
        let mut nodal_work_j = 0.0;
        for (dof, force) in nodal_force.iter().enumerate() {
            if dof % 128 == 0 { poll(gate)?; }
            let increment: f64 = self.basis.modes().iter().enumerate().map(|(i,mode)|
                mode.displacement[dof]*(native.x[i]-state[i])).sum();
            nodal_work_j = finite(nodal_work_j+force*increment)?;
        }
        let stored_before_j = self.system.hamiltonian(state);
        let stored_after_j = self.system.hamiltonian(&native.x);
        let loss: f64 = material_port_dissipation_j.iter().chain(&fluid_member_dissipation_j).sum();
        let energy_residual_j = finite(stored_after_j-stored_before_j+loss-nodal_work_j)?;
        let port_residual = finite(native.supplied-nodal_work_j)?;
        if energy_residual_j.abs() > energy_tolerance_j || port_residual.abs() > energy_tolerance_j {
            return Err(FrameFluidError::Balance {
                residual_j: if energy_residual_j.abs() >= port_residual.abs() { energy_residual_j } else { port_residual },
                tolerance_j: energy_tolerance_j,
            });
        }
        poll(gate)?;
        Ok(FrameFluidStep { native, stored_before_j, stored_after_j, material_port_dissipation_j,
            fluid_member_dissipation_j, nodal_work_j, energy_residual_j })
    }

    fn admit_state(&self, state: &[f64]) -> Result<()> {
        if state.len() != self.system.state_dim() || state.iter().any(|v| !v.is_finite()) {
            return Err(invalid("complete finite frame/fluid state required"));
        }
        finite(self.system.hamiltonian(state))?;
        if self.system.effort(state).iter().any(|v| !v.is_finite()) {
            return Err(invalid("unrepresentable frame/fluid effort"));
        }
        Ok(())
    }
}

fn project(basis: &FrameVibrationBasis, dofs: &[usize; 6], weights: &[f64; 6]) -> Result<Vec<f64>> {
    basis.modes().iter().map(|mode| finite(dofs.iter().zip(weights)
        .map(|(&dof,&weight)| weight*mode.displacement[dof]).sum())).collect()
}
fn add_outer(matrix: &mut [f64], n: usize, offset: usize, shape: &[f64], coefficient: f64) -> Result<()> {
    for i in 0..shape.len() {
        for k in 0..=i {
            let value = finite(coefficient*shape[i]*shape[k])?;
            let a = (offset+i)*n+offset+k;
            matrix[a] = finite(matrix[a]+value)?;
            if k != i { matrix[(offset+k)*n+offset+i] = matrix[a]; }
        }
    }
    Ok(())
}
fn inverse_mass(mass: &[f64], n: usize, tolerance: f64, gate: &CancelGate) -> Result<(Vec<f64>, f64)> {
    let factor = cholesky(mass, n).map_err(FrameFluidError::Factor)?;
    let mut inverse = vec![0.0; n*n];
    for column in 0..n {
        poll(gate)?;
        let mut rhs = vec![0.0; n];
        rhs[column] = 1.0;
        factor.solve(&mut rhs);
        for row in 0..n { inverse[row*n+column] = finite(rhs[row])?; }
    }
    // Retain an exactly symmetric Hamiltonian; check the actual averaged
    // inverse against the ORIGINAL full mass, rather than assuming symmetry
    // or silently accepting a poorly resolved factorization.
    for i in 0..n {
        for k in 0..i {
            let value = 0.5*inverse[i*n+k]+0.5*inverse[k*n+i];
            inverse[i*n+k] = value; inverse[k*n+i] = value;
        }
    }
    let mut residual = 0.0_f64;
    for i in 0..n {
        poll(gate)?;
        let mut row_sum = 0.0;
        for k in 0..n {
            let product: f64 = (0..n).map(|j| mass[i*n+j]*inverse[j*n+k]).sum();
            row_sum += finite(product-if i == k { 1.0 } else { 0.0 })?.abs();
        }
        residual = residual.max(finite(row_sum)?);
    }
    if residual > tolerance { return Err(invalid("loaded modal mass inverse exceeds the declared residual allowance")); }
    Ok((inverse, residual))
}

#[cfg(test)]
mod tests;
