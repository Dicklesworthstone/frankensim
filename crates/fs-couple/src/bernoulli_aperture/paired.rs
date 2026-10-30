//! Two independently moving blades sharing one fluid slit and contact interface.
//!
//! Each blade retains its own opening-positive coordinate, mass, stiffness,
//! loss and pressure area. No identical-blade or symmetric-motion assumption.
//! The common gap is H+x0+x1; pressure forces are -Ai*dp and swept flow is
//! -sum(Ai*vi). This is the same work-conjugate pressure/flow convention as
//! DynamicAperture. fs-dcontact owns the ONE gap-contact potential and loss;
//! moving owns the continuous Bernoulli law. Each blade's linear response over
//! a step, midpoint or exact zero-order hold ([`BladePropagator`]), is condensed
//! to an affine mean velocity `free + mobility * F` in the held force. That
//! feeds a scalar gap equation for the uniform image or a coupled
//! two-velocity equation for the spatial image, never a fictitious oscillator.
//! Antisymmetric motion therefore remains a real state. The exact blocks do not
//! warp frequency; the contact/Bernoulli forces are still held over each step
//! and solved at the midpoint opening.
//!
//! Fixed-state one-mode blades with a uniform slit or an explicitly selected
//! matching-planar spatial interface. The latter retains separate nodal gaps
//! and contact forces and solves both midpoint blade velocities. This does
//! not add curved-shell modes, spatial material evolution, staple pressure
//! recovery, or measured oboe/bassoon accuracy. Numerical root selection is
//! deterministic, not a uniqueness theorem. Midpoint accuracy is second order
//! only on smooth branches. The existing symmetric single-coordinate image
//! remains available and unchanged. Optional independent regional Maxwell laws
//! use the same fs-material/fs-plate/fs-phs projection and update as the scalar
//! plate valve. Material history participates in the selected linear block and
//! publishes only with both blades and the connected acoustic network. The exact
//! block integrates each blade's own Maxwell force, storage and dissipation; the
//! pressure/contact forces are still held over a step and require time refinement.

/// Spatial slit/contact coupled to the original independent blade mechanics.
mod spatial;

use super::plate::closure::paired::PairedPlateClosure;
use super::dynamic::{ApertureDrive, ApertureFrame, ApertureState};
use super::dynamic::relaxation::{
    ApertureRelaxation, InitialApertureMemory, MemoryTrial, PlateRelaxationSpec,
};
use super::plate::{PlateApertureMotion, PlateApertureReduction};
use crate::acoustic_realize::AcousticRealizeError;
use crate::unilateral_contact::distributed::ApertureContactStep;
use fs_dcontact::{ContactStorage, Obstacle};
use fs_material::gas::GasState;
use fs_phs::Storage;
use std::rc::Rc;

/// One effective blade in its own opening-positive mean-edge coordinate.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ApertureBlade {
    /// Effective mass [kg].
    pub mass_kg: f64,
    /// Linear restoring stiffness [N/m].
    pub stiffness_n_m: f64,
    /// Explicit viscous damping coefficient [N s/m].
    pub damping_n_s_m: f64,
    /// Positive uniform-pressure area, also used for swept volume [m²].
    pub pressure_area_m2: f64,
}

/// Physical displacement from the supplied rest shape, not half the total gap.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct BladeState {
    /// Displacement in this blade's opening-positive coordinate [m].
    pub displacement_m: f64,
    /// Velocity in that same coordinate [m/s].
    pub velocity_m_s: f64,
}

/// Linear blade propagator, always chosen explicitly by the caller.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BladePropagator {
    /// Implicit midpoint on each blade's linear block. It warps frequency, so it
    /// is admitted only inside the one-cent band (DEF-04), and it carries the
    /// optional Maxwell material memory.
    Midpoint,
    /// Exact zero-order-hold exponential of the mass/damping/stiffness block,
    /// with the step's contact, pressure and external forces held over the
    /// step. No frequency warping at any omega*dt; the nonlinear contact/
    /// Bernoulli junction solve is unchanged. Each original blade's optional
    /// Maxwell arms evolve inside its own block, without a split memory update.
    ExactZoh,
}

/// One Maxwell arm on a scalar coordinate x: energy `(c x - z)^2 / 2`,
/// memory `dz/dt = (c x - z) / tau`, restoring force `-c (c x - z)`.
#[derive(Debug, Clone, Copy, PartialEq)]
pub(crate) struct ExactArm {
    /// `sqrt(k_arm) * projection` [sqrt(N/m)].
    pub(crate) coupling: f64,
    /// Relaxation time [s].
    pub(crate) relaxation_time_s: f64,
}

/// Prepared exact maps of one blade over one step. The scaled state
/// `[omega x, v, F/(m omega), w_1..w_n, omega int v, omega^2 int x]` with
/// `w_j = omega z_j / c_j` keeps every generator entry at omega*dt, c*dt/m,
/// dt/tau or (k_j/k) omega*dt. The `int v` row gives the displacement
/// increment directly, without subtracting nearly equal values. With arms, one
/// extra row integrates their restoring force divided by mass, avoiding an
/// endpoint subtraction when reporting the step-mean material force.
#[derive(Debug, Clone, PartialEq)]
pub(crate) struct ExactBlade {
    omega: f64,
    mass: f64,
    arms: Vec<ExactArm>,
    // Square map, side 5 + arms plus one material-force integral when nonempty.
    map: Vec<f64>,
    // Exact c*int v^2 and arm int (c x - z)^2/tau over the step, as quadratic
    // forms on the leading 3 + arms (state and held force) coordinates.
    viscous: Vec<f64>,
    arm_loss: Vec<f64>,
}

impl ExactBlade {
    fn new(blade: &ApertureBlade, dt: f64) -> Result<Self, AcousticRealizeError> {
        Self::from_mechanics(blade.mass_kg, blade.damping_n_s_m, blade.stiffness_n_m, dt)
    }
    /// One mass/viscous-damper/spring coordinate with positive mass and stiffness.
    pub(crate) fn from_mechanics(mass: f64, damping_n_s_m: f64, stiffness_n_m: f64, dt: f64)
        -> Result<Self, AcousticRealizeError>
    {
        Self::with_arms(mass, damping_n_s_m, stiffness_n_m, &[], dt)
    }
    /// The same coordinate with its Maxwell arms inside the exact block.
    pub(crate) fn with_arms(mass: f64, damping_n_s_m: f64, stiffness_n_m: f64, arms: &[ExactArm], dt: f64)
        -> Result<Self, AcousticRealizeError>
    {
        let refuse = |_| invalid("exact blade propagator could not be prepared");
        if !(mass > 0.0 && stiffness_n_m > 0.0 && damping_n_s_m >= 0.0 && dt > 0.0)
            || ![mass, damping_n_s_m, stiffness_n_m, dt].iter().all(|x| x.is_finite())
            || arms.iter().any(|a| !(a.coupling.is_finite() && a.coupling != 0.0
                && a.relaxation_time_s.is_finite() && a.relaxation_time_s > 0.0)) {
            return Err(invalid("exact blade propagator needs positive mass, stiffness, step and arms"));
        }
        let omega = fs_math::det::sqrt(stiffness_n_m / mass);
        let (w, g) = (omega * dt, damping_n_s_m / mass * dt);
        let a = arms.len();
        // Preserve the original no-memory matrix and its arithmetic exactly.
        let (core, n) = (3 + a, 5 + a + usize::from(a != 0));
        let mut generator = vec![0.0; n * n];
        generator[1] = w;
        generator[n] = -w;
        generator[n + 1] = -g;
        generator[n + 2] = w;
        for (j, arm) in arms.iter().enumerate() {
            let beta = arm.coupling * arm.coupling / (mass * omega) * dt;
            let rate = dt / arm.relaxation_time_s;
            generator[n] -= beta;
            generator[n + 3 + j] = beta;
            generator[(3 + j) * n] = rate;
            generator[(3 + j) * n + 3 + j] = -rate;
            // dI/dt = sum(c*(c*x-z))/m. The stored scaled coordinates are
            // omega*x and omega*z/c, so beta is the same momentum coefficient.
            generator[(5 + a) * n] += beta;
            generator[(5 + a) * n + 3 + j] = -beta;
        }
        generator[(3 + a) * n + 1] = w;
        generator[(4 + a) * n] = w;
        let map = fs_phs::linear_exact::exponential(&generator, n).map_err(refuse)?;
        // Van Loan: exp([[-A^T, W], [0, A]] dt) = [[., F12], [0, F22]], and
        // int_0^dt e^{A^T s} W e^{A s} ds = F22^T F12 on the closed core block.
        let quadratic = |weights: &[f64]| -> Result<Vec<f64>, AcousticRealizeError> {
            let side = 2 * core;
            let mut van_loan = vec![0.0; side * side];
            for i in 0..core {
                for j in 0..core {
                    van_loan[i * side + j] = -generator[j * n + i];
                    van_loan[(i + core) * side + j + core] = generator[i * n + j];
                    van_loan[i * side + j + core] = weights[i * core + j] * dt;
                }
            }
            let blocks = fs_phs::linear_exact::exponential(&van_loan, side).map_err(refuse)?;
            let mut form = vec![0.0; core * core];
            for i in 0..core {
                for j in 0..core {
                    form[i * core + j] = (0..core)
                        .map(|l| blocks[(l + core) * side + i + core] * blocks[l * side + j + core])
                        .sum();
                }
            }
            Ok(form)
        };
        let mut viscous_weights = vec![0.0; core * core];
        viscous_weights[core + 1] = damping_n_s_m;
        let mut arm_weights = vec![0.0; core * core];
        for (j, arm) in arms.iter().enumerate() {
            // (c x - z)^2 / tau = (c/omega)^2 (omega x - w)^2 / tau.
            let k = arm.coupling * arm.coupling / (omega * omega * arm.relaxation_time_s);
            arm_weights[0] += k;
            arm_weights[3 + j] -= k;
            arm_weights[(3 + j) * core] -= k;
            arm_weights[(3 + j) * core + 3 + j] += k;
        }
        let (viscous, arm_loss) = (quadratic(&viscous_weights)?, quadratic(&arm_weights)?);
        if !(omega.is_finite() && omega > 0.0)
            || !map.iter().chain(&viscous).chain(&arm_loss).all(|x| x.is_finite()) {
            return Err(invalid("exact blade propagator is not representable"));
        }
        Ok(Self { omega, mass, arms: arms.to_vec(), map, viscous, arm_loss })
    }
    fn side(&self) -> usize { 5 + self.arms.len() + usize::from(!self.arms.is_empty()) }
    // Leading state: scaled displacement, velocity, held force and arm memories.
    fn scaled(&self, state: BladeState, force_n: f64, memory: &[f64]) -> Vec<f64> {
        assert_eq!(memory.len(), self.arms.len(), "exact blade memory must name every arm");
        let mut z = vec![self.omega * state.displacement_m, state.velocity_m_s,
            force_n / (self.mass * self.omega)];
        z.extend(self.arms.iter().zip(memory).map(|(arm, m)| self.omega * m / arm.coupling));
        z
    }
    fn row(&self, row: usize, z: &[f64]) -> f64 {
        let n = self.side();
        z.iter().enumerate().map(|(j, value)| self.map[row * n + j] * value).sum()
    }
    // Read the same scaled coordinates without allocating a temporary state.
    // The contact island evaluates these rows at every trial velocity.
    fn coordinate(&self, state: BladeState, force_n: f64, memory: &[f64], j: usize) -> f64 {
        match j {
            0 => self.omega * state.displacement_m,
            1 => state.velocity_m_s,
            2 => force_n / (self.mass * self.omega),
            _ => self.omega * memory[j - 3] / self.arms[j - 3].coupling,
        }
    }
    fn row_state(&self, row: usize, state: BladeState, force_n: f64, memory: &[f64]) -> f64 {
        assert_eq!(memory.len(), self.arms.len(), "exact blade memory must name every arm");
        let n = self.side();
        (0..3 + self.arms.len()).map(|j|
            self.map[row * n + j] * self.coordinate(state, force_n, memory, j)).sum()
    }
    fn form_state(&self, matrix: &[f64], state: BladeState, force_n: f64, memory: &[f64]) -> f64 {
        assert_eq!(memory.len(), self.arms.len(), "exact blade memory must name every arm");
        let core = 3 + self.arms.len();
        (0..core).map(|i| {
            let left = self.coordinate(state, force_n, memory, i);
            (0..core).map(|j| left * matrix[i * core + j]
                * self.coordinate(state, force_n, memory, j)).sum::<f64>()
        }).sum()
    }
    fn form(&self, matrix: &[f64], z: &[f64]) -> f64 {
        let core = z.len();
        (0..core).map(|i| (0..core).map(|j| z[i] * matrix[i * core + j] * z[j]).sum::<f64>()).sum()
    }
    /// Mean velocity `dx/dt` over the step is `free + mobility * F` for held F.
    pub(crate) fn free_and_mobility(&self, state: BladeState, external_n: f64, dt: f64, memory: &[f64])
        -> (f64, f64)
    {
        let eta = 3 + self.arms.len();
        let free = self.row_state(eta, state, external_n, memory) / (self.omega * dt);
        let mobility = self.map[eta * self.side() + 2] / (self.mass * self.omega * self.omega * dt);
        (free, mobility)
    }
    pub(crate) fn end_velocity(&self, state: BladeState, held_n: f64, memory: &[f64]) -> f64 {
        self.row_state(1, state, held_n, memory)
    }
    pub(crate) fn mean_displacement(&self, state: BladeState, held_n: f64, dt: f64, memory: &[f64]) -> f64 {
        self.row(4 + self.arms.len(), &self.scaled(state, held_n, memory)) / (self.omega * self.omega * dt)
    }
    pub(crate) fn viscous_loss_j(&self, state: BladeState, held_n: f64, memory: &[f64]) -> f64 {
        self.form_state(&self.viscous, state, held_n, memory)
    }
    /// End-of-step arm memories z_j, energy-normalized like the arms' own owner.
    pub(crate) fn end_memory(&self, state: BladeState, held_n: f64, memory: &[f64]) -> Vec<f64> {
        let z = self.scaled(state, held_n, memory);
        self.arms.iter().enumerate().map(|(j, arm)| arm.coupling * self.row(3 + j, &z) / self.omega).collect()
    }
    /// Step-mean Maxwell restoring force [N], integrated in the same block.
    /// This is not the endpoint or midpoint force and does not subtract nearly
    /// equal end memories. Empty spectra preserve the no-memory arithmetic.
    fn mean_memory_force_n(&self, state: BladeState, held_n: f64, dt: f64, memory: &[f64]) -> f64 {
        if self.arms.is_empty() { return 0.0; }
        self.mass / dt * self.row(5 + self.arms.len(), &self.scaled(state, held_n, memory))
    }
    /// Exact arm dissipation over the step [J].
    pub(crate) fn arm_loss_j(&self, state: BladeState, held_n: f64, memory: &[f64]) -> f64 {
        self.form(&self.arm_loss, &self.scaled(state, held_n, memory))
    }
}

/// Explicit physical/solver allowances; no response is clipped to these bounds.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairedApertureLimits {
    /// Bound on each blade's displacement from its own rest shape [m].
    pub maximum_abs_displacement_m: f64,
    /// Bound on each blade's endpoint speed [m/s].
    pub maximum_abs_velocity_m_s: f64,
    /// Bound on combined mechanical, material-memory and contact storage [J].
    pub maximum_energy_j: f64,
    /// Total nonlinear trial evaluations, including bracketing, in 1..=512.
    pub max_evaluations: usize,
    /// Absolute residual allowance for uniform total-gap velocity or EACH
    /// independently solved spatial blade midpoint velocity [m/s].
    pub velocity_tolerance_m_s: f64,
    /// Relative residual allowance, in (0,1).
    pub relative_tolerance: f64,
}

/// Shared geometry/fluid/clock inputs and independent mechanical data.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PairedApertureSpec {
    /// Total rest clearance of the single shared slit [m].
    pub rest_opening_m: f64,
    /// Width of the shared slit [m], not the sum of both blade widths.
    pub width_m: f64,
    /// Independent blade parameters; authored scalars carry no source authority.
    pub blades: [ApertureBlade; 2],
    /// Fluid density [kg/m³].
    pub density_kg_m3: f64,
    /// Bore pressure/volume-flow impedance [Pa s/m³].
    pub impedance_pa_s_m3: f64,
    /// Positive fixed mechanical interval [s].
    pub time_step_s: f64,
    /// Total accepted samples, at most 2^53, including resumed windows.
    pub max_steps: u64,
    /// Caller-selected physical and numerical limits.
    pub limits: PairedApertureLimits,
}

/// Same-clock observations without erasing either blade's state.
#[derive(Debug, Clone, Copy, Default, PartialEq)]
pub struct PairedApertureFrame {
    /// Common fluid/contact work account. State is the reference H+x0+x1
    /// coordinate; on a spatial interface use its probe for actual local gaps.
    pub aperture: ApertureFrame,
    /// Actual individual endpoint displacements and velocities.
    pub blades: [BladeState; 2],
    /// Midpoint upstream-minus-bore pressure [Pa].
    pub pressure_drop_pa: f64,
    /// Shared uniform force, or total physical normal resultant of a spatial
    /// interface [N]. A spatial resultant is not either generalized blade force.
    pub contact_force_n: f64,
    /// Actual work-conjugate contact forces on the independent blade modes [N].
    pub blade_contact_forces_n: [f64; 2],
    /// Independent blade momentum residuals in force units [N].
    pub momentum_residual_n: [f64; 2],
}

/// Stateful paired mechanical source; trials never advance accepted state.
/// Clones copy both blades and their complete material histories; immutable
/// contact potentials and reference interface geometry are shared.
#[derive(Clone)]
pub struct PairedAperture {
    spec: PairedApertureSpec,
    blades: [BladeState; 2],
    contact: Obstacle,
    contact_storage: Rc<ContactStorage>,
    plates: Option<Box<[PlateApertureReduction; 2]>>,
    gas: Option<GasState>,
    viscous_slit: Option<super::viscous::ViscousSlit>,
    relaxation: [Option<ApertureRelaxation>; 2],
    closure: Option<Rc<PairedPlateClosure>>,
    exact: Option<[ExactBlade; 2]>,
    samples: u64,
}

fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn nonlinear(e: impl core::fmt::Display) -> AcousticRealizeError {
    AcousticRealizeError::Nonlinear(e.to_string())
}
struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 { 0.0 }
    fn gradient(&self, _: &[f64], out: &mut [f64]) { out.fill(0.0); }
}
#[derive(Clone, Copy)]
struct Trial { velocity: [f64; 2], dp: f64, residual: f64, scale: f64 }

// Complete immutable candidate. Public observations remain Copy, but material
// histories cannot be accidentally discarded when a network accepts a frame.
pub(super) struct PairedTrial {
    pub(super) frame: PairedApertureFrame,
    memory: [Option<MemoryTrial>; 2],
}

impl PairedAperture {
    /// Admit explicit mechanics and ONE scalar gap obstacle. For opposed blades
    /// its physical law must describe their interface, not two separate lays.
    /// Initial gap penetration is allowed only as stored compliant contact.
    /// Zero contact stiffness is an explicit no-contact comparison.
    pub fn new(spec: PairedApertureSpec, blades: [BladeState; 2], contact: Obstacle)
        -> Result<Self, AcousticRealizeError>
    {
        Self::new_with_propagator(spec, blades, contact, BladePropagator::Midpoint)
    }

    /// [`Self::new`] with an explicitly selected linear blade propagator.
    /// `Midpoint` keeps the one-cent band refusal; `ExactZoh` admits any
    /// omega*dt because the elastic blade blocks are propagated exactly.
    pub fn new_with_propagator(spec: PairedApertureSpec, blades: [BladeState; 2], contact: Obstacle,
        propagator: BladePropagator) -> Result<Self, AcousticRealizeError>
    {
        let l = spec.limits;
        if ![spec.rest_opening_m, spec.width_m, spec.density_kg_m3,
            spec.impedance_pa_s_m3, spec.time_step_s, l.maximum_abs_displacement_m,
            l.maximum_abs_velocity_m_s, l.maximum_energy_j, l.velocity_tolerance_m_s,
            l.relative_tolerance].iter().all(|x| x.is_finite() && *x > 0.0)
            || l.relative_tolerance >= 1.0 || !(1..=512).contains(&l.max_evaluations)
            || spec.max_steps == 0 || spec.max_steps > (1_u64 << 53)
            || !(spec.time_step_s * spec.max_steps as f64).is_finite()
            || contact.provenance().trim().is_empty()
        { return Err(invalid("paired aperture needs explicit finite physical, time and solve limits")); }
        for b in &spec.blades {
            if ![b.mass_kg, b.stiffness_n_m, b.pressure_area_m2].iter().all(|x| x.is_finite() && *x > 0.0)
                || !b.damping_n_s_m.is_finite() || b.damping_n_s_m < 0.0 {
                return Err(invalid("each aperture blade needs positive mass, stiffness, area and nonnegative viscous loss"));
            }
        }
        let gap = spec.rest_opening_m + blades[0].displacement_m + blades[1].displacement_m;
        ApertureContactStep::new(&contact, gap)
            .and_then(|c| c.response(gap, blades[0].velocity_m_s + blades[1].velocity_m_s))
            .map_err(nonlinear)?;
        let contact_storage = ContactStorage::new(Box::new(ZeroStorage), 1, vec![contact.clone()])
            .map_err(nonlinear)?;
        let exact = match propagator {
            BladePropagator::Midpoint => None,
            BladePropagator::ExactZoh => Some([ExactBlade::new(&spec.blades[0], spec.time_step_s)?,
                ExactBlade::new(&spec.blades[1], spec.time_step_s)?]),
        };
        let model = Self { spec, blades, contact, contact_storage: contact_storage.into(), plates: None, gas: None, viscous_slit: None, relaxation: [None, None], closure: None, exact, samples: 0 };
        model.admit_linear_clock()?;
        model.validate_state(blades)?;
        Ok(model)
    }

    /// Compose two independently reduced source plates. No coefficient is
    /// averaged, inferred from the other blade, or selected by a material name.
    /// Each reduction must be a single blade. Both slit widths and declared
    /// total rest gaps must agree; they are not doubled. Material receipts,
    /// source shapes and individual linear-slope limits remain attached.
    #[allow(clippy::too_many_arguments)]
    pub fn from_plates(plates: [PlateApertureReduction; 2], gas: GasState,
        impedance_pa_s_m3: f64, time_step_s: f64, max_steps: u64,
        limits: PairedApertureLimits, states: [BladeState; 2], contact: Obstacle)
        -> Result<Self, AcousticRealizeError>
    {
        Self::from_plates_with_propagator(plates, gas, impedance_pa_s_m3, time_step_s,
            max_steps, limits, states, contact, BladePropagator::Midpoint)
    }

    /// Retain both source plates with an explicitly chosen structural propagator.
    /// Exact propagation carries each attached material spectrum inside its own
    /// block; it does not remove the source band, time-resolution, geometry or
    /// nonlinear junction limits. No plate reduction or material law is replaced.
    #[allow(clippy::too_many_arguments)]
    pub fn from_plates_with_propagator(plates: [PlateApertureReduction; 2], gas: GasState,
        impedance_pa_s_m3: f64, time_step_s: f64, max_steps: u64,
        limits: PairedApertureLimits, states: [BladeState; 2], contact: Obstacle,
        propagator: BladePropagator) -> Result<Self, AcousticRealizeError>
    {
        gas.validate_coherence().map_err(|_| invalid("paired aperture requires one coherent source gas"))?;
        let h = plates[0].options().rest_opening_m;
        let width = plates[0].width_m();
        if plates.iter().any(|p| p.moving_blade_count() != 1)
            || plates[1].options().rest_opening_m.to_bits() != h.to_bits()
            || (plates[1].width_m() - width).abs() > 64.0 * f64::EPSILON * width {
            return Err(invalid("independent blade pairing requires two single-blade charts with the same slit width and total rest gap"));
        }
        let blade = |p: &PlateApertureReduction| ApertureBlade {
            mass_kg: p.mass_kg(), stiffness_n_m: p.stiffness_n_m(),
            damping_n_s_m: 2.0 * p.damping_ratio()
                * fs_math::det::sqrt(p.mass_kg() * p.stiffness_n_m()),
            pressure_area_m2: p.pressure_area_m2(),
        };
        let spec = PairedApertureSpec { rest_opening_m: h, width_m: width,
            blades: [blade(&plates[0]), blade(&plates[1])], density_kg_m3: gas.density,
            impedance_pa_s_m3, time_step_s, max_steps, limits };
        let mut model = Self::new_with_propagator(spec, states, contact, propagator)?;
        model.plates = Some(Box::new(plates));
        model.gas = Some(gas);
        model.validate_state(states)?;
        Ok(model)
    }

    /// Select one common viscous channel before the first accepted step. Both
    /// independent blades share its pressure and flow; width is counted once.
    /// A source-bound law retains its full gas for subsequent network admission;
    /// it must agree with any gas already supplied by the material construction.
    pub fn set_viscous_slit(&mut self, law: super::viscous::ViscousSlit)
        -> Result<(), AcousticRealizeError>
    {
        if self.samples != 0 {
            return Err(invalid("paired viscous slit must be selected before the first accepted sample"));
        }
        match &self.closure {
            Some(closure) => { closure.viscous_flow(law, self.blades.map(|s| s.displacement_m),
                self.spec.density_kg_m3, 0.0)?; }
            None => { law.evaluate(self.spec.width_m, self.gap(self.blades).max(0.0), self.spec.density_kg_m3, 0.0)?; }
        }
        // Preserve the complete source binding for connected network admission.
        // Density alone cannot identify viscosity, heat capacity or temperature.
        if let Some(gas) = law.gas().copied() {
            if self.gas.is_some_and(|retained| retained != gas) {
                return Err(invalid("paired viscous slit must use the complete retained source gas"));
            }
            self.gas = Some(gas);
        }
        self.viscous_slit = Some(law);
        Ok(())
    }

    /// Derive viscosity from the retained gas, preserving its material identity.
    /// The shared streamwise channel length is supplied geometry, not a preset.
    pub fn set_material_viscous_slit(&mut self, length_m: f64) -> Result<(), AcousticRealizeError> {
        let gas = self.gas.as_ref().ok_or_else(|| invalid("paired channel requires its retained material gas"))?;
        self.set_viscous_slit(super::viscous::ViscousSlit::from_gas(gas, length_m)?)
    }

    pub(super) fn rebind_exact_gas_load(&mut self, gas: &GasState, impedance_pa_s_m3: f64)
        -> Result<(), AcousticRealizeError>
    {
        gas.validate_coherence().map_err(|_| invalid("live paired aperture requires one coherent source gas"))?;
        if self.exact.is_none() || !impedance_pa_s_m3.is_finite() || impedance_pa_s_m3 <= 0.0 {
            return Err(invalid("live paired gas refresh requires ExactZoh mechanics and positive finite impedance"));
        }
        let law = self.viscous_slit.map(|law| law.rebind_source_gas(gas)).transpose()?;
        if let Some(law) = law {
            match &self.closure {
                Some(closure) => { closure.viscous_flow(law, self.blades.map(|s| s.displacement_m), gas.density, 0.0)?; }
                None => { law.evaluate(self.spec.width_m, self.gap(self.blades).max(0.0), gas.density, 0.0)?; }
            }
        }
        // Both blades and all solid/contact memories keep their original
        // coordinates; the quasisteady gas law has no stored energy to remap.
        self.spec.density_kg_m3 = gas.density;
        self.spec.impedance_pa_s_m3 = impedance_pa_s_m3;
        self.viscous_slit = law;
        self.gas = Some(*gas);
        Ok(())
    }

    /// Original explicit or material-bound channel law, shared by both blades.
    #[must_use]
    pub const fn viscous_slit(&self) -> Option<super::viscous::ViscousSlit> { self.viscous_slit }

    /// Attach regional hereditary bending to ONE original source blade, before
    /// any accepted step. Call separately for the other blade; its material,
    /// history and even branch count may differ, or it may remain elastic.
    ///
    /// The existing plate projection checks equilibrium sections, triangle
    /// coverage, material bands and time resolution. Only bending is relaxed;
    /// mass, prestress, pressure area and the shared contact law are unchanged.
    /// This blade must have zero separate viscous damping to avoid counting the
    /// same loss twice. No spectra or history are inferred from a material name.
    /// Initial viscous displacements use this blade's OWN coordinate in metres.
    ///
    /// # Errors
    /// Invalid blade index, missing source plate, repeated/late attachment,
    /// separate damping, incompatible material/clock/history, or combined energy
    /// above the explicit budget. No active history is ever reset or transferred.
    pub fn with_blade_relaxation(mut self, blade: usize, spec: PlateRelaxationSpec,
        initial: InitialApertureMemory) -> Result<Self, AcousticRealizeError>
    {
        if blade>=2 || self.samples!=0 || self.relaxation[blade].is_some()
            || self.spec.blades[blade].damping_n_s_m!=0.0 {
            return Err(invalid("blade relaxation requires a valid cold blade, one attachment and zero separate damping"));
        }
        let plate=self.plates.as_ref().ok_or_else(||invalid("blade relaxation requires its original physical plate"))?;
        let memory=ApertureRelaxation::from_plate(&plate[blade],self.spec.time_step_s,
            self.blades[blade].displacement_m,spec,initial)?;
        self.validate_energy(self.stored_energy_j()+memory.stored_energy_j())?;
        if let Some(exact) = &mut self.exact {
            let b = self.spec.blades[blade];
            exact[blade] = ExactBlade::with_arms(b.mass_kg, b.damping_n_s_m,
                b.stiffness_n_m, &memory.exact_arms()?, self.spec.time_step_s)?;
        }
        self.relaxation[blade]=Some(memory);
        Ok(self)
    }

    /// Accepted material history and original region laws for one blade.
    /// An elastic blade or out-of-range index has no material-memory state.
    #[must_use]
    pub fn blade_relaxation(&self, blade: usize) -> Option<&ApertureRelaxation> {
        self.relaxation.get(blade).and_then(Option::as_ref)
    }

    /// Original mechanical/fluid/clock choices.
    #[must_use]
    pub const fn spec(&self) -> &PairedApertureSpec { &self.spec }
    /// The linear blade propagator selected at construction.
    #[must_use]
    pub const fn blade_propagator(&self) -> BladePropagator {
        if self.exact.is_some() { BladePropagator::ExactZoh } else { BladePropagator::Midpoint }
    }
    /// Individual accepted physical states, in original blade order.
    #[must_use]
    pub const fn blade_states(&self) -> &[BladeState; 2] { &self.blades }
    /// Exact independent plate sources, absent on the authored-scalar path.
    #[must_use]
    pub fn plates(&self) -> Option<&[PlateApertureReduction; 2]> { self.plates.as_deref() }
    /// Complete source gas retained by material geometry or a source-bound channel.
    #[must_use]
    pub const fn gas(&self) -> Option<&GasState> { self.gas.as_ref() }
    /// Accepted sample count; refusals never advance it.
    #[must_use]
    pub const fn accepted_steps(&self) -> u64 { self.samples }
    /// Combined mechanical, independent material-memory and shared-interface energy [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.energy(self.blades)+self.relaxation.iter().flatten()
            .map(ApertureRelaxation::stored_energy_j).sum::<f64>()
    }
    fn gap(&self, state: [BladeState; 2]) -> f64 {
        self.spec.rest_opening_m + state[0].displacement_m + state[1].displacement_m
    }
    fn energy(&self, state: [BladeState; 2]) -> f64 {
        // Public storage observes accepted states. During a numerical trial an
        // unrepresentable spatial potential remains nonfinite and is refused
        // by validate_energy; it is never silently replaced with zero storage.
        let mut energy = match &self.closure {
            Some(closure) => closure.stored_energy_j(state.map(|s|s.displacement_m)).unwrap_or(f64::NAN),
            None => self.contact_storage.hamiltonian(&[self.gap(state), 0.0]),
        };
        for (s, b) in state.iter().zip(&self.spec.blades) {
            energy += 0.5*b.mass_kg*s.velocity_m_s*s.velocity_m_s
                + 0.5*b.stiffness_n_m*s.displacement_m*s.displacement_m;
        }
        energy
    }
    fn validate_state(&self, state: [BladeState; 2]) -> Result<(), AcousticRealizeError> {
        for (i, s) in state.iter().enumerate() {
            if !s.displacement_m.is_finite() || !s.velocity_m_s.is_finite()
                || s.displacement_m.abs() > self.spec.limits.maximum_abs_displacement_m
                || s.velocity_m_s.abs() > self.spec.limits.maximum_abs_velocity_m_s {
                return Err(invalid("paired blade state exceeds its physical bounds"));
            }
            if let Some(plates) = &self.plates {
                plates[i].validate_opening(plates[i].options().rest_opening_m + s.displacement_m)?;
            }
        }
        if let Some(closure) = &self.closure {
            closure.validate_displacements(state.map(|s|s.displacement_m))?;
        }
        self.validate_energy(self.energy(state))
    }
    fn validate_energy(&self, energy: f64) -> Result<(), AcousticRealizeError> {
        if !energy.is_finite() || energy < 0.0 || energy > self.spec.limits.maximum_energy_j {
            return Err(invalid("paired aperture storage exceeds its finite energy budget"));
        }
        Ok(())
    }
    /// Observe either blade on its original mesh, without inventing symmetry.
    pub fn nodal_motion(&self, blade: usize, node: usize) -> Result<PlateApertureMotion, AcousticRealizeError> {
        let plate = self.plates.as_ref().and_then(|p| p.get(blade))
            .ok_or_else(|| invalid("paired nodal motion requires a retained source blade"))?;
        let state = self.blades[blade];
        plate.nodal_motion(node, plate.options().rest_opening_m + state.displacement_m, state.velocity_m_s)
    }
    /// Increase the horizon without changing either blade or contact state.
    pub fn extend_step_budget(&mut self, total: u64) -> Result<(), AcousticRealizeError> {
        if total <= self.spec.max_steps || total > (1_u64 << 53)
            || !(total as f64*self.spec.time_step_s).is_finite() {
            return Err(invalid("paired aperture budget extension must increase a finite exact-step horizon"));
        }
        self.spec.max_steps = total;
        Ok(())
    }
    /// Advance pressure, both blades and their shared contact as one transaction.
    pub fn step(&mut self, drive: ApertureDrive) -> Result<PairedApertureFrame, AcousticRealizeError> {
        self.step_with_forces(drive, [0.0; 2])
    }
    /// Independent opening-positive held forces [N]. Each supplies its OWN
    /// F_i*v_i*dt work; a load is not silently duplicated onto the other blade.
    pub fn step_with_forces(&mut self, drive: ApertureDrive, forces_n: [f64; 2])
        -> Result<PairedApertureFrame, AcousticRealizeError>
    {
        let trial = self.preview(drive, forces_n)?;
        let frame = trial.frame;
        self.accept(trial);
        Ok(frame)
    }
    pub(super) fn accept(&mut self, trial: PairedTrial) {
        for (memory, next) in self.relaxation.iter_mut().zip(trial.memory) {
            if let (Some(memory), Some(next)) = (memory, next) { memory.accept(next); }
        }
        self.blades = trial.frame.blades;
        self.samples = trial.frame.aperture.step;
    }

    pub(super) fn preview(&self, drive: ApertureDrive, forces: [f64; 2])
        -> Result<PairedTrial, AcousticRealizeError>
    {
        self.admit_linear_clock()?;
        if self.samples >= self.spec.max_steps { return Err(invalid("paired aperture step budget exhausted")); }
        if ![drive.upstream_pressure_pa, drive.incoming_pressure_pa, drive.body_flow_m3_s,
            forces[0], forces[1]].iter().all(|x| x.is_finite()) {
            return Err(invalid("paired aperture needs finite pressure, flow and independent force inputs"));
        }
        let dt = self.spec.time_step_s;
        let h = 0.5*dt;
        let b = self.spec.blades;
        let old = self.blades;
        let gap = self.gap(old);
        let contact = ApertureContactStep::new(&self.contact, gap).map_err(nonlinear)?;
        let mut free = [0.0; 2];
        let mut mobility = [0.0; 2];
        let arm_memory: [&[f64]; 2] = std::array::from_fn(|i|
            self.relaxation[i].as_ref().map_or(&[][..], ApertureRelaxation::memory_sqrt_j));
        for i in 0..2 {
            if let Some(exact) = &self.exact {
                // Exact held-force map: the mean velocity dx/dt stays affine in the
                // step's force, so the junction solve below is unchanged.
                (free[i], mobility[i]) = exact[i].free_and_mobility(old[i], forces[i], dt, arm_memory[i]);
                if !free[i].is_finite() || !mobility[i].is_finite() || mobility[i] <= 0.0 {
                    return Err(invalid("paired blade exact response is not representable"));
                }
                continue;
            }
            let mut d = b[i].mass_kg + h*b[i].damping_n_s_m + h*h*b[i].stiffness_n_m;
            let mut numerator = b[i].mass_kg*old[i].velocity_m_s
                - h*b[i].stiffness_n_m*old[i].displacement_m + h*forces[i];
            if let Some(memory) = self.relaxation[i].as_ref().filter(|m| !m.branches().is_empty()) {
                let (bias, resistance) = memory.midpoint_response(old[i].displacement_m, dt)?;
                d += h*resistance;
                numerator -= h*bias;
            }
            mobility[i] = h/d;
            free[i] = numerator/d;
            if !d.is_finite() || d <= 0.0 || !mobility[i].is_finite() || mobility[i] <= 0.0
                || !free[i].is_finite() { return Err(invalid("paired blade midpoint response is not representable")); }
        }
        let z = self.spec.impedance_pa_s_m3;
        let trial = if let Some(closure) = &self.closure {
            self.solve_profiled(drive, free, mobility, closure)?
        } else {
            let pressure_mobility = (0..2).map(|i| mobility[i]*b[i].pressure_area_m2*b[i].pressure_area_m2).sum::<f64>();
            let denominator = 1.0 + z*pressure_mobility;
            if !denominator.is_finite() { return Err(invalid("paired swept-flow compliance overflowed")); }
            let mut evaluations = 0;
            let mut evaluate = |velocity: f64| -> Result<Trial, AcousticRealizeError> {
                evaluations += 1;
                if evaluations > self.spec.limits.max_evaluations { return Err(invalid("paired aperture nonlinear work budget exhausted")); }
                let next_gap = gap + dt*velocity;
                let force = contact.response(next_gap, velocity).map_err(nonlinear)?.force;
                let base = [free[0]+mobility[0]*force, free[1]+mobility[1]*force];
                let moving = drive.body_flow_m3_s - b[0].pressure_area_m2*base[0] - b[1].pressure_area_m2*base[1];
                let source = (drive.upstream_pressure_pa - 2.0*drive.incoming_pressure_pa - z*moving)/denominator;
                let opening = f64::midpoint(gap, next_gap).max(0.0);
                // Eliminate the independent pressure/swept-motion response first.
                // This is the ORIGINAL continuous Bernoulli junction at its reduced
                // impedance, not an extra source or a fixed-opening approximation.
                let dp = match self.viscous_slit {
                    Some(law) => law.fixed_opening_pressure_drop(self.spec.width_m, opening,
                        self.spec.density_kg_m3, z/denominator, source)?,
                    None => super::moving::fixed_opening_pressure_drop(self.spec.width_m, opening,
                        self.spec.density_kg_m3, z/denominator, source)?,
                };
                let velocities = [base[0]-mobility[0]*b[0].pressure_area_m2*dp,
                    base[1]-mobility[1]*b[1].pressure_area_m2*dp];
                let residual = velocity - velocities[0] - velocities[1];
                let scale = velocity.abs() + velocities[0].abs() + velocities[1].abs()
                    + free[0].abs() + free[1].abs();
                if ![velocity, next_gap, dp, residual, scale].iter().all(|x| x.is_finite()) {
                    return Err(invalid("paired aperture trial left the finite physical domain"));
                }
                Ok(Trial { velocity: velocities, dp, residual, scale })
            };
            let resolved = |t: Trial| t.residual.abs() <= self.spec.limits.velocity_tolerance_m_s
                + self.spec.limits.relative_tolerance*t.scale;
            let initial = evaluate(free[0]+free[1])?;
            let trial = if resolved(initial) { initial } else {
                let mut span = (free[0].abs()+free[1].abs()).max(0.01);
                let (mut lo, mut hi) = (-span, span);
                let (mut left, mut right) = (evaluate(lo)?, evaluate(hi)?);
                while left.residual > 0.0 || right.residual < 0.0 {
                    span *= 2.0;
                    if !span.is_finite() { return Err(invalid("paired aperture root has no finite bracket")); }
                    if left.residual > 0.0 { lo = -span; left = evaluate(lo)?; }
                    if right.residual < 0.0 { hi = span; right = evaluate(hi)?; }
                }
                loop {
                    if resolved(left) { break left; }
                    if resolved(right) { break right; }
                    let mid = f64::midpoint(lo,hi);
                    if mid == lo || mid == hi { return Err(invalid("paired aperture unresolved gap momentum balance")); }
                    let next = evaluate(mid)?;
                    if resolved(next) { break next; }
                    if next.residual < 0.0 { lo=mid; left=next; } else { hi=mid; right=next; }
                }
            };
            trial
        };
        let after = [old[0].displacement_m+dt*trial.velocity[0], old[1].displacement_m+dt*trial.velocity[1]];
        let next_gap = self.spec.rest_opening_m + after[0] + after[1];
        let vm = trial.velocity[0]+trial.velocity[1];
        let (contact_forces, contact_normal, contact_loss, opening) = if let Some(closure) = &self.closure {
            let before=old.map(|s|s.displacement_m);
            let response=closure.response(before,after,trial.velocity)?;
            let middle=[f64::midpoint(before[0],after[0]),f64::midpoint(before[1],after[1])];
            (response.forces_n,response.normal_force_n,response.dissipated_power_w,
                closure.open_area_m2(middle)?/self.spec.width_m)
        } else {
            let response=contact.response(next_gap,vm).map_err(nonlinear)?;
            ([response.force;2],response.force,response.dissipated_power,
                f64::midpoint(gap,next_gap).max(0.0))
        };
        // Total force each blade holds over the step: its own external load, the
        // work-conjugate contact force and the uniform pressure drop.
        let held = [forces[0]+contact_forces[0]-b[0].pressure_area_m2*trial.dp,
            forces[1]+contact_forces[1]-b[1].pressure_area_m2*trial.dp];
        let next: [BladeState; 2] = std::array::from_fn(|i| BladeState {
            displacement_m: after[i],
            velocity_m_s: match &self.exact {
                Some(exact) => exact[i].end_velocity(old[i], held[i], arm_memory[i]),
                None => 2.0*trial.velocity[i]-old[i].velocity_m_s,
            },
        });
        self.validate_state(next)?;
        let jet = match self.viscous_slit {
            Some(law) => match &self.closure {
                Some(closure) => {
                    let middle = [f64::midpoint(old[0].displacement_m, after[0]),
                        f64::midpoint(old[1].displacement_m, after[1])];
                    closure.viscous_flow(law, middle, self.spec.density_kg_m3, trial.dp)?.flow_m3_s
                }
                None => law.evaluate(self.spec.width_m, opening, self.spec.density_kg_m3, trial.dp)?.flow_m3_s,
            },
            None => super::moving::volume_flow(self.spec.width_m,opening,trial.dp,self.spec.density_kg_m3),
        };
        let swept = -b[0].pressure_area_m2*trial.velocity[0]-b[1].pressure_area_m2*trial.velocity[1];
        let no_drop = drive.upstream_pressure_pa - drive.incoming_pressure_pa;
        let load_pressure = z*(jet+swept+drive.body_flow_m3_s);
        let outgoing = if no_drop.abs().max(trial.dp.abs())
            <= drive.incoming_pressure_pa.abs().max(load_pressure.abs()) {
            no_drop-trial.dp
        } else { drive.incoming_pressure_pa+load_pressure };
        let bore = outgoing+drive.incoming_pressure_pa;
        let wave = (outgoing-drive.incoming_pressure_pa)/z;
        let viscous = (0..2).map(|i| b[i].damping_n_s_m*trial.velocity[i]*trial.velocity[i]).sum::<f64>();
        let mechanical = dt*(forces[0]*trial.velocity[0]+forces[1]*trial.velocity[1]);
        let mut memory: [Option<MemoryTrial>; 2] = [None, None];
        let (mut memory_energy, mut memory_loss) = (0.0, 0.0);
        for i in 0..2 {
            if let Some(material) = &self.relaxation[i] {
                let result = match &self.exact {
                    Some(exact) => material.exact_trial(
                        exact[i].end_memory(old[i], held[i], arm_memory[i]),
                        next[i].displacement_m,
                        exact[i].arm_loss_j(old[i], held[i], arm_memory[i]))?,
                    None => material.preview(old[i].displacement_m, next[i].displacement_m, dt)?,
                };
                memory_energy += result.energy;
                memory_loss += result.loss;
                memory[i] = Some(result);
            }
        }
        let energy = self.energy(next)+memory_energy;
        self.validate_energy(energy)?;
        let aperture = ApertureFrame {
            step: self.samples+1, time_s: (self.samples+1) as f64*dt,
            state: ApertureState { opening_m: next_gap, opening_velocity_m_s: next[0].velocity_m_s+next[1].velocity_m_s },
            midpoint_opening_m: opening, outgoing_pressure_pa: outgoing, bore_pressure_pa: bore,
            jet_flow_m3_s: jet, swept_flow_m3_s: swept, bore_flow_m3_s: wave,
            flow_residual_m3_s: jet+swept+drive.body_flow_m3_s-wave,
            stored_energy_j: energy, storage_change_j: energy-self.stored_energy_j(),
            dissipated_energy_j: match &self.exact {
                // The exact blocks integrate c*v^2 along their own trajectories.
                Some(exact) => dt*(trial.dp*jet+contact_loss)+memory_loss
                    + exact[0].viscous_loss_j(old[0], held[0], arm_memory[0])
                    + exact[1].viscous_loss_j(old[1], held[1], arm_memory[1]),
                None => dt*(trial.dp*jet+viscous+contact_loss)+memory_loss,
            },
            pressure_work_j: dt*trial.dp*(wave-drive.body_flow_m3_s), mechanical_work_j: mechanical,
            relaxation_energy_j: memory_energy, relaxation_loss_j: memory_loss,
        };
        let mut momentum = [0.0; 2];
        for i in 0..2 {
            // Step-mean displacement: exact time average, or the midpoint value.
            let mean_x = match &self.exact {
                Some(exact) => exact[i].mean_displacement(old[i], held[i], dt, arm_memory[i]),
                None => f64::midpoint(old[i].displacement_m,next[i].displacement_m),
            };
            momentum[i] = b[i].mass_kg*(next[i].velocity_m_s-old[i].velocity_m_s)/dt
                + b[i].stiffness_n_m*mean_x
                + b[i].damping_n_s_m*trial.velocity[i]+b[i].pressure_area_m2*trial.dp-contact_forces[i]-forces[i];
            if let Some(material) = &self.relaxation[i] {
                momentum[i] += match &self.exact {
                    Some(exact) => exact[i].mean_memory_force_n(old[i], held[i], dt, arm_memory[i]),
                    None => material.restoring_force(old[i].displacement_m, next[i].displacement_m, dt)?,
                };
            }
        }
        let flow_scale = jet.abs()+swept.abs()+drive.body_flow_m3_s.abs()+wave.abs();
        let scale = self.stored_energy_j()+energy+aperture.pressure_work_j.abs()+mechanical.abs()
            + aperture.dissipated_energy_j.abs();
        let allowance = 128.0*self.spec.limits.relative_tolerance.max(f64::EPSILON);
        if ![outgoing,bore,jet,swept,wave,aperture.flow_residual_m3_s,
            aperture.balance_residual_j(),aperture.dissipated_energy_j,momentum[0],momentum[1],scale,flow_scale]
            .iter().all(|x| x.is_finite())
            || aperture.balance_residual_j().abs() > allowance*scale.max(f64::MIN_POSITIVE)
            || aperture.dissipated_energy_j < -allowance*scale
            || aperture.flow_residual_m3_s.abs() > allowance*flow_scale.max(f64::MIN_POSITIVE)
        { return Err(invalid("paired aperture refused an unresolved physical work or flow balance")); }
        Ok(PairedTrial { frame: PairedApertureFrame { aperture, blades: next, pressure_drop_pa: trial.dp,
            contact_force_n: contact_normal, blade_contact_forces_n: contact_forces, momentum_residual_n: momentum }, memory })
    }
}

impl PairedAperture {
    fn admit_linear_clock(&self) -> Result<(), AcousticRealizeError> {
        // Exact blocks do not warp; the midpoint band applies only to Midpoint.
        if self.exact.is_some() { return Ok(()); }
        let blades = self.spec.blades;
        let weights = blades.map(|b| b.pressure_area_m2 / fs_math::det::sqrt(b.mass_kg));
        let weight_sum = weights.iter().map(|a| a.abs()).sum::<f64>();
        let mut rate = 0.0_f64;
        for (i,b) in blades.iter().enumerate() {
            let loss = b.damping_n_s_m / b.mass_kg
                + self.spec.impedance_pa_s_m3 * weights[i].abs() * weight_sum;
            let mut local = fs_math::det::sqrt(b.stiffness_n_m / b.mass_kg) + loss;
            if let Some(memory) = self.relaxation[i].as_ref().filter(|m| !m.branches().is_empty()) {
                let instantaneous = b.stiffness_n_m + memory.branches().iter().map(|arm| arm.stiffness).sum::<f64>();
                local = 2.0 * fs_math::det::sqrt(instantaneous / b.mass_kg) + loss;
                for arm in memory.branches() { local = local.max(2.0 / arm.relaxation_time_s); }
            }
            rate = rate.max(local);
        }
        fs_phs::linear_exact::admit_midpoint_band(rate, self.spec.time_step_s)
            .map_err(AcousticRealizeError::TimeResolution)
    }
}

#[cfg(test)]
mod tests;

#[cfg(test)]
mod relaxation_tests;

#[cfg(test)]
mod spatial_tests;

#[cfg(test)]
mod exact_memory_tests;
