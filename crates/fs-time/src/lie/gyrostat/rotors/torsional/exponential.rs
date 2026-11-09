//! Exact linear shaft stages, with the original canonical rotor/carrier storage.
//!
//! At fixed attitude, x=(q-q_rest,p,L) has Hamiltonian
//! H=sum(K*q²+p²/J)/2+(L-A*p)^T I^-1 (L-A*p)/2. Canonical q/p pairs,
//! resistance on p, and physical shaft/external torque ports give the original
//! coupled equations. The existing port-Hamiltonian exact-ZOH owner propagates
//! ALL these coordinates and independently integrates positive physical heat.

// These indices are the explicit original coordinate and port ordering.
#![allow(clippy::needless_range_loop, clippy::cast_precision_loss)]

use super::super::{
    ExponentialViscousRotorStepper, HeldTorque, LieStepError, RotorBank, RotorImpulse, RotorState,
    components, dd, dot_dd, increment, invalid, product, scalar, vector,
};
use crate::lie::compensated;
use fs_ga::Vec3;
use fs_math::dd::Dd;
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_phs::{PhsError, PortHamiltonian, QuadraticStorage, ResistanceWitness};

// torsional.rs loads this file through #[path], so a bare `mod` would resolve beside it.
#[path = "exponential/relaxation.rs"]
mod relaxation;
pub use relaxation::{
    MAX_TORSIONAL_RELAXATION_BRANCHES, TorsionalRelaxationBranch, TorsionalRelaxationStep,
    TorsionalRelaxationStepper,
};

/// Hard bound on the complete dense spring-bearing preparation: 27 original
/// coordinates and at most 110 augmented exponential rows. No shaft is dropped.
/// The all-zero-stiffness path retains the viscous owner's 64-rotor bound.
pub const MAX_EXPONENTIAL_TORSIONAL_ROTORS: usize = 12;

/// Native mechanics, bounded linear propagation or cooperative cancellation.
#[derive(Debug, Clone, PartialEq)]
pub enum TorsionalExponentialError {
    /// Original Lie-group, impulse, recoil or local-work admission failed.
    Native(LieStepError),
    /// The shared full-state exponential/positive-heat owner refused.
    Linear(PhsError),
    /// The caller cancelled before a complete result was returned.
    Cancelled,
}

impl From<LieStepError> for TorsionalExponentialError {
    fn from(value: LieStepError) -> Self {
        Self::Native(value)
    }
}
impl From<PhsError> for TorsionalExponentialError {
    fn from(value: PhsError) -> Self {
        Self::Linear(value)
    }
}
impl core::fmt::Display for TorsionalExponentialError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Native(error) => error.fmt(f),
            Self::Linear(error) => error.fmt(f),
            Self::Cancelled => f.write_str("torsional exponential step cancelled"),
        }
    }
}
impl std::error::Error for TorsionalExponentialError {}

type Result<T> = core::result::Result<T, TorsionalExponentialError>;

struct Coupled {
    flow: PreparedExactZoh<'static>,
    initial: Vec<f64>,
    input: Vec<f64>,
    endpoint: Vec<f64>,
    output: Vec<f64>,
    drive_input: Vec<f64>,
    external_input: Vec<f64>,
    memory: Vec<MemoryBranch>,
}

struct MemoryBranch {
    projection: Vec<f64>,
    stiffness: f64,
    root_stiffness: f64,
}

enum Cache<const N: usize> {
    Viscous(Box<ExponentialViscousRotorStepper<N>>),
    Coupled(Box<Coupled>),
}

/// Prepared full coupled spring/drag dynamics for one fixed positive clock.
///
/// The frozen-attitude shaft stages use the native exact-ZOH/positive-heat
/// owner, retaining every elastic, damping, inertial and carrier-recoil entry.
/// Zero stiffness/damping directions remain in the complete original system.
/// No fitted modal damping or independent-shaft approximation is substituted.
/// Preparation is bounded dense O(N³); cached shaft stages are O(N²) and
/// allocate nothing, factor nothing and evaluate no exponentials.
///
/// Two exact shaft half-stages surround the native gyrostat drift, so the
/// complete trajectory remains second order. Refine for attitude/precession
/// and spatially held torque reorientation even when shaft relaxation is stiff.
/// The native dense owner may refuse unrepresentable/ill-conditioned dynamics;
/// arbitrary stiffness, interval certification and discrete adjoints are not
/// implied. Original coefficients and rest phases remain fixed.
pub struct ExponentialTorsionalRotorStepper<const N: usize> {
    bank: RotorBank<N>,
    stiffness: [f64; N],
    rest: [f64; N],
    damping: [f64; N],
    duration: f64,
    cache: Cache<N>,
}

impl<const N: usize> core::fmt::Debug for ExponentialTorsionalRotorStepper<N> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("ExponentialTorsionalRotorStepper")
            .field("bank", &self.bank)
            .field("stiffness", &self.stiffness)
            .field("rest", &self.rest)
            .field("damping", &self.damping)
            .field("duration", &self.duration)
            .finish_non_exhaustive()
    }
}

/// Accepted mechanics and separate original physical work/storage/heat accounts.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct ExponentialTorsionalRotorStep<const N: usize> {
    /// Complete canonical state for the next native rotor step.
    pub state_after: RotorState<N>,
    /// Actual accepted changes of original axial momenta [N m s].
    pub axial_momentum_change_nms: [f64; N],
    /// Actual external angular impulse in inertial coordinates [N m s].
    pub external_angular_impulse_world: Vec3,
    /// Signed held drive torque times its continuous integrated relative angle.
    pub drive_work_j: [f64; N],
    /// Independently integrated original viscous and optional Maxwell losses,
    /// never an energy gap. No per-shaft allocation of aggregate heat is supplied.
    pub dissipated_energy_j: f64,
    /// Polarized spring storage changes of the actual accepted unwrapped phases.
    pub elastic_energy_change_j: [f64; N],
    /// Maxwell internal-storage change [J]. Zero for the memoryless entry.
    pub internal_energy_change_j: f64,
    /// Drive work minus spring/internal storage change minus physical heat.
    pub net_shaft_work_j: f64,
    /// Signed continuous external torque work from the same shaft trajectories.
    pub external_work_j: f64,
    /// Polarized total kinetic change of the accepted rotor/carrier state [J].
    pub kinetic_energy_change_j: f64,
    /// The native free gyrostat's numerical defect, kept distinct from heat.
    pub drift_energy_defect_j: f64,
    /// Delta kinetic + Delta elastic + Delta internal - drive - external
    /// + heat - native drift defect.
    pub work_energy_residual_j: f64,
    /// Accepted world momentum change minus the actual external impulse.
    pub angular_momentum_residual_world: Vec3,
}

fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(TorsionalExponentialError::Cancelled)
    } else {
        Ok(())
    }
}

impl<const N: usize> RotorBank<N> {
    /// Prepare full coupled, stiff spring/drag shaft subflows at one fixed clock.
    ///
    /// Springs and viscous coefficients are finite/nonnegative; rest phases are
    /// finite. Spring-bearing systems admit at most 12 original rotors. All-zero
    /// stiffness delegates exactly to `viscous_exponential_stepper`, including
    /// its larger inventory cap and zero-damping native held-torque behavior.
    pub fn torsional_exponential_stepper(
        &self,
        stiffness_nm_per_rad: [f64; N],
        rest_angle_rad: [f64; N],
        damping_nms: [f64; N],
        duration_s: f64,
    ) -> Result<ExponentialTorsionalRotorStepper<N>> {
        self.torsional_exponential_stepper_with_cancel(
            stiffness_nm_per_rad,
            rest_angle_rad,
            damping_nms,
            duration_s,
            || false,
        )
    }

    /// Same preparation with polling before/after bounded native dense work.
    /// Cancellation returns no partially prepared owner. Native exponential and
    /// factor calls do not poll internally; the explicit inventory cap bounds
    /// their work and storage.
    pub fn torsional_exponential_stepper_with_cancel(
        &self,
        stiffness_nm_per_rad: [f64; N],
        rest_angle_rad: [f64; N],
        damping_nms: [f64; N],
        duration_s: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExponentialTorsionalRotorStepper<N>> {
        poll(&mut cancelled)?;
        scalar(duration_s, "exponential torsional duration")?;
        if duration_s <= 0.0 || 0.5 * duration_s == 0.0 {
            return Err(invalid("positive resolved exponential torsional duration").into());
        }
        let mut has_spring = false;
        for i in 0..N {
            let k = scalar(stiffness_nm_per_rad[i], "exponential torsional stiffness")?;
            let d = scalar(damping_nms[i], "exponential torsional damping")?;
            scalar(rest_angle_rad[i], "exponential spring rest phase")?;
            if k < 0.0 || d < 0.0 {
                return Err(invalid("nonnegative exponential spring and damping").into());
            }
            has_spring |= k > 0.0;
        }
        let cache = if has_spring {
            if N > MAX_EXPONENTIAL_TORSIONAL_ROTORS {
                return Err(invalid("bounded complete exponential torsional inventory").into());
            }
            Cache::Coupled(Box::new(prepare(
                self,
                stiffness_nm_per_rad,
                damping_nms,
                &[],
                0.5 * duration_s,
                &mut cancelled,
            )?))
        } else {
            Cache::Viscous(Box::new(
                self.viscous_exponential_stepper(damping_nms, duration_s)?,
            ))
        };
        poll(&mut cancelled)?;
        Ok(ExponentialTorsionalRotorStepper {
            bank: *self,
            stiffness: stiffness_nm_per_rad,
            rest: rest_angle_rad,
            damping: damping_nms,
            duration: duration_s,
            cache,
        })
    }
}

fn prepare<const N: usize>(
    bank: &RotorBank<N>,
    stiffness: [f64; N],
    damping: [f64; N],
    branches: &[TorsionalRelaxationBranch<N>],
    half_s: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Coupled> {
    let n = 2 * N + 3;
    // Physical ports: N shaft drives, three external torques. N additional
    // zero-input observers integrate spring efforts for the original impulse law.
    let m = 2 * N + 3;
    let mut q = vec![0.0; n * n];
    let mut j = vec![0.0; n * n];
    let mut r = vec![0.0; n * n];
    let mut root = vec![0.0; n * n];
    let mut g = vec![0.0; n * m];
    let inverse = [
        bank.inertia.velocity(Vec3::new(1.0, 0.0, 0.0))?,
        bank.inertia.velocity(Vec3::new(0.0, 1.0, 0.0))?,
        bank.inertia.velocity(Vec3::new(0.0, 0.0, 1.0))?,
    ];
    for row in 0..3 {
        for col in row..3 {
            let entry = scalar(
                0.5 * components(inverse[col])[row] + 0.5 * components(inverse[row])[col],
                "exponential carrier inverse inertia",
            )?;
            q[(2 * N + row) * n + 2 * N + col] = entry;
            q[(2 * N + col) * n + 2 * N + row] = entry;
        }
        g[(2 * N + row) * m + N + row] = 1.0;
    }
    for a in 0..N {
        poll(cancelled)?;
        q[a * n + a] = stiffness[a];
        j[a * n + N + a] = 1.0;
        j[(N + a) * n + a] = -1.0;
        r[(N + a) * n + N + a] = damping[a];
        root[(N + a) * n + N + a] = fs_math::det::sqrt(damping[a]);
        g[(N + a) * m + a] = 1.0;
        g[a * m + N + 3 + a] = 1.0;
        let axis = bank.rotors[a].axis_body;
        let mobility = bank.inertia.velocity(axis)?;
        for b in a..N {
            let entry = scalar(
                (dot_dd(bank.rotors[b].axis_body, mobility)
                    + if a == b {
                        dd(1.0 / bank.rotors[a].axial_inertia_kg_m2)
                    } else {
                        Dd::ZERO
                    })
                .to_f64(),
                "exponential shaft mobility",
            )?;
            q[(N + a) * n + N + b] = entry;
            q[(N + b) * n + N + a] = entry;
        }
        for row in 0..3 {
            let entry = -components(mobility)[row];
            q[(N + a) * n + 2 * N + row] = entry;
            q[(2 * N + row) * n + N + a] = entry;
        }
    }
    let witness = ResistanceWitness::try_from_source_rows(n, r.clone(), root)?;
    let storage = QuadraticStorage::new(q, n)?;
    let system =
        PortHamiltonian::new(n, m, j, r, g, Box::new(storage))?.with_resistance_witness(witness)?;
    let mut memory = Vec::with_capacity(branches.len());
    let mut laws = Vec::with_capacity(branches.len());
    for branch in branches {
        poll(cancelled)?;
        let mut projection = vec![0.0; n];
        projection[..N].copy_from_slice(&branch.projection);
        laws.push(fs_phs::RelaxationBranch {
            projection,
            stiffness: branch.stiffness_nm_per_rad,
            relaxation_time_s: branch.relaxation_time_s,
        });
        memory.push(MemoryBranch {
            projection: branch.projection.to_vec(),
            stiffness: branch.stiffness_nm_per_rad,
            root_stiffness: fs_math::det::sqrt(branch.stiffness_nm_per_rad),
        });
    }
    let system = system.with_relaxation_branches(laws)?;
    let n = system.state_dim();
    poll(cancelled)?;
    let flow = PreparedExactZoh::new_owned(
        system,
        half_s,
        LinearSplitBudget {
            max_modes: n,
            max_omega_dt: 1.0,
        },
    )?;
    poll(cancelled)?;
    Ok(Coupled {
        flow,
        initial: vec![0.0; n],
        input: vec![0.0; m],
        endpoint: vec![0.0; n],
        output: vec![0.0; m],
        drive_input: vec![0.0; m],
        external_input: vec![0.0; m],
        memory,
    })
}

impl<const N: usize> ExponentialTorsionalRotorStepper<N> {
    /// Prepared positive full interval [s].
    #[must_use]
    pub const fn duration_s(&self) -> f64 {
        self.duration
    }

    /// Original spring potential of the actual unwrapped relative phases [J].
    pub fn potential_energy(&self, state: &RotorState<N>) -> Result<f64> {
        self.bank.admit_state(state)?;
        let mut value = Dd::ZERO;
        for i in 0..N {
            if self.stiffness[i] == 0.0 {
                continue;
            }
            let strain = dd(state.relative_angle_rad[i]) - dd(self.rest[i]);
            value = value + dd(0.5) * dd(self.stiffness[i]) * strain * strain;
        }
        Ok(scalar(value.to_f64(), "exponential torsional potential")?)
    }

    /// Advance held drives and body- or space-held external torque, preserving
    /// caller state on any refusal. Only reusable private scratch is mutable.
    pub fn step(
        &mut self,
        state: &RotorState<N>,
        drive_torque_nm: [f64; N],
        external: HeldTorque,
    ) -> Result<ExponentialTorsionalRotorStep<N>> {
        self.step_with_cancel(state, drive_torque_nm, external, || false)
    }

    /// Poll cancellation before/after each bounded shaft stage and native drift.
    /// A cancelled call can be retried with the same original state and cache.
    pub fn step_with_cancel(
        &mut self,
        state: &RotorState<N>,
        drive: [f64; N],
        external: HeldTorque,
        cancelled: impl FnMut() -> bool,
    ) -> Result<ExponentialTorsionalRotorStep<N>> {
        self.step_with_memory(state, &[], &mut [], &mut [], drive, external, cancelled)
    }

    #[allow(clippy::too_many_arguments)]
    fn step_with_memory(
        &mut self,
        state: &RotorState<N>,
        memory_before: &[f64],
        memory_midpoint: &mut [f64],
        memory_after: &mut [f64],
        drive: [f64; N],
        external: HeldTorque,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<ExponentialTorsionalRotorStep<N>> {
        poll(&mut cancelled)?;
        if let Cache::Viscous(viscous) = &self.cache {
            let step = viscous.step(state, drive, external)?;
            poll(&mut cancelled)?;
            return Ok(ExponentialTorsionalRotorStep {
                state_after: step.state_after,
                axial_momentum_change_nms: step.axial_momentum_change_nms,
                external_angular_impulse_world: step.external_angular_impulse_world,
                drive_work_j: step.drive_work_j,
                dissipated_energy_j: step.dissipated_energy_j,
                elastic_energy_change_j: [0.0; N],
                internal_energy_change_j: 0.0,
                net_shaft_work_j: step.net_shaft_work_j,
                external_work_j: step.external_work_j,
                kinetic_energy_change_j: step.kinetic_energy_change_j,
                drift_energy_defect_j: step.drift_energy_defect_j,
                work_energy_residual_j: step.work_energy_residual_j,
                angular_momentum_residual_world: step.angular_momentum_residual_world,
            });
        }
        let Cache::Coupled(cache) = &mut self.cache else {
            unreachable!()
        };
        let first = stage(
            &self.bank,
            self.stiffness,
            self.rest,
            self.damping,
            cache,
            state,
            drive,
            external,
            memory_before,
            memory_midpoint,
        )?;
        poll(&mut cancelled)?;
        let mut drifted = first.impulse.state_after;
        drifted.relative_angle_rad = first.phase;
        (drifted.rotation, drifted.total_momentum_body) = self.bank.inertia.momentum_step(
            drifted.rotation,
            drifted.total_momentum_body,
            self.bank.rotor_momentum(&drifted.axial_momentum_nms)?,
            self.duration,
        )?;
        poll(&mut cancelled)?;
        let second = stage(
            &self.bank,
            self.stiffness,
            self.rest,
            self.damping,
            cache,
            &drifted,
            drive,
            external,
            memory_midpoint,
            memory_after,
        )?;
        poll(&mut cancelled)?;
        let native = self.bank.finish_step(
            state,
            &first.impulse,
            &drifted,
            &second.impulse,
            self.duration,
            Some(second.phase),
        )?;
        let elastic = potential_changes(
            self.stiffness,
            self.rest,
            state.relative_angle_rad,
            second.phase,
        )?;
        let mut drive_work = [0.0; N];
        let heat = dd(first.heat) + dd(second.heat);
        let external_work = dd(first.external_work) + dd(second.external_work);
        let internal_change = memory_energy_change(
            &cache.memory,
            self.rest,
            state.relative_angle_rad,
            memory_before,
            second.phase,
            memory_after,
        )?;
        let mut net = -heat;
        if !cache.memory.is_empty() {
            net = net - dd(internal_change);
        }
        for i in 0..N {
            drive_work[i] = scalar(
                (dd(first.drive_work[i]) + dd(second.drive_work[i])).to_f64(),
                "exponential torsional drive work",
            )?;
            net = net + dd(drive_work[i]) - dd(elastic[i]);
        }
        let residual = dd(native.kinetic_energy_change_j)
            - dd(native.drift_energy_defect_j)
            - external_work
            - net;
        Ok(ExponentialTorsionalRotorStep {
            state_after: native.state_after,
            axial_momentum_change_nms: native.axial_momentum_change_nms,
            external_angular_impulse_world: native.external_angular_impulse_world,
            drive_work_j: drive_work,
            dissipated_energy_j: scalar(heat.to_f64(), "exponential torsional heat")?,
            elastic_energy_change_j: elastic,
            internal_energy_change_j: internal_change,
            net_shaft_work_j: scalar(net.to_f64(), "exponential torsional net work")?,
            external_work_j: scalar(external_work.to_f64(), "exponential external torque work")?,
            kinetic_energy_change_j: native.kinetic_energy_change_j,
            drift_energy_defect_j: native.drift_energy_defect_j,
            work_energy_residual_j: scalar(residual.to_f64(), "exponential torsional balance")?,
            angular_momentum_residual_world: native.angular_momentum_residual_world,
        })
    }
}

fn potential_changes<const N: usize>(
    stiffness: [f64; N],
    rest: [f64; N],
    before: [f64; N],
    after: [f64; N],
) -> Result<[f64; N]> {
    let mut change = [0.0; N];
    for i in 0..N {
        if stiffness[i] == 0.0 {
            continue;
        }
        let mean = (dd(before[i]) - dd(rest[i])) * dd(0.5) + (dd(after[i]) - dd(rest[i])) * dd(0.5);
        change[i] = scalar(
            (dd(stiffness[i]) * mean * (dd(after[i]) - dd(before[i]))).to_f64(),
            "exponential actual spring storage increment",
        )?;
    }
    Ok(change)
}

struct Stage<const N: usize> {
    impulse: RotorImpulse<N>,
    phase: [f64; N],
    drive_work: [f64; N],
    external_work: f64,
    heat: f64,
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn stage<const N: usize>(
    bank: &RotorBank<N>,
    stiffness: [f64; N],
    rest: [f64; N],
    damping: [f64; N],
    cache: &mut Coupled,
    state: &RotorState<N>,
    drive: [f64; N],
    external: HeldTorque,
    memory_before: &[f64],
    memory_after: &mut [f64],
) -> Result<Stage<N>> {
    bank.admit_state(state)?;
    if memory_before.len() != cache.memory.len() || memory_after.len() != cache.memory.len() {
        return Err(invalid("complete torsional material memory state").into());
    }
    let h = cache.flow.duration();
    let torque = match external {
        HeldTorque::Body(value) => vector(value, "exponential external body torque")?,
        HeldTorque::Space(value) => state
            .rotation
            .inverse()
            .rotate(value)
            .map_err(LieStepError::from)?,
    };
    let requested = Vec3::new(
        product(torque.x, h, "torsional external impulse")?,
        product(torque.y, h, "torsional external impulse")?,
        product(torque.z, h, "torsional external impulse")?,
    );
    let next_total = crate::lie::torque_kick(state.total_momentum_body, requested, 1.0)?;
    let actual_external = next_total - state.total_momentum_body;
    if !cache.memory.is_empty()
        && components(requested).iter().all(|value| *value == 0.0)
        && components(actual_external).iter().all(|value| *value == 0.0)
    {
        // Exact physical equilibrium of this frozen-attitude shaft stage.
        // A dense exponential can perturb a nonbinary nullspace by one ULP;
        // do not manufacture material loss or unresolved impulses there.
        // Both DD words must vanish: there is no tolerance or clipping.
        let mut equilibrium = bank.relative_rates(state)?.iter().all(|rate| *rate == 0.0);
        for i in 0..N {
            scalar(drive[i], "exponential shaft drive")?;
            let effort = dd(stiffness[i])
                * (dd(state.relative_angle_rad[i]) - dd(rest[i]))
                - dd(drive[i]);
            equilibrium &= effort == Dd::ZERO;
        }
        for (arm, branch) in cache.memory.iter().enumerate() {
            scalar(memory_before[arm], "finite torsional viscous strain")?;
            equilibrium &= memory_strain(
                branch,
                rest,
                state.relative_angle_rad,
                memory_before[arm],
            ) == Dd::ZERO;
        }
        if equilibrium {
            let impulse = bank.apply_impulse(state, [0.0; N], Vec3::default())?;
            memory_after.copy_from_slice(memory_before);
            return Ok(Stage {
                impulse,
                phase: state.relative_angle_rad,
                drive_work: [0.0; N],
                external_work: 0.0,
                heat: 0.0,
            });
        }
    }
    cache.input.fill(0.0);
    cache.drive_input.fill(0.0);
    cache.external_input.fill(0.0);
    for i in 0..N {
        scalar(drive[i], "exponential shaft drive")?;
        // Translate each spring by its held-force equilibrium. This is the
        // SAME affine trajectory, with physical drive work observed separately.
        // It preserves zero/small oscillations about a large preload without
        // cancelling two prepared full-state/source maps. Free angles have no
        // spring origin and retain their real drive as a source.
        cache.initial[i] = if stiffness[i] == 0.0 {
            if cache
                .memory
                .iter()
                .any(|branch| branch.projection[i] != 0.0)
            {
                scalar(
                    (dd(state.relative_angle_rad[i]) - dd(rest[i])).to_f64(),
                    "original relaxing shaft strain",
                )?
            } else {
                0.0
            }
        } else {
            scalar(
                compensated::divide(
                    dd(stiffness[i]) * (dd(state.relative_angle_rad[i]) - dd(rest[i]))
                        - dd(drive[i]),
                    dd(stiffness[i]),
                )?
                .to_f64(),
                "exponential force-relative initial spring strain",
            )?
        };
        cache.initial[N + i] = state.axial_momentum_nms[i];
        cache.input[i] = if stiffness[i] == 0.0 { drive[i] } else { 0.0 };
        cache.drive_input[i] = drive[i];
    }
    for i in 0..3 {
        cache.initial[2 * N + i] = components(state.total_momentum_body)[i];
        let source = scalar(
            compensated::divide(dd(components(actual_external)[i]), dd(h))?.to_f64(),
            "represented continuous external torque",
        )?;
        cache.input[N + i] = source;
        cache.external_input[N + i] = source;
    }
    for (arm, branch) in cache.memory.iter().enumerate() {
        scalar(memory_before[arm], "finite torsional viscous strain")?;
        // The spring equilibrium translation above must translate EVERY
        // material memory coordinate too. Its physical strain and stress stay
        // unchanged; a held drive must never reset a relaxation history.
        let shifted = dd(memory_before[arm]) - memory_translation(branch, stiffness, drive)?;
        cache.initial[2 * N + 3 + arm] = scalar(
            (dd(branch.root_stiffness) * shifted).to_f64(),
            "normalized torsional memory strain",
        )?;
    }
    let (record, observed_work) = cache.flow.step_with_work_inputs_into(
        &cache.initial,
        &cache.input,
        [&cache.drive_input, &cache.external_input],
        &mut cache.endpoint,
        &mut cache.output,
    )?;
    let impulses = core::array::from_fn(|i| {
        (dd(cache.endpoint[N + i]) - dd(state.axial_momentum_nms[i])).to_f64()
    });
    let impulse = bank.apply_impulse(state, impulses, requested)?;
    let mut phase = state.relative_angle_rad;
    let mut drive_work = [0.0; N];
    let allowance = 32768.0 * f64::EPSILON * (2 * N + 4) as f64;
    let mut drive_sum = Dd::ZERO;
    for (arm, branch) in cache.memory.iter().enumerate() {
        let index = 2 * N + 3 + arm;
        // Advance the PHYSICAL strain from the normalized increment. Rebuilding
        // an absolute strain by multiply/divide can move a stationary history
        // by one ULP, particularly for a nonbinary sqrt(stiffness).
        let normalized_change = dd(cache.endpoint[index]) - dd(cache.initial[index]);
        let change = compensated::divide(normalized_change, dd(branch.root_stiffness))?;
        let before = dd(memory_before[arm]);
        memory_after[arm] = scalar(
            (before + change).to_f64(),
            "accepted torsional viscous strain",
        )?;
        let actual_change = dd(memory_after[arm]) - before;
        // This is a bound on public f64 storage rounding in radians, not an
        // arbitrary relative-to-increment tolerance. Loss/work checks below
        // still use the actual accepted physical history. Erased nonzero
        // history increments refuse instead of silently freezing the material.
        let scale = memory_before[arm]
            .abs()
            .max(memory_after[arm].abs())
            .max(change.to_f64().abs());
        if (normalized_change.to_f64() != 0.0 && actual_change.to_f64() == 0.0)
            || (actual_change - change).to_f64().abs() > 8.0 * f64::EPSILON * scale
        {
            return Err(invalid("resolved torsional material memory increment").into());
        }
    }
    for i in 0..N {
        let integrated = product(h, cache.output[i], "exponential torsional phase increment")?;
        phase[i] = increment(
            phase[i],
            integrated,
            "exponential unwrapped torsional phase",
        )?;
        let actual = dd(phase[i]) - dd(state.relative_angle_rad[i]);
        if scalar(
            (actual - dd(integrated)).to_f64(),
            "torsional phase reconstruction",
        )?
        .abs()
            > allowance * actual.to_f64().abs().max(integrated.abs())
        {
            return Err(invalid("resolved exponential torsional phase").into());
        }
        let authored = dd(h) * dd(drive[i]);
        let drag = dd(damping[i]) * dd(integrated);
        let spring = dd(h)
            * (dd(cache.output[N + 3 + i])
                + if stiffness[i] == 0.0 {
                    Dd::ZERO
                } else {
                    dd(drive[i])
                });
        let residual = dd(impulse.axial_momentum_change_nms[i]) - authored + drag + spring;
        let scale = impulse.axial_momentum_change_nms[i]
            .abs()
            .max(authored.to_f64().abs())
            .max(drag.to_f64().abs())
            .max(spring.to_f64().abs());
        if scalar(residual.to_f64(), "original integrated spring/drag law")?.abs()
            > allowance * scale
        {
            return Err(invalid("resolved original exponential shaft effort").into());
        }
        drive_work[i] = scalar(
            (dd(drive[i]) * dd(integrated)).to_f64(),
            "torsional drive work",
        )?;
        drive_sum = drive_sum + dd(drive_work[i]);
    }
    let elastic = potential_changes(stiffness, rest, state.relative_angle_rad, phase)?;
    let internal_change = memory_energy_change(
        &cache.memory,
        rest,
        state.relative_angle_rad,
        memory_before,
        phase,
        memory_after,
    )?;
    let mut balance = dd(impulse.kinetic_energy_change_j) - dd(observed_work[1]) - drive_sum
        + dd(record.dissipated);
    if !cache.memory.is_empty() {
        balance = balance + dd(internal_change);
    }
    let mut scale = impulse
        .kinetic_energy_change_j
        .abs()
        .max(observed_work[1].abs())
        .max(observed_work[0].abs())
        .max(record.dissipated.abs())
        .max(internal_change.abs());
    // Include actual opposing carrier/rotor transfers, never unchanged energy.
    let carrier_delta =
        bank.carrier_momentum(&impulse.state_after)? - bank.carrier_momentum(state)?;
    let mean_omega = bank.angular_velocity(state)?.scale(0.5)
        + bank.angular_velocity(&impulse.state_after)?.scale(0.5);
    for i in 0..3 {
        scale = scale.max(
            scalar(
                (dd(components(carrier_delta)[i]) * dd(components(mean_omega)[i])).to_f64(),
                "exponential carrier transfer work",
            )?
            .abs(),
        );
    }
    for i in 0..N {
        let mean_spin = (0.5 * state.axial_momentum_nms[i]
            + 0.5 * impulse.state_after.axial_momentum_nms[i])
            / bank.rotors[i].axial_inertia_kg_m2;
        scale = scale
            .max(
                scalar(
                    (dd(mean_spin) * dd(impulse.axial_momentum_change_nms[i])).to_f64(),
                    "exponential rotor transfer work",
                )?
                .abs(),
            )
            .max(elastic[i].abs())
            .max(drive_work[i].abs());
        balance = balance + dd(elastic[i]);
    }
    if record.dissipated < 0.0
        || scalar(balance.to_f64(), "exponential local torsional work balance")?.abs()
            > allowance * scale
        || scalar(
            (drive_sum - dd(observed_work[0])).to_f64(),
            "continuous drive work reconstruction",
        )?
        .abs()
            > allowance * scale
    {
        return Err(invalid("resolved exponential physical torsional work").into());
    }
    Ok(Stage {
        impulse,
        phase,
        drive_work,
        external_work: observed_work[1],
        heat: record.dissipated,
    })
}

fn memory_translation<const N: usize>(
    branch: &MemoryBranch,
    stiffness: [f64; N],
    drive: [f64; N],
) -> Result<Dd> {
    let mut translation = Dd::ZERO;
    for i in 0..N {
        if stiffness[i] != 0.0 && branch.projection[i] != 0.0 {
            translation = translation
                + dd(branch.projection[i]) * compensated::divide(dd(drive[i]), dd(stiffness[i]))?;
        }
    }
    Ok(translation)
}

fn memory_strain<const N: usize>(
    branch: &MemoryBranch,
    rest: [f64; N],
    phase: [f64; N],
    memory: f64,
) -> Dd {
    let mut strain = -dd(memory);
    for i in 0..N {
        strain = strain + dd(branch.projection[i]) * (dd(phase[i]) - dd(rest[i]));
    }
    strain
}

fn memory_energy_change<const N: usize>(
    branches: &[MemoryBranch],
    rest: [f64; N],
    phase_before: [f64; N],
    before: &[f64],
    phase_after: [f64; N],
    after: &[f64],
) -> Result<f64> {
    let mut change = Dd::ZERO;
    for (i, branch) in branches.iter().enumerate() {
        let a = memory_strain(branch, rest, phase_before, before[i]);
        let b = memory_strain(branch, rest, phase_after, after[i]);
        change = change + dd(0.5) * dd(branch.stiffness) * (a + b) * (b - a);
    }
    Ok(scalar(
        change.to_f64(),
        "original torsional internal storage change",
    )?)
}
