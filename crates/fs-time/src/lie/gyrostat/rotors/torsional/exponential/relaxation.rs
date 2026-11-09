//! Maxwell material memory on the original complete rotor/carrier trajectory.

use super::{
    Cache, ExponentialTorsionalRotorStep, ExponentialTorsionalRotorStepper, HeldTorque,
    MAX_EXPONENTIAL_TORSIONAL_ROTORS, Result, RotorBank, RotorState, dd, invalid, memory_strain,
    poll, prepare, scalar,
};
use fs_math::dd::Dd;

/// At most 32 retained material arms, in addition to at most 12 original rotors.
/// Every arm keeps its own state; no fitted common time or truncation is used.
pub const MAX_TORSIONAL_RELAXATION_BRANCHES: usize = 32;

/// One passive Maxwell shaft arm acting on a linear combination of relative angles.
///
/// `epsilon = sum projection[i] * (q_i - rest_i)`, with dimensionless projection;
/// `stress = k * (epsilon - z)` [N m], `z_dot = (epsilon - z)/tau` [rad/s].
/// The arm stores `k * (epsilon-z)^2/2` [J] and dissipates
/// `k * (epsilon-z)^2/tau` [W]. It is added to the equilibrium stiffness and
/// memoryless drag supplied separately. Several arms give a generalized Maxwell
/// law; no material coefficients or applicability band are inferred.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionalRelaxationBranch<const N: usize> {
    pub projection: [f64; N],
    pub stiffness_nm_per_rad: f64,
    pub relaxation_time_s: f64,
}

/// Prepared causal, frequency-dependent shaft loss on the complete native bank.
///
/// Each frozen-attitude stage is the existing full-state exponential with the
/// canonical `fs_phs::RelaxationBranch` storage and resistance. Their independently
/// integrated heat never comes from an energy difference. Two exact half-stages
/// surround the original gyrostat drift; the complete split trajectory is second
/// order and still needs refinement for tumbling or spatial-torque reorientation.
///
/// The cold dense preparation is bounded by the public rotor/arm caps. Warm
/// steps reuse all matrices and private scratch, allocate nothing, and return
/// original physical viscous strains for explicit checkpoint/resume. Rates,
/// geometry, rest angles and material parameters stay fixed. Thermal evolution,
/// fitted loss factors, nonlinear shaft strain and discrete adjoints are absent.
pub struct TorsionalRelaxationStepper<const N: usize, const B: usize> {
    native: ExponentialTorsionalRotorStepper<N>,
}

impl<const N: usize, const B: usize> core::fmt::Debug for TorsionalRelaxationStepper<N, B> {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        f.debug_struct("TorsionalRelaxationStepper")
            .field("native", &self.native)
            .field("material_arms", &B)
            .finish()
    }
}

/// Complete accepted mechanical and material state with separate storage/heat.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct TorsionalRelaxationStep<const N: usize, const B: usize> {
    /// Native mechanics; `internal_energy_change_j` is the arms' storage change,
    /// and `dissipated_energy_j` includes their loss plus the original drag once.
    pub motion: ExponentialTorsionalRotorStep<N>,
    /// Physical viscous strain z [rad], in the original supplied branch order.
    /// Pass these together with `motion.state_after` to the next call.
    pub viscous_angle_after_rad: [f64; B],
}

impl<const N: usize> RotorBank<N> {
    /// Prepare equilibrium springs, Kelvin drag and retained Maxwell shaft arms.
    ///
    /// Positive branch stiffness and time and finite projections are required.
    /// The spring parameters are the EQUILIBRIUM stiffness; adding the arms to
    /// an already unrelaxed stiffness would double-count their storage.
    /// `B=0` delegates exactly to the existing memoryless exponential owner.
    pub fn torsional_relaxation_stepper<const B: usize>(
        &self,
        equilibrium_stiffness_nm_per_rad: [f64; N],
        rest_angle_rad: [f64; N],
        damping_nms: [f64; N],
        branches: [TorsionalRelaxationBranch<N>; B],
        duration_s: f64,
    ) -> Result<TorsionalRelaxationStepper<N, B>> {
        self.torsional_relaxation_stepper_with_cancel(
            equilibrium_stiffness_nm_per_rad,
            rest_angle_rad,
            damping_nms,
            branches,
            duration_s,
            || false,
        )
    }

    /// Same cold preparation with cancellation at bounded native work boundaries.
    #[allow(clippy::too_many_arguments)]
    pub fn torsional_relaxation_stepper_with_cancel<const B: usize>(
        &self,
        stiffness: [f64; N],
        rest: [f64; N],
        damping: [f64; N],
        branches: [TorsionalRelaxationBranch<N>; B],
        duration: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<TorsionalRelaxationStepper<N, B>> {
        poll(&mut cancelled)?;
        if B == 0 {
            return Ok(TorsionalRelaxationStepper {
                native: self.torsional_exponential_stepper_with_cancel(
                    stiffness, rest, damping, duration, cancelled,
                )?,
            });
        }
        if B > MAX_TORSIONAL_RELAXATION_BRANCHES || N > MAX_EXPONENTIAL_TORSIONAL_ROTORS {
            return Err(invalid("bounded complete torsional material inventory").into());
        }
        scalar(duration, "torsional relaxation duration")?;
        if duration <= 0.0 || duration * 0.5 == 0.0 {
            return Err(invalid("positive resolved torsional relaxation duration").into());
        }
        for i in 0..N {
            scalar(stiffness[i], "equilibrium torsional stiffness")?;
            scalar(damping[i], "torsional memoryless damping")?;
            scalar(rest[i], "unwrapped torsional rest phase")?;
            if stiffness[i] < 0.0 || damping[i] < 0.0 {
                return Err(invalid("nonnegative equilibrium stiffness and damping").into());
            }
        }
        let cache = Cache::Coupled(Box::new(prepare(
            self, stiffness, damping, &branches, 0.5 * duration, &mut cancelled,
        )?));
        poll(&mut cancelled)?;
        Ok(TorsionalRelaxationStepper {
            native: ExponentialTorsionalRotorStepper {
                bank: *self,
                stiffness,
                rest,
                damping,
                duration,
                cache,
            },
        })
    }
}

impl<const N: usize, const B: usize> TorsionalRelaxationStepper<N, B> {
    pub const fn duration_s(&self) -> f64 {
        self.native.duration_s()
    }

    /// Fully relaxed material initialization at the actual original shaft angles.
    /// An instantaneous strain from an unstrained material instead starts z=0.
    /// This method is initialization only; never call it to advance a history.
    pub fn relaxed_memory(&self, state: &RotorState<N>) -> Result<[f64; B]> {
        self.native.bank.admit_state(state)?;
        let mut memory = [0.0; B];
        if let Cache::Coupled(cache) = &self.native.cache {
            for (i, branch) in cache.memory.iter().enumerate() {
                memory[i] = scalar(
                    memory_strain(branch, self.native.rest, state.relative_angle_rad, 0.0).to_f64(),
                    "relaxed original torsional strain",
                )?;
            }
        }
        Ok(memory)
    }

    /// The arms' retained physical energy only; equilibrium springs are separate.
    pub fn internal_energy_j(&self, state: &RotorState<N>, memory: &[f64; B]) -> Result<f64> {
        self.native.bank.admit_state(state)?;
        let mut energy = Dd::ZERO;
        if let Cache::Coupled(cache) = &self.native.cache {
            for (i, branch) in cache.memory.iter().enumerate() {
                scalar(memory[i], "finite torsional material memory")?;
                let strain = memory_strain(
                    branch, self.native.rest, state.relative_angle_rad, memory[i],
                );
                energy = energy + dd(0.5) * dd(branch.stiffness) * strain * strain;
            }
        }
        Ok(scalar(energy.to_f64(), "torsional material stored energy")?)
    }

    pub fn step(
        &mut self,
        state: &RotorState<N>,
        memory: &[f64; B],
        drive: [f64; N],
        external: HeldTorque,
    ) -> Result<TorsionalRelaxationStep<N, B>> {
        self.step_with_cancel(state, memory, drive, external, || false)
    }

    /// Atomic mechanical/material advance. Cancellation or refusal publishes
    /// neither; scratch may change, and retry uses the original supplied pair.
    pub fn step_with_cancel(
        &mut self,
        state: &RotorState<N>,
        memory: &[f64; B],
        drive: [f64; N],
        external: HeldTorque,
        cancelled: impl FnMut() -> bool,
    ) -> Result<TorsionalRelaxationStep<N, B>> {
        let mut midpoint = [0.0; B];
        let mut endpoint = [0.0; B];
        let motion = self.native.step_with_memory(
            state, memory, &mut midpoint, &mut endpoint, drive, external, cancelled,
        )?;
        Ok(TorsionalRelaxationStep {
            motion,
            viscous_angle_after_rad: endpoint,
        })
    }
}
