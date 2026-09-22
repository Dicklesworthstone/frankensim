//! Implicit driven dry friction on the existing two-way modal network.
//!
//! The attachment coordinate is x=B^T q (left minus optional right). A held
//! traction F acts through B; the prescribed surface advances by u*dt. The
//! SAME displacement increment defines slip and work: d=u*dt-B^T delta_q.
//! Sticking solves d=0 inside the static cone. Sliding solves F=mu(|d|/dt) N
//! with the sign of d, using fs-tribo's kinetic law, not a velocity regularizer.
//! The network supplies its complete bilateral held-force response, so every
//! connected body reacts back during this solve. Its original time stepper,
//! component/connection limits and state transaction remain in charge.
//!
//! This is step-averaged sticking, not an exact endpoint velocity constraint.
//! One tangential coordinate and a prescribed nonnegative normal load are
//! supported. No normal collision, finite contact patch, thermal/wear evolution,
//! unique Stribeck branch, experimental validation or hard-real-time claim is
//! implied. Input law coefficients retain caller authority. A nonzero drive
//! speed is an explicit work source; zero drive gives a passive internal brake.

use fs_exec::CancelGate;
use fs_tribo::{FrictionLaw, TriboError};
use crate::modal_acoustic_time::{ModalAcousticTimeModel, advance_exact_zoh};
use super::effective_compliance;
use super::super::{CoupledModalSystem, ModalAttachment, ModalCouplingError,
    dot, finite, invalid, limit, poll};

/// Signed tangential attachment; a missing right side is a prescribed surface.
#[derive(Clone, Debug)]
pub struct ModalFrictionPort {
    /// Positive side of x, receiving +B_left F.
    pub left: ModalAttachment,
    /// Optional negative side, receiving -B_right F in the SAME solve.
    pub right: Option<ModalAttachment>,
    /// Unmodified caller-declared fs-tribo law, without inferred source authority.
    pub law: FrictionLaw,
}

/// External controls held over one accepted mechanical sample.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionDrive {
    /// Signed prescribed relative surface speed [m/s]. Zero is an undriven brake.
    pub speed_m_s: f64,
    /// Prescribed compressive load [N]. Zero releases friction, not vibration.
    pub normal_force_n: f64,
}

/// Caller-selected physical limits and nonlinear-work budget; no clipping.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ModalFrictionConfig {
    /// At most this many bisections, in 1..=128.
    pub max_iterations: usize,
    /// Maximum absolute accepted tangential traction [N].
    pub maximum_force_n: f64,
    /// Maximum absolute accepted step-averaged slip speed [m/s].
    pub maximum_slip_speed_m_s: f64,
    /// Absolute force-law residual tolerance [N].
    pub force_absolute_tolerance_n: f64,
    /// Relative force-law residual tolerance, in (0,1).
    pub force_relative_tolerance: f64,
    /// Absolute displacement mismatch allowed on the sticking branch [m].
    pub sticking_tolerance_m: f64,
}

/// Resolved discrete branch, never inferred from a small velocity threshold.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum ModalFrictionRegime {
    /// The normal load was explicitly zero; the network still evolves.
    Released,
    /// The displacement constraint is satisfied within its declared tolerance.
    Sticking,
    /// The fs-tribo kinetic law is satisfied at actual step-averaged slip.
    Sliding,
}

/// Network, constitutive and step-averaged sticking failures stay distinct.
#[derive(Debug)]
pub enum ModalFrictionError {
    /// An original numerical owner or whole-system budget refused the trial.
    Network(ModalCouplingError),
    /// The unchanged fs-tribo law refused its parameters or slip speed.
    Law(TriboError),
    /// Actual endpoint displacements failed the sticking constraint.
    Sticking { /// Signed mismatch [m].
        residual_m: f64, /// Explicit tolerance [m].
        tolerance_m: f64 },
}
impl From<ModalCouplingError> for ModalFrictionError {
    fn from(value: ModalCouplingError) -> Self { Self::Network(value) }
}
impl From<TriboError> for ModalFrictionError {
    fn from(value: TriboError) -> Self { Self::Law(value) }
}
impl core::fmt::Display for ModalFrictionError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Network(e) => write!(f, "friction network: {e}"),
            Self::Law(e) => write!(f, "friction law: {e}"),
            Self::Sticking { residual_m, tolerance_m } => write!(f,
                "sticking mismatch {residual_m:e} m exceeds {tolerance_m:e} m"),
        }
    }
}
impl std::error::Error for ModalFrictionError {}

/// A complete accepted energy/work window, in SI units.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FrictionModalFrame {
    /// Original network's one-based accepted sample clock.
    pub sample: u64,
    /// Sum of the original read-only pressure observations [Pa].
    pub observer_pressure_pa: f64,
    /// Modal plus bilateral-spring stored energy [J].
    pub network_energy_j: f64,
    /// Work of external modal forces, excluding friction [J].
    pub external_work_j: f64,
    /// Work of the prescribed moving surface, F*u*dt [J]. May be negative.
    pub drive_work_j: f64,
    /// Original component and bilateral-dashpot loss [J].
    pub network_dissipation_j: f64,
    /// Nonnegative F*d on sliding; zero for sticking/released [J].
    pub friction_dissipation_j: f64,
    /// Storage change + losses - external and surface work [J].
    pub energy_residual_j: f64,
    /// Original network's scaled energy tolerance [J].
    pub energy_tolerance_j: f64,
    /// Signed traction on the positive attachment [N].
    pub traction_n: f64,
    /// Actual u*dt - B^T delta_q [m], including sticking roundoff.
    pub slip_distance_m: f64,
    /// Actual signed kinetic residual [N]; zero on other branches.
    pub constitutive_residual_n: f64,
    /// Resolved branch for this step.
    pub regime: ModalFrictionRegime,
    /// Bisection refinements used; zero for a closed-form branch or endpoint.
    pub iterations: usize,
}

/// Retained physical network with one implicitly solved tangential port.
/// Admission adds no stored energy and can preserve an already vibrating network.
/// Refusal/cancellation publishes neither candidate state nor a new frame.
pub struct FrictionModalSystem {
    network: CoupledModalSystem,
    port: ModalFrictionPort,
    config: ModalFrictionConfig,
    column: Vec<f64>,
    compliance: f64,
    forces: Vec<f64>,
    last: Option<FrictionModalFrame>,
}

impl FrictionModalSystem {
    /// Condense the admitted bilateral response once; never substitute a fixed
    /// receiver or an uncoupled effective mass for the real mechanical network.
    pub fn new(network: CoupledModalSystem, port: ModalFrictionPort,
        config: ModalFrictionConfig, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if !(1..=128).contains(&config.max_iterations)
            || [config.maximum_force_n, config.maximum_slip_speed_m_s,
                config.force_absolute_tolerance_n, config.force_relative_tolerance,
                config.sticking_tolerance_m].iter().any(|x| !x.is_finite() || *x <= 0.0)
            || config.force_relative_tolerance >= 1.0 {
            return Err(invalid("friction requires explicit positive finite physical and numerical budgets").into());
        }
        // This public owner entry validates all parameters before its zero-slip
        // early return. No regularized traction is used in the actual dynamics.
        port.law.regularized_traction_1d(0.0, 0.0, 1.0)?;
        let mut column = vec![0.0; network.mode_count()];
        for (attachment, sign) in std::iter::once((&port.left, 1.0))
            .chain(port.right.as_ref().map(|right| (right, -1.0))) {
            let model = network.models.get(attachment.component)
                .ok_or_else(|| invalid("friction attachment names an unknown component"))?;
            if attachment.shapes.len() != model.modes().len()
                || attachment.shapes.iter().any(|x| !x.is_finite()) {
                return Err(invalid("friction attachment must match its finite mass-normalized basis").into());
            }
            for (k, b) in attachment.shapes.iter().enumerate() {
                let i = network.offsets[attachment.component] + k;
                column[i] = finite(column[i] + sign*b)?;
            }
        }
        let compliance = effective_compliance(&network, &column, gate)?;
        poll(Some(gate))?;
        Ok(Self { forces: vec![0.0; network.mode_count()], network, port,
            config, column, compliance, last: None })
    }

    /// Accepted components, with their existing basis, budgets and state.
    #[must_use]
    pub fn components(&self) -> &[ModalAcousticTimeModel] { self.network.components() }
    /// Flattened external-force dimension, component then mode order.
    #[must_use]
    pub fn mode_count(&self) -> usize { self.network.mode_count() }
    /// Original time step [s].
    #[must_use]
    pub fn sample_period_s(&self) -> f64 { self.network.sample_period_s() }
    /// Original accepted sample clock; refused trials never advance it.
    #[must_use]
    pub fn samples_rendered(&self) -> u64 { self.network.samples_rendered() }
    /// Current modal plus bilateral connection energy [J].
    pub fn total_energy_j(&self) -> Result<f64, ModalCouplingError> { self.network.total_energy_j() }
    /// Last fully accepted friction-inclusive frame.
    #[must_use]
    pub const fn last_frame(&self) -> Option<&FrictionModalFrame> { self.last.as_ref() }
    /// Unmodified interface law and attachment data.
    #[must_use]
    pub const fn port(&self) -> &ModalFrictionPort { &self.port }

    /// Advance one held-input physical sample without a cancellation request.
    pub fn step(&mut self, external: &[f64], drive: FrictionDrive)
        -> Result<&FrictionModalFrame, ModalFrictionError> { self.step_inner(external, drive, None) }
    /// Same transaction with bounded cancellation polls at prediction/root/stage boundaries.
    pub fn step_under_gate(&mut self, external: &[f64], drive: FrictionDrive, gate: &CancelGate)
        -> Result<&FrictionModalFrame, ModalFrictionError> { self.step_inner(external, drive, Some(gate)) }

    fn step_inner(&mut self, external: &[f64], drive: FrictionDrive, gate: Option<&CancelGate>)
        -> Result<&FrictionModalFrame, ModalFrictionError>
    {
        poll(gate)?;
        if !drive.speed_m_s.is_finite() || !drive.normal_force_n.is_finite() || drive.normal_force_n < 0.0 {
            return Err(invalid("friction drive requires finite signed speed and nonnegative normal load").into());
        }
        let before = self.network.total_energy_j()?;
        self.network.prepare_forces(external, gate, false)?;
        let mut free_increment = 0.0;
        let mut i = 0;
        for model in &self.network.models {
            poll(gate)?;
            for (&mode, &state) in model.modes().iter().zip(model.states()) {
                let q1 = advance_exact_zoh(mode, state, self.network.forces[i], self.network.dt)
                    .displacement_m_sqrt_kg;
                free_increment = finite(free_increment + self.column[i]*(q1-self.network.old_q[i]))?;
                i += 1;
            }
        }
        let drive_distance = finite(drive.speed_m_s*self.network.dt)?;
        let free_slip = finite(drive_distance-free_increment)?;
        let (traction, regime, iterations) = self.solve(free_slip, drive.normal_force_n, gate)?;
        self.forces.copy_from_slice(external);
        if traction != 0.0 {
            for (force, b) in self.forces.iter_mut().zip(&self.column) { *force = finite(*force+b*traction)?; }
        }
        self.network.stage_inner(&self.forces, gate)?;
        let slip = finite(drive_distance-dot(&self.column, &self.network.free_delta)?)?;
        let speed = finite(slip/self.network.dt)?;
        let mut constitutive_residual_n = 0.0;
        let loss = match regime {
            ModalFrictionRegime::Released => 0.0,
            ModalFrictionRegime::Sticking => {
                if slip.abs() > self.config.sticking_tolerance_m {
                    return Err(ModalFrictionError::Sticking { residual_m: slip,
                        tolerance_m: self.config.sticking_tolerance_m });
                }
                0.0
            }
            ModalFrictionRegime::Sliding => {
                limit("friction slip speed", speed.abs(), self.config.maximum_slip_speed_m_s)?;
                if speed == 0.0 || (traction != 0.0 && traction.signum() != speed.signum()) {
                    return Err(invalid("actual sliding force must oppose nonzero relative slip").into());
                }
                let expected = finite(self.port.law.kinetic_coefficient(speed.abs())?
                    *drive.normal_force_n*speed.signum())?;
                constitutive_residual_n = finite(traction-expected)?;
                let tolerance = self.force_tolerance(traction, expected)?;
                if constitutive_residual_n.abs() > tolerance {
                    return Err(ModalCouplingError::ContactSolve { residual_n: constitutive_residual_n,
                        tolerance_n: tolerance, iterations }.into());
                }
                // Applied traction, not a nearby law value: this closes the
                // actual discrete work exchanged with the original stepper.
                finite(traction*slip)?
            }
        };
        let frame = &self.network.staged_frame;
        let after = finite(frame.modal_energy_j+frame.connection_energy_j)?;
        let network_loss = finite(frame.component_dissipation_j+frame.connection_dissipation_j)?;
        let external_work = dot(external, &self.network.free_delta)?;
        let drive_work = finite(traction*drive_distance)?;
        let residual = finite((after-before)+network_loss+loss-external_work-drive_work)?;
        let scale = before.max(after).max(network_loss).max(loss).max(external_work.abs()).max(drive_work.abs());
        let tolerance = finite(self.network.config.energy_absolute_tolerance_j
            +self.network.config.energy_relative_tolerance*scale)?;
        if loss < 0.0 || residual.abs() > tolerance {
            return Err(ModalCouplingError::EnergyBalance { residual_j: residual, tolerance_j: tolerance }.into());
        }
        let accepted = FrictionModalFrame {
            sample: frame.sample, observer_pressure_pa: frame.observer_pressure_pa,
            network_energy_j: after, external_work_j: external_work, drive_work_j: drive_work,
            network_dissipation_j: network_loss, friction_dissipation_j: loss,
            energy_residual_j: residual, energy_tolerance_j: tolerance, traction_n: traction,
            slip_distance_m: slip, constitutive_residual_n, regime, iterations,
        };
        poll(gate)?;
        self.network.publish_staged();
        Ok(self.last.insert(accepted))
    }

    fn force_tolerance(&self, applied: f64, expected: f64) -> Result<f64, ModalCouplingError> {
        finite(self.config.force_absolute_tolerance_n
            +self.config.force_relative_tolerance*applied.abs().max(expected.abs()))
    }

    fn solve(&self, free_slip: f64, normal: f64, gate: Option<&CancelGate>)
        -> Result<(f64, ModalFrictionRegime, usize), ModalFrictionError>
    {
        if normal == 0.0 { return Ok((0.0, ModalFrictionRegime::Released, 0)); }
        let static_mu = match self.port.law {
            FrictionLaw::Coulomb { static_mu, .. } | FrictionLaw::Stribeck { static_mu, .. }
                | FrictionLaw::VelocityDependent { static_mu, .. } => static_mu,
        };
        let cap = finite(static_mu*normal)?;
        let required = finite(free_slip/self.compliance)?;
        if required.abs() <= cap {
            limit("friction traction", required.abs(), self.config.maximum_force_n)?;
            return Ok((required, ModalFrictionRegime::Sticking, 0));
        }
        let sign = free_slip.signum();
        let evaluate = |magnitude: f64| -> Result<(f64, f64), ModalFrictionError> {
            poll(gate)?;
            let slip = finite(free_slip.abs()-self.compliance*magnitude)?.max(0.0);
            // Use the one-sided kinetic limit at the sticking endpoint; never
            // feed zero into the owner's strictly-positive kinetic query.
            let mu = if slip == 0.0 {
                match self.port.law {
                    FrictionLaw::Coulomb { kinetic_mu, .. } => kinetic_mu,
                    FrictionLaw::Stribeck { static_mu, .. } => static_mu,
                    FrictionLaw::VelocityDependent { mu_zero, .. } => mu_zero,
                }
            } else { self.port.law.kinetic_coefficient(finite(slip/self.network.dt)?)? };
            let expected = finite(mu*normal)?;
            Ok((finite(magnitude-expected)?, self.force_tolerance(magnitude, expected)?))
        };
        let (zero, _) = evaluate(0.0)?;
        if zero == 0.0 { return Ok((0.0, ModalFrictionRegime::Sliding, 0)); }
        let mut lo = 0.0;
        let mut hi = required.abs().min(self.config.maximum_force_n);
        let (mut residual, mut tolerance) = evaluate(hi)?;
        if residual.abs() <= tolerance { return Ok((sign*hi, ModalFrictionRegime::Sliding, 0)); }
        if residual < 0.0 {
            return Err(invalid("friction root is not bracketed within the declared force ceiling").into());
        }
        let mut used = 0;
        for iteration in 1..=self.config.max_iterations {
            used = iteration;
            let mid = f64::midpoint(lo, hi);
            (residual, tolerance) = evaluate(mid)?;
            if residual.abs() <= tolerance { return Ok((sign*mid, ModalFrictionRegime::Sliding, iteration)); }
            if mid == lo || mid == hi { break; }
            if residual < 0.0 { lo = mid; } else { hi = mid; }
        }
        Err(ModalCouplingError::ContactSolve { residual_n: residual, tolerance_n: tolerance, iterations: used }.into())
    }
}
