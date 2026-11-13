//! Finite-amplitude internal resistance through the original physical owners.
//!
//! xdot = (J-R) grad H(x) - D(x,grad H(x)) + G u(t). The SAME nonlinear
//! resisting law used by fs-phs time stepping and equilibrium linearization is
//! evaluated at every original harmonic/time sample. Its differential slope
//! is not replaced by D/e, and its work is not reported as external supply.

use super::{PeriodicPortInput, PhsStorageHessianAction, tangent::StorageTangent};
use core::cell::RefCell;
use fs_math::dd::Dd;
use fs_orbit::{
    AdaptiveHbOrbit, HbAnchor, HbBudget, HbRefinementBudget, HbScales, OrbitError, OrbitProblem,
    ScaledHbOrbit, ShootBudget, VariationalShootControls, VariationalShootError,
    VariationalShootOrbit, solve_hb_adaptive_scaled_seeded, solve_hb_scaled_seeded,
    solve_shooting_variational,
};
pub use fs_phs::harmonic::equilibrium::NonlinearDissipation;
use fs_phs::{PortHamiltonian, harmonic::equilibrium::LinearizationError};

/// A complete source or native solver refusal. No partial result is published.
#[derive(Clone, Debug, PartialEq)]
pub enum DissipativeError {
    Invalid(&'static str),
    SourceRefused,
    Dissipation(LinearizationError),
    Orbit(OrbitError),
    Shooting(VariationalShootError),
}
impl core::fmt::Display for DissipativeError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "nonlinear dissipative periodic source refused: {self:?}")
    }
}
impl std::error::Error for DissipativeError {}
type Result<T> = core::result::Result<T, DissipativeError>;

/// Separate instantaneous physical accounts of the original source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DissipativePower {
    pub energy: f64,
    /// Signed external input e^T G u.
    pub port_power: f64,
    /// Original linear resistance e^T R e, preserving raw signed R.
    pub resistance_power: f64,
    /// Original nonlinear loss e^T D, not an incremental/tangent power.
    /// Only roundoff-sized negative values can pass the source admission.
    pub nonlinear_power: f64,
    /// Independently contracted grad H dot xdot.
    pub energy_rate: f64,
}

struct Scratch {
    effort: Vec<f64>,
    input: Vec<f64>,
    resistance: Vec<f64>,
    state_direction: Vec<f64>,
    effort_direction: Vec<f64>,
    resistance_direction: Vec<f64>,
    flow: Vec<f64>,
}

/// Prepared full-state nonlinear resistance and periodic driving.
///
/// All constitutive memory remains in the original states. The original H, J,
/// R and G are borrowed unchanged. D must be a read-only smooth law with its
/// complete action D_x dx + D_e de; hidden history and one-sided corner slopes
/// are unsupported. Pointwise nonnegative e dot D does not imply nonnegative
/// incremental damping. A real negative differential slope is retained.
///
/// Preparation admits at most 256 states and 256 physical ports and retains
/// (6*n + m) scalar scratch entries plus the existing storage-tangent owner.
/// The supplied max_tangent_actions must cover all n original directions.
/// Native HB/shooting budgets separately bound their numerical work; callbacks'
/// own storage/work are not included. Warm source evaluation allocates nothing
/// in this adapter. Source callbacks and original storage may have their own
/// allocation behavior. The prepared value is intentionally not shared mutable
/// state across threads.
pub struct DissipativeHamiltonian<'a> {
    system: &'a PortHamiltonian,
    input: &'a PeriodicPortInput,
    dissipation: NonlinearDissipation<'a>,
    tangent: StorageTangent<'a>,
    scratch: RefCell<Scratch>,
    failure: RefCell<Option<DissipativeError>>,
}

impl<'a> DissipativeHamiltonian<'a> {
    pub fn new(
        system: &'a PortHamiltonian,
        input: &'a PeriodicPortInput,
        dissipation: NonlinearDissipation<'a>,
    ) -> Result<Self> {
        Self::prepare(system, input, dissipation, None)
    }

    /// Supply an exact COMPLETE Hessian action when the original storage does
    /// not declare its own. The action includes every material-memory row.
    pub fn new_analytic(
        system: &'a PortHamiltonian,
        input: &'a PeriodicPortInput,
        dissipation: NonlinearDissipation<'a>,
        hessian: &'a PhsStorageHessianAction<'a>,
    ) -> Result<Self> {
        Self::prepare(system, input, dissipation, Some(hessian))
    }

    fn prepare(
        system: &'a PortHamiltonian,
        input: &'a PeriodicPortInput,
        dissipation: NonlinearDissipation<'a>,
        hessian: Option<&'a PhsStorageHessianAction<'a>>,
    ) -> Result<Self> {
        let n = system.state_dim();
        let m = system.port_dim();
        let period = core::f64::consts::TAU / input.omega_rad_s();
        if !(1..=256).contains(&n)
            || m > 256
            || input.port_count() != m
            || dissipation.max_tangent_actions < n
            || !period.is_finite()
            || period <= 0.0
        {
            return Err(DissipativeError::Invalid(
                "bounded original states, complete ports, tangent actions and representable period required",
            ));
        }
        let tangent = match hessian {
            Some(action) => StorageTangent::hessian(n, action),
            None => StorageTangent::from_storage(system),
        };
        if !tangent.available() {
            return Err(DissipativeError::Invalid(
                "complete original storage tangent required",
            ));
        }
        Ok(Self {
            system,
            input,
            dissipation,
            tangent,
            scratch: RefCell::new(Scratch {
                effort: vec![0.0; n],
                input: vec![0.0; m],
                resistance: vec![0.0; n],
                state_direction: vec![0.0; n],
                effort_direction: vec![0.0; n],
                resistance_direction: vec![0.0; n],
                flow: vec![0.0; n],
            }),
            failure: RefCell::new(None),
        })
    }

    fn evaluate(&self, time: f64, state: &[f64], scratch: &mut Scratch) -> Result<f64> {
        if !time.is_finite()
            || state.len() != self.system.state_dim()
            || state.iter().any(|v| !v.is_finite())
        {
            return Err(DissipativeError::Invalid(
                "complete finite original state and time required",
            ));
        }
        scratch.effort.fill(f64::NAN);
        self.system.storage().gradient(state, &mut scratch.effort);
        if scratch.effort.iter().any(|v| !v.is_finite())
            || !self.input.sample_into(time, &mut scratch.input)
        {
            return Err(DissipativeError::SourceRefused);
        }
        self.dissipation
            .evaluate_into(state, &scratch.effort, &mut scratch.resistance, || false)
            .map_err(DissipativeError::Dissipation)
    }

    /// Original full nonlinear field; all four terms share one final rounding.
    pub fn vector_field_into(&self, time: f64, state: &[f64], out: &mut [f64]) -> Result<()> {
        let mut scratch = self.scratch.borrow_mut();
        self.evaluate(time, state, &mut scratch)?;
        self.system
            .vector_field_with_dissipation_from_effort_into(
                &scratch.effort,
                &scratch.input,
                &scratch.resistance,
                out,
            )
            .map_err(|_| DissipativeError::SourceRefused)
    }

    /// Full continuous tangent (J-R) Hxx - D_x - D_e Hxx. State and effort
    /// directions are independent callback arguments and their chain is added
    /// exactly once. No secant coefficient, symmetrization or finite difference.
    pub fn state_jacobian_into(&self, time: f64, state: &[f64], out: &mut [f64]) -> Result<()> {
        let n = self.system.state_dim();
        if out.len() != n * n {
            return Err(DissipativeError::Invalid(
                "complete row-major state Jacobian required",
            ));
        }
        let mut scratch = self.scratch.borrow_mut();
        self.evaluate(time, state, &mut scratch)?;
        let Scratch {
            effort,
            input,
            state_direction,
            effort_direction,
            resistance_direction,
            flow,
            ..
        } = &mut *scratch;
        input.fill(0.0);
        state_direction.fill(0.0);
        for column in 0..n {
            state_direction[column] = 1.0;
            if !self
                .tangent
                .apply_hessian(state, state_direction, effort_direction)
            {
                return Err(DissipativeError::SourceRefused);
            }
            resistance_direction.fill(f64::NAN);
            if !(self.dissipation.tangent)(
                state,
                effort,
                state_direction,
                effort_direction,
                resistance_direction,
            ) || resistance_direction.iter().any(|v| !v.is_finite())
            {
                return Err(DissipativeError::Dissipation(
                    LinearizationError::DissipationRefused {
                        column: Some(column),
                    },
                ));
            }
            self.system
                .vector_field_with_dissipation_from_effort_into(
                    effort_direction,
                    input,
                    resistance_direction,
                    flow,
                )
                .map_err(|_| DissipativeError::SourceRefused)?;
            for row in 0..n {
                out[row * n + column] = flow[row];
            }
            state_direction[column] = 0.0;
        }
        Ok(())
    }

    /// Independent work-conjugate powers of the SAME complete finite-amplitude
    /// source. This is an instantaneous identity, not a cycle-energy, temporal
    /// convergence, passivity or stability certificate.
    pub fn power(&self, time: f64, state: &[f64]) -> Result<DissipativePower> {
        let mut scratch = self.scratch.borrow_mut();
        let nonlinear_power = self.evaluate(time, state, &mut scratch)?;
        let Scratch {
            effort,
            input,
            resistance,
            flow,
            ..
        } = &mut *scratch;
        self.system
            .vector_field_with_dissipation_from_effort_into(effort, input, resistance, flow)
            .map_err(|_| DissipativeError::SourceRefused)?;
        let energy = self.system.hamiltonian(state);
        let port_power = self
            .system
            .supplied_power_from_effort(effort, input)
            .map_err(|_| DissipativeError::SourceRefused)?;
        let resistance_power = self
            .system
            .resistance_power_from_effort(effort)
            .map_err(|_| DissipativeError::SourceRefused)?;
        let energy_rate = effort
            .iter()
            .zip(flow)
            .fold(Dd::ZERO, |s, (&e, d)| {
                s + Dd::from_f64(e) * Dd::from_f64(*d)
            })
            .to_f64();
        if [
            energy,
            port_power,
            resistance_power,
            nonlinear_power,
            energy_rate,
        ]
        .iter()
        .any(|v| !v.is_finite())
        {
            return Err(DissipativeError::SourceRefused);
        }
        Ok(DissipativePower {
            energy,
            port_power,
            resistance_power,
            nonlinear_power,
            energy_rate,
        })
    }

    fn admit_order(&self, harmonics: usize) -> Result<()> {
        if self.input.highest_nonzero_harmonic() > harmonics {
            Err(DissipativeError::Invalid(
                "every nonzero drive harmonic must fit retained order",
            ))
        } else {
            Ok(())
        }
    }

    fn finish<T>(&self, result: Result<T>) -> Result<T> {
        match self.failure.borrow_mut().take() {
            Some(failure) => Err(failure),
            None => result,
        }
    }

    /// Forced nonlinear harmonic balance of every original coordinate. The
    /// original clock fixes phase; all nonlinear response harmonics are solved
    /// together. Fixed-order sampled convergence does not bound omitted modes.
    pub fn solve_harmonic_balance_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        scales: &HbScales,
    ) -> Result<ScaledHbOrbit> {
        self.admit_order(budget.harmonics)?;
        *self.failure.borrow_mut() = None;
        self.finish(
            solve_hb_scaled_seeded(
                &Adapter(self),
                HbAnchor::Forced {
                    omega: self.input.omega_rad_s(),
                },
                seed,
                budget,
                scales,
            )
            .map_err(DissipativeError::Orbit),
        )
    }

    /// Existing full-equation/successive-waveform harmonic ladder, with the
    /// complete original nonlinear loss law evaluated at every resolved sample.
    pub fn solve_harmonic_balance_adaptive_scaled_seeded(
        &self,
        seed: Vec<f64>,
        budget: &HbBudget,
        refinement: &HbRefinementBudget,
        scales: &HbScales,
    ) -> Result<AdaptiveHbOrbit> {
        self.admit_order(budget.harmonics)?;
        *self.failure.borrow_mut() = None;
        self.finish(
            solve_hb_adaptive_scaled_seeded(
                &Adapter(self),
                HbAnchor::Forced {
                    omega: self.input.omega_rad_s(),
                },
                seed,
                budget,
                refinement,
                scales,
            )
            .map_err(DissipativeError::Orbit),
        )
    }

    /// Independent original-ODE shooting and its exact discrete RK4 monodromy.
    /// Refine the time grid separately. All nonlinear-loss/memory derivatives
    /// participate; no effective damping or HB tangent is substituted.
    pub fn solve_shooting_variational(
        &self,
        seed: &[f64],
        budget: &ShootBudget,
        controls: &VariationalShootControls,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<VariationalShootOrbit> {
        if cancelled() {
            return Err(DissipativeError::Shooting(VariationalShootError::Cancelled));
        }
        if budget.steps_per_period == 0
            || self.input.highest_nonzero_harmonic() > (budget.steps_per_period - 1) / 2
        {
            return Err(DissipativeError::Invalid(
                "shooting grid must resolve every drive harmonic",
            ));
        }
        *self.failure.borrow_mut() = None;
        let result = solve_shooting_variational(
            &Adapter(self),
            seed,
            core::f64::consts::TAU / self.input.omega_rad_s(),
            budget,
            controls,
            &mut cancelled,
        );
        if matches!(result, Err(VariationalShootError::Cancelled)) {
            return Err(DissipativeError::Shooting(VariationalShootError::Cancelled));
        }
        self.finish(result.map_err(DissipativeError::Shooting))
    }
}

struct Adapter<'s, 'a>(&'s DissipativeHamiltonian<'a>);
impl Adapter<'_, '_> {
    fn record(&self, result: Result<()>) -> bool {
        if let Err(error) = result {
            let mut failure = self.0.failure.borrow_mut();
            if failure.is_none() {
                *failure = Some(error);
            }
            false
        } else {
            true
        }
    }
    fn jacobian(&self, time: f64, state: &[f64], out: &mut [f64]) -> bool {
        self.record(self.0.state_jacobian_into(time, state, out))
    }
}
impl OrbitProblem for Adapter<'_, '_> {
    fn dim(&self) -> usize {
        self.0.system.state_dim()
    }
    fn autonomous(&self) -> bool {
        false
    }
    fn island(&self, time: f64, state: &[f64], out: &mut [f64]) {
        if !self.record(self.0.vector_field_into(time, state, out)) {
            out.fill(f64::NAN);
        }
    }
    fn island_state_jacobian(&self) -> Option<fs_orbit::IslandStateJacobian<Self>> {
        Some(Self::jacobian)
    }
    fn time_domain_rhs(&self) -> Option<fs_orbit::TimeDomainRhs<Self>> {
        Some(Self::island)
    }
    fn time_domain_state_jacobian(&self) -> Option<fs_orbit::TimeDomainStateJacobian<Self>> {
        Some(Self::jacobian)
    }
}

#[cfg(test)]
mod tests;
