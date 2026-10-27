//! One chemical, caloric and inertial-wall endpoint on the original gas owner.
//!
//! The endpoint gas pressure drives the native constant-load impulse and
//! midpoint displacement. Its actual p dV work enters the original implicit
//! reaction equations at the same volume and total internal energy. This is
//! first order in the gas and chemistry; no midpoint gas receipt is fabricated.

use super::boundary::{BoundaryGasReactionCandidate, BoundaryGasReactionStep, GasReactionBoundary};
use super::{
    GasPiston, GasReaction, GasReactionError, GasReactionSolve, PistonError, PistonInput,
    PistonSolve,
};
use crate::cavity_phs::piston::{MechanicalTrial, PistonSupport, thermal::feedback};
use crate::port_adoption::{PortTransfer, StepAccount};

/// Original mechanical bounds plus the endpoint pressure equation's bracket.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoupledReactivePistonSolve {
    /// Volume, clock, gas-source and physical residual bounds. Its existing
    /// maximum_iterations bounds scalar pressure trials (at most cap + 2).
    pub piston: PistonSolve,
    /// Positive absolute pressures [Pa], containing the accepted initial p.
    /// Every queried pressure must produce an admitted volume and gas state.
    pub pressure_interval_pa: [f64; 2],
    /// Independent maximum |p_EOS,end - p_applied| [Pa].
    pub pressure_tolerance_pa: f64,
}

/// Original source failures and the actual coupled pressure-root refusal.
#[derive(Clone, Debug, PartialEq)]
pub enum CoupledReactivePistonError {
    Invalid(&'static str),
    Reaction(GasReactionError),
    Piston(PistonError),
    OutsidePressureBracket {
        lower_residual_pa: f64,
        upper_residual_pa: f64,
    },
    NotConverged {
        iterations: usize,
        residual_pa: f64,
    },
    Cancelled,
}
impl From<GasReactionError> for CoupledReactivePistonError {
    fn from(error: GasReactionError) -> Self {
        if error == GasReactionError::Cancelled {
            Self::Cancelled
        } else {
            Self::Reaction(error)
        }
    }
}
impl From<PistonError> for CoupledReactivePistonError {
    fn from(error: PistonError) -> Self {
        if error == PistonError::Cancelled {
            Self::Cancelled
        } else {
            Self::Piston(error)
        }
    }
}
impl core::fmt::Display for CoupledReactivePistonError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "coupled reactive piston: {self:?}")
    }
}
impl std::error::Error for CoupledReactivePistonError {}

/// Actual same-interval chemical, mechanical and extensive energy balances.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledReactivePistonStep {
    pub reaction: BoundaryGasReactionStep,
    /// Pressure held during the native wall impulse and work [Pa]. At
    /// acceptance it agrees with the original endpoint EOS to declared bounds.
    pub held_pressure_pa: f64,
    /// p_EOS,end - held_pressure_pa [Pa].
    pub pressure_residual_pa: f64,
    /// Pressure mismatch times actual volume change [J], independently bounded.
    pub pressure_eos_work_residual_j: f64,
    pub displacement_m: f64,
    pub velocity_m_s: [f64; 2],
    pub kinetic_energy_change_j: f64,
    pub exterior_pressure_work_j: f64,
    pub external_force_work_j: f64,
    pub support_impulse_n_s: f64,
    pub support_work_j: f64,
    pub momentum_residual_n_s: f64,
    pub mechanical_energy_residual_j: f64,
    /// Delta U_gas + Delta K + exterior work - force/support work - heat [J].
    /// Formation energy is already in U_gas; no chemical heat is added again.
    pub total_energy_residual_j: f64,
    /// Actual pressure proposals, including the initial and bracket probes.
    pub trials: usize,
    /// Original implicit mechanism evaluations over ALL pressure proposals,
    /// including each nested Jacobian and rejected chemical trial.
    pub chemical_evaluations: usize,
}

/// Complete native gas, volume, velocity and port history; no input is mutated.
#[derive(Clone, Debug, PartialEq)]
pub struct CoupledReactivePistonCandidate<P = GasPiston> {
    piston: P,
    step: CoupledReactivePistonStep,
}
impl<P> CoupledReactivePistonCandidate<P> {
    pub const fn piston(&self) -> &P {
        &self.piston
    }
    pub const fn step(&self) -> &CoupledReactivePistonStep {
        &self.step
    }
    pub fn into_parts(self) -> (P, CoupledReactivePistonStep) {
        (self.piston, self.step)
    }
}
impl CoupledReactivePistonCandidate {
    pub(in crate::cavity_phs::piston) fn map_piston<P>(
        self,
        bind: impl FnOnce(GasPiston) -> P,
    ) -> CoupledReactivePistonCandidate<P> {
        CoupledReactivePistonCandidate {
            piston: bind(self.piston),
            step: self.step,
        }
    }
}

type Error = CoupledReactivePistonError;
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<(), Error> {
    if cancelled() {
        Err(Error::Cancelled)
    } else {
        Ok(())
    }
}
fn feedback_error(error: feedback::ScalarError<Error>) -> Error {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => Error::Invalid(message),
        feedback::ScalarError::Cancelled => Error::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => Error::OutsidePressureBracket {
            lower_residual_pa: lower,
            upper_residual_pa: upper,
        },
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => Error::NotConverged {
            iterations,
            residual_pa: residual,
        },
    }
}

struct Trial {
    chemical: BoundaryGasReactionCandidate,
    motion: MechanicalTrial,
    pressure_pa: f64,
    residual_pa: f64,
    trials: usize,
}

impl GasPiston {
    /// Solve reaction rates, gas heat/work, EOS pressure and native momentum
    /// together. No reaction-then-motion stage or second gas inventory exists.
    ///
    /// At each proposed p, the original inertial wall supplies its constant-load
    /// impulse and midpoint displacement. The implicit chemical owner then uses
    /// V1 and U1=U0+Q-p(V1-V0). The shared scalar solver closes p=p_EOS,end.
    /// Every rate and source evaluates the actual final composition/T/pressure.
    /// Fixed support retains its original volume and zero wall velocity, with
    /// an explicit balancing support impulse and exactly zero support work.
    ///
    /// Requires `GasReaction::with_implicit_mechanism`, one shared gas-inversion
    /// policy and no chemistry-only subcycling. Refine the whole interval to
    /// assess this first-order method. Pressure bracketing is not a uniqueness,
    /// monotonicity or unconditional convergence claim. Source and volume
    /// bounds apply to every query, not merely to the accepted endpoint.
    ///
    /// Independent pressure, impulse, work, total-energy and signed heat-entropy
    /// admission precede return. Neither numerical residual nor reaction
    /// enthalpy is deposited as heat. No midpoint gas, analytic constitutive
    /// tangent, spatial flame/shock or entropy-production isolation is claimed.
    pub fn coupled_reactive_candidate(
        &self,
        reaction: &GasReaction,
        input: PistonInput,
        chemistry: GasReactionSolve,
        controls: CoupledReactivePistonSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<CoupledReactivePistonCandidate, Error> {
        poll(&mut cancelled)?;
        let native = controls.piston;
        let p0 = self.gas.state().state().pressure;
        let v0 = self.gas.volume_m3();
        let [pmin, pmax] = controls.pressure_interval_pa;
        let [vmin, vmax] = native.volume_bracket_m3;
        if !positive(input.duration_s)
            || !input.heat_j.is_finite()
            || !input.exterior_pressure_pa.is_finite()
            || input.exterior_pressure_pa < 0.0
            || !input.external_force_n.is_finite()
            || !positive(pmin)
            || !positive(pmax)
            || pmin >= pmax
            || p0 < pmin
            || p0 > pmax
            || !positive(vmin)
            || !positive(vmax)
            || vmin > vmax
            || (vmin == vmax && self.support == PistonSupport::Free)
            || v0 < vmin
            || v0 > vmax
            || !positive(controls.pressure_tolerance_pa)
            || !positive(native.momentum_tolerance_n_s)
            || !positive(native.energy_tolerance_j)
            || !(1..=128).contains(&native.maximum_iterations)
            || chemistry.gas != native.gas
        {
            return Err(Error::Invalid(
                "bounded common pressure, volume, clock and reaction/gas controls required",
            ));
        }
        chemistry.validate(input.duration_s)?;
        self.admit_linear_clock(self.gas.state(), v0, input.duration_s)?;
        let mut chemical_evaluations = 0usize;
        let trial = feedback::solve_scalar(
            [pmin - p0, pmax - p0],
            controls.pressure_tolerance_pa,
            native.maximum_iterations,
            |delta, trials, cancelled| {
                poll(cancelled)?;
                let pressure_pa = p0 + delta;
                if !positive(pressure_pa) || pressure_pa < pmin || pressure_pa > pmax {
                    return Err(Error::Invalid("coupled pressure proposal left its bracket"));
                }
                let volume = if self.support == PistonSupport::Fixed {
                    v0
                } else {
                    let force = self.area_m2 * (pressure_pa - input.exterior_pressure_pa)
                        + input.external_force_n;
                    let velocity = self.velocity_m_s + input.duration_s * force / self.mass_kg;
                    let displacement =
                        input.duration_s * f64::midpoint(self.velocity_m_s, velocity);
                    v0 + self.area_m2 * displacement
                };
                if !volume.is_finite() || volume < vmin || volume > vmax {
                    return Err(Error::Invalid(
                        "native reacting wall motion left its volume bracket",
                    ));
                }
                let motion = self.mechanical_trial(volume, pressure_pa, input)?;
                let work_by_gas_j = pressure_pa * (volume - v0);
                let chemical = reaction.mechanical_implicit_trial(
                    &self.gas,
                    input.duration_s,
                    GasReactionBoundary {
                        final_volume_m3: volume,
                        heat_j: input.heat_j,
                        work_by_gas_j,
                    },
                    chemistry,
                    &mut *cancelled,
                )?;
                chemical_evaluations = chemical_evaluations
                    .checked_add(
                        chemical
                            .step()
                            .reaction
                            .implicit_mechanism
                            .as_ref()
                            .expect("original implicit endpoint retains its work report")
                            .evaluations,
                    )
                    .ok_or(Error::Invalid("coupled chemical evaluation count overflow"))?;
                let residual_pa = chemical.gas().state().state().pressure - pressure_pa;
                if !residual_pa.is_finite() {
                    return Err(Error::Invalid(
                        "coupled endpoint pressure residual is unrepresentable",
                    ));
                }
                // The shared bracket owner uses applied - EOS pressure.
                Ok((
                    Trial {
                        chemical,
                        motion,
                        pressure_pa,
                        residual_pa,
                        trials,
                    },
                    -residual_pa,
                ))
            },
            &mut cancelled,
        )
        .map_err(feedback_error)?;
        poll(&mut cancelled)?;
        let Trial {
            chemical,
            motion,
            pressure_pa,
            residual_pa,
            trials,
        } = trial;
        let (gas, reaction) = chemical.into_parts();
        self.admit_linear_clock(gas.state(), gas.volume_m3(), input.duration_s)?;
        let before = self.velocity_m_s;
        let kinetic =
            (0.5 * self.mass_kg * (motion.velocity - before)) * (motion.velocity + before);
        let exterior = input.exterior_pressure_pa * (gas.volume_m3() - v0);
        let external = input.external_force_n * motion.displacement;
        let support_work = motion.support_impulse * f64::midpoint(before, motion.velocity);
        let mechanical =
            kinetic - reaction.boundary.work_by_gas_j + exterior - external - support_work;
        let total = reaction.reaction.stored_energy_change_j + kinetic + exterior
            - external
            - support_work
            - input.heat_j;
        let pressure_work = residual_pa * (gas.volume_m3() - v0);
        if residual_pa.abs() > controls.pressure_tolerance_pa
            || motion.residual.abs() > native.momentum_tolerance_n_s
            || [
                kinetic,
                exterior,
                external,
                support_work,
                mechanical,
                total,
                pressure_work,
            ]
            .iter()
            .any(|v| !v.is_finite())
            || mechanical.abs() > native.energy_tolerance_j
            || total.abs() > native.energy_tolerance_j
            || pressure_work.abs() > native.energy_tolerance_j
        {
            return Err(Error::Invalid(
                "coupled pressure, native impulse or actual gas/wall work did not close",
            ));
        }
        if reaction.entropy_residual_j_k < -chemistry.entropy_tolerance_j_k {
            return Err(GasReactionError::EntropyDecrease {
                change_j_k: reaction.entropy_residual_j_k,
                tolerance_j_k: chemistry.entropy_tolerance_j_k,
            }
            .into());
        }
        let step = CoupledReactivePistonStep {
            reaction,
            held_pressure_pa: pressure_pa,
            pressure_residual_pa: residual_pa,
            pressure_eos_work_residual_j: pressure_work,
            displacement_m: motion.displacement,
            velocity_m_s: [before, motion.velocity],
            kinetic_energy_change_j: kinetic,
            exterior_pressure_work_j: exterior,
            external_force_work_j: external,
            support_impulse_n_s: motion.support_impulse,
            support_work_j: support_work,
            momentum_residual_n_s: motion.residual,
            mechanical_energy_residual_j: mechanical,
            total_energy_residual_j: total,
            trials,
            chemical_evaluations,
        };
        let mut next = self.clone();
        next.gas = gas;
        next.velocity_m_s = motion.velocity;
        next.ports.accept(
            StepAccount {
                tick: next.ports.next_tick(),
                duration_s: input.duration_s,
                storage_change_j: Some(step.reaction.reaction.stored_energy_change_j + kinetic),
                loss_j: Some(0.0),
                included_work_j: 0.0,
                entropy_change_j_k: step.reaction.reaction.entropy_change_j_k,
                tolerance_j: Some(native.energy_tolerance_j),
                unavailable: None,
            },
            [
                PortTransfer::work(external - exterior + support_work, true),
                PortTransfer::heat(input.heat_j, next.gas.state().state().temperature),
            ],
        );
        poll(&mut cancelled)?;
        Ok(CoupledReactivePistonCandidate { piston: next, step })
    }
}
