//! A sliding interface under a prescribed slip rate, on the composition
//! contract.
//!
//! A [`TangentialLaw`] owns the traction, the solved slip, the regime, any
//! reversible storage and the dissipation it exports. This adapter drives
//! one with a released slip increment of `rate * dt` at zero compliance
//! (prescribed motion) and with the interface temperature its heat port
//! receives, and hands the dissipation to that port. It adds no heat
//! fraction, no heat capacity and no flash temperature. Where the heat goes
//! is decided by what the port is connected to: a thermal node between two
//! bodies divides it by their conductances and their temperatures.
//!
//! Port `heat` receives the interface temperature [K] and returns the heat
//! rate into the interface [W], which is minus the dissipation rate. Port
//! `traction` is a signal output of two values, the traction entering the
//! interface [N] along the two tangent axes.
//!
//! The law at a temperature. Friction that depends on temperature enters as
//! the caller's source `law_at(T)`, which returns the law to use at the
//! interface temperature of a step; a law that is itself bounded in
//! temperature also receives that temperature as its thermal boundary. The
//! source is called at every evaluation of a step and must be a pure
//! function.
//!
//! Work. No mechanical port is offered. A dry law is set-valued in the slip
//! rate, and a kernel that differences smooth port maps cannot solve for the
//! slip through a port. The motion is prescribed, so the work it does
//! against the traction is the component's own account, as a prescribed
//! force's is: the stored energy is the law's reversible energy less the
//! work received so far, and a step's receipt books the law's closure
//! residual `storage + dissipation - work` as numerical loss.
//!
//! Entropy generated: the dissipation over the interface temperature.
//!
//! Class: `Smooth` in its one input, as far as the source and the law are
//! smooth in temperature at a prescribed slip. No tangent is offered; a
//! kernel differences the map.
//!
//! Refused: an interface temperature or a step that is not positive and
//! finite; whatever the source or the law refuses, unchanged; a step whose
//! work balance is off by more than the declared tolerance, or whose
//! reported stored change disagrees with the law's stored energies by more
//! than it (the admission of the hand-written friction-heat coupling's
//! `tangential_heat_candidate`); a source whose law gives the held state
//! another stored energy than the state was left with. A reversible law
//! whose stiffness depends on temperature needs a thermoelastic account that
//! this adapter does not have.
//!
//! One call does a bounded amount of work (the law's own root budget), and
//! [`TriboError`] has no cancellation value, so `cx` is not polled here.

use core::fmt;
use core::marker::PhantomData;
use core::num::NonZeroUsize;

use fs_exec::Cx;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortRole, PortSpec, Proposal, Receipt};

use super::{TangentialLaw, TangentialSolveConfig, TangentialStepInput};
use crate::TriboError;

/// The load and motion a [`SlidingInterface`] is held to.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PrescribedSlip {
    /// Compressive normal load [N], non-negative.
    pub normal_n: f64,
    /// Relative slip rate along the two tangent axes [m/s].
    pub slip_rate_m_s: [f64; 2],
    /// Nominal pressure [Pa], for a law that is bounded in it.
    pub nominal_pressure_pa: Option<f64>,
    /// The law's numerical and effort budgets.
    pub solve: TangentialSolveConfig,
    /// Largest admitted defect of a step's work balance [J].
    pub work_tolerance_j: f64,
}

/// State of a [`SlidingInterface`]: the law's own checkpoint and the two
/// terms of the stored energy.
#[derive(Debug, Clone, PartialEq)]
pub struct SlidingState<S> {
    tangential: S,
    /// The law's reversible energy at `tangential` [J].
    reversible_j: f64,
    /// Work the prescribed motion has done on the interface so far [J].
    work_j: f64,
}

impl<S> SlidingState<S> {
    /// The law's checkpoint: its slip, its thermal boundary, its layers.
    #[must_use]
    pub const fn tangential(&self) -> &S {
        &self.tangential
    }

    /// Work the prescribed motion has done on the interface so far [J].
    #[must_use]
    pub const fn work_j(&self) -> f64 {
        self.work_j
    }
}

/// A sliding interface whose slip rate and normal load are prescribed, with
/// the law `law_at` returns for the interface temperature of each step.
pub struct SlidingInterface<L, F> {
    law_at: F,
    drive: PrescribedSlip,
    ports: [PortSpec; 2],
    law: PhantomData<fn() -> L>,
}

impl<L, F> fmt::Debug for SlidingInterface<L, F> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("SlidingInterface")
            .field("drive", &self.drive)
            .finish_non_exhaustive()
    }
}

fn refused(field: &'static str) -> TriboError {
    TriboError::InvalidInput { field }
}

fn finite(value: f64, field: &'static str) -> Result<f64, TriboError> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(refused(field))
    }
}

/// One evaluated step: the candidate state, the heat rate into the
/// interface [W], the traction [N] and the receipt.
type Stepped<S> = (SlidingState<S>, f64, [f64; 2], Receipt);

impl<L, F> SlidingInterface<L, F>
where
    L: TangentialLaw,
    F: Fn(f64) -> Result<L, TriboError>,
{
    /// An interface held to `drive`, with the law `law_at` gives at an
    /// interface temperature [K].
    ///
    /// # Errors
    /// A normal load or nominal pressure that is negative or not finite, a
    /// slip rate that is not finite, a work tolerance that is negative or
    /// not finite, or solve budgets the law's owner refuses.
    pub fn new(law_at: F, drive: PrescribedSlip) -> Result<Self, TriboError> {
        drive.solve.validate()?;
        crate::nonnegative_finite(drive.normal_n, "normal_force")?;
        if let Some(pressure_pa) = drive.nominal_pressure_pa {
            crate::nonnegative_finite(pressure_pa, "contact_pressure_pa")?;
        }
        for rate in drive.slip_rate_m_s {
            finite(rate, "slip_rate_m_s")?;
        }
        crate::nonnegative_finite(drive.work_tolerance_j, "work_tolerance_j")?;
        let name = |id: &str| StableId::new(id).map_err(|_| refused("port_identifier"));
        Ok(Self {
            law_at,
            drive,
            ports: [
                PortSpec::heat(name("heat")?, PortRole::EffortIn)
                    .map_err(|_| refused("heat_port"))?,
                PortSpec::signal_out(
                    name("traction")?,
                    NonZeroUsize::new(2).ok_or(refused("traction_port"))?,
                ),
            ],
            law: PhantomData,
        })
    }

    /// The load and motion this interface is held to.
    #[must_use]
    pub const fn drive(&self) -> &PrescribedSlip {
        &self.drive
    }

    /// The unloaded state of the law at `interface_temperature_k`, with no
    /// work received: the state to hand to a kernel with this component.
    ///
    /// # Errors
    /// A temperature that is not positive and finite, the source's refusal,
    /// or a reversible energy that is not finite.
    pub fn start(&self, interface_temperature_k: f64) -> Result<SlidingState<L::State>, TriboError> {
        let law = self.law_at_k(interface_temperature_k)?;
        let tangential = law.initial_state();
        let reversible_j = finite(law.stored_energy_j(&tangential), "stored_energy_j")?;
        Ok(SlidingState {
            tangential,
            reversible_j,
            work_j: 0.0,
        })
    }

    fn law_at_k(&self, interface_temperature_k: f64) -> Result<L, TriboError> {
        if !(interface_temperature_k.is_finite() && interface_temperature_k > 0.0) {
            return Err(refused("interface_temperature_k"));
        }
        (self.law_at)(interface_temperature_k)
    }

    /// The law's step from `old` over `dt_s` at the interface temperature
    /// `inputs[0]`, admitted as the friction-heat coupling admits one.
    fn step(
        &self,
        old: &SlidingState<L::State>,
        dt_s: f64,
        inputs: &[f64],
    ) -> Result<Stepped<L::State>, TriboError> {
        let &[kelvin] = inputs else {
            return Err(refused("sliding_interface_inputs"));
        };
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(refused("dt_s"));
        }
        let law = self.law_at_k(kelvin)?;
        let tolerance_j = self.drive.work_tolerance_j;
        // The held state must be worth to this law what it was left worth:
        // otherwise changing the law has created or destroyed energy.
        let before_j = law.stored_energy_j(&old.tangential);
        if !((before_j - old.reversible_j).abs() <= tolerance_j) {
            return Err(refused("sliding_interface_law_changed_stored_energy"));
        }
        let step = law.advance(
            &old.tangential,
            &TangentialStepInput {
                normal_n: self.drive.normal_n,
                slip_increment_m: self.drive.slip_rate_m_s.map(|rate| rate * dt_s),
                dt_s,
                contact_temperature_k: Some(kelvin),
                nominal_pressure_pa: self.drive.nominal_pressure_pa,
                // Prescribed motion: the bodies do not give under traction.
                compliance_m_per_n: 0.0,
                solve: self.drive.solve,
            },
        )?;
        // Work is traction times the solved slip, never the released motion.
        let work_j = step.traction_n[0] * step.slip_increment_m[0]
            + step.traction_n[1] * step.slip_increment_m[1];
        let after_j = law.stored_energy_j(&step.state);
        let storage_j = after_j - before_j;
        let closure_j = storage_j + step.dissipated_j - work_j;
        let work_total_j = old.work_j + work_j;
        let heat_w = -step.dissipated_j / dt_s;
        if !(work_j.is_finite()
            && storage_j.is_finite()
            && closure_j.is_finite()
            && work_total_j.is_finite()
            && heat_w.is_finite()
            && step.dissipated_j >= 0.0
            && step.stored_change_j.is_finite()
            && step.closure_residual_j.is_finite()
            && closure_j.abs() <= tolerance_j
            && (storage_j - step.stored_change_j).abs() <= tolerance_j)
        {
            return Err(refused("sliding_interface_work_balance"));
        }
        let scale = step.dissipated_j
            + work_j.abs()
            + before_j.abs()
            + after_j.abs()
            + old.work_j.abs()
            + work_total_j.abs();
        let receipt = Receipt {
            stored_change_j: storage_j - work_j,
            dissipated_j: step.dissipated_j,
            // What the law's own closure leaves of `storage + heat - work`.
            numerical_loss_j: -closure_j,
            entropy_generated_j_per_k: step.dissipated_j / kelvin,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE,
            residual: closure_j.abs(),
            tolerance: tolerance_j,
        };
        let state = SlidingState {
            tangential: step.state,
            reversible_j: after_j,
            work_j: work_total_j,
        };
        Ok((state, heat_w, step.traction_n, receipt))
    }
}

impl<L, F> Component for SlidingInterface<L, F>
where
    L: TangentialLaw,
    F: Fn(f64) -> Result<L, TriboError>,
{
    type State = SlidingState<L::State>;
    type Error = TriboError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &Self::State) -> Result<f64, TriboError> {
        finite(state.reversible_j - state.work_j, "stored_energy_j")
    }

    fn propose(
        &self,
        old: &Self::State,
        dt_s: f64,
        inputs: &[f64],
        new: &mut Self::State,
        outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<Proposal, TriboError> {
        let [heat_w, traction_n @ ..] = outputs else {
            return Err(refused("sliding_interface_outputs"));
        };
        let (state, heat, traction, receipt) = self.step(old, dt_s, inputs)?;
        if traction_n.len() != traction.len() {
            return Err(refused("sliding_interface_outputs"));
        }
        *new = state;
        *heat_w = heat;
        traction_n.copy_from_slice(&traction);
        Ok(Proposal::Step(receipt))
    }
}
