//! An equivalent-circuit battery whose two resistances heat what their heat
//! ports are joined to.
//!
//! [`EquivalentCircuitBattery`] puts `fs_phs::circuit::battery` behind
//! `fs_ports::component::Component`: an OCV/SOC curve owning the chemical
//! energy, a series resistance, one polarization RC pair and a bus capacitor.
//! Nothing is integrated here: a step is the owner's own
//! `PreparedBatterySource::step_controlled` (its componentwise discrete
//! gradient, exact across OCV knots), on the battery built with the
//! resistances of that step.
//!
//! Ports, in this order:
//!
//! - `terminal`: electrical, [`PortRole::FlowIn`]. It receives the current
//!   entering the bus \[A\] (positive charges the battery, negative draws
//!   from it) and returns the owner's discrete bus voltage \[V\], the mean
//!   of the bus capacitor's voltages at the two ends of the step.
//! - `heat_series`, then `heat_polarization`: heat, [`PortRole::EffortIn`].
//!   Each receives the temperature \[K\] where that resistance's heat leaves
//!   and returns minus that heat over the step divided by the step \[W\].
//!
//! State: `[remaining chemical charge, polarization charge, bus charge]`
//! \[C\], the owner's own coordinates. Its stored energy is the owner's
//! Hamiltonian, which no resistance changes.
//!
//! Heat: the owner's own two accounts, `R_s I_cell^2 dt` with the cell
//! current from the step's discrete efforts and `V_p^2 / R_p dt` with the
//! polarization voltage at the midpoint charge. Supplied work is the owner's
//! terminal work. Their difference from the stored change is booked as
//! numerical loss and refused beyond the rounding of the step's arithmetic
//! plus the declared energy tolerance, as is any difference between these
//! accounts and the owner's native step record.
//!
//! Resistance: a resistance with no law keeps the value in the parameters.
//! One following a [`ResistanceLaw`] is given `law.ohms(T)` at the
//! temperature `T` its heat port receives in that step. With a backward-Euler
//! thermal timeline returning end-of-step temperatures, that is the endpoint
//! material resistance of `fs_couple::circuit_enthalpy::battery`. A law with
//! no answer refuses the step; it is never clamped. A law that follows the
//! operating point is refused when it is given.
//!
//! Class: [`Class::Smooth`]. The chemical energy is the integral of a
//! continuous piecewise-linear OCV, so it is continuously differentiable and
//! so is its discrete gradient; no tangent is claimed, so a kernel
//! differences the port map.
//!
//! Not claimed: OCV that depends on temperature, reversible entropic heat,
//! capacitances that drift, ageing or protection events. The SOC is bounded:
//! a step that would leave the curve's support is the owner's refusal.

use core::fmt;
use core::num::NonZeroUsize;
use std::sync::Arc;

use fs_exec::Cx;
use fs_exec::solver::snapshot_v2::StateIdentityCharterV2;
use fs_phs::circuit::battery::{BatteryOcvCurve, BatteryParameters, BatterySource};
use fs_phs::{PhsError, PreparedStepError, SolveBudget, arithmetic_rounding_floor};
use fs_ports::component::{
    Class, Component, ContractViolation, PortLayout, PortRole, PortSpec, Proposal, Receipt,
};
use fs_ports::{CoupleError, PortKind, StableId};

use crate::circuit::ResistanceLaw;

/// One of the battery's two loss channels.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum BatteryLoss {
    /// The series resistance, carrying the cell current.
    Series,
    /// The polarization resistance, across the polarization capacitor.
    Polarization,
}

impl BatteryLoss {
    const fn index(self) -> usize {
        match self {
            Self::Series => 0,
            Self::Polarization => 1,
        }
    }
}

/// Refusal of the battery owner or of its Component boundary.
#[derive(Debug)]
pub enum BatteryComponentError {
    /// The owner refused the curve, the parameters or a state.
    Model(PhsError),
    /// The owner refused the step, unchanged.
    Step(PreparedStepError),
    /// A resistance law had no answer, or answered a resistance that is not
    /// positive and finite.
    Law {
        /// The channel whose law it is.
        loss: BatteryLoss,
        /// The law's reason.
        reason: String,
    },
    /// Invalid or repeated port metadata.
    Port(ContractViolation),
    /// An invalid port identifier.
    Identifier(CoupleError),
    /// A shape, a value or a construction that is refused.
    Invalid {
        /// What was refused.
        what: &'static str,
    },
    /// Work, heat and storage do not close within the declared tolerance and
    /// rounding, or the owner's two accounts of the step disagree.
    Balance {
        /// The largest defect \[J\].
        defect_j: f64,
        /// The allowance \[J\].
        tolerance_j: f64,
    },
    /// Cancellation was requested.
    Cancelled,
}

impl fmt::Display for BatteryComponentError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::Model(error) => write!(f, "battery owner refused: {error:?}"),
            Self::Step(error) => write!(f, "battery step refused: {error:?}"),
            Self::Law { loss, reason } => {
                write!(
                    f,
                    "resistance law of the {loss:?} loss has no answer: {reason}"
                )
            }
            Self::Port(error) => write!(f, "battery port refused: {error:?}"),
            Self::Identifier(error) => write!(f, "battery port identifier refused: {error:?}"),
            Self::Invalid { what } => write!(f, "battery component refused: {what}"),
            Self::Balance {
                defect_j,
                tolerance_j,
            } => write!(
                f,
                "battery balance defect {defect_j} J exceeds {tolerance_j} J"
            ),
            Self::Cancelled => write!(f, "cancelled"),
        }
    }
}

impl std::error::Error for BatteryComponentError {}

fn invalid(what: &'static str) -> BatteryComponentError {
    BatteryComponentError::Invalid { what }
}

fn step_error(error: PreparedStepError) -> BatteryComponentError {
    match error {
        PreparedStepError::Cancelled => BatteryComponentError::Cancelled,
        other @ PreparedStepError::Solver(_) => BatteryComponentError::Step(other),
    }
}

/// An equivalent-circuit battery of `fs_phs::circuit::battery` as one
/// Component.
///
/// See the module documentation for the ports, the state, the heat and the
/// resistances.
pub struct EquivalentCircuitBattery {
    capacity_c: f64,
    knots: Vec<(f64, f64)>,
    curve: BatteryOcvCurve,
    parameters: BatteryParameters,
    laws: [Option<Arc<dyn ResistanceLaw>>; 2],
    solve: SolveBudget,
    backtracking_halvings: u8,
    /// The owner at the parameters' resistances, for states and energies,
    /// neither of which depends on a resistance.
    reference: BatterySource,
    ports: Vec<PortSpec>,
    energy_tolerance_j: f64,
}

impl fmt::Debug for EquivalentCircuitBattery {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("EquivalentCircuitBattery")
            .field("capacity_c", &self.capacity_c)
            .field("knots", &self.knots)
            .field("parameters", &self.parameters)
            .field("laws", &self.laws)
            .field("solve", &self.solve)
            .field("backtracking_halvings", &self.backtracking_halvings)
            .field("energy_tolerance_j", &self.energy_tolerance_j)
            .finish_non_exhaustive()
    }
}

impl EquivalentCircuitBattery {
    /// The battery whose OCV is the piecewise-linear curve through
    /// `soc_voltage` (state of charge, volts) over `capacity_c` \[C\], with
    /// `parameters`, solved by the owner within `solve` and at most
    /// `backtracking_halvings` step halvings of its Newton update.
    ///
    /// `energy_tolerance_j` \[J\] is the largest balance defect a step may
    /// leave beyond the rounding of its own arithmetic: the owner's solve
    /// residual times its efforts, which this boundary cannot bound by
    /// itself. A defect within it is booked as numerical loss.
    ///
    /// # Errors
    /// The owner's refusal of the curve or the parameters
    /// ([`BatteryComponentError::Model`]); a tolerance that is not positive
    /// and finite.
    pub fn new(
        capacity_c: f64,
        soc_voltage: &[(f64, f64)],
        parameters: BatteryParameters,
        solve: SolveBudget,
        backtracking_halvings: u8,
        energy_tolerance_j: f64,
    ) -> Result<Self, BatteryComponentError> {
        if !(energy_tolerance_j.is_finite() && energy_tolerance_j > 0.0) {
            return Err(invalid("energy tolerance must be positive and finite"));
        }
        let curve =
            BatteryOcvCurve::new(capacity_c, soc_voltage).map_err(BatteryComponentError::Model)?;
        let reference =
            BatterySource::new(curve.clone(), parameters).map_err(BatteryComponentError::Model)?;
        let mut ports = Vec::with_capacity(3);
        ports.push(
            PortSpec::power(
                StableId::new("terminal").map_err(BatteryComponentError::Identifier)?,
                PortKind::ElectricalVoltageCurrent,
                PortRole::FlowIn,
                NonZeroUsize::MIN,
            )
            .map_err(BatteryComponentError::Port)?,
        );
        for name in ["heat_series", "heat_polarization"] {
            let id = StableId::new(name).map_err(BatteryComponentError::Identifier)?;
            ports
                .push(PortSpec::heat(id, PortRole::EffortIn).map_err(BatteryComponentError::Port)?);
        }
        PortLayout::new(&ports).map_err(BatteryComponentError::Port)?;
        Ok(Self {
            capacity_c,
            knots: soc_voltage.to_vec(),
            curve,
            parameters,
            laws: [None, None],
            solve,
            backtracking_halvings,
            reference,
            ports,
            energy_tolerance_j,
        })
    }

    /// The same battery with the resistance of `loss` following `law` at the
    /// temperature its heat port receives, in every proposal.
    ///
    /// # Errors
    /// That resistance already follows a law, or the law follows the
    /// operating point
    /// ([`ResistanceLaw::follows_operating_point`]), which the battery does
    /// not solve: it asks its laws at temperature alone.
    pub fn following(
        mut self,
        loss: BatteryLoss,
        law: Arc<dyn ResistanceLaw>,
    ) -> Result<Self, BatteryComponentError> {
        if law.follows_operating_point() {
            return Err(invalid(
                "a battery resistance follows temperature, not the operating point",
            ));
        }
        if self.laws[loss.index()].replace(law).is_some() {
            return Err(invalid("a resistance follows one law only"));
        }
        Ok(self)
    }

    /// The parameters the battery was built with.
    #[must_use]
    pub const fn parameters(&self) -> BatteryParameters {
        self.parameters
    }

    /// A state at state of charge `soc`, polarization voltage
    /// `polarization_v` \[V\] and bus voltage `bus_v` \[V\].
    ///
    /// # Errors
    /// The owner's refusal: a value that is not finite or a state of charge
    /// outside the curve.
    pub fn state(
        &self,
        soc: f64,
        polarization_v: f64,
        bus_v: f64,
    ) -> Result<[f64; 3], BatteryComponentError> {
        self.reference
            .state(soc, polarization_v, bus_v)
            .map_err(BatteryComponentError::Model)
    }

    fn energy(&self, state: &[f64; 3]) -> Result<f64, BatteryComponentError> {
        let energy = self.reference.system().hamiltonian(state);
        if energy.is_finite() {
            Ok(energy)
        } else {
            Err(invalid(
                "the battery's state is outside its curve or its energy is not finite",
            ))
        }
    }

    /// The parameters of a step whose heat ports received `kelvin`.
    fn parameters_at(&self, kelvin: [f64; 2]) -> Result<BatteryParameters, BatteryComponentError> {
        let mut ohms = [
            self.parameters.series_resistance_ohm,
            self.parameters.polarization_resistance_ohm,
        ];
        for (loss, (law, temperature)) in [BatteryLoss::Series, BatteryLoss::Polarization]
            .into_iter()
            .zip(self.laws.iter().zip(kelvin))
        {
            let Some(law) = law else { continue };
            let answer = law
                .ohms(temperature)
                .map_err(|reason| BatteryComponentError::Law { loss, reason })?;
            if !(answer.is_finite() && answer > 0.0) {
                return Err(BatteryComponentError::Law {
                    loss,
                    reason: format!("{answer} ohm at {temperature} K is not positive and finite"),
                });
            }
            ohms[loss.index()] = answer;
        }
        Ok(BatteryParameters {
            series_resistance_ohm: ohms[0],
            polarization_resistance_ohm: ohms[1],
            ..self.parameters
        })
    }
}

impl Component for EquivalentCircuitBattery {
    type State = [f64; 3];
    type Error = BatteryComponentError;

    const STATE_CHARTER: Option<StateIdentityCharterV2> = Some(StateIdentityCharterV2 {
        owner: "fs-system::battery::EquivalentCircuitBattery",
        state_family: "fs_phs::circuit::battery charges",
        schema_grammar: "[f64; 3]: remaining chemical charge, polarization charge, bus charge, in coulombs",
        codec_grammar: "one word per coordinate, the f64 bits, in schema order",
        codec_version: 1,
    });

    fn encode_state(&self, state: &[f64; 3], words: &mut Vec<u64>) -> bool {
        words.extend(state.iter().map(|value| value.to_bits()));
        true
    }

    fn decode_state(&self, words: &[u64]) -> Option<[f64; 3]> {
        let &[a, b, c] = words else { return None };
        let state = [a, b, c].map(f64::from_bits);
        state.iter().all(|value| value.is_finite()).then_some(state)
    }

    fn identity_words(&self, words: &mut Vec<u64>) {
        words.push(self.capacity_c.to_bits());
        words.push(self.knots.len() as u64);
        for &(soc, volts) in &self.knots {
            words.extend([soc.to_bits(), volts.to_bits()]);
        }
        let p = self.parameters;
        words.extend([
            p.series_resistance_ohm.to_bits(),
            p.polarization_resistance_ohm.to_bits(),
            p.polarization_capacitance_f.to_bits(),
            p.bus_capacitance_f.to_bits(),
        ]);
        words.extend([
            self.solve.relative.to_bits(),
            self.solve.absolute.to_bits(),
            self.solve.max_iterations as u64,
            u64::from(self.backtracking_halvings),
        ]);
        // Each law's words are counted, so that one law's words cannot read
        // as another's.
        for law in &self.laws {
            match law {
                None => words.push(0),
                Some(law) => {
                    words.push(1);
                    let count = words.len();
                    words.push(0);
                    law.identity_words(words);
                    words[count] = (words.len() - count - 1) as u64;
                }
            }
        }
        words.push(self.energy_tolerance_j.to_bits());
    }

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::Smooth
    }

    fn stored_energy_j(&self, state: &[f64; 3]) -> Result<f64, BatteryComponentError> {
        self.energy(state)
    }

    fn propose(
        &self,
        old: &[f64; 3],
        dt_s: f64,
        inputs: &[f64],
        new: &mut [f64; 3],
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, BatteryComponentError> {
        cx.checkpoint()
            .map_err(|_| BatteryComponentError::Cancelled)?;
        if !(dt_s.is_finite() && dt_s > 0.0)
            || inputs.len() != 3
            || outputs.len() != 3
            || inputs.iter().any(|value| !value.is_finite())
            || inputs[1..].iter().any(|&kelvin| kelvin <= 0.0)
        {
            return Err(invalid(
                "positive duration, complete finite ports and positive temperatures required",
            ));
        }
        let (current_a, kelvin) = (inputs[0], [inputs[1], inputs[2]]);
        let parameters = self.parameters_at(kelvin)?;
        let source = BatterySource::new(self.curve.clone(), parameters)
            .map_err(BatteryComponentError::Model)?;
        let mut solver = source
            .prepare_controlled(self.solve, self.backtracking_halvings, || {
                cx.is_cancel_requested()
            })
            .map_err(step_error)?;
        let mut next = [0.0; 3];
        let step = solver
            .step_controlled(old, current_a, dt_s, &mut next, || cx.is_cancel_requested())
            .map_err(step_error)?;
        let account = step.battery;
        let heats = [account.series_heat_j, account.polarization_heat_j];
        if heats.iter().any(|heat| !(heat.is_finite() && *heat >= 0.0)) {
            return Err(invalid("a resistance's heat is not finite and nonnegative"));
        }
        let rates = heats.map(|heat| heat / dt_s);
        if rates
            .iter()
            .zip(heats)
            .any(|(rate, heat)| !rate.is_finite() || (heat != 0.0 && *rate == 0.0))
        {
            return Err(invalid("a resistance's heat rate is not representable"));
        }

        let stored_change_j = account.storage_change_j();
        let heat_j = account.heat_j();
        let supplied_j = account.terminal_work_j;
        let numerical_loss_j = -account.balance_residual_j();
        let energy_old = self.energy(old)?;
        let energy_new = self.energy(&next)?;
        let magnitude = energy_old.abs()
            + energy_new.abs()
            + supplied_j.abs()
            + heat_j
            + step.native.delta_h.abs()
            + step.native.supplied.abs()
            + step.native.dissipated.abs();
        let rounding = arithmetic_rounding_floor(64, magnitude);
        let allowance = rounding + self.energy_tolerance_j;
        // The owner's component accounts against its native step record.
        let defect_j = [
            numerical_loss_j,
            stored_change_j - step.native.delta_h,
            heat_j - step.native.dissipated,
            supplied_j - step.native.supplied,
        ]
        .into_iter()
        .map(f64::abs)
        .fold(0.0, f64::max);
        let entropy = heats[0] / kelvin[0] + heats[1] / kelvin[1];
        if ![
            rounding,
            defect_j,
            stored_change_j,
            numerical_loss_j,
            entropy,
        ]
        .iter()
        .all(|value| value.is_finite())
        {
            return Err(invalid("the battery's energy account is not finite"));
        }
        if defect_j > allowance {
            return Err(BatteryComponentError::Balance {
                defect_j,
                tolerance_j: allowance,
            });
        }
        let certificate = step
            .native
            .solver_certificate
            .ok_or_else(|| invalid("the owner's step carries no solve certificate"))?;
        cx.checkpoint()
            .map_err(|_| BatteryComponentError::Cancelled)?;
        outputs[0] = account.discrete_terminal_voltage_v;
        outputs[1] = -rates[0];
        outputs[2] = -rates[1];
        *new = next;
        Ok(Proposal::Step(Receipt {
            stored_change_j,
            dissipated_j: heat_j,
            numerical_loss_j,
            entropy_generated_j_per_k: entropy,
            energy_tolerance_j: rounding,
            residual: certificate.residual,
            tolerance: certificate.effective_budget,
        }))
    }
}
