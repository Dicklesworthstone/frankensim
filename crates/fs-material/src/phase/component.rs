//! A uniform body on the original equilibrium enthalpy chart, on the
//! composition contract.
//!
//! The chart owns `h -> (T, phase)` and the entropy-conjugate temperature of
//! a finite enthalpy step. This adapter gives them one heat port. It adds no
//! heat capacity, no apparent latent heat, no transport law and no geometry.
//!
//! Port `heat` receives the heat rate [W] into the body and returns
//! `theta = delta_h / delta_s` [K] of the step that heat makes: the
//! temperature at which the body's entropy change is exactly `Q / theta`.
//! Over an interval of constant heat capacity that is the logarithmic mean
//! of the two end temperatures, and on a latent plateau it is the transition
//! temperature. A conductance that is given this temperature therefore
//! produces `G (theta_a - theta_b)^2 / (theta_a theta_b) dt` of entropy and
//! the body none, whatever the size of the step.
//!
//! The body has a fixed mass and stores enthalpy with negligible `p V` work,
//! as its chart does. Its volume, its surface and the law of any conductance
//! that touches it belong to whoever connects it; a phase change moves
//! nothing mechanical here.

use fs_exec::Cx;
use fs_ports::StableId;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use super::{
    EnthalpyEntropyMean, EquilibriumEnthalpyPhaseCurve, EquilibriumPhaseState, PhaseStateError,
};

/// A fixed mass of the material of one enthalpy chart, at one enthalpy.
///
/// Its [`Component::State`] is the chart's own [`EquilibriumPhaseState`].
#[derive(Clone, Debug)]
pub struct UniformEnthalpyBody {
    curve: EquilibriumEnthalpyPhaseCurve,
    mass_kg: f64,
    ports: [PortSpec; 1],
}

fn refused(what: &'static str) -> PhaseStateError {
    PhaseStateError::InvalidIsobaricHeatStep { what }
}

fn poll(cx: &Cx<'_>) -> Result<(), PhaseStateError> {
    cx.checkpoint()
        .map_err(|_| PhaseStateError::EntropyEvaluationCancelled)
}

impl UniformEnthalpyBody {
    /// `mass_kg` of the material of `curve`. Pass a state of the same curve
    /// as the initial state when adding the component.
    ///
    /// # Errors
    /// A mass that is not positive and finite.
    pub fn new(
        curve: EquilibriumEnthalpyPhaseCurve,
        mass_kg: f64,
    ) -> Result<Self, PhaseStateError> {
        if !(mass_kg.is_finite() && mass_kg > 0.0) {
            return Err(refused(
                "a uniform enthalpy body requires a positive finite mass",
            ));
        }
        Ok(Self {
            curve,
            mass_kg,
            ports: [PortSpec::heat(
                StableId::new("heat").expect("static port name"),
                PortRole::FlowIn,
            )
            .expect("heat input port")],
        })
    }

    /// The chart every state of this body belongs to.
    #[must_use]
    pub const fn curve(&self) -> &EquilibriumEnthalpyPhaseCurve {
        &self.curve
    }

    /// Mass of the body [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> f64 {
        self.mass_kg
    }

    /// What the heat rate `inputs[0]` makes of `old` in `dt_s`: the state
    /// at the new enthalpy, the chart's mean over the step, and the heat
    /// [J]. Both come from the chart's own calls at the two enthalpies.
    fn step(
        &self,
        old: &EquilibriumPhaseState,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(EquilibriumPhaseState, EnthalpyEntropyMean, f64), PhaseStateError> {
        poll(cx)?;
        let [heat_w] = inputs else {
            return Err(refused("a uniform enthalpy body requires one heat rate"));
        };
        if !heat_w.is_finite() {
            return Err(PhaseStateError::NonFiniteSpecificEnergyIncrement);
        }
        if !(dt_s.is_finite() && dt_s > 0.0) {
            return Err(refused(
                "a uniform enthalpy body requires a positive finite timestep",
            ));
        }
        if old.phase_curve_identity() != self.curve.identity() {
            return Err(PhaseStateError::StateCurveMismatch);
        }
        let heat_j = dt_s * heat_w;
        let initial = old.specific_enthalpy_j_kg();
        // Heat too small to move the enthalpy leaves it where it is; the
        // receipt books what rounding kept and what it lost.
        let reached = initial + heat_j / self.mass_kg;
        let endpoint = self.curve.state_at_specific_enthalpy(reached)?;
        let mean = self.curve.thermodynamic_mean_temperature(
            initial,
            reached,
            self.curve.knots().len() - 1,
            || cx.is_cancel_requested(),
        )?;
        Ok((endpoint, mean, heat_j))
    }

    /// Slope of the returned temperature in the heat rate [K/W] at the step
    /// [`Self::step`] takes: the chart's derivative of its mean in the final
    /// enthalpy, with its one-sided convention at a corner of the chart and
    /// zero across a plateau.
    fn slope(
        &self,
        old: &EquilibriumPhaseState,
        dt_s: f64,
        inputs: &[f64],
        cx: &Cx<'_>,
    ) -> Result<f64, PhaseStateError> {
        let (_, mean, _) = self.step(old, dt_s, inputs, cx)?;
        let slope = mean.d_temperature_d_final_enthalpy_kg_k_j * (dt_s / self.mass_kg);
        if slope.is_finite() {
            Ok(slope)
        } else {
            Err(PhaseStateError::UnrepresentableTemperatureDerivative)
        }
    }

    /// `slope * seed` into `slot`: the tangent in either direction, the map
    /// being one by one.
    fn tangent(
        &self,
        old: &EquilibriumPhaseState,
        dt_s: f64,
        inputs: &[f64],
        seed: &[f64],
        slot: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), PhaseStateError> {
        let ([seed], [slot]) = (seed, slot) else {
            return Err(refused(
                "a uniform enthalpy body has one tangent value a side",
            ));
        };
        let value = self.slope(old, dt_s, inputs, cx)? * seed;
        if !value.is_finite() {
            return Err(PhaseStateError::UnrepresentableTemperatureDerivative);
        }
        poll(cx)?;
        *slot = value;
        Ok(())
    }
}

impl Component for UniformEnthalpyBody {
    type State = EquilibriumPhaseState;
    type Error = PhaseStateError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        // One input and one output. The mean never falls as the heat rises,
        // and it is flat across a latent plateau.
        Class::Monotone1d
    }

    fn stored_energy_j(&self, state: &EquilibriumPhaseState) -> Result<f64, PhaseStateError> {
        if state.phase_curve_identity() != self.curve.identity() {
            return Err(PhaseStateError::StateCurveMismatch);
        }
        // Enthalpy in the chart's own datum.
        Ok(self.mass_kg * state.specific_enthalpy_j_kg())
    }

    fn propose(
        &self,
        old: &EquilibriumPhaseState,
        dt_s: f64,
        inputs: &[f64],
        new: &mut EquilibriumPhaseState,
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, PhaseStateError> {
        let [returned] = outputs else {
            return Err(refused(
                "a uniform enthalpy body requires one scalar port output",
            ));
        };
        let (endpoint, mean, heat_j) = self.step(old, dt_s, inputs, cx)?;
        let (initial, reached) = (
            old.specific_enthalpy_j_kg(),
            endpoint.specific_enthalpy_j_kg(),
        );
        // The represented change, as the chart's other owners form it.
        let stored_change_j = self.mass_kg * (reached - initial);
        let scale = self.mass_kg * (initial.abs() + reached.abs()) + heat_j.abs();
        let receipt = Receipt {
            stored_change_j,
            dissipated_j: 0.0,
            // What forming the new enthalpy rounded away of the heat.
            numerical_loss_j: heat_j - stored_change_j,
            // The body's entropy change is the heat over the returned
            // temperature exactly, so it generates none.
            entropy_generated_j_per_k: 0.0,
            energy_tolerance_j: 64.0 * f64::EPSILON * scale + f64::MIN_POSITIVE,
            residual: 0.0,
            tolerance: 0.0,
        };
        if !(stored_change_j.is_finite()
            && receipt.numerical_loss_j.is_finite()
            && receipt.energy_tolerance_j.is_finite())
        {
            return Err(refused(
                "a uniform enthalpy body's energy account is unrepresentable",
            ));
        }
        poll(cx)?;
        *new = endpoint;
        *returned = mean.temperature_k;
        Ok(Proposal::Step(receipt))
    }
}

impl PortTangent for UniformEnthalpyBody {
    fn jvp(
        &self,
        old: &EquilibriumPhaseState,
        dt_s: f64,
        inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), PhaseStateError> {
        self.tangent(old, dt_s, inputs, d_inputs, d_outputs, cx)
    }

    fn vjp(
        &self,
        old: &EquilibriumPhaseState,
        dt_s: f64,
        inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<(), PhaseStateError> {
        self.tangent(old, dt_s, inputs, bar_outputs, bar_inputs, cx)
    }
}
