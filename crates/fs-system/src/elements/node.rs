//! A thermal node of no heat capacity, joined to `n` bodies by conductances,
//! with one port for a heat source at the node.
//!
//! It is [`super::Conductance`] with a source in the middle: two bodies
//! joined through it without a source exchange the heat of the two
//! conductances in series. With a source it is the plane a heat is generated
//! on between bodies (a sliding interface, a thin heater, a glue line): the
//! heat divides by where the bodies' temperatures put the node, not by a
//! prescribed fraction.
//!
//! The node stores nothing, so its temperature is fixed by its balance,
//! `P + sum G_k (T_k - T_n) = 0`:
//! `T_n = (P + sum G_k T_k) / sum G_k`.
//! Each conductance carries heat between `T_k` and `T_n` and generates
//! `G_k (T_k - T_n)^2 / (T_k T_n)` of entropy per second.

use fs_exec::Cx;
use fs_ports::component::{Class, Component, PortRole, PortSpec, PortTangent, Proposal, Receipt};

use super::{ElementError, heat, inputs_finite, positive, roundoff, temperature};

/// A massless thermal node with port `source` and ports `body0`, `body1`,
/// ... in the order of its conductances.
///
/// `source` receives the heat rate [W] put into the node and returns the
/// node's temperature [K]. `body<k>` receives that body's temperature [K]
/// and returns the heat rate [W] that enters the node from it,
/// `G_k (T_k - T_n)`; the body, on the other side of the wire, receives the
/// opposite.
#[derive(Debug, Clone)]
pub struct HeatNode {
    conductances_w_k: Vec<f64>,
    /// Their sum.
    total_w_k: f64,
    ports: Vec<PortSpec>,
}

impl HeatNode {
    /// A node joined to one body per entry of `conductances_w_k` [W/K].
    ///
    /// # Errors
    /// No conductance, or one that is not positive and finite, or a sum
    /// that is not finite.
    pub fn new(conductances_w_k: &[f64]) -> Result<Self, ElementError> {
        if conductances_w_k.is_empty() {
            return Err(ElementError::Parameter {
                name: "conductances_w_k (count)",
                value: 0.0,
            });
        }
        let mut ports = vec![heat("source", PortRole::FlowIn)];
        let mut total_w_k = 0.0;
        for (index, &conductance) in conductances_w_k.iter().enumerate() {
            total_w_k += positive("conductances_w_k", conductance)?;
            ports.push(heat(&format!("body{index}"), PortRole::EffortIn));
        }
        Ok(Self {
            conductances_w_k: conductances_w_k.to_vec(),
            total_w_k: positive("conductances_w_k (sum)", total_w_k)?,
            ports,
        })
    }

    /// The node's temperature [K] for the source rate `inputs[0]` and the
    /// bodies' temperatures after it.
    fn node_k(&self, inputs: &[f64]) -> Result<f64, ElementError> {
        inputs_finite(inputs)?;
        let mut weighted_w = inputs[0];
        for (&conductance, &kelvin) in self.conductances_w_k.iter().zip(&inputs[1..]) {
            weighted_w += conductance * temperature(kelvin)?;
        }
        // A source that takes more heat than the conductances can bring at
        // any positive temperature leaves no temperature for the node.
        temperature(weighted_w / self.total_w_k)
    }
}

impl Component for HeatNode {
    type State = ();
    type Error = ElementError;

    fn ports(&self) -> &[PortSpec] {
        &self.ports
    }

    fn class(&self) -> Class {
        Class::ExactLinear
    }

    fn stored_energy_j(&self, _state: &()) -> Result<f64, ElementError> {
        Ok(0.0)
    }

    fn propose(
        &self,
        _old: &(),
        dt_s: f64,
        inputs: &[f64],
        _new: &mut (),
        outputs: &mut [f64],
        cx: &Cx<'_>,
    ) -> Result<Proposal, ElementError> {
        cx.checkpoint()?;
        let node_k = self.node_k(inputs)?;
        let (mut terms_w, mut entropy_w_k) = (inputs[0].abs(), 0.0);
        for (index, (&conductance, &kelvin)) in
            self.conductances_w_k.iter().zip(&inputs[1..]).enumerate()
        {
            let rate_w = conductance * (kelvin - node_k);
            outputs[1 + index] = rate_w;
            terms_w += conductance * (kelvin + node_k);
            entropy_w_k += rate_w * (kelvin - node_k) / (kelvin * node_k);
        }
        outputs[0] = node_k;
        Ok(Proposal::Step(Receipt {
            stored_change_j: 0.0,
            dissipated_j: 0.0,
            numerical_loss_j: 0.0,
            entropy_generated_j_per_k: dt_s * entropy_w_k,
            // The rates sum to zero by the node's balance, to the rounding
            // of its terms: each rate is a difference of `G_k T_k` and
            // `G_k T_n`, which do not shrink when the rate does (bodies
            // near one temperature).
            energy_tolerance_j: roundoff(dt_s * terms_w),
            residual: 0.0,
            tolerance: 0.0,
        }))
    }
}

impl PortTangent for HeatNode {
    fn jvp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        d_inputs: &[f64],
        d_outputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        let mut weighted = d_inputs[0];
        for (&conductance, &d_kelvin) in self.conductances_w_k.iter().zip(&d_inputs[1..]) {
            weighted += conductance * d_kelvin;
        }
        let d_node = weighted / self.total_w_k;
        d_outputs[0] = d_node;
        for (index, (&conductance, &d_kelvin)) in
            self.conductances_w_k.iter().zip(&d_inputs[1..]).enumerate()
        {
            d_outputs[1 + index] = conductance * (d_kelvin - d_node);
        }
        Ok(())
    }

    fn vjp(
        &self,
        _old: &(),
        _dt_s: f64,
        _inputs: &[f64],
        bar_outputs: &[f64],
        bar_inputs: &mut [f64],
        _cx: &Cx<'_>,
    ) -> Result<(), ElementError> {
        // With w = (bar_node - sum G_k bar_k) / sum G: the source's entry
        // is w and body j's is G_j (w + bar_j).
        let mut weighted = bar_outputs[0];
        for (&conductance, &bar) in self.conductances_w_k.iter().zip(&bar_outputs[1..]) {
            weighted -= conductance * bar;
        }
        let w = weighted / self.total_w_k;
        bar_inputs[0] = w;
        for (index, (&conductance, &bar)) in
            self.conductances_w_k.iter().zip(&bar_outputs[1..]).enumerate()
        {
            bar_inputs[1 + index] = conductance * (w + bar);
        }
        Ok(())
    }
}
