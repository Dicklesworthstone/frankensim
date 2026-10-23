//! Fixed nonlinear terminal losses sharing the original wave/load transaction.
use super::{NetworkNode, WaveguideNetwork};
use crate::{quadratic_loss::QuadraticLossSpec, waveguide::WaveguideError};

impl WaveguideNetwork {
    /// Attach `dp = K U|U|` in series with original RLC/relaxation terminals.
    /// Every law uses the same solved current flow as the retained load states;
    /// no previous-flow resistance, artificial lag or discarded stored energy.
    /// Loss coefficients occupy the node payload admitted by `new`.
    ///
    /// This local law does not support a mutually coupled MIMO terminal group.
    /// Such a group requires a simultaneous nonlinear multivariate solve.
    ///
    /// # Errors
    /// Advanced graph, duplicate/repeated attachment, nonterminal, memoryless
    /// reflection or MIMO port. Admission is atomic and never changes history.
    pub fn install_terminal_quadratic_losses(&mut self,
        assignments: &[(usize, QuadraticLossSpec)]) -> Result<(), WaveguideError> {
        if self.accepted_steps != 0 || assignments.len() > self.nodes.len() {
            return Err(WaveguideError("quadratic terminal losses require a cold bounded graph"));
        }
        for (index, &(node, _)) in assignments.iter().enumerate() {
            if !matches!(self.nodes.get(node), Some(NetworkNode::Impedance { .. } | NetworkNode::Relaxation { .. }))
                || assignments[..index].iter().any(|(previous, _)| *previous == node)
                || self.quadratic_losses[node].is_some()
                || self.mimo.as_ref().is_some_and(|m| m.index(node).is_some()) {
                return Err(WaveguideError("quadratic loss needs a distinct original local impedance terminal"));
            }
        }
        for &(node, loss) in assignments { self.quadratic_losses[node] = Some(loss); }
        Ok(())
    }

    /// Original explicit nonlinear law, separate from the linear node load.
    #[must_use]
    pub fn terminal_quadratic_loss(&self, node: usize) -> Option<QuadraticLossSpec> {
        self.quadratic_losses.get(node).copied().flatten()
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::{impedance::SeriesImpedanceSpec,
        relaxation::{ComplianceRelaxationTerm, RelaxationImpedanceSpec, RelaxationTerm},
        waveguide::network::NetworkSegment};

    fn graph(k: Option<f64>, limit: f64) -> WaveguideNetwork {
        let load = RelaxationImpedanceSpec::new(SeriesImpedanceSpec {
            resistance_pa_s_m3: 200.0, inertance_pa_s2_m3: 0.2,
            compliance_m3_pa: Some(2e-5),
        }, &[RelaxationTerm { resistance_pa_s_m3: 300.0, rate_per_s: 400.0 }]).unwrap()
            .with_compliance_terms(&[ComplianceRelaxationTerm {
                resistance_pa_s_m3: 100.0, rate_per_s: 100.0,
            }]).unwrap();
        let mut graph = WaveguideNetwork::new(&[NetworkNode::Inlet, NetworkNode::Relaxation { load }],
            &[NetworkSegment { nodes: [0, 1], one_way_samples: 3, impedance_pa_s_m3: 1000.0 }],
            1.0/48000.0, 100000).unwrap();
        if let Some(k) = k {
            graph.install_terminal_quadratic_losses(&[(1, QuadraticLossSpec::new(k, limit).unwrap())]).unwrap();
        }
        graph
    }

    #[test]
    fn g0_quadratic_terminal_preserves_all_reactive_storage_and_closes_work() {
        let mut positive = graph(Some(1e7), 1.0);
        let mut negative = graph(Some(1e7), 1.0);
        let mut accumulated_loss = 0.0;
        for index in 0..4096 {
            let input = if index < 2000 { 3.0*(0.027*index as f64).sin() } else { 0.0 };
            let a = positive.step(input).unwrap();
            let b = negative.step(-input).unwrap();
            let node = positive.node_frame(1).unwrap();
            let opposite = negative.node_frame(1).unwrap();
            let flow = node.net_flow_into_node_m3_s;
            let scale = a.stored_energy_j + a.terminal_loss_j + a.inlet_work_j.abs();
            assert!(a.balance_residual_j().abs() < 2e-13*scale.max(1e-20));
            assert!((node.quadratic_pressure_drop_pa - 1e7*flow*flow.abs()).abs()
                < 2e-12*node.quadratic_pressure_drop_pa.abs().max(1e-20));
            assert!(node.quadratic_loss_j >= 0.0);
            assert_eq!(a.stored_energy_j, b.stored_energy_j);
            assert_eq!(node.pressure_pa, -opposite.pressure_pa);
            accumulated_loss += node.quadratic_loss_j;
        }
        assert!(accumulated_loss > 0.0);
        assert!(positive.terminal_relaxation_flows(1).unwrap()[0] != 0.0);
        assert!(positive.load_compliance_pressures(1).unwrap()[0] != 0.0);
    }

    #[test]
    fn g3_zero_coefficient_is_bitwise_linear_and_refusal_preserves_histories() {
        let mut linear = graph(None, 1.0);
        let mut zero = graph(Some(0.0), 1.0);
        for index in 0..128 {
            let input = (index as f64*0.1).sin();
            assert_eq!(linear.step(input).unwrap(), zero.step(input).unwrap());
            assert_eq!(linear.node_frame(1), zero.node_frame(1));
        }
        let mut limited = graph(Some(1e7), 1e-8);
        for _ in 0..3 { limited.step(2.0).unwrap(); }
        let state = limited.terminal_state(1);
        let steps = limited.accepted_steps();
        let storage = limited.stored_energy_j();
        let observation = limited.node_frame(1).copied();
        assert!(limited.step(2.0).is_err());
        assert_eq!(limited.accepted_steps(), steps);
        assert_eq!(limited.terminal_state(1), state);
        assert_eq!(limited.stored_energy_j(), storage);
        assert_eq!(limited.node_frame(1).copied(), observation);
    }
}
