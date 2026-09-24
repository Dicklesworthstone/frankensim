//! Eliminate interior passive elements into the EXISTING one-port load solver.
//!
//! Shunt: a_eq = sum(a_i/Z_i)/sum(1/Z_i), Z_eq = 1/sum(1/Z_i).
//! Series: a_eq = a_left-a_right, Z_eq = Z_left+Z_right.
//! These are the parallel/series wave-adaptor relations (J. O. Smith, PASP,
//! "Adaptors for Wave Digital Elements"). No delay or extra integrator is used.

use super::{NetworkNode, NetworkSegment, NodeFrame, WaveguideError, WaveguideNetwork};

pub(super) fn load_impedance(
    kind: NetworkNode, adjacent: &[usize], segments: &[NetworkSegment],
) -> Result<f64, WaveguideError> {
    // Node degree and section indices have already been admitted by the owner.
    let first = segments[adjacent[0] / 2].impedance_pa_s_m3;
    let z = match kind {
        NetworkNode::Series { .. } => first + segments[adjacent[1] / 2].impedance_pa_s_m3,
        NetworkNode::Shunt { .. } | NetworkNode::ShuntAdmittance { .. } => {
            let min_z = adjacent.iter().map(|&p| segments[p / 2].impedance_pa_s_m3)
                .fold(f64::INFINITY, f64::min);
            let sum = adjacent.iter().fold(0.0, |s, &p| s + min_z / segments[p / 2].impedance_pa_s_m3);
            min_z / sum
        },
        _ => first,
    };
    if !z.is_finite() || z <= 0.0 {
        return Err(WaveguideError("interior equivalent impedance is not representable"));
    }
    Ok(z)
}

impl WaveguideNetwork {
    /// Stage one interior load and all its departing waves. Return only the
    /// interface work defect; never relabel load storage as junction dissipation.
    pub(super) fn preview_interior(&mut self, node: usize) -> Result<f64, WaveguideError> {
        let range = self.offsets[node]..self.offsets[node + 1];
        let first = self.ports[range.start];
        let series = matches!(self.nodes[node], NetworkNode::Series { .. });
        let incident = if series {
            self.arriving(first) - self.arriving(self.ports[range.start + 1])
        } else {
            range.clone().fold(0.0, |sum, i| {
                let port = self.ports[i];
                sum + self.weights[port] * self.arriving(port)
            })
        };
        let trial = self.loads[node].as_ref().expect("admitted interior load").preview_step(incident)?;
        let load = trial.port();
        let mut observation = NodeFrame {
            pressure_pa: load.pressure_pa,
            load_pressure_pa: load.pressure_pa,
            load_flow_m3_s: load.flow_m3_s,
            stored_energy_j: load.stored_energy_j,
            storage_change_j: load.storage_change_j,
            absorbed_energy_j: load.dissipated_energy_j,
            ..NodeFrame::default()
        };
        let mut work = 0.0;
        for index in range {
            let port = self.ports[index];
            let a = self.arriving(port);
            let z = self.lines[port / 2].spec.impedance_pa_s_m3;
            let b = if series {
                // The two section flows into this element are equal/opposite.
                let sign = if port == first { 1.0 } else { -1.0 };
                a - sign * (z * load.flow_m3_s)
            } else { load.pressure_pa - a };
            let pressure = a + b;
            let flow = (a - b) / z;
            if ![b, pressure, flow].iter().all(|x| x.is_finite()) {
                return Err(WaveguideError("interior scattering left the finite set"));
            }
            if series {
                if port == first { observation.pressure_pa = pressure; }
                else { observation.series_other_pressure_pa = Some(pressure); }
            }
            observation.net_flow_into_node_m3_s += flow;
            work += pressure * (flow * self.time_step_s);
            self.departing[port] = b;
        }
        let residual = work - load.supplied_work_j;
        if ![observation.net_flow_into_node_m3_s, residual].iter().all(|x| x.is_finite()) {
            return Err(WaveguideError("interior interface work left the finite set"));
        }
        self.load_candidate[node] = trial;
        self.candidate[node] = observation;
        Ok(residual)
    }
}
