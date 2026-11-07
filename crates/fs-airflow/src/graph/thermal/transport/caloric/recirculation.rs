//! Imposed return streams with original source enthalpy and an analytic Jacobian.
//!
//! At each supply, h_mix = (1-r) h_fresh + sum r_j h_return,j. The native
//! DAG transports that same h; the reduced Newton system acts only on the
//! returned supplies. Passive channel gains are d h_out / d h_in, not frozen
//! temperature attenuation. All trial states and the final outer watt ledger
//! remain subject to the original source support and transport budgets.

use super::super::recirculation::{
    RecirculatedStream, RecirculatedSupply, RecirculationLink, RecirculationReport, ReturnTopology,
    admit_returns, solve_feedback,
};
use super::*;

/// Bounded nonlinear work and temperature accuracy for source-caloric returns.
/// The network's original absolute watt tolerance also gates all mixers.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricRecirculationConfig {
    /// Positive upper bound on each mixture-temperature residual, K.
    pub temperature_tolerance_k: f64,
    /// Maximum accepted residual evaluations, including the initial state.
    /// Must lie in 1..=256. Each Newton step has at most 32 damped trials.
    pub max_iterations: usize,
}

#[derive(Debug)]
pub(super) struct CaloricFeedback {
    links: Vec<RecirculationLink>,
    sources: Vec<usize>,
    fractions: Vec<f64>,
    remaining_exhaust: Vec<f64>,
    minimum_cp: f64,
    fresh_enthalpies: Vec<Option<f64>>,
    config: CaloricRecirculationConfig,
}

struct Trial {
    inlets: Vec<Option<f64>>,
    inlet_enthalpies: Vec<Option<f64>>,
    node_enthalpies: Vec<Option<f64>>,
    branch_enthalpy_gains: Vec<f64>,
    march: TransportMarch,
    mismatch: Vec<f64>,
    temperature_residual: f64,
    heat_residual: f64,
}

impl CaloricTransportNetwork<'_> {
    /// Bind prescribed adiabatic return fractions to this original source law.
    /// Inlet declarations become fresh-makeup temperatures. The return path's
    /// pressure reset, fan work, residence time and hydraulic losses remain
    /// outside the model, as in the native constant-Cp return owner.
    ///
    /// Fractions are receiving-supply MASS fractions. Positive fresh makeup,
    /// exhaust availability, deterministic ordering, native resource caps and
    /// cancellation are checked before execution. No Cp average substitutes
    /// for the source integral, and no source extrapolation is introduced.
    pub fn with_recirculation(
        mut self,
        cx: &Cx<'_>,
        links: Vec<RecirculationLink>,
        config: CaloricRecirculationConfig,
    ) -> Result<Self, TransportError> {
        poll(cx)?;
        positive(
            config.temperature_tolerance_k,
            "recirculation temperature tolerance",
        )?;
        if self.circulation.is_some() {
            return Err(TransportError::InvalidInput(
                "imposed return links cannot be added to a closed source-caloric circulation",
            ));
        }
        if self.feedback.is_some() || !(1..=256).contains(&config.max_iterations) {
            return Err(TransportError::InvalidInput(
                "source-caloric returns already bound or iteration budget outside 1..=256",
            ));
        }
        // Admit actual MASS before any Cp-scaled capacity arithmetic can
        // reject a complete exhaust draw because of reference-Cp roundoff.
        // Share native ownership/caps without constructing affine matrices
        // that the source-caloric Newton solve does not use.
        let ReturnTopology {
            links,
            sources,
            fractions,
            remaining_exhaust,
        } = admit_returns(cx, links, &self.external_mass)?;
        if links.is_empty() {
            return Ok(self);
        }
        let minimum_cp = self
            .caloric
            .source_line()
            .knots
            .iter()
            .map(|(_, answer)| answer.evidence.value.value)
            .fold(f64::INFINITY, f64::min);
        positive(minimum_cp, "original minimum source Cp")?;
        let mut fresh_enthalpies = Vec::with_capacity(self.external_mass.len());
        for &temperature in &self.topology.inlet_temperatures {
            poll(cx)?;
            fresh_enthalpies.push(
                temperature
                    .map(|temperature| {
                        self.context()
                            .enthalpy(self.topology.reference_k, temperature)
                    })
                    .transpose()?,
            );
        }
        self.feedback = Some(CaloricFeedback {
            links,
            sources,
            fractions,
            remaining_exhaust,
            minimum_cp,
            fresh_enthalpies,
            config,
        });
        Ok(self)
    }

    /// Original imposed return topology, in deterministic supply/return order.
    #[must_use]
    pub fn recirculation_links(&self) -> &[RecirculationLink] {
        self.feedback
            .as_ref()
            .map_or(&[], |feedback| &feedback.links)
    }
}

impl CaloricFeedback {
    pub(super) fn march(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        walls: Option<&[f64]>,
    ) -> Result<TransportMarch, TransportError> {
        let mut current = self.trial(network, cx, walls, self.fresh_enthalpies.clone())?;
        for iteration in 0..self.config.max_iterations {
            poll(cx)?;
            if self.converged(network, &current) {
                let report = self.audit(network, cx, &current)?;
                current.march.recirculation = Some(report);
                return Ok(current.march);
            }
            if iteration + 1 == self.config.max_iterations {
                return Err(self.no_convergence(iteration + 1, &current));
            }
            let matrix = self.jacobian(network, cx, &current.branch_enthalpy_gains)?;
            let correction = solve_feedback(cx, &matrix, &current.mismatch, false)?;
            let mut scale = 1.0;
            let mut accepted = None;
            for _ in 0..32 {
                poll(cx)?;
                let mut enthalpies = current.inlet_enthalpies.clone();
                let mut supported = true;
                for (&node, &delta_h) in self.sources.iter().zip(&correction) {
                    poll(cx)?;
                    let start = current.inlet_enthalpies[node].expect("admitted return supply");
                    match finite(start + scale * delta_h, "trial supply enthalpy") {
                        Ok(enthalpy) => enthalpies[node] = Some(enthalpy),
                        Err(_) => {
                            supported = false;
                            break;
                        }
                    }
                }
                if supported {
                    match self.trial(network, cx, walls, enthalpies) {
                        Ok(trial)
                            if self.converged(network, &trial)
                                || trial.heat_residual
                                    < current.heat_residual * (1.0 - 1e-4 * scale) =>
                        {
                            accepted = Some(trial);
                            break;
                        }
                        Err(TransportError::Interrupted) => {
                            return Err(TransportError::Interrupted);
                        }
                        // A candidate is never published unless every branch,
                        // source query and native watt gate admits that state.
                        _ => {}
                    }
                }
                scale *= 0.5;
            }
            let Some(next) = accepted else {
                return Err(self.no_convergence(iteration + 1, &current));
            };
            current = next;
        }
        unreachable!("a positive bounded iteration count returns above")
    }

    fn converged(&self, network: &CaloricTransportNetwork<'_>, trial: &Trial) -> bool {
        trial.temperature_residual <= self.config.temperature_tolerance_k
            && trial.heat_residual <= network.topology.config.absolute_heat_tolerance_w
    }

    fn no_convergence(&self, iterations: usize, trial: &Trial) -> TransportError {
        TransportError::CaloricRecirculationNotConverged {
            iterations,
            max_mixing_residual_k: trial.temperature_residual,
            mixing_heat_residual_w: trial.heat_residual,
        }
    }

    fn trial(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        walls: Option<&[f64]>,
        inlet_enthalpies: Vec<Option<f64>>,
    ) -> Result<Trial, TransportError> {
        let reference = network.topology.reference_k;
        let mut inlets = Vec::with_capacity(inlet_enthalpies.len());
        for &enthalpy in &inlet_enthalpies {
            poll(cx)?;
            inlets.push(
                enthalpy
                    .map(|h| network.context().temperature(reference, h))
                    .transpose()?,
            );
        }
        let mut source = network.context();
        source.reference_k = Some(reference);
        source.inlet_enthalpies = Some(&inlet_enthalpies);
        let mut node_enthalpies = vec![None; inlet_enthalpies.len()];
        let mut branch_enthalpy_gains = vec![1.0; network.topology.models.len()];
        let march = network.topology.march_once_with_source(
            cx,
            walls,
            &inlets,
            None,
            Some(source),
            Some(CaloricCapture {
                node_enthalpies: &mut node_enthalpies,
                branch_enthalpy_gains: &mut branch_enthalpy_gains,
                circulation_residuals_w: None,
            }),
        )?;
        let mut mismatch = vec![0.0; self.sources.len()];
        let mut corrections = vec![0.0; self.sources.len()];
        for (i, &node) in self.sources.iter().enumerate() {
            poll(cx)?;
            mismatch[i] = finite(
                (1.0 - self.fractions[i])
                    * (self.fresh_enthalpies[node].expect("fresh supply")
                        - inlet_enthalpies[node].expect("mixed supply")),
                "fresh mixing enthalpy",
            )?;
        }
        for link in &self.links {
            poll(cx)?;
            let i = self
                .sources
                .binary_search(&link.supply_node)
                .expect("admitted supply");
            let returned = node_enthalpies[link.return_node].ok_or(
                TransportError::InvalidInput("returned exhaust has no retained enthalpy"),
            )?;
            let change = returned - inlet_enthalpies[link.supply_node].expect("mixed supply");
            add(
                &mut mismatch[i],
                &mut corrections[i],
                link.fraction * change,
            )?;
        }
        let (mut heat_residual, mut temperature_residual) = (0.0, 0.0_f64);
        for (i, &node) in self.sources.iter().enumerate() {
            heat_residual = finite(
                heat_residual + network.external_mass[node] * mismatch[i].abs(),
                "mixing heat residual",
            )?;
            // The positive source minimum bounds the inverse's temperature
            // error even across knots. It also avoids trying to invert a
            // sub-ulp residual merely to decide that an accurate solve is done.
            temperature_residual = temperature_residual.max(finite(
                mismatch[i].abs() / self.minimum_cp,
                "mixing temperature bound",
            )?);
        }
        Ok(Trial {
            inlets,
            inlet_enthalpies,
            node_enthalpies,
            branch_enthalpy_gains,
            march,
            mismatch,
            temperature_residual,
            heat_residual,
        })
    }

    /// Eliminate the native DAG using the gains retained by the accepted
    /// primal channel pass, then share them across all supply columns. No
    /// duplicate channel march or finite-difference network probe is needed.
    fn jacobian(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        gains: &[f64],
    ) -> Result<Vec<Vec<f64>>, TransportError> {
        let topology = &network.topology;
        let n = self.sources.len();
        let mut matrix = vec![vec![0.0; n]; n];
        for (i, row) in matrix.iter_mut().enumerate() {
            row[i] = 1.0;
        }
        for (column, &source) in self.sources.iter().enumerate() {
            poll(cx)?;
            let mut mass: Vec<_> = network
                .external_mass
                .iter()
                .map(|flow| flow.max(0.0))
                .collect();
            let mut values = vec![0.0; mass.len()];
            values[source] = 1.0;
            for &node in &topology.order {
                poll(cx)?;
                for &branch in &topology.outgoing[node] {
                    poll(cx)?;
                    let down = topology.downstream[branch];
                    let total = positive(
                        mass[down] + topology.mass_flows[branch],
                        "response mixing mass",
                    )?;
                    let arriving = values[node] * gains[branch];
                    values[down] = finite(
                        values[down]
                            + (arriving - values[down]) * (topology.mass_flows[branch] / total),
                        "return enthalpy response",
                    )?;
                    mass[down] = total;
                }
            }
            for link in &self.links {
                let row = self
                    .sources
                    .binary_search(&link.supply_node)
                    .expect("admitted supply");
                matrix[row][column] = finite(
                    matrix[row][column] - link.fraction * values[link.return_node],
                    "return enthalpy Jacobian",
                )?;
            }
        }
        Ok(matrix)
    }

    fn audit(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        trial: &Trial,
    ) -> Result<RecirculationReport, TransportError> {
        let mut report = RecirculationReport {
            caloric_reference_temperature_k: Some(network.reference_temperature.value()),
            supplies: Vec::with_capacity(self.sources.len()),
            streams: Vec::with_capacity(self.links.len()),
            external_heat_gain_w: 0.0,
            heat_imbalance_w: 0.0,
            max_mixing_residual_k: trial.temperature_residual,
        };
        let cp = network.topology.air.specific_heat_j_kg_k;
        for (i, &node) in self.sources.iter().enumerate() {
            poll(cx)?;
            let fresh = (1.0 - self.fractions[i]) * network.external_mass[node];
            let returned = self.fractions[i] * network.external_mass[node];
            report.supplies.push(RecirculatedSupply {
                node,
                fresh_temperature_k: network.topology.inlet_temperatures[node]
                    .expect("fresh supply"),
                mixed_temperature_k: trial.inlets[node].expect("mixed supply"),
                fresh_mass_flow_kg_per_s: fresh,
                return_mass_flow_kg_per_s: returned,
                fresh_capacity_w_per_k: finite(fresh * cp, "reference fresh capacity")?,
                return_capacity_w_per_k: finite(returned * cp, "reference return capacity")?,
            });
        }
        for link in &self.links {
            poll(cx)?;
            let mass = link.fraction * network.external_mass[link.supply_node];
            report.streams.push(RecirculatedStream {
                supply_node: link.supply_node,
                return_node: link.return_node,
                fraction: link.fraction,
                mass_flow_kg_per_s: mass,
                capacity_w_per_k: finite(mass * cp, "reference returned capacity")?,
                temperature_k: trial.march.node_temperatures_k[link.return_node]
                    .expect("accepted return"),
            });
        }
        let (mut gross, mut correction) = (0.0, 0.0);
        for (node, &external) in network.external_mass.iter().enumerate() {
            poll(cx)?;
            let heat = if external > 0.0 {
                let fraction = self
                    .sources
                    .binary_search(&node)
                    .map_or(0.0, |i| self.fractions[i]);
                -external * (1.0 - fraction) * self.fresh_enthalpies[node].expect("fresh supply")
            } else if self.remaining_exhaust[node] > 0.0 {
                self.remaining_exhaust[node]
                    * trial.node_enthalpies[node].ok_or(TransportError::InvalidInput(
                        "remaining exhaust has no retained enthalpy",
                    ))?
            } else {
                0.0
            };
            add(&mut report.external_heat_gain_w, &mut correction, heat)?;
            gross = finite(gross + heat.abs(), "gross source return heat")?;
        }
        let mut gross_wall = 0.0;
        for branch in &trial.march.branches {
            if let Some(path) = &branch.march {
                for segment in &path.segments {
                    poll(cx)?;
                    gross_wall = finite(
                        gross_wall + segment.heat_rate_w.abs(),
                        "gross source wall heat",
                    )?;
                }
            }
        }
        report.heat_imbalance_w = finite(
            report.external_heat_gain_w - trial.march.wall_heat_rate_w,
            "outer source return heat imbalance",
        )?;
        let tolerance = finite(
            network.topology.config.absolute_heat_tolerance_w
                + network.topology.config.relative_heat_tolerance * gross.max(gross_wall),
            "source return watt tolerance",
        )?;
        if report.heat_imbalance_w.abs() > tolerance {
            return Err(TransportError::HeatImbalance {
                residual_w: report.heat_imbalance_w,
                tolerance_w: tolerance,
                hydraulic_defect_w: trial.march.hydraulic_energy_defect_w,
            });
        }
        poll(cx)?;
        Ok(report)
    }
}

fn add(sum: &mut f64, correction: &mut f64, amount: f64) -> Result<(), TransportError> {
    let adjusted = finite(amount - *correction, "source return sum increment")?;
    let next = finite(*sum + adjusted, "source return sum")?;
    *correction = (next - *sum) - adjusted;
    *sum = next;
    Ok(())
}
