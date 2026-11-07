//! Finite exhaust-to-fresh heat recovery on the native source-caloric DAG.
//!
//! An exchanger transfers heat between two external streams without returning
//! any mass. The receiving inlet enthalpies close a reduced nonlinear system;
//! the native DAG still owns branch heat, mixing and its inner watt gate.
//! A second ledger uses the ORIGINAL fresh supplies and POST-exchanger exhaust.
//! One retained signed Q appears at both exchanger ports exactly once.

use super::super::recirculation::{MAX_RECIRCULATION_SUPPLIES, solve_feedback};
use super::*;
use crate::conjugate::caloric::exchanger::{
    CaloricHeatExchanger, ExchangerConfig, ExchangerError, ExchangerInletResponse,
    ExchangerResponseConfig, ExchangerState, ExchangerStream, FlowArrangement,
};
use fs_math::dd::Dd;

/// One whole external exhaust exchanging heat with one whole fresh supply.
/// Each supply and each exhaust can belong to at most one recuperator.
/// The first exchanger stream is always exhaust; the second is fresh supply.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricHeatRecoveryLink {
    pub supply_node: usize,
    pub exhaust_node: usize,
    pub arrangement: FlowArrangement,
    pub conductance_w_per_k: f64,
}

/// Native reduced Newton limits plus the original exchanger's own accuracy
/// and work limits. Exchanger allowances do not relax either network watt gate.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricHeatRecoveryConfig {
    /// Positive bound on the final Newton correction divided by minimum
    /// original Cp. This is a local convergence test, not a global error proof.
    pub temperature_tolerance_k: f64,
    /// Accepted iterates, including the initial trial, in 1..=256. Each
    /// correction has at most 32 damped trials. A failed fresh seed may use one
    /// zero-wall-heat seed solve with the same bounded work policy.
    pub max_iterations: usize,
    pub exchanger: ExchangerConfig,
    pub response: ExchangerResponseConfig,
}

impl Default for CaloricHeatRecoveryConfig {
    fn default() -> Self {
        Self {
            temperature_tolerance_k: 1e-8,
            max_iterations: 64,
            exchanger: ExchangerConfig::default(),
            response: ExchangerResponseConfig::default(),
        }
    }
}

/// One accepted finite exchanger, including its original physical audits.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricRecoveredStream {
    pub link: CaloricHeatRecoveryLink,
    pub supply_mass_flow_kg_per_s: f64,
    pub exhaust_mass_flow_kg_per_s: f64,
    pub fresh_temperature_k: f64,
    /// Actual external supply marched by the DAG. Its difference from the
    /// exchanger's second outlet is bounded by the independently gated closure.
    pub supply_temperature_k: f64,
    pub exhaust_inlet_temperature_k: f64,
    /// Positive heat in exchanger means recovered heat enters the supply.
    /// Negative heat precools a hotter fresh supply and warms the exhaust.
    pub exchanger: ExchangerState,
    pub response: ExchangerInletResponse,
    /// m_supply*(h_fresh-h_supply)+Q, recomputed in relative enthalpy.
    pub supply_heat_residual_w: f64,
    /// Original calorimetry at the represented exhaust inlet minus its retained
    /// DAG enthalpy coordinate, W. This exposes the temperature handoff rounding.
    pub exhaust_inlet_energy_residual_w: f64,
    /// Numerical allowance added ONLY to this endpoint handoff audit. It is
    /// never added to the interface, inner DAG, or outer boundary watt budgets.
    pub exhaust_inlet_energy_roundoff_allowance_w: f64,
}

/// Independent recovered-boundary balance. TransportMarch's original inner
/// balance remains unchanged: its supplies are conditioned and its exhausts
/// precede recovery. This report counts fresh supplies and recovered exhausts.
#[derive(Clone, Debug, PartialEq)]
pub struct CaloricHeatRecoveryReport {
    pub streams: Vec<CaloricRecoveredStream>,
    pub enthalpy_reference_k: f64,
    /// POST-exchanger exhaust enthalpy minus ORIGINAL fresh inlet enthalpy, W.
    pub external_heat_gain_w: f64,
    pub heat_imbalance_w: f64,
    /// Sum of absolute supply/exchanger closure residuals, W.
    pub supply_heat_residual_w: f64,
    pub temperature_correction_bound_k: f64,
    /// Sum over recuperators only; wall exchange and native junction mixing
    /// are not relabeled as a complete-network entropy account.
    pub exchanger_entropy_generation_w_per_k: f64,
    pub exchanger_entropy_roundoff_allowance_w_per_k: f64,
    pub iterations: usize,
}

#[derive(Debug)]
pub(super) struct CaloricRecovery {
    links: Vec<CaloricHeatRecoveryLink>,
    exhaust_link: Vec<Option<usize>>,
    fresh_enthalpies: Vec<Option<f64>>,
    minimum_cp: f64,
    config: CaloricHeatRecoveryConfig,
}

struct Trial {
    inlet_enthalpies: Vec<Option<f64>>,
    node_enthalpies: Vec<Option<f64>>,
    branch_enthalpy_gains: Vec<f64>,
    march: TransportMarch,
    streams: Vec<CaloricRecoveredStream>,
    mismatch: Vec<f64>,
    heat_residual: f64,
}

impl CaloricTransportNetwork<'_> {
    /// Recover finite exhaust heat into fresh external supplies on this DAG.
    ///
    /// Actual native external MASS flows determine both exchanger streams.
    /// The source isobar, density snapshot, hydraulic graph, branch directions,
    /// wall regions and shared-solid driver remain native. No return mass,
    /// pressure-reset path, fan work, new reservoir or hydraulic edge is added.
    ///
    /// One-to-one whole-port ownership and the native 64-supply feedback cap
    /// bound this first composition. Mass-return links and closed circulations
    /// remain distinct modes and cannot be combined with recovery. Individual
    /// exchangers retain their original source, pinch, entropy and UA admission.
    pub fn with_heat_recovery(
        mut self,
        cx: &Cx<'_>,
        mut links: Vec<CaloricHeatRecoveryLink>,
        config: CaloricHeatRecoveryConfig,
    ) -> Result<Self, TransportError> {
        poll(cx)?;
        positive(
            config.temperature_tolerance_k,
            "heat recovery temperature tolerance",
        )?;
        if !(1..=256).contains(&config.max_iterations) || links.len() > MAX_RECIRCULATION_SUPPLIES {
            return Err(TransportError::InvalidInput(
                "heat recovery iteration or supply cap",
            ));
        }
        if self.feedback.is_some() || self.circulation.is_some() || self.recovery.is_some() {
            return Err(TransportError::InvalidInput(
                "finite heat recovery requires an unbound source-caloric DAG",
            ));
        }
        links.sort_by_key(|link| (link.supply_node, link.exhaust_node));
        let mut previous_supply = None;
        let mut exhaust_link = vec![None; self.external_mass.len()];
        for (index, link) in links.iter().enumerate() {
            poll(cx)?;
            if link.supply_node >= self.external_mass.len()
                || link.exhaust_node >= self.external_mass.len()
                || self.external_mass[link.supply_node] <= 0.0
                || self.external_mass[link.exhaust_node] >= 0.0
                || previous_supply == Some(link.supply_node)
                || exhaust_link[link.exhaust_node].is_some()
            {
                return Err(TransportError::InvalidInput(
                    "heat recovery needs uniquely owned external supply and exhaust ports",
                ));
            }
            finite(link.conductance_w_per_k, "heat recovery conductance")?;
            if link.conductance_w_per_k < 0.0 {
                return Err(TransportError::InvalidInput(
                    "negative heat recovery conductance",
                ));
            }
            positive(self.external_mass[link.supply_node], "recovery supply mass")?;
            positive(
                -self.external_mass[link.exhaust_node],
                "recovery exhaust mass",
            )?;
            previous_supply = Some(link.supply_node);
            exhaust_link[link.exhaust_node] = Some(index);
        }
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
        positive(minimum_cp, "minimum recovery source Cp")?;
        let mut fresh_enthalpies = Vec::with_capacity(self.external_mass.len());
        for &temperature in &self.topology.inlet_temperatures {
            poll(cx)?;
            fresh_enthalpies.push(
                temperature
                    .map(|t| self.context().enthalpy(self.topology.reference_k, t))
                    .transpose()?,
            );
        }
        self.recovery = Some(CaloricRecovery {
            links,
            exhaust_link,
            fresh_enthalpies,
            minimum_cp,
            config,
        });
        Ok(self)
    }

    #[must_use]
    pub fn heat_recovery_links(&self) -> &[CaloricHeatRecoveryLink] {
        self.recovery
            .as_ref()
            .map_or(&[], |recovery| recovery.links.as_slice())
    }
}

impl CaloricRecovery {
    pub(super) fn march(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        walls: Option<&[f64]>,
    ) -> Result<TransportMarch, TransportError> {
        let (trial, iterations, temperature_bound) = self.solve(network, cx, walls)?;
        self.audit(network, cx, trial, iterations, temperature_bound)
    }

    fn solve(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        walls: Option<&[f64]>,
    ) -> Result<(Trial, usize, f64), TransportError> {
        let mut current = match self.trial(network, cx, walls, self.fresh_enthalpies.clone()) {
            Ok(trial) => trial,
            Err(TransportError::Interrupted) => return Err(TransportError::Interrupted),
            // Cross-component recovery can make the unconditioned fresh seed
            // infeasible even when the final wall path is supported. The native
            // zero-wall-heat composition supplies one bounded supported retry.
            Err(_) if walls.is_some() => {
                let (seed, _, _) = self.solve(network, cx, None)?;
                self.trial(network, cx, walls, seed.inlet_enthalpies)?
            }
            Err(error) => return Err(error),
        };
        for iteration in 0..self.config.max_iterations {
            poll(cx)?;
            let correction = self.correction(network, cx, &current)?;
            let mut temperature_bound = 0.0_f64;
            for &delta in &correction {
                temperature_bound = temperature_bound.max(finite(
                    delta.abs() / self.minimum_cp,
                    "heat recovery temperature correction",
                )?);
            }
            if current.heat_residual <= network.topology.config.absolute_heat_tolerance_w
                && temperature_bound <= self.config.temperature_tolerance_k
            {
                return Ok((current, iteration + 1, temperature_bound));
            }
            if iteration + 1 == self.config.max_iterations {
                return Err(self.no_convergence(iteration + 1, temperature_bound, &current));
            }
            let mut scale = 1.0;
            let mut accepted = None;
            for _ in 0..32 {
                poll(cx)?;
                let mut enthalpies = current.inlet_enthalpies.clone();
                let mut supported = true;
                for (link, &delta) in self.links.iter().zip(&correction) {
                    let old = current.inlet_enthalpies[link.supply_node].expect("recovery supply");
                    match finite(old + scale * delta, "trial recovery inlet enthalpy") {
                        Ok(value) => enthalpies[link.supply_node] = Some(value),
                        Err(_) => {
                            supported = false;
                            break;
                        }
                    }
                }
                if supported {
                    match self.trial(network, cx, walls, enthalpies) {
                        Ok(trial)
                            if trial.heat_residual
                                < current.heat_residual * (1.0 - 1e-4 * scale) =>
                        {
                            accepted = Some(trial);
                            break;
                        }
                        Err(TransportError::Interrupted) => {
                            return Err(TransportError::Interrupted);
                        }
                        _ => {}
                    }
                }
                scale *= 0.5;
            }
            current = accepted
                .ok_or_else(|| self.no_convergence(iteration + 1, temperature_bound, &current))?;
        }
        unreachable!("positive bounded recovery iteration count returns above")
    }

    fn trial(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        walls: Option<&[f64]>,
        inlet_enthalpies: Vec<Option<f64>>,
    ) -> Result<Trial, TransportError> {
        poll(cx)?;
        let reference = network.topology.reference_k;
        let mut inlets = Vec::with_capacity(inlet_enthalpies.len());
        for (node, &enthalpy) in inlet_enthalpies.iter().enumerate() {
            poll(cx)?;
            inlets.push(match enthalpy {
                Some(h) if Some(h) == self.fresh_enthalpies[node] => {
                    network.topology.inlet_temperatures[node]
                }
                Some(h) => Some(network.context().temperature(reference, h)?),
                None => None,
            });
        }
        let mut source = network.context();
        source.reference_k = Some(reference);
        source.inlet_enthalpies = Some(&inlet_enthalpies);
        let mut node_enthalpies = vec![None; inlets.len()];
        let mut branch_enthalpy_gains = vec![1.0; network.topology.models.len()];
        // This is the ordinary admitted DAG pass, including its independent
        // inner watt gate. It never uses the closed-loop trial-gate bypass.
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
        let mut streams = Vec::with_capacity(self.links.len());
        let mut mismatch = Vec::with_capacity(self.links.len());
        let mut total_residual = Dd::ZERO;
        for link in &self.links {
            poll(cx)?;
            let supply_mass = network.external_mass[link.supply_node];
            let exhaust_mass = -network.external_mass[link.exhaust_node];
            let fresh_temperature = network.topology.inlet_temperatures[link.supply_node]
                .expect("original fresh supply");
            let exhaust_temperature = march.node_temperatures_k[link.exhaust_node].ok_or(
                TransportError::InvalidInput("recovery exhaust has no native temperature"),
            )?;
            let exhaust_enthalpy = node_enthalpies[link.exhaust_node].ok_or(
                TransportError::InvalidInput("recovery exhaust has no retained enthalpy"),
            )?;
            let (rounding_residual, roundoff) =
                exhaust_handoff(network, exhaust_mass, exhaust_temperature, exhaust_enthalpy)?;
            let exchanger = CaloricHeatExchanger::new(
                ExchangerStream::new(exhaust_temperature, exhaust_mass, network.caloric.clone())
                    .map_err(exchanger_error)?,
                ExchangerStream::new(fresh_temperature, supply_mass, network.caloric.clone())
                    .map_err(exchanger_error)?,
                link.arrangement,
                link.conductance_w_per_k,
            )
            .map_err(exchanger_error)?;
            let (state, response) = exchanger
                .solve_with_inlet_response(cx, self.config.exchanger, self.config.response)
                .map_err(exchanger_error)?;
            let fresh_h = self.fresh_enthalpies[link.supply_node].expect("fresh enthalpy");
            let supply_h = inlet_enthalpies[link.supply_node].expect("conditioned supply enthalpy");
            let residual = Dd::from_f64(supply_mass)
                * (Dd::from_f64(fresh_h) - Dd::from_f64(supply_h))
                + Dd::from_f64(state.heat_rate_first_to_second_w);
            let watts = finite(residual.to_f64(), "recovery supply heat residual")?;
            total_residual = total_residual + residual.abs();
            mismatch.push(finite(
                watts / supply_mass,
                "recovery specific-enthalpy residual",
            )?);
            streams.push(CaloricRecoveredStream {
                link: *link,
                supply_mass_flow_kg_per_s: supply_mass,
                exhaust_mass_flow_kg_per_s: exhaust_mass,
                fresh_temperature_k: fresh_temperature,
                supply_temperature_k: inlets[link.supply_node].expect("conditioned supply"),
                exhaust_inlet_temperature_k: exhaust_temperature,
                exchanger: state,
                response,
                supply_heat_residual_w: watts,
                exhaust_inlet_energy_residual_w: rounding_residual,
                exhaust_inlet_energy_roundoff_allowance_w: roundoff,
            });
        }
        Ok(Trial {
            inlet_enthalpies,
            node_enthalpies,
            branch_enthalpy_gains,
            march,
            streams,
            mismatch,
            heat_residual: finite(total_residual.to_f64(), "total recovery supply residual")?,
        })
    }

    fn correction(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        trial: &Trial,
    ) -> Result<Vec<f64>, TransportError> {
        let n = self.links.len();
        let mut matrix = vec![vec![0.0; n]; n];
        for (column, link) in self.links.iter().enumerate() {
            poll(cx)?;
            let values = network.inlet_enthalpy_response(
                cx,
                link.supply_node,
                &trial.branch_enthalpy_gains,
            )?;
            for (row, stream) in trial.streams.iter().enumerate() {
                let response = stream.response;
                let propagation = Dd::from_f64(values[stream.link.exhaust_node])
                    * Dd::from_f64(stream.exhaust_mass_flow_kg_per_s)
                    / Dd::from_f64(stream.supply_mass_flow_kg_per_s);
                // The smaller retained component supplies the accurate loss.
                // A near-unit transfer never erases the finite-UA anchor when
                // the reduced matrix subtracts that feedback from identity.
                let coefficient = if response.first_outlet_enthalpy_gain <= 0.5 {
                    propagation - propagation * Dd::from_f64(response.first_outlet_enthalpy_gain)
                } else {
                    propagation * Dd::from_f64(response.heat_from_first_inlet_enthalpy_kg_per_s)
                        / Dd::from_f64(stream.exhaust_mass_flow_kg_per_s)
                };
                matrix[row][column] = finite(
                    (Dd::from_f64(if row == column { 1.0 } else { 0.0 }) - coefficient).to_f64(),
                    "heat recovery enthalpy Jacobian",
                )?;
            }
        }
        solve_feedback(cx, &matrix, &trial.mismatch, false)
    }

    fn audit(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        trial: Trial,
        iterations: usize,
        temperature_bound: f64,
    ) -> Result<TransportMarch, TransportError> {
        let mut external = Dd::ZERO;
        let mut gross = Dd::ZERO;
        for (node, &mass) in network.external_mass.iter().enumerate() {
            poll(cx)?;
            let term = if mass > 0.0 {
                -Dd::from_f64(mass)
                    * Dd::from_f64(self.fresh_enthalpies[node].expect("original fresh inlet"))
            } else if mass < 0.0 {
                let before = Dd::from_f64(-mass)
                    * Dd::from_f64(trial.node_enthalpies[node].ok_or(
                        TransportError::InvalidInput(
                            "external exhaust lacks retained source enthalpy",
                        ),
                    )?);
                match self.exhaust_link[node] {
                    Some(index) => {
                        before
                            - Dd::from_f64(
                                trial.streams[index].exchanger.heat_rate_first_to_second_w,
                            )
                    }
                    None => before,
                }
            } else {
                Dd::ZERO
            };
            finite(term.to_f64(), "recovered external heat term")?;
            external = external + term;
            gross = gross + term.abs();
        }
        let mut gross_wall = 0.0;
        for branch in &trial.march.branches {
            if let Some(path) = &branch.march {
                for segment in &path.segments {
                    poll(cx)?;
                    gross_wall = finite(
                        gross_wall + segment.heat_rate_w.abs(),
                        "recovery gross wall heat",
                    )?;
                }
            }
        }
        let external_heat = finite(external.to_f64(), "recovered external heat gain")?;
        let imbalance = finite(
            (external - Dd::from_f64(trial.march.wall_heat_rate_w)).to_f64(),
            "outer recovered heat imbalance",
        )?;
        let gross = finite(gross.to_f64(), "gross recovered external heat")?;
        let tolerance = finite(
            network.topology.config.absolute_heat_tolerance_w
                + network.topology.config.relative_heat_tolerance * gross.max(gross_wall),
            "outer recovered heat tolerance",
        )?;
        if imbalance.abs() > tolerance {
            return Err(TransportError::HeatImbalance {
                residual_w: imbalance,
                tolerance_w: tolerance,
                hydraulic_defect_w: trial.march.hydraulic_energy_defect_w,
            });
        }
        let mut entropy = Dd::ZERO;
        let mut entropy_allowance = Dd::ZERO;
        for stream in &trial.streams {
            poll(cx)?;
            entropy = entropy + Dd::from_f64(stream.exchanger.entropy_generation_w_per_k);
            entropy_allowance = entropy_allowance
                + Dd::from_f64(stream.exchanger.entropy_roundoff_allowance_w_per_k);
        }
        let entropy = finite(entropy.to_f64(), "recovery exchanger entropy")?;
        let entropy_allowance = finite(entropy_allowance.to_f64(), "recovery entropy allowance")?;
        if entropy < -entropy_allowance {
            return Err(TransportError::InvalidInput(
                "recovery exchanger entropy balance",
            ));
        }
        let report = CaloricHeatRecoveryReport {
            streams: trial.streams,
            enthalpy_reference_k: trial.march.enthalpy_reference_k,
            external_heat_gain_w: external_heat,
            heat_imbalance_w: imbalance,
            supply_heat_residual_w: trial.heat_residual,
            temperature_correction_bound_k: temperature_bound,
            exchanger_entropy_generation_w_per_k: entropy,
            exchanger_entropy_roundoff_allowance_w_per_k: entropy_allowance,
            iterations,
        };
        let mut march = trial.march;
        march.heat_recovery = Some(report);
        poll(cx)?;
        Ok(march)
    }

    fn no_convergence(&self, iterations: usize, bound: f64, trial: &Trial) -> TransportError {
        TransportError::CaloricHeatRecoveryNotConverged {
            iterations,
            temperature_correction_bound_k: bound,
            supply_heat_residual_w: trial.heat_residual,
        }
    }
}

fn exhaust_handoff(
    network: &CaloricTransportNetwork<'_>,
    mass: f64,
    temperature: f64,
    retained_h: f64,
) -> Result<(f64, f64), TransportError> {
    let reference = network.topology.reference_k;
    let represented = network.context().enthalpy(reference, temperature)?;
    let residual =
        (Dd::from_f64(mass) * (Dd::from_f64(represented) - Dd::from_f64(retained_h))).to_f64();
    let inlet_cp = network
        .caloric
        .specific_heat_cp_j_kg_k(temperature)
        .map_err(source_error)?;
    let reference_cp = network
        .caloric
        .specific_heat_cp_j_kg_k(reference)
        .map_err(source_error)?;
    let allowance = finite(
        64.0 * f64::EPSILON
            * mass
            * (retained_h.abs() + inlet_cp * temperature + reference_cp * reference),
        "exhaust enthalpy handoff allowance",
    )?;
    finite(residual, "exhaust enthalpy handoff residual")?;
    if residual.abs() > network.topology.config.absolute_heat_tolerance_w + allowance {
        return Err(TransportError::InvalidInput(
            "original exhaust enthalpy handoff did not close",
        ));
    }
    Ok((residual, allowance))
}

fn exchanger_error(error: ExchangerError) -> TransportError {
    match error {
        ExchangerError::Interrupted => TransportError::Interrupted,
        error => TransportError::Exchanger(error),
    }
}
