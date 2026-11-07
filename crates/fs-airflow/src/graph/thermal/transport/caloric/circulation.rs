//! Closed hydraulic circulation in original relative specific enthalpy.
//!
//! The nonlinear unknowns are the native node enthalpies. Every trial uses
//! original source channel heat and dh_out/dh_in; mass-weighted mixing supplies
//! the residual and analytic Jacobian. A small watt residual alone does not
//! establish temperature accuracy when wall anchors are weak, so admission
//! also bounds the Newton correction using the original source minimum Cp.

use super::super::recirculation::solve_feedback;
use super::*;

/// Bounded nonlinear work for closed source-caloric circulation.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CaloricCirculationConfig {
    /// Positive upper bound on the final Newton correction divided by the
    /// minimum source Cp, K. This is a local convergence test, not a certified
    /// bound on nonlinear, hydraulic, spatial or material-model error.
    pub temperature_tolerance_k: f64,
    /// Maximum nonlinear residual evaluations along accepted Newton iterates,
    /// including the initial state, in 1..=256. Each has at most 32 damped trials.
    pub max_iterations: usize,
}

#[derive(Debug)]
pub(super) struct CaloricCirculation {
    incoming_mass: Vec<f64>,
    minimum_cp: f64,
    config: CaloricCirculationConfig,
}

struct Trial {
    enthalpies: Vec<f64>,
    gains: Vec<f64>,
    mismatch: Vec<f64>,
    heat_residual: f64,
}

impl<'a> CaloricTransportNetwork<'a> {
    /// Admit a closed hydraulic circulation with the original source Cp(T).
    ///
    /// Pressure boundaries are datums, never external fluid supplies. Native
    /// admission checks mass balance at every node, original branch direction,
    /// reference density/Cp, the 64-node cap and a resolved wall-heat anchor in
    /// each active component. A wholly adiabatic active loop has no unique
    /// steady temperature and refuses.
    ///
    /// The fixed-composition source isobar owns every branch and mixer. The
    /// physical reference temperature supplies the native material snapshot
    /// and initial solid guess; it is independent of the source's energy datum.
    /// It does not add a reservoir. Final node mixing and enclosing wall watts
    /// are independently re-evaluated through the native admission gates.
    ///
    /// Flow, density and film coefficients remain prescribed. This does not
    /// infer hydraulic dissipation as heat, pressure work, phase change,
    /// composition mixing, buoyancy feedback, fluid residence time or adjoints.
    #[allow(clippy::too_many_arguments)]
    pub fn new_circulating(
        cx: &Cx<'_>,
        flow: &'a GraphSolution,
        reference_density: Density,
        reference_temperature: Temperature,
        caloric: Arc<IsobaricHeatCapacityLaw>,
        models: Vec<BranchThermalModel>,
        config: TransportConfig,
        circulation: CaloricCirculationConfig,
    ) -> Result<Self, TransportError> {
        Self::new_inner(
            cx,
            flow,
            reference_density,
            reference_temperature,
            caloric,
            models,
            &[],
            config,
            Some(circulation),
        )
    }
}

impl CaloricCirculation {
    pub(super) fn new(
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        config: CaloricCirculationConfig,
    ) -> Result<Self, TransportError> {
        poll(cx)?;
        positive(
            config.temperature_tolerance_k,
            "circulation temperature tolerance",
        )?;
        if !(1..=256).contains(&config.max_iterations) {
            return Err(TransportError::InvalidInput(
                "source-caloric circulation iteration budget must lie in 1..=256",
            ));
        }
        let mut incoming_mass = vec![0.0; network.topology.order.len()];
        for branches in &network.topology.outgoing {
            for &branch in branches {
                poll(cx)?;
                let down = network.topology.downstream[branch];
                incoming_mass[down] = finite(
                    incoming_mass[down] + network.topology.mass_flows[branch],
                    "source circulation incoming mass",
                )?;
            }
        }
        let minimum_cp = network
            .caloric
            .source_line()
            .knots
            .iter()
            .map(|(_, answer)| answer.evidence.value.value)
            .fold(f64::INFINITY, f64::min);
        positive(minimum_cp, "source circulation minimum Cp")?;
        Ok(Self {
            incoming_mass,
            minimum_cp,
            config,
        })
    }

    pub(super) fn march(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        walls: &[f64],
    ) -> Result<TransportMarch, TransportError> {
        poll(cx)?;
        let reference = network.topology.reference_k;
        let knots = &network.caloric.source_line().knots;
        let lower = knots.first().expect("admitted source support").0;
        let upper = knots.last().expect("admitted source support").0;
        // Start in the supported part of the wall-temperature hull. This is
        // only a Newton seed: all accepted states still solve the same source
        // equation and the original declared reference owns energy accounting.
        let seed = if walls.is_empty() {
            reference
        } else {
            let min_wall = walls.iter().copied().fold(f64::INFINITY, f64::min);
            let max_wall = walls.iter().copied().fold(0.0_f64, f64::max);
            (min_wall + 0.5 * (max_wall - min_wall)).clamp(lower, upper)
        };
        let h = network.context().enthalpy(reference, seed)?;
        let mut current = match self.trial(network, cx, walls, vec![h; self.incoming_mass.len()]) {
            Ok(trial) => trial,
            Err(TransportError::Interrupted) => return Err(TransportError::Interrupted),
            // A bath can lie outside the source span while its fluid path is
            // wholly supported. Clipping the preferred seed to that span's
            // edge must not prevent the admitted physical reference from
            // supplying a feasible initial trajectory. Try it at most once.
            Err(_) if seed != reference => {
                self.trial(network, cx, walls, vec![0.0; self.incoming_mass.len()])?
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
                    "source circulation temperature correction",
                )?);
            }
            if temperature_bound <= self.config.temperature_tolerance_k
                && current.heat_residual <= network.topology.config.absolute_heat_tolerance_w
            {
                // Re-march the accepted states without the private trial
                // capture. Both native mixer and global watt gates now apply;
                // no residual, hydraulic defect or rejected heat is corrected.
                let mut source = network.context();
                source.reference_k = Some(reference);
                return network.topology.march_once_with_source(
                    cx,
                    Some(walls),
                    &network.topology.inlet_temperatures,
                    Some(&current.enthalpies),
                    Some(source),
                    None,
                );
            }
            if iteration + 1 == self.config.max_iterations {
                return Err(self.no_convergence(iteration + 1, temperature_bound, &current));
            }
            let mut scale = 1.0;
            let mut accepted = None;
            for _ in 0..32 {
                poll(cx)?;
                let candidate: Option<Vec<f64>> = current
                    .enthalpies
                    .iter()
                    .zip(&correction)
                    .map(|(&h, &delta)| {
                        let value = h + scale * delta;
                        value.is_finite().then_some(value)
                    })
                    .collect();
                if let Some(enthalpies) = candidate {
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
                        // All material states and each physical branch must be
                        // supported; no extrapolation or partial field escapes.
                        _ => {}
                    }
                }
                scale *= 0.5;
            }
            current = accepted
                .ok_or_else(|| self.no_convergence(iteration + 1, temperature_bound, &current))?;
        }
        unreachable!("positive bounded circulation iteration count returns above")
    }

    fn trial(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        walls: &[f64],
        enthalpies: Vec<f64>,
    ) -> Result<Trial, TransportError> {
        let mut node_enthalpies = vec![None; enthalpies.len()];
        let mut gains = vec![1.0; network.topology.models.len()];
        let mut mismatch = vec![0.0; enthalpies.len()];
        let mut source = network.context();
        source.reference_k = Some(network.topology.reference_k);
        network.topology.march_once_with_source(
            cx,
            Some(walls),
            &network.topology.inlet_temperatures,
            Some(&enthalpies),
            Some(source),
            Some(CaloricCapture {
                node_enthalpies: &mut node_enthalpies,
                branch_enthalpy_gains: &mut gains,
                circulation_residuals_w: Some(&mut mismatch),
            }),
        )?;
        let mut heat_residual = 0.0;
        for (value, &mass) in mismatch.iter_mut().zip(&self.incoming_mass) {
            poll(cx)?;
            heat_residual = finite(
                heat_residual + value.abs(),
                "source circulation mixer residual",
            )?;
            if mass > 0.0 {
                *value = finite(*value / mass, "source circulation enthalpy residual")?;
            }
        }
        Ok(Trial {
            enthalpies,
            gains,
            mismatch,
            heat_residual,
        })
    }

    fn correction(
        &self,
        network: &CaloricTransportNetwork<'_>,
        cx: &Cx<'_>,
        trial: &Trial,
    ) -> Result<Vec<f64>, TransportError> {
        let topology = &network.topology;
        let n = self.incoming_mass.len();
        let mut matrix = vec![vec![0.0; n]; n];
        for (node, row) in matrix.iter_mut().enumerate() {
            row[node] = 1.0;
        }
        for (node, branches) in topology.outgoing.iter().enumerate() {
            for &branch in branches {
                poll(cx)?;
                let down = topology.downstream[branch];
                let weight = topology.mass_flows[branch] / self.incoming_mass[down];
                matrix[down][node] = finite(
                    matrix[down][node] - weight * trial.gains[branch],
                    "source circulation enthalpy Jacobian",
                )?;
            }
        }
        solve_feedback(cx, &matrix, &trial.mismatch, false)
    }

    fn no_convergence(
        &self,
        iterations: usize,
        temperature_bound: f64,
        trial: &Trial,
    ) -> TransportError {
        TransportError::CaloricCirculationNotConverged {
            iterations,
            temperature_correction_bound_k: temperature_bound,
            mixing_heat_residual_w: trial.heat_residual,
        }
    }
}
