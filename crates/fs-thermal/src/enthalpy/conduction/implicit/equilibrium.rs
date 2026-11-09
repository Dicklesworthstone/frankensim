//! Stationary operating points of the original nonlinear conduction graph.
//!
//! Solve P(E)=0 with E=m*h and the SAME heat rates and analytic -dP/dE used by
//! native implicit evolution and thermal impedance. No artificial timestep,
//! heat-capacity fit, finite difference or temperature clamp selects the root.
//! Driven steady states may carry heat and produce entropy; they are not
//! thermodynamic equilibrium of an isolated inventory.
//!
//! Every currently conducting component must reach a temperature-dependent
//! external bath. Insulated inventories need an energy constraint, so this
//! owner refuses them instead of selecting an arbitrary null-space temperature.
//! Clocks, authored events and prescribed geometry stay frozen at the queried
//! instant. The returned point is a read-only field, never a time trajectory or
//! a heat/work credit. Existence, uniqueness, stability, phase allocation and
//! continuum accuracy are not inferred from a small stationary residual.

use super::harmonic::{ThermalImpedanceBudget, ThermalImpedanceError, ThermalImpedanceModel};
use super::*;

/// Conserved original storage in closed, disconnected conducting components.
pub mod isolated;

/// Bounded nonlinear work and physical heat-rate acceptance, in watts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalOperatingControls {
    pub maximum_newton_iterations: usize,
    pub maximum_backtracks: usize,
    /// Includes the initial query and all rejected constitutive trials.
    pub maximum_residual_evaluations: usize,
    pub cell_tolerance_w: f64,
    pub total_tolerance_w: f64,
}

/// Work counts and the actual original stationary balance at the returned point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalOperatingReport {
    pub nonlinear_iterations: usize,
    pub linear_iterations: usize,
    pub residual_evaluations: usize,
    pub maximum_cell_residual_w: f64,
    pub total_residual_w: f64,
}

#[derive(Clone, Debug, PartialEq)]
pub enum ThermalOperatingError {
    Original(ImplicitConductionError),
    /// This component has no currently active temperature-dependent bath.
    Unanchored {
        cell: CellId,
    },
    NotConverged {
        maximum_residual_w: f64,
    },
}
impl From<ImplicitConductionError> for ThermalOperatingError {
    fn from(error: ImplicitConductionError) -> Self {
        Self::Original(error)
    }
}
impl From<DepositionError> for ThermalOperatingError {
    fn from(error: DepositionError) -> Self {
        Self::Original(error.into())
    }
}
impl From<TridiagonalError> for ThermalOperatingError {
    fn from(error: TridiagonalError) -> Self {
        Self::Original(error.into())
    }
}
impl core::fmt::Display for ThermalOperatingError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "stationary thermal balance: {self:?}")
    }
}
impl std::error::Error for ThermalOperatingError {}

/// A stationary field bound to the exact original materials, geometry, sources
/// and baths. The private constitutive snapshot cannot be mistaken for an
/// evolved timeline. Use its fields to author initial data explicitly.
pub struct ThermalOperatingPoint<'s, 'a, 'm> {
    owner: &'s ImplicitConductionChain<'a>,
    state: DepositionSimulation<'m>,
    enthalpies: Vec<f64>,
    temperatures: Vec<f64>,
    net_powers: Vec<f64>,
    report: ThermalOperatingReport,
}
impl ThermalOperatingPoint<'_, '_, '_> {
    pub fn cell_ids(&self) -> &[CellId] {
        &self.owner.ids
    }
    pub fn specific_enthalpies_j_kg(&self) -> &[f64] {
        &self.enthalpies
    }
    pub fn temperatures_k(&self) -> &[f64] {
        &self.temperatures
    }
    /// Net inward watts in cell_ids() order, not only a summed external balance.
    pub fn net_powers_w(&self) -> &[f64] {
        &self.net_powers
    }
    pub fn operating_time_s(&self) -> f64 {
        self.state.time_s
    }
    pub fn report(&self) -> ThermalOperatingReport {
        self.report
    }

    /// Reuse this exact source-bound operating point for small-signal response.
    /// The harmonic owner's support, resource and residual gates remain active.
    pub fn prepare_thermal_impedance(
        &self,
        budget: ThermalImpedanceBudget,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThermalImpedanceModel, ThermalImpedanceError> {
        self.owner
            .prepare_thermal_impedance(&self.state, budget, checkpoint)
    }
}

impl<'a> ImplicitConductionChain<'a> {
    /// Find a stationary point using the original analytic energy tangent.
    /// The caller supplies an initial guess within every original source domain.
    /// Line search, cancellation, singular tangents and work exhaustion publish
    /// no field. Advection and isolated energy-constrained equilibria are outside
    /// this conduction-only lane. A latent plateau is never regularized.
    pub fn solve_operating_point<'s, 'm>(
        &'s self,
        initial_guess: &DepositionSimulation<'m>,
        controls: ThermalOperatingControls,
        workspace: &mut ImplicitConductionWorkspace,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ThermalOperatingPoint<'s, 'a, 'm>, ThermalOperatingError> {
        checkpoint()?;
        if workspace.cells() != self.ids.len()
            || controls.maximum_residual_evaluations == 0
            || !nonnegative(controls.cell_tolerance_w)
            || !nonnegative(controls.total_tolerance_w)
            || self.advection.is_some()
        {
            return Err(ImplicitConductionError::Invalid(
                "invalid stationary controls/workspace or unsupported advection",
            )
            .into());
        }
        self.admit(initial_guess, &mut checkpoint)?;
        let mut state = initial_guess.clone();
        let mut candidate = initial_guess.clone();
        let mut work = empty_report();
        let (mut endpoint, mut powers) =
            self.operating_trial(&state, controls, &mut work, &mut checkpoint)?;
        self.admit_operating_anchors(&state, &endpoint, &mut checkpoint)?;
        for iteration in 0..=controls.maximum_newton_iterations {
            checkpoint()?;
            let (maximum, total) = residuals(&powers, &mut checkpoint)?;
            if maximum <= controls.cell_tolerance_w && total.abs() <= controls.total_tolerance_w {
                let mut enthalpies = Vec::with_capacity(self.ids.len());
                let mut temperatures = Vec::with_capacity(self.ids.len());
                for (i, &id) in self.ids.iter().enumerate() {
                    if i % 256 == 0 {
                        checkpoint()?;
                    }
                    enthalpies.push(state.cells[&id].specific_h_j_kg);
                    temperatures.push(endpoint.states[i].temperature_k);
                }
                checkpoint()?;
                return Ok(ThermalOperatingPoint {
                    owner: self,
                    state,
                    enthalpies,
                    temperatures,
                    net_powers: powers,
                    report: ThermalOperatingReport {
                        nonlinear_iterations: work.nonlinear_iterations,
                        linear_iterations: work.linear_iterations,
                        residual_evaluations: work.residual_evaluations,
                        maximum_cell_residual_w: maximum,
                        total_residual_w: total,
                    },
                });
            }
            if iteration == controls.maximum_newton_iterations {
                return Err(ThermalOperatingError::NotConverged {
                    maximum_residual_w: maximum,
                });
            }
            // -dP/dE * delta E = P. The zero identity shift is essential:
            // time-stepping to a large dt would perturb the stationary problem.
            self.assemble_tangent(&state, &endpoint, 1.0, 0.0, workspace, &mut checkpoint)?;
            workspace.rhs.copy_from_slice(&powers);
            self.solve_direction(workspace, &mut work, &mut checkpoint)?;
            work.nonlinear_iterations += 1;
            let mut fraction = 1.0;
            let mut accepted = None;
            for _ in 0..=controls.maximum_backtracks {
                checkpoint()?;
                let mut changed = false;
                let mut in_domain = true;
                for (i, &id) in self.ids.iter().enumerate() {
                    checkpoint()?;
                    let cell = &state.cells[&id];
                    let h =
                        cell.specific_h_j_kg + fraction * (workspace.direction[i] / cell.mass_kg);
                    let [lo, hi] = self.laws[&cell.material_id].enthalpy_bounds_j_kg();
                    in_domain &= h.is_finite() && (lo..=hi).contains(&h);
                    changed |= h != cell.specific_h_j_kg;
                    candidate
                        .cells
                        .get_mut(&id)
                        .ok_or(DepositionError::Arithmetic)?
                        .specific_h_j_kg = h;
                }
                if in_domain && changed {
                    match self.operating_trial(&candidate, controls, &mut work, &mut checkpoint) {
                        Ok((trial_endpoint, trial_powers)) => {
                            let (next_maximum, next_total) =
                                residuals(&trial_powers, &mut checkpoint)?;
                            if next_maximum <= (1.0 - 1e-4 * fraction) * maximum
                                || (next_maximum <= controls.cell_tolerance_w
                                    && next_total.abs() <= controls.total_tolerance_w)
                            {
                                match self.admit_operating_anchors(
                                    &candidate,
                                    &trial_endpoint,
                                    &mut checkpoint,
                                ) {
                                    Ok(()) => {
                                        accepted = Some((trial_endpoint, trial_powers));
                                        break;
                                    }
                                    Err(ThermalOperatingError::Unanchored { .. }) => {}
                                    Err(error) => return Err(error),
                                }
                            }
                        }
                        Err(ThermalOperatingError::Original(ImplicitConductionError::Thermal(
                            DepositionError::MaterialDomain | DepositionError::Arithmetic,
                        ))) => {}
                        Err(error) => return Err(error),
                    }
                }
                fraction *= 0.5;
            }
            let Some((next_endpoint, next_powers)) = accepted else {
                return Err(ThermalOperatingError::NotConverged {
                    maximum_residual_w: maximum,
                });
            };
            std::mem::swap(&mut state, &mut candidate);
            endpoint = next_endpoint;
            powers = next_powers;
        }
        unreachable!("bounded stationary loop returns at its limit")
    }

    fn operating_trial(
        &self,
        state: &DepositionSimulation<'_>,
        controls: ThermalOperatingControls,
        work: &mut ImplicitConductionReport,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(Endpoint, Vec<f64>), ThermalOperatingError> {
        checkpoint()?;
        if work.residual_evaluations >= controls.maximum_residual_evaluations {
            return Err(DepositionError::Budget("stationary thermal residual evaluations").into());
        }
        work.residual_evaluations += 1;
        let endpoint = self.endpoint(state, checkpoint)?;
        let powers = self.stationary_heat_rates(state, &endpoint, checkpoint)?;
        Ok((endpoint, powers))
    }

    /// Shared with harmonic preparation; source feedback is credited once on
    /// top of the original paired conduction/boundary heat-rate owner.
    pub(super) fn stationary_heat_rates(
        &self,
        thermal: &DepositionSimulation<'_>,
        endpoint: &Endpoint,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<f64>, ImplicitConductionError> {
        let mut rates = thermal.heat_rates_with_coefficients(
            f64::MAX,
            |id, cell| {
                Ok((
                    thermal.temperature_k(id)?,
                    thermal.material_temperature_slopes[&cell.material_id],
                ))
            },
            &endpoint.boundary_values,
            &endpoint.exchange_values,
            Some(&endpoint.link_values),
            checkpoint,
        )?;
        for (&id, source) in &endpoint.heat_source_states {
            checkpoint()?;
            rates
                .watts
                .get_mut(&id)
                .ok_or(DepositionError::Arithmetic)?
                .add(source.power_w)?;
        }
        if let (Some(source), Some(state)) = (
            &self.coupled_heat_source,
            &endpoint.coupled_heat_source_state,
        ) {
            for (&id, &power) in source.cells.iter().zip(&state.powers_w) {
                checkpoint()?;
                rates
                    .watts
                    .get_mut(&id)
                    .ok_or(DepositionError::Arithmetic)?
                    .add(power)?;
            }
        }
        self.ids
            .iter()
            .map(|id| {
                checkpoint()?;
                rates.watts[id].value()
            })
            .collect::<Result<_, _>>()
            .map_err(Into::into)
    }

    fn admit_operating_anchors(
        &self,
        thermal: &DepositionSimulation<'_>,
        endpoint: &Endpoint,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), ThermalOperatingError> {
        let n = self.ids.len();
        let mut adjacency = vec![Vec::new(); n];
        let mut covered: Vec<f64> = self
            .ids
            .iter()
            .map(|id| {
                thermal.cells[id]
                    .boundary
                    .reservoir
                    .map_or(0.0, |r| r.area_m2)
            })
            .collect();
        for (&id, link) in &thermal.links {
            checkpoint()?;
            let a = self.index[&link.a];
            let b = self.index[&link.b];
            let areas = thermal.contact_areas(link);
            covered[a] = finite(covered[a] + areas[0])?;
            covered[b] = finite(covered[b] + areas[1])?;
            if endpoint.link_values[&id] > 0.0 {
                adjacency[a].push(b);
                adjacency[b].push(a);
            }
        }
        let mut anchored = vec![false; n];
        for (i, &id) in self.ids.iter().enumerate() {
            checkpoint()?;
            let cell = &thermal.cells[&id];
            let area = finite(cell.geometry.surface_area_m2 - covered[i])?;
            if area < 0.0 {
                return Err(DepositionError::Arithmetic.into());
            }
            let ambient = if let Some(exchange) = endpoint.ambient_exchange_states.get(&id) {
                exchange.power_enthalpy_derivative_w_kg_j != 0.0
            } else {
                endpoint
                    .boundary_values
                    .ambient(id, area, cell.boundary.htc_w_m2_k)?
                    > 0.0
            };
            let reservoir = if let Some(exchange) = endpoint.reservoir_exchange_states.get(&id) {
                exchange.power_enthalpy_derivative_w_kg_j != 0.0
            } else {
                endpoint
                    .boundary_values
                    .reservoir_w_k
                    .get(&id)
                    .is_some_and(|&g| g > 0.0)
            };
            anchored[i] = ambient || reservoir;
        }
        let mut seen = vec![false; n];
        let mut pending = Vec::with_capacity(n);
        for start in 0..n {
            checkpoint()?;
            if seen[start] {
                continue;
            }
            pending.clear();
            pending.push(start);
            seen[start] = true;
            let mut has_bath = false;
            while let Some(i) = pending.pop() {
                checkpoint()?;
                has_bath |= anchored[i];
                for &j in &adjacency[i] {
                    if !seen[j] {
                        seen[j] = true;
                        pending.push(j);
                    }
                }
            }
            if !has_bath {
                return Err(ThermalOperatingError::Unanchored {
                    cell: self.ids[start],
                });
            }
        }
        Ok(())
    }
}

fn residuals(
    powers: &[f64],
    checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
) -> Result<(f64, f64), DepositionError> {
    let mut total = Sum::default();
    let mut maximum = 0.0_f64;
    for (i, &power) in powers.iter().enumerate() {
        if i % 256 == 0 {
            checkpoint()?;
        }
        maximum = maximum.max(power.abs());
        total.add(power)?;
    }
    Ok((maximum, total.total()?))
}
