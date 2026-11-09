//! Conserved-storage equilibrium of closed native conduction components.
//!
//! Each connected component retains its own original sum m*e. The historical
//! h-labelled state may store enthalpy, fixed-volume internal energy, or a
//! declared sensible-storage coordinate; this owner never converts those bases.
//! Common temperatures use the original monotone caloric laws and sharp latent
//! intervals. No bath, heat source, advection, elapsed time or work is invented.

pub mod caloric;
use super::*;
pub use caloric::{CaloricEquilibriumControls, CaloricEquilibriumError};
use caloric::{CaloricInventory, CaloricState, equilibrate_caloric_inventory};

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct IsolatedThermalControls {
    /// This evaluation budget is shared by all disconnected components.
    pub caloric: CaloricEquilibriumControls,
    /// Independent final original heat-rate tolerance in EVERY cell, W.
    pub cell_power_tolerance_w: f64,
}
impl Default for IsolatedThermalControls {
    fn default() -> Self {
        Self {
            caloric: CaloricEquilibriumControls::default(),
            cell_power_tolerance_w: 1e-9,
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
pub enum IsolatedThermalError {
    Original(ImplicitConductionError),
    Caloric(CaloricEquilibriumError<DepositionError>),
    Unsupported(&'static str),
    /// An initially conducting connection became insulating, or a connection
    /// supplied no conducting path from which to define a closed inventory.
    ConductingTopology {
        link: u64,
    },
    Residual {
        cell: CellId,
        power_w: f64,
    },
}
impl From<ImplicitConductionError> for IsolatedThermalError {
    fn from(value: ImplicitConductionError) -> Self {
        Self::Original(value)
    }
}
impl From<DepositionError> for IsolatedThermalError {
    fn from(value: DepositionError) -> Self {
        Self::Original(value.into())
    }
}
impl core::fmt::Display for IsolatedThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "closed thermal equilibrium: {self:?}")
    }
}
impl std::error::Error for IsolatedThermalError {}

#[derive(Clone, Debug, PartialEq)]
pub struct IsolatedComponentReport {
    /// Stable original cells; no inventory passes between distinct components.
    pub cells: Vec<CellId>,
    pub energy_residual_j: f64,
    pub temperature_spread_k: f64,
    pub temperature_iterations: usize,
    pub caloric_evaluations: usize,
    /// True includes an originally equilibrated latent allocation, kept exactly.
    pub retained_initial_state: bool,
}

/// Read-only equilibrium field and a separate conservation row per component.
pub struct IsolatedThermalEquilibrium<'s, 'a, 'm> {
    field: ThermalOperatingPoint<'s, 'a, 'm>,
    components: Vec<IsolatedComponentReport>,
    caloric_evaluations: usize,
}
impl<'s, 'a, 'm> IsolatedThermalEquilibrium<'s, 'a, 'm> {
    /// Same immutable material/geometry-bound field as driven operating points.
    /// Its original harmonic preparation may still refuse an unanchored DC mode.
    pub fn field(&self) -> &ThermalOperatingPoint<'s, 'a, 'm> {
        &self.field
    }
    pub fn components(&self) -> &[IsolatedComponentReport] {
        &self.components
    }
    pub fn caloric_evaluations(&self) -> usize {
        self.caloric_evaluations
    }
}

impl<'a> ImplicitConductionChain<'a> {
    /// Solve closed components without marching through a fictitious time.
    ///
    /// Sources, baths, advection and custom flux links are outside this lane.
    /// Native geometric Fourier links and passive contact resistances remain
    /// original. Every authored link must conduct initially and at the result;
    /// changing insulating topology requires an explicit physical owner. Truly
    /// disconnected cells/components are supported by the graph entry point.
    ///
    /// An already equilibrated component keeps every original specific energy,
    /// including its latent allocations. Otherwise one unresolved latent parcel
    /// is determined by its energy row; several independent latent allocations
    /// refuse. No phase distribution, entropy production amount or dynamical
    /// reachability is inferred from these equilibrium conditions.
    pub fn solve_isolated_equilibrium<'s, 'm>(
        &'s self,
        initial: &DepositionSimulation<'m>,
        controls: IsolatedThermalControls,
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<IsolatedThermalEquilibrium<'s, 'a, 'm>, IsolatedThermalError> {
        checkpoint()?;
        if !nonnegative(controls.cell_power_tolerance_w) {
            return Err(IsolatedThermalError::Unsupported(
                "finite nonnegative watt tolerance required",
            ));
        }
        self.admit(initial, &mut checkpoint)?;
        if self.advection.is_some()
            || !self.heat_sources.is_empty()
            || self.coupled_heat_source.is_some()
            || !self.link_laws.is_empty()
            || !self.reservoir_laws.is_empty()
            || !self.ambient_laws.is_empty()
            || !self.reservoir_exchanges.is_empty()
            || !self.ambient_exchanges.is_empty()
            || !self.reservoir_histories.is_empty()
            || !self.ambient_histories.is_empty()
            || !self.conduction.reservoirs.is_empty()
            || !self.conduction.ambient_faces.is_empty()
        {
            return Err(IsolatedThermalError::Unsupported(
                "closed equilibrium requires no bath, source, advection or custom flux binding",
            ));
        }
        for cell in initial.cells.values() {
            checkpoint()?;
            if cell.boundary.source_w != 0.0
                || cell.boundary.htc_w_m2_k != 0.0
                || cell.boundary.reservoir.is_some()
            {
                return Err(IsolatedThermalError::Unsupported(
                    "closed equilibrium requires insulated source-free cells",
                ));
            }
        }
        let initial_endpoint = self.endpoint(initial, &mut checkpoint)?;
        let components = self.closed_components(initial, &initial_endpoint, &mut checkpoint)?;
        let mut state = initial.clone();
        let mut reports = Vec::with_capacity(components.len());
        let mut used = 0usize;
        let mut iterations = 0usize;
        for component in &components {
            checkpoint()?;
            let mut parcels = Vec::with_capacity(component.len());
            for &index in component {
                let cell = &initial.cells[&self.ids[index]];
                parcels.push(CaloricInventory {
                    mass_kg: cell.mass_kg,
                    initial_specific_energy_j_kg: cell.specific_h_j_kg,
                    specific_energy_bounds_j_kg: self.laws[&cell.material_id]
                        .enthalpy_bounds_j_kg(),
                });
            }
            let remaining = controls
                .caloric
                .maximum_caloric_evaluations
                .checked_sub(used)
                .filter(|&remaining| remaining != 0)
                .ok_or(DepositionError::Budget("closed caloric evaluations"))?;
            let solved = equilibrate_caloric_inventory(
                &parcels,
                CaloricEquilibriumControls {
                    maximum_caloric_evaluations: remaining,
                    ..controls.caloric
                },
                |i, energy| {
                    let cell = &initial.cells[&self.ids[component[i]]];
                    let law = self.laws[&cell.material_id];
                    // The exact original caloric owner is authoritative. Its
                    // analytic tangent accelerates only the scalar temperature
                    // solve, and must agree with the source's admitted slope.
                    let temperature = law.caloric_law().temperature_k(energy)?;
                    let tangent = law.state_at_specific_enthalpy(energy)?;
                    if temperature != tangent.temperature_k
                        || tangent.temperature_slope_kg_k_j
                            > initial.material_temperature_slopes[&cell.material_id]
                    {
                        return Err(DepositionError::MaterialDomain);
                    }
                    Ok(CaloricState {
                        temperature_k: temperature,
                        temperature_slope_kg_k_j: tangent.temperature_slope_kg_k_j,
                    })
                },
                &mut checkpoint,
            )
            .map_err(IsolatedThermalError::Caloric)?;
            used = used
                .checked_add(solved.caloric_evaluations)
                .ok_or(DepositionError::Arithmetic)?;
            iterations = iterations
                .checked_add(solved.temperature_iterations)
                .ok_or(DepositionError::Arithmetic)?;
            let ids: Vec<_> = component.iter().map(|&i| self.ids[i]).collect();
            for (&id, &energy) in ids.iter().zip(&solved.specific_energies_j_kg) {
                checkpoint()?;
                let cell = state
                    .cells
                    .get_mut(&id)
                    .ok_or(DepositionError::Arithmetic)?;
                cell.specific_h_j_kg = energy;
                self.laws[&cell.material_id]
                    .caloric_law()
                    .validate_cell(cell)?;
            }
            reports.push(IsolatedComponentReport {
                cells: ids,
                energy_residual_j: solved.energy_residual_j,
                temperature_spread_k: solved.temperature_spread_k,
                temperature_iterations: solved.temperature_iterations,
                caloric_evaluations: solved.caloric_evaluations,
                retained_initial_state: solved.retained_initial_state,
            });
        }
        self.admit(&state, &mut checkpoint)?;
        let endpoint = self.endpoint(&state, &mut checkpoint)?;
        if self.closed_components(&state, &endpoint, &mut checkpoint)? != components {
            return Err(IsolatedThermalError::Unsupported(
                "closed conducting component inventory changed",
            ));
        }
        let powers = self.stationary_heat_rates(&state, &endpoint, &mut checkpoint)?;
        for (&id, &power) in self.ids.iter().zip(&powers) {
            checkpoint()?;
            if power.abs() > controls.cell_power_tolerance_w {
                return Err(IsolatedThermalError::Residual {
                    cell: id,
                    power_w: power,
                });
            }
        }
        let (maximum, total) = residuals(&powers, &mut checkpoint)?;
        let enthalpies = self
            .ids
            .iter()
            .map(|id| state.cells[id].specific_h_j_kg)
            .collect();
        let temperatures = endpoint
            .states
            .iter()
            .map(|value| value.temperature_k)
            .collect();
        checkpoint()?;
        Ok(IsolatedThermalEquilibrium {
            field: ThermalOperatingPoint {
                owner: self,
                state,
                enthalpies,
                temperatures,
                net_powers: powers,
                report: ThermalOperatingReport {
                    nonlinear_iterations: iterations,
                    linear_iterations: 0,
                    residual_evaluations: 1,
                    maximum_cell_residual_w: maximum,
                    total_residual_w: total,
                },
            },
            components: reports,
            caloric_evaluations: used,
        })
    }

    fn closed_components(
        &self,
        state: &DepositionSimulation<'_>,
        endpoint: &Endpoint,
        checkpoint: &mut impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<Vec<usize>>, IsolatedThermalError> {
        let mut adjacency = vec![Vec::new(); self.ids.len()];
        for (&id, link) in &state.links {
            checkpoint()?;
            if !endpoint.link_values.get(&id).is_some_and(|&g| g > 0.0) {
                return Err(IsolatedThermalError::ConductingTopology { link: id });
            }
            let a = self.index[&link.a];
            let b = self.index[&link.b];
            adjacency[a].push(b);
            adjacency[b].push(a);
        }
        let mut seen = vec![false; self.ids.len()];
        let mut components = Vec::new();
        for start in 0..self.ids.len() {
            checkpoint()?;
            if seen[start] {
                continue;
            }
            let mut pending = vec![start];
            let mut component = Vec::new();
            seen[start] = true;
            while let Some(index) = pending.pop() {
                checkpoint()?;
                component.push(index);
                for &neighbor in &adjacency[index] {
                    if !seen[neighbor] {
                        seen[neighbor] = true;
                        pending.push(neighbor);
                    }
                }
            }
            component.sort_unstable();
            components.push(component);
        }
        Ok(components)
    }
}
