//! Source-backed walls on the existing steady gas enthalpy and acoustic field.
//!
//! Wall conductivity, geometry and an explicitly declared gas/contact film
//! determine each thermal substep's reservoir flux. The same gas solution
//! supplies density, calorics, viscosity and acoustics. This does not solve
//! finite wall storage or infer a convection correlation from a material name.

use super::network::{SteadyGasNetwork, SteadyGasNetworkSolution};
use super::{
    MAX_STEADY_GAS_CELLS, SteadyAxialGas, SteadyGasError, SteadyGasSolution, SteadyGasSolve, bad,
    poll,
};
use crate::thermal::conduction::{
    LayerConductionConfig, ReservoirLayerGeometry, ReservoirLayerResponse,
    source::ResolvedConductionLayer,
};
use fs_thermal::enthalpy::{
    CellId, DepositionError, DepositionObserver, DepositionSimulation, Event,
};
use std::collections::BTreeMap;

/// One original conductivity source bound to one gas control volume.
/// Geometry and source remain immutable during the steady solve.
#[derive(Clone, Copy, Debug)]
pub struct MaterialGasWall<'a, 'm> {
    cell: usize,
    layer: &'a ResolvedConductionLayer<'m>,
    geometry: ReservoirLayerGeometry,
    film_resistance_k_w: f64,
}

impl<'a, 'm> MaterialGasWall<'a, 'm> {
    /// The external material face uses the cell's prescribed wall temperature.
    /// Film resistance is TOTAL K/W between well-mixed gas and inner material
    /// face. Zero explicitly assumes perfect contact. A measured or separately
    /// resolved film remains independent data; this API invents no convection.
    ///
    /// Binding replaces this cell's supplied wall conductance, never adds a
    /// second copy of its heat. Other cells retain their original wall inputs.
    pub fn new(
        cell: usize,
        layer: &'a ResolvedConductionLayer<'m>,
        geometry: ReservoirLayerGeometry,
        film_resistance_k_w: f64,
    ) -> Result<Self, SteadyGasError> {
        if !film_resistance_k_w.is_finite() || film_resistance_k_w < 0.0 {
            return Err(bad(
                "material gas wall film resistance must be finite and nonnegative",
            ));
        }
        Ok(Self {
            cell,
            layer,
            geometry,
            film_resistance_k_w,
        })
    }
    pub const fn cell(&self) -> usize {
        self.cell
    }
    /// Retains all original material/source-line and pure-phase admission data.
    pub const fn layer(&self) -> &'a ResolvedConductionLayer<'m> {
        self.layer
    }

    /// Actual nonlinear wall flux at the stated gas/reservoir temperatures.
    pub fn evaluate(
        &self,
        temperatures_k: [f64; 2],
        config: LayerConductionConfig,
        checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<ReservoirLayerResponse, SteadyGasError> {
        Ok(self.geometry.evaluate(
            self.layer.conductivity_curve(),
            temperatures_k,
            self.film_resistance_k_w,
            config,
            checkpoint,
        )?)
    }
}

#[derive(Clone)]
struct PreparedWall<'a, 'm> {
    wall: MaterialGasWall<'a, 'm>,
    maximum_conductance_w_k: f64,
}

#[derive(Clone)]
struct WallObserver<'a, 'm> {
    walls: BTreeMap<CellId, PreparedWall<'a, 'm>>,
    controls: LayerConductionConfig,
    extra_conductance_w_k: BTreeMap<CellId, f64>,
}

impl<'a, 'm> WallObserver<'a, 'm> {
    fn new(
        network: &SteadyGasNetwork,
        walls: &[MaterialGasWall<'a, 'm>],
        controls: LayerConductionConfig,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, SteadyGasError> {
        poll(&mut cancelled)?;
        if walls.len() > network.cells.len() || walls.len() > MAX_STEADY_GAS_CELLS {
            return Err(bad("material gas walls exceed the cell inventory"));
        }
        let mut prepared = BTreeMap::new();
        for wall in walls {
            poll(&mut cancelled)?;
            let cell = network
                .cells
                .get(wall.cell)
                .ok_or(bad("material gas wall names an inactive cell"))?;
            // A smaller wall leaves the rest of the existing surface insulated.
            // The tolerance covers only multiplication-order roundoff in geometry.
            if wall.geometry.surface_area_m2()
                > cell.geometry.surface_area_m2 * (1.0 + 16.0 * f64::EPSILON)
            {
                return Err(bad("material gas wall exceeds the cell's exposed surface"));
            }
            wall.evaluate(
                [cell.initial_temperature_k, cell.wall_temperature_k],
                controls,
                || poll(&mut cancelled),
            )?;
            let hottest_conductivity_knot = wall
                .layer
                .conductivity_curve()
                .knots()
                .iter()
                .max_by(|a, b| a[1].total_cmp(&b[1]))
                .expect("admitted conductivity chart");
            let maximum_conductance_w_k = wall
                .evaluate([hottest_conductivity_knot[0]; 2], controls, || {
                    poll(&mut cancelled)
                })?
                .conductance_w_k;
            if prepared
                .insert(
                    CellId(wall.cell as u64),
                    PreparedWall {
                        wall: *wall,
                        maximum_conductance_w_k,
                    },
                )
                .is_some()
            {
                return Err(bad("duplicate material gas wall cell"));
            }
        }
        Ok(Self {
            walls: prepared,
            controls,
            extra_conductance_w_k: BTreeMap::new(),
        })
    }

    fn validate_endpoint(
        &self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), SteadyGasError> {
        for (&id, prepared) in &self.walls {
            checkpoint()?;
            // Birth events build the initial inventory one cell at a time.
            if let Some(cell) = thermal.cells().get(&id) {
                let reservoir = cell
                    .boundary
                    .reservoir
                    .ok_or(bad("material gas wall lost its physical reservoir"))?;
                prepared.wall.evaluate(
                    [thermal.temperature_k(id)?, reservoir.temperature_k],
                    self.controls,
                    &mut *checkpoint,
                )?;
            }
        }
        Ok(())
    }
}

impl DepositionObserver for WallObserver<'_, '_> {
    type Error = SteadyGasError;

    fn reservoir_conductances(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<BTreeMap<CellId, f64>, Self::Error> {
        self.extra_conductance_w_k.clear();
        let mut conductances = BTreeMap::new();
        for (&id, prepared) in &self.walls {
            checkpoint()?;
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(bad("material gas wall names an inactive thermal cell"))?;
            let reservoir = cell
                .boundary
                .reservoir
                .ok_or(bad("material gas wall needs the existing reservoir"))?;
            let response = prepared.wall.evaluate(
                [thermal.temperature_k(id)?, reservoir.temperature_k],
                self.controls,
                &mut *checkpoint,
            )?;
            conductances.insert(id, response.conductance_w_k);
            // The native CFL already counts the secant G. Add only the excess
            // of the maximum possible material tangent across the whole source
            // span; advection then composes its own rate on this same bound.
            self.extra_conductance_w_k.insert(
                id,
                (prepared.maximum_conductance_w_k - response.conductance_w_k).max(0.0),
            );
        }
        Ok(conductances)
    }

    fn limit_timestep_s(
        &mut self,
        thermal: &DepositionSimulation<'_>,
        proposed_s: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<f64, Self::Error> {
        Ok(thermal.limit_timestep_for_additional_conductances(
            proposed_s,
            &self.extra_conductance_w_k,
            checkpoint,
        )?)
    }

    fn event(&mut self, _: &Event, thermal: &DepositionSimulation<'_>) -> Result<(), Self::Error> {
        self.validate_endpoint(thermal, &mut || Ok(()))
    }

    fn step(
        &mut self,
        _: f64,
        _: f64,
        _: &BTreeMap<CellId, f64>,
        thermal: &DepositionSimulation<'_>,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Self::Error> {
        self.validate_endpoint(thermal, checkpoint)
    }
}

impl SteadyGasNetwork {
    /// Solve the SAME network with source-bound nonlinear wall heat. Geometry,
    /// species, pressure, calorics, advection and the energy audit keep their
    /// original owners. No wall flux is applied as a second external source.
    /// Failure/cancellation publishes no gas candidate and changes no input.
    pub fn solve_with_material_walls(
        &self,
        walls: &[MaterialGasWall<'_, '_>],
        wall_controls: LayerConductionConfig,
        controls: SteadyGasSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SteadyGasNetworkSolution, SteadyGasError> {
        let observer = WallObserver::new(self, walls, wall_controls, &mut cancelled)?;
        let mut solution = self.solve_with_observer(controls, observer, &mut cancelled)?;
        // For the monotone film/material residual, dF/dQ >= 1. Its absolute
        // residual therefore bounds the remaining physical flux error. Include
        // it in the published balance; loose wall tolerances cannot masquerade
        // as a converged gas field.
        for wall in walls {
            let response = wall.evaluate(
                [
                    solution.states()[wall.cell].state().temperature,
                    self.cells[wall.cell].wall_temperature_k,
                ],
                wall_controls,
                || poll(&mut cancelled),
            )?;
            let bound =
                solution.cell_power_residual_bounds_w[wall.cell] + response.heat_residual_w.abs();
            if !bound.is_finite() || bound > controls.steady_power_tolerance_w {
                return Err(SteadyGasError::PublishedBalance {
                    cell: wall.cell,
                    residual_bound_w: bound,
                });
            }
            solution.cell_power_residual_bounds_w[wall.cell] = bound;
        }
        poll(&mut cancelled)?;
        Ok(solution)
    }
}

impl SteadyAxialGas {
    /// Preserve the original axial-to-network lowering and publish the SAME
    /// AxialGasField consumed by regional ducts, reeds and radiation. Changing
    /// only a wall source changes heat transfer, local gas state and acoustics.
    pub fn solve_with_material_walls(
        &self,
        walls: &[MaterialGasWall<'_, '_>],
        wall_controls: LayerConductionConfig,
        controls: SteadyGasSolve,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<SteadyGasSolution, SteadyGasError> {
        poll(&mut cancelled)?;
        let network = self.to_network()?;
        let solution =
            network.solve_with_material_walls(walls, wall_controls, controls, &mut cancelled)?;
        self.finish_network_solution(&network, solution, controls, cancelled)
    }
}

#[cfg(test)]
mod tests;
