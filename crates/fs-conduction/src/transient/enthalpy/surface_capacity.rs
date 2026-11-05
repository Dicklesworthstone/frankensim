//! Sensible, perfectly bonded film storage on the original host enthalpy.
//!
//! Each row stores E_i(h_i)=m_i*h_i+C_i*T_i(h_i), plus each nonlinear film's
//! original integral of C(T) through that row's temperature. Constant C_i is
//! the exact original P1 surface-capacity row moment.
//! The host chart still owns all latent heat. A fixed host temperature stores
//! no additional sensible film heat, even while the host phase fraction moves.
//! This is not a film latent-heat law or a separately meshed material inventory.

use std::collections::BTreeSet;

use fs_math::dd::Dd;

use super::*;
use crate::material::TemperatureSpan;
use crate::transient::{MaterialHeatCapacity, backward_euler::SurfaceHeatCapacity};

mod nonlinear;
pub use nonlinear::SurfaceHeatStorage;

#[derive(Debug)]
pub(super) struct SurfaceStorage<'c> {
    nodal_j_per_k: Vec<f64>,
    temperature_spans: Vec<TemperatureSpan>,
    materials: Vec<MaterialHeatCapacity>,
    nonlinear: Option<nonlinear::NonlinearSurfaceStorage<'c>>,
}

impl SurfaceStorage<'_> {
    fn new(vertices: usize) -> Self {
        Self {
            nodal_j_per_k: vec![0.0; vertices],
            temperature_spans: vec![TemperatureSpan::Unbounded; vertices],
            materials: Vec::new(),
            nonlinear: None,
        }
    }
}

impl EnthalpyBackwardEuler<'_, '_> {
    /// Add the existing perfectly bonded, sensible film capacity to the
    /// original enthalpy storage. No host reference mass or caloric chart changes.
    ///
    /// The shared row moments give C_i=integral(C_s*N_i dA), in J/K.
    /// Every physical balance uses m_i*(h_new-h_old)+C_i*(T_new-T_old);
    /// Newton uses m_i+C_i*dT/dh. The host's sharp latent plateau remains
    /// exact, with zero sensible film increment at constant temperature.
    ///
    /// Backward Euler, SDIRK2, adaptive audits, local contact/radiation owners,
    /// spectral heating and existing history/endpoint adjoints retain this
    /// storage. Film coefficients, geometry and material sources stay fixed;
    /// their control derivatives and film phase transitions are not introduced.
    ///
    /// A face may occur once per call. Further calls add explicitly declared
    /// layers. Empty/zero layers add no storage or source restriction. Avoid
    /// double-counting a coating already represented by volume elements.
    /// Original material source domains constrain every receiving vertex at
    /// both history and candidate temperatures. A refusal publishes no owner.
    pub fn with_surface_heat_capacity(
        mut self,
        cx: &Cx<'_>,
        surfaces: &[SurfaceHeatCapacity<'_>],
    ) -> Result<Self, EnthalpyError> {
        poll(cx, 0)?;
        let n = self.masses.len();
        let mut seen = BTreeSet::new();
        for (index, &surface) in surfaces.iter().enumerate() {
            if index % ASSEMBLY_TILE == 0 {
                poll(cx, index)?;
            }
            let (vertices, moments) = surface.nodal_moments(self.mesh)?;
            if !seen.insert(surface.boundary_face()) {
                return Err(EnthalpyError::InvalidInput(
                    "surface capacity face is declared more than once",
                ));
            }
            if moments.iter().all(|&value| value == 0.0) {
                continue;
            }
            // Build support before mutably borrowing the retained film owner.
            // A source need only overlap the host chart; every actual state is
            // checked later, without widening either original validity domain.
            let host_spans = vertices.map(|vertex| {
                let knots = self.curve_for_vertex(vertex).knots();
                TemperatureSpan::Sampled {
                    low: knots[0].temperature_k,
                    high: knots[knots.len() - 1].temperature_k,
                }
            });
            let prepared = self.surface_capacity.get_or_insert_with(|| SurfaceStorage::new(n));
            for (local, &vertex) in vertices.iter().enumerate() {
                prepared.nodal_j_per_k[vertex] =
                    finite(prepared.nodal_j_per_k[vertex] + moments[local])?;
                if let Some(material) = surface.material() {
                    let span = prepared.temperature_spans[vertex]
                        .intersect(material.temperature_span())
                        .intersect(host_spans[local]);
                    if matches!(span, TemperatureSpan::Sampled { low, high } if low > high) {
                        return Err(EnthalpyError::InvalidInput(
                            "host chart and film capacity sources have disjoint temperature support",
                        ));
                    }
                    prepared.temperature_spans[vertex] = span;
                }
            }
            if let Some(material) = surface.material()
                && !prepared.materials.contains(material)
            {
                prepared.materials.push(material.clone());
            }
        }
        poll(cx, surfaces.len())?;
        Ok(self)
    }

    /// Constant film capacities in original vertex order, J/K. Nonlinear
    /// inventories returned by `surface_heat_storage` are additional storage.
    /// None means no active film and preserves the original storage path.
    #[must_use]
    pub fn surface_nodal_heat_capacity_j_per_k(&self) -> Option<&[f64]> {
        self.surface_capacity
            .as_ref()
            .map(|film| film.nodal_j_per_k.as_slice())
    }

    /// Original validity-wide constant density/Cp film bundles, without
    /// assigning sourced authority to explicitly declared coefficients.
    #[must_use]
    pub fn surface_material_capacities(&self) -> &[MaterialHeatCapacity] {
        self.surface_capacity
            .as_ref()
            .map_or(&[], |film| film.materials.as_slice())
    }

    pub(super) fn has_film_storage(&self, vertex: usize) -> bool {
        self.surface_capacity.as_ref().is_some_and(|film| {
            film.nodal_j_per_k[vertex] > 0.0
                || film.nonlinear.as_ref().is_some_and(|laws| laws.at_vertex(vertex).next().is_some())
        })
    }

    // A support edge is exact; a tolerance would erase a valid neighborhood.
    #[allow(clippy::float_cmp)]
    pub(super) fn film_temperature_endpoint(&self, vertex: usize, temperature: f64) -> bool {
        self.surface_capacity.as_ref().is_some_and(|film| {
            matches!(film.temperature_spans[vertex], TemperatureSpan::Sampled { low, high }
                if temperature == low || temperature == high)
        })
    }

    pub(super) fn check_film_temperature(
        &self,
        vertex: usize,
        temperature: f64,
    ) -> Result<(), EnthalpyError> {
        if let Some(film) = &self.surface_capacity {
            film.temperature_spans[vertex].check(temperature)?;
        }
        Ok(())
    }

    pub(super) fn storage_derivative(
        &self,
        vertex: usize,
        temperature: f64,
        temperature_derivative: f64,
    ) -> Result<f64, EnthalpyError> {
        if !self.has_film_storage(vertex) {
            return Ok(self.masses[vertex]);
        }
        self.check_film_temperature(vertex, temperature)?;
        // An interior host plateau has exactly zero sensible response. The
        // source temperature is still checked, but no capacity can add latent heat.
        if temperature_derivative == 0.0 { return Ok(self.masses[vertex]); }
        let film = self.surface_capacity.as_ref().expect("active film storage");
        let mut derivative = finite(self.masses[vertex]
            + finite(film.nodal_j_per_k[vertex] * temperature_derivative)?)?;
        if let Some(laws) = &film.nonlinear {
            for law in laws.at_vertex(vertex) {
                let capacity = law.heat_capacity_j_k(temperature)
                    .map_err(|source| EnthalpyError::FilmStorage { vertex, source })?;
                derivative = finite(derivative + finite((capacity / 3.0) * temperature_derivative)?)?;
            }
        }
        Ok(derivative)
    }

    fn film_temperature_change(
        &self,
        vertex: usize,
        old_h: f64,
        temperature: f64,
    ) -> Result<(f64, f64), EnthalpyError> {
        let old_temperature = self
            .curve_for_vertex(vertex)
            .values_at_specific_enthalpy(old_h)
            .map_err(|source| EnthalpyError::Phase { vertex, source })?
            .temperature_k;
        self.check_film_temperature(vertex, old_temperature)?;
        self.check_film_temperature(vertex, temperature)?;
        Ok((old_temperature, temperature))
    }

    pub(super) fn storage_change_j(
        &self,
        vertex: usize,
        new_h: f64,
        old_h: f64,
        temperature: f64,
    ) -> Result<f64, EnthalpyError> {
        let bulk = finite(self.masses[vertex] * (new_h - old_h))?;
        if !self.has_film_storage(vertex) {
            return Ok(bulk);
        }
        let (old_temperature, temperature) =
            self.film_temperature_change(vertex, old_h, temperature)?;
        let film = self.surface_capacity.as_ref().expect("active film storage");
        let mut change = finite(bulk + finite(film.nodal_j_per_k[vertex]
            * (temperature - old_temperature))?)?;
        if let Some(laws) = &film.nonlinear {
            for law in laws.at_vertex(vertex) {
                let energy = law.energy_change_j(old_temperature, temperature)
                    .map_err(|source| EnthalpyError::FilmStorage { vertex, source })?;
                change = finite(change + energy / 3.0)?;
            }
        }
        Ok(change)
    }

    // The genuine SDIRK carry and spectral-source accounts retain their
    // existing double-double subtraction, including the additional film term.
    pub(super) fn storage_change_dd(
        &self,
        vertex: usize,
        new_h: f64,
        old_h: f64,
        temperature: f64,
    ) -> Result<Dd, EnthalpyError> {
        use super::interference::{checked_dd, product};
        let bulk = product(
            Dd::from_f64(self.masses[vertex]),
            checked_dd(Dd::from_f64(new_h) - Dd::from_f64(old_h))?,
        )?;
        if !self.has_film_storage(vertex) {
            return Ok(bulk);
        }
        let (old_temperature, temperature) =
            self.film_temperature_change(vertex, old_h, temperature)?;
        let film = self.surface_capacity.as_ref().expect("active film storage");
        let mut change = checked_dd(bulk + product(
            Dd::from_f64(film.nodal_j_per_k[vertex]),
            checked_dd(Dd::from_f64(temperature) - Dd::from_f64(old_temperature))?,
        )?)?;
        if let Some(laws) = &film.nonlinear {
            for law in laws.at_vertex(vertex) {
                // The native relative integral avoids subtracting large stored
                // primitives; retain its result in the existing stage accumulator.
                let energy = law.energy_change_j(old_temperature, temperature)
                    .map_err(|source| EnthalpyError::FilmStorage { vertex, source })?;
                change = checked_dd(change + Dd::from_f64(energy) / Dd::from_f64(3.0))?;
            }
        }
        Ok(change)
    }
}
