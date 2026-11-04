//! Generalized-Newtonian film pressure in the existing joint time equation.
//!
//! fs-tribo owns the complete shear-dependent pressure solve and its analytic
//! tangent. This binding supplies the original configuration and actual
//! discrete-gradient velocity, and adds the resulting work-conjugate force
//! once. The original fs-phs owner advances motion and accounts dissipation.
//! Laws and numerical budgets remain immutable throughout mechanical trials.

use fs_tribo::resistive_film::MAX_CELLS;
pub use fs_tribo::resistive_film::power_law::{
    PowerLawControls, PowerLawFilm, PowerLawFilmReport, RegularizedPowerLaw,
};

use super::{ImpactError, ImpactSystem, MAX_IMPACT_MODES, invalid};

impl ImpactSystem {
    /// Number of original generalized-Newtonian pressure islands.
    /// Insertion indices are independent of normal and sliding film indices.
    pub fn power_law_film_count(&self) -> usize {
        self.power_law_films.len()
    }

    /// Complete accepted constitutive owner, including its SI channel laws
    /// and explicit nonlinear flow/work admission budgets.
    pub fn power_law_film(&self, index: usize) -> Option<&PowerLawFilm> {
        self.power_law_films.get(index)
    }

    /// Attach shear-thinning or shear-thickening pressure resistance to the
    /// original mechanical force, directional tangent and acceptance path.
    ///
    /// The film covers every mechanical mode. All resistance-only fluid
    /// families share four pressure islands; compressible-gas images remain
    /// exclusive. Mapped support memory retains its original recomposition
    /// restriction. Neither body state nor material history is reset.
    ///
    /// Reference and prepared execution use the complete nonlinear owner,
    /// never its Newtonian-core geometry view. Each solve keeps its explicit
    /// bounded work and closure controls. Original mechanical cancellation
    /// and final endpoint admission remain authoritative.
    ///
    /// Its physical work can be routed through the existing material heating
    /// owner with `PowerLawFilmHeatRoute`. This stationary-wall, full-film
    /// reduction retains frozen laws; it adds no sliding entrainment,
    /// temperature feedback, cavitation or fluid inertia.
    pub fn with_power_law_film(mut self, film: PowerLawFilm) -> Result<Self, ImpactError> {
        self.require_recomposition_without_mapped_memory()?;
        let geometry = film.pressure_operator();
        if self.gas_film.is_some()
            || geometry.port_count() != self.modes
            || self.incompressible_film_count() >= 4
        {
            return Err(invalid(
                "power-law film mode count, duplicate gas image or four-island ceiling failed",
            ));
        }
        let mut q = [0.0; MAX_IMPACT_MODES];
        for (i, value) in q[..self.modes].iter_mut().enumerate() {
            *value = self.x[2 * i];
        }
        let mut force = [0.0; MAX_IMPACT_MODES];
        let mut pressure = [0.0; MAX_CELLS];
        film.evaluate_into(
            &q[..self.modes],
            &[0.0; MAX_IMPACT_MODES][..self.modes],
            &mut force[..self.modes],
            &mut pressure[..geometry.cell_count()],
        )
        .map_err(|error| invalid(error.0))?;
        self.power_law_films.push(film);
        Ok(self)
    }
}

#[cfg(test)]
mod tests;
