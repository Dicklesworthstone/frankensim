//! Pressure-fed wall friction in the original joint mechanical equation.
//!
//! The film owns pressure, entrainment, traction and nonnegative viscous heat.
//! Reservoir work is an external supply, including return flow into a supply.
//! The same fs-phs step advances mechanics once and accounts for both powers;
//! a pressurized film is never passed off as a passive resistance.
//! Reservoir pressures, wall rows and layer profiles are fixed during a step.

use super::{ImpactError, ImpactSystem, MAX_IMPACT_MODES, invalid};
use fs_phs::NonlinearPortPower;
pub use fs_tribo::resistive_film::shear::PressureDrivenShearedFilm;
use fs_tribo::resistive_film::{MAX_CELLS, MAX_CHANNELS};

impl ImpactSystem {
    /// Number of pressure-fed sliding films, in insertion order.
    /// These indices are independent of passive sliding and squeeze indices.
    pub fn driven_sliding_film_count(&self) -> usize {
        self.driven_sliding_films.len()
    }

    /// Original pressure-fed operator, retaining its fixed pressure supplies.
    pub fn driven_sliding_film(&self, index: usize) -> Option<&PressureDrivenShearedFilm> {
        self.driven_sliding_films.get(index)
    }

    /// Attach pressure supply, wall entrainment and complete shear traction to
    /// the existing mechanical clock without resetting accepted motion/history.
    ///
    /// The same four-island limit covers passive and supplied incompressible
    /// films. Every velocity row covers all mechanical modes. Gas-film images
    /// and mapped support memory keep their existing exclusion. The film's
    /// fixed pressure supplies are external work; their inventories and thermal
    /// destinations remain with the caller. Stationary and resonance analysis
    /// currently refuse this supplied image explicitly.
    pub fn with_driven_sliding_film(
        mut self,
        film: PressureDrivenShearedFilm,
    ) -> Result<Self, ImpactError> {
        self.require_recomposition_without_mapped_memory()?;
        if self.gas_film.is_some()
            || film.port_count() != self.modes
            || self.incompressible_film_count() >= 4
        {
            return Err(invalid(
                "pressure-fed sliding film mode count, gas image or four-island ceiling failed",
            ));
        }
        let mut q = [0.0; MAX_IMPACT_MODES];
        for (i, value) in q[..self.modes].iter_mut().enumerate() {
            *value = self.x[2 * i];
        }
        let mut force = [0.0; MAX_IMPACT_MODES];
        let mut pressure = [0.0; MAX_CELLS];
        film.evaluate_force_into(
            &q[..self.modes],
            &[0.0; MAX_IMPACT_MODES][..self.modes],
            &mut force[..self.modes],
            &mut pressure[..film.cell_count()],
        )
        .map_err(|error| invalid(error.0))?;
        self.driven_sliding_films.push(film);
        Ok(self)
    }

    // Query the existing passive owners first, then add supplied-film forces.
    // Do not reconstruct a film's small positive heat by subtracting large
    // opposing mechanical and hydraulic powers.
    pub(in crate::render::plate::impact) fn supplied_film_flow_into(
        &self,
        x: &[f64],
        e: &[f64],
        out: &mut [f64],
    ) -> Option<NonlinearPortPower> {
        if !self.dissipative_flow_into(x, e, out) {
            return None;
        }
        let mut loss = 0.0_f64;
        let mut absolute = 0.0_f64;
        for (&effort, &force) in e.iter().zip(out.iter()) {
            let product = effort * force;
            loss += product;
            absolute += product.abs();
        }
        // Same arithmetic-only allowance as fs-phs's passive port helper.
        let tolerance = 8.0 * f64::EPSILON * e.len() as f64 * absolute;
        if ![loss, absolute, tolerance].iter().all(|v| v.is_finite()) || loss < -tolerance {
            return None;
        }
        let mut power = NonlinearPortPower {
            supplied_power_w: 0.0,
            dissipated_power_w: loss.max(0.0),
        };
        let mut q = [0.0; MAX_IMPACT_MODES];
        let mut velocity = [0.0; MAX_IMPACT_MODES];
        let mut force = [0.0; MAX_IMPACT_MODES];
        let mut pressure = [0.0; MAX_CELLS];
        for i in 0..self.modes {
            q[i] = x[2 * i];
            velocity[i] = e[2 * i + 1];
        }
        for film in &self.driven_sliding_films {
            let report = film
                .evaluate_force_into(
                    &q[..self.modes],
                    &velocity[..self.modes],
                    &mut force[..self.modes],
                    &mut pressure[..film.cell_count()],
                )
                .ok()?;
            for i in 0..self.modes {
                out[2 * i + 1] += force[i];
            }
            power.dissipated_power_w += report.film.dissipated_power_w;
            power.supplied_power_w += report.reservoir_power_w;
        }
        if !out.iter().all(|v| v.is_finite())
            || !power.dissipated_power_w.is_finite()
            || !power.supplied_power_w.is_finite()
        {
            return None;
        }
        Some(power)
    }

    // Differentiate the complete resisting force at the independent state and
    // discrete-effort directions, with reservoir pressures held fixed.
    pub(in crate::render::plate::impact) fn supplied_film_tangent_into(
        &self,
        x: &[f64],
        e: &[f64],
        dx: &[f64],
        de: &[f64],
        out: &mut [f64],
    ) -> bool {
        if !self.dissipative_flow_tangent_into(x, e, dx, de, out) {
            return false;
        }
        let mut q = [0.0; MAX_IMPACT_MODES];
        let mut velocity = [0.0; MAX_IMPACT_MODES];
        let mut dq = [0.0; MAX_IMPACT_MODES];
        let mut dv = [0.0; MAX_IMPACT_MODES];
        let mut force = [0.0; MAX_IMPACT_MODES];
        let reservoir_directions = [0.0; MAX_CHANNELS];
        for i in 0..self.modes {
            q[i] = x[2 * i];
            velocity[i] = e[2 * i + 1];
            dq[i] = dx[2 * i];
            dv[i] = de[2 * i + 1];
        }
        for film in &self.driven_sliding_films {
            if film
                .tangent_into(
                    &q[..self.modes],
                    &velocity[..self.modes],
                    &dq[..self.modes],
                    &dv[..self.modes],
                    &reservoir_directions[..film.channel_count()],
                    &mut force[..self.modes],
                )
                .is_err()
            {
                return false;
            }
            for i in 0..self.modes {
                out[2 * i + 1] += force[i];
            }
        }
        out.iter().all(|v| v.is_finite())
    }
}

#[cfg(test)]
mod tests;
