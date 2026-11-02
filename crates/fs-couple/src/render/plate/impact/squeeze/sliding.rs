//! Sliding film tractions in the original joint mechanical time equation.
//!
//! The native Reynolds owner computes entrainment, pressure and wall shear.
//! Its complete work-conjugate force acts through the same generalized rates
//! as contacts, supports and normal squeeze films. The existing fs-phs owner
//! advances the system once and accounts for this resistance once.
//! Wall velocity rows and layer profiles are immutable during a solve.

pub use fs_tribo::resistive_film::shear::{
    ShearChannelReport, ShearedFilm, WallVelocityPorts,
};
use fs_tribo::resistive_film::MAX_CELLS;
use super::{ImpactError, ImpactSystem, MAX_IMPACT_MODES, invalid};

impl ImpactSystem {
    /// Number of original sliding pressure islands, in insertion order.
    /// These indices are independent of stationary-wall squeeze-film indices.
    pub fn sliding_film_count(&self) -> usize { self.sliding_films.len() }

    /// Complete accepted pressure/wall operator, retaining original wall rows
    /// and, where supplied, its original through-gap viscosity profiles.
    pub fn sliding_film(&self, index: usize) -> Option<&ShearedFilm> {
        self.sliding_films.get(index)
    }

    /// Attach simultaneous wall entrainment, pressure traction and relative
    /// shear without resetting motion/history or adding a time integrator.
    ///
    /// The film rows cover every mechanical mode; zeros represent unrelated
    /// bodies. Both finite-difference and analytic prepared execution consume
    /// the same native force and tangent at the actual discrete-gradient rate.
    /// Pressure and sliding islands share the existing four-island ceiling.
    /// Compressible-gas images and mapped support memory retain their existing
    /// exclusivity; no overlapping fluid domain or unbounded recomposition is
    /// inferred. All endpoints retain the film's gap and pressure admission.
    pub fn with_sliding_film(mut self, film: ShearedFilm) -> Result<Self, ImpactError> {
        self.require_recomposition_without_mapped_memory()?;
        let pressure_operator=film.pressure_operator();
        if self.gas_film.is_some() || pressure_operator.port_count()!=self.modes
            || self.squeeze_films.len()+self.sliding_films.len()>=4 {
            return Err(invalid("sliding film mode count, duplicate gas image or four-island ceiling failed"));
        }
        let mut q=[0.0;MAX_IMPACT_MODES];
        for (i,value) in q[..self.modes].iter_mut().enumerate() {*value=self.x[2*i];}
        let mut force=[0.0;MAX_IMPACT_MODES];
        let mut pressure=[0.0;MAX_CELLS];
        film.evaluate_into(&q[..self.modes],&[0.0;MAX_IMPACT_MODES][..self.modes],
            &mut force[..self.modes],&mut pressure[..pressure_operator.cell_count()])
            .map_err(|error|invalid(error.0))?;
        self.sliding_films.push(film);
        Ok(self)
    }
}

#[cfg(test)]
mod tests;
