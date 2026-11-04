//! Prescribed pressure supplies with simultaneous sliding-wall entrainment.
//!
//! The original Reynolds graph owns pressure, normal force, Couette flow and
//! wall traction. Supply power uses TOTAL boundary flow, while irreversible
//! heat is pressure-driven flow times pressure drop plus relative-wall drag:
//! `F dot v - sum_boundary(p_b * Q_total,out) = viscous_heat`.
//! Supplies can drive the walls, and moving walls can pump against a supply.
//!
//! Pressures retain the original mechanical and viscosity gauge datum. Barus
//! pressure dependence and retained layer profiles use the same pressure solve;
//! this adds no reservoir inventory, thermal time owner, inertia or cavitation.

use super::{
    FilmError, PressureDrivenReport, ReservoirPressures, ShearChannelReport, ShearedFilm,
};

/// Pressure-fed version of one original sliding-wall Reynolds film.
///
/// Every original boundary channel requires `Some(gauge_pressure_pa)`, and
/// every internal channel requires `None`. Supply potentials are compiled once.
/// Evaluation and tangents reuse the native factorization and bounded stack
/// storage, with no heap allocation. Refused evaluations leave outputs intact.
#[derive(Clone, Debug)]
pub struct PressureDrivenShearedFilm {
    film: ShearedFilm,
    reservoirs: ReservoirPressures,
}

impl PressureDrivenShearedFilm {
    /// Bind prescribed pressures in original channel order.
    pub fn new(
        film: ShearedFilm,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        let reservoirs = ReservoirPressures::new(&film.film, pressures_pa)?;
        Ok(Self { film, reservoirs })
    }

    /// Original geometry, wall rows and retained viscosity profiles.
    #[must_use]
    pub const fn sheared_film(&self) -> &ShearedFilm {
        &self.film
    }

    /// Prescribed gauge pressures in original channel order.
    #[must_use]
    pub fn reservoir_pressures_pa(&self) -> &[Option<f64>] {
        &self.reservoirs.reservoirs_pa
    }

    /// Rebind supply pressures without changing geometry, profiles or walls.
    pub fn with_reservoir_pressures(
        &self,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        Self::new(self.film.clone(), pressures_pa)
    }

    /// Refresh layer viscosities while retaining the prescribed supplies.
    ///
    /// Values follow the original sliding film's flat layer order. This is a
    /// constitutive rebind; it introduces no thermal state or time integration.
    pub fn with_layer_viscosities(
        &self,
        viscosities_pa_s: &[f64],
    ) -> Result<Self, FilmError> {
        Self::new(
            self.film.with_layer_viscosities(viscosities_pa_s)?,
            &self.reservoirs.reservoirs_pa,
        )
    }

    /// Evaluate simultaneous force, cell pressures, total flows and work.
    ///
    /// Mechanical power and reservoir power may each have either sign. Their
    /// sum equals the nonnegative irreversible heat in the returned film report.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearChannelReport],
    ) -> Result<PressureDrivenReport, FilmError> {
        self.film.evaluate_with_reservoirs::<true, false>(
            q,
            velocity,
            force,
            pressure,
            channels,
            &mut [],
            Some(&self.reservoirs),
        )
    }

    /// Also distribute viscous work through the original retained gap layers.
    ///
    /// Local pressure/shear cross terms are retained. Layer powers sum to each
    /// channel's irreversible heat, rather than to its total pressure transport.
    #[allow(clippy::too_many_arguments)]
    pub fn evaluate_with_layers_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [ShearChannelReport],
        layer_power_w: &mut [f64],
    ) -> Result<PressureDrivenReport, FilmError> {
        self.film.evaluate_with_reservoirs::<true, true>(
            q,
            velocity,
            force,
            pressure,
            channels,
            layer_power_w,
            Some(&self.reservoirs),
        )
    }

    /// Exact force derivative in configuration, rate and supply-pressure directions.
    ///
    /// Supply directions use original channel order and must be zero on internal
    /// channels. Geometry, pressure-dependent viscosity, wall traction and the
    /// boundary contribution are differentiated together in the native solve.
    #[allow(clippy::too_many_arguments)]
    pub fn tangent_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        reservoir_directions_pa: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        self.film.tangent_with_reservoirs(
            q,
            velocity,
            dq,
            dv,
            out,
            Some(&self.reservoirs),
            reservoir_directions_pa,
        )
    }
}
