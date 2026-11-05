//! Pressure supplies for the original regularized power-law Reynolds film.
//!
//! Each passage retains its exact no-slip constitutive flux. The same bounded
//! Newton solve enforces cell continuity with the authored boundary pressures.
//! Mechanical input plus `-sum_boundary(p_b Q_out)` equals viscous heat
//! `sum_channels(Q delta_p)`. Either input may be negative; heat remains
//! nonnegative. Pressure offsets affect normal loading but do not change this
//! pressure-independent rheology. Relative pressure coordinates keep a large
//! common supply bias out of small hydraulic drops.

use super::{
    FilmChannelReport, FilmError, PowerLawFilm, PowerLawFilmReport, RegularizedPowerLaw,
    ReservoirPressures, fail,
};

/// Nonlinear film diagnostics with the explicit mechanical and supply ports.
/// These are instantaneous powers, not integrated energy or stored heat.
#[derive(Debug, Clone, Copy)]
pub struct PressureDrivenPowerLawReport {
    /// Original pressure, flow, work and bounded-solve diagnostics. Its
    /// `work_defect_w` includes both mechanical and reservoir input.
    pub power_law: PowerLawFilmReport,
    /// `F dot v`, positive into the film; may be negative under a supply load.
    pub mechanical_power_w: f64,
    /// `-sum_boundary(p_b Q_out)`, positive into the film from supplies.
    pub reservoir_power_w: f64,
}

/// Pressure-fed generalized-Newtonian film with explicit active power ports.
///
/// Every original boundary passage requires `Some(gauge_pressure_pa)`, including
/// a shut passage; every internal passage requires `None`. All pressures use
/// the original mechanical gauge datum and validity limits. This separate type
/// cannot silently substitute for the passive `PowerLawFilm` resistance.
///
/// Geometry, laws and supplies are frozen within a solve. Evaluation and the
/// full tangent use bounded stack storage without per-query heap allocation.
/// Fluid inertia, sliding walls, cavitation, supply storage and a thermal time
/// owner are not supplied by this instantaneous Reynolds reduction.
#[derive(Debug, Clone)]
pub struct PressureDrivenPowerLawFilm {
    film: PowerLawFilm,
    reservoirs: ReservoirPressures,
    ambient: bool,
}

impl PressureDrivenPowerLawFilm {
    /// Bind the supplied pressures to the original channel order and domain.
    pub fn new(film: PowerLawFilm, pressures_pa: &[Option<f64>]) -> Result<Self, FilmError> {
        let reservoirs = ReservoirPressures::new(&film.film, pressures_pa)?;
        let ambient = pressures_pa.iter().flatten().all(|p| *p == 0.0);
        Ok(Self {
            film,
            reservoirs,
            ambient,
        })
    }

    /// Original nonlinear constitutive owner. Its own evaluations omit the
    /// supplies; use this binding's methods for the complete supplied film.
    #[must_use]
    pub const fn power_law_film(&self) -> &PowerLawFilm {
        &self.film
    }

    /// Prescribed gauge pressures in the original channel order.
    #[must_use]
    pub fn reservoir_pressures_pa(&self) -> &[Option<f64>] {
        &self.reservoirs.reservoirs_pa
    }

    /// Rebind supplies without changing the geometry or constitutive laws.
    /// Refusals leave this original immutable binding intact.
    pub fn with_reservoir_pressures(
        &self,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        Self::new(self.film.clone(), pressures_pa)
    }

    /// Refresh the channel laws while retaining geometry, controls and supplies.
    /// This introduces no thermal state, material authority or elapsed time.
    pub fn with_laws(&self, laws: Vec<RegularizedPowerLaw>) -> Result<Self, FilmError> {
        Self::new(
            PowerLawFilm::new(&self.film.film, laws, self.film.controls)?,
            self.reservoir_pressures_pa(),
        )
    }

    fn evaluation_reservoirs(&self) -> Option<&ReservoirPressures> {
        if self.ambient {
            None
        } else {
            Some(&self.reservoirs)
        }
    }

    /// Solve force and physical pressure with complete mechanical/supply work.
    pub fn evaluate_force_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        self.film.evaluate_with_reservoirs::<false>(
            q,
            velocity,
            force,
            pressure,
            &mut [],
            self.evaluation_reservoirs(),
            &mut || Ok(()),
        )
    }

    /// Publish force, cell pressure and original-order passage flow/heat together.
    /// A surrounding integrator must account for the signed reservoir power.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        self.evaluate_with_check_into(q, velocity, force, pressure, channels, || Ok(()))
    }

    /// Same bounded nonlinear solve with original caller-owned checks between
    /// channel batches, Newton/line-search stages and final publication. Any
    /// refusal, including the callback, leaves every output slice untouched.
    pub fn evaluate_with_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        self.film.evaluate_with_reservoirs::<true>(
            q,
            velocity,
            force,
            pressure,
            channels,
            self.evaluation_reservoirs(),
            &mut check,
        )
    }

    /// Differentiate geometry, normal velocity and reservoir pressure together.
    ///
    /// Supply directions follow original channel order and are zero on internal
    /// passages. The accepted nonlinear differential conductances give
    /// `J delta_p = A B delta_v - div(Q_h delta_h) + sum_boundary(Q_p delta_p_b)`.
    /// Its same factorization is reused; no pressure finite differences are used.
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
        self.tangent_with_check_into(q, velocity, dq, dv, reservoir_directions_pa, out, || Ok(()))
    }

    /// Same analytic tangent and atomic output with caller-owned bounded checks.
    #[allow(clippy::too_many_arguments)]
    pub fn tangent_with_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        reservoir_directions_pa: &[f64],
        out: &mut [f64],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(), FilmError> {
        if reservoir_directions_pa.len() != self.film.film.channel_count() {
            return fail("invalid pressure-driven power-law supply direction count");
        }
        if self.ambient && reservoir_directions_pa.iter().all(|d| *d == 0.0) {
            return self
                .film
                .tangent_with_check_into(q, velocity, dq, dv, out, check);
        }
        self.film.tangent_with_reservoirs(
            q,
            velocity,
            dq,
            dv,
            out,
            Some(&self.reservoirs),
            reservoir_directions_pa,
            &mut check,
        )
    }
}
