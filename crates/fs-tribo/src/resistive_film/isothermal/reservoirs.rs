//! Pressure-fed compressible films on the original isothermal mass owner.
//!
//! Fixed reservoir pressures contribute chemical free-energy work. They are
//! not a passive damping law: both their signed supply and nonnegative loss
//! must accompany the same discrete-effort mass flux in a coupled time owner.

use super::{GasObservation, IsothermalFilm, log_ratio};
use crate::resistive_film::{FilmError, fail};
use std::sync::Arc;

#[derive(Clone, Copy, Debug)]
pub(super) struct Reservoir {
    pub(super) pressure_pa: f64,
    /// Chemical potential divided by RT, relative to the unchanged ambient.
    pub(super) effort: f64,
}

/// One original passage's signed mass transfer and free-energy powers.
#[derive(Clone, Copy, Debug, Default)]
pub struct GasChannelFlow {
    /// Positive from the channel's `from` cell toward its neighbor or reservoir.
    pub mass_flow_kg_s: f64,
    /// Nonnegative isothermal free-energy dissipation, in watts.
    pub dissipated_power_w: f64,
    /// Positive into the film from this reservoir; zero on internal channels.
    pub reservoir_power_w: f64,
}

/// Power and mass balance at one supplied discrete-gradient effort.
///
/// `mass_storage_power + dissipated_power = reservoir_power`. The storage term
/// includes mass exchange only. Pressure work on the moving walls remains in
/// the original storage gradient and cancels its conjugate mechanical work.
/// These are rates; a joint integrator must use this SAME effort to integrate
/// work. This is an isothermal free-energy account, not an isolated-gas heat or
/// enthalpy balance, and adds no thermal destination or finite supply inventory.
#[derive(Clone, Copy, Debug, Default)]
pub struct GasFilmFlowReport {
    pub dissipated_power_w: f64,
    pub reservoir_power_w: f64,
    pub mass_storage_power_w: f64,
    /// Net gas mass entering the film from all external passages, in kg/s.
    pub reservoir_mass_rate_kg_s: f64,
    pub power_residual_w: f64,
}

/// Independently pressurized reservoirs on an existing compressible film.
///
/// The original `IsothermalFilm` still owns every cell mass, volume, ideal-gas
/// pressure, reciprocal normal force and free-energy derivative. The common
/// temperature, gas constant, opposing mechanical pressure and free-energy
/// reference remain fixed. Reservoirs share that temperature and gas species.
/// Only the exterior Poiseuille endpoints change; gas viscosity remains
/// pressure independent, no-slip and laminar with negligible channel inertia.
///
/// This is a distinct type because a pressure-fed film can deliver power.
/// It must not be inserted into a passive-only dissipative callback. There is
/// no time integrator here; mass histories stay in the caller's joint state.
/// Construction caches reservoir logarithms. Flow and complete directional
/// derivatives use the original mobility, bounded scratch and no allocations.
#[derive(Clone, Debug)]
pub struct PressureDrivenGasFilm {
    film: IsothermalFilm,
    reservoirs: Arc<[Option<Reservoir>]>,
}

impl PressureDrivenGasFilm {
    /// Bind one absolute pressure per original exterior channel, including
    /// currently closed passages. Every internal channel must contain `None`.
    /// Reservoirs obey the original film's positive absolute-pressure and gauge
    /// validity limits. No initial cell mass, scale or pressure is reset.
    pub fn new(film: IsothermalFilm, pressures_pa: &[Option<f64>]) -> Result<Self, FilmError> {
        if pressures_pa.len() != film.geometry.channel_count() {
            return fail("gas reservoir pressure count must match original channels");
        }
        let mut reservoirs = Vec::with_capacity(pressures_pa.len());
        for (channel, &pressure) in film.geometry.channels.iter().zip(pressures_pa) {
            if channel.to.is_none() != pressure.is_some() {
                return fail("gas reservoir pressure topology");
            }
            let bound = if let Some(p) = pressure {
                if !p.is_finite()
                    || p <= 0.0
                    || !p.powi(2).is_finite()
                    || (p - film.gas.pressure_pa).abs() > film.geometry.limits.maximum_pressure_pa
                {
                    return fail("gas reservoir absolute pressure or gauge domain");
                }
                let effort = log_ratio(p, film.gas.pressure_pa);
                if !effort.is_finite() {
                    return fail("gas reservoir free-energy effort is unrepresentable");
                }
                Some(Reservoir {
                    pressure_pa: p,
                    effort,
                })
            } else {
                None
            };
            reservoirs.push(bound);
        }
        Ok(Self {
            film,
            reservoirs: reservoirs.into(),
        })
    }

    /// Original channel order; absolute reservoir pressure, or no internal
    /// boundary. This immutable view allocates nothing.
    pub fn reservoir_pressures_pa(&self) -> impl ExactSizeIterator<Item = Option<f64>> + '_ {
        self.reservoirs.iter().map(|b| b.map(|b| b.pressure_pa))
    }

    /// Change supplies without resetting mass or changing the reference EOS,
    /// free energy or mechanical back pressure. Hold the new pressures fixed
    /// within a step and account their work from the returned flow report.
    pub fn with_reservoir_pressures(
        &self,
        pressures_pa: &[Option<f64>],
    ) -> Result<Self, FilmError> {
        Self::new(self.film.clone(), pressures_pa)
    }

    /// Refresh uniform transport while preserving reservoirs and gas inventory.
    pub fn with_viscosity(&self, viscosity_pa_s: f64) -> Result<Self, FilmError> {
        Ok(Self {
            film: self.film.with_viscosity(viscosity_pa_s)?,
            reservoirs: Arc::clone(&self.reservoirs),
        })
    }

    /// Refresh one viscosity per original passage; storage and supplies persist.
    pub fn with_channel_viscosities(&self, viscosities_pa_s: &[f64]) -> Result<Self, FilmError> {
        Ok(Self {
            film: self.film.with_channel_viscosities(viscosities_pa_s)?,
            reservoirs: Arc::clone(&self.reservoirs),
        })
    }

    pub fn port_count(&self) -> usize {
        self.film.port_count()
    }
    pub fn cell_count(&self) -> usize {
        self.film.cell_count()
    }
    pub fn channel_count(&self) -> usize {
        self.film.geometry.channel_count()
    }
    /// The original cold mass coordinates; binding supplies does not refill cells.
    pub fn initial_state(&self) -> &[f64] {
        self.film.initial_state()
    }
    pub fn observe(&self, q: &[f64], z: &[f64]) -> Result<GasObservation, FilmError> {
        self.film.observe(q, z)
    }
    pub fn gradient_into(
        &self,
        q: &[f64],
        z: &[f64],
        gq: &mut [f64],
        gz: &mut [f64],
    ) -> Result<(), FilmError> {
        self.film.gradient_into(q, z, gq, gz)
    }
    #[allow(clippy::too_many_arguments)]
    pub fn hessian_into(
        &self,
        q: &[f64],
        z: &[f64],
        dq: &[f64],
        dz: &[f64],
        hq: &mut [f64],
        hz: &mut [f64],
    ) -> Result<(), FilmError> {
        self.film.hessian_into(q, z, dq, dz, hq, hz)
    }

    /// Positive RESISTING mass-state flow, so the caller advances `z_dot=-out`.
    ///
    /// At the physical storage gradient each passage recovers compressible
    /// Poiseuille flow `G*(p_from²-p_to²)/(2*R*T)`. At an arbitrary joint discrete
    /// effort, the SAME logarithmic-mean mobility preserves internal mass and
    /// the exact supply/loss identity. Closed passages return exact zero.
    /// Both output slices publish only after shape, state and work admission.
    pub fn flow_into(
        &self,
        q: &[f64],
        z: &[f64],
        effort: &[f64],
        out: &mut [f64],
        channels: &mut [GasChannelFlow],
    ) -> Result<GasFilmFlowReport, FilmError> {
        self.film
            .flow_with_reservoirs::<true>(q, z, effort, &self.reservoirs, out, channels)
    }

    /// Complete derivative in geometry, retained gas masses, discrete effort,
    /// and prescribed supply pressure. `reservoir_direction_pa` has the SAME
    /// channel mapping as the supplies; use `Some(0.0)` for a held reservoir.
    /// Both logarithmic-mean mobility and reservoir chemical effort vary with
    /// supply pressure. Failure leaves `out` untouched.
    #[allow(clippy::too_many_arguments)]
    pub fn flow_tangent_into(
        &self,
        q: &[f64],
        z: &[f64],
        effort: &[f64],
        dq: &[f64],
        dz: &[f64],
        de: &[f64],
        reservoir_direction_pa: &[Option<f64>],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        if reservoir_direction_pa.len() != self.channel_count()
            || self
                .reservoirs
                .iter()
                .zip(reservoir_direction_pa)
                .any(|(b, d)| b.is_some() != d.is_some() || d.is_some_and(|d| !d.is_finite()))
        {
            return fail("gas reservoir pressure direction shape or topology");
        }
        self.film.flow_tangent_with_reservoirs::<true>(
            q,
            z,
            effort,
            dq,
            dz,
            de,
            &self.reservoirs,
            reservoir_direction_pa,
            out,
        )
    }
}

#[cfg(test)]
mod tests;
