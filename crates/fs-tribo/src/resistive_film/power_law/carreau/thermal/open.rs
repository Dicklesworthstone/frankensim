//! Open, uniformly mixed Carreau films with changing liquid inventory.
//!
//! One original single-liquid caloric chart supplies density and every sensible
//! enthalpy difference. All its density knots must agree: the Reynolds graph is
//! incompressible. Retained mass is rho times resolved volume-cell storage;
//! hydraulic passages carry flow but have no independent stored volume.
//!
//! At the accepted midpoint flow, a reservoir is the donor only when its flow
//! enters the film. Outflow carries film-midpoint enthalpy. Pressure/squeeze/
//! sliding work is already converted to native viscous heat; pressure flow work
//! must not be added again to these fixed-reference sensible enthalpy fluxes.
//!
//! The temperature root shares the original bounded kernel. With dM=M1-M0,
//! dh=h(T1)-h(T0), dh_mid=h(Tmid)-h(T0), it solves
//! M1*dh-dM*dh_mid = viscous/boundary heat
//!     + dt*sum_in(mdot*(h(T_supply)-h(Tmid))).
//! With mass residual dM+sum_out(m)-sum_in(m), the full step-relative and
//! retained-reference energy defects add that residual times h(Tmid)-h(T0)
//! and h(Tmid)-h(T_reference), respectively. Both are energy-gated; neither
//! numerical defect is turned into heat.
//! Native enthalpy differences keep this independent of the source's datum.
//!
//! This is full-film, one-fluid, uniform-temperature mixing at fixed density,
//! using the chart's fixed pressure/composition caloric convention. No local
//! temperature field, latent/cavitated heat, thermal expansion, compressibility,
//! finite supply-reservoir evolution or new material calibration is inferred.

use fs_material::phase::ResolvedHeatCapacityCurve;

use super::{
    Balance, CarreauFilm, CarreauFilmChannelReport, Exchange, FilmError, ThermalFilmBoundary,
    ThermalFilmControls, ThermalFilmReport, Volumes, advance_kernel, dot, fail, finite,
    nonnegative, positive, product,
};
use crate::resistive_film::ResistiveFilm;

#[cfg(test)]
mod tests;

/// Accepted full-film geometry, original liquid source and mixed temperature.
#[derive(Clone, Debug)]
pub struct OpenThermalCarreauFilm<'a> {
    film: CarreauFilm,
    caloric: &'a ResolvedHeatCapacityCurve,
    configuration: Vec<f64>,
    temperature_k: f64,
    time_s: f64,
    temperature_span_k: [f64; 2],
    density_kg_m3: f64,
    mass_kg: f64,
    reference_temperature_k: f64,
}

/// One physical reservoir channel, with positive transfer INTO the film.
/// Energies use the owner's retained reference temperature, not absolute h.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReservoirEnthalpy {
    pub channel: usize,
    pub mass_into_film_kg: f64,
    pub sensible_enthalpy_into_film_j: f64,
    /// None only for exactly zero transfer.
    pub donor_temperature_k: Option<f64>,
}

#[derive(Clone, Debug)]
pub struct OpenThermalFilmReport {
    /// Native midpoint mechanics, volume, heat and full sensible energy audit.
    pub thermal: ThermalFilmReport,
    pub reference_temperature_k: f64,
    pub initial_liquid_mass_kg: f64,
    pub final_liquid_mass_kg: f64,
    pub supplied_liquid_mass_kg: f64,
    pub returned_liquid_mass_kg: f64,
    /// Final minus original retained mass plus net outward transported mass.
    /// Includes represented endpoint-geometry roundoff.
    pub mass_balance_residual_kg: f64,
    /// Original channel order, exterior channels only.
    pub reservoirs: Vec<ReservoirEnthalpy>,
}

#[derive(Debug)]
pub struct OpenThermalFilmStep<'a> {
    pub next: OpenThermalCarreauFilm<'a>,
    pub report: OpenThermalFilmReport,
}

impl<'a> OpenThermalCarreauFilm<'a> {
    /// Bind the same liquid source to storage and all future supply/return
    /// enthalpy. Constant density and an entirely liquid chart are mandatory.
    /// Source receipts and the original Cp(T) integral remain with the chart.
    pub fn new(
        film: CarreauFilm,
        caloric: &'a ResolvedHeatCapacityCurve,
        configuration: &[f64],
        time_s: f64,
    ) -> Result<Self, FilmError> {
        if !nonnegative(time_s) {
            return fail("open Carreau thermal clock must be finite and nonnegative");
        }
        film.film.validate_configuration(configuration)?;
        let knots = caloric.curve().knots();
        let density = knots[0].bulk_density_kg_m3;
        if !positive(density)
            || knots
                .iter()
                .any(|k| k.liquid_mass_fraction != 1.0 || k.bulk_density_kg_m3 != density)
        {
            return fail("open Carreau film needs one original constant-density liquid chart");
        }
        let temperature_k = film.passages[0].temperature_k();
        let mut span = [knots[0].temperature_k, knots[knots.len() - 1].temperature_k];
        for passage in &film.passages {
            if passage.temperature_k() != temperature_k {
                return fail("open mixed Carreau film needs one initial material temperature");
            }
            let local = passage.law().temperature_shift().validity_k();
            span[0] = span[0].max(local[0]);
            span[1] = span[1].min(local[1]);
        }
        if span[0] >= span[1] || temperature_k < span[0] || temperature_k > span[1] {
            return fail("open Carreau caloric and rheology temperature supports do not overlap");
        }
        caloric
            .secant_heat_capacity(temperature_k, temperature_k)
            .map_err(|_| {
                FilmError("original open-film caloric source refused its initial temperature")
            })?;
        let mut volume = 0.0;
        for cell in film.film.cells.iter() {
            volume = finite(
                volume
                    + product(
                        cell.area_m2,
                        cell.gap.reference_m - dot(&cell.gap.closure, configuration),
                    )?,
            )?;
        }
        let mass_kg = product(density, volume)?;
        if !positive(mass_kg) {
            return fail("open Carreau liquid inventory must be positive and representable");
        }
        Ok(Self {
            film,
            caloric,
            configuration: configuration.to_vec(),
            temperature_k,
            time_s,
            temperature_span_k: span,
            density_kg_m3: density,
            mass_kg,
            reference_temperature_k: temperature_k,
        })
    }

    pub const fn film(&self) -> &CarreauFilm {
        &self.film
    }
    pub const fn caloric(&self) -> &ResolvedHeatCapacityCurve {
        self.caloric
    }
    pub fn configuration(&self) -> &[f64] {
        &self.configuration
    }
    pub const fn temperature_k(&self) -> f64 {
        self.temperature_k
    }
    pub const fn time_s(&self) -> f64 {
        self.time_s
    }
    pub const fn density_kg_m3(&self) -> f64 {
        self.density_kg_m3
    }
    pub const fn liquid_mass_kg(&self) -> f64 {
        self.mass_kg
    }
    pub const fn reference_temperature_k(&self) -> f64 {
        self.reference_temperature_k
    }
    pub const fn temperature_span_k(&self) -> [f64; 2] {
        self.temperature_span_k
    }

    /// Advance affine configuration, retained liquid mass and temperature.
    ///
    /// Supply temperatures follow original channel order: Some(T) for EVERY
    /// exterior channel and None internally, even for a presently outward port
    /// because its direction can reverse. They are held during this interval.
    /// Only actual inflow uses the supply value; outflow uses film temperature.
    /// The same source law describes all transported liquid.
    ///
    /// The mechanical outputs are the accepted midpoint values. Returned
    /// per-reservoir mass/enthalpy transfers allow a surrounding finite-reservoir
    /// owner to propose its own matching update. No reservoir is mutated here.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_into(
        &self,
        velocity: &[f64],
        duration_s: f64,
        supply_temperatures_k: &[Option<f64>],
        boundary: ThermalFilmBoundary,
        controls: ThermalFilmControls,
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
    ) -> Result<OpenThermalFilmStep<'a>, FilmError> {
        self.advance_with_check_into(
            velocity,
            duration_s,
            supply_temperatures_k,
            boundary,
            controls,
            force,
            pressure,
            channels,
            || Ok(()),
        )
    }

    /// Bounded source/pressure/thermal work and transactional publication, using
    /// the existing thermal controls and original material refusals.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_check_into(
        &self,
        velocity: &[f64],
        duration_s: f64,
        supply_temperatures_k: &[Option<f64>],
        boundary: ThermalFilmBoundary,
        controls: ThermalFilmControls,
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<OpenThermalFilmStep<'a>, FilmError> {
        check()?;
        if supply_temperatures_k.len() != self.film.film.channel_count() {
            return fail("open Carreau film needs source temperatures in original channel order");
        }
        for (edge, temperature) in self.film.film.channels.iter().zip(supply_temperatures_k) {
            check()?;
            match (edge.to, temperature) {
                (None, Some(t)) => {
                    self.caloric
                        .specific_enthalpy_change(self.temperature_k, *t)
                        .map_err(|_| {
                            FilmError(
                                "open Carreau supply leaves the original liquid caloric support",
                            )
                        })?;
                }
                (Some(_), None) => {}
                _ => return fail("open Carreau supply temperatures must match exterior topology"),
            }
        }
        let liquid = LiquidBalance {
            caloric: self.caloric,
            density: self.density_kg_m3,
            initial_mass: self.mass_kg,
            reference_temperature: self.reference_temperature_k,
            supplies: supply_temperatures_k,
        };
        let advance = advance_kernel::<true>(
            &self.film,
            &self.configuration,
            self.temperature_k,
            self.time_s,
            self.temperature_span_k,
            Balance::Open(liquid),
            velocity,
            duration_s,
            boundary,
            controls,
            force,
            pressure,
            channels,
            &mut check,
        )?;
        let report = OpenThermalFilmReport {
            thermal: advance.report,
            reference_temperature_k: self.reference_temperature_k,
            initial_liquid_mass_kg: advance.advection.initial_mass,
            final_liquid_mass_kg: advance.advection.final_mass,
            supplied_liquid_mass_kg: advance.advection.incoming_mass,
            returned_liquid_mass_kg: advance.advection.outgoing_mass,
            mass_balance_residual_kg: advance.advection.mass_residual,
            reservoirs: advance.advection.reservoirs,
        };
        let next = Self {
            film: advance.film,
            caloric: self.caloric,
            configuration: advance.configuration,
            temperature_k: report.thermal.final_temperature_k,
            time_s: advance.time_s,
            temperature_span_k: self.temperature_span_k,
            density_kg_m3: self.density_kg_m3,
            mass_kg: report.final_liquid_mass_kg,
            reference_temperature_k: self.reference_temperature_k,
        };
        Ok(OpenThermalFilmStep { next, report })
    }
}

#[derive(Clone, Copy)]
pub(super) struct LiquidBalance<'a> {
    caloric: &'a ResolvedHeatCapacityCurve,
    density: f64,
    initial_mass: f64,
    reference_temperature: f64,
    supplies: &'a [Option<f64>],
}

#[derive(Default)]
pub(super) struct Advection {
    pub(super) initial_mass: f64,
    pub(super) final_mass: f64,
    pub(super) incoming_mass: f64,
    pub(super) outgoing_mass: f64,
    pub(super) mass_residual: f64,
    pub(super) reservoirs: Vec<ReservoirEnthalpy>,
}

#[derive(Default)]
struct Sum {
    value: f64,
    correction: f64,
}
impl Sum {
    fn add(&mut self, value: f64) -> Result<(), FilmError> {
        let corrected = value - self.correction;
        let next = finite(self.value + corrected)?;
        self.correction = (next - self.value) - corrected;
        self.value = next;
        Ok(())
    }
}

impl LiquidBalance<'_> {
    fn dh(self, from: f64, to: f64) -> Result<f64, FilmError> {
        self.caloric
            .specific_enthalpy_change(from, to)
            .map_err(|_| FilmError("original open-film sensible enthalpy refused its temperature"))
    }

    fn masses(self, volumes: Volumes) -> Result<(f64, f64), FilmError> {
        let represented_initial = product(self.density, volumes.initial_m3)?;
        let final_mass = product(self.density, volumes.final_m3)?;
        if represented_initial != self.initial_mass || !positive(final_mass) {
            return fail("open Carreau retained mass differs from original full-film geometry");
        }
        Ok((self.initial_mass, final_mass))
    }

    pub(super) fn capacity(
        self,
        temperature: f64,
        volumes: Volumes,
        incoming_mass: f64,
    ) -> Result<f64, FilmError> {
        let (m0, m1) = self.masses(volumes)?;
        let cp = self
            .caloric
            .secant_heat_capacity(temperature, temperature)
            .map_err(|_| {
                FilmError("original open-film specific capacity refused its temperature")
            })?;
        product(finite(0.5 * m0 + 0.5 * m1 + 0.5 * incoming_mass)?, cp)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn evaluate(
        self,
        initial: f64,
        endpoint: f64,
        duration: f64,
        volumes: Volumes,
        geometry: &ResistiveFilm,
        channels: &[CarreauFilmChannelReport],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<Exchange, FilmError> {
        check()?;
        let (m0, m1) = self.masses(volumes)?;
        let change = m1 - m0;
        let midpoint = initial + 0.5 * (endpoint - initial);
        let dh = self.dh(initial, endpoint)?;
        let dh_mid = self.dh(initial, midpoint)?;
        let h0_ref = self.dh(self.reference_temperature, initial)?;
        let hmid_ref = self.dh(self.reference_temperature, midpoint)?;
        let mut mixing = Sum::default();
        let mut incoming = Sum::default();
        let mut outgoing = Sum::default();
        let mut net_out = Sum::default();
        for (index, (edge, channel)) in geometry.channels.iter().zip(channels).enumerate() {
            check()?;
            if edge.to.is_some() {
                continue;
            }
            let mass_out = product(duration, product(self.density, channel.volume_flow_m3_s)?)?;
            net_out.add(mass_out)?;
            if mass_out < 0.0 {
                let supply = self.supplies[index].ok_or(FilmError(
                    "open Carreau inflow has no original supply temperature",
                ))?;
                incoming.add(-mass_out)?;
                mixing.add(product(-mass_out, self.dh(midpoint, supply)?)?)?;
            } else {
                outgoing.add(mass_out)?;
            }
        }
        let stored = finite(product(m1, dh)? + product(change, h0_ref)?)?;
        let reduced = finite(product(m1, dh)? - product(change, dh_mid)?)?;
        let mass_residual = finite(change + net_out.value)?;
        let continuity = product(mass_residual, hmid_ref)?;
        let advected = finite(mixing.value - product(net_out.value, hmid_ref)?)?;
        Ok(Exchange {
            stored_j: stored,
            reduced_storage_j: reduced,
            mixing_j: mixing.value,
            advected_j: advected,
            continuity_defect_j: continuity,
            step_continuity_defect_j: product(mass_residual, dh_mid)?,
            initial_mass_kg: m0,
            final_mass_kg: m1,
            incoming_mass_kg: incoming.value,
            outgoing_mass_kg: outgoing.value,
            mass_residual_kg: mass_residual,
        })
    }

    pub(super) fn mixing_power(
        self,
        temperature: f64,
        geometry: &ResistiveFilm,
        channels: &[CarreauFilmChannelReport],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<f64, FilmError> {
        let mut power = Sum::default();
        for (index, (edge, channel)) in geometry.channels.iter().zip(channels).enumerate() {
            check()?;
            if edge.to.is_none() && channel.volume_flow_m3_s < 0.0 {
                let supply = self.supplies[index].ok_or(FilmError(
                    "open Carreau inflow has no original supply temperature",
                ))?;
                power.add(product(
                    product(-self.density, channel.volume_flow_m3_s)?,
                    self.dh(temperature, supply)?,
                )?)?;
            }
        }
        Ok(power.value)
    }

    #[allow(clippy::too_many_arguments)]
    pub(super) fn report(
        self,
        initial: f64,
        endpoint: f64,
        duration: f64,
        exchange: Exchange,
        geometry: &ResistiveFilm,
        channels: &[CarreauFilmChannelReport],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<Advection, FilmError> {
        let midpoint = initial + 0.5 * (endpoint - initial);
        let mut reservoirs =
            Vec::with_capacity(geometry.channels.iter().filter(|e| e.to.is_none()).count());
        for (index, (edge, channel)) in geometry.channels.iter().zip(channels).enumerate() {
            check()?;
            if edge.to.is_some() {
                continue;
            }
            let mass_into = product(-duration, product(self.density, channel.volume_flow_m3_s)?)?;
            let donor = if mass_into > 0.0 {
                self.supplies[index]
            } else if mass_into < 0.0 {
                Some(midpoint)
            } else {
                None
            };
            let energy = if let Some(t) = donor {
                product(mass_into, self.dh(self.reference_temperature, t)?)?
            } else {
                0.0
            };
            reservoirs.push(ReservoirEnthalpy {
                channel: index,
                mass_into_film_kg: mass_into,
                sensible_enthalpy_into_film_j: energy,
                donor_temperature_k: donor,
            });
        }
        Ok(Advection {
            initial_mass: exchange.initial_mass_kg,
            final_mass: exchange.final_mass_kg,
            incoming_mass: exchange.incoming_mass_kg,
            outgoing_mass: exchange.outgoing_mass_kg,
            mass_residual: exchange.mass_residual_kg,
            reservoirs,
        })
    }
}
