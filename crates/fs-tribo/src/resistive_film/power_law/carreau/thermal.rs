//! Original Carreau dissipation coupled to sensible heat storage and cooling.
//!
//! The thermal body is an explicitly uniform-temperature recirculating assembly:
//! its fixed inventory includes the film AND its returning supply fluid. All
//! pressure ports return at the same body temperature. The authored total C(T)
//! belongs to that inventory; through-flow is not silently treated as a closed
//! isolated film. Open supplies, local advected enthalpy, transverse thermal
//! gradients, thermal expansion and phase changes require another owner.
//!
//! The fixed-geometry advance admits tangential motion. The moving advance
//! also evolves the original affine gaps with held generalized velocity,
//! `q(t)=q0+t*v`. Squeeze and sliding share the original pressure graph at
//! `q_mid=q0+dt*v/2`, with source admission at both physical endpoints.
//! Film-cell volume change and signed reservoir outflow close a separate
//! continuity audit. Liquid displaced from the film stays in the declared
//! returning thermal inventory; its caloric table is never rescaled by gap.
//!
//! Held reservoir pressures, prescribed heat and bath data give
//! `E(T1)-E(T0) = dt * [D(q_mid,(T0+T1)/2) + Q - H*((T0+T1)/2-Tbath)]`.
//! Every D evaluation re-prepares the original material laws and solves their
//! complete pressure/squeeze/sliding graph. This is implicit midpoint geometry
//! and temperature quadrature, not an exact trajectory or an unconditional
//! uniqueness claim. A bounded scalar solve refuses unresolved intervals.
//! The net-power sign guard uses the SAME midpoint geometry at both thermal
//! endpoints, rejecting detected frozen-load equilibrium overshoots without
//! confusing time-varying squeeze power with a stationary thermal equilibrium.
//! It does not isolate every interior equilibrium of a nonmonotone load or
//! certify source support throughout an arbitrary nonlinear trajectory. The
//! caller owns temporal refinement and the uniform-temperature admission.

use fs_material::thermal::transport::storage::TabulatedHeatStorage;

use super::{
    CarreauFilm, CarreauFilmChannelReport, CarreauFilmReport, FilmError, LaminarCarreauSlit,
};
use crate::resistive_film::{MAX_CELLS, MAX_CHANNELS, MAX_PORTS, dot, fail};

/// Hard ceiling on original pressure solves in one proposed thermal interval.
pub const MAX_THERMAL_EVALUATIONS: usize = 256;

/// Explicit numerical budgets; none certifies a material or a lumped model.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalFilmControls {
    /// Includes bracket search, root refinement and endpoint source admission.
    pub maximum_evaluations: usize,
    /// Width of the final scalar-temperature bracket [K].
    pub temperature_tolerance_k: f64,
    pub absolute_energy_tolerance_j: f64,
    pub relative_energy_tolerance: f64,
    /// Search only this far from the accepted temperature, inside source support.
    pub maximum_temperature_change_k: f64,
}

impl Default for ThermalFilmControls {
    fn default() -> Self {
        Self {
            maximum_evaluations: 96,
            temperature_tolerance_k: 1e-8,
            absolute_energy_tolerance_j: 1e-10,
            relative_energy_tolerance: 1e-9,
            maximum_temperature_change_k: 5.0,
        }
    }
}

impl ThermalFilmControls {
    fn validate(self) -> Result<(), FilmError> {
        if !(2..=MAX_THERMAL_EVALUATIONS).contains(&self.maximum_evaluations)
            || !positive(self.temperature_tolerance_k)
            || !positive(self.maximum_temperature_change_k)
            || !nonnegative(self.absolute_energy_tolerance_j)
            || !nonnegative(self.relative_energy_tolerance)
            || self.relative_energy_tolerance >= 1.0
            || self.absolute_energy_tolerance_j == 0.0 && self.relative_energy_tolerance == 0.0
        {
            return fail("invalid Carreau thermal solve controls");
        }
        Ok(())
    }
}

/// Held external heat input and reciprocal exchange with one prescribed bath.
///
/// Positive prescribed power enters the assembly. Positive bath conductance is
/// a caller-supplied whole-assembly coefficient [W/K], not an inferred fluid or
/// surface property. All native viscous heat is already deposited once.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ThermalFilmBoundary {
    pub prescribed_power_w: f64,
    pub bath_conductance_w_k: f64,
    pub bath_temperature_k: f64,
}

impl ThermalFilmBoundary {
    fn validate(self) -> Result<(), FilmError> {
        if !self.prescribed_power_w.is_finite()
            || !nonnegative(self.bath_conductance_w_k)
            || !positive(self.bath_temperature_k)
        {
            return fail("invalid Carreau thermal power or bath data");
        }
        Ok(())
    }
}

/// Accepted original film and the caloric state of its returning inventory.
///
/// The exact declared caloric table is borrowed, never refit or extended. The
/// initial passage temperatures must agree exactly. The original affine graph
/// and fixed thermal inventory are distinct declarations: this type does not infer
/// mass, density, heat capacity, supply-loop losses or source-card provenance.
#[derive(Clone, Debug)]
pub struct ThermalCarreauFilm<'a> {
    film: CarreauFilm,
    storage: &'a TabulatedHeatStorage,
    configuration: Vec<f64>,
    temperature_k: f64,
    time_s: f64,
    temperature_span_k: [f64; 2],
}

/// Accepted midpoint mechanics and separate physical/numerical energy accounts.
#[derive(Clone, Copy, Debug)]
pub struct ThermalFilmReport {
    pub initial_temperature_k: f64,
    pub final_temperature_k: f64,
    pub midpoint_temperature_k: f64,
    pub duration_s: f64,
    /// Resolved volume cells, excluding the zero-storage hydraulic passages.
    pub initial_film_cell_volume_m3: f64,
    pub final_film_cell_volume_m3: f64,
    /// Affine volume increment from the original cell-area/closure rows.
    /// Evaluated from rates, avoiding subtraction of nearly equal volumes.
    pub film_cell_volume_change_m3: f64,
    /// Signed midpoint outflow from the film into its returning inventory.
    pub reservoir_outflow_volume_m3: f64,
    /// Cell volume change plus reservoir outflow; a numerical mass-balance
    /// defect, not an extra fluid inventory or heat source.
    pub volume_balance_residual_m3: f64,
    pub stored_energy_change_j: f64,
    pub dissipated_heat_j: f64,
    /// Signed imposed heat into the fixed thermal inventory [J].
    pub prescribed_heat_j: f64,
    /// Signed bath heat into the fixed thermal inventory [J].
    pub bath_heat_j: f64,
    pub mechanical_work_j: f64,
    pub reservoir_work_j: f64,
    /// Storage change minus physical viscous and boundary heat [J].
    pub thermal_residual_j: f64,
    /// Native midpoint mechanical/supply work minus viscous heat [J].
    /// It is a numerical residual and is never deposited as heat.
    pub mechanical_work_residual_j: f64,
    /// Storage change minus imposed mechanical/supply and boundary work [J].
    pub combined_energy_residual_j: f64,
    pub thermal_residual_allowance_j: f64,
    pub temperature_bracket_width_k: f64,
    /// Original complete mechanical report at midpoint temperature and geometry.
    pub midpoint_mechanics: CarreauFilmReport,
    /// Includes the final original-law endpoint mechanical support check.
    pub evaluations: usize,
}

/// Complete immutable proposal. All mutable output slices publish with it.
#[derive(Debug)]
pub struct ThermalFilmStep<'a> {
    pub next: ThermalCarreauFilm<'a>,
    pub report: ThermalFilmReport,
}

impl<'a> ThermalCarreauFilm<'a> {
    pub fn new(
        film: CarreauFilm,
        storage: &'a TabulatedHeatStorage,
        configuration: &[f64],
        time_s: f64,
    ) -> Result<Self, FilmError> {
        if !nonnegative(time_s) {
            return fail("Carreau thermal clock must be finite and nonnegative");
        }
        film.film.validate_configuration(configuration)?;
        let temperature_k = film.passages[0].temperature_k();
        if film
            .passages
            .iter()
            .any(|p| p.temperature_k() != temperature_k)
        {
            return fail("uniform Carreau thermal inventory needs identical initial temperatures");
        }
        let mut span = storage.temperature_span_k();
        for passage in &film.passages {
            let local = passage.law().temperature_shift().validity_k();
            span[0] = span[0].max(local[0]);
            span[1] = span[1].min(local[1]);
        }
        if span[0] >= span[1] || temperature_k < span[0] || temperature_k > span[1] {
            return fail("Carreau material and caloric temperature supports do not overlap");
        }
        storage.heat_capacity_j_k(temperature_k).map_err(|_| {
            FilmError("original Carreau assembly heat storage refused its initial state")
        })?;
        Ok(Self {
            film,
            storage,
            configuration: configuration.to_vec(),
            temperature_k,
            time_s,
            temperature_span_k: span,
        })
    }

    pub const fn film(&self) -> &CarreauFilm {
        &self.film
    }
    pub const fn storage(&self) -> &TabulatedHeatStorage {
        self.storage
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
    pub const fn temperature_span_k(&self) -> [f64; 2] {
        self.temperature_span_k
    }

    /// Advance with held tangential velocities and boundary data. Returned
    /// forces, pressures and channel reports are the actual MIDPOINT mechanics;
    /// the returned owner's film is prepared at the accepted endpoint instead.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_into(
        &self,
        velocity: &[f64],
        duration_s: f64,
        boundary: ThermalFilmBoundary,
        controls: ThermalFilmControls,
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
    ) -> Result<ThermalFilmStep<'a>, FilmError> {
        self.advance_with_check_into(
            velocity,
            duration_s,
            boundary,
            controls,
            force,
            pressure,
            channels,
            || Ok(()),
        )
    }

    /// Advance configuration and temperature together under held generalized
    /// velocity: `q_next=q+duration*velocity`. Original cell/channel affine gaps
    /// are admitted at the midpoint and both endpoints. Native midpoint normal
    /// pressure work and tangential work feed the SAME viscous heat calculation.
    ///
    /// The caloric inventory remains the supplied whole recirculating body.
    /// Fluid squeezed out of the resolved film stays in that body at its uniform
    /// temperature. No open-supply enthalpy, local heat partition, density change,
    /// inertia, cavitation or new thermal storage is inferred.
    ///
    /// Returned forces/pressures are midpoint values. The next owner retains
    /// endpoint configuration and source temperature, ready for continuation.
    /// Time/geometry quadrature needs refinement for the actual gap evolution.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_moving_into(
        &self,
        velocity: &[f64],
        duration_s: f64,
        boundary: ThermalFilmBoundary,
        controls: ThermalFilmControls,
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
    ) -> Result<ThermalFilmStep<'a>, FilmError> {
        self.advance_moving_with_check_into(
            velocity,
            duration_s,
            boundary,
            controls,
            force,
            pressure,
            channels,
            || Ok(()),
        )
    }

    /// Checked moving-gap advance. Initial and final physical source admission
    /// count against the SAME pressure-solve budget as the thermal root. Refusal
    /// preserves the accepted configuration, temperature, clock and all outputs.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_moving_with_check_into(
        &self,
        velocity: &[f64],
        duration_s: f64,
        boundary: ThermalFilmBoundary,
        controls: ThermalFilmControls,
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
        check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<ThermalFilmStep<'a>, FilmError> {
        self.advance_kernel::<true>(
            velocity, duration_s, boundary, controls, force, pressure, channels, check,
        )
    }

    /// Refusal/cancellation leaves self and every caller output unchanged.
    /// One candidate film is allocated per interval; thermal root iterations
    /// re-use its prepared passage/coefficient arrays and bounded stack scratch.
    #[allow(clippy::too_many_arguments)]
    pub fn advance_with_check_into(
        &self,
        velocity: &[f64],
        duration_s: f64,
        boundary: ThermalFilmBoundary,
        controls: ThermalFilmControls,
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
        check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<ThermalFilmStep<'a>, FilmError> {
        self.advance_kernel::<false>(
            velocity, duration_s, boundary, controls, force, pressure, channels, check,
        )
    }

    #[allow(clippy::too_many_arguments)]
    fn advance_kernel<const MOVING: bool>(
        &self,
        velocity: &[f64],
        duration_s: f64,
        boundary: ThermalFilmBoundary,
        controls: ThermalFilmControls,
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [CarreauFilmChannelReport],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<ThermalFilmStep<'a>, FilmError> {
        check()?;
        controls.validate()?;
        boundary.validate()?;
        let final_time = self.time_s + duration_s;
        let geometry = &self.film.film;
        if !positive(duration_s)
            || !final_time.is_finite()
            || final_time <= self.time_s
            || velocity.len() != geometry.port_count()
            || velocity.iter().any(|v| !v.is_finite())
            || force.len() != geometry.port_count()
            || pressure.len() != geometry.cell_count()
            || channels.len() != geometry.channel_count()
        {
            return fail("invalid Carreau thermal interval, velocity or output shape");
        }
        let moving_gaps = geometry
            .cells
            .iter()
            .any(|cell| dot(&cell.gap.closure, velocity) != 0.0)
            || geometry
                .channels
                .iter()
                .any(|edge| dot(&edge.gap.closure, velocity) != 0.0);
        if !MOVING && moving_gaps {
            return fail("Carreau thermal fixed geometry refuses gap-changing velocity");
        }
        let mut midpoint_q = [0.0; MAX_PORTS];
        let mut final_q = [0.0; MAX_PORTS];
        for (index, (&position, &rate)) in self.configuration.iter().zip(velocity).enumerate() {
            midpoint_q[index] = if MOVING {
                finite((0.5 * duration_s).mul_add(rate, position))?
            } else {
                position
            };
            final_q[index] = if MOVING {
                finite(duration_s.mul_add(rate, position))?
            } else {
                position
            };
            if MOVING && rate != 0.0 && final_q[index] == position {
                return fail("Carreau moving thermal configuration increment is unrepresentable");
            }
        }
        let midpoint_q = &midpoint_q[..geometry.port_count()];
        let final_q = &final_q[..geometry.port_count()];
        // For affine gaps, endpoint admission retains the entire geometric
        // domain; admit the rounded midpoint independently as well.
        geometry.validate_configuration(midpoint_q)?;
        geometry.validate_configuration(final_q)?;
        let final_checks = if moving_gaps { 2 } else { 1 };
        if controls.maximum_evaluations <= final_checks + usize::from(moving_gaps) {
            return fail("Carreau thermal pressure-solve budget exhausted");
        }
        let mut candidate = self.film.clone();
        let mut scratch = Scratch::new();
        let mut evaluations = 0;
        if moving_gaps {
            // A mid-gap trial cannot conceal an unsupported starting pressure,
            // rate or wall state. This scratch is replaced by accepted midpoint
            // mechanics before publication.
            scratch.evaluate(&self.film, &self.configuration, velocity, &mut check)?;
            evaluations += 1;
        }
        let mut solve = |endpoint: f64, candidate: &mut CarreauFilm, scratch: &mut Scratch| {
            if evaluations + final_checks >= controls.maximum_evaluations {
                return fail("Carreau thermal pressure-solve budget exhausted");
            }
            evaluations += 1;
            check()?;
            let midpoint = self.temperature_k + 0.5 * (endpoint - self.temperature_k);
            prepare_temperature(candidate, midpoint, &mut check)?;
            let mechanics = scratch.evaluate(candidate, midpoint_q, velocity, &mut check)?;
            let stored = self
                .storage
                .energy_change_j(self.temperature_k, endpoint)
                .map_err(|_| {
                    FilmError("original Carreau assembly heat storage refused its endpoint")
                })?;
            let heat = product(duration_s, mechanics.film.dissipated_power_w)?;
            let imposed = product(duration_s, boundary.prescribed_power_w)?;
            let bath = product(duration_s, bath_power(boundary, midpoint)?)?;
            let residual = finite((stored - heat) - (imposed + bath))?;
            let scale = stored.abs().max(heat).max(imposed.abs()).max(bath.abs());
            let allowance = finite(
                controls.absolute_energy_tolerance_j + controls.relative_energy_tolerance * scale,
            )?;
            Ok(Trial {
                temperature: endpoint,
                stored,
                heat,
                imposed,
                bath,
                residual,
                allowance,
                mechanics,
            })
        };

        let initial = solve(self.temperature_k, &mut candidate, &mut scratch)?;
        let (accepted, bracket_width) = if initial.residual == 0.0 {
            (initial, 0.0)
        } else {
            let direction = -initial.residual.signum();
            let bound = if direction > 0.0 {
                (self.temperature_k + controls.maximum_temperature_change_k)
                    .min(self.temperature_span_k[1])
            } else {
                (self.temperature_k - controls.maximum_temperature_change_k)
                    .max(self.temperature_span_k[0])
            };
            let limit = (bound - self.temperature_k).abs();
            if !positive(limit) {
                return fail("Carreau thermal step has no supported local temperature interval");
            }
            let capacity = self
                .storage
                .heat_capacity_j_k(self.temperature_k)
                .map_err(|_| {
                    FilmError("original Carreau heat capacity refused its accepted temperature")
                })?;
            let scale =
                finite(capacity + product(0.5 * duration_s, boundary.bath_conductance_w_k)?)?;
            let mut distance = (initial.residual.abs() / scale).min(limit);
            if !positive(distance) {
                return fail("Carreau thermal bracket increment is unrepresentable");
            }
            let outer = loop {
                let temperature = if distance == limit {
                    bound
                } else {
                    self.temperature_k + direction * distance
                };
                if temperature == self.temperature_k {
                    return fail("Carreau thermal temperature change is unrepresentable");
                }
                let trial = solve(temperature, &mut candidate, &mut scratch)?;
                if trial.residual == 0.0 || trial.residual.signum() != initial.residual.signum() {
                    break trial;
                }
                if distance == limit {
                    return fail(
                        "Carreau thermal step has no bracket inside its local source limits",
                    );
                }
                distance = (2.0 * distance).min(limit);
            };
            if outer.residual == 0.0 {
                (outer, 0.0)
            } else {
                let (mut low, mut high) = if direction > 0.0 {
                    (initial, outer)
                } else {
                    (outer, initial)
                };
                loop {
                    let width = high.temperature - low.temperature;
                    // Ratio arithmetic avoids overflowing the sum of residual
                    // magnitudes. The interior safeguard guarantees progress.
                    let a = low.residual.abs();
                    let b = high.residual.abs();
                    let fraction = if a >= b {
                        1.0 / (1.0 + b / a)
                    } else {
                        let ratio = a / b;
                        ratio / (1.0 + ratio)
                    };
                    let fraction = if (0.1..=0.9).contains(&fraction) {
                        fraction
                    } else {
                        0.5
                    };
                    let temperature = low.temperature + fraction * width;
                    if temperature <= low.temperature || temperature >= high.temperature {
                        return fail("Carreau thermal bracket exhausted floating-point resolution");
                    }
                    let trial = solve(temperature, &mut candidate, &mut scratch)?;
                    if trial.residual == 0.0
                        || width <= controls.temperature_tolerance_k
                            && trial.residual.abs() <= trial.allowance
                    {
                        break (trial, if trial.residual == 0.0 { 0.0 } else { width });
                    }
                    if trial.residual.signum() == low.residual.signum() {
                        low = trial;
                    } else {
                        high = trial;
                    }
                }
            }
        };
        // Release the trial closure's mutable counters/check before the final
        // endpoint solve. Its outputs must not overwrite accepted midpoint work.
        drop(solve);
        if accepted.residual.abs() > accepted.allowance {
            return fail("Carreau thermal candidate missed its caloric energy gate");
        }
        prepare_temperature(&mut candidate, accepted.temperature, &mut check)?;
        let mut endpoint_scratch = Scratch::new();
        let endpoint = endpoint_scratch.evaluate(&candidate, midpoint_q, velocity, &mut check)?;
        evaluations += 1;
        if moving_gaps {
            // Admit the endpoint geometry with its endpoint material. Keep the
            // midpoint-geometry thermal endpoint above for the frozen-load
            // equilibrium guard; nonautonomous squeeze power may itself turn.
            endpoint_scratch.evaluate(&candidate, final_q, velocity, &mut check)?;
            evaluations += 1;
        }
        let endpoint_net = finite(
            endpoint.film.dissipated_power_w
                + boundary.prescribed_power_w
                + bath_power(boundary, accepted.temperature)?,
        )?;
        let endpoint_heat = product(duration_s, endpoint_net)?;
        // A reversed thermal-endpoint power at the SAME midpoint geometry
        // detects a crossed stationary temperature of this frozen-load solve,
        // including stiff cooling. Physical geometry changes may move that
        // equilibrium; their endpoint source check is separate above.
        if initial.residual < 0.0 && endpoint_heat < -accepted.allowance
            || initial.residual > 0.0 && endpoint_heat > accepted.allowance
        {
            return fail("Carreau thermal step crosses a stationary temperature; refine duration");
        }
        let mechanical_work_j = product(duration_s, accepted.mechanics.mechanical_power_w)?;
        let reservoir_work_j = product(duration_s, accepted.mechanics.reservoir_power_w)?;
        let mechanical_work_residual_j = product(duration_s, accepted.mechanics.work_defect_w)?;
        let mut initial_volume = 0.0;
        let mut final_volume = 0.0;
        let mut volume_rate = 0.0;
        for (index, cell) in geometry.cells.iter().enumerate() {
            if index.is_multiple_of(32) {
                check()?;
            }
            initial_volume = finite(
                initial_volume
                    + product(
                        cell.area_m2,
                        cell.gap.reference_m - dot(&cell.gap.closure, &self.configuration),
                    )?,
            )?;
            final_volume = finite(
                final_volume
                    + product(
                        cell.area_m2,
                        cell.gap.reference_m - dot(&cell.gap.closure, final_q),
                    )?,
            )?;
            volume_rate =
                finite(volume_rate - product(cell.area_m2, dot(&cell.gap.closure, velocity))?)?;
        }
        let volume_change = if MOVING {
            product(duration_s, volume_rate)?
        } else {
            0.0
        };
        let mut outflow = 0.0;
        let mut compensation = 0.0;
        for (index, (edge, channel)) in geometry.channels.iter().zip(&scratch.channels).enumerate()
        {
            if index.is_multiple_of(32) {
                check()?;
            }
            if edge.to.is_none() {
                let corrected = channel.volume_flow_m3_s - compensation;
                let next = finite(outflow + corrected)?;
                compensation = (next - outflow) - corrected;
                outflow = next;
            }
        }
        let displaced_volume = product(duration_s, outflow)?;
        let report = ThermalFilmReport {
            initial_temperature_k: self.temperature_k,
            final_temperature_k: accepted.temperature,
            midpoint_temperature_k: self.temperature_k
                + 0.5 * (accepted.temperature - self.temperature_k),
            duration_s,
            initial_film_cell_volume_m3: initial_volume,
            final_film_cell_volume_m3: final_volume,
            film_cell_volume_change_m3: volume_change,
            reservoir_outflow_volume_m3: displaced_volume,
            volume_balance_residual_m3: finite(volume_change + displaced_volume)?,
            stored_energy_change_j: accepted.stored,
            dissipated_heat_j: accepted.heat,
            prescribed_heat_j: accepted.imposed,
            bath_heat_j: accepted.bath,
            mechanical_work_j,
            reservoir_work_j,
            thermal_residual_j: accepted.residual,
            mechanical_work_residual_j,
            combined_energy_residual_j: finite(accepted.residual - mechanical_work_residual_j)?,
            thermal_residual_allowance_j: accepted.allowance,
            temperature_bracket_width_k: bracket_width,
            midpoint_mechanics: accepted.mechanics,
            evaluations,
        };
        let next = Self {
            film: candidate,
            storage: self.storage,
            configuration: final_q.to_vec(),
            temperature_k: accepted.temperature,
            time_s: final_time,
            temperature_span_k: self.temperature_span_k,
        };
        check()?;
        force.copy_from_slice(&scratch.force[..geometry.port_count()]);
        pressure.copy_from_slice(&scratch.pressure[..geometry.cell_count()]);
        channels.copy_from_slice(&scratch.channels[..geometry.channel_count()]);
        Ok(ThermalFilmStep { next, report })
    }
}

#[derive(Clone, Copy)]
struct Trial {
    temperature: f64,
    stored: f64,
    heat: f64,
    imposed: f64,
    bath: f64,
    residual: f64,
    allowance: f64,
    mechanics: CarreauFilmReport,
}

struct Scratch {
    force: [f64; MAX_PORTS],
    pressure: [f64; MAX_CELLS],
    channels: [CarreauFilmChannelReport; MAX_CHANNELS],
}

impl Scratch {
    fn new() -> Self {
        Self {
            force: [0.0; MAX_PORTS],
            pressure: [0.0; MAX_CELLS],
            channels: [CarreauFilmChannelReport::default(); MAX_CHANNELS],
        }
    }

    fn evaluate(
        &mut self,
        film: &CarreauFilm,
        q: &[f64],
        velocity: &[f64],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<CarreauFilmReport, FilmError> {
        film.evaluate_with_check_into(
            q,
            velocity,
            &mut self.force[..film.film.port_count()],
            &mut self.pressure[..film.film.cell_count()],
            &mut self.channels[..film.film.channel_count()],
            check,
        )
    }
}

/// The candidate owns these arrays; accepted preparations are never touched.
fn prepare_temperature(
    film: &mut CarreauFilm,
    temperature: f64,
    check: &mut impl FnMut() -> Result<(), FilmError>,
) -> Result<(), FilmError> {
    for (index, passage) in film.passages.iter_mut().enumerate() {
        check()?;
        let next = LaminarCarreauSlit::new(
            *passage.law(),
            temperature,
            passage.width_m(),
            passage.gap_m(),
            passage.length_m(),
            passage.maximum_laminar_shear_rate_per_s(),
        )
        .map_err(|_| FilmError("original Carreau source refused the thermal trial"))?;
        let viscosity = next.zero_shear_viscosity_pa_s();
        // Same native Newtonian geometry seed as with_channel_viscosities;
        // thermal trials retain already allocated, originally ordered arrays.
        let edge = &film.film.channels[index];
        let factor = (edge.width_m / edge.length_m) / (12.0 * viscosity);
        if !positive(factor) {
            return fail("Carreau thermal conductance seed is unrepresentable");
        }
        film.film.factors[index] = factor;
        film.film.channel_viscosities_pa_s[index] = viscosity;
        *passage = next;
    }
    Ok(())
}

fn bath_power(boundary: ThermalFilmBoundary, temperature: f64) -> Result<f64, FilmError> {
    product(
        boundary.bath_conductance_w_k,
        boundary.bath_temperature_k - temperature,
    )
}
fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}
fn nonnegative(value: f64) -> bool {
    value.is_finite() && value >= 0.0
}
fn finite(value: f64) -> Result<f64, FilmError> {
    if value.is_finite() {
        Ok(value)
    } else {
        fail("unrepresentable Carreau thermal arithmetic")
    }
}
fn product(a: f64, b: f64) -> Result<f64, FilmError> {
    let value = finite(a * b)?;
    if a != 0.0 && b != 0.0 && value == 0.0 {
        fail("nonzero Carreau thermal work underflowed")
    } else {
        Ok(value)
    }
}

#[cfg(test)]
mod tests;
