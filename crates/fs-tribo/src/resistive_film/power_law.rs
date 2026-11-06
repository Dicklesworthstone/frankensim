//! Regularized power-law pressure flow on the native Reynolds film graph.
//!
//! Each passage has the SI generalized-Newtonian law
//! mu(gamma) = K * max(gamma, gamma0)^(n-1), with n>0 and gamma0>0.
//! gamma0 is a declared physical Newtonian-core transition, not a numerical
//! strain clamp. The stress law is continuous and strictly increasing. There
//! is no yield stress, elastic memory or inferred material/source authority.
//!
//! For half-gap a=h/2 and tau_w=|delta_p|*a/length, integration of
//! Q=2*width*integral_0^a(y*gamma(|dp/dx|*y) dy) gives the exact no-slip
//! channel flux. Below tau0=K*gamma0^n it is ordinary Poiseuille flow.
//! Above tau0, with m=2+1/n and R=(tau0/tau_w)^m,
//! |Q|=2*width*a^2*gamma_w*((1-R)/m+R/3).
//! The positive blends used below give its pressure and aperture derivatives
//! without subtracting nearly equal terms at the core transition.
//!
//! A bounded damped Newton solve enforces the original cell-volume balance.
//! It reuses the native equilibrated pressure factorization; force remains
//! B^T*A*p and irreversible power remains sum(Q*delta_p). The directional
//! tangent uses the SAME accepted differential conductances. The n=1 image
//! retains the original Newtonian pressure solve and tangent arithmetic.
//!
//! Geometry/laws are fixed during a call. This is an incompressible,
//! stationary-wall, full-film thin-gap reduction. The original owner uses
//! zero-gauge reservoirs; [`pressure::PressureDrivenPowerLawFilm`] adds explicit
//! pressure supplies and their work without changing the constitutive solver.
//! [`shear::ShearedPowerLawFilm`] supplies the combined nonlinear sliding and
//! pressure profile with work-conjugate wall tractions. Pressure-dependent
//! viscosity, cavitation, fluid inertia, layer/thermal evolution and shear-history
//! models are not inferred.

use fs_math::det;

mod nonlinear;
use nonlinear::ReynoldsLaw;
/// Original thermal Carreau source on the shared Reynolds pressure solver.
pub mod carreau;

/// Prescribed pressure supplies on the same non-Newtonian pressure operator.
pub mod pressure;
pub use pressure::{PressureDrivenPowerLawFilm, PressureDrivenPowerLawReport};
/// Moving walls on the same pressure-fed nonlinear film.
pub mod shear;
pub use shear::{ShearedPowerLawChannelReport, ShearedPowerLawFilm};

use super::pressure::ReservoirPressures;
use super::{
    FilmChannelReport, FilmError, FilmReport, FlowBalance, MAX_CELLS, MAX_CHANNELS, MAX_PORTS,
    ResistiveFilm, Solution, dot, fail, solve_factored,
};

/// Explicit SI law for one original hydraulic channel.
///
/// Shear rate is in 1/s and dynamic viscosity is in Pa s. The consistency's
/// units depend on n. For n=1, consistency is exactly the Newtonian viscosity
/// and the positive core-transition input has no numerical effect.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct RegularizedPowerLaw {
    /// Positive consistency K [Pa s^n].
    pub consistency_pa_s_n: f64,
    /// Positive stress power n; below one thins, above one thickens.
    pub flow_index: f64,
    /// Positive physical transition gamma0 [1/s].
    pub minimum_shear_rate_per_s: f64,
}

/// Explicit numerical work and acceptance budgets for the pressure island.
///
/// Absolute flow tolerances have units m^3/s and apply to the worst cell.
/// Work defects have units W; they are numerical residuals, never extra heat.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct PowerLawControls {
    /// At most 128 pressure corrections. A moving-wall binding also uses this
    /// ceiling for each local scalar stress-offset solve. Zero refuses.
    pub maximum_iterations: usize,
    /// At most 64 candidate trials per correction. Zero refuses.
    pub maximum_line_search_steps: usize,
    /// Nonnegative absolute cell-volume-flow defect [m^3/s].
    pub absolute_flow_tolerance_m3_s: f64,
    /// Relative cell defect against the largest original swept-volume rate.
    /// A pressure-fed film also includes its largest actual passage flow.
    /// Moving-wall films use reduced pressure forcing and pressure-only flux,
    /// retaining small pressure response under large divergence-free transport.
    pub relative_flow_tolerance: f64,
    /// Nonnegative absolute force/flow work defect [W].
    pub absolute_power_tolerance_w: f64,
    /// Relative work defect against max(port work, channel work). A supplied
    /// film includes the magnitudes of mechanical and reservoir input powers.
    pub relative_power_tolerance: f64,
}

impl PowerLawControls {
    fn validate(self) -> Result<(), FilmError> {
        if self.maximum_iterations == 0
            || self.maximum_iterations > 128
            || self.maximum_line_search_steps == 0
            || self.maximum_line_search_steps > 64
            || !self.absolute_flow_tolerance_m3_s.is_finite()
            || self.absolute_flow_tolerance_m3_s < 0.0
            || !self.absolute_power_tolerance_w.is_finite()
            || self.absolute_power_tolerance_w < 0.0
            || !self.relative_flow_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_flow_tolerance)
            || !self.relative_power_tolerance.is_finite()
            || !(0.0..1.0).contains(&self.relative_power_tolerance)
        {
            return fail("invalid power-law film work or residual budgets");
        }
        Ok(())
    }

    fn admits_work(self, work: f64, power: f64) -> bool {
        let allowance =
            self.absolute_power_tolerance_w + self.relative_power_tolerance * work.abs().max(power);
        allowance.is_finite()
            && work.is_finite()
            && work >= 0.0
            && power.is_finite()
            && power >= 0.0
            && (work - power).abs() <= allowance
    }

    fn admits_supplied_work(self, work: f64, state: &State, supplied: bool) -> bool {
        if !supplied {
            return self.admits_work(work, state.power);
        }
        let magnitude = (work.abs() + state.supply_magnitude).max(state.power);
        let residual = work + state.reservoir_power - state.power;
        let allowance = self.absolute_power_tolerance_w + self.relative_power_tolerance * magnitude;
        work.is_finite()
            && magnitude.is_finite()
            && residual.is_finite()
            && allowance.is_finite()
            && state.power.is_finite()
            && state.power >= 0.0
            && residual.abs() <= allowance
    }
}

/// Original film diagnostics plus the actual nonlinear acceptance residuals.
#[derive(Debug, Clone, Copy)]
pub struct PowerLawFilmReport {
    /// Native gap, pressure and physical dissipation diagnostics.
    pub film: FilmReport,
    /// Number of nonlinear corrections; zero for the native n=1 image.
    pub iterations: usize,
    /// Rejected line-search candidates across those corrections.
    pub backtracks: usize,
    /// Largest absolute div(Q)-swept-volume defect over the original cells.
    pub maximum_flow_defect_m3_s: f64,
    /// Signed mechanical plus supplied reservoir work minus dissipation [W].
    /// Reservoir work is zero for the original passive owner.
    pub work_defect_w: f64,
}

#[derive(Debug, Clone, Copy)]
struct CompiledLaw {
    inverse_n: f64,
    exponent: f64,
    transition_stress_pa: f64,
    log_transition_stress: f64,
    log_minimum_shear: f64,
    newtonian: bool,
}

impl RegularizedPowerLaw {
    fn compile(self) -> Result<(CompiledLaw, f64), FilmError> {
        if [
            self.consistency_pa_s_n,
            self.flow_index,
            self.minimum_shear_rate_per_s,
        ]
        .iter()
        .any(|x| !x.is_finite() || *x <= 0.0)
        {
            return fail(
                "power-law film needs positive finite SI consistency, index and core shear",
            );
        }
        if self.flow_index == 1.0 {
            return Ok((
                CompiledLaw {
                    inverse_n: 1.0,
                    exponent: 3.0,
                    transition_stress_pa: 0.0,
                    log_transition_stress: 0.0,
                    log_minimum_shear: 0.0,
                    newtonian: true,
                },
                self.consistency_pa_s_n,
            ));
        }
        let inverse_n = 1.0 / self.flow_index;
        let exponent = 2.0 + inverse_n;
        if !inverse_n.is_finite() || inverse_n <= 0.0 || !exponent.is_finite() {
            return fail("power-law stress exponent is unrepresentable");
        }
        let log_minimum_shear = det::ln(self.minimum_shear_rate_per_s);
        let core_viscosity = positive_or_log(
            self.consistency_pa_s_n
                * det::pow(self.minimum_shear_rate_per_s, self.flow_index - 1.0),
            || det::ln(self.consistency_pa_s_n) + (self.flow_index - 1.0) * log_minimum_shear,
        )?;
        let transition_stress_pa = core_viscosity * self.minimum_shear_rate_per_s;
        if !transition_stress_pa.is_finite() || transition_stress_pa <= 0.0 {
            return fail("power-law Newtonian-core stress is unrepresentable");
        }
        Ok((
            CompiledLaw {
                inverse_n,
                exponent,
                transition_stress_pa,
                log_transition_stress: det::ln(transition_stress_pa),
                log_minimum_shear,
                newtonian: false,
            },
            core_viscosity,
        ))
    }
}

/// Generalized-Newtonian pressure resistance with the original film geometry.
///
/// Compilation shares the admitted native geometry. Each evaluation uses
/// bounded stack scratch and performs no heap allocation. Nonlinear work is
/// capped by explicit controls. Checked entry points let the caller poll
/// cancellation between channel batches, pressure solves and trial states;
/// the unchecked conveniences remain bounded synchronous scalar leaves.
///
/// Every output is staged until pressure, local flow and mechanical work
/// admission succeeds. Laws, geometry and prior caller outputs never mutate
/// on a failed candidate.
#[derive(Debug, Clone)]
pub struct PowerLawFilm {
    film: ResistiveFilm,
    laws: Vec<RegularizedPowerLaw>,
    compiled: Vec<CompiledLaw>,
    controls: PowerLawControls,
    all_newtonian: bool,
}

impl PowerLawFilm {
    /// Bind one SI constitutive law to every ORIGINAL native channel.
    ///
    /// Existing pressure-dependent coefficients refuse. This binding replaces
    /// each uniform channel coefficient with the law's Newtonian-core value,
    /// sharing geometry and preserving the original gap/pressure domains.
    pub fn new(
        film: &ResistiveFilm,
        laws: Vec<RegularizedPowerLaw>,
        controls: PowerLawControls,
    ) -> Result<Self, FilmError> {
        controls.validate()?;
        if film.pressure_viscosity_coefficient_per_pa() != 0.0 {
            return fail("power-law film requires pressure-independent constitutive laws");
        }
        if laws.len() != film.channel_count() {
            return fail("power-law film needs one law per original channel");
        }
        let mut compiled = Vec::with_capacity(laws.len());
        let mut core_viscosities = Vec::with_capacity(laws.len());
        for law in &laws {
            let (prepared, viscosity) = law.compile()?;
            compiled.push(prepared);
            core_viscosities.push(viscosity);
        }
        let all_newtonian = compiled.iter().all(|law| law.newtonian);
        Ok(Self {
            film: film.with_channel_viscosities(&core_viscosities)?,
            laws,
            compiled,
            controls,
            all_newtonian,
        })
    }

    /// Exact original SI laws in channel order.
    #[must_use]
    pub fn laws(&self) -> &[RegularizedPowerLaw] {
        &self.laws
    }

    /// Numerical acceptance and work budgets, held fixed for every evaluation.
    #[must_use]
    pub const fn controls(&self) -> PowerLawControls {
        self.controls
    }

    /// Original geometry with each channel's NEWTONIAN-CORE viscosity.
    ///
    /// This view is insufficient for non-Newtonian mechanics. Use this
    /// binding's evaluations and tangent for the actual shear-dependent law.
    #[must_use]
    pub const fn pressure_operator(&self) -> &ResistiveFilm {
        &self.film
    }

    // Native mass-conserving films share this constitutive leaf; pressure-only
    // flow remains separate from transported liquid content and mean motion.
    pub(super) fn passage_pressure_response(
        &self,
        index: usize,
        gap: f64,
        core_conductance: f64,
        pressure_drop: f64,
        relative_wall_speed: f64,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(f64, f64), FilmError> {
        let channel = if relative_wall_speed == 0.0 {
            self.channel(index, gap, core_conductance, pressure_drop)?
        } else {
            if gap <= 0.0 {
                return fail("moving power-law passage needs a positive gap");
            }
            shear::channel::evaluate(
                self,
                index,
                gap,
                core_conductance,
                pressure_drop,
                shear::Motion {
                    mean: 0.0,
                    relative: relative_wall_speed,
                    entrainment: 0.0,
                },
                check,
            )?
        };
        Ok((channel.flow, channel.differential))
    }

    pub(super) const fn is_newtonian(&self) -> bool {
        self.all_newtonian
    }

    /// Evaluate the complete nonlinear pressure resistance at one configuration.
    pub fn evaluate_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
    ) -> Result<PowerLawFilmReport, FilmError> {
        self.evaluate::<false>(q, velocity, force, pressure, &mut [], &mut || Ok(()))
    }

    /// Return original-order channel flow and physical work from the same solve.
    pub fn evaluate_with_channels_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
    ) -> Result<PowerLawFilmReport, FilmError> {
        self.evaluate::<true>(q, velocity, force, pressure, channels, &mut || Ok(()))
    }

    /// Same complete solve, with caller-owned cancellation/work checks.
    /// A callback refusal, including the final prepublication check, leaves
    /// force, pressure and all channel reports untouched.
    pub fn evaluate_with_channels_and_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<PowerLawFilmReport, FilmError> {
        self.evaluate::<true>(q, velocity, force, pressure, channels, &mut check)
    }

    fn evaluate<const REPORT_CHANNELS: bool>(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<PowerLawFilmReport, FilmError> {
        self.evaluate_with_reservoirs::<REPORT_CHANNELS>(
            q, velocity, force, pressure, channels, None, check,
        )
        .map(|report| report.power_law)
    }

    #[allow(clippy::too_many_arguments)]
    fn evaluate_with_reservoirs<const REPORT_CHANNELS: bool>(
        &self,
        q: &[f64],
        velocity: &[f64],
        force: &mut [f64],
        pressure: &mut [f64],
        channels: &mut [FilmChannelReport],
        reservoirs: Option<&ReservoirPressures>,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<PressureDrivenPowerLawReport, FilmError> {
        if force.len() != self.film.ports
            || pressure.len() != self.film.cell_count()
            || (REPORT_CHANNELS && channels.len() != self.film.channel_count())
        {
            return fail("invalid power-law film output shape");
        }
        let state = self.solve_with_reservoirs(q, velocity, reservoirs, check)?;
        let mut channel_candidate = [FilmChannelReport::default(); MAX_CHANNELS];
        let work = state.port_work;
        if REPORT_CHANNELS {
            for (k, edge) in self.film.channels.iter().enumerate() {
                let dp = state.system.potential[edge.from]
                    - edge.to.map_or_else(
                        || reservoirs.map_or(0.0, |r| r.potential[k]),
                        |j| state.system.potential[j],
                    );
                channel_candidate[k] = FilmChannelReport {
                    volume_flow_m3_s: state.channels[k].flow,
                    dissipated_power_w: state.channels[k].flow * dp,
                };
            }
        }
        let report = PowerLawFilmReport {
            film: FilmReport {
                minimum_gap_m: state.system.minimum_gap,
                peak_pressure_pa: state.system.pressure[..self.film.cell_count()]
                    .iter()
                    .fold(0.0_f64, |m, p| m.max(p.abs())),
                dissipated_power_w: state.power,
            },
            iterations: state.iterations,
            backtracks: state.backtracks,
            maximum_flow_defect_m3_s: state.maximum_residual,
            work_defect_w: work + state.reservoir_power - state.power,
        };
        check()?;
        force.copy_from_slice(&state.force[..self.film.ports]);
        pressure.copy_from_slice(&state.system.pressure[..self.film.cell_count()]);
        if REPORT_CHANNELS {
            channels.copy_from_slice(&channel_candidate[..self.film.channel_count()]);
        }
        Ok(PressureDrivenPowerLawReport {
            power_law: report,
            mechanical_power_w: work,
            reservoir_power_w: state.reservoir_power,
        })
    }

    /// Analytic directional action for independently supplied q and v directions.
    ///
    /// If J is the graph of accepted positive dQ/d(delta_p), then
    /// J*delta_p=A*B*delta_v-div((dQ/dh)*delta_h). The already factored J is
    /// reused. At exact channel closure the finite-core h_+^3 derivative is
    /// zero. The stationary constitutive profile has no additional state.
    pub fn tangent_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        out: &mut [f64],
    ) -> Result<(), FilmError> {
        self.tangent_with_check_into(q, velocity, dq, dv, out, || Ok(()))
    }

    /// The analytic tangent with caller checks and atomic output publication.
    pub fn tangent_with_check_into(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        out: &mut [f64],
        mut check: impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(), FilmError> {
        self.tangent_with_reservoirs(q, velocity, dq, dv, out, None, &[], &mut check)
    }

    #[allow(clippy::too_many_arguments)]
    fn tangent_with_reservoirs(
        &self,
        q: &[f64],
        velocity: &[f64],
        dq: &[f64],
        dv: &[f64],
        out: &mut [f64],
        reservoirs: Option<&ReservoirPressures>,
        reservoir_directions_pa: &[f64],
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<(), FilmError> {
        if dq.len() != self.film.ports
            || dv.len() != self.film.ports
            || out.len() != self.film.ports
            || dq.iter().chain(dv).any(|x| !x.is_finite())
        {
            return fail("invalid power-law film tangent shape or direction");
        }
        if reservoirs.is_some()
            && (reservoir_directions_pa.len() != self.film.channel_count()
                || reservoir_directions_pa.iter().any(|x| !x.is_finite())
                || self
                    .film
                    .channels
                    .iter()
                    .zip(reservoir_directions_pa)
                    .any(|(edge, &dp)| edge.to.is_some() && dp != 0.0))
        {
            return fail("invalid power-law film reservoir pressure direction");
        }
        check()?;
        let state = self.solve_with_reservoirs(q, velocity, reservoirs, check)?;
        let mut candidate = [0.0; MAX_PORTS];
        if self.all_newtonian && reservoirs.is_none() {
            self.film.tangent_from_solution(
                q,
                &state.system,
                dq,
                dv,
                &mut candidate[..self.film.ports],
            )?;
        } else {
            let mut rhs = [0.0; MAX_CELLS];
            for (i, cell) in self.film.cells.iter().enumerate() {
                rhs[i] = cell.area_m2 * dot(&cell.gap.closure, dv);
            }
            for (k, edge) in self.film.channels.iter().enumerate() {
                if k.is_multiple_of(32) {
                    check()?;
                }
                let dh = -dot(&edge.gap.closure, dq);
                let change = state.channels[k].aperture_derivative * dh;
                rhs[edge.from] -= change;
                if let Some(j) = edge.to {
                    rhs[j] += change;
                } else if reservoirs.is_some() {
                    rhs[edge.from] += state.channels[k].differential * reservoir_directions_pa[k];
                }
            }
            solve_factored(
                &state.system.lower,
                &state.system.scale,
                self.film.cell_count(),
                &mut rhs,
            )?;
            self.film.project(&rhs, &mut candidate);
        }
        if candidate[..self.film.ports].iter().any(|x| !x.is_finite()) {
            return fail("nonfinite power-law film tangent");
        }
        check()?;
        out.copy_from_slice(&candidate[..self.film.ports]);
        Ok(())
    }

    fn channel(
        &self,
        index: usize,
        gap: f64,
        core_conductance: f64,
        dp: f64,
    ) -> Result<Channel, FilmError> {
        if gap <= 0.0 {
            return Ok(Channel::default());
        }
        let prepared = self.compiled[index];
        let edge = &self.film.channels[index];
        let magnitude = dp.abs();
        let mut wall_stress = 0.0;
        if !prepared.newtonian {
            wall_stress = magnitude * (0.5 * gap / edge.length_m);
            if !wall_stress.is_finite() || (magnitude > 0.0 && wall_stress == 0.0) {
                wall_stress = det::exp(
                    det::ln(magnitude) + det::ln(gap)
                        - std::f64::consts::LN_2
                        - det::ln(edge.length_m),
                );
            }
        }
        if prepared.newtonian || wall_stress <= prepared.transition_stress_pa {
            let flow = core_conductance * dp;
            let aperture_derivative = 3.0 * self.film.factors[index] * gap * gap * dp;
            if !flow.is_finite() || !aperture_derivative.is_finite() {
                return fail("power-law core flow or gap derivative is unrepresentable");
            }
            return Ok(Channel {
                flow,
                differential: core_conductance,
                aperture_derivative,
                ..Channel::default()
            });
        }
        if !wall_stress.is_finite() || wall_stress <= 0.0 {
            return fail("power-law wall stress is unrepresentable");
        }
        let ratio = wall_stress / prepared.transition_stress_pa;
        let log_ratio = if ratio.is_finite() {
            det::ln(ratio)
        } else {
            det::ln(wall_stress) - prepared.log_transition_stress
        };
        let law = self.laws[index];
        let rate_exponent = prepared.inverse_n * log_ratio;
        let wall_rate = positive_or_log(
            law.minimum_shear_rate_per_s * det::exp(rate_exponent),
            || prepared.log_minimum_shear + rate_exponent,
        )?;
        let core_exponent = -prepared.exponent * log_ratio;
        let core = det::exp(core_exponent);
        let outside = -det::expm1(core_exponent);
        let flux_fraction = outside / prepared.exponent + core / 3.0;
        let differential_fraction = outside * (prepared.inverse_n / prepared.exponent) + core / 3.0;
        let base = (0.5 * edge.width_m * gap * gap) * wall_rate;
        let log_base = || {
            -std::f64::consts::LN_2
                + det::ln(edge.width_m)
                + 2.0 * det::ln(gap)
                + det::ln(wall_rate)
        };
        let flow = positive_or_log(base * flux_fraction, || log_base() + det::ln(flux_fraction))?;
        let differential = positive_or_log((base / magnitude) * differential_fraction, || {
            log_base() - det::ln(magnitude) + det::ln(differential_fraction)
        })?;
        let aperture_derivative = positive_or_log(base / gap, || log_base() - det::ln(gap))?;
        Ok(Channel {
            flow: flow.copysign(dp),
            differential,
            aperture_derivative: aperture_derivative.copysign(dp),
            ..Channel::default()
        })
    }
}

impl ReynoldsLaw for PowerLawFilm {
    fn pressure_operator(&self) -> &ResistiveFilm {
        &self.film
    }
    fn pressure_controls(&self) -> PowerLawControls {
        self.controls
    }
    fn native_newtonian(&self) -> bool {
        self.all_newtonian
    }
    fn core_requires_positive_conductance(&self, index: usize) -> bool {
        !self.compiled[index].newtonian
    }
    fn channel_response(
        &self,
        index: usize,
        gap: f64,
        core_conductance: f64,
        pressure_drop: f64,
        motion: shear::Motion,
        moving: bool,
        check: &mut impl FnMut() -> Result<(), FilmError>,
    ) -> Result<Channel, FilmError> {
        if moving {
            shear::channel::evaluate(
                self,
                index,
                gap,
                core_conductance,
                pressure_drop,
                motion,
                check,
            )
        } else {
            self.channel(index, gap, core_conductance, pressure_drop)
        }
    }
}

#[derive(Clone, Copy)]
struct ChannelWork {
    dissipated_power_w: f64,
    constitutive_work_defect_w: f64,
}

#[derive(Clone, Copy, Default)]
struct Channel {
    flow: f64,
    differential: f64,
    aperture_derivative: f64,
    wall: shear::WallResponse,
    work: Option<ChannelWork>,
}

struct Geometry {
    gaps: [f64; MAX_CHANNELS],
    core_conductance: [f64; MAX_CHANNELS],
    source: [f64; MAX_CELLS],
    pressure_source: [f64; MAX_CELLS],
    source_scale: f64,
    minimum_gap: f64,
    motion: [shear::Motion; MAX_CHANNELS],
}

struct State {
    system: Solution,
    channels: [Channel; MAX_CHANNELS],
    residual: [f64; MAX_CELLS],
    maximum_residual: f64,
    source_work: f64,
    wall_work: f64,
    reservoir_power: f64,
    supply_magnitude: f64,
    flow_scale: f64,
    power: f64,
    force: [f64; MAX_PORTS],
    port_work: f64,
    iterations: usize,
    backtracks: usize,
}

impl State {
    fn new(n: usize, minimum_gap: f64) -> Self {
        Self {
            system: Solution::new(n, minimum_gap),
            channels: [Channel::default(); MAX_CHANNELS],
            residual: [0.0; MAX_CELLS],
            maximum_residual: 0.0,
            source_work: 0.0,
            wall_work: 0.0,
            reservoir_power: 0.0,
            supply_magnitude: 0.0,
            flow_scale: 0.0,
            power: 0.0,
            force: [0.0; MAX_PORTS],
            port_work: 0.0,
            iterations: 0,
            backtracks: 0,
        }
    }
}

/// Fast ordinary arithmetic is retained when its positive result is finite.
/// Only an exceptional intermediate uses the same deterministic logarithms.
fn positive_or_log(value: f64, logarithm: impl FnOnce() -> f64) -> Result<f64, FilmError> {
    if value.is_finite() && value > 0.0 {
        return Ok(value);
    }
    let result = det::exp(logarithm());
    if !result.is_finite() || result <= 0.0 {
        return fail("power-law film coefficient, rate or differential is unrepresentable");
    }
    Ok(result)
}

#[cfg(test)]
mod tests;
