//! Finite cooling-wall storage with a zero-capable, fully-wetting liquid film.
//!
//! An outer original-wall enthalpy trial and inner sensible-heat trial reuse
//! the shared bounded scalar solver. Each gas trial receives sensible heat,
//! then the original wall-film owner advances phase against that heated gas.
//! Its endpoint affinity guard therefore sees the ACTUAL final gas. The original
//! lumped wall consumes the film/phase heat port and evaluates convection,
//! radiation, source conductivity and Biot admission at the accepted endpoint.
//!
//! Wall material mass and film mass remain distinct. Film h is supplied by its
//! own original chart, never by the wall's heat capacity or caloric datum.
//! Dry inventory has no temperature or thermal storage. The full requested
//! interval is retained when frozen-rate evaporation reaches exact dry-out.

use crate::cavity_phs::piston::thermal::{HeatExchangeSolve, feedback};
use fs_conduction::{
    ConductionError,
    lumped::{LumpedEnthalpyMarch, LumpedThermalEnvironment},
};
use fs_exec::Cx;
use fs_material::{
    evaporation::gas::{
        GasEvaporationError, GasEvaporationSolve, GasLiquidBinding,
        wall_film::{
            GasWallFilmInventory, GasWallFilmStep, GasWallFilmSurface, advance_wall_film_with_gas,
        },
    },
    phase::{UniformEnthalpyStep, UniformEnthalpyStepInput},
    thermochemical_gas::{
        ThermochemicalGasModel, energy::GasHeatStep, energy::parcel::ClosedGasParcel,
    },
};
use fs_thermal::evaporation::{EvaporationError, EvaporationLaw};

/// Original liquid sources for a measured patch on an isothermal finite wall.
/// Both callbacks must describe the SAME original pure-liquid caloric branch.
#[derive(Clone, Copy)]
pub struct WallFilmThermalBoundary<'a> {
    pub law: &'a dyn EvaporationLaw,
    pub binding: &'a GasLiquidBinding,
    pub area_m2: f64,
    pub mass_transfer_coefficient_m_s: f64,
    pub liquid_enthalpy_at_temperature: &'a dyn Fn(f64) -> Result<f64, EvaporationError>,
    pub liquid_temperature_from_enthalpy: &'a dyn Fn(f64) -> Result<f64, EvaporationError>,
}

/// Original source controls and two separately bounded physical feedback roots.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct WallFilmThermalSolve {
    pub phase: GasEvaporationSolve,
    /// Signed sensible heat into gas; every trial must lie inside source support.
    pub sensible: HeatExchangeSolve,
    /// Bracket is original WALL storage change m_wall*(h_trial-h_initial), J.
    /// Its energy allowance bounds the complete wall+film+gas transaction.
    pub wall: HeatExchangeSolve,
    /// Original liquid h(T) round trip and accepted film/wall mismatch [K].
    pub temperature_tolerance_k: f64,
}

/// Repeatable proposal; neither gas, film nor accepted wall is mutated.
#[derive(Clone)]
pub struct WallFilmThermalEnvironment<'a> {
    pub gas: &'a ClosedGasParcel,
    pub catalog: &'a ThermochemicalGasModel,
    pub film: GasWallFilmInventory,
    pub boundary: WallFilmThermalBoundary<'a>,
    /// Its convection temperature is replaced by the actual final gas. Its
    /// explicit h, radiation enclosure, conductivity and Biot policy survive.
    pub body: &'a LumpedThermalEnvironment,
    /// Independent signed heat into the gas over the full requested interval.
    pub gas_heat_j: f64,
    pub solve: WallFilmThermalSolve,
}

#[derive(Debug)]
pub enum WallFilmThermalError {
    Phase(GasEvaporationError),
    Thermal(super::Error),
    Invalid(&'static str),
    MassBalance {
        residual_kg: f64,
        tolerance_kg: f64,
    },
    EnergyBalance {
        residual_j: f64,
        tolerance_j: f64,
    },
    /// Sensible residuals are joules; wall residuals are dimensionless maxima
    /// of the declared energy and temperature consistency allowances.
    OutsideBracket {
        stage: &'static str,
        lower_residual: f64,
        upper_residual: f64,
    },
    NotConverged {
        stage: &'static str,
        iterations: usize,
        residual: f64,
    },
    Cancelled,
}
impl core::fmt::Display for WallFilmThermalError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite wall-film thermal exchange: {self:?}")
    }
}
impl std::error::Error for WallFilmThermalError {}
impl From<GasEvaporationError> for WallFilmThermalError {
    fn from(error: GasEvaporationError) -> Self {
        if error == GasEvaporationError::Cancelled {
            Self::Cancelled
        } else {
            Self::Phase(error)
        }
    }
}
impl From<super::Error> for WallFilmThermalError {
    fn from(error: super::Error) -> Self {
        if matches!(error, super::Error::Cancelled) {
            Self::Cancelled
        } else {
            Self::Thermal(error)
        }
    }
}

/// Original physical receipts and the gas/film candidate to publish with wall.
#[derive(Clone, Debug, PartialEq)]
pub struct WallFilmThermalReport {
    pub next_gas: ClosedGasParcel,
    pub next_film: GasWallFilmInventory,
    /// Phase-only receipt, starting from the gas after sensible heating.
    pub phase: GasWallFilmStep,
    /// Original frozen-composition gas heating, before phase transfer.
    pub gas_heat: GasHeatStep,
    /// Original wall solve; its internal source includes the signed film debit.
    pub body: LumpedEnthalpyMarch,
    pub solid_to_gas_heat_j: f64,
    /// Positive heat leaves the finite wall for film plus phase-changing gas.
    pub wall_heat_into_phase_system_j: f64,
    pub external_gas_heat_j: f64,
    /// External plus sensible gas heat, excluding phase-carried vapor energy.
    pub net_gas_heat_j: f64,
    /// Stored-minus-source defect from combining external and sensible gas heat.
    pub gas_source_combination_residual_j: f64,
    /// Stored-minus-source defect from combining wall and phase source ports.
    pub wall_source_combination_residual_j: f64,
    /// Roundoff between separate boundary ports and returned net wall heat.
    pub external_heat_roundoff_j: f64,
    pub radiation_heat_j: f64,
    /// q plus the original wall owner's actual convective heat into the wall.
    pub source_residual_j: f64,
    /// Actual wall storage minus all declared internal and boundary heat.
    pub wall_energy_residual_j: f64,
    /// Original wall trial minus its independently solved storage change.
    pub wall_trial_energy_residual_j: f64,
    pub wall_temperature_residual_k: f64,
    /// Source-relative gas + film + wall closure; phase reference rounding
    /// remains independently mass-budgeted as in the original species owner.
    pub combined_energy_residual_j: f64,
    /// Absolute phase, gas-heating, wall-solve, source and rounding defects.
    pub sum_absolute_energy_residual_j: f64,
    pub mass_residual_kg: f64,
    pub volatile_mass_residual_kg: f64,
    pub sum_absolute_mass_residual_kg: f64,
    pub trials: usize,
    /// Includes all rejected outer wall proposals.
    pub total_sensible_trials: usize,
}

type Error = WallFilmThermalError;

fn scalar_error(stage: &'static str, error: feedback::ScalarError<Error>) -> Error {
    match error {
        feedback::ScalarError::Evaluation(error) => error,
        feedback::ScalarError::Invalid(message) => Error::Invalid(message),
        feedback::ScalarError::Cancelled => Error::Cancelled,
        feedback::ScalarError::OutsideBracket { lower, upper } => Error::OutsideBracket {
            stage,
            lower_residual: lower,
            upper_residual: upper,
        },
        feedback::ScalarError::NotConverged {
            iterations,
            residual,
        } => Error::NotConverged {
            stage,
            iterations,
            residual,
        },
    }
}
fn finite(value: f64) -> Result<f64, Error> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(Error::Invalid(
            "finite wall-film account is unrepresentable",
        ))
    }
}
fn combine(a: f64, b: f64) -> Result<(f64, f64), Error> {
    let (sum, roundoff) = fs_math::eft::two_sum(a, b);
    let sum = finite(sum)?;
    let roundoff = finite(roundoff)?;
    if (a != 0.0 && sum == b) || (b != 0.0 && sum == a) {
        return Err(Error::Invalid(
            "a nonzero wall-film heat source rounds away",
        ));
    }
    Ok((sum, roundoff))
}

struct SensibleCandidate {
    gas: ClosedGasParcel,
    film: GasWallFilmInventory,
    phase: GasWallFilmStep,
    gas_heat: GasHeatStep,
    gas_source_roundoff: f64,
    q: f64,
    absolute_mass_defect: f64,
}

impl WallFilmThermalEnvironment<'_> {
    /// Propose one full physical interval of finite wall, gas and retained film.
    ///
    /// Heat is applied to the original finite gas before the original phase
    /// owner. Its full-interval frozen-rate closure reports any exact dry-out
    /// time and applies zero evaporation afterward. The same interval drives
    /// the original wall thermal owner; no clock is shortened or advanced twice.
    ///
    /// The sensible coefficient is additional to all vapor-carried latent and
    /// sensible energy. Whole-wall convection uses input.surface_area_m2; the
    /// declared phase patch cannot exceed it. Radiation belongs to its separate
    /// enclosing field. Phase heat is homogenized by the original uniform-wall
    /// model; Biot admission alone does not justify rapid localized heating.
    ///
    /// This is a first-order heat/phase treatment with endpoint wall feedback,
    /// requiring physical-timestep refinement. Brackets must support EVERY
    /// queried gas, liquid and wall state. A trial source refusal is not a proof
    /// that an interior physical solution does not exist. The fully-wetting
    /// planar closure adds no nucleation barrier, moving film geometry, pressure
    /// work, runoff, spray, or independent entropy theorem.
    pub fn advance(
        &self,
        cx: &Cx<'_>,
        input: UniformEnthalpyStepInput<'_>,
    ) -> Result<UniformEnthalpyStep<WallFilmThermalReport>, Error> {
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        feedback::admit::<ConductionError>(self.solve.wall)?;
        feedback::admit::<ConductionError>(self.solve.sensible)?;
        if !self.gas_heat_j.is_finite()
            || !input.mass_kg.is_finite()
            || input.mass_kg <= 0.0
            || !input.duration_s.is_finite()
            || input.duration_s <= 0.0
            || !input.surface_area_m2.is_finite()
            || input.surface_area_m2 <= 0.0
            || !input.internal_heat_j.is_finite()
            || !self.boundary.area_m2.is_finite()
            || self.boundary.area_m2 < 0.0
            || self.boundary.area_m2 > input.surface_area_m2
            || !self.body.convection_w_per_m2_k.is_finite()
            || self.body.convection_w_per_m2_k < 0.0
            || !self.solve.temperature_tolerance_k.is_finite()
            || self.solve.temperature_tolerance_k <= 0.0
            || input.initial.phase_curve_identity() != input.curve.identity()
        {
            return Err(Error::Invalid(
                "invalid finite wall-film geometry, source or temperature controls",
            ));
        }
        let initial_h = input.initial.specific_enthalpy_j_kg();
        let conductance_time =
            finite(self.body.convection_w_per_m2_k * input.surface_area_m2 * input.duration_s)?;
        if self.body.convection_w_per_m2_k > 0.0 && conductance_time == 0.0 {
            return Err(Error::Invalid(
                "nonzero wall convection exposure underflows",
            ));
        }
        let mut total_sensible_trials = 0usize;
        let (candidate, body, wall_trial_energy_residual_j, wall_temperature_residual_k, trials) =
            feedback::solve_scalar(
                self.solve.wall.heat_bracket_j,
                1.0,
                self.solve.wall.maximum_iterations,
                |wall_energy, trials, cancelled| {
                    let trial_h = finite(initial_h + wall_energy / input.mass_kg)?;
                    if wall_energy != 0.0 && trial_h == initial_h {
                        return Err(Error::Invalid(
                            "wall trial energy cannot change its caloric coordinate",
                        ));
                    }
                    let trial_wall = input
                        .curve
                        .state_at_specific_enthalpy(trial_h)
                        .map_err(super::Error::from)?;
                    let wall_t = trial_wall.temperature_k();
                    let liquid_h = (self.boundary.liquid_enthalpy_at_temperature)(wall_t)
                        .map_err(GasEvaporationError::from)?;
                    let liquid_t = (self.boundary.liquid_temperature_from_enthalpy)(liquid_h)
                        .map_err(GasEvaporationError::from)?;
                    if !liquid_h.is_finite()
                        || !liquid_t.is_finite()
                        || (liquid_t - wall_t).abs() > self.solve.temperature_tolerance_k
                    {
                        return Err(Error::Invalid(
                            "original liquid h(T) and T(h) do not match the wall trial",
                        ));
                    }
                    let candidate = feedback::solve_scalar(
                        self.solve.sensible.heat_bracket_j,
                        self.solve.sensible.source_tolerance_j,
                        self.solve.sensible.maximum_iterations,
                        |q, _, cancelled| {
                            total_sensible_trials = total_sensible_trials
                                .checked_add(1)
                                .ok_or(Error::Invalid("wall-film feedback work overflow"))?;
                            let (net_heat, gas_source_roundoff) = combine(self.gas_heat_j, q)?;
                            let mut heated_gas = self.gas.clone();
                            let gas_heat = heated_gas
                                .add_heat(net_heat, self.solve.phase.gas.gas, &mut *cancelled)
                                .map_err(super::gas_error)?;
                            if gas_heat.residual_j.abs() > self.solve.sensible.energy_tolerance_j {
                                return Err(Error::EnergyBalance {
                                    residual_j: gas_heat.residual_j,
                                    tolerance_j: self.solve.sensible.energy_tolerance_j,
                                });
                            }
                            let phase = advance_wall_film_with_gas(
                                GasWallFilmSurface {
                                    law: self.boundary.law,
                                    binding: self.boundary.binding,
                                    film: self.film,
                                    area_m2: self.boundary.area_m2,
                                    interface_liquid_enthalpy_j_kg: liquid_h,
                                    temperature_from_enthalpy: self
                                        .boundary
                                        .liquid_temperature_from_enthalpy,
                                    mass_transfer_coefficient_m_s: self
                                        .boundary
                                        .mass_transfer_coefficient_m_s,
                                },
                                &heated_gas,
                                self.catalog,
                                input.duration_s,
                                self.solve.phase,
                                &mut *cancelled,
                            )?;
                            let absolute_mass_defect = finite(
                                phase.step().liquid_mass_residual_kg.abs()
                                    + phase.gas_exchange().step().mass_residual_kg.abs()
                                    + phase.gas_exchange().step().maximum_species_mass_residual_kg,
                            )?;
                            let residual = finite(
                                q - conductance_time
                                    * (wall_t - phase.gas().state().state().temperature),
                            )?;
                            let (film, gas, _, phase) = phase.into_parts();
                            Ok((
                                SensibleCandidate {
                                    gas,
                                    film,
                                    phase,
                                    gas_heat,
                                    gas_source_roundoff,
                                    q,
                                    absolute_mass_defect,
                                },
                                residual,
                            ))
                        },
                        &mut *cancelled,
                    )
                    .map_err(|error| scalar_error("sensible heat", error))?;
                    let (internal_heat, _) = combine(
                        input.internal_heat_j,
                        -candidate.phase.wall_heat_into_system_j,
                    )?;
                    let mut environment = self.body.clone();
                    environment.temperature_k = candidate.gas.state().state().temperature;
                    let body = environment
                        .advance(
                            cx,
                            UniformEnthalpyStepInput {
                                internal_heat_j: internal_heat,
                                ..input
                            },
                        )
                        .map_err(|error| Error::from(super::Error::Interface(error)))?;
                    let residual =
                        finite(input.mass_kg * (trial_h - body.state.specific_enthalpy_j_kg()))?;
                    let temperature_residual =
                        finite(candidate.phase.wall_temperature_k - body.state.temperature_k())?;
                    let endpoint = body.report.samples().last().ok_or(Error::Invalid(
                        "original thermal owner returned no wall endpoint",
                    ))?;
                    let source_residual =
                        finite(candidate.q + input.duration_s * endpoint.convection_into_body_w)?;
                    let scaled = finite(
                        (residual.abs() / self.solve.wall.source_tolerance_j)
                            .max(temperature_residual.abs() / self.solve.temperature_tolerance_k)
                            .max(source_residual.abs() / self.solve.sensible.source_tolerance_j),
                    )?;
                    let signed = scaled.copysign(residual);
                    Ok((
                        (candidate, body, residual, temperature_residual, trials),
                        signed,
                    ))
                },
                &mut || cx.checkpoint().is_err(),
            )
            .map_err(|error| scalar_error("wall storage and temperature", error))?;

        let wall_t = body.state.temperature_k();
        let gas_t = candidate.gas.state().state().temperature;
        let maximum_increment = self.solve.phase.maximum_temperature_increment_k;
        if (wall_t - input.initial.temperature_k()).abs() > maximum_increment
            || (gas_t - self.gas.state().state().temperature).abs() > maximum_increment
        {
            return Err(GasEvaporationError::StepSizeLimit {
                maximum_duration_s: (input.duration_s * 0.5).next_down().max(0.0),
                reason: "complete finite wall-film temperature departure",
            }
            .into());
        }
        if (candidate.q > 0.0 && wall_t < gas_t) || (candidate.q < 0.0 && wall_t > gas_t) {
            return Err(Error::Invalid(
                "sensible heat disagrees with the actual endpoint temperature ordering",
            ));
        }
        let sample = body.report.samples().last().ok_or(Error::Invalid(
            "original thermal owner returned no wall endpoint",
        ))?;
        let convection_heat = finite(input.duration_s * sample.convection_into_body_w)?;
        let radiation_heat_j = finite(input.duration_s * sample.radiation_into_body_w)?;
        let source_residual_j = finite(candidate.q + convection_heat)?;
        if source_residual_j.abs() > self.solve.sensible.source_tolerance_j {
            return Err(Error::EnergyBalance {
                residual_j: source_residual_j,
                tolerance_j: self.solve.sensible.source_tolerance_j,
            });
        }
        let wall_heat_into_phase_system_j = candidate.phase.wall_heat_into_system_j;
        let (external_sensible, external_roundoff_a) =
            fs_math::eft::two_sum(radiation_heat_j, -candidate.q);
        let (external_heat_j, external_roundoff_b) =
            fs_math::eft::two_sum(external_sensible, -wall_heat_into_phase_system_j);
        let external_heat_j = finite(external_heat_j)?;
        let external_heat_roundoff_j = finite(external_roundoff_a + external_roundoff_b)?;
        let body_storage =
            finite(input.mass_kg * (body.state.specific_enthalpy_j_kg() - initial_h))?;
        let (internal_heat, wall_source_roundoff) =
            combine(input.internal_heat_j, -wall_heat_into_phase_system_j)?;
        let body_residual = finite(body_storage - internal_heat - body.external_heat_j)?;
        let gas_source_combination_residual_j = -candidate.gas_source_roundoff;
        let wall_source_combination_residual_j = -wall_source_roundoff;
        let wall_energy_residual_j =
            finite(body_residual + source_residual_j + wall_source_combination_residual_j)?;
        let combined_energy_residual_j = finite(
            wall_energy_residual_j
                + candidate.phase.energy_residual_j
                + candidate.gas_heat.residual_j
                + gas_source_combination_residual_j,
        )?;
        let sum_absolute_energy_residual_j = finite(
            body_residual.abs()
                + source_residual_j.abs()
                + candidate.phase.sum_absolute_energy_residual_j
                + candidate.gas_heat.residual_j.abs()
                + gas_source_combination_residual_j.abs()
                + wall_source_combination_residual_j.abs()
                + external_heat_roundoff_j.abs(),
        )?;
        let liquid_change = candidate.film.mass_kg() - self.film.mass_kg();
        let mass_residual_kg =
            finite(liquid_change + (candidate.gas.mass_kg() - self.gas.mass_kg()))?;
        let volatile = &self.boundary.binding.species;
        let volatile_mass_residual_kg = finite(
            liquid_change
                + candidate
                    .gas
                    .species_mass_kg(volatile)
                    .map_err(GasEvaporationError::from)?
                - self
                    .gas
                    .species_mass_kg(volatile)
                    .map_err(GasEvaporationError::from)?,
        )?;
        let worst_mass = mass_residual_kg
            .abs()
            .max(volatile_mass_residual_kg.abs())
            .max(candidate.absolute_mass_defect);
        if worst_mass > self.solve.phase.gas.mass_tolerance_kg {
            return Err(Error::MassBalance {
                residual_kg: worst_mass,
                tolerance_kg: self.solve.phase.gas.mass_tolerance_kg,
            });
        }
        let worst_energy = combined_energy_residual_j
            .abs()
            .max(wall_energy_residual_j.abs())
            .max(sum_absolute_energy_residual_j);
        if worst_energy > self.solve.wall.energy_tolerance_j {
            return Err(Error::EnergyBalance {
                residual_j: worst_energy,
                tolerance_j: self.solve.wall.energy_tolerance_j,
            });
        }
        cx.checkpoint().map_err(|_| Error::Cancelled)?;
        Ok(UniformEnthalpyStep {
            state: body.state,
            external_heat_j,
            energy_residual_tolerance_j: self.solve.wall.energy_tolerance_j,
            report: WallFilmThermalReport {
                next_gas: candidate.gas,
                next_film: candidate.film,
                phase: candidate.phase,
                gas_heat: candidate.gas_heat,
                body: body.report,
                solid_to_gas_heat_j: candidate.q,
                wall_heat_into_phase_system_j,
                external_gas_heat_j: self.gas_heat_j,
                net_gas_heat_j: candidate.gas_heat.supplied_heat_j,
                gas_source_combination_residual_j,
                wall_source_combination_residual_j,
                external_heat_roundoff_j,
                radiation_heat_j,
                source_residual_j,
                wall_energy_residual_j,
                wall_trial_energy_residual_j,
                wall_temperature_residual_k,
                combined_energy_residual_j,
                sum_absolute_energy_residual_j,
                mass_residual_kg,
                volatile_mass_residual_kg,
                sum_absolute_mass_residual_kg: candidate.absolute_mass_defect,
                trials,
                total_sensible_trials,
            },
        })
    }
}
