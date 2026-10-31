//! Source-backed squeeze pressure on the original articulated/thermal timeline.
//!
//! The native Reynolds graph supplies BOTH `B^T A p` and its nonnegative power
//! at each actual generalized position/rate. Its positive resistance is
//! subtracted from the native ABA efforts. Geometry is the explicitly authored
//! linearized gap map in the mechanism's full scalar-coordinate ordering;
//! density/inertia, pressure and viscous loss are not replaced by fitted damping.
//!
//! The pressure core is normally embedded implicitly. This adapter instead
//! admits its nonstiff ODE image through the existing bounded adaptive RK45
//! stages, integrated loss coordinate and independent mechanical-energy gate.
//! It makes no stiff-film stability claim: unresolved dynamics, collapsed gaps,
//! trapped pockets and pressure/source limits refuse without a candidate.
//!
//! One finite receiving wall sets the massless, isothermal film temperature and
//! receives its loss once. Original coefficients may depend only on that T.
//! Their pressure support must cover the declared small gauge-pressure range.
//! Absolute pressure must exceed the SAME bundle's saturation pressure at every
//! stage and refreshed endpoint. This conservative full-liquid gate is not a
//! cavitation-inception prediction or the separate Elrod--Adams content solver.
//! Liquid inertia/storage, compressibility, phase change, entrainment, evolving
//! surface geometry and non-Newtonian response are outside this adapter.

use super::{Error, finite, nonnegative, positive};
use fs_matdb::{ClaimSelection, QueryPoint};
use fs_material::liquid::ResolvedLiquidState;
use fs_material::liquid_temperature::LiquidTemperatureLaw;
use fs_material::state_point::MaterialStatePointError;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Pressure, QuantitySpec, Temperature};
use fs_thermal::enthalpy::{CellId, DepositionError, DepositionSimulation};
use fs_tribo::resistive_film::{FilmReport, ResistiveFilm};

/// Physical pressure envelope and independent coefficient/integration controls.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SqueezeLimits {
    /// Maximum `abs(gauge_pressure) / ambient_pressure`, in (0, 0.1].
    /// All original property claims must support this complete pressure range.
    pub maximum_relative_gauge_pressure: f64,
    /// Maximum fractional viscosity change in one accepted thermal substep.
    pub maximum_relative_viscosity_change: f64,
    /// Absolute tolerance on this film's integrated native viscous loss [J].
    pub absolute_loss_tolerance_j: f64,
}

/// One original pressure graph and source liquid, with one finite heat receiver.
#[derive(Clone, Debug)]
pub struct SqueezeFilm<'a> {
    /// Existing finite caloric cell; all film heat and its temperature belong here.
    pub heat_cell: CellId,
    /// Immutable cells/passages and their work-conjugate generalized gap map.
    /// Its authored viscosity is replaced by the selected original liquid value.
    pub geometry: ResistiveFilm,
    /// Original full source bundle, including the explicitly named p_sat claim.
    pub liquid: LiquidTemperatureLaw<'a>,
    /// Existing semantic Pressure/Static coordinate, on an absolute reference.
    /// Its retained value is the pressure of every exterior film reservoir.
    pub pressure_axis: &'a str,
    /// Pressure-valued equilibrium vapor-pressure claim in the SAME bundle.
    /// Construct that bundle with `LiquidTemperatureLaw::try_from_material`.
    pub saturation_pressure_property: &'a str,
    pub limits: SqueezeLimits,
}

/// Accepted source, pressure and force observations in the original port order.
#[derive(Clone, Debug, PartialEq)]
pub struct SqueezeFilmState {
    pub heat_cell: CellId,
    pub liquid: ResolvedLiquidState,
    /// Positive resisting efforts; the native mechanism receives their negatives.
    /// Units follow each original joint: N for translation, N m for rotation.
    pub generalized_resistance: Vec<f64>,
    pub gauge_pressure_pa: Vec<f64>,
    pub ambient_pressure_pa: f64,
    pub saturation_pressure_pa: f64,
    pub minimum_gap_m: f64,
    pub dissipated_power_w: f64,
    /// Cumulative accepted work already deposited into the named finite cell.
    pub dissipated_heat_j: f64,
}

#[derive(Clone)]
pub(super) struct BoundSqueeze<'a> {
    pub definition: SqueezeFilm<'a>,
    pub state: SqueezeFilmState,
    film: ResistiveFilm,
}

impl<'a> BoundSqueeze<'a> {
    pub fn new(
        definition: SqueezeFilm<'a>,
        thermal: &DepositionSimulation<'_>,
        position: &[f64],
        velocity: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        let limits = definition.limits;
        if definition.geometry.port_count() != position.len()
            || position.len() != velocity.len()
            || !positive(limits.maximum_relative_gauge_pressure)
            || limits.maximum_relative_gauge_pressure > 0.1
            || !positive(limits.maximum_relative_viscosity_change)
            || limits.maximum_relative_viscosity_change > 1.0
            || !positive(limits.absolute_loss_tolerance_j)
        {
            return Err(Error::Invalid(
                "squeeze film needs the full joint port map and positive bounded pressure/viscosity/loss controls",
            ));
        }
        let reference = definition.liquid.reference_state();
        let axis = definition.pressure_axis;
        let point = reference.material().query();
        let pressure =
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static));
        if axis == definition.liquid.temperature_axis()
            || point.axis_quantities().get(axis) != Some(&pressure)
            || reference.material().properties().iter().any(|property| {
                !property.is_constant_at_fixed_temperature()
                    || !property.has_only_fixed_external_conditions(
                        Some(definition.liquid.temperature_axis()),
                        false,
                    )
            })
        {
            return Err(Error::Invalid(
                "squeeze film requires absolute reservoir pressure and constant or single-temperature Newtonian source claims",
            ));
        }
        let ambient = *point.axes().get(axis).ok_or(Error::Invalid(
            "squeeze film reservoir pressure coordinate is absent",
        ))?;
        if !positive(ambient) {
            return Err(Error::Invalid(
                "squeeze film reservoir pressure must be positive absolute Pa",
            ));
        }
        let liquid = definition.liquid.at_temperature(Temperature::new(
            thermal.temperature_k(definition.heat_cell)?,
        ))?;
        temperature_support(&definition.liquid, reference, &liquid, checkpoint)?;
        pressure_support(&definition, &liquid, ambient, checkpoint)?;
        let saturation = saturation_pressure(&definition, &liquid)?;
        let film = definition
            .geometry
            .with_viscosity(liquid.dynamic_viscosity_pa_s())?;
        let mut result = Self {
            state: SqueezeFilmState {
                heat_cell: definition.heat_cell,
                liquid,
                generalized_resistance: vec![0.0; film.port_count()],
                gauge_pressure_pa: vec![0.0; film.cell_count()],
                ambient_pressure_pa: ambient,
                saturation_pressure_pa: saturation,
                minimum_gap_m: 0.0,
                dissipated_power_w: 0.0,
                dissipated_heat_j: 0.0,
            },
            definition,
            film,
        };
        result.update_motion(position, velocity)?;
        checkpoint()?;
        Ok(result)
    }

    /// One native pressure solution provides both the mechanical port and heat.
    pub fn response(
        &self,
        position: &[f64],
        velocity: &[f64],
    ) -> Result<(Vec<f64>, Vec<f64>, FilmReport), Error> {
        let mut force = vec![0.0; self.film.port_count()];
        let mut pressure = vec![0.0; self.film.cell_count()];
        let report = self
            .film
            .evaluate_into(position, velocity, &mut force, &mut pressure)?;
        let ambient = self.state.ambient_pressure_pa;
        let saturation = self.state.saturation_pressure_pa;
        if ambient <= saturation {
            return Err(Error::Invalid(
                "squeeze film reservoir is outside the source full-liquid regime",
            ));
        }
        for &gauge in &pressure {
            let absolute = finite(ambient + gauge)?;
            if !positive(absolute) || absolute <= saturation {
                return Err(Error::Invalid(
                    "squeeze film reaches saturation; a mass-conserving cavitation owner is required",
                ));
            }
            if gauge.abs() > self.definition.limits.maximum_relative_gauge_pressure * ambient {
                return Err(Error::Budget(
                    "squeeze film exceeds its supported low-pressure envelope",
                ));
            }
        }
        Ok((force, pressure, report))
    }

    fn update_motion(&mut self, position: &[f64], velocity: &[f64]) -> Result<(), Error> {
        let (force, pressure, report) = self.response(position, velocity)?;
        self.state.generalized_resistance = force;
        self.state.gauge_pressure_pa = pressure;
        self.state.minimum_gap_m = report.minimum_gap_m;
        self.state.dissipated_power_w = report.dissipated_power_w;
        Ok(())
    }

    pub fn accept(
        &mut self,
        loss_j: f64,
        thermal: &DepositionSimulation<'_>,
        position: &[f64],
        velocity: &[f64],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        checkpoint()?;
        if !nonnegative(loss_j) {
            return Err(Error::Budget(
                "squeeze film integrated loss must be nonnegative",
            ));
        }
        let liquid = self.definition.liquid.at_temperature(Temperature::new(
            thermal.temperature_k(self.definition.heat_cell)?,
        ))?;
        temperature_support(
            &self.definition.liquid,
            &self.state.liquid,
            &liquid,
            checkpoint,
        )?;
        pressure_support(
            &self.definition,
            &liquid,
            self.state.ambient_pressure_pa,
            checkpoint,
        )?;
        let old_mu = self.state.liquid.dynamic_viscosity_pa_s();
        let change = finite((liquid.dynamic_viscosity_pa_s() - old_mu).abs() / old_mu)?;
        if change > self.definition.limits.maximum_relative_viscosity_change {
            return Err(Error::Budget(
                "squeeze viscosity change requires a smaller thermal substep",
            ));
        }
        self.film = self
            .definition
            .geometry
            .with_viscosity(liquid.dynamic_viscosity_pa_s())?;
        self.state.saturation_pressure_pa = saturation_pressure(&self.definition, &liquid)?;
        self.state.liquid = liquid;
        self.update_motion(position, velocity)?;
        self.state.dissipated_heat_j = finite(self.state.dissipated_heat_j + loss_j)?;
        checkpoint()?;
        Ok(())
    }
}

fn saturation_pressure(
    definition: &SqueezeFilm<'_>,
    liquid: &ResolvedLiquidState,
) -> Result<f64, Error> {
    let property = liquid.material().property(definition.saturation_pressure_property)
        .filter(|property| property.requirement().quantity() == QuantitySpec::dimensional(Pressure::DIMS))
        .ok_or(Error::Invalid("squeeze film requires its original pressure-valued saturation claim in the liquid bundle"))?;
    let pressure = property.value_si();
    if !nonnegative(pressure) {
        return Err(Error::Invalid(
            "squeeze saturation pressure must be finite and nonnegative",
        ));
    }
    Ok(pressure)
}

fn temperature_support(
    law: &LiquidTemperatureLaw<'_>,
    from: &ResolvedLiquidState,
    to: &ResolvedLiquidState,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), Error> {
    let a = from.material().query();
    let b = to.material().query();
    let axis = law.temperature_axis();
    let (low, high) = if a.axes()[axis] <= b.axes()[axis] {
        (&a, &b)
    } else {
        (&b, &a)
    };
    source_line(law, axis, low, high, checkpoint)
}

fn pressure_support(
    definition: &SqueezeFilm<'_>,
    liquid: &ResolvedLiquidState,
    ambient: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), Error> {
    let excursion = definition.limits.maximum_relative_gauge_pressure * ambient;
    let point = liquid.material().query();
    let lower = point
        .clone()
        .with(definition.pressure_axis, ambient - excursion)
        .map_err(|_| Error::Invalid("squeeze lower pressure support is unrepresentable"))?;
    let upper = point
        .with(definition.pressure_axis, finite(ambient + excursion)?)
        .map_err(|_| Error::Invalid("squeeze upper pressure support is unrepresentable"))?;
    source_line(
        &definition.liquid,
        definition.pressure_axis,
        &lower,
        &upper,
        checkpoint,
    )
}

fn source_line(
    law: &LiquidTemperatureLaw<'_>,
    axis: &str,
    lower: &QueryPoint,
    upper: &QueryPoint,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
) -> Result<(), Error> {
    for property in law.reference_state().material().properties() {
        checkpoint()?;
        let selected = property.answer().receipt.selected;
        let claim = law
            .card()
            .claims()
            .claim(selected)
            .ok_or(Error::Invalid("squeeze selected liquid claim is absent"))?;
        let line = law
            .card()
            .claims()
            .query_line_typed(
                &claim.key,
                axis,
                lower,
                upper,
                ClaimSelection::Pinned(selected),
            )
            .map_err(|source| MaterialStatePointError::Support {
                property: property.requirement().name().into(),
                source: Box::new(source),
            })?;
        // Reuse the original complete bundle's numerical admission at every
        // source breakpoint; signed optional roles keep their own requirements.
        for (coordinate, _) in line.knots {
            checkpoint()?;
            let point = lower
                .clone()
                .with(axis, coordinate)
                .map_err(|_| Error::Invalid("squeeze source support coordinate is invalid"))?;
            law.reference_state()
                .material()
                .at_query_point(law.card(), &point)?;
        }
    }
    Ok(())
}
