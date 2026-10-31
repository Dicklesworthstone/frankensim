//! Original T/shear source lines drive the existing journal torque/heat path.
//!
//! The native finite-gap owner compiles the exact piecewise-affine viscosity
//! line at each accepted thermal endpoint. Mechanical RK stages use that
//! immutable line; they do not re-query the material database or approximate
//! the cylindrical profile with one selected apparent viscosity.

use super::{Error, JournalBearing, JournalGeometry, ThermalJournalBinding, positive};
use fs_matdb::ClaimSelection;
use fs_material::liquid::{ResolvedLiquidState, DYNAMIC_VISCOSITY_PROPERTY, KINEMATIC_VISCOSITY_PROPERTY};
use fs_material::liquid::surface::ResolvedLiquidSurfaceState;
use fs_material::liquid_temperature::{LiquidTemperatureLaw, LiquidViscosityBasis};
use fs_material::state_point::MaterialStatePointError;
use fs_thermal::enthalpy::DepositionError;
use fs_tribo::journal::{GeneralizedNewtonianJournal, JournalShearResponse,
    JournalViscosityKnot, MAX_JOURNAL_VISCOSITY_KNOTS};
pub use fs_tribo::journal::JournalShearSolve;

/// The existing journal with a complete monotone generalized-Newtonian source.
#[derive(Clone, Debug)]
pub struct RheologicalJournalBearing<'a> {
    /// Its original liquid reference point must declare zero physical shear.
    pub bearing: JournalBearing<'a>,
    /// Explicit dimensional 1/s coordinate, not cyclic/angular frequency.
    pub shear_rate_axis: &'a str,
    /// Original continuous source support is required over [0, maximum].
    pub maximum_shear_rate_per_s: f64,
    pub solve: JournalShearSolve,
    /// Optional SAME source-owned cylindrical expansion as the mechanism.
    pub expansion: Option<ThermalJournalBinding>,
}

/// Accepted wall states retain the complete original source bundle separately.
#[derive(Clone, Debug, PartialEq)]
pub struct JournalRheologyState {
    pub inner_liquid: ResolvedLiquidState,
    pub outer_liquid: ResolvedLiquidState,
    pub shear_rate_per_s: [f64; 2],
    pub angular_rate_defect_per_s: f64,
    pub iterations: usize,
}

#[derive(Clone, Copy)]
pub(super) struct RheologyBinding<'a> {
    pub rate_axis: &'a str,
    pub maximum_rate: f64,
    pub solve: JournalShearSolve,
}

impl RheologyBinding<'_> {
    pub fn admit(&self, law: &LiquidTemperatureLaw<'_>) -> Result<(), Error> {
        if !positive(self.maximum_rate)
            || law.reference_point().axes().get(self.rate_axis) != Some(&0.0)
        {
            return Err(Error::Invalid("journal rheology needs a zero-shear reference and positive rate support"));
        }
        ResolvedLiquidSurfaceState::try_for_shear_rate_evolution(
            law.reference_state().material().clone(), law.viscosity_basis(),
            law.temperature_axis(), self.rate_axis,
        )?;
        if ["density", "specific_heat_capacity"].iter().any(|name|
            law.reference_state().material().property(name).is_some_and(|property|
                property.source_axes().iter().any(|(axis, _)| axis == self.rate_axis)))
        {
            return Err(Error::Invalid("journal rheology cannot own shear-dependent density or caloric storage"));
        }
        Ok(())
    }
}

#[derive(Clone)]
pub(super) struct BoundRheology<'a> {
    pub binding: RheologyBinding<'a>,
    pub native: GeneralizedNewtonianJournal,
}

impl<'a> BoundRheology<'a> {
    pub fn compile(binding: RheologyBinding<'a>, law: &LiquidTemperatureLaw<'_>,
        liquid: &ResolvedLiquidState, geometry: JournalGeometry,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<Self, Error> {
        checkpoint()?;
        let lower = liquid.material().query();
        let upper = lower.clone().with(binding.rate_axis, binding.maximum_rate)
            .map_err(|_| Error::Invalid("invalid journal maximum shear coordinate"))?;
        let viscosity_name = match law.viscosity_basis() {
            LiquidViscosityBasis::Dynamic => DYNAMIC_VISCOSITY_PROPERTY,
            LiquidViscosityBasis::Kinematic => KINEMATIC_VISCOSITY_PROPERTY,
        };
        let mut viscosity = None;
        for property in liquid.material().properties() {
            checkpoint()?;
            let selected = property.answer().receipt.selected;
            let claim = law.card().claims().claim(selected)
                .ok_or(Error::Invalid("journal selected rheology claim is missing"))?;
            let line = law.card().claims().query_line_typed(&claim.key,
                binding.rate_axis, &lower, &upper, ClaimSelection::Pinned(selected))
                .map_err(|source| MaterialStatePointError::Support {
                    property: property.requirement().name().to_owned(), source: Box::new(source),
                })?;
            if line.knots.len() > MAX_JOURNAL_VISCOSITY_KNOTS {
                return Err(Error::Budget("journal source line exceeds its native breakpoint budget"));
            }
            let mut knots = Vec::with_capacity(line.knots.len());
            for (rate, _) in line.knots {
                checkpoint()?;
                let point = lower.clone().with(binding.rate_axis, rate)
                    .map_err(|_| Error::Invalid("invalid journal source shear coordinate"))?;
                // Every retained role keeps its own quantity and admissibility.
                let sampled = liquid.at_query_point(law.card(), &point)?;
                if property.requirement().name() == viscosity_name {
                    knots.push(JournalViscosityKnot { shear_rate_per_s: rate,
                        dynamic_viscosity_pa_s: sampled.dynamic_viscosity_pa_s() });
                }
            }
            if property.requirement().name() == viscosity_name { viscosity = Some(knots); }
        }
        let native = GeneralizedNewtonianJournal::new(geometry,
            viscosity.ok_or(Error::Invalid("journal original viscosity line is missing"))?,
            binding.solve)?;
        checkpoint()?;
        Ok(Self { binding, native })
    }

    pub fn at_temperature(&self, law: &LiquidTemperatureLaw<'_>,
        from: &ResolvedLiquidState, to: &ResolvedLiquidState, geometry: JournalGeometry,
        previous_wall_rates: [f64; 2], maximum_relative_change: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<Self, Error> {
        // Source lines follow the declared thermal split at held shear. The
        // next rate sweep then uses the full new-temperature line. This does
        // not pretend corner samples establish arbitrary 2-D box coverage.
        for rate in self.native.knots().iter().map(|knot| knot.shear_rate_per_s)
            .chain(previous_wall_rates)
        {
            checkpoint()?;
            temperature_path(law, from, to, self.binding.rate_axis, rate, checkpoint)?;
        }
        let next = Self::compile(self.binding, law, to, geometry, checkpoint)?;
        // A ratio of positive affine laws is monotone within each shared
        // segment. The union of original breakpoints therefore gives its
        // exact maximum, including temperature-induced breakpoint movement.
        for point in self.native.knots().iter().chain(next.native.knots()) {
            checkpoint()?;
            let old = self.native.viscosity_at_shear_rate(point.shear_rate_per_s)?;
            let new = next.native.viscosity_at_shear_rate(point.shear_rate_per_s)?;
            let change = (new - old).abs() / old;
            if !change.is_finite() || change > maximum_relative_change {
                return Err(Error::Budget("journal rheology change requires a smaller thermal substep"));
            }
        }
        Ok(next)
    }

    pub fn sample(&self, law: &LiquidTemperatureLaw<'_>, liquid: &ResolvedLiquidState,
        response: JournalShearResponse) -> Result<JournalRheologyState, Error> {
        let at = |rate| {
            let point = liquid.material().query().with(self.binding.rate_axis, rate)
                .map_err(|_| Error::Invalid("invalid journal endpoint shear coordinate"))?;
            Ok::<_, Error>(liquid.at_query_point(law.card(), &point)?)
        };
        Ok(JournalRheologyState {
            inner_liquid: at(response.shear_rate_per_s[0])?,
            outer_liquid: at(response.shear_rate_per_s[1])?,
            shear_rate_per_s: response.shear_rate_per_s,
            angular_rate_defect_per_s: response.angular_rate_defect_per_s,
            iterations: response.iterations,
        })
    }
}

fn temperature_path(law: &LiquidTemperatureLaw<'_>, from: &ResolvedLiquidState,
    to: &ResolvedLiquidState, rate_axis: &str, rate: f64,
    checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>) -> Result<(), Error> {
    let axis = law.temperature_axis();
    let a = from.material().query().with(rate_axis, rate)
        .map_err(|_| Error::Invalid("invalid journal thermal source coordinate"))?;
    let b = to.material().query().with(rate_axis, rate)
        .map_err(|_| Error::Invalid("invalid journal thermal source coordinate"))?;
    let (lo, hi) = if a.axes()[axis] <= b.axes()[axis] { (&a, &b) } else { (&b, &a) };
    for property in from.material().properties() {
        checkpoint()?;
        let selected = property.answer().receipt.selected;
        let claim = law.card().claims().claim(selected)
            .ok_or(Error::Invalid("journal thermal rheology claim is missing"))?;
        let line = law.card().claims().query_line_typed(&claim.key, axis, lo, hi,
            ClaimSelection::Pinned(selected)).map_err(|source| MaterialStatePointError::Support {
                property: property.requirement().name().into(), source: Box::new(source),
            })?;
        if line.knots.len() > MAX_JOURNAL_VISCOSITY_KNOTS {
            return Err(Error::Budget("journal thermal source path exceeds its breakpoint budget"));
        }
        for (temperature, _) in line.knots {
            checkpoint()?;
            let point = lo.clone().with(axis, temperature)
                .map_err(|_| Error::Invalid("invalid journal thermal breakpoint"))?;
            from.at_query_point(law.card(), &point)?;
        }
    }
    Ok(())
}
