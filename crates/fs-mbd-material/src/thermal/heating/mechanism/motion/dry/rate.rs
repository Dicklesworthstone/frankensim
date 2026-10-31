//! Original speed/temperature friction sources on the existing dry-brake path.
//!
//! Only thermal endpoints query source data. Each mechanical stage evaluates
//! the exact affine line between retained source breakpoints. The native
//! regularization samples at max(|v|, v_reg), then ramps traction below v_reg.
//! Its conjugate mechanical power enters the existing loss accounting.
//! This owns instantaneous speed dependence, not frictional memory or sticking.
use super::{
    DryJointBrake, Error, FrictionLaw, KINETIC_FRICTION_PROPERTY, ResolvedInterfaceStatePoint,
    STATIC_FRICTION_PROPERTY, law, nonnegative, positive,
};
use fs_matdb::ClaimSelection;
use fs_material::state_point::MaterialStatePointError;
use fs_qty::Velocity;
use fs_thermal::enthalpy::DepositionError;

/// Bound for the union of both original coefficient-line breakpoint sets.
pub const MAX_FRICTION_RATE_KNOTS: usize = 256;

/// Ordered dry brake with one explicitly owned physical slip-speed coordinate.
#[derive(Clone, Debug)]
pub struct RateDependentDryJointBrake<'a> {
    /// The reference interface context must contain zero speed at the named axis.
    pub brake: DryJointBrake<'a>,
    /// A distinct, typed m/s coordinate. Direction remains the sign of slip.
    pub slip_speed_axis: &'a str,
    /// Continuous original-source support is required over [0, maximum].
    pub maximum_slip_speed_m_s: f64,
}

/// Endpoint coefficients and original source receipts at the actual |v| and T.
#[derive(Clone, Debug, PartialEq)]
pub struct DryBrakeRateState {
    pub resolved: ResolvedInterfaceStatePoint,
    pub law: FrictionLaw,
    /// Explicit native regularization sample max(|v|, v_reg), distinct from
    /// the actual physical speed retained in `resolved`.
    pub traction_speed_m_s: f64,
    pub traction_law: FrictionLaw,
}

#[derive(Clone, Copy)]
pub(super) struct RateBinding<'a> {
    pub axis: &'a str,
    pub maximum_speed: f64,
}

impl RateBinding<'_> {
    pub fn admit(&self, definition: &DryJointBrake<'_>) -> Result<(), Error> {
        let point = definition.point;
        if self.axis.trim().is_empty()
            || !positive(self.maximum_speed)
            || self.maximum_speed < definition.limits.regularization_m_s
            || point.context.axes().get(self.axis) != Some(&0.0)
            || point
                .context
                .axis_quantities()
                .get(self.axis)
                .is_none_or(|quantity| quantity.dims() != Velocity::DIMS)
            || point.temperature_axes.contains(&self.axis)
            || point.liquid_fraction_axes.contains(&Some(self.axis))
            || point.mean_temperature_axis == Some(self.axis)
            || definition.pressure_axis == self.axis
        {
            return Err(Error::Invalid(
                "rate-dependent dry brake needs a distinct typed zero-speed axis and positive speed support",
            ));
        }
        Ok(())
    }
}

#[derive(Clone, Copy)]
struct Knot {
    speed: f64,
    coefficients: [f64; 2],
}

#[derive(Clone)]
pub(super) struct BoundRate<'a> {
    binding: RateBinding<'a>,
    knots: Vec<Knot>,
}

impl<'a> BoundRate<'a> {
    pub fn compile(
        binding: RateBinding<'a>,
        definition: &DryJointBrake<'_>,
        reference: &ResolvedInterfaceStatePoint,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        checkpoint()?;
        let lower = reference.query();
        let upper = lower
            .clone()
            .with(binding.axis, binding.maximum_speed)
            .map_err(|_| Error::Invalid("invalid dry brake speed support endpoint"))?;
        let claims = definition.point.interface.claims();
        let mut coordinates = Vec::new();
        for property in reference.properties() {
            checkpoint()?;
            let selected = property.answer().receipt.selected;
            let claim = claims
                .claim(selected)
                .ok_or(Error::Invalid("missing original friction source"))?;
            let line = claims
                .query_line_typed(
                    &claim.key,
                    binding.axis,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|source| MaterialStatePointError::Support {
                    property: claim.key.name().to_owned(),
                    source: Box::new(source),
                })?;
            if line.knots.len() > MAX_FRICTION_RATE_KNOTS {
                return Err(Error::Budget(
                    "dry friction source exceeds its speed breakpoint budget",
                ));
            }
            coordinates.extend(line.knots.into_iter().map(|(speed, _)| speed));
        }
        coordinates.sort_by(f64::total_cmp);
        coordinates.dedup();
        if coordinates.len() > MAX_FRICTION_RATE_KNOTS
            || coordinates.first() != Some(&0.0)
            || coordinates.last() != Some(&binding.maximum_speed)
        {
            return Err(Error::Budget(
                "dry friction needs a complete bounded speed source line",
            ));
        }
        let mut knots = Vec::with_capacity(coordinates.len());
        for speed in coordinates {
            checkpoint()?;
            let point = lower
                .clone()
                .with(binding.axis, speed)
                .map_err(|_| Error::Invalid("invalid dry friction speed breakpoint"))?;
            let resolved = reference.at_query_point(definition.point.interface, &point)?;
            // Both sources are affine between their union of breakpoints.
            // Their nonnegative ordering therefore holds on whole segments.
            law(&resolved)?;
            knots.push(Knot {
                speed,
                coefficients: coefficients(&resolved),
            });
        }
        checkpoint()?;
        Ok(Self { binding, knots })
    }

    fn coefficients_at_speed(&self, speed: f64) -> Result<[f64; 2], Error> {
        if !nonnegative(speed) || speed > self.binding.maximum_speed {
            return Err(Error::Invalid(
                "dry brake slip speed exceeds original source support",
            ));
        }
        let upper = self.knots.partition_point(|knot| knot.speed < speed);
        let next = self.knots[upper];
        if upper == 0 || next.speed == speed {
            return Ok(next.coefficients);
        }
        let previous = self.knots[upper - 1];
        let t = (speed - previous.speed) / (next.speed - previous.speed);
        let result = core::array::from_fn(|i| {
            (1.0 - t) * previous.coefficients[i] + t * next.coefficients[i]
        });
        if result.into_iter().any(|value| !nonnegative(value)) {
            return Err(Error::Invalid(
                "dry friction coefficient interpolation overflow",
            ));
        }
        Ok(result)
    }

    /// Allocation-free source interpolation; no material queries in RK stages.
    pub fn law_at_speed(&self, speed: f64) -> Result<FrictionLaw, Error> {
        let [static_mu, kinetic_mu] = self.coefficients_at_speed(speed)?;
        let law = FrictionLaw::Coulomb {
            static_mu,
            kinetic_mu,
        };
        law.validate()?;
        Ok(law)
    }

    pub fn sample(
        &self,
        definition: &DryJointBrake<'_>,
        reference: &ResolvedInterfaceStatePoint,
        speed: f64,
    ) -> Result<DryBrakeRateState, Error> {
        self.law_at_speed(speed)?;
        let point = reference
            .query()
            .with(self.binding.axis, speed)
            .map_err(|_| Error::Invalid("invalid dry friction endpoint speed"))?;
        let resolved = reference.at_query_point(definition.point.interface, &point)?;
        let law = law(&resolved)?;
        let traction_speed_m_s = speed.max(definition.limits.regularization_m_s);
        let traction_law = self.law_at_speed(traction_speed_m_s)?;
        Ok(DryBrakeRateState {
            resolved,
            law,
            traction_speed_m_s,
            traction_law,
        })
    }

    pub fn at_temperature(
        &self,
        definition: &DryJointBrake<'_>,
        from: &ResolvedInterfaceStatePoint,
        to: &ResolvedInterfaceStatePoint,
        temperature_axis: Option<&str>,
        endpoint_speeds: [f64; 2],
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, Error> {
        let next = Self::compile(self.binding, definition, to, checkpoint)?;
        // This is the existing thermal split: coefficient updates at fixed
        // speed, then mechanical motion over the complete new source line.
        // Do not infer arbitrary two-dimensional support from corner samples.
        for speed in self
            .knots
            .iter()
            .chain(&next.knots)
            .map(|knot| knot.speed)
            .chain(endpoint_speeds)
        {
            checkpoint()?;
            self.law_at_speed(speed)?;
            next.law_at_speed(speed)?;
            if let Some(axis) = temperature_axis {
                self.temperature_path(definition, from, to, axis, speed, checkpoint)?;
            }
        }
        // A difference of affine source lines is affine on every common
        // segment. The union of knots bounds the entire thermal change.
        for knot in self.knots.iter().chain(&next.knots) {
            checkpoint()?;
            let old = self.coefficients_at_speed(knot.speed)?;
            let new = next.coefficients_at_speed(knot.speed)?;
            if old
                .into_iter()
                .zip(new)
                .any(|(old, new)| (new - old).abs() > definition.limits.maximum_coefficient_change)
            {
                return Err(Error::Budget(
                    "dry friction speed law changed too far; refine thermal substeps",
                ));
            }
        }
        Ok(next)
    }

    fn temperature_path(
        &self,
        definition: &DryJointBrake<'_>,
        from: &ResolvedInterfaceStatePoint,
        to: &ResolvedInterfaceStatePoint,
        axis: &str,
        speed: f64,
        checkpoint: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<(), Error> {
        let old = from.query();
        let a = old.axes()[axis];
        let b = to.query().axes()[axis];
        let lower = old
            .clone()
            .with(self.binding.axis, speed)
            .and_then(|point| point.with(axis, a.min(b)))
            .map_err(|_| Error::Invalid("invalid dry friction thermal path"))?;
        let upper = lower
            .clone()
            .with(axis, a.max(b))
            .map_err(|_| Error::Invalid("invalid dry friction thermal endpoint"))?;
        let claims = definition.point.interface.claims();
        let mut coordinates = Vec::new();
        for property in from.properties() {
            checkpoint()?;
            let selected = property.answer().receipt.selected;
            let claim = claims
                .claim(selected)
                .ok_or(Error::Invalid("missing original friction thermal source"))?;
            let line = claims
                .query_line_typed(
                    &claim.key,
                    axis,
                    &lower,
                    &upper,
                    ClaimSelection::Pinned(selected),
                )
                .map_err(|source| MaterialStatePointError::Support {
                    property: claim.key.name().to_owned(),
                    source: Box::new(source),
                })?;
            if line.knots.len() > MAX_FRICTION_RATE_KNOTS {
                return Err(Error::Budget(
                    "dry friction thermal path exceeds its breakpoint budget",
                ));
            }
            coordinates.extend(line.knots.into_iter().map(|(temperature, _)| temperature));
        }
        coordinates.sort_by(f64::total_cmp);
        coordinates.dedup();
        if coordinates.len() > MAX_FRICTION_RATE_KNOTS {
            return Err(Error::Budget(
                "dry friction thermal union exceeds its breakpoint budget",
            ));
        }
        for temperature in coordinates {
            checkpoint()?;
            let point = lower
                .clone()
                .with(axis, temperature)
                .map_err(|_| Error::Invalid("invalid dry friction thermal breakpoint"))?;
            law(&from.at_query_point(definition.point.interface, &point)?)?;
        }
        Ok(())
    }
}

fn coefficients(resolved: &ResolvedInterfaceStatePoint) -> [f64; 2] {
    [STATIC_FRICTION_PROPERTY, KINETIC_FRICTION_PROPERTY].map(|name| {
        resolved
            .property(name)
            .expect("resolved friction requirement")
            .value_si()
    })
}
