//! Pinned source contact resistance in the original implicit thermal graph.
//!
//! Scalar constants and continuously supported linear source curves are compiled
//! through the canonical source-line owner. Their response coordinate may be
//! either side's temperature/fraction, the declared mean temperature, or a fixed
//! external coordinate. Actual two-sided applicability is checked at every trial.
//! A sampled surface needs a canonical source-cell gradient owner and is refused;
//! this module does not infer gradients from endpoint values or source differences.
//!
//! The original caloric owners provide numeric states and exact local tangents.
//! Successful runtime queries allocate nothing and perform no source query,
//! identity hash, source selection, thermal solve or physical-state mutation.

use super::{ContactInterfaceState, MaterialContactLink, resolve_contact};
use crate::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    ConductionError,
};
use fs_matdb::{
    ClaimSelection, InterpolationPolicy, PropertyKey, PropertyLineSlice, PropertySupportError,
    PropertyValue, QueryPoint,
};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, PhaseStateError};
use fs_material::state_point::MaterialPropertySelection;
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_material::thermal::source::ThermalMaterialError;
use fs_thermal::enthalpy::conduction::implicit::contact::{
    ImplicitContactLaw, ImplicitContactState,
};
use fs_thermal::enthalpy::{DepositionError, DepositionSimulation, ThermalLaw};

/// Original preparation failures and allocation-free runtime domain diagnostics.
#[derive(Debug)]
pub enum PreparedContactError {
    Invalid(&'static str),
    Contact(ConductionError),
    Material(ThermalMaterialError),
    Phase(PhaseStateError),
    Source(PropertySupportError),
    Thermal(DepositionError),
    OutsideEnthalpyDomain {
        side: usize,
        enthalpy_j_kg: f64,
        bounds_j_kg: [f64; 2],
    },
    /// Index into `source_bounds`, whose original axis name remains available.
    OutsideSourceDomain {
        axis: usize,
        value: f64,
        bounds: [f64; 2],
    },
    Unrepresentable,
}
impl core::fmt::Display for PreparedContactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "implicit material contact: {self:?}")
    }
}
impl std::error::Error for PreparedContactError {}
impl From<ConductionError> for PreparedContactError {
    fn from(error: ConductionError) -> Self {
        Self::Contact(error)
    }
}
impl From<ThermalMaterialError> for PreparedContactError {
    fn from(error: ThermalMaterialError) -> Self {
        Self::Material(error)
    }
}
impl From<PhaseStateError> for PreparedContactError {
    fn from(error: PhaseStateError) -> Self {
        Self::Phase(error)
    }
}
impl From<PropertySupportError> for PreparedContactError {
    fn from(error: PropertySupportError) -> Self {
        Self::Source(error)
    }
}
impl From<DepositionError> for PreparedContactError {
    fn from(error: DepositionError) -> Self {
        Self::Thermal(error)
    }
}

#[derive(Clone, Copy, Debug)]
enum Coordinate {
    Temperature(usize),
    Fraction(usize),
    MeanTemperature,
    Fixed(f64),
}
impl Coordinate {
    fn from_axis(
        point: ThermalInterfacePoint<'_>,
        query: &QueryPoint,
        name: &str,
    ) -> Result<Self, PreparedContactError> {
        if let Some(side) = point.temperature_axes.iter().position(|&axis| axis == name) {
            Ok(Self::Temperature(side))
        } else if let Some(side) = point
            .liquid_fraction_axes
            .iter()
            .position(|&axis| axis == Some(name))
        {
            Ok(Self::Fraction(side))
        } else if point.mean_temperature_axis == Some(name) {
            Ok(Self::MeanTemperature)
        } else {
            query
                .axes()
                .get(name)
                .copied()
                .map(Self::Fixed)
                .ok_or(PreparedContactError::Invalid(
                    "source applicability axis is absent",
                ))
        }
    }
    fn value(self, states: [EnthalpyPhaseKnot; 2]) -> f64 {
        match self {
            Self::Temperature(side) => states[side].temperature_k,
            Self::Fraction(side) => states[side].liquid_mass_fraction,
            Self::MeanTemperature => 0.5 * states[0].temperature_k + 0.5 * states[1].temperature_k,
            Self::Fixed(value) => value,
        }
    }
    fn derivatives(
        self,
        curves: [&EquilibriumEnthalpyPhaseCurve; 2],
        at: [f64; 2],
    ) -> Result<[f64; 2], PhaseStateError> {
        let mut result = [0.0; 2];
        match self {
            Self::Temperature(side) => {
                result[side] = curves[side].temperature_derivative_at_specific_enthalpy(at[side])?
            }
            Self::Fraction(side) => {
                result[side] =
                    curves[side].liquid_mass_fraction_derivative_at_specific_enthalpy(at[side])?
            }
            Self::MeanTemperature => {
                for side in 0..2 {
                    result[side] =
                        0.5 * curves[side].temperature_derivative_at_specific_enthalpy(at[side])?;
                }
            }
            Self::Fixed(_) => {}
        }
        Ok(result)
    }
}

#[derive(Clone, Debug)]
struct Constraint {
    axis: String,
    coordinate: Coordinate,
    bounds: [f64; 2],
}

/// An original ordered interface claim on explicit, borrowed caloric domains.
/// The complete source line and initial physical source result retain authority.
/// Numerical trial samples are compiled responses, not fresh query receipts.
#[derive(Clone, Debug)]
pub struct PreparedMaterialContact<'a> {
    curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    bounds: [[f64; 2]; 2],
    coordinate: Coordinate,
    source: PropertyLineSlice,
    slopes: Vec<f64>,
    constraints: Vec<Constraint>,
    reference: ContactInterfaceState,
}

impl<'a> PreparedMaterialContact<'a> {
    /// Pin the source chosen by the existing material-contact owner at the
    /// actual initial state. Compile its full scalar response interval, while
    /// retaining every other applicability coordinate for actual trial checks.
    ///
    /// The canonical line's knot/work bound remains authoritative. Original
    /// source knots are retained; no empirical mixing rule or numerical source
    /// derivative is inferred. Unsupported nonlinear/multiaxis payloads refuse.
    /// Evolving phase requires its explicitly declared fraction coordinate.
    pub fn new(
        thermal: &DepositionSimulation<'_>,
        binding: MaterialContactLink<'a>,
        enthalpy_bounds_j_kg: [[f64; 2]; 2],
        mut checkpoint: impl FnMut() -> Result<(), DepositionError>,
    ) -> Result<Self, PreparedContactError> {
        checkpoint()?;
        let reference = resolve_contact(thermal, &binding)?;
        let query = binding.point.at(reference.phases)?;
        let curves = binding.sides.map(|side| side.curve);
        let mut lower = [
            curves[0].values_at_specific_enthalpy(enthalpy_bounds_j_kg[0][0])?,
            curves[1].values_at_specific_enthalpy(enthalpy_bounds_j_kg[1][0])?,
        ];
        let mut upper = lower;
        for side in 0..2 {
            checkpoint()?;
            let [lo, hi] = enthalpy_bounds_j_kg[side];
            if !lo.is_finite() || !hi.is_finite() || lo >= hi {
                return Err(PreparedContactError::Invalid(
                    "finite increasing contact enthalpy bounds required",
                ));
            }
            lower[side] = curves[side].values_at_specific_enthalpy(lo)?;
            upper[side] = curves[side].values_at_specific_enthalpy(hi)?;
            if lower[side].liquid_mass_fraction != upper[side].liquid_mass_fraction
                && binding.point.liquid_fraction_axes[side].is_none()
            {
                return Err(PreparedContactError::Invalid(
                    "phase-changing contact needs its evolving fraction axis",
                ));
            }
        }
        let selected = reference.interface.receipt().selected;
        // Reuse the common static-response admission, including frequency and
        // explicitly owned/fixed response coordinates, before compiling a line.
        binding.point.resolve_resistance(
            reference.phases,
            MaterialPropertySelection::PinnedByProperty(vec![(RP.to_owned(), selected)]),
        )?;
        let claim = binding.point.interface.claims().claim(selected).ok_or(
            PreparedContactError::Invalid("selected contact source is absent"),
        )?;
        let axis = match (&claim.value, claim.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {
                binding.point.temperature_axes[0]
            }
            (PropertyValue::Curve { abscissa, .. }, InterpolationPolicy::LinearInside) => {
                abscissa.as_str()
            }
            _ => {
                return Err(PreparedContactError::Invalid(
                    "implicit contact requires a constant scalar or continuously supported linear source curve",
                ));
            }
        };
        let coordinate = Coordinate::from_axis(binding.point, &query, axis)?;
        let endpoints = [coordinate.value(lower), coordinate.value(upper)];
        let line_lower = query
            .clone()
            .with(axis, endpoints[0])
            .map_err(|_| PreparedContactError::Invalid("invalid contact source coordinate"))?;
        let line_upper = query
            .clone()
            .with(axis, endpoints[1])
            .map_err(|_| PreparedContactError::Invalid("invalid contact source coordinate"))?;
        let source = binding.point.interface.claims().query_line_typed(
            &PropertyKey::new(RP, RD),
            axis,
            &line_lower,
            &line_upper,
            ClaimSelection::Pinned(selected),
        )?;
        checkpoint()?;
        let mut slopes = Vec::with_capacity(source.knots.len().saturating_sub(1));
        for (index, (_, answer)) in source.knots.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint()?;
            }
            let value = answer.evidence.value.value;
            if !value.is_finite() || value <= 0.0 {
                return Err(PreparedContactError::Invalid(
                    "contact resistance must stay positive on every source interval",
                ));
            }
            if index > 0 {
                let (x0, old) = &source.knots[index - 1];
                let difference = value - old.evidence.value.value;
                let width = source.knots[index].0 - x0;
                let slope = difference / width;
                if !width.is_finite()
                    || width <= 0.0
                    || !slope.is_finite()
                    || difference != 0.0 && slope == 0.0
                {
                    return Err(PreparedContactError::Unrepresentable);
                }
                slopes.push(slope);
            }
        }
        let mut constraints = Vec::with_capacity(claim.validity.bounds().len());
        for (axis, &(lo, hi)) in claim.validity.bounds() {
            checkpoint()?;
            constraints.push(Constraint {
                axis: axis.clone(),
                coordinate: Coordinate::from_axis(binding.point, &query, axis)?,
                bounds: [lo, hi],
            });
        }
        let prepared = Self {
            curves,
            bounds: enthalpy_bounds_j_kg,
            coordinate,
            source,
            slopes,
            constraints,
            reference,
        };
        prepared.state_at_specific_enthalpies(
            prepared
                .reference
                .phases
                .map(|s| s.specific_enthalpy_j_kg()),
        )?;
        checkpoint()?;
        Ok(prepared)
    }

    pub const fn caloric_curves(&self) -> [&'a EquilibriumEnthalpyPhaseCurve; 2] {
        self.curves
    }
    pub const fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        self.bounds
    }
    pub fn source_line(&self) -> &PropertyLineSlice {
        &self.source
    }
    /// Original physical reference query, including its exact interface pin.
    pub fn reference(&self) -> &ContactInterfaceState {
        &self.reference
    }
    /// Original applicability bounds in the order used by runtime errors.
    pub fn source_bounds(&self) -> impl ExactSizeIterator<Item = (&str, [f64; 2])> {
        self.constraints
            .iter()
            .map(|bound| (bound.axis.as_str(), bound.bounds))
    }

    /// Exact source response and analytic chain-rule derivatives at these h values.
    /// Source corners select right slope; the compiled upper coordinate selects
    /// left. Caloric tangents likewise select left at each admitted upper h.
    pub fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<ImplicitContactState, PreparedContactError> {
        for side in 0..2 {
            if !h[side].is_finite()
                || !(self.bounds[side][0]..=self.bounds[side][1]).contains(&h[side])
            {
                return Err(PreparedContactError::OutsideEnthalpyDomain {
                    side,
                    enthalpy_j_kg: h[side],
                    bounds_j_kg: self.bounds[side],
                });
            }
        }
        let values = [
            self.curves[0].values_at_specific_enthalpy(h[0])?,
            self.curves[1].values_at_specific_enthalpy(h[1])?,
        ];
        for (axis, constraint) in self.constraints.iter().enumerate() {
            let value = constraint.coordinate.value(values);
            if !value.is_finite() || !(constraint.bounds[0]..=constraint.bounds[1]).contains(&value)
            {
                return Err(PreparedContactError::OutsideSourceDomain {
                    axis,
                    value,
                    bounds: constraint.bounds,
                });
            }
        }
        let x = self.coordinate.value(values);
        let knots = &self.source.knots;
        let first = knots
            .first()
            .ok_or(PreparedContactError::Invalid("empty compiled source line"))?;
        let last = &knots[knots.len() - 1];
        if !x.is_finite() || !(first.0..=last.0).contains(&x) {
            return Err(PreparedContactError::Invalid(
                "contact left its compiled source interval",
            ));
        }
        let (resistance, slope) = if knots.len() == 1 {
            (first.1.evidence.value.value, 0.0)
        } else {
            let index = knots.partition_point(|k| k.0 <= x).min(knots.len() - 1) - 1;
            let (a, ka) = &knots[index];
            let (b, kb) = &knots[index + 1];
            let ra = ka.evidence.value.value;
            let rb = kb.evidence.value.value;
            let fraction = (x - a) / (b - a);
            let value = if x == *a {
                ra
            } else if x == *b {
                rb
            } else if fraction <= 0.5 {
                (rb - ra).mul_add(fraction, ra)
            } else {
                (ra - rb).mul_add(1.0 - fraction, rb)
            };
            (value, self.slopes[index])
        };
        let mut at = h;
        for side in 0..2 {
            at[side] = if h[side] == self.bounds[side][1] {
                h[side].next_down()
            } else {
                h[side]
            };
        }
        let coordinate_derivatives = if slope == 0.0 {
            [0.0; 2]
        } else {
            self.coordinate.derivatives(self.curves, at)?
        };
        let derivatives = coordinate_derivatives.map(|value| value * slope);
        if !resistance.is_finite()
            || resistance <= 0.0
            || derivatives.iter().any(|slope| !slope.is_finite())
            || derivatives
                .iter()
                .zip(coordinate_derivatives)
                .any(|(&value, coordinate)| value == 0.0 && coordinate != 0.0 && slope != 0.0)
        {
            return Err(PreparedContactError::Unrepresentable);
        }
        Ok(ImplicitContactState {
            resistance_m2_k_w: resistance,
            resistance_enthalpy_derivatives_m2_kg_k_w_j: derivatives,
        })
    }
}

impl ImplicitContactLaw for PreparedMaterialContact<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        [self.curves[0], self.curves[1]]
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        self.bounds
    }
    fn state_at_specific_enthalpies(
        &self,
        h: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        self.state_at_specific_enthalpies(h)
            .map_err(|_| DepositionError::MaterialDomain)
    }
}
