//! Explicit physical coordinates for source-pinned friction queries.
use super::*;

/// Physical meaning of an additional, explicitly evolved source coordinate.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum FrictionQueryAxisKind {
    /// Absolute temperature of an original surface or declared thermal mean [K].
    Temperature,
    /// An explicitly supplied nominal pressure [Pa].
    NominalPressure,
    /// Nonnegative actual relative sliding speed [m/s].
    SlidingSpeed,
    /// Nonnegative accumulated actual sliding arc length [m].
    AccumulatedDistance,
}

impl FrictionQueryAxisKind {
    pub(super) const fn identity_tag(self) -> u8 {
        match self {
            Self::Temperature => 0,
            Self::NominalPressure => 1,
            Self::SlidingSpeed => 2,
            Self::AccumulatedDistance => 3,
        }
    }

    fn accepts(self, quantity: QuantitySpec, value: f64) -> bool {
        if !value.is_finite() || value < 0.0 {
            return false;
        }
        match self {
            Self::Temperature => {
                value > 0.0
                    && (quantity == QuantitySpec::dimensional(Temperature::DIMS)
                        || quantity.semantic_type().is_some_and(|q| {
                            q.kind() == QuantityKind::AbsoluteTemperature
                                && q.form() == ValueForm::Static
                        }))
            }
            Self::NominalPressure => {
                quantity == QuantitySpec::dimensional(fs_qty::Pressure::DIMS)
                    || quantity.semantic_type().is_some_and(|q| {
                        q.kind() == QuantityKind::Pressure && q.form() == ValueForm::Static
                    })
            }
            Self::SlidingSpeed => quantity == QuantitySpec::dimensional(Dims([1, 0, -1, 0, 0, 0])),
            Self::AccumulatedDistance => {
                quantity == QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]))
            }
        }
    }
}

/// One existing typed query axis whose value is supplied by the mechanics or
/// thermal owner. This declaration does not itself evolve speed or history.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct FrictionQueryAxis {
    pub name: String,
    pub kind: FrictionQueryAxisKind,
}

pub(super) fn validate_additional_axes(
    point: &QueryPoint,
    primary: &ThermalFrictionAxes,
    axes: &[FrictionQueryAxis],
) -> Result<(), SourceThermalCoulombError> {
    for (i, axis) in axes.iter().enumerate() {
        if axis.name.trim().is_empty()
            || axis.name == primary.temperature
            || primary.nominal_pressure.as_ref() == Some(&axis.name)
            || axes[..i].iter().any(|previous| previous.name == axis.name)
            || !point.axis_quantities().get(&axis.name).is_some_and(|q| {
                point
                    .axes()
                    .get(&axis.name)
                    .is_some_and(|&value| axis.kind.accepts(*q, value))
            })
        {
            return Err(SourceThermalCoulombError::Invalid(
                "owned friction axes require distinct existing coordinates with physical quantities and ranges",
            ));
        }
    }
    Ok(())
}

pub(super) fn validate_refresh(
    initial: &QueryPoint,
    point: &QueryPoint,
    primary: &ThermalFrictionAxes,
    axes: &[FrictionQueryAxis],
    temperature_k: Option<f64>,
    pressure_pa: Option<f64>,
) -> Result<(), SourceThermalCoulombError> {
    if initial.axis_quantities() != point.axis_quantities()
        || initial.axes().keys().ne(point.axes().keys())
    {
        return Err(SourceThermalCoulombError::Invalid(
            "friction refresh requires the original query coordinate schema",
        ));
    }
    validate_additional_axes(point, primary, axes)?;
    for (name, value) in point.axes() {
        let owned = name == &primary.temperature
            || primary.nominal_pressure.as_ref() == Some(name)
            || axes.iter().any(|axis| &axis.name == name);
        if !owned && initial.axes().get(name) != Some(value) {
            return Err(SourceThermalCoulombError::Invalid(
                "friction refresh cannot change an unowned external coordinate",
            ));
        }
    }
    if !temperature_k.is_some_and(|t| {
        t.is_finite() && t > 0.0 && point.axes().get(&primary.temperature) == Some(&t)
    }) {
        return Err(SourceThermalCoulombError::Invalid(
            "friction query temperature must equal the supplied constitutive temperature",
        ));
    }
    match (&primary.nominal_pressure, pressure_pa) {
        (Some(axis), Some(p))
            if p.is_finite() && p >= 0.0 && point.axes().get(axis) == Some(&p) => {}
        (None, None) => {}
        _ => {
            return Err(SourceThermalCoulombError::Invalid(
                "friction query pressure must equal the supplied constitutive pressure",
            ));
        }
    }
    Ok(())
}
