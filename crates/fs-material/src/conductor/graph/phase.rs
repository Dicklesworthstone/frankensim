//! Two-coordinate derivatives of the original pinned Ohmic source.
use super::*;
use fs_matdb::{ClaimSelection, InterpolationPolicy, PreparedSurfaceResponse, PropertyValue};
use fs_qty::QuantitySpec;

#[derive(Clone, Debug)]
enum Gradient {
    Constant,
    Curve {
        coordinate: usize,
        knots: Vec<(f64, f64)>,
    },
    Surface {
        response: PreparedSurfaceResponse,
        coordinates: [usize; 2],
    },
}

#[derive(Clone, Debug)]
pub(super) struct PreparedOhmicPhase {
    pub axis: String,
    pub reference: f64,
    pub range: [f64; 2],
    gradient: Gradient,
}

impl PreparedOhmicPhase {
    pub fn new(
        response: &ScalarConductorResponse,
        temperature_axis: &str,
        temperature_range: [f64; 2],
        range: [f64; 2],
    ) -> Result<Self, ConductorError> {
        let original = response.original();
        let axis = original.liquid_mass_fraction_axis().ok_or_else(|| {
            invalid("phase Ohmic node requires the original explicit fraction axis")
        })?;
        let query = original.material().query();
        let reference = *query
            .axes()
            .get(axis)
            .ok_or_else(|| invalid("phase Ohmic source fraction is absent"))?;
        if axis == temperature_axis
            || !range[0].is_finite()
            || !range[1].is_finite()
            || range[0] < 0.0
            || range[1] > 1.0
            || range[0] >= range[1]
            || !(range[0]..=range[1]).contains(&reference)
        {
            return Err(invalid(
                "phase Ohmic interval must contain its original mass fraction in [0,1]",
            ));
        }
        let owned = [(axis, QuantitySpec::dimensional(Dims::NONE))];
        if original.material().properties().iter().any(|property| {
            !property.has_only_owned_coordinates_and_fixed_external_conditions(&owned)
        }) {
            return Err(invalid(
                "phase Ohmic source has an unowned response coordinate",
            ));
        }
        let property = original.source_property();
        let selected = property.answer().receipt.selected;
        let claim = response
            .card()
            .claims()
            .claim(selected)
            .ok_or_else(|| invalid("original phase Ohmic source is absent"))?;
        let coordinate = |name: &str| -> Result<usize, ConductorError> {
            if name == temperature_axis {
                Ok(0)
            } else if name == axis {
                Ok(1)
            } else {
                Err(invalid(
                    "phase Ohmic response must use its original T/f coordinates",
                ))
            }
        };
        let gradient = match (&claim.value, claim.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {
                Gradient::Constant
            }
            (PropertyValue::Curve { abscissa, .. }, InterpolationPolicy::LinearInside) => {
                let coordinate = coordinate(abscissa)?;
                let limits = if coordinate == 0 {
                    temperature_range
                } else {
                    range
                };
                let lower = query
                    .clone()
                    .with(abscissa, limits[0])
                    .map_err(|_| invalid("invalid phase Ohmic lower source coordinate"))?;
                let upper = query
                    .clone()
                    .with(abscissa, limits[1])
                    .map_err(|_| invalid("invalid phase Ohmic upper source coordinate"))?;
                let line = response
                    .card()
                    .claims()
                    .query_line_typed(
                        &claim.key,
                        abscissa,
                        &lower,
                        &upper,
                        ClaimSelection::Pinned(selected),
                    )
                    .map_err(|source| MaterialStatePointError::Support {
                        property: claim.key.name().to_owned(),
                        source: Box::new(source),
                    })?;
                let knots: Vec<_> = line
                    .knots
                    .into_iter()
                    .map(|(x, answer)| (x, answer.evidence.value.value))
                    .collect();
                if knots.len() < 2
                    || knots.first().map(|p| p.0) != Some(limits[0])
                    || knots.last().map(|p| p.0) != Some(limits[1])
                    || knots.windows(2).any(|p| {
                        let width = p[1].0 - p[0].0;
                        let difference = p[1].1 - p[0].1;
                        let slope = difference / width;
                        !width.is_finite()
                            || width <= 0.0
                            || !slope.is_finite()
                            || difference != 0.0 && slope == 0.0
                    })
                {
                    return Err(invalid(
                        "phase Ohmic source curve or derivative is unrepresentable",
                    ));
                }
                Gradient::Curve { coordinate, knots }
            }
            (PropertyValue::Surface2D { .. }, InterpolationPolicy::LinearInside) => {
                let response = claim.value.prepare_surface_response().map_err(|_| {
                    invalid("phase Ohmic canonical source surface cannot be prepared")
                })?;
                let coordinates = [
                    coordinate(&response.axes()[0])?,
                    coordinate(&response.axes()[1])?,
                ];
                if coordinates[0] == coordinates[1] {
                    return Err(invalid(
                        "phase Ohmic source surface requires distinct T/f axes",
                    ));
                }
                Gradient::Surface {
                    response,
                    coordinates,
                }
            }
            _ => {
                return Err(invalid(
                    "phase Ohmic source needs constant or continuously interpolated response",
                ));
            }
        };
        Ok(Self {
            axis: axis.into(),
            reference,
            range,
            gradient,
        })
    }

    pub fn admit_fraction(&self, fraction: f64) -> Result<(), ConductorError> {
        if !fraction.is_finite() || !(self.range[0]..=self.range[1]).contains(&fraction) {
            return Err(invalid(
                "Ohmic fraction left its declared material interval",
            ));
        }
        Ok(())
    }

    pub fn response_at(
        &self,
        conductor: &ResolvedConductor,
        temperature: f64,
        fraction: f64,
        current: f64,
    ) -> Result<([f64; 2], [f64; 6]), GraphError> {
        let point = [temperature, fraction];
        let mut source_slopes = [0.0; 2];
        match &self.gradient {
            Gradient::Constant => {}
            Gradient::Curve { coordinate, knots } => {
                let index = knots
                    .partition_point(|p| p.0 <= point[*coordinate])
                    .saturating_sub(1)
                    .min(knots.len() - 2);
                let (a, va) = knots[index];
                let (b, vb) = knots[index + 1];
                source_slopes[*coordinate] = (vb - va) / (b - a);
            }
            Gradient::Surface {
                response,
                coordinates,
            } => {
                let value = response
                    .evaluate(coordinates.map(|i| point[i]))
                    .map_err(|_| {
                        graph_error("phase Ohmic source surface lacks actual support or gradient")
                    })?;
                for i in 0..2 {
                    source_slopes[coordinates[i]] = value.gradient[i];
                }
            }
        }
        let resistance = conductor.resistance_ohm();
        let sign = match conductor.source_basis() {
            ElectricalPropertyBasis::Resistivity => 1.0,
            ElectricalPropertyBasis::Conductivity => -1.0,
        };
        let derivatives = source_slopes
            .map(|slope| resistance * (slope / conductor.source_property().value_si()) * sign);
        let voltage = current * resistance;
        let heat = conductor
            .joule_power_w(current)
            .map_err(|_| graph_error("phase Ohmic Joule power is unrepresentable"))?;
        let vt = current * derivatives[0];
        let vf = current * derivatives[1];
        // Row-major [d/dT, d/dI, d/df] for [V, Qdot]. The reciprocal source
        // basis changes BOTH derivatives; interpolating resistivity would not.
        let tangent = [
            vt,
            resistance,
            vf,
            current * vt,
            2.0 * voltage,
            current * vf,
        ];
        if !voltage.is_finite()
            || current != 0.0 && voltage == 0.0
            || tangent.iter().any(|v| !v.is_finite())
            || derivatives
                .iter()
                .zip(source_slopes)
                .any(|(&d, s)| !d.is_finite() || s != 0.0 && d == 0.0)
        {
            return Err(graph_error(
                "phase Ohmic response or derivative is unrepresentable",
            ));
        }
        Ok(([voltage, heat], tangent))
    }
}
