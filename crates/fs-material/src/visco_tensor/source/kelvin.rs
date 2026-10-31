//! Original Kelvin viscosity claims on the same tensor material and query.
//! Every symmetric Mandel entry is explicit. Thermal admission checks all
//! knots in their common refinement; the PSD cone is convex between those
//! knots because every admitted coefficient is affine there.
use super::*;
use fs_matdb::{InterpolationPolicy, PropertyValue};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{DynViscosity, Temperature};

/// Full symmetric Kelvin tensor, in Mandel order and Pa s. Upper-triangle
/// entries are C11,C12,...,C16,C22,...,C66, including every explicit zero.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct TensorKelvinViscosityProperties {
    pub viscosity_upper: [String; 21],
}

/// Same-card viscosity receipts; their physical query follows the complete
/// tensor material atomically, while no Kelvin history or energy is invented.
#[derive(Clone, Debug, PartialEq)]
pub struct ResolvedKelvinViscosity {
    properties: TensorKelvinViscosityProperties,
    state: ResolvedMaterialStatePoint,
}
impl ResolvedKelvinViscosity {
    pub const fn properties(&self) -> &TensorKelvinViscosityProperties {
        &self.properties
    }
    pub const fn state(&self) -> &ResolvedMaterialStatePoint {
        &self.state
    }

    pub(super) fn at_query_point(&self, card: &MaterialCard, point: &QueryPoint) -> Result<Self> {
        Ok(Self {
            properties: self.properties.clone(),
            state: self.state.at_query_point(card, point)?,
        })
    }

    fn apply(&self, law: TensorProny) -> Result<TensorProny> {
        let mut band = law.band_rad_s();
        let mut matrix = [[0.0; 6]; 6];
        let mut entry = 0;
        for row in 0..6 {
            for column in row..6 {
                let name = &self.properties.viscosity_upper[entry];
                let property = self
                    .state
                    .property(name)
                    .ok_or_else(|| invalid("missing original Kelvin coefficient"))?;
                if property.requirement().quantity()
                    != QuantitySpec::dimensional(DynViscosity::DIMS)
                    || !property.is_constant_at_fixed_environment()
                    || !property.value_si().is_finite()
                {
                    return Err(invalid(
                        "Kelvin coefficients require Pa s and fixed-environment, frequency-independent values",
                    ));
                }
                if let Some((low, high)) = property.angular_frequency_band_rad_s()? {
                    band = (band.0.max(low), band.1.min(high));
                }
                matrix[row][column] = property.value_si();
                matrix[column][row] = property.value_si();
                entry += 1;
            }
        }
        law.with_kelvin_viscosity_mandel(matrix)
            .and_then(|law| law.with_band(band))
            .map_err(|e| invalid(&e.to_string()))
    }
}

impl ResolvedTensorMaterial {
    /// Resolve all 21 viscosity coefficients on the original immutable card,
    /// physical query and external context. Missing coefficients never become
    /// zero. Adding this law creates no hidden strain or memory coordinates.
    pub fn with_kelvin_viscosity(
        self,
        card: &MaterialCard,
        properties: TensorKelvinViscosityProperties,
        selection: MaterialPropertySelection,
    ) -> Result<Self> {
        if self.kelvin.is_some() || card.content_hash() != self.state.card_identity() {
            return Err(invalid(
                "Kelvin binding requires the original card and no prior viscosity owner",
            ));
        }
        let mut requirements = Vec::with_capacity(21);
        for name in &properties.viscosity_upper {
            if name.len() > 128 {
                return Err(invalid("Kelvin coefficient name exceeds 128 bytes"));
            }
            if requirements.iter().any(|r: &Requirement| r.name() == name) {
                continue;
            }
            requirements.push(Requirement::try_new(
                name,
                DynViscosity::DIMS,
                Admit::Finite,
            )?);
        }
        let state =
            resolve_material_state_point(card, &self.state.query(), &requirements, selection)?;
        self.bind_kelvin(ResolvedKelvinViscosity { properties, state })
    }

    pub(super) fn bind_kelvin(mut self, source: ResolvedKelvinViscosity) -> Result<Self> {
        if self.kelvin.is_some()
            || source.state.card_identity() != self.state.card_identity()
            || source.state.query_point() != self.state.query_point()
            || source.state.axis_quantities() != self.state.axis_quantities()
        {
            return Err(invalid(
                "Kelvin and elastic sources must share the exact material, query and axis meanings",
            ));
        }
        self.law = source.apply(self.law)?;
        self.kelvin = Some(source);
        Ok(self)
    }

    /// Original Kelvin data, receipts and source selectors, when bound.
    pub fn kelvin_source(&self) -> Option<&ResolvedKelvinViscosity> {
        self.kelvin.as_ref()
    }

    /// Admit the complete viscosity interval under the original tensor owner's
    /// PSD arithmetic policy. Only validity-wide constants or piecewise-linear
    /// temperature curves are accepted. At most 4096 source knots are inspected;
    /// all union knots are checked, so unequal coefficient grids cannot hide an
    /// indefinite interior matrix. This is not a certificate for nonlinear
    /// interpolation, extrapolation, or a temperature integrator.
    pub fn admit_kelvin_thermal_interval(
        &self,
        card: &MaterialCard,
        axis: &str,
        interval_k: [f64; 2],
    ) -> Result<()> {
        let Some(source) = &self.kelvin else {
            return Ok(());
        };
        let [low, high] = interval_k;
        let original = source.state.state_coordinate(axis).ok_or_else(|| {
            invalid("Kelvin thermal admission needs the original temperature coordinate")
        })?;
        let absolute = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let quantity = source.state.axis_quantities().get(axis);
        if !low.is_finite()
            || !high.is_finite()
            || low <= 0.0
            || high <= low
            || !(low..=high).contains(&original)
            || card.content_hash() != self.state.card_identity()
            || quantity.is_some_and(|q| {
                *q != absolute && *q != QuantitySpec::dimensional(Temperature::DIMS)
            })
            || quantity.is_none() && axis != "T"
        {
            return Err(invalid(
                "Kelvin thermal interval must contain the original state on its original card",
            ));
        }
        let mut temperatures = vec![low, high];
        let mut knot_count = 0usize;
        for property in source.state.properties() {
            if !property.has_only_fixed_external_conditions(Some(axis), true) {
                return Err(invalid(
                    "Kelvin thermal feedback cannot freeze unowned source coordinates",
                ));
            }
            let claim = card
                .claims()
                .claim(property.answer().receipt.selected)
                .ok_or_else(|| invalid("original Kelvin source claim is absent"))?;
            match (&claim.value, claim.interpolation) {
                (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {}
                (
                    PropertyValue::Curve {
                        abscissa,
                        abscissa_dims,
                        knots,
                        ..
                    },
                    InterpolationPolicy::LinearInside,
                ) if abscissa == axis && *abscissa_dims == Temperature::DIMS => {
                    knot_count = knot_count
                        .checked_add(knots.len())
                        .ok_or_else(|| invalid("Kelvin thermal knot count overflow"))?;
                    if knot_count > 4096
                        || knots.len() < 2
                        || knots[0].0 > low
                        || knots[knots.len() - 1].0 < high
                        || knots.iter().any(|&(t, v)| !t.is_finite() || !v.is_finite())
                        || knots.windows(2).any(|p| p[0].0 >= p[1].0)
                    {
                        return Err(invalid(
                            "Kelvin thermal curves require a bounded, complete, strictly ordered source grid",
                        ));
                    }
                    temperatures.extend(
                        knots
                            .iter()
                            .filter_map(|&(t, _)| (low <= t && t <= high).then_some(t)),
                    );
                }
                _ => {
                    return Err(invalid(
                        "Kelvin thermal feedback requires constant or piecewise-linear temperature source laws",
                    ));
                }
            }
        }
        temperatures.sort_by(f64::total_cmp);
        temperatures.dedup();
        for temperature in temperatures {
            let point = self
                .state
                .query()
                .with(axis, temperature)
                .map_err(|e| invalid(&e.to_string()))?;
            // Same selected claims, density, all elastic tensors and the full
            // Kelvin matrix are independently admitted at every union knot.
            self.at_query_point(card, &point)?;
        }
        Ok(())
    }
}

fn invalid(reason: &str) -> TensorSourceError {
    TensorSourceError::Invalid(reason.into())
}
