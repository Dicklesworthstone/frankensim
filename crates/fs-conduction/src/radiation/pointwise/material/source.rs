//! One pinned material emissivity query shared by spatial radiation owners.

use super::*;
use fs_matdb::{ClaimId, InterpolationPolicy, PreparedSurfaceResponse, PropertyValue};

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(in crate::radiation) enum DerivativeMode {
    None,
    Primal,
    Strict,
}

/// A source or phase branch has no admitted two-sided derivative.
#[derive(Debug, Clone, PartialEq)]
pub enum MaterialRadiationSmoothnessError {
    ChartEndpoint {
        vertex: usize,
        specific_enthalpy_j_kg: f64,
    },
    ChartKink {
        vertex: usize,
        specific_enthalpy_j_kg: f64,
    },
    SourceEndpoint {
        vertex: usize,
        axis: String,
        coordinate: f64,
    },
    SourceKink {
        vertex: usize,
        axis: String,
        coordinate: f64,
    },
    SourceSurface {
        vertex: usize,
        axes: [String; 2],
        coordinates: [f64; 2],
        upstream: String,
    },
    ValidityEndpoint {
        vertex: usize,
        axis: String,
        coordinate: f64,
    },
}
impl std::fmt::Display for MaterialRadiationSmoothnessError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        write!(
            f,
            "material radiation two-sided derivative refused: {self:?}"
        )
    }
}
impl std::error::Error for MaterialRadiationSmoothnessError {}

#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub(super) enum Coordinate {
    Temperature,
    LiquidFraction,
}

impl Coordinate {
    pub(super) fn value(self, phase: EquilibriumPhaseState) -> f64 {
        match self {
            Self::Temperature => phase.temperature_k(),
            Self::LiquidFraction => phase.liquid_mass_fraction(),
        }
    }

    fn slope(self, dt_dh: f64, df_dh: f64) -> f64 {
        match self {
            Self::Temperature => dt_dh,
            Self::LiquidFraction => df_dh,
        }
    }
}

#[derive(Debug, Clone)]
pub(super) enum Response<'a> {
    Constant,
    Curve {
        coordinate: Coordinate,
        knots: &'a [(f64, f64)],
    },
    Surface {
        coordinates: [Coordinate; 2],
        response: PreparedSurfaceResponse,
    },
}

#[derive(Debug, Clone)]
pub(in crate::radiation) struct MaterialEmissivitySource<'a> {
    pub(super) region: String,
    pub(super) snapshot: SurfaceEmissivity,
    pub(super) curve: &'a EquilibriumEnthalpyPhaseCurve,
    pub(super) material: ThermalMaterialPoint<'a>,
    pub(super) claim: ClaimId,
    pub(super) response: Response<'a>,
}

#[derive(Debug, Clone)]
pub(in crate::radiation) struct MaterialNodeResponse {
    pub state: MaterialRadiationNodeState,
    pub dt_dh: f64,
    pub depsilon_dh: f64,
}

impl<'a> MaterialEmissivitySource<'a> {
    pub(in crate::radiation) fn prepare(
        region: &str,
        curve: &'a EquilibriumEnthalpyPhaseCurve,
        material: ThermalMaterialPoint<'a>,
        selection: MaterialPropertySelection,
    ) -> Result<Self, ConductionError> {
        if curve.material_card_identity() != material.material.content_hash() {
            return Err(invalid(
                "surface emissivity and enthalpy require the same original material card",
            ));
        }
        let (snapshot, _) = SurfaceEmissivity::from_condition(region, material, selection)?;
        let claim = snapshot.receipt().selected;
        let original = material
            .material
            .claims()
            .claim(claim)
            .ok_or_else(|| invalid("prepared emissivity claim disappeared"))?;
        let coordinate = |axis: &str| {
            if axis == material.temperature_axis {
                Ok(Coordinate::Temperature)
            } else if Some(axis) == material.liquid_fraction_axis {
                Ok(Coordinate::LiquidFraction)
            } else {
                Err(invalid(
                    "material emissivity requires the declared temperature or liquid mass fraction axis",
                ))
            }
        };
        let response = match (&original.value, &original.interpolation) {
            (PropertyValue::Scalar { .. }, InterpolationPolicy::ConstantWithinValidity) => {
                Response::Constant
            }
            (
                PropertyValue::Curve {
                    abscissa, knots, ..
                },
                InterpolationPolicy::LinearInside,
            ) => {
                let coordinate = coordinate(abscissa)?;
                if knots.len() < 2 {
                    return Err(invalid(
                        "material emissivity requires at least two source curve knots",
                    ));
                }
                Response::Curve { coordinate, knots }
            }
            (PropertyValue::Surface2D { axes, .. }, InterpolationPolicy::LinearInside) => {
                let coordinates = [coordinate(&axes[0])?, coordinate(&axes[1])?];
                if coordinates[0] == coordinates[1] {
                    return Err(invalid(
                        "material emissivity surfaces require both temperature and liquid mass fraction",
                    ));
                }
                Response::Surface {
                    coordinates,
                    response: original.value.prepare_surface_response()
                        .map_err(|error| invalid(&error.to_string()))?,
                }
            }
            _ => {
                return Err(invalid(
                    "material radiation supports static scalars, linear thermal curves and temperature/fraction source surfaces",
                ));
            }
        };
        Ok(Self {
            region: region.into(),
            snapshot,
            curve,
            material,
            claim,
            response,
        })
    }

    pub(in crate::radiation) fn snapshot(&self) -> &SurfaceEmissivity {
        &self.snapshot
    }
    pub(in crate::radiation) fn curve(&self) -> &'a EquilibriumEnthalpyPhaseCurve {
        self.curve
    }
    pub(in crate::radiation) fn claim(&self) -> ClaimId {
        self.claim
    }

    pub(in crate::radiation) fn resolve(
        &self,
        vertex: usize,
        h: f64,
        supplied_temperature: Option<f64>,
        mode: DerivativeMode,
    ) -> Result<MaterialNodeResponse, ConductionError> {
        let failure = |upstream| ConductionError::MaterialQuery {
            property: SURFACE_EMISSIVITY_PROPERTY.into(),
            temperature: supplied_temperature.unwrap_or(f64::NAN),
            upstream,
        };
        let phase = self
            .curve
            .state_at_specific_enthalpy(h)
            .map_err(|e| failure(e.to_string()))?;
        let temperature = phase.temperature_k();
        if supplied_temperature.is_some_and(|value| value != temperature) {
            return Err(invalid(
                "material radiation temperature differs from its original caloric state",
            ));
        }
        let (emissivity, material) = SurfaceEmissivity::from_thermal_state(
            &self.region,
            self.material,
            phase,
            MaterialPropertySelection::PinnedByProperty(vec![(
                SURFACE_EMISSIVITY_PROPERTY.into(),
                self.claim,
            )]),
        )
        .map_err(|error| match error {
            ConductionError::MaterialQuery { .. } => error,
            _ => ConductionError::MaterialQuery {
                property: SURFACE_EMISSIVITY_PROPERTY.into(),
                temperature,
                upstream: error.to_string(),
            },
        })?;
        if mode == DerivativeMode::Strict {
            self.admit_smooth(vertex, h, phase)
                .map_err(ConductionError::RadiationMaterialSmoothness)?;
        }
        let mut dt_dh = 0.0;
        let mut depsilon_dh = 0.0;
        if mode != DerivativeMode::None {
            dt_dh = self
                .curve
                .temperature_derivative_at_specific_enthalpy(h)
                .map_err(|e| failure(e.to_string()))?;
            if let Response::Curve { coordinate, knots } = &self.response {
                let (value, coordinate_slope) = match coordinate {
                    Coordinate::Temperature => (temperature, dt_dh),
                    Coordinate::LiquidFraction => (
                        phase.liquid_mass_fraction(),
                        self.curve
                            .liquid_mass_fraction_derivative_at_specific_enthalpy(h)
                            .map_err(|e| failure(e.to_string()))?,
                    ),
                };
                if value < knots[0].0 || value > knots[knots.len() - 1].0 {
                    return Err(failure(
                        "emissivity derivative is outside its original source segments".into(),
                    ));
                }
                let upper = knots
                    .partition_point(|knot| knot.0 <= value)
                    .min(knots.len() - 1);
                let a = knots[upper - 1];
                let b = knots[upper];
                let slope = finite((b.1 - a.1) / (b.0 - a.0))?;
                if a.1 != b.1 && slope == 0.0 {
                    return Err(invalid(
                        "nonzero emissivity source derivative is unrepresentable",
                    ));
                }
                depsilon_dh = finite(slope * coordinate_slope)?;
                if slope != 0.0 && coordinate_slope != 0.0 && depsilon_dh == 0.0 {
                    return Err(invalid(
                        "nonzero emissivity enthalpy derivative is unrepresentable",
                    ));
                }
            }
            if let Response::Surface { coordinates, response } = &self.response {
                let df_dh = self.curve
                    .liquid_mass_fraction_derivative_at_specific_enthalpy(h)
                    .map_err(|e| failure(e.to_string()))?;
                let point = coordinates.map(|axis| axis.value(phase));
                let slopes = coordinates.map(|axis| axis.slope(dt_dh, df_dh));
                let varying = slopes.map(|slope| slope != 0.0);
                // A selected cell may only touch a held-coordinate slice at a
                // vertex. Use its unique active-slice gradient when admitted,
                // so the forward h tangent equals every strict adjoint tangent.
                // Nonsmooth primal trials retain the canonical one-sided cell.
                let gradient = match response.gradient_for_axes(point, varying) {
                    Ok(gradient) => gradient,
                    Err(_) if mode == DerivativeMode::Primal => response.evaluate(point)
                        .map_err(|e| failure(e.to_string()))?.gradient,
                    Err(error) => return Err(ConductionError::RadiationMaterialSmoothness(
                        MaterialRadiationSmoothnessError::SourceSurface {
                            vertex,
                            axes: response.axes().clone(),
                            coordinates: point,
                            upstream: error.to_string(),
                        },
                    )),
                };
                for axis in 0..2 {
                    let contribution = finite(gradient[axis] * slopes[axis])?;
                    if gradient[axis] != 0.0 && slopes[axis] != 0.0 && contribution == 0.0 {
                        return Err(invalid(
                            "nonzero emissivity enthalpy derivative is unrepresentable",
                        ));
                    }
                    depsilon_dh = finite(depsilon_dh + contribution)?;
                }
            }
        }
        Ok(MaterialNodeResponse {
            state: MaterialRadiationNodeState {
                vertex,
                phase,
                emissivity,
                material,
            },
            dt_dh,
            depsilon_dh,
        })
    }

    pub(in crate::radiation) fn admit_smooth(
        &self,
        vertex: usize,
        h: f64,
        phase: EquilibriumPhaseState,
    ) -> Result<(), MaterialRadiationSmoothnessError> {
        use MaterialRadiationSmoothnessError as E;
        let knots = self.curve.knots();
        if h <= knots[0].specific_enthalpy_j_kg
            || h >= knots[knots.len() - 1].specific_enthalpy_j_kg
        {
            return Err(E::ChartEndpoint {
                vertex,
                specific_enthalpy_j_kg: h,
            });
        }
        let upper = knots
            .partition_point(|knot| knot.specific_enthalpy_j_kg <= h)
            .min(knots.len() - 1);
        let slopes = |a: usize, b: usize| {
            let dh = knots[b].specific_enthalpy_j_kg - knots[a].specific_enthalpy_j_kg;
            (
                (knots[b].temperature_k - knots[a].temperature_k) / dh,
                (knots[b].liquid_mass_fraction - knots[a].liquid_mass_fraction) / dh,
            )
        };
        let (dt, df) = slopes(upper - 1, upper);
        if h == knots[upper - 1].specific_enthalpy_j_kg && upper > 1 {
            let before = slopes(upper - 2, upper - 1);
            if before.0 != dt || self.material.liquid_fraction_axis.is_some() && before.1 != df {
                return Err(E::ChartKink {
                    vertex,
                    specific_enthalpy_j_kg: h,
                });
            }
        }
        let original = self
            .material
            .material
            .claims()
            .claim(self.claim)
            .expect("pinned immutable source");
        for (axis, value, slope) in [
            (
                Some(self.material.temperature_axis),
                phase.temperature_k(),
                dt,
            ),
            (
                self.material.liquid_fraction_axis,
                phase.liquid_mass_fraction(),
                df,
            ),
        ] {
            if let Some(axis) = axis
                && let Some(&(low, high)) = original.validity.bounds().get(axis)
                && slope != 0.0
                && (value <= low || value >= high)
            {
                return Err(E::ValidityEndpoint {
                    vertex,
                    axis: axis.into(),
                    coordinate: value,
                });
            }
        }
        if let Response::Curve { coordinate, knots } = &self.response {
            let (axis, value, slope) = match coordinate {
                Coordinate::Temperature => (self.material.temperature_axis, phase.temperature_k(), dt),
                Coordinate::LiquidFraction => (
                    self.material
                        .liquid_fraction_axis
                        .expect("prepared fraction axis"),
                    phase.liquid_mass_fraction(),
                    df,
                ),
            };
            // A strictly smooth chart may hold this coordinate constant on an
            // open h neighborhood: latent T or pure-phase fraction. Its source
            // composition then has zero h slope even at an original source
            // knot or endpoint; ordinate controls still use that exact knot.
            if slope == 0.0 { return Ok(()); }
            if value <= knots[0].0 || value >= knots[knots.len() - 1].0 {
                return Err(E::SourceEndpoint {
                    vertex,
                    axis: axis.into(),
                    coordinate: value,
                });
            }
            if let Ok(i) = knots.binary_search_by(|knot| knot.0.total_cmp(&value)) {
                let left = (knots[i].1 - knots[i - 1].1) / (knots[i].0 - knots[i - 1].0);
                let right = (knots[i + 1].1 - knots[i].1) / (knots[i + 1].0 - knots[i].0);
                if left != right {
                    return Err(E::SourceKink {
                        vertex,
                        axis: axis.into(),
                        coordinate: value,
                    });
                }
            }
        }
        if let Response::Surface { coordinates, response } = &self.response {
            let point = coordinates.map(|axis| axis.value(phase));
            let varying = coordinates.map(|axis| axis.slope(dt, df) != 0.0);
            response.gradient_for_axes(point, varying).map_err(|error| E::SourceSurface {
                vertex,
                axes: response.axes().clone(),
                coordinates: point,
                upstream: error.to_string(),
            })?;
        }
        Ok(())
    }
}
