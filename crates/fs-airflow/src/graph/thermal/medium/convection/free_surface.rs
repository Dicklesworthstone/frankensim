//! Source-backed plate convection on the original moving-liquid heat owner.
//!
//! The liquid solver supplies its actual trial interface T, original gas at
//! bulk T, and the same source claims at mean-film T. The existing correlation
//! owner computes every group and h from that one film state. The original
//! liquid-side resistance, radiation and finite-gas energy owner remain active.

use super::*;
use fs_lbm::material::free_surface::enthalpy::FreeSurfaceEnthalpyError;
use fs_lbm::material::free_surface::enthalpy::conduction::surface::convection::{
    FreeSurfaceConvectionQuery, FreeSurfaceGasConvection, FreeSurfaceHeatTransferCoefficient,
};
use fs_material::fluid_transport::{AidingVerticalPlateFrame, ThermalFluidTransportState};

/// Explicitly selected average flat-plate forced-convection correlation.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum PlanarForcedConvectionRegime {
    /// Original laminar flat-plate average, including its Reynolds/Prandtl limits.
    Laminar,
    /// Original turbulent average including its laminar leading-edge correction.
    Turbulent,
}

/// One declared external-flow model, held throughout a thermal/native-flow split.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum PlanarFreeSurfaceHeatLaw {
    Forced {
        regime: PlanarForcedConvectionRegime,
        /// Physical gas velocity RELATIVE to the liquid surface [m/s]. Must be
        /// tangent to the plate. The caller supplies the actual relative field;
        /// this bounded model does not infer gas momentum or moving-wall drag.
        relative_velocity_m_s: [f64; 3],
    },
    /// Quiescent gas and physical gravity; only vertical/horizontal plates.
    /// Source-owned thermal expansion selects hot/cold stable/unstable buoyancy.
    Natural { gravity_m_s2: [f64; 3] },
    /// Laminar aiding natural/forced flow on an isothermal vertical plate.
    ///
    /// All properties come from the actual mean-film state, while buoyancy
    /// uses the full interface/bulk temperature contrast. This frozen-property
    /// approximation is an adapter policy; Manning–Qureshi assumes constant
    /// properties. The density bound does not establish small variation of the
    /// other properties or extend the original correlation's domain.
    AidingVerticalPlate {
        /// Physical plate-up direction, tangent to the actual face. Need not
        /// be normalized. Nonzero gravity must be opposite this direction.
        plate_up: [f64; 3],
        /// Actual gas velocity relative to the moving liquid [m/s], parallel
        /// to plate-up. Its sign must aid the source-owned buoyant response.
        relative_velocity_m_s: [f64; 3],
        /// Physical gravity [m/s²]. Zero retains the pure forced limit.
        gravity_m_s2: [f64; 3],
        /// Caller-admitted bound on |beta * (T_surface - T_bulk)| in [0, 1).
        maximum_relative_density_departure: f64,
    },
}

/// A declared physical isothermal plate patch in the moving liquid's world frame.
///
/// `characteristic_length` is the WHOLE physical plate's streamwise length for
/// forced flow, height for a vertical natural/mixed plate, or area/perimeter for a
/// horizontal natural plate. It must not become a fresh thermal entrance at
/// each mesh face. A provider can be shared by several faces of that same plate;
/// every face retains its physical integration area and original liquid owner.
///
/// Applying a plate-average relation locally assumes the declared patch is
/// approximately isothermal with locally held film properties. Spatially varying
/// wall temperatures require refinement/model assessment; this is not a resolved
/// gas boundary layer. Unsupported regimes fail instead of switching models.
/// The authored unit normal must match the actual liquid-to-gas normal, including
/// reconstructed planar interfaces. Unit length, direction agreement and
/// tangency admit only `64 * f64::EPSILON` normalization roundoff; this is not
/// an empirical angular allowance. Gas speed and gravity are never projected
/// onto an otherwise unsupported orientation.
///
/// Curved, enclosed and evaporating surfaces are outside this adapter. Natural
/// convection admits only vertical/horizontal plates relative to actual gravity.
/// Mixed convection admits only the declared aiding vertical-plate model, not
/// opposing flow or crossflow. Geometry and area remain with the liquid owner.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanarFreeSurfaceConvection {
    /// Physical unit normal, outward from the liquid into gas.
    pub outward_normal: [f64; 3],
    pub characteristic_length: Length,
    pub law: PlanarFreeSurfaceHeatLaw,
}

fn invalid(what: &'static str) -> FreeSurfaceEnthalpyError {
    FreeSurfaceEnthalpyError::Convection(what.into())
}

fn original(error: impl core::fmt::Display) -> FreeSurfaceEnthalpyError {
    FreeSurfaceEnthalpyError::Convection(error.to_string())
}

fn thermal_film(
    query: &FreeSurfaceConvectionQuery<'_>,
) -> Result<ThermalFluidTransportState, FreeSurfaceEnthalpyError> {
    let thermal = query.film.thermal_transport_state().map_err(original)?;
    let transport = query.film.transport_state().map_err(original)?;
    if thermal.transport() != transport || thermal.temperature_k() != query.film.state().temperature
    {
        return Err(invalid(
            "film thermal and momentum views must describe the same source state",
        ));
    }
    Ok(thermal)
}

// One roundoff allowance for dimensionless unit length and orientation tests.
// Comparing vector differences, rather than 1 - dot(a,b), does not admit a
// sqrt(epsilon)-sized tilt when nearly parallel dot products round to one.
const DIRECTION_ROUNDOFF: f64 = 64.0 * f64::EPSILON;

fn direction(vector: [f64; 3]) -> Result<([f64; 3], f64), FreeSurfaceEnthalpyError> {
    if vector.into_iter().any(|x| !x.is_finite()) {
        return Err(invalid(
            "plate convection requires a finite physical vector",
        ));
    }
    let scale = vector.into_iter().map(f64::abs).fold(0.0, f64::max);
    if scale == 0.0 {
        return Err(invalid(
            "plate convection requires a nonzero finite physical vector",
        ));
    }
    let scaled = vector.map(|x| x / scale);
    let length = fs_math::det::sqrt(scaled.into_iter().map(|x| x * x).sum());
    // Plate-up is only a direction and can be normalized even if its raw
    // magnitude overflows. Physical speed/gravity and unit normals separately
    // admit their returned magnitude before using it.
    Ok((scaled.map(|x| x / length), scale * length))
}

fn physical_vector(vector: [f64; 3]) -> Result<([f64; 3], f64), FreeSurfaceEnthalpyError> {
    let (unit, magnitude) = direction(vector)?;
    if !magnitude.is_finite() || magnitude <= 0.0 {
        return Err(invalid(
            "plate convection requires a nonzero finite physical vector",
        ));
    }
    Ok((unit, magnitude))
}

fn unit_normal(vector: [f64; 3]) -> Result<[f64; 3], FreeSurfaceEnthalpyError> {
    let (unit, magnitude) = physical_vector(vector)?;
    if (magnitude - 1.0).abs() > DIRECTION_ROUNDOFF {
        return Err(invalid("plate convection requires physical unit normals"));
    }
    Ok(unit)
}

fn aligned(a: [f64; 3], b: [f64; 3]) -> bool {
    a.into_iter()
        .zip(b)
        .map(|(x, y)| (x - y) * (x - y))
        .sum::<f64>()
        <= DIRECTION_ROUNDOFF * DIRECTION_ROUNDOFF
}

fn tangent(normal: [f64; 3], direction: [f64; 3]) -> bool {
    normal
        .into_iter()
        .zip(direction)
        .map(|(x, y)| x * y)
        .sum::<f64>()
        .abs()
        <= DIRECTION_ROUNDOFF
}

impl FreeSurfaceGasConvection for PlanarFreeSurfaceConvection {
    fn coefficient(
        &self,
        query: FreeSurfaceConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>
    {
        let normal = unit_normal(self.outward_normal)?;
        let actual_normal = unit_normal(query.outward_normal)?;
        if !aligned(normal, actual_normal) {
            return Err(invalid(
                "the declared plate normal must match its actual liquid/gas surface",
            ));
        }
        if !self.characteristic_length.value().is_finite()
            || self.characteristic_length.value() <= 0.0
        {
            return Err(invalid(
                "plate convection requires a positive physical characteristic length",
            ));
        }
        let wall = query.surface_temperature_k;
        let bulk = query.bulk.state().temperature;
        if [wall, bulk].into_iter().any(|x| !x.is_finite() || x <= 0.0)
            || query.film.state().temperature != 0.5 * wall + 0.5 * bulk
        {
            return Err(invalid(
                "plate convection requires original bulk and actual mean-film temperatures",
            ));
        }
        let evaluated = match self.law {
            PlanarFreeSurfaceHeatLaw::Forced {
                regime,
                relative_velocity_m_s,
            } => {
                let (flow, speed) = physical_vector(relative_velocity_m_s)?;
                if !tangent(actual_normal, flow) {
                    return Err(invalid(
                        "forced plate convection requires relative gas flow tangent to the surface",
                    ));
                }
                let correlation = match regime {
                    PlanarForcedConvectionRegime::Laminar => CorrelationId::FlatPlateLaminarAverage,
                    PlanarForcedConvectionRegime::Turbulent => {
                        CorrelationId::FlatPlateTurbulentAverage
                    }
                };
                evaluate_forced_from_material(
                    correlation,
                    query.film,
                    speed,
                    self.characteristic_length,
                    ForcedConvectionSettings {
                        direction: if wall < bulk {
                            ThermalDirection::CoolingFluid
                        } else {
                            ThermalDirection::HeatingFluid
                        },
                        ..ForcedConvectionSettings::default()
                    },
                )
                .map_err(original)?
            }
            PlanarFreeSurfaceHeatLaw::Natural { gravity_m_s2 } => {
                let (down, gravity) = physical_vector(gravity_m_s2)?;
                let geometry = if tangent(actual_normal, down) {
                    NaturalConvectionGeometry::VerticalPlate
                } else if aligned(actual_normal, down) {
                    NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Down)
                } else if aligned(actual_normal.map(|x| -x), down) {
                    NaturalConvectionGeometry::HorizontalPlate(HorizontalPlateFacing::Up)
                } else {
                    return Err(invalid(
                        "natural plate convection supports only gravity-aligned horizontal or vertical plates",
                    ));
                };
                let thermal = thermal_film(&query)?;
                let transport = thermal.transport();
                let beta = thermal.volumetric_expansion_per_k().ok_or_else(|| {
                    invalid("natural plate convection requires source-owned thermal expansion")
                })?;
                // Admit all physical driving inputs even at exact equilibrium.
                // No finite-Ra correlation value or uncertainty is invented there.
                natural_rayleigh(
                    transport,
                    gravity,
                    beta,
                    Temperature::new(wall),
                    Temperature::new(bulk),
                    self.characteristic_length,
                )
                .map_err(original)?;
                if wall == bulk {
                    return Ok(None);
                }
                evaluate_natural_from_material(
                    geometry,
                    query.film,
                    gravity,
                    beta,
                    Temperature::new(wall),
                    Temperature::new(bulk),
                    self.characteristic_length,
                )
                .map_err(original)?
            }
            PlanarFreeSurfaceHeatLaw::AidingVerticalPlate {
                plate_up,
                relative_velocity_m_s,
                gravity_m_s2,
                maximum_relative_density_departure,
            } => {
                let (up, _) = direction(plate_up)?;
                let tangent_flow = relative_velocity_m_s == [0.0; 3]
                    || tangent(actual_normal, direction(relative_velocity_m_s)?.0);
                if !tangent(actual_normal, up) || !tangent_flow {
                    return Err(invalid(
                        "aiding vertical-plate axis and relative gas flow must be tangent to the actual face",
                    ));
                }
                let settings = AidingVerticalPlateFrame::try_new(plate_up, gravity_m_s2)
                    .and_then(|frame| {
                        frame.settings(relative_velocity_m_s, maximum_relative_density_departure)
                    })
                    .map_err(original)?;
                let thermal = thermal_film(&query)?;
                let transport = thermal.transport();
                let beta = thermal.volumetric_expansion_per_k().ok_or_else(|| {
                    invalid("mixed plate convection requires source-owned thermal expansion")
                })?;
                let groups = settings
                    .groups(
                        transport,
                        beta,
                        Temperature::new(wall),
                        Temperature::new(bulk),
                        self.characteristic_length,
                    )
                    .map_err(original)?;
                // Admit frame, source and density policy even at equilibrium.
                // No coefficient is invented at the card's Re = Ra = 0 point.
                if wall == bulk && groups.reynolds == 0.0 {
                    return Ok(None);
                }
                finish(
                    CorrelationId::ManningQureshiAidingVerticalPlate,
                    query.film,
                    transport,
                    self.characteristic_length,
                    CorrelationInputs::mixed(groups.reynolds, groups.rayleigh, groups.prandtl),
                )
                .map_err(original)?
            }
        };
        let (_, _, coefficient) = evaluated.into_parts();
        Ok(Some(coefficient))
    }

    fn admits_uniform_boundary_evolution(&self) -> bool {
        // One authored plate frame, whole-plate length and external-flow law
        // apply independently of old patch indices. Newly reconstructed faces
        // must still pass every actual-normal and material-state check above.
        true
    }
}

#[cfg(test)]
mod tests;
