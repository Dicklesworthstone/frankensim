//! A conducting solid face and its original finite gas share one heat rate.

use super::*;
use fs_evidence::Evidence;
use fs_matdb::PropertyLineSlice;
use fs_material::thermal::conduction::{
    ConvectiveLayerError, LayerConductionConfig, ReservoirLayerGeometry,
};
use fs_material::thermal::conduction::source::ResolvedConductionLayer;
use std::sync::Arc;

use super::super::super::surface::convection::FreeSurfaceHeatTransferCoefficient;
use super::super::super::radiation::{self, FreeSurfaceBodyRadiationSample};

/// Actual signed stationary wall geometry and original gas states at one trial.
#[derive(Clone, Copy, Debug)]
pub struct FreeSurfaceBodyConvectionQuery<'a> {
    pub body_id: u64,
    pub face: FreeSurfaceBodyGasFace,
    /// Unit normal from the original wall cell into its adjacent dry gas cell.
    pub outward_normal: [f64; 3],
    pub surface_temperature_k: f64,
    /// The original finite atmosphere at its current EOS pressure.
    pub bulk: &'a ResolvedThermochemicalGasState,
    /// The SAME gas claims at (surface T + bulk T)/2 and current EOS pressure.
    pub film: &'a ResolvedThermochemicalGasState,
}

/// Repeatable gas-side constitutive law, evaluated at the actual solid surface.
///
/// The provider owns correlation geometry, physical boundary-layer length and
/// regime admission, including equilibrium. `None` is admitted only at exact
/// surface/bulk equilibrium; every driven coefficient must be finite and positive.
/// Calls include rejected trials and transactional clones. The provider must
/// not advance physical history or integrate heat. A refusal propagates through
/// the original solid/liquid/gas rollback; no stagnant-film fallback occurs.
pub trait FreeSurfaceBodyGasConvection: core::fmt::Debug {
    fn coefficient(
        &self,
        query: FreeSurfaceBodyConvectionQuery<'_>,
    ) -> Result<Option<Evidence<FreeSurfaceHeatTransferCoefficient>>, FreeSurfaceEnthalpyError>;
}

/// The exact prepared solid source and converged gas-side evidence for a face.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceBodyConvectionSample {
    pub surface_temperature_k: f64,
    /// Signed material integral-coordinate heat, positive from cell to surface.
    /// Retained directly to avoid subtracting rounded, nearly equal temperatures.
    pub solid_power_to_surface_w: f64,
    /// Material integral power minus total gas/radiant surface power [W].
    pub surface_heat_residual_w: f64,
    pub iterations: usize,
    pub source_line: Arc<PropertyLineSlice>,
    /// None denotes only exact equilibrium admitted by the original provider.
    pub coefficient: Option<Evidence<FreeSurfaceHeatTransferCoefficient>>,
    /// Gas/radiant-bath pair after eliminating the shared material surface.
    pub gas_radiation_conductance_w_k: f64,
    /// The SAME solved surface's radiant rate and original emissivity. The
    /// radiation timeline consumes this sample once, without another root.
    pub radiation: Option<FreeSurfaceBodyRadiationSample>,
}

#[derive(Clone, Debug)]
pub(in super::super::super) struct PreparedBodyGasConvection<'m> {
    pub(in super::super::super) provider: &'m dyn FreeSurfaceBodyGasConvection,
    pub(in super::super::super) layer: Arc<ResolvedConductionLayer<'m>>,
    pub(in super::super::super) config: LayerConductionConfig,
}

impl PreparedBodyGasConvection<'_> {
    pub(in super::super::super) fn sample(
        &self,
        body: &FreeSurfaceThermalBody<'_>,
        face: FreeSurfaceBodyGasFace,
        radiation_face: Option<FreeSurfaceBodyRadiationFace>,
        solid: ResolvedMaterialStatePoint,
        specific_h_j_kg: f64,
        parcel: &IsochoricGasParcel,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<FreeSurfaceBodyGasSample, Error> {
        let body_id = body.id();
        let solid_temperature_k = self.layer.temperature_at_enthalpy(specific_h_j_kg)?;
        let bulk = parcel.state();
        let gas_temperature_k = bulk.state().temperature;
        // Admission has established a unit Cartesian adjacency. Signed wall
        // coordinates remain intact, including the outer-domain wall at -1.
        let outward_normal = std::array::from_fn(|axis| {
            if face.wall[axis] < face.gas[axis] as i64 {
                1.0
            } else if face.wall[axis] > face.gas[axis] as i64 {
                -1.0
            } else {
                0.0
            }
        });
        let geometry = ReservoirLayerGeometry::planar(face.area_m2, face.distances_m[0])?;
        let gas_film = |surface_temperature_k| {
                let point = parcel.point().clone()
                    .with("temperature", 0.5 * surface_temperature_k + 0.5 * gas_temperature_k)
                    .map_err(|_| Error::Invalid("invalid conducting solid gas film temperature"))?;
                let film = parcel.model().evaluate(&point).map_err(MaterialLbmError::from)?;
                if !environment::same_sources(bulk, &film) {
                    return Err(Error::Invalid("conducting solid film left its finite gas source claims"));
                }
                let coefficient = self.provider.coefficient(FreeSurfaceBodyConvectionQuery {
                    body_id,
                    face,
                    outward_normal,
                    surface_temperature_k,
                    bulk,
                    film: &film,
                })?;
                let conductance = match coefficient.as_ref() {
                    Some(value) if value.value.value().is_finite() && value.value.value() > 0.0 => {
                        let conductance = environment::finite(face.area_m2 * value.value.value())?;
                        if conductance <= 0.0 {
                            return Err(DepositionError::Arithmetic.into());
                        }
                        conductance
                    }
                    None if surface_temperature_k == gas_temperature_k => 0.0,
                    _ => return Err(Error::Invalid("solid convection needs positive finite h or exact equilibrium")),
                };
                let power = environment::finite(conductance * (surface_temperature_k - gas_temperature_k))?;
                Ok((power, conductance, (film, coefficient)))
            };
        if let Some(radiation_face) = radiation_face {
            let response = geometry.evaluate_material_to_two_baths(
                self.layer.conductivity_curve(),
                solid_temperature_k,
                [gas_temperature_k, radiation_face.surroundings_k],
                self.layer.conductivity_curve().temperature_range_k(),
                self.config,
                |surface_temperature_k| {
                    let (power, conductance, (film, coefficient)) = gas_film(surface_temperature_k)?;
                    let emissivity_source = body.emissivity_at_surface(
                        specific_h_j_kg, surface_temperature_k,
                    )?;
                    let radiant = radiation::radiant_conductance(
                        emissivity_source.value_si(), radiation_face.area_m2,
                        surface_temperature_k, radiation_face.surroundings_k,
                    )?;
                    let radiant_power = environment::finite(
                        radiant * (surface_temperature_k - radiation_face.surroundings_k),
                    )?;
                    Ok(([power, radiant_power], [conductance, radiant],
                        (film, coefficient, emissivity_source)))
                },
                check,
            ).map_err(|error| balance_error(body_id, error))?;
            let (gas, coefficient, emissivity_source) = response.film;
            let radiation = FreeSurfaceBodyRadiationSample {
                body_id,
                face: radiation_face,
                solid_temperature_k,
                surface_temperature_k: response.surface_temperature_k,
                emissivity: emissivity_source.value_si(),
                emissivity_source,
                conductance_w_k: response.material_bath_conductances_w_k[1],
                heat_into_solid_w: -response.powers_to_baths_w[1],
            };
            return Ok(FreeSurfaceBodyGasSample {
                body_id,
                face,
                solid_temperature_k,
                gas_temperature_k,
                solid,
                gas,
                conductance_w_k: response.material_bath_conductances_w_k[0],
                gas_decay_conductance_w_k: environment::finite(
                    response.material_bath_conductances_w_k[0] + response.bath_pair_conductance_w_k,
                )?,
                heat_into_solid_w: -response.powers_to_baths_w[0],
                convection: Some(FreeSurfaceBodyConvectionSample {
                    surface_temperature_k: response.surface_temperature_k,
                    solid_power_to_surface_w: response.solid_power_to_surface_w,
                    surface_heat_residual_w: response.heat_residual_w,
                    iterations: response.iterations,
                    source_line: Arc::clone(self.layer.source_line()),
                    coefficient,
                    gas_radiation_conductance_w_k: response.bath_pair_conductance_w_k,
                    radiation: Some(radiation),
                }),
            });
        }
        let response = geometry.evaluate_convective_material_to_bath(
            self.layer.conductivity_curve(),
            [solid_temperature_k, gas_temperature_k],
            self.config,
            gas_film,
            check,
        ).map_err(|error| balance_error(body_id, error))?;
        Ok(FreeSurfaceBodyGasSample {
            body_id,
            face,
            solid_temperature_k,
            gas_temperature_k,
            solid,
            gas: response.film.0,
            conductance_w_k: response.conductance_w_k,
            gas_decay_conductance_w_k: response.conductance_w_k,
            heat_into_solid_w: -response.power_to_bath_w,
            convection: Some(FreeSurfaceBodyConvectionSample {
                surface_temperature_k: response.surface_temperature_k,
                solid_power_to_surface_w: response.solid_power_to_bath_w,
                surface_heat_residual_w: response.heat_residual_w,
                iterations: response.iterations,
                source_line: Arc::clone(self.layer.source_line()),
                coefficient: response.film.1,
                gas_radiation_conductance_w_k: 0.0,
                radiation: None,
            }),
        })
    }
}

fn balance_error(body_id: u64, error: ConvectiveLayerError<Error>) -> Error {
    match error {
        ConvectiveLayerError::Thermal(error) => Error::from(error),
        ConvectiveLayerError::Film(error) => error,
        ConvectiveLayerError::SurfaceBalance { residual_w, iterations } => {
            Error::BodyConvectionBalance { body_id, residual_w, iterations }
        }
    }
}
