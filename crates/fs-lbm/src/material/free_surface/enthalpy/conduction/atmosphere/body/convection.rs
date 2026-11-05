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
    /// Signed material integral-coordinate heat, positive from solid to gas.
    /// Retained directly to avoid subtracting rounded, nearly equal temperatures.
    pub solid_power_to_gas_w: f64,
    /// Material integral power minus actual gas-side power [W].
    pub surface_heat_residual_w: f64,
    pub iterations: usize,
    pub source_line: Arc<PropertyLineSlice>,
    /// None denotes only exact equilibrium admitted by the original provider.
    pub coefficient: Option<Evidence<FreeSurfaceHeatTransferCoefficient>>,
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
        body_id: u64,
        face: FreeSurfaceBodyGasFace,
        solid: ResolvedMaterialStatePoint,
        specific_h_j_kg: f64,
        parcel: &IsochoricGasParcel,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<FreeSurfaceBodyGasSample, Error> {
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
        let response = geometry.evaluate_convective_material_to_bath(
            self.layer.conductivity_curve(),
            [solid_temperature_k, gas_temperature_k],
            self.config,
            |surface_temperature_k| {
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
            },
            check,
        ).map_err(|error| match error {
            ConvectiveLayerError::Thermal(error) => Error::from(error),
            ConvectiveLayerError::Film(error) => error,
            ConvectiveLayerError::SurfaceBalance { residual_w, iterations } => {
                Error::BodyConvectionBalance { body_id, residual_w, iterations }
            }
        })?;
        Ok(FreeSurfaceBodyGasSample {
            body_id,
            face,
            solid_temperature_k,
            gas_temperature_k,
            solid,
            gas: response.film.0,
            conductance_w_k: response.conductance_w_k,
            heat_into_solid_w: -response.power_to_bath_w,
            convection: Some(FreeSurfaceBodyConvectionSample {
                surface_temperature_k: response.surface_temperature_k,
                solid_power_to_gas_w: response.solid_power_to_bath_w,
                surface_heat_residual_w: response.heat_residual_w,
                iterations: response.iterations,
                source_line: Arc::clone(self.layer.source_line()),
                coefficient: response.film.1,
            }),
        })
    }
}
