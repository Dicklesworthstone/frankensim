//! Solid and liquid share one finite source-backed gas inventory.
//!
//! Film geometry belongs to the specimen; solid conductivity and gas transport
//! retain their original owners. The existing thermal timeline advances every
//! finite store, and its one aggregate gas debit includes both kinds of patch.

use super::*;
use fs_material::thermochemical_gas::ResolvedThermochemicalGasState;
use super::super::radiation::FreeSurfaceBodyRadiationFace;

mod convection;
pub(in super::super) use convection::PreparedBodyGasConvection;
pub use convection::{
    FreeSurfaceBodyConvectionQuery, FreeSurfaceBodyConvectionSample, FreeSurfaceBodyGasConvection,
};

/// Exposed stationary solid/gas face, attached to its finite solid cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceBodyGasFace {
    /// Native stationary wall cell, including a signed outer-domain wall.
    pub wall: [i64; 3],
    /// Face-adjacent dry native cell at the accepted flow configuration.
    pub gas: [usize; 3],
    /// Patch area [m2], at most one Cartesian face; zero insulates a prescribed
    /// gas film. A conducting convection surface requires positive area.
    pub area_m2: f64,
    /// Solid-centroid distance and explicit continuum gas-film thickness [m].
    /// Both are positive for a prescribed film. Source convection requires a
    /// positive solid distance and ZERO gas thickness, replacing that film.
    pub distances_m: [f64; 2],
}

/// Original sources and rates actually used at one solid/gas patch.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceBodyGasSample {
    pub body_id: u64,
    pub face: FreeSurfaceBodyGasFace,
    pub solid_temperature_k: f64,
    pub gas_temperature_k: f64,
    pub solid: ResolvedMaterialStatePoint,
    /// Gas transport at the actual surface/bulk mean film temperature and
    /// current EOS pressure. A prescribed film uses the solid cell as surface.
    pub gas: ResolvedThermochemicalGasState,
    /// Material/gas pair conductance after eliminating the common surface.
    pub conductance_w_k: f64,
    /// Gas row includes both material/gas and gas/radiant-bath exchange.
    /// This is a frozen secant scale, not a derivative of the original source.
    pub gas_decay_conductance_w_k: f64,
    /// Positive into the solid; its opposite enters the one gas account [W].
    pub heat_into_solid_w: f64,
    /// Present only for an explicitly selected conducting convection surface.
    pub convection: Option<FreeSurfaceBodyConvectionSample>,
}

/// Direct gas/solid transfer, separate from the solid/liquid contact account.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceBodyGasExchange {
    pub heat_to_solids_j: f64,
    pub last_applied: Vec<FreeSurfaceBodyGasSample>,
    pub thermal_endpoint: Vec<FreeSurfaceBodyGasSample>,
}

#[derive(Clone, Default)]
pub(in super::super) struct BodyAtmosphere {
    pub(in super::super) faces: Vec<(usize, FreeSurfaceBodyGasFace)>,
    radiation: Vec<Option<FreeSurfaceBodyRadiationFace>>,
    pub(in super::super) reference: Vec<FreeSurfaceBodyGasSample>,
    pub(in super::super) last: Vec<FreeSurfaceBodyGasSample>,
    pub(in super::super) endpoint: Vec<FreeSurfaceBodyGasSample>,
    pub(in super::super) heat_j: f64,
    pub(in super::super) absolute_heat_j: f64,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(in super::super) fn admit_body_atmosphere(
        &self,
        bodies: &[FreeSurfaceThermalBody<'_>],
        finite: bool,
        dx: f64,
        gate: &CancelGate,
    ) -> Result<BodyAtmosphere, Error> {
        let mut timeline = BodyAtmosphere::default();
        let mut faces = BTreeSet::new();
        let mut radiating = BTreeMap::new();
        if bodies.iter().any(|body| body.gas_convection().is_some()) {
            for (index, body) in bodies.iter().enumerate() {
                for face in body.radiation_faces() {
                    poll(gate)?;
                    if face.area_m2 > 0.0
                        && radiating.insert((face.wall, face.gas), (index, *face)).is_some()
                    {
                        return Err(Error::Invalid("a radiating native wall face needs one owner"));
                    }
                }
            }
        }
        for (index, body) in bodies.iter().enumerate() {
            let convective = body.gas_convection().is_some();
            for &face in body.gas_faces() {
                poll(gate)?;
                let [Ok(x), Ok(y), Ok(z)] = face.gas.map(i64::try_from) else {
                    return Err(Error::Invalid("solid gas coordinate is unrepresentable"));
                };
                let gas = [x, y, z];
                let distance = (0..3).fold(0u64, |sum, axis| {
                    sum.saturating_add(face.wall[axis].abs_diff(gas[axis]))
                });
                if !finite
                    || distance != 1
                    || !self.flow.lattice.transport_wall_at(face.wall)
                    || self.flow.lattice.transport_wall_at(gas)
                    || self
                        .flow
                        .lattice
                        .sample_at(face.gas[0], face.gas[1], face.gas[2])
                        .is_some_and(|cell| cell.fill_fraction != 0.0)
                    || !faces.insert((face.wall, face.gas))
                    || !face.area_m2.is_finite()
                    || face.area_m2 < 0.0
                    || face.area_m2 > dx * dx
                    || !face.distances_m[0].is_finite()
                    || face.distances_m[0] <= 0.0
                    || !face.distances_m[1].is_finite()
                    || if convective {
                        face.area_m2 <= 0.0 || face.distances_m[1] != 0.0
                    } else {
                        face.distances_m[1] <= 0.0
                    }
                {
                    return Err(Error::Invalid(
                        "solid film needs the shared finite gas and a unique adjacent exposed wall face",
                    ));
                }
                if convective
                    && let Some((radiating_body, radiation)) = radiating.get(&(face.wall, face.gas))
                    && (*radiating_body != index || radiation.area_m2 != face.area_m2)
                {
                    return Err(Error::Invalid(
                        "a shared conducting surface needs the same body and exact convection/radiation area",
                    ));
                }
                timeline.faces.push((index, face));
            }
        }
        timeline
            .faces
            .sort_by_key(|(index, face)| (bodies[*index].id(), face.wall, face.gas));
        for (index, face) in &timeline.faces {
            poll(gate)?;
            timeline.radiation.push(if bodies[*index].gas_convection().is_some() {
                radiating.get(&(face.wall, face.gas)).map(|(_, face)| *face)
            } else {
                None
            });
        }
        Ok(timeline)
    }
}

impl Conduction<'_, '_> {
    pub(in super::super) fn sample_body_atmosphere(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreeSurfaceBodyGasSample>, Error> {
        if self.body_atmosphere.faces.is_empty() {
            return Ok(Vec::new());
        }
        let atmosphere = self.finite_atmosphere.as_ref().ok_or(Error::Invalid(
            "solid gas faces require their shared finite atmosphere",
        ))?;
        let mut samples = Vec::with_capacity(self.body_atmosphere.faces.len());
        let mut covered = BTreeMap::new();
        for (sample_index, &(body_index, face)) in self.body_atmosphere.faces.iter().enumerate() {
            check()?;
            let body = &self.bodies[body_index];
            let id = self.body_ids[body_index];
            let area = covered.entry(id).or_insert(0.0);
            *area = environment::finite(*area + face.area_m2)?;
            // The original thermal owner subtracts all solid/solid and
            // solid/liquid contact coverage, including insulating links.
            if *area > thermal.exposed_area_m2(id)? {
                return Err(Error::Invalid(
                    "solid gas films exceed the uncovered solid area",
                ));
            }
            let solid = body.resolve(thermal.cells()[&id].specific_h_j_kg)?;
            if let Some(convection) = body.gas_convection() {
                let sample = convection.sample(body, face, self.body_atmosphere.radiation[sample_index], solid,
                    thermal.cells()[&id].specific_h_j_kg, &atmosphere.parcel, check)?;
                if self.body_atmosphere.reference.get(sample_index)
                    .is_some_and(|start| !environment::same_sources(&start.gas, &sample.gas))
                {
                    return Err(Error::Invalid("solid convection left its original source claims"));
                }
                samples.push(sample);
                continue;
            }
            let solid_temperature_k = thermal.temperature_k(id)?;
            let gas_temperature_k = atmosphere.parcel.state().state().temperature;
            let film = 0.5 * solid_temperature_k + 0.5 * gas_temperature_k;
            let point = atmosphere
                .parcel
                .point()
                .clone()
                .with("temperature", film)
                .map_err(|_| Error::Invalid("invalid solid gas film temperature"))?;
            let gas = atmosphere
                .parcel
                .model()
                .evaluate(&point)
                .map_err(MaterialLbmError::from)?;
            if self
                .body_atmosphere
                .reference
                .get(sample_index)
                .is_some_and(|start| !environment::same_sources(&start.gas, &gas))
            {
                return Err(Error::Invalid(
                    "solid gas film left its original source claims",
                ));
            }
            let conductance_w_k = series_conductance_w_k(
                face.area_m2,
                face.distances_m,
                [body.conductivity(&solid), gas.state().thermal_conductivity],
            )?;
            samples.push(FreeSurfaceBodyGasSample {
                body_id: body.id(),
                face,
                solid_temperature_k,
                gas_temperature_k,
                solid,
                gas,
                conductance_w_k,
                gas_decay_conductance_w_k: conductance_w_k,
                heat_into_solid_w: environment::finite(
                    conductance_w_k * (gas_temperature_k - solid_temperature_k),
                )?,
                convection: None,
            });
        }
        Ok(samples)
    }
}
