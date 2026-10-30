//! Original solid emissivity on the same finite wall/liquid/gas heat timeline.

use super::*;

/// An exposed stationary solid patch facing black isothermal surroundings.
///
/// The source-bound total hemispherical emissivity belongs to its finite body.
/// This is the same opaque diffuse-gray, unit-view-factor approximation as the
/// liquid boundary. Surface T is the body's resolved cell T; callers must
/// resolve the boundary layer or justify that uniform-temperature approximation.
/// No gas absorption, finite radiative reservoir or inferred visibility is added.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FreeSurfaceBodyRadiationFace {
    pub wall: [i64; 3],
    /// Face-adjacent dry native cell; wet cells and permanent walls refuse.
    pub gas: [usize; 3],
    /// Radiating area [m2], no greater than one Cartesian face.
    pub area_m2: f64,
    /// Prescribed surrounding temperature [K], positive.
    pub surroundings_k: f64,
}

/// Actual original-source coefficient and physical heat rate at one wall patch.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceBodyRadiationSample {
    pub body_id: u64,
    pub face: FreeSurfaceBodyRadiationFace,
    pub solid_temperature_k: f64,
    pub emissivity: f64,
    pub emissivity_source: ResolvedScalarProperty,
    pub conductance_w_k: f64,
    /// Positive into the finite solid; never debited from the finite gas [W].
    pub heat_into_solid_w: f64,
}

/// External radiant supply, separate from internal wall/liquid and gas exchange.
#[derive(Clone, Debug, PartialEq)]
pub struct FreeSurfaceBodyRadiationExchange {
    pub heat_to_solids_j: f64,
    pub last_applied: Vec<FreeSurfaceBodyRadiationSample>,
    pub thermal_endpoint: Vec<FreeSurfaceBodyRadiationSample>,
}

#[derive(Clone, Default)]
pub(in super::super) struct BodyRadiation {
    pub(in super::super) faces: Vec<(usize, FreeSurfaceBodyRadiationFace)>,
    pub(in super::super) last: Vec<FreeSurfaceBodyRadiationSample>,
    pub(in super::super) endpoint: Vec<FreeSurfaceBodyRadiationSample>,
    pub(in super::super) heat_j: f64,
    pub(in super::super) absolute_heat_j: f64,
}

impl MaterialFreeSurfaceEnthalpy3<'_, '_> {
    pub(in super::super) fn admit_body_radiation(
        &self,
        bodies: &[FreeSurfaceThermalBody<'_>],
        dx: f64,
        gate: &CancelGate,
    ) -> Result<BodyRadiation, Error> {
        let mut timeline = BodyRadiation::default();
        let mut faces = BTreeSet::new();
        for (index, body) in bodies.iter().enumerate() {
            for &face in body.radiation_faces() {
                poll(gate)?;
                let [Ok(x), Ok(y), Ok(z)] = face.gas.map(i64::try_from) else {
                    return Err(Error::Invalid(
                        "radiating wall gas coordinate is unrepresentable",
                    ));
                };
                let gas = [x, y, z];
                let distance = (0..3).fold(0u64, |sum, axis| {
                    sum.saturating_add(face.wall[axis].abs_diff(gas[axis]))
                });
                if distance != 1
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
                    || !face.surroundings_k.is_finite()
                    || face.surroundings_k <= 0.0
                {
                    return Err(Error::Invalid(
                        "solid radiation needs a unique exposed wall face and finite area/positive kelvin bath",
                    ));
                }
                timeline.faces.push((index, face));
            }
        }
        timeline
            .faces
            .sort_by_key(|(index, face)| (bodies[*index].id(), face.wall, face.gas));
        Ok(timeline)
    }
}

impl Conduction<'_, '_> {
    pub(in super::super) fn sample_body_radiation(
        &self,
        thermal: &DepositionSimulation<'_>,
        check: &mut dyn FnMut() -> Result<(), DepositionError>,
    ) -> Result<Vec<FreeSurfaceBodyRadiationSample>, Error> {
        let mut samples = Vec::with_capacity(self.body_radiation.faces.len());
        let mut covered = BTreeMap::new();
        for &(body_index, face) in &self.body_radiation.faces {
            check()?;
            let body = &self.bodies[body_index];
            let id = self.body_ids[body_index];
            let area = covered.entry(id).or_insert(0.0);
            *area = environment::finite(*area + face.area_m2)?;
            // Contact coverage is subtracted by the original thermal owner.
            // Gas-film and radiative exchange may share the uncovered patch.
            if *area > thermal.exposed_area_m2(id)? {
                return Err(Error::Invalid(
                    "solid radiation exceeds uncovered solid area",
                ));
            }
            let material = body.resolve(thermal.cells()[&id].specific_h_j_kg)?;
            let emissivity_source =
                material
                    .property(SURFACE_EMISSIVITY_PROPERTY)
                    .ok_or(Error::Invalid(
                        "radiating wall lost its original emissivity",
                    ))?;
            let emissivity = emissivity_source.value_si();
            let solid_temperature_k = thermal.temperature_k(id)?;
            let conductance_w_k = radiant_conductance(
                emissivity,
                face.area_m2,
                solid_temperature_k,
                face.surroundings_k,
            )?;
            samples.push(FreeSurfaceBodyRadiationSample {
                body_id: body.id(),
                face,
                solid_temperature_k,
                emissivity,
                emissivity_source: emissivity_source.clone(),
                conductance_w_k,
                heat_into_solid_w: environment::finite(
                    conductance_w_k * (face.surroundings_k - solid_temperature_k),
                )?,
            });
        }
        Ok(samples)
    }
}
