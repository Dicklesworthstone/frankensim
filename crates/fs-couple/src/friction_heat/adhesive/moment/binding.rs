//! Original caloric material, ordered real face and native marker registration.

use super::{AdhesiveMomentHeat, AdhesiveMomentHeatError};
use fs_math::det;
use fs_mbd::{MassProperties, RigidBodyState, Vec3};
use fs_thermal::enthalpy::DepositionSimulation;

impl AdhesiveMomentHeat<'_> {
    pub(super) fn validate_receivers<E>(
        &self,
        thermal: &DepositionSimulation<'_>,
    ) -> Result<(), AdhesiveMomentHeatError<E>> {
        self.route.validate(thermal)?;
        let elastic = self.native.normal_source().normal_state().elastic();
        let identities = [
            elastic.surface_a_state().resolved().card_identity(),
            elastic.surface_b_state().resolved().card_identity(),
        ];
        let area = self.native.normal_response().contact_area_m2;
        if self.body_ids[0] == self.body_ids[1]
            || (area - self.route.contact_area_m2()).abs()
                > 64.0 * f64::EPSILON * area.max(self.route.contact_area_m2())
        {
            return Err(AdhesiveMomentHeatError::Invalid(
                "distinct original bodies and actual attached JKR area required",
            ));
        }
        for (side, id) in self.route.cells().into_iter().enumerate() {
            let cell = thermal
                .cells()
                .get(&id)
                .ok_or(AdhesiveMomentHeatError::Invalid(
                    "missing original finite caloric receiver",
                ))?;
            if cell.body_id != self.body_ids[side]
                || self.caloric[side].material_card_identity() != identities[side]
                || !thermal.uses_material_law(cell.material_id, self.caloric[side])
            {
                return Err(AdhesiveMomentHeatError::Invalid(
                    "original ordered native bodies and registered source caloric laws required",
                ));
            }
        }
        Ok(())
    }

    pub(super) fn validate_geometry<E>(
        &self,
        bodies: &[RigidBodyState],
        properties: &[MassProperties],
        tolerance: f64,
    ) -> Result<(), AdhesiveMomentHeatError<E>> {
        let geometry = self.route.geometry;
        let centre = [geometry[0].start_m, geometry[0].end_m]
            .into_iter()
            .find(|p| *p == geometry[1].start_m || *p == geometry[1].end_m)
            .ok_or(AdhesiveMomentHeatError::Invalid(
                "original heat face no longer adjoins",
            ))?;
        let binding = self.native.binding();
        let indices = [binding.body_a, binding.body_b];
        let arms = [binding.arm_a_body, binding.arm_b_body];
        let axis = self.native.frame().normal();
        let normal = Vec3::new(axis[0], axis[1], axis[2]);
        for side in 0..2 {
            let body = bodies
                .get(indices[side])
                .ok_or(AdhesiveMomentHeatError::Invalid("missing native body"))?;
            let property =
                properties
                    .get(indices[side])
                    .ok_or(AdhesiveMomentHeatError::Invalid(
                        "missing native mass properties",
                    ))?;
            let point = body
                .point_kinematics(*property, arms[side])
                .map_err(|_| AdhesiveMomentHeatError::Invalid("native marker kinematics refused"))?
                .point_world;
            let g = geometry[side];
            let axis = Vec3::new(
                g.end_m[0] - g.start_m[0],
                g.end_m[1] - g.start_m[1],
                g.end_m[2] - g.start_m[2],
            );
            let interior = if g.start_m == centre {
                g.end_m
            } else {
                g.start_m
            };
            let inward = Vec3::new(
                0.5 * (interior[0] - centre[0]),
                0.5 * (interior[1] - centre[1]),
                0.5 * (interior[2] - centre[2]),
            );
            let depth = inward.dot(normal) * if side == 0 { -1.0 } else { 1.0 };
            if magnitude(point.sub(Vec3::new(centre[0], centre[1], centre[2]))) > tolerance
                || magnitude(axis.cross(normal)) > tolerance
                || !depth.is_finite()
                || depth <= tolerance
            {
                return Err(AdhesiveMomentHeatError::Invalid(
                    "native marker and ordered A-to-B normal must coincide with original heat face",
                ));
            }
        }
        Ok(())
    }
}
fn magnitude(v: Vec3) -> f64 {
    det::hypot(det::hypot(v.x, v.y), v.z)
}
