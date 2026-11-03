//! Curved closest points become the unchanged native work-conjugate trace.

use std::collections::BTreeMap;

use super::*;
use crate::inelastic_impact::multiple::{FiniteSolidContactCounterpart, RigidImpactPoint};
use fs_exec::Cx;
use fs_mbd::RigidBodyState;
use fs_query::{TriangleFeature, triangle::quadratic::QuadraticTrianglePoint};
use fs_solid::linear3::tet10::{tet10_point, tet10_shapes};

impl QuadraticSolidContactSurface {
    /// Refit all actual six-node patches once, then query bounded sphere geometry
    /// in deterministic ID order. `positions` follows the complete original node
    /// ordering. Explicit coordinate images do not assert an accepted native
    /// material state; the native helpers bind that state and reference identity.
    pub fn sphere_separations_at_positions(
        &self,
        positions: &[[f64; 3]],
        rigid_bodies: &[RigidBodyState],
        spheres: &[FiniteSolidSphere],
        controls: QuadraticSolidSurfaceControls,
        cx: &Cx<'_>,
    ) -> Result<QuadraticSolidSurfaceQuery, QuadraticSolidSurfaceError> {
        checkpoint(cx)?;
        let admission = controls.surface;
        if positions.len() != self.reference_nodes_m.len()
            || spheres.len() > admission.maximum_spheres
            || !admission.activation_tolerance_m.is_finite()
            || admission.activation_tolerance_m < 0.0
            || !admission.minimum_scaled_jacobian.is_finite()
            || admission.minimum_scaled_jacobian <= 0.0
            || !controls.maximum_witness_span_m.is_finite()
            || controls.maximum_witness_span_m <= 0.0
            || !controls.query.distance_tolerance.is_finite()
            || controls.query.distance_tolerance <= 0.0
            || !controls.query.witness_tolerance.is_finite()
            || controls.query.witness_tolerance <= 0.0
            || controls.query.maximum_cells == 0
            || controls.query.maximum_active_cells == 0
        {
            return Err(invalid("matching bounded curved geometry and positive finite tolerances required"));
        }
        let mut query_controls = controls.query;
        // Leave room for the original native trace's independent floating-point
        // reconstruction. Its actual error is checked below, never assumed zero.
        query_controls.witness_tolerance = query_controls
            .witness_tolerance
            .min(0.5 * controls.maximum_witness_span_m);
        if query_controls.witness_tolerance <= 0.0 {
            return Err(invalid("quadratic witness tolerance is not representable"));
        }
        for (index, position) in positions.iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            vector(*position)?;
        }
        let surface = self.triangle_surface.refit(positions, cx)?;
        let mut order: Vec<_> = (0..spheres.len()).collect();
        order.sort_unstable_by_key(|&index| spheres[index].id);
        let mut report = QuadraticSolidSurfaceQuery {
            closest: Vec::with_capacity(spheres.len()),
            contacts: Vec::new(),
            unresolved_activation: None,
            visited_nodes: 0,
            tested_triangles: 0,
            examined_cells: 0,
            candidate_triangles: 0,
        };
        let mut candidates: BTreeMap<u64, QuadraticSolidSphereWitness> = BTreeMap::new();
        for (slot, &index) in order.iter().enumerate() {
            checkpoint(cx)?;
            let sphere = spheres[index];
            if sphere.rigid_body >= rigid_bodies.len()
                || !sphere.center_body_m.is_finite()
                || !sphere.radius_m.is_finite()
                || sphere.radius_m <= 0.0
                || !sphere.friction_coefficient.is_finite()
                || sphere.friction_coefficient < 0.0
                || (slot > 0 && spheres[order[slot - 1]].id == sphere.id)
            {
                return Err(invalid("distinct sphere IDs, valid rigid body, positive radius and nonnegative friction required"));
            }
            sphere
                .id
                .checked_mul(self.feature_count)
                .and_then(|base| base.checked_add(self.feature_count - 1))
                .ok_or(FiniteSolidSurfaceError::ContactIdOverflow {
                    sphere_id: sphere.id,
                })?;
            let rigid = rigid_bodies[sphere.rigid_body];
            let center = rigid.pose().point_world_from_body(sphere.center_body_m)?;
            let within = finite(
                finite(sphere.radius_m + admission.activation_tolerance_m)?.next_up(),
            )?;
            query_controls.maximum_cells = controls
                .query
                .maximum_cells
                .saturating_sub(report.examined_cells);
            if query_controls.maximum_cells == 0 {
                return Err(QuadraticQueryError::CellBudgetExhausted {
                    maximum_cells: controls.query.maximum_cells,
                }
                .into());
            }
            query_controls.maximum_active_cells = controls
                .query
                .maximum_active_cells
                .min(query_controls.maximum_cells);
            let query = surface.point_query(
                array(center),
                within,
                admission.maximum_query_visits.saturating_sub(
                    report.visited_nodes.saturating_add(report.tested_triangles),
                ),
                admission
                    .maximum_candidates
                    .saturating_sub(report.candidate_triangles),
                query_controls,
                cx,
            )?;
            report.visited_nodes += query.visited_nodes;
            report.tested_triangles += query.tested_triangles;
            report.examined_cells += query.examined_cells;
            report.candidate_triangles += query.nearby.len();
            if report.visited_nodes.saturating_add(report.tested_triangles)
                > admission.maximum_query_visits
            {
                return Err(FiniteSolidSurfaceError::WorkBudgetExhausted {
                    maximum_query_visits: admission.maximum_query_visits,
                }
                .into());
            }
            // Enclose the outward side at EVERY patch still eligible as globally
            // nearest. A sampled normal on one approximate minimizer cannot
            // establish exterior containment in a thin or concave curved shell.
            if query.nearest_oriented_side_bounds[0] <= 0.0 {
                return Err(FiniteSolidSurfaceError::UnresolvedExterior {
                    sphere_id: sphere.id,
                }
                .into());
            }
            let closest = self.witness(
                positions,
                sphere,
                rigid,
                center,
                query.nearest_triangle,
                &query.nearest,
                [query.distance_lower, query.distance_upper],
                controls,
            )?;
            report.closest.push(closest);
            for (face, point) in query.nearby {
                checkpoint(cx)?;
                // A patch wholly facing away cannot supply an exterior contact.
                // A sign-straddling possible active patch remains unresolved.
                if point.oriented_side_bounds[1] < 0.0 {
                    continue;
                }
                if point.oriented_side_bounds[0] <= 0.0 {
                    return Err(FiniteSolidSurfaceError::UnresolvedExterior {
                        sphere_id: sphere.id,
                    }
                    .into());
                }
                let witness = self.witness(
                    positions,
                    sphere,
                    rigid,
                    center,
                    face,
                    &point,
                    [point.distance_lower, point.distance_upper],
                    controls,
                )?;
                if witness.gap_lower_m > admission.activation_tolerance_m {
                    continue;
                }
                let candidate_count = candidates.len();
                match candidates.entry(witness.contact.id) {
                    std::collections::btree_map::Entry::Vacant(entry) => {
                        if candidate_count == admission.maximum_contacts {
                            return Err(FiniteSolidSurfaceError::ContactBudgetExhausted {
                                maximum_contacts: admission.maximum_contacts,
                            }
                            .into());
                        }
                        entry.insert(witness);
                    }
                    std::collections::btree_map::Entry::Occupied(mut entry) => {
                        let previous = entry.get_mut();
                        let lower = previous.gap_lower_m.min(witness.gap_lower_m);
                        let upper = previous.gap_upper_m.max(witness.gap_upper_m);
                        if witness.face < previous.face {
                            *previous = witness;
                        }
                        previous.gap_lower_m = lower;
                        previous.gap_upper_m = upper;
                    }
                }
            }
        }
        for witness in candidates.into_values() {
            checkpoint(cx)?;
            if witness.gap_upper_m <= admission.activation_tolerance_m
                && witness.gap_m <= admission.activation_tolerance_m
            {
                report.contacts.push(witness);
            } else if report.unresolved_activation.is_none() {
                report.unresolved_activation = Some(witness);
            }
        }
        checkpoint(cx)?;
        Ok(report)
    }

    #[allow(clippy::too_many_arguments)]
    fn witness(
        &self,
        positions: &[[f64; 3]],
        sphere: FiniteSolidSphere,
        rigid: RigidBodyState,
        center: Vec3,
        face: usize,
        point: &QuadraticTrianglePoint,
        distance_bounds: [f64; 2],
        controls: QuadraticSolidSurfaceControls,
    ) -> Result<QuadraticSolidSphereWitness, QuadraticSolidSurfaceError> {
        let material = &self.material_faces[face];
        let mut barycentric = [0.0; 4];
        for (local, weight) in material.local_corners.into_iter().zip(point.barycentric) {
            barycentric[local] = weight;
        }
        let sum: f64 = barycentric.iter().sum();
        if barycentric
            .iter()
            .any(|weight| !weight.is_finite() || !(0.0..=1.0).contains(weight))
            || !sum.is_finite()
            || (sum - 1.0).abs() > 64.0 * f64::EPSILON
        {
            return Err(invalid("quadratic query returned inadmissible native barycentrics"));
        }
        let solid_point = Tet10ImpulsePoint {
            tetrahedron: material.tetrahedron,
            barycentric,
        };
        // Store original natural coordinates, then perform exactly the same
        // one-time normalization and signed shape evaluation as native inertia.
        let lambda = barycentric.map(|weight| weight / sum);
        let nodes = self.tetrahedra[material.tetrahedron];
        for coordinates in [&self.reference_nodes_m[..], positions] {
            tet10_point(
                nodes.map(|node| coordinates[node]),
                lambda,
                controls.surface.minimum_scaled_jacobian,
                material.tetrahedron,
                face,
            )?;
        }
        let mut position = [0.0; 3];
        for (node, weight) in nodes.into_iter().zip(tet10_shapes(lambda).0) {
            for axis in 0..3 {
                position[axis] = finite(weight.mul_add(positions[node][axis], position[axis]))?;
            }
        }
        let solid_point_world_m = vec(vector(position)?);
        let trace_error_m = magnitude(solid_point_world_m.sub(vec(point.point)))?;
        let witness_error = finite((point.witness_error_bound + trace_error_m).next_up())?;
        if witness_error > controls.maximum_witness_span_m {
            return Err(QuadraticSolidSurfaceError::AmbiguousWitness {
                sphere_id: sphere.id,
                face,
                error_bound_m: witness_error,
                maximum_error_m: controls.maximum_witness_span_m,
            });
        }
        let radial = center.sub(solid_point_world_m);
        let distance = magnitude(radial)?;
        if distance == 0.0 {
            return Err(FiniteSolidSurfaceError::UnresolvedExterior {
                sphere_id: sphere.id,
            }
            .into());
        }
        let normal_world = Vec3::new(radial.x / distance, radial.y / distance, radial.z / distance);
        let arm_body = sphere.center_body_m.sub(
            rigid
                .pose()
                .orientation()
                .rotate_world_to_body(normal_world.scale(sphere.radius_m)),
        );
        let rigid_point_world_m = rigid.pose().point_world_from_body(arm_body)?;
        let gap_m = finite(normal_world.dot(rigid_point_world_m.sub(solid_point_world_m)))?;
        // The rigid material arm and native trace are the actual mechanical
        // witnesses. Preserve their rounding differences in the geometric gap
        // interval instead of moving either point to manufacture zero separation.
        let rigid_error = finite((gap_m - finite(distance - sphere.radius_m)?).abs())?;
        let reconstruction_error = finite((trace_error_m + rigid_error).next_up())?;
        let gap_lower_m = finite(
            ((distance_bounds[0] - sphere.radius_m).next_down() - reconstruction_error).next_down(),
        )?
        .min(gap_m);
        let gap_upper_m = finite(
            ((distance_bounds[1] - sphere.radius_m).next_up() + reconstruction_error).next_up(),
        )?
        .max(gap_m);
        let feature = match point.feature {
            TriangleFeature::Vertex(local) => usize::try_from(self.faces()[face][local])
                .map_err(|_| invalid("quadratic feature node index is not representable"))?,
            TriangleFeature::Edge(a, b) => {
                let local_edge = match (a, b) {
                    (0, 1) => 0,
                    (0, 2) => 1,
                    (1, 2) => 2,
                    _ => return Err(invalid("unknown quadratic triangle edge")),
                };
                self.reference_nodes_m.len() + material.edge_ids[local_edge]
            }
            TriangleFeature::Face => self.reference_nodes_m.len() + self.edge_count + face,
        };
        let feature = u64::try_from(feature).map_err(|_| {
            FiniteSolidSurfaceError::ContactIdOverflow {
                sphere_id: sphere.id,
            }
        })?;
        let id = sphere
            .id
            .checked_mul(self.feature_count)
            .and_then(|base| base.checked_add(feature))
            .ok_or(FiniteSolidSurfaceError::ContactIdOverflow {
                sphere_id: sphere.id,
            })?;
        let components = [normal_world.x.abs(), normal_world.y.abs(), normal_world.z.abs()];
        let axis = if components[0] <= components[1] && components[0] <= components[2] {
            Vec3::new(1.0, 0.0, 0.0)
        } else if components[1] <= components[2] {
            Vec3::new(0.0, 1.0, 0.0)
        } else {
            Vec3::new(0.0, 0.0, 1.0)
        };
        let tangent = normal_world.cross(axis);
        let tangent_world = tangent.scale(magnitude(tangent)?.recip());
        Ok(QuadraticSolidSphereWitness {
            sphere_id: sphere.id,
            face,
            rigid_point_world_m,
            solid_point_world_m,
            gap_m,
            gap_lower_m,
            gap_upper_m,
            contact: InelasticMixedContact {
                id,
                point: RigidImpactPoint {
                    body: sphere.rigid_body,
                    arm_body,
                },
                counterpart: FiniteSolidContactCounterpart::Solid(solid_point),
                normal_world,
                tangent_world,
                activation_tolerance_m: controls.surface.activation_tolerance_m,
                friction_coefficient: sphere.friction_coefficient,
                normal_law: sphere.normal_law,
            },
        })
    }
}

fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn vec(value: [f64; 3]) -> Vec3 {
    Vec3::new(value[0], value[1], value[2])
}
fn magnitude(value: Vec3) -> Result<f64, QuadraticSolidSurfaceError> {
    finite(value.x.hypot(value.y).hypot(value.z))
}
