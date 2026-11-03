//! Whole-path sphere clearance against a prepared quadratic solid boundary.
//!
//! Endpoint poses locate each sphere centre; the declared path then interpolates
//! those centres and every solid node affinely. This does not prescribe rigid
//! rotational interpolation, integrate a solid trajectory, locate an impact
//! root, or certify volume injectivity along the path. Initial exterior and
//! actual-witness gap checks reuse the original curved geometry owner. The
//! complete path proof can resolve an initially uncertain lower gap at touch.

use super::{
    FiniteSolidSphere, FiniteSolidSurfaceError, QuadraticSolidContactSurface,
    QuadraticSolidSurfaceControls, QuadraticSolidSurfaceError, checkpoint, finite, invalid,
};
use fs_exec::Cx;
use fs_mbd::RigidBodyState;
use fs_query::triangle::{
    TriangleQueryError,
    quadratic::{
        QuadraticQueryError,
        sweep::{QuadraticSweepControls, linear_sphere_quadratic_triangle},
    },
    sweep::LinearSweepVerdict,
};

/// Initial geometry admission and aggregate whole-path work limits.
#[derive(Clone, Copy, Debug, Default)]
pub struct QuadraticSolidPathControls {
    /// Reuses sphere identity, local Jacobian and exterior witness admission.
    /// Query visits and candidate limits include both initial geometry and
    /// subsequent swept broad-phase traversal across all spheres. The finite
    /// activation tolerance also pads swept queries: a culled face has at least
    /// that much clearance. Activation-band ambiguity alone does not refuse a
    /// path with a nonnegative initial actual gap. An uncertain initial lower
    /// gap must still be resolved by the complete path proof.
    pub initial: QuadraticSolidSurfaceControls,
    /// Window and total distance-cell caps are shared across all face paths.
    /// The total distance-cell cap also includes initial geometry cells. The
    /// nested distance controls still cap each individual global distance solve.
    pub sweep: QuadraticSweepControls,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum QuadraticSolidPathVerdict {
    /// Every sphere is nonpenetrating throughout [0,1]. This rigorous lower
    /// bound is capped by the authored broad-phase padding, including when no
    /// pairs remain. A zero lower bound does not assert positive separation.
    Clear { minimum_clearance_m: f64 },
    /// Earliest unresolved window across all sphere/face paths. This is neither
    /// an asserted collision nor an impact-time root. Equal start fractions
    /// select the lowest sphere ID, then canonical face index.
    PossibleContact {
        sphere_id: u64,
        face: usize,
        fraction_start: f64,
        fraction_end: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticSolidPathReport {
    pub verdict: QuadraticSolidPathVerdict,
    /// Initial and swept BVH visits, excluding the mesh-bounded linear refits.
    pub visited_nodes: usize,
    /// Initial narrow-phase patches plus swept broad-phase patch box tests.
    pub tested_triangles: usize,
    /// Initial near candidates plus whole-path candidate pairs, before any
    /// geometric feature deduplication. Each sphere/face pair is counted once
    /// in the swept traversal, even when initial geometry also examined it.
    pub candidate_triangles: usize,
    pub tested_paths: usize,
    pub examined_windows: usize,
    /// Initial global distance work plus all narrow-phase path distance work.
    pub examined_distance_cells: usize,
}

impl QuadraticSolidContactSurface {
    /// Certify the complete declared affine sphere-centre and six-node paths.
    ///
    /// Node order and rigid-body indices agree at both endpoints. Sphere centres
    /// are the endpoint poses applied to each original body-frame centre; the
    /// straight segment between those two points is the declared centre path.
    /// Rigid velocities, spins and the solid's constitutive history are neither
    /// advanced nor modified. Volume validity remains a separate admission.
    ///
    /// The initial query must resolve exterior containment and return no
    /// negative actual witness gap. An initially uncertain lower gap may proceed
    /// to the complete path proof, permitting exact touching/receding paths
    /// without an epsilon allowance. Initially contained spheres, negative
    /// actual gaps and ambiguous exterior witnesses refuse. Unresolved activation
    /// alone is irrelevant to path clearance.
    #[allow(clippy::too_many_arguments)]
    pub fn affine_sphere_paths_at_positions(
        &self,
        positions_start: &[[f64; 3]],
        positions_end: &[[f64; 3]],
        rigid_start: &[RigidBodyState],
        rigid_end: &[RigidBodyState],
        spheres: &[FiniteSolidSphere],
        controls: QuadraticSolidPathControls,
        cx: &Cx<'_>,
    ) -> Result<QuadraticSolidPathReport, QuadraticSolidSurfaceError> {
        checkpoint(cx)?;
        controls.sweep.validate()?;
        // Validate authored initial query limits before clipping them to the
        // shared budget, including when there are no spheres to query.
        let mut initial_limits = controls.sweep;
        initial_limits.distance = controls.initial.query;
        initial_limits.validate()?;
        if positions_start.len() != self.reference_nodes_m.len()
            || positions_end.len() != self.reference_nodes_m.len()
            || rigid_start.len() != rigid_end.len()
        {
            return Err(invalid(
                "matching complete node and rigid-body endpoint images required",
            ));
        }
        let mut initial_controls = controls.initial;
        initial_controls.query.maximum_cells = initial_controls
            .query
            .maximum_cells
            .min(controls.sweep.maximum_distance_cells);
        initial_controls.query.maximum_active_cells = initial_controls
            .query
            .maximum_active_cells
            .min(initial_controls.query.maximum_cells);
        let initial = self.sphere_separations_at_positions(
            positions_start,
            rigid_start,
            spheres,
            initial_controls,
            cx,
        )?;
        for witness in &initial.closest {
            checkpoint(cx)?;
            if witness.gap_m < 0.0 {
                return Err(FiniteSolidSurfaceError::Penetration {
                    sphere_id: witness.sphere_id,
                    gap_m: witness.gap_m,
                }
                .into());
            }
            // A straddling interval is not a penetration allowance. Every face
            // relevant at t=0 remains in the complete swept broad phase; only
            // the ensuing whole-path proof can admit the path as Clear.
        }
        let padding = controls.initial.surface.activation_tolerance_m;
        let mut report = QuadraticSolidPathReport {
            verdict: QuadraticSolidPathVerdict::Clear {
                minimum_clearance_m: padding,
            },
            visited_nodes: initial.visited_nodes,
            tested_triangles: initial.tested_triangles,
            candidate_triangles: initial.candidate_triangles,
            tested_paths: 0,
            examined_windows: 0,
            examined_distance_cells: initial.examined_cells,
        };
        drop(initial);
        // This validates every end coordinate even for an empty sphere set.
        // Endpoint Bernstein hull unions cover every intermediate nodal image.
        let swept = self
            .triangle_surface
            .refit_swept(positions_start, positions_end, cx)?;
        let mut order: Vec<_> = (0..spheres.len()).collect();
        order.sort_unstable_by_key(|&index| spheres[index].id);
        let mut minimum = padding;
        let mut possible: Option<(f64, f64, u64, usize)> = None;
        for index in order {
            checkpoint(cx)?;
            // Shared initial admission has already checked IDs, body indices,
            // body-frame centres, original radii and descriptor material data.
            let sphere = spheres[index];
            let a = rigid_start[sphere.rigid_body]
                .pose()
                .point_world_from_body(sphere.center_body_m)?;
            let b = rigid_end[sphere.rigid_body]
                .pose()
                .point_world_from_body(sphere.center_body_m)?;
            let center_start = [a.x, a.y, a.z];
            let center_end = [b.x, b.y, b.z];
            let search_radius = if padding == 0.0 {
                sphere.radius_m
            } else {
                // Original geometry owner uses this outward bound for the
                // single finite radius-plus-band addition as well.
                finite(finite(sphere.radius_m + padding)?.next_up())?
            };
            let candidates = swept.sphere_candidates(
                center_start,
                center_end,
                search_radius,
                controls
                    .initial
                    .surface
                    .maximum_query_visits
                    .saturating_sub(report.visited_nodes.saturating_add(report.tested_triangles)),
                controls
                    .initial
                    .surface
                    .maximum_candidates
                    .saturating_sub(report.candidate_triangles),
                cx,
            )?;
            report.visited_nodes += candidates.visited_nodes;
            report.tested_triangles += candidates.tested_triangles;
            report.candidate_triangles += candidates.triangles.len();
            for face in candidates.triangles {
                checkpoint(cx)?;
                let remaining_windows = controls
                    .sweep
                    .time
                    .maximum_windows
                    .saturating_sub(report.examined_windows);
                if remaining_windows == 0 {
                    return Err(QuadraticQueryError::Geometry(
                        TriangleQueryError::WorkBudgetExhausted {
                            maximum_visits: controls.sweep.time.maximum_windows,
                        },
                    )
                    .into());
                }
                let remaining_cells = controls
                    .sweep
                    .maximum_distance_cells
                    .saturating_sub(report.examined_distance_cells);
                if remaining_cells == 0 {
                    return Err(QuadraticQueryError::CellBudgetExhausted {
                        maximum_cells: controls.sweep.maximum_distance_cells,
                    }
                    .into());
                }
                let mut limits = controls.sweep;
                limits.time.maximum_windows = remaining_windows;
                limits.maximum_distance_cells = remaining_cells;
                let nodes = self.triangle_surface.triangles()[face];
                let path = linear_sphere_quadratic_triangle(
                    center_start,
                    center_end,
                    sphere.radius_m,
                    nodes.map(|node| positions_start[node as usize]),
                    nodes.map(|node| positions_end[node as usize]),
                    limits,
                    cx,
                )?;
                report.tested_paths += 1;
                report.examined_windows += path.path.examined_windows;
                report.examined_distance_cells += path.examined_distance_cells;
                match path.path.verdict {
                    LinearSweepVerdict::Clear { minimum_clearance } => {
                        minimum = minimum.min(minimum_clearance);
                    }
                    LinearSweepVerdict::PossibleContact {
                        fraction_start,
                        fraction_end,
                    } => {
                        // Canonical traversal breaks equal-start ties. All
                        // remaining faces still run so an earlier window wins.
                        if possible.is_none_or(|earlier| fraction_start < earlier.0) {
                            possible = Some((fraction_start, fraction_end, sphere.id, face));
                        }
                    }
                }
            }
        }
        report.verdict = match possible {
            Some((fraction_start, fraction_end, sphere_id, face)) => {
                QuadraticSolidPathVerdict::PossibleContact {
                    sphere_id,
                    face,
                    fraction_start,
                    fraction_end,
                }
            }
            None => QuadraticSolidPathVerdict::Clear {
                minimum_clearance_m: minimum,
            },
        };
        checkpoint(cx)?;
        Ok(report)
    }
}

#[cfg(test)]
mod tests;
