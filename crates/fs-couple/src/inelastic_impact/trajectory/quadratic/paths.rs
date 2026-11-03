//! Remaining rigid-pair paths share the authored boundary work allowance.

use super::{
    FiniteSolidSphere, FiniteSolidSphereTrajectoryError, MixedViscoSpherePathReport,
    MixedViscoSphereTrajectoryError, QuadraticSolidPathControls, QuadraticSolidPathReport,
    QuadraticSolidVolumePathReport, poll,
};
use fs_exec::Cx;
use fs_mbd::{RigidBodyState, Vec3};
use fs_query::triangle::sweep::{LinearSweepVerdict, linear_spheres};

#[allow(clippy::too_many_arguments)]
pub(super) fn admit_pairs(
    volume: QuadraticSolidVolumePathReport,
    boundary: QuadraticSolidPathReport,
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    spheres: &[FiniteSolidSphere],
    controls: QuadraticSolidPathControls,
    cx: &Cx<'_>,
) -> Result<MixedViscoSpherePathReport, MixedViscoSphereTrajectoryError> {
    use MixedViscoSphereTrajectoryError as E;
    let mut report = MixedViscoSpherePathReport {
        volume,
        boundary,
        sphere_pairs_checked: 0,
        sphere_pair_windows: 0,
        minimum_pair_clearance_m: None,
    };
    for (index, a) in spheres.iter().enumerate() {
        for b in &spheres[index + 1..] {
            poll(cx)?;
            let visits = boundary
                .visited_nodes
                .saturating_add(boundary.tested_triangles)
                .saturating_add(report.sphere_pairs_checked);
            let windows = boundary
                .examined_windows
                .saturating_add(report.sphere_pair_windows);
            if visits >= controls.initial.surface.maximum_query_visits
                || windows >= controls.sweep.time.maximum_windows
            {
                return Err(E::Path(FiniteSolidSphereTrajectoryError::PathBudget));
            }
            let mut time = controls.sweep.time;
            time.maximum_windows -= windows;
            let pair = linear_spheres(
                array(before[a.rigid_body].pose().position_world()),
                array(after[a.rigid_body].pose().position_world()),
                a.radius_m,
                array(before[b.rigid_body].pose().position_world()),
                array(after[b.rigid_body].pose().position_world()),
                b.radius_m,
                time,
                cx,
            )
            .map_err(|e| E::Path(FiniteSolidSphereTrajectoryError::Query(e)))?;
            report.sphere_pairs_checked += 1;
            report.sphere_pair_windows += pair.examined_windows;
            match pair.verdict {
                LinearSweepVerdict::Clear { minimum_clearance } => {
                    report.minimum_pair_clearance_m = Some(
                        report
                            .minimum_pair_clearance_m
                            .map_or(minimum_clearance, |previous| {
                                previous.min(minimum_clearance)
                            }),
                    );
                }
                LinearSweepVerdict::PossibleContact { fraction_start, .. } => {
                    return Err(E::Path(
                        FiniteSolidSphereTrajectoryError::UnsupportedSpherePair {
                            first: a.id,
                            second: b.id,
                            fraction: fraction_start,
                        },
                    ));
                }
            }
        }
    }
    Ok(report)
}

fn array(v: Vec3) -> [f64; 3] {
    [v.x, v.y, v.z]
}
