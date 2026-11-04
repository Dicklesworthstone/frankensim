//! Remaining rigid-pair paths share the authored boundary work allowance.

use super::{
    FiniteSolidSphere, FiniteSolidSphereTrajectoryError, MixedViscoSpherePathRefusal,
    MixedViscoSpherePathReport, MixedViscoSphereTrajectoryError, QuadraticSolidContactSurface,
    QuadraticSolidPathControls, QuadraticSolidPathReport, QuadraticSolidPathVerdict,
    QuadraticSolidVolumePathReport, QuadraticSphereScenePathVerdict, TrajectoryScene, poll,
};
use fs_exec::Cx;
use fs_mbd::{RigidBodyState, Vec3};
use fs_query::triangle::sweep::{LinearSweepVerdict, linear_spheres};

/// Both entry points share the same native motion and physical accounts. A
/// prepared scene additionally authorizes its explicitly declared pair laws.
#[allow(clippy::too_many_arguments)]
pub(super) fn admit_path(
    scene: TrajectoryScene<'_>,
    surface: &QuadraticSolidContactSurface,
    volume: QuadraticSolidVolumePathReport,
    start: &[[f64; 3]],
    end: &[[f64; 3]],
    before: &[RigidBodyState],
    after: &[RigidBodyState],
    spheres: &[FiniteSolidSphere],
    controls: QuadraticSolidPathControls,
    cx: &Cx<'_>,
) -> Result<
    Result<MixedViscoSpherePathReport, MixedViscoSpherePathRefusal>,
    MixedViscoSphereTrajectoryError,
> {
    use MixedViscoSphereTrajectoryError as E;
    if let TrajectoryScene::Complete(scene) = scene {
        let path = scene
            .affine_quadratic_sphere_paths_at_positions(
                surface, start, end, before, after, controls, cx,
            )
            .map_err(E::Scene)?;
        if matches!(
            path.verdict,
            QuadraticSphereScenePathVerdict::PossibleContact { .. }
        ) {
            return Ok(Err(MixedViscoSpherePathRefusal::ScenePath(path)));
        }
        return Ok(Ok(MixedViscoSpherePathReport {
            contact_keys: Vec::new(),
            volume,
            boundary: path.boundary,
            sphere_pair_comparisons: path.sphere_pair_comparisons,
            sphere_pair_candidates: path.sphere_pair_candidates,
            sphere_pairs_checked: path.sphere_pairs_checked,
            sphere_pair_windows: path.sphere_pair_windows,
            minimum_pair_clearance_m: path.minimum_pair_clearance_m,
        }));
    }
    let boundary = surface
        .affine_sphere_paths_at_positions(start, end, before, after, spheres, controls, cx)
        .map_err(E::Surface)?;
    if matches!(
        boundary.verdict,
        QuadraticSolidPathVerdict::PossibleContact { .. }
    ) {
        return Ok(Err(MixedViscoSpherePathRefusal::Boundary(boundary)));
    }
    admit_pairs(volume, boundary, before, after, spheres, controls, cx).map(Ok)
}

#[allow(clippy::too_many_arguments)]
fn admit_pairs(
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
        contact_keys: Vec::new(),
        volume,
        boundary,
        sphere_pair_comparisons: 0,
        sphere_pair_candidates: 0,
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
            report.sphere_pair_comparisons += 1;
            report.sphere_pair_candidates += 1;
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
