//! Complete affine paths for the original prepared scene's rigid spheres.
//!
//! Swept closed boxes prune distant pairs before the original directed sphere
//! query. Pair laws remain bound by the scene constructor; this path owner
//! supplies no contact response or guessed interface coefficients.

use super::{
    FiniteSolidSphereScene, FiniteSolidSphereSceneError, SweepConfig, SweepItem,
    sphere_pair_witness, sweep_candidates, swept_sphere_aabb,
};
use fs_exec::Cx;
use fs_mbd::{RigidBodyState, Vec3};
use fs_query::triangle::{
    TriangleQueryError,
    sweep::{LinearSweepControls, LinearSweepVerdict, linear_spheres},
};

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) enum SpherePairPathVerdict {
    Clear,
    PossibleContact {
        first: u64,
        second: u64,
        fraction_start: f64,
        fraction_end: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SpherePairPathReport {
    pub verdict: SpherePairPathVerdict,
    pub compared_pairs: usize,
    pub candidate_pairs: usize,
    pub tested_pairs: usize,
    pub examined_windows: usize,
    /// A lower bound for all pairs, including broad-phase exclusions. Absent
    /// when fewer than two spheres exist or a path remains unresolved.
    pub minimum_clearance_m: Option<f64>,
}

#[derive(Debug)]
pub(crate) enum SpherePairPathError {
    Scene(FiniteSolidSphereSceneError),
    Query(TriangleQueryError),
    PathBudget,
    Cancelled,
}

impl FiniteSolidSphereScene {
    /// Remaining work limits belong to the enclosing complete scene path.
    /// Zero remaining capacity is permitted only if no corresponding work is
    /// needed. Setup remains linear storage in the admitted scene size.
    #[allow(clippy::too_many_arguments)]
    pub(crate) fn affine_pair_paths(
        &self,
        before: &[RigidBodyState],
        after: &[RigidBodyState],
        padding_m: f64,
        maximum_visits: usize,
        maximum_candidates: usize,
        controls: LinearSweepControls,
        cx: &Cx<'_>,
    ) -> Result<SpherePairPathReport, SpherePairPathError> {
        use SpherePairPathError as E;
        poll(cx)?;
        if before.len() != self.spheres.len()
            || after.len() != self.spheres.len()
            || !padding_m.is_finite()
            || padding_m < 0.0
            || controls.maximum_windows > 1 << 20
            || !controls.fraction_tolerance.is_finite()
            || controls.fraction_tolerance <= 0.0
            || controls.fraction_tolerance >= 1.0
        {
            return Err(E::Scene(FiniteSolidSphereSceneError::Invalid(
                "matching sphere endpoints and bounded affine path controls required",
            )));
        }
        let mut items = Vec::with_capacity(self.spheres.len());
        for sphere in &self.spheres {
            poll(cx)?;
            let bounds = swept_sphere_aabb(
                array(before[sphere.rigid_body].pose().position_world()),
                array(after[sphere.rigid_body].pose().position_world()),
                sphere.radius_m,
                padding_m,
            )
            .map_err(|e| E::Scene(e.into()))?;
            items.push(SweepItem {
                id: sphere.id,
                bounds,
            });
        }
        let broad = sweep_candidates(
            &items,
            SweepConfig {
                max_candidates: self.maximum_pairs.min(maximum_candidates),
                max_comparisons: self.maximum_pairs.min(maximum_visits),
            },
            || cx.checkpoint().is_err(),
        )
        .map_err(|e| E::Scene(e.into()))?;
        let mut report = SpherePairPathReport {
            verdict: SpherePairPathVerdict::Clear,
            compared_pairs: broad.checked_pairs,
            candidate_pairs: broad.pairs.len(),
            tested_pairs: 0,
            examined_windows: 0,
            // Each excluded pair has a separating coordinate beyond both
            // outward radius+padding boxes. Capping at one padding retains a
            // finite conservative lower bound even if every pair was culled.
            minimum_clearance_m: (self.spheres.len() >= 2).then_some(padding_m),
        };
        for pair in &broad.pairs {
            poll(cx)?;
            if report.compared_pairs.saturating_add(report.tested_pairs) >= maximum_visits
                || report.examined_windows >= controls.maximum_windows
            {
                return Err(E::PathBudget);
            }
            let a = self.spheres[pair.index_a];
            let b = self.spheres[pair.index_b];
            // The original paired material witness admits the initial state.
            // A path proof is never permission to begin inside another sphere.
            let initial = sphere_pair_witness(
                a.rigid_body,
                before[a.rigid_body],
                a.radius_m,
                b.rigid_body,
                before[b.rigid_body],
                b.radius_m,
            )
            .map_err(|e| E::Scene(e.into()))?;
            let gap_m = initial.centre_gap_m.min(initial.material_gap_m);
            if gap_m < 0.0 {
                return Err(E::Scene(FiniteSolidSphereSceneError::Penetration {
                    first: a.id,
                    second: b.id,
                    gap_m,
                }));
            }
            let mut remaining = controls;
            remaining.maximum_windows -= report.examined_windows;
            let path = linear_spheres(
                array(before[a.rigid_body].pose().position_world()),
                array(after[a.rigid_body].pose().position_world()),
                a.radius_m,
                array(before[b.rigid_body].pose().position_world()),
                array(after[b.rigid_body].pose().position_world()),
                b.radius_m,
                remaining,
                cx,
            )
            .map_err(E::Query)?;
            report.tested_pairs += 1;
            report.examined_windows += path.examined_windows;
            match path.verdict {
                LinearSweepVerdict::Clear { minimum_clearance } => {
                    report.minimum_clearance_m = report
                        .minimum_clearance_m
                        .map(|previous| previous.min(minimum_clearance));
                }
                LinearSweepVerdict::PossibleContact {
                    fraction_start,
                    fraction_end,
                } => {
                    let earlier = match report.verdict {
                        SpherePairPathVerdict::Clear => true,
                        SpherePairPathVerdict::PossibleContact {
                            first,
                            second,
                            fraction_start: previous,
                            ..
                        } => {
                            fraction_start < previous
                                || (fraction_start == previous && (a.id, b.id) < (first, second))
                        }
                    };
                    if earlier {
                        report.verdict = SpherePairPathVerdict::PossibleContact {
                            first: a.id,
                            second: b.id,
                            fraction_start,
                            fraction_end,
                        };
                    }
                }
            }
        }
        if !matches!(report.verdict, SpherePairPathVerdict::Clear) {
            report.minimum_clearance_m = None;
        }
        poll(cx)?;
        Ok(report)
    }
}

fn array(value: Vec3) -> [f64; 3] {
    [value.x, value.y, value.z]
}
fn poll(cx: &Cx<'_>) -> Result<(), SpherePairPathError> {
    cx.checkpoint().map_err(|_| SpherePairPathError::Cancelled)
}
