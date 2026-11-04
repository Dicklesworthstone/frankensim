//! Compose complete curved-boundary and original sphere-pair affine paths.

use super::{FiniteSolidSphereScene, QuadraticSphereSceneError, validate_bodies};
use crate::finite_solid::contact::scene::path::{SpherePairPathError, SpherePairPathVerdict};
use crate::inelastic_impact::surface::quadratic::{
    QuadraticSolidContactSurface,
    path::{QuadraticSolidPathControls, QuadraticSolidPathReport, QuadraticSolidPathVerdict},
};
use fs_exec::Cx;
use fs_mbd::RigidBodyState;

/// A swept face does not yet have a resolved point-feature contact ID. These
/// original face/pair keys identify its path without inventing such an ID.
/// Equal-time ties order solid faces first, then the canonical sphere IDs.
#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord)]
pub enum QuadraticSphereScenePathKey {
    Solid { sphere: u64, face: usize },
    Spheres { first: u64, second: u64 },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub enum QuadraticSphereScenePathVerdict {
    /// All declared centre and curved boundary paths are nonpenetrating. The
    /// lower bound includes culled pairs and is capped by the authored padding.
    Clear { minimum_clearance_m: f64 },
    /// Earliest uncertified path window; not an asserted collision time.
    PossibleContact {
        key: QuadraticSphereScenePathKey,
        fraction_start: f64,
        fraction_end: f64,
    },
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub struct QuadraticSphereScenePathReport {
    pub verdict: QuadraticSphereScenePathVerdict,
    pub boundary: QuadraticSolidPathReport,
    pub sphere_pair_comparisons: usize,
    pub sphere_pair_candidates: usize,
    pub sphere_pairs_checked: usize,
    pub sphere_pair_windows: usize,
    /// Covers every pair, including broad-phase exclusions. Absent with fewer
    /// than two spheres or an unresolved pair path.
    pub minimum_pair_clearance_m: Option<f64>,
}

impl FiniteSolidSphereScene {
    /// Admit complete curved-solid and sphere-pair paths in one scene. Original
    /// scene preparation already binds every pair law; no interface is inferred.
    /// Initial pair penetration refuses, including when endpoints later separate.
    ///
    /// Query visits, candidate storage and time windows are shared with the
    /// complete boundary query. Pair visits count broad-phase comparisons plus
    /// narrow-phase paths. Linear refit/setup storage is bounded separately by
    /// original mesh and prepared-scene size. All centres and solid nodes move
    /// affinely between the supplied endpoints; volume regularity and mechanics
    /// remain separate owners. There is no continuum-flow or self-contact claim.
    #[allow(clippy::too_many_arguments)]
    pub fn affine_quadratic_sphere_paths_at_positions(
        &self,
        surface: &QuadraticSolidContactSurface,
        positions_start: &[[f64; 3]],
        positions_end: &[[f64; 3]],
        before: &[RigidBodyState],
        after: &[RigidBodyState],
        controls: QuadraticSolidPathControls,
        cx: &Cx<'_>,
    ) -> Result<QuadraticSphereScenePathReport, QuadraticSphereSceneError> {
        use QuadraticSphereSceneError as E;
        validate_bodies(self, before, cx)?;
        validate_bodies(self, after, cx)?;
        let boundary = surface.affine_sphere_paths_at_positions(
            positions_start,
            positions_end,
            before,
            after,
            self.spheres(),
            controls,
            cx,
        )?;
        let mut time = controls.sweep.time;
        time.maximum_windows = time
            .maximum_windows
            .saturating_sub(boundary.examined_windows);
        let pairs = self
            .affine_pair_paths(
                before,
                after,
                controls.initial.surface.activation_tolerance_m,
                controls
                    .initial
                    .surface
                    .maximum_query_visits
                    .saturating_sub(
                        boundary
                            .visited_nodes
                            .saturating_add(boundary.tested_triangles),
                    ),
                controls
                    .initial
                    .surface
                    .maximum_candidates
                    .saturating_sub(boundary.candidate_triangles),
                time,
                cx,
            )
            .map_err(|error| match error {
                SpherePairPathError::Scene(error) => E::Scene(error),
                SpherePairPathError::Query(error) => E::Query(error),
                SpherePairPathError::PathBudget => E::PathBudget,
                SpherePairPathError::Cancelled => E::Cancelled,
            })?;
        let mut verdict = match boundary.verdict {
            QuadraticSolidPathVerdict::Clear {
                minimum_clearance_m,
            } => QuadraticSphereScenePathVerdict::Clear {
                minimum_clearance_m: pairs
                    .minimum_clearance_m
                    .map_or(minimum_clearance_m, |pair| minimum_clearance_m.min(pair)),
            },
            QuadraticSolidPathVerdict::PossibleContact {
                sphere_id,
                face,
                fraction_start,
                fraction_end,
            } => QuadraticSphereScenePathVerdict::PossibleContact {
                key: QuadraticSphereScenePathKey::Solid {
                    sphere: sphere_id,
                    face,
                },
                fraction_start,
                fraction_end,
            },
        };
        if let SpherePairPathVerdict::PossibleContact {
            first,
            second,
            fraction_start,
            fraction_end,
        } = pairs.verdict
        {
            let key = QuadraticSphereScenePathKey::Spheres { first, second };
            let earlier = match verdict {
                QuadraticSphereScenePathVerdict::Clear { .. } => true,
                QuadraticSphereScenePathVerdict::PossibleContact {
                    key: previous_key,
                    fraction_start: previous,
                    ..
                } => {
                    fraction_start < previous || (fraction_start == previous && key < previous_key)
                }
            };
            if earlier {
                verdict = QuadraticSphereScenePathVerdict::PossibleContact {
                    key,
                    fraction_start,
                    fraction_end,
                };
            }
        }
        cx.checkpoint().map_err(|_| E::Cancelled)?;
        Ok(QuadraticSphereScenePathReport {
            verdict,
            boundary,
            sphere_pair_comparisons: pairs.compared_pairs,
            sphere_pair_candidates: pairs.candidate_pairs,
            sphere_pairs_checked: pairs.tested_pairs,
            sphere_pair_windows: pairs.examined_windows,
            minimum_pair_clearance_m: pairs.minimum_clearance_m,
        })
    }
}
