//! Authority-graded converter mesh → SDF (plan §7.3 edge 1): point-triangle
//! distance + winding sign sampled onto fs-rep-sdf grids.
//!
//! Certificate honesty (the bead's core requirement): point-triangle distance
//! is exact up to fp rounding, but the winding SIGN needs component orientation,
//! nesting, vertex-link manifoldness, and self-intersection certificates that
//! this converter does not yet have. Edge-use and aggregate-volume checks are
//! useful diagnostics, but cannot establish that theorem. Generic mesh inputs
//! therefore produce at best an Estimate payload and receipt whose model
//! evidence names the winding-sign heuristic; weak/NoClaim sampled payloads can
//! never be promoted by mesh quality.
//!
//! The incremental path ([`IncrementalMeshSdf`]) re-samples only tiles
//! touched by an edit region and is BIT-IDENTICAL to full regeneration
//! (samples are recomputed at exactly the original positions — the G5
//! law, proven in rmesh-007).

use crate::chart::MeshChart;
use crate::winding::Soup;
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, NumericalKind, ProvenanceHash,
    SensitivitySummary, StatisticalCertificate, ValidityDomain,
};
use fs_exec::Cx;
use fs_geom::{Aabb, Chart, ChartSample, Differentiability, Point3, TraceStepClaim};
use fs_rep_sdf::{SdfBuildError, TiledSdf};
use std::collections::BTreeMap;

/// Deterministic input-quality diagnostics for the winding-sign model.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MeshQuality {
    /// Triangles containing at least one out-of-range vertex reference.
    pub invalid_triangles: usize,
    /// Edges traversed exactly once (open boundary).
    pub boundary_edges: usize,
    /// Edges traversed more than twice (non-manifold fins).
    pub nonmanifold_edges: usize,
    /// Two-use edges traversed in the same direction by both incident faces.
    pub orientation_conflicts: usize,
    /// Whether the soup's aggregate signed volume is finite and positive under
    /// the winding convention used by `MeshChart`.
    pub outward_oriented: bool,
}

impl MeshQuality {
    /// Whether basic edge-use and aggregate-orientation checks pass.
    ///
    /// This is deliberately not a sign certificate. Disconnected components
    /// can have opposing orientations while their aggregate signed volume is
    /// positive, and edge counts do not prove vertex-link manifoldness,
    /// nesting, or freedom from self-intersection.
    #[must_use]
    pub fn passes_basic_orientation_checks(&self) -> bool {
        self.invalid_triangles == 0
            && self.boundary_edges == 0
            && self.nonmanifold_edges == 0
            && self.orientation_conflicts == 0
            && self.outward_oriented
    }

    /// Compatibility alias for [`Self::passes_basic_orientation_checks`].
    ///
    /// Despite the historical name, this result is only a diagnostic screen
    /// and must never authorize an `ExactDistance` theorem or rigorous receipt.
    #[deprecated(
        since = "0.1.0",
        note = "not a sign certificate; use passes_basic_orientation_checks"
    )]
    #[must_use]
    pub fn sign_certified(&self) -> bool {
        self.passes_basic_orientation_checks()
    }
}

/// Assess a soup's edge usage (deterministic).
#[must_use]
pub fn assess_quality(soup: &Soup) -> MeshQuality {
    let mut edge_use: BTreeMap<[u32; 2], (u32, i32)> = BTreeMap::new();
    for tri in &soup.triangles {
        for c in 0..3 {
            let (a, b) = (tri[c], tri[(c + 1) % 3]);
            let key = if a < b { [a, b] } else { [b, a] };
            let (count, orientation) = edge_use.entry(key).or_insert((0, 0));
            *count += 1;
            *orientation += match a.cmp(&b) {
                core::cmp::Ordering::Less => 1,
                core::cmp::Ordering::Greater => -1,
                core::cmp::Ordering::Equal => 0,
            };
        }
    }
    let reference = soup
        .positions
        .first()
        .copied()
        .unwrap_or(Point3::new(0.0, 0.0, 0.0));
    let mut signed_volume_six = 0.0f64;
    let mut invalid_triangles = 0usize;
    for &[ia, ib, ic] in &soup.triangles {
        let vertex = |index: u32| {
            usize::try_from(index)
                .ok()
                .and_then(|index| soup.positions.get(index))
                .copied()
        };
        let (Some(a), Some(b), Some(c)) = (vertex(ia), vertex(ib), vertex(ic)) else {
            invalid_triangles += 1;
            continue;
        };
        let a = a.delta_from(reference);
        let b = b.delta_from(reference);
        let c = c.delta_from(reference);
        signed_volume_six += a.x * (b.y * c.z - b.z * c.y)
            + a.y * (b.z * c.x - b.x * c.z)
            + a.z * (b.x * c.y - b.y * c.x);
    }
    MeshQuality {
        invalid_triangles,
        boundary_edges: edge_use.values().filter(|&&(count, _)| count == 1).count(),
        nonmanifold_edges: edge_use.values().filter(|&&(count, _)| count > 2).count(),
        orientation_conflicts: edge_use
            .values()
            .filter(|&&(count, orientation)| count == 2 && orientation != 0)
            .count(),
        outward_oriented: signed_volume_six.is_finite() && signed_volume_six > 0.0,
    }
}

/// Conversion failure (the underlying dense-build refusals pass through).
pub type MeshSdfError = SdfBuildError;

/// Private adapter for a non-rigorous sampled-field approximation.
///
/// The unsigned distance-to-triangle-set has unit slope, which is useful as a
/// nominal reconstruction scale. The generalized-winding sign is not globally
/// certified, so this adapter explicitly retains `TraceStepClaim::NoClaim` and
/// forwards the raw chart's Estimate/NoClaim numerical authority. It must never
/// be exposed as a general mesh chart or used to authorize sphere tracing.
struct MeshSamplingEstimate<'a> {
    chart: &'a MeshChart,
}

impl Chart for MeshSamplingEstimate<'_> {
    fn eval(&self, x: Point3, cx: &Cx<'_>) -> ChartSample {
        let mut sample = self.chart.eval(x, cx);
        sample.lipschitz = sample.signed_distance.is_finite().then_some(1.0);
        sample
    }

    fn support(&self) -> Aabb {
        self.chart.support()
    }

    fn trace_step_claim(&self) -> TraceStepClaim {
        TraceStepClaim::NoClaim
    }

    fn name(&self) -> &'static str {
        self.chart.name()
    }

    fn differentiability(&self) -> Differentiability {
        self.chart.differentiability()
    }
}

/// Convert a mesh chart to a dense tiled SDF with an honesty-graded
/// receipt: at best Estimate + named winding heuristic for generic mesh input,
/// and NoClaim whenever the sampled payload has no abstract signed-distance
/// authority (see module docs). The QoI is the total abstract-distance estimate
/// bound when available, else the finite nominal-field reconstruction bound.
///
/// # Errors
/// [`SdfBuildError`] refusals from the dense sampler (teaching text).
pub fn mesh_to_sdf(
    chart: &MeshChart,
    target_h: f64,
    cx: &Cx<'_>,
) -> Result<Evidence<TiledSdf>, MeshSdfError> {
    let quality = assess_quality(chart.soup());
    let sdf = TiledSdf::build(&MeshSamplingEstimate { chart }, target_h, cx)?;
    Ok(mesh_sdf_evidence(chart, quality, sdf))
}

fn cap_generic_mesh_authority(sdf: &mut TiledSdf) {
    let _ = sdf.downgrade_abstract_distance_authority(NumericalKind::Estimate);
}

fn mesh_sdf_evidence(
    chart: &MeshChart,
    quality: MeshQuality,
    mut sdf: TiledSdf,
) -> Evidence<TiledSdf> {
    cap_generic_mesh_authority(&mut sdf);
    let payload_kind = sdf.abstract_distance_kind();
    let receipt_kind = payload_kind.max(NumericalKind::Estimate);
    let qoi = sdf
        .abstract_distance_bound()
        .unwrap_or(sdf.nominal_field_bound());
    let provenance = ProvenanceHash::chain(
        "convert/mesh-to-sdf",
        &[ProvenanceHash::of_bytes(chart.name().as_bytes())],
    );
    let numerical = match receipt_kind {
        NumericalKind::Exact | NumericalKind::Enclosure => {
            NumericalCertificate::enclosure(0.0, qoi)
        }
        NumericalKind::Estimate => NumericalCertificate::estimate(0.0, qoi),
        NumericalKind::NoClaim => NumericalCertificate::no_claim(),
    };
    let model = if quality.passes_basic_orientation_checks() {
        ModelEvidence {
            cards: vec!["winding-sign-heuristic".to_string()],
            assumptions: vec![
                "basic edge-use and aggregate-orientation checks pass, but component orientation, \
                 nesting, vertex-link manifoldness, and self-intersection are uncertified; the \
                 winding sign is therefore a model heuristic"
                    .to_string(),
            ],
            validity: ValidityDomain::unconstrained(),
            discrepancy_rel: 0.0,
            in_domain: true,
        }
    } else {
        ModelEvidence {
            cards: vec!["winding-sign-heuristic".to_string()],
            assumptions: vec![format!(
                "input fails the basic edge-use/orientation screen ({} invalid triangles, {} \
                 boundary edges, {} non-manifold edges, {} orientation conflicts, aggregate outward={}): component \
                 orientation, nesting, vertex-link manifoldness, and self-intersection remain \
                 uncertified, so the winding sign is a model heuristic",
                quality.invalid_triangles,
                quality.boundary_edges,
                quality.nonmanifold_edges,
                quality.orientation_conflicts,
                quality.outward_oriented
            )],
            validity: ValidityDomain::unconstrained(),
            discrepancy_rel: 0.0,
            in_domain: true,
        }
    };
    Evidence {
        value: sdf,
        qoi,
        numerical,
        statistical: StatisticalCertificate::None,
        model,
        sensitivity: SensitivitySummary::default(),
        provenance,
        adjoint_ref: None,
    }
}

/// The optimization-loop path: a mesh-backed SDF that re-samples only the
/// tiles an edit touched.
pub struct IncrementalMeshSdf {
    chart: MeshChart,
    sdf: TiledSdf,
    /// Samples refreshed by the last update (dirty-work evidence).
    pub last_update_samples: u64,
}

impl IncrementalMeshSdf {
    /// Build the initial field and cap generic-mesh authority at Estimate
    /// before exposing the raw sampled payload.
    ///
    /// # Errors
    /// [`SdfBuildError`] from the dense sampler.
    pub fn build(chart: MeshChart, target_h: f64, cx: &Cx<'_>) -> Result<Self, MeshSdfError> {
        let mut sdf = TiledSdf::build(&MeshSamplingEstimate { chart: &chart }, target_h, cx)?;
        cap_generic_mesh_authority(&mut sdf);
        Ok(IncrementalMeshSdf {
            chart,
            sdf,
            last_update_samples: 0,
        })
    }

    /// The current field.
    #[must_use]
    pub fn sdf(&self) -> &TiledSdf {
        &self.sdf
    }

    /// The current mesh chart.
    #[must_use]
    pub fn chart(&self) -> &MeshChart {
        &self.chart
    }

    /// Replace the mesh with an edited version and refresh only samples
    /// inside `dirty` (the union box of everything the edit moved,
    /// inflated by the old/new geometry's reach). BIT-IDENTICAL to a full
    /// rebuild when `dirty` covers the true change support (rmesh-007's
    /// G5 law); a too-small `dirty` box is the CALLER's bug — this type
    /// records what it refreshed so audits can catch it.
    ///
    /// # Errors
    /// [`SdfBuildError`] when the dirty box is inadmissible or cancellation is
    /// observed mid-refresh. The update is transactional: on error the prior
    /// chart, field, authority, and refreshed-sample count remain paired and
    /// unchanged. After success, generic-mesh authority is capped at Estimate
    /// before the edited chart is published.
    pub fn update(
        &mut self,
        edited: MeshChart,
        dirty: Aabb,
        cx: &Cx<'_>,
    ) -> Result<(), MeshSdfError> {
        // Resampling stages its values transactionally. Commit the chart and
        // audit count only after the field refresh has succeeded, so callers
        // never observe a new chart paired with an old or partial field.
        let last_update_samples =
            self.sdf
                .resample_box(&MeshSamplingEstimate { chart: &edited }, dirty, cx)?;
        cap_generic_mesh_authority(&mut self.sdf);
        self.chart = edited;
        self.last_update_samples = last_update_samples;
        Ok(())
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::shapes;
    use asupersync::types::Budget;
    use fs_exec::{CancelGate, ExecMode, StreamKey};
    use fs_geom::Point3;

    fn with_gate_cx<R>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                gate,
                arena,
                StreamKey {
                    seed: 0xC0,
                    kernel_id: 1,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            );
            f(&cx)
        })
    }

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        with_gate_cx(&gate, f)
    }

    #[test]
    fn quality_assessment_distinguishes_closed_from_soup() {
        let closed = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        let q = assess_quality(&closed);
        assert!(q.passes_basic_orientation_checks(), "{q:?}");
        let open = shapes::corrupt(closed, 0, 0, 0..0, Some(3));
        let q = assess_quality(&open);
        assert_eq!(q.boundary_edges, 3);
        assert!(!q.passes_basic_orientation_checks());

        let flipped = shapes::corrupt(
            shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
            0,
            0,
            0..1,
            None,
        );
        let q = assess_quality(&flipped);
        assert!(q.orientation_conflicts > 0, "{q:?}");
        assert!(!q.passes_basic_orientation_checks());

        let cube = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        let face_count = cube.triangles.len();
        let inward = shapes::corrupt(cube, 0, 0, 0..face_count, None);
        let q = assess_quality(&inward);
        assert_eq!(
            q.orientation_conflicts, 0,
            "globally reversed is consistent"
        );
        assert!(!q.outward_oriented, "{q:?}");
        assert!(!q.passes_basic_orientation_checks());

        let invalid = Soup {
            positions: vec![Point3::new(0.0, 0.0, 0.0), Point3::new(1.0, 0.0, 0.0)],
            triangles: vec![[0, 1, u32::MAX]],
        };
        let q = assess_quality(&invalid);
        assert_eq!(q.invalid_triangles, 1);
        assert!(!q.passes_basic_orientation_checks());
    }

    #[test]
    fn generic_mesh_inputs_are_estimates_and_name_the_sign_heuristic() {
        with_cx(|cx| {
            let clean = MeshChart::new(shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0));
            let adapter = MeshSamplingEstimate { chart: &clean };
            assert_eq!(adapter.trace_step_claim(), TraceStepClaim::NoClaim);
            let adapter_sample = adapter.eval(Point3::new(2.0, 0.0, 0.0), cx);
            assert_eq!(adapter_sample.lipschitz, Some(1.0));
            assert_eq!(adapter_sample.error.kind, NumericalKind::Estimate);
            let receipt = mesh_to_sdf(&clean, 0.2, cx).expect("build");
            assert_eq!(
                receipt.numerical.kind,
                NumericalKind::Estimate,
                "basic mesh diagnostics are not a global sign certificate"
            );
            assert_eq!(
                receipt.value.abstract_distance_kind(),
                NumericalKind::Estimate
            );
            assert_eq!(
                receipt.qoi.to_bits(),
                receipt
                    .value
                    .abstract_distance_bound()
                    .expect("mesh estimate carries a finite bound")
                    .to_bits()
            );
            assert!(
                receipt
                    .model
                    .cards
                    .contains(&"winding-sign-heuristic".to_string())
            );
            assert!(
                receipt.certified().is_err(),
                "generic mesh must not certify"
            );

            let soup = MeshChart::new(shapes::corrupt(
                shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
                0,
                0,
                0..0,
                Some(3),
            ));
            let receipt = mesh_to_sdf(&soup, 0.2, cx).expect("build");
            assert_eq!(
                receipt.value.abstract_distance_kind(),
                NumericalKind::Estimate,
                "heuristic winding sign weakens the sampled payload itself"
            );
            let err = receipt.certified().expect_err("open soup must not certify");
            assert!(err.to_string().contains("rigorous"), "{err}");
        });
    }

    #[test]
    fn aggregate_orientation_screen_does_not_certify_mixed_components() {
        with_cx(|cx| {
            let mut mixed = shapes::cube(Point3::new(0.0, 0.0, 0.0), 2.0);
            let inner = shapes::cube(Point3::new(0.0, 0.0, 0.0), 0.5);
            let inner_face_count = inner.triangles.len();
            let inward = shapes::corrupt(inner, 0, 0, 0..inner_face_count, None);
            let offset = u32::try_from(mixed.positions.len()).expect("fixture index range");
            mixed.positions.extend(inward.positions);
            mixed
                .triangles
                .extend(inward.triangles.into_iter().map(|[a, b, c]| {
                    [
                        a.checked_add(offset).expect("fixture index range"),
                        b.checked_add(offset).expect("fixture index range"),
                        c.checked_add(offset).expect("fixture index range"),
                    ]
                }));

            let quality = assess_quality(&mixed);
            assert!(
                quality.passes_basic_orientation_checks(),
                "aggregate volume masks the inward component: {quality:?}"
            );
            let receipt = mesh_to_sdf(&MeshChart::new(mixed), 0.5, cx).expect("build");
            assert_eq!(receipt.numerical.kind, NumericalKind::Estimate);
            assert_eq!(
                receipt.value.abstract_distance_kind(),
                NumericalKind::Estimate
            );
            assert!(
                receipt
                    .model
                    .assumptions
                    .iter()
                    .any(|assumption| assumption.contains("component orientation")),
                "model evidence names the missing certificate"
            );
            assert!(receipt.certified().is_err());
        });
    }

    #[test]
    fn incremental_generic_mesh_payload_stays_estimate() {
        with_cx(|cx| {
            let clean = MeshChart::new(shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0));
            let mut incremental = IncrementalMeshSdf::build(clean, 0.5, cx).expect("clean build");
            assert_eq!(
                incremental.sdf().abstract_distance_kind(),
                NumericalKind::Estimate
            );
            let dirty = incremental.sdf().support();
            let open = MeshChart::new(shapes::corrupt(
                shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
                0,
                0,
                0..0,
                Some(3),
            ));
            incremental.update(open, dirty, cx).expect("open update");
            assert_eq!(
                incremental.sdf().abstract_distance_kind(),
                NumericalKind::Estimate,
                "incremental refresh cannot promote generic mesh authority"
            );
        });
    }

    #[test]
    fn clean_mesh_quality_cannot_promote_a_no_claim_sampled_payload() {
        with_cx(|cx| {
            let clean = MeshChart::new(shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0));
            let quality = assess_quality(clean.soup());
            assert!(quality.passes_basic_orientation_checks());
            let clip = Aabb::new(Point3::new(-0.5, -0.5, -0.5), Point3::new(0.5, 0.5, 0.5));
            let sdf =
                TiledSdf::build_clipped(&MeshSamplingEstimate { chart: &clean }, 0.25, clip, cx)
                    .expect("clipped nominal field still samples");
            assert_eq!(sdf.abstract_distance_kind(), NumericalKind::NoClaim);
            let receipt = mesh_sdf_evidence(&clean, quality, sdf);
            assert_eq!(receipt.numerical.kind, NumericalKind::NoClaim);
            assert!(receipt.value.nominal_field_bound().is_finite());
            assert!(receipt.value.abstract_distance_bound().is_none());
            assert!(
                receipt.certified().is_err(),
                "clean input quality cannot launder a NoClaim payload"
            );
        });
    }

    #[test]
    fn cancelled_incremental_update_preserves_paired_state() {
        with_cx(|cx| {
            let mut incremental = IncrementalMeshSdf::build(
                MeshChart::new(shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0)),
                0.5,
                cx,
            )
            .expect("initial build");
            incremental.last_update_samples = 17;
            let before_positions = incremental.chart().soup().positions.clone();
            let probe = Point3::new(0.2, 0.1, -0.1);
            let before_sample = incremental.sdf().eval(probe, cx).signed_distance.to_bits();
            let before_kind = incremental.sdf().abstract_distance_kind();
            let before_bound = incremental.sdf().abstract_distance_bound();
            let dirty = incremental.sdf().support();

            let cancel_gate = CancelGate::new();
            cancel_gate.request();
            with_gate_cx(&cancel_gate, |cancel_cx| {
                let error = incremental
                    .update(
                        MeshChart::new(shapes::cube(Point3::new(0.25, 0.0, 0.0), 1.0)),
                        dirty,
                        cancel_cx,
                    )
                    .expect_err("cancelled refresh");
                assert_eq!(error, SdfBuildError::Cancelled);
            });

            assert_eq!(incremental.chart().soup().positions, before_positions);
            assert_eq!(incremental.last_update_samples, 17);
            assert_eq!(incremental.sdf().abstract_distance_kind(), before_kind);
            assert_eq!(incremental.sdf().abstract_distance_bound(), before_bound);
            assert_eq!(
                incremental.sdf().eval(probe, cx).signed_distance.to_bits(),
                before_sample
            );
        });
    }
}
