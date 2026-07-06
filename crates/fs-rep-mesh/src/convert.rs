//! Certified converter mesh → SDF (plan §7.3 edge 1): exact point-triangle
//! distance + winding sign, sampled onto fs-rep-sdf grids — "the strongest
//! edge in the conversion graph".
//!
//! Certificate honesty (the bead's core requirement): each SAMPLE is exact
//! distance up to fp rounding, and the SIGN is exact ONLY when the input
//! is closed and consistently oriented. The converter therefore assesses
//! input quality first: clean meshes get an enclosure-grade `Certified`
//! receipt; soups with boundary edges or non-manifold fins get an HONESTLY
//! DOWNGRADED Estimate receipt whose model evidence names the winding-sign
//! heuristic and the defect counts (Evidence reflects input quality —
//! nothing is silently promoted).
//!
//! The incremental path ([`IncrementalMeshSdf`]) re-samples only tiles
//! touched by an edit region and is BIT-IDENTICAL to full regeneration
//! (samples are recomputed at exactly the original positions — the G5
//! law, proven in rmesh-007).

use crate::chart::MeshChart;
use crate::winding::Soup;
use fs_evidence::{
    Evidence, ModelEvidence, NumericalCertificate, ProvenanceHash, SensitivitySummary,
    StatisticalCertificate, ValidityDomain,
};
use fs_exec::{Cancelled, Cx};
use fs_geom::{Aabb, Chart};
use fs_rep_sdf::{SdfBuildError, TiledSdf};
use std::collections::BTreeMap;

/// Input-quality assessment: what the sign certificate may claim.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MeshQuality {
    /// Edges traversed exactly once (open boundary).
    pub boundary_edges: usize,
    /// Edges traversed more than twice (non-manifold fins).
    pub nonmanifold_edges: usize,
}

impl MeshQuality {
    /// True when the winding sign is exact (closed, 2-manifold edge use).
    #[must_use]
    pub fn sign_certified(&self) -> bool {
        self.boundary_edges == 0 && self.nonmanifold_edges == 0
    }
}

/// Assess a soup's edge usage (deterministic).
#[must_use]
pub fn assess_quality(soup: &Soup) -> MeshQuality {
    let mut edge_use: BTreeMap<[u32; 2], u32> = BTreeMap::new();
    for tri in &soup.triangles {
        for c in 0..3 {
            let (a, b) = (tri[c], tri[(c + 1) % 3]);
            let key = if a < b { [a, b] } else { [b, a] };
            *edge_use.entry(key).or_insert(0) += 1;
        }
    }
    MeshQuality {
        boundary_edges: edge_use.values().filter(|&&c| c == 1).count(),
        nonmanifold_edges: edge_use.values().filter(|&&c| c > 2).count(),
    }
}

/// Conversion failure (the underlying dense-build refusals pass through).
pub type MeshSdfError = SdfBuildError;

/// Convert a mesh chart to a dense tiled SDF with an honesty-graded
/// receipt: `Certified`-grade enclosure for clean input, Estimate + named
/// heuristic for soup (see module docs). The QoI is the field's declared
/// error bound.
///
/// # Errors
/// [`SdfBuildError`] refusals from the dense sampler (teaching text).
pub fn mesh_to_sdf(
    chart: &MeshChart,
    target_h: f64,
    cx: &Cx<'_>,
) -> Result<Evidence<TiledSdf>, MeshSdfError> {
    let quality = assess_quality(chart.soup());
    let sdf = TiledSdf::build(chart, target_h, cx)?;
    let bound = sdf.bound();
    let provenance = ProvenanceHash::chain(
        "convert/mesh-to-sdf",
        &[ProvenanceHash::of_bytes(chart.name().as_bytes())],
    );
    let receipt = if quality.sign_certified() {
        Evidence {
            value: sdf,
            qoi: bound,
            numerical: NumericalCertificate::enclosure(0.0, bound),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence::none(),
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        }
    } else {
        // Honest downgrade: the sign is a winding heuristic on this input.
        Evidence {
            value: sdf,
            qoi: bound,
            numerical: NumericalCertificate::estimate(0.0, bound),
            statistical: StatisticalCertificate::None,
            model: ModelEvidence {
                cards: vec!["winding-sign-heuristic".to_string()],
                assumptions: vec![format!(
                    "input is not a closed 2-manifold ({} boundary edges, {} non-manifold \
                     edges): the sign near defects is a generalized-winding vote, not exact",
                    quality.boundary_edges, quality.nonmanifold_edges
                )],
                validity: ValidityDomain::unconstrained(),
                discrepancy_rel: 0.0,
                in_domain: true,
            },
            sensitivity: SensitivitySummary::default(),
            provenance,
            adjoint_ref: None,
        }
    };
    Ok(receipt)
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
    /// Build the initial field.
    ///
    /// # Errors
    /// [`SdfBuildError`] from the dense sampler.
    pub fn build(chart: MeshChart, target_h: f64, cx: &Cx<'_>) -> Result<Self, MeshSdfError> {
        let sdf = TiledSdf::build(&chart, target_h, cx)?;
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
    /// [`Cancelled`] mid-refresh (samples written so far are complete).
    pub fn update(
        &mut self,
        edited: MeshChart,
        dirty: Aabb,
        cx: &Cx<'_>,
    ) -> Result<(), Cancelled> {
        self.chart = edited;
        self.last_update_samples = self.sdf.resample_box(&self.chart, dirty, cx)?;
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

    fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
        let gate = CancelGate::new();
        let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
        pool.scope(|arena| {
            let cx = Cx::new(
                &gate,
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

    #[test]
    fn quality_assessment_distinguishes_closed_from_soup() {
        let closed = shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0);
        let q = assess_quality(&closed);
        assert!(q.sign_certified(), "{q:?}");
        let open = shapes::corrupt(closed, 0, 0, 0..0, Some(3));
        let q = assess_quality(&open);
        assert_eq!(q.boundary_edges, 3);
        assert!(!q.sign_certified());
    }

    #[test]
    fn clean_input_certifies_and_soup_downgrades() {
        with_cx(|cx| {
            let clean = MeshChart::new(shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0));
            let receipt = mesh_to_sdf(&clean, 0.2, cx).expect("build");
            assert!(receipt.certified().is_ok(), "clean cube certifies");

            let soup = MeshChart::new(shapes::corrupt(
                shapes::cube(Point3::new(0.0, 0.0, 0.0), 1.0),
                0,
                0,
                0..0,
                Some(3),
            ));
            let receipt = mesh_to_sdf(&soup, 0.2, cx).expect("build");
            let err = receipt.certified().expect_err("open soup must not certify");
            assert!(err.to_string().contains("rigorous"), "{err}");
        });
    }
}
