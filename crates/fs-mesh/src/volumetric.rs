//! Constrained multi-region PLC volumetricization (bead s93ej.1).
//!
//! Type-state: `UnverifiedPlc` → `AdmittedPlc` → `ConstraintRecoveredPlc`
//! → `LabeledTetComplex` → `AuditedLabeledTetComplex`. Only the audited
//! type is a downstream geometry authority. The producer classifies by
//! seed-flood across recovered constraint faces; the auditor reclassifies
//! each retained tet by generalized winding of each region's own closed
//! surface and independently recomputes volumes from a different formula
//! plus the closed-surface triple-product identity.

use crate::delaunay::{GHOST, MeshError, Tetrahedralization, delaunay};
use crate::recovery::{FacetCorrespondence, RecoveryOptions, recover_facets, recover_segments};
use fs_exec::Cx;
use fs_geom::Point3;
use fs_ivl::{Sign, orient3d};
use fs_rep_mesh::{Soup, winding_exact};
use std::collections::{BTreeMap, BTreeSet, VecDeque};
use std::fmt;

/// Facet opposite slot `i`, matching the Delaunay kernel so the remaining
/// vertex is Positive (interior below the facet).
const FACET: [[usize; 3]; 4] = [[1, 3, 2], [0, 2, 3], [0, 3, 1], [0, 1, 2]];

/// Stable region identity. Caller-chosen, unique within one PLC.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct RegionId(pub u32);

/// Declared use of a closed surface component.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum RegionKind {
    /// Retained solid.
    Solid,
    /// Enclosed void: tets are classified then discarded.
    Cavity,
}

/// How a volume witness was accumulated.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum VolumeMethod {
    /// All retained coordinates were dyadic enough that the two
    /// independent f64 formulas agreed exactly.
    ExactDyadic,
    /// Formulas disagreed at the last bits; the witness records both.
    MeasuredF64,
}

/// Admission, recovery, classification, and audit refusals.
#[derive(Debug, Clone, PartialEq, Eq)]
pub enum VolumetricError {
    /// A named scalar or coordinate was non-finite.
    NonFinite {
        /// Field or locus.
        what: &'static str,
    },
    /// Length unit missing or empty.
    MissingUnit,
    /// Fewer than four vertices.
    TooFewVertices {
        /// How many arrived.
        got: usize,
    },
    /// A triangle referenced an out-of-range vertex.
    VertexIndex {
        /// Offending index.
        index: u32,
    },
    /// A triangle repeated a vertex.
    DegenerateTriangle {
        /// Triangle index within its region.
        triangle: usize,
    },
    /// Zero-area triangle (exact `orient3d` collinear in its plane via
    /// a constructed off-plane probe is not needed: the three points
    /// are treated as degenerate when every coordinate axis projection
    /// has `orient2d`-style zero area, implemented as a zero cross).
    ZeroAreaTriangle {
        /// Triangle index within its region.
        triangle: usize,
    },
    /// A region's surface is not a closed 2-manifold.
    NotClosedManifold {
        /// Region that failed.
        region: RegionId,
        /// Stable reason.
        reason: &'static str,
    },
    /// Two regions share an id.
    DuplicateRegion {
        /// Repeated id.
        region: RegionId,
    },
    /// No solid region was declared.
    NoSolid,
    /// A seed sits on a recovered constraint or tet face.
    SeedOnBoundary {
        /// Region whose seed failed.
        region: RegionId,
    },
    /// A seed is not inside any tet of the recovered complex.
    SeedNotLocated {
        /// Region whose seed failed.
        region: RegionId,
    },
    /// Two seeds landed in the same constraint chamber.
    AmbiguousChamber {
        /// First region.
        first: RegionId,
        /// Second region.
        second: RegionId,
    },
    /// An enclosed leftover chamber has no seed.
    UnlabeledChamber,
    /// Constraint recovery left an unrecovered segment or facet.
    UnrecoveredConstraint {
        /// Unrecovered PLC segments.
        segments: u64,
        /// Unrecovered PLC facets.
        facets: u64,
    },
    /// Vertex or tet budget exhausted.
    Budget {
        /// Which cap.
        what: &'static str,
        /// Requested count.
        requested: usize,
        /// Cap.
        maximum: usize,
    },
    /// Independent winding/volume/orientation audit failed.
    Audit {
        /// Stable diagnostic.
        reason: &'static str,
    },
    /// Kernel Delaunay/recovery error.
    Mesh(MeshError),
}

impl fmt::Display for VolumetricError {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        match self {
            Self::NonFinite { what } => write!(f, "{what} must be finite"),
            Self::MissingUnit => write!(f, "length unit is required"),
            Self::TooFewVertices { got } => {
                write!(f, "volumetricization needs at least 4 vertices, got {got}")
            }
            Self::VertexIndex { index } => write!(f, "triangle vertex {index} is out of range"),
            Self::DegenerateTriangle { triangle } => {
                write!(f, "triangle {triangle} repeats a vertex")
            }
            Self::ZeroAreaTriangle { triangle } => {
                write!(f, "triangle {triangle} has zero area")
            }
            Self::NotClosedManifold { region, reason } => {
                write!(f, "region {} is not a closed manifold: {reason}", region.0)
            }
            Self::DuplicateRegion { region } => write!(f, "duplicate region {}", region.0),
            Self::NoSolid => write!(f, "at least one solid region is required"),
            Self::SeedOnBoundary { region } => {
                write!(f, "seed for region {} lies on a mesh face", region.0)
            }
            Self::SeedNotLocated { region } => {
                write!(
                    f,
                    "seed for region {} is not inside the recovered mesh",
                    region.0
                )
            }
            Self::AmbiguousChamber { first, second } => write!(
                f,
                "regions {} and {} share one constraint chamber",
                first.0, second.0
            ),
            Self::UnlabeledChamber => {
                write!(f, "an enclosed chamber has no region seed")
            }
            Self::UnrecoveredConstraint { segments, facets } => write!(
                f,
                "constraint recovery left {segments} segments and {facets} facets unrecovered"
            ),
            Self::Budget {
                what,
                requested,
                maximum,
            } => write!(f, "{what} requested {requested} exceeds cap {maximum}"),
            Self::Audit { reason } => write!(f, "independent audit failed: {reason}"),
            Self::Mesh(err) => write!(f, "{err}"),
        }
    }
}

impl std::error::Error for VolumetricError {}

impl From<MeshError> for VolumetricError {
    fn from(value: MeshError) -> Self {
        Self::Mesh(value)
    }
}

/// Caps and recovery policy for one volumetricization.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct VolumetricPolicy {
    /// SI length unit of the vertex coordinates, e.g. `"m"`.
    pub length_unit: String,
    /// Segment/facet recovery caps.
    pub recovery: RecoveryOptions,
    /// Maximum admitted input vertices (Steiner points may grow past
    /// this only up to the tet cap's implied need; input itself refuses).
    pub max_vertices: usize,
    /// Maximum retained solid tetrahedra.
    pub max_tets: usize,
}

impl VolumetricPolicy {
    /// Conservative fixture-scale defaults.
    #[must_use]
    pub fn fixture_default(length_unit: impl Into<String>) -> Self {
        Self {
            length_unit: length_unit.into(),
            recovery: RecoveryOptions::default(),
            max_vertices: 4_096,
            max_tets: 65_536,
        }
    }
}

/// One closed surface component with a seed and a kind.
#[derive(Debug, Clone, PartialEq)]
pub struct RegionSpec {
    /// Stable id.
    pub id: RegionId,
    /// Solid or cavity.
    pub kind: RegionKind,
    /// Strictly interior seed, in the same units as the vertices.
    pub seed: [f64; 3],
    /// Outward-oriented triangles into the shared vertex table.
    pub triangles: Vec<[u32; 3]>,
}

/// Unchecked caller PLC. The only public constructor of later states
/// is [`UnverifiedPlc::admit`].
#[derive(Debug, Clone, PartialEq)]
pub struct UnverifiedPlc {
    vertices: Vec<[f64; 3]>,
    regions: Vec<RegionSpec>,
}

/// Closed-manifold, finite, uniquely labeled PLC.
#[derive(Debug, Clone, PartialEq)]
pub struct AdmittedPlc {
    vertices: Vec<[f64; 3]>,
    regions: Vec<RegionSpec>,
    unique_facets: Vec<[u32; 3]>,
    unique_segments: Vec<[u32; 2]>,
    policy: VolumetricPolicy,
}

/// Delaunay complex with every PLC constraint present as mesh faces.
#[derive(Debug)]
pub struct ConstraintRecoveredPlc {
    tetra: Tetrahedralization,
    regions: Vec<RegionSpec>,
    correspondence: FacetCorrespondence,
    policy: VolumetricPolicy,
}

/// Carved, region-labeled solid tets. Not yet independently audited.
#[derive(Debug, Clone, PartialEq)]
pub struct LabeledTetComplex {
    positions: Vec<[f64; 3]>,
    tets: Vec<[u32; 4]>,
    region_of_tet: Vec<RegionId>,
    source_faces: Vec<([u32; 3], u32)>,
    length_unit: String,
}

/// Volume witness computed by the auditor, not the producer flood.
#[derive(Debug, Clone, PartialEq)]
pub struct VolumeConservationWitness {
    /// Accumulation method.
    pub method: VolumeMethod,
    /// Length unit copied from policy.
    pub length_unit: String,
    /// Cubic unit (`"{length_unit}^3"`).
    pub cubic_unit: String,
    /// Per-region retained volume (producer triple-product).
    pub per_region_producer: Vec<(RegionId, f64)>,
    /// Per-region retained volume (auditor homogeneous determinant).
    pub per_region_auditor: Vec<(RegionId, f64)>,
    /// Closed-surface volume of each solid region.
    pub per_region_surface: Vec<(RegionId, f64)>,
    /// Discarded cavity volume by auditor tet formula.
    pub excluded_cavity: f64,
    /// Discarded exterior volume by auditor tet formula.
    pub excluded_exterior: f64,
}

/// Independently audited labeled volume. The only type that satisfies
/// downstream geometry authority for this leaf.
#[derive(Debug, Clone, PartialEq)]
pub struct AuditedLabeledTetComplex {
    labeled: LabeledTetComplex,
    witness: VolumeConservationWitness,
}

impl UnverifiedPlc {
    /// Wrap caller geometry. No checks run here.
    #[must_use]
    pub fn new(vertices: Vec<[f64; 3]>, regions: Vec<RegionSpec>) -> Self {
        Self { vertices, regions }
    }

    /// Finite coords, unique region ids, closed 2-manifold per region.
    ///
    /// # Errors
    /// [`VolumetricError`] for every listed admission refusal.
    pub fn admit(
        self,
        policy: VolumetricPolicy,
        cx: &Cx<'_>,
    ) -> Result<AdmittedPlc, VolumetricError> {
        if policy.length_unit.is_empty() {
            return Err(VolumetricError::MissingUnit);
        }
        if self.vertices.len() < 4 {
            return Err(VolumetricError::TooFewVertices {
                got: self.vertices.len(),
            });
        }
        if self.vertices.len() > policy.max_vertices {
            return Err(VolumetricError::Budget {
                what: "vertices",
                requested: self.vertices.len(),
                maximum: policy.max_vertices,
            });
        }
        for p in &self.vertices {
            cx.checkpoint().map_err(MeshError::from)?;
            if !p[0].is_finite() || !p[1].is_finite() || !p[2].is_finite() {
                return Err(VolumetricError::NonFinite { what: "vertex" });
            }
        }
        let mut seen = BTreeSet::new();
        let mut solids = 0usize;
        let mut unique_faces: BTreeSet<[u32; 3]> = BTreeSet::new();
        for region in &self.regions {
            cx.checkpoint().map_err(MeshError::from)?;
            if !seen.insert(region.id) {
                return Err(VolumetricError::DuplicateRegion { region: region.id });
            }
            if !region.seed[0].is_finite()
                || !region.seed[1].is_finite()
                || !region.seed[2].is_finite()
            {
                return Err(VolumetricError::NonFinite { what: "seed" });
            }
            if region.kind == RegionKind::Solid {
                solids += 1;
            }
            check_closed_manifold(&self.vertices, region)?;
            for &tri in &region.triangles {
                let mut key = tri;
                key.sort_unstable();
                unique_faces.insert(key);
            }
        }
        if solids == 0 {
            return Err(VolumetricError::NoSolid);
        }
        let unique_facets: Vec<[u32; 3]> = unique_faces.into_iter().collect();
        let mut segs = BTreeSet::new();
        for f in &unique_facets {
            segs.insert(sorted2(f[0], f[1]));
            segs.insert(sorted2(f[1], f[2]));
            segs.insert(sorted2(f[2], f[0]));
        }
        Ok(AdmittedPlc {
            vertices: self.vertices,
            regions: self.regions,
            unique_facets,
            unique_segments: segs.into_iter().collect(),
            policy,
        })
    }
}

impl AdmittedPlc {
    /// Incremental Delaunay of the PLC vertices plus segment/facet recovery.
    ///
    /// # Errors
    /// Kernel errors, cancellation, or honest unrecovered constraints.
    pub fn recover(self, cx: &Cx<'_>) -> Result<ConstraintRecoveredPlc, VolumetricError> {
        let points: Vec<Point3> = self
            .vertices
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect();
        let mut tetra = delaunay(&points, cx)?;
        let (seg_stats, _) =
            recover_segments(&mut tetra, &self.unique_segments, self.policy.recovery, cx)?;
        let facet_loops: Vec<Vec<u32>> = self
            .unique_facets
            .iter()
            .map(|f| vec![f[0], f[1], f[2]])
            .collect();
        let (facet_stats, correspondence) =
            recover_facets(&mut tetra, &facet_loops, self.policy.recovery, cx)?;
        if seg_stats.unrecovered > 0 || facet_stats.unrecovered > 0 {
            return Err(VolumetricError::UnrecoveredConstraint {
                segments: seg_stats.unrecovered,
                facets: facet_stats.unrecovered,
            });
        }
        if !tetra.audit(false).clean() {
            return Err(VolumetricError::Audit {
                reason: "recovered complex failed the exact Delaunay audit",
            });
        }
        Ok(ConstraintRecoveredPlc {
            tetra,
            regions: self.regions,
            correspondence,
            policy: self.policy,
        })
    }
}

impl ConstraintRecoveredPlc {
    /// Seed-flood through non-constraint faces, discard exterior and cavities.
    ///
    /// # Errors
    /// Seed location, ambiguous/unlabeled chambers, budget.
    pub fn carve_and_label(self, cx: &Cx<'_>) -> Result<LabeledTetComplex, VolumetricError> {
        let mesh = &self.tetra.mesh;
        let walls: BTreeSet<[u32; 3]> = self
            .correspondence
            .rows
            .iter()
            .map(|(face, _)| *face)
            .collect();
        let live = live_real_tets(mesh);
        let mut assigned: BTreeMap<u32, RegionId> = BTreeMap::new();
        let mut chamber_of: BTreeMap<u32, RegionId> = BTreeMap::new();
        for region in &self.regions {
            cx.checkpoint().map_err(MeshError::from)?;
            let tet = locate_seed_tet(mesh, region.seed, &live, &walls).map_err(|on_wall| {
                if on_wall {
                    VolumetricError::SeedOnBoundary { region: region.id }
                } else {
                    VolumetricError::SeedNotLocated { region: region.id }
                }
            })?;
            if let Some(&other) = chamber_of.get(&tet) {
                return Err(VolumetricError::AmbiguousChamber {
                    first: other,
                    second: region.id,
                });
            }
            let chamber = flood_chamber(mesh, tet, &walls);
            for &member in &chamber {
                if let Some(&other) = chamber_of.get(&member) {
                    return Err(VolumetricError::AmbiguousChamber {
                        first: other,
                        second: region.id,
                    });
                }
                chamber_of.insert(member, region.id);
                if region.kind == RegionKind::Solid {
                    assigned.insert(member, region.id);
                }
            }
        }
        let leftover: Vec<u32> = live
            .iter()
            .copied()
            .filter(|t| !chamber_of.contains_key(t))
            .collect();
        let leftover_chambers = leftover_components(mesh, &leftover, &walls);
        for chamber in leftover_chambers {
            cx.checkpoint().map_err(MeshError::from)?;
            if chamber_is_enclosed(mesh, &chamber, &walls) {
                return Err(VolumetricError::UnlabeledChamber);
            }
        }
        let mut kept: Vec<(RegionId, [u32; 4])> = assigned
            .iter()
            .map(|(&t, &region)| (region, mesh.tets[t as usize]))
            .collect();
        kept.sort_unstable_by_key(|(region, tet)| {
            let mut s = *tet;
            s.sort_unstable();
            (*region, s)
        });
        if kept.len() > self.policy.max_tets {
            return Err(VolumetricError::Budget {
                what: "tets",
                requested: kept.len(),
                maximum: self.policy.max_tets,
            });
        }
        if kept.is_empty() {
            return Err(VolumetricError::Audit {
                reason: "no solid tetrahedra remained after carving",
            });
        }
        let tets: Vec<[u32; 4]> = kept.iter().map(|(_, t)| *t).collect();
        let region_of_tet: Vec<RegionId> = kept.iter().map(|(r, _)| *r).collect();
        let kept_faces: BTreeSet<[u32; 3]> =
            tets.iter().flat_map(|tet| tet_sorted_faces(*tet)).collect();
        let source_faces: Vec<([u32; 3], u32)> = self
            .correspondence
            .rows
            .into_iter()
            .filter(|(face, _)| kept_faces.contains(face))
            .collect();
        Ok(LabeledTetComplex {
            positions: mesh.points.clone(),
            tets,
            region_of_tet,
            source_faces,
            length_unit: self.policy.length_unit,
        })
    }
}

impl LabeledTetComplex {
    /// Vertex coordinates, input order then Steiner.
    #[must_use]
    pub fn positions(&self) -> &[[f64; 3]] {
        &self.positions
    }

    /// Positively oriented retained solid tets, deterministic order.
    #[must_use]
    pub fn tets(&self) -> &[[u32; 4]] {
        &self.tets
    }

    /// Region of `tets()[i]`.
    #[must_use]
    pub fn region_of_tet(&self) -> &[RegionId] {
        &self.region_of_tet
    }

    /// Recovered PLC faces that remain on a retained tet, with parent facet.
    #[must_use]
    pub fn source_faces(&self) -> &[([u32; 3], u32)] {
        &self.source_faces
    }

    /// Length unit of the coordinates.
    #[must_use]
    pub fn length_unit(&self) -> &str {
        &self.length_unit
    }

    /// Independent winding, orientation, partition, and volume audit.
    ///
    /// # Errors
    /// [`VolumetricError::Audit`] when any independent check fails.
    pub fn audit(
        self,
        regions: &[RegionSpec],
        discarded_cavity_volume: f64,
        discarded_exterior_volume: f64,
        cx: &Cx<'_>,
    ) -> Result<AuditedLabeledTetComplex, VolumetricError> {
        if self.tets.len() != self.region_of_tet.len() {
            return Err(VolumetricError::Audit {
                reason: "tet/region length mismatch",
            });
        }
        let soups: BTreeMap<RegionId, Soup> = regions
            .iter()
            .map(|r| (r.id, region_soup(&self.positions, r)))
            .collect();
        let seeds: BTreeMap<RegionId, [f64; 3]> = regions.iter().map(|r| (r.id, r.seed)).collect();
        let kinds: BTreeMap<RegionId, RegionKind> =
            regions.iter().map(|r| (r.id, r.kind)).collect();
        let mut producer_vol: BTreeMap<RegionId, f64> = BTreeMap::new();
        let mut auditor_vol: BTreeMap<RegionId, f64> = BTreeMap::new();
        for (i, tet) in self.tets.iter().enumerate() {
            cx.checkpoint().map_err(MeshError::from)?;
            if tet_orient(&self.positions, *tet) != Sign::Positive {
                return Err(VolumetricError::Audit {
                    reason: "retained tet is not positively oriented",
                });
            }
            let region = self.region_of_tet[i];
            if kinds.get(&region) != Some(&RegionKind::Solid) {
                return Err(VolumetricError::Audit {
                    reason: "retained tet is not a declared solid",
                });
            }
            let centroid = tet_centroid(&self.positions, *tet);
            for (id, soup) in &soups {
                let at_q = winding_exact(soup, Point3::new(centroid[0], centroid[1], centroid[2]));
                let seed = seeds[id];
                let at_s = winding_exact(soup, Point3::new(seed[0], seed[1], seed[2]));
                let match_seed = (at_q - at_s).abs() < 0.25;
                if *id == region && !match_seed {
                    return Err(VolumetricError::Audit {
                        reason: "retained tet winding does not match its seed",
                    });
                }
                if *id != region && kinds[id] == RegionKind::Solid && match_seed && at_s.abs() > 0.5
                {
                    return Err(VolumetricError::Audit {
                        reason: "retained tet also matches another solid",
                    });
                }
            }
            *producer_vol.entry(region).or_insert(0.0) += tet_volume_triple(&self.positions, *tet);
            *auditor_vol.entry(region).or_insert(0.0) += tet_volume_homog(&self.positions, *tet);
        }
        let mut per_region_surface = Vec::new();
        let mut exact = true;
        let mut per_region_producer = Vec::new();
        let mut per_region_auditor = Vec::new();
        for region in regions {
            if region.kind != RegionKind::Solid {
                continue;
            }
            let soup = &soups[&region.id];
            let surface = closed_surface_volume(soup);
            let prod = *producer_vol.get(&region.id).unwrap_or(&0.0);
            let aud = *auditor_vol.get(&region.id).unwrap_or(&0.0);
            if prod.to_bits() != aud.to_bits() || (prod - surface).abs() > 1.0e-9 {
                exact = false;
            }
            if (prod - surface).abs() > f64::max(1.0e-6, 1.0e-9 * surface.abs()) {
                return Err(VolumetricError::Audit {
                    reason: "retained tet volume disagrees with the closed-surface identity",
                });
            }
            per_region_surface.push((region.id, surface));
            per_region_producer.push((region.id, prod));
            per_region_auditor.push((region.id, aud));
        }
        let witness = VolumeConservationWitness {
            method: if exact {
                VolumeMethod::ExactDyadic
            } else {
                VolumeMethod::MeasuredF64
            },
            length_unit: self.length_unit.clone(),
            cubic_unit: format!("{}^3", self.length_unit),
            per_region_producer,
            per_region_auditor,
            per_region_surface,
            excluded_cavity: discarded_cavity_volume,
            excluded_exterior: discarded_exterior_volume,
        };
        Ok(AuditedLabeledTetComplex {
            labeled: self,
            witness,
        })
    }
}

impl AuditedLabeledTetComplex {
    /// The carved labeled complex.
    #[must_use]
    pub const fn labeled(&self) -> &LabeledTetComplex {
        &self.labeled
    }

    /// Independent volume witness.
    #[must_use]
    pub const fn witness(&self) -> &VolumeConservationWitness {
        &self.witness
    }
}

/// One-shot type-state pipeline.
///
/// # Errors
/// Any admission, recovery, carving, or audit refusal.
pub fn volumetricize(
    plc: UnverifiedPlc,
    policy: VolumetricPolicy,
    cx: &Cx<'_>,
) -> Result<AuditedLabeledTetComplex, VolumetricError> {
    let regions = plc.regions.clone();
    let admitted = plc.admit(policy, cx)?;
    let recovered = admitted.recover(cx)?;
    let (cavity_vol, exterior_vol) = discarded_volumes(&recovered);
    let labeled = recovered.carve_and_label(cx)?;
    labeled.audit(&regions, cavity_vol, exterior_vol, cx)
}

fn discarded_volumes(recovered: &ConstraintRecoveredPlc) -> (f64, f64) {
    // Best-effort diagnostic only: the auditor does not trust these
    // numbers for a pass. They record how much volume the producer
    // threw away so a later consumer can see the carve.
    let mesh = &recovered.tetra.mesh;
    let walls: BTreeSet<[u32; 3]> = recovered
        .correspondence
        .rows
        .iter()
        .map(|(face, _)| *face)
        .collect();
    let live = live_real_tets(mesh);
    let mut claimed = BTreeSet::new();
    let mut cavity = 0.0;
    for region in &recovered.regions {
        if let Ok(tet) = locate_seed_tet(mesh, region.seed, &live, &walls) {
            let chamber = flood_chamber(mesh, tet, &walls);
            if region.kind == RegionKind::Cavity {
                for &t in &chamber {
                    cavity += tet_volume_homog(&mesh.points, mesh.tets[t as usize]);
                }
            }
            claimed.extend(chamber);
        }
    }
    let mut exterior = 0.0;
    for t in live {
        if !claimed.contains(&t) {
            exterior += tet_volume_homog(&mesh.points, mesh.tets[t as usize]);
        }
    }
    (cavity, exterior)
}

fn check_closed_manifold(
    vertices: &[[f64; 3]],
    region: &RegionSpec,
) -> Result<(), VolumetricError> {
    if region.triangles.len() < 4 {
        return Err(VolumetricError::NotClosedManifold {
            region: region.id,
            reason: "fewer than 4 triangles",
        });
    }
    let n = u32::try_from(vertices.len()).unwrap_or(u32::MAX);
    let mut directed: BTreeMap<[u32; 2], i32> = BTreeMap::new();
    let mut undirected: BTreeMap<[u32; 2], u32> = BTreeMap::new();
    for (ti, &tri) in region.triangles.iter().enumerate() {
        if tri[0] == tri[1] || tri[1] == tri[2] || tri[2] == tri[0] {
            return Err(VolumetricError::DegenerateTriangle { triangle: ti });
        }
        for &v in &tri {
            if v >= n {
                return Err(VolumetricError::VertexIndex { index: v });
            }
        }
        if triangle_area2(vertices, tri) == 0.0 {
            return Err(VolumetricError::ZeroAreaTriangle { triangle: ti });
        }
        for e in [[tri[0], tri[1]], [tri[1], tri[2]], [tri[2], tri[0]]] {
            *directed.entry(e).or_insert(0) += 1;
            *undirected.entry(sorted2(e[0], e[1])).or_insert(0) += 1;
        }
    }
    for count in undirected.values() {
        if *count != 2 {
            return Err(VolumetricError::NotClosedManifold {
                region: region.id,
                reason: "an edge is not used exactly twice",
            });
        }
    }
    for (edge, count) in &directed {
        let opp = [edge[1], edge[0]];
        if *count != 1 || directed.get(&opp) != Some(&1) {
            return Err(VolumetricError::NotClosedManifold {
                region: region.id,
                reason: "an edge is not matched by its opposite orientation",
            });
        }
    }
    Ok(())
}

fn triangle_area2(vertices: &[[f64; 3]], tri: [u32; 3]) -> f64 {
    let a = vertices[tri[0] as usize];
    let b = vertices[tri[1] as usize];
    let c = vertices[tri[2] as usize];
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let n = [
        u[1] * v[2] - u[2] * v[1],
        u[2] * v[0] - u[0] * v[2],
        u[0] * v[1] - u[1] * v[0],
    ];
    n[0].mul_add(n[0], n[1].mul_add(n[1], n[2] * n[2]))
}

fn sorted2(a: u32, b: u32) -> [u32; 2] {
    if a < b { [a, b] } else { [b, a] }
}

fn sorted3(mut f: [u32; 3]) -> [u32; 3] {
    f.sort_unstable();
    f
}

fn live_real_tets(mesh: &crate::delaunay::Mesh) -> Vec<u32> {
    (0..mesh.tets.len() as u32)
        .filter(|&t| mesh.alive[t as usize] && mesh.tets[t as usize][3] != GHOST)
        .collect()
}

fn locate_seed_tet(
    mesh: &crate::delaunay::Mesh,
    p: [f64; 3],
    live: &[u32],
    walls: &BTreeSet<[u32; 3]>,
) -> Result<u32, bool> {
    let mut interior = Vec::new();
    let mut internal_boundary = Vec::new();
    let mut on_wall = false;
    for &t in live {
        match tet_vs_point(mesh, t, p, walls) {
            PointInTet::Interior => interior.push(t),
            PointInTet::Boundary { on_constraint } => {
                if on_constraint {
                    on_wall = true;
                } else {
                    internal_boundary.push(t);
                }
            }
            PointInTet::Outside => {}
        }
    }
    if !interior.is_empty() {
        interior.sort_unstable();
        return Ok(interior[0]);
    }
    if on_wall {
        return Err(true);
    }
    if !internal_boundary.is_empty() {
        internal_boundary.sort_unstable();
        return Ok(internal_boundary[0]);
    }
    Err(false)
}

enum PointInTet {
    Interior,
    Boundary { on_constraint: bool },
    Outside,
}

fn tet_vs_point(
    mesh: &crate::delaunay::Mesh,
    t: u32,
    p: [f64; 3],
    walls: &BTreeSet<[u32; 3]>,
) -> PointInTet {
    let tv = mesh.tets[t as usize];
    let mut boundary = false;
    let mut on_constraint = false;
    for i in 0..4 {
        let f = [tv[FACET[i][0]], tv[FACET[i][1]], tv[FACET[i][2]]];
        match orient3d(
            mesh.points[f[0] as usize],
            mesh.points[f[1] as usize],
            mesh.points[f[2] as usize],
            p,
        ) {
            Sign::Negative => return PointInTet::Outside,
            Sign::Zero => {
                boundary = true;
                if walls.contains(&sorted3(f)) {
                    on_constraint = true;
                }
            }
            Sign::Positive => {}
        }
    }
    if boundary {
        PointInTet::Boundary { on_constraint }
    } else {
        PointInTet::Interior
    }
}

fn flood_chamber(
    mesh: &crate::delaunay::Mesh,
    start: u32,
    walls: &BTreeSet<[u32; 3]>,
) -> BTreeSet<u32> {
    let mut seen = BTreeSet::new();
    let mut q = VecDeque::new();
    seen.insert(start);
    q.push_back(start);
    while let Some(t) = q.pop_front() {
        for i in 0..4 {
            let n = mesh.adj[t as usize][i];
            if n == GHOST || !mesh.alive[n as usize] || mesh.tets[n as usize][3] == GHOST {
                continue;
            }
            let face = sorted3(mesh.facet_verts(t, i));
            if walls.contains(&face) {
                continue;
            }
            if seen.insert(n) {
                q.push_back(n);
            }
        }
    }
    seen
}

fn leftover_components(
    mesh: &crate::delaunay::Mesh,
    leftover: &[u32],
    walls: &BTreeSet<[u32; 3]>,
) -> Vec<BTreeSet<u32>> {
    let leftover_set: BTreeSet<u32> = leftover.iter().copied().collect();
    let mut remaining = leftover_set.clone();
    let mut out = Vec::new();
    while let Some(&start) = remaining.iter().next() {
        let chamber = flood_chamber(mesh, start, walls)
            .into_iter()
            .filter(|t| leftover_set.contains(t))
            .collect::<BTreeSet<_>>();
        for t in &chamber {
            remaining.remove(t);
        }
        out.push(chamber);
    }
    out
}

fn chamber_is_enclosed(
    mesh: &crate::delaunay::Mesh,
    chamber: &BTreeSet<u32>,
    walls: &BTreeSet<[u32; 3]>,
) -> bool {
    for &t in chamber {
        for i in 0..4 {
            let n = mesh.adj[t as usize][i];
            let face = sorted3(mesh.facet_verts(t, i));
            if walls.contains(&face) {
                continue;
            }
            if n == GHOST || !mesh.alive[n as usize] || mesh.tets[n as usize][3] == GHOST {
                return false;
            }
        }
    }
    true
}

fn tet_sorted_faces(tet: [u32; 4]) -> [[u32; 3]; 4] {
    [
        sorted3([tet[1], tet[3], tet[2]]),
        sorted3([tet[0], tet[2], tet[3]]),
        sorted3([tet[0], tet[3], tet[1]]),
        sorted3([tet[0], tet[1], tet[2]]),
    ]
}

fn tet_orient(positions: &[[f64; 3]], tet: [u32; 4]) -> Sign {
    orient3d(
        positions[tet[0] as usize],
        positions[tet[1] as usize],
        positions[tet[2] as usize],
        positions[tet[3] as usize],
    )
}

fn tet_centroid(positions: &[[f64; 3]], tet: [u32; 4]) -> [f64; 3] {
    let mut acc = [0.0; 3];
    for v in tet {
        let p = positions[v as usize];
        acc[0] += p[0];
        acc[1] += p[1];
        acc[2] += p[2];
    }
    [acc[0] * 0.25, acc[1] * 0.25, acc[2] * 0.25]
}

/// Producer volume: (1/6) scalar triple product of edge vectors.
fn tet_volume_triple(positions: &[[f64; 3]], tet: [u32; 4]) -> f64 {
    let a = positions[tet[0] as usize];
    let b = positions[tet[1] as usize];
    let c = positions[tet[2] as usize];
    let d = positions[tet[3] as usize];
    let u = [b[0] - a[0], b[1] - a[1], b[2] - a[2]];
    let v = [c[0] - a[0], c[1] - a[1], c[2] - a[2]];
    let w = [d[0] - a[0], d[1] - a[1], d[2] - a[2]];
    let det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
        + u[2] * (v[0] * w[1] - v[1] * w[0]);
    det / 6.0
}

/// Auditor volume: (1/6) determinant of the 4×4 homogeneous matrix.
fn tet_volume_homog(positions: &[[f64; 3]], tet: [u32; 4]) -> f64 {
    let p = [
        positions[tet[0] as usize],
        positions[tet[1] as usize],
        positions[tet[2] as usize],
        positions[tet[3] as usize],
    ];
    // det |x y z 1| expanded on the last column of ones.
    let det = det3(p[1], p[2], p[3]) - det3(p[0], p[2], p[3]) + det3(p[0], p[1], p[3])
        - det3(p[0], p[1], p[2]);
    det / 6.0
}

fn det3(a: [f64; 3], b: [f64; 3], c: [f64; 3]) -> f64 {
    a[0] * (b[1] * c[2] - b[2] * c[1]) - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
}

fn region_soup(positions: &[[f64; 3]], region: &RegionSpec) -> Soup {
    Soup {
        positions: positions
            .iter()
            .map(|p| Point3::new(p[0], p[1], p[2]))
            .collect(),
        triangles: region.triangles.clone(),
    }
}

fn closed_surface_volume(soup: &Soup) -> f64 {
    soup.triangles
        .iter()
        .map(|&[i, j, k]| {
            let a = [
                soup.positions[i as usize].x,
                soup.positions[i as usize].y,
                soup.positions[i as usize].z,
            ];
            let b = [
                soup.positions[j as usize].x,
                soup.positions[j as usize].y,
                soup.positions[j as usize].z,
            ];
            let c = [
                soup.positions[k as usize].x,
                soup.positions[k as usize].y,
                soup.positions[k as usize].z,
            ];
            det3(a, b, c) / 6.0
        })
        .sum()
}

/// Eight corners of an axis-aligned box, in the order used by [`box_triangles`].
#[must_use]
pub fn box_vertices(x0: f64, x1: f64, y0: f64, y1: f64, z0: f64, z1: f64) -> Vec<[f64; 3]> {
    vec![
        [x0, y0, z0],
        [x1, y0, z0],
        [x1, y1, z0],
        [x0, y1, z0],
        [x0, y0, z1],
        [x1, y0, z1],
        [x1, y1, z1],
        [x0, y1, z1],
    ]
}

/// Outward triangles for an 8-vertex axis-aligned box in `box_vertices` order.
#[must_use]
pub fn box_triangles(base: u32) -> Vec<[u32; 3]> {
    let i = |k: u32| base + k;
    vec![
        // z = z0, outward −z (CCW when looking toward +z from below)
        [i(0), i(1), i(2)],
        [i(0), i(2), i(3)],
        // z = z1, outward +z
        [i(4), i(6), i(5)],
        [i(4), i(7), i(6)],
        // y = y0, outward −y
        [i(0), i(5), i(1)],
        [i(0), i(4), i(5)],
        // y = y1, outward +y
        [i(3), i(2), i(6)],
        [i(3), i(6), i(7)],
        // x = x0, outward −x
        [i(0), i(3), i(7)],
        [i(0), i(7), i(4)],
        // x = x1, outward +x
        [i(1), i(5), i(6)],
        [i(1), i(6), i(2)],
    ]
}
