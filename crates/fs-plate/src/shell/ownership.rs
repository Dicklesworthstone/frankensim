//! Explicit acoustic-fluid boundary ownership audit and cavity skin closure (`M1.2-W2b`).
//!
//! Every two-sided shell component (`ShellRadiationSurface`) exposes two skin faces per
//! midsurface triangle (`ShellFace::Positive` and `ShellFace::Negative`). In a multi-component
//! tied shell assembly (such as a 6-panel guitar/violin corpus or closed resonator box):
//! - Every exposed skin face must be owned by **exactly one** acoustic domain
//!   ([`AcousticDomainId::InteriorCavity`] or [`AcousticDomainId::Exterior`]).
//! - Faces covered by 2-D glue overlaps ([`FaceTie`]) and the finite-thickness glue-overlap
//!   strips at 1-D tied edges/corners ([`EdgeTie`]) are classified as solid-interior and
//!   owned by **no** fluid.
//! - Open apertures on an interior cavity must be explicitly declared as closed boundary loops
//!   ([`ApertureLoop`]); any undeclared open boundary edge on an interior cavity refuses.
//! - A flipped panel or inconsistent skin-side assignment across a tied edge refuses without
//!   silently reorienting.

use super::junction::{EdgeTie, FaceTie, ShellComponent};
use super::reduction::radiation::{ShellFace, ShellRadiationSurface};
use super::{cross, dot, sub};
use crate::PlateError;
use fs_math::det;
use std::collections::{BTreeMap, BTreeSet};

/// Acoustic fluid domain that owns an exposed shell skin face.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub enum AcousticDomainId {
    /// Enclosed (or aperture-vented) interior acoustic cavity with canonical identifier `0, 1, ...`.
    InteriorCavity(usize),
    /// Unbounded exterior acoustic fluid domain.
    Exterior,
}

/// Canonical `(component, facet, side)` address of one extruded shell skin triangle.
///
/// Ordered deterministically by `(component, facet, side)`.
#[derive(Debug, Clone, Copy, PartialEq, Eq, PartialOrd, Ord, Hash)]
pub struct SkinFaceId {
    /// Component index in `0..components.len()`.
    pub component: usize,
    /// Midsurface triangle index in `0..components[component].mesh.tris.len()`.
    pub facet: usize,
    /// Extruded skin side (`Positive` or `Negative`).
    pub side: ShellFace,
}

/// Explicitly declared open aperture boundary loop (e.g., soundhole or f-hole) on an interior cavity.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct ApertureLoop {
    /// Target interior cavity identifier matching [`AcousticDomainId::InteriorCavity`].
    pub cavity_id: usize,
    /// Component index carrying this aperture boundary loop.
    pub component: usize,
    /// Ordered sequence of boundary node indices forming a closed cycle (`len >= 3`).
    pub nodes: Vec<usize>,
}

/// Complete fluid-boundary ownership assignment and tie topology input to [`FluidBoundaryOwnership::audit`].
#[derive(Debug, Clone)]
pub struct FluidBoundaryAssignments {
    /// Explicit `(SkinFaceId, AcousticDomainId)` assignments for fluid-wetted skin faces.
    pub face_owners: Vec<(SkinFaceId, AcousticDomainId)>,
    /// 1-D edge ties joining component boundary polylines.
    pub edge_ties: Vec<EdgeTie>,
    /// 2-D overlapping face ties whose mating skin facets are solid-interior glue overlaps.
    pub face_ties: Vec<FaceTie>,
    /// Explicitly declared open aperture loops on interior cavities.
    pub open_apertures: Vec<ApertureLoop>,
}

impl FluidBoundaryAssignments {
    /// Build deterministic sorted `(component, facet, side)` ownership assignments from per-component
    /// `(positive_domain, negative_domain)` rules, skipping any facets covered by `face_ties` glue overlaps.
    #[must_use]
    pub fn from_component_sides(
        components: &[ShellComponent<'_>],
        component_domains: &[(AcousticDomainId, AcousticDomainId)],
        edge_ties: Vec<EdgeTie>,
        face_ties: Vec<FaceTie>,
        open_apertures: Vec<ApertureLoop>,
    ) -> Self {
        let solid_faces = collect_face_tie_solid_faces(&face_ties);
        let mut face_owners = Vec::new();
        for (c, comp) in components.iter().enumerate() {
            let Some(&(pos_dom, neg_dom)) = component_domains.get(c) else {
                continue;
            };
            for f in 0..comp.mesh.tris.len() {
                let pos_id = SkinFaceId {
                    component: c,
                    facet: f,
                    side: ShellFace::Positive,
                };
                if !solid_faces.contains(&pos_id) {
                    face_owners.push((pos_id, pos_dom));
                }
                let neg_id = SkinFaceId {
                    component: c,
                    facet: f,
                    side: ShellFace::Negative,
                };
                if !solid_faces.contains(&neg_id) {
                    face_owners.push((neg_id, neg_dom));
                }
            }
        }
        Self {
            face_owners,
            edge_ties,
            face_ties,
            open_apertures,
        }
    }
}

/// Per-domain geometry, area, and enclosed volume audit summary in [`OwnershipReceipt`].
#[derive(Debug, Clone)]
pub struct DomainOwnershipSummary {
    /// Acoustic domain identifier.
    pub domain: AcousticDomainId,
    /// Owned skin faces in deterministic sorted `(component, facet, side)` order.
    pub faces: Vec<SkinFaceId>,
    /// Sum of raw extruded skin triangle areas [$\text{m}^2$] before glue-overlap trimming.
    pub raw_skin_area_m2: f64,
    /// Area [$\text{m}^2$] of glue-overlap strips at tied edges and corners subtracted from this domain.
    pub glue_overlap_area_m2: f64,
    /// Net fluid-wetted surface area [$\text{m}^2$] (`raw_skin_area_m2 - glue_overlap_area_m2`).
    pub wetted_area_m2: f64,
    /// Oriented closure triangles built at tied edges, corners, and declared apertures.
    pub closure_triangles: Vec<[[f64; 3]; 3]>,
    /// Enclosed fluid volume [$\text{m}^3$] for [`AcousticDomainId::InteriorCavity`] domains
    /// (`None` for the unbounded [`AcousticDomainId::Exterior`] domain).
    pub enclosed_volume_m3: Option<f64>,
}

/// Auditable receipt returned by [`FluidBoundaryOwnership::audit`].
#[derive(Debug, Clone)]
pub struct OwnershipReceipt {
    /// Per-domain ownership summaries keyed by [`AcousticDomainId`].
    pub domains: BTreeMap<AcousticDomainId, DomainOwnershipSummary>,
    /// Skin faces classified as solid-interior due to 2-D [`FaceTie`] glue overlaps.
    pub solid_interior_faces: Vec<SkinFaceId>,
    /// Number of tied-edge glue-overlap closure triangles classified as solid-interior.
    pub solid_interior_closure_triangles: usize,
    /// Total solid-interior glue-overlap area [$\text{m}^2$] owned by no fluid.
    pub solid_interior_area_m2: f64,
    /// Number of validated open aperture loops.
    pub open_aperture_loops: usize,
}

/// Stateless fluid-boundary ownership auditor.
pub struct FluidBoundaryOwnership;

fn norm3(v: [f64; 3]) -> f64 {
    det::sqrt(dot(v, v))
}

fn add3(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [a[0] + b[0], a[1] + b[1], a[2] + b[2]]
}

fn scale3(v: [f64; 3], s: f64) -> [f64; 3] {
    [v[0] * s, v[1] * s, v[2] * s]
}

fn collect_face_tie_solid_faces(face_ties: &[FaceTie]) -> BTreeSet<SkinFaceId> {
    let mut set = BTreeSet::new();
    for ftie in face_ties {
        let m_side = match ftie.offset {
            super::junction::TieOffset::Skin { side: -1 } => ShellFace::Negative,
            _ => ShellFace::Positive,
        };
        let s_side = match m_side {
            ShellFace::Positive => ShellFace::Negative,
            ShellFace::Negative => ShellFace::Positive,
        };
        for &mt in &ftie.master_tris {
            set.insert(SkinFaceId {
                component: ftie.master_component,
                facet: mt,
                side: m_side,
            });
        }
        for &st in &ftie.slave_tris {
            set.insert(SkinFaceId {
                component: ftie.slave_component,
                facet: st,
                side: s_side,
            });
        }
    }
    set
}

fn compute_component_directors(comp: &ShellComponent<'_>) -> Result<Vec<[f64; 3]>, PlateError> {
    let n_nodes = comp.mesh.nodes.len();
    if comp.nodal_thickness_m.len() != n_nodes
        || comp
            .nodal_thickness_m
            .iter()
            .any(|&h| !h.is_finite() || h <= 0.0)
    {
        return Err(PlateError::BadSection {
            what: "ownership audit requires positive finite nodal_thickness_m at every node",
        });
    }
    let mut out = vec![[0.0; 3]; n_nodes];
    for (e, tri) in comp.mesh.tris.iter().enumerate() {
        let g = comp.mesh.facet(e)?;
        let area_normal = scale3(g.frame[2], 2.0 * g.area_m2);
        for &i in tri {
            out[i] = add3(out[i], area_normal);
        }
    }
    for (i, d) in out.iter_mut().enumerate() {
        let len = norm3(*d);
        if !len.is_finite() || len <= 0.0 {
            return Err(PlateError::BadSection {
                what: "ownership audit encountered an undefined vertex director",
            });
        }
        *d = scale3(*d, 0.5 * comp.nodal_thickness_m[i] / len);
    }
    Ok(out)
}

impl FluidBoundaryOwnership {
    /// Audit fluid-boundary ownership across all component skins, ties, and declared apertures.
    ///
    /// # Errors
    /// Returns [`PlateError::BadConstraint`] or [`PlateError::BadSection`] if:
    /// - component and skin slices mismatch;
    /// - any non-glue-overlap skin face is missing an owner;
    /// - any skin face is assigned more than once or a glue-overlap face is assigned to a fluid;
    /// - a panel has a flipped normal or inconsistent interior/exterior orientation across a tied edge;
    /// - an interior cavity has an untied boundary edge not covered by a declared [`ApertureLoop`].
    pub fn audit(
        components: &[ShellComponent<'_>],
        skins: &[ShellRadiationSurface],
        assignments: &FluidBoundaryAssignments,
    ) -> Result<OwnershipReceipt, PlateError> {
        if components.is_empty() || components.len() != skins.len() {
            return Err(PlateError::BadConstraint {
                what: "ownership audit requires non-empty matching components and skins slices",
            });
        }

        let mut directors_per_comp: Vec<Vec<[f64; 3]>> = Vec::with_capacity(components.len());
        for (c, (comp, skin)) in components.iter().zip(skins.iter()).enumerate() {
            if comp.id != c {
                return Err(PlateError::BadConstraint {
                    what: "ShellComponent.id must equal its slice index in ownership audit",
                });
            }
            let n_tris = comp.mesh.tris.len();
            if n_tris == 0 || skin.triangles().len() < 2 * n_tris {
                return Err(PlateError::BadConstraint {
                    what: "ShellRadiationSurface panel count does not match component triangles",
                });
            }
            directors_per_comp.push(compute_component_directors(comp)?);
        }

        // Validate FaceTie indices and collect 2-D glue-overlap faces (owned by no fluid).
        for ftie in &assignments.face_ties {
            if ftie.master_component >= components.len()
                || ftie.slave_component >= components.len()
                || ftie.master_component == ftie.slave_component
            {
                return Err(PlateError::BadConstraint {
                    what: "FaceTie component index out of range in ownership audit",
                });
            }
            for &mt in &ftie.master_tris {
                if mt >= components[ftie.master_component].mesh.tris.len() {
                    return Err(PlateError::BadConstraint {
                        what: "FaceTie master triangle index out of range in ownership audit",
                    });
                }
            }
            for &st in &ftie.slave_tris {
                if st >= components[ftie.slave_component].mesh.tris.len() {
                    return Err(PlateError::BadConstraint {
                        what: "FaceTie slave triangle index out of range in ownership audit",
                    });
                }
            }
        }
        let solid_faces_set = collect_face_tie_solid_faces(&assignments.face_ties);

        // Parse and validate face_owners: reject out-of-range, glue-overlap assignment, or duplicate owners.
        let mut owner_map: BTreeMap<SkinFaceId, AcousticDomainId> = BTreeMap::new();
        for &(face_id, domain) in &assignments.face_owners {
            if face_id.component >= components.len()
                || face_id.facet >= components[face_id.component].mesh.tris.len()
            {
                return Err(PlateError::BadConstraint {
                    what: "SkinFaceId component or facet index out of range",
                });
            }
            if solid_faces_set.contains(&face_id) {
                return Err(PlateError::BadConstraint {
                    what: "double owner: glue-overlap face is solid-interior and cannot be owned by a fluid",
                });
            }
            if owner_map.insert(face_id, domain).is_some() {
                return Err(PlateError::BadConstraint {
                    what: "double owner: skin face is assigned to more than one acoustic domain",
                });
            }
        }

        // Deterministic sorted check over all (component, facet, side): every non-glue-overlap face must have an owner,
        // and the two faces of a facet cannot both belong to the same domain.
        let mut domain_faces: BTreeMap<AcousticDomainId, Vec<SkinFaceId>> = BTreeMap::new();
        let mut domain_raw_area: BTreeMap<AcousticDomainId, f64> = BTreeMap::new();
        let mut solid_interior_faces: Vec<SkinFaceId> = Vec::new();
        let mut solid_interior_area_m2 = 0.0_f64;

        for (c, comp) in components.iter().enumerate() {
            for f in 0..comp.mesh.tris.len() {
                let pos_id = SkinFaceId {
                    component: c,
                    facet: f,
                    side: ShellFace::Positive,
                };
                let neg_id = SkinFaceId {
                    component: c,
                    facet: f,
                    side: ShellFace::Negative,
                };
                let pos_owner = if solid_faces_set.contains(&pos_id) {
                    solid_interior_faces.push(pos_id);
                    solid_interior_area_m2 += skins[c].areas_m2()[2 * f];
                    None
                } else {
                    let dom = *owner_map.get(&pos_id).ok_or(PlateError::BadConstraint {
                        what: "missing owner: exposed Positive skin face has no acoustic domain assignment",
                    })?;
                    domain_faces.entry(dom).or_default().push(pos_id);
                    *domain_raw_area.entry(dom).or_insert(0.0) += skins[c].areas_m2()[2 * f];
                    Some(dom)
                };

                let neg_owner = if solid_faces_set.contains(&neg_id) {
                    solid_interior_faces.push(neg_id);
                    solid_interior_area_m2 += skins[c].areas_m2()[2 * f + 1];
                    None
                } else {
                    let dom = *owner_map.get(&neg_id).ok_or(PlateError::BadConstraint {
                        what: "missing owner: exposed Negative skin face has no acoustic domain assignment",
                    })?;
                    domain_faces.entry(dom).or_default().push(neg_id);
                    *domain_raw_area.entry(dom).or_insert(0.0) += skins[c].areas_m2()[2 * f + 1];
                    Some(dom)
                };

                if let (Some(dp), Some(dn)) = (pos_owner, neg_owner)
                    && dp == dn
                {
                    return Err(PlateError::BadConstraint {
                        what: "orientation conflict: both Positive and Negative faces of a facet are assigned to the same domain",
                    });
                }
            }
        }

        // Collect allinterior cavity IDs present in the assignments.
        let cavity_ids: Vec<usize> = domain_faces
            .keys()
            .filter_map(|d| match *d {
                AcousticDomainId::InteriorCavity(id) => Some(id),
                AcousticDomainId::Exterior => None,
            })
            .collect();

        let mut summaries: BTreeMap<AcousticDomainId, DomainOwnershipSummary> = BTreeMap::new();
        let mut total_closure_tris = 0usize;

        // Audit each InteriorCavity for consistent panel orientation, boundary closure, and enclosed volume.
        for cav_id in cavity_ids {
            let dom = AcousticDomainId::InteriorCavity(cav_id);
            let faces = domain_faces.get(&dom).cloned().unwrap_or_default();
            let raw_area = domain_raw_area.get(&dom).copied().unwrap_or(0.0);

            let (closure_tris, overlap_area, volume_m3) = audit_interior_cavity(
                cav_id,
                components,
                skins,
                &directors_per_comp,
                &faces,
                &assignments.edge_ties,
                &assignments.open_apertures,
            )?;

            total_closure_tris += closure_tris.len();
            solid_interior_area_m2 += overlap_area;
            summaries.insert(
                dom,
                DomainOwnershipSummary {
                    domain: dom,
                    faces,
                    raw_skin_area_m2: raw_area,
                    glue_overlap_area_m2: overlap_area,
                    wetted_area_m2: (raw_area - overlap_area).max(0.0),
                    closure_triangles: closure_tris,
                    enclosed_volume_m3: Some(volume_m3),
                },
            );
        }

        // Populate Exterior domain summary if present.
        if let Some(ext_faces) = domain_faces.get(&AcousticDomainId::Exterior) {
            let raw_area = domain_raw_area
                .get(&AcousticDomainId::Exterior)
                .copied()
                .unwrap_or(0.0);
            summaries.insert(
                AcousticDomainId::Exterior,
                DomainOwnershipSummary {
                    domain: AcousticDomainId::Exterior,
                    faces: ext_faces.clone(),
                    raw_skin_area_m2: raw_area,
                    glue_overlap_area_m2: 0.0,
                    wetted_area_m2: raw_area,
                    closure_triangles: Vec::new(),
                    enclosed_volume_m3: None,
                },
            );
        }

        Ok(OwnershipReceipt {
            domains: summaries,
            solid_interior_faces,
            solid_interior_closure_triangles: total_closure_tris,
            solid_interior_area_m2,
            open_aperture_loops: assignments.open_apertures.len(),
        })
    }
}

#[allow(clippy::too_many_lines)]
fn audit_interior_cavity(
    cav_id: usize,
    components: &[ShellComponent<'_>],
    skins: &[ShellRadiationSurface],
    directors_per_comp: &[Vec<[f64; 3]>],
    cavity_faces: &[SkinFaceId],
    edge_ties: &[EdgeTie],
    open_apertures: &[ApertureLoop],
) -> Result<(Vec<[[f64; 3]; 3]>, f64, f64), PlateError> {
    // Map each component participating in this cavity to its uniform side (`Positive` or `Negative`),
    // refusing if a single connected component flips sides mid-panel.
    let mut comp_cavity_side: BTreeMap<usize, ShellFace> = BTreeMap::new();
    let mut comp_facet_count: BTreeMap<usize, usize> = BTreeMap::new();
    for &face in cavity_faces {
        if let Some(&existing) = comp_cavity_side.get(&face.component) {
            if existing != face.side {
                return Err(PlateError::BadConstraint {
                    what: "orientation conflict: component switches cavity skin side across facets",
                });
            }
        } else {
            comp_cavity_side.insert(face.component, face.side);
        }
        *comp_facet_count.entry(face.component).or_insert(0) += 1;
    }

    // For each participating component, collect its directed boundary edges (u -> v) in the
    // cavity-outward triangle orientation.
    // Note: in ShellRadiationSurface, every skin triangle [v0, v1, v2] has normal pointing from
    // the solid wall into the fluid (cavity-inward). Reversing it to [v0, v2, v1] gives the
    // cavity-outward triangle!
    // In terms of midsurface vertex indices:
    // - If side == ShellFace::Negative, radiation.rs built [tri[2], tri[1], tri[0]], so reversing
    //   gives [tri[2], tri[0], tri[1]], which has the SAME cyclic winding as [tri[0], tri[1], tri[2]]!
    // - If side == ShellFace::Positive, radiation.rs built [tri[0], tri[1], tri[2]], so reversing
    //   gives [tri[0], tri[2], tri[1]], which has the OPPOSITE cyclic winding to [tri[0], tri[1], tri[2]]!
    let mut comp_boundary_edges: BTreeMap<usize, Vec<(usize, usize)>> = BTreeMap::new();
    let mut comp_boundary_lookup: BTreeMap<(usize, usize, usize), (usize, usize)> = BTreeMap::new();

    for (&c, &side) in &comp_cavity_side {
        let mesh = components[c].mesh;
        let mut edge_uses: BTreeMap<(usize, usize), Vec<(usize, usize)>> = BTreeMap::new();
        for tri in &mesh.tris {
            let cyc = match side {
                ShellFace::Negative => [tri[0], tri[1], tri[2]],
                ShellFace::Positive => [tri[0], tri[2], tri[1]],
            };
            for k in 0..3 {
                let u = cyc[k];
                let v = cyc[(k + 1) % 3];
                edge_uses.entry((u.min(v), u.max(v))).or_default().push((u, v));
            }
        }
        let mut bdy = Vec::new();
        for (&(lo, hi), uses) in &edge_uses {
            if uses.len() == 1 {
                let dir_edge = uses[0];
                bdy.push(dir_edge);
                comp_boundary_lookup.insert((c, lo, hi), dir_edge);
            }
        }
        comp_boundary_edges.insert(c, bdy);
    }

    // Track which unordered boundary edges (c, u.min(v), u.max(v)) are covered by an EdgeTie or ApertureLoop.
    let mut covered_boundary_edges: BTreeSet<(usize, usize, usize)> = BTreeSet::new();

    // Also accumulate, for each boundary vertex (c, v), the inward cavity offset vectors d_{c'}^{cav}
    // of the transverse panels tied to (c, v), keyed by quantized unit normal direction so a 3-way
    // corner collects each transverse wall offset once.
    let mut vertex_transverse_offsets: BTreeMap<(usize, usize), Vec<[f64; 3]>> = BTreeMap::new();

    let side_sign_fn = |side: ShellFace| -> f64 {
        match side {
            ShellFace::Positive => 1.0,
            ShellFace::Negative => -1.0,
        }
    };

    // Audit every EdgeTie touching any component of this cavity.
    for tie in edge_ties {
        let c_a = tie.a.component;
        let c_b = tie.b.component;
        let a_in = comp_cavity_side.get(&c_a).copied();
        let b_in = comp_cavity_side.get(&c_b).copied();
        match (a_in, b_in) {
            (None, None) => continue,
            (Some(_), None) | (None, Some(_)) => {
                return Err(PlateError::BadConstraint {
                    what: "unclosed interior: EdgeTie connects a cavity panel to a panel not in the cavity",
                });
            }
            (Some(side_a), Some(side_b)) => {
                let s_a = side_sign_fn(side_a);
                let s_b = side_sign_fn(side_b);

                // Verify that every segment of tie.a and tie.b is a boundary edge of c_a and c_b,
                // and check manifold orientation compatibility across the seam.
                let n_a = tie.a.nodes.len();
                let n_b = tie.b.nodes.len();
                let chord_a = sub(
                    components[c_a].mesh.nodes[tie.a.nodes[n_a - 1]],
                    components[c_a].mesh.nodes[tie.a.nodes[0]],
                );
                let chord_b = sub(
                    components[c_b].mesh.nodes[tie.b.nodes[n_b - 1]],
                    components[c_b].mesh.nodes[tie.b.nodes[0]],
                );
                if dot(chord_a, chord_b) <= 0.0 {
                    return Err(PlateError::BadConstraint {
                        what: "orientation conflict: EdgeTie polylines run in opposite directions",
                    });
                }

                // Check first segment boundary winding on A and B: in a consistently oriented closed
                // surface, the cavity-outward boundary half-edges on A and B must run in OPPOSITE
                // 3-D directions along the shared seam!
                let key_a0 = (
                    c_a,
                    tie.a.nodes[0].min(tie.a.nodes[1]),
                    tie.a.nodes[0].max(tie.a.nodes[1]),
                );
                let key_b0 = (
                    c_b,
                    tie.b.nodes[0].min(tie.b.nodes[1]),
                    tie.b.nodes[0].max(tie.b.nodes[1]),
                );
                let &(ua, va) = comp_boundary_lookup.get(&key_a0).ok_or(
                    PlateError::BadConstraint {
                        what: "EdgeTie master segment is not a boundary edge of the component mesh",
                    },
                )?;
                let &(ub, vb) = comp_boundary_lookup.get(&key_b0).ok_or(
                    PlateError::BadConstraint {
                        what: "EdgeTie slave segment is not a boundary edge of the component mesh",
                    },
                )?;
                let dir_a = sub(components[c_a].mesh.nodes[va], components[c_a].mesh.nodes[ua]);
                let dir_b = sub(components[c_b].mesh.nodes[vb], components[c_b].mesh.nodes[ub]);
                if dot(dir_a, dir_b) >= 0.0 {
                    return Err(PlateError::BadConstraint {
                        what: "orientation conflict: tied panels have conflicting cavity-outward windings across seam",
                    });
                }

                // Mark all segments on both sides as covered by this tie.
                for w in tie.a.nodes.windows(2) {
                    let key = (c_a, w[0].min(w[1]), w[0].max(w[1]));
                    if !comp_boundary_lookup.contains_key(&key) {
                        return Err(PlateError::BadConstraint {
                            what: "EdgeTie master segment is not on the component boundary",
                        });
                    }
                    covered_boundary_edges.insert(key);
                }
                for w in tie.b.nodes.windows(2) {
                    let key = (c_b, w[0].min(w[1]), w[0].max(w[1]));
                    if !comp_boundary_lookup.contains_key(&key) {
                        return Err(PlateError::BadConstraint {
                            what: "EdgeTie slave segment is not on the component boundary",
                        });
                    }
                    covered_boundary_edges.insert(key);
                }

                // Representative cavity-inward offset vectors for c_a and c_b along this seam.
                let d_a_rep = scale3(directors_per_comp[c_a][tie.a.nodes[0]], s_a);
                let d_b_rep = scale3(directors_per_comp[c_b][tie.b.nodes[0]], s_b);
                for &na in &tie.a.nodes {
                    push_unique_transverse_offset(
                        vertex_transverse_offsets.entry((c_a, na)).or_default(),
                        directors_per_comp[c_a][na],
                        d_b_rep,
                    );
                }
                for &nb in &tie.b.nodes {
                    push_unique_transverse_offset(
                        vertex_transverse_offsets.entry((c_b, nb)).or_default(),
                        directors_per_comp[c_b][nb],
                        d_a_rep,
                    );
                }
            }
        }
    }

    // Propagate transverse offsets across corner node clusters so that at a 3-way corner where
    // A-B and A-C are tied, B and C also know about each other's transverse wall offset even if
    // the third tie B-C wasn't the one that first visited them.
    propagate_corner_transverse_offsets(
        components,
        &comp_cavity_side,
        directors_per_comp,
        edge_ties,
        &mut vertex_transverse_offsets,
    );

    // Validate declared open apertures for this cavity.
    let mut aperture_closure_tris: Vec<[[f64; 3]; 3]> = Vec::new();
    for ap in open_apertures.iter().filter(|a| a.cavity_id == cav_id) {
        if !comp_cavity_side.contains_key(&ap.component) || ap.nodes.len() < 3 {
            return Err(PlateError::BadConstraint {
                what: "ApertureLoop references a component not in the cavity or has fewer than 3 nodes",
            });
        }
        let c = ap.component;
        let s_c = side_sign_fn(comp_cavity_side[&c]);
        let n_loop = ap.nodes.len();
        let mut loop_pts = Vec::with_capacity(n_loop);
        let mut center = [0.0; 3];
        for i in 0..n_loop {
            let u = ap.nodes[i];
            let v = ap.nodes[(i + 1) % n_loop];
            let key = (c, u.min(v), u.max(v));
            let &(b_u, b_v) = comp_boundary_lookup.get(&key).ok_or(PlateError::BadConstraint {
                what: "ApertureLoop segment is not an untied boundary edge of the component",
            })?;
            if covered_boundary_edges.contains(&key) {
                return Err(PlateError::BadConstraint {
                    what: "ApertureLoop segment is already tied or declared twice",
                });
            }
            covered_boundary_edges.insert(key);
            let p_u = add3(
                components[c].mesh.nodes[b_u],
                scale3(directors_per_comp[c][b_u], s_c),
            );
            let p_v = add3(
                components[c].mesh.nodes[b_v],
                scale3(directors_per_comp[c][b_v], s_c),
            );
            loop_pts.push((p_u, p_v));
            center = add3(center, scale3(p_u, 1.0 / (n_loop as f64)));
        }
        let n_cav_out = scale3(directors_per_comp[c][ap.nodes[0]], -s_c);
        for (p_u, p_v) in loop_pts {
            let cand = [p_u, p_v, center];
            let area_n = cross(sub(cand[1], cand[0]), sub(cand[2], cand[0]));
            if dot(area_n, n_cav_out) >= 0.0 {
                aperture_closure_tris.push(cand);
            } else {
                aperture_closure_tris.push([p_v, p_u, center]);
            }
        }
    }

    // Verify that EVERY boundary edge of EVERY component in this cavity is covered!
    for (&c, bdy_edges) in &comp_boundary_edges {
        for &(u, v) in bdy_edges {
            let key = (c, u.min(v), u.max(v));
            if !covered_boundary_edges.contains(&key) {
                return Err(PlateError::BadConstraint {
                    what: "unclosed interior: cavity has an untied boundary edge without a declared ApertureLoop",
                });
            }
        }
    }

    // Build tied-edge/corner closure triangles and compute enclosed cavity volume.
    let mut closure_tris: Vec<[[f64; 3]; 3]> = Vec::new();
    let mut overlap_area_m2 = 0.0_f64;

    for (&c, bdy_edges) in &comp_boundary_edges {
        let s_c = side_sign_fn(comp_cavity_side[&c]);
        for &(a, b) in bdy_edges {
            // Only build glue-overlap closure quads on tied edges (aperture edges already have cap triangles).
            let is_aperture_edge = open_apertures.iter().any(|ap| {
                ap.cavity_id == cav_id
                    && ap.component == c
                    && (0..ap.nodes.len()).any(|k| {
                        let u = ap.nodes[k];
                        let v = ap.nodes[(k + 1) % ap.nodes.len()];
                        (u.min(v), u.max(v)) == (a.min(b), a.max(b))
                    })
            });
            if is_aperture_edge {
                continue;
            }
            let p_a_raw = add3(
                components[c].mesh.nodes[a],
                scale3(directors_per_comp[c][a], s_c),
            );
            let p_b_raw = add3(
                components[c].mesh.nodes[b],
                scale3(directors_per_comp[c][b], s_c),
            );
            let mut p_a_junc = p_a_raw;
            if let Some(offs) = vertex_transverse_offsets.get(&(c, a)) {
                for &off in offs {
                    p_a_junc = add3(p_a_junc, off);
                }
            }
            let mut p_b_junc = p_b_raw;
            if let Some(offs) = vertex_transverse_offsets.get(&(c, b)) {
                for &off in offs {
                    p_b_junc = add3(p_b_junc, off);
                }
            }

            let tri1 = [p_b_raw, p_a_raw, p_a_junc];
            let tri2 = [p_b_raw, p_a_junc, p_b_junc];
            let area1 = 0.5 * norm3(cross(sub(tri1[1], tri1[0]), sub(tri1[2], tri1[0])));
            let area2 = 0.5 * norm3(cross(sub(tri2[1], tri2[0]), sub(tri2[2], tri2[0])));
            overlap_area_m2 += area1 + area2;
            closure_tris.push(tri1);
            closure_tris.push(tri2);
        }
    }
    closure_tris.extend(aperture_closure_tris);

    // Compute signed enclosed cavity volume from (cavity-outward raw skin triangles + closure triangles).
    let origin = components[0].mesh.nodes[0];
    let mut volume_m3 = 0.0_f64;
    for &face in cavity_faces {
        let panel_idx = match face.side {
            ShellFace::Positive => 2 * face.facet,
            ShellFace::Negative => 2 * face.facet + 1,
        };
        let raw_tri = skins[face.component].triangles()[panel_idx];
        // Reverse winding from wall-outward (cavity-inward) to cavity-outward.
        let cav_tri = [raw_tri[0], raw_tri[2], raw_tri[1]];
        let area_normal = cross(sub(cav_tri[1], cav_tri[0]), sub(cav_tri[2], cav_tri[0]));
        volume_m3 += dot(sub(cav_tri[0], origin), area_normal) / 6.0;
    }
    for tri in &closure_tris {
        let area_normal = cross(sub(tri[1], tri[0]), sub(tri[2], tri[0]));
        volume_m3 += dot(sub(tri[0], origin), area_normal) / 6.0;
    }

    if !volume_m3.is_finite() || volume_m3 <= 0.0 {
        return Err(PlateError::BadConstraint {
            what: "orientation conflict: interior cavity enclosed volume is non-positive",
        });
    }

    Ok((closure_tris, overlap_area_m2, volume_m3))
}

fn push_unique_transverse_offset(
    list: &mut Vec<[f64; 3]>,
    own_director: [f64; 3],
    other_offset: [f64; 3],
) {
    let own_len = norm3(own_director);
    let other_len = norm3(other_offset);
    if own_len <= 1e-15 || other_len <= 1e-15 {
        return;
    }
    let own_unit = scale3(own_director, 1.0 / own_len);
    let other_unit = scale3(other_offset, 1.0 / other_len);
    // Ignore offsets parallel to the panel's own normal.
    if dot(own_unit, other_unit).abs() > 0.99 {
        return;
    }
    for existing in list.iter() {
        let ex_len = norm3(*existing);
        if ex_len > 1e-15 {
            let ex_unit = scale3(*existing, 1.0 / ex_len);
            if dot(ex_unit, other_unit).abs() > 0.99 {
                return;
            }
        }
    }
    list.push(other_offset);
}

fn propagate_corner_transverse_offsets(
    components: &[ShellComponent<'_>],
    comp_cavity_side: &BTreeMap<usize, ShellFace>,
    directors_per_comp: &[Vec<[f64; 3]>],
    edge_ties: &[EdgeTie],
    vertex_transverse_offsets: &mut BTreeMap<(usize, usize), Vec<[f64; 3]>>,
) {
    // Run a fixed number of passes across EdgeTie endpoints to share transverse wall offsets at 3-way corners.
    for _ in 0..components.len() {
        for tie in edge_ties {
            let c_a = tie.a.component;
            let c_b = tie.b.component;
            if !comp_cavity_side.contains_key(&c_a) || !comp_cavity_side.contains_key(&c_b) {
                continue;
            }
            let ends = [
                (tie.a.nodes[0], tie.b.nodes[0]),
                (
                    tie.a.nodes[tie.a.nodes.len() - 1],
                    tie.b.nodes[tie.b.nodes.len() - 1],
                ),
            ];
            for (na, nb) in ends {
                let offs_a = vertex_transverse_offsets
                    .get(&(c_a, na))
                    .cloned()
                    .unwrap_or_default();
                let offs_b = vertex_transverse_offsets
                    .get(&(c_b, nb))
                    .cloned()
                    .unwrap_or_default();
                for off in offs_b {
                    push_unique_transverse_offset(
                        vertex_transverse_offsets.entry((c_a, na)).or_default(),
                        directors_per_comp[c_a][na],
                        off,
                    );
                }
                for off in offs_a {
                    push_unique_transverse_offset(
                        vertex_transverse_offsets.entry((c_b, nb)).or_default(),
                        directors_per_comp[c_b][nb],
                        off,
                    );
                }
            }
        }
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::PlateSection;
    use crate::shell::junction::{RotationCoupling, TieKind, TieOffset, TieSide};
    use crate::shell::reduction::radiation::RadiationSurfaceBudget;
    use crate::shell::reduction::{ReductionBudget, ShellReduction};
    use crate::shell::{ShellMesh, ShellSupport, assemble_shell, modes_shell};
    use fs_modal::SliceOptions;

    struct ClosedBoxFixture {
        meshes: Vec<ShellMesh>,
        models: Vec< crate::shell::ShellModel>,
        thicknesses: Vec<Vec<f64>>,
        skins: Vec<ShellRadiationSurface>,
        ties: Vec<EdgeTie>,
    }

    fn build_closed_box_fixture(a: f64, b: f64, c: f64, h: f64, flipped_panel: Option<usize>) -> ClosedBoxFixture {
        let sec = PlateSection::isotropic(10e9, 0.3, h, 500.0).unwrap();
        // 6 panels with outward normals:
        // 0: Bottom (z=0, normal -z)
        // 1: Top (z=c, normal +z)
        // 2: Front (y=0, normal -y)
        // 3: Back (y=b, normal +y)
        // 4: Left (x=0, normal -x)
        // 5: Right (x=a, normal +x)
        let raw_panels: [(Vec<[f64; 3]>, Vec<[usize; 3]>); 6] = [
            // 0: Bottom z=0, outward normal (0,0,-1)
            (
                vec![[0.0, 0.0, 0.0], [a, 0.0, 0.0], [a, b, 0.0], [0.0, b, 0.0]],
                vec![[0, 2, 1], [0, 3, 2]],
            ),
            // 1: Top z=c, outward normal (0,0,+1)
            (
                vec![[0.0, 0.0, c], [a, 0.0, c], [a, b, c], [0.0, b, c]],
                vec![[0, 1, 2], [0, 2, 3]],
            ),
            // 2: Front y=0, outward normal (0,-1,0)
            (
                vec![[0.0, 0.0, 0.0], [a, 0.0, 0.0], [a, 0.0, c], [0.0, 0.0, c]],
                vec![[0, 1, 2], [0, 2, 3]],
            ),
            // 3: Back y=b, outward normal (0,+1,0)
            (
                vec![[0.0, b, 0.0], [a, b, 0.0], [a, b, c], [0.0, b, c]],
                vec![[0, 2, 1], [0, 3, 2]],
            ),
            // 4: Left x=0, outward normal (-1,0,0)
            (
                vec![[0.0, 0.0, 0.0], [0.0, b, 0.0], [0.0, b, c], [0.0, 0.0, c]],
                vec![[0, 2, 1], [0, 3, 2]],
            ),
            // 5: Right x=a, outward normal (+1,0,0)
            (
                vec![[a, 0.0, 0.0], [a, b, 0.0], [a, b, c], [a, 0.0, c]],
                vec![[0, 1, 2], [0, 2, 3]],
            ),
        ];

        let mut meshes = Vec::new();
        let mut models = Vec::new();
        let mut thicknesses = Vec::new();
        let mut skins = Vec::new();
        let rad_budget = RadiationSurfaceBudget {
            max_panels: 256,
            max_panel_modes: 1024,
        };

        for (idx, (nodes, mut tris)) in raw_panels.into_iter().enumerate() {
            if flipped_panel == Some(idx) {
                for t in &mut tris {
                    t.swap(1, 2);
                }
            }
            let mesh = ShellMesh::new(nodes, tris).unwrap();
            let model = assemble_shell(&mesh, &sec, &[], ShellSupport::Free).unwrap();
            let h_vec = vec![h; mesh.nodes.len()];
            let report = modes_shell(
                &model,
                (10.0, 1.0e10),
                &SliceOptions {
                    ritz_tol: 1e-12,
                    max_lanczos: 64,
                    max_restarts: 32,
                    mass_diagonal_equilibration: true,
                    ..SliceOptions::default()
                },
            )
            .unwrap();
            let reduction = ShellReduction::new(
                &mesh,
                &vec![sec; mesh.tris.len()],
                &model,
                &report.modes[..1],
                ReductionBudget {
                    max_modes: 8,
                    max_facet_modes: 4096,
                    relative_tolerance: 1e-6,
                },
            )
            .unwrap();
            let skin = reduction.radiation_surface(&h_vec, rad_budget).unwrap();
            meshes.push(mesh);
            models.push(model);
            thicknesses.push(h_vec);
            skins.push(skin);
        }

        let tie = |ca: usize, na: [usize; 2], cb: usize, nb: [usize; 2]| EdgeTie {
            a: TieSide {
                component: ca,
                nodes: na.to_vec(),
            },
            b: TieSide {
                component: cb,
                nodes: nb.to_vec(),
            },
            offset: TieOffset::Midsurface,
            rotations: RotationCoupling::Full,
            kind: TieKind::Rigid,
        };

        let ties = vec![
            // Bottom (0) to Front (2), Back (3), Left (4), Right (5)
            tie(0, [0, 1], 2, [0, 1]),
            tie(0, [3, 2], 3, [0, 1]),
            tie(0, [0, 3], 4, [0, 1]),
            tie(0, [1, 2], 5, [0, 1]),
            // Top (1) to Front (2), Back (3), Left (4), Right (5)
            tie(1, [0, 1], 2, [3, 2]),
            tie(1, [3, 2], 3, [3, 2]),
            tie(1, [0, 3], 4, [3, 2]),
            tie(1, [1, 2], 5, [3, 2]),
            // Vertical edges: Front (2) & Back (3) to Left (4) & Right (5)
            tie(2, [0, 3], 4, [0, 3]),
            tie(2, [1, 2], 5, [0, 3]),
            tie(3, [0, 3], 4, [1, 2]),
            tie(3, [1, 2], 5, [1, 2]),
        ];

        ClosedBoxFixture {
            meshes,
            models,
            thicknesses,
            skins,
            ties,
        }
    }

    fn components_from_fixture(f: &ClosedBoxFixture) -> Vec<ShellComponent<'_>> {
        (0..f.meshes.len())
            .map(|i| ShellComponent {
                id: i,
                mesh: &f.meshes[i],
                model: &f.models[i],
                nodal_thickness_m: &f.thicknesses[i],
            })
            .collect()
    }

    #[test]
    fn g0_missing_owner_refuses() {
        let f = build_closed_box_fixture(0.3, 0.2, 0.15, 0.004, None);
        let comps = components_from_fixture(&f);
        let doms = vec![(AcousticDomainId::Exterior, AcousticDomainId::InteriorCavity(0)); 6];
        let mut assignments = FluidBoundaryAssignments::from_component_sides(
            &comps,
            &doms,
            f.ties.clone(),
            Vec::new(),
            Vec::new(),
        );
        // Remove one face's ownership assignment.
        assignments.face_owners.pop();
        let err = FluidBoundaryOwnership::audit(&comps, &f.skins, &assignments).unwrap_err();
        assert!(matches!(err, PlateError::BadConstraint { .. }));
    }

    #[test]
    fn g0_double_owner_refuses() {
        let f = build_closed_box_fixture(0.3, 0.2, 0.15, 0.004, None);
        let comps = components_from_fixture(&f);
        let doms = vec![(AcousticDomainId::Exterior, AcousticDomainId::InteriorCavity(0)); 6];
        let mut assignments = FluidBoundaryAssignments::from_component_sides(
            &comps,
            &doms,
            f.ties.clone(),
            Vec::new(),
            Vec::new(),
        );
        // Duplicate an existing face ownership entry.
        let dup = assignments.face_owners[0];
        assignments.face_owners.push((dup.0, AcousticDomainId::InteriorCavity(1)));
        let err = FluidBoundaryOwnership::audit(&comps, &f.skins, &assignments).unwrap_err();
        assert!(matches!(err, PlateError::BadConstraint { .. }));
    }

    #[test]
    fn g0_flipped_panel_refuses() {
        // Flip panel 2's triangle winding without updating its domain assignment -> orientation conflict across seam.
        let f = build_closed_box_fixture(0.3, 0.2, 0.15, 0.004, Some(2));
        let comps = components_from_fixture(&f);
        let doms = vec![(AcousticDomainId::Exterior, AcousticDomainId::InteriorCavity(0)); 6];
        let assignments = FluidBoundaryAssignments::from_component_sides(
            &comps,
            &doms,
            f.ties.clone(),
            Vec::new(),
            Vec::new(),
        );
        let err = FluidBoundaryOwnership::audit(&comps, &f.skins, &assignments).unwrap_err();
        assert!(matches!(err, PlateError::BadConstraint { .. }));
    }

    #[test]
    fn g0_open_box_without_declared_aperture_refuses() {
        let f = build_closed_box_fixture(0.3, 0.2, 0.15, 0.004, None);
        let comps = components_from_fixture(&f);
        let doms = vec![(AcousticDomainId::Exterior, AcousticDomainId::InteriorCavity(0)); 6];
        // Omit one edge tie so the interior cavity has an undeclared open seam.
        let mut incomplete_ties = f.ties.clone();
        incomplete_ties.pop();
        let assignments = FluidBoundaryAssignments::from_component_sides(
            &comps,
            &doms,
            incomplete_ties,
            Vec::new(),
            Vec::new(),
        );
        let err = FluidBoundaryOwnership::audit(&comps, &f.skins, &assignments).unwrap_err();
        assert!(matches!(err, PlateError::BadConstraint { .. }));
    }
}
