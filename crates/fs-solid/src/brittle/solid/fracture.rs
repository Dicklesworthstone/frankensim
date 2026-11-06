//! Native continuation after an explicitly completed, traction-free fracture.
//!
//! This changes connectivity, not constitutive representation. Both incident
//! bands must already have failed on the selected face normal, and neither may
//! carry compressive contact. Original element strains, source owners, stored
//! energy and cumulative fracture loss survive. No new cohesive spring, fresh
//! surface work, artificial impulse or reset of the damaged bulk is introduced.
//! Each incident volume keeps its own previously accumulated V/A fracture loss;
//! the selected face does not replace that inventory with one fresh Gc*A.
//!
//! Released faces retain a fixed-reference paired opening constraint. Its
//! affine P1 minimum is attained at a vertex, so all three vertex pairs are
//! checked on every accepted mechanical endpoint. Closing endpoints refuse
//! transactionally until a contact response is supplied. This is not global
//! collision detection, finite-sliding contact, or insertion of a partially
//! damaged loaded surface. Explicitly selected conforming faces must align with
//! the completed bands; the original V/A regularization is not recalibrated.

use std::collections::BTreeSet;

use fs_exec::Cx;

use super::{BrittleTetError, BrittleTetSolid, checkpoint, energy_gate};
use crate::brittle::TetCrackBandGeometry;
use crate::brittle::topology::{FractureTopologyError, TetFractureTopology};
use crate::cohesive::p1::frame;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::SmallStrainTetSolution;

/// Explicit numerical gates; none supplies a new strength or fracture law.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetFractureControls {
    /// Chord distance between unit band/face normals, independent of sign.
    pub normal_alignment_tolerance: f64,
    /// Euclidean norm of actual physical traction on each incident side, Pa.
    pub traction_tolerance_pa: f64,
    /// Separate absolute storage and cumulative-fracture transfer gates, J.
    pub energy_tolerance_j: f64,
    /// Roundoff allowance for fixed-reference paired opening, metres.
    pub opening_tolerance_m: f64,
}

/// A previously released face bound to the current element/node numbering.
/// Its normal points from A to B in the unchanged reference configuration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReleasedBrittleFace {
    pub id: u64,
    pub elements: [usize; 2],
    pub paired_nodes: [[usize; 3]; 2],
    pub normal: [f64; 3],
    pub opening_tolerance_m: f64,
}

/// Representation-change accounts. The operation does no physical work.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BrittleTetFractureEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub fracture_dissipation_before_j: f64,
    pub fracture_dissipation_after_j: f64,
    pub storage_defect_j: f64,
    pub fracture_history_defect_j: f64,
}

/// Complete native solid on the admitted split geometry.
#[derive(Debug)]
pub struct BrittleTetFractureStep<'a> {
    pub solid: BrittleTetSolid<'a>,
    pub energy: BrittleTetFractureEnergy,
}

#[derive(Clone, Debug, PartialEq)]
pub enum BrittleTetFractureError {
    Invalid(&'static str),
    Topology(FractureTopologyError),
    Solid(BrittleTetError),
}
impl core::fmt::Display for BrittleTetFractureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "completed brittle fracture transfer refused: {self:?}")
    }
}
impl std::error::Error for BrittleTetFractureError {}
impl From<BrittleTetError> for BrittleTetFractureError {
    fn from(value: BrittleTetError) -> Self {
        Self::Solid(value)
    }
}
impl From<FractureTopologyError> for BrittleTetFractureError {
    fn from(value: FractureTopologyError) -> Self {
        Self::Topology(value)
    }
}

impl<'a> BrittleTetSolid<'a> {
    /// Retained released traces for observation and an explicit contact owner.
    #[must_use]
    pub fn released_faces(&self) -> &[ReleasedBrittleFace] {
        &self.released_faces
    }

    /// Release explicitly selected, completed fractures without losing history.
    ///
    /// Construct `topology` from this solid's CURRENT reference mesh, then admit
    /// its output with the original load-geometry owner. Element/local-corner
    /// order and reference coordinates must be unchanged. Both incident native
    /// material histories must admit the face as failed and traction-free.
    /// Previous releases are remapped by their owning element's local corners;
    /// seam IDs remain unique through subsequent cuts.
    ///
    /// Material/source ownership, strain, stress and irreversible memory stay
    /// exact. The latest equilibrium record is cleared because its old global
    /// DOF/reaction numbering no longer describes the receiving mesh. Boundary
    /// conditions and applied loads must be supplied in the new node numbering.
    pub fn split_fracture<'next>(
        &self,
        topology: &TetFractureTopology<'_>,
        target_geometry: &'next TetLoadGeometry3<'next>,
        controls: BrittleTetFractureControls,
        cx: &Cx<'_>,
    ) -> Result<BrittleTetFractureStep<'next>, BrittleTetFractureError>
    where
        'a: 'next,
    {
        checkpoint(cx)?;
        if !controls.normal_alignment_tolerance.is_finite()
            || !(0.0..1.0).contains(&controls.normal_alignment_tolerance)
            || [
                controls.traction_tolerance_pa,
                controls.energy_tolerance_j,
                controls.opening_tolerance_m,
            ]
            .iter()
            .any(|x| !x.is_finite() || *x < 0.0)
        {
            return Err(BrittleTetFractureError::Invalid(
                "finite nonnegative release gates required",
            ));
        }
        self.require_split_geometry(topology, target_geometry, cx)?;
        self.require_released_faces(&self.solution.displacement_m, cx)?;
        let face_count = self
            .released_faces
            .len()
            .checked_add(topology.seams().len())
            .ok_or(BrittleTetFractureError::Invalid(
                "released face count overflow",
            ))?;
        if topology.seams().is_empty() || face_count > self.materials.len().saturating_mul(2) {
            return Err(BrittleTetFractureError::Invalid(
                "nonempty bounded new fracture faces required",
            ));
        }
        let mut released_faces = Vec::with_capacity(face_count);
        let mut ids = BTreeSet::new();
        for old in self.released_faces.iter() {
            checkpoint(cx)?;
            let mut next = *old;
            for side in 0..2 {
                let element = old.elements[side];
                for vertex in 0..3 {
                    let corner = self.geometry.tetrahedra()[element]
                        .iter()
                        .position(|node| *node == old.paired_nodes[side][vertex])
                        .ok_or(BrittleTetFractureError::Invalid(
                            "retained face lost its element corner",
                        ))?;
                    next.paired_nodes[side][vertex] = target_geometry.tetrahedra()[element][corner];
                }
            }
            ids.insert(next.id);
            released_faces.push(next);
        }
        for seam in topology.seams() {
            checkpoint(cx)?;
            if !ids.insert(seam.id) {
                return Err(BrittleTetFractureError::Invalid(
                    "fracture face ID was already released",
                ));
            }
            let opposite = |side: usize, nodes: [usize; 3]| {
                target_geometry.tetrahedra()[seam.elements[side]]
                    .iter()
                    .copied()
                    .find(|node| !nodes.contains(node))
                    .ok_or(BrittleTetFractureError::Invalid(
                        "missing opposite fracture vertex",
                    ))
            };
            let (axes, _) = frame(
                target_geometry.nodes_m(),
                seam.side_a,
                opposite(0, seam.side_a)?,
                opposite(1, seam.side_b)?,
                0.0,
            )
            .map_err(|_| BrittleTetFractureError::Invalid("invalid native fracture face frame"))?;
            for element in seam.elements {
                let released = self.materials[element]
                    .law()
                    .release_failed_plane(
                        &self.solution.element_states[element],
                        axes[0],
                        controls.normal_alignment_tolerance,
                        controls.traction_tolerance_pa,
                    )
                    .map_err(|source| BrittleTetError::Material { element, source })?;
                if released.stress_pa != self.solution.element_stress[element] {
                    return Err(BrittleTetFractureError::Invalid(
                        "release changed accepted native stress",
                    ));
                }
            }
            released_faces.push(ReleasedBrittleFace {
                id: seam.id,
                elements: seam.elements,
                paired_nodes: [seam.side_a, seam.side_b],
                normal: axes[0],
                opening_tolerance_m: controls.opening_tolerance_m,
            });
        }
        let mut displacement_m = Vec::with_capacity(target_geometry.nodes_m().len());
        for (index, &node) in topology.original_nodes().iter().enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            displacement_m.push(self.solution.displacement_m[node]);
        }
        let mut bands = Vec::with_capacity(self.bands.len());
        let mut stored = 0.0;
        let mut loss = 0.0;
        for (element, source) in self.materials.iter().enumerate() {
            checkpoint(cx)?;
            let band =
                TetCrackBandGeometry::from_element(target_geometry, element, source.law().normal())
                    .map_err(BrittleTetError::Geometry)?;
            let old = self.bands[element];
            if band.reference_volume_m3() != old.reference_volume_m3()
                || band.projected_area_m2() != old.projected_area_m2()
                || band.band_width_m() != old.band_width_m()
            {
                return Err(BrittleTetFractureError::Invalid(
                    "release changed the original band geometry",
                ));
            }
            let state = &self.solution.element_states[element];
            stored = band
                .reference_volume_m3()
                .mul_add(state.stored_energy_j_m3(), stored);
            loss = band
                .reference_volume_m3()
                .mul_add(state.dissipated_j_m3(), loss);
            bands.push(band);
        }
        let energy = BrittleTetFractureEnergy {
            stored_before_j: self.stored_energy_j,
            stored_after_j: stored,
            fracture_dissipation_before_j: self.fracture_dissipation_j,
            fracture_dissipation_after_j: loss,
            storage_defect_j: stored - self.stored_energy_j,
            fracture_history_defect_j: loss - self.fracture_dissipation_j,
        };
        energy_gate(
            "fracture topology storage",
            energy.storage_defect_j,
            controls.energy_tolerance_j,
        )?;
        energy_gate(
            "fracture topology irreversible history",
            energy.fracture_history_defect_j,
            controls.energy_tolerance_j,
        )?;
        let solid = BrittleTetSolid {
            geometry: target_geometry,
            materials: self.materials,
            bands: bands.into(),
            released_faces: released_faces.into(),
            budget: self.budget,
            solution: SmallStrainTetSolution {
                displacement_m,
                element_strain: self.solution.element_strain.clone(),
                element_stress: self.solution.element_stress.clone(),
                element_states: self.solution.element_states.clone(),
                increments: Vec::new(),
            },
            stored_energy_j: stored,
            fracture_dissipation_j: loss,
        };
        solid.require_released_faces(&solid.solution.displacement_m, cx)?;
        checkpoint(cx)?;
        Ok(BrittleTetFractureStep { solid, energy })
    }

    fn require_split_geometry(
        &self,
        topology: &TetFractureTopology<'_>,
        target: &TetLoadGeometry3<'_>,
        cx: &Cx<'_>,
    ) -> Result<(), BrittleTetFractureError> {
        let old_nodes = self.geometry.nodes_m();
        let map = topology.original_nodes();
        if target.nodes_m().len() > self.budget.maximum_nodes
            || target.nodes_m().len().checked_mul(3).is_none()
            || target.tetrahedra().len() != self.materials.len()
            || target.nodes_m().len() != topology.nodes_m().len()
            || target.tetrahedra().len() != topology.tetrahedra().len()
            || map.len() != target.nodes_m().len()
            || map.len() < old_nodes.len()
        {
            return Err(BrittleTetFractureError::Invalid(
                "receiving mesh must be the bounded native split mesh",
            ));
        }
        for (next, &previous) in map.iter().enumerate() {
            if next % 256 == 0 {
                checkpoint(cx)?;
            }
            if previous >= old_nodes.len()
                || (next < old_nodes.len() && next != previous)
                || target.nodes_m()[next] != topology.nodes_m()[next]
                || target.nodes_m()[next] != old_nodes[previous]
            {
                return Err(BrittleTetFractureError::Invalid(
                    "split changed reference nodes or origin map",
                ));
            }
        }
        for (element, tet) in target.tetrahedra().iter().enumerate() {
            checkpoint(cx)?;
            if *tet != topology.tetrahedra()[element]
                || tet.map(|node| map[node]) != self.geometry.tetrahedra()[element]
                || target.element_volumes_m3()[element]
                    != self.geometry.element_volumes_m3()[element]
            {
                return Err(BrittleTetFractureError::Invalid(
                    "split changed element order, local corners or volume",
                ));
            }
        }
        Ok(())
    }

    /// Fixed-reference matching trace guard, shared by every native integrator.
    pub(crate) fn require_released_faces(
        &self,
        displacement_m: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<(), BrittleTetError> {
        if displacement_m.len() != self.geometry.nodes_m().len() {
            return Err(BrittleTetError::Invalid(
                "released-face displacement dimension",
            ));
        }
        for face in self.released_faces.iter() {
            checkpoint(cx)?;
            for vertex in 0..3 {
                let a = displacement_m[face.paired_nodes[0][vertex]];
                let b = displacement_m[face.paired_nodes[1][vertex]];
                let mut opening = 0.0;
                for component in 0..3 {
                    opening = face.normal[component].mul_add(b[component] - a[component], opening);
                }
                if !opening.is_finite() || opening < -face.opening_tolerance_m {
                    return Err(BrittleTetError::ReleasedFaceClosure {
                        face: face.id,
                        vertex,
                        opening_m: opening,
                        tolerance_m: face.opening_tolerance_m,
                    });
                }
            }
        }
        Ok(())
    }
}
