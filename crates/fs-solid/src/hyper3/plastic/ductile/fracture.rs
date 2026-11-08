//! Native finite ductile failure-to-separation on explicitly selected P1 faces.
//!
//! This operation only duplicates original nodal stars. It preserves each
//! material owner, multiplicative metric, initiation and work coordinate,
//! deformation, physical stress, directional length, loss history and mass.
//! Both incident points must reach an explicitly finite-support terminal law;
//! rounded damage or exponential underflow cannot authorize a crack.
//! Newly released faces carry no compressive storage and no accepted traction.
//! No extra surface energy, impulse, constitutive return or physical time is
//! introduced. The resulting deformable fragments retain native finite inertia.
//!
//! Contact uses normals frozen from the CURRENT release geometry. Continuation
//! outside the admitted face-normal cone refuses; this is a matched-trace seam
//! owner, not arbitrary rotating-fragment collision detection.

use core::fmt::Debug;
use std::collections::BTreeSet;
use std::sync::Arc;

use fs_exec::Cx;
use fs_material::ductile::evolution::DifferentiableDuctileInitiationLaw;
use fs_material::plastic::finite::ductile::FiniteDuctileEnergy;
use fs_math::det;

use super::dynamics::{
    FiniteDuctileMotionControls, FiniteDuctileMotionStep, FiniteDuctileTetMotion,
};
use super::{
    FiniteDuctileTetEnergy, FiniteDuctileTetError, FiniteDuctileTetSolid, FinitePlasticTetState,
    FrozenPoint, HyperTetProblem, HyperelasticLawField, checkpoint,
};
use crate::brittle::TetCrackBandGeometry;
use crate::brittle::topology::{FractureTopologyError, TetFractureTopology};
use crate::ductile::DuctileStepConditions;
use crate::inertia::InertialEnergy;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::{IncrementRecord, IncrementSettings, SmallStrainTetError};

mod transfer;

/// Independent absolute release gates; none supplies a material calibration.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctileFractureControls {
    /// Maximum accepted physical traction magnitude on either current trace.
    pub traction_tolerance_pa: f64,
    /// Maximum absolute native versus integrated kinetic/storage roundoff, J.
    pub energy_tolerance_j: f64,
    /// Permitted paired-opening roundoff for contact-free continuation, m.
    pub opening_tolerance_m: f64,
    /// Maximum 1-cos(angle) for authored reference-band and current trace normals.
    pub normal_alignment_tolerance: f64,
}

/// Original face identity, retained incident cells and paired local witnesses.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct ReleasedDuctileFace {
    pub id: u64,
    pub elements: [usize; 2],
    /// Corresponding original vertices, oriented A-to-B in the release geometry.
    pub paired_nodes: [[usize; 3]; 2],
    /// Unit normal of the CURRENT face when released, fixed during continuation.
    pub normal: [f64; 3],
}

/// Consistent P1 integrals for constant-cell or retained P1 nodal density.
/// First and angular moments use current positions about the coordinate origin.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct FiniteDuctileFractureInventory {
    pub mass_kg: f64,
    pub first_moment_kg_m: [f64; 3],
    pub linear_momentum_kg_m_s: [f64; 3],
    pub angular_momentum_kg_m2_s: [f64; 3],
    pub kinetic_energy_j: f64,
}

/// Native work histories before and after a zero-duration connectivity change.
/// These cumulative physical losses are retained, never credited a second time.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct FiniteDuctileFractureEnergy {
    pub stored_before_j: f64,
    pub stored_after_j: f64,
    pub plastic_dissipation_before_j: f64,
    pub plastic_dissipation_after_j: f64,
    pub damage_dissipation_before_j: f64,
    pub damage_dissipation_after_j: f64,
    pub endpoint_defect_before_j: f64,
    pub endpoint_defect_after_j: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    pub kinetic_defect_j: f64,
}

/// Native deformable motion with persistent released-face identity and guards.
#[derive(Debug)]
pub struct FiniteDuctileFractureMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub(super) motion: FiniteDuctileTetMotion<'g, 'm, L>,
    pub(super) released_faces: Arc<[ReleasedDuctileFace]>,
    pub(super) opening_tolerance_m: f64,
    pub(super) normal_alignment_tolerance: f64,
}

/// Atomic topology handoff, including independently integrated matter inventory.
#[derive(Debug)]
pub struct FiniteDuctileFractureMotionStep<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub motion: FiniteDuctileFractureMotion<'g, 'm, L>,
    pub inventory_before: FiniteDuctileFractureInventory,
    pub inventory_after: FiniteDuctileFractureInventory,
    /// Reference-element connected components in the topology owner's order.
    pub fragments: Vec<FiniteDuctileFractureInventory>,
    pub energy: FiniteDuctileFractureEnergy,
}

/// One accepted native physical interval retaining the released-face owner.
#[derive(Debug)]
pub struct FiniteDuctileFractureAdvance<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub motion: FiniteDuctileFractureMotion<'g, 'm, L>,
    pub equilibrium: IncrementRecord,
    pub material_energy: FiniteDuctileTetEnergy,
    pub element_energy: Vec<FiniteDuctileEnergy>,
    pub inertia: InertialEnergy,
    pub balance_error_j: f64,
}

#[derive(Debug)]
pub enum FiniteDuctileFractureError {
    Invalid(&'static str),
    Topology(FractureTopologyError),
    Native(FiniteDuctileTetError),
    Face {
        face: u64,
        element: usize,
        reason: &'static str,
    },
    Energy {
        quantity: &'static str,
        value_j: f64,
        maximum_abs_j: f64,
    },
}
impl core::fmt::Display for FiniteDuctileFractureError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "finite ductile separation refused: {self:?}")
    }
}
impl std::error::Error for FiniteDuctileFractureError {}
impl From<FiniteDuctileTetError> for FiniteDuctileFractureError {
    fn from(error: FiniteDuctileTetError) -> Self {
        Self::Native(error)
    }
}
impl From<SmallStrainTetError> for FiniteDuctileFractureError {
    fn from(error: SmallStrainTetError) -> Self {
        Self::Native(error.into())
    }
}
impl From<FractureTopologyError> for FiniteDuctileFractureError {
    fn from(error: FractureTopologyError) -> Self {
        Self::Topology(error)
    }
}

impl<'g, 'm, L> FiniteDuctileTetMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    /// Split explicitly selected fully failed conforming faces, without taking a
    /// time step. Target geometry and bands must be exact node reindexings of the
    /// same original cells; caller-owned laws and directional lengths stay fixed.
    #[allow(clippy::too_many_arguments)]
    pub fn split_fracture<'next>(
        &self,
        topology: &TetFractureTopology<'_>,
        target_geometry: &'next TetLoadGeometry3<'next>,
        target_bands: &'next [TetCrackBandGeometry],
        controls: FiniteDuctileFractureControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileFractureMotionStep<'next, 'm, L>, FiniteDuctileFractureError> {
        split_motion(
            self,
            &[],
            topology,
            target_geometry,
            target_bands,
            controls,
            cx,
        )
    }
}

impl<'g, 'm, L> FiniteDuctileFractureMotion<'g, 'm, L>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    pub fn motion(&self) -> &FiniteDuctileTetMotion<'g, 'm, L> {
        &self.motion
    }
    pub fn released_faces(&self) -> &[ReleasedDuctileFace] {
        &self.released_faces
    }
    pub const fn opening_tolerance_m(&self) -> f64 {
        self.opening_tolerance_m
    }
    pub const fn normal_alignment_tolerance(&self) -> f64 {
        self.normal_alignment_tolerance
    }

    /// Add new explicit cuts while preserving all already released face IDs,
    /// current contact normals and original incident-cell/local-corner lineage.
    #[allow(clippy::too_many_arguments)]
    pub fn split_fracture<'next>(
        &self,
        topology: &TetFractureTopology<'_>,
        target_geometry: &'next TetLoadGeometry3<'next>,
        target_bands: &'next [TetCrackBandGeometry],
        controls: FiniteDuctileFractureControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileFractureMotionStep<'next, 'm, L>, FiniteDuctileFractureError> {
        self.admit_endpoint(self.motion.solid.displacement_m(), cx)?;
        if controls.opening_tolerance_m > self.opening_tolerance_m
            || controls.normal_alignment_tolerance > self.normal_alignment_tolerance
        {
            return Err(FiniteDuctileFractureError::Invalid(
                "additional cuts cannot relax retained released-face gates",
            ));
        }
        split_motion(
            &self.motion,
            &self.released_faces,
            topology,
            target_geometry,
            target_bands,
            controls,
            cx,
        )
    }

    /// Native finite backward Euler without applied seam contact. Every
    /// published endpoint retains all trace-opening and current-normal gates.
    #[allow(clippy::too_many_arguments)]
    pub fn advance(
        &self,
        prescribed_m: &[(usize, f64)],
        nodal_forces_n: &[(usize, f64)],
        body_force_n_m3: [f64; 3],
        conditions: DuctileStepConditions<'_>,
        settings: IncrementSettings,
        controls: FiniteDuctileMotionControls,
        cx: &Cx<'_>,
    ) -> Result<FiniteDuctileFractureAdvance<'g, 'm, L>, FiniteDuctileFractureError> {
        self.admit_endpoint(self.motion.solid.displacement_m(), cx)?;
        let step = self.motion.advance(
            prescribed_m,
            nodal_forces_n,
            body_force_n_m3,
            conditions,
            settings,
            controls,
            cx,
        )?;
        self.admit_endpoint(step.motion.solid.displacement_m(), cx)?;
        checkpoint(cx)?;
        Ok(self.retain_step(step))
    }

    fn retain_step(
        &self,
        step: FiniteDuctileMotionStep<'g, 'm, L>,
    ) -> FiniteDuctileFractureAdvance<'g, 'm, L> {
        FiniteDuctileFractureAdvance {
            motion: Self {
                motion: step.motion,
                released_faces: Arc::clone(&self.released_faces),
                opening_tolerance_m: self.opening_tolerance_m,
                normal_alignment_tolerance: self.normal_alignment_tolerance,
            },
            equilibrium: step.equilibrium,
            material_energy: step.material_energy,
            element_energy: step.element_energy,
            inertia: step.inertia,
            balance_error_j: step.balance_error_j,
        }
    }

    fn admit_endpoint(
        &self,
        displacement_m: &[[f64; 3]],
        cx: &Cx<'_>,
    ) -> Result<(), FiniteDuctileFractureError> {
        admit_faces(
            self.motion.solid.reference_nodes_m(),
            &self.released_faces,
            displacement_m,
            self.opening_tolerance_m,
            self.normal_alignment_tolerance,
            cx,
        )
    }
}

#[allow(clippy::too_many_arguments, clippy::too_many_lines)]
fn split_motion<'next, 'm, L>(
    old: &FiniteDuctileTetMotion<'_, 'm, L>,
    retained_faces: &[ReleasedDuctileFace],
    topology: &TetFractureTopology<'_>,
    target: &'next TetLoadGeometry3<'next>,
    bands: &'next [TetCrackBandGeometry],
    controls: FiniteDuctileFractureControls,
    cx: &Cx<'_>,
) -> Result<FiniteDuctileFractureMotionStep<'next, 'm, L>, FiniteDuctileFractureError>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    use FiniteDuctileFractureError as E;
    checkpoint(cx)?;
    if [
        controls.traction_tolerance_pa,
        controls.energy_tolerance_j,
        controls.opening_tolerance_m,
    ]
    .iter()
    .any(|x| !x.is_finite() || *x < 0.0)
        || !controls.normal_alignment_tolerance.is_finite()
        || !(0.0..1.0).contains(&controls.normal_alignment_tolerance)
    {
        return Err(E::Invalid(
            "finite nonnegative release gates and a bounded normal cone required",
        ));
    }
    let solid = &old.solid;
    let map = topology.original_nodes();
    if topology.seams().is_empty()
        || target.nodes_m().len() > solid.state.budget.maximum_nodes
        || target.nodes_m().len().checked_mul(3).is_none()
        || target.nodes_m() != topology.nodes_m()
        || target.tetrahedra() != topology.tetrahedra()
        || target.tetrahedra().len() != solid.state.elements.len()
        || bands.len() != solid.bands.len()
        || map.len() != target.nodes_m().len()
        || map.len() < solid.state.nodes_m.len()
    {
        return Err(E::Invalid(
            "a nonempty bounded original-cell split and complete target bands required",
        ));
    }
    for (next, &previous) in map.iter().enumerate() {
        if next % 256 == 0 {
            checkpoint(cx)?;
        }
        if previous >= solid.state.nodes_m.len()
            || (next < solid.state.nodes_m.len() && next != previous)
            || target.nodes_m()[next] != solid.state.nodes_m[previous]
        {
            return Err(E::Invalid(
                "split changed reference coordinates or original node identity",
            ));
        }
    }
    for (element, tet) in target.tetrahedra().iter().enumerate() {
        checkpoint(cx)?;
        if tet.map(|node| map[node]) != solid.state.tetrahedra[element]
            || target.element_volumes_m3()[element] != solid.state.elements[element].volume
        {
            return Err(E::Invalid(
                "split changed element order, local corners or reference volume",
            ));
        }
        let previous = solid.bands[element];
        let band = bands[element];
        let admitted = TetCrackBandGeometry::from_element(target, element, previous.normal())
            .map_err(FiniteDuctileTetError::Band)?;
        if band.element() != element
            || band.nodes() != *tet
            || band.normal() != previous.normal()
            || band.reference_volume_m3() != previous.reference_volume_m3()
            || band.projected_area_m2() != previous.projected_area_m2()
            || band.band_width_m() != previous.band_width_m()
            || admitted.nodes() != band.nodes()
            || admitted.reference_volume_m3() != band.reference_volume_m3()
            || (admitted.projected_area_m2() / band.projected_area_m2() - 1.0).abs()
                > 64.0 * f64::EPSILON
            || solid.materials[element].settings().characteristic_length_m != band.band_width_m()
        {
            return Err(E::Invalid(
                "target bands must retain the exact original direction, width and work owner",
            ));
        }
        let state = &solid.state.element_states[element];
        if !solid.materials[element].owns_state(state)
            || !state.is_valid()
            || state.deformation() != &solid.state.element_deformation[element]
            || state.first_piola_pa() != &solid.state.element_first_piola_pa[element]
        {
            return Err(E::Invalid(
                "original accepted finite history and physical stress must stay bound",
            ));
        }
    }
    let displacement: Vec<_> = map
        .iter()
        .map(|&node| solid.state.displacement_m[node])
        .collect();
    let velocity: Vec<_> = map.iter().map(|&node| old.velocity_mps[node]).collect();
    let mut faces = Vec::with_capacity(retained_faces.len() + topology.seams().len());
    let mut ids = BTreeSet::new();
    for face in retained_faces {
        checkpoint(cx)?;
        let mut next = *face;
        for side in 0..2 {
            let element = face.elements[side];
            for vertex in 0..3 {
                let corner = solid.state.tetrahedra[element]
                    .iter()
                    .position(|node| *node == face.paired_nodes[side][vertex])
                    .ok_or(E::Invalid(
                        "retained face lost its original incident-cell corner",
                    ))?;
                next.paired_nodes[side][vertex] = target.tetrahedra()[element][corner];
            }
        }
        if !ids.insert(next.id) {
            return Err(E::Invalid("duplicate retained face ID"));
        }
        faces.push(next);
    }
    for seam in topology.seams() {
        checkpoint(cx)?;
        if !ids.insert(seam.id) {
            return Err(E::Invalid("a released face ID cannot be reused"));
        }
        let mut paired = [seam.side_a, seam.side_b];
        let opposite = |side: usize| {
            target.tetrahedra()[seam.elements[side]]
                .iter()
                .copied()
                .find(|node| !paired[side].contains(node))
                .ok_or(E::Invalid(
                    "selected face must retain its incident opposite vertex",
                ))
        };
        let across = difference(target.nodes_m(), &displacement, opposite(0)?, opposite(1)?);
        let (mut normal, area) = face_normal(target.nodes_m(), &displacement, paired[0])?;
        if dot(normal, across) < 0.0 {
            paired[0].swap(1, 2);
            paired[1].swap(1, 2);
            normal = normal.map(|x| -x);
        }
        if dot(normal, across) <= 0.0 {
            return Err(E::Invalid(
                "current facet must separate its original incident cells",
            ));
        }
        let reference_edges = [
            subtract(
                target.nodes_m()[paired[0][1]],
                target.nodes_m()[paired[0][0]],
            ),
            subtract(
                target.nodes_m()[paired[0][2]],
                target.nodes_m()[paired[0][0]],
            ),
        ];
        let reference_cross = cross(reference_edges[0], reference_edges[1]);
        let reference_norm = magnitude(reference_cross)?;
        let reference_normal = reference_cross.map(|x| x / reference_norm);
        let reference_area = 0.5 * reference_norm;
        for element in seam.elements {
            let fail = |reason| E::Face {
                face: seam.id,
                element,
                reason,
            };
            let state = &solid.state.element_states[element];
            if !state.is_fully_failed() {
                return Err(fail(
                    "explicit finite-support work completion is required on both incident cells",
                ));
            }
            if state.compression_stored_energy_j_m3() != 0.0 {
                return Err(fail(
                    "a newly released face cannot discard retained compressive support",
                ));
            }
            if 1.0 - dot(reference_normal, solid.bands[element].normal()).abs()
                > controls.normal_alignment_tolerance
            {
                return Err(fail(
                    "selected reference facet differs from the authored crack-band direction",
                ));
            }
            let j = determinant(state.deformation());
            if !j.is_finite() || j <= 0.0 {
                return Err(fail("positive accepted finite volume ratio required"));
            }
            let traction: [f64; 3] = core::array::from_fn(|i| {
                (0..3).fold(0.0, |s, k| {
                    state.kirchhoff_pa()[3 * i + k].mul_add(normal[k], s)
                }) / j
            });
            let piola_traction: [f64; 3] = core::array::from_fn(|i| {
                (0..3).fold(0.0, |s, k| {
                    state.first_piola_pa()[3 * i + k].mul_add(reference_normal[k], s)
                }) * (reference_area / area)
            });
            if magnitude(traction)? > controls.traction_tolerance_pa
                || magnitude(piola_traction)? > controls.traction_tolerance_pa
            {
                return Err(fail("accepted physical facet traction is not released"));
            }
        }
        faces.push(ReleasedDuctileFace {
            id: seam.id,
            elements: seam.elements,
            paired_nodes: paired,
            normal,
        });
    }
    let frozen: Vec<_> = solid
        .materials
        .iter()
        .zip(&solid.state.element_states)
        .map(|(model, state)| FrozenPoint {
            law: model.native(),
            accepted: state.native_state(),
        })
        .collect();
    let elements = FinitePlasticTetState::reference_elements(
        target.nodes_m(),
        target.tetrahedra(),
        HyperelasticLawField::PerElement(&frozen),
        solid.state.budget,
        cx,
    )
    .map_err(FiniteDuctileTetError::from)?;
    for (element, (new, previous)) in elements.iter().zip(solid.state.elements.iter()).enumerate() {
        checkpoint(cx)?;
        if new.volume != previous.volume
            || new.gradients != previous.gradients
            || HyperTetProblem::<FrozenPoint<'_>>::deformation_gradient(
                new,
                displacement.as_flattened(),
            ) != solid.state.element_deformation[element]
        {
            return Err(E::Invalid(
                "reindexed native kinematics must equal the accepted original deformation",
            ));
        }
    }
    let transferred = transfer::motion(
        old,
        topology,
        target,
        &displacement,
        &velocity,
        controls.energy_tolerance_j,
        cx,
    )?;
    let before = histories(solid)?;
    let next_solid = FiniteDuctileTetSolid {
        materials: solid.materials,
        bands,
        state: FinitePlasticTetState {
            nodes_m: target.nodes_m(),
            tetrahedra: target.tetrahedra(),
            elements,
            budget: solid.state.budget,
            displacement_m: displacement,
            element_deformation: solid.state.element_deformation.clone(),
            element_states: solid.state.element_states.clone(),
            element_first_piola_pa: solid.state.element_first_piola_pa.clone(),
            stored_energy_j: solid.state.stored_energy_j,
        },
    };
    let after = histories(&next_solid)?;
    if before != after {
        return Err(E::Invalid(
            "a connectivity change cannot alter any native material work history",
        ));
    }
    energy_gate(
        "accepted storage inventory",
        before[0] - solid.state.stored_energy_j,
        controls.energy_tolerance_j,
    )?;
    let energy = FiniteDuctileFractureEnergy {
        stored_before_j: solid.state.stored_energy_j,
        stored_after_j: next_solid.state.stored_energy_j,
        plastic_dissipation_before_j: before[1],
        plastic_dissipation_after_j: after[1],
        damage_dissipation_before_j: before[2],
        damage_dissipation_after_j: after[2],
        endpoint_defect_before_j: before[3],
        endpoint_defect_after_j: after[3],
        kinetic_before_j: transferred.kinetic_before_j,
        kinetic_after_j: transferred.kinetic_after_j,
        kinetic_defect_j: transferred.kinetic_after_j - transferred.kinetic_before_j,
    };
    let motion = FiniteDuctileFractureMotion {
        motion: FiniteDuctileTetMotion {
            solid: next_solid,
            velocity_mps: velocity,
            reference_density_kg_m3: Arc::clone(&old.reference_density_kg_m3),
            reference_nodal_density_kg_m3: transferred.nodal_density,
            mass: Arc::new(transferred.mass),
        },
        released_faces: faces.into(),
        opening_tolerance_m: controls.opening_tolerance_m,
        normal_alignment_tolerance: controls.normal_alignment_tolerance,
    };
    motion.admit_endpoint(motion.motion.solid.displacement_m(), cx)?;
    checkpoint(cx)?;
    Ok(FiniteDuctileFractureMotionStep {
        motion,
        inventory_before: transferred.before,
        inventory_after: transferred.after,
        fragments: transferred.fragments,
        energy,
    })
}

fn histories<L>(
    solid: &FiniteDuctileTetSolid<'_, '_, L>,
) -> Result<[f64; 4], FiniteDuctileFractureError>
where
    L: DifferentiableDuctileInitiationLaw + Clone + Debug + PartialEq,
{
    let mut values = [0.0; 4];
    for (element, state) in solid.state.elements.iter().zip(&solid.state.element_states) {
        for (sum, value) in values.iter_mut().zip([
            state.stored_energy_j_m3(),
            state.dissipated_plastic_j_m3(),
            state.dissipated_damage_j_m3(),
            state.endpoint_work_defect_j_m3(),
        ]) {
            *sum = element.volume.mul_add(value, *sum);
        }
    }
    if values.iter().any(|v| !v.is_finite()) {
        return Err(FiniteDuctileFractureError::Invalid(
            "finite retained material inventories required",
        ));
    }
    Ok(values)
}

fn admit_faces(
    nodes: &[[f64; 3]],
    faces: &[ReleasedDuctileFace],
    u: &[[f64; 3]],
    opening_tolerance_m: f64,
    alignment_tolerance: f64,
    cx: &Cx<'_>,
) -> Result<(), FiniteDuctileFractureError> {
    if u.len() != nodes.len() || u.iter().flatten().any(|v| !v.is_finite()) {
        return Err(FiniteDuctileFractureError::Invalid(
            "complete finite released-face displacement required",
        ));
    }
    for face in faces {
        checkpoint(cx)?;
        for side in 0..2 {
            let (normal, _) = face_normal(nodes, u, face.paired_nodes[side])?;
            if 1.0 - dot(normal, face.normal) > alignment_tolerance {
                return Err(FiniteDuctileFractureError::Face {
                    face: face.id,
                    element: face.elements[side],
                    reason: "current trace left its admitted fixed release-normal cone",
                });
            }
        }
        for vertex in 0..3 {
            let [a, b] = [face.paired_nodes[0][vertex], face.paired_nodes[1][vertex]];
            let gap = dot(face.normal, subtract(u[b], u[a]));
            if !gap.is_finite() || gap < -opening_tolerance_m {
                return Err(FiniteDuctileFractureError::Face {
                    face: face.id,
                    element: face.elements[0],
                    reason: "released traces interpenetrate without an admitted contact response",
                });
            }
        }
    }
    Ok(())
}

fn energy_gate(
    quantity: &'static str,
    value_j: f64,
    maximum_abs_j: f64,
) -> Result<(), FiniteDuctileFractureError> {
    if !value_j.is_finite() || value_j.abs() > maximum_abs_j {
        Err(FiniteDuctileFractureError::Energy {
            quantity,
            value_j,
            maximum_abs_j,
        })
    } else {
        Ok(())
    }
}
fn difference(nodes: &[[f64; 3]], u: &[[f64; 3]], a: usize, b: usize) -> [f64; 3] {
    core::array::from_fn(|i| (nodes[b][i] - nodes[a][i]) + (u[b][i] - u[a][i]))
}
fn subtract(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    core::array::from_fn(|i| a[i] - b[i])
}
fn dot(a: [f64; 3], b: [f64; 3]) -> f64 {
    a[0].mul_add(b[0], a[1].mul_add(b[1], a[2] * b[2]))
}
fn cross(a: [f64; 3], b: [f64; 3]) -> [f64; 3] {
    [
        a[1].mul_add(b[2], -a[2] * b[1]),
        a[2].mul_add(b[0], -a[0] * b[2]),
        a[0].mul_add(b[1], -a[1] * b[0]),
    ]
}
fn magnitude(x: [f64; 3]) -> Result<f64, FiniteDuctileFractureError> {
    if x.iter().any(|v| !v.is_finite()) {
        return Err(FiniteDuctileFractureError::Invalid(
            "finite physical vector required",
        ));
    }
    let scale = x.into_iter().fold(0.0_f64, |s, v| s.max(v.abs()));
    if scale == 0.0 {
        return Ok(0.0);
    }
    let normalized = x.map(|v| v / scale);
    let value = scale * det::sqrt(dot(normalized, normalized));
    if !value.is_finite() {
        return Err(FiniteDuctileFractureError::Invalid(
            "unrepresentable physical vector norm",
        ));
    }
    Ok(value)
}
fn face_normal(
    nodes: &[[f64; 3]],
    u: &[[f64; 3]],
    face: [usize; 3],
) -> Result<([f64; 3], f64), FiniteDuctileFractureError> {
    let normal = cross(
        difference(nodes, u, face[0], face[1]),
        difference(nodes, u, face[0], face[2]),
    );
    let length = magnitude(normal)?;
    if length <= 0.0 || 0.5 * length <= 0.0 {
        return Err(FiniteDuctileFractureError::Invalid(
            "positive current released-face area required",
        ));
    }
    Ok((normal.map(|v| v / length), 0.5 * length))
}
fn determinant(f: &[f64; 9]) -> f64 {
    f[0].mul_add(
        f[4].mul_add(f[8], -f[5] * f[7]),
        (-f[1]).mul_add(
            f[3].mul_add(f[8], -f[5] * f[6]),
            f[2] * f[3].mul_add(f[7], -f[4] * f[6]),
        ),
    )
}
