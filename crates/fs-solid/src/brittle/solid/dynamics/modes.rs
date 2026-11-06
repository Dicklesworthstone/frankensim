//! Vibration pencils from the original accepted brittle history and inertia.
//!
//! Native frozen constitutive storage, the existing P1 tangent assembly and
//! the EXACT retained consistent mass supply `K` and `M`. Only explicitly fixed
//! incremental DOFs are removed. Fracture does not delete mass, substitute the
//! softening loading tangent, or add artificial stiffness to zero modes.
//!
//! This is an instantaneous, conditional linearization: coefficients, damage
//! and the opening/compression branch stay fixed. It is not a finite-amplitude
//! vibration guarantee or the stability tangent of evolving fracture. A state
//! on the damage envelope can evolve under an arbitrarily small outward load.
//! The caller establishes equilibrium if natural frequencies about equilibrium
//! are required. Geometric initial-stress stiffness and active contact
//! constraints are absent. Released paired faces must be strictly open in the
//! original fixed-reference normal; tied crack-front vertices remain tied.
//! This check does not establish a contact-free vibration amplitude or discover
//! arbitrary fragment collisions. The result enters the existing `fs-modal`
//! pencil consumers without a new eigensolver.

use fs_evidence::{Ambition, ModelCard, ValidityDomain};
use fs_exec::Cx;
use fs_material::brittle::{BrittleCrackBand, BrittleState};
use fs_material::{MaterialAdmissibility, SmallStrainLaw, Tangent6, Voigt};
use fs_sparse::{Coo, Csr};

use super::BrittleTetMotion;
use crate::brittle::solid::fracture::{ReleasedBrittleFace, require_released_faces};
use crate::brittle::solid::{BrittleTetError, checkpoint};
use crate::inertia::P1ReferenceMass;
use crate::linear3::TetAssemblyBudget;
use crate::loads3::TetLoadGeometry3;
use crate::plastic3::{IncrementSettings, SmallStrainLawField, SmallStrainTetProblem};

/// Explicit bounds for a read-only sparse extraction; no dense operator is built.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct FrozenBrittleModesControls {
    /// Maximum full displacement DOFs, before incremental support elimination.
    pub maximum_dofs: usize,
    /// Applied separately to full stiffness staging (144 entries per P1 tet)
    /// and each reduced COO matrix. The retained mass is not rebuilt.
    pub maximum_staged_entries: usize,
}

/// Immutable frozen-history operators in the original reference coordinates.
///
/// Stiffness may be positive SEMIDEFINITE: rigid motions and completely failed
/// transfer modes remain zero. The mass is the exact supported submatrix of
/// the motion's retained reference mass. No modal truncation, extra damping,
/// fracture heat, or physical time advancement is hidden in this observation.
#[derive(Clone, Debug)]
pub struct FrozenBrittlePencil {
    stiffness: Csr,
    mass: Csr,
    free_dofs: Vec<usize>,
    full_dofs: usize,
    reference_mass_kg: f64,
    element_damage: Vec<f64>,
}

impl FrozenBrittlePencil {
    /// Work-conjugate frozen stiffness [N/m], including physical zero modes.
    #[must_use]
    pub const fn stiffness_n_m(&self) -> &Csr {
        &self.stiffness
    }

    /// Original consistent reference mass after support elimination [kg].
    #[must_use]
    pub const fn mass_kg(&self) -> &Csr {
        &self.mass
    }

    /// Ascending full node-major xyz DOF for each reduced row/column.
    #[must_use]
    pub fn free_dofs(&self) -> &[usize] {
        &self.free_dofs
    }

    /// Size of the original coordinate space, before support elimination.
    #[must_use]
    pub const fn full_dofs(&self) -> usize {
        self.full_dofs
    }

    /// Complete retained physical mass, including supported coordinates [kg].
    #[must_use]
    pub const fn reference_mass_kg(&self) -> f64 {
        self.reference_mass_kg
    }

    /// Original accepted damage, in connectivity order. This is an observation;
    /// surviving stiffness uses the native integrity, not `1 - damage`.
    #[must_use]
    pub fn element_damage(&self) -> &[f64] {
        &self.element_damage
    }
}

impl BrittleTetMotion<'_> {
    /// Extract a conditional vibration pencil without changing this motion.
    ///
    /// `fixed_dofs` declare zero incremental displacement in original node-major
    /// xyz coordinates; they need not equal the supports of a previous load
    /// step. At least one coordinate must remain free. Damaged zero opening
    /// refuses because its compression and opening stiffnesses differ. Released
    /// faces must be strictly open beyond their original geometric tolerance;
    /// touching or closed faces need an explicit contact linearization.
    /// All original constitutive, geometry, count and cancellation refusals
    /// propagate before a complete candidate is returned.
    pub fn frozen_vibration_pencil(
        &self,
        fixed_dofs: &[usize],
        controls: FrozenBrittleModesControls,
        cx: &Cx<'_>,
    ) -> Result<FrozenBrittlePencil, BrittleTetError> {
        require_open_released_faces(
            self.solid.geometry.nodes_m().len(),
            &self.solid.released_faces,
            &self.solid.solution.displacement_m,
            cx,
        )?;
        assemble_frozen_pencil(
            self.solid.geometry,
            self.solid.budget,
            self.solid.materials.iter().map(|source| source.law()),
            &self.solid.solution.element_states,
            &self.solid.solution.displacement_m,
            &self.mass,
            fixed_dofs,
            controls,
            cx,
        )
    }
}

/// Use the native paired-face geometry and penetration guard, then require a
/// strict opening margin wherever the topology introduced independent traces.
/// A shared crack-front node has no relative DOF and needs no contact tangent.
pub(crate) fn require_open_released_faces(
    node_count: usize,
    released_faces: &[ReleasedBrittleFace],
    displacement_m: &[[f64; 3]],
    cx: &Cx<'_>,
) -> Result<(), BrittleTetError> {
    checkpoint(cx)?;
    require_released_faces(node_count, released_faces, displacement_m, cx)?;
    for face in released_faces {
        checkpoint(cx)?;
        for vertex in 0..3 {
            let a = face.paired_nodes[0][vertex];
            let b = face.paired_nodes[1][vertex];
            if a == b {
                continue;
            }
            let mut opening = 0.0;
            for component in 0..3 {
                opening = face.normal[component].mul_add(
                    displacement_m[b][component] - displacement_m[a][component],
                    opening,
                );
            }
            if opening <= face.opening_tolerance_m {
                return Err(BrittleTetError::Invalid(
                    "released brittle faces require an explicit contact linearization",
                ));
            }
        }
    }
    Ok(())
}

// Private operator adapter: never exposed as a replacement evolution law.
// Its derivative contract is exact, but only tangent_matrix consumes it.
struct ConstantFrozenTangent(Tangent6);

impl SmallStrainLaw for ConstantFrozenTangent {
    type State = ();
    fn initial_state(&self) {}
    fn stress(&self, strain: &Voigt, (): &()) -> Voigt {
        self.0
            .map(|row| row.iter().zip(strain).map(|(a, b)| a * b).sum())
    }
    fn tangent(&self, _strain: &Voigt, (): &()) -> Tangent6 {
        self.0
    }
    fn update_state(&self, _strain: &Voigt, (): &()) {}
    fn admissibility(&self) -> MaterialAdmissibility {
        MaterialAdmissibility {
            has_stored_energy: true,
            dissipation_nonnegative: true,
            polyconvex: None,
            tangent_symmetric: true,
            failure_envelope: "private frozen-history assembly; no constitutive evolution",
        }
    }
    fn card(&self) -> ModelCard {
        ModelCard::new(
            "solid.frozen-brittle-storage-assembly",
            "0.1.0",
            Ambition::Solid,
            vec!["original damage, coefficients and opening/compression branch held fixed".into()],
            ValidityDomain::unconstrained(),
            vec!["no evolving-fracture stability or finite-amplitude vibration claim".into()],
            f64::INFINITY,
        )
    }
}

/// Shared by mechanical and thermal motion owners. The latter supplies the
/// ORIGINAL mechanical-strain histories and its SAME retained reference mass;
/// thermal eigenstrain is not reinterpreted as a new damage-driving strain.
#[allow(clippy::too_many_arguments)]
pub(crate) fn assemble_frozen_pencil<'m>(
    geometry: &TetLoadGeometry3<'_>,
    budget: TetAssemblyBudget,
    materials: impl ExactSizeIterator<Item = &'m BrittleCrackBand>,
    states: &[BrittleState<'_>],
    displacement_m: &[[f64; 3]],
    mass: &P1ReferenceMass,
    fixed_dofs: &[usize],
    controls: FrozenBrittleModesControls,
    cx: &Cx<'_>,
) -> Result<FrozenBrittlePencil, BrittleTetError> {
    checkpoint(cx)?;
    let full_dofs = geometry
        .nodes_m()
        .len()
        .checked_mul(3)
        .ok_or(BrittleTetError::Invalid(
            "frozen vibration DOF count overflow",
        ))?;
    let elements = geometry.tetrahedra().len();
    let stiffness_entries = elements.checked_mul(144).ok_or(BrittleTetError::Invalid(
        "frozen vibration stiffness count overflow",
    ))?;
    if full_dofs == 0
        || full_dofs > controls.maximum_dofs
        || stiffness_entries > controls.maximum_staged_entries
        || mass.matrix().nnz() > controls.maximum_staged_entries
    {
        return Err(BrittleTetError::Invalid(
            "frozen vibration extraction exceeds its count budget",
        ));
    }
    if states.len() != elements
        || materials.len() != elements
        || displacement_m.len() != geometry.nodes_m().len()
        || mass.ndofs() != full_dofs
        || fixed_dofs.len() > full_dofs
    {
        return Err(BrittleTetError::Invalid(
            "inconsistent frozen vibration source arrays",
        ));
    }
    let mut fixed = vec![false; full_dofs];
    let mut prescribed = Vec::with_capacity(fixed_dofs.len());
    for (index, &dof) in fixed_dofs.iter().enumerate() {
        if index % 256 == 0 {
            checkpoint(cx)?;
        }
        if dof >= full_dofs || fixed[dof] {
            return Err(BrittleTetError::Invalid(
                "incremental supports are out of range or repeated",
            ));
        }
        fixed[dof] = true;
        prescribed.push((dof, 0.0));
    }
    let mut free_dofs = Vec::with_capacity(full_dofs - fixed_dofs.len());
    let mut reduced_of = vec![usize::MAX; full_dofs];
    for (dof, is_fixed) in fixed.into_iter().enumerate() {
        if dof % 256 == 0 {
            checkpoint(cx)?;
        }
        if !is_fixed {
            reduced_of[dof] = free_dofs.len();
            free_dofs.push(dof);
        }
    }
    if free_dofs.is_empty() {
        return Err(BrittleTetError::Invalid(
            "frozen vibration needs at least one free coordinate",
        ));
    }
    let mut laws = Vec::with_capacity(elements);
    let mut element_damage = Vec::with_capacity(elements);
    for (element, (material, state)) in materials.zip(states).enumerate() {
        checkpoint(cx)?;
        let tangent = material
            .frozen_tangent(state)
            .map_err(|source| BrittleTetError::Material { element, source })?;
        laws.push(ConstantFrozenTangent(tangent));
        element_damage.push(state.damage());
    }
    let problem = SmallStrainTetProblem {
        nodes_m: geometry.nodes_m(),
        tetrahedra: geometry.tetrahedra(),
        laws: SmallStrainLawField::PerElement(&laws),
        prescribed_m: &prescribed,
        nodal_forces_n: &[],
        load_path: &[1.0],
        budget,
        settings: IncrementSettings::default(),
    };
    let full_stiffness =
        problem.tangent_matrix(displacement_m.as_flattened(), &vec![(); elements], cx)?;
    let stiffness = reduce(&full_stiffness, &free_dofs, &reduced_of, cx)?;
    let reduced_mass = reduce(mass.matrix(), &free_dofs, &reduced_of, cx)?;
    checkpoint(cx)?;
    Ok(FrozenBrittlePencil {
        stiffness,
        mass: reduced_mass,
        free_dofs,
        full_dofs,
        reference_mass_kg: mass.reference_mass_kg(),
        element_damage,
    })
}

fn reduce(
    full: &Csr,
    free: &[usize],
    reduced_of: &[usize],
    cx: &Cx<'_>,
) -> Result<Csr, BrittleTetError> {
    let mut entries = Coo::new(free.len(), free.len());
    for (row, &dof) in free.iter().enumerate() {
        checkpoint(cx)?;
        let (columns, values) = full.row(dof);
        for (index, (&column, &value)) in columns.iter().zip(values).enumerate() {
            if index % 256 == 0 {
                checkpoint(cx)?;
            }
            let reduced_column = reduced_of[column];
            if reduced_column != usize::MAX {
                entries.push(row, reduced_column, value);
            }
        }
    }
    checkpoint(cx)?;
    Ok(entries.assemble())
}
