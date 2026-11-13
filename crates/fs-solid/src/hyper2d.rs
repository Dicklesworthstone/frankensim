//! Finite-strain plane hyperelasticity: kinematics, weak forms,
//! and the globalized Newton loop. Constitutive behavior comes from
//! fs-material cards (exact AD Piola stress and consistent 9×9
//! tangent); the default embeds the 2D displacement gradient into the
//! 3D deformation gradient (F₃₃ = 1), integrates the residual
//! `∫ P : ∇δu`, and assembles the consistent tangent
//! `∫ ∇δu : A : ∇Δu`.
//!
//! Globalization: backtracking (Armijo) line search on the total
//! potential energy with dead loads, with material refusals
//! (det F ≤ 0 on a trial state) absorbed by step halving, plus load
//! stepping — the battery drives a buckling-adjacent fixture through
//! exactly this recipe.
//! Prescribed motion enters the same Newton direction through `K_fp du_p`;
//! free and prescribed updates backtrack together. While boundary motion is
//! pending, an admissible step need not decrease energy: the moving support
//! does work. Force convergence is accepted only at the exact target supports.
//! Equilibrium and the Newton equation use relative residuals, so changing
//! the common scale of stiffness and dead traction cannot erase deformation.
//! Reference and current P1/Q1 maps must preserve orientation at every corner;
//! positive constitutive Jacobians at quadrature sites alone cannot admit a
//! folded quadrilateral. The same geometry refusal governs energy line searches,
//! force probes, tangents and fully prescribed equilibria.
//!
//! [`mean_dilatation`] adds an explicit Q1 condensed-volumetric formulation
//! for quasi-incompressible finite compression, sharing this Newton driver.
//! Its energy and stress recovery differ from the default pointwise law.
//! [`plane_stress`] instead relaxes all three transverse tractions for thin
//! continua and uses the corresponding reduced energy and consistent tangent.

use crate::SolidError;
use crate::mesh2::{
    Mesh2, Patch, checked_reference_element, positive_deformed_orientation, quad_points, shapes_at,
};
use fs_exec::Cx;
use fs_material::hyper::Hyperelastic;
use fs_material::hyper::cards::{
    FiniteStrainPotential, HyperelasticMaterial, ResolvedHyperelastic,
};
use fs_material::hyper::incompressible_plane_stress::IncompressiblePlaneStressPotential;
use fs_material::hyper::reference::ReferencedHyperelastic;
use fs_material::visco_hyper::FrozenViscoHyperelastic;
use fs_solver::krylov::MinresState;
use fs_solver::op::{CsrOp, LinearOp};
use fs_sparse::Coo;
use std::collections::BTreeMap;

/// Native reduced-plane equilibrium paths, including prescribed support ramps.
pub mod continuation;
/// Opt-in Q1 mean-dilatation finite-strain compression.
pub mod mean_dilatation;
/// Traction-free thickness relaxation on the original hyperelastic energy.
pub mod plane_stress;
/// Retained finite-plastic P1/Q1 sheet forming with solved transverse thickness.
pub mod plastic;
/// Loaded equilibria and their original consistent-mass tangent spectra.
pub mod vibration;
/// Persistent finite-strain viscoelastic P1/Q1 histories and physical intervals.
pub mod visco;

/// One reference-area-weighted material point of the native P1/Q1 continuum.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlaneMaterialPoint {
    /// Element index in the original mesh connectivity.
    pub element: usize,
    /// Stable local quadrature index: centroid zero for P1, original Gauss
    /// order for Q1. History owners must retain this ordering.
    pub local_point: usize,
    /// Positive quadrature area in square metres; the physical mesh coordinates
    /// must be in metres. Multiply by original thickness to integrate an
    /// original-volume material density.
    pub reference_area_m2: f64,
    /// Compatible total in-plane deformation in `[F00,F01,F10,F11]` order.
    pub in_plane: [f64; 4],
}

/// Recover native material-point kinematics without evaluating a constitutive
/// law. Geometry, orientation admission, shape gradients and material-point
/// ordering are shared with the plane-stress assembly. The borrowed mesh and
/// displacement remain unchanged on refusal or cancellation.
pub fn material_point_kinematics(
    mesh: &Mesh2,
    displacement: &[[f64; 2]],
    cx: &Cx<'_>,
) -> Result<Vec<PlaneMaterialPoint>, SolidError> {
    plane_stress::material_point_kinematics(mesh, displacement, cx)
}

fn material_quadrature(arity: usize) -> Vec<(f64, f64, f64)> {
    if arity == 3 {
        // A P1 displacement has constant F. Its single retained material point
        // is shared by plane strain and both plane-stress reductions.
        vec![(1.0 / 3.0, 1.0 / 3.0, 0.5)]
    } else {
        quad_points(4)
    }
}

#[derive(Clone, Copy)]
enum Reduction<'a, 'cx> {
    MeanDilatation(&'a mean_dilatation::Prepared, &'a Cx<'cx>),
    PlaneStress(&'a plane_stress::Prepared<'a>, &'a Cx<'cx>),
}

type Projection<'a, 'cx> = Option<Reduction<'a, 'cx>>;

fn checkpoint(projection: Projection<'_, '_>) -> Result<(), SolidError> {
    if let Some(Reduction::MeanDilatation(_, cx) | Reduction::PlaneStress(_, cx)) = projection {
        cx.checkpoint().map_err(|_| SolidError::Cancelled)?;
    }
    Ok(())
}

/// Newton controls.
#[derive(Debug, Clone, Copy)]
pub struct NewtonSettings {
    /// Relative free-force residual tolerance, strictly between zero and one.
    /// Each static load step uses the largest of its initial free residual,
    /// applied free load and linearized prescribed-motion force norms and the
    /// force its starting deformation carries (stiffness times motion relative
    /// to each unsupported row's node, as hyper3's static gate).
    /// There is no implicit force-unit floor.
    pub tol: f64,
    /// Iteration cap per load step.
    pub max_iters: usize,
    /// Load steps (the Dirichlet data and tractions ramp linearly).
    pub load_steps: usize,
    /// Backtracking halvings before declaring a stall.
    pub max_backtracks: usize,
}

impl Default for NewtonSettings {
    fn default() -> Self {
        NewtonSettings {
            tol: 1e-10,
            max_iters: 30,
            load_steps: 1,
            max_backtracks: 25,
        }
    }
}

/// Convergence evidence for one solve (all load steps).
#[derive(Debug, Clone, Default)]
pub struct NewtonReport {
    /// Relative free-force residuals per iteration, per load step.
    pub histories: Vec<Vec<f64>>,
    /// Fixed force normalization for each corresponding history, per unit
    /// reference thickness. Zero denotes an exactly unloaded equilibrium.
    pub force_scales: Vec<f64>,
    /// Total line-search backtracks taken.
    pub backtracks: usize,
}

/// A displacement/traction field on a patch.
pub type VecField<'m> = &'m dyn Fn(f64, f64) -> [f64; 2];

/// An immutable plane-continuum potential and its exact volume constraint.
///
/// Reference-deformation and frozen-history adapters retain their own energy,
/// stress, tangent and relaxed-volume convention through every spatial trial.
/// A frozen inelastic potential is an incremental potential; equilibrium does
/// not commit its history or turn its value into physical stored energy/heat.
pub trait PlaneMaterial: FiniteStrainPotential + IncompressiblePlaneStressPotential {
    /// Admit public parameters once before spatial assembly. Actual deformation
    /// and history admission remains with the constitutive evaluation at each
    /// material point; identity need not lie in a prestrained law's domain.
    fn validate_plane_material(&self) -> Result<(), fs_material::MaterialError> {
        Ok(())
    }

    /// Return a bare native law only when it is the complete potential owner.
    ///
    /// This admits the existing scalar-generic mean-dilatation formulation and
    /// preserves optimized native transverse reductions. Reference/history
    /// adapters must keep the default `None`: extracting their elastic backbone
    /// would discard the prescribed strain or material memory.
    fn native_hyperelastic(&self) -> Option<&Hyperelastic> {
        None
    }

    /// Original stored-energy owner, when this potential is elastic.
    ///
    /// Vibration uses this only to admit physical elastic storage; assembly
    /// still calls this complete plane material, including its reference map.
    /// Frozen inelastic incremental potentials must retain the default `None`:
    /// their algorithmic Hessian is not an elastic vibration stiffness.
    fn elastic_potential(&self) -> Option<&dyn HyperelasticMaterial> {
        None
    }
}

fn validate_native(material: &Hyperelastic) -> Result<(), fs_material::MaterialError> {
    Hyperelastic::new(material.model.clone(), material.stretch_limit).map(|_| ())
}

impl PlaneMaterial for Hyperelastic {
    fn validate_plane_material(&self) -> Result<(), fs_material::MaterialError> {
        validate_native(self)
    }

    fn native_hyperelastic(&self) -> Option<&Hyperelastic> {
        Some(self)
    }

    fn elastic_potential(&self) -> Option<&dyn HyperelasticMaterial> {
        Some(self)
    }
}

impl PlaneMaterial for ResolvedHyperelastic {
    fn native_hyperelastic(&self) -> Option<&Hyperelastic> {
        Some(self.hyperelastic())
    }

    fn elastic_potential(&self) -> Option<&dyn HyperelasticMaterial> {
        Some(self)
    }
}

impl<L: HyperelasticMaterial + IncompressiblePlaneStressPotential + ?Sized> PlaneMaterial
    for ReferencedHyperelastic<'_, L>
{
    fn validate_plane_material(&self) -> Result<(), fs_material::MaterialError> {
        validate_native(self.material().hyperelastic())
    }

    fn elastic_potential(&self) -> Option<&dyn HyperelasticMaterial> {
        Some(self)
    }
}

impl PlaneMaterial for FrozenViscoHyperelastic<'_> {}

/// Original finite-strain potential owners in mesh-element order.
///
/// A single law is uniform. Arrays, vectors and borrowed slices assign exactly
/// one original law per element; no coefficients or fibre directions are
/// averaged across an interface. The selected law is borrowed before the
/// element's quadrature/AD work, without cloning materials in the assembly loop.
///
/// Implementations must keep the mapping fixed while a problem is borrowed.
/// `None` from `material_count` declares the same uniform law at every index.
/// A point field declares an exact element count and exact point counts, and
/// returns each immutable owner from `material_at_point`. Its element-only
/// access may return `None`; consumers never replace it with the first point.
pub trait HyperelasticMaterials {
    /// Exact number of element laws, or `None` for a uniform law.
    fn material_count(&self) -> Option<usize>;
    /// Original law at this element, or `None` for missing coverage.
    fn material_at(&self, element: usize) -> Option<&dyn PlaneMaterial>;

    /// Exact number of independent material points in this element. `None`
    /// means the element's one law is shared by all its quadrature points.
    /// Point fields require one centroid for P1 or four Gauss points for Q1.
    fn material_point_count(&self, _element: usize) -> Option<usize> {
        None
    }

    /// Original owner at this element and local quadrature point. A point
    /// field overrides this and `material_point_count` together; neither law
    /// coefficients nor accepted internal histories are averaged.
    fn material_at_point(&self, element: usize, _local_point: usize) -> Option<&dyn PlaneMaterial> {
        self.material_at(element)
    }

    /// Optional accepted transverse `[F02,F12,F22]` used only to initialize
    /// local plane-stress relaxation. Exact elastic volume still determines
    /// its own new thickness; the seed never imposes a transverse constraint.
    fn transverse_seed_at(&self, _element: usize, _local_point: usize) -> Option<[f64; 3]> {
        None
    }
}

impl<L: PlaneMaterial> HyperelasticMaterials for L {
    fn material_count(&self) -> Option<usize> {
        None
    }

    fn material_at(&self, _: usize) -> Option<&dyn PlaneMaterial> {
        Some(self)
    }
}

impl<L: PlaneMaterial> HyperelasticMaterials for Vec<L> {
    fn material_count(&self) -> Option<usize> {
        Some(self.len())
    }

    fn material_at(&self, element: usize) -> Option<&dyn PlaneMaterial> {
        self.get(element).map(|law| law as &dyn PlaneMaterial)
    }
}

impl<L: PlaneMaterial, const N: usize> HyperelasticMaterials for [L; N] {
    fn material_count(&self) -> Option<usize> {
        Some(N)
    }

    fn material_at(&self, element: usize) -> Option<&dyn PlaneMaterial> {
        self.get(element).map(|law| law as &dyn PlaneMaterial)
    }
}

impl<L: PlaneMaterial> HyperelasticMaterials for &[L] {
    fn material_count(&self) -> Option<usize> {
        Some(self.len())
    }

    fn material_at(&self, element: usize) -> Option<&dyn PlaneMaterial> {
        self.get(element).map(|law| law as &dyn PlaneMaterial)
    }
}

/// A finite-strain problem on a body-fitted mesh (plane strain by default).
/// [`Self::solve_plane_stress`] admits traction-free thickness relaxation.
pub struct HyperProblem<'m> {
    /// The mesh.
    pub mesh: &'m Mesh2,
    /// A uniform original law, or exactly one law per mesh element.
    /// Native/source-bound, reference-deformation and frozen-history owners are
    /// retained. Pass a vector, array or borrowed slice reference to assign one
    /// owner per element. A frozen-history solve returns equilibrium only; its
    /// material owner must subsequently accept the converged constitutive step.
    pub material: &'m dyn HyperelasticMaterials,
    /// Dirichlet patches (ramped by the load factor).
    pub dirichlet: Vec<(Patch, VecField<'m>)>,
    /// Dead-load traction patches (ramped by the load factor).
    pub traction: Vec<(Patch, VecField<'m>)>,
    /// Newton controls.
    pub settings: NewtonSettings,
}

impl HyperProblem<'_> {
    fn element_material(&self, element: usize) -> Result<&dyn PlaneMaterial, SolidError> {
        self.material
            .material_at(element)
            .ok_or_else(|| SolidError::InvalidInput {
                what: format!("hyperelastic material missing at element {element}"),
            })
    }

    fn point_material(
        &self,
        element: usize,
        local_point: usize,
    ) -> Result<&dyn PlaneMaterial, SolidError> {
        self.material
            .material_at_point(element, local_point)
            .ok_or_else(|| SolidError::InvalidInput {
                what: format!(
                    "finite-strain material missing at element {element}, point {local_point}"
                ),
            })
    }

    fn mean_dilatation_material(&self, element: usize) -> Result<&Hyperelastic, SolidError> {
        if self.material.material_point_count(element).is_some() {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "mean dilatation has no admitted condensation for independent point materials at element {element}"
                ),
            });
        }
        self.element_material(element)?
            .native_hyperelastic()
            .ok_or_else(|| SolidError::InvalidInput {
                what: format!(
                    "mean dilatation requires a bare native potential at element {element}; \
                     reference-deformation and frozen-history adapters need their own condensation"
                ),
            })
    }

    fn validate_materials(&self) -> Result<(), SolidError> {
        self.validate_materials_with(|| Ok(()))
    }

    fn validate_materials_with(
        &self,
        mut checkpoint: impl FnMut() -> Result<(), SolidError>,
    ) -> Result<(), SolidError> {
        let count = self.material.material_count();
        if let Some(count) = count
            && count != self.mesh.elems.len()
        {
            return Err(SolidError::InvalidInput {
                what: format!(
                    "hyperelastic field has {count} laws for {} elements",
                    self.mesh.elems.len()
                ),
            });
        }
        // Exact point coverage is admitted before any assembly. Uniform and
        // per-element owners still validate parameters only once per owner;
        // independent point histories each retain their own admission.
        for (element, conn) in self.mesh.elems.iter().enumerate() {
            checkpoint()?;
            let expected = match conn.len() {
                3 => 1,
                4 => 4,
                _ => {
                    return Err(SolidError::InvalidInput {
                        what: format!("finite-strain element {element} must be P1 or Q1"),
                    });
                }
            };
            if checked_reference_element(&self.mesh.nodes, conn).is_none() {
                return Err(SolidError::InvalidInput {
                    what: format!(
                        "finite-strain reference element {element} has invalid connectivity or a folded, reversed or nonfinite map"
                    ),
                });
            }
            let points = self.material.material_point_count(element);
            if let Some(points) = points
                && (count.is_none() || points != expected)
            {
                return Err(SolidError::InvalidInput {
                    what: format!(
                        "finite-strain element {element} needs {expected} material points and an explicit element count; got {points} points"
                    ),
                });
            }
            if count.is_none() && element != 0 {
                continue;
            }
            for local_point in 0..points.unwrap_or(1) {
                checkpoint()?;
                self.point_material(element, local_point)?
                    .validate_plane_material()
                    .map_err(|error| SolidError::MaterialRefused {
                        what: format!(
                            "finite-strain element {element}, point {local_point}: {error}"
                        ),
                    })?;
            }
        }
        Ok(())
    }

    fn validate_displacement(&self, displacement: &[f64]) -> Result<(), SolidError> {
        if self.mesh.nodes.len().checked_mul(2) != Some(displacement.len())
            || displacement.iter().any(|value| !value.is_finite())
        {
            return Err(SolidError::InvalidInput {
                what: "finite-strain displacement must contain one finite xy pair per mesh node"
                    .into(),
            });
        }
        Ok(())
    }

    fn admit_current_element(
        &self,
        element: usize,
        conn: &[usize],
        displacement: &[f64],
    ) -> Result<(), SolidError> {
        if !positive_deformed_orientation(&self.mesh.nodes, conn, displacement) {
            return Err(SolidError::MaterialRefused {
                what: format!(
                    "finite-strain current element {element} is folded, reversed or unrepresentable"
                ),
            });
        }
        Ok(())
    }

    /// Solve by load stepping + globalized Newton; returns nodal
    /// displacements and the convergence evidence.
    ///
    /// # Errors
    /// [`SolidError::NewtonStalled`], [`SolidError::SolveFailed`],
    /// [`SolidError::MaterialRefused`] (only if globalization cannot
    /// absorb a refusal), [`SolidError::UnknownPatch`].
    pub fn solve(&self) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        self.solve_with_projection(None)
    }

    /// Solve from an admitted nodal displacement in the original mesh frame.
    ///
    /// This is a numerical starting state, not a new reference configuration
    /// or a constitutive history commit. The original absolute boundary data
    /// and load-step ramp are retained; use one load step for a changed-reference
    /// endpoint with already deformed initial geometry. Each prescribed value
    /// is lifted exactly at its load step. The borrowed input is never changed.
    ///
    /// # Errors
    /// Refuses incomplete/nonfinite initial data or an inadmissible initial
    /// deformation before boundary lifting, as well as the failures of [`Self::solve`].
    pub fn solve_from(
        &self,
        initial_displacement: &[[f64; 2]],
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        self.solve_with_projection_from(None, Some(initial_displacement))
    }

    fn solve_with_projection(
        &self,
        projection: Projection<'_, '_>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        self.solve_with_projection_from(projection, None)
    }

    fn solve_with_projection_from(
        &self,
        projection: Projection<'_, '_>,
        initial_displacement: Option<&[[f64; 2]]>,
    ) -> Result<(Vec<[f64; 2]>, NewtonReport), SolidError> {
        checkpoint(projection)?;
        // Reduced formulations already admitted every law during preparation.
        if projection.is_none() {
            self.validate_materials()?;
        }
        if !(self.settings.tol.is_finite() && self.settings.tol > 0.0 && self.settings.tol < 1.0)
            || self.settings.max_iters == 0
            || self.settings.load_steps == 0
            || self.settings.max_backtracks == 0
        {
            return Err(SolidError::InvalidInput {
                what: "Newton needs a relative tolerance in (0, 1) and positive iteration/load/backtrack limits".into(),
            });
        }
        let n = self.mesh.node_count();
        let mut u = if let Some(initial) = initial_displacement {
            if initial.len() != n {
                return Err(SolidError::InvalidInput {
                    what: "initial finite-strain displacement must cover every mesh node".into(),
                });
            }
            let mut u = Vec::with_capacity(2 * n);
            for (node, value) in initial.iter().enumerate() {
                if node % 256 == 0 {
                    checkpoint(projection)?;
                }
                if value.iter().any(|v| !v.is_finite()) {
                    return Err(SolidError::InvalidInput {
                        what: "initial finite-strain displacement must be finite".into(),
                    });
                }
                u.extend_from_slice(value);
            }
            // Use the same actual-state reduction as Newton. No identity
            // substitute, reference reset or target-boundary overwrite may
            // authorize an inadmissible supplied starting state.
            if !self
                .energy_with_projection(&u, 0.0, projection)?
                .is_finite()
            {
                return Err(SolidError::MaterialRefused {
                    what: "initial finite-strain deformation is outside the material domain".into(),
                });
            }
            u
        } else {
            vec![0.0f64; 2 * n]
        };
        let mut report = NewtonReport::default();
        for step in 1..=self.settings.load_steps {
            checkpoint(projection)?;
            #[allow(clippy::cast_precision_loss)]
            let load = step as f64 / self.settings.load_steps as f64;
            let fixed = self.fixed_dofs(load)?;
            // With no free motion, the target itself is the complete state;
            // admit it directly as before. Otherwise a boundary-only jump can
            // invert the first element layer even for a valid affine squeeze.
            if fixed.len() == u.len() {
                for (&dof, &val) in &fixed {
                    u[dof] = val;
                }
            }
            let mut applied = vec![0.0; u.len()];
            self.for_traction(load, |dof, value| {
                if !fixed.contains_key(&dof) {
                    applied[dof] += value;
                }
            })?;
            if let Some(Reduction::PlaneStress(prepared, cx)) = projection {
                prepared.complete_applied_force(&mut applied, &fixed, cx)?;
            }
            let free: Vec<usize> = match projection {
                Some(Reduction::PlaneStress(prepared, _)) if prepared.is_dynamic() => (0..u.len())
                    .filter(|dof| !fixed.contains_key(dof))
                    .collect(),
                _ => Vec::new(),
            };
            let applied_norm = norm(&applied);
            if !applied_norm.is_finite() {
                return Err(SolidError::InvalidInput {
                    what: "Newton free-load norm is not finite".into(),
                });
            }
            let mut force_scale = applied_norm;
            let mut history = Vec::new();
            let mut converged = false;
            for _ in 0..self.settings.max_iters {
                checkpoint(projection)?;
                let pending_any = fixed.iter().any(|(&dof, &target)| u[dof] != target);
                let pending: BTreeMap<_, _> = if pending_any {
                    fixed
                        .iter()
                        .map(|(&dof, &target)| (dof, target - u[dof]))
                        .collect()
                } else {
                    BTreeMap::new()
                };
                if pending.values().any(|value| !value.is_finite()) {
                    return Err(SolidError::InvalidInput {
                        what: "prescribed finite-strain motion must be finite and representable"
                            .into(),
                    });
                }
                let no_constraints = BTreeMap::new();
                let (mut resid, tangent) = self.assemble_with_projection(
                    &u,
                    load,
                    if pending_any { &no_constraints } else { &fixed },
                    projection,
                )?;
                for &dof in fixed.keys() {
                    resid[dof] = 0.0;
                }
                // The current free residual is the convergence/Armijo quantity;
                // its boundary linearization belongs only to the Newton RHS.
                let mut lifted_residual = pending_any.then(|| resid.clone());
                let tangent = if let Some(rhs) = &mut lifted_residual {
                    eliminate_prescribed_increment(tangent, rhs, &pending)
                } else {
                    tangent
                };
                let newton_residual = lifted_residual.as_deref().unwrap_or(&resid);
                let rnorm = norm(&resid);
                let drive_norm = if pending_any {
                    norm(newton_residual)
                } else {
                    rnorm
                };
                if !rnorm.is_finite() || !drive_norm.is_finite() {
                    return Err(SolidError::InvalidInput {
                        what: "Newton free-force residual or prescribed-motion drive is not finite"
                            .into(),
                    });
                }
                let momentum_scale = match projection {
                    Some(Reduction::PlaneStress(prepared, cx)) => {
                        prepared.dynamic_force_scale(&u, &resid, &applied, &free, cx)?
                    }
                    _ => None,
                };
                if let Some(scale) = momentum_scale {
                    force_scale = scale;
                } else if history.is_empty() {
                    force_scale = force_scale
                        .max(rnorm)
                        .max(drive_norm)
                        .max(self.deformation_force(&u, load, &fixed, projection)?);
                }
                let relative = if rnorm == 0.0 {
                    0.0
                } else {
                    rnorm / force_scale
                };
                history.push(relative);
                if !pending_any && relative < self.settings.tol {
                    converged = true;
                    break;
                }
                // Solve K δ = −R (MINRES: symmetric, indefinite-safe
                // near buckling). Normalize BOTH the operator and RHS before
                // Krylov dot products. Identity rows for supports are not a
                // physical stiffness scale; their iterates remain exactly zero.
                let mut stiffness_scale = 0.0_f64;
                for row in 0..resid.len() {
                    if !fixed.contains_key(&row) {
                        for &value in tangent.row(row).1 {
                            if !value.is_finite() {
                                return Err(SolidError::SolveFailed {
                                    iters: 0,
                                    rel_residual: f64::INFINITY,
                                });
                            }
                            stiffness_scale = stiffness_scale.max(value.abs());
                        }
                    }
                }
                if stiffness_scale == 0.0 && drive_norm != 0.0 {
                    return Err(SolidError::SolveFailed {
                        iters: 0,
                        rel_residual: 1.0,
                    });
                }
                let mut direction = vec![0.0; u.len()];
                if drive_norm != 0.0 {
                    let b: Vec<f64> = newton_residual.iter().map(|r| -r / drive_norm).collect();
                    let op = normalized_newton_operator(tangent, stiffness_scale, &fixed);
                    let mut st = MinresState::new(&op, &b);
                    while st.iters < 20_000 && st.rel_residual() >= 1e-10 {
                        checkpoint(projection)?;
                        let before = st.iters;
                        let allowance = (20_000 - before).min(32);
                        let _ = st.run(&op, 1e-10, allowance);
                        if st.iters == before {
                            break;
                        }
                    }
                    // Recurrence convergence alone cannot admit the actual
                    // lifted Newton equation, including its support coupling.
                    let mut defect = vec![0.0; b.len()];
                    op.apply(&st.x, &mut defect);
                    for (value, rhs) in defect.iter_mut().zip(&b) {
                        *value -= rhs;
                    }
                    let linear_relative = norm(&defect) / norm(&b);
                    if !linear_relative.is_finite() || linear_relative > 1e-6 {
                        return Err(SolidError::SolveFailed {
                            iters: st.iters,
                            rel_residual: linear_relative,
                        });
                    }
                    let displacement_scale = drive_norm / stiffness_scale;
                    for (value, &x) in direction.iter_mut().zip(&st.x) {
                        *value = x * displacement_scale;
                    }
                }
                for (&dof, &value) in &pending {
                    direction[dof] = value;
                }
                // Armijo backtracking on the potential energy; material
                // refusals (det F ≤ 0) halve the step.
                let slope: f64 = resid.iter().zip(&direction).map(|(r, d)| r * d).sum();
                let e0 = self.energy_with_projection(&u, load, projection)?;
                if !e0.is_finite() || !slope.is_finite() {
                    return Err(SolidError::NewtonStalled { history });
                }
                let mut alpha = 1.0f64;
                let mut taken = false;
                for _ in 0..self.settings.max_backtracks {
                    checkpoint(projection)?;
                    let mut trial: Vec<f64> = u
                        .iter()
                        .zip(&direction)
                        .map(|(ui, di)| ui + alpha * di)
                        .collect();
                    // A full accepted step installs the exact target bits;
                    // rounded subtraction/addition cannot leave a tiny pending
                    // support increment and keep the solver in the work phase.
                    if alpha == 1.0 {
                        for (&dof, &target) in &fixed {
                            trial[dof] = target;
                        }
                    }
                    if let Some(e1) = self.try_energy_with_projection(&trial, load, projection)?
                        && e1.is_finite()
                    {
                        let mut acceptable =
                            pending_any || (slope < 0.0 && e1 <= e0 + 1e-4 * alpha * slope);
                        // At an indefinite tangent use residual globalization;
                        // near energy roundoff also demand actual force progress.
                        // Neither gate introduces a dimensional unit floor.
                        let roundoff =
                            64.0 * f64::EPSILON * e0.abs().max(e1.abs()).max((alpha * slope).abs());
                        let force_merit = slope >= 0.0 || e1 <= e0 + roundoff;
                        // Once the predicted decrease itself is below energy
                        // resolution, cancellation in a constitutive potential
                        // can exceed this net-energy estimate. Admit that case
                        // only if the actual candidate force already satisfies
                        // the requested nonlinear tolerance. A small predicted
                        // decrease alone must never authorize an unconverged
                        // uphill step or loosen the force acceptance criterion.
                        let unresolved_decrease = (alpha * slope).abs() <= roundoff;
                        if !acceptable && (force_merit || unresolved_decrease) {
                            match self.residual_with_projection(&trial, load, &fixed, projection) {
                                Ok(trial_residual) => {
                                    let trial_norm = norm(&trial_residual);
                                    let trial_scale = match projection {
                                        Some(Reduction::PlaneStress(prepared, cx)) => prepared
                                            .dynamic_force_scale(
                                                &trial,
                                                &trial_residual,
                                                &applied,
                                                &free,
                                                cx,
                                            )?
                                            .unwrap_or(force_scale),
                                        _ => force_scale,
                                    };
                                    acceptable = trial_norm.is_finite()
                                        && trial_norm < (1.0 - 1e-4 * alpha) * rnorm
                                        && (force_merit
                                            || trial_norm / trial_scale < self.settings.tol);
                                }
                                Err(SolidError::MaterialRefused { .. }) => {}
                                Err(error) => return Err(error),
                            }
                        }
                        if acceptable {
                            u = trial;
                            taken = true;
                            break;
                        }
                    }
                    alpha *= 0.5;
                    report.backtracks += 1;
                }
                if !taken {
                    return Err(SolidError::NewtonStalled { history });
                }
            }
            if !converged {
                return Err(SolidError::NewtonStalled { history });
            }
            report.histories.push(history);
            report.force_scales.push(force_scale);
        }
        checkpoint(projection)?;
        Ok(((0..n).map(|i| [u[2 * i], u[2 * i + 1]]).collect(), report))
    }

    fn residual_with_projection(
        &self,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
        projection: Projection<'_, '_>,
    ) -> Result<Vec<f64>, SolidError> {
        let mut residual = match projection {
            Some(Reduction::PlaneStress(prepared, cx)) => {
                prepared.force_residual(self, u, load, cx)?
            }
            Some(Reduction::MeanDilatation(_, _)) => {
                self.assemble_with_projection(u, load, fixed, projection)?.0
            }
            None => self.force_residual(u, load)?,
        };
        for &dof in fixed.keys() {
            residual[dof] = 0.0;
        }
        Ok(residual)
    }

    fn assemble_with_projection(
        &self,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
        projection: Projection<'_, '_>,
    ) -> Result<(Vec<f64>, fs_sparse::Csr), SolidError> {
        match projection {
            Some(Reduction::MeanDilatation(prepared, cx)) => {
                prepared.assemble(self, u, load, fixed, cx)
            }
            Some(Reduction::PlaneStress(prepared, cx)) => {
                prepared.assemble(self, u, load, fixed, cx)
            }
            None => self.assemble(u, load, fixed),
        }
    }

    /// The force a deformed state already carries: on each unsupported row,
    /// stiffness times motion relative to the row's node, from the
    /// unconstrained tangent (the deformation term of hyper3's
    /// `static_force_scale`). A state at rest carries none and is not assembled
    /// again. Without it an equilibrium start under no load and no pending
    /// motion had only its own rounding-level residual as its scale.
    fn deformation_force(
        &self,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
        projection: Projection<'_, '_>,
    ) -> Result<f64, SolidError> {
        if u.iter().all(|&value| value == 0.0) {
            return Ok(0.0);
        }
        let (_, tangent) = self.assemble_with_projection(u, load, &BTreeMap::new(), projection)?;
        let mut total = 0.0_f64;
        for row in 0..u.len() {
            if row % 256 == 0 {
                checkpoint(projection)?;
            }
            if fixed.contains_key(&row) {
                continue;
            }
            let origin = 2 * (row / 2);
            let (columns, coefficients) = tangent.row(row);
            let mut row_force = 0.0_f64;
            for (&column, &coefficient) in columns.iter().zip(coefficients) {
                let relative = u[column] - u[origin + column % 2];
                row_force = coefficient.abs().mul_add(relative.abs(), row_force);
            }
            total = total.hypot(row_force);
        }
        if !total.is_finite() {
            return Err(SolidError::InvalidInput {
                what: "connected deformation force scale is not representable".into(),
            });
        }
        Ok(total)
    }

    fn try_energy_with_projection(
        &self,
        u: &[f64],
        load: f64,
        projection: Projection<'_, '_>,
    ) -> Result<Option<f64>, SolidError> {
        match projection {
            Some(reduction) => {
                let energy = match reduction {
                    Reduction::MeanDilatation(prepared, cx) => {
                        prepared.potential_energy(self, u, load, cx)
                    }
                    Reduction::PlaneStress(prepared, cx) => {
                        prepared.potential_energy(self, u, load, cx)
                    }
                };
                match energy {
                    Ok(energy) => Ok(Some(energy)),
                    Err(SolidError::MaterialRefused { .. }) => Ok(None),
                    Err(error) => Err(error),
                }
            }
            None => Ok(self.try_energy(u, load)),
        }
    }

    fn energy_with_projection(
        &self,
        u: &[f64],
        load: f64,
        projection: Projection<'_, '_>,
    ) -> Result<f64, SolidError> {
        self.try_energy_with_projection(u, load, projection)
            .map(|value| value.unwrap_or(f64::INFINITY))
    }

    /// Public probe for the consistency merge gate: residual and
    /// consistent tangent at an unconstrained state (no pinned DOFs).
    ///
    /// # Errors
    /// [`SolidError::MaterialRefused`] when det F ≤ 0 at `u`.
    pub fn residual_and_tangent(
        &self,
        u: &[f64],
        load: f64,
    ) -> Result<(Vec<f64>, fs_sparse::Csr), SolidError> {
        self.validate_materials()?;
        self.assemble(u, load, &BTreeMap::new())
    }

    /// Unconstrained force residual, without material tangents or sparse assembly.
    ///
    /// Contact line searches and refreshed friction lags need only forces.
    /// This uses the same Piola stresses, quadrature order, and dead tractions
    /// as [`Self::residual_and_tangent`], avoiding nine nested-dual sweeps per
    /// quadrature point when no Newton matrix will be consumed.
    ///
    /// # Errors
    /// [`SolidError::MaterialRefused`] when the material rejects `u`, or
    /// [`SolidError::UnknownPatch`] for an unknown traction patch.
    pub fn unconstrained_residual(&self, u: &[f64], load: f64) -> Result<Vec<f64>, SolidError> {
        self.validate_materials()?;
        self.force_residual(u, load)
    }

    fn force_residual(&self, u: &[f64], load: f64) -> Result<Vec<f64>, SolidError> {
        self.validate_displacement(u)?;
        let mut resid = vec![0.0; u.len()];
        for (element, conn) in self.mesh.elems.iter().enumerate() {
            self.admit_current_element(element, conn, u)?;
            for (point, &(xi, eta, w)) in material_quadrature(conn.len()).iter().enumerate() {
                let material = self.point_material(element, point)?;
                let (_, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                let wq = w * det;
                let f = Self::def_grad(&grads, conn, u);
                let p = material
                    .first_piola(&f)
                    .map_err(|e| SolidError::MaterialRefused {
                        what: format!("{e:?}"),
                    })?;
                for (a, &na) in conn.iter().enumerate() {
                    for i in 0..2 {
                        resid[2 * na + i] +=
                            wq * (p[3 * i] * grads[a][0] + p[3 * i + 1] * grads[a][1]);
                    }
                }
            }
        }
        self.for_traction(load, |dof, v| resid[dof] -= v)?;
        Ok(resid)
    }

    /// Public probe: the total potential energy at a state (`None`
    /// when the material refuses it).
    #[must_use]
    pub fn potential_energy(&self, u: &[f64], load: f64) -> Option<f64> {
        self.validate_materials().ok()?;
        self.try_energy(u, load)
    }

    fn fixed_dofs(&self, load: f64) -> Result<BTreeMap<usize, f64>, SolidError> {
        let mut fixed = BTreeMap::new();
        for (patch, g) in &self.dirichlet {
            if self.mesh.patch_edges(*patch).is_none() {
                return Err(SolidError::UnknownPatch { patch: *patch });
            }
            for node in self.mesh.patch_nodes(*patch) {
                let p = self.mesh.nodes[node];
                let val = g(p[0], p[1]);
                fixed.insert(2 * node, load * val[0]);
                fixed.insert(2 * node + 1, load * val[1]);
            }
        }
        Ok(fixed)
    }

    /// The 3D deformation gradient at a quadrature point (row-major,
    /// F₃₃ = 1 for plane strain).
    fn def_grad(grads: &[[f64; 2]], conn: &[usize], u: &[f64]) -> [f64; 9] {
        let mut h = [[0.0f64; 2]; 2];
        // Partition of unity removes the common displacement. Form relative
        // motion first so a large rigid translation cannot create strain.
        for (a, &node) in conn.iter().enumerate().skip(1) {
            for i in 0..2 {
                let ui = u[2 * node + i] - u[2 * conn[0] + i];
                h[i][0] += ui * grads[a][0];
                h[i][1] += ui * grads[a][1];
            }
        }
        [
            1.0 + h[0][0],
            h[0][1],
            0.0,
            h[1][0],
            1.0 + h[1][1],
            0.0,
            0.0,
            0.0,
            1.0,
        ]
    }

    /// Residual and consistent tangent on the free DOFs (fixed DOFs
    /// pinned with identity rows and zero residual).
    #[allow(clippy::type_complexity)]
    fn assemble(
        &self,
        u: &[f64],
        load: f64,
        fixed: &BTreeMap<usize, f64>,
    ) -> Result<(Vec<f64>, fs_sparse::Csr), SolidError> {
        self.validate_displacement(u)?;
        let ndof = u.len();
        let mut resid = vec![0.0f64; ndof];
        let mut coo = Coo::new(ndof, ndof);
        for (element, conn) in self.mesh.elems.iter().enumerate() {
            self.admit_current_element(element, conn, u)?;
            let nn = conn.len();
            for (point, &(xi, eta, w)) in material_quadrature(nn).iter().enumerate() {
                let material = self.point_material(element, point)?;
                let (_, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                let wq = w * det;
                let f = Self::def_grad(&grads, conn, u);
                let p = material
                    .first_piola(&f)
                    .map_err(|e| SolidError::MaterialRefused {
                        what: format!("{e:?}"),
                    })?;
                let a_mat =
                    material
                        .consistent_tangent(&f)
                        .map_err(|e| SolidError::MaterialRefused {
                            what: format!("{e:?}"),
                        })?;
                for (a, &na) in conn.iter().enumerate() {
                    for i in 0..2 {
                        let ia = 2 * na + i;
                        let ra = wq * (p[3 * i] * grads[a][0] + p[3 * i + 1] * grads[a][1]);
                        if !fixed.contains_key(&ia) {
                            resid[ia] += ra;
                        }
                        for (b, &nb) in conn.iter().enumerate() {
                            for k in 0..2 {
                                let ib = 2 * nb + k;
                                if fixed.contains_key(&ia) || fixed.contains_key(&ib) {
                                    continue;
                                }
                                let mut v = 0.0;
                                for jj in 0..2 {
                                    for ll in 0..2 {
                                        v += a_mat[3 * i + jj][3 * k + ll]
                                            * grads[a][jj]
                                            * grads[b][ll];
                                    }
                                }
                                coo.push(ia, ib, wq * v);
                            }
                        }
                    }
                }
            }
        }
        // Dead-load tractions (ramped).
        self.for_traction(load, |dof, v| {
            if !fixed.contains_key(&dof) {
                resid[dof] -= v;
            }
        })?;
        for &dof in fixed.keys() {
            coo.push(dof, dof, 1.0);
        }
        Ok((resid, coo.assemble()))
    }

    fn for_traction(&self, load: f64, mut sink: impl FnMut(usize, f64)) -> Result<(), SolidError> {
        for (patch, t) in &self.traction {
            let edges = self
                .mesh
                .patch_edges(*patch)
                .ok_or(SolidError::UnknownPatch { patch: *patch })?;
            for &(a, b) in edges {
                let (pa, pb) = (self.mesh.nodes[a], self.mesh.nodes[b]);
                let len = ((pb[0] - pa[0]).powi(2) + (pb[1] - pa[1]).powi(2)).sqrt();
                let g = 0.5 / 3.0f64.sqrt();
                for s in [0.5 - g, 0.5 + g] {
                    let p = [pa[0] + s * (pb[0] - pa[0]), pa[1] + s * (pb[1] - pa[1])];
                    let tv = t(p[0], p[1]);
                    let w = 0.5 * len * load;
                    for (node, shape) in [(a, 1.0 - s), (b, s)] {
                        for (c, tc) in tv.iter().enumerate() {
                            sink(2 * node + c, w * shape * tc);
                        }
                    }
                }
            }
        }
        Ok(())
    }

    /// Total potential energy (None when the material refuses the
    /// state — the line search treats that as "worse").
    fn try_energy(&self, u: &[f64], load: f64) -> Option<f64> {
        self.validate_displacement(u).ok()?;
        let mut internal = 0.0f64;
        for (element, conn) in self.mesh.elems.iter().enumerate() {
            self.admit_current_element(element, conn, u).ok()?;
            for (point, &(xi, eta, w)) in material_quadrature(conn.len()).iter().enumerate() {
                let material = self.point_material(element, point).ok()?;
                let (_, grads, det) = shapes_at(&self.mesh.nodes, conn, xi, eta);
                let f = Self::def_grad(&grads, conn, u);
                let j = f[0] * f[4] - f[1] * f[3];
                if !(j.is_finite() && j > 0.0) {
                    return None;
                }
                // Each original law owns its stretch/domain admission. A
                // permissive neighbour cannot authorize this trial state.
                let potential = if let Some(native) = material.native_hyperelastic() {
                    native.check_state(&f).ok()?;
                    native.energy(&f)
                } else {
                    // The combined operation lets frozen inelastic owners
                    // prepare one return map for both admission and potential.
                    material.stress_and_potential(&f).ok()?.1
                };
                internal += w * det * potential;
            }
        }
        let mut external = 0.0f64;
        let mut ok = true;
        self.for_traction(load, |dof, v| {
            if u.get(dof).is_some() {
                external += v * u[dof];
            } else {
                ok = false;
            }
        })
        .ok()?;
        let potential = internal - external;
        (ok && potential.is_finite()).then_some(potential)
    }
}

/// Continuation plumbing: a [`HyperProblem`] with HOMOGENEOUS
/// Dirichlet data (clamps) is a [`crate::continuation::PathResidual`]
/// — the load parameter scales the dead tractions.
impl crate::continuation::PathResidual for HyperProblem<'_> {
    fn ndof(&self) -> usize {
        2 * self.mesh.node_count()
    }

    fn residual(&self, u: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.validate_materials()?;
        let fixed = self.fixed_dofs(lambda)?;
        Ok(self.assemble(u, lambda, &fixed)?.0)
    }

    fn tangent(&self, u: &[f64], lambda: f64) -> Result<fs_sparse::Csr, SolidError> {
        self.validate_materials()?;
        let fixed = self.fixed_dofs(lambda)?;
        Ok(self.assemble(u, lambda, &fixed)?.1)
    }

    fn load_vector(&self) -> Vec<f64> {
        let mut f = vec![0.0f64; 2 * self.mesh.node_count()];
        let fixed = self.fixed_dofs(1.0).unwrap_or_default();
        let _ = self.for_traction(1.0, |dof, v| {
            if !fixed.contains_key(&dof) {
                f[dof] += v;
            }
        });
        f
    }
}

// Eliminate supports only after retaining their consistent contribution to
// the free equations: K_ff du_f = -(r_f + K_fp du_p). The original full
// tangent is needed only while a boundary increment remains pending.
fn eliminate_prescribed_increment(
    tangent: fs_sparse::Csr,
    residual: &mut [f64],
    pending: &BTreeMap<usize, f64>,
) -> fs_sparse::Csr {
    let mut constrained = Coo::new(tangent.nrows(), tangent.ncols());
    for (row, rhs) in residual.iter_mut().enumerate() {
        if pending.contains_key(&row) {
            *rhs = 0.0;
            constrained.push(row, row, 1.0);
            continue;
        }
        let (columns, values) = tangent.row(row);
        for (&column, &coefficient) in columns.iter().zip(values) {
            if let Some(&increment) = pending.get(&column) {
                *rhs = coefficient.mul_add(increment, *rhs);
            } else {
                constrained.push(row, column, coefficient);
            }
        }
    }
    constrained.assemble()
}

fn normalized_newton_operator(
    mut tangent: fs_sparse::Csr,
    stiffness_scale: f64,
    fixed: &BTreeMap<usize, f64>,
) -> CsrOp {
    // Scale the coefficients BEFORE multiplication. Dividing K*x afterward
    // cannot recover a finite normalized result once the unscaled product or
    // row sum has overflowed (or underflowed). The matrix is already owned by
    // this Newton step, so no second sparse allocation or per-apply division
    // is needed. Support rows are dimensionless identities, not stiffness;
    // keep them at one even when the free coefficients are subnormal.
    for row in 0..tangent.nrows() {
        if !fixed.contains_key(&row) {
            for value in tangent.row_values_mut(row) {
                *value /= stiffness_scale;
            }
        }
    }
    CsrOp::symmetric(tangent)
}

fn norm(v: &[f64]) -> f64 {
    if v.iter().any(|value| !value.is_finite()) {
        return f64::NAN;
    }
    let scale = v.iter().map(|value| value.abs()).fold(0.0_f64, f64::max);
    if scale == 0.0 {
        return 0.0;
    }
    let sum: f64 = v.iter().map(|value| (value / scale).powi(2)).sum();
    scale * fs_math::det::sqrt(sum)
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::hyper::{FibreInvariants, HyperelasticModel, OgdenTerm, Volumetric};

    /// G0: the known support increment enters the SAME tangent with its actual
    /// sign; reaction rows cannot become extra applied forces or solve rows.
    #[test]
    fn prescribed_increment_elimination_retains_free_force_and_symmetric_operator() {
        let mut coo = Coo::new(3, 3);
        for (row, values) in [[4.0, -2.0, 1.0], [-2.0, 6.0, -3.0], [1.0, -3.0, 5.0]]
            .into_iter()
            .enumerate()
        {
            for (column, value) in values.into_iter().enumerate() {
                coo.push(row, column, value);
            }
        }
        let pending = BTreeMap::from([(1, -0.25)]);
        let mut residual = [2.0, 100.0, -1.0];
        let tangent = eliminate_prescribed_increment(coo.assemble(), &mut residual, &pending);
        assert_eq!(residual, [2.5, 0.0, -0.25]);
        let mut action = [0.0; 3];
        CsrOp::symmetric(tangent).apply(&[2.0, -8.0, 3.0], &mut action);
        assert_eq!(action, [11.0, -8.0, 17.0]);
    }

    /// G0/G1: quadrature-positive elasticity must not authorize a folded body,
    /// including when every displacement is prescribed and the free residual
    /// would otherwise vanish without evaluating an energy line search.
    #[test]
    fn plane_strain_rejects_fold_hidden_between_gauss_points_and_boundary() {
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 5.0,
            },
            100.0,
        )
        .unwrap();
        let prescribed = |x: f64, y: f64| [1.2 * (1.0 - x) * y, 0.0];
        let displacement: Vec<f64> = mesh
            .nodes
            .iter()
            .flat_map(|&[x, y]| prescribed(x, y))
            .collect();
        let problem = HyperProblem {
            mesh: &mesh,
            material: &material,
            dirichlet: vec![
                (Patch::Left, &prescribed),
                (Patch::Right, &prescribed),
                (Patch::Top, &prescribed),
                (Patch::Bottom, &prescribed),
            ],
            traction: vec![],
            settings: NewtonSettings::default(),
        };
        // Independent map: x' = x + 1.2 (1-x)y, y'=y, det F=1-1.2y.
        // Both Gauss y coordinates are below 5/6; the entire top edge is
        // reversed. The original pointwise material calls all admit this.
        for (xi, eta, _) in quad_points(4) {
            let (_, gradients, _) = shapes_at(&mesh.nodes, &mesh.elems[0], xi, eta);
            let f = HyperProblem::def_grad(&gradients, &mesh.elems[0], &displacement);
            let expected = 1.0 - 1.2 * (0.5 * (eta + 1.0));
            assert!(expected > 0.0);
            assert!((f[0] * f[4] - f[1] * f[3] - expected).abs() < 1e-14);
            // Both imported traits name first_piola; the FiniteStrainPotential blanket impl
            // forwards to HyperelasticMaterial's, so naming the trait changes no call.
            FiniteStrainPotential::first_piola(&material, &f).unwrap();
        }
        assert!(matches!(
            problem.unconstrained_residual(&displacement, 1.0),
            Err(SolidError::MaterialRefused { .. })
        ));
        assert!(matches!(
            problem.residual_and_tangent(&displacement, 1.0),
            Err(SolidError::MaterialRefused { .. })
        ));
        assert!(problem.potential_energy(&displacement, 1.0).is_none());
        assert!(matches!(
            problem.solve(),
            Err(SolidError::MaterialRefused { .. })
        ));
    }

    /// G0: original reference charts cannot introduce negative elastic storage
    /// or panic through connectivity before the current material is evaluated.
    #[test]
    fn plane_strain_rejects_invalid_reference_charts_and_incomplete_motion() {
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 5.0,
            },
            3.0,
        )
        .unwrap();
        for kind in 0..3 {
            let mut mesh = Mesh2::quads(1.0, 1.0, 1, 1);
            match kind {
                0 => mesh.elems[0].reverse(),
                1 => mesh.elems[0][0] = mesh.nodes.len(),
                _ => mesh.nodes[2][0] = 1.2,
            }
            let problem = HyperProblem {
                mesh: &mesh,
                material: &material,
                dirichlet: vec![],
                traction: vec![],
                settings: NewtonSettings::default(),
            };
            let u = vec![0.0; 2 * mesh.node_count()];
            assert!(matches!(
                problem.residual_and_tangent(&u, 0.0),
                Err(SolidError::InvalidInput { .. })
            ));
            assert!(problem.potential_energy(&u, 0.0).is_none());
        }
        let mesh = Mesh2::quads(1.0, 1.0, 1, 1);
        let problem = HyperProblem {
            mesh: &mesh,
            material: &material,
            dirichlet: vec![],
            traction: vec![],
            settings: NewtonSettings::default(),
        };
        for u in [vec![0.0; 7], vec![0.0; 9], vec![f64::NAN; 8]] {
            assert!(matches!(
                problem.unconstrained_residual(&u, 0.0),
                Err(SolidError::InvalidInput { .. })
            ));
            assert!(matches!(
                problem.residual_and_tangent(&u, 0.0),
                Err(SolidError::InvalidInput { .. })
            ));
            assert!(problem.potential_energy(&u, 0.0).is_none());
        }
    }

    /// G3: force-unit normalization must precede products, including the true
    /// residual action on a recovered Krylov iterate larger than one.
    #[test]
    fn normalized_newton_action_preserves_finite_products_and_supports() {
        for scale in [1e-310, 1.0, 1e308] {
            let mut tangent = Coo::new(3, 3);
            tangent.push(0, 0, scale);
            tangent.push(0, 1, -scale);
            tangent.push(1, 0, -scale);
            tangent.push(1, 1, scale);
            tangent.push(2, 2, 1.0);
            let fixed = BTreeMap::from([(2, 0.0)]);
            let op = normalized_newton_operator(tangent.assemble(), scale, &fixed);
            let mut actual = [0.0; 3];
            op.apply(&[2.0, -2.0, 0.0], &mut actual);
            assert_eq!(actual, [4.0, -4.0, 0.0]);
            // Tiny products would round to zero in the physical matrix.
            op.apply(&[1e-20, 0.0, 3.0], &mut actual);
            assert_eq!(actual, [1e-20, -1e-20, 3.0]);
            let mut transpose = [0.0; 3];
            op.apply_transpose(&[1e-20, 0.0, 3.0], &mut transpose);
            assert_eq!(actual, transpose);
        }
    }

    /// A symmetric indefinite tangent remains usable near a limit point.
    /// Both the Krylov path and the final true residual use the same scaling.
    #[test]
    fn normalized_newton_krylov_solves_across_extreme_force_units() {
        for scale in [1e-300, 1.0, 1e300] {
            let mut tangent = Coo::new(3, 3);
            tangent.push(0, 0, scale);
            tangent.push(0, 1, scale);
            tangent.push(1, 0, scale);
            tangent.push(1, 1, -scale);
            tangent.push(2, 2, 1.0);
            let fixed = BTreeMap::from([(2, 0.0)]);
            let op = normalized_newton_operator(tangent.assemble(), scale, &fixed);
            let rhs = [0.6, 0.8, 0.0];
            let mut state = MinresState::new(&op, &rhs);
            state.run(&op, 1e-12, 8);
            let mut actual = [0.0; 3];
            op.apply(&state.x, &mut actual);
            for (a, b) in actual.iter().zip(rhs) {
                assert!((a - b).abs() < 2e-14, "scale={scale:e}, action={actual:?}");
            }
            assert!((state.x[0] - 0.7).abs() < 2e-14);
            assert!((state.x[1] + 0.1).abs() < 2e-14);
            assert_eq!(state.x[2], 0.0);
        }
    }

    /// The cheaper force probe must preserve the fully assembled residual bits.
    #[test]
    fn residual_only_matches_all_materials_on_triangles_and_quads() {
        let models = [
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 5.0,
            },
            HyperelasticModel::MooneyRivlin {
                c10: 0.7,
                c01: 0.3,
                kappa: 8.0,
            },
            HyperelasticModel::Ogden {
                terms: vec![OgdenTerm {
                    mu_pa: 2.0,
                    alpha: 2.0,
                }],
                volumetric: Volumetric::LogSquared { kappa: 8.0 },
            },
            HyperelasticModel::HolzapfelGasserOgden {
                c10: 0.7,
                k1: 1.2,
                k2: 2.0,
                dispersion: 0.1,
                fibres: vec![[1.0, 0.0, 0.0]],
                invariants: FibreInvariants::Full,
                volumetric: Volumetric::Quadratic { kappa: 8.0 },
            },
        ];
        let traction = |x: f64, y: f64| [0.2 + 0.1 * y, -0.15 * x];
        let prescribed = |_: f64, _: f64| [0.03, -0.01];
        for mesh in [
            Mesh2::quads(1.0, 0.7, 2, 2),
            Mesh2::triangles(1.0, 0.7, 2, 2),
        ] {
            for model in &models {
                let material = Hyperelastic::new(model.clone(), 3.0).unwrap();
                let problem = HyperProblem {
                    mesh: &mesh,
                    material: &material,
                    // Both public probes intentionally leave Dirichlet DOFs free.
                    dirichlet: vec![(Patch::Left, &prescribed)],
                    traction: vec![(Patch::Top, &traction), (Patch::Right, &traction)],
                    settings: NewtonSettings::default(),
                };
                for scale in [0.0, 0.1, -0.15] {
                    let u: Vec<f64> = mesh
                        .nodes
                        .iter()
                        .flat_map(|&[x, y]| {
                            [
                                scale * (0.2 * x + 0.1 * y + 0.15 * x * y),
                                scale * (-0.1 * x + 0.3 * y - 0.1 * x * y),
                            ]
                        })
                        .collect();
                    for load in [0.0, 0.6] {
                        let actual = problem.unconstrained_residual(&u, load).unwrap();
                        let expected = problem.residual_and_tangent(&u, load).unwrap().0;
                        assert_eq!(
                            actual.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                            expected.iter().map(|v| v.to_bits()).collect::<Vec<_>>()
                        );
                    }
                }
                let inverted: Vec<f64> = mesh
                    .nodes
                    .iter()
                    .flat_map(|&[x, _]| [-2.0 * x, 0.0])
                    .collect();
                assert!(matches!(
                    problem.unconstrained_residual(&inverted, 0.6),
                    Err(SolidError::MaterialRefused { .. })
                ));
                assert!(matches!(
                    problem.residual_and_tangent(&inverted, 0.6),
                    Err(SolidError::MaterialRefused { .. })
                ));
            }
        }
    }

    /// G3: rigid translations leave every public constitutive probe unchanged.
    #[test]
    fn skew_elements_preserve_force_energy_and_tangent_under_large_translation() {
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2.0,
                lambda: 5.0,
            },
            3.0,
        )
        .unwrap();

        for mut mesh in [
            Mesh2::quads(1.0, 1.0, 1, 1),
            Mesh2::triangles(1.0, 1.0, 1, 1),
        ] {
            // Dyadic geometry and motion survive the translation exactly, while
            // the inverse element Jacobian has non-dyadic shape gradients.
            for x in &mut mesh.nodes {
                *x = [1.125 * x[0] + 0.375 * x[1], 0.25 * x[0] + 0.875 * x[1]];
            }
            let problem = HyperProblem {
                mesh: &mesh,
                material: &material,
                dirichlet: vec![],
                traction: vec![],
                settings: NewtonSettings::default(),
            };

            for scale in [0.0, 1.0] {
                let u: Vec<f64> = mesh
                    .nodes
                    .iter()
                    .flat_map(|&[x, y]| {
                        [
                            scale * (0.125 * x + 0.25 * y),
                            scale * (-0.25 * x + 0.125 * y),
                        ]
                    })
                    .collect();
                let translation = [1_099_511_627_776.0, -2_199_023_255_552.0];
                let translated: Vec<f64> = u
                    .iter()
                    .enumerate()
                    .map(|(i, value)| value + translation[i % 2])
                    .collect();

                let (force, tangent) = problem.residual_and_tangent(&u, 0.0).unwrap();
                let (shifted_force, shifted_tangent) =
                    problem.residual_and_tangent(&translated, 0.0).unwrap();
                assert_eq!(force, shifted_force);
                assert_eq!(
                    problem.unconstrained_residual(&translated, 0.0).unwrap(),
                    force
                );
                assert_eq!(
                    problem.potential_energy(&u, 0.0).unwrap().to_bits(),
                    problem
                        .potential_energy(&translated, 0.0)
                        .unwrap()
                        .to_bits()
                );
                for i in 0..u.len() {
                    for j in 0..u.len() {
                        assert_eq!(
                            tangent.get(i, j).to_bits(),
                            shifted_tangent.get(i, j).to_bits()
                        );
                    }
                }
            }
        }
    }
}
