//! Native reduced-plane equilibrium paths through load-control limit points.
//!
//! The existing bordered continuation owner follows the SAME finite-bulk
//! plane-stress, exact-volume plane-stress or Q1 mean-dilatation potential.
//! Prepared reference geometry, material-point owners, thickness relaxation and
//! consistent condensed Hessians are retained; no second constitutive law or
//! prestress matrix is introduced.
//!
//! Free path coordinates are q=u_free/L. Prescribed coordinates remain the
//! original absolute ramp u_fixed=lambda*g, while dead tractions scale by lambda.
//! Thus R=h*(f_int-lambda*f_dead)/F and -R_lambda=h*(f_dead-K_fc*g)/F.
//! The support term is evaluated from the current COMPLETE native tangent.
//! Supports and reaction rows never enter the free arclength metric.
//!
//! Inputs and callbacks describe one held elastic problem. Temperature,
//! reference maps and material history do not evolve along the path; frozen
//! incremental inelastic potentials refuse. Negative tangents are permitted.
//! Path events are numerical diagnostics, not certified bifurcation locations,
//! global injectivity, continuum stability, bending or wrinkling claims.

use std::cell::RefCell;
use std::collections::BTreeMap;

use super::HyperProblem;
use super::mean_dilatation;
use super::plane_stress::{self, Constraint};
use crate::SolidError;
use crate::continuation::{self, ArcSettings, PathResidual, PathState};
use fs_exec::Cx;
use fs_material::hyper::plane_stress::PlaneStressSettings;
use fs_sparse::Csr;

/// Original native reduction used at every predictor and corrector state.
#[derive(Clone, Copy, Debug)]
pub enum PlanePathFormulation {
    PlaneStress(PlaneStressSettings),
    IncompressiblePlaneStress(PlaneStressSettings),
    /// Native Q1 plane-strain mean dilatation and its original material limits.
    MeanDilatation,
}

/// Positive physical scales of the dimensionless residual and free coordinates.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct PlanePathScales {
    pub displacement_m: f64,
    pub force_n: f64,
}

/// Physical response of one path state. The residual remains explicit even for
/// a caller-edited state; observing a state does not certify equilibrium.
#[derive(Clone, Debug, PartialEq)]
pub struct PlanePathResponse {
    pub load_factor: f64,
    pub displacement_m: Vec<[f64; 2]>,
    /// Internal minus applied force on prescribed DOFs, in N.
    pub reactions_n: Vec<(usize, f64)>,
    pub free_force_residual_n: Vec<f64>,
    /// Original material stored energy, including its authored energy datum.
    pub strain_energy_j: f64,
    /// Stored energy minus the current dead-load work.
    /// Work performed by moving supports is not a dead-load potential.
    pub potential_energy_j: f64,
}

enum Prepared {
    Plane(plane_stress::Prepared<'static>),
    Mean(mean_dilatation::Prepared),
}

impl Prepared {
    fn assemble(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        load: f64,
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, Csr), SolidError> {
        match self {
            Self::Plane(prepared) => prepared.assemble(problem, u, load, &BTreeMap::new(), cx),
            Self::Mean(prepared) => prepared.assemble(problem, u, load, &BTreeMap::new(), cx),
        }
    }
    fn energy(
        &self,
        problem: &HyperProblem<'_>,
        u: &[f64],
        cx: &Cx<'_>,
    ) -> Result<f64, SolidError> {
        match self {
            Self::Plane(prepared) => prepared.potential_energy(problem, u, 0.0, cx),
            Self::Mean(prepared) => prepared.potential_energy(problem, u, 0.0, cx),
        }
    }
}

/// Immutable native source binding and original prescribed-motion ramp.
pub struct PlaneContinuation<'p, 'm> {
    problem: &'p HyperProblem<'m>,
    formulation: PlanePathFormulation,
    prepared: Prepared,
    thickness_m: f64,
    scales: PlanePathScales,
    force_scale: f64,
    tangent_scale: f64,
    prescribed: BTreeMap<usize, f64>,
    free: Vec<usize>,
    reduced: Vec<usize>,
    reference_force: Vec<f64>,
    reference_force_n: Vec<f64>,
}

fn invalid(what: &str) -> SolidError {
    SolidError::InvalidInput { what: what.into() }
}
fn checkpoint(cx: &Cx<'_>) -> Result<(), SolidError> {
    cx.checkpoint().map_err(|_| SolidError::Cancelled)
}
fn scaled(value: f64, scale: f64) -> Result<f64, SolidError> {
    let result = value * scale;
    if !result.is_finite() || (value != 0.0 && scale != 0.0 && result == 0.0) {
        Err(invalid(
            "unrepresentable plane continuation force, tangent or displacement",
        ))
    } else {
        Ok(result)
    }
}
fn norm(values: &[f64]) -> f64 {
    values
        .iter()
        .fold(0.0_f64, |norm, &value| norm.hypot(value))
}

struct Linearization {
    residual: Vec<f64>,
    tangent: Csr,
    load_column: Vec<f64>,
}

impl<'m> HyperProblem<'m> {
    /// Prepare native reduced-plane continuation with original support/traction
    /// ramps. All boundary and material callbacks must remain pure while held.
    ///
    /// The reference thickness is physical and fixed. Scales and their h/F and
    /// h*L/F ratios must be finite, positive and representable. The explicit
    /// free-DOF cap is checked before reference-geometry preparation; native
    /// point/element work and the existing bordered solver retain their limits.
    /// There must be free coordinates and a nonzero dead-force or support ramp.
    /// This does not infer a material history path or reference-state update.
    pub fn prepare_continuation<'p>(
        &'p self,
        formulation: PlanePathFormulation,
        reference_thickness_m: f64,
        scales: PlanePathScales,
        maximum_free_dofs: usize,
        cx: &Cx<'_>,
    ) -> Result<PlaneContinuation<'p, 'm>, SolidError> {
        checkpoint(cx)?;
        let force_scale = reference_thickness_m / scales.force_n;
        let tangent_scale = force_scale * scales.displacement_m;
        if [
            reference_thickness_m,
            scales.displacement_m,
            scales.force_n,
            force_scale,
            tangent_scale,
        ]
        .iter()
        .any(|value| !value.is_finite() || *value <= 0.0)
        {
            return Err(invalid(
                "positive representable plane path thickness and scales required",
            ));
        }
        let ndof = self
            .mesh
            .nodes
            .len()
            .checked_mul(2)
            .ok_or_else(|| invalid("plane path displacement count overflows"))?;
        if self
            .mesh
            .patches
            .iter()
            .flat_map(|(_, edges)| edges)
            .any(|&(a, b)| a >= self.mesh.nodes.len() || b >= self.mesh.nodes.len())
        {
            return Err(invalid("plane path boundary references an absent node"));
        }
        let prescribed = self.fixed_dofs(1.0)?;
        if prescribed.values().any(|value| !value.is_finite()) {
            return Err(invalid("finite prescribed plane motion required"));
        }
        let free: Vec<_> = (0..ndof)
            .filter(|dof| !prescribed.contains_key(dof))
            .collect();
        if free.is_empty() || free.len() > maximum_free_dofs {
            return Err(invalid(
                "nonempty plane path free coordinates exceed admission",
            ));
        }
        let prepared = match formulation {
            PlanePathFormulation::PlaneStress(local) => Prepared::Plane(
                plane_stress::Prepared::new(self, local, Constraint::FiniteBulk, cx)?,
            ),
            PlanePathFormulation::IncompressiblePlaneStress(local) => Prepared::Plane(
                plane_stress::Prepared::new(self, local, Constraint::Incompressible, cx)?,
            ),
            PlanePathFormulation::MeanDilatation => {
                Prepared::Mean(mean_dilatation::Prepared::new(self, cx)?)
            }
        };
        for (element, _) in self.mesh.elems.iter().enumerate() {
            checkpoint(cx)?;
            for point in 0..self.material.material_point_count(element).unwrap_or(1) {
                if self
                    .point_material(element, point)?
                    .elastic_potential()
                    .is_none()
                {
                    return Err(invalid(
                        "plane continuation requires elastic stored energy, not a frozen incremental potential",
                    ));
                }
            }
        }
        let mut reduced = vec![usize::MAX; ndof];
        for (index, &dof) in free.iter().enumerate() {
            reduced[dof] = index;
        }
        let mut reference_force = vec![0.0; ndof];
        self.for_traction(1.0, |dof, value| reference_force[dof] += value)?;
        let reference_force_n = reference_force
            .iter()
            .map(|&force| scaled(force, reference_thickness_m))
            .collect::<Result<Vec<_>, _>>()?;
        if !free.iter().any(|&dof| reference_force[dof] != 0.0)
            && !prescribed.values().any(|&value| value != 0.0)
        {
            return Err(invalid(
                "plane continuation needs a nonzero free dead force or prescribed motion ramp",
            ));
        }
        checkpoint(cx)?;
        Ok(PlaneContinuation {
            problem: self,
            formulation,
            prepared,
            thickness_m: reference_thickness_m,
            scales,
            force_scale,
            tangent_scale,
            prescribed,
            free,
            reduced,
            reference_force,
            reference_force_n,
        })
    }
}

impl<'p, 'm> PlaneContinuation<'p, 'm> {
    pub fn problem(&self) -> &HyperProblem<'m> {
        self.problem
    }
    pub const fn formulation(&self) -> PlanePathFormulation {
        self.formulation
    }
    pub const fn scales(&self) -> PlanePathScales {
        self.scales
    }
    pub const fn reference_thickness_m(&self) -> f64 {
        self.thickness_m
    }
    /// Increasing reduced-to-full node-major xy map; fixed rows are absent.
    pub fn free_dofs(&self) -> &[usize] {
        &self.free
    }
    /// Original complete reference dead force, including force at supports.
    pub fn reference_force_n(&self) -> &[f64] {
        &self.reference_force_n
    }
    /// Original complete support targets at lambda=1, in m.
    pub fn prescribed_motion_m(&self) -> &BTreeMap<usize, f64> {
        &self.prescribed
    }

    fn displacement(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        if q.len() != self.free.len() || !lambda.is_finite() {
            return Err(invalid(
                "finite dimension-matched plane path state required",
            ));
        }
        let mut u = vec![0.0; self.reduced.len()];
        for (&dof, &value) in self.free.iter().zip(q) {
            u[dof] = scaled(value, self.scales.displacement_m)?;
        }
        for (&dof, &target) in &self.prescribed {
            u[dof] = scaled(target, lambda)?;
        }
        for &force in &self.reference_force_n {
            scaled(force, lambda)?;
        }
        Ok(u)
    }

    fn linearization(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Linearization, SolidError> {
        checkpoint(cx)?;
        let u = self.displacement(q, lambda)?;
        let (full_residual, full_tangent) = self.prepared.assemble(self.problem, &u, lambda, cx)?;
        let n = self.free.len();
        let mut residual = Vec::with_capacity(n);
        let mut load_column = Vec::with_capacity(n);
        let mut rows = Vec::with_capacity(n + 1);
        let mut columns = Vec::new();
        let mut values = Vec::new();
        rows.push(0);
        for &row in &self.free {
            checkpoint(cx)?;
            residual.push(scaled(full_residual[row], self.force_scale)?);
            let mut load = self.reference_force[row];
            let (native_columns, native_values) = full_tangent.row(row);
            for (&column, &value) in native_columns.iter().zip(native_values) {
                let index = self.reduced[column];
                if index == usize::MAX {
                    // q holds free displacement fixed while lambda moves supports.
                    // Native K_fc is retained BEFORE eliminating prescribed DOFs.
                    load = (-self.prescribed[&column]).mul_add(value, load);
                } else {
                    columns.push(index);
                    values.push(scaled(value, self.tangent_scale)?);
                }
            }
            load_column.push(scaled(load, self.force_scale)?);
            rows.push(values.len());
        }
        let tangent =
            Csr::try_from_parts_with_checkpoint(n, n, rows, columns, values, || checkpoint(cx))?
                .ok_or_else(|| invalid("noncanonical native plane path restriction"))?;
        checkpoint(cx)?;
        Ok(Linearization {
            residual,
            tangent,
            load_column,
        })
    }

    /// Dimensionless native free-force residual at current q and lambda.
    pub fn residual(&self, q: &[f64], lambda: f64, cx: &Cx<'_>) -> Result<Vec<f64>, SolidError> {
        Ok(self.linearization(q, lambda, cx)?.residual)
    }
    /// Exact q derivative, including the original transverse/volumetric condensation.
    pub fn tangent(&self, q: &[f64], lambda: f64, cx: &Cx<'_>) -> Result<Csr, SolidError> {
        Ok(self.linearization(q, lambda, cx)?.tangent)
    }
    /// Exact -dR/dlambda at fixed q, including the original moving supports.
    pub fn load_vector_at(
        &self,
        q: &[f64],
        lambda: f64,
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, SolidError> {
        Ok(self.linearization(q, lambda, cx)?.load_column)
    }

    /// Start at zero displacement/load only if that is an admitted equilibrium.
    /// A held prestrain/reference may instead require start_from.
    pub fn start(&self, settings: &ArcSettings, cx: &Cx<'_>) -> Result<PathState, SolidError> {
        self.start_from(
            &vec![[0.0; 2]; self.problem.mesh.nodes.len()],
            0.0,
            settings,
            cx,
        )
    }

    /// Start from a complete physical equilibrium without changing its supports.
    /// Prescribed entries must equal lambda*g exactly. The initial path tangent
    /// points toward increasing lambda; resumed states retain their own orientation.
    pub fn start_from(
        &self,
        displacement_m: &[[f64; 2]],
        lambda: f64,
        settings: &ArcSettings,
        cx: &Cx<'_>,
    ) -> Result<PathState, SolidError> {
        checkpoint(cx)?;
        let u = displacement_m.as_flattened();
        if u.len() != self.reduced.len()
            || !lambda.is_finite()
            || u.iter().any(|value| !value.is_finite())
            || self
                .prescribed
                .iter()
                .any(|(&dof, &target)| u[dof] != lambda * target)
        {
            return Err(invalid(
                "complete finite path displacement with exact ramped supports required",
            ));
        }
        let mut state = PathState::start(self.free.len(), settings);
        state.lambda = lambda;
        state.u = self
            .free
            .iter()
            .map(|&dof| {
                let q = u[dof] / self.scales.displacement_m;
                if !q.is_finite() || (u[dof] != 0.0 && q == 0.0) {
                    Err(invalid("unrepresentable initial plane path coordinate"))
                } else {
                    Ok(q)
                }
            })
            .collect::<Result<Vec<_>, _>>()?;
        state.trace = vec![(lambda, state.u.iter().map(|v| v.abs()).fold(0.0, f64::max))];
        self.advance(&state, settings, 0, cx)
    }

    /// Return a complete candidate; failed steps/cancellation leave previous
    /// unchanged. Resumption first rechecks the original physical force equation.
    pub fn advance(
        &self,
        previous: &PathState,
        settings: &ArcSettings,
        steps: usize,
        cx: &Cx<'_>,
    ) -> Result<PathState, SolidError> {
        checkpoint(cx)?;
        let adapter = PathAdapter {
            owner: self,
            cx,
            cache: RefCell::new(None),
        };
        let mut candidate = previous.clone();
        continuation::advance(&adapter, &mut candidate, settings, 0)?;
        if norm(&adapter.residual(&previous.u, previous.lambda)?) > settings.tol {
            return Err(invalid(
                "initial plane path state does not satisfy its original force equation",
            ));
        }
        continuation::advance(&adapter, &mut candidate, settings, steps)?;
        checkpoint(cx)?;
        Ok(candidate)
    }

    /// Native physical forces and elastic energy, with no stability relabeling.
    pub fn response(
        &self,
        state: &PathState,
        cx: &Cx<'_>,
    ) -> Result<PlanePathResponse, SolidError> {
        checkpoint(cx)?;
        let u = self.displacement(&state.u, state.lambda)?;
        let (residual, _) = self.prepared.assemble(self.problem, &u, state.lambda, cx)?;
        let strain_energy_j = scaled(
            self.prepared.energy(self.problem, &u, cx)?,
            self.thickness_m,
        )?;
        let mut dead_work = 0.0;
        for (&force, &displacement) in self.reference_force_n.iter().zip(&u) {
            dead_work = force.mul_add(displacement, dead_work);
        }
        let potential_energy_j = (-state.lambda).mul_add(dead_work, strain_energy_j);
        if !potential_energy_j.is_finite() {
            return Err(invalid("plane path potential overflow"));
        }
        let reactions_n = self
            .prescribed
            .keys()
            .map(|&dof| Ok((dof, scaled(residual[dof], self.thickness_m)?)))
            .collect::<Result<Vec<_>, SolidError>>()?;
        let free_force_residual_n = self
            .free
            .iter()
            .map(|&dof| scaled(residual[dof], self.thickness_m))
            .collect::<Result<Vec<_>, _>>()?;
        checkpoint(cx)?;
        Ok(PlanePathResponse {
            load_factor: state.lambda,
            displacement_m: u.chunks_exact(2).map(|u| [u[0], u[1]]).collect(),
            reactions_n,
            free_force_residual_n,
            strain_energy_j,
            potential_energy_j,
        })
    }
}

struct Cached {
    q: Vec<f64>,
    lambda: f64,
    value: Linearization,
}
struct PathAdapter<'s, 'p, 'm, 'c, 'cx> {
    owner: &'s PlaneContinuation<'p, 'm>,
    cx: &'c Cx<'cx>,
    cache: RefCell<Option<Cached>>,
}
impl PathAdapter<'_, '_, '_, '_, '_> {
    // The bordered owner asks for residual, K and parameter column at the same
    // state. Reuse ONE native constitutive reduction/assembly per exact state.
    fn with_linearization<T>(
        &self,
        q: &[f64],
        lambda: f64,
        read: impl FnOnce(&Linearization) -> T,
    ) -> Result<T, SolidError> {
        checkpoint(self.cx)?;
        let hit = self.cache.borrow().as_ref().is_some_and(|cache| {
            cache.lambda.to_bits() == lambda.to_bits()
                && cache.q.len() == q.len()
                && cache
                    .q
                    .iter()
                    .zip(q)
                    .all(|(a, b)| a.to_bits() == b.to_bits())
        });
        if !hit {
            let value = self.owner.linearization(q, lambda, self.cx)?;
            *self.cache.borrow_mut() = Some(Cached {
                q: q.to_vec(),
                lambda,
                value,
            });
        }
        Ok(read(
            &self
                .cache
                .borrow()
                .as_ref()
                .expect("prepared exact state")
                .value,
        ))
    }
}
impl PathResidual for PathAdapter<'_, '_, '_, '_, '_> {
    fn checkpoint(&self) -> Result<(), SolidError> {
        checkpoint(self.cx)
    }
    fn ndof(&self) -> usize {
        self.owner.free.len()
    }
    fn residual(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.with_linearization(q, lambda, |value| value.residual.clone())
    }
    fn tangent(&self, q: &[f64], lambda: f64) -> Result<Csr, SolidError> {
        self.with_linearization(q, lambda, |value| value.tangent.clone())
    }
    fn load_vector(&self) -> Vec<f64> {
        self.owner
            .free
            .iter()
            .map(|&dof| self.owner.reference_force[dof] * self.owner.force_scale)
            .collect()
    }
    fn load_vector_at(&self, q: &[f64], lambda: f64) -> Result<Vec<f64>, SolidError> {
        self.with_linearization(q, lambda, |value| value.load_column.clone())
    }
}
