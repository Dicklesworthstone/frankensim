//! Stationary native gravity on the evolving finite solid.
//!
//! The same reference mass supplies inertia and current-position quadrature.
//! The shared field supplies potential, acceleration and its analytic spatial
//! derivative. Original generalized-alpha/Newton owns the motion; this module
//! only binds that existing force/tangent/potential into its existing operator.

/// Reciprocal tidal impulses on the original solid and a finite point source.
pub mod reciprocal;

use std::fmt;
use std::sync::Arc;

use fs_math::dd::Dd;
use fs_mbd::gravity::PointGravityField;
use fs_mbd_material::density_field::TetrahedralDensityBody;
use fs_solid::hyper3::gravity::TetGravityLoad3;
use fs_solid::hyper3::pressure::{ClosedTetPressureSurface, TetPressure};
use fs_sparse::{Coo, Csr};

use super::{
    Cx, DynamicOperator, FiniteSolidBody, FiniteSolidError, FiniteSolidInitialization,
    FiniteSolidMaterial, HyperTetError, HyperTetProblem, HyperelasticLawField,
    PreparedHyperTetProblem, admit_reference_density_field, checkpoint,
};

/// One borrowed stationary conservative field in the solid's world frame.
///
/// The original field owns its domain and derivatives. It must remain
/// stationary and deterministic for the lifetime of this binding; a mutable
/// source trajectory needs its own external-work law. No Earth preset,
/// rotating-frame force, source recoil or solid self-gravity is inferred.
#[derive(Clone, Copy)]
pub struct StationaryGravity<'a> {
    field: &'a dyn PointGravityField,
    frame: u64,
}

impl<'a> StationaryGravity<'a> {
    /// Bind the existing field; actual material samples are admitted at preparation.
    pub const fn new(field: &'a dyn PointGravityField, frame: u64) -> Self {
        Self { field, frame }
    }

    /// Explicit inertial world frame used by every current-position query.
    pub const fn frame(self) -> u64 {
        self.frame
    }

    /// The original field owner; no copied or frozen acceleration field.
    pub const fn field(self) -> &'a dyn PointGravityField {
        self.field
    }
}

impl fmt::Debug for StationaryGravity<'_> {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        f.debug_struct("StationaryGravity")
            .field("frame", &self.frame)
            .finish_non_exhaustive()
    }
}

/// The original field's actual potential and midpoint force over one step.
#[derive(Clone, Debug)]
pub struct GravityStep {
    pub frame: u64,
    pub potential_before_j: f64,
    pub potential_after_j: f64,
    /// Minus the change of the original quadrature potential. Included exactly
    /// once in the parent source-work account; do not add it again as storage.
    pub work_j: f64,
    /// Actual full nodal force used by the native midpoint equation, including
    /// supported coordinates, for source-impulse and support-reaction consumers.
    pub midpoint_force_n: Vec<f64>,
    /// Actual midpoint force paired with the accepted displacement increment.
    pub midpoint_work_j: f64,
    /// Midpoint work minus potential work: signed time-discretization/rounding
    /// defect. It is zero to roundoff for affine gravity and never physical heat.
    pub quadrature_defect_j: f64,
}

#[derive(Clone, Debug)]
enum ReferenceDensity {
    PerElement(Arc<[f64]>),
    Nodal(Arc<[f64]>),
}

#[derive(Clone, Debug)]
pub(super) struct BoundGravity<'a> {
    source: StationaryGravity<'a>,
    density: ReferenceDensity,
}

impl<'a> BoundGravity<'a> {
    pub(super) fn new(
        problem: &HyperTetProblem<'a, FiniteSolidMaterial>,
        nodal: Option<&TetrahedralDensityBody>,
        source: StationaryGravity<'a>,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteSolidError> {
        checkpoint(cx)?;
        let density = match nodal {
            Some(field) => {
                ReferenceDensity::Nodal(admit_reference_density_field(problem, field, cx)?.into())
            }
            None => {
                let values: Vec<f64> = match problem.material {
                    HyperelasticLawField::Uniform(material) => {
                        vec![material.reference_density_kg_m3(); problem.tetrahedra.len()]
                    }
                    HyperelasticLawField::PerElement(materials) => materials
                        .iter()
                        .map(FiniteSolidMaterial::reference_density_kg_m3)
                        .collect(),
                };
                ReferenceDensity::PerElement(values.into())
            }
        };
        Ok(Self { source, density })
    }

    pub(super) fn evaluate(
        &self,
        problem: &PreparedHyperTetProblem<'a, FiniteSolidMaterial>,
        displacement: &[f64],
        cx: &Cx<'_>,
    ) -> Result<TetGravityLoad3, HyperTetError> {
        let field = crate::solid_gravity::bind_field(self.source.field, self.source.frame);
        match &self.density {
            ReferenceDensity::PerElement(values) => {
                problem.gravity_load(values, displacement, &field, cx)
            }
            ReferenceDensity::Nodal(values) => {
                problem.gravity_load_with_nodal_density(values, displacement, &field, cx)
            }
        }
    }

    fn force_and_potential(
        &self,
        problem: &PreparedHyperTetProblem<'a, FiniteSolidMaterial>,
        displacement: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, f64), HyperTetError> {
        let field = crate::solid_gravity::bind_field(self.source.field, self.source.frame);
        match &self.density {
            ReferenceDensity::PerElement(values) => {
                problem.gravity_force_and_potential(values, displacement, &field, cx)
            }
            ReferenceDensity::Nodal(values) => problem
                .gravity_force_and_potential_with_nodal_density(values, displacement, &field, cx),
        }
    }
}

impl<'a> FiniteSolidBody<'a> {
    /// Initialize the original finite-solid motion in one stationary gravity field.
    ///
    /// Every later ordinary, pressure or energy/momentum candidate automatically
    /// resamples this same field on its current material configuration. Source
    /// density, support reactions and initial acceleration use the original mass.
    /// The native midpoint method is exact in its gravity work for affine fields;
    /// nonlinear gravity retains its measured potential-work defect and timestep
    /// refinement requirement. A field-domain failure publishes no candidate.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_gravity(
        problem: HyperTetProblem<'a, FiniteSolidMaterial>,
        gravity: StationaryGravity<'a>,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        time_s: f64,
        controls: FiniteSolidInitialization,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteSolidError> {
        Self::new_with_gravity_and_sources(
            problem,
            None,
            None,
            gravity,
            displacement_m,
            velocity_m_s,
            time_s,
            controls,
            cx,
        )
    }

    /// Combine the same field with an explicit original P1 density and/or pressure.
    ///
    /// Optional reference matter has exactly the existing density-field admission
    /// and replaces the cellwise inertial density. An optional exterior and its
    /// initial pressure have the existing current-area pressure semantics; every
    /// later pressure interval must still be declared through a pressure entry.
    /// Neither pressure nor gravity is absorbed into a frozen body-force density.
    #[allow(clippy::too_many_arguments)]
    pub fn new_with_gravity_and_sources<'p>(
        problem: HyperTetProblem<'a, FiniteSolidMaterial>,
        reference_density_field: Option<Arc<TetrahedralDensityBody>>,
        pressure: Option<(Arc<ClosedTetPressureSurface>, TetPressure<'p>)>,
        gravity: StationaryGravity<'a>,
        displacement_m: &[[f64; 3]],
        velocity_m_s: &[[f64; 3]],
        time_s: f64,
        controls: FiniteSolidInitialization,
        cx: &Cx<'_>,
    ) -> Result<Self, FiniteSolidError> {
        let (surface, pressure) = match pressure {
            Some((surface, pressure)) => (Some(surface), pressure),
            None => (None, TetPressure::Uniform(0.0)),
        };
        Self::initialize(
            problem,
            reference_density_field,
            Some(gravity),
            surface,
            pressure,
            displacement_m,
            velocity_m_s,
            time_s,
            controls,
            cx,
        )
    }

    /// Retained source and inertial frame, shared unchanged by impulse candidates.
    pub fn stationary_gravity(&self) -> Option<StationaryGravity<'a>> {
        self.gravity.as_ref().map(|gravity| gravity.source)
    }
}

impl DynamicOperator<'_, '_, '_> {
    pub(super) fn ensure_gravity(&self, displacement: &[f64]) -> Result<(), HyperTetError> {
        let Some(gravity) = self.gravity else {
            return Ok(());
        };
        let rebuild = self
            .gravity_at
            .borrow()
            .as_ref()
            .is_none_or(|(at, _)| at.as_slice() != displacement);
        if rebuild {
            let load = gravity.evaluate(self.problem, displacement, self.cx)?;
            self.gravity_at.replace(Some((displacement.to_vec(), load)));
        }
        Ok(())
    }

    pub(super) fn subtract_gravity_force(
        &self,
        full: &mut [f64],
        displacement: &[f64],
    ) -> Result<(), HyperTetError> {
        self.ensure_gravity(displacement)?;
        if let Some((_, load)) = self.gravity_at.borrow().as_ref() {
            for (value, gravity) in full.iter_mut().zip(&load.force_n) {
                *value -= gravity;
            }
        }
        Ok(())
    }

    pub(super) fn retain_gravity_force(&self, displacement: &[f64]) -> Result<(), HyperTetError> {
        if self.gravity.is_none() {
            return Ok(());
        }
        let cached = self.gravity_at.borrow();
        let Some((at, load)) = cached.as_ref() else {
            return Err(invalid("native gravity force stage was not evaluated"));
        };
        if at.as_slice() != displacement || load.force_n.len() != displacement.len() {
            return Err(invalid(
                "native gravity force does not match its evaluated stage",
            ));
        }
        self.gravity_force_at
            .replace(Some((at.clone(), load.force_n.clone())));
        Ok(())
    }

    pub(super) fn combine_gravity_tangent(
        &self,
        full: Csr,
        displacement: &[f64],
    ) -> Result<Csr, HyperTetError> {
        if self.gravity.is_none() {
            return Ok(full);
        }
        self.ensure_gravity(displacement)?;
        let gravity = self.gravity_at.borrow();
        let load = &gravity.as_ref().expect("admitted gravity load").1;
        // The discrete material/pressure gradient is differentiated with
        // respect to the endpoint. Its native midpoint reconstruction supplies
        // the factor two. Gravity is evaluated directly at the midpoint, so
        // its analytic tangent must NOT receive that endpoint chain factor.
        let material_scale = if self.energy_momentum_start.is_some() {
            2.0
        } else {
            1.0
        };
        let mut combined = Coo::new(full.nrows(), full.ncols());
        for row in 0..full.nrows() {
            checkpoint(self.cx)?;
            for (matrix, scale) in [(&full, material_scale), (&load.force_tangent_n_m, -1.0)] {
                let (columns, values) = matrix.row(row);
                for (&column, &value) in columns.iter().zip(values) {
                    combined.push(row, column, scale * value);
                }
            }
        }
        let combined = combined.assemble();
        for row in 0..combined.nrows() {
            checkpoint(self.cx)?;
            if combined.row(row).1.iter().any(|value| !value.is_finite()) {
                return Err(invalid(
                    "native solid minus gravity tangent is not representable",
                ));
            }
        }
        Ok(combined)
    }

    pub(super) fn gravity_balance(
        &self,
        before: &[f64],
        after: &[f64],
    ) -> Result<Option<GravityStep>, HyperTetError> {
        let Some(gravity) = self.gravity else {
            return Ok(None);
        };
        // Native Newton retains the residual of its accepted point. Its last
        // force callback therefore belongs to the converged stage, including
        // zero-iteration convergence. Trust-region acceptance subsequently
        // applies an OLD-point Jacobian; endpoint admission also samples another
        // state. Neither is allowed to overwrite this force-only snapshot.
        let cached = self.gravity_force_at.borrow();
        let Some((stage, force_n)) = cached.as_ref() else {
            return Err(invalid("accepted native gravity force stage is missing"));
        };
        if before.len() != after.len()
            || stage.len() != before.len()
            || force_n.len() != before.len()
        {
            return Err(invalid(
                "accepted native gravity stage shape is inconsistent",
            ));
        }
        for (dof, ((&old, &new), &at)) in before.iter().zip(after).zip(stage).enumerate() {
            if dof % 4096 == 0 {
                checkpoint(self.cx)?;
            }
            // fs-time forms q_stage = q_before + 0.5*dq with FMA, before the
            // accepted endpoint rounds q_before+dq. Check its relation to the
            // published endpoints without replacing it by a reconstructed point.
            let midpoint = 0.5 * old + 0.5 * new;
            let scale = old
                .abs()
                .max(new.abs())
                .max(at.abs())
                .max(f64::MIN_POSITIVE);
            if !at.is_finite() || (at - midpoint).abs() > 8.0 * f64::EPSILON * scale {
                return Err(invalid(
                    "gravity force stage does not match the accepted interval",
                ));
            }
        }
        let first = gravity.evaluate(self.problem, before, self.cx)?;
        let last = gravity.evaluate(self.problem, after, self.cx)?;
        let work = dd(first.potential_energy_j) - dd(last.potential_energy_j);
        let mut midpoint_work = dd(0.0);
        for (dof, ((&old, &new), &force)) in before.iter().zip(after).zip(force_n).enumerate() {
            if dof % 4096 == 0 {
                checkpoint(self.cx)?;
            }
            midpoint_work = midpoint_work + dd(force) * (dd(new) - dd(old));
        }
        let defect = midpoint_work - work;
        if [work.to_f64(), midpoint_work.to_f64(), defect.to_f64()]
            .iter()
            .any(|value| !value.is_finite())
        {
            return Err(invalid("gravity work is not representable"));
        }
        Ok(Some(GravityStep {
            frame: gravity.source.frame,
            potential_before_j: first.potential_energy_j,
            potential_after_j: last.potential_energy_j,
            work_j: work.to_f64(),
            midpoint_force_n: force_n.clone(),
            midpoint_work_j: midpoint_work.to_f64(),
            quadrature_defect_j: defect.to_f64(),
        }))
    }
}

fn dd(value: f64) -> Dd {
    Dd::from_f64(value)
}
fn invalid(what: &str) -> HyperTetError {
    HyperTetError::InvalidInput { what: what.into() }
}
