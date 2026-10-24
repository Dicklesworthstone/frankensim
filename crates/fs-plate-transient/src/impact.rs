//! Reciprocal contact between a finite normal mass and a full spatial plate.
//!
//! Both coordinates use the same positive transverse direction. Relative closure
//! is z - b.q, where b is the existing finite-area plate interpolation. The
//! supplied scalar contact law acts on that closure; its original quadrature
//! weights, gaps, stiffness, exponent and Hunt--Crossley coefficient are retained.
//! The footprint distributes the total reaction, so its physical area is not
//! multiplied into the constitutive law a second time.
//!
//! Every residual and tangent uses the existing fs-dcontact potential and loss
//! port at the same implicit displacement/velocity. The original generalized-
//! alpha/Newton--Krylov owner advances all plate DOFs and the finite mass together.
//! This is a fixed-footprint, normal-translation model, not contact detection,
//! moving-patch geometry, local pressure resolution, or a rotational rigid body.
//! Generalized-alpha has its own numerical dissipation and is not asserted to
//! preserve the nonlinear contact energy exactly.

use crate::{PlateDynamics, PlateLoad, PlateMassPreconditioner};
use fs_dcontact::{ContactStorage, DContactError, Obstacle};
use fs_phs::Storage;
use fs_plate::loading::{PlateLoadError, PlatePatchStencil, PlatePointSample};
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::second_order_adjoint::SecondOrderVjp;
use fs_time::galpha::second_order_adjoint::trajectory::StructuralTrajectoryModel;
use fs_time::galpha::{SecondOrderOperatorWeights, SecondOrderProblem};

/// Admission bounds for the combined spatial model and original contact law.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateImpactBudget {
    /// Combined plate plus indenter DOFs, including the final indenter coordinate.
    pub max_free_dofs: usize,
    /// Maximum original scalar contact-law quadrature points.
    pub max_contact_points: usize,
}

/// A malformed model, exceeded cap, owner refusal, or unrepresentable arithmetic.
#[derive(Clone, Debug, PartialEq)]
pub enum PlateImpactError {
    /// A physical coefficient or state violated the adapter's contract.
    Invalid(&'static str),
    /// A declared size bound or checked layout was exceeded.
    Budget(&'static str),
    /// Scratch allocation was refused before model or state publication.
    Allocation,
    /// The original normal-law owner refused its input.
    Contact(DContactError),
    /// The original finite-patch owner refused its input.
    Patch(PlateLoadError),
}
impl core::fmt::Display for PlateImpactError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        write!(f, "spatial plate impact refused: {self:?}")
    }
}
impl std::error::Error for PlateImpactError {}
impl From<DContactError> for PlateImpactError {
    fn from(error: DContactError) -> Self {
        Self::Contact(error)
    }
}
impl From<PlateLoadError> for PlateImpactError {
    fn from(error: PlateLoadError) -> Self {
        Self::Patch(error)
    }
}

fn zeros(n: usize) -> Result<Vec<f64>, PlateImpactError> {
    let mut values = Vec::new();
    values
        .try_reserve_exact(n)
        .map_err(|_| PlateImpactError::Allocation)?;
    values.resize(n, 0.0);
    Ok(values)
}
fn finite(values: &[f64]) -> bool {
    values.iter().all(|v| v.is_finite())
}
fn dot(left: &[f64], right: &[f64]) -> f64 {
    left.iter().zip(right).map(|(a, b)| a * b).sum()
}

struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self, _: &[f64]) -> f64 {
        0.0
    }
    fn gradient(&self, _: &[f64], output: &mut [f64]) {
        output.fill(0.0);
    }
}

/// Physical observations of one state, not an integrated energy-balance receipt.
#[derive(Clone, Debug, PartialEq)]
pub struct PlateImpactProbe {
    /// Original scaled spatial plate kinetic energy [J].
    pub plate_kinetic_energy_j: f64,
    /// Original scaled spatial plate elastic energy [J].
    pub plate_elastic_energy_j: f64,
    /// Finite indenter kinetic energy [J].
    pub indenter_kinetic_energy_j: f64,
    /// Original fs-dcontact potential, counted exactly once [J].
    pub contact_energy_j: f64,
    /// Sum of the four represented storage terms [J].
    pub total_energy_j: f64,
    /// Indenter displacement minus area-averaged plate displacement [m].
    pub closure_m: f64,
    /// The same signed relative material velocity [m/s].
    pub relative_velocity_m_per_s: f64,
    /// Maximum positive penetration among the original law points [m].
    pub maximum_penetration_m: f64,
    /// Total compressive reaction, including the non-attractive loss law [N].
    pub normal_force_n: f64,
    /// Original Hunt--Crossley irreversible power [W]; numerical loss is excluded.
    pub contact_dissipated_power_w: f64,
    /// Original plate Rayleigh power v.C.v [W]. Nonnegativity requires the
    /// original plate's stated stiffness/PSD assumptions.
    pub plate_dissipated_power_w: f64,
    /// Area-averaged plate displacement, slopes and material velocity.
    pub patch: PlatePointSample,
    /// Physical indenter displacement in the positive transverse direction [m].
    pub indenter_displacement_m: f64,
    /// Physical indenter velocity [m/s].
    pub indenter_velocity_m_per_s: f64,
    /// Original contact-law quadrature points with positive penetration.
    pub active_contact_points: usize,
}

/// Original spatial plate plus one finite-mass, normal-translation coordinate.
///
/// The final component of every q/v/a vector belongs to the indenter. All
/// contact points use the same patch-average closure with their own original
/// gaps and weights: this is a declared uniform-pressure interface reduction,
/// not a resolved distribution of local penetration. The indenter mass, fixed
/// footprint, and contact-law parameters are held fixed by the adjoint; its
/// parameter order is exactly the original plate's four parameters plus loads.
pub struct PlateImpact<'a, L: ?Sized> {
    plate: PlateDynamics<'a, L>,
    patch: PlatePatchStencil,
    indenter_mass_kg: f64,
    contact: ContactStorage,
    closure: Vec<f64>,
    dimension: usize,
}

#[derive(Clone, Copy)]
enum ContactPart {
    Elastic,
    Total,
    ElasticTangent,
    DisplacementTangent,
    VelocityTangent,
}

impl<'a, L: PlateLoad + ?Sized> PlateImpact<'a, L> {
    /// Bind a fixed footprint and an original scalar positive-closure law.
    ///
    /// Every original collocation row must be [1]. Raw-parts obstacles are
    /// rebuilt through Obstacle::new/with_internal_loss before any evaluation;
    /// no physical coefficient is inferred or rescaled. The patch support map
    /// must match the supplied plate. Both retain the caller's unchanged mesh.
    ///
    /// Work is bounded by the DOF and contact-point caps; contact callbacks use
    /// O(DOFs) fallible scratch and O(DOFs * contact points) arithmetic. They
    /// allocate and are not a hard-real-time API. The existing integrator owns
    /// cancellation between bounded callbacks and transactional time-step state.
    pub fn new(
        plate: PlateDynamics<'a, L>,
        patch: PlatePatchStencil,
        indenter_mass_kg: f64,
        law: Obstacle,
        budget: PlateImpactBudget,
    ) -> Result<Self, PlateImpactError> {
        let dimension = plate
            .dimension()
            .checked_add(1)
            .ok_or(PlateImpactError::Budget("combined dimension overflow"))?;
        if dimension > budget.max_free_dofs
            || law.n_points() > budget.max_contact_points
            || dimension.checked_mul(2).is_none()
            || dimension.checked_mul(law.n_points()).is_none()
        {
            return Err(PlateImpactError::Budget("combined DOF/contact-point cap"));
        }
        if !indenter_mass_kg.is_finite()
            || indenter_mass_kg <= 0.0
            || !indenter_mass_kg.recip().is_finite()
        {
            return Err(PlateImpactError::Invalid(
                "positive finite invertible indenter mass required",
            ));
        }
        if law.n_points() == 0
            || law.collocation().len() != law.n_points()
            || law.gaps().len() != law.n_points()
            || law.weights().len() != law.n_points()
            || law.collocation().iter().any(|b| *b != 1.0)
        {
            return Err(PlateImpactError::Invalid(
                "scalar positive-closure rows [1] required",
            ));
        }
        patch.validate_model(plate.model())?;
        let admitted = Obstacle::new(
            law.collocation().to_vec(),
            law.n_points(),
            1,
            law.gaps().to_vec(),
            law.weights().to_vec(),
            law.stiffness(),
            law.alpha(),
            law.provenance().to_owned(),
        )?
        .with_internal_loss(law.internal_loss())?;
        let n = dimension - 1;
        let mut closure = zeros(dimension)?;
        patch.add_load([1.0, 0.0, 0.0], &mut closure[..n])?;
        for b in &mut closure[..n] {
            *b = -*b;
        }
        closure[n] = 1.0;
        let mut collocation = zeros(dimension * admitted.n_points())?;
        for row in collocation.chunks_exact_mut(dimension) {
            row.copy_from_slice(&closure);
        }
        let lifted = Obstacle::new(
            collocation,
            admitted.n_points(),
            dimension,
            admitted.gaps().to_vec(),
            admitted.weights().to_vec(),
            admitted.stiffness(),
            admitted.alpha(),
            admitted.provenance().to_owned(),
        )?
        .with_internal_loss(admitted.internal_loss())?;
        let contact = ContactStorage::new(Box::new(ZeroStorage), dimension, vec![lifted])?;
        Ok(Self {
            plate,
            patch,
            indenter_mass_kg,
            contact,
            closure,
            dimension,
        })
    }

    /// The original spatial model, physical parameters and external load owner.
    #[must_use]
    pub fn plate(&self) -> &PlateDynamics<'a, L> {
        &self.plate
    }
    /// Number of plate DOFs preceding the indenter coordinate.
    #[must_use]
    pub fn plate_dofs(&self) -> usize {
        self.dimension - 1
    }
    /// Combined displacement/velocity/acceleration dimension.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.dimension
    }
    /// Original finite normal mass [kg].
    #[must_use]
    pub fn indenter_mass_kg(&self) -> f64 {
        self.indenter_mass_kg
    }

    /// Exact block-diagonal inverse mass for consistent initialization and its
    /// transposed solve. This does not approximate the nonlinear effective Jacobian.
    #[must_use]
    pub fn mass_preconditioner(&self) -> PlateImpactMassPreconditioner<'_> {
        PlateImpactMassPreconditioner {
            plate: self.plate.mass_preconditioner(),
            plate_dofs: self.plate_dofs(),
            indenter_mass_kg: self.indenter_mass_kg,
        }
    }

    fn check(&self, fields: &[&[f64]]) -> Result<(), PlateImpactError> {
        if fields
            .iter()
            .any(|x| x.len() != self.dimension || !finite(x))
        {
            Err(PlateImpactError::Invalid(
                "finite combined-dimensional states required",
            ))
        } else {
            Ok(())
        }
    }
    fn io(&self, fields: &[&[f64]], output: &mut [f64]) -> bool {
        if output.len() != self.dimension || self.check(fields).is_err() {
            output.fill(f64::NAN);
            false
        } else {
            true
        }
    }
    fn pack(&self, q: &[f64]) -> Result<Vec<f64>, PlateImpactError> {
        let mut x = zeros(2 * self.dimension)?;
        for (pair, value) in x.chunks_exact_mut(2).zip(q) {
            pair[0] = *value;
        }
        Ok(x)
    }
    fn contact_action(
        &self,
        q: &[f64],
        v: Option<&[f64]>,
        direction: Option<&[f64]>,
        part: ContactPart,
    ) -> Result<Vec<f64>, PlateImpactError> {
        self.check(&[q])?;
        if let Some(v) = v {
            self.check(&[v])?;
        }
        if let Some(d) = direction {
            self.check(&[d])?;
        }
        let x = self.pack(q)?;
        let mut conservative = zeros(x.len())?;
        let mut dissipative = zeros(x.len())?;
        match part {
            ContactPart::Elastic | ContactPart::Total => {
                self.contact.gradient(&x, &mut conservative)
            }
            ContactPart::ElasticTangent | ContactPart::DisplacementTangent => {
                let dx =
                    self.pack(direction.ok_or(PlateImpactError::Invalid("missing direction"))?)?;
                if !self
                    .contact
                    .hessian_vector_with(&x, &dx, &mut conservative, |_, _, out| {
                        out.fill(0.0);
                        true
                    })
                {
                    return Err(PlateImpactError::Invalid("contact Hessian action"));
                }
            }
            ContactPart::VelocityTangent => {}
        }
        if matches!(
            part,
            ContactPart::Total | ContactPart::DisplacementTangent | ContactPart::VelocityTangent
        ) {
            let velocity = v.ok_or(PlateImpactError::Invalid("missing contact velocity"))?;
            let mut e = zeros(x.len())?;
            for (pair, value) in e.chunks_exact_mut(2).zip(velocity) {
                pair[1] = *value;
            }
            let valid = if matches!(part, ContactPart::Total) {
                self.contact.dissipative_flow_into(&x, &e, &mut dissipative)
            } else {
                let d = direction.ok_or(PlateImpactError::Invalid("missing direction"))?;
                let mut dx = zeros(x.len())?;
                let mut de = zeros(x.len())?;
                for (i, value) in d.iter().enumerate() {
                    if matches!(part, ContactPart::DisplacementTangent) {
                        dx[2 * i] = *value;
                    } else {
                        de[2 * i + 1] = *value;
                    }
                }
                self.contact
                    .dissipative_flow_tangent_into(&x, &e, &dx, &de, &mut dissipative)
            };
            if !valid {
                return Err(PlateImpactError::Invalid("contact loss action"));
            }
        }
        let mut out = zeros(self.dimension)?;
        for (i, value) in out.iter_mut().enumerate() {
            *value = conservative[2 * i] + dissipative[2 * i + 1];
        }
        if !finite(&out) {
            return Err(PlateImpactError::Invalid("nonfinite contact action"));
        }
        Ok(out)
    }

    fn add_contact(
        output: &mut [f64],
        contact: Result<Vec<f64>, PlateImpactError>,
    ) -> Result<(), String> {
        let values = contact.map_err(|e| e.to_string())?;
        for (out, value) in output.iter_mut().zip(values) {
            *out += value;
        }
        if !finite(output) {
            return Err("nonfinite combined plate/contact action".into());
        }
        Ok(())
    }

    /// Observe the actual physical state without changing history or solver state.
    pub fn probe(&self, q: &[f64], v: &[f64]) -> Result<PlateImpactProbe, PlateImpactError> {
        self.check(&[q, v])?;
        let n = self.plate_dofs();
        let mut work = zeros(n)?;
        self.plate.mass_apply(&v[..n], &mut work);
        let plate_kinetic_energy_j = 0.5 * dot(&v[..n], &work);
        self.plate.internal_force(&q[..n], &mut work);
        let plate_elastic_energy_j = 0.5 * dot(&q[..n], &work);
        self.plate.damping_apply(&v[..n], &mut work);
        let plate_dissipated_power_w = dot(&v[..n], &work);
        let indenter_kinetic_energy_j = 0.5 * v[n] * (self.indenter_mass_kg * v[n]);
        let x = self.pack(q)?;
        let contact = self.contact.probe(&x);
        let reaction = self.contact_action(q, Some(v), None, ContactPart::Total)?;
        let mut loss_power = [0.0];
        if !self
            .contact
            .dissipated_power_by_obstacle_into(&x, v, &mut loss_power)
        {
            return Err(PlateImpactError::Invalid("contact dissipation observation"));
        }
        let closure_m = dot(&self.closure, q);
        let relative_velocity_m_per_s = dot(&self.closure, v);
        let total_energy_j = plate_kinetic_energy_j
            + plate_elastic_energy_j
            + indenter_kinetic_energy_j
            + contact.contact_energy;
        if !finite(&[
            plate_kinetic_energy_j,
            plate_elastic_energy_j,
            plate_dissipated_power_w,
            indenter_kinetic_energy_j,
            contact.contact_energy,
            contact.max_penetration,
            reaction[n],
            loss_power[0],
            closure_m,
            relative_velocity_m_per_s,
            total_energy_j,
        ]) {
            return Err(PlateImpactError::Invalid(
                "nonfinite physical impact observation",
            ));
        }
        Ok(PlateImpactProbe {
            plate_kinetic_energy_j,
            plate_elastic_energy_j,
            indenter_kinetic_energy_j,
            contact_energy_j: contact.contact_energy,
            total_energy_j,
            closure_m,
            relative_velocity_m_per_s,
            maximum_penetration_m: contact.max_penetration,
            normal_force_n: reaction[n],
            contact_dissipated_power_w: loss_power[0],
            plate_dissipated_power_w,
            patch: self.patch.sample(&q[..n], &v[..n])?,
            indenter_displacement_m: q[n],
            indenter_velocity_m_per_s: v[n],
            active_contact_points: contact.active_points,
        })
    }
}

impl<L: PlateLoad + ?Sized> SecondOrderProblem for PlateImpact<'_, L> {
    fn dimension(&self) -> usize {
        self.dimension
    }
    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        if !self.io(&[input], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate.mass_apply(&input[..n], &mut output[..n]);
        output[n] = self.indenter_mass_kg * input[n];
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }
    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        if !self.io(&[input], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate.damping_apply(&input[..n], &mut output[..n]);
        output[n] = 0.0;
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }
    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        if !self.io(&[q], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate.internal_force(&q[..n], &mut output[..n]);
        output[n] = 0.0;
        if Self::add_contact(
            output,
            self.contact_action(q, None, None, ContactPart::Elastic),
        )
        .is_err()
        {
            output.fill(f64::NAN);
        }
    }
    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        if !self.io(&[q, direction], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate
            .tangent_apply(&q[..n], &direction[..n], &mut output[..n]);
        output[n] = 0.0;
        if Self::add_contact(
            output,
            self.contact_action(q, None, Some(direction), ContactPart::ElasticTangent),
        )
        .is_err()
        {
            output.fill(f64::NAN);
        }
    }
    fn resistance(&self, q: &[f64], v: &[f64], output: &mut [f64]) {
        if !self.io(&[q, v], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate.resistance(&q[..n], &v[..n], &mut output[..n]);
        output[n] = 0.0;
        if Self::add_contact(
            output,
            self.contact_action(q, Some(v), None, ContactPart::Total),
        )
        .is_err()
        {
            output.fill(f64::NAN);
        }
    }
    fn resistance_displacement_apply(
        &self,
        q: &[f64],
        v: &[f64],
        direction: &[f64],
        output: &mut [f64],
    ) {
        if !self.io(&[q, v, direction], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate.resistance_displacement_apply(
            &q[..n],
            &v[..n],
            &direction[..n],
            &mut output[..n],
        );
        output[n] = 0.0;
        if Self::add_contact(
            output,
            self.contact_action(
                q,
                Some(v),
                Some(direction),
                ContactPart::DisplacementTangent,
            ),
        )
        .is_err()
        {
            output.fill(f64::NAN);
        }
    }
    fn resistance_velocity_apply(
        &self,
        q: &[f64],
        v: &[f64],
        direction: &[f64],
        output: &mut [f64],
    ) {
        if !self.io(&[q, v, direction], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate
            .resistance_velocity_apply(&q[..n], &v[..n], &direction[..n], &mut output[..n]);
        output[n] = 0.0;
        if Self::add_contact(
            output,
            self.contact_action(q, Some(v), Some(direction), ContactPart::VelocityTangent),
        )
        .is_err()
        {
            output.fill(f64::NAN);
        }
    }
    fn preconditioner_apply(
        &self,
        q: &[f64],
        weights: SecondOrderOperatorWeights,
        outer: usize,
        inner: usize,
        residual: &[f64],
        output: &mut [f64],
    ) {
        if !self.io(&[q, residual], output) {
            return;
        }
        let n = self.plate_dofs();
        self.plate.preconditioner_apply(
            &q[..n],
            weights,
            outer,
            inner,
            &residual[..n],
            &mut output[..n],
        );
        let diagonal = weights.mass * self.indenter_mass_kg;
        if !diagonal.is_finite() || diagonal <= 0.0 || !diagonal.recip().is_finite() {
            output.fill(f64::NAN);
            return;
        }
        // A block Jacobi approximation only: the actual Newton operator retains
        // all state-dependent contact displacement and velocity derivatives.
        output[n] = residual[n] / diagonal;
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }
}

impl<L: PlateLoad + ?Sized> SecondOrderVjp for PlateImpact<'_, L> {
    fn parameter_count(&self) -> usize {
        self.plate.parameter_count()
    }
    fn mass_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        if !self.io(&[seed], output) {
            return Err("invalid impact mass transpose".into());
        }
        let n = self.plate_dofs();
        self.plate
            .mass_transpose_apply(&seed[..n], &mut output[..n])?;
        output[n] = self.indenter_mass_kg * seed[n];
        if !finite(output) {
            return Err("nonfinite impact mass transpose".into());
        }
        Ok(())
    }
    fn damping_transpose_apply(&self, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        if !self.io(&[seed], output) {
            return Err("invalid impact damping transpose".into());
        }
        let n = self.plate_dofs();
        self.plate
            .damping_transpose_apply(&seed[..n], &mut output[..n])?;
        output[n] = 0.0;
        Ok(())
    }
    fn tangent_transpose_apply(
        &self,
        q: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        if !self.io(&[q, seed], output) {
            return Err("invalid impact elastic transpose".into());
        }
        let n = self.plate_dofs();
        self.plate
            .tangent_transpose_apply(&q[..n], &seed[..n], &mut output[..n])?;
        output[n] = 0.0;
        Self::add_contact(
            output,
            self.contact_action(q, None, Some(seed), ContactPart::ElasticTangent),
        )
    }
    fn resistance_displacement_transpose_apply(
        &self,
        q: &[f64],
        v: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        if !self.io(&[q, v, seed], output) {
            return Err("invalid impact displacement transpose".into());
        }
        let n = self.plate_dofs();
        self.plate.resistance_displacement_transpose_apply(
            &q[..n],
            &v[..n],
            &seed[..n],
            &mut output[..n],
        )?;
        output[n] = 0.0;
        // Each contact contribution is a scalar times B^T B, including the
        // unsaturated or capped loss derivative. Only this local part is symmetric.
        Self::add_contact(
            output,
            self.contact_action(q, Some(v), Some(seed), ContactPart::DisplacementTangent),
        )
    }
    fn resistance_velocity_transpose_apply(
        &self,
        q: &[f64],
        v: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        if !self.io(&[q, v, seed], output) {
            return Err("invalid impact velocity transpose".into());
        }
        let n = self.plate_dofs();
        self.plate.resistance_velocity_transpose_apply(
            &q[..n],
            &v[..n],
            &seed[..n],
            &mut output[..n],
        )?;
        output[n] = 0.0;
        Self::add_contact(
            output,
            self.contact_action(q, Some(v), Some(seed), ContactPart::VelocityTangent),
        )
    }
    fn residual_parameter_vjp(
        &self,
        q: &[f64],
        v: &[f64],
        a: &[f64],
        seed: &[f64],
        output: &mut [f64],
    ) -> Result<(), String> {
        self.check(&[q, v, a, seed]).map_err(|e| e.to_string())?;
        let n = self.plate_dofs();
        self.plate
            .residual_parameter_vjp(&q[..n], &v[..n], &a[..n], &seed[..n], output)
    }
}

impl<L: PlateLoad + ?Sized> StructuralTrajectoryModel for PlateImpact<'_, L> {
    fn forcing(&self, time: f64, output: &mut [f64]) -> Result<(), String> {
        if output.len() != self.dimension {
            return Err("invalid impact forcing dimension".into());
        }
        let n = self.plate_dofs();
        self.plate.forcing(time, &mut output[..n])?;
        output[n] = 0.0;
        Ok(())
    }
    fn forcing_vjp(&self, time: f64, seed: &[f64], output: &mut [f64]) -> Result<(), String> {
        self.check(&[seed]).map_err(|e| e.to_string())?;
        self.plate
            .forcing_vjp(time, &seed[..self.plate_dofs()], output)
    }
}

/// Exact plate/indenter mass inverse, with the same fixed physical masses as
/// PlateImpact. Suitable for the shared consistent-acceleration forward and
/// transpose solves; the nonlinear effective operator has a separate hook.
pub struct PlateImpactMassPreconditioner<'a> {
    plate: PlateMassPreconditioner<'a>,
    plate_dofs: usize,
    indenter_mass_kg: f64,
}
impl FlexiblePreconditioner for PlateImpactMassPreconditioner<'_> {
    fn apply(&self, iteration: usize, residual: &[f64], output: &mut [f64]) {
        let n = self.plate_dofs;
        if residual.len() != n + 1 || output.len() != n + 1 || !finite(residual) {
            output.fill(f64::NAN);
            return;
        }
        self.plate
            .apply(iteration, &residual[..n], &mut output[..n]);
        output[n] = residual[n] / self.indenter_mass_kg;
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }
}
