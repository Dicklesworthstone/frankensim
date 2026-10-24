//! Simultaneous finite bodies on one shared spatial plate and physical clock.
//!
//! State order is `[plate DOFs, body_0, body_1, ...]`. Each original contact
//! sees the same trial plate motion and its own normal body coordinate. The
//! plate operators act once; only contact reactions are gathered and scattered.
//! Thus the original implicit solve retains the complete interaction through
//! the spatial plate without another contact solver or duplicated plate storage.

use super::{
    ContactPart, PlateDynamics, PlateImpact, PlateImpactBodyProbe, PlateImpactError, PlateLoad,
    dot, finite, zeros,
};
use core::ops::Range;
use fs_solver::FlexiblePreconditioner;
use fs_time::galpha::second_order_adjoint::SecondOrderVjp;
use fs_time::galpha::{SecondOrderOperatorWeights, SecondOrderProblem};

mod adjoint;

/// Aggregate admission ceilings, in addition to each original impact's bounds.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PlateImpactSetBudget {
    /// Maximum simultaneous finite bodies; at least one must be supplied.
    pub max_bodies: usize,
    /// Maximum plate plus body displacement coordinates.
    pub max_free_dofs: usize,
    /// Maximum sum of original contact-law quadrature points over all bodies.
    pub max_contact_points: usize,
    /// Maximum shared plate/load, body calibration and body-force parameters.
    pub max_parameters: usize,
}

/// Physical observation of one joint state, with the shared plate counted once.
/// These are instantaneous storage/power observations, not a step energy audit.
#[derive(Clone, Debug, PartialEq)]
pub struct PlateImpactSetProbe {
    /// Shared spatial plate kinetic energy [J].
    pub plate_kinetic_energy_j: f64,
    /// Shared spatial plate elastic energy [J].
    pub plate_elastic_energy_j: f64,
    /// Shared original Rayleigh damping power [W].
    pub plate_dissipated_power_w: f64,
    /// Sum of finite-body kinetic energies [J].
    pub indenter_kinetic_energy_j: f64,
    /// Sum of original contact potentials [J].
    pub contact_energy_j: f64,
    /// Sum of original contact irreversible powers [W].
    pub contact_dissipated_power_w: f64,
    /// Shared plate plus every body and contact storage term [J].
    pub total_energy_j: f64,
    /// Per-body observations in the original constructor order.
    pub bodies: Vec<PlateImpactBodyProbe>,
}

/// One spatial plate with independent finite masses at fixed normal footprints.
///
/// Prepare each original impact from a copy of the SAME [`PlateDynamics`].
/// Admission checks model and load identity, physical parameters and the cached
/// load parameter count. Equal-sized different plates or loads cannot be mixed.
/// All original authored/source owners and source limits remain attached.
///
/// The common plate/load parameter prefix appears once. Each body's existing
/// optional calibration suffix follows in body order; source-bound and fixed
/// authored bodies append no parameters. Optional applied body-force controls
/// append after those suffixes. Contact geometry, finite footprints,
/// and body identities are fixed. This adds neither collision detection nor
/// tangential/rotational body motion or body-to-body contacts.
pub struct PlateImpactSet<'a, L: ?Sized> {
    bodies: Vec<PlateImpact<'a, L>>,
    plate_dofs: usize,
    dimension: usize,
    parameters: usize,
    parameter_ranges: Vec<Range<usize>>,
    parameter_cap: usize,
    body_forcing: Option<(&'a dyn PlateLoad, usize)>,
}

impl<'a, L: PlateLoad + ?Sized> PlateImpactSet<'a, L> {
    /// Admit a nonempty set before allocating any aggregate scratch or layout.
    pub fn new(
        bodies: Vec<PlateImpact<'a, L>>,
        budget: PlateImpactSetBudget,
    ) -> Result<Self, PlateImpactError> {
        if bodies.is_empty() {
            return Err(PlateImpactError::Invalid(
                "at least one impact body required",
            ));
        }
        if bodies.len() > budget.max_bodies {
            return Err(PlateImpactError::Budget("simultaneous impact body count"));
        }
        let plate = bodies[0].plate();
        let plate_dofs = plate.dimension();
        let dimension = plate_dofs
            .checked_add(bodies.len())
            .filter(|n| *n <= budget.max_free_dofs)
            .ok_or(PlateImpactError::Budget("combined plate/body dimension"))?;
        let mut points = 0usize;
        let mut parameters = plate.parameter_count();
        for body in &bodies {
            if !plate.same_model_and_load(body.plate()) {
                return Err(PlateImpactError::Invalid(
                    "impact bodies must share the same plate, load and physical parameters",
                ));
            }
            for law in body.contact.obstacles() {
                points = points
                    .checked_add(law.n_points())
                    .ok_or(PlateImpactError::Budget(
                        "aggregate contact quadrature point count",
                    ))?;
            }
            parameters = parameters
                .checked_add(body.contact_parameters)
                .ok_or(PlateImpactError::Budget("aggregate impact parameter count"))?;
        }
        if points > budget.max_contact_points || parameters > budget.max_parameters {
            return Err(PlateImpactError::Budget(
                "aggregate impact contact points or parameters",
            ));
        }
        let mut parameter_ranges = Vec::new();
        parameter_ranges
            .try_reserve_exact(bodies.len())
            .map_err(|_| PlateImpactError::Allocation)?;
        let mut cursor = plate.parameter_count();
        for body in &bodies {
            let end = cursor + body.contact_parameters;
            parameter_ranges.push(cursor..end);
            cursor = end;
        }
        Ok(Self {
            bodies,
            plate_dofs,
            dimension,
            parameters,
            parameter_ranges,
            parameter_cap: budget.max_parameters,
            body_forcing: None,
        })
    }

    /// Bind prescribed normal forces on the ordered finite-body coordinates.
    ///
    /// Reuses the existing pure time-dependent [`PlateLoad`] protocol with an
    /// output slice of length body_count(), in constructor order. These are
    /// external forces in newtons; state-dependent contact remains in the
    /// coupled residual. Positive force follows the positive body coordinate.
    ///
    /// Force-history parameters append after all existing plate/contact
    /// coordinates, including when contact coefficients are source-owned. Their
    /// residual partials are zero; the trajectory forcing pullback owns their
    /// contribution. The original set parameter cap remains binding. Rebinding
    /// replaces this history and its parameter suffix, preserving every prior
    /// contact range. Rebuild at each new immutable physical parameter point.
    pub fn with_body_forcing(mut self, load: &'a dyn PlateLoad) -> Result<Self, PlateImpactError> {
        let old_count = self.body_forcing.map_or(0, |(_, count)| count);
        let fixed = self.parameters - old_count;
        let count = load.parameter_count();
        let parameters = fixed
            .checked_add(count)
            .filter(|n| *n <= self.parameter_cap)
            .ok_or(PlateImpactError::Budget("body-force parameter count"))?;
        self.body_forcing = Some((load, count));
        self.parameters = parameters;
        Ok(self)
    }

    /// Ordered body-force history coordinates in the complete parameter vector.
    /// None denotes no bound history; a prescribed parameter-free history has
    /// a present empty range.
    #[must_use]
    pub fn body_forcing_parameter_range(&self) -> Option<Range<usize>> {
        self.body_forcing
            .map(|(_, count)| self.parameters - count..self.parameters)
    }

    fn checked_body_forcing(&self) -> Result<Option<(&dyn PlateLoad, Range<usize>)>, String> {
        let Some((load, count)) = self.body_forcing else {
            return Ok(None);
        };
        if load.parameter_count() != count {
            return Err("body-force history parameter count changed".into());
        }
        Ok(Some((load, self.parameters - count..self.parameters)))
    }

    /// Number of original finite bodies.
    #[must_use]
    pub fn body_count(&self) -> usize {
        self.bodies.len()
    }

    /// Number of shared reduced plate displacement/slope coordinates.
    #[must_use]
    pub fn plate_dofs(&self) -> usize {
        self.plate_dofs
    }

    /// Shared plate coordinates plus one normal coordinate per body.
    #[must_use]
    pub fn dimension(&self) -> usize {
        self.dimension
    }

    /// Original shared spatial matrices, physical parameter point and load.
    #[must_use]
    pub fn plate(&self) -> &PlateDynamics<'a, L> {
        self.bodies[0].plate()
    }

    /// Original body/contact owner in constructor order. Its single-body APIs
    /// use local vectors `[plate DOFs, this body]`, not the set's full layout.
    #[must_use]
    pub fn body(&self, index: usize) -> Option<&PlateImpact<'a, L>> {
        self.bodies.get(index)
    }

    /// Optional calibration suffix in the full parameter vector. A present but
    /// empty range means this body's physical parameters are fixed; None means
    /// the body index is outside the admitted set.
    #[must_use]
    pub fn contact_parameter_range(&self, body: usize) -> Option<Range<usize>> {
        self.parameter_ranges.get(body).cloned()
    }

    /// Exact original plate and finite-body mass inverse for initialization and
    /// transposed mass solves. The effective nonlinear solve has its own hook.
    #[must_use]
    pub fn mass_preconditioner(&self) -> PlateImpactSetMassPreconditioner<'_, 'a, L> {
        PlateImpactSetMassPreconditioner { system: self }
    }

    fn check(&self, fields: &[&[f64]]) -> Result<(), PlateImpactError> {
        if fields
            .iter()
            .any(|field| field.len() != self.dimension || !finite(field))
        {
            return Err(PlateImpactError::Invalid(
                "finite complete impact-set states required",
            ));
        }
        Ok(())
    }

    fn io(&self, fields: &[&[f64]], output: &mut [f64]) -> bool {
        if output.len() != self.dimension || self.check(fields).is_err() {
            output.fill(f64::NAN);
            return false;
        }
        true
    }

    fn local(&self, field: &[f64], body: usize) -> Result<Vec<f64>, PlateImpactError> {
        self.check(&[field])?;
        if body >= self.bodies.len() {
            return Err(PlateImpactError::Invalid("impact body index outside set"));
        }
        let n = self.plate_dofs;
        let mut local = zeros(n + 1)?;
        local[..n].copy_from_slice(&field[..n]);
        local[n] = field[n + body];
        Ok(local)
    }

    // Only contact contributions: the caller has already applied the one
    // original plate action. Each local contact derivative is symmetric;
    // transposed global actions still use the plate's explicit transpose.
    fn add_contacts(
        &self,
        q: &[f64],
        v: Option<&[f64]>,
        direction: Option<&[f64]>,
        part: ContactPart,
        output: &mut [f64],
    ) -> Result<(), String> {
        self.check(&[q]).map_err(|error| error.to_string())?;
        if let Some(v) = v {
            self.check(&[v]).map_err(|error| error.to_string())?;
        }
        if let Some(d) = direction {
            self.check(&[d]).map_err(|error| error.to_string())?;
        }
        if output.len() != self.dimension || !finite(output) {
            return Err("invalid impact-set contact action output".into());
        }
        let n = self.plate_dofs;
        for (i, body) in self.bodies.iter().enumerate() {
            let q_local = self.local(q, i).map_err(|error| error.to_string())?;
            let v_local = v
                .map(|v| self.local(v, i))
                .transpose()
                .map_err(|error| error.to_string())?;
            let d_local = direction
                .map(|d| self.local(d, i))
                .transpose()
                .map_err(|error| error.to_string())?;
            let contact = body
                .contact_action(&q_local, v_local.as_deref(), d_local.as_deref(), part)
                .map_err(|error| format!("impact body {i}: {error}"))?;
            for (out, force) in output[..n].iter_mut().zip(&contact[..n]) {
                *out += force;
            }
            output[n + i] += contact[n];
        }
        if !finite(output) {
            return Err("nonfinite combined impact-set contact action".into());
        }
        Ok(())
    }

    /// Observe all contacts at one actual state, with plate storage/loss once.
    /// Refusal returns no partially published observation and changes no state.
    pub fn probe(&self, q: &[f64], v: &[f64]) -> Result<PlateImpactSetProbe, PlateImpactError> {
        self.check(&[q, v])?;
        let mut bodies = Vec::new();
        bodies
            .try_reserve_exact(self.bodies.len())
            .map_err(|_| PlateImpactError::Allocation)?;
        let mut indenter_kinetic_energy_j = 0.0;
        let mut contact_energy_j = 0.0;
        let mut contact_dissipated_power_w = 0.0;
        for (i, body) in self.bodies.iter().enumerate() {
            let observation = body.body_probe(&self.local(q, i)?, &self.local(v, i)?)?;
            indenter_kinetic_energy_j += observation.indenter_kinetic_energy_j;
            contact_energy_j += observation.contact_energy_j;
            contact_dissipated_power_w += observation.contact_dissipated_power_w;
            bodies.push(observation);
        }
        let n = self.plate_dofs;
        let mut work = zeros(n)?;
        self.plate().mass_apply(&v[..n], &mut work);
        let plate_kinetic_energy_j = 0.5 * dot(&v[..n], &work);
        self.plate().internal_force(&q[..n], &mut work);
        let plate_elastic_energy_j = 0.5 * dot(&q[..n], &work);
        self.plate().damping_apply(&v[..n], &mut work);
        let plate_dissipated_power_w = dot(&v[..n], &work);
        let total_energy_j = plate_kinetic_energy_j
            + plate_elastic_energy_j
            + indenter_kinetic_energy_j
            + contact_energy_j;
        if !finite(&[
            plate_kinetic_energy_j,
            plate_elastic_energy_j,
            plate_dissipated_power_w,
            indenter_kinetic_energy_j,
            contact_energy_j,
            contact_dissipated_power_w,
            total_energy_j,
        ]) {
            return Err(PlateImpactError::Invalid(
                "nonfinite physical impact-set observation",
            ));
        }
        Ok(PlateImpactSetProbe {
            plate_kinetic_energy_j,
            plate_elastic_energy_j,
            plate_dissipated_power_w,
            indenter_kinetic_energy_j,
            contact_energy_j,
            contact_dissipated_power_w,
            total_energy_j,
            bodies,
        })
    }
}

impl<L: PlateLoad + ?Sized> SecondOrderProblem for PlateImpactSet<'_, L> {
    fn dimension(&self) -> usize {
        self.dimension
    }

    fn validate_step(&self, step_s: f64) -> Result<(), String> {
        self.plate().validate_step(step_s)?;
        for (i, body) in self.bodies.iter().enumerate() {
            if let Some(source) = body.source() {
                source
                    .validate_kinematics(0.0, 0.0, step_s)
                    .map_err(|error| {
                        format!("impact body {i}: {}", PlateImpactError::Source(error))
                    })?;
            }
        }
        Ok(())
    }

    fn validate_state(&self, q: &[f64], v: &[f64]) -> Result<(), String> {
        self.check(&[q, v]).map_err(|error| error.to_string())?;
        self.plate()
            .validate_state(&q[..self.plate_dofs], &v[..self.plate_dofs])?;
        for (i, body) in self.bodies.iter().enumerate() {
            let local_q = self.local(q, i).map_err(|error| error.to_string())?;
            let local_v = self.local(v, i).map_err(|error| error.to_string())?;
            body.validate_source_kinematics(&local_q, Some(&local_v))
                .map_err(|error| format!("impact body {i}: {error}"))?;
        }
        Ok(())
    }

    fn mass_apply(&self, input: &[f64], output: &mut [f64]) {
        if !self.io(&[input], output) {
            return;
        }
        let n = self.plate_dofs;
        self.plate().mass_apply(&input[..n], &mut output[..n]);
        for (i, body) in self.bodies.iter().enumerate() {
            output[n + i] = body.indenter_mass_kg * input[n + i];
        }
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }

    fn damping_apply(&self, input: &[f64], output: &mut [f64]) {
        if !self.io(&[input], output) {
            return;
        }
        let n = self.plate_dofs;
        self.plate().damping_apply(&input[..n], &mut output[..n]);
        output[n..].fill(0.0);
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }

    fn internal_force(&self, q: &[f64], output: &mut [f64]) {
        if !self.io(&[q], output) {
            return;
        }
        let n = self.plate_dofs;
        self.plate().internal_force(&q[..n], &mut output[..n]);
        output[n..].fill(0.0);
        if self
            .add_contacts(q, None, None, ContactPart::Elastic, output)
            .is_err()
        {
            output.fill(f64::NAN);
        }
    }

    fn tangent_apply(&self, q: &[f64], direction: &[f64], output: &mut [f64]) {
        if !self.io(&[q, direction], output) {
            return;
        }
        let n = self.plate_dofs;
        self.plate()
            .tangent_apply(&q[..n], &direction[..n], &mut output[..n]);
        output[n..].fill(0.0);
        if self
            .add_contacts(
                q,
                None,
                Some(direction),
                ContactPart::ElasticTangent,
                output,
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
        let n = self.plate_dofs;
        self.plate().resistance(&q[..n], &v[..n], &mut output[..n]);
        output[n..].fill(0.0);
        if self
            .add_contacts(q, Some(v), None, ContactPart::Total, output)
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
        let n = self.plate_dofs;
        self.plate().resistance_displacement_apply(
            &q[..n],
            &v[..n],
            &direction[..n],
            &mut output[..n],
        );
        output[n..].fill(0.0);
        if self
            .add_contacts(
                q,
                Some(v),
                Some(direction),
                ContactPart::DisplacementTangent,
                output,
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
        let n = self.plate_dofs;
        self.plate()
            .resistance_velocity_apply(&q[..n], &v[..n], &direction[..n], &mut output[..n]);
        output[n..].fill(0.0);
        if self
            .add_contacts(
                q,
                Some(v),
                Some(direction),
                ContactPart::VelocityTangent,
                output,
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
        let n = self.plate_dofs;
        self.plate().preconditioner_apply(
            &q[..n],
            weights,
            outer,
            inner,
            &residual[..n],
            &mut output[..n],
        );
        for (i, body) in self.bodies.iter().enumerate() {
            let diagonal = weights.mass * body.indenter_mass_kg;
            if !diagonal.is_finite() || diagonal <= 0.0 || !diagonal.recip().is_finite() {
                output.fill(f64::NAN);
                return;
            }
            output[n + i] = residual[n + i] / diagonal;
        }
        // Block Jacobi only: all coupled contact terms remain in the true
        // Newton operator and its separately implemented transpose.
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }
}

/// Exact shared plate and body diagonal mass inverse without stored scratch.
pub struct PlateImpactSetMassPreconditioner<'set, 'model, L: ?Sized> {
    system: &'set PlateImpactSet<'model, L>,
}

impl<L: PlateLoad + ?Sized> FlexiblePreconditioner for PlateImpactSetMassPreconditioner<'_, '_, L> {
    fn apply(&self, iteration: usize, residual: &[f64], output: &mut [f64]) {
        if !self.system.io(&[residual], output) {
            return;
        }
        let n = self.system.plate_dofs;
        self.system.plate().mass_preconditioner().apply(
            iteration,
            &residual[..n],
            &mut output[..n],
        );
        for (i, body) in self.system.bodies.iter().enumerate() {
            output[n + i] = residual[n + i] / body.indenter_mass_kg;
        }
        if !finite(output) {
            output.fill(f64::NAN);
        }
    }
}
