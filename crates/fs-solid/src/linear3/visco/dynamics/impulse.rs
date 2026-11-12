//! Instantaneous momentum changes on the original tensor-Prony trajectory.

use super::*;
use fs_math::dd::Dd;
use fs_solver::LinearOp;

/// A zero-duration event in full node-major coordinates. The material history,
/// displacement, physical clock and accepted interval count are unchanged.
#[derive(Clone, Debug)]
pub struct ViscoDynamicImpulse {
    pub time_s: f64,
    pub iterations: usize,
    /// Original full consistent mass times the actual represented delta-v.
    /// Rounded for reporting; the residual retains compensated endpoint differences.
    pub momentum_change_kg_m_s: Vec<f64>,
    /// Impulse ON the solid from supports; zero on free coordinates.
    pub support_impulse_kg_m_s: Vec<f64>,
    /// M delta-v minus applied and support impulse, in full nodal order.
    pub impulse_residual_kg_m_s: Vec<f64>,
    pub maximum_impulse_residual_kg_m_s: f64,
    pub kinetic_before_j: f64,
    pub kinetic_after_j: f64,
    /// Average velocity dotted with M delta-v, without subtracting large K.
    pub kinetic_change_j: f64,
    pub external_work_j: f64,
    pub support_work_j: f64,
    /// Average velocity dotted with the actual impulse residual; never heat.
    pub balance_residual_j: f64,
}

/// Admitted impulse staging with exclusive access to the accepted trajectory.
/// Outer contact checks can inspect both endpoints without cloning material
/// histories. Dropping this handle leaves the original trajectory unchanged.
#[must_use]
pub struct ViscoDynamicImpulseProposal<'s, 'p, 'a> {
    original: &'s mut ViscoDynamicState<'p, 'a>,
    velocity_mps: Vec<f64>,
    external_work_j: f64,
    support_work_j: f64,
    report: ViscoDynamicImpulse,
}

impl<'p, 'a> ViscoDynamicImpulseProposal<'_, 'p, 'a> {
    pub fn original(&self) -> &ViscoDynamicState<'p, 'a> {
        self.original
    }

    pub fn velocity_mps(&self) -> &[f64] {
        &self.velocity_mps
    }

    pub fn report(&self) -> &ViscoDynamicImpulse {
        &self.report
    }

    /// Publish the already checked velocity and work, then release the borrow.
    /// This performs no solve, allocation, history copy or cancellation check.
    /// A coupled owner completes all partner checks before publishing either.
    pub fn publish(self) -> ViscoDynamicImpulse {
        self.original.velocity_mps = self.velocity_mps;
        self.original.kinetic_j = self.report.kinetic_after_j;
        self.original.external_work_j = self.external_work_j;
        self.original.support_work_j = self.support_work_j;
        self.report
    }
}

// Reuse the native sparse solve while borrowing its original mass. No mass
// assembly, matrix copy, dense inverse or constitutive evaluation is needed.
struct BorrowedMass<'a>(&'a Csr);
impl LinearOp for BorrowedMass<'_> {
    fn n(&self) -> usize {
        self.0.nrows()
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        self.0.spmv(x, y);
    }
    // The consistent mass matrix is symmetric.
    fn apply_transpose(&self, x: &[f64], y: &mut [f64]) {
        self.apply(x, y);
    }
}

fn finite_dd(value: Dd) -> Result<Dd, ViscoElementError> {
    finite(value.hi)?;
    finite(value.lo)?;
    finite(value.to_f64())?;
    Ok(value)
}

fn work(before: &[f64], after: &[f64], impulse: &[f64]) -> Result<f64, ViscoElementError> {
    // Retain a midpoint between adjacent f64 velocities before pairing it with
    // impulse. Rounding that midpoint first can erase quadratic kinetic work.
    let sum =
        before
            .iter()
            .zip(after)
            .zip(impulse)
            .fold(Dd::ZERO, |sum, ((&old, &new), &applied)| {
                let mean = Dd::from_f64(0.5 * old) + Dd::from_f64(0.5 * new);
                sum + mean * Dd::from_f64(applied)
            });
    Ok(finite_dd(sum)?.to_f64())
}

impl<'p, 'a> ViscoDynamicState<'p, 'a> {
    /// Original admitted nodes; current positions add `material_history()`'s
    /// full displacement. Impulse geometry never comes from unrelated nodes.
    pub fn reference_nodes_m(&self) -> &[[f64; 3]] {
        self.prepared.system.reference_nodes_m()
    }

    /// Full nodal velocity response to a world impulse [kg m/s], with every
    /// support velocity held. This is the native M_ff^-1 action needed by a
    /// contact mobility, lifted to full coordinates with zero supported rows.
    /// All off-diagonal consistent-mass terms are retained. No state changes.
    ///
    /// The same full reference mass required by `apply_momentum_impulse` must
    /// be retained. Interpolation traces and the contact law remain with their
    /// physical owner; this method does not infer a point from a nodal vector.
    pub fn momentum_response(
        &self,
        nodal_impulse_kg_m_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<Vec<f64>, ViscoElementError> {
        self.momentum_response_impl(nodal_impulse_kg_m_s, cx)
            .map(|(velocity, _)| velocity)
    }

    fn momentum_response_impl(
        &self,
        nodal_impulse_kg_m_s: &[f64],
        cx: &Cx<'_>,
    ) -> Result<(Vec<f64>, usize), ViscoElementError> {
        checkpoint(cx)?;
        let system = self.prepared.system;
        if system.full_reference_mass().is_none()
            || nodal_impulse_kg_m_s.len() != system.full_dofs
            || nodal_impulse_kg_m_s.iter().any(|j| !j.is_finite())
        {
            return Err(invalid(
                "impulses require retained full reference mass and complete finite nodal components",
            ));
        }
        let free = &system.equilibrium.free_dofs;
        let rhs: Vec<_> = free.iter().map(|&d| nodal_impulse_kg_m_s[d]).collect();
        let scale = finite(euclidean_norm(rhs.iter().copied()))?;
        let mut response = vec![0.0; system.full_dofs];
        if scale == 0.0 {
            checkpoint(cx)?;
            return Ok((response, 0));
        }
        let mass = &system.equilibrium.mass;
        let (delta, iterations) = solve_linear(
            &BorrowedMass(mass),
            &Jacobi::new(mass),
            &rhs,
            scale,
            self.prepared.config.solve,
            cx,
        )?;
        for (i, (&d, value)) in free.iter().zip(delta).enumerate() {
            if i % 256 == 0 {
                checkpoint(cx)?;
            }
            response[d] = finite(value)?;
        }
        checkpoint(cx)?;
        Ok((response, iterations))
    }

    /// Apply an instantaneous nodal impulse through the original reference
    /// mass, retaining the complete accepted material/thermal history and clock.
    /// Supports keep their actual velocities, including nonzero moving rates;
    /// full mass coupling supplies their reaction impulses and signed work.
    /// Applied impulses on supported rows enter those reactions as well.
    ///
    /// Free bodies already retain full mass. Constrained sources must use
    /// `ViscoElements::new_with_full_reference_mass`. Both force and impulse
    /// responses then use the same inventory, without remeshing or lumping.
    /// The caller owns any restitution/friction law and partner impulse.
    ///
    /// No material interval passes: Maxwell memory, Kelvin heat, thermal strain,
    /// displacement, source work, time and lifetime step budget all survive.
    /// The next native interval starts from the returned velocity. Solve,
    /// actual-velocity momentum, work and cancellation checks precede mutation.
    pub fn apply_momentum_impulse(
        &mut self,
        nodal_impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicImpulse, ViscoElementError> {
        Ok(self
            .prepare_momentum_impulse(nodal_impulse_kg_m_s, maximum_impulse_error_kg_m_s, cx)?
            .publish())
    }

    /// Prepare the same native impulse without changing accepted motion.
    /// The exclusive borrow prevents state changes while outer contact or
    /// partner audits inspect the original and proposed endpoints. Every
    /// native solve, momentum, work and cancellation gate finishes here.
    /// Dropping the proposal discards it; publishing is infallible.
    pub fn prepare_momentum_impulse<'s>(
        &'s mut self,
        nodal_impulse_kg_m_s: &[f64],
        maximum_impulse_error_kg_m_s: f64,
        cx: &Cx<'_>,
    ) -> Result<ViscoDynamicImpulseProposal<'s, 'p, 'a>, ViscoElementError> {
        checkpoint(cx)?;
        if !maximum_impulse_error_kg_m_s.is_finite() || maximum_impulse_error_kg_m_s < 0.0 {
            return Err(invalid("finite nonnegative impulse allowance required"));
        }
        let (response, iterations) = self.momentum_response_impl(nodal_impulse_kg_m_s, cx)?;
        let system = self.prepared.system;
        let n = system.full_dofs;
        let mut velocity = self.velocity_mps.clone();
        let mut delta = vec![Dd::ZERO; n];
        for d in 0..n {
            if d % 256 == 0 {
                checkpoint(cx)?;
            }
            if response[d] != 0.0 {
                velocity[d] = finite(velocity[d] + response[d])?;
            }
            // A large accepted velocity can absorb the old velocity in f64
            // subtraction. Audit the two represented endpoints before rounding.
            delta[d] = finite_dd(Dd::from_f64(velocity[d]) - Dd::from_f64(self.velocity_mps[d]))?;
        }
        let mass = system
            .full_reference_mass()
            .ok_or_else(|| invalid("impulses require retained full reference mass"))?;
        let mut momentum_change_kg_m_s = vec![0.0; n];
        let mut support_impulse_kg_m_s = vec![0.0; n];
        let mut impulse_residual_kg_m_s = vec![0.0; n];
        let mut maximum = Dd::ZERO;
        let mut kinetic_change = Dd::ZERO;
        let mut balance_residual = Dd::ZERO;
        for d in 0..n {
            if d % 256 == 0 {
                checkpoint(cx)?;
            }
            let (columns, values) = mass.row(d);
            let momentum = finite_dd(
                columns
                    .iter()
                    .zip(values)
                    .fold(Dd::ZERO, |sum, (&column, &value)| {
                        sum + Dd::from_f64(value) * delta[column]
                    }),
            )?;
            momentum_change_kg_m_s[d] = momentum.to_f64();
            let difference = finite_dd(momentum - Dd::from_f64(nodal_impulse_kg_m_s[d]))?;
            if self.prepared.reduced[d].is_none() {
                support_impulse_kg_m_s[d] = finite(difference.to_f64())?;
            }
            // Check the returned reaction too: rounding a large opposing
            // support impulse must not erase a smaller physical M delta-v.
            let residual = finite_dd(difference - Dd::from_f64(support_impulse_kg_m_s[d]))?;
            impulse_residual_kg_m_s[d] = residual.to_f64();
            if maximum.lt(residual.abs()) {
                maximum = residual.abs();
            }
            // Pair the compensated mass action/residual with the endpoint
            // midpoint before rounding either quantity to the public report.
            let mean = Dd::from_f64(0.5 * self.velocity_mps[d]) + Dd::from_f64(0.5 * velocity[d]);
            kinetic_change = kinetic_change + mean * momentum;
            balance_residual = balance_residual + mean * residual;
        }
        if Dd::from_f64(maximum_impulse_error_kg_m_s).lt(maximum) {
            return Err(ViscoElementError::ImpulseMomentum {
                residual_kg_m_s: maximum.to_f64(),
                allowance_kg_m_s: maximum_impulse_error_kg_m_s,
            });
        }
        let kinetic_after_j = kinetic(system, &velocity)?;
        let kinetic_change_j = finite_dd(kinetic_change)?.to_f64();
        let external_work_j = work(&self.velocity_mps, &velocity, nodal_impulse_kg_m_s)?;
        let support_work_j = work(&self.velocity_mps, &velocity, &support_impulse_kg_m_s)?;
        let balance_residual = finite_dd(balance_residual)?;
        let balance_residual_j = balance_residual.to_f64();
        let config = self.prepared.config;
        let allowance_j = finite(
            config.absolute_energy_tolerance_j
                + config.relative_energy_tolerance
                    * (external_work_j.abs() + support_work_j.abs() + kinetic_change_j.abs()),
        )?;
        if Dd::from_f64(allowance_j).lt(balance_residual.abs()) {
            return Err(ViscoElementError::DynamicEnergy {
                residual_j: balance_residual_j,
                allowance_j,
            });
        }
        let external_work = finite(self.external_work_j + external_work_j)?;
        let support_work = finite(self.support_work_j + support_work_j)?;
        let report = ViscoDynamicImpulse {
            time_s: self.time_s,
            iterations,
            momentum_change_kg_m_s,
            support_impulse_kg_m_s,
            impulse_residual_kg_m_s,
            maximum_impulse_residual_kg_m_s: maximum.to_f64(),
            kinetic_before_j: self.kinetic_j,
            kinetic_after_j,
            kinetic_change_j,
            external_work_j,
            support_work_j,
            balance_residual_j,
        };
        checkpoint(cx)?;
        Ok(ViscoDynamicImpulseProposal {
            original: self,
            velocity_mps: velocity,
            external_work_j: external_work,
            support_work_j: support_work,
            report,
        })
    }
}
