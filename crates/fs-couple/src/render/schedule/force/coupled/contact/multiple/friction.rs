//! One authored tangential coordinate per compliant modal contact.
//!
//! The constitutive owner is fs-tribo, reached through StribeckFriction with
//! equal static/dynamic coefficients. This selects its continuous, monotone
//! regularized Coulomb rung, NOT velocity weakening or set-valued sticking.
//! With y = C^T q, the signed reaction T is applied as -C T and evaluated at
//! the finite-step slip (y1-y0)/dt and the SAME solved normal reaction R.
//! The complete normal/tangent displacement response is retained in both
//! directions; no normal force is frozen and no body is treated as a driver.
//!
//! Fixed authored shapes are not contact discovery, a 2-D friction cone,
//! finite-patch partial slip, a rigid-impact law, or an admitted material card.
//! Source labels remain caller declarations. Loss is mechanical dissipation,
//! not a claimed heat partition. Refinement and convergence remain necessary.

use super::*;
use crate::stribeck_friction::StribeckFriction;

/// Explicit 1-D regularized Coulomb attachment on an existing normal pair.
/// Shapes use that pair's unchanged component indices and modal bases.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalFriction {
    /// Tangential displacement participation on the normal contact's left body.
    /// Units: 1/sqrt(kg); one entry per retained mode.
    pub left_shapes: Vec<f64>,
    /// Tangential participation on its right body; subtracted from the left.
    pub right_shapes: Vec<f64>,
    /// Finite nonnegative coefficient, used for BOTH static and dynamic mu.
    pub coefficient: f64,
    /// Positive ramp speed [m/s], an explicit model parameter, not solver epsilon.
    pub regularization_speed_m_s: f64,
    /// Positive finite applied tangential-force ceiling [N], never a target.
    pub maximum_force_n: f64,
    /// Nonblank authored/source identity. No experimental authority is inferred.
    pub source: String,
}

/// Tangential diagnostics evaluated on the same accepted endpoints as the normals.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct TangentialContactFrame {
    /// Signed reaction T [N]; applied modal forces are -C T.
    pub reaction_n: f64,
    /// Applied reaction minus the shared law on actual motion and normal load [N].
    pub constitutive_residual_n: f64,
    /// Original normal contact's absolute-plus-relative force allowance [N].
    pub force_tolerance_n: f64,
    /// Relative displacement during the complete sample divided by its duration.
    /// This is a step-average velocity [m/s], not an instantaneous endpoint value.
    pub slip_velocity_m_s: f64,
    /// Actual removed work T * delta_y [J], checked nonnegative before publication.
    pub dissipation_j: f64,
}

impl MultiContactModalSystem {
    /// Opt into tangential friction before the first sample, once only.
    ///
    /// Entries correspond exactly to normal-contact construction order; None
    /// explicitly leaves a contact frictionless. Existing normal-only callers
    /// retain their arithmetic. Each tangent inherits its normal contact's root
    /// iteration/residual budgets; max_sweeps now bounds their JOINT solve.
    /// max_setup_terms covers original normal setup plus this extension with
    /// the conservative bound 3*p*(n*(k+2)+(k+1)^2)+4*n*p^2.
    ///
    /// Refuses malformed shapes/coefficients, zero tangent compliance, repeated
    /// or mid-run insertion, cancellation, and insufficient setup work budget.
    pub fn with_friction(
        mut self,
        specifications: Vec<Option<ModalFriction>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || self.friction.is_some() {
            return Err(invalid("friction requires a sample-zero network and may be attached only once"));
        }
        let friction = RegularizedSet::new(
            &self.network, &self.points, specifications, self.config.max_setup_terms, gate,
        )?;
        poll(Some(gate))?;
        self.frame.friction = vec![None; self.points.len()];
        self.friction = Some(TangentialSet::Regularized(friction));
        Ok(self)
    }

    /// Retained authored friction input; None means absent or an unknown index.
    #[must_use]
    pub fn friction_law(&self, index: usize) -> Option<&ModalFriction> {
        self.friction.as_ref()?.regularized_law(index)
    }
}

struct TangentialPoint {
    spec: ModalFriction,
    law: StribeckFriction,
    column: Vec<f64>,
}

pub(super) struct RegularizedSet {
    points: Vec<Option<TangentialPoint>>,
    // Row-major normal-from-tangent, tangent-from-normal and tangent-from-tangent.
    // Both mixed directions are evaluated, not assumed numerically reciprocal.
    nt: Vec<f64>,
    tn: Vec<f64>,
    tt: Vec<f64>,
    old_y: Vec<f64>,
    free_y: Vec<f64>,
    reactions: Vec<f64>,
    staged: Vec<Option<TangentialContactFrame>>,
}

impl RegularizedSet {
    fn new(
        network: &CoupledModalSystem,
        normals: &[ContactPoint],
        specifications: Vec<Option<ModalFriction>>,
        max_setup_terms: usize,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        let p = normals.len();
        if specifications.len() != p {
            return Err(invalid("friction entries must match normal-contact construction order"));
        }
        let n = network.mode_count();
        let k = network.columns.len();
        let terms = n.checked_mul(k+2).and_then(|v| v.checked_add((k+1)*(k+1)))
            .and_then(|v| v.checked_mul(3*p))
            .and_then(|v| n.checked_mul(4*p*p).and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("normal/tangential setup work overflow"))?;
        if terms > max_setup_terms {
            return Err(invalid("normal/tangential setup exceeds max_setup_terms"));
        }
        let mut points = Vec::with_capacity(p);
        for (normal, spec) in normals.iter().zip(specifications) {
            poll(Some(gate))?;
            let Some(spec) = spec else { points.push(None); continue; };
            if !spec.maximum_force_n.is_finite() || spec.maximum_force_n <= 0.0
                || spec.source.trim().is_empty() {
                return Err(invalid("friction requires a positive finite force ceiling and a source label"));
            }
            let law = StribeckFriction::try_new(
                spec.coefficient, spec.coefficient, spec.regularization_speed_m_s,
            ).map_err(invalid)?;
            let mut column = vec![0.0; n];
            for (attachment, shapes, sign) in [
                (&normal.contact.left, &spec.left_shapes, 1.0),
                (&normal.contact.right, &spec.right_shapes, -1.0),
            ] {
                if shapes.len() != attachment.shapes.len() || shapes.iter().any(|x| !x.is_finite()) {
                    return Err(invalid("tangential shapes must match the normal pair's finite modal bases"));
                }
                for (i, b) in shapes.iter().enumerate() {
                    let at = network.offsets[attachment.component] + i;
                    column[at] = finite(column[at] + sign*b)?;
                }
            }
            points.push(Some(TangentialPoint { spec, law, column }));
        }
        let mut nt = vec![0.0; p*p];
        let mut tn = vec![0.0; p*p];
        let mut tt = vec![0.0; p*p];
        for j in 0..p {
            poll(Some(gate))?;
            if let Some(point) = &points[j] {
                let response = network_response(network, &point.column, gate)?;
                for i in 0..p {
                    poll(Some(gate))?;
                    nt[i*p+j] = dot(&normals[i].column, &response)?;
                    if let Some(other) = &points[i] {
                        tt[i*p+j] = dot(&other.column, &response)?;
                    }
                }
                if tt[j*p+j] <= 0.0 {
                    return Err(invalid("each tangent requires positive representable network compliance"));
                }
            }
            let response = network_response(network, &normals[j].column, gate)?;
            for i in 0..p {
                poll(Some(gate))?;
                if let Some(point) = &points[i] { tn[i*p+j] = dot(&point.column, &response)?; }
            }
        }
        Ok(Self { points, nt, tn, tt, old_y: vec![0.0; p], free_y: vec![0.0; p],
            reactions: vec![0.0; p], staged: vec![None; p] })
    }

    pub(super) fn prepare(&mut self, old_q: &[f64], free_q: &[f64], gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        self.reactions.fill(0.0);
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            if let Some(point) = point {
                self.old_y[i] = dot(&point.column, old_q)?;
                self.free_y[i] = dot(&point.column, free_q)?;
            }
        }
        Ok(())
    }

    pub(super) fn normal_endpoint(&self, i: usize, mut x: f64) -> Result<f64, ModalCouplingError> {
        let p = self.points.len();
        for j in 0..p {
            if self.reactions[j] != 0.0 { x = finite(x - self.nt[i*p+j]*self.reactions[j])?; }
        }
        Ok(x)
    }

    fn delta(&self, i: usize, normals: &[f64], exclude_self: bool) -> Result<f64, ModalCouplingError> {
        let p = self.points.len();
        let mut delta = finite(self.free_y[i] - self.old_y[i])?;
        for j in 0..p {
            if normals[j] != 0.0 { delta = finite(delta - self.tn[i*p+j]*normals[j])?; }
            if !(exclude_self && i == j) && self.reactions[j] != 0.0 {
                delta = finite(delta - self.tt[i*p+j]*self.reactions[j])?;
            }
        }
        Ok(delta)
    }

    pub(super) fn solve_coordinate(&mut self, i: usize, normals: &[f64], dt: f64,
        config: ModalContactConfig, gate: Option<&CancelGate>) -> Result<(), ModalCouplingError>
    {
        let Some(point) = &self.points[i] else { return Ok(()); };
        let delta = self.delta(i, normals, true)?;
        self.reactions[i] = solve_tangent(point, normals[i], delta,
            self.tt[i*self.points.len()+i], dt, config, gate)?;
        Ok(())
    }

    pub(super) fn residuals(&self, normals: &[f64], points: &[ContactPoint], dt: f64,
        worst: &mut (f64, f64), gate: Option<&CancelGate>) -> Result<(), ModalCouplingError>
    {
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            if let Some(point) = point {
                let velocity = finite(self.delta(i, normals, false)?/dt)?;
                let expected = point.law.traction(velocity, normals[i]).map_err(invalid)?;
                let residual = finite(self.reactions[i] - expected)?;
                let tolerance = force_tolerance(self.reactions[i], expected, points[i].config)?;
                if residual.abs()/tolerance > worst.0.abs()/worst.1 { *worst = (residual, tolerance); }
            }
        }
        Ok(())
    }

    pub(super) fn add_forces(&self, forces: &mut [f64], gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            if self.reactions[i] != 0.0 {
                if let Some(point) = point {
                    for (f, b) in forces.iter_mut().zip(&point.column) { *f = finite(*f - b*self.reactions[i])?; }
                }
            }
        }
        Ok(())
    }

    pub(super) fn stage(&mut self, network: &CoupledModalSystem, normals: &[f64],
        points: &[ContactPoint], sweeps: usize, gate: Option<&CancelGate>)
        -> Result<f64, ModalCouplingError>
    {
        let mut total = 0.0;
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            let Some(point) = point else { continue; };
            // Use the actual component displacements, the same dual pairing as
            // external work, not a free prediction or endpoint instantaneous speed.
            let delta = dot(&point.column, &network.free_delta)?;
            let velocity = finite(delta/network.dt)?;
            let applied = self.reactions[i];
            let expected = point.law.traction(velocity, normals[i]).map_err(invalid)?;
            let tolerance = force_tolerance(applied, expected, points[i].config)?;
            let residual = finite(applied - expected)?;
            if residual.abs() > tolerance {
                return Err(ModalCouplingError::ContactSolve {
                    residual_n: residual, tolerance_n: tolerance, iterations: sweeps,
                });
            }
            limit("contact tangential force", applied.abs(), point.spec.maximum_force_n)?;
            let loss = finite(applied*delta)?;
            if loss < 0.0 { return Err(invalid("tangential contact work must not create energy")); }
            total = finite(total + loss)?;
            self.staged[i] = Some(TangentialContactFrame {
                reaction_n: applied, constitutive_residual_n: residual, force_tolerance_n: tolerance,
                slip_velocity_m_s: velocity, dissipation_j: loss,
            });
        }
        Ok(total)
    }

    pub(super) fn publish(&mut self, frames: &mut Vec<Option<TangentialContactFrame>>) {
        std::mem::swap(frames, &mut self.staged);
    }
}

// Positive scalar compliance and the constant-mu regularized law make this
// coordinate residual monotone. The whole normal/tangential problem need not
// be contractive: only JOINT residual acceptance can publish a sample.
fn solve_tangent(point: &TangentialPoint, normal: f64, free_delta: f64,
    compliance: f64, dt: f64, config: ModalContactConfig, gate: Option<&CancelGate>)
    -> Result<f64, ModalCouplingError>
{
    let free_force = point.law.traction(finite(free_delta/dt)?, normal).map_err(invalid)?;
    if free_force == 0.0 { return Ok(0.0); }
    let sign = free_force.signum();
    let evaluate = |magnitude: f64| -> Result<(f64, f64), ModalCouplingError> {
        poll(gate)?;
        let applied = sign*magnitude;
        let velocity = finite(finite(free_delta - compliance*applied)?/dt)?;
        let expected = point.law.traction(velocity, normal).map_err(invalid)?;
        Ok((finite(sign*(applied - expected))?, force_tolerance(applied, expected, config)?))
    };
    let mut lo = 0.0;
    let mut hi = free_force.abs().min(point.spec.maximum_force_n);
    let (mut residual, mut tolerance) = evaluate(hi)?;
    if residual.abs() <= tolerance || residual < 0.0 {
        // A ceiling-limited coordinate remains a TRIAL. Other contacts may
        // relieve it; otherwise the joint constitutive check refuses it.
        return Ok(sign*hi);
    }
    let mut used = 0;
    for iteration in 1..=config.max_iterations {
        used = iteration;
        let mid = f64::midpoint(lo, hi);
        (residual, tolerance) = evaluate(mid)?;
        if residual.abs() <= tolerance { return Ok(sign*mid); }
        if mid == lo || mid == hi { break; }
        if residual < 0.0 { lo = mid; } else { hi = mid; }
    }
    Err(ModalCouplingError::ContactSolve { residual_n: sign*residual, tolerance_n: tolerance, iterations: used })
}

// Both physical models share the existing normal-contact transaction. Dispatch
// leaves the v5 regularized arithmetic untouched and never changes models at run time.
pub(super) enum TangentialSet {
    Regularized(RegularizedSet),
    Coulomb(super::coulomb::CoulombSet),
}
impl TangentialSet {
    pub(super) fn regularized_law(&self, i: usize) -> Option<&ModalFriction> {
        match self { Self::Regularized(s) => s.points.get(i)?.as_ref().map(|p| &p.spec), Self::Coulomb(_) => None }
    }
    pub(super) fn coulomb_law(&self, i: usize) -> Option<&super::coulomb::ModalCoulombFriction> {
        match self { Self::Coulomb(s) => s.law(i), Self::Regularized(_) => None }
    }
    pub(super) fn prepare(&mut self, old: &[f64], free: &[f64], gate: Option<&CancelGate>) -> Result<(), ModalCouplingError> {
        match self { Self::Regularized(s) => s.prepare(old, free, gate), Self::Coulomb(s) => s.prepare(old, free, gate) }
    }
    pub(super) fn normal_endpoint(&self, i: usize, x: f64) -> Result<f64, ModalCouplingError> {
        match self { Self::Regularized(s) => s.normal_endpoint(i, x), Self::Coulomb(s) => s.normal_endpoint(i, x) }
    }
    pub(super) fn solve_coordinate(&mut self, i: usize, normals: &[f64], dt: f64, config: ModalContactConfig,
        gate: Option<&CancelGate>) -> Result<(), ModalCouplingError>
    {
        match self { Self::Regularized(s) => s.solve_coordinate(i, normals, dt, config, gate),
            Self::Coulomb(s) => s.solve_coordinate(i, normals, dt, config, gate) }
    }
    pub(super) fn residuals(&self, normals: &[f64], points: &[ContactPoint], dt: f64, worst: &mut (f64, f64),
        gate: Option<&CancelGate>) -> Result<(), ModalCouplingError>
    {
        match self { Self::Regularized(s) => s.residuals(normals, points, dt, worst, gate),
            Self::Coulomb(s) => s.residuals(normals, points, dt, worst, gate) }
    }
    pub(super) fn add_forces(&self, forces: &mut [f64], gate: Option<&CancelGate>) -> Result<(), ModalCouplingError> {
        match self { Self::Regularized(s) => s.add_forces(forces, gate), Self::Coulomb(s) => s.add_forces(forces, gate) }
    }
    pub(super) fn stage(&mut self, network: &CoupledModalSystem, normals: &[f64], points: &[ContactPoint], sweeps: usize,
        gate: Option<&CancelGate>) -> Result<f64, ModalCouplingError>
    {
        match self { Self::Regularized(s) => s.stage(network, normals, points, sweeps, gate),
            Self::Coulomb(s) => s.stage(network, normals, points, sweeps, gate) }
    }
    pub(super) fn publish(&mut self, regularized: &mut Vec<Option<TangentialContactFrame>>,
        coulomb: &mut Vec<Option<super::coulomb::CoulombContactFrame>>)
    {
        match self { Self::Regularized(s) => s.publish(regularized), Self::Coulomb(s) => s.publish(coulomb) }
    }
}
