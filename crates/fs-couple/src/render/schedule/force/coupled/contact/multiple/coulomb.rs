//! Set-valued Coulomb friction in a caller-declared orthonormal tangent plane.
//!
//! This is a finite-step contact-port embedding, not a rigid-impact stepper.
//! The fs-tribo Coulomb owner supplies the constant-mu capacity. The solve uses
//! T = project_disk(T + delta_y / L, mu*R), with L the trace of the local
//! tangential displacement response. Inside the disk this requires delta_y=0;
//! on its boundary T opposes the applied slip force (-C T). Normal forces and
//! BOTH tangent directions participate in the existing joint contact sweeps.
//! No velocity regularization, friction pyramid, or independent axis clamps.
//!
//! A projected block iteration is a bounded baseline, not a general conic
//! solver or a global convergence theorem. It refuses unresolvable blocks.
//! Sticking refers to step-average slip within the explicit velocity tolerance,
//! not a claim about every instant between endpoints. Geometry and tangent
//! orthonormality are authored; no discovery, rotating-frame transport,
//! partial-slip memory, unequal static/kinetic coefficients or heat partition.

use super::*;
use fs_tribo::{ContactFrame, FrictionLaw, InterfaceSystemRef, TangentialSlip};

/// A dry constant-coefficient law on two tangent directions of a normal pair.
#[derive(Clone, Debug, PartialEq)]
pub struct ModalCoulombFriction {
    /// Two displacement/force maps [1/sqrt(kg)] on the normal pair's left body.
    pub left_shapes: [Vec<f64>; 2],
    /// Corresponding maps on its right body; subtracted from the left maps.
    pub right_shapes: [Vec<f64>; 2],
    /// Common static and kinetic coefficient. Nonnegative and finite.
    pub coefficient: f64,
    /// Positive finite norm ceiling [N], not a replacement Coulomb capacity.
    pub maximum_force_n: f64,
    /// Positive finite step-average sticking tolerance [m/s]; not a creep law.
    pub velocity_tolerance_m_s: f64,
    /// Ordered dry-interface, history and source declarations retained verbatim.
    /// The fs-tribo owner checks these without upgrading caller authority.
    pub interface: InterfaceSystemRef,
}

/// Regime of the accepted discrete graph, not of unobserved intra-step motion.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoulombContactRegime {
    /// Zero friction capacity (zero load or zero coefficient).
    Inactive,
    /// Step-average slip norm is within the declared velocity tolerance.
    Sticking,
    /// Nonzero resolved slip with saturated, aligned reaction.
    Sliding,
}

/// Diagnostics from the SAME accepted endpoints as the normal-contact report.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct CoulombContactFrame {
    /// Reaction T [N]; forces on the network are -C T.
    pub reaction_n: [f64; 2],
    /// Actual relative displacement divided by the sample duration [m/s].
    pub slip_velocity_m_s: [f64; 2],
    /// Physical disk radius mu*R [N], independent of the force budget.
    pub capacity_n: f64,
    /// Maximum projection, cone-excess and resolved sliding-direction residual [N].
    pub graph_residual_n: f64,
    /// Original normal contact's absolute-plus-relative force allowance [N].
    pub force_tolerance_n: f64,
    /// Regime only after graph, cone, velocity and work checks pass.
    pub regime: CoulombContactRegime,
    /// Nonnegative graph dissipation mu*R*|delta_y| [J].
    pub dissipation_j: f64,
    /// Signed ACTUAL removed work T dot delta_y [J].
    pub removed_work_j: f64,
    /// Actual removed work minus graph dissipation [J], never hidden in heat.
    pub work_residual_j: f64,
    /// Work allowance derived from the admitted force/velocity tolerances [J].
    pub work_tolerance_j: f64,
}

impl MultiContactModalSystem {
    /// Select two-direction set-valued Coulomb friction before the first sample.
    ///
    /// One Some(law) or explicit None per normal contact. This is a distinct
    /// model from with_friction's regularized 1-D law; the APIs cannot be mixed
    /// or attached twice. The existing root/sweep budgets bound local projected
    /// iterations and their joint normal/tangential solve. Setup is bounded by
    /// 5*p*(n*(k+2)+(k+1)^2)+9*n*p^2, including the original normal setup.
    /// Each local tangent block must have a positive-definite displacement
    /// response; redundancy BETWEEN contacts remains legal.
    pub fn with_coulomb_friction(
        mut self,
        specifications: Vec<Option<ModalCoulombFriction>>,
        gate: &CancelGate,
    ) -> Result<Self, ModalCouplingError> {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || self.friction.is_some() {
            return Err(invalid("Coulomb friction requires sample zero and no previously attached friction set"));
        }
        let set = CoulombSet::new(&self.network, &self.points, specifications,
            self.config.max_setup_terms, gate)?;
        poll(Some(gate))?;
        self.frame.coulomb_friction = vec![None; self.points.len()];
        self.friction = Some(friction::TangentialSet::Coulomb(set));
        Ok(self)
    }

    /// Retained dry-interface and tangent declarations, without material inference.
    #[must_use]
    pub fn coulomb_friction_law(&self, index: usize) -> Option<&ModalCoulombFriction> {
        self.friction.as_ref()?.coulomb_law(index)
    }
}

struct CoulombPoint {
    spec: ModalCoulombFriction,
    columns: [Vec<f64>; 2],
    unit_capacity: f64,
    response_trace: f64,
}
impl CoulombPoint {
    fn capacity(&self, normal: f64) -> Result<f64, ModalCouplingError> {
        if !normal.is_finite() || normal < 0.0 { return Err(invalid("Coulomb load must be finite and nonnegative")); }
        finite(self.unit_capacity * normal)
    }

    // The SAME graph checks govern local, joint and actual-state acceptance.
    // Near stick/slide transition, a small projected residual need not imply a
    // small angular force error; resolved sliding is checked separately.
    fn graph(&self, x: [f64; 2], delta: [f64; 2], radius: f64, dt: f64,
        config: ModalContactConfig) -> Result<GraphEvaluation, ModalCouplingError>
    {
        let projected = projected_step(x, delta, self.response_trace, radius)?;
        let magnitude = norm(x)?;
        let mut residual = norm([finite(x[0]-projected[0])?, finite(x[1]-projected[1])?])?
            .max(magnitude-radius);
        let tolerance = force_tolerance(magnitude, norm(projected)?, config)?;
        let velocity = [finite(delta[0]/dt)?, finite(delta[1]/dt)?];
        let speed = norm(velocity)?;
        let regime = if radius == 0.0 { CoulombContactRegime::Inactive }
            else if speed <= self.spec.velocity_tolerance_m_s { CoulombContactRegime::Sticking }
            else {
                let expected = [finite(radius*(velocity[0]/speed))?, finite(radius*(velocity[1]/speed))?];
                residual = residual.max(norm([finite(x[0]-expected[0])?, finite(x[1]-expected[1])?])?);
                CoulombContactRegime::Sliding
            };
        Ok(GraphEvaluation { projected, residual, tolerance, velocity, regime })
    }
}
struct GraphEvaluation {
    projected: [f64; 2],
    residual: f64,
    tolerance: f64,
    velocity: [f64; 2],
    regime: CoulombContactRegime,
}

pub(super) struct CoulombSet {
    points: Vec<Option<CoulombPoint>>,
    // p by 2p, 2p by p, and 2p by 2p. Mixed directions are both measured.
    nt: Vec<f64>,
    tn: Vec<f64>,
    tt: Vec<f64>,
    old_y: Vec<f64>,
    free_y: Vec<f64>,
    reactions: Vec<f64>,
    staged: Vec<Option<CoulombContactFrame>>,
}

impl CoulombSet {
    fn new(network: &CoupledModalSystem, normals: &[ContactPoint],
        specifications: Vec<Option<ModalCoulombFriction>>, max_setup_terms: usize, gate: &CancelGate)
        -> Result<Self, ModalCouplingError>
    {
        let p = normals.len();
        if specifications.len() != p { return Err(invalid("Coulomb entries must match the normal-contact set")); }
        let n = network.mode_count();
        let k = network.columns.len();
        let terms = n.checked_mul(k+2).and_then(|v| v.checked_add((k+1)*(k+1)))
            .and_then(|v| v.checked_mul(5*p))
            .and_then(|v| n.checked_mul(9*p*p).and_then(|w| v.checked_add(w)))
            .ok_or_else(|| invalid("Coulomb setup work overflow"))?;
        if terms > max_setup_terms { return Err(invalid("Coulomb setup exceeds max_setup_terms")); }
        let frame = ContactFrame::new([0.0, 0.0, 1.0]).map_err(|_| invalid("Coulomb tangent frame refused"))?;
        let at_rest = TangentialSlip::new(&frame, [0.0; 3]).map_err(|_| invalid("Coulomb rest slip refused"))?;
        let mut points = Vec::with_capacity(p);
        for (normal, spec) in normals.iter().zip(specifications) {
            poll(Some(gate))?;
            let Some(spec) = spec else { points.push(None); continue; };
            if !spec.maximum_force_n.is_finite() || spec.maximum_force_n <= 0.0
                || !spec.velocity_tolerance_m_s.is_finite() || spec.velocity_tolerance_m_s <= 0.0 {
                return Err(invalid("Coulomb force and velocity budgets must be positive and finite"));
            }
            // The constitutive owner admits coefficient/interface and supplies
            // the unit-load capacity. Constant mu then scales linearly with R.
            let unit_capacity = FrictionLaw::Coulomb { static_mu: spec.coefficient, kinetic_mu: spec.coefficient }
                .evaluate(&spec.interface, 1.0, at_rest)
                .map_err(|_| invalid("fs-tribo refused the declared dry Coulomb law/interface"))?.static_limit;
            let mut columns = [vec![0.0; n], vec![0.0; n]];
            for axis in 0..2 {
                for (attachment, shapes, sign) in [
                    (&normal.contact.left, &spec.left_shapes[axis], 1.0),
                    (&normal.contact.right, &spec.right_shapes[axis], -1.0),
                ] {
                    if shapes.len() != attachment.shapes.len() || shapes.iter().any(|x| !x.is_finite()) {
                        return Err(invalid("Coulomb tangent maps must match the normal pair's finite modal bases"));
                    }
                    for (mode, b) in shapes.iter().enumerate() {
                        let at = network.offsets[attachment.component]+mode;
                        columns[axis][at] = finite(columns[axis][at]+sign*b)?;
                    }
                }
            }
            points.push(Some(CoulombPoint { spec, columns, unit_capacity, response_trace: 0.0 }));
        }
        let t = 2*p;
        let mut nt = vec![0.0; p*t];
        let mut tn = vec![0.0; t*p];
        let mut tt = vec![0.0; t*t];
        for j in 0..p {
            poll(Some(gate))?;
            let response = network_response(network, &normals[j].column, gate)?;
            for (i, point) in points.iter().enumerate() {
                if let Some(point) = point {
                    for axis in 0..2 { tn[(2*i+axis)*p+j] = dot(&point.columns[axis], &response)?; }
                }
            }
            if let Some(point) = &points[j] {
                for axis in 0..2 {
                    let col = 2*j+axis;
                    let response = network_response(network, &point.columns[axis], gate)?;
                    for i in 0..p {
                        poll(Some(gate))?;
                        nt[i*t+col] = dot(&normals[i].column, &response)?;
                        if let Some(other) = &points[i] {
                            for row_axis in 0..2 { tt[(2*i+row_axis)*t+col] = dot(&other.columns[row_axis], &response)?; }
                        }
                    }
                }
            }
        }
        for (i, point) in points.iter_mut().enumerate() {
            if let Some(point) = point {
                let u = 2*i;
                let a = tt[u*t+u];
                let b = tt[u*t+u+1];
                let c = tt[(u+1)*t+u];
                let d = tt[(u+1)*t+u+1];
                let scale = a.abs().max(b.abs()).max(c.abs()).max(d.abs());
                if scale <= 0.0 || a <= 0.0 || d <= 0.0
                    || (b/scale-c/scale).abs() > network.config.solve_relative_tolerance
                    || (a/scale)*(d/scale)-(b/scale)*(c/scale) <= 0.0 {
                    return Err(invalid("each Coulomb tangent block needs a finite symmetric positive-definite response"));
                }
                point.response_trace = finite(a+d)?;
                if point.response_trace <= 0.0 { return Err(invalid("Coulomb response trace must be positive")); }
            }
        }
        Ok(Self { points, nt, tn, tt, old_y: vec![0.0; t], free_y: vec![0.0; t],
            reactions: vec![0.0; t], staged: vec![None; p] })
    }

    pub(super) fn law(&self, index: usize) -> Option<&ModalCoulombFriction> {
        self.points.get(index)?.as_ref().map(|p| &p.spec)
    }

    pub(super) fn prepare(&mut self, old: &[f64], free: &[f64], gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        self.reactions.fill(0.0);
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            if let Some(point) = point {
                for axis in 0..2 {
                    self.old_y[2*i+axis] = dot(&point.columns[axis], old)?;
                    self.free_y[2*i+axis] = dot(&point.columns[axis], free)?;
                }
            }
        }
        Ok(())
    }

    pub(super) fn normal_endpoint(&self, i: usize, mut x: f64) -> Result<f64, ModalCouplingError> {
        for j in 0..self.reactions.len() {
            if self.reactions[j] != 0.0 { x = finite(x-self.nt[i*self.reactions.len()+j]*self.reactions[j])?; }
        }
        Ok(x)
    }

    fn delta(&self, i: usize, normals: &[f64], exclude_block: bool) -> Result<[f64; 2], ModalCouplingError> {
        let p = self.points.len();
        let t = self.reactions.len();
        let mut delta = [0.0; 2];
        for (axis, value) in delta.iter_mut().enumerate() {
            let row = 2*i+axis;
            *value = finite(self.free_y[row]-self.old_y[row])?;
            for (j, r) in normals.iter().enumerate() {
                if *r != 0.0 { *value = finite(*value-self.tn[row*p+j]*r)?; }
            }
            for j in 0..t {
                if !(exclude_block && j/2 == i) && self.reactions[j] != 0.0 {
                    *value = finite(*value-self.tt[row*t+j]*self.reactions[j])?;
                }
            }
        }
        Ok(delta)
    }

    pub(super) fn solve_coordinate(&mut self, i: usize, normals: &[f64], dt: f64,
        config: ModalContactConfig, gate: Option<&CancelGate>) -> Result<(), ModalCouplingError>
    {
        let Some(point) = &self.points[i] else { return Ok(()); };
        let capacity = point.capacity(normals[i])?;
        // The work ceiling may bound a TRIAL but never changes the physical
        // graph checked jointly or on actual candidate motion.
        let radius = capacity.min(point.spec.maximum_force_n);
        let free = self.delta(i, normals, true)?;
        let t = self.reactions.len();
        let u = 2*i;
        let mut x = project([self.reactions[u], self.reactions[u+1]], radius)?;
        let mut last = (0.0, 0.0);
        for iteration in 0..=config.max_iterations {
            poll(gate)?;
            let delta = [finite(free[0]-self.tt[u*t+u]*x[0]-self.tt[u*t+u+1]*x[1])?,
                finite(free[1]-self.tt[(u+1)*t+u]*x[0]-self.tt[(u+1)*t+u+1]*x[1])?];
            let graph = point.graph(x, delta, radius, dt, config)?;
            last = (graph.residual, graph.tolerance);
            if graph.residual <= graph.tolerance {
                self.reactions[u] = x[0]; self.reactions[u+1] = x[1];
                return Ok(());
            }
            if iteration < config.max_iterations { x = graph.projected; }
        }
        Err(ModalCouplingError::ContactSolve { residual_n: last.0, tolerance_n: last.1, iterations: config.max_iterations })
    }

    pub(super) fn residuals(&self, normals: &[f64], points: &[ContactPoint], dt: f64,
        worst: &mut (f64, f64), gate: Option<&CancelGate>) -> Result<(), ModalCouplingError>
    {
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            if let Some(point) = point {
                let x = [self.reactions[2*i], self.reactions[2*i+1]];
                let graph = point.graph(x, self.delta(i, normals, false)?, point.capacity(normals[i])?, dt, points[i].config)?;
                if graph.residual/graph.tolerance > worst.0.abs()/worst.1 { *worst = (graph.residual, graph.tolerance); }
            }
        }
        Ok(())
    }

    pub(super) fn add_forces(&self, forces: &mut [f64], gate: Option<&CancelGate>)
        -> Result<(), ModalCouplingError>
    {
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            if let Some(point) = point {
                for axis in 0..2 {
                    let reaction = self.reactions[2*i+axis];
                    if reaction != 0.0 {
                        for (f, b) in forces.iter_mut().zip(&point.columns[axis]) { *f = finite(*f-b*reaction)?; }
                    }
                }
            }
        }
        Ok(())
    }

    pub(super) fn stage(&mut self, network: &CoupledModalSystem, normals: &[f64],
        points: &[ContactPoint], sweeps: usize, gate: Option<&CancelGate>) -> Result<f64, ModalCouplingError>
    {
        let mut total = 0.0;
        for (i, point) in self.points.iter().enumerate() {
            poll(gate)?;
            let Some(point) = point else { continue; };
            let delta = [dot(&point.columns[0], &network.free_delta)?, dot(&point.columns[1], &network.free_delta)?];
            let x = [self.reactions[2*i], self.reactions[2*i+1]];
            let radius = point.capacity(normals[i])?;
            let graph = point.graph(x, delta, radius, network.dt, points[i].config)?;
            if graph.residual > graph.tolerance {
                return Err(ModalCouplingError::ContactSolve {
                    residual_n: graph.residual, tolerance_n: graph.tolerance, iterations: sweeps });
            }
            limit("contact Coulomb force norm", norm(x)?, point.spec.maximum_force_n)?;
            let distance = norm(delta)?;
            let loss = finite(radius*distance)?;
            let work = finite(finite(x[0]*delta[0])?+finite(x[1]*delta[1])?)?;
            let work_residual = finite(work-loss)?;
            let work_tolerance = finite(graph.tolerance*distance + 2.0*radius*network.dt*point.spec.velocity_tolerance_m_s)?;
            if work_residual.abs() > work_tolerance {
                return Err(ModalCouplingError::EnergyBalance { residual_j: work_residual, tolerance_j: work_tolerance });
            }
            total = finite(total+loss)?;
            self.staged[i] = Some(CoulombContactFrame { reaction_n: x, slip_velocity_m_s: graph.velocity,
                capacity_n: radius, graph_residual_n: graph.residual, force_tolerance_n: graph.tolerance,
                regime: graph.regime, dissipation_j: loss, removed_work_j: work, work_residual_j: work_residual,
                work_tolerance_j: work_tolerance });
        }
        Ok(total)
    }

    pub(super) fn publish(&mut self, frames: &mut Vec<Option<CoulombContactFrame>>) {
        std::mem::swap(frames, &mut self.staged);
    }
}

// Euclidean disk, not a componentwise friction pyramid. hypot avoids squaring
// extreme raw inputs; a nonrepresentable norm or arithmetic operation refuses.
fn norm(v: [f64; 2]) -> Result<f64, ModalCouplingError> { finite(v[0].hypot(v[1])) }
fn project(x: [f64; 2], radius: f64) -> Result<[f64; 2], ModalCouplingError> {
    let length = norm(x)?;
    if radius == 0.0 { return Ok([0.0; 2]); }
    if length <= radius { return Ok(x); }
    Ok([finite(radius*(x[0]/length))?, finite(radius*(x[1]/length))?])
}
fn projected_step(x: [f64; 2], delta: [f64; 2], trace: f64, radius: f64)
    -> Result<[f64; 2], ModalCouplingError>
{
    project([finite(x[0]+delta[0]/trace)?, finite(x[1]+delta[1]/trace)?], radius)
}
