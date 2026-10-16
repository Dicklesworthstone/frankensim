//! Fixed-activity IFT derivatives of a settled mechanical/contact network.
//!
//! R(q,p)=K0*q + sum b_i*r_i(b_i^T*q-gap_i) - f = 0 (rest loads included).
//! The tangent is K0 + sum k_i*b_i*b_i^T. Reuse the SAME static inverse K0^-1
//! and fs-la Cholesky for its contact-space Woodbury correction; do not
//! differentiate nonlinear sweeps, time integration, or a fictitious warm-up.
//! One transpose solve serves any scalar objective's parameter pullback.
//!
//! This is local sensitivity of the admitted reduced static equations, not a
//! derivative across activity switches, mode-basis recomputation, friction,
//! acoustic radiation, dynamic trajectories or uncertain material identification.
use super::*;
use super::super::contact::{ModalContact, ModalContactConfig, contact_column, force_tolerance};
use fs_dcontact::{OpeningContactStep, ContactStorage};
use fs_phs::Storage;

/// Explicit derivative preparation/query work screens, not wall-clock promises.
#[derive(Clone, Copy, Debug)]
pub struct SensitivityBudget {
    /// Maximum normal contacts, also hard-capped at 32; zero is legal.
    pub max_contacts: usize,
    /// Screen n*(k+p+1)^2+(k+p+1)^3; K0 preparation has its original caps too.
    pub max_setup_terms: usize,
    /// Per solve/pullback screen n*(k+p+1)+(k+p+1)^2.
    pub max_query_terms: usize,
    /// Refuse when |closure-gap| <= this finite nonnegative distance [m].
    pub minimum_contact_margin_m: f64,
}

/// Evidence of the actual force-balance check at the frozen initial state.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct EquilibriumLinearizationReport {
    /// Maximum |R_j| / (original relative allowance + mapped contact allowance).
    pub primal_residual_to_tolerance: f64,
    /// Smallest distance from any activity switch [m], absent without contacts.
    pub minimum_contact_margin_m: Option<f64>,
    /// Number of strictly penetrating contacts (not inferred from force alone).
    pub active_contacts: usize,
}

/// Original tangent solution; the physical network is never changed by a query.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumTangentSolution {
    /// Component-major, mode-major solution of the requested tangent system.
    pub values: Vec<f64>,
    /// Recomputed full-coordinate residual normalized by absolute term sums.
    pub relative_residual: f64,
}

/// Parameter partials for one declared bilateral spring, in construction order.
#[derive(Clone, Debug, PartialEq)]
pub struct SpringParameterPullback {
    /// Partial with respect to spring stiffness, not its dashpot coefficient.
    pub stiffness: f64,
    /// Partial with respect to the stress-free relative displacement [m].
    pub rest_extension: f64,
    /// Partial with respect to B=left-right in the full retained modal layout.
    /// Left attachment entries take these values; right entries take negatives.
    pub column: Vec<f64>,
}

/// Parameter partials for one declared normal contact, in construction order.
#[derive(Clone, Debug, PartialEq)]
pub struct ContactParameterPullback {
    /// Partial with respect to the original power-law stiffness coefficient.
    pub stiffness: f64,
    /// Partial with respect to the original gap [m].
    pub gap: f64,
    /// Partial with respect to the supplied quadrature weight.
    pub weight: f64,
    /// Partial with respect to B=left-right, including force AND closure mapping.
    pub column: Vec<f64>,
}

/// (dR/dp)^T lambda, NOT the total objective gradient. For a goal J use
/// dJ/dp = explicit_dJ/dp - this value after solving K^T lambda = dJ/dq.
/// Damping has zero static partial and is not identifiable from these goals.
#[derive(Clone, Debug, PartialEq)]
pub struct EquilibriumParameterPullback {
    /// Generalized held force partials; equal to -lambda.
    pub external_forces: Vec<f64>,
    /// Natural angular frequency partials. Explicit free coordinates are None:
    /// a physical mass derivative must include its coordinate/shape chain rule.
    pub angular_frequencies: Vec<Option<f64>>,
    /// One entry per bilateral connection.
    pub springs: Vec<SpringParameterPullback>,
    /// One entry per normal contact.
    pub contacts: Vec<ContactParameterPullback>,
}

struct Point {
    column: Vec<f64>,
    force: f64,
    tangent: f64,
    force_per_stiffness: f64,
    force_per_weight: f64,
}

/// Immutable borrow of an ACTUAL equilibrium. Keeping it alive prevents source
/// state/parameter mutation through this API. Every query reuses one prepared
/// contact factor; no n-by-n modal matrix and no per-parameter primal solve.
pub struct EquilibriumLinearization<'a> {
    network: &'a CoupledModalSystem,
    base: StaticResponse,
    q: Vec<f64>,
    points: Vec<Point>,
    responses: Vec<Vec<f64>>,
    roots: Vec<f64>,
    matrix: Vec<f64>,
    factor: Cholesky,
    budget: SensitivityBudget,
    report: EquilibriumLinearizationReport,
}

impl<'a> EquilibriumLinearization<'a> {
    /// Prepare at sample zero with exactly zero velocity, checking the complete
    /// original force balance against `external` and the supplied contact set.
    /// Supported free coordinates are legal under K0's existing support screen.
    /// Missing forces/contacts, moving states, switching margins, nonfinite data,
    /// exceeded caps and cancellation all refuse without mutating the network.
    pub fn new(network: &'a CoupledModalSystem, external: &[f64],
        contacts: &[(ModalContact, ModalContactConfig)], budget: SensitivityBudget, gate: &CancelGate)
        -> Result<Self, ModalCouplingError>
    {
        poll(Some(gate))?;
        let (n,p,k)=(network.mode_count(),contacts.len(),network.columns.len());
        if budget.max_contacts>32 || p>budget.max_contacts || !budget.minimum_contact_margin_m.is_finite()
            || budget.minimum_contact_margin_m<0.0 || network.samples_rendered()!=0
            || external.len()!=n || external.iter().any(|x|!x.is_finite())
            || network.models.iter().flat_map(|m|m.states()).any(|s|s.velocity_m_sqrt_kg_per_s!=0.0) {
            return Err(invalid("static sensitivity requires sample-zero stationary states, complete forces and bounded contacts/margins"));
        }
        let width=k+p+1;
        let setup=n.checked_mul(width*width).and_then(|v|v.checked_add(width*width*width))
            .ok_or_else(||invalid("static sensitivity setup overflow"))?;
        if setup>budget.max_setup_terms {return Err(invalid("static sensitivity setup exceeds max_setup_terms"));}
        let q:Vec<f64>=network.models.iter().flat_map(|m|m.states()).map(|s|s.displacement_m_sqrt_kg).collect();
        let mut residual:Vec<f64>=external.iter().map(|f|-f).collect();
        let mut scales:Vec<f64>=external.iter().map(|f|f.abs()).collect();
        let mut allowance=vec![0.0;n];
        for (j,mode) in network.models.iter().flat_map(|m|m.modes()).enumerate() {
            add_term(&mut residual[j],&mut scales[j],finite(mode.angular_frequency_rad_s*mode.angular_frequency_rad_s*q[j])?)?;
        }
        for (column,link) in network.columns.iter().zip(&network.connections) {
            poll(Some(gate))?;
            let force=finite(link.stiffness_n_m*finite(dot(column,&q)?-link.rest_extension_m)?)?;
            limit("static sensitivity spring force",force.abs(),network.config.maximum_abs_connection_force_n)?;
            for j in 0..n {add_term(&mut residual[j],&mut scales[j],finite(column[j]*force)?)?;}
        }
        let mut points=Vec::with_capacity(p);
        let mut report=EquilibriumLinearizationReport {primal_residual_to_tolerance:0.0,minimum_contact_margin_m:None,active_contacts:0};
        let mut energy=network.total_energy_j()?;
        for (contact,c) in contacts {
            poll(Some(gate))?;
            let column=contact_column(network,contact,*c)?;
            let x=dot(&column,&q)?;
            let d=OpeningContactStep::new(&contact.law,-x).and_then(|s|s.static_differential(budget.minimum_contact_margin_m))
                .map_err(ModalCouplingError::ContactLaw)?;
            limit("static sensitivity contact force",d.force_n,c.maximum_force_n)?;
            limit("static sensitivity penetration",d.signed_penetration_m.max(0.0),c.maximum_penetration_m)?;
            let margin=d.signed_penetration_m.abs();
            report.minimum_contact_margin_m=Some(report.minimum_contact_margin_m.map_or(margin,|m|m.min(margin)));
            report.active_contacts+=usize::from(d.signed_penetration_m>0.0);
            let storage=ContactStorage::new(Box::new(ZeroStorage),1,vec![contact.law.clone()])
                .map_err(ModalCouplingError::ContactLaw)?;
            energy=finite(energy+finite(storage.probe(&[-x,0.0]).contact_energy)?)?;
            let force_allowance=force_tolerance(d.force_n,d.force_n,*c)?;
            for j in 0..n {
                add_term(&mut residual[j],&mut scales[j],finite(column[j]*d.force_n)?)?;
                allowance[j]=finite(allowance[j]+column[j].abs()*force_allowance)?;
            }
            points.push(Point {column,force:d.force_n,tangent:d.closure_stiffness_n_m,
                force_per_stiffness:d.force_per_stiffness,force_per_weight:d.force_per_weight});
        }
        limit("static sensitivity stored energy",energy,network.config.maximum_total_energy_j)?;
        for j in 0..n {
            let tolerance=finite(network.config.solve_relative_tolerance*scales[j]+allowance[j])?;
            if residual[j].abs()>tolerance {return Err(invalid("static sensitivity state does not satisfy the declared full force balance"));}
            if tolerance>0.0 {report.primal_residual_to_tolerance=report.primal_residual_to_tolerance.max(residual[j].abs()/tolerance);}
        }
        let base=StaticResponse::new(network,gate)?;
        let mut responses=Vec::with_capacity(p);
        for point in &points {responses.push(base.solve(network,&point.column,false,false,gate)?);}
        let roots:Vec<f64>=points.iter().map(|p|p.tangent.sqrt()).collect();
        let mut matrix=vec![0.0;p*p];
        for i in 0..p {
            poll(Some(gate))?;
            for j in 0..=i {
                let a=finite(roots[i]*dot(&points[i].column,&responses[j])?*roots[j])?;
                let b=finite(roots[j]*dot(&points[j].column,&responses[i])?*roots[i])?;
                if (a-b).abs()>network.config.solve_relative_tolerance*a.abs().max(b.abs()).max(1.0) {
                    return Err(invalid("static contact tangent lost numerical symmetry"));
                }
                let value=finite(f64::midpoint(a,b)+if i==j {1.0}else{0.0})?;
                matrix[i*p+j]=value;matrix[j*p+i]=value;
            }
        }
        let factor=cholesky(&matrix,p).map_err(ModalCouplingError::Factor)?;
        poll(Some(gate))?;
        Ok(Self {network,base,q,points,responses,roots,matrix,factor,budget,report})
    }

    /// Actual equilibrium admission and activity separation, not a certificate.
    #[must_use]
    pub const fn report(&self)->EquilibriumLinearizationReport {self.report}
    /// Frozen displacement in the original retained basis.
    #[must_use]
    pub fn displacement(&self)->&[f64] {&self.q}

    /// Tangent action; symmetric, so also its transpose. No state is advanced.
    pub fn apply(&self, direction:&[f64], gate:&CancelGate)->Result<Vec<f64>,ModalCouplingError> {
        self.admit_query(direction,gate)?;
        Ok(self.action(direction,gate)?.0)
    }

    /// Solve a load tangent OR a scalar-goal adjoint using the same symmetric
    /// tangent. Rest extensions belong only in the primal, never this inverse.
    /// A small contact-space residual is insufficient: recompute all modal rows.
    pub fn solve(&self, rhs:&[f64], gate:&CancelGate)->Result<EquilibriumTangentSolution,ModalCouplingError> {
        self.admit_query(rhs,gate)?;
        let mut values=self.base.solve(self.network,rhs,false,false,gate)?;
        let mut reduced:Vec<f64>=self.points.iter().zip(&self.roots)
            .map(|(p,r)|finite(r*dot(&p.column,&values)?)).collect::<Result<_,_>>()?;
        let original=reduced.clone();self.factor.solve(&mut reduced);
        check_solve(&self.matrix,&reduced,&original,self.network.config.solve_relative_tolerance)?;
        for (i,response) in self.responses.iter().enumerate() {
            poll(Some(gate))?;
            for (x,y) in values.iter_mut().zip(response) {*x=finite(*x-y*self.roots[i]*reduced[i])?;}
        }
        let (applied,scales)=self.action(&values,gate)?;
        let mut relative=0.0_f64;
        for j in 0..rhs.len() {
            let denominator=finite(scales[j]+rhs[j].abs())?;
            let residual=finite(applied[j]-rhs[j])?.abs();
            if denominator>0.0 {relative=relative.max(residual/denominator);}
            else if residual!=0.0 {return Err(invalid("static tangent residual has zero scale"));}
        }
        if relative>self.network.config.solve_relative_tolerance {
            return Err(ModalCouplingError::SolveResidual {relative,tolerance:self.network.config.solve_relative_tolerance});
        }
        poll(Some(gate))?;
        Ok(EquilibriumTangentSolution {values,relative_residual:relative})
    }

    /// Analytic (dR/dp)^T lambda for original modal frequencies, forces, spring
    /// parameters and contact parameters/shapes. Both left/right shape effects
    /// are retained through the signed conjugate column; no finite differences.
    pub fn parameter_pullback(&self, lambda:&[f64], gate:&CancelGate)
        ->Result<EquilibriumParameterPullback,ModalCouplingError>
    {
        self.admit_query(lambda,gate)?;
        let angular_frequencies=self.network.models.iter().flat_map(|m|m.modes()).enumerate()
            .map(|(j,m)|if m.angular_frequency_rad_s==0.0 {Ok(None)} else {
                finite(2.0*m.angular_frequency_rad_s*self.q[j]*lambda[j]).map(Some)
            }).collect::<Result<Vec<_>,_>>()?;
        let mut springs=Vec::with_capacity(self.network.columns.len());
        for (b,link) in self.network.columns.iter().zip(&self.network.connections) {
            poll(Some(gate))?;
            let s=dot(b,lambda)?;let x=finite(dot(b,&self.q)?-link.rest_extension_m)?;
            let column=self.q.iter().zip(lambda).map(|(q,l)|finite(link.stiffness_n_m*finite(x*l+s*q)?))
                .collect::<Result<_,_>>()?;
            springs.push(SpringParameterPullback {stiffness:finite(x*s)?,rest_extension:finite(-link.stiffness_n_m*s)?,column});
        }
        let mut contacts=Vec::with_capacity(self.points.len());
        for point in &self.points {
            poll(Some(gate))?;
            let s=dot(&point.column,lambda)?;
            let column=self.q.iter().zip(lambda).map(|(q,l)|finite(point.force*l+point.tangent*s*q))
                .collect::<Result<_,_>>()?;
            contacts.push(ContactParameterPullback {stiffness:finite(point.force_per_stiffness*s)?,gap:finite(-point.tangent*s)?,
                weight:finite(point.force_per_weight*s)?,column});
        }
        poll(Some(gate))?;
        Ok(EquilibriumParameterPullback {external_forces:lambda.iter().map(|l|-l).collect(),angular_frequencies,springs,contacts})
    }

    fn admit_query(&self,x:&[f64],gate:&CancelGate)->Result<(),ModalCouplingError> {
        poll(Some(gate))?;
        if x.len()!=self.q.len() || x.iter().any(|v|!v.is_finite()) {return Err(invalid("static derivative vector must match every finite modal coordinate"));}
        let width=self.network.columns.len()+self.points.len()+1;
        let terms=self.q.len().checked_mul(width).and_then(|v|v.checked_add(width*width))
            .ok_or_else(||invalid("static derivative query work overflow"))?;
        if terms>self.budget.max_query_terms {return Err(invalid("static derivative query exceeds max_query_terms"));}
        Ok(())
    }
    fn action(&self,x:&[f64],gate:&CancelGate)->Result<(Vec<f64>,Vec<f64>),ModalCouplingError> {
        let mut out=vec![0.0;x.len()];let mut scale=vec![0.0;x.len()];
        for (j,m) in self.network.models.iter().flat_map(|m|m.modes()).enumerate() {
            add_term(&mut out[j],&mut scale[j],finite(m.angular_frequency_rad_s*m.angular_frequency_rad_s*x[j])?)?;
        }
        for (b,k) in self.network.columns.iter().zip(&self.network.connections).map(|(b,l)|(b,l.stiffness_n_m))
            .chain(self.points.iter().map(|p|(&p.column,p.tangent))) {
            poll(Some(gate))?;
            let force=finite(k*dot(b,x)?)?;
            for j in 0..x.len() {add_term(&mut out[j],&mut scale[j],finite(b[j]*force)?)?;}
        }
        Ok((out,scale))
    }
}
fn add_term(value:&mut f64,scale:&mut f64,term:f64)->Result<(),ModalCouplingError> {
    *value=finite(*value+term)?;*scale=finite(*scale+term.abs())?;Ok(())
}

struct ZeroStorage;
impl Storage for ZeroStorage {
    fn hamiltonian(&self,_:&[f64])->f64 {0.0}
    fn gradient(&self,_:&[f64],out:&mut [f64]) {out.fill(0.0);}
}

/// Physical displacement objectives and correctly projected actuator gradients.
pub mod objective;
