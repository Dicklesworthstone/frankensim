//! Omitted modes at fixed friction ports, in the original joint load line.
//!
//! The producer supplies the static compliance G after any support condensation.
//! Contact displacement is b^T q + G F, so a step-held traction gives
//! slip = u dt - b^T delta_q - G (F - F_old). Its work splits exactly into
//! delta(F^T G F / 2) + delta_F^T G delta_F / 2. The second term is the
//! backward-Euler truncation's NUMERICAL loss, never a constitutive loss.
//! Refining the retained basis drives G and this loss to zero. This is a
//! quasi-static tail, not a dynamic replacement for omitted resonances.

use super::*;
use surfaces::geometry::FrictionGeometry;

pub(super) struct FrictionResidual {
    compliance: Vec<f64>,
    accepted_force: Vec<f64>,
}
impl FrictionResidual {
    fn new(compliance: Vec<f64>, count: usize) -> Result<Self, ModalFrictionError> {
        if compliance.len() != count*count || compliance.iter().any(|g| !g.is_finite()) {
            return Err(invalid("friction residual compliance needs one finite entry per port pair").into());
        }
        let mut largest = 0.0_f64;
        for i in 0..count {
            if compliance[i*count+i] < 0.0 {
                return Err(invalid("friction residual compliance must be positive semidefinite").into());
            }
            largest = largest.max(compliance[i*count+i]);
            for j in 0..i {
                if compliance[i*count+j].to_bits() != compliance[j*count+i].to_bits() {
                    return Err(invalid("friction residual compliance must be exactly symmetric").into());
                }
            }
        }
        let (values, _) = fs_la::eigen::jacobi_eigh(&compliance, count);
        if values.iter().any(|&v| !(v >= -1e-12*largest)) {
            return Err(invalid("friction residual compliance must be positive semidefinite").into());
        }
        Ok(Self { compliance, accepted_force: vec![0.0; count] })
    }
    pub(super) fn compliance(&self) -> &[f64] { &self.compliance }
    pub(super) fn displacement(&self, port: usize, force: &[f64]) -> Result<f64, ModalCouplingError> {
        let p = self.accepted_force.len();
        dot(&self.compliance[port*p..(port+1)*p], force)
    }
    pub(super) fn old_displacement(&self, port: usize) -> Result<f64, ModalCouplingError> {
        self.displacement(port, &self.accepted_force)
    }
    pub(super) fn increment(&self, port: usize, force: &[f64]) -> Result<f64, ModalCouplingError> {
        finite(self.displacement(port, force)? - self.old_displacement(port)?)
    }
    fn quadratic(&self, force: &[f64], difference: bool) -> Result<f64, ModalCouplingError> {
        let p = self.accepted_force.len();
        let (mut sum, mut scale) = (0.0, 0.0);
        for i in 0..p {
            let a = force[i] - if difference { self.accepted_force[i] } else { 0.0 };
            for j in 0..p {
                let b = force[j] - if difference { self.accepted_force[j] } else { 0.0 };
                let term = finite(0.5*a*self.compliance[i*p+j]*b)?;
                sum = finite(sum + term)?;
                scale = finite(scale + term.abs())?;
            }
        }
        if sum < -1e-12*scale {
            return Err(invalid("friction residual quadratic form must be nonnegative"));
        }
        Ok(sum.max(0.0))
    }
    pub(super) fn energy_j(&self) -> Result<f64, ModalCouplingError> {
        self.quadratic(&self.accepted_force, false)
    }
    pub(super) fn staged_energy_loss_j(&self, force: &[f64]) -> Result<(f64,f64), ModalCouplingError> {
        Ok((self.quadratic(force, false)?, self.quadratic(force, true)?))
    }
    pub(super) fn publish(&mut self, force: &[f64]) { self.accepted_force.copy_from_slice(force); }
}

impl MultiFrictionModalSystem {
    /// Bind the omitted modes' static response at fixed friction ports, once
    /// before playback. `geometry` is the retained force/displacement map after
    /// the producer's support condensation; laws and port order stay unchanged.
    /// `compliance` is its symmetric positive-semidefinite row-major matrix in
    /// m/N. Every off-diagonal reaction participates in the simultaneous solve.
    ///
    /// The same producer must derive maps and compliance from the omitted
    /// spectrum. This owner stores F G F / 2 and reports dF G dF / 2 separately
    /// as numerical loss. It does not infer compliance from friction or tune
    /// retained eigenfrequencies. Modal and time refinement remain necessary.
    ///
    /// # Errors
    /// Playback/rebinding, moving maps, normal contacts, nonlinear storage or
    /// source-bound interface histories (whose omitted endpoint velocities are
    /// not reconstructed), invalid maps/matrix, budget or cancellation.
    pub fn with_friction_residual_ports(mut self, geometry: Vec<FrictionGeometry>,
        compliance: Vec<f64>, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || self.friction_residual.is_some() {
            return Err(invalid("friction residual compliance binds once before playback").into());
        }
        self.admit_friction_residual()?;
        let residual = FrictionResidual::new(compliance, self.port_count())?;
        if geometry.len() != self.port_count() {
            return Err(invalid("friction residual needs one condensed map per original port").into());
        }
        let ports = geometry.into_iter().zip(&self.points).map(|(map, old)| ModalFrictionPort {
            left: map.left, right: map.right, law: old.port.law.clone(),
        }).collect();
        let (points, compliance) = self.prepare_surface_maps(ports, gate)?;
        poll(Some(gate))?;
        self.points = points;
        self.compliance = compliance;
        self.core.port = self.points[0].port.clone();
        self.core.column = self.points[0].column.clone();
        self.core.compliance = self.compliance[0];
        self.friction_residual = Some(residual);
        Ok(self)
    }
    /// Admitted omitted-mode compliance in the unchanged friction-port order.
    #[must_use]
    pub fn friction_residual_compliance(&self) -> Option<&[f64]> {
        self.friction_residual.as_ref().map(FrictionResidual::compliance)
    }
    pub(super) fn admit_friction_residual(&self) -> Result<(), ModalFrictionError> {
        if self.obstacles.is_some() || self.core.has_nonlinear_potentials()
            || self.endpoint_controls.is_some()
            || self.shear_layers.iter().any(Option::is_some)
            || self.mindlin_contacts.iter().any(Option::is_some)
            || (0..self.port_count()).any(|i| self.interface(i).is_some()) {
            return Err(invalid("friction residual compliance requires fixed dry ports without normal contacts, nonlinear storage or interface histories").into());
        }
        Ok(())
    }
    pub(super) fn friction_compliance(&self, i: usize, j: usize) -> Result<f64, ModalCouplingError> {
        let index = i*self.port_count()+j;
        match &self.friction_residual {
            Some(residual) => finite(self.compliance[index]+residual.compliance[index]),
            None => Ok(self.compliance[index]),
        }
    }
}

