//! Omitted modes at fixed friction ports, in the original joint load line.
//!
//! The producer supplies the static compliance G after any support condensation.
//! Contact displacement is b^T q + G F, so a step-held traction gives
//! slip = u dt - b^T delta_q - G (F - F_old). Its work splits exactly into
//! delta(F^T G F / 2) + delta_F^T G delta_F / 2. The second term is the
//! backward-Euler truncation's NUMERICAL loss, never a constitutive loss.
//! Refining the retained basis drives G and this loss to zero. This is a
//! quasi-static tail, not a dynamic replacement for omitted resonances.
//!
//! A material film uses that same corrected slip, including every cross-port
//! reaction, for its elastic/plastic/Maxwell return map. Its storage and heat
//! remain separate from the tail: G is in series with the physical contact,
//! not an added film modulus or plastic viscosity. Finite-cell thermal trials
//! repeat this joint solve from the accepted film AND residual-force history.

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

// Sparse signed support rows enter the SAME complete force/displacement column
// for tangential and normal contacts. One row per component, in original index
// order, bounds both validation and accumulation by the retained mode count.
pub(super) fn add_support_corrections(network: &CoupledModalSystem, column: &mut [f64],
    corrections: &[ModalAttachment], gate: &CancelGate) -> Result<(), ModalCouplingError>
{
    if corrections.len() > network.models.len() {
        return Err(invalid("residual support maps exceed the retained component count"));
    }
    let mut previous = None;
    for correction in corrections {
        poll(Some(gate))?;
        if previous.is_some_and(|component| component >= correction.component) {
            return Err(invalid("residual support maps require unique components in ascending order"));
        }
        let model = network.models.get(correction.component)
            .ok_or_else(|| invalid("residual support map names an unknown retained component"))?;
        if correction.shapes.len() != model.modes().len()
            || correction.shapes.iter().any(|shape| !shape.is_finite()) {
            return Err(invalid("residual support map must match its finite retained basis"));
        }
        let offset = network.offsets[correction.component];
        for (local, shape) in correction.shapes.iter().enumerate() {
            column[offset+local] = finite(column[offset+local]+shape)?;
        }
        previous = Some(correction.component);
    }
    Ok(())
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
    /// Fixed material shear films and series compliant drives share the actual
    /// slip increment on this load line. Film plastic/relaxation loss may enter
    /// finite heat cells; the residual's numerical loss never does.
    ///
    /// # Errors
    /// Playback/rebinding, moving maps, normal contacts, nonlinear storage or
    /// Mindlin or dry-interface source histories (whose omitted endpoint
    /// velocities are not reconstructed), invalid maps/matrix, budget or
    /// cancellation. Fixed material films use interval slip and are admitted.
    pub fn with_friction_residual_ports(self, geometry: Vec<FrictionGeometry>,
        compliance: Vec<f64>, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        self.bind_friction_residual_ports(geometry, None, compliance, gate)
    }

    /// Bind the same omitted-mode tail while retaining both physical contact
    /// surfaces and a sparse collection of signed support maps per port. A condensed
    /// string against resolved fibers, for example, still applies the opposite
    /// traction to the actual fiber while its bridge reacts on the body.
    ///
    /// The complete retained map is `b_left - b_right + sum(b_correction)`.
    /// It supplies BOTH displacement and force in the original joint solver;
    /// every cross-port term is recomputed from that complete map. The producer
    /// must derive this correction and `compliance` from the SAME support Schur
    /// complement. A correction has no independent force law, state or heat.
    /// The local contact attachments returned by `port` retain their original
    /// participant identities. All fixed-map restrictions above still apply.
    /// Supply one collection per original port, including an empty collection
    /// when it has no support correction. Within each collection, original
    /// component indices must be unique and ascending; combine same-body rows
    /// before binding. Thus support setup stays bounded by the original modes.
    pub fn with_friction_residual_ports_and_corrections(self,
        geometry: Vec<FrictionGeometry>, corrections: Vec<Vec<ModalAttachment>>,
        compliance: Vec<f64>, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        self.bind_friction_residual_ports(geometry, Some(corrections), compliance, gate)
    }

    /// Bind the same omitted-mode tail at the friction ports beside already bound
    /// FIXED normal contacts that carry their own residual compliance
    /// (`with_obstacle_residual_ports_and_corrections`), with the omitted modes'
    /// contact/port block. A finger stop on a bowed string, for example, deflects
    /// the bow point through the modes the string does not retain, and the bow
    /// deflects the stop's closure.
    ///
    /// `obstacle_cross` is row-major q x p [m/N]: contact k's closure under a unit
    /// traction of port j. The producer derives it, `compliance` and the contacts'
    /// own block from ONE omitted spectrum. With f = [t; -R], the joint matrix
    /// [[G_pp, G_sp^T], [G_sp, G_ss]] must be symmetric positive semidefinite; it is
    /// checked here. Closure becomes N^T q - G_ss R + G_sp t and port slip gains
    /// G_sp^T dR. Both enter the same joint solve, its active-set rescue on the
    /// condensed stops, and its ledger: the contacts book the cross part
    /// -R^T G_sp t of the joint storage f^T G f / 2 and -dR^T G_sp dt of its
    /// declared numerical loss. Moving contacts (normal translation or attachment
    /// motion) and contact-loaded ports refuse, in either binding order.
    ///
    /// # Errors
    /// As `with_friction_residual_ports_and_corrections`; no bound contacts,
    /// contacts without their own residual or moving ones, a contact-loaded port,
    /// a block of the wrong size, or a joint matrix that is not PSD.
    pub fn with_friction_residual_ports_beside_obstacles(self,
        geometry: Vec<FrictionGeometry>, corrections: Vec<Vec<ModalAttachment>>,
        compliance: Vec<f64>, obstacle_cross: Vec<f64>, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        let (p, q) = (self.port_count(), self.obstacle_count());
        let stops = self.obstacles.as_ref()
            .ok_or_else(|| invalid("a joint stop/port residual needs bound normal contacts"))?;
        let stop_block = self.obstacle_residual_compliance()
            .ok_or_else(|| invalid("a joint stop/port residual needs the contacts' own residual compliance"))?;
        if !stops.fixed_residual() || stop_block.len() != q*q || obstacle_cross.len() != q*p
            || compliance.len() != p*p {
            return Err(invalid("a joint stop/port residual needs fixed contacts and one entry per contact and port").into());
        }
        // The joint matrix over [ports; contacts], f = [t; -R].
        let s = p+q;
        let mut joint = vec![0.0; s*s];
        for i in 0..s {
            for j in 0..s {
                joint[i*s+j] = match (i < p, j < p) {
                    (true, true) => compliance[i*p+j],
                    (true, false) => obstacle_cross[(j-p)*p+i],
                    (false, true) => obstacle_cross[(i-p)*p+j],
                    (false, false) => stop_block[(i-p)*q+(j-p)],
                };
            }
        }
        if joint.iter().any(|g| !g.is_finite()) {
            return Err(invalid("a joint stop/port residual needs finite entries").into());
        }
        let largest = (0..s).map(|i| joint[i*s+i]).fold(0.0_f64, f64::max);
        let (values, _) = fs_la::eigen::jacobi_eigh(&joint, s);
        if values.iter().any(|&v| !(v >= -1e-12*largest)) {
            return Err(invalid("a joint stop/port residual compliance must be positive semidefinite").into());
        }
        self.bind_ports(geometry, Some(corrections), compliance, Some(obstacle_cross), gate)
    }

    fn bind_friction_residual_ports(self, geometry: Vec<FrictionGeometry>,
        corrections: Option<Vec<Vec<ModalAttachment>>>, compliance: Vec<f64>, gate: &CancelGate)
        -> Result<Self, ModalFrictionError>
    {
        self.bind_ports(geometry, corrections, compliance, None, gate)
    }

    // With bound contacts, their stop/port block is bound on the ports' new maps
    // in the same transaction; without contacts it must be absent.
    fn bind_ports(mut self, geometry: Vec<FrictionGeometry>,
        corrections: Option<Vec<Vec<ModalAttachment>>>, compliance: Vec<f64>,
        obstacle_cross: Option<Vec<f64>>, gate: &CancelGate) -> Result<Self, ModalFrictionError>
    {
        poll(Some(gate))?;
        if self.samples_rendered() != 0 || self.friction_residual.is_some() {
            return Err(invalid("friction residual compliance binds once before playback").into());
        }
        if self.obstacles.is_some() != obstacle_cross.is_some() {
            return Err(invalid("friction residual compliance beside normal contacts binds with their joint cross-compliance").into());
        }
        self.admit_friction_residual()?;
        let residual = FrictionResidual::new(compliance, self.port_count())?;
        if geometry.len() != self.port_count() {
            return Err(invalid("friction residual needs one condensed map per original port").into());
        }
        let ports = geometry.into_iter().zip(&self.points).map(|(map, old)| ModalFrictionPort {
            left: map.left, right: map.right, law: old.port.law.clone(),
        }).collect();
        let (points, compliance, endpoint) =
            self.prepare_surface_maps_with_corrections(ports, corrections.as_deref(), gate)?;
        poll(Some(gate))?;
        if let (Some(stops), Some(cross)) = (self.obstacles.as_mut(), obstacle_cross) {
            stops.bind_port_cross(&self.core.network, &points, cross, gate)?;
        }
        self.points = points;
        self.compliance = compliance;
        self.endpoint = endpoint;
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
        if self.core.has_nonlinear_potentials()
            || self.endpoint_controls.is_some()
            || self.mindlin_contacts.iter().any(Option::is_some)
            || (0..self.port_count()).any(|i| self.interface(i).is_some()) {
            return Err(invalid("friction residual compliance requires fixed ports without nonlinear storage, Mindlin or dry-interface histories").into());
        }
        // Normal contacts only as fixed ones carrying their own residual, joined
        // to the ports by the stop/port block (`with_friction_residual_ports_beside_obstacles`).
        if let Some(stops) = &self.obstacles {
            let joined = self.friction_residual.is_none() || stops.has_port_cross();
            if !stops.fixed_residual() || !joined || self.contact_loads.iter().any(Option::is_some) {
                return Err(invalid("friction residual compliance beside normal contacts requires fixed contacts with their own residual compliance, their joint cross-compliance, and no contact-loaded port").into());
            }
        }
        Ok(())
    }
    // Port i's row of its own line, with the omitted modes when they are bound.
    pub(super) fn friction_compliance(&self, i: usize, j: usize) -> Result<f64, ModalCouplingError> {
        let line = self.line(i)[j];
        match &self.friction_residual {
            Some(residual) => finite(line+residual.compliance[i*self.port_count()+j]),
            None => Ok(line),
        }
    }
}

