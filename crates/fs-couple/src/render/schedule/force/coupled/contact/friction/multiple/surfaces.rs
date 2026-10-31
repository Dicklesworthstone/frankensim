//! Two-sided friction without replacing the retained mechanical state.
//!
//! Slip is `u dt - (b_left - b_right)^T delta_q`. The same column applies
//! equal and opposite contact tractions, so the surface's inertia, compliance
//! and material loss participate in the existing simultaneous force solve.
//! Prescribed speed is additional actuator motion, not the resolved surface
//! velocity; zero speed makes the contact an internal passive interface.

use super::*;

/// Atomic held-map changes during playback, using this same response compiler.
pub mod geometry;

/// A physical receiving surface in the original complete component basis.
#[derive(Clone, Debug, PartialEq)]
pub struct FrictionSurfaceAttachment {
    /// Existing friction port, in the original drive order.
    pub port: usize,
    /// Displacement along the same signed tangent as that port's left side.
    /// The existing component receives the opposite generalized traction.
    pub surface: ModalAttachment,
}

impl MultiFrictionModalSystem {
    /// Bind resolved surface motion to previously prescribed friction ports.
    ///
    /// The actual receiving components must already exist in the network, with
    /// their geometry-derived displacement maps. Their complete linear and
    /// nonlinear dynamics, material memory and initial state remain unchanged.
    /// This recomputes the full contact compliance matrix, including off-diagonal
    /// reactions between contacts sharing a compliant surface.
    ///
    /// Bind before playback and before adding normal obstacles (which retain
    /// cross-compliances to the friction columns). A port binds at most once.
    /// Every source interface, pressure observer and numerical budget survives.
    /// This does not infer contact geometry, normal load, bow hair or rosin.
    /// Nonzero drive speed remains an explicit external work source.
    pub fn with_surface_attachments(
        mut self,
        attachments: Vec<FrictionSurfaceAttachment>,
        gate: &CancelGate,
    ) -> Result<Self, ModalFrictionError> {
        poll(Some(gate))?;
        if attachments.is_empty() {
            return Ok(self);
        }
        let p = self.points.len();
        if self.samples_rendered() != 0 || self.obstacles.is_some() || attachments.len() > p {
            return Err(invalid("friction surfaces must bind before playback and normal obstacles within the existing port budget").into());
        }
        let mut ports: Vec<_> = self.points.iter().map(|point| point.port.clone()).collect();
        for attachment in attachments {
            poll(Some(gate))?;
            let port = ports
                .get_mut(attachment.port)
                .ok_or_else(|| invalid("friction surface names an absent port"))?;
            if port.right.is_some()
                || attachment.surface.shapes.is_empty()
                || attachment.surface.shapes.iter().any(|v| !v.is_finite())
                || attachment.surface.shapes.iter().all(|v| *v == 0.0)
            {
                return Err(invalid("each friction surface needs one nonzero finite receiver and cannot replace an existing attachment").into());
            }
            port.right = Some(attachment.surface);
        }
        let (points, compliance) = self.prepare_surface_maps(ports, gate)?;
        self.points = points;
        self.compliance = compliance;
        // The one-port fast path must use exactly the same new physical map.
        self.core.port = self.points[0].port.clone();
        self.core.column = self.points[0].column.clone();
        self.core.compliance = self.compliance[0];
        Ok(self)
    }

    // Shared by cold surface binding and live held-map changes. No mechanical,
    // interface or constitutive history is rebuilt by geometric preparation.
    pub(super) fn prepare_surface_maps(&self, ports: Vec<ModalFrictionPort>, gate: &CancelGate)
        -> Result<(Vec<Point>, Vec<f64>), ModalFrictionError>
    {
        if self.friction_residual.is_some() {
            return Err(invalid("friction residual compliance requires its original fixed attachment maps").into());
        }
        let p = self.points.len();
        if ports.len() != p {
            return Err(invalid("surface maps must retain every original friction port").into());
        }
        let mut columns = Vec::with_capacity(p);
        for (port, point) in ports.iter().zip(&self.points) {
            poll(Some(gate))?;
            columns.push(scalar::admit_port(&self.core.network, port, point.config)?);
        }
        // The original p,n,k setup bound still applies: no modes or ports are
        // added. Recompute every cross-term; diagonal-only updates would leave
        // contacts sharing the receiving body mechanically disconnected.
        let mut compliance = vec![0.0; p * p];
        for j in 0..p {
            let response = network_response(&self.core.network, &columns[j], gate)?;
            for i in 0..p {
                compliance[i * p + j] = dot(&columns[i], &response)?;
            }
            if compliance[j * p + j] <= 0.0 {
                return Err(invalid(
                    "each two-sided friction port needs positive network compliance",
                )
                .into());
            }
        }
        poll(Some(gate))?;
        let points = self.points.iter().zip(ports).zip(columns)
            .map(|((old, port), column)| Point { port, column, config: old.config }).collect();
        Ok((points, compliance))
    }
}
