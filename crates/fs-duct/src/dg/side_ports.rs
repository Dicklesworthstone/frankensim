//! Point side junctions using the original alternating DG flux and port work.
//!
//! At an interior face the original flux takes pressure from the downstream
//! element and through-flow from upstream. Adding a held inflow `u` to that
//! downstream pressure equation gives `q_hat_right = q_hat_left + u`, with
//! the SAME pressure conjugate on both through-flow faces and the side port.
//! The energy input is therefore exactly `u * integral(p_hat dt)`.
//!
//! This is a zero-volume plane-wave junction. A finite chimney, aperture area,
//! radiation, junction inertance, higher modes or moving pad must be supplied
//! by its own coupled physical owner; none is inferred from the bore area.

use super::*;

/// Maximum attached point side ports, additional to the two boundary ports.
pub const MAX_DG_SIDE_PORTS: usize = 32;

/// A side-flow junction at an EXISTING interior element face.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct DgSidePort {
    /// Unique port identity across both endpoints and every side junction.
    pub id: PortId,
    /// Element whose left face supplies the original alternating pressure
    /// trace. Must be in `1..mesh.elements`; no station is rounded or snapped.
    pub downstream_element: usize,
}

/// Geometry and complete local material state at the original junction face.
#[derive(Debug, Clone, Copy, PartialEq)]
pub struct DgSidePortSource {
    /// Original identity and mesh-face declaration.
    pub port: DgSidePort,
    /// Axial station of that exact face [m].
    pub station_m: f64,
    /// Downstream BORE cross-sectional area [m2], not a side-aperture area.
    pub bore_area_m2: f64,
    /// Original gas or actual evaluation of its primitive field at the face.
    pub gas: GasState,
}

fn invalid(what: &'static str) -> DuctError { DuctError::BadParameter { what } }

impl<M: DgMedium> DgDuct<M> {
    /// Attach independent held-inflow ports at existing interior faces.
    /// Each port is conjugate to the SAME pressure trace used by the original
    /// alternating through-flow flux. A load consumes the opposite flow sign.
    ///
    /// Only input/output columns change. Acoustic and wall operators, physical
    /// coordinates, storage and the original exponential are unchanged. This
    /// may precede or follow either viscothermal attachment. No direct wall
    /// forcing or duplicate inviscid storage is introduced.
    ///
    /// All unknown loads must be solved together through multiport_response;
    /// independent scalar responses would omit their same-step interaction.
    /// An empty list preserves the original owner and identity verbatim.
    ///
    /// # Errors
    /// Invalid/interior-face index, repeated face or PortId, or more than
    /// MAX_DG_SIDE_PORTS total side ports. A boundary endpoint uses its existing
    /// DgBoundary::Port instead. Existing states cannot be reused after changing
    /// the port basis: the original operator identity binds these columns too.
    pub fn with_side_ports(mut self, ports: &[DgSidePort]) -> Result<Self, DuctError> {
        if ports.is_empty() { return Ok(self); }
        if ports.len() > MAX_DG_SIDE_PORTS - self.side_ports.len() {
            return Err(invalid("DG side-port count exceeds MAX_DG_SIDE_PORTS"));
        }
        for (index, port) in ports.iter().enumerate() {
            if port.downstream_element == 0 || port.downstream_element >= self.mesh.elements {
                return Err(invalid("DG side port must name an existing interior element face"));
            }
            if self.ports.iter().any(|(id, _)| *id == port.id)
                || ports[..index].iter().any(|other| other.id == port.id)
            {
                return Err(invalid("DG side-port IDs must be distinct from every declared port"));
            }
            if self.side_ports.iter().chain(&ports[..index]).any(|other|
                other.downstream_element == port.downstream_element)
            {
                return Err(invalid("DG side ports require distinct interior faces"));
            }
        }
        let n_local = self.mesh.order + 1;
        for port in ports {
            let element = &self.elements[port.downstream_element];
            let base = port.downstream_element * n_local;
            let mut column = vec![0.0; self.total_dofs()];
            for i in 0..n_local {
                column[base + i] = element.mp_inv_sqrt[i * n_local];
            }
            self.ports.push((port.id, column));
            self.side_ports.push(*port);
        }
        self.operator_identity = hash_operator(
            self.mesh, self.inlet, self.outlet, &self.elements, &self.l_tilde, &self.ports,
        );
        Ok(self)
    }

    /// Point junction declarations in attachment order, after boundary ports.
    #[must_use]
    pub fn side_ports(&self) -> &[DgSidePort] { &self.side_ports }

}

impl DgDuct {
    /// Read the actual face geometry and evaluate its original material source.
    /// Source/phase refusals propagate without substituting inlet or mean gas.
    /// Unknown IDs refuse; boundary ports retain their existing endpoint APIs.
    pub fn side_port_source(&self, id: PortId) -> Result<DgSidePortSource, GasFieldDuctError> {
        let port = *self.side_ports.iter().find(|port| port.id == id)
            .ok_or_else(|| invalid("DG side-port source references an undeclared side port"))?;
        let element = &self.elements[port.downstream_element];
        let station_m = element.x_left;
        let gas = match &self.gas_field {
            Some(field) => *field.state_at(station_m)?.state(),
            None => self.gas,
        };
        Ok(DgSidePortSource {
            port, station_m, bore_area_m2: element.areas[0], gas,
        })
    }
}

#[cfg(test)]
mod tests;
