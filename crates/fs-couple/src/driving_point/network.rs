//! Explicit frequency-domain driving point of a generic duct graph. This is
//! the port seam for player-airway/bore composition; no termination is defaulted.
//! It is not a characteristic-line time-domain network realization.

use fs_duct::network::{
    CoupledTermination, DuctNetwork, NetworkDriveResponse, NodeBoundary, PortTermination,
    network_port_impedance,
};
use fs_duct::{DuctError, LossModel};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_ports::{PortKind, PortOrientation, PortSchema, PortValueShape};

/// A selected active port and explicit loads for every other graph port.
#[derive(Debug, Clone)]
pub struct NetworkDrivingPoint {
    network: DuctNetwork,
    active_port: String,
    terminations: Vec<PortTermination>,
    schema: PortSchema,
    coupled_ports: Vec<String>,
}

/// Complex pressure and owner-outward volume flow bound to an explicit schema.
#[derive(Debug, Clone)]
pub struct HarmonicNetworkPort {
    /// Angular frequency [rad/s] bound to these phasors.
    pub omega: f64,
    /// Caller-declared pressure/volume-flow schema and reference frame.
    pub schema: PortSchema,
    /// Pressure phasor [Pa], exp(-i omega t).
    pub pressure: C64,
    /// Volume flow leaving the network [m³/s]; minus the injected flow.
    pub outward_volume_flow: C64,
}

impl NetworkDrivingPoint {
    /// Bind all port labels at construction. Empty/missing/duplicate/unknown
    /// assignments refuse; medium, frequency and loss remain explicit per solve.
    pub fn new(
        network: DuctNetwork,
        active_port: String,
        terminations: Vec<PortTermination>,
        schema: PortSchema,
    ) -> Result<Self, DuctError> {
        Self::new_with_coupled_ports(network, active_port, terminations, Vec::new(), schema)
    }

    /// Declare one coupled boundary's port order, plus independent loads on every
    /// remaining port. The frequency-specific matrix is supplied to
    /// `response_coupled`; it is never replaced with independent diagonal loads.
    /// An empty coupled list has exactly the original constructor's semantics.
    pub fn new_with_coupled_ports(
        network: DuctNetwork,
        active_port: String,
        terminations: Vec<PortTermination>,
        coupled_ports: Vec<String>,
        schema: PortSchema,
    ) -> Result<Self, DuctError> {
        network.validate()?;
        if schema.kind() != PortKind::FluidPressureFlux
            || schema.shape() != PortValueShape::Scalar
            || schema.coordinates().orientation() != PortOrientation::OutwardFromOwner
        {
            return Err(DuctError::BadParameter {
                what: "network driving point needs an explicit scalar pressure/volume-flow schema with owner-outward orientation",
            });
        }
        let ports: Vec<&str> = network
            .nodes
            .iter()
            .filter_map(|n| match &n.boundary {
                NodeBoundary::Port { label } => Some(label.as_str()),
                _ => None,
            })
            .collect();
        let bad = || DuctError::BadParameter {
            what: "network driving point requires an active port and exactly one explicit load for every other port",
        };
        if !ports.contains(&active_port.as_str()) || terminations.len().checked_add(coupled_ports.len()).and_then(|n| n.checked_add(1)) != Some(ports.len()) {
            return Err(bad());
        }
        let mut labels = vec![active_port.as_str()];
        for label in terminations.iter().map(|t| t.label.as_str())
            .chain(coupled_ports.iter().map(String::as_str)) {
            if !ports.contains(&label) || labels.contains(&label) {
                return Err(bad());
            }
            labels.push(label);
        }
        Ok(Self {
            network,
            active_port,
            terminations,
            schema,
            coupled_ports,
        })
    }

    /// Complete pressure/volume-flow response for unit flow into the active port.
    pub fn response(
        &self,
        gas: &GasState,
        omega: f64,
        loss: LossModel,
    ) -> Result<NetworkDriveResponse, DuctError> {
        if !self.coupled_ports.is_empty() {
            return Err(DuctError::BadParameter {
                what: "declared coupled ports require their complete boundary matrix",
            });
        }
        network_port_impedance(&self.network, gas, omega, loss)?.drive(
            &self.active_port,
            &self.terminations,
            gas,
        )
    }

    /// Driving-point impedance [Pa s/m³], exp(-i omega t).
    pub fn impedance(&self, gas: &GasState, omega: f64, loss: LossModel) -> Result<C64, DuctError> {
        Ok(self.response(gas, omega, loss)?.impedance)
    }

    /// Declared graph geometry and port labels.
    #[must_use]
    pub fn network(&self) -> &DuctNetwork {
        &self.network
    }

    /// Actual harmonic transfer, with the graph's inward current pulled back
    /// to the schema's owner-outward flow. No window-energy audit is implied.
    pub fn harmonic_port(
        &self,
        gas: &GasState,
        omega: f64,
        loss: LossModel,
        injected_volume_flow: C64,
    ) -> Result<HarmonicNetworkPort, DuctError> {
        if !injected_volume_flow.re.is_finite() || !injected_volume_flow.im.is_finite() {
            return Err(DuctError::BadParameter {
                what: "network driving flow phasor must be finite",
            });
        }
        let pressure = self.impedance(gas, omega, loss)? * injected_volume_flow;
        if !pressure.re.is_finite() || !pressure.im.is_finite() {
            return Err(DuctError::BadParameter {
                what: "network harmonic pressure overflowed",
            });
        }
        Ok(HarmonicNetworkPort {
            omega,
            schema: self.schema.clone(),
            pressure,
            outward_volume_flow: -injected_volume_flow,
        })
    }
}

impl NetworkDrivingPoint {
    /// Original complete matrix order, not the network node order. A producer
    /// must use these labels in this order; transposition/permutation is never
    /// guessed from dimensions. The independent loads are disjoint by admission.
    #[must_use]
    pub fn coupled_ports(&self) -> &[String] { &self.coupled_ports }

    /// Solve one complete physical boundary. The original duct owner checks the
    /// exact frequency, source identity, finite entries and full Hermitian
    /// passivity; this adapter neither symmetrizes nor clips its matrix.
    pub fn response_coupled(
        &self, gas: &GasState, omega: f64, loss: LossModel,
        boundary: &CoupledTermination,
    ) -> Result<NetworkDriveResponse, DuctError> {
        if self.coupled_ports.is_empty() || boundary.labels != self.coupled_ports {
            return Err(DuctError::BadParameter {
                what: "coupled boundary must use the declared nonempty port order",
            });
        }
        network_port_impedance(&self.network, gas, omega, loss)?.drive_coupled(
            &self.active_port, &self.terminations, boundary, gas,
        )
    }

    /// The same coupled solution at the original dimensioned active-port seam.
    /// No second boundary solve or independent-port load is substituted.
    pub fn harmonic_port_coupled(
        &self, gas: &GasState, omega: f64, loss: LossModel,
        boundary: &CoupledTermination, injected_volume_flow: C64,
    ) -> Result<HarmonicNetworkPort, DuctError> {
        if !injected_volume_flow.re.is_finite() || !injected_volume_flow.im.is_finite() {
            return Err(DuctError::BadParameter { what: "network driving flow phasor must be finite" });
        }
        let pressure = self.response_coupled(gas, omega, loss, boundary)?.impedance * injected_volume_flow;
        if !pressure.re.is_finite() || !pressure.im.is_finite() {
            return Err(DuctError::BadParameter { what: "network harmonic pressure overflowed" });
        }
        Ok(HarmonicNetworkPort { omega, schema: self.schema.clone(), pressure,
            outward_volume_flow: -injected_volume_flow })
    }
}

#[cfg(test)]
mod tests;

/// Complete exterior opening feedback and coherent microphone fields.
pub mod radiation;
