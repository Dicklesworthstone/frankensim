//! Explicit frequency-domain driving point of a generic duct graph. This is
//! the port seam for player-airway/bore composition; no termination is defaulted.
//! It is not a characteristic-line time-domain network realization.

use fs_duct::network::{
    DuctNetwork, NetworkDriveResponse, NodeBoundary, PortTermination, network_port_impedance,
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
        if !ports.contains(&active_port.as_str()) || terminations.len() + 1 != ports.len() {
            return Err(bad());
        }
        let mut labels = vec![active_port.as_str()];
        for t in &terminations {
            if !ports.contains(&t.label.as_str()) || labels.contains(&t.label.as_str()) {
                return Err(bad());
            }
            labels.push(t.label.as_str());
        }
        Ok(Self {
            network,
            active_port,
            terminations,
            schema,
        })
    }

    /// Complete pressure/volume-flow response for unit flow into the active port.
    pub fn response(
        &self,
        gas: &GasState,
        omega: f64,
        loss: LossModel,
    ) -> Result<NetworkDriveResponse, DuctError> {
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
