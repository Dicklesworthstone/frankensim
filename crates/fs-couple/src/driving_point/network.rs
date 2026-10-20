//! Explicit frequency-domain driving point of a generic duct graph. This is
//! the port seam for player-airway/bore composition; no termination is defaulted.
//! [`NetworkDrivingPoint::realize_line`] realizes it as a characteristic line
//! on the same ladder as a single duct.

use super::{
    DrivingPointError, LineLadder, LineRealization, bulk_delay_pole_state,
    realize_characteristic_response,
};
use fs_duct::network::{
    CoupledTermination, DuctNetwork, NetworkDriveResponse, NodeBoundary, PortTermination,
    network_port_impedance,
};
use fs_duct::{DuctError, LossModel, Segment, Termination};
use fs_material::gas::GasState;
use fs_math::c64::C64;
use fs_ports::{PortKind, PortOrientation, PortSchema, PortValueShape};
use fs_vfit::discretize::{DelayedFilter, reflectance};

/// A selected active port and explicit loads for every other graph port.
#[derive(Debug, Clone)]
pub struct NetworkDrivingPoint {
    network: DuctNetwork,
    active_port: String,
    terminations: Vec<PortTermination>,
    coupled_labels: Vec<String>,
    schema: PortSchema,
}

/// Per-frequency source of the coupled load, for example a multi-aperture
/// radiation matrix evaluated at each sampled frequency.
pub type CoupledSource<'a> = &'a mut dyn FnMut(f64) -> Result<CoupledTermination, DuctError>;

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
        Self::new_coupled(network, active_port, terminations, Vec::new(), schema)
    }

    /// [`Self::new`] with the ports in `coupled_labels` terminated together by
    /// one [`CoupledTermination`] per solve; `terminations` covers every other
    /// non-active port.
    pub fn new_coupled(
        network: DuctNetwork,
        active_port: String,
        terminations: Vec<PortTermination>,
        coupled_labels: Vec<String>,
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
        if !ports.contains(&active_port.as_str())
            || terminations.len() + coupled_labels.len() + 1 != ports.len()
        {
            return Err(bad());
        }
        let mut labels = vec![active_port.as_str()];
        for label in terminations
            .iter()
            .map(|t| t.label.as_str())
            .chain(coupled_labels.iter().map(String::as_str))
        {
            if !ports.contains(&label) || labels.contains(&label) {
                return Err(bad());
            }
            labels.push(label);
        }
        Ok(Self {
            network,
            active_port,
            terminations,
            coupled_labels,
            schema,
        })
    }

    /// Complete pressure/volume-flow response for unit flow into the active port.
    /// Declared coupled ports refuse here; they need [`Self::response_coupled`].
    pub fn response(
        &self,
        gas: &GasState,
        omega: f64,
        loss: LossModel,
    ) -> Result<NetworkDriveResponse, DuctError> {
        if !self.coupled_labels.is_empty() {
            return Err(DuctError::BadParameter {
                what: "network driving point with coupled ports needs their coupled termination",
            });
        }
        network_port_impedance(&self.network, gas, omega, loss)?.drive(
            &self.active_port,
            &self.terminations,
            gas,
        )
    }

    /// [`Self::response`] with the declared coupled ports terminated through
    /// `coupled`, whose labels must be exactly those ports.
    pub fn response_coupled(
        &self,
        gas: &GasState,
        omega: f64,
        loss: LossModel,
        coupled: &CoupledTermination,
    ) -> Result<NetworkDriveResponse, DuctError> {
        if coupled.labels.len() != self.coupled_labels.len()
            || coupled
                .labels
                .iter()
                .any(|l| !self.coupled_labels.contains(l))
        {
            return Err(DuctError::BadParameter {
                what: "coupled termination labels differ from the declared coupled ports",
            });
        }
        network_port_impedance(&self.network, gas, omega, loss)?.drive_coupled(
            &self.active_port,
            &self.terminations,
            coupled,
            gas,
        )
    }

    /// Realize the driving point as a causal characteristic line with the
    /// ladder, Fourier convention and passivity gate of [`super::realize_line`].
    /// `R = (Z - zc)/(Z + zc)` is sampled from the network solve at every
    /// ladder frequency and at the realizer's near-DC point (1e-3 rad/s).
    /// A network solve may refuse that close to DC as ill-conditioned
    /// (`Singular`). Below the largest rung's first bin, such a refusal is
    /// answered from the lowest frequency in decades up to that bin at which
    /// the solve is admitted. The realizer's own fallback, `R = -1`, ignores
    /// the Poiseuille resistance of a lossy bore. Its DC jump (`R(0+) = -0.945`
    /// on a 2 mm bore) left a tail that only halved per rung.
    /// `coupled_at` supplies the coupled termination at each sampled
    /// frequency, and is required exactly when coupled ports are declared.
    /// The nominal round trip is twice the longest shortest axial path from
    /// the active port to any node.
    ///
    /// # Errors
    /// A missing or unexpected coupled source, network or coupled-load
    /// refusals, ladder non-convergence, or a non-passive converged line.
    pub fn realize_line(
        &self,
        gas: &GasState,
        loss: LossModel,
        sample_rate_hz: u32,
        zc: f64,
        ladder: &LineLadder,
        mut coupled_at: Option<CoupledSource<'_>>,
    ) -> Result<(DelayedFilter, LineRealization), DrivingPointError> {
        if self.coupled_labels.is_empty() != coupled_at.is_none() {
            return Err(DrivingPointError::Invalid {
                what: "network line needs a coupled source exactly when coupled ports are declared",
            });
        }
        // The single-duct arithmetic, so a one-edge network keeps its ladder.
        let dt = 1.0 / f64::from(sample_rate_hz);
        let geo_delay = 2.0 * self.longest_path()? / gas.sound_speed / dt;
        let first_bin = core::f64::consts::TAU * f64::from(sample_rate_hz) / ladder.max_taps as f64;
        let mut solve = |omega: f64| -> Result<C64, DuctError> {
            let response = match coupled_at.as_mut() {
                Some(source) => self.response_coupled(gas, omega, loss, &source(omega)?)?,
                None => self.response(gas, omega, loss)?,
            };
            Ok(reflectance(response.impedance, zc))
        };
        realize_characteristic_response(
            geo_delay,
            // Selects the realizer's near-DC sample of R for the DC bin; the
            // network's own loads decide its value.
            Termination::IdealOpen,
            sample_rate_hz,
            zc,
            vec![bulk_delay_pole_state(geo_delay)],
            ladder,
            |omega| match solve(omega) {
                Err(DuctError::Singular) if omega < first_bin => {
                    let mut w = omega;
                    loop {
                        w = (10.0 * w).min(first_bin);
                        match solve(w) {
                            Err(DuctError::Singular) if w < first_bin => {}
                            other => return other,
                        }
                    }
                }
                other => other,
            },
        )
    }

    /// Radius of the active port's plane: the matching end of its one edge.
    ///
    /// # Errors
    /// An active port with other than one incident edge.
    pub fn active_port_radius(&self) -> Result<f64, DuctError> {
        let node = self.active_node();
        let mut ends = self.network.edges.iter().filter_map(|e| {
            if e.from == node {
                e.duct.segments.first().map(Segment::inlet_radius)
            } else if e.to == node {
                e.duct.segments.last().map(Segment::outlet_radius)
            } else {
                None
            }
        });
        match (ends.next(), ends.next()) {
            (Some(radius), None) => Ok(radius),
            _ => Err(DuctError::BadParameter {
                what: "network active port needs exactly one incident edge",
            }),
        }
    }

    fn active_node(&self) -> usize {
        self.network
            .nodes
            .iter()
            .position(|n| match &n.boundary {
                NodeBoundary::Port { label } => *label == self.active_port,
                _ => false,
            })
            .expect("active port admitted at construction")
    }

    // Bellman-Ford over axial edge lengths; validated graphs are connected.
    fn longest_path(&self) -> Result<f64, DrivingPointError> {
        let mut distance = vec![f64::INFINITY; self.network.nodes.len()];
        distance[self.active_node()] = 0.0;
        for _ in 0..self.network.nodes.len() {
            for e in &self.network.edges {
                let length: f64 = e.duct.segments.iter().map(Segment::axial_length).sum();
                let (a, b) = (distance[e.from], distance[e.to]);
                distance[e.to] = distance[e.to].min(a + length);
                distance[e.from] = distance[e.from].min(b + length);
            }
        }
        let longest = distance.iter().copied().fold(0.0, f64::max);
        if longest > 0.0 && longest.is_finite() {
            Ok(longest)
        } else {
            Err(DrivingPointError::Invalid {
                what: "network line needs a positive finite axial path",
            })
        }
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
