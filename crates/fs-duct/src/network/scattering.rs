//! Coherent reflection and transmission at real-reference characteristic ports.
//!
//! With flow directed INTO the graph, p=a+b and U=(a-b)/R. Therefore the
//! original pressure-source boundary p+R U=2a supplies the scattering problem
//! directly, including frequencies where an open-circuit impedance has a pole.
//! One prepared boundary factorization serves every incident-wave column.
//!
//! Pressure waves are peak Pa phasors on the network's exp(-i omega t) clock.
//! Power-normalized waves are a/sqrt(R), so mean power is |a/sqrt(R)|²/2.
//! This is the real-reference normalization in J. O. Smith, Physical Audio
//! Signal Processing, "Normalized Scattering":
//! <https://www.dsprelated.com/freebooks/pasp/Normalized_Scattering.html>.
//! Reference impedances are explicit coordinates, not inferred radiation or
//! constitutive laws. Complex reference impedances and time-varying ports are
//! outside this API; original propagation, load and condition bounds remain.

use super::{
    CoupledTermination, DuctError, DuctNetwork, GasState, LossModel, MAX_NETWORK_PORTS,
    NetworkMedia, PortExcitation, PortSource, PortTermination, PreparedNetworkSources, bad,
    coherent_sum, finite, network_prepare_sources, network_prepare_sources_regional,
};
use crate::fluid::FluidDuctError;
use fs_material::fluid_acoustics::FluidAcoustics;
use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluid;
use fs_math::{c64::C64, det};

/// A labelled characteristic-wave port with a real positive reference.
#[derive(Clone, Debug, PartialEq)]
pub struct WavePort {
    /// Exact external port label in the original graph.
    pub label: String,
    /// Real reference impedance R [Pa s/m³], used in p=a+b, U=(a-b)/R.
    pub reference_impedance_pa_s_m3: f64,
}

/// Fixed-frequency N-port scattering with every omitted port explicitly loaded.
/// Both matrix axes and incident/output vectors use lexicographic label order.
/// Preparation allocates bounded matrices; applying a new coherent incident
/// vector performs no allocation, factorization or material reevaluation.
#[derive(Clone, Debug)]
pub struct PreparedNetworkScattering {
    sources: PreparedNetworkSources,
    ports: Vec<WavePort>,
    pressure: Vec<C64>,
    power: Vec<C64>,
}

impl PreparedNetworkScattering {
    /// Prepare waves in a uniform gas. Wave ports and the scalar/coupled loads
    /// must partition all external ports exactly once. Unexcited wave ports
    /// retain their reference impedance; they are not closed ends.
    pub fn new(
        net: &DuctNetwork,
        gas: &GasState,
        omega: f64,
        loss: LossModel,
        ports: &[WavePort],
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<Self, DuctError> {
        let sources = wave_sources(ports)?;
        let prepared =
            network_prepare_sources(net, gas, omega, loss, &sources, None, loads, coupled)?;
        Self::from_prepared(prepared, ports)
    }

    /// The same scattering problem with the original per-edge gases and
    /// independent exterior material. Reference R never replaces these laws.
    pub fn regional(
        net: &DuctNetwork,
        media: &NetworkMedia,
        omega: f64,
        loss: LossModel,
        ports: &[WavePort],
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<Self, DuctError> {
        let sources = wave_sources(ports)?;
        let prepared = network_prepare_sources_regional(
            net, media, omega, loss, &sources, None, loads, coupled,
        )?;
        Self::from_prepared(prepared, ports)
    }

    /// Scattering in one original gas/liquid material. Its actual rho and c
    /// drive the existing branched pressure/flow graph and same-fluid loads;
    /// the explicit real wave references never substitute for that material.
    /// Prepared application retains the same allocation-free matrices.
    ///
    /// # Errors
    /// Original material, graph, source/load and wave-reference refusals from
    /// [`super::fluid::prepare_sources`]. Open-circuit poles need no inverse.
    pub fn from_material<M: FluidAcoustics + ?Sized>(
        net: &DuctNetwork,
        material: &M,
        omega: f64,
        ports: &[WavePort],
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<Self, FluidDuctError> {
        let sources = wave_sources(ports)?;
        let prepared =
            super::fluid::prepare_sources(net, material, omega, &sources, None, loads, coupled)?;
        Ok(Self::from_prepared(prepared, ports)?)
    }

    /// Material-driven scattering with original transport/caloric data and
    /// an explicit lossy-cone refinement axis. Radiation and porous boundaries
    /// use the same single complete snapshot as every propagation edge.
    ///
    /// # Errors
    /// As [`super::fluid::prepare_sources_with_loss`], plus wave-reference
    /// normalization refusals. No clipping or passive projection is applied.
    #[allow(clippy::too_many_arguments)]
    pub fn from_material_with_loss<M: ThermoviscousFluid + ?Sized>(
        net: &DuctNetwork,
        material: &M,
        omega: f64,
        loss: LossModel,
        cone_slices: crate::ConeSlices,
        ports: &[WavePort],
        loads: &[PortTermination],
        coupled: Option<&CoupledTermination>,
    ) -> Result<Self, FluidDuctError> {
        let sources = wave_sources(ports)?;
        let prepared = super::fluid::prepare_sources_with_loss(
            net,
            material,
            omega,
            loss,
            cone_slices,
            &sources,
            None,
            loads,
            coupled,
        )?;
        Ok(Self::from_prepared(prepared, ports)?)
    }

    fn from_prepared(
        sources: PreparedNetworkSources,
        ports: &[WavePort],
    ) -> Result<Self, DuctError> {
        let ports: Vec<_> = sources
            .source_coordinates()
            .iter()
            .map(|source| {
                ports
                    .iter()
                    .find(|port| port.label == source.label)
                    .expect("prepared sources retain admitted wave labels")
                    .clone()
            })
            .collect();
        let n = ports.len();
        let mut pressure = vec![C64::ZERO; n * n];
        let mut power = pressure.clone();
        for (row, &port_index) in sources.source_port_indices().iter().enumerate() {
            let reference = ports[row].reference_impedance_pa_s_m3;
            for column in 0..n {
                let index = port_index * n + column;
                // Each retained column was driven by p_open=1. For a=1,
                // b=(p-R U)/2 receives the factor two from p_open=2a.
                let value = sources.pressure_transfer()[index]
                    - sources.flow_transfer()[index].scale(reference);
                // sqrt(Rj)/sqrt(Ri) avoids overflowing Rj/Ri first.
                let scale =
                    det::sqrt(ports[column].reference_impedance_pa_s_m3) / det::sqrt(reference);
                let normalized = value.scale(scale);
                if !finite(value) || !finite(normalized) || !scale.is_finite() || scale <= 0.0 {
                    return Err(bad("network wave normalization is unrepresentable"));
                }
                pressure[row * n + column] = value;
                power[row * n + column] = normalized;
            }
        }
        Ok(Self {
            sources,
            ports,
            pressure,
            power,
        })
    }

    /// Fixed physical frequency [rad/s], inherited from the original graph.
    #[must_use]
    pub fn omega(&self) -> f64 {
        self.sources.omega()
    }

    /// Canonical matrix/vector order and original reference impedances.
    #[must_use]
    pub fn ports(&self) -> &[WavePort] {
        &self.ports
    }

    /// Row-major b=S a for complex peak pressure waves [Pa]. Diagonal entries
    /// reflect; off-diagonal entries transmit. Their phases are retained.
    #[must_use]
    pub fn pressure_matrix(&self) -> &[C64] {
        &self.pressure
    }

    /// Row-major scattering of a/sqrt(R). A reciprocal lossless graph with
    /// no absorbing loads has a symmetric unitary matrix in these coordinates.
    /// No clipping, passivity projection or discarded interference is applied.
    #[must_use]
    pub fn power_normalized_matrix(&self) -> &[C64] {
        &self.power
    }

    /// Original pressure/flow transfer and work/statistics owner. Its inputs
    /// are p_open=2a in the same canonical order, NOT incident pressure a.
    #[must_use]
    pub fn prepared_sources(&self) -> &PreparedNetworkSources {
        &self.sources
    }

    /// Apply coherent incident peak pressures to reflected/transmitted peak
    /// pressures. Invalid/nonfinite input or output preserves the output slice.
    pub fn scatter_into(&self, incident: &[C64], outgoing: &mut [C64]) -> Result<(), DuctError> {
        self.apply(&self.pressure, incident, outgoing)
    }

    /// Apply power-normalized incident waves; mean incident/outgoing power is
    /// respectively sum(|wave|²)/2. This retains coherent cross-port terms.
    pub fn scatter_power_normalized_into(
        &self,
        incident: &[C64],
        outgoing: &mut [C64],
    ) -> Result<(), DuctError> {
        self.apply(&self.power, incident, outgoing)
    }

    fn apply(
        &self,
        matrix: &[C64],
        incident: &[C64],
        outgoing: &mut [C64],
    ) -> Result<(), DuctError> {
        let n = self.ports.len();
        if incident.len() != n || outgoing.len() != n || incident.iter().any(|&v| !finite(v)) {
            return Err(bad(
                "network scattering needs complete finite canonical wave vectors",
            ));
        }
        let mut result = [C64::ZERO; MAX_NETWORK_PORTS];
        for row in 0..n {
            result[row] =
                coherent_sum((0..n).map(|column| matrix[row * n + column] * incident[column]));
            if !finite(result[row]) {
                return Err(bad("network scattered wave is unrepresentable"));
            }
        }
        outgoing.copy_from_slice(&result[..n]);
        Ok(())
    }
}

fn wave_sources(ports: &[WavePort]) -> Result<Vec<PortSource>, DuctError> {
    if ports.is_empty() || ports.len() > MAX_NETWORK_PORTS {
        return Err(bad("network scattering needs 1..=32 wave ports"));
    }
    let mut sources = Vec::with_capacity(ports.len());
    for (index, port) in ports.iter().enumerate() {
        let reference = port.reference_impedance_pa_s_m3;
        if port.label.trim().is_empty()
            || !reference.is_finite()
            || reference <= 0.0
            || ports[..index].iter().any(|prior| prior.label == port.label)
        {
            return Err(bad(
                "network wave ports need unique labels and finite positive real references",
            ));
        }
        sources.push(PortSource {
            label: port.label.clone(),
            excitation: PortExcitation::Pressure {
                open_circuit_pressure_pa: C64::ZERO,
                impedance_pa_s_m3: C64::from_re(reference),
            },
        });
    }
    Ok(sources)
}

#[cfg(test)]
mod tests;
