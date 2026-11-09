//! Initialized passive storage behind an explicit pressure source.
//!
//! The original `fs_vfit::impedance::SeriesImpedance` owns R, L and C,
//! their midpoint update, physical heat and quadratic energy. This adapter
//! reverses its flow into the source: `p = P_s + p_load(-U_out)`.
//! Thus `delta H + heat + W_out = W_source`, with
//! `W_source = P_s * U_out * dt`. Initially charged compliance emits pressure
//! even when P_s is zero; initial storage is never counted as supplied work.
//!
//! A valve can solve the original `UpstreamPort` and then commit its flow.
//! A characteristic duct can instead attach this source to its sole additional
//! `Driven` endpoint. That path uses the original simultaneous wave/load
//! solve and accepts source plus network atomically, with no FIR approximation
//! of the compliance pole or extra sample of boundary delay.
//!
//! Coefficients and topology are fixed. This is the declared linear R-L-C
//! model, not a finite thermodynamic gas inventory, moving lung volume, inferred
//! tissue parameters or a certificate for a fitted airway FIR. In particular,
//! it must not be superposed on an airway already terminated by the same source.

use super::UpstreamPort;
use fs_vfit::impedance::{ImpedanceFrame, SeriesImpedance};
use fs_vfit::waveguide::{
    WaveguideError,
    network::{DrivenBoundary, DrivenWave, NetworkFrame, WaveguideNetwork},
};

/// One accepted source sample. Positive flow and outgoing work leave the source.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StoredSourceFrame {
    /// Original passive owner's state, storage, irreversible heat and work.
    /// Its flow is INTO the source, opposite to `flow_out_m3_s`.
    pub load: ImpedanceFrame,
    /// Step-mean pressure at the actual source/duct or source/valve port [Pa].
    pub pressure_pa: f64,
    /// Step-mean volume flow leaving this source [m^3/s].
    pub flow_out_m3_s: f64,
    /// Signed work supplied by the explicit ideal pressure source [J].
    pub source_work_j: f64,
    /// Signed work leaving through the physical coupled port [J].
    pub outgoing_work_j: f64,
}

impl StoredSourceFrame {
    /// Local balance [J], with initial storage and reversible release retained.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.load.storage_change_j + self.load.dissipated_energy_j + self.outgoing_work_j
            - self.source_work_j
    }
}

/// One simultaneous source and original wave-network step.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct StoredSourceNetworkFrame {
    pub network: NetworkFrame,
    pub source: StoredSourceFrame,
}

impl StoredSourceNetworkFrame {
    /// Independently evaluated endpoint storage of both owners [J].
    #[must_use]
    pub fn stored_energy_j(&self) -> f64 {
        self.network.stored_energy_j + self.source.load.stored_energy_j
    }

    /// Combined balance [J]. Source/duct work is internal and counted once.
    #[must_use]
    pub fn balance_residual_j(&self) -> f64 {
        self.network.storage_change_j
            + self.source.load.storage_change_j
            + self.network.terminal_loss_j
            + self.network.interior_loss_j
            + self.network.coupled_terminal_work_j
            + self.source.load.dissipated_energy_j
            - self.network.inlet_work_j
            - self.source.source_work_j
    }
}

/// A pressure source with retained, possibly initialized R-L-C storage.
///
/// Construct the original `SeriesImpedance` with `new_with_state` to supply
/// initial compliance pressure or inertive flow. Its passive flow coordinate
/// points INTO this source. Construction and stepping allocate no heap storage.
#[derive(Clone, Copy, Debug)]
pub struct StoredPressureSource {
    load: SeriesImpedance,
    pending: Option<f64>,
}

impl StoredPressureSource {
    /// Retain the original admitted impedance, physical history and timestep.
    #[must_use]
    pub const fn new(load: SeriesImpedance) -> Self {
        Self {
            load,
            pending: None,
        }
    }

    /// Immutable original storage owner; no replacement or resetting of history.
    #[must_use]
    pub const fn impedance(&self) -> &SeriesImpedance {
        &self.load
    }

    /// Prepare the original valve-generic port with this step's ideal source
    /// pressure [Pa]. Repeated preparation replaces only the pending source,
    /// never advances storage. Solve flow using this exact port, then commit.
    pub fn port(&mut self, source_pa: f64) -> Result<UpstreamPort, WaveguideError> {
        let (history, impedance) = self.load.affine_pressure()?;
        let open_circuit_pa = source_pa + history;
        if !source_pa.is_finite() || !open_circuit_pa.is_finite() {
            return Err(WaveguideError("stored pressure source port must be finite"));
        }
        self.pending = Some(source_pa);
        Ok(UpstreamPort {
            open_circuit_pa,
            impedance,
        })
    }

    fn frame(
        &self,
        load: ImpedanceFrame,
        source_pa: f64,
    ) -> Result<StoredSourceFrame, WaveguideError> {
        let pressure_pa = source_pa + load.pressure_pa;
        let flow_out_m3_s = -load.flow_m3_s;
        let volume = flow_out_m3_s * self.load.time_step_s();
        let frame = StoredSourceFrame {
            load,
            pressure_pa,
            flow_out_m3_s,
            source_work_j: source_pa * volume,
            outgoing_work_j: pressure_pa * volume,
        };
        if ![
            source_pa,
            pressure_pa,
            flow_out_m3_s,
            frame.source_work_j,
            frame.outgoing_work_j,
            frame.balance_residual_j(),
        ]
        .iter()
        .all(|x| x.is_finite())
        {
            return Err(WaveguideError(
                "stored pressure source work left the finite set",
            ));
        }
        Ok(frame)
    }

    /// Observe the next source state and original work at the solved outgoing
    /// flow. A preview or refusal never advances or consumes the pending port.
    pub fn preview_commit(&self, flow_out_m3_s: f64) -> Result<StoredSourceFrame, WaveguideError> {
        let source_pa = self.pending.ok_or(WaveguideError(
            "stored pressure source commit needs this sample's prepared port",
        ))?;
        self.frame(self.load.preview_at_flow(-flow_out_m3_s)?, source_pa)
    }

    /// Publish the solved source step only after all original work is finite.
    /// A refusal retains storage and the pending port, permitting a retry.
    pub fn commit(&mut self, flow_out_m3_s: f64) -> Result<StoredSourceFrame, WaveguideError> {
        let frame = self.preview_commit(flow_out_m3_s)?;
        // Copy only this bounded two-coordinate state. No stale or caller-forged
        // frame can be committed, and a failed source cannot partly advance.
        let mut next = self.load;
        next.step_at_flow(-flow_out_m3_s)?;
        self.load = next;
        self.pending = None;
        Ok(frame)
    }

    /// Shared physical boundary candidate for larger native transactions.
    /// Only this bounded source state is copied; the caller retains and checks
    /// the original network clock and publishes the candidate after all peers.
    pub(crate) fn preview_boundary(
        &self,
        boundary: DrivenBoundary,
        source_pa: f64,
    ) -> Result<(Self, StoredSourceFrame, DrivenWave), WaveguideError> {
        if boundary.impedance_pa_s_m3.to_bits() != self.load.port_impedance().to_bits() {
            return Err(WaveguideError(
                "stored source and adjacent characteristic references must match",
            ));
        }
        // p_load = p_node - P_s; subtract half the series source from both
        // characteristic waves. The native boundary solves its full history.
        let mut next = *self;
        let load = next
            .load
            .step(boundary.incoming_pressure_pa - 0.5 * source_pa)?;
        let frame = self.frame(load, source_pa)?;
        let outgoing_pressure_pa = load.reflected_pressure_pa + 0.5 * source_pa;
        if !outgoing_pressure_pa.is_finite() {
            return Err(WaveguideError(
                "stored pressure source outgoing wave must be finite",
            ));
        }
        next.pending = None;
        Ok((
            next,
            frame,
            DrivenWave {
                node: boundary.node,
                outgoing_pressure_pa,
            },
        ))
    }

    /// Couple to the sole additional driven endpoint of an existing native
    /// characteristic network. `outgoing_at_inlet_pa` is the independently
    /// solved inlet wave; all endpoints read the same old traveling state.
    ///
    /// Source and adjacent wave references and clocks must match exactly.
    /// The network keeps its own admitted memory/work bounds. Cancellation is
    /// polled before preparation and after the original native preview, before
    /// either owner publishes. Interruption inside that preview is not added.
    ///
    /// # Errors
    /// Wrong topology/reference/clock, cancellation, or either owner's native
    /// refusal. Neither accepted clock nor physical state moves on failure.
    pub fn step_network(
        &mut self,
        network: &mut WaveguideNetwork,
        source_node: usize,
        outgoing_at_inlet_pa: f64,
        source_pa: f64,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<StoredSourceNetworkFrame, WaveguideError> {
        if cancelled() {
            return Err(WaveguideError(
                "stored pressure source network step cancelled",
            ));
        }
        if network.time_step_s().to_bits() != self.load.time_step_s().to_bits()
            || !network.driven_nodes().eq(core::iter::once(source_node))
        {
            return Err(WaveguideError(
                "stored source requires the same clock and sole driven endpoint",
            ));
        }
        let boundary = network.driven_boundary(source_node)?;
        let (next, source, wave) = self.preview_boundary(boundary, source_pa)?;
        let frame = network.step_with_driven_checked(outgoing_at_inlet_pa, &[wave], |network| {
            let combined = StoredSourceNetworkFrame {
                network: *network,
                source,
            };
            if !combined.stored_energy_j().is_finite() || !combined.balance_residual_j().is_finite()
            {
                return Err(WaveguideError(
                    "combined source/network energy is not representable",
                ));
            }
            if cancelled() {
                return Err(WaveguideError(
                    "stored pressure source network step cancelled",
                ));
            }
            Ok(())
        })?;
        *self = next;
        Ok(StoredSourceNetworkFrame {
            network: frame,
            source,
        })
    }
}

#[cfg(test)]
mod tests;
