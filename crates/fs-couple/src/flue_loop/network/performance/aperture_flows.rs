//! Read-only record of every observed opening's solved volume flow.
//!
//! On every mechanical frame the tap copies the same signed flows that the
//! baffled receiver consumes: each selected outlet's net flow into its terminal
//! node and, when a mouth is attached, the flow into the source node. Both are
//! outward volume flows [m^3/s] (see `MouthPressure::push`). Nothing is fed
//! back: the source, its loads and the microphone signal are unchanged, so an
//! exterior observer of another geometry can replay exactly this source.
//!
//! File layout: a text header ending in a `data` line, then one little-endian
//! f64 per aperture per mechanical frame, in header order. The sink is flushed
//! and checked when the declared window completes; a short file is an error,
//! never a silent truncation.
//!
//! The labium variant (schema v2) appends two f64 to every frame: the jet
//! displacement and labium offset that the frame's accepted cut used. Its
//! header states the cut half-width. A v1 reader refuses it by its first line.
use super::{FlueOutletPerformance, invalid};
use crate::flue_loop::{FlueError, LabiumSample};
use crate::flue_loop::network::FlueNetwork;
use crate::render::schedule::GestureCompileError;
use std::io::Write;

/// First line of an aperture-flow file.
pub const APERTURE_FLOWS_SCHEMA: &str = "frankensim-aperture-flows-v1";
/// First line of an aperture-flow file that also carries the labium displacement.
pub const APERTURE_FLOWS_LABIUM_SCHEMA: &str = "frankensim-aperture-flows-v2";
/// Bytes buffered before a write to the sink.
const CHUNK_BYTES: usize = 1 << 20;

pub(super) struct ApertureFlowTap {
    nodes: Vec<usize>,
    labium: bool,
    sink: Box<dyn Write + Send>,
    buffer: Vec<u8>,
    frames: u64,
    total_frames: u64,
}

impl ApertureFlowTap {
    /// Append one mechanical frame; flush and check the sink after the last one.
    pub(super) fn record(
        &mut self,
        network: &FlueNetwork,
        labium: LabiumSample,
    ) -> Result<(), FlueError> {
        if self.frames == self.total_frames {
            return Err(FlueError::Invalid {
                what: "aperture flow tap ran past its declared window",
            });
        }
        for &node in &self.nodes {
            let flow = network
                .node_frame(node)
                .ok_or(FlueError::Invalid {
                    what: "tapped aperture has no flow",
                })?
                .net_flow_into_node_m3_s;
            self.buffer.extend_from_slice(&flow.to_le_bytes());
        }
        if self.labium {
            self.buffer.extend_from_slice(&labium.eta_m.to_le_bytes());
            self.buffer.extend_from_slice(&labium.offset_m.to_le_bytes());
        }
        self.frames += 1;
        if self.buffer.len() >= CHUNK_BYTES || self.frames == self.total_frames {
            self.sink
                .write_all(&self.buffer)
                .map_err(|_| FlueError::Invalid {
                    what: "aperture flow tap write failed",
                })?;
            self.buffer.clear();
        }
        if self.frames == self.total_frames {
            self.sink.flush().map_err(|_| FlueError::Invalid {
                what: "aperture flow tap flush failed",
            })?;
        }
        Ok(())
    }
}

impl FlueOutletPerformance {
    /// Copy every observed opening's outward volume flow, on every mechanical
    /// frame of the declared window, to `sink`. The header names each aperture
    /// by role, network node and baffle-plane center, in data order: the
    /// selected outlets as `outlet_nodes()` orders them, then the mouth.
    ///
    /// Observation only: the rendered pressure is bit-identical with or without
    /// the tap. Attach once, before the first sample.
    ///
    /// # Errors
    /// A rendered, poisoned or already tapped performance, a window too long to
    /// count, or a failed header write.
    pub fn with_aperture_flow_tap(
        self,
        sink: Box<dyn Write + Send>,
    ) -> Result<Self, GestureCompileError> {
        self.attach_flow_tap(sink, false)
    }

    /// The same tap, schema v2: after each frame's apertures it writes the jet
    /// displacement `eta_m` and labium offset `offset_m` [m] that the frame's
    /// accepted cut `0.5 q_max (1 + tanh((eta_m - offset_m)/b))` used, b being
    /// the header's cut half-width. The header also gives the forced jet's
    /// exit half-width, or `none`.
    ///
    /// # Errors
    /// As [`Self::with_aperture_flow_tap`].
    pub fn with_aperture_flow_and_labium_tap(
        self,
        sink: Box<dyn Write + Send>,
    ) -> Result<Self, GestureCompileError> {
        self.attach_flow_tap(sink, true)
    }

    fn attach_flow_tap(
        mut self,
        mut sink: Box<dyn Write + Send>,
        labium: bool,
    ) -> Result<Self, GestureCompileError> {
        if self.poisoned || self.voice.samples_rendered() != 0 || self.flow_tap.is_some() {
            return Err(invalid(
                "an aperture flow tap is attached once, before rendering",
            ));
        }
        let total_frames = self
            .config
            .samples
            .checked_mul(self.mechanics_per_output as u64)
            .ok_or_else(|| invalid("aperture flow tap window overflows"))?;
        let centers = std::iter::once([0.0; 2]).chain(
            self.radiation
                .iter()
                .flat_map(|r| r.additional_outlets.iter().map(|o| o.center_m)),
        );
        let mut nodes = Vec::new();
        let schema = if labium { APERTURE_FLOWS_LABIUM_SCHEMA } else { APERTURE_FLOWS_SCHEMA };
        let mut header = format!(
            "{schema}\nmechanics_rate_hz {}\noutput_rate_hz {}\nwindow_frames {total_frames}\n",
            self.mechanics_rate_hz, self.config.sample_rate_hz
        );
        for (&node, center) in self.outlet_nodes.iter().zip(centers) {
            header += &format!(
                "aperture outlet node {node} center_m {:e} {:e}\n",
                center[0], center[1]
            );
            nodes.push(node);
        }
        if let Some(mouth) = &self.mouth {
            let center = mouth.input.center_m;
            header += &format!(
                "aperture mouth node {} center_m {:e} {:e}\n",
                mouth.inlet_node, center[0], center[1]
            );
            nodes.push(mouth.inlet_node);
        }
        header += "values outward_volume_flow_m3_s f64le aperture_order_per_mechanical_frame\n";
        if labium {
            let forced = self.voice.forced_jet().map_or_else(
                || "none".to_owned(),
                |jet| format!("{:e}", jet.half_width_m()),
            );
            let cut_half_width_m = 0.5 * self.voice.geometry.flue_height_m;
            header += "labium values eta_m offset_m f64le after_apertures_per_mechanical_frame\n";
            header += &format!("labium cut_half_width_m {cut_half_width_m:e}\n");
            header += &format!("labium forced_jet_half_width_m {forced}\n");
        }
        header += "data\n";
        sink.write_all(header.as_bytes())
            .map_err(|_| invalid("aperture flow tap header write failed"))?;
        self.flow_tap = Some(ApertureFlowTap {
            nodes,
            labium,
            sink,
            buffer: Vec::with_capacity(CHUNK_BYTES + 64),
            frames: 0,
            total_frames,
        });
        Ok(self)
    }

    /// Mechanical frames recorded by an attached aperture flow tap.
    #[must_use]
    pub fn aperture_flow_frames(&self) -> Option<u64> {
        self.flow_tap.as_ref().map(|t| t.frames)
    }

    /// Write every buffered frame and flush the sink, for a run that stops
    /// before its declared window ends. The complete window flushes itself.
    ///
    /// # Errors
    /// No tap is attached, or the sink refuses the write or flush.
    pub fn flush_aperture_flow_tap(&mut self) -> Result<(), GestureCompileError> {
        let tap = self
            .flow_tap
            .as_mut()
            .ok_or_else(|| invalid("no aperture flow tap is attached"))?;
        tap.sink
            .write_all(&tap.buffer)
            .and_then(|()| tap.sink.flush())
            .map_err(|_| invalid("aperture flow tap flush failed"))?;
        tap.buffer.clear();
        Ok(())
    }
}
