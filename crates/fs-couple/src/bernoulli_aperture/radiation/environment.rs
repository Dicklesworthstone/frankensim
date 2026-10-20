//! One exterior material state for radiation reaction AND observed pressure.
//!
//! Cold composition only: DuctPressure owns port/medium/receiver admission;
//! BaffledRadiationLoad owns the compact positive-real boundary law. Selected
//! memoryless endpoints become those physical loads before propagation starts.
//! No reaction is inferred merely from attaching a one-way pressure observer.
use super::BaffledRadiationLoad;
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec};
use crate::pcm_wav::baffled::openings::duct::{
    DuctOpening, DuctPressure, DuctRadiationSpec, DuctView,
};
use fs_exec::CancelGate;
use fs_material::gas::GasState;

/// Original port and medium plus the exact boundary model installed there.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct BaffledOutletBinding {
    opening: DuctOpening,
    interior: GasState,
    previous_reflection: f64,
    load: BaffledRadiationLoad,
}
impl BaffledOutletBinding {
    /// Verify a retained load against the actual graph and observation medium.
    /// A caller's separately supplied gas label cannot override installed physics.
    pub fn validate_receiver(
        &self,
        graph: &TubeNetworkSpec,
        radiation: &DuctRadiationSpec,
        exterior: GasState,
    ) -> Result<(), crate::pcm_wav::observation::ObservationError> {
        use crate::pcm_wav::baffled::RayleighMedium;
        use crate::pcm_wav::observation::{ObservationError, validate_reaction_medium};
        validate_reaction_medium(
            RayleighMedium {
                density: self.load.density_kg_m3(),
                sound_speed: self.load.sound_speed_m_s(),
            },
            RayleighMedium {
                density: exterior.density,
                sound_speed: exterior.sound_speed,
            },
        )?;
        let matches_node = graph.nodes.get(self.opening.node) == Some(&self.load.termination());
        let mut adjacent = graph
            .sections
            .iter()
            .filter(|s| s.nodes.contains(&self.opening.node));
        let section = adjacent.next();
        if !matches_node
            || section.is_none_or(|s| s.radius_m != self.load.radius_m())
            || adjacent.next().is_some()
            || !radiation.openings.contains(&self.opening)
            || self.load.time_step_s().to_bits()
                != (1.0 / f64::from(radiation.reconstruction.mechanics_rate_hz())).to_bits()
            || radiation.receiver.maximum_frequency_hz > self.load.maximum_frequency_hz()
        {
            return Err(ObservationError::ReactionBinding {
                what: "terminal, radius, clock or band differs from the installed load",
            });
        }
        Ok(())
    }
    /// Network/node identity and original physical location in the baffle.
    pub const fn opening(&self) -> DuctOpening { self.opening }
    /// Actual adjacent bore state. It is NOT replaced with the exterior gas.
    pub const fn interior_gas(&self) -> &GasState { &self.interior }
    /// Explicit memoryless condition replaced by this cold construction.
    pub const fn previous_reflection(&self) -> f64 { self.previous_reflection }
    /// Complete original-owner load, including band and approximation checks.
    pub const fn load(&self) -> &BaffledRadiationLoad { &self.load }
}

/// A common exterior selection and its installed independent terminal loads.
/// The networks retain all other geometry, gases, wall laws and source controls.
#[derive(Clone, Debug)]
pub struct BaffledEnvironment {
    radiation: DuctRadiationSpec,
    bindings: Vec<BaffledOutletBinding>,
}
impl BaffledEnvironment {
    /// Attach this same declaration to a pressure performer after constructing
    /// the physical networks. Its exterior is resolved explicitly, including
    /// the legacy one-port local-gas choice, so reaction and observation agree.
    pub const fn radiation(&self) -> &DuctRadiationSpec { &self.radiation }
    /// Boundary models in the same order as the observed openings.
    pub fn bindings(&self) -> &[BaffledOutletBinding] { &self.bindings }
    /// Separate the performer declaration from immutable boundary diagnostics.
    pub fn into_parts(self) -> (DuctRadiationSpec, Vec<BaffledOutletBinding>) {
        (self.radiation, self.bindings)
    }
}
fn invalid(what: &'static str) -> AcousticRealizeError {
    AcousticRealizeError::InvalidDescription { what }
}
fn poll(gate: &CancelGate) -> Result<(), AcousticRealizeError> {
    if gate.is_requested() { Err(AcousticRealizeError::Cancelled) } else { Ok(()) }
}

/// Derive every selected load from its ACTUAL radius and the SAME exterior gas
/// used by the shared microphone. The local bore impedance and propagation gas
/// remain unchanged. This is the stationary pressure-matched two-fluid acoustic
/// interface: volume flow is shared; pressure reacts through the exterior load.
///
/// Only explicitly selected memoryless terminals can be replaced. Existing
/// reactive loads, enclosed wall/cavity nodes, and interior junctions refuse;
/// this function must not silently discard a separately authored physical law.
/// No waveform/state exists in these graph descriptions. Every load, observer,
/// domain and cancellation check passes BEFORE any graph node is changed.
///
/// The boundary image is the existing independent compact baffled piston:
/// ka <= 0.5 and the original analogue/discrete 5% checks still apply. The
/// receiver band is also the load band; a broader observer cannot certify this
/// compact load. Multiple disks interfere at the microphone, but this image
/// omits mutual radiation reaction between them. No free-field instrument-body
/// scattering, mean gas flow/mixing, live atmosphere switch or measured timbre
/// claim follows. Existing one-way observers remain unchanged and available.
///
/// # Errors
/// Mismatched graph/gas lists, original observer/medium/clock refusal, a selected
/// non-memoryless terminal, compact-load applicability, or cancellation. Every
/// input graph remains unchanged on error, including a late load refusal.
pub fn bind_baffled_environment(
    graphs: &mut [TubeNetworkSpec], gases: &[&[GasState]],
    mut radiation: DuctRadiationSpec, gate: &CancelGate,
) -> Result<BaffledEnvironment, AcousticRealizeError> {
    poll(gate)?;
    if graphs.len() != gases.len() {
        return Err(invalid("baffled environment needs one gas mapping per network"));
    }
    let views: Vec<_> = graphs.iter().zip(gases).map(|(graph, gases)|
        DuctView { graph, gases }).collect();
    // Reuse the real receiver's full geometry, static-pressure, duplicate-port,
    // clock and reconstruction admission; do not maintain a second validator.
    let receiver = DuctPressure::new(&views, radiation.clone())
        .map_err(AcousticRealizeError::Nonlinear)?;
    let exterior = *receiver.exterior_gas();
    drop(receiver); // cold admission only; do not keep an unused pressure history
    drop(views);
    radiation.exterior = Some(exterior);
    let dt = 1.0 / f64::from(radiation.reconstruction.mechanics_rate_hz());
    let band = radiation.receiver.maximum_frequency_hz;
    let mut bindings = Vec::new();
    bindings.try_reserve_exact(radiation.openings.len())
        .map_err(|_| invalid("baffled boundary allocation failed"))?;
    for &opening in &radiation.openings {
        poll(gate)?;
        let graph = &graphs[opening.network]; // shared receiver admitted addresses
        let NetworkNode::Termination { reflection } = graph.nodes[opening.node] else {
            return Err(invalid("baffled environment only replaces explicitly selected memoryless terminals"));
        };
        if !reflection.is_finite() || reflection.abs() > 1.0 {
            return Err(invalid("original terminal reflection must be finite and passive"));
        }
        let (index, section) = graph.sections.iter().enumerate()
            .find(|(_, s)| s.nodes.contains(&opening.node)).expect("admitted physical outlet");
        let load = BaffledRadiationLoad::new(section.radius_m, exterior.density,
            exterior.sound_speed, dt, band, gate)?;
        bindings.push(BaffledOutletBinding { opening, interior: gases[opening.network][index],
            previous_reflection: reflection, load });
    }
    poll(gate)?;
    // Infallible publication. No accepted mechanical or observer state is reset.
    for binding in &bindings {
        graphs[binding.opening.network].nodes[binding.opening.node] = binding.load.termination();
    }
    Ok(BaffledEnvironment { radiation, bindings })
}

mod instruments;

/// Source-bound mutual reaction through the original exact-ZOH multiport owner.
pub mod multiport;
