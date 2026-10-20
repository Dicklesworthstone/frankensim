//! Source-bound mutual terminal reaction through the shared MIMO ZOH runtime.
//!
//! This adapter consumes an externally produced impedance model. It does not
//! invent a radiation approximation, refit a material substitution, normalize
//! a pressure signal or promote a sampled passivity check into certification.
//! The first observation binding is the existing coplanar baffled receiver.

/// Bounded native source-artifact admission shared by wind instrument loaders.
pub mod input;
/// Geometry- and gas-derived impedance through the existing BEM and VF owners.
pub mod bake;
/// Shared selection of an imported operator or a computed baffled exterior.
pub mod source;

use super::{invalid, poll};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::network::{NetworkNode, TubeNetworkSpec};
use crate::pcm_wav::baffled::openings::CircularOpening;
use crate::pcm_wav::baffled::openings::duct::{DuctPressure, DuctRadiationSpec, DuctView};
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_vfit::mimo::MimoModel;
use fs_vfit::waveguide::network::multiport::MimoTerminalSpec;

/// Original offline material/geometry declaration and its impedance realization.
/// Uniform signed volume flows enter the apertures in this exact matrix order.
#[derive(Clone, Debug)]
pub struct MimoRadiationModel {
    /// Caller-retained identity of the geometry/material radiation declaration.
    /// This binding preserves the identity; it does not verify an offline bake.
    pub declaration_identity: [u8; 32],
    /// Caller-retained identity of this rational realization and its source fit.
    pub realization_identity: [u8; 32],
    /// Complete gas state at which the source radiation operator was produced.
    pub exterior: GasState,
    /// Source aperture centers/radii in the shared z=0 infinite rigid baffle.
    /// Arbitrary unbaffled surfaces require their matching observer adapter.
    pub apertures: Vec<CircularOpening>,
    /// Complete physical band [Hz], covering the source mechanics Nyquist limit.
    pub band_hz: [f64; 2],
    /// Proper undelayed reciprocal `s=+i omega` impedance, in Pa s/m^3.
    /// Acoustic `exp(-i omega t)` data must already have been conjugated.
    pub impedance: MimoModel,
}

/// Explicit work/storage limits for the original dense rational terminal host.
#[derive(Clone, Copy, Debug)]
pub struct MimoRadiationLimits {
    /// Real state-dimension cap, at most 256.
    pub max_states: usize,
    /// Separate retained terminal payload allowance [bytes].
    pub max_memory_bytes: usize,
    /// Actual sampled-response passivity grid, including DC/Nyquist; 17..=65537.
    pub sampled_passivity_points: usize,
}

/// One retained source, port order and exterior for reaction and microphone.
#[derive(Clone, Debug)]
pub struct MimoEnvironment {
    source: MimoRadiationModel,
    radiation: DuctRadiationSpec,
}

impl MimoEnvironment {
    /// Original source material, aperture basis, identities and rational matrix.
    #[must_use]
    pub const fn source(&self) -> &MimoRadiationModel {
        &self.source
    }
    /// Resolved exterior and ordered receiver declaration used at construction.
    #[must_use]
    pub const fn radiation(&self) -> &DuctRadiationSpec {
        &self.radiation
    }

    pub(crate) fn prepare(
        graph: &TubeNetworkSpec,
        gases: &[GasState],
        source: MimoRadiationModel,
        mut radiation: DuctRadiationSpec,
        limits: MimoRadiationLimits,
        gate: &CancelGate,
    ) -> Result<(Self, MimoTerminalSpec), AcousticRealizeError> {
        poll(gate)?;
        graph.validate_section_gases(gases)?;
        let observer = DuctPressure::new(&[DuctView { graph, gases }], radiation.clone())
            .map_err(AcousticRealizeError::Nonlinear)?;
        let exterior = *observer.exterior_gas();
        radiation.exterior = Some(exterior);
        source
            .exterior
            .validate_coherence()
            .map_err(|_| invalid("MIMO radiation source needs a coherent original gas state"))?;
        if source.exterior != exterior {
            return Err(invalid(
                "MIMO reaction and receiver must use the source operator's exact exterior gas",
            ));
        }
        if source.declaration_identity == [0; 32]
            || source.realization_identity == [0; 32]
            || !source.impedance.reciprocal
            || source.apertures.len() != source.impedance.ports
            || source.apertures.len() != radiation.openings.len()
            || !source.band_hz.iter().all(|f| f.is_finite() && *f > 0.0)
            || source.band_hz[0] >= source.band_hz[1]
            || source.band_hz[1] < 0.5 * f64::from(radiation.reconstruction.mechanics_rate_hz())
        {
            return Err(invalid(
                "MIMO radiation needs source identities, a reciprocal ordered basis and full clock bandwidth",
            ));
        }
        let n = source.impedance.ports;
        // The low-level owner admits complete shapes before indexing matrices;
        // symmetry is checked here because this exterior is reciprocal physics.
        let symmetric_real = |m: &[f64]| {
            m.len() == n * n
                && (0..n * n).all(|i| m[i].to_bits() == m[(i % n) * n + i / n].to_bits())
        };
        let symmetric_complex = |m: &[fs_math::c64::C64]| {
            m.len() == n * n && (0..n * n).all(|i| m[i] == m[(i % n) * n + i / n])
        };
        if !symmetric_real(&source.impedance.d)
            || !symmetric_real(&source.impedance.e)
            || source
                .impedance
                .residues
                .iter()
                .any(|r| !symmetric_complex(r))
        {
            return Err(invalid(
                "MIMO radiation's reciprocal flag must match every original matrix entry",
            ));
        }
        for (aperture, opening) in source.apertures.iter().zip(&radiation.openings) {
            poll(gate)?;
            let section = graph
                .sections
                .iter()
                .find(|s| s.nodes.contains(&opening.node))
                .ok_or_else(|| invalid("MIMO radiation aperture has no actual physical section"))?;
            if opening.network != 0
                || aperture.center_m != opening.center_m
                || aperture.radius_m.to_bits() != section.radius_m.to_bits()
                || !matches!(graph.nodes[opening.node], NetworkNode::Termination { .. })
            {
                return Err(invalid(
                    "MIMO source aperture geometry/order must match selected memoryless outlets",
                ));
            }
        }
        poll(gate)?;
        let spec = MimoTerminalSpec {
            nodes: radiation.openings.iter().map(|o| o.node).collect(),
            model: source.impedance.clone(),
            band_rad_s: source.band_hz.map(|f| 2.0 * core::f64::consts::PI * f),
            max_states: limits.max_states,
            max_memory_bytes: limits.max_memory_bytes,
            sampled_passivity_points: limits.sampled_passivity_points,
        };
        Ok((Self { source, radiation }, spec))
    }

    /// Ensure a later shared pressure performer uses the installed operator's
    /// exact aperture order, geometry, exterior and physical clock. Receiver
    /// position/reconstruction/output rate may change within their own owner;
    /// the source mechanics and material may not be silently relabeled.
    pub fn validate_receiver(
        &self,
        graph: &TubeNetworkSpec,
        radiation: &DuctRadiationSpec,
        exterior: GasState,
    ) -> Result<(), AcousticRealizeError> {
        if exterior != self.source.exterior
            || radiation.openings != self.radiation.openings
            || radiation.reconstruction.mechanics_rate_hz()
                != self.radiation.reconstruction.mechanics_rate_hz()
            || radiation.receiver.maximum_frequency_hz > self.source.band_hz[1]
        {
            return Err(invalid(
                "MIMO reaction/receiver gas, port order, clock or band mismatch",
            ));
        }
        for (aperture, opening) in self.source.apertures.iter().zip(&radiation.openings) {
            let mut adjacent = graph
                .sections
                .iter()
                .filter(|s| s.nodes.contains(&opening.node));
            if adjacent
                .next()
                .is_none_or(|s| s.radius_m.to_bits() != aperture.radius_m.to_bits())
                || adjacent.next().is_some()
            {
                return Err(invalid(
                    "MIMO receiver geometry no longer matches the installed aperture basis",
                ));
            }
        }
        Ok(())
    }
}
