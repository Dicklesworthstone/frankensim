//! Compose the shared boundary/environment selection with existing instruments.
use super::{BaffledEnvironment, BaffledOutletBinding, bind_baffled_environment, poll};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::fractional::FractionalTubePolicy;
use crate::bernoulli_aperture::two_sided::{
    ApertureDuctSpec, LoadedAperture, TwoSidedAperture, TwoSidedLimits,
};
use crate::bernoulli_aperture::two_sided::observer::TwoSidedPressure;
use crate::bernoulli_aperture::wall::material::MaterialWallInput;
use crate::flue_loop::network::FlueNetworkSpec;
use crate::pcm_wav::baffled::openings::duct::DuctRadiationSpec;
use crate::render::RenderError;
use fs_exec::CancelGate;

impl TwoSidedPressure {
    /// Build physical radiation reaction and the microphone from one exterior
    /// declaration. The supplied unadvanced single/paired valve keeps all initial
    /// mechanical/contact/material energy. Its impedance must still be the sum
    /// of the two ORIGINAL inlet impedances; radiation does not retune that port.
    ///
    /// Networks 0/1 are upstream/downstream. Selected memoryless endpoints become
    /// the original compact baffled loads. Independent material-wall assignments
    /// and optional fractional/viscothermal propagation then use their existing
    /// constructors. All original node addresses survive lowering. The same
    /// resolved exterior goes to the shared multirate pressure performer; neither
    /// the source clock nor pressure gestures are resampled or reset here.
    ///
    /// Returned boundary records contain the exact load models used by the
    /// mechanical networks, not separately fitted microphone coefficients.
    /// Existing one-way `with_radiation` remains a different explicit selection.
    /// This fixed-state compact independent-outlet model has the applicability
    /// limits documented by `bind_baffled_environment`, including no mutual
    /// radiation reaction or live mixing/temperature evolution.
    ///
    /// # Errors
    /// Original valve, boundary, material, propagation or observation refusal;
    /// inconsistent clocks, an already advanced valve, or cancellation. No
    /// partially running system or pressure output is returned.
    #[allow(clippy::too_many_arguments)]
    pub fn from_baffled_environment(
        valve: LoadedAperture,
        sources: [ApertureDuctSpec; 2],
        walls: [Option<MaterialWallInput>; 2],
        limits: TwoSidedLimits,
        policies: Option<[FractionalTubePolicy; 2]>,
        radiation: DuctRadiationSpec,
        max_block: usize,
        gate: &CancelGate,
    ) -> Result<(Self, Vec<BaffledOutletBinding>), RenderError> {
        Self::build_baffled_environment(valve, sources, walls, limits, policies,
            radiation, max_block, gate, None)
    }

    /// Original material reed, both material ducts, a fixed no-slip supply
    /// passage and one shared exterior reaction/microphone declaration.
    /// Construct the valve load with combined_inlet_state_with_supply. Its
    /// single/independent/modal source and all original material histories remain
    /// in the same mechanics. The passage is additional supplied geometry, with
    /// the StaticPressure law and the identical inlet gas and physical clock.
    #[allow(clippy::too_many_arguments)]
    pub fn from_baffled_environment_with_supply(
        valve: LoadedAperture,
        sources: [ApertureDuctSpec; 2],
        walls: [Option<MaterialWallInput>; 2],
        channel: impl Into<crate::bernoulli_aperture::slit_channel::FixedChannel>,
        limits: TwoSidedLimits,
        policies: Option<[FractionalTubePolicy; 2]>,
        radiation: DuctRadiationSpec,
        max_block: usize,
        gate: &CancelGate,
    ) -> Result<(Self, Vec<BaffledOutletBinding>), RenderError> {
        Self::build_baffled_environment(valve, sources, walls, limits, policies,
            radiation, max_block, gate, Some(channel.into()))
    }

    #[allow(clippy::too_many_arguments)]
    fn build_baffled_environment(
        valve: LoadedAperture,
        sources: [ApertureDuctSpec; 2],
        walls: [Option<MaterialWallInput>; 2],
        limits: TwoSidedLimits,
        policies: Option<[FractionalTubePolicy; 2]>,
        radiation: DuctRadiationSpec,
        max_block: usize,
        gate: &CancelGate,
        supply: Option<crate::bernoulli_aperture::slit_channel::FixedChannel>,
    ) -> Result<(Self, Vec<BaffledOutletBinding>), RenderError> {
        poll(gate).map_err(RenderError::Voice)?;
        let rate = radiation.reconstruction.mechanics_rate_hz();
        if rate == 0 || max_block == 0 || valve.accepted_steps() != 0
            || valve.time_step_s().to_bits() != (1.0/f64::from(rate)).to_bits() {
            return Err(RenderError::Sizing { what: "baffled reed environment requires an unadvanced valve and its actual mechanical clock" });
        }
        if policies.is_some_and(|ps| ps.iter().any(|p|
            p.maximum_frequency_hz < radiation.receiver.maximum_frequency_hz)) {
            return Err(RenderError::Sizing { what: "fractional propagation must cover the exterior radiation band" });
        }
        let [up, down] = sources;
        let mut graphs = [up.graph, down.graph];
        let environment = bind_baffled_environment(&mut graphs, &[&up.gases, &down.gases],
            radiation, gate).map_err(RenderError::Voice)?;
        let [up_graph, down_graph] = graphs;
        let sources = [ApertureDuctSpec { graph: up_graph, ..up },
            ApertureDuctSpec { graph: down_graph, ..down }];
        let system = match supply {
            Some(channel) => TwoSidedAperture::with_material_walls_and_supply_channel(
                valve, sources, walls, channel, limits, policies, gate),
            None => TwoSidedAperture::with_material_walls(
                valve, sources, walls, limits, policies, gate),
        }.map_err(RenderError::Voice)?;
        let (radiation, bindings) = environment.into_parts();
        let performer = Self::with_observation_policy(
            system,
            radiation,
            max_block,
            crate::pcm_wav::observation::ObservationPolicy::Acceptance,
            &bindings,
        )?;
        poll(gate).map_err(RenderError::Voice)?;
        Ok((performer, bindings))
    }
}

impl FlueNetworkSpec {
    /// Install exterior-derived terminal reaction before handing this source to
    /// the existing integer, fractional or material-wall jet constructor. Network
    /// indices in the shared declaration must be zero for this one graph.
    /// Return the SAME resolved observation declaration alongside the unchanged
    /// jet/gas/loss selections. Its receiver can consume actual terminal flows
    /// through DuctPressure, or configure the existing outlet performance wrapper.
    ///
    /// The original jet card, Reynolds validity, optional contraction loss and
    /// wall-material histories are not replaced. Reaction is computed in the
    /// source network, not simulated by an audio gain or a second jet equation.
    ///
    /// # Errors
    /// The shared compact boundary, geometry, exterior-gas or cancellation refusal.
    pub fn with_baffled_environment(mut self, radiation: DuctRadiationSpec,
        gate: &CancelGate) -> Result<(Self, BaffledEnvironment), AcousticRealizeError>
    {
        let environment = bind_baffled_environment(core::slice::from_mut(&mut self.graph),
            &[&self.gases], radiation, gate)?;
        Ok((self, environment))
    }
}

impl crate::flue_loop::FlueVoice {
    /// Construct the existing jet-driven instrument with one exterior declaration
    /// controlling BOTH terminal reaction and subsequent microphone propagation.
    /// The mechanical rate is read from that declaration, never supplied twice.
    /// Optional material walls and fractional loss/travel use their original
    /// constructors. The returned environment carries the matching observation
    /// declaration and the actual compact boundary models for all selected ports.
    ///
    /// Original source/card, compact-radiation and material applicability limits
    /// remain in force. This does not introduce turbulent jet noise or a second
    /// integrator; missing cards still refuse before acoustic construction.
    #[allow(clippy::too_many_arguments)]
    pub fn from_baffled_environment(
        card: Option<&fs_aeroac::jetcard::JetCard>, displacement_gain: f64,
        seed_amplitude_rel: f64, geometry: crate::flue_loop::FlueGeometry,
        source: FlueNetworkSpec, walls: Option<MaterialWallInput>,
        policy: Option<FractionalTubePolicy>, radiation: DuctRadiationSpec,
        voice_seed: u64, voice_index: u32, gate: &CancelGate,
    ) -> Result<(Self, BaffledEnvironment), crate::flue_loop::FlueError> {
        use crate::flue_loop::FlueError;
        poll(gate).map_err(FlueError::Physics)?;
        let card = card.ok_or(FlueError::NoCard)?;
        let rate = radiation.reconstruction.mechanics_rate_hz();
        if policy.is_some_and(|p| p.maximum_frequency_hz < radiation.receiver.maximum_frequency_hz) {
            return Err(FlueError::Physics(super::invalid("fractional propagation must cover the exterior radiation band")));
        }
        let (source, environment) = source.with_baffled_environment(radiation, gate)
            .map_err(FlueError::Physics)?;
        let voice = match (policy, walls) {
            (Some(policy), walls) => Self::from_fractional_network(Some(card),
                displacement_gain, seed_amplitude_rel, geometry, source, walls,
                policy, rate, voice_seed, voice_index, gate)?,
            (None, Some(walls)) => Self::from_material_network(Some(card),
                displacement_gain, seed_amplitude_rel, geometry, source, walls,
                rate, voice_seed, voice_index, gate)?,
            (None, None) => Self::from_network(Some(card), displacement_gain,
                seed_amplitude_rel, geometry, source, rate, voice_seed, voice_index, gate)?,
        };
        poll(gate).map_err(FlueError::Physics)?;
        Ok((voice, environment))
    }
}
