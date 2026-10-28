//! A fixed physical supply passage in the SAME reed/tract/bore transaction.
//!
//! StaticPressure SlitChannel has the discrete law dp = Z_c Q - p_history.
//! It therefore adds Z_c to the original valve load and p_history to its
//! effective upstream pressure. No second nonlinear solve or feedback delay is
//! needed. The channel is additional supplied geometry, not a resized reed gap.
//! Its fluid storage, wall heat and backward-Euler numerical loss remain separate
//! and enter the complete account once. This is the existing small-signal
//! through-flow approximation, not moving-wall squeeze flow or an anatomical
//! vocal tract model.

use super::{ApertureDuctSpec, FractionalTubePolicy, LoadedAperture,
    TwoSidedAperture, TwoSidedLimits, combined_inlet_state, invalid};
use crate::acoustic_realize::AcousticRealizeError;
use crate::bernoulli_aperture::slit_channel::{FixedChannel, SlitChannelStep};
use crate::bernoulli_aperture::wall::material::MaterialWallInput;
use fs_exec::CancelGate;
use fs_material::gas::GasState;

/// Accepted supply-port exchange, kept separate from the moving reed gap.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SupplyChannelFrame {
    /// Actual pressure lost between source-side and upstream reed face [Pa].
    /// It may be negative when stored channel momentum returns work.
    pub pressure_drop_pa: f64,
    /// Original no-slip field's flow, storage, physical and numerical loss.
    pub channel: SlitChannelStep,
}

/// Resolve the actual gas and complete valve load INCLUDING the fixed passage.
/// Construct any original single, independent-pair or modal material valve with
/// this impedance before composing it with the same channel. Both duct inlets
/// and the channel must retain the identical coherent material gas.
///
/// The channel's geometry and step period determine its discrete impedance.
/// Its history pressure is state, never an additional authored blowing command.
/// A nonlinear kinetic-head terminal refuses rather than silently losing its law.
pub fn combined_inlet_state_with_supply(
    sides: &[ApertureDuctSpec; 2], channel: impl Into<FixedChannel>,
) -> Result<(GasState, f64), AcousticRealizeError> {
    let channel = channel.into();
    let (gas, z) = combined_inlet_state(sides)?;
    if channel.gas() != &gas {
        return Err(invalid("reed supply passage and both physical inlets must share the complete material gas"));
    }
    let total = z + channel.pressure_load()?.impedance_pa_s_m3;
    if !total.is_finite() || total <= z {
        return Err(invalid("reed supply passage impedance cannot be represented beside the actual duct load"));
    }
    Ok((gas, total))
}

impl TwoSidedAperture {
    /// Couple one fixed no-slip passage between the series pressure source and
    /// the upstream reed face. The original upstream duct remains a separate
    /// physical part; its geometry must not include this same passage twice.
    ///
    /// The channel, either acoustic side and every original valve image share
    /// the SAME through-flow and sample transaction. Pressure feedback includes
    /// both the channel's current profile and the original tube reflections.
    /// Existing channel motion is retained as explicit initial fluid energy.
    /// The valve must be unadvanced and use combined_inlet_state_with_supply;
    /// the channel time step must exactly match the valve clock.
    ///
    /// The supplied channel must have fixed walls and StaticPressure terminals.
    /// Its declared planar or circular, incompressible, fully developed scope remains; this
    /// does not attach the fixed mesh to a closing or spatially varying reed gap.
    #[allow(clippy::too_many_arguments)]
    pub fn with_supply_channel(
        valve: LoadedAperture, sources: [ApertureDuctSpec; 2], channel: impl Into<FixedChannel>,
        limits: TwoSidedLimits, policies: Option<[FractionalTubePolicy; 2]>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        Self::build(valve, sources, limits, gate, policies, Some(channel.into()))
    }

    /// The same supply passage together with original regional duct materials,
    /// Maxwell wall histories, optional fractional travel and viscothermal loss.
    /// Every original material receipt and physical body/node index survives.
    /// This constructor delegates to the same material and coupled state owners.
    #[allow(clippy::too_many_arguments)]
    pub fn with_material_walls_and_supply_channel(
        valve: LoadedAperture, sources: [ApertureDuctSpec; 2],
        walls: [Option<MaterialWallInput>; 2], channel: impl Into<FixedChannel>,
        limits: TwoSidedLimits, policies: Option<[FractionalTubePolicy; 2]>,
        gate: &CancelGate,
    ) -> Result<Self, AcousticRealizeError> {
        Self::build_material_walls(valve, sources, walls, limits, policies, gate, Some(channel.into()))
    }

    /// Actual accepted transverse profile, source gas, geometry and energy ledger.
    #[must_use]
    pub const fn supply_channel(&self) -> Option<&FixedChannel> {
        self.supply_channel.as_ref()
    }
}
