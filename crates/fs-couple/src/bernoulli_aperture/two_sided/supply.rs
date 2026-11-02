//! A fixed physical supply passage in the SAME reed/tract/bore transaction.
//!
//! StaticPressure SlitChannel has the discrete law dp = Z_c Q - p_history.
//! It therefore adds Z_c to the original valve load and p_history to its
//! effective upstream pressure. No second nonlinear solve or feedback delay is
//! needed. The channel is additional supplied geometry, not a resized reed gap.
//! Z_c = rho L / (A dt) + wall friction. The first part is ADDED MASS (inertance
//! per step), not loss, so the valve's clock admission counts only the friction
//! as a rate; the passage steps on the reed's step-mean flow, which leaves one
//! parasitic root -mu/(mu+2) that `admit_parasitic_root` bounds.
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
use fs_phs::linear_exact::MIDPOINT_ONE_CENT_ANGULAR_STEP;

/// Accepted supply-port exchange, kept separate from the moving reed gap.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct SupplyChannelFrame {
    /// Actual pressure lost between source-side and upstream reed face [Pa].
    /// It may be negative when stored channel momentum returns work.
    pub pressure_drop_pa: f64,
    /// Original no-slip field's flow, storage, physical and numerical loss.
    pub channel: SlitChannelStep,
}

/// The decay the coupling's parasitic mode must reach within one period of the fastest admitted physical mode:
/// -60 dB, the reverberation-decay convention (T60, ISO 3382-1). THIS IS A CHOICE: a cited convention for "has
/// died away", not a measured or fitted number, and no source says the parasitic mode must meet it. Everything below
/// is derived from it.
pub const SUPPLY_PARASITIC_DECAY: f64 = 1.0e-3;

/// Largest admitted magnitude of the supply coupling's parasitic root, rho_max = 0.956977 (admitted added-mass ratio
/// mu <= 2 rho_max / (1 - rho_max) = 44.5). DERIVED, given the choice [`SUPPLY_PARASITIC_DECAY`], from the one-cent
/// band so the two cannot drift apart: the fastest admitted mode has `MIDPOINT_ONE_CENT_ANGULAR_STEP` = 0.04 radians
/// per step, so one of its periods is `TAU / 0.04` = 157.08 steps, and the parasitic root, which decays as rho^n,
/// must reach `SUPPLY_PARASITIC_DECAY` over exactly that many: rho_max^(TAU / step) = decay, the `pow` below.
#[must_use]
pub fn maximum_parasitic_root() -> f64 {
    fs_math::det::pow(SUPPLY_PARASITIC_DECAY, MIDPOINT_ONE_CENT_ANGULAR_STEP / core::f64::consts::TAU)
}

/// How the passage's discrete impedance Z (dp = Z Q - H) splits for the valve. The part rho L / (A dt) is ADDED MASS:
/// it stores kinetic energy and returns it, and only lowers every mechanical rate. The excess is wall friction.
#[derive(Clone, Copy, Debug, PartialEq)]
pub(crate) struct SupplyLoad {
    pub(crate) gas: GasState,
    /// Duct load plus wall friction: what the moving valve is declared (and band-admitted) with.
    pub(crate) resistive_pa_s_m3: f64,
    /// The passage's added-mass impedance, installed on the valve on top of the declared load.
    pub(crate) inertance_pa_s_m3: f64,
}

pub(crate) fn supply_load(sides: &[ApertureDuctSpec; 2], channel: &FixedChannel)
    -> Result<SupplyLoad, AcousticRealizeError>
{
    let (gas, z) = combined_inlet_state(sides)?;
    if channel.gas() != &gas {
        return Err(invalid("reed supply passage and both physical inlets must share the complete material gas"));
    }
    let impedance = channel.pressure_load()?.impedance_pa_s_m3;
    let total = z + impedance;
    if !total.is_finite() || total <= z {
        return Err(invalid("reed supply passage impedance cannot be represented beside the actual duct load"));
    }
    let inertance = channel.inertial_impedance_pa_s_m3();
    // Viscous diffusion only raises Z above rho L / (A dt); rounding may leave a -1 ulp residue.
    let resistive = z + (impedance - inertance).max(0.0);
    if !(inertance.is_finite() && inertance > 0.0 && resistive.is_finite() && resistive > 0.0) {
        return Err(invalid("reed supply passage inertance cannot be represented beside the actual duct load"));
    }
    Ok(SupplyLoad { gas, resistive_pa_s_m3: resistive, inertance_pa_s_m3: inertance })
}

/// Admit the parasitic mode of the supply coupling. The passage steps on the reed's step-mean flow, so the coupled
/// recursion has three levels and, besides the physical pair, one real root -mu/(mu+2) per step (mu = the passage's
/// added mass against the moving masses, (rho L / A_t) sum(A_i^2 / m_i)). It lies strictly inside the unit circle for
/// every finite mu, so the scheme is stable, but it tends to -1 as the passage gets heavy and then rings at the
/// Nyquist rate for many steps. Its magnitude must not exceed [`maximum_parasitic_root`]. This holds for every
/// propagator: the exact block has the same root and no band admission.
pub(crate) fn admit_parasitic_root(valve: &LoadedAperture, inertance_pa_s_m3: f64)
    -> Result<(), AcousticRealizeError>
{
    let mass_ratio = inertance_pa_s_m3 * valve.time_step_s() * valve.pressure_area_squared_over_mass();
    let root = mass_ratio / (mass_ratio + 2.0);
    let maximum_root = maximum_parasitic_root();
    if mass_ratio.is_finite() && root <= maximum_root {
        Ok(())
    } else {
        Err(AcousticRealizeError::SupplyParasiticRoot { mass_ratio, root, maximum_root })
    }
}

/// Resolve the actual gas and the valve's declared load for a fixed passage: the duct load plus the passage's wall
/// friction. Construct any original single, independent-pair or modal material valve with this impedance before
/// composing it with the same channel; [`TwoSidedAperture::with_supply_channel`] then adds the passage's added mass
/// to the valve's load itself, so the clock admission counts it as inertance and not as loss, and admits the
/// coupling's parasitic mode ([`admit_parasitic_root`]). Both duct inlets and the channel must retain the identical
/// coherent material gas.
///
/// The channel's geometry and step period determine its discrete impedance.
/// Its history pressure is state, never an additional authored blowing command.
/// A nonlinear kinetic-head terminal refuses rather than silently losing its law.
pub fn combined_inlet_state_with_supply(
    sides: &[ApertureDuctSpec; 2], channel: impl Into<FixedChannel>,
) -> Result<(GasState, f64), AcousticRealizeError> {
    let load = supply_load(sides, &channel.into())?;
    Ok((load.gas, load.resistive_pa_s_m3))
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
