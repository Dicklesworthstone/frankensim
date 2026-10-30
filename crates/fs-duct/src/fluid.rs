//! Small-signal duct propagation from the shared gas/liquid acoustic response.
//!
//! Density and adiabatic sound speed come from one original [`FluidAcoustics`]
//! evaluation. Lossy entry points use [`ThermoviscousFluid`] for the same
//! state's viscosity, conductivity, Cp and actual Cp/Cv ratio. The material
//! may be a sourced liquid or a real gas; no ideal-gas relation, viscosity or
//! atmospheric default is manufactured to fit a gas-shaped API.
//!
//! Cylinders and rigid cones use the existing analytic wave bases. Smooth
//! profiles use the existing fourth-order Magnus integrator, and complete
//! chains use the existing stable scattering/composition owner. [`WallModel`]
//! supplies the same local, material-owned impedance as the gas path, at each
//! actual radius and surface slope. Thus liquid compressibility, wall storage
//! and wall dissipation all enter the same physical pressure/volume-flow port.
//!
//! Lossless calls need only rho and c. Thermoviscous calls use the existing
//! no-slip, isothermal-boundary Zwikker–Kosten model, with full Bessel or
//! explicitly admitted wide-tube response. A separate compliant wall or
//! terminal load may also dissipate energy. Mean flow, molecular/bulk-viscous
//! absorption, cavitation, nonlinear water hammer, multiphase acoustics and
//! shell dynamics are outside this model. Gas-only tone-hole radiation is not
//! inferred for a liquid: those branches need their own pressure/flow boundary.
//! The caller retains the source material and its original query receipts.

use crate::compose::{self, ChainResult, ComposeOptions, PreparedChain, TwoPort};
use crate::{ConeSlices, Duct, DuctError, LossModel, Segment, SegmentWave, WallModel};
use fs_material::MaterialError;
use fs_material::fluid_acoustics::thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState};
use fs_material::fluid_acoustics::{FluidAcousticState, FluidAcoustics};
use fs_math::c64::C64;

/// Original material failure or the existing geometry/propagation refusal.
#[derive(Debug, Clone, PartialEq)]
pub enum FluidDuctError {
    /// The material owner could not provide a coherent acoustic state.
    Material(MaterialError),
    /// The native duct owner refused geometry, work, conditioning or a wall.
    Duct(DuctError),
}

impl core::fmt::Display for FluidDuctError {
    fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result {
        match self {
            Self::Material(error) => write!(f, "FS-DUCT-FLUID-MATERIAL: {error}"),
            Self::Duct(error) => write!(f, "{error}"),
        }
    }
}

impl std::error::Error for FluidDuctError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            Self::Material(error) => Some(error),
            Self::Duct(error) => Some(error),
        }
    }
}

impl From<MaterialError> for FluidDuctError {
    fn from(error: MaterialError) -> Self {
        Self::Material(error)
    }
}

impl From<DuctError> for FluidDuctError {
    fn from(error: DuctError) -> Self {
        Self::Duct(error)
    }
}

fn bad(what: &'static str) -> DuctError {
    DuctError::BadParameter { what }
}

fn positive(value: f64) -> bool {
    value.is_finite() && value > 0.0
}

/// Shared numerical owner of the lossless wave, including the existing gas
/// entry point. The snapshot has already admitted rho, c and adiabatic K.
pub(crate) fn wave_from_state(
    state: &FluidAcousticState,
    radius: f64,
    omega: f64,
) -> Result<SegmentWave, DuctError> {
    if !positive(radius) || !positive(omega) {
        return Err(bad(
            "fluid wave requires positive finite radius and frequency",
        ));
    }
    let area = core::f64::consts::PI * radius * radius;
    let k = omega / state.sound_speed_m_s();
    let specific = state.characteristic_impedance_pa_s_m();
    let characteristic = specific / area;
    if !positive(area) || !positive(k) || !positive(characteristic) {
        return Err(bad(
            "fluid wave area, wavenumber and impedance must be representable",
        ));
    }
    Ok(SegmentWave {
        wavenumber: C64::from_re(k),
        characteristic_impedance: C64::from_re(characteristic),
        specific_impedance: C64::from_re(specific),
        shear_number: f64::INFINITY,
    })
}

/// Lossless wave at one actual material state, in SI and exp(-i omega t).
/// A supplied [`fs_material::gas::GasState`] still undergoes its full original
/// coherence admission; a liquid never needs invented gas transport data.
pub fn segment_wave<M: FluidAcoustics + ?Sized>(
    material: &M,
    radius: f64,
    omega: f64,
) -> Result<SegmentWave, FluidDuctError> {
    let state = material.acoustic_state()?;
    Ok(wave_from_state(&state, radius, omega)?)
}

/// Full Bessel or first-order wide-tube wall diffusion from one coherent
/// material state. Cp/Cv is the actual caloric ratio, including for a liquid;
/// no relation between that ratio, pressure and sound speed is imposed.
/// The wide-tube path requires BOTH viscous and coupled thermal boundary
/// layers to be thin, with a small coupled thermal correction. Use
/// [`LossModel::Bessel`] outside that approximation.
pub fn segment_wave_with_loss<M: ThermoviscousFluid + ?Sized>(
    material: &M,
    radius: f64,
    omega: f64,
    loss: LossModel,
) -> Result<SegmentWave, FluidDuctError> {
    let state = material.thermoviscous_state()?;
    Ok(thermoviscous_wave_from_state(&state, radius, omega, loss)?)
}

fn thermoviscous_wave_from_state(
    state: &ThermoviscousFluidState,
    radius: f64,
    omega: f64,
    loss: LossModel,
) -> Result<SegmentWave, DuctError> {
    let acoustic = state.acoustic();
    let base = wave_from_state(&acoustic, radius, omega)?;
    if loss == LossModel::Lossless {
        return Ok(base);
    }
    let transport = state.transport();
    let inputs = crate::WallDiffusionInputs {
        density: acoustic.density_kg_m3(),
        sound_speed: acoustic.sound_speed_m_s(),
        dynamic_viscosity: transport.dynamic_viscosity_pa_s(),
        gamma: state.heat_capacity_ratio(),
        prandtl: transport.prandtl(),
    };
    let wave = match loss {
        LossModel::WideTube => crate::wide_tube_wave_from_inputs(
            &inputs,
            radius,
            omega,
            base.wavenumber.re,
            base.characteristic_impedance.re,
            true,
        )?,
        LossModel::Bessel => crate::bessel_wave_from_inputs(&inputs, radius, omega)?,
        LossModel::Lossless => unreachable!("lossless response returned above"),
    };
    if [
        wave.wavenumber,
        wave.characteristic_impedance,
        wave.specific_impedance,
    ]
    .iter()
    .any(|value| !value.re.is_finite() || !value.im.is_finite())
    {
        return Err(bad("thermoviscous fluid wave is not representable"));
    }
    Ok(wave)
}

fn poll(cancelled: &mut impl FnMut() -> bool, completed_steps: usize) -> Result<(), DuctError> {
    if cancelled() {
        Err(DuctError::Cancelled { completed_steps })
    } else {
        Ok(())
    }
}

fn admit_blocks(current: usize, additional: usize, maximum: usize) -> Result<(), DuctError> {
    if additional > maximum.saturating_sub(current) {
        Err(bad("fluid duct composition block work cap exceeded"))
    } else {
        Ok(())
    }
}

/// Prepare the original pressure/volume-flow scattering blocks from one
/// coherent, uniform material state. The prepared two-port can be connected
/// to an explicitly supplied load,
/// another chain, or the existing graph/port algebra.
///
/// The wall is sampled at each true local radius. `cone_slices` is the
/// refinable axis for compliant cones; a rigid, lossless cone is exact in one
/// block. Smooth profiles retain their authored [`crate::RiccatiSteps`].
/// Cancellation is polled before material admission and every segment/slice
/// or profile step. No partial chain is returned on any failure.
#[allow(clippy::too_many_arguments)]
pub fn duct_chain<M: FluidAcoustics + ?Sized>(
    duct: &Duct,
    material: &M,
    omega: f64,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedChain, FluidDuctError> {
    admit_geometry(duct, options, &mut cancelled)?;
    let state = material.acoustic_state()?;
    Ok(chain_from_waves(
        duct,
        &state,
        omega,
        wall,
        cone_slices,
        options,
        true,
        &mut cancelled,
        &mut |radius| wave_from_state(&state, radius, omega),
    )?)
}

/// The same geometry and stable composer with material-driven viscosity and
/// heat diffusion. Lossy cones use the caller's `cone_slices` even with rigid
/// walls, because boundary-layer response varies with local radius. Profiles
/// evaluate the same material law at both original Magnus quadrature radii.
#[allow(clippy::too_many_arguments)]
pub fn duct_chain_with_loss<M: ThermoviscousFluid + ?Sized>(
    duct: &Duct,
    material: &M,
    omega: f64,
    loss: LossModel,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<PreparedChain, FluidDuctError> {
    admit_geometry(duct, options, &mut cancelled)?;
    let state = material.thermoviscous_state()?;
    Ok(chain_from_waves(
        duct,
        &state.acoustic(),
        omega,
        wall,
        cone_slices,
        options,
        loss == LossModel::Lossless,
        &mut cancelled,
        &mut |radius| thermoviscous_wave_from_state(&state, radius, omega, loss),
    )?)
}

fn admit_geometry(
    duct: &Duct,
    options: ComposeOptions,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<(), DuctError> {
    options.validate()?;
    poll(cancelled, 0)?;
    if duct.segments.is_empty() {
        return Err(DuctError::EmptyDuct);
    }
    admit_blocks(0, duct.segments.len(), options.maximum_blocks)?;
    for segment in &duct.segments {
        segment.validate()?;
        if matches!(segment, Segment::ToneHole { .. }) {
            return Err(bad(
                "generic fluid branches need an explicit pressure/flow load; gas tone-hole radiation is not implied",
            ));
        }
    }
    Ok(())
}

#[allow(clippy::too_many_arguments)]
fn chain_from_waves(
    duct: &Duct,
    state: &FluidAcousticState,
    omega: f64,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    lossless_bulk: bool,
    cancelled: &mut impl FnMut() -> bool,
    wave_at: &mut impl FnMut(f64) -> Result<SegmentWave, DuctError>,
) -> Result<PreparedChain, DuctError> {
    let first = duct.segments.first().ok_or(DuctError::EmptyDuct)?;
    let reference = wave_from_state(state, first.inlet_radius(), omega)?
        .characteristic_impedance
        .re;
    let mut blocks = Vec::new();
    let mut min_shear_number = f64::INFINITY;
    for (index, segment) in duct.segments.iter().enumerate() {
        poll(cancelled, index)?;
        match *segment {
            Segment::Profile { ref profile, steps } => {
                admit_blocks(blocks.len(), steps.0, options.maximum_blocks)?;
                let profile_reference = wave_from_state(state, profile.inlet_radius(), omega)?
                    .characteristic_impedance
                    .re;
                let (local, shear) = crate::profile::profile_blocks_with_waves_controlled(
                    profile,
                    steps,
                    profile_reference,
                    !wall.is_rigid(),
                    cancelled,
                    &mut |radius, slant| wall.apply(wave_at(radius)?, radius, omega, slant),
                )?;
                min_shear_number = min_shear_number.min(shear);
                for matrix in local {
                    blocks.push(TwoPort::from_abcd(matrix, reference)?);
                }
            }
            Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } if !lossless_bulk || !wall.is_rigid() => {
                let count = cone_slices.count(inlet_radius, outlet_radius)?;
                admit_blocks(blocks.len(), count, options.maximum_blocks)?;
                let dx = length / count as f64;
                for i in 0..count {
                    poll(cancelled, i)?;
                    let ra =
                        inlet_radius + (outlet_radius - inlet_radius) * (i as f64 / count as f64);
                    let rb = inlet_radius
                        + (outlet_radius - inlet_radius) * ((i + 1) as f64 / count as f64);
                    let radius = f64::midpoint(ra, rb);
                    let slope = (rb - ra) / dx;
                    let wave = wall.apply(
                        wave_at(radius)?,
                        radius,
                        omega,
                        (1.0 + slope * slope).sqrt(),
                    )?;
                    min_shear_number = min_shear_number.min(wave.shear_number);
                    blocks.push(TwoPort::from_wave(
                        &Segment::Cone {
                            inlet_radius: ra,
                            outlet_radius: rb,
                            length: dx,
                        },
                        &wave,
                        reference,
                    )?);
                }
            }
            Segment::Cylinder { .. } | Segment::Cone { .. } => {
                admit_blocks(blocks.len(), 1, options.maximum_blocks)?;
                let radius = segment.mean_radius();
                let wave = wall.apply(wave_at(radius)?, radius, omega, 1.0)?;
                min_shear_number = min_shear_number.min(wave.shear_number);
                blocks.push(TwoPort::from_wave(segment, &wave, reference)?);
            }
            Segment::ToneHole { .. } => unreachable!("branch geometry refused before assembly"),
        }
    }
    Ok(PreparedChain {
        blocks,
        reference_impedance: reference,
        min_shear_number,
    })
}

/// Input impedance and actual outlet pressure/volume-flow transfer for an
/// explicitly supplied passive load, using the existing stable composer.
/// `None` is a rigid cap and `Some(C64::ZERO)` is ideal pressure release.
/// An impedance may come directly from a material-driven BEM/radiation owner
/// or [`crate::TabulatedLoad::z_at`], without converting the medium to a gas.
///
/// A closed outlet returns its actual pressure transfer and zero volume flow,
/// not a radiation sentinel.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_load<M: FluidAcoustics + ?Sized>(
    duct: &Duct,
    material: &M,
    omega: f64,
    load: Option<C64>,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ChainResult, FluidDuctError> {
    admit_load(load)?;
    let chain = duct_chain(
        duct,
        material,
        omega,
        wall,
        cone_slices,
        options,
        &mut cancelled,
    )?;
    Ok(compose::compose_chain_controlled(
        &chain.blocks,
        load,
        options,
        &mut cancelled,
    )?)
}

/// Native pressure/flow transmission with full thermoviscous wall diffusion
/// and an explicit passive terminal load. Loss ownership is unchanged:
/// Zwikker–Kosten supplies fluid boundary diffusion, `wall` supplies its own
/// mechanical dissipation, and the load supplies terminal absorption.
#[allow(clippy::too_many_arguments)]
pub fn input_impedance_load_with_loss<M: ThermoviscousFluid + ?Sized>(
    duct: &Duct,
    material: &M,
    omega: f64,
    loss: LossModel,
    load: Option<C64>,
    wall: WallModel<'_>,
    cone_slices: ConeSlices,
    options: ComposeOptions,
    mut cancelled: impl FnMut() -> bool,
) -> Result<ChainResult, FluidDuctError> {
    admit_load(load)?;
    let chain = duct_chain_with_loss(
        duct,
        material,
        omega,
        loss,
        wall,
        cone_slices,
        options,
        &mut cancelled,
    )?;
    Ok(compose::compose_chain_controlled(
        &chain.blocks,
        load,
        options,
        &mut cancelled,
    )?)
}

fn admit_load(load: Option<C64>) -> Result<(), DuctError> {
    if load.is_some_and(|z| !z.re.is_finite() || !z.im.is_finite() || z.re < 0.0) {
        return Err(bad(
            "fluid duct terminal impedance must be finite and passive",
        ));
    }
    Ok(())
}
