//! Small-signal duct propagation from the shared gas/liquid acoustic response.
//!
//! Density and adiabatic sound speed come from one original [`FluidAcoustics`]
//! evaluation. The actual material may be a sourced liquid, an IF97 water
//! state, or a real gas; no gas heat-capacity ratio, viscosity or atmospheric
//! default is manufactured to fit that material into a gas-shaped API.
//!
//! Cylinders and rigid cones use the existing analytic wave bases. Smooth
//! profiles use the existing fourth-order Magnus integrator, and complete
//! chains use the existing stable scattering/composition owner. [`WallModel`]
//! supplies the same local, material-owned impedance as the gas path, at each
//! actual radius and surface slope. Thus liquid compressibility, wall storage
//! and wall dissipation all enter the same physical pressure/volume-flow port.
//!
//! Bulk propagation here is explicitly lossless and has no mean flow. A wall
//! or supplied terminal impedance may dissipate energy. This does not assert
//! liquid viscothermal losses, cavitation, nonlinear water hammer, multiphase
//! acoustics or a shell model. Gas-only tone-hole radiation is not inferred for
//! a liquid: those branches need their own explicit pressure/flow boundary.
//! The caller retains the source material and its original query receipts.

use crate::compose::{self, ChainResult, ComposeOptions, PreparedChain, TwoPort};
use crate::{ConeSlices, Duct, DuctError, Segment, SegmentWave, WallModel};
use fs_material::MaterialError;
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
    options.validate()?;
    poll(&mut cancelled, 0)?;
    let first = duct.segments.first().ok_or(DuctError::EmptyDuct)?;
    admit_blocks(0, duct.segments.len(), options.maximum_blocks)?;
    for segment in &duct.segments {
        segment.validate()?;
        if matches!(segment, Segment::ToneHole { .. }) {
            return Err(bad("generic fluid branches need an explicit pressure/flow load; gas tone-hole radiation is not implied").into());
        }
    }
    let state = material.acoustic_state()?;
    let reference = wave_from_state(&state, first.inlet_radius(), omega)?
        .characteristic_impedance
        .re;
    let mut blocks = Vec::new();
    for (index, segment) in duct.segments.iter().enumerate() {
        poll(&mut cancelled, index)?;
        match *segment {
            Segment::Profile { ref profile, steps } => {
                admit_blocks(blocks.len(), steps.0, options.maximum_blocks)?;
                let profile_reference = wave_from_state(&state, profile.inlet_radius(), omega)?
                    .characteristic_impedance
                    .re;
                let (local, _) = crate::profile::profile_blocks_with_waves_controlled(
                    profile,
                    steps,
                    profile_reference,
                    !wall.is_rigid(),
                    &mut cancelled,
                    &mut |radius, slant| {
                        wall.apply(
                            wave_from_state(&state, radius, omega)?,
                            radius,
                            omega,
                            slant,
                        )
                    },
                )?;
                for matrix in local {
                    blocks.push(TwoPort::from_abcd(matrix, reference)?);
                }
            }
            Segment::Cone {
                inlet_radius,
                outlet_radius,
                length,
            } if !wall.is_rigid() => {
                let count = cone_slices.count(inlet_radius, outlet_radius)?;
                admit_blocks(blocks.len(), count, options.maximum_blocks)?;
                let dx = length / count as f64;
                for i in 0..count {
                    poll(&mut cancelled, i)?;
                    let ra =
                        inlet_radius + (outlet_radius - inlet_radius) * (i as f64 / count as f64);
                    let rb = inlet_radius
                        + (outlet_radius - inlet_radius) * ((i + 1) as f64 / count as f64);
                    let radius = f64::midpoint(ra, rb);
                    let slope = (rb - ra) / dx;
                    let wave = wall.apply(
                        wave_from_state(&state, radius, omega)?,
                        radius,
                        omega,
                        (1.0 + slope * slope).sqrt(),
                    )?;
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
                let wave =
                    wall.apply(wave_from_state(&state, radius, omega)?, radius, omega, 1.0)?;
                blocks.push(TwoPort::from_wave(segment, &wave, reference)?);
            }
            Segment::ToneHole { .. } => unreachable!("branch geometry refused before assembly"),
        }
    }
    Ok(PreparedChain {
        blocks,
        reference_impedance: reference,
        min_shear_number: f64::INFINITY,
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
    if load.is_some_and(|z| !z.re.is_finite() || !z.im.is_finite() || z.re < 0.0) {
        return Err(bad("fluid duct terminal impedance must be finite and passive").into());
    }
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
