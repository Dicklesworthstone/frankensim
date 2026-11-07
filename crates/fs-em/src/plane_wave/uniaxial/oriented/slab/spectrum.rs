//! Sparse finite-slab illumination, reusing the original coherent-field sum.
use super::super::super::kernel::unit_axis;
use super::super::kernel::FieldSum;
use super::super::spectrum::{UniaxialIncidentMode, UniaxialSpectrumField, UniaxialSpectrumLimits};
use super::super::{Incidence, ReflectionError};
use super::{OrientedUniaxialSlab, SlabRegion, UniaxialSlabResponse, kernel};
use crate::plane_wave::angular_spectrum::AngularSpectrumError;
use fs_math::c64::C64;

type Result<T> = core::result::Result<T, AngularSpectrumError>;

struct Mode {
    axis: [f64; 2],
    response: UniaxialSlabResponse,
}

/// One source-bound coherent angular spectrum through the complete finite
/// crystal. Original spectral laws are evaluated once, and each occupied mode
/// gets one two-face solve against the SAME global optic axis. Queries retain
/// four-wave internal interference and all intermode field cross terms.
pub struct PreparedUniaxialSlabSpectrum<'s, 'i, 'm> {
    source: &'s OrientedUniaxialSlab<'i, 'm>,
    omega_rad_s: f64,
    modes: Vec<Mode>,
    workspace_bytes: usize,
}

impl<'s, 'i, 'm> PreparedUniaxialSlabSpectrum<'s, 'i, 'm> {
    /// Same authored mode/phase convention as the half-space spectrum. No
    /// nonzero coefficient is filtered; original support and all field/power
    /// gates must pass before the complete preparation is returned.
    pub fn prepare(
        source: &'s OrientedUniaxialSlab<'i, 'm>,
        omega_rad_s: f64,
        input: &[UniaxialIncidentMode],
        limits: UniaxialSpectrumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let workspace_bytes = admit(input.len(), limits)?;
        let data = source.resolve_frequency(omega_rad_s)?;
        let modes = prepare_modes(data, input, &mut cancelled)?;
        Ok(Self {
            source,
            omega_rad_s: data.crystal.omega,
            modes,
            workspace_bytes,
        })
    }

    pub const fn source(&self) -> &'s OrientedUniaxialSlab<'i, 'm> {
        self.source
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    pub fn retained_mode_count(&self) -> usize {
        self.modes.len()
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Physical [x,y,z] [m] relative to the original parallel faces. This
    /// allocates nothing and performs no material or boundary reevaluation.
    pub fn field_at(
        &self,
        region: SlabRegion,
        position_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<UniaxialSpectrumField> {
        field_at(
            &self.modes,
            self.source.thickness_m(),
            region,
            position_m,
            &mut cancelled,
        )
    }
}

fn bad(message: &'static str) -> AngularSpectrumError {
    AngularSpectrumError::Invalid(message)
}
fn poll(cancelled: &mut impl FnMut() -> bool) -> Result<()> {
    if cancelled() {
        Err(AngularSpectrumError::Cancelled)
    } else {
        Ok(())
    }
}
fn admit(count: usize, limits: UniaxialSpectrumLimits) -> Result<usize> {
    if count > limits.maximum_modes {
        return Err(bad("uniaxial slab exceeds its mode/work budget"));
    }
    count
        .checked_mul(core::mem::size_of::<Mode>())
        .and_then(|v| v.checked_add(core::mem::size_of::<FieldSum>()))
        .and_then(|v| v.checked_add(core::mem::size_of::<UniaxialSpectrumField>()))
        .filter(|&v| v <= limits.maximum_workspace_bytes)
        .ok_or_else(|| bad("uniaxial slab exceeds its workspace budget"))
}
fn prepare_modes(
    data: kernel::Input,
    input: &[UniaxialIncidentMode],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Mode>> {
    let mut modes = Vec::new();
    modes
        .try_reserve_exact(input.len())
        .map_err(|_| bad("uniaxial slab mode allocation failed"))?;
    for authored in input {
        poll(cancelled)?;
        if !authored.tangential_phase_index.is_finite()
            || authored
                .electric_v_m
                .iter()
                .any(|v| !v.re.is_finite() || !v.im.is_finite())
        {
            return Err(bad("uniaxial slab incident modes must be finite"));
        }
        let axis = unit_axis(authored.incidence_plane_axis).map_err(bad)?;
        if authored.electric_v_m == [C64::ZERO; 2] {
            continue;
        }
        let [x, y, z] = data.axis;
        let local = [x * axis[0] + y * axis[1], y * axis[0] - x * axis[1], z];
        let [ex, ey] = authored.electric_v_m;
        let electric = [
            ex.scale(axis[0]) + ey.scale(axis[1]),
            ey.scale(axis[0]) - ex.scale(axis[1]),
        ];
        let mut response = kernel::solve(
            kernel::Input {
                axis: local,
                ..data
            },
            Incidence::TangentialIndex(authored.tangential_phase_index),
            electric,
        )
        .map_err(|message| AngularSpectrumError::Source(ReflectionError::Invalid(message)))?;
        response.rotate(axis);
        modes.push(Mode { axis, response });
    }
    Ok(modes)
}
fn field_at(
    modes: &[Mode],
    depth: f64,
    region: SlabRegion,
    point: [f64; 3],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<UniaxialSpectrumField> {
    poll(cancelled)?;
    kernel::admit_point(region, point, depth)?;
    let mut sum = FieldSum::at_point(point).map_err(bad)?;
    for mode in modes {
        poll(cancelled)?;
        mode.response.accumulate(&mut sum, region, mode.axis)?;
    }
    let field = sum.finish().map_err(bad)?;
    Ok(UniaxialSpectrumField {
        electric_v_m: field.electric_v_m,
        magnetic_a_m: field.magnetic_a_m,
        poynting_w_m2: field.poynting_w_m2,
    })
}

#[cfg(test)]
mod tests;
