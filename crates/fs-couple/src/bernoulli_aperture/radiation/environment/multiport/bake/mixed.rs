//! One reciprocal Rayleigh image of rectangular driven mouths and circular
//! terminals. Leading channels follow the declared rectangle order; remaining
//! channels retain the original terminal order. Each field is normal velocity per physical outward
//! volume flow, so its work pairing is pressure times volume flow.
//!
//! Geometry uses the EXACT existing receiver quadratures, and every matrix entry
//! uses the original power-consistent Rayleigh owner. The unchanged common bake
//! retains source causality, three orders, held-out errors and both passivity
//! admissions. These finite surface/frequency checks do not certify continuum
//! error or turn the original jet source into resolved compressible flow.

use super::{
    BaffledMimoBakeError, BaffledMimoBakeOptions, BaffledMimoBakeReport,
    MimoRadiationLimits, admit_apertures, admit_common, declaration, fit_scratch_bytes,
    fit_source, invalid, poll, realization, scalar,
};
use crate::pcm_wav::baffled::{BaffledPressure, CircularOutletReceiver};
use crate::pcm_wav::baffled::openings::{CircularOpening, RectangularOutletReceiver};
use fs_bem::helmholtz::{Medium, aperture::reciprocity_and_passivity};
use fs_bem::helmholtz::rayleigh_modal::{
    AreaSample, PreparedSurface, MAX_MODAL_RAYLEIGH_CELLS, rectangular_static_inertance,
};
use fs_blake3::DomainHasher;
use fs_exec::CancelGate;
use fs_material::gas::GasState;
use fs_math::{c64::C64, det};
use fs_vfit::mimo::MimoModel;

/// Bound the complete cold source sweep, not just one matrix evaluation.
pub const MAX_MIXED_BAFFLED_SWEEP_WORK: u128 = 20_000_000_000;

/// Actual rectangular opening and explicit quadrature of every source channel.
/// Geometry is fixed in the same z=0 infinite rigid baffle as the terminals.
/// Source admission resolves the FULL mechanical band, independently of the
/// microphone's narrower observation band.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct DrivenMouthSurface {
    pub center_m: [f64; 2],
    pub size_m: [f64; 2],
    pub cells: [usize; 2],
    pub terminal_radial_rings: usize,
    pub terminal_angular_points: usize,
}

struct MixedSurface {
    surface: PreparedSurface,
    ports: usize,
    scratch_bytes: usize,
}

impl MixedSurface {
    fn new(
        mouth: DrivenMouthSurface, terminals: &[CircularOpening], gas: &GasState,
        band_hz: [f64; 2], options: BaffledMimoBakeOptions, gate: &CancelGate,
    ) -> Result<Self, BaffledMimoBakeError> {
        Self::new_many(&[mouth], terminals, gas, band_hz, options, gate)
    }

    fn new_many(
        mouths: &[DrivenMouthSurface], terminals: &[CircularOpening], gas: &GasState,
        band_hz: [f64; 2], options: BaffledMimoBakeOptions, gate: &CancelGate,
    ) -> Result<Self, BaffledMimoBakeError> {
        poll(gate)?;
        let ports = mouths.len().checked_add(terminals.len())
            .ok_or_else(|| invalid("mixed source channel count overflows"))?;
        if mouths.is_empty() || terminals.is_empty() || ports > super::MAX_BAFFLED_MIMO_PORTS
            || mouths.iter().any(|mouth|
                mouth.center_m.iter().any(|value| !value.is_finite())
                || mouth.size_m.iter().any(|value| !value.is_finite() || *value <= 0.0)
                || mouth.cells.iter().any(|count| !(1..=512).contains(count))
                || !(1..=128).contains(&mouth.terminal_radial_rings)
                || !(8..=512).contains(&mouth.terminal_angular_points))
        {
            return Err(invalid("mixed source requires bounded rectangles and original circular quadrature"));
        }
        let quadrature = mouths[0];
        if mouths.iter().any(|mouth|
            mouth.terminal_radial_rings != quadrature.terminal_radial_rings
                || mouth.terminal_angular_points != quadrature.terminal_angular_points)
        {
            return Err(invalid("joint mouths must declare the same shared circular quadrature"));
        }
        admit_apertures(terminals)?;
        for (index, mouth) in mouths.iter().enumerate() {
            poll(gate)?;
            for previous in &mouths[..index] {
                let dx = (mouth.center_m[0] - previous.center_m[0]).abs()
                    - (0.5 * mouth.size_m[0] + 0.5 * previous.size_m[0]);
                let dy = (mouth.center_m[1] - previous.center_m[1]).abs()
                    - (0.5 * mouth.size_m[1] + 0.5 * previous.size_m[1]);
                if !dx.is_finite() || !dy.is_finite() || (dx <= 0.0 && dy <= 0.0) {
                    return Err(invalid("joint rectangular mouths must be disjoint and separated"));
                }
            }
            for opening in terminals {
                let dx = ((opening.center_m[0] - mouth.center_m[0]).abs()
                    - 0.5 * mouth.size_m[0]).max(0.0);
                let dy = ((opening.center_m[1] - mouth.center_m[1]).abs()
                    - 0.5 * mouth.size_m[1]).max(0.0);
                if !dx.is_finite() || !dy.is_finite() || det::hypot(dx, dy) <= opening.radius_m {
                    return Err(invalid("mixed mouth and terminal openings must be disjoint and separated"));
                }
            }
        }
        let rectangle_cells = mouths.iter().try_fold(0usize, |sum, mouth|
            sum.checked_add(mouth.cells[0] * mouth.cells[1]))
            .ok_or_else(|| invalid("mixed rectangle surface count overflows"))?;
        let circle_cells = quadrature.terminal_radial_rings * quadrature.terminal_angular_points;
        let count = rectangle_cells.checked_add(circle_cells * terminals.len())
            .ok_or_else(|| invalid("mixed source surface count overflows"))?;
        if count > MAX_MODAL_RAYLEIGH_CELLS {
            return Err(invalid("mixed source exceeds the original Rayleigh surface-cell allowance"));
        }
        let cells = count as u128;
        let fields = ports as u128;
        // Training, withheld and causality rows plus the independent static row.
        let work = 3 * options.frequency_samples as u128
            * cells * fields * (cells + fields);
        if work > MAX_MIXED_BAFFLED_SWEEP_WORK {
            return Err(invalid("mixed Rayleigh source exceeds its complete sweep-work allowance"));
        }
        // Pair distances, complete sample fields, original pressure scratch and
        // the shared fitter coexist. Check their combined declared allowance
        // before either receiver surface allocates its first cell.
        let surface_bytes = 8 * (cells * (cells - 1) / 2 + cells * (12 + 8 * fields));
        let required = fit_scratch_bytes(ports, options) as u128 + surface_bytes;
        if required > options.max_fit_memory_bytes as u128 {
            return Err(BaffledMimoBakeError::Memory { kind: "mixed surface and fit scratch",
                required: usize::try_from(required).unwrap_or(usize::MAX),
                allowed: options.max_fit_memory_bytes });
        }
        let medium = gas.into();
        // Receiver position is not used by these surface-only owners. The real
        // microphone remains independently constructed at its authored position.
        let mut samples = Vec::with_capacity(count);
        for (index, mouth) in mouths.iter().enumerate() {
            poll(gate)?;
            let rectangle = BaffledPressure::rectangular_surface(mouth.size_m,
                RectangularOutletReceiver { center_m: mouth.center_m, position_m: [0.0; 3],
                    cells: mouth.cells, maximum_frequency_hz: band_hz[1] }, medium)
                .map_err(BaffledMimoBakeError::Surface)?;
            for cell in rectangle {
                let mut shape = vec![0.0; ports];
                shape[index] = cell.mode_shape[0];
                samples.push(AreaSample { position_m: [cell.position_m[0], cell.position_m[1]],
                    area_m2: cell.area_m2, mode_shape: shape });
            }
        }
        for (index, opening) in terminals.iter().enumerate() {
            poll(gate)?;
            let circle = BaffledPressure::circular_surface(opening.radius_m,
                CircularOutletReceiver { position_m: [0.0; 3],
                    radial_rings: quadrature.terminal_radial_rings,
                    angular_points: quadrature.terminal_angular_points,
                    maximum_frequency_hz: band_hz[1] }, medium)
                .map_err(BaffledMimoBakeError::Surface)?;
            for cell in circle {
                let mut shape = vec![0.0; ports];
                shape[index + mouths.len()] = cell.mode_shape[0];
                samples.push(AreaSample {
                    position_m: [cell.position_m[0] + opening.center_m[0],
                        cell.position_m[1] + opening.center_m[1]],
                    area_m2: cell.area_m2, mode_shape: shape,
                });
            }
        }
        let surface = PreparedSurface::new(&samples)
            .map_err(|error| BaffledMimoBakeError::Surface(error.to_string()))?;
        // Retain the compact mouth's independent five-percent static-shape
        // check, and apply the original circular piston limit to every disk.
        // No self entry is repaired or replaced by this independent reference.
        let omega_small = core::f64::consts::TAU * band_hz[1] / 1_048_576.0;
        let image = surface.power_consistent_impedance(omega_small, Medium::from(gas))
            .map_err(|source| BaffledMimoBakeError::Sample {
                frequency_hz: omega_small / core::f64::consts::TAU, source })?;
        for port in 0..ports {
            let exact = if port < mouths.len() {
                rectangular_static_inertance(mouths[port].size_m, gas.density)
                    .map_err(|error| BaffledMimoBakeError::Surface(error.to_string()))?
            } else {
                8.0 * gas.density / (3.0 * core::f64::consts::PI.powi(2)
                    * terminals[port - mouths.len()].radius_m)
            };
            let relative = (-image[port * ports + port].im / omega_small / exact - 1.0).abs();
            if !relative.is_finite() || relative > 0.05 {
                return Err(invalid("mixed source exceeds five-percent exact static aperture-inertia error; refine its cells"));
            }
        }
        Ok(Self { surface, ports, scratch_bytes: surface_bytes as usize })
    }

    fn rows(&self, gas: &GasState, hz: &[f64], gate: &CancelGate)
        -> Result<Vec<Vec<C64>>, BaffledMimoBakeError>
    {
        let mut rows = Vec::with_capacity(hz.len());
        for &frequency_hz in hz {
            poll(gate)?;
            let row: Vec<_> = self.surface.power_consistent_impedance(
                core::f64::consts::TAU * frequency_hz, Medium::from(gas))
                .map_err(|source| BaffledMimoBakeError::Sample { frequency_hz, source })?
                .iter().map(|value| value.conj()).collect();
            // The surface fields already include 1/area once on each side.
            // Only the exp(-iwt) to s=+iw convention changes here.
            let (reciprocity, margin) = reciprocity_and_passivity(&row, self.ports);
            if row.iter().any(|value| !value.re.is_finite() || !value.im.is_finite())
                || !reciprocity.is_finite() || !margin.is_finite()
                || reciprocity > 1e-12 || margin < -1e-12
            {
                return Err(invalid("original mixed Rayleigh matrix fails reciprocity or sampled passivity"));
            }
            rows.push(row);
        }
        Ok(rows)
    }
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn bake_driven_mouth(
    exterior: &GasState, mouth: DrivenMouthSurface, terminals: &[CircularOpening],
    mechanics_rate_hz: u32, band_hz: [f64; 2], options: BaffledMimoBakeOptions,
    limits: MimoRadiationLimits, gate: &CancelGate,
) -> Result<(MimoModel, BaffledMimoBakeReport, [u8; 32], [u8; 32]), BaffledMimoBakeError> {
    bake_driven_mouths(exterior, &[mouth], terminals, mechanics_rate_hz,
        band_hz, options, limits, gate)
}

#[allow(clippy::too_many_arguments)]
pub(crate) fn bake_driven_mouths(
    exterior: &GasState, mouths: &[DrivenMouthSurface], terminals: &[CircularOpening],
    mechanics_rate_hz: u32, band_hz: [f64; 2], options: BaffledMimoBakeOptions,
    limits: MimoRadiationLimits, gate: &CancelGate,
) -> Result<(MimoModel, BaffledMimoBakeReport, [u8; 32], [u8; 32]), BaffledMimoBakeError> {
    poll(gate)?;
    let ports = mouths.len().checked_add(terminals.len())
        .ok_or_else(|| invalid("mixed source channel count overflows"))?;
    admit_common(exterior, ports, mechanics_rate_hz, band_hz, options, limits)?;
    let surface = if mouths.len() == 1 {
        MixedSurface::new(mouths[0], terminals, exterior, band_hz, options, gate)?
    } else {
        MixedSurface::new_many(mouths, terminals, exterior, band_hz, options, gate)?
    };
    // The retained surface coexists with both fitting and residue repair.
    // Keep its allocation out of the shared owner's remaining allowance.
    let fit_options = BaffledMimoBakeOptions {
        max_fit_memory_bytes: options.max_fit_memory_bytes - surface.scratch_bytes,
        ..options
    };
    let (model, report) = fit_source(ports, mechanics_rate_hz, band_hz,
        fit_options, limits, gate, |hz, gate| surface.rows(exterior, hz, gate))?;
    // Preserve the original singleton declaration and realization identities.
    let mut hash = DomainHasher::new(if mouths.len() == 1 {
        "org.frankensim.fs-couple.mixed-baffled-flow-source.v1"
    } else { "org.frankensim.fs-couple.multiple-mixed-baffled-flow-source.v1" });
    hash.update(&declaration(exterior, terminals, mechanics_rate_hz, band_hz));
    if mouths.len() != 1 { hash.update(&(mouths.len() as u64).to_le_bytes()); }
    for mouth in mouths {
        for value in mouth.center_m.into_iter().chain(mouth.size_m) { scalar(&mut hash, value); }
        for count in [mouth.cells[0], mouth.cells[1], mouth.terminal_radial_rings,
            mouth.terminal_angular_points] { hash.update(&(count as u64).to_le_bytes()); }
    }
    let source_identity = hash.finalize().0;
    let realization_identity = realization(source_identity, &model, options);
    Ok((model, report, source_identity, realization_identity))
}

#[cfg(test)]
mod tests;
