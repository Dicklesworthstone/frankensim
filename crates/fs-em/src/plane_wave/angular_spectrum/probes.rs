//! Direct Maxwell fields at nonuniform physical observations.
//!
//! The original modal amplitudes are summed before any quadratic observable.
//! No FFT interpolation, paraxial approximation or sampled-intensity mixture
//! is used. The periodic envelope and its physical carrier remain explicit.

use super::*;
use fs_math::dd::Dd;

type ComplexDd = [Dd; 2];
const ZERO: ComplexDd = [Dd::ZERO; 2];

/// Bounds a nonuniform query before output allocation. For A exactly nonzero
/// incident modes and P probes, work is A*P native Maxwell mode evaluations
/// and phase evaluations, plus two linear scans of the original mode array
/// and O(P) carrier/observable work.
/// Storage includes the owner's conservative workspace, A indices, P outputs
/// and one six-component double-double accumulator. Borrowed observations,
/// previously retained outputs and allocator overhead are excluded.
#[derive(Clone, Copy, Debug)]
pub struct AngularSpectrumProbeLimits {
    pub maximum_points: usize,
    pub maximum_mode_evaluations: usize,
    pub maximum_workspace_bytes: usize,
}

/// Complete complex peak fields at one observation. The complex cycle-mean
/// Poynting vector is E cross H* / 2: its real part is signed energy transport,
/// and its imaginary part retains reactive near-field flow. A point has no
/// implicit detector area and these values are not integrated watts.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumProbeField {
    pub electric_v_per_m: [C64; 3],
    pub magnetic_a_per_m: [C64; 3],
    pub complex_poynting_w_m2: [C64; 3],
}

/// Position inside a prepared original film. Depth is measured from that
/// film's first face and includes either boundary. x/y use the original
/// periodic envelope coordinates and may lie outside the original cell.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumLayerObservation {
    pub depth_m: f64,
    pub position_m: [f64; 2],
}

/// Maxwell field and positive local constitutive heat inside the selected
/// film. Heat includes longitudinal fields and coherent spatial interference.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct AngularSpectrumLayerProbe {
    pub field: AngularSpectrumProbeField,
    pub absorbed_power_density_w_m3: f64,
}

impl PlanarAngularSpectrum<'_, '_> {
    /// Resolve arbitrary exterior probes, including irregular or tilted sensor
    /// layouts and different depths/sides in one query. Each result corresponds
    /// to the observation at the same index. All six fields use the same native
    /// mode reconstruction as `field`, retaining evanescent transmitted decay,
    /// exact-critical exits, reflected signs and standing-wave interference.
    ///
    /// Sparse authored spectra visit only exactly nonzero modes. No tolerance
    /// removes a source coefficient. Materials and multilayer responses are
    /// not reevaluated. This direct sum is intended for sparse receivers; the
    /// existing FFT `field` remains faster for a complete uniform plane.
    pub fn probe_fields(
        &self,
        observations: &[AngularSpectrumObservation],
        limits: AngularSpectrumProbeLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Vec<AngularSpectrumProbeField>, AngularSpectrumError> {
        poll(&mut cancelled)?;
        if observations.len() > limits.maximum_points {
            return Err(bad("probe count exceeds its budget"));
        }
        for (index, observation) in observations.iter().enumerate() {
            if index % 256 == 0 {
                poll(&mut cancelled)?;
            }
            admit_position(observation.position_m)?;
            if !observation.distance_from_face_m.is_finite()
                || observation.distance_from_face_m < 0.0
            {
                return Err(bad("probe depth must be finite and nonnegative"));
            }
        }
        let active = admit::<AngularSpectrumProbeField>(
            self,
            observations.len(),
            self.workspace_bytes,
            limits,
            &mut cancelled,
        )?;
        let mut output = Vec::with_capacity(observations.len());
        for observation in observations {
            output.push(synthesize(
                self,
                &active,
                observation.position_m,
                |index| {
                    self.mode_field(
                        &self.modes[index],
                        observation.side,
                        observation.distance_from_face_m,
                    )
                },
                &mut cancelled,
            )?);
        }
        poll(&mut cancelled)?;
        Ok(output)
    }
}

pub(super) fn admit_position(position: [f64; 2]) -> Result<(), AngularSpectrumError> {
    if position.iter().any(|v| !v.is_finite()) {
        return Err(bad("probe transverse coordinates must be finite"));
    }
    Ok(())
}

pub(super) fn admit<T>(
    beam: &PlanarAngularSpectrum<'_, '_>,
    points: usize,
    owner_workspace_bytes: usize,
    limits: AngularSpectrumProbeLimits,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<usize>, AngularSpectrumError> {
    poll(cancelled)?;
    let mut count = 0usize;
    for (index, mode) in beam.modes.iter().enumerate() {
        if index % 256 == 0 {
            poll(cancelled)?;
        }
        if mode.incident_e != [C64::ZERO; 2] {
            count += 1;
        }
    }
    let evaluations = count
        .checked_mul(points)
        .ok_or_else(|| bad("probe evaluation count overflow"))?;
    let bytes = count
        .checked_mul(core::mem::size_of::<usize>())
        .and_then(|v| {
            points
                .checked_mul(core::mem::size_of::<T>())
                .and_then(|p| v.checked_add(p))
        })
        .and_then(|v| v.checked_add(core::mem::size_of::<[ComplexDd; 6]>()))
        .and_then(|v| v.checked_add(owner_workspace_bytes))
        .ok_or_else(|| bad("probe workspace extent overflow"))?;
    if points > limits.maximum_points
        || evaluations > limits.maximum_mode_evaluations
        || bytes > limits.maximum_workspace_bytes
    {
        return Err(bad("probes exceed their work or workspace budget"));
    }
    let mut active = Vec::with_capacity(count);
    for (index, mode) in beam.modes.iter().enumerate() {
        if index % 256 == 0 {
            poll(cancelled)?;
        }
        if mode.incident_e != [C64::ZERO; 2] {
            active.push(index);
        }
    }
    Ok(active)
}

pub(super) fn synthesize(
    beam: &PlanarAngularSpectrum<'_, '_>,
    active: &[usize],
    position: [f64; 2],
    mut field: impl FnMut(usize) -> Result<([C64; 3], [C64; 3]), AngularSpectrumError>,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AngularSpectrumProbeField, AngularSpectrumError> {
    poll(cancelled)?;
    let mut sum = [ZERO; 6];
    for (ordinal, &index) in active.iter().enumerate() {
        if ordinal % 32 == 0 {
            poll(cancelled)?;
        }
        let (e, h) = field(index)?;
        // Empty outgoing fields need no spatial phase even at a distant probe.
        if e == [C64::ZERO; 3] && h == [C64::ZERO; 3] {
            continue;
        }
        let weight = phase(
            [
                beam.envelope_wavenumber(0, index % beam.grid.width),
                beam.envelope_wavenumber(1, index / beam.grid.width),
            ],
            position,
        )?;
        for (component, value) in sum.iter_mut().zip(e.into_iter().chain(h)) {
            let product = multiply([Dd::from_f64(value.re), Dd::from_f64(value.im)], weight);
            component[0] = component[0] + product[0];
            component[1] = component[1] + product[1];
        }
    }
    // A common carrier cancels from E cross H*, so form the complex flux from
    // the complete envelope before applying its phase once to the six fields.
    let mut complex_poynting_w_m2 = [C64::ZERO; 3];
    for (axis, (a, b)) in [(1, 2), (2, 0), (0, 1)].into_iter().enumerate() {
        let x = multiply(sum[a], [sum[b + 3][0], -sum[b + 3][1]]);
        let y = multiply(sum[b], [sum[a + 3][0], -sum[a + 3][1]]);
        complex_poynting_w_m2[axis] = rounded([
            (x[0] - y[0]) * Dd::from_f64(0.5),
            (x[1] - y[1]) * Dd::from_f64(0.5),
        ])?;
    }
    if beam.carrier_rad_m != [0.0; 2] && sum != [ZERO; 6] {
        let carrier = phase(beam.carrier_rad_m, position)?;
        for value in &mut sum {
            *value = multiply(*value, carrier);
        }
    }
    let mut fields = [C64::ZERO; 6];
    for (value, total) in fields.iter_mut().zip(sum) {
        *value = rounded(total)?;
    }
    poll(cancelled)?;
    Ok(AngularSpectrumProbeField {
        electric_v_per_m: [fields[0], fields[1], fields[2]],
        magnetic_a_per_m: [fields[3], fields[4], fields[5]],
        complex_poynting_w_m2,
    })
}

fn phase(k: [f64; 2], position: [f64; 2]) -> Result<ComplexDd, AngularSpectrumError> {
    let argument = Dd::from_f64(k[0]) * Dd::from_f64(position[0])
        + Dd::from_f64(k[1]) * Dd::from_f64(position[1]);
    if !argument.hi.is_finite() || !argument.lo.is_finite() {
        return Err(bad("probe spatial phase is unrepresentable"));
    }
    let (sine, cosine) = det::sin_cos_dd(argument);
    Ok([cosine, -sine])
}

fn multiply(a: ComplexDd, b: ComplexDd) -> ComplexDd {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn rounded(value: ComplexDd) -> Result<C64, AngularSpectrumError> {
    let result = C64::new(value[0].to_f64(), value[1].to_f64());
    if !finite(result) {
        return Err(bad("probe field or Poynting vector is unrepresentable"));
    }
    Ok(result)
}
