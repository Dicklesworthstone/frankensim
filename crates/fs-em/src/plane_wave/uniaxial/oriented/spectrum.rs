//! Coherent vector angular spectra at one source-supported positive frequency.
//!
//! Sparse authored modes retain their amplitudes and phases; there is no FFT,
//! paraxial approximation, interpolation, clipping, filtering or power rescale.
//! Every mode sees the SAME fixed global crystal axis. Its ordinary and
//! extraordinary fields propagate separately, retaining birefringent phase,
//! spatial walk-off, evanescence and interference between different directions.
//! This is a monochromatic, planar, homogeneous transparent half-space, not a
//! finite crystal, absorbing anisotropy or a causal time-domain reservoir.

use super::{Incidence, OrientedCrystalResponse, OrientedUniaxialInterface, ResolvedFrequency};
use crate::plane_wave::angular_spectrum::{
    AngularSpectrumError, AngularSpectrumSide, modes::propagation,
};
use fs_math::{c64::C64, dd::Dd, det};

type Result<T> = core::result::Result<T, AngularSpectrumError>;

/// One Fourier coefficient, expressed in a fixed Cartesian interface frame.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniaxialIncidentMode {
    /// Signed k_t/k0, k0=omega/c0. Its magnitude must be less than the original
    /// incident index for every nonzero coefficient. Keeping this coordinate
    /// avoids an angle round-trip at grazing or crystal critical incidence.
    pub tangential_phase_index: f64,
    /// Finite nonzero [x,y] direction of the incidence plane, normalized on
    /// preparation. The physical tangential wavevector is k_t times this axis.
    /// Its sign and the signed phase index together specify the direction.
    pub incidence_plane_axis: [f64; 2],
    /// Global Cartesian [Ex,Ey] complex peak coefficient [V/m] at the origin.
    /// It multiplies exp(+i omega t-i k.r), with no implicit 1/N factor.
    pub electric_v_m: [C64; 2],
}

/// Checked before allocating modes or evaluating any original material law.
#[derive(Clone, Copy, Debug)]
pub struct UniaxialSpectrumLimits {
    /// Bound both preparation solves and O(N) work per point. Each query uses
    /// at most two propagated plane-wave fields per retained input mode in a
    /// half-space, or four inside a finite crystal slab.
    pub maximum_modes: usize,
    /// Logical retained mode storage plus one field and DD accumulator.
    /// Caller input/source storage, stack temporaries and allocator overhead
    /// are outside this bound. Queries allocate nothing.
    pub maximum_workspace_bytes: usize,
}

/// Actual coherent fields at the requested physical point, in global axes.
/// Real fields are Re(E exp(+i omega t)) and Re(H exp(+i omega t)).
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniaxialSpectrumField {
    pub electric_v_m: [C64; 3],
    pub magnetic_a_m: [C64; 3],
    /// Cycle-mean E cross H*/2, including every intermode cross term.
    pub poynting_w_m2: [f64; 3],
}

struct Mode {
    axis: [f64; 2],
    response: OrientedCrystalResponse,
}

/// Prepared coherent incident spectrum tied to its original interface owner.
/// Material resolution and vector Maxwell boundary solves happen only once.
pub struct PreparedUniaxialAngularSpectrum<'s, 'm> {
    source: &'s OrientedUniaxialInterface<'m>,
    omega_rad_s: f64,
    modes: Vec<Mode>,
    workspace_bytes: usize,
}

impl<'s, 'm> PreparedUniaxialAngularSpectrum<'s, 'm> {
    pub fn prepare(
        source: &'s OrientedUniaxialInterface<'m>,
        omega_rad_s: f64,
        input: &[UniaxialIncidentMode],
        limits: UniaxialSpectrumLimits,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self> {
        poll(&mut cancelled)?;
        let workspace_bytes = admit(input.len(), limits)?;
        let frequency = source.resolve_frequency(omega_rad_s)?;
        let modes = prepare_modes(frequency, source.optic_axis(), input, &mut cancelled)?;
        Ok(Self {
            source,
            omega_rad_s: frequency.omega,
            modes,
            workspace_bytes,
        })
    }

    pub const fn source(&self) -> &'s OrientedUniaxialInterface<'m> {
        self.source
    }
    pub const fn angular_frequency_rad_s(&self) -> f64 {
        self.omega_rad_s
    }
    /// Exactly zero coefficients are omitted; no nonzero coefficient is cut.
    pub fn retained_mode_count(&self) -> usize {
        self.modes.len()
    }
    pub const fn workspace_bytes(&self) -> usize {
        self.workspace_bytes
    }

    /// Full spatial fields at [x,y,z] [m], with interface z=0. Transmitted
    /// queries require z>=0; all incident-side queries require z<=0.
    /// TotalIncidentSide retains incident/reflected standing-wave terms.
    /// All modes and both crystal eigenmodes are summed before computing flux.
    pub fn field_at(
        &self,
        side: AngularSpectrumSide,
        position_m: [f64; 3],
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<UniaxialSpectrumField> {
        field_at(&self.modes, side, position_m, &mut cancelled)
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
fn finite(z: C64) -> bool {
    z.re.is_finite() && z.im.is_finite()
}
fn admit(count: usize, limits: UniaxialSpectrumLimits) -> Result<usize> {
    if count > limits.maximum_modes {
        return Err(bad("uniaxial spectrum exceeds its mode/work budget"));
    }
    count
        .checked_mul(core::mem::size_of::<Mode>())
        .and_then(|bytes| bytes.checked_add(core::mem::size_of::<FieldSum>()))
        .and_then(|bytes| bytes.checked_add(core::mem::size_of::<UniaxialSpectrumField>()))
        .filter(|&bytes| bytes <= limits.maximum_workspace_bytes)
        .ok_or_else(|| bad("uniaxial spectrum exceeds its workspace budget"))
}

fn prepare_modes(
    frequency: ResolvedFrequency,
    crystal_axis: [f64; 3],
    input: &[UniaxialIncidentMode],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<Vec<Mode>> {
    let mut modes = Vec::new();
    modes
        .try_reserve_exact(input.len())
        .map_err(|_| bad("uniaxial spectrum mode allocation failed"))?;
    for authored in input {
        poll(cancelled)?;
        if !authored.tangential_phase_index.is_finite()
            || authored.electric_v_m.iter().any(|&z| !finite(z))
        {
            return Err(bad("uniaxial spectrum modes must be finite"));
        }
        let p = super::super::kernel::unit_axis(authored.incidence_plane_axis).map_err(bad)?;
        if authored.electric_v_m == [C64::ZERO; 2] {
            continue;
        }
        // Rotation is passive: the material axis remains fixed in global
        // space while each incidence plane gets the original x-z solve.
        let [x, y, z] = crystal_axis;
        let axis = [x * p[0] + y * p[1], y * p[0] - x * p[1], z];
        let [ex, ey] = authored.electric_v_m;
        let electric = [
            ex.scale(p[0]) + ey.scale(p[1]),
            ey.scale(p[0]) - ex.scale(p[1]),
        ];
        let mut r = frequency.solve(
            axis,
            Incidence::TangentialIndex(authored.tangential_phase_index),
            electric,
        )?;
        r.incident_electric_v_m = rotate(r.incident_electric_v_m, p);
        r.incident_magnetic_a_m = rotate(r.incident_magnetic_a_m, p);
        r.reflected_electric_v_m = rotate(r.reflected_electric_v_m, p);
        r.reflected_magnetic_a_m = rotate(r.reflected_magnetic_a_m, p);
        for mode in &mut r.transmitted {
            mode.electric_v_m = rotate(mode.electric_v_m, p);
            mode.magnetic_a_m = rotate(mode.magnetic_a_m, p);
            mode.phase_direction = mode.phase_direction.map(|v| rotate_real(v, p));
            mode.energy_direction = mode.energy_direction.map(|v| rotate_real(v, p));
        }
        modes.push(Mode {
            axis: p,
            response: r,
        });
    }
    Ok(modes)
}

pub(super) fn rotate(v: [C64; 3], p: [f64; 2]) -> [C64; 3] {
    [
        v[0].scale(p[0]) - v[1].scale(p[1]),
        v[0].scale(p[1]) + v[1].scale(p[0]),
        v[2],
    ]
}
fn rotate_real(v: [f64; 3], p: [f64; 2]) -> [f64; 3] {
    [v[0] * p[0] - v[1] * p[1], v[0] * p[1] + v[1] * p[0], v[2]]
}

pub(super) struct FieldSum {
    electric: [[Dd; 2]; 3],
    magnetic: [[Dd; 2]; 3],
}
impl FieldSum {
    pub(super) const fn new() -> Self {
        Self {
            electric: [[Dd::ZERO; 2]; 3],
            magnetic: [[Dd::ZERO; 2]; 3],
        }
    }
    pub(super) fn add(&mut self, electric: [C64; 3], magnetic: [C64; 3], phase: C64) {
        let [pr, pi] = [phase.re, phase.im].map(Dd::from_f64);
        for (sum, vector) in [
            (&mut self.electric, electric),
            (&mut self.magnetic, magnetic),
        ] {
            for j in 0..3 {
                let [vr, vi] = [vector[j].re, vector[j].im].map(Dd::from_f64);
                sum[j][0] = sum[j][0] + vr * pr - vi * pi;
                sum[j][1] = sum[j][1] + vr * pi + vi * pr;
            }
        }
    }
    pub(super) fn finish(self) -> Result<UniaxialSpectrumField> {
        let electric = self
            .electric
            .map(|v| C64::new(v[0].to_f64(), v[1].to_f64()));
        let magnetic = self
            .magnetic
            .map(|v| C64::new(v[0].to_f64(), v[1].to_f64()));
        let pair = |j: usize, k: usize| {
            self.electric[j][0] * self.magnetic[k][0] + self.electric[j][1] * self.magnetic[k][1]
        };
        let flux = [(1, 2), (2, 0), (0, 1)]
            .map(|(j, k)| ((pair(j, k) - pair(k, j)) * Dd::from_f64(0.5)).to_f64());
        if electric.iter().chain(&magnetic).any(|&z| !finite(z))
            || flux.iter().any(|v| !v.is_finite())
        {
            return Err(bad("coherent uniaxial field is unrepresentable"));
        }
        Ok(UniaxialSpectrumField {
            electric_v_m: electric,
            magnetic_a_m: magnetic,
            poynting_w_m2: flux,
        })
    }
}

pub(super) fn transverse_phase(k: f64, axis: [f64; 2], point: [f64; 3]) -> Result<C64> {
    // Preserve the authored binary64 spatial products until trigonometric
    // reduction, just as the shared longitudinal propagation owner does.
    let phase = Dd::from_f64(k)
        * (Dd::from_f64(axis[0]) * Dd::from_f64(point[0])
            + Dd::from_f64(axis[1]) * Dd::from_f64(point[1]));
    if !phase.hi.is_finite() || !phase.lo.is_finite() {
        return Err(bad("uniaxial transverse phase is unrepresentable"));
    }
    let (sine, cosine) = det::sin_cos_dd(phase);
    Ok(C64::new(cosine.to_f64(), -sine.to_f64()))
}

fn field_at(
    modes: &[Mode],
    side: AngularSpectrumSide,
    point: [f64; 3],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<UniaxialSpectrumField> {
    poll(cancelled)?;
    if point.iter().any(|v| !v.is_finite())
        || if side == AngularSpectrumSide::Transmitted {
            point[2] < 0.0
        } else {
            point[2] > 0.0
        }
    {
        return Err(bad(
            "uniaxial field point must lie in its selected half-space",
        ));
    }
    let mut sum = FieldSum::new();
    for mode in modes {
        poll(cancelled)?;
        let r = &mode.response;
        let transverse = transverse_phase(r.tangential_wavenumber_rad_m, mode.axis, point)?;
        let mut add = |e, h, kz| -> Result<()> {
            let longitudinal = propagation(kz, point[2])
                .ok_or_else(|| bad("uniaxial propagation is unrepresentable"))?;
            sum.add(e, h, transverse * longitudinal);
            Ok(())
        };
        if matches!(
            side,
            AngularSpectrumSide::Incident | AngularSpectrumSide::TotalIncidentSide
        ) {
            add(
                r.incident_electric_v_m,
                r.incident_magnetic_a_m,
                C64::from_re(r.incident_normal_wavenumber_rad_m),
            )?;
        }
        if matches!(
            side,
            AngularSpectrumSide::Reflected | AngularSpectrumSide::TotalIncidentSide
        ) {
            add(
                r.reflected_electric_v_m,
                r.reflected_magnetic_a_m,
                C64::from_re(-r.incident_normal_wavenumber_rad_m),
            )?;
        }
        if side == AngularSpectrumSide::Transmitted {
            for wave in &r.transmitted {
                add(
                    wave.electric_v_m,
                    wave.magnetic_a_m,
                    wave.normal_wavenumber_rad_m,
                )?;
            }
        }
    }
    sum.finish()
}

#[cfg(test)]
mod tests;
