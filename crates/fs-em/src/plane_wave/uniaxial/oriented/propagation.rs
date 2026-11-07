//! Spatial phasors from the original solved Maxwell boundary fields.
use super::{OrientedCrystalResponse, Result, Vector};
use fs_math::{c64::C64, dd::Dd, decay, det};

/// Components selected at one physical point. The interface is z=0, with the
/// original incident medium at z<0 and the transparent crystal at z>0.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum OrientedCrystalFieldSide {
    Incident,
    Reflected,
    /// Coherent incident plus reflected fields, including standing waves.
    TotalIncidentSide,
    /// The transmitted ordinary eigenmode by itself.
    Ordinary,
    /// The transmitted extraordinary eigenmode by itself.
    Extraordinary,
    /// Both transmitted eigenmodes, summed before power is evaluated.
    Transmitted,
}

/// Full Cartesian peak phasors at the requested physical point.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedCrystalField {
    pub electric_v_m: Vector,
    pub magnetic_a_m: Vector,
    /// Re(E cross conj(H))/2 [W/m²], with the physical Cartesian sign.
    /// Coherent products are accumulated before rounding the result. This
    /// includes cross-mode interference, not a sum of individual mode powers.
    pub poynting_w_m2: [f64; 3],
}

type Complex = [Dd; 2];
type Fields = [Complex; 6];
const ZERO: Complex = [Dd::ZERO; 2];

#[derive(Clone, Copy)]
struct Wave<'a> {
    normal: C64,
    electric: &'a Vector,
    magnetic: &'a Vector,
}

/// Geometry of the first active wave, retained as a shared phase reference.
#[derive(Clone, Copy)]
struct Carrier {
    tangential: f64,
    axis: [f64; 2],
    normal: f64,
}

/// One physical-point accumulator shared by single waves and sparse spectra.
/// The fields are in the caller's global Cartesian frame. No source law,
/// field rotation, allocation, budget or cancellation policy lives here.
pub(crate) struct FieldSum {
    fields: Fields,
    side: OrientedCrystalFieldSide,
    point: [f64; 3],
    reference: Option<Carrier>,
}

impl FieldSum {
    /// Accumulate fields whose spatial/face-anchored phasors are already
    /// owned by the caller. The zero carrier adds no phase or admission rule.
    pub(crate) const fn from_phasors() -> Self {
        Self {
            fields: [ZERO; 6],
            side: OrientedCrystalFieldSide::Transmitted,
            point: [0.0; 3],
            reference: Some(Carrier {
                tangential: 0.0,
                axis: [0.0; 2],
                normal: 0.0,
            }),
        }
    }

    /// Infallible accumulation for existing callers with admitted phasors.
    /// Nonfinite products are rejected by finish before a field is published.
    pub(crate) fn add_phasor(&mut self, electric: Vector, magnetic: Vector, phase: C64) {
        let [pr, pi] = [phase.re, phase.im].map(Dd::from_f64);
        for (j, sum) in self.fields.iter_mut().enumerate() {
            let value = if j < 3 { electric[j] } else { magnetic[j - 3] };
            let [vr, vi] = [value.re, value.im].map(Dd::from_f64);
            // Keep the existing slab accumulator's arithmetic association.
            sum[0] = sum[0] + vr * pr - vi * pi;
            sum[1] = sum[1] + vr * pi + vi * pr;
        }
    }

    pub(crate) fn new(side: OrientedCrystalFieldSide, point: [f64; 3]) -> Result<Self> {
        use OrientedCrystalFieldSide::{Extraordinary, Ordinary, Transmitted};
        if point.iter().any(|v| !v.is_finite()) {
            return Err("crystal observation coordinates must be finite metres");
        }
        let crystal_side = matches!(side, Ordinary | Extraordinary | Transmitted);
        if crystal_side && point[2] < 0.0 || !crystal_side && point[2] > 0.0 {
            return Err("crystal field observation is outside the selected half-space");
        }
        Ok(Self { fields: [ZERO; 6], side, point, reference: None })
    }

    /// Add one already solved incident direction. Its response fields have
    /// already been rotated to global axes; axis is its normalized [x,y]
    /// incidence-plane direction. The scalar tangential index keeps its sign.
    pub(crate) fn add_response(
        &mut self,
        boundary: &OrientedCrystalResponse,
        axis: [f64; 2],
    ) -> Result<()> {
        let Some((values, normal)) = pair_field(boundary, self.side, self.point[2])? else {
            return Ok(());
        };
        let carrier = Carrier {
            tangential: boundary.tangential_wavenumber_rad_m,
            axis,
            normal,
        };
        if let Some(reference) = self.reference {
            let rotation = relative_rotation(carrier, reference, self.point)?;
            for (sum, value) in self.fields.iter_mut().zip(values) {
                let value = multiply(value, rotation);
                *sum = [sum[0] + value[0], sum[1] + value[1]];
            }
        } else {
            self.reference = Some(carrier);
            self.fields = values;
        }
        Ok(())
    }

    pub(crate) fn finish(self) -> Result<OrientedCrystalField> {
        if self.fields.iter().all(|&v| v == ZERO) {
            // Empty, cancelled or fully underflowed physical fields need no
            // absolute phase, even when that phase would overflow.
            return Ok(OrientedCrystalField {
                electric_v_m: [C64::ZERO; 3],
                magnetic_a_m: [C64::ZERO; 3],
                poynting_w_m2: [0.0; 3],
            });
        }
        let reference = self.reference.ok_or("crystal field has no phase reference")?;
        // All directions and eigenmodes interfere in the SAME DD sum. A
        // common unit phase cannot change the physical energy transport.
        let poynting_w_m2 = core::array::from_fn(|j| {
            let a = (j + 1) % 3;
            let b = (j + 2) % 3;
            ((real_product(self.fields[a], self.fields[b + 3])
                - real_product(self.fields[b], self.fields[a + 3]))
                * Dd::from_f64(0.5))
            .to_f64()
        });
        let rotation = relative_rotation(
            reference,
            Carrier { tangential: 0.0, axis: [0.0; 2], normal: 0.0 },
            self.point,
        )?;
        let rounded = self.fields.map(|value| {
            let value = multiply(value, rotation);
            C64::new(value[0].to_f64(), value[1].to_f64())
        });
        if rounded.iter().any(|&z| !super::finite(z))
            || poynting_w_m2.iter().any(|v| !v.is_finite())
        {
            return Err("propagated crystal field or Poynting vector is unrepresentable");
        }
        Ok(OrientedCrystalField {
            electric_v_m: [rounded[0], rounded[1], rounded[2]],
            magnetic_a_m: [rounded[3], rounded[4], rounded[5]],
            poynting_w_m2,
        })
    }
}

pub(crate) fn field(
    boundary: &OrientedCrystalResponse,
    side: OrientedCrystalFieldSide,
    point: [f64; 3],
) -> Result<OrientedCrystalField> {
    let mut sum = FieldSum::new(side, point)?;
    sum.add_response(boundary, [1.0, 0.0])?;
    sum.finish()
}

fn pair_field(
    boundary: &OrientedCrystalResponse,
    side: OrientedCrystalFieldSide,
    z: f64,
) -> Result<Option<(Fields, f64)>> {
    use OrientedCrystalFieldSide::{
        Extraordinary, Incident, Ordinary, Reflected, TotalIncidentSide, Transmitted,
    };
    let incident = Wave {
        normal: C64::from_re(boundary.incident_normal_wavenumber_rad_m),
        electric: &boundary.incident_electric_v_m,
        magnetic: &boundary.incident_magnetic_a_m,
    };
    let reflected = Wave {
        normal: -incident.normal,
        electric: &boundary.reflected_electric_v_m,
        magnetic: &boundary.reflected_magnetic_a_m,
    };
    let transmitted = |index: usize| Wave {
        normal: boundary.transmitted[index].normal_wavenumber_rad_m,
        electric: &boundary.transmitted[index].electric_v_m,
        magnetic: &boundary.transmitted[index].magnetic_a_m,
    };
    let waves = match side {
        Incident => [Some(incident), None],
        Reflected => [Some(reflected), None],
        TotalIncidentSide => [Some(incident), Some(reflected)],
        Ordinary => [Some(transmitted(0)), None],
        Extraordinary => [Some(transmitted(1)), None],
        Transmitted => [Some(transmitted(0)), Some(transmitted(1))],
    };
    let mut sum = [ZERO; 6];
    let mut reference = None;
    for wave in waves.into_iter().flatten() {
        let values = attenuated(wave, z)?;
        if values.iter().all(|&v| v == ZERO) {
            continue;
        }
        if let Some(normal) = reference {
            // Remove the shared carrier before the two waves interfere.
            // A large absolute phase must not erase weak retardance or the
            // relative phase of an incident/reflected standing wave.
            let relative = (Dd::from_f64(wave.normal.re) - Dd::from_f64(normal))
                * Dd::from_f64(z);
            let rotation = phasor(relative)?;
            for (total, value) in sum.iter_mut().zip(values) {
                let value = multiply(value, rotation);
                *total = [total[0] + value[0], total[1] + value[1]];
            }
        } else {
            reference = Some(wave.normal.re);
            sum = values;
        }
    }
    Ok(reference.map(|normal| (sum, normal)))
}

fn relative_rotation(carrier: Carrier, reference: Carrier, point: [f64; 3]) -> Result<Complex> {
    let kt = Dd::from_f64(carrier.tangential);
    let reference_kt = Dd::from_f64(reference.tangential);
    let transverse = |dimension: usize| {
        (kt * Dd::from_f64(carrier.axis[dimension])
            - reference_kt * Dd::from_f64(reference.axis[dimension]))
            * Dd::from_f64(point[dimension])
    };
    // Reduce each spatial product separately. A huge x carrier must not
    // round away a physically significant y or z phase; close directions
    // retain their difference before multiplication by a large coordinate.
    let x = phasor(transverse(0))?;
    let y = phasor(transverse(1))?;
    let z = phasor(
        (Dd::from_f64(carrier.normal) - Dd::from_f64(reference.normal))
            * Dd::from_f64(point[2]),
    )?;
    Ok(multiply(multiply(x, y), z))
}

fn attenuated(wave: Wave<'_>, z: f64) -> Result<Fields> {
    let exposure = -wave.normal.im * z;
    if !(exposure >= 0.0) {
        return Err("crystal evanescent field must decay into its half-space");
    }
    let remaining = decay::remaining_fraction(exposure);
    // expm1 retains weak attenuation in cancellation-sensitive fields; the
    // shared decay owner retains a late tail whose fraction underflows first.
    let factor = if exposure < 0.5 {
        Dd::ONE + Dd::from_f64(det::expm1(-exposure))
    } else {
        Dd::from_f64(remaining)
    };
    Ok(core::array::from_fn(|j| {
        let value = if j < 3 {
            wave.electric[j]
        } else {
            wave.magnetic[j - 3]
        };
        [value.re, value.im].map(|v| {
            if remaining >= f64::MIN_POSITIVE {
                Dd::from_f64(v) * factor
            } else {
                Dd::from_f64(decay::attenuate(v, exposure))
            }
        })
    }))
}

fn phasor(phase: Dd) -> Result<Complex> {
    if !phase.hi.is_finite() || !phase.lo.is_finite() {
        return Err("crystal propagation phase is unrepresentable");
    }
    let (sine, cosine) = det::sin_cos_dd(phase);
    Ok([cosine, -sine])
}

fn multiply(a: Complex, b: Complex) -> Complex {
    [a[0] * b[0] - a[1] * b[1], a[0] * b[1] + a[1] * b[0]]
}

fn real_product(a: Complex, b: Complex) -> Dd {
    a[0] * b[0] + a[1] * b[1]
}

#[cfg(test)]
#[path = "propagation_tests.rs"]
mod tests;
