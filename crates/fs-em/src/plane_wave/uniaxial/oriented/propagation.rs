//! Spatial phasors from the original solved Maxwell boundary fields.
use super::{OrientedCrystalResponse, Result, Vector};
use fs_math::{c64::C64, dd::Dd, decay, det};

#[path = "propagation/frequency.rs"]
mod frequency;
pub(in crate::plane_wave::uniaxial::oriented) use frequency::{
    FrequencyWave, field_frequency_derivatives,
};

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

/// An already solved plane wave, with fields anchored at one physical z face.
/// The normal wavevector keeps its global sign: backward evanescence has
/// Im(kz)>0 and decays toward z<anchor_z_m. Transverse origins remain x=y=0.
#[derive(Clone, Copy)]
pub(crate) struct WaveGeometry {
    pub(crate) tangential_wavenumber_rad_m: f64,
    pub(crate) incidence_plane_axis: [f64; 2],
    pub(crate) normal_wavenumber_rad_m: C64,
    pub(crate) anchor_z_m: f64,
}

/// Geometry of the first active wave, retained as a shared phase reference.
#[derive(Clone, Copy)]
struct Carrier {
    tangential: f64,
    axis: [f64; 2],
    normal: f64,
    anchor_z: f64,
}

/// One physical-point accumulator shared by single waves and sparse spectra.
/// The fields are in the caller's global Cartesian frame. No source law,
/// field rotation, allocation, budget or cancellation policy lives here.
pub(crate) struct FieldSum {
    fields: Fields,
    side: OrientedCrystalFieldSide,
    point: [f64; 2],
    depth: Dd,
    reference: Option<Carrier>,
}

impl FieldSum {
    /// Accumulate fields whose spatial/face-anchored phasors are already
    /// owned by the caller. The zero carrier adds no phase or admission rule.
    pub(crate) const fn from_phasors() -> Self {
        Self {
            fields: [ZERO; 6],
            side: OrientedCrystalFieldSide::Transmitted,
            point: [0.0; 2],
            depth: Dd::ZERO,
            reference: Some(Carrier {
                tangential: 0.0,
                axis: [0.0; 2],
                normal: 0.0,
                anchor_z: 0.0,
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

    /// Construct a physical-point sum after the caller's region admission.
    /// Generic anchored waves need no single-interface half-space selection.
    pub(crate) fn at_point(point: [f64; 3]) -> Result<Self> {
        Self::at_local_point([point[0], point[1]], Dd::from_f64(point[2]))
    }

    /// Retain a compensated axial coordinate in the caller's physical frame.
    /// Layer owners can supply exact global-coordinate subtraction internally,
    /// and their original authored thickness at an admitted exit face. Region
    /// admission and the choice of axial origin remain entirely with them.
    pub(crate) fn at_local_point(point: [f64; 2], depth: Dd) -> Result<Self> {
        if point.iter().any(|v| !v.is_finite())
            || !depth.hi.is_finite()
            || !depth.lo.is_finite()
            || !depth.to_f64().is_finite()
        {
            return Err("crystal observation coordinates must be finite metres");
        }
        Ok(Self {
            fields: [ZERO; 6],
            side: OrientedCrystalFieldSide::Transmitted,
            point,
            depth,
            reference: None,
        })
    }

    pub(crate) fn new(side: OrientedCrystalFieldSide, point: [f64; 3]) -> Result<Self> {
        use OrientedCrystalFieldSide::{Extraordinary, Ordinary, Transmitted};
        let mut sum = Self::at_point(point)?;
        let crystal_side = matches!(side, Ordinary | Extraordinary | Transmitted);
        if crystal_side && point[2] < 0.0 || !crystal_side && point[2] > 0.0 {
            return Err("crystal field observation is outside the selected half-space");
        }
        sum.side = side;
        Ok(sum)
    }

    /// Add full global E/H fields before rounding their propagation factor.
    /// Decay uses the authored face distance, including its subtraction
    /// residual; phase is relative to the first active wave across all faces.
    pub(crate) fn add_wave(
        &mut self,
        electric: &Vector,
        magnetic: &Vector,
        geometry: WaveGeometry,
    ) -> Result<()> {
        if !geometry.anchor_z_m.is_finite() {
            return Err("crystal wave anchor must be finite metres");
        }
        let values = attenuated(
            Wave {
                normal: geometry.normal_wavenumber_rad_m,
                electric,
                magnetic,
            },
            self.depth - Dd::from_f64(geometry.anchor_z_m),
        )?;
        self.add_fields(
            values,
            Carrier {
                tangential: geometry.tangential_wavenumber_rad_m,
                axis: geometry.incidence_plane_axis,
                normal: geometry.normal_wavenumber_rad_m.re,
                anchor_z: geometry.anchor_z_m,
            },
        )
    }

    /// Add one already solved incident direction. Its response fields have
    /// already been rotated to global axes; axis is its normalized [x,y]
    /// incidence-plane direction. The scalar tangential index keeps its sign.
    pub(crate) fn add_response(
        &mut self,
        boundary: &OrientedCrystalResponse,
        axis: [f64; 2],
    ) -> Result<()> {
        let Some((values, normal)) = pair_field(boundary, self.side, self.depth)? else {
            return Ok(());
        };
        let carrier = Carrier {
            tangential: boundary.tangential_wavenumber_rad_m,
            axis,
            normal,
            anchor_z: 0.0,
        };
        self.add_fields(values, carrier)
    }

    fn add_fields(&mut self, values: Fields, carrier: Carrier) -> Result<()> {
        if values.iter().all(|&v| v == ZERO) {
            return Ok(());
        }
        if let Some(reference) = self.reference {
            let rotation = relative_rotation(carrier, reference, self.point, self.depth)?;
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
        let reference = self
            .reference
            .ok_or("crystal field has no phase reference")?;
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
            Carrier {
                tangential: 0.0,
                axis: [0.0; 2],
                normal: 0.0,
                anchor_z: 0.0,
            },
            self.point,
            self.depth,
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
    z: Dd,
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
            let relative = (Dd::from_f64(wave.normal.re) - Dd::from_f64(normal)) * z;
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

fn relative_rotation(
    carrier: Carrier,
    reference: Carrier,
    point: [f64; 2],
    depth: Dd,
) -> Result<Complex> {
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
    let k = Dd::from_f64(carrier.normal);
    let reference_k = Dd::from_f64(reference.normal);
    let anchor = Dd::from_f64(carrier.anchor_z);
    // k*(z-a) - kr*(z-ar) = (k-kr)*(z-a) + kr*(ar-a).
    // Difference the carriers and face coordinates BEFORE their products.
    // Equal k at different faces keeps its physical anchor phase; equal
    // faces retain weak retardance without subtracting huge absolute phases.
    let z = phasor(
        (k - reference_k) * (depth - anchor)
            + reference_k * (Dd::from_f64(reference.anchor_z) - anchor),
    )?;
    Ok(multiply(multiply(x, y), z))
}

fn attenuated(wave: Wave<'_>, distance: Dd) -> Result<Fields> {
    let exposure = attenuation_exposure(wave.normal.im, distance)?;
    let values = core::array::from_fn(|j| {
        let value = if j < 3 {
            wave.electric[j]
        } else {
            wave.magnetic[j - 3]
        };
        [value.re, value.im].map(Dd::from_f64)
    });
    Ok(attenuated_fields(values, exposure))
}

fn attenuation_exposure(normal_imaginary: f64, distance: Dd) -> Result<f64> {
    let leading = -normal_imaginary * distance.hi;
    // Positive overflow means complete extinction, not an invalid DD product.
    // Otherwise retain the authored distance residual until exposure rounds.
    let exposure = if leading.is_infinite() {
        leading
    } else {
        (Dd::from_f64(-normal_imaginary) * distance).to_f64()
    };
    if !(exposure >= 0.0) {
        return Err("crystal evanescent field must decay into its half-space");
    }
    Ok(exposure)
}

// Shared by native carriers and weighted spatial frequency polynomials.
// Keep DD interference residuals through weak attenuation; late tails use the
// same original signed decay owner independently on the retained DD parts.
fn attenuated_fields(values: Fields, exposure: f64) -> Fields {
    let remaining = decay::remaining_fraction(exposure);
    let factor = if exposure < 0.5 {
        Dd::ONE + Dd::from_f64(det::expm1(-exposure))
    } else {
        Dd::from_f64(remaining)
    };
    values.map(|value| {
        value.map(|v| {
            if remaining >= f64::MIN_POSITIVE {
                v * factor
            } else {
                let leading = Dd::from_f64(decay::attenuate(v.hi, exposure));
                if v.lo == 0.0 {
                    leading
                } else {
                    leading + Dd::from_f64(decay::attenuate(v.lo, exposure))
                }
            }
        })
    })
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
