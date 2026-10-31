//! Harmonic vector fields of the existing prescribed-flow Rayleigh surface.
//!
//! With exp(-i omega t), a surface volume-velocity element Q gives
//! p = -i rho omega Q exp(ikR)/(2 pi R),
//! grad(p) = p (ik - 1/R) r_hat, and u = grad(p)/(i rho omega).
//! The velocity therefore retains the reactive 1/R² term. Coherent signed
//! surface contributions are summed before applying the existing `FirstOrder`
//! microphone; a pressure null need not be a velocity or directional null.
//!
//! This is the continuous-frequency image of `BaffledPressure`'s same infinite
//! rigid, flat baffle and caller-supplied area/velocity quadrature. It adds no
//! time filter or frequency warping. Spatial quadrature must still be refined
//! at the requested frequency and receiver clearance: no continuum error bound,
//! finite-baffle scattering, reaction load or molecular absorption is claimed.
//! The source convention is the conjugate of the Rayleigh formula in
//! <https://euphonics.org/4-3-2-the-rayleigh-integral-and-the-baffled-piston/>.

use super::{MAX_RAYLEIGH_INPUTS, RayleighMedium, SurfaceSample};
use fs_bem::helmholtz::Medium;
use fs_bem::near_field::FirstOrder;
use fs_math::{c64::C64, det};
use std::f64::consts::TAU;

/// Same receiver-array bound as the existing exterior BEM observer.
pub const MAX_RECEIVERS: usize = fs_bem::near_field::MAX_RECEIVERS;
/// Hard cap on (surface sample, input coordinate, receiver) contributions.
pub const MAX_KERNEL_TERMS: usize = 1 << 24;

/// Cold preparation budget and separation from the singular baffle plane.
#[derive(Clone, Copy, Debug)]
pub struct Options {
    /// Positive receiver height [m]. The default matches `BaffledPressure`.
    /// A smaller clearance requires correspondingly finer surface quadrature.
    pub minimum_clearance_m: f64,
    /// Complete quadrature work, no partial matrix on budget exhaustion.
    pub maximum_kernel_terms: usize,
}

impl Default for Options {
    fn default() -> Self {
        Self {
            minimum_clearance_m: 0.05,
            maximum_kernel_terms: 1 << 21,
        }
    }
}

/// Co-located SI fields with the common exp(-i omega t) convention.
#[derive(Clone, Copy, Debug, Default, PartialEq)]
pub struct Field {
    /// Acoustic pressure [Pa].
    pub pressure_pa: C64,
    /// Spatial derivative with respect to the receiver coordinates [Pa/m].
    pub pressure_gradient_pa_m: [C64; 3],
    /// Fluid particle velocity [m/s], including its reactive near field.
    pub particle_velocity_m_s: [C64; 3],
}

#[derive(Clone, Copy, Default)]
struct Column {
    pressure: C64,
    velocity: [C64; 3],
}

impl Column {
    fn field(self, omega_rho: f64) -> Field {
        Field {
            pressure_pa: self.pressure,
            pressure_gradient_pa_m: self.velocity.map(|u| u * C64::new(0.0, omega_rho)),
            particle_velocity_m_s: self.velocity,
        }
    }
}

#[derive(Clone, Copy, Default)]
struct Sum {
    value: C64,
    correction: C64,
}

impl Sum {
    fn add(&mut self, value: C64) {
        let adjusted = value - self.correction;
        let next = self.value + adjusted;
        self.correction = (next - self.value) - adjusted;
        self.value = next;
    }
}

fn finite(value: C64) -> bool {
    value.re.is_finite() && value.im.is_finite()
}
fn finite_field(value: &Field) -> bool {
    finite(value.pressure_pa)
        && value.pressure_gradient_pa_m.iter().copied().all(finite)
        && value.particle_velocity_m_s.iter().copied().all(finite)
}

/// One frequency's receiver-major transfer matrix, reusable across source phasors.
/// Preparation is O(samples * inputs * receivers); evaluation is O(inputs *
/// receivers), allocation-free, with compensated complex sums. Unit input
/// coordinates use the original signed `SurfaceSample::mode_shape`; a 1/area
/// shape therefore still means a volume-flow input, not a velocity input.
pub struct Prepared {
    columns: Vec<Column>,
    inputs: usize,
    receivers: usize,
    frequency_hz: f64,
    omega_rho: f64,
    medium: RayleighMedium,
}

impl Prepared {
    /// Prepare exact Green kernels at the supplied positive-area quadrature
    /// points, without changing their source basis, normal orientation or phase.
    /// Cancellation is polled before each receiver and every 64 surface samples.
    /// There is no sample-rate/Nyquist restriction on a continuous phasor solve.
    ///
    /// # Errors
    /// Invalid complete flat-surface basis, medium, positive frequency, clearance,
    /// receiver/input count, total work budget, cancellation or nonfinite fields.
    pub fn new(
        surface: &[SurfaceSample],
        inputs: usize,
        positions_m: &[[f64; 3]],
        frequency_hz: f64,
        medium: RayleighMedium,
        options: Options,
        mut cancelled: impl FnMut() -> bool,
    ) -> Result<Self, String> {
        let terms = surface
            .len()
            .checked_mul(inputs)
            .and_then(|n| n.checked_mul(positions_m.len()));
        if !(1..=MAX_RAYLEIGH_INPUTS).contains(&inputs)
            || !(1..=MAX_RECEIVERS).contains(&positions_m.len())
            || surface.is_empty()
            || surface.len() > 120_000
            || !(1..=MAX_KERNEL_TERMS).contains(&options.maximum_kernel_terms)
            || terms.is_none_or(|n| n > options.maximum_kernel_terms)
        {
            return Err(
                "harmonic Rayleigh input, receiver or complete kernel-work budget exceeded".into(),
            );
        }
        let omega = TAU * frequency_hz;
        let omega_rho = omega * medium.density;
        let k = omega / medium.sound_speed;
        if !frequency_hz.is_finite()
            || frequency_hz <= 0.0
            || !medium.density.is_finite()
            || medium.density <= 0.0
            || !medium.sound_speed.is_finite()
            || medium.sound_speed <= 0.0
            || !omega_rho.is_finite()
            || omega_rho <= 0.0
            || !k.is_finite()
            || k <= 0.0
            || !options.minimum_clearance_m.is_finite()
            || options.minimum_clearance_m <= 0.0
            || positions_m
                .iter()
                .any(|p| p.iter().any(|v| !v.is_finite()) || p[2] < options.minimum_clearance_m)
        {
            return Err("harmonic Rayleigh needs a finite positive medium/frequency and admitted baffle clearance".into());
        }
        for (index, sample) in surface.iter().enumerate() {
            if index % 64 == 0 && cancelled() {
                return Err("harmonic Rayleigh preparation cancelled".into());
            }
            if !sample.area_m2.is_finite()
                || sample.area_m2 <= 0.0
                || sample.position_m.iter().any(|x| !x.is_finite())
                || sample.position_m[2] != 0.0
                || sample.mode_shape.len() != inputs
                || sample.mode_shape.iter().any(|x| !x.is_finite())
            {
                return Err(
                    "Rayleigh sample must be a finite flat-panel area with all modal shapes".into(),
                );
            }
        }
        let mut columns = Vec::with_capacity(inputs * positions_m.len());
        let mut sums = vec![[Sum::default(); 4]; inputs];
        for position in positions_m {
            if cancelled() {
                return Err("harmonic Rayleigh preparation cancelled".into());
            }
            sums.fill([Sum::default(); 4]);
            for (index, sample) in surface.iter().enumerate() {
                if index % 64 == 0 && cancelled() {
                    return Err("harmonic Rayleigh preparation cancelled".into());
                }
                let delta: [f64; 3] = std::array::from_fn(|i| position[i] - sample.position_m[i]);
                let radius = det::sqrt(delta.iter().map(|x| x * x).sum());
                let phase = k * radius;
                if !radius.is_finite() || radius < options.minimum_clearance_m || !phase.is_finite()
                {
                    return Err(
                        "harmonic Rayleigh separation or propagation phase is unrepresentable"
                            .into(),
                    );
                }
                let green = C64::new(det::cos(phase), det::sin(phase))
                    .scale(sample.area_m2 / (TAU * radius));
                let pressure = green * C64::new(0.0, -omega_rho);
                // Evaluate u directly: division by omega would lose the regular
                // low-frequency volume-flow limit and its 1/R² reactive term.
                let radial_velocity = green * C64::new(1.0 / radius, -k);
                let velocity = delta.map(|x| radial_velocity.scale(x / radius));
                for (sum, &shape) in sums.iter_mut().zip(&sample.mode_shape) {
                    sum[0].add(pressure.scale(shape));
                    for c in 0..3 {
                        sum[c + 1].add(velocity[c].scale(shape));
                    }
                }
            }
            for sum in &sums {
                let column = Column {
                    pressure: sum[0].value,
                    velocity: std::array::from_fn(|c| sum[c + 1].value),
                };
                if !finite_field(&column.field(omega_rho)) {
                    return Err("harmonic Rayleigh transfer overflow".into());
                }
                columns.push(column);
            }
        }
        Ok(Self {
            columns,
            inputs,
            receivers: positions_m.len(),
            frequency_hz,
            omega_rho,
            medium,
        })
    }

    /// Number of original source coordinates.
    pub fn input_count(&self) -> usize {
        self.inputs
    }
    /// Number of receivers, in the caller's original order.
    pub fn receiver_count(&self) -> usize {
        self.receivers
    }
    /// Continuous frequency used by every field and transfer column [Hz].
    pub fn frequency_hz(&self) -> f64 {
        self.frequency_hz
    }
    /// Exact medium used at preparation and by directional observation.
    pub fn medium(&self) -> RayleighMedium {
        self.medium
    }

    /// Field per unit input coordinate, suitable for frequency-transfer sampling.
    pub fn transfer(&self, receiver: usize, input: usize) -> Option<Field> {
        if receiver >= self.receivers || input >= self.inputs {
            return None;
        }
        Some(self.columns[receiver * self.inputs + input].field(self.omega_rho))
    }

    /// Apply all coherent source phasors to all receivers without allocation.
    /// Output is unchanged on any refusal, including an overflowing source sum.
    ///
    /// # Errors
    /// Incomplete/nonfinite source phasors, output count mismatch or overflow.
    pub fn evaluate_into(
        &self,
        source_velocity: &[C64],
        output: &mut [Field],
    ) -> Result<(), String> {
        if source_velocity.len() != self.inputs
            || output.len() != self.receivers
            || source_velocity.iter().copied().any(|v| !finite(v))
        {
            return Err(
                "harmonic Rayleigh requires complete finite source phasors and receiver outputs"
                    .into(),
            );
        }
        let mut pending = [Field::default(); MAX_RECEIVERS];
        for (receiver, row) in self.columns.chunks_exact(self.inputs).enumerate() {
            let mut sums = [Sum::default(); 4];
            for (column, &source) in row.iter().zip(source_velocity) {
                sums[0].add(column.pressure * source);
                for c in 0..3 {
                    sums[c + 1].add(column.velocity[c] * source);
                }
            }
            pending[receiver] = Column {
                pressure: sums[0].value,
                velocity: std::array::from_fn(|c| sums[c + 1].value),
            }
            .field(self.omega_rho);
            if !finite_field(&pending[receiver]) {
                return Err("harmonic Rayleigh field overflow".into());
            }
        }
        output.copy_from_slice(&pending[..self.receivers]);
        Ok(())
    }

    /// Pa-equivalent ideal microphone outputs using the original BEM `FirstOrder`
    /// pattern and this same medium. No capsule response or backreaction is added.
    /// Output is unchanged on refusal; source and receiver phases stay coherent.
    ///
    /// # Errors
    /// Mismatched pattern/output counts or any field/projection refusal.
    pub fn observe_into(
        &self,
        source_velocity: &[C64],
        patterns: &[FirstOrder],
        output: &mut [C64],
    ) -> Result<(), String> {
        if patterns.len() != self.receivers || output.len() != self.receivers {
            return Err("harmonic Rayleigh needs one pattern and output per receiver".into());
        }
        let mut fields = [Field::default(); MAX_RECEIVERS];
        self.evaluate_into(source_velocity, &mut fields[..self.receivers])?;
        let medium = Medium {
            density: self.medium.density,
            sound_speed: self.medium.sound_speed,
        };
        let mut pending = [C64::ZERO; MAX_RECEIVERS];
        for (i, (&pattern, field)) in patterns.iter().zip(&fields).enumerate() {
            pending[i] = pattern
                .observe(field.pressure_pa, field.particle_velocity_m_s, medium)
                .map_err(|error| error.to_string())?;
        }
        output.copy_from_slice(&pending[..self.receivers]);
        Ok(())
    }
}

#[cfg(test)]
#[path = "baffled_pressure_harmonic_tests.rs"]
mod tests;
