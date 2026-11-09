//! Original coherent coating/etalon frequency derivatives, without differencing.
//!
//! The native slab stores a bounded matrix M and scale s with physical T=M/s.
//! At the carrier we HOLD s fixed while differentiating T; the same applies to
//! each homogeneous state normalization. No derivative of an arbitrary max
//! normalization, growing exponential, or nearly-unit reflection product is
//! needed. Original M anchors every carrier value and the native power solve
//! remains authoritative. The entire functions cos(sqrt(w)), sinc(sqrt(w))
//! retain the regular finite-film limit at w=(k_z d)^2=0.

use super::{LayeredReflection, PlanarDielectricLayers, Slab, bad, divided, maximum};
use crate::plane_wave::{
    Polarization, ReflectionError, WaveMedium,
    dispersion::{
        FresnelAmplitudeDispersion, IncidenceFrequencyConstraint, normal_squared_derivatives,
        root_derivatives,
    },
    finite, positive,
};
use fs_math::c64::C64;

type Result<T> = core::result::Result<T, ReflectionError>;

/// Coherent frequency derivatives attached to the original material/geometry
/// stack. Reflection is referenced to its entrance and transmission to its
/// physical exit face, including the original film propagation phase.
#[derive(Clone, Copy)]
pub struct LayeredFrequencyDispersion<'s, 'm> {
    source: &'s PlanarDielectricLayers<'m>,
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
    carrier: LayeredReflection,
    reflected: FresnelAmplitudeDispersion,
    transmitted: Option<[FresnelAmplitudeDispersion; 2]>,
}

impl<'s, 'm> LayeredFrequencyDispersion<'s, 'm> {
    pub const fn source(self) -> &'s PlanarDielectricLayers<'m> {
        self.source
    }
    pub const fn incidence_cosine(self) -> f64 {
        self.cosine
    }
    pub const fn polarization(self) -> Polarization {
        self.polarization
    }
    pub const fn frequency_constraint(self) -> IncidenceFrequencyConstraint {
        self.constraint
    }
    pub const fn carrier(self) -> LayeredReflection {
        self.carrier
    }
    pub const fn reflected_tangential_e(self) -> FresnelAmplitudeDispersion {
        self.reflected
    }
    /// None after propagation underflow: a rounded zero carrier does not prove
    /// that its first or second frequency derivative is zero.
    pub fn transmitted_tangential_e(self) -> Option<FresnelAmplitudeDispersion> {
        self.transmitted.map(|value| value[0])
    }
    /// The same propagation-underflow refusal applies to both exit fields.
    pub fn transmitted_conjugate_h(self) -> Option<FresnelAmplitudeDispersion> {
        self.transmitted.map(|value| value[1])
    }
}

impl<'m> PlanarDielectricLayers<'m> {
    /// Analytic reflection/transmission jets and local delay/GDD of the exact
    /// original coherent stack at its prepared carrier. Both material dispersion
    /// and finite-film propagation/interference are retained. Work is O(layers)
    /// with the existing 64-layer/pole bounds and no allocation in this query.
    ///
    /// Internal critical films use entire-function derivatives; an exactly
    /// critical exit, nonsmooth/out-of-band source jet or non-lossless incident
    /// jet refuses. The native carrier response remains independently usable.
    /// The original attenuation-underflow limit retains its reflected
    /// half-space jet; transmitted jets become unavailable. The existing
    /// first-order `frequency_dispersion` query retains its cheaper tangent
    /// path and broader first-derivative representability domain.
    /// No incoherent phase, broadband pulse approximation, absorption derivative
    /// or measured-material validity beyond the original sources is inferred.
    pub fn second_order_frequency_dispersion(
        &self,
        incidence_cosine: f64,
        polarization: Polarization,
        constraint: IncidenceFrequencyConstraint,
    ) -> Result<LayeredFrequencyDispersion<'_, 'm>> {
        let carrier = self.response(incidence_cosine, polarization)?;
        let incident = self.incident_response.bulk_dispersion(self.omega)?;
        let ei = incident
            .scaled_permittivity_derivatives_f_m()
            .map(|value| divided(value, self.incident_epsilon_f_m));
        if ei.iter().any(|value| value.im != 0.0) {
            return Err(bad("incident frequency jet must remain lossless"));
        }
        let exit = self.exit.bulk_dispersion(self.omega)?;
        let exit_et = exit
            .scaled_permittivity_derivatives_f_m()
            .map(|value| divided(value, self.incident_epsilon_f_m));
        let mut state = admittance_jets(
            self.exit_wave,
            incidence_cosine,
            polarization,
            constraint,
            ei,
            exit_et,
        )?;
        let sum = state[0] + state[1];
        let mut transmitted = Some([state[0].scale(2.0) / sum, state[1].scale(2.0) / sum]);
        for layer in self.layers.iter().rev() {
            // Even an identity/hidden film retains its authored source domain.
            let bulk = layer.source.response.bulk_dispersion(self.omega)?;
            let et = bulk
                .scaled_permittivity_derivatives_f_m()
                .map(|value| divided(value, self.incident_epsilon_f_m));
            if layer.source.thickness_m == 0.0 {
                continue;
            }
            let slab = Slab::new(
                layer.wave,
                self.wave_number,
                incidence_cosine,
                polarization,
                layer.source.thickness_m,
            )?;
            if slab.scale == C64::ZERO {
                state = admittance_jets(
                    layer.wave,
                    incidence_cosine,
                    polarization,
                    constraint,
                    ei,
                    et,
                )?;
                transmitted = None;
                continue;
            }
            let matrix = slab_jets(
                slab,
                layer.wave,
                self.wave_number * layer.source.thickness_m,
                incidence_cosine,
                polarization,
                constraint,
                ei,
                et,
            )?;
            let mut head = [
                matrix[0][0] * state[0] + matrix[0][1] * state[1],
                matrix[1][0] * state[0] + matrix[1][1] * state[1],
            ];
            let normalization = positive(maximum(head[0].value).max(maximum(head[1].value)))?;
            head = head.map(|value| value.divide_real(normalization));
            if head.into_iter().any(|value| !value.is_finite()) {
                return Err(bad("coherent film frequency jet is unrepresentable"));
            }
            if let Some(tail) = transmitted {
                let gain = Jet::constant(divided(slab.scale, normalization))
                    * ((state[0] + state[1]) / (head[0] + head[1]));
                let next = tail.map(|value| gain * value);
                if gain.value == C64::ZERO
                    || tail
                        .into_iter()
                        .zip(next)
                        .any(|(old, new)| old.value != C64::ZERO && new.value == C64::ZERO)
                {
                    transmitted = None;
                } else {
                    if next.into_iter().any(|value| !value.is_finite()) {
                        return Err(bad("coherent film transmission jet is unrepresentable"));
                    }
                    transmitted = Some(next);
                }
            }
            state = head;
        }
        let reflection = reflection_jet(state);
        let reflected = FresnelAmplitudeDispersion::try_from_scaled(
            carrier.reflected_tangential_e,
            self.omega,
            reflection.derivatives(),
        )?;
        let transmitted = transmitted
            .filter(|_| {
                carrier.transmitted_tangential_e != C64::ZERO
                    && carrier.transmitted_conjugate_h != C64::ZERO
            })
            .map(|jets| -> Result<_> {
                Ok([
                    FresnelAmplitudeDispersion::try_from_scaled(
                        carrier.transmitted_tangential_e,
                        self.omega,
                        jets[0].derivatives(),
                    )?,
                    FresnelAmplitudeDispersion::try_from_scaled(
                        carrier.transmitted_conjugate_h,
                        self.omega,
                        jets[1].derivatives(),
                    )?,
                ])
            })
            .transpose()?;
        Ok(LayeredFrequencyDispersion {
            source: self,
            cosine: incidence_cosine,
            polarization,
            constraint,
            carrier,
            reflected,
            transmitted,
        })
    }
}

// Frequency-local scalar algebra. First/second are omega*H' and omega²*H'',
// ordinary derivatives, not factorial-scaled Taylor coefficients.
#[derive(Clone, Copy)]
struct Jet {
    value: C64,
    first: C64,
    second: C64,
}
impl Jet {
    const fn constant(value: C64) -> Self {
        Self {
            value,
            first: C64::ZERO,
            second: C64::ZERO,
        }
    }
    fn new(value: C64, derivatives: [C64; 2]) -> Self {
        Self {
            value,
            first: derivatives[0],
            second: derivatives[1],
        }
    }
    fn derivatives(self) -> [C64; 2] {
        [self.first, self.second]
    }
    fn scale(self, value: f64) -> Self {
        Self {
            value: self.value.scale(value),
            first: self.first.scale(value),
            second: self.second.scale(value),
        }
    }
    fn divide_real(self, value: f64) -> Self {
        Self {
            value: divided(self.value, value),
            first: divided(self.first, value),
            second: divided(self.second, value),
        }
    }
    fn root_at(self, root: C64) -> Self {
        Self::new(root, root_derivatives(root, self.derivatives()))
    }
    fn is_finite(self) -> bool {
        [self.value, self.first, self.second]
            .into_iter()
            .all(finite)
    }
}
impl core::ops::Add for Jet {
    type Output = Self;
    fn add(self, other: Self) -> Self {
        Self {
            value: self.value + other.value,
            first: self.first + other.first,
            second: self.second + other.second,
        }
    }
}
impl core::ops::Mul for Jet {
    type Output = Self;
    fn mul(self, other: Self) -> Self {
        Self {
            value: self.value * other.value,
            first: self.first * other.value + self.value * other.first,
            second: self.second * other.value
                + (self.first * other.first).scale(2.0)
                + self.value * other.second,
        }
    }
}
impl core::ops::Div for Jet {
    type Output = Self;
    fn div(self, other: Self) -> Self {
        let value = self.value / other.value;
        let first = (self.first - value * other.first) / other.value;
        Self {
            value,
            first,
            second: (self.second - value * other.second - (first * other.first).scale(2.0))
                / other.value,
        }
    }
}

fn reduced_jets(
    medium: WaveMedium,
    cosine: f64,
    constraint: IncidenceFrequencyConstraint,
    ei: [C64; 2],
    et: [C64; 2],
) -> (Jet, Jet) {
    let (di, dt) = normal_squared_derivatives(medium, cosine, constraint, ei, et);
    let qi = Jet::new(C64::from_re(cosine * cosine), di).root_at(C64::from_re(cosine));
    (qi, Jet::new(medium.normal_squared(cosine), dt))
}

fn admittance_jets(
    medium: WaveMedium,
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
    ei: [C64; 2],
    et: [C64; 2],
) -> Result<[Jet; 2]> {
    let (qi, q_squared) = reduced_jets(medium, cosine, constraint, ei, et);
    let q = medium.normal(cosine)?;
    if q == C64::ZERO {
        return Err(bad("critical exit has no finite frequency jet"));
    }
    let q = q_squared.root_at(q);
    let (a, b) = match polarization {
        Polarization::Te => (qi.scale(medium.mu_ratio), q),
        Polarization::Tm => (
            Jet::new(C64::ONE, ei) * q,
            Jet::new(medium.epsilon_ratio, et) * qi,
        ),
    };
    let normalization = positive(maximum(a.value).max(maximum(b.value)))?;
    let state = [a.divide_real(normalization), b.divide_real(normalization)];
    if state.into_iter().any(|value| !value.is_finite()) {
        return Err(bad("exterior admittance frequency jet is unrepresentable"));
    }
    Ok(state)
}

fn reflection_jet([a, b]: [Jet; 2]) -> Jet {
    let denominator = a.value + b.value;
    let first = ((a.first * b.value - a.value * b.first).scale(2.0) / denominator) / denominator;
    let second = ((a.second * b.value - a.value * b.second).scale(2.0) / denominator
        - (first * (a.first + b.first)).scale(2.0))
        / denominator;
    Jet::new((a.value - b.value) / denominator, [first, second])
}

fn slab_jets(
    slab: Slab,
    medium: WaveMedium,
    depth: f64,
    cosine: f64,
    polarization: Polarization,
    constraint: IncidenceFrequencyConstraint,
    ei: [C64; 2],
    et: [C64; 2],
) -> Result<[[Jet; 2]; 2]> {
    let (qi, q_squared) = reduced_jets(medium, cosine, constraint, ei, et);
    let a = match polarization {
        Polarization::Te => qi.scale(medium.mu_ratio),
        Polarization::Tm => Jet::new(medium.epsilon_ratio, et) * qi / Jet::new(C64::ONE, ei),
    };
    let depth = Jet::new(C64::from_re(depth), [C64::from_re(depth), C64::ZERO]);
    let imaginary_depth = Jet::constant(C64::new(0.0, 1.0)) * depth;
    let u = imaginary_depth * a;
    let v = imaginary_depth * q_squared / a;
    let w = depth * depth * q_squared;
    if !w.is_finite() || !u.is_finite() || !v.is_finite() || u.value == C64::ZERO {
        return Err(bad("film phase frequency jet is unrepresentable"));
    }
    let original_u = match polarization {
        Polarization::Te => slab.matrix[0][1],
        Polarization::Tm => slab.matrix[1][0],
    };
    let (sinc, sinc_w, sinc_ww) = if w.value.abs() <= 1.0 {
        let [s, first, second] = sinc_series(w.value);
        (slab.scale * s, slab.scale * first, slab.scale * second)
    } else {
        let s = original_u / u.value;
        let first = ((slab.matrix[0][0] - s) / w.value).scale(0.5);
        let second = ((s.scale(-0.5) - first.scale(3.0)) / w.value).scale(0.5);
        (s, first, second)
    };
    let cosine = Jet::new(
        slab.matrix[0][0],
        [
            sinc.scale(-0.5) * w.first,
            (sinc_w * w.first * w.first + sinc * w.second).scale(-0.5),
        ],
    );
    let sinc = Jet::new(
        sinc,
        [
            sinc_w * w.first,
            sinc_ww * w.first * w.first + sinc_w * w.second,
        ],
    );
    let (upper, lower) = match polarization {
        Polarization::Te => (u * sinc, v * sinc),
        Polarization::Tm => (v * sinc, u * sinc),
    };
    let mut matrix = [[cosine, upper], [lower, cosine]];
    // Preserve the exact native carrier, including its DD phase reduction and
    // reciprocal branch. Only derivatives use the entire-function identities.
    for row in 0..2 {
        for column in 0..2 {
            matrix[row][column].value = slab.matrix[row][column];
        }
    }
    if matrix.iter().flatten().any(|value| !value.is_finite()) {
        return Err(bad("film matrix frequency jet is unrepresentable"));
    }
    Ok(matrix)
}

// Entire sinc(sqrt(w)) and its first two w derivatives. The fixed polynomial
// includes w^10; at |w|<=1 its omitted second-derivative term is <5e-21.
#[inline]
pub(super) fn sinc_series(w: C64) -> [C64; 3] {
    const COEFFICIENTS: [f64; 11] = [
        1.0,
        -1.0 / 6.0,
        1.0 / 120.0,
        -1.0 / 5040.0,
        1.0 / 362880.0,
        -1.0 / 39916800.0,
        1.0 / 6227020800.0,
        -1.0 / 1307674368000.0,
        1.0 / 355687428096000.0,
        -1.0 / 121645100408832000.0,
        1.0 / 51090942171709440000.0,
    ];
    let mut value = C64::from_re(COEFFICIENTS[10]);
    let mut first = C64::ZERO;
    let mut second = C64::ZERO;
    for &coefficient in COEFFICIENTS[..10].iter().rev() {
        second = second * w + first.scale(2.0);
        first = first * w + value;
        value = value * w + C64::from_re(coefficient);
    }
    [value, first, second]
}
