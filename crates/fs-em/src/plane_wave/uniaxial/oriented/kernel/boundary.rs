//! Analytic vector boundary jets in the native carrier's frozen local gauge.
//!
//! The original normalized basis and solved coefficients are authoritative at
//! the carrier. Its normalization scale is held constant while differentiating
//! the raw Maxwell vectors; coefficient times basis is gauge invariant. The
//! two derivative solves reuse the same fixed-size carrier matrix.

use super::{
    BoundaryWitness, Input as NativeInput, OrientedCrystalResponse, Result, Vector, inverse_tensor,
    positive, raw_eigenvectors, unit,
};
use crate::plane_wave::dispersion::IncidenceFrequencyConstraint;
use crate::plane_wave::uniaxial::dispersion::kernel::{
    Input as NormalInput, Z, normal_index_derivatives,
};
use crate::plane_wave::uniaxial::{UniaxialMode, UniaxialWaveKind};
use fs_math::{c64::C64, dd::Dd};

const TWO: Dd = Dd { hi: 2.0, lo: 0.0 };
const HALF: Dd = Dd { hi: 0.5, lo: 0.0 };
const C0: f64 = 299_792_458.0;
type JVector = [Jet; 3];

pub(in crate::plane_wave::uniaxial::oriented) struct Input<'a> {
    pub original: NativeInput,
    pub boundary: &'a OrientedCrystalResponse,
    pub witness: &'a BoundaryWitness,
    /// [omega*epsilon'/epsilon, omega²*epsilon''/epsilon] of the three originals.
    pub relative_permittivity_derivatives: [[f64; 2]; 3],
    /// The exact original [K', omega*K''] of each principal material.
    pub principal_wavenumber_derivatives: [[f64; 2]; 2],
    pub constraint: IncidenceFrequencyConstraint,
}

#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::plane_wave::uniaxial::oriented) struct FieldDerivatives {
    /// [omega*E', omega²*E''] in original physical electric-field units.
    pub electric: [Vector; 2],
    /// [omega*H', omega²*H''] in original physical magnetic-field units.
    pub magnetic: [Vector; 2],
}
#[derive(Clone, Copy, Debug, PartialEq)]
pub(in crate::plane_wave::uniaxial::oriented) struct BoundaryDerivatives {
    pub incident: FieldDerivatives,
    pub reflected: FieldDerivatives,
    pub transmitted: [FieldDerivatives; 2],
}

pub(in crate::plane_wave::uniaxial::oriented) fn derivatives(
    data: Input<'_>,
) -> Result<BoundaryDerivatives> {
    let original = data.original;
    let basis = data.witness.basis;
    let axis = unit(original.axis)?;
    let relative = data.relative_permittivity_derivatives;
    if relative
        .iter()
        .flatten()
        .chain(data.principal_wavenumber_derivatives.iter().flatten())
        .any(|value| !value.is_finite())
    {
        return Err("finite original vector boundary frequency jets required");
    }
    if basis
        .modes
        .iter()
        .any(|mode| mode.kind == UniaxialWaveKind::Critical)
    {
        return Err("critical crystal boundaries have no admitted vector frequency jet");
    }
    if original.no == original.ne
        && axis[1] != 0.0
        && (relative[1] != relative[2]
            || data.principal_wavenumber_derivatives[0] != data.principal_wavenumber_derivatives[1])
    {
        return Err(
            "a degenerate out-of-plane crystal carrier requires matching principal frequency jets",
        );
    }

    // A common local index scaling preserves every electric boundary vector
    // and scales every magnetic vector by that same index law. Keep the exact
    // identity even when tiny slopes disappear from the original bulk K'.
    let fixed_direction = data.constraint == IncidenceFrequencyConstraint::FixedAngle
        || (relative[0] == [-2.0, 6.0] && data.principal_wavenumber_derivatives == [[0.0; 2]; 2]);
    if fixed_direction && relative[0] == relative[1] && relative[1] == relative[2] {
        let [r1, r2] = relative[0].map(Dd::from_f64);
        let factors = [r1 * HALF, r2 * HALF - r1 * r1 * Dd::from_f64(0.25)];
        let fields = |magnetic: Vector| -> Result<FieldDerivatives> {
            let mut result = FieldDerivatives {
                electric: [[C64::ZERO; 3]; 2],
                magnetic: [[C64::ZERO; 3]; 2],
            };
            for order in 0..2 {
                for component in 0..3 {
                    result.magnetic[order][component] = Z::from_c64(magnetic[component])
                        .scale(factors[order])
                        .finish()?;
                }
            }
            Ok(result)
        };
        return Ok(BoundaryDerivatives {
            incident: fields(data.boundary.incident_magnetic_a_m)?,
            reflected: fields(data.boundary.reflected_magnetic_a_m)?,
            transmitted: [
                fields(data.boundary.transmitted[0].magnetic_a_m)?,
                fields(data.boundary.transmitted[1].magnetic_a_m)?,
            ],
        });
    }

    let [i, o, e] = [
        basis.i,
        original.no / basis.scale,
        original.ne / basis.scale,
    ];
    let squares = [i, o, e].map(|index| Dd::from_f64(index) * Dd::from_f64(index));
    let epsilon: [Jet; 3] = core::array::from_fn(|j| {
        Jet::real_jet(
            squares[j],
            relative[j].map(|value| squares[j] * Dd::from_f64(value)),
        )
    });
    let p0 = Dd::from_f64(basis.p);
    let qi0 = Dd::from_f64(basis.qi);
    let [ri1, ri2] = relative[0].map(Dd::from_f64);
    let (p, qi) = if data.constraint == IncidenceFrequencyConstraint::FixedAngle {
        let second = ri2 * HALF - ri1 * ri1 * Dd::from_f64(0.25);
        (
            Jet::real_jet(p0, [p0 * ri1 * HALF, p0 * second]),
            Jet::real_jet(qi0, [qi0 * ri1 * HALF, qi0 * second]),
        )
    } else {
        // Factor against the original normal coordinate. In particular
        // n proportional to 1/omega gives qi'=-qi and qi''=2qi even when
        // an authored cos² is too small to survive subtraction from one.
        let first = squares[0] * (ri1 + TWO) * HALF / qi0 - qi0;
        let second = squares[0] * (ri2 - Dd::from_f64(6.0)) * HALF / qi0 + qi0 * Dd::from_f64(3.0)
            - first * (first / qi0);
        (
            Jet::real_jet(p0, [-p0, p0 * TWO]),
            Jet::real_jet(qi0, [first, second]),
        )
    };
    let matched_isotropic = original.ni == original.no
        && original.no == original.ne
        && relative[0] == relative[1]
        && relative[1] == relative[2]
        && data.principal_wavenumber_derivatives[0] == data.principal_wavenumber_derivatives[1];
    let matched_ordinary =
        axis[1] == 0.0 && original.ni == original.no && relative[0] == relative[1];
    if matched_isotropic
        || (matched_ordinary && (original.electric[0] == C64::ZERO || axis[2] == 0.0))
    {
        // Locally identical phase laws retain the same angle through the
        // interface. Their TE/TM reflection is the constant permeability
        // contrast, even under fixed kx. The same identity independently
        // applies to ordinary TE in the exact principal plane.
        // Do not recover its derivative by
        // subtracting A' times a rounded native coefficient from b': at
        // grazing that residual would be amplified as 1/cos².
        let ratio = Dd::from_f64(basis.ratio_mu);
        let denominator = ratio + Dd::ONE;
        let reflection = (ratio - Dd::ONE) / denominator;
        let electric_transmission = (ratio / denominator) * TWO;
        let magnetic_transmission = TWO / denominator;
        let amplitude = data.witness.amplitude;
        let input = original
            .electric
            .map(|z| C64::new(z.re / amplitude, z.im / amplitude));
        let ex = Jet::constant(input[0]);
        let ey = Jet::constant(input[1]);
        let te_e = [Jet::zero(), ey, Jet::zero()];
        let te_h = [
            qi.mul(ey).scale(-ratio),
            Jet::zero(),
            p.mul(ey).scale(ratio),
        ];
        // An unexcited TM field is identically zero. Its unused longitudinal
        // quotient can have an unrepresentable second derivative while the
        // physical TE derivatives are still finite at authored grazing.
        let (tm_e, tm_h, reflected_tm_e, reflected_tm_h, transmitted_tm_e, transmitted_tm_h) =
            if input[0] == C64::ZERO {
                (
                    [Jet::zero(); 3],
                    [Jet::zero(); 3],
                    [Jet::zero(); 3],
                    [Jet::zero(); 3],
                    [Jet::zero(); 3],
                    [Jet::zero(); 3],
                )
            } else {
                let (tm_reflection, tm_transmission, index_ratio) = if matched_isotropic {
                    (
                        Jet::real_jet(reflection, [Dd::ZERO; 2]),
                        Jet::real_jet(electric_transmission, [Dd::ZERO; 2]),
                        Jet::one(),
                    )
                } else {
                    // For an optic axis along x and a locally matched ordinary
                    // law, q_e=(ne/no)*qi. Cancel this common cosine BEFORE
                    // differentiating the transverse Fresnel amplitudes.
                    // Otherwise finite material slopes are lost between
                    // normal-root derivatives growing as 1/cos.
                    let index = |value: Dd, r: [f64; 2]| {
                        let [r1, r2] = r.map(Dd::from_f64);
                        Jet::real_jet(
                            value,
                            [
                                value * r1 * HALF,
                                value * (r2 * HALF - r1 * r1 * Dd::from_f64(0.25)),
                            ],
                        )
                    };
                    let ordinary = index(Dd::from_f64(i), relative[0]);
                    let extraordinary = index(Dd::from_f64(e), relative[2]);
                    let numerator = ordinary.scale(ratio);
                    let denominator = numerator.add(extraordinary);
                    (
                        numerator.sub(extraordinary).divide(denominator)?,
                        numerator.scale(TWO).divide(denominator)?,
                        extraordinary.divide(ordinary)?,
                    )
                };
                // Weight the longitudinal quotients by the actual excitation
                // before dividing, retaining representable weak TM fields.
                let longitudinal = p.mul(ex).divide(qi)?.neg();
                let admittance = epsilon[0].mul(ex).divide(qi)?;
                let e = [ex, Jet::zero(), longitudinal];
                let h = [Jet::zero(), admittance.scale(ratio), Jet::zero()];
                let reflected_e = [
                    ex.mul(tm_reflection),
                    Jet::zero(),
                    longitudinal.mul(tm_reflection).neg(),
                ];
                let reflected_h = [Jet::zero(), h[1].mul(tm_reflection).neg(), Jet::zero()];
                let transmitted_e = [
                    ex.mul(tm_transmission),
                    Jet::zero(),
                    longitudinal.mul(index_ratio).mul(tm_transmission),
                ];
                let transmitted_h = [
                    Jet::zero(),
                    admittance.mul(index_ratio).mul(tm_transmission),
                    Jet::zero(),
                ];
                (e, h, reflected_e, reflected_h, transmitted_e, transmitted_h)
            };
        let incident_e = core::array::from_fn(|j| te_e[j].add(tm_e[j]));
        let incident_h = core::array::from_fn(|j| te_h[j].add(tm_h[j]));
        let reflected_e = [reflected_tm_e[0], ey.scale(reflection), reflected_tm_e[2]];
        let reflected_h = [
            te_h[0].scale(-reflection),
            reflected_tm_h[1],
            te_h[2].scale(reflection),
        ];
        let hscale = positive(basis.scale / (original.mut_ * C0))?;
        return Ok(BoundaryDerivatives {
            incident: finish_fields(incident_e, incident_h, amplitude, hscale)?,
            reflected: finish_fields(reflected_e, reflected_h, amplitude, hscale)?,
            transmitted: [
                finish_fields(
                    te_e.map(|v| v.scale(electric_transmission)),
                    te_h.map(|v| v.scale(magnetic_transmission)),
                    amplitude,
                    hscale,
                )?,
                finish_fields(transmitted_tm_e, transmitted_tm_h, amplitude, hscale)?,
            ],
        });
    }

    let inverse_o = Jet::one().divide(epsilon[1])?;
    let inverse_e = Jet::one().divide(epsilon[2])?;
    let mut electric = [[Jet::zero(); 3]; 2];
    let mut magnetic = [[Jet::zero(); 3]; 2];
    for mode in 0..2 {
        let native = basis.modes[mode];
        let normal = normal_index_derivatives(
            NormalInput {
                omega: original.omega,
                indices: [original.ni, original.no, original.ne],
                relative_permittivity_derivatives: relative,
                principal_wavenumber_derivatives: data.principal_wavenumber_derivatives,
                incidence: original.incidence,
                axis: original.axis,
                mode: if mode == 0 {
                    UniaxialMode::Ordinary
                } else {
                    UniaxialMode::Extraordinary
                },
                kind: native.kind,
                normal_wavenumber: data.boundary.transmitted[mode].normal_wavenumber_rad_m,
                constraint: data.constraint,
            },
            native.z,
        )?;
        let normal = if original.ni == original.no
            && relative[0] == relative[1]
            && (mode == 0 || (original.no == original.ne && relative[1] == relative[2]))
        {
            // Exact matched scalar roots share the incident coordinate even
            // at authored grazing; do not subtract independent rounded jets.
            qi.with_value(native.z)
        } else {
            Jet {
                value: Z::from_c64(native.z),
                first: normal[0],
                second: normal[1],
            }
        };
        let q = [p, Jet::zero(), normal];
        let (mut ej, mut hj) = eigenvectors(
            axis,
            q,
            mode == 0,
            original.no == original.ne,
            inverse_o,
            inverse_e,
        );
        let (raw_e, raw_h) = raw_eigenvectors(
            axis,
            [C64::from_re(basis.p), C64::ZERO, native.z],
            mode == 0,
            original.no == original.ne,
            |v| inverse_tensor(v, axis, o * o, e * e),
        );
        let gauge = positive(
            raw_e
                .iter()
                .chain(&raw_h)
                .fold(0.0_f64, |value, z| value.max(z.abs())),
        )?;
        let inverse_gauge = Dd::from_f64(1.0 / gauge);
        for j in 0..3 {
            // Zeroth order is exactly the original normalized native basis,
            // while the gauge is frozen for the two derivative orders.
            ej[j] = ej[j].scale(inverse_gauge).with_value(native.e[j]);
            hj[j] = hj[j].scale(inverse_gauge).with_value(native.h[j]);
        }
        electric[mode] = ej;
        magnetic[mode] = hj;
    }

    let ratio = Dd::from_f64(basis.ratio_mu);
    let mut matrix = [[Jet::zero(); 2]; 2];
    for mode in 0..2 {
        matrix[0][mode] = magnetic[mode][1]
            .mul(qi)
            .add(electric[mode][0].mul(epsilon[0]).scale(ratio))
            .with_value(
                basis.modes[mode].h[1].scale(basis.qi)
                    + basis.modes[mode].e[0].scale(basis.ratio_mu * basis.i * basis.i),
            );
        matrix[1][mode] = magnetic[mode][0]
            .neg()
            .add(electric[mode][1].mul(qi).scale(ratio))
            .with_value(
                -basis.modes[mode].h[0] + basis.modes[mode].e[1].scale(basis.ratio_mu * basis.qi),
            );
    }
    let amplitude = data.witness.amplitude;
    let incident = original
        .electric
        .map(|z| C64::new(z.re / amplitude, z.im / amplitude));
    let mut rhs = [
        epsilon[0]
            .mul(Jet::constant(incident[0]))
            .scale(ratio * TWO),
        qi.mul(Jet::constant(incident[1])).scale(ratio * TWO),
    ];
    if matched_ordinary {
        // In the principal plane the ordinary row contains only the TE
        // coefficient. Cancel its exact common normal root analytically:
        // that coefficient is constant in the frozen electric-field gauge.
        // This also applies to the TE part of mixed incident polarization.
        for entry in &mut matrix[1] {
            entry.first = Z::real(Dd::ZERO);
            entry.second = Z::real(Dd::ZERO);
        }
        rhs[1].first = Z::real(Dd::ZERO);
        rhs[1].second = Z::real(Dd::ZERO);
    }
    for row in 0..2 {
        let row_scale = positive(
            matrix[row][0]
                .value
                .finish()?
                .abs()
                .max(matrix[row][1].value.finish()?.abs()),
        )?;
        let scale = Dd::from_f64(1.0 / row_scale);
        matrix[row] = matrix[row].map(|entry| entry.scale(scale));
        rhs[row] = rhs[row].scale(scale);
    }
    let carrier = data.witness.coefficients.map(Z::from_c64);
    let a0 = matrix.map(|row| row.map(|entry| entry.value));
    let first_rhs = core::array::from_fn(|row| {
        rhs[row]
            .first
            .sub(matrix[row][0].first.mul(carrier[0]))
            .sub(matrix[row][1].first.mul(carrier[1]))
    });
    let first = solve(a0, first_rhs)?;
    let second_rhs = core::array::from_fn(|row| {
        rhs[row]
            .second
            .sub(matrix[row][0].second.mul(carrier[0]))
            .sub(matrix[row][1].second.mul(carrier[1]))
            .sub(matrix[row][0].first.mul(first[0]).scale(TWO))
            .sub(matrix[row][1].first.mul(first[1]).scale(TWO))
    });
    let second = solve(a0, second_rhs)?;
    for mode in 0..2 {
        let coefficient = Jet {
            value: carrier[mode],
            first: first[mode],
            second: second[mode],
        };
        for j in 0..3 {
            electric[mode][j] = electric[mode][j]
                .mul(coefficient)
                .with_value(basis.modes[mode].e[j] * data.witness.coefficients[mode]);
            magnetic[mode][j] = magnetic[mode][j]
                .mul(coefficient)
                .with_value(basis.modes[mode].h[j] * data.witness.coefficients[mode]);
        }
    }

    if matched_ordinary {
        let denominator = ratio + Dd::ONE;
        let ey = Jet::constant(incident[1]);
        electric[0] = [
            Jet::zero(),
            ey.scale((ratio / denominator) * TWO),
            Jet::zero(),
        ];
        let transmission = (ratio / denominator) * TWO;
        magnetic[0] = [
            qi.mul(ey).scale(-transmission),
            Jet::zero(),
            p.mul(ey).scale(transmission),
        ];
    }

    let incident_e = [
        Jet::constant(incident[0]),
        Jet::constant(incident[1]),
        p.divide(qi)?
            .mul(Jet::constant(incident[0]))
            .neg()
            .with_value(incident[0].scale(-basis.p / basis.qi)),
    ];
    let incident_h = cross([p, Jet::zero(), qi], incident_e).map(|value| value.scale(ratio));
    let reflected_carrier = core::array::from_fn::<_, 2, _>(|j| {
        basis.modes[0].e[j] * data.witness.coefficients[0]
            + basis.modes[1].e[j] * data.witness.coefficients[1]
            - incident[j]
    });
    let reflected_x = electric[0][0]
        .add(electric[1][0])
        .sub(incident_e[0])
        .with_value(reflected_carrier[0]);
    let reflected_y = if matched_ordinary {
        // The exact constant TE reflection must not inherit a rounded
        // normalized-basis/coefficient product residual from the carrier.
        Jet::constant(incident[1]).scale((ratio - Dd::ONE) / (ratio + Dd::ONE))
    } else {
        electric[0][1]
            .add(electric[1][1])
            .sub(incident_e[1])
            .with_value(reflected_carrier[1])
    };
    let reflected_e = [
        reflected_x,
        reflected_y,
        p.divide(qi)?
            .mul(reflected_x)
            .with_value(reflected_carrier[0].scale(basis.p / basis.qi)),
    ];
    let reflected_h =
        cross([p, Jet::zero(), qi.neg()], reflected_e).map(|value| value.scale(ratio));
    let hscale = positive(basis.scale / (original.mut_ * C0))?;
    let fields = |e: JVector, h: JVector| finish_fields(e, h, amplitude, hscale);
    Ok(BoundaryDerivatives {
        incident: fields(incident_e, incident_h)?,
        reflected: fields(reflected_e, reflected_h)?,
        transmitted: [
            fields(electric[0], magnetic[0])?,
            fields(electric[1], magnetic[1])?,
        ],
    })
}

fn eigenvectors(
    axis: [f64; 3],
    q: JVector,
    ordinary: bool,
    isotropic: bool,
    inverse_o: Jet,
    inverse_e: Jet,
) -> (JVector, JVector) {
    let a = axis.map(|value| Jet::constant(C64::from_re(value)));
    let inverse = |v: JVector| {
        let projection = v[0]
            .scale(Dd::from_f64(axis[0]))
            .add(v[1].scale(Dd::from_f64(axis[1])))
            .add(v[2].scale(Dd::from_f64(axis[2])));
        core::array::from_fn(|j| {
            v[j].mul(inverse_o).add(
                projection
                    .mul(inverse_e.sub(inverse_o))
                    .scale(Dd::from_f64(axis[j])),
            )
        })
    };
    if isotropic || axis[1] == 0.0 {
        if ordinary {
            let e = [Jet::zero(), Jet::one(), Jet::zero()];
            (e, cross(q, e))
        } else {
            let h = [Jet::zero(), Jet::one(), Jet::zero()];
            (inverse(cross(q, h).map(Jet::neg)), h)
        }
    } else if ordinary {
        let e = cross(q, a);
        (e, cross(q, e))
    } else {
        let h = cross(q, a);
        (inverse(cross(q, h).map(Jet::neg)), h)
    }
}

fn solve(matrix: [[Z; 2]; 2], rhs: [Z; 2]) -> Result<[Z; 2]> {
    let determinant = matrix[0][0]
        .mul(matrix[1][1])
        .sub(matrix[0][1].mul(matrix[1][0]));
    Ok([
        rhs[0]
            .mul(matrix[1][1])
            .sub(rhs[1].mul(matrix[0][1]))
            .divide(determinant)?,
        matrix[0][0]
            .mul(rhs[1])
            .sub(matrix[1][0].mul(rhs[0]))
            .divide(determinant)?,
    ])
}

fn cross(a: JVector, b: JVector) -> JVector {
    [
        a[1].mul(b[2]).sub(a[2].mul(b[1])),
        a[2].mul(b[0]).sub(a[0].mul(b[2])),
        a[0].mul(b[1]).sub(a[1].mul(b[0])),
    ]
}

fn finish_fields(e: JVector, h: JVector, amplitude: f64, hscale: f64) -> Result<FieldDerivatives> {
    let mut electric = [[C64::ZERO; 3]; 2];
    let mut magnetic = [[C64::ZERO; 3]; 2];
    for j in 0..3 {
        for order in 0..2 {
            electric[order][j] = [e[j].first, e[j].second][order]
                .scale(Dd::from_f64(amplitude))
                .finish()?;
            magnetic[order][j] = [h[j].first, h[j].second][order]
                .scale(Dd::from_f64(amplitude))
                .scale(Dd::from_f64(hscale))
                .finish()?;
        }
    }
    Ok(FieldDerivatives { electric, magnetic })
}

/// Ordinary derivatives, scaled by the one fixed carrier frequency.
#[derive(Clone, Copy)]
struct Jet {
    value: Z,
    first: Z,
    second: Z,
}
impl Jet {
    fn constant(value: C64) -> Self {
        Self {
            value: Z::from_c64(value),
            first: Z::real(Dd::ZERO),
            second: Z::real(Dd::ZERO),
        }
    }
    fn zero() -> Self {
        Self::constant(C64::ZERO)
    }
    fn one() -> Self {
        Self::constant(C64::ONE)
    }
    fn real_jet(value: Dd, derivatives: [Dd; 2]) -> Self {
        Self {
            value: Z::real(value),
            first: Z::real(derivatives[0]),
            second: Z::real(derivatives[1]),
        }
    }
    fn with_value(self, value: C64) -> Self {
        Self {
            value: Z::from_c64(value),
            ..self
        }
    }
    fn add(self, other: Self) -> Self {
        Self {
            value: self.value.add(other.value),
            first: self.first.add(other.first),
            second: self.second.add(other.second),
        }
    }
    fn sub(self, other: Self) -> Self {
        Self {
            value: self.value.sub(other.value),
            first: self.first.sub(other.first),
            second: self.second.sub(other.second),
        }
    }
    fn scale(self, scale: Dd) -> Self {
        Self {
            value: self.value.scale(scale),
            first: self.first.scale(scale),
            second: self.second.scale(scale),
        }
    }
    fn neg(self) -> Self {
        self.scale(-Dd::ONE)
    }
    fn mul(self, other: Self) -> Self {
        Self {
            value: self.value.mul(other.value),
            first: self.first.mul(other.value).add(self.value.mul(other.first)),
            second: self
                .second
                .mul(other.value)
                .add(self.first.mul(other.first).scale(TWO))
                .add(self.value.mul(other.second)),
        }
    }
    fn divide(self, other: Self) -> Result<Self> {
        let value = self.value.divide(other.value)?;
        let first = self.first.sub(value.mul(other.first)).divide(other.value)?;
        let second = self
            .second
            .sub(value.mul(other.second))
            .sub(first.mul(other.first).scale(TWO))
            .divide(other.value)?;
        Ok(Self {
            value,
            first,
            second,
        })
    }
}

#[cfg(test)]
mod tests;
