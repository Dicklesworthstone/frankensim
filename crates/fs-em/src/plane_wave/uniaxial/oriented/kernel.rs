//! Fixed-size vector Maxwell boundary kernel; no source resolution or allocation.
use super::super::kernel::{
    Incidence, UniaxialWaveKind, binary_scale, dispersive_group_denominator, outgoing_root,
};
use fs_math::{c64::C64, dd::Dd, det};

const C0: f64 = 299_792_458.0;
type Vector = [C64; 3];
type Result<T> = core::result::Result<T, &'static str>;

pub(super) mod boundary;

#[path = "propagation.rs"]
mod propagation;
pub(super) use propagation::{FieldSum, WaveGeometry, field};
pub use propagation::{OrientedCrystalField, OrientedCrystalFieldSide};

/// Actual transmitted field of one outgoing eigenmode at the interface.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedCrystalMode {
    pub wave_kind: UniaxialWaveKind,
    pub normal_wavenumber_rad_m: C64,
    /// Full Cartesian electric field at z=0 [V/m], with the incident phase.
    pub electric_v_m: Vector,
    /// Full Cartesian magnetic field at z=0 [A/m].
    pub magnetic_a_m: Vector,
    pub normal_power_w_m2: f64,
    pub phase_direction: Option<[f64; 3]>,
    pub energy_direction: Option<[f64; 3]>,
}

/// Complete coherent boundary fields for the caller's incident polarization.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct OrientedCrystalResponse {
    pub tangential_wavenumber_rad_m: f64,
    pub incident_normal_wavenumber_rad_m: f64,
    pub incident_electric_v_m: Vector,
    pub incident_magnetic_a_m: Vector,
    pub reflected_electric_v_m: Vector,
    pub reflected_magnetic_a_m: Vector,
    /// Ordinary then extraordinary. In the exactly isotropic limit these
    /// degenerate modes use the deterministic TE then TM basis.
    pub transmitted: [OrientedCrystalMode; 2],
    pub incident_normal_power_w_m2: f64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
}

#[derive(Clone, Copy)]
pub(super) struct Input {
    pub(super) omega: f64,
    pub(super) ni: f64,
    pub(super) no: f64,
    pub(super) ne: f64,
    pub(super) mui: f64,
    pub(super) mut_: f64,
    pub(super) axis: [f64; 3],
    pub(super) incidence: Incidence,
    pub(super) electric: [C64; 2],
}
#[derive(Clone, Copy)]
pub(super) struct Mode {
    pub(super) z: C64,
    pub(super) kind: UniaxialWaveKind,
    pub(super) e: Vector,
    pub(super) h: Vector,
}

fn invalid(message: &'static str) -> &'static str {
    message
}
fn positive(v: f64) -> Result<f64> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err(invalid(
            "oriented crystal requires positive representable wave data",
        ))
    }
}
fn finite(v: C64) -> bool {
    v.re.is_finite() && v.im.is_finite()
}
pub(super) fn unit(v: [f64; 3]) -> Result<[f64; 3]> {
    if v.iter().any(|x| !x.is_finite()) {
        return Err(invalid("finite nonzero optic axis required"));
    }
    let scale = positive(v.iter().fold(0.0_f64, |m, x| m.max(x.abs())))?;
    let v = v.map(|x| x / scale);
    let length = det::hypot(det::hypot(v[0], v[1]), v[2]);
    Ok(v.map(|x| if x == 0.0 { 0.0 } else { x / length }))
}
pub(super) fn cross(a: Vector, b: Vector) -> Vector {
    [
        a[1] * b[2] - a[2] * b[1],
        a[2] * b[0] - a[0] * b[2],
        a[0] * b[1] - a[1] * b[0],
    ]
}
fn power(e: Vector, h: Vector) -> [f64; 3] {
    cross(e, h.map(C64::conj)).map(|z| 0.5 * z.re)
}
/// The original two outgoing eigenvectors before interface amplitudes.
/// Shared by the half-space and finite-slab boundary owners.
#[derive(Clone, Copy)]
pub(super) struct Eigenbasis {
    pub(super) i: f64,
    pub(super) p: f64,
    pub(super) qi: f64,
    pub(super) scale: f64,
    pub(super) ratio_mu: f64,
    pub(super) modes: [Mode; 2],
}

/// Exact native carrier gauge retained for analytic boundary derivatives.
/// Fixed-size only: no material resolution, root reconstruction or field solve.
#[derive(Clone, Copy)]
pub(super) struct BoundaryWitness {
    pub(super) basis: Eigenbasis,
    pub(super) coefficients: [C64; 2],
    pub(super) amplitude: f64,
    pub(super) incident_electric: [C64; 2],
}

pub(super) fn eigenbasis(data: Input) -> Result<Eigenbasis> {
    let axis = unit(data.axis)?;
    let scale = binary_scale(positive(data.ni.max(data.no).max(data.ne))?);
    let (i, o, e) = (
        positive(data.ni / scale)?,
        positive(data.no / scale)?,
        positive(data.ne / scale)?,
    );
    let ratio_mu = positive(data.mut_ / data.mui)?;
    let od = Dd::from_f64(o);
    let ed = Dd::from_f64(e);
    let o2 = od * od;
    let e2 = ed * ed;
    let (p, p2, qi) = incident_coordinates(data.ni, scale, i, data.incidence)?;
    let (zo, ko) = if data.ni == data.no {
        (C64::from_re(qi), UniaxialWaveKind::Propagating)
    } else {
        outgoing_root(o2 - p2).map_err(invalid)?
    };
    let (ze, ke) = if data.no == data.ne {
        (zo, ko)
    } else {
        let [ax, ay, az] = axis.map(Dd::from_f64);
        let axis2 = ax * ax + ay * ay + az * az;
        // Normalize direction weights BEFORE the positive constitutive sums.
        // Dividing a rounded weighted numerator afterwards can leave a tiny
        // signed residual at an exactly critical symmetric principal section.
        // This matches the original principal-plane kernel without a cutoff
        // tolerance or cancellation of a small index against a large one.
        let xx = ax * ax / axis2;
        let yy = ay * ay / axis2;
        let zz = az * az / axis2;
        let xz = ax * az / axis2;
        let f = o2 * (xx + yy) + e2 * zz;
        let g = (e2 - o2) * xz;
        let h = o2 * (yy + zz) + e2 * xx;
        let gp = g * Dd::from_f64(p);
        let constant = h * p2 - o2 * e2;
        // Factor f*h-g*g before applying p². This avoids subtracting
        // nearly equal shifted-root terms and keeps one conserved p².
        let determinant = o2 * (o2 * yy + e2 * (xx + zz));
        let discriminant = if axis[1] == 0.0 {
            o2 * e2 * (f - p2)
        } else {
            f * o2 * e2 - determinant * p2
        };
        let id = Dd::from_f64(i);
        let i2 = id * id;
        let matched_cutoff = if axis[1] == 0.0 {
            f == i2
        } else {
            f * o2 * e2 == determinant * i2
        };
        let (s, k) = if matched_cutoff {
            // The original incident normal coordinate remains positive even
            // when 1-cos(theta)^2 rounds to one. At an exactly matched cutoff,
            // D=det*qi^2; retain sqrt(D)=sqrt(det)*qi without squaring the
            // grazing coordinate or inventing a coalesced critical mode.
            let factor = if axis[1] == 0.0 {
                od * ed
            } else {
                determinant.sqrt()
            };
            let root = positive((factor * Dd::from_f64(qi)).to_f64())?;
            (C64::from_re(root), UniaxialWaveKind::Propagating)
        } else {
            outgoing_root(discriminant).map_err(invalid)?
        };
        let real = if k == UniaxialWaveKind::Propagating && Dd::ZERO.lt(gp) {
            (constant / (-gp - Dd::from_f64(s.re))).to_f64()
        } else {
            ((-gp + Dd::from_f64(s.re)) / f).to_f64()
        };
        (C64::new(real, s.im / f.to_f64()), k)
    };
    let mut modes = [
        Mode {
            z: zo,
            kind: ko,
            e: [C64::ZERO; 3],
            h: [C64::ZERO; 3],
        },
        Mode {
            z: ze,
            kind: ke,
            e: [C64::ZERO; 3],
            h: [C64::ZERO; 3],
        },
    ];
    for (m, mode) in modes.iter_mut().enumerate() {
        (mode.e, mode.h) = eigenvectors(
            axis,
            [C64::from_re(p), C64::ZERO, mode.z],
            m == 0,
            data.no == data.ne,
            |v| inverse_tensor(v, axis, o * o, e * e),
        )?;
    }
    Ok(Eigenbasis {
        i,
        p,
        qi,
        scale,
        ratio_mu,
        modes,
    })
}

/// Preserve authored grazing cosines and signed tangential indices. Shared
/// by the real and complex constitutive branches before their dispersion.
pub(super) fn incident_coordinates(
    incident_index: f64,
    scale: f64,
    i: f64,
    incidence: Incidence,
) -> Result<(f64, Dd, f64)> {
    let id = Dd::from_f64(i);
    Ok(match incidence {
        Incidence::Cosine(c) => {
            positive(c)?;
            if c > 1.0 {
                return Err(invalid("crystal incidence cosine must lie in (0,1]"));
            }
            let cd = Dd::from_f64(c);
            let s2 = (Dd::ONE - cd) * (Dd::ONE + cd);
            (i * det::sqrt(s2.to_f64()), id * id * s2, positive(i * c)?)
        }
        Incidence::TangentialIndex(value) => {
            if !value.is_finite() || value.abs() >= incident_index {
                return Err(invalid(
                    "tangential index needs a propagating incident port",
                ));
            }
            let p = value / scale;
            if value != 0.0 && p == 0.0 {
                return Err(invalid("crystal tangential index is unrepresentable"));
            }
            let pd = Dd::from_f64(p);
            let p2 = pd * pd;
            (p, p2, positive((id * id - p2).sqrt().to_f64())?)
        }
    })
}

/// Native Maxwell eigenvectors, shared across transparent and passive
/// complex principal tensors. The caller owns only the inverse constitutive
/// tensor; exact principal-section/isotropic degeneracies use TE then TM.
pub(super) fn eigenvectors(
    axis: [f64; 3],
    q: Vector,
    ordinary: bool,
    isotropic: bool,
    inverse: impl Fn(Vector) -> Vector,
) -> Result<(Vector, Vector)> {
    let (e, h) = raw_eigenvectors(axis, q, ordinary, isotropic, inverse);
    let basis_scale = positive(e.iter().chain(&h).fold(0.0_f64, |v, z| v.max(z.abs())))?;
    Ok((
        e.map(|z| z.scale(1.0 / basis_scale)),
        h.map(|z| z.scale(1.0 / basis_scale)),
    ))
}

fn raw_eigenvectors(
    axis: [f64; 3],
    q: Vector,
    ordinary: bool,
    isotropic: bool,
    inverse: impl Fn(Vector) -> Vector,
) -> (Vector, Vector) {
    let a = axis.map(C64::from_re);
    if isotropic || axis[1] == 0.0 {
        if ordinary {
            let e = [C64::ZERO, C64::ONE, C64::ZERO];
            (e, cross(q, e))
        } else {
            let h = [C64::ZERO, C64::ONE, C64::ZERO];
            (inverse(cross(q, h).map(|z| -z)), h)
        }
    } else if ordinary {
        let e = cross(q, a);
        (e, cross(q, e))
    } else {
        let h = cross(q, a);
        (inverse(cross(q, h).map(|z| -z)), h)
    }
}

pub(super) fn solve(data: Input) -> Result<OrientedCrystalResponse> {
    solve_into(data, None)
}

pub(super) fn solve_with_witness(
    data: Input,
) -> Result<(OrientedCrystalResponse, BoundaryWitness)> {
    let mut witness = None;
    let boundary = solve_into(data, Some(&mut witness))?;
    Ok((
        boundary,
        witness.ok_or("native crystal boundary witness is missing")?,
    ))
}

fn solve_into(
    data: Input,
    capture: Option<&mut Option<BoundaryWitness>>,
) -> Result<OrientedCrystalResponse> {
    if data.electric.iter().any(|z| !finite(*z)) {
        return Err(invalid("finite crystal incident field required"));
    }
    let amplitude = positive(
        data.electric
            .iter()
            .fold(0.0_f64, |m, z| m.max(z.re.abs()).max(z.im.abs())),
    )?;
    let ei = data
        .electric
        .map(|z| C64::new(z.re / amplitude, z.im / amplitude));
    let basis = eigenbasis(data)?;
    let Eigenbasis {
        i,
        p,
        qi,
        scale,
        ratio_mu,
        mut modes,
    } = basis;
    let ko = modes[0].kind;
    // Eliminate the reflected field. The two equations below are the exact
    // tangential Maxwell conditions with the incident admittances multiplied
    // through, so no 1/qi or critical transmitted admittance is formed.
    let mut matrix = [[C64::ZERO; 2]; 2];
    for (m, mode) in modes.iter().enumerate() {
        matrix[0][m] = mode.h[1].scale(qi) + mode.e[0].scale(ratio_mu * i * i);
        matrix[1][m] = -mode.h[0] + mode.e[1].scale(ratio_mu * qi);
    }
    let mut rhs = [
        ei[0].scale(2.0 * ratio_mu * i * i),
        ei[1].scale(2.0 * ratio_mu * qi),
    ];
    for row in 0..2 {
        let row_scale = positive(matrix[row][0].abs().max(matrix[row][1].abs()))?;
        matrix[row] = matrix[row].map(|z| z.scale(1.0 / row_scale));
        rhs[row] = rhs[row].scale(1.0 / row_scale);
    }
    let determinant = matrix[0][0] * matrix[1][1] - matrix[0][1] * matrix[1][0];
    positive(determinant.abs())?;
    let coefficients = [
        (rhs[0] * matrix[1][1] - rhs[1] * matrix[0][1]) / determinant,
        (matrix[0][0] * rhs[1] - matrix[1][0] * rhs[0]) / determinant,
    ];
    for m in 0..2 {
        modes[m].e = modes[m].e.map(|z| z * coefficients[m]);
        modes[m].h = modes[m].h.map(|z| z * coefficients[m]);
    }
    let te = core::array::from_fn(|j| modes[0].e[j] + modes[1].e[j]);
    let th = core::array::from_fn(|j| modes[0].h[j] + modes[1].h[j]);
    let incident_e = [ei[0], ei[1], ei[0].scale(-p / qi)];
    let incident_h = cross([C64::from_re(p), C64::ZERO, C64::from_re(qi)], incident_e)
        .map(|z| z.scale(ratio_mu));
    let reflected_e = [te[0] - ei[0], te[1] - ei[1], (te[0] - ei[0]).scale(p / qi)];
    let reflected_h = cross([C64::from_re(p), C64::ZERO, C64::from_re(-qi)], reflected_e)
        .map(|z| z.scale(ratio_mu));
    let incoming = positive(power(incident_e, incident_h)[2])?;
    let reflected = -power(reflected_e, reflected_h)[2] / incoming;
    let transmitted = power(te, th)[2] / incoming;
    if !reflected.is_finite()
        || !transmitted.is_finite()
        || transmitted < -1e-12
        || (reflected + transmitted - 1.0).abs() > 4096.0 * f64::EPSILON
    {
        return Err(invalid(
            "oriented crystal Maxwell normal-power balance failed",
        ));
    }
    let hscale = positive(scale / (data.mut_ * C0))?;
    let kscale = positive((data.omega / C0) * scale)?;
    let mut result = OrientedCrystalResponse {
        tangential_wavenumber_rad_m: kscale * p,
        incident_normal_wavenumber_rad_m: kscale * qi,
        incident_electric_v_m: incident_e.map(|z| z.scale(amplitude)),
        incident_magnetic_a_m: incident_h.map(|z| z.scale(amplitude * hscale)),
        reflected_electric_v_m: reflected_e.map(|z| z.scale(amplitude)),
        reflected_magnetic_a_m: reflected_h.map(|z| z.scale(amplitude * hscale)),
        transmitted: [OrientedCrystalMode {
            wave_kind: ko,
            normal_wavenumber_rad_m: C64::ZERO,
            electric_v_m: [C64::ZERO; 3],
            magnetic_a_m: [C64::ZERO; 3],
            normal_power_w_m2: 0.0,
            phase_direction: None,
            energy_direction: None,
        }; 2],
        incident_normal_power_w_m2: incoming * amplitude * amplitude * hscale,
        reflected_power_fraction: reflected,
        transmitted_power_fraction: transmitted.max(0.0),
    };
    for (m, mode) in modes.iter().enumerate() {
        let s = power(mode.e, mode.h);
        result.transmitted[m] = OrientedCrystalMode {
            wave_kind: mode.kind,
            normal_wavenumber_rad_m: mode.z.scale(kscale),
            electric_v_m: mode.e.map(|z| z.scale(amplitude)),
            magnetic_a_m: mode.h.map(|z| z.scale(amplitude * hscale)),
            normal_power_w_m2: if mode.kind == UniaxialWaveKind::Propagating {
                s[2].max(0.0) * amplitude * amplitude * hscale
            } else {
                0.0
            },
            phase_direction: if mode.kind == UniaxialWaveKind::Evanescent {
                None
            } else {
                Some(unit([p, 0.0, mode.z.re])?)
            },
            energy_direction: if mode.kind == UniaxialWaveKind::Evanescent
                || coefficients[m] == C64::ZERO
            {
                None
            } else {
                Some(unit(s)?)
            },
        };
    }
    if !result.incident_normal_power_w_m2.is_finite()
        || result.incident_normal_power_w_m2 <= 0.0
        || result
            .incident_electric_v_m
            .iter()
            .chain(&result.incident_magnetic_a_m)
            .chain(&result.reflected_electric_v_m)
            .chain(&result.reflected_magnetic_a_m)
            .any(|z| !finite(*z))
        || result.transmitted.iter().any(|m| {
            !finite(m.normal_wavenumber_rad_m)
                || !m.normal_power_w_m2.is_finite()
                || m.electric_v_m
                    .iter()
                    .chain(&m.magnetic_a_m)
                    .any(|z| !finite(*z))
        })
    {
        return Err(invalid("oriented crystal fields are unrepresentable"));
    }
    if let Some(capture) = capture {
        *capture = Some(BoundaryWitness {
            basis,
            coefficients,
            amplitude,
            incident_electric: data.electric,
        });
    }
    Ok(result)
}

fn inverse_tensor(
    v: Vector,
    axis: [f64; 3],
    ordinary_squared: f64,
    extraordinary_squared: f64,
) -> Vector {
    let projection = v[0].scale(axis[0]) + v[1].scale(axis[1]) + v[2].scale(axis[2]);
    core::array::from_fn(|j| {
        v[j].scale(1.0 / ordinary_squared)
            + projection.scale(axis[j] * (1.0 / extraordinary_squared - 1.0 / ordinary_squared))
    })
}

/// Local bulk derivatives of the original native outgoing modes. No boundary
/// solve, source replacement or coherent-superposition velocity is introduced.
pub(super) fn mode_velocities(
    boundary: &OrientedCrystalResponse,
    omega_rad_s: f64,
    indices: [f64; 2],
    group_indices: [f64; 2],
    optic_axis: [f64; 3],
) -> Result<[Option<[f64; 3]>; 2]> {
    let k0 = positive(positive(omega_rad_s)? / C0)?;
    let [no, ne] = indices.map(positive);
    let (no, ne) = (no?, ne?);
    if group_indices.iter().any(|value| !value.is_finite()) {
        return Err(invalid("uniaxial group indices must be finite"));
    }
    let axis = unit(optic_axis)?;
    // At exact isotropy the native kernel deliberately chooses TE/TM. When
    // ay!=0 and the original slopes split, that basis mixes spectral branches.
    // Preserve the carrier but refuse an unsupported per-mode interpretation.
    if no == ne && group_indices[0] != group_indices[1] && axis[1] != 0.0 {
        return Err(invalid(
            "equal principal phase indices with unequal group slopes need a spectral eigenbasis for an out-of-plane optic axis",
        ));
    }
    let mut velocities = [None; 2];
    for (m, mode) in boundary.transmitted.iter().enumerate() {
        let Some(mut direction) = mode.energy_direction else {
            continue;
        };
        if mode.wave_kind == UniaxialWaveKind::Evanescent {
            continue;
        }
        if mode.wave_kind == UniaxialWaveKind::Critical {
            // The native field cross product may retain roundoff normal flux.
            // An exactly critical branch has only tangent group transport.
            direction[2] = 0.0;
            direction = unit(direction)?;
        }
        let speed = if m == 0 {
            if group_indices[0] == 0.0 {
                continue;
            }
            C0 / group_indices[0]
        } else {
            let qx = boundary.tangential_wavenumber_rad_m / k0;
            let qz = mode.normal_wavenumber_rad_m.re / k0;
            let along = (axis[0] * qx + axis[2] * qz) / no;
            // |q cross a| equals the perpendicular wavevector magnitude.
            // This preserves small transverse components without subtracting
            // two nearly equal squared norms.
            let across = [
                (-axis[1] * qz) / ne,
                (axis[0] * qz - axis[2] * qx) / ne,
                (axis[1] * qx) / ne,
            ];
            if ![qx, qz, along]
                .into_iter()
                .chain(across)
                .all(f64::is_finite)
            {
                return Err(invalid(
                    "uniaxial group-wave coordinates are unrepresentable",
                ));
            }
            let weights = [
                Dd::from_f64(along) * Dd::from_f64(along),
                across.into_iter().fold(Dd::ZERO, |sum, value| {
                    let value = Dd::from_f64(value);
                    sum + value * value
                }),
            ];
            let Some(denominator) = dispersive_group_denominator(weights, [no, ne], group_indices)?
            else {
                continue;
            };
            let perpendicular = det::hypot(det::hypot(across[0], across[1]), across[2]);
            let metric_magnitude = positive(det::hypot(along / no, perpendicular / ne))?;
            C0 * (metric_magnitude / denominator)
        };
        let velocity = direction.map(|value| value * speed);
        if !speed.is_finite()
            || speed == 0.0
            || velocity.iter().any(|value| !value.is_finite())
            || velocity == [0.0; 3]
        {
            return Err(invalid("uniaxial group velocity is unrepresentable"));
        }
        velocities[m] = Some(velocity);
    }
    Ok(velocities)
}

#[cfg(test)]
#[path = "group_tests.rs"]
mod group_tests;

#[cfg(test)]
#[path = "tests.rs"]
mod tests;
