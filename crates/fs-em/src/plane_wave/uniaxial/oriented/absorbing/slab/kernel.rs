//! Fixed-size matching of the original complex Maxwell modes at both faces.
use super::super::super::{
    kernel::{cross, unit},
    slab::kernel::{exit_basis, finite, flux, positive, solve_boundary, tangential},
};
use super::super::kernel::{Input, eigenbasis};
use super::{AbsorbingSlabResponse, AbsorbingSlabWave, C0, Vector};
use crate::plane_wave::{angular_spectrum::modes::propagation, uniaxial::kernel::UniaxialWaveKind};
use fs_math::{c64::C64, dd::Dd};
type Result<T> = core::result::Result<T, &'static str>;
type Wave = AbsorbingSlabWave;

fn combine(waves: [Wave; 2], c: [C64; 2]) -> (Vector, Vector) {
    (
        core::array::from_fn(|j| waves[0].electric_v_m[j] * c[0] + waves[1].electric_v_m[j] * c[1]),
        core::array::from_fn(|j| waves[0].magnetic_a_m[j] * c[0] + waves[1].magnetic_a_m[j] * c[1]),
    )
}

pub(in crate::plane_wave::uniaxial::oriented::absorbing) fn solve(
    data: Input,
    exit_index: f64,
    exit_mu: f64,
    depth: f64,
) -> Result<AbsorbingSlabResponse> {
    solve_from_both_sides(data, exit_index, exit_mu, depth, [C64::ZERO; 2])
}

pub(super) fn solve_from_both_sides(
    data: Input,
    exit_index: f64,
    exit_mu: f64,
    depth: f64,
    right_electric: [C64; 2],
) -> Result<AbsorbingSlabResponse> {
    positive(depth)?;
    positive(exit_index)?;
    positive(exit_mu)?;
    if data
        .electric
        .iter()
        .chain(&right_electric)
        .any(|&v| !finite(v))
        || data.loss_f_m.iter().any(|v| !v.is_finite() || *v < 0.0)
        || data
            .indices
            .iter()
            .zip(data.loss_f_m)
            .any(|(n, l)| (n.im == 0.0) != (l == 0.0))
    {
        return Err("absorbing film needs finite fields and passive original dielectric losses");
    }
    let right_is_driven = right_electric.iter().any(|&v| v != C64::ZERO);
    let amplitude = data
        .electric
        .iter()
        .chain(&right_electric)
        .fold(0.0_f64, |m, v| m.max(v.re.abs()).max(v.im.abs()));
    if amplitude == 0.0 {
        return Err("absorbing film requires at least one nonzero incident port field");
    }
    let amplitude = positive(amplitude)?;
    let electric = data
        .electric
        .map(|v| C64::new(v.re / amplitude, v.im / amplitude));
    let right_electric = right_electric.map(|v| C64::new(v.re / amplitude, v.im / amplitude));
    let basis = eigenbasis(data)?;
    // Reflect the tensor through z=0, solve its outgoing modes, and transform
    // polar E / axial H back. Negating the forward kz alone is wrong for a
    // tilted anisotropic tensor, including its complex extraordinary root.
    let reversed = eigenbasis(Input {
        axis: [data.axis[0], data.axis[1], -data.axis[2]],
        ..data
    })?;
    let kscale = positive((data.omega / C0) * basis.scale)?;
    let hscale = positive(basis.scale / (data.mut_ * C0))?;
    let mut forward = basis.modes.map(|m| Wave {
        normal_wavenumber_rad_m: m.z.scale(kscale),
        electric_v_m: m.e,
        magnetic_a_m: m.h.map(|v| v.scale(hscale)),
    });
    let mut backward = reversed.modes.map(|m| Wave {
        normal_wavenumber_rad_m: -m.z.scale(kscale),
        electric_v_m: [m.e[0], m.e[1], -m.e[2]],
        magnetic_a_m: [
            -m.h[0].scale(hscale),
            -m.h[1].scale(hscale),
            m.h[2].scale(hscale),
        ],
    });
    for (original, wave) in basis
        .modes
        .iter()
        .zip(&forward)
        .chain(reversed.modes.iter().zip(&backward))
    {
        let k = wave.normal_wavenumber_rad_m;
        if !finite(k) || original.z.re != 0.0 && k.re == 0.0 || original.z.im != 0.0 && k.im == 0.0
        {
            return Err("absorbing film wavenumber or attenuation is unrepresentable");
        }
    }
    let (p, qi) = (basis.p, basis.qi);
    let incident_hscale = positive(basis.scale / (data.mui * C0))?;
    let ie = [electric[0], electric[1], electric[0].scale(-p / qi)];
    let ih =
        cross([C64::from_re(p), C64::ZERO, C64::from_re(qi)], ie).map(|v| v.scale(incident_hscale));
    let incoming_left = flux(ie, ih);
    let port_p = p * basis.scale;
    let reflected = [
        Wave {
            normal_wavenumber_rad_m: C64::from_re(-kscale * qi),
            electric_v_m: [C64::ZERO, C64::ONE, C64::ZERO],
            magnetic_a_m: [
                C64::from_re(qi * incident_hscale),
                C64::ZERO,
                C64::from_re(p * incident_hscale),
            ],
        },
        Wave {
            normal_wavenumber_rad_m: C64::from_re(-kscale * qi),
            electric_v_m: [
                C64::from_re(qi / basis.i),
                C64::ZERO,
                C64::from_re(p / basis.i),
            ],
            magnetic_a_m: [
                C64::ZERO,
                C64::from_re(-basis.i * incident_hscale),
                C64::ZERO,
            ],
        },
    ];
    let (exit, exit_kind) = if exit_index == data.ni {
        // Preserve a positive original cosine when the conserved tangential
        // index has rounded to ni. Equal phase indices have equal normal k.
        let h = positive(basis.scale / (exit_mu * C0))?;
        let kz = C64::from_re(kscale * qi);
        (
            [
                Wave {
                    normal_wavenumber_rad_m: kz,
                    electric_v_m: [C64::ZERO, C64::ONE, C64::ZERO],
                    magnetic_a_m: [-C64::from_re(qi * h), C64::ZERO, C64::from_re(p * h)],
                },
                Wave {
                    normal_wavenumber_rad_m: kz,
                    electric_v_m: [
                        C64::from_re(qi / basis.i),
                        C64::ZERO,
                        C64::from_re(-p / basis.i),
                    ],
                    magnetic_a_m: [C64::ZERO, C64::from_re(basis.i * h), C64::ZERO],
                },
            ],
            UniaxialWaveKind::Propagating,
        )
    } else {
        let basis = exit_basis(data.omega, port_p, exit_index, exit_mu)?;
        let kind = basis[0].wave_kind;
        (
            basis.map(|w| Wave {
                normal_wavenumber_rad_m: w.normal_wavenumber_rad_m,
                electric_v_m: w.electric_v_m,
                magnetic_a_m: w.magnetic_a_m,
            }),
            kind,
        )
    };
    let (right_e, right_h) = if right_is_driven {
        if exit_kind != UniaxialWaveKind::Propagating {
            return Err("a driven exit port must be a propagating lossless channel");
        }
        // The exit basis is outgoing. Reflect its polar E and axial H to
        // obtain the incoming port with the same global Ex/Ey at z=depth.
        let coefficients = [
            right_electric[1],
            right_electric[0] / exit[1].electric_v_m[0],
        ];
        let (e, h) = combine(exit, coefficients);
        (
            [right_electric[0], right_electric[1], -e[2]],
            [-h[0], -h[1], h[2]],
        )
    } else {
        ([C64::ZERO; 3], [C64::ZERO; 3])
    };
    let incoming_right = -flux(right_e, right_h);
    if !incoming_left.is_finite()
        || incoming_left < 0.0
        || !incoming_right.is_finite()
        || incoming_right < 0.0
    {
        return Err("absorbing film incident port power is unrepresentable");
    }
    let incoming = positive((Dd::from_f64(incoming_left) + Dd::from_f64(incoming_right)).to_f64())?;
    let mut pf = [C64::ZERO; 2];
    let mut pb = [C64::ZERO; 2];
    for j in 0..2 {
        pf[j] = propagation(forward[j].normal_wavenumber_rad_m, depth)
            .ok_or("absorbing film forward phase is unrepresentable")?;
        pb[j] = propagation(-backward[j].normal_wavenumber_rad_m, depth)
            .ok_or("absorbing film backward phase is unrepresentable")?;
    }
    let mut matrix = [[C64::ZERO; 8]; 8];
    let mut rhs = [C64::ZERO; 8];
    let incident = tangential(ie, ih);
    let right_incident = tangential(right_e, right_h);
    for row in 0..4 {
        rhs[row] = -incident[row];
        rhs[row + 4] = right_incident[row];
        for m in 0..2 {
            let r = tangential(reflected[m].electric_v_m, reflected[m].magnetic_a_m);
            let a = tangential(forward[m].electric_v_m, forward[m].magnetic_a_m);
            let b = tangential(backward[m].electric_v_m, backward[m].magnetic_a_m);
            let t = tangential(exit[m].electric_v_m, exit[m].magnetic_a_m);
            matrix[row][m] = r[row];
            matrix[row][m + 2] = -a[row];
            matrix[row][m + 4] = -b[row] * pb[m];
            matrix[row + 4][m + 2] = a[row] * pf[m];
            matrix[row + 4][m + 4] = b[row];
            matrix[row + 4][m + 6] = -t[row];
        }
    }
    let c = solve_boundary(matrix, rhs)?;
    let (re, rh) = combine(reflected, [c[0], c[1]]);
    let (te, th) = combine(exit, [c[6], c[7]]);
    for m in 0..2 {
        forward[m].electric_v_m = forward[m].electric_v_m.map(|v| v * c[m + 2]);
        forward[m].magnetic_a_m = forward[m].magnetic_a_m.map(|v| v * c[m + 2]);
        backward[m].electric_v_m = backward[m].electric_v_m.map(|v| v * c[m + 4]);
        backward[m].magnetic_a_m = backward[m].magnetic_a_m.map(|v| v * c[m + 4]);
    }
    // Check the solved tangential Maxwell equations independently of the
    // elimination. Local energy density is deliberately not evaluated here.
    let magnitude = |v: Vector| v.iter().fold(0.0_f64, |s, z| s.max(z.abs()));
    for (factors_f, factors_b, exterior_e, exterior_h, exterior_scale) in [
        (
            [C64::ONE; 2],
            pb,
            core::array::from_fn::<_, 3, _>(|j| ie[j] + re[j]),
            core::array::from_fn::<_, 3, _>(|j| ih[j] + rh[j]),
            [
                magnitude(ie).max(magnitude(re)),
                magnitude(ih).max(magnitude(rh)),
            ],
        ),
        (
            pf,
            [C64::ONE; 2],
            core::array::from_fn::<_, 3, _>(|j| te[j] + right_e[j]),
            core::array::from_fn::<_, 3, _>(|j| th[j] + right_h[j]),
            [
                magnitude(te).max(magnitude(right_e)),
                magnitude(th).max(magnitude(right_h)),
            ],
        ),
    ] {
        let (ef, hf) = combine(forward, factors_f);
        let (eb, hb) = combine(backward, factors_b);
        for (a, b, c, exterior_bound) in [
            (exterior_e, ef, eb, exterior_scale[0]),
            (exterior_h, hf, hb, exterior_scale[1]),
        ] {
            let scale = a
                .iter()
                .chain(&b)
                .chain(&c)
                .fold(exterior_bound, |s, v| s.max(v.abs()));
            if !scale.is_finite()
                || (0..2).any(|j| (a[j] - b[j] - c[j]).abs() > 8192.0 * f64::EPSILON * scale)
            {
                return Err("absorbing film tangential Maxwell continuity is unresolved");
            }
        }
    }
    let r = -flux(re, rh) / incoming;
    let mut t = flux(te, th) / incoming;
    if exit_kind != UniaxialWaveKind::Propagating {
        if t.abs() > 8192.0 * f64::EPSILON {
            return Err("absorbing film nonpropagating exit carries normal power");
        }
        t = 0.0;
    }
    let a = (Dd::ONE - Dd::from_f64(r) - Dd::from_f64(t)).to_f64();
    let tolerance = 8192.0 * f64::EPSILON;
    if !r.is_finite()
        || !t.is_finite()
        || !a.is_finite()
        || r < 0.0
        || t < 0.0
        || a < -tolerance
        || data.loss_f_m == [0.0; 2] && a.abs() > tolerance
    {
        return Err("absorbing film violates passive normal-power balance");
    }
    let tangential = kscale * p;
    if !tangential.is_finite() || p != 0.0 && tangential == 0.0 {
        return Err("absorbing film tangential wavenumber is unrepresentable");
    }
    for wave in forward.iter_mut().chain(&mut backward) {
        wave.electric_v_m = wave.electric_v_m.map(|v| v.scale(amplitude));
        wave.magnetic_a_m = wave.magnetic_a_m.map(|v| v.scale(amplitude));
    }
    let response = AbsorbingSlabResponse {
        thickness_m: depth,
        tangential_wavenumber_rad_m: tangential,
        incident_normal_wavenumber_rad_m: positive(kscale * qi)?,
        exit_normal_wavenumber_rad_m: exit[0].normal_wavenumber_rad_m,
        incident_electric_v_m: ie.map(|v| v.scale(amplitude)),
        incident_magnetic_a_m: ih.map(|v| v.scale(amplitude)),
        reflected_electric_v_m: re.map(|v| v.scale(amplitude)),
        reflected_magnetic_a_m: rh.map(|v| v.scale(amplitude)),
        transmitted_electric_v_m: te.map(|v| v.scale(amplitude)),
        transmitted_magnetic_a_m: th.map(|v| v.scale(amplitude)),
        forward,
        backward,
        incident_normal_power_w_m2: (incoming_left * amplitude) * amplitude,
        right_incident_electric_v_m: right_e.map(|v| v.scale(amplitude)),
        right_incident_magnetic_a_m: right_h.map(|v| v.scale(amplitude)),
        right_incident_normal_power_w_m2: (incoming_right * amplitude) * amplitude,
        reflected_power_fraction: r,
        transmitted_power_fraction: t,
        absorbed_power_fraction: if data.loss_f_m == [0.0; 2] {
            0.0
        } else {
            a.max(0.0)
        },
        omega: data.omega,
        axis: unit(data.axis)?,
        loss_f_m: data.loss_f_m,
        plane: [1.0, 0.0],
    };
    positive(response.total_incident_normal_power_w_m2())?;
    response.validate_fields()?;
    Ok(response)
}

impl AbsorbingSlabResponse {
    fn validate_fields(&self) -> Result<()> {
        let ports = [
            self.incident_electric_v_m,
            self.incident_magnetic_a_m,
            self.reflected_electric_v_m,
            self.reflected_magnetic_a_m,
            self.transmitted_electric_v_m,
            self.transmitted_magnetic_a_m,
            self.right_incident_electric_v_m,
            self.right_incident_magnetic_a_m,
        ];
        if ports.iter().flatten().any(|&v| !finite(v))
            || self.forward.iter().chain(&self.backward).any(|w| {
                w.electric_v_m
                    .iter()
                    .chain(&w.magnetic_a_m)
                    .any(|&v| !finite(v))
            })
        {
            return Err("absorbing film dimensional fields are unrepresentable");
        }
        Ok(())
    }

    pub(super) fn rotate(&mut self, plane: [f64; 2], global_axis: [f64; 3]) -> Result<()> {
        let rotate = |v: Vector| {
            [
                v[0].scale(plane[0]) - v[1].scale(plane[1]),
                v[0].scale(plane[1]) + v[1].scale(plane[0]),
                v[2],
            ]
        };
        self.incident_electric_v_m = rotate(self.incident_electric_v_m);
        self.incident_magnetic_a_m = rotate(self.incident_magnetic_a_m);
        self.reflected_electric_v_m = rotate(self.reflected_electric_v_m);
        self.reflected_magnetic_a_m = rotate(self.reflected_magnetic_a_m);
        self.transmitted_electric_v_m = rotate(self.transmitted_electric_v_m);
        self.transmitted_magnetic_a_m = rotate(self.transmitted_magnetic_a_m);
        self.right_incident_electric_v_m = rotate(self.right_incident_electric_v_m);
        self.right_incident_magnetic_a_m = rotate(self.right_incident_magnetic_a_m);
        for wave in self.forward.iter_mut().chain(&mut self.backward) {
            wave.electric_v_m = rotate(wave.electric_v_m);
            wave.magnetic_a_m = rotate(wave.magnetic_a_m);
        }
        self.axis = global_axis;
        self.plane = plane;
        self.validate_fields()
    }
}

#[cfg(test)]
mod tests;
