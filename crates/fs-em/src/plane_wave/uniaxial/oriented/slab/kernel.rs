//! Fixed-size two-face Maxwell matching with bounded propagation factors.
use super::super::super::kernel::{UniaxialWaveKind, binary_scale, outgoing_root};
use super::super::kernel::{Eigenbasis, FieldSum, Input as ModeInput, WaveGeometry, eigenbasis};
use super::super::spectrum::{UniaxialSpectrumField, rotate};
use super::super::{Incidence, ResolvedFrequency};
use crate::plane_wave::angular_spectrum::{AngularSpectrumError, modes::propagation};
use fs_math::{c64::C64, dd::Dd};

const C0: f64 = 299_792_458.0;
type Vector = [C64; 3];
type Result<T> = core::result::Result<T, &'static str>;

/// Physical region for a field point, with entrance z=0 and exit z=d.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum SlabRegion {
    Incident,
    Reflected,
    TotalIncidentSide,
    Crystal,
    Transmitted,
}

/// One actual internal eigenmode. Forward fields are anchored at the selected
/// layer's entrance; backward fields at its exit. The wavevector is global:
/// a backward evanescent mode has Im(kz)>0 and decays as z decreases.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct SlabWave {
    pub wave_kind: UniaxialWaveKind,
    pub normal_wavenumber_rad_m: C64,
    pub electric_v_m: Vector,
    pub magnetic_a_m: Vector,
}

/// Source-evaluated coherent fields and independently computed normal power.
/// Transmitted E/H are at the actual exit z=d, without a vacuum phase removal.
#[derive(Clone, Copy, Debug, PartialEq)]
pub struct UniaxialSlabResponse {
    pub thickness_m: f64,
    pub tangential_wavenumber_rad_m: f64,
    pub incident_normal_wavenumber_rad_m: f64,
    pub exit_normal_wavenumber_rad_m: C64,
    pub exit_wave_kind: UniaxialWaveKind,
    pub incident_electric_v_m: Vector,
    pub incident_magnetic_a_m: Vector,
    pub reflected_electric_v_m: Vector,
    pub reflected_magnetic_a_m: Vector,
    pub transmitted_electric_v_m: Vector,
    pub transmitted_magnetic_a_m: Vector,
    /// Ordinary then extraordinary, with positive outgoing normal energy.
    pub forward: [SlabWave; 2],
    /// Ordinary then extraordinary, with negative outgoing normal energy.
    pub backward: [SlabWave; 2],
    pub incident_normal_power_w_m2: f64,
    pub reflected_power_fraction: f64,
    pub transmitted_power_fraction: f64,
}

impl UniaxialSlabResponse {
    /// Actual complex peak E/H and cycle-mean E cross H*/2 at [x,y,z] [m].
    /// Single-mode responses use the original x-z incidence plane; prepared
    /// spectra additionally rotate every mode into its global azimuth.
    pub fn field_at(
        &self,
        region: SlabRegion,
        position_m: [f64; 3],
    ) -> core::result::Result<UniaxialSpectrumField, AngularSpectrumError> {
        admit_point(region, position_m, self.thickness_m)?;
        let mut sum = FieldSum::at_point(position_m).map_err(AngularSpectrumError::Invalid)?;
        self.accumulate(&mut sum, region, [1.0, 0.0])?;
        let field = sum.finish().map_err(AngularSpectrumError::Invalid)?;
        Ok(UniaxialSpectrumField {
            electric_v_m: field.electric_v_m,
            magnetic_a_m: field.magnetic_a_m,
            poynting_w_m2: field.poynting_w_m2,
        })
    }

    pub(super) fn rotate(&mut self, axis: [f64; 2]) {
        self.incident_electric_v_m = rotate(self.incident_electric_v_m, axis);
        self.incident_magnetic_a_m = rotate(self.incident_magnetic_a_m, axis);
        self.reflected_electric_v_m = rotate(self.reflected_electric_v_m, axis);
        self.reflected_magnetic_a_m = rotate(self.reflected_magnetic_a_m, axis);
        self.transmitted_electric_v_m = rotate(self.transmitted_electric_v_m, axis);
        self.transmitted_magnetic_a_m = rotate(self.transmitted_magnetic_a_m, axis);
        for wave in self.forward.iter_mut().chain(&mut self.backward) {
            wave.electric_v_m = rotate(wave.electric_v_m, axis);
            wave.magnetic_a_m = rotate(wave.magnetic_a_m, axis);
        }
    }

    /// The single response or prepared spectrum admits the point once before
    /// adding waves. Every prepared mode belongs to the same slab thickness.
    pub(super) fn accumulate(
        &self,
        sum: &mut FieldSum,
        region: SlabRegion,
        axis: [f64; 2],
    ) -> core::result::Result<(), AngularSpectrumError> {
        let mut add = |e, h, k, anchor| {
            sum.add_wave(
                e,
                h,
                WaveGeometry {
                    tangential_wavenumber_rad_m: self.tangential_wavenumber_rad_m,
                    incidence_plane_axis: axis,
                    normal_wavenumber_rad_m: k,
                    anchor_z_m: anchor,
                },
            )
            .map_err(AngularSpectrumError::Invalid)
        };
        if matches!(region, SlabRegion::Incident | SlabRegion::TotalIncidentSide) {
            add(
                &self.incident_electric_v_m,
                &self.incident_magnetic_a_m,
                C64::from_re(self.incident_normal_wavenumber_rad_m),
                0.0,
            )?;
        }
        if matches!(
            region,
            SlabRegion::Reflected | SlabRegion::TotalIncidentSide
        ) {
            add(
                &self.reflected_electric_v_m,
                &self.reflected_magnetic_a_m,
                C64::from_re(-self.incident_normal_wavenumber_rad_m),
                0.0,
            )?;
        }
        if region == SlabRegion::Crystal {
            for wave in &self.forward {
                add(
                    &wave.electric_v_m,
                    &wave.magnetic_a_m,
                    wave.normal_wavenumber_rad_m,
                    0.0,
                )?;
            }
            for wave in &self.backward {
                // Keep the actual backward root and exit-face anchor. The
                // common kernel retains z-d before phase/decay multiplication;
                // this is exp(-i*(-k_minus)*(d-z)) without rounded distances.
                add(
                    &wave.electric_v_m,
                    &wave.magnetic_a_m,
                    wave.normal_wavenumber_rad_m,
                    self.thickness_m,
                )?;
            }
        }
        if region == SlabRegion::Transmitted {
            add(
                &self.transmitted_electric_v_m,
                &self.transmitted_magnetic_a_m,
                self.exit_normal_wavenumber_rad_m,
                self.thickness_m,
            )?;
        }
        Ok(())
    }
}

pub(super) fn admit_point(
    region: SlabRegion,
    point: [f64; 3],
    depth: f64,
) -> core::result::Result<(), AngularSpectrumError> {
    let supported = match region {
        SlabRegion::Crystal => (0.0..=depth).contains(&point[2]),
        SlabRegion::Transmitted => point[2] >= depth,
        _ => point[2] <= 0.0,
    };
    if point.iter().any(|v| !v.is_finite()) || !supported {
        return Err(AngularSpectrumError::Invalid(
            "uniaxial slab point is outside its selected region",
        ));
    }
    Ok(())
}

#[derive(Clone, Copy)]
pub(super) struct Input {
    pub(super) crystal: ResolvedFrequency,
    pub(super) exit_index: f64,
    pub(super) exit_mu_h_m: f64,
    pub(super) thickness_m: f64,
    pub(super) axis: [f64; 3],
}

pub(in super::super) fn positive(v: f64) -> Result<f64> {
    if v.is_finite() && v > 0.0 {
        Ok(v)
    } else {
        Err("uniaxial slab wave data must be positive and representable")
    }
}
pub(in super::super) fn finite(v: C64) -> bool {
    v.re.is_finite() && v.im.is_finite()
}
pub(in super::super) fn tangential(e: Vector, h: Vector) -> [C64; 4] {
    [e[0], e[1], h[0], h[1]]
}
pub(in super::super) fn flux(e: Vector, h: Vector) -> f64 {
    let product = |j: usize, k: usize| {
        Dd::from_f64(e[j].re) * Dd::from_f64(h[k].re)
            + Dd::from_f64(e[j].im) * Dd::from_f64(h[k].im)
    };
    ((product(0, 1) - product(1, 0)) * Dd::from_f64(0.5)).to_f64()
}
pub(in super::super) fn combine(waves: [SlabWave; 2], coefficients: [C64; 2]) -> (Vector, Vector) {
    let e = core::array::from_fn(|j| {
        waves[0].electric_v_m[j] * coefficients[0] + waves[1].electric_v_m[j] * coefficients[1]
    });
    let h = core::array::from_fn(|j| {
        waves[0].magnetic_a_m[j] * coefficients[0] + waves[1].magnetic_a_m[j] * coefficients[1]
    });
    (e, h)
}

pub(super) fn solve(
    data: Input,
    incidence: Incidence,
    electric: [C64; 2],
) -> Result<UniaxialSlabResponse> {
    let depth = positive(data.thickness_m)?;
    if electric.iter().any(|&v| !finite(v)) {
        return Err("finite incident slab fields required");
    }
    let amplitude = positive(
        electric
            .iter()
            .fold(0.0_f64, |m, z| m.max(z.re.abs()).max(z.im.abs())),
    )?;
    let electric = electric.map(|z| C64::new(z.re / amplitude, z.im / amplitude));
    let f = data.crystal;
    let original = ModeInput {
        omega: f.omega,
        ni: f.ni,
        no: f.no,
        ne: f.ne,
        mui: f.mui,
        mut_: f.mut_,
        axis: data.axis,
        incidence,
        electric,
    };
    let (basis, mut forward, mut backward) = crystal_basis(original)?;
    let kscale = positive((f.omega / C0) * basis.scale)?;
    let incident_hscale = positive(basis.scale / (f.mui * C0))?;
    let (p, qi, i) = (basis.p, basis.qi, basis.i);
    let incident_k = positive(kscale * qi)?;
    let tangential_k = kscale * p;
    if !tangential_k.is_finite() || p != 0.0 && tangential_k == 0.0 {
        return Err("uniaxial slab tangential wavenumber is unrepresentable");
    }
    let incident_e = [
        electric[0],
        electric[1],
        -electric[0].scale(p) / C64::from_re(qi),
    ];
    let incident_h = [
        -incident_e[1].scale(qi * incident_hscale),
        (incident_e[0].scale(qi) - incident_e[2].scale(p)).scale(incident_hscale),
        incident_e[1].scale(p * incident_hscale),
    ];
    let incoming = positive(flux(incident_e, incident_h))?;
    let reflected = [
        SlabWave {
            wave_kind: UniaxialWaveKind::Propagating,
            normal_wavenumber_rad_m: C64::from_re(-kscale * qi),
            electric_v_m: [C64::ZERO, C64::ONE, C64::ZERO],
            magnetic_a_m: [
                C64::from_re(qi * incident_hscale),
                C64::ZERO,
                C64::from_re(p * incident_hscale),
            ],
        },
        SlabWave {
            wave_kind: UniaxialWaveKind::Propagating,
            normal_wavenumber_rad_m: C64::from_re(-kscale * qi),
            electric_v_m: [C64::from_re(qi / i), C64::ZERO, C64::from_re(p / i)],
            magnetic_a_m: [C64::ZERO, C64::from_re(-i * incident_hscale), C64::ZERO],
        },
    ];
    let exit = exit_basis(
        f.omega,
        p * basis.scale,
        positive(data.exit_index)?,
        positive(data.exit_mu_h_m)?,
    )?;
    let mut phase_forward = [C64::ZERO; 2];
    let mut phase_backward = [C64::ZERO; 2];
    for j in 0..2 {
        phase_forward[j] = propagation(forward[j].normal_wavenumber_rad_m, depth)
            .ok_or("forward slab phase is unrepresentable")?;
        phase_backward[j] = propagation(-backward[j].normal_wavenumber_rad_m, depth)
            .ok_or("backward slab phase is unrepresentable")?;
    }
    let mut matrix = [[C64::ZERO; 8]; 8];
    let mut rhs = [C64::ZERO; 8];
    let incident = tangential(incident_e, incident_h);
    for row in 0..4 {
        rhs[row] = -incident[row];
        for mode in 0..2 {
            let r = tangential(reflected[mode].electric_v_m, reflected[mode].magnetic_a_m);
            let a = tangential(forward[mode].electric_v_m, forward[mode].magnetic_a_m);
            let b = tangential(backward[mode].electric_v_m, backward[mode].magnetic_a_m);
            let t = tangential(exit[mode].electric_v_m, exit[mode].magnetic_a_m);
            matrix[row][mode] = r[row];
            matrix[row][mode + 2] = -a[row];
            matrix[row][mode + 4] = -b[row] * phase_backward[mode];
            matrix[row + 4][mode + 2] = a[row] * phase_forward[mode];
            matrix[row + 4][mode + 4] = b[row];
            matrix[row + 4][mode + 6] = -t[row];
        }
    }
    let coefficients = solve_boundary(matrix, rhs)?;
    let (re, rh) = combine(reflected, [coefficients[0], coefficients[1]]);
    let (te, th) = combine(exit, [coefficients[6], coefficients[7]]);
    for m in 0..2 {
        forward[m].electric_v_m = forward[m].electric_v_m.map(|v| v * coefficients[m + 2]);
        forward[m].magnetic_a_m = forward[m].magnetic_a_m.map(|v| v * coefficients[m + 2]);
        backward[m].electric_v_m = backward[m].electric_v_m.map(|v| v * coefficients[m + 4]);
        backward[m].magnetic_a_m = backward[m].magnetic_a_m.map(|v| v * coefficients[m + 4]);
    }
    let reflected_power = -flux(re, rh) / incoming;
    let transmitted_flux = flux(te, th) / incoming;
    let transmitted_power = if exit[0].wave_kind == UniaxialWaveKind::Propagating {
        positive(transmitted_flux)?
    } else {
        if transmitted_flux.abs() > 4096.0 * f64::EPSILON {
            return Err("nonpropagating slab exit has unresolved normal power");
        }
        0.0
    };
    if !reflected_power.is_finite()
        || reflected_power < 0.0
        || (reflected_power + transmitted_power - 1.0).abs() > 4096.0 * f64::EPSILON
    {
        return Err("uniaxial slab Maxwell normal-power balance failed");
    }
    let mut result = UniaxialSlabResponse {
        thickness_m: depth,
        tangential_wavenumber_rad_m: tangential_k,
        incident_normal_wavenumber_rad_m: incident_k,
        exit_normal_wavenumber_rad_m: exit[0].normal_wavenumber_rad_m,
        exit_wave_kind: exit[0].wave_kind,
        incident_electric_v_m: incident_e,
        incident_magnetic_a_m: incident_h,
        reflected_electric_v_m: re,
        reflected_magnetic_a_m: rh,
        transmitted_electric_v_m: te,
        transmitted_magnetic_a_m: th,
        forward,
        backward,
        incident_normal_power_w_m2: incoming,
        reflected_power_fraction: reflected_power,
        transmitted_power_fraction: transmitted_power,
    };
    boundary_gate(&result, positive(i * incident_hscale)?)?;
    result.incident_electric_v_m = result.incident_electric_v_m.map(|v| v.scale(amplitude));
    result.incident_magnetic_a_m = result.incident_magnetic_a_m.map(|v| v.scale(amplitude));
    result.reflected_electric_v_m = result.reflected_electric_v_m.map(|v| v.scale(amplitude));
    result.reflected_magnetic_a_m = result.reflected_magnetic_a_m.map(|v| v.scale(amplitude));
    result.transmitted_electric_v_m = result.transmitted_electric_v_m.map(|v| v.scale(amplitude));
    result.transmitted_magnetic_a_m = result.transmitted_magnetic_a_m.map(|v| v.scale(amplitude));
    for wave in result.forward.iter_mut().chain(&mut result.backward) {
        wave.electric_v_m = wave.electric_v_m.map(|v| v.scale(amplitude));
        wave.magnetic_a_m = wave.magnetic_a_m.map(|v| v.scale(amplitude));
    }
    result.incident_normal_power_w_m2 = positive((incoming * amplitude) * amplitude)?;
    let ports = [
        result.incident_electric_v_m,
        result.incident_magnetic_a_m,
        result.reflected_electric_v_m,
        result.reflected_magnetic_a_m,
        result.transmitted_electric_v_m,
        result.transmitted_magnetic_a_m,
    ];
    if ports.iter().flatten().any(|&v| !finite(v))
        || result.forward.iter().chain(&result.backward).any(|w| {
            !finite(w.normal_wavenumber_rad_m)
                || w.electric_v_m
                    .iter()
                    .chain(&w.magnetic_a_m)
                    .any(|&v| !finite(v))
        })
    {
        return Err("uniaxial slab dimensional fields are unrepresentable");
    }
    Ok(result)
}

/// Original dimensional forward/backward basis, shared by slabs and stacks.
pub(in super::super) fn crystal_basis(
    original: ModeInput,
) -> Result<(Eigenbasis, [SlabWave; 2], [SlabWave; 2])> {
    let basis = eigenbasis(original)?;
    // Under z'=-z, E is polar, H axial: E=(Ex',Ey',-Ez'),
    // H=(-Hx',-Hy',Hz'), k=(kx',ky',-kz'). The optic axis transforms as E.
    // This obtains the OTHER extraordinary root even for a tilted tensor;
    // simply negating the forward kz would be physically incorrect.
    let reversed = eigenbasis(ModeInput {
        axis: [original.axis[0], original.axis[1], -original.axis[2]],
        ..original
    })?;
    if basis
        .modes
        .iter()
        .chain(&reversed.modes)
        .any(|m| m.kind == UniaxialWaveKind::Critical)
    {
        return Err("coalesced critical crystal modes require a generalized slab basis");
    }
    let kscale = positive((original.omega / C0) * basis.scale)?;
    let hscale = positive(basis.scale / (original.mut_ * C0))?;
    let forward = basis.modes.map(|m| SlabWave {
        wave_kind: m.kind,
        normal_wavenumber_rad_m: m.z.scale(kscale),
        electric_v_m: m.e,
        magnetic_a_m: m.h.map(|v| v.scale(hscale)),
    });
    let backward = reversed.modes.map(|m| SlabWave {
        wave_kind: m.kind,
        normal_wavenumber_rad_m: -m.z.scale(kscale),
        electric_v_m: [m.e[0], m.e[1], -m.e[2]],
        magnetic_a_m: [
            -m.h[0].scale(hscale),
            -m.h[1].scale(hscale),
            m.h[2].scale(hscale),
        ],
    });
    for (basis_mode, wave) in basis
        .modes
        .iter()
        .zip(&forward)
        .chain(reversed.modes.iter().zip(&backward))
    {
        let k = wave.normal_wavenumber_rad_m;
        if !finite(k)
            || basis_mode.z.re != 0.0 && k.re == 0.0
            || basis_mode.z.im != 0.0 && k.im == 0.0
        {
            return Err("uniaxial slab internal wavenumber is unrepresentable");
        }
    }
    Ok((basis, forward, backward))
}

pub(in super::super) fn exit_basis(omega: f64, p: f64, index: f64, mu: f64) -> Result<[SlabWave; 2]> {
    let scale = binary_scale(index.max(p.abs()));
    let n = positive(index / scale)?;
    if p != 0.0 && p / scale == 0.0 {
        return Err("uniaxial slab exit tangential index is unrepresentable");
    }
    let pd = Dd::from_f64(p / scale);
    let nd = Dd::from_f64(n);
    let (q, kind) = outgoing_root(nd * nd - pd * pd)?;
    let h = positive(scale / (mu * C0))?;
    let kz = q.scale(positive((omega / C0) * scale)?);
    if !finite(kz) || q.re != 0.0 && kz.re == 0.0 || q.im != 0.0 && kz.im == 0.0 {
        return Err("uniaxial slab exit wavenumber is unrepresentable");
    }
    Ok([
        SlabWave {
            wave_kind: kind,
            normal_wavenumber_rad_m: kz,
            electric_v_m: [C64::ZERO, C64::ONE, C64::ZERO],
            magnetic_a_m: [-q.scale(h), C64::ZERO, C64::from_re((p / scale) * h)],
        },
        SlabWave {
            wave_kind: kind,
            normal_wavenumber_rad_m: kz,
            electric_v_m: [q.scale(1.0 / n), C64::ZERO, C64::from_re(-(p / scale) / n)],
            magnetic_a_m: [C64::ZERO, C64::from_re(n * h), C64::ZERO],
        },
    ])
}

/// Row scaling and partial pivoting; an unresolved pivot is a refusal, never
/// a fabricated zero mode or a regularized material. All work is fixed-size.
pub(in super::super) fn solve_boundary<const N: usize>(
    mut a: [[C64; N]; N], mut b: [C64; N],
) -> Result<[C64; N]> {
    for row in 0..N {
        let scale = positive(a[row].iter().fold(0.0_f64, |v, z| v.max(z.abs())))?;
        if a[row].iter().any(|&z| !finite(z)) || !finite(b[row]) {
            return Err("finite slab boundary data required");
        }
        a[row] = a[row].map(|v| C64::new(v.re / scale, v.im / scale));
        b[row] = C64::new(b[row].re / scale, b[row].im / scale);
    }
    for column in 0..N {
        let mut pivot = column;
        for row in column + 1..N {
            if a[row][column].abs() > a[pivot][column].abs() {
                pivot = row;
            }
        }
        if !finite(a[pivot][column]) || a[pivot][column].abs() <= 64.0 * f64::EPSILON {
            return Err("uniaxial slab boundary system is numerically unresolved");
        }
        a.swap(pivot, column);
        b.swap(pivot, column);
        for row in column + 1..N {
            let factor = a[row][column] / a[column][column];
            a[row][column] = C64::ZERO;
            for j in column + 1..N {
                a[row][j] = a[row][j] - factor * a[column][j];
            }
            b[row] = b[row] - factor * b[column];
        }
    }
    let mut x = [C64::ZERO; N];
    for row in (0..N).rev() {
        let mut rhs = b[row];
        for j in row + 1..N {
            rhs = rhs - a[row][j] * x[j];
        }
        x[row] = rhs / a[row][row];
        if !finite(x[row]) {
            return Err("uniaxial slab amplitudes are unrepresentable");
        }
    }
    Ok(x)
}

fn boundary_gate(r: &UniaxialSlabResponse, admittance: f64) -> Result<()> {
    let get = |region, z| {
        r.field_at(region, [0.0, 0.0, z])
            .map_err(|_| "uniaxial slab boundary fields are unrepresentable")
    };
    let left = get(SlabRegion::TotalIncidentSide, 0.0)?;
    let front = get(SlabRegion::Crystal, 0.0)?;
    let back = get(SlabRegion::Crystal, r.thickness_m)?;
    let right = get(SlabRegion::Transmitted, r.thickness_m)?;
    // Scale by coherent boundary fields, never by large cancelling waves.
    // The unit-port floor admits ordinary floating-point zeros, while a
    // near-degenerate basis that loses finite transmitted fields refuses.
    let e_scale = [&left, &front, &back, &right].iter().fold(1.0_f64, |m, f| {
        m.max(f.electric_v_m[0].abs()).max(f.electric_v_m[1].abs())
    });
    let h_scale = [&left, &front, &back, &right]
        .iter()
        .fold(admittance, |m, f| {
            m.max(f.magnetic_a_m[0].abs()).max(f.magnetic_a_m[1].abs())
        });
    for (a, b) in [(left, front), (back, right)] {
        for j in 0..2 {
            if (a.electric_v_m[j] - b.electric_v_m[j]).abs() > 4096.0 * f64::EPSILON * e_scale
                || (a.magnetic_a_m[j] - b.magnetic_a_m[j]).abs() > 4096.0 * f64::EPSILON * h_scale
            {
                return Err("uniaxial slab tangential Maxwell continuity is unresolved");
            }
        }
    }
    Ok(())
}

#[cfg(test)]
#[path = "kernel_tests.rs"]
mod tests;
