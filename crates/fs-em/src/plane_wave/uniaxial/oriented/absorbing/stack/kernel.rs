//! Passive complex Maxwell modes in the existing O(layers) Schur sweep.
//!
//! At each right face, b = R a and a_next = T a. At its left face the
//! loaded tangential basis is F + B Q R P: P propagates forward to the right,
//! Q propagates backward to the left. Neither decaying factor is inverted.
//! The order matters for independently oriented anisotropic layers.

use super::super::super::{
    Incidence, OrientedCrystalField,
    kernel::FieldSum,
    slab::kernel::{exit_basis_from_incidence, finite, flux, positive, solve_boundary, tangential},
};
use super::super::{
    AbsorbingCrystalField,
    kernel::{Eigenbasis, Input, eigenbasis},
    slab::AbsorbingSlabWave,
};
use super::{
    AbsorbingStackLayerField, AbsorbingStackResponse, Frequency, ResolvedLayer, Result,
    StackRegion, bad, poll,
};
use crate::plane_wave::{angular_spectrum::modes::propagation, uniaxial::kernel::UniaxialWaveKind};
use fs_math::{c64::C64, dd::Dd};

const C0: f64 = 299_792_458.0;
type Vector = [C64; 3];
type Wave = AbsorbingSlabWave;
type Matrix = [[C64; 2]; 2];
type Trace = [[C64; 4]; 2]; // [mode][Ex,Ey,Hx,Hy]

struct LayerWork {
    field: AbsorbingStackLayerField,
    forward_phase: [C64; 2],
    backward_phase: [C64; 2],
    reflection: Matrix,
    through: Matrix,
}

/// Retained responses plus one frequency-resolution/elimination workspace.
pub(super) fn workspace_bytes(layers: usize, modes: usize) -> Result<usize> {
    let retained = layers
        .checked_mul(core::mem::size_of::<AbsorbingStackLayerField>())
        .and_then(|n| n.checked_add(core::mem::size_of::<AbsorbingStackResponse>()))
        .and_then(|n| n.checked_mul(modes));
    let scratch = layers
        .checked_mul(core::mem::size_of::<LayerWork>() + core::mem::size_of::<ResolvedLayer>());
    retained
        .and_then(|n| n.checked_add(scratch?))
        .and_then(|n| n.checked_add(core::mem::size_of::<Frequency>()))
        .and_then(|n| {
            n.checked_add(
                core::mem::size_of::<FieldSum>() + core::mem::size_of::<AbsorbingCrystalField>(),
            )
        })
        .ok_or_else(|| bad("absorbing stack workspace size overflow"))
}

fn apply(matrix: Matrix, vector: [C64; 2]) -> [C64; 2] {
    matrix.map(|row| row[0] * vector[0] + row[1] * vector[1])
}

fn trace(waves: [Wave; 2]) -> Trace {
    waves.map(|wave| tangential(wave.electric_v_m, wave.magnetic_a_m))
}

fn combine(waves: [Wave; 2], coefficients: [C64; 2]) -> (Vector, Vector) {
    (
        core::array::from_fn(|j| {
            waves[0].electric_v_m[j] * coefficients[0] + waves[1].electric_v_m[j] * coefficients[1]
        }),
        core::array::from_fn(|j| {
            waves[0].magnetic_a_m[j] * coefficients[0] + waves[1].magnetic_a_m[j] * coefficients[1]
        }),
    )
}

fn interface(backward: [Wave; 2], next: Trace, forward: [Wave; 2]) -> Result<(Matrix, Matrix)> {
    let backward = trace(backward);
    let matrix = core::array::from_fn(|row| {
        [
            backward[0][row],
            backward[1][row],
            -next[0][row],
            -next[1][row],
        ]
    });
    let mut reflection = [[C64::ZERO; 2]; 2];
    let mut through = reflection;
    for (column, rhs) in trace(forward).into_iter().enumerate() {
        let solution = solve_boundary(matrix, rhs.map(|v| -v)).map_err(bad)?;
        for row in 0..2 {
            reflection[row][column] = solution[row];
            through[row][column] = solution[row + 2];
        }
    }
    Ok((reflection, through))
}

fn front(work: &LayerWork) -> Trace {
    let forward = trace(work.field.forward);
    let backward = trace(work.field.backward);
    core::array::from_fn(|column| {
        core::array::from_fn(|row| {
            forward[column][row]
                + backward[0][row]
                    * work.backward_phase[0]
                    * work.reflection[0][column]
                    * work.forward_phase[column]
                + backward[1][row]
                    * work.backward_phase[1]
                    * work.reflection[1][column]
                    * work.forward_phase[column]
        })
    })
}

fn weight(waves: [Wave; 2], coefficients: [C64; 2]) -> [Wave; 2] {
    core::array::from_fn(|m| Wave {
        electric_v_m: waves[m].electric_v_m.map(|v| v * coefficients[m]),
        magnetic_a_m: waves[m].magnetic_a_m.map(|v| v * coefficients[m]),
        ..waves[m]
    })
}

fn crystal_basis(data: Input) -> Result<(Eigenbasis, [Wave; 2], [Wave; 2])> {
    let basis = eigenbasis(data).map_err(bad)?;
    // Reflection transforms the constitutive tensor as well as the waves.
    // A tilted extraordinary backward root is not minus its forward root.
    let reversed = eigenbasis(Input {
        axis: [data.axis[0], data.axis[1], -data.axis[2]],
        ..data
    })
    .map_err(bad)?;
    let kscale = positive((data.omega / C0) * basis.scale).map_err(bad)?;
    let hscale = positive(basis.scale / (data.mut_ * C0)).map_err(bad)?;
    let forward = basis.modes.map(|mode| Wave {
        normal_wavenumber_rad_m: mode.z.scale(kscale),
        electric_v_m: mode.e,
        magnetic_a_m: mode.h.map(|v| v.scale(hscale)),
    });
    let backward = reversed.modes.map(|mode| Wave {
        normal_wavenumber_rad_m: -mode.z.scale(kscale),
        electric_v_m: [mode.e[0], mode.e[1], -mode.e[2]],
        magnetic_a_m: [
            -mode.h[0].scale(hscale),
            -mode.h[1].scale(hscale),
            mode.h[2].scale(hscale),
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
            return Err(bad(
                "absorbing stack wavenumber or attenuation is unrepresentable",
            ));
        }
    }
    if forward
        .iter()
        .any(|wave| wave.normal_wavenumber_rad_m.im > 0.0)
        || backward
            .iter()
            .any(|wave| wave.normal_wavenumber_rad_m.im < 0.0)
    {
        return Err(bad(
            "absorbing stack requires decaying-side passive eigenwaves",
        ));
    }
    Ok((basis, forward, backward))
}

pub(super) fn solve(
    data: &Frequency,
    incidence: Incidence,
    electric: [C64; 2],
    plane: [f64; 2],
    cancelled: &mut impl FnMut() -> bool,
) -> Result<AbsorbingStackResponse> {
    poll(cancelled)?;
    if data.layers.is_empty() || electric.iter().any(|&v| !finite(v)) {
        return Err(bad(
            "finite nonzero absorbing stack illumination and layers required",
        ));
    }
    positive(data.nt).map_err(bad)?;
    positive(data.mut_).map_err(bad)?;
    let amplitude = positive(
        electric
            .iter()
            .fold(0.0_f64, |m, v| m.max(v.re.abs()).max(v.im.abs())),
    )
    .map_err(bad)?;
    let electric = electric.map(|v| C64::new(v.re / amplitude, v.im / amplitude));
    let electric = [
        electric[0].scale(plane[0]) + electric[1].scale(plane[1]),
        electric[1].scale(plane[0]) - electric[0].scale(plane[1]),
    ];
    let mut work = Vec::new();
    work.try_reserve_exact(data.layers.len())
        .map_err(|_| bad("absorbing stack workspace allocation failed"))?;
    let mut first = None;
    for layer in &data.layers {
        poll(cancelled)?;
        positive(layer.thickness).map_err(bad)?;
        if layer.loss_f_m.iter().any(|v| !v.is_finite() || *v < 0.0)
            || layer
                .indices
                .iter()
                .zip(layer.loss_f_m)
                .any(|(n, loss)| (n.im == 0.0) != (loss == 0.0))
        {
            return Err(bad(
                "absorbing stack requires passive original dielectric losses",
            ));
        }
        let [x, y, z] = layer.axis;
        let (basis, forward, backward) = crystal_basis(Input {
            omega: data.omega,
            ni: data.ni,
            indices: layer.indices,
            mui: data.mui,
            mut_: layer.mu,
            axis: [x * plane[0] + y * plane[1], y * plane[0] - x * plane[1], z],
            // Preserve the authored cosine in EVERY layer. Converting it to a
            // rounded transverse index can erase a positive grazing normal k.
            incidence,
            electric,
            loss_f_m: layer.loss_f_m,
        })?;
        if first.is_none() {
            first = Some(basis);
        }
        let mut phases = [[C64::ZERO; 2]; 2];
        for m in 0..2 {
            phases[0][m] = propagation(forward[m].normal_wavenumber_rad_m, layer.thickness)
                .ok_or_else(|| bad("absorbing stack forward phase is unrepresentable"))?;
            phases[1][m] = propagation(-backward[m].normal_wavenumber_rad_m, layer.thickness)
                .ok_or_else(|| bad("absorbing stack backward phase is unrepresentable"))?;
        }
        work.push(LayerWork {
            field: AbsorbingStackLayerField {
                entrance_z_m: layer.entrance,
                exit_z_m: layer.exit,
                thickness_m: layer.thickness,
                forward,
                backward,
                omega: data.omega,
                axis: layer.axis,
                loss_f_m: layer.loss_f_m,
            },
            forward_phase: phases[0],
            backward_phase: phases[1],
            reflection: [[C64::ZERO; 2]; 2],
            through: [[C64::ZERO; 2]; 2],
        });
    }
    let basis = first.ok_or_else(|| bad("empty absorbing stack"))?;
    let kscale = positive((data.omega / C0) * basis.scale).map_err(bad)?;
    let hscale = positive(basis.scale / (data.mui * C0)).map_err(bad)?;
    let (p, q, n) = (basis.p, basis.qi, basis.i);
    let normal = positive(kscale * q).map_err(bad)?;
    let tangential_k = kscale * p;
    if !tangential_k.is_finite() || p != 0.0 && tangential_k == 0.0 {
        return Err(bad(
            "absorbing stack tangential wavenumber is unrepresentable",
        ));
    }
    let incident_e = [
        electric[0],
        electric[1],
        -electric[0].scale(p) / C64::from_re(q),
    ];
    let incident_h = [
        -incident_e[1].scale(q * hscale),
        (incident_e[0].scale(q) - incident_e[2].scale(p)).scale(hscale),
        incident_e[1].scale(p * hscale),
    ];
    let incoming = positive(flux(incident_e, incident_h)).map_err(bad)?;
    let reflected = [
        Wave {
            normal_wavenumber_rad_m: C64::from_re(-normal),
            electric_v_m: [C64::ZERO, C64::ONE, C64::ZERO],
            magnetic_a_m: [
                C64::from_re(q * hscale),
                C64::ZERO,
                C64::from_re(p * hscale),
            ],
        },
        Wave {
            normal_wavenumber_rad_m: C64::from_re(-normal),
            electric_v_m: [C64::from_re(q / n), C64::ZERO, C64::from_re(p / n)],
            magnetic_a_m: [C64::ZERO, C64::from_re(-n * hscale), C64::ZERO],
        },
    ];
    let (exit, exit_kind) = if data.nt == data.ni {
        // Equal-index ports share the preserved normal coordinate even when
        // n*sin(theta) has rounded to n. Permeabilities may still differ.
        let h = positive(basis.scale / (data.mut_ * C0)).map_err(bad)?;
        (
            [
                Wave {
                    normal_wavenumber_rad_m: C64::from_re(normal),
                    electric_v_m: [C64::ZERO, C64::ONE, C64::ZERO],
                    magnetic_a_m: [-C64::from_re(q * h), C64::ZERO, C64::from_re(p * h)],
                },
                Wave {
                    normal_wavenumber_rad_m: C64::from_re(normal),
                    electric_v_m: [C64::from_re(q / n), C64::ZERO, C64::from_re(-p / n)],
                    magnetic_a_m: [C64::ZERO, C64::from_re(n * h), C64::ZERO],
                },
            ],
            UniaxialWaveKind::Propagating,
        )
    } else {
        let exit =
            exit_basis_from_incidence(data.omega, data.ni, incidence, data.nt, data.mut_).map_err(bad)?;
        let kind = exit[0].wave_kind;
        (
            exit.map(|wave| Wave {
                normal_wavenumber_rad_m: wave.normal_wavenumber_rad_m,
                electric_v_m: wave.electric_v_m,
                magnetic_a_m: wave.magnetic_a_m,
            }),
            kind,
        )
    };
    let mut next = trace(exit);
    for layer in work.iter_mut().rev() {
        poll(cancelled)?;
        (layer.reflection, layer.through) =
            interface(layer.field.backward, next, layer.field.forward)?;
        next = front(layer);
    }
    let reflected_trace = trace(reflected);
    let entrance = core::array::from_fn(|row| {
        [
            reflected_trace[0][row],
            reflected_trace[1][row],
            -next[0][row],
            -next[1][row],
        ]
    });
    let coefficients =
        solve_boundary(entrance, tangential(incident_e, incident_h).map(|v| -v)).map_err(bad)?;
    let (re, rh) = combine(reflected, [coefficients[0], coefficients[1]]);
    let mut forward = [coefficients[2], coefficients[3]];
    let mut layers = Vec::new();
    layers
        .try_reserve_exact(work.len())
        .map_err(|_| bad("absorbing stack field allocation failed"))?;
    for layer in work {
        poll(cancelled)?;
        let arriving = core::array::from_fn(|m| forward[m] * layer.forward_phase[m]);
        let backward = apply(layer.reflection, arriving);
        layers.push(AbsorbingStackLayerField {
            forward: weight(layer.field.forward, forward),
            backward: weight(layer.field.backward, backward),
            ..layer.field
        });
        forward = apply(layer.through, arriving);
    }
    let (te, th) = combine(exit, forward);
    let reflected_power = -flux(re, rh) / incoming;
    let mut transmitted_power = flux(te, th) / incoming;
    let tolerance = 8192.0 * f64::EPSILON;
    if exit_kind != UniaxialWaveKind::Propagating {
        if !transmitted_power.is_finite() || transmitted_power.abs() > tolerance {
            return Err(bad(
                "absorbing stack nonpropagating exit carries unresolved normal power",
            ));
        }
        transmitted_power = 0.0;
    }
    let absorbed =
        (Dd::ONE - Dd::from_f64(reflected_power) - Dd::from_f64(transmitted_power)).to_f64();
    let lossless = data.layers.iter().all(|layer| layer.loss_f_m == [0.0; 2]);
    if !reflected_power.is_finite()
        || !transmitted_power.is_finite()
        || !absorbed.is_finite()
        || reflected_power < 0.0
        || transmitted_power < 0.0
        || absorbed < -tolerance
        || lossless && absorbed.abs() > tolerance
    {
        return Err(bad(
            "absorbing stack violates passive exterior normal-power balance",
        ));
    }
    let mut result = AbsorbingStackResponse {
        thickness_m: data
            .layers
            .last()
            .ok_or_else(|| bad("empty absorbing stack"))?
            .exit,
        tangential_wavenumber_rad_m: tangential_k,
        incident_normal_wavenumber_rad_m: normal,
        exit_normal_wavenumber_rad_m: exit[0].normal_wavenumber_rad_m,
        exit_wave_kind: exit_kind,
        incident_electric_v_m: incident_e,
        incident_magnetic_a_m: incident_h,
        reflected_electric_v_m: re,
        reflected_magnetic_a_m: rh,
        transmitted_electric_v_m: te,
        transmitted_magnetic_a_m: th,
        layers,
        incident_normal_power_w_m2: incoming,
        reflected_power_fraction: reflected_power,
        transmitted_power_fraction: transmitted_power,
        absorbed_power_fraction: if lossless { 0.0 } else { absorbed.max(0.0) },
        omega: data.omega,
        plane: [1.0, 0.0],
    };
    boundary_gate(&result, positive(n * hscale).map_err(bad)?, cancelled)?;
    for vector in [
        &mut result.incident_electric_v_m,
        &mut result.incident_magnetic_a_m,
        &mut result.reflected_electric_v_m,
        &mut result.reflected_magnetic_a_m,
        &mut result.transmitted_electric_v_m,
        &mut result.transmitted_magnetic_a_m,
    ] {
        *vector = vector.map(|v| v.scale(amplitude));
        if vector.iter().any(|&v| !finite(v)) {
            return Err(bad(
                "absorbing stack dimensional fields are unrepresentable",
            ));
        }
    }
    for layer in &mut result.layers {
        poll(cancelled)?;
        for wave in layer.forward.iter_mut().chain(&mut layer.backward) {
            wave.electric_v_m = wave.electric_v_m.map(|v| v.scale(amplitude));
            wave.magnetic_a_m = wave.magnetic_a_m.map(|v| v.scale(amplitude));
            if wave
                .electric_v_m
                .iter()
                .chain(&wave.magnetic_a_m)
                .any(|&v| !finite(v))
            {
                return Err(bad(
                    "absorbing stack dimensional fields are unrepresentable",
                ));
            }
        }
    }
    result.incident_normal_power_w_m2 =
        positive((incoming * amplitude) * amplitude).map_err(bad)?;
    result.rotate(plane)?;
    poll(cancelled)?;
    Ok(result)
}

fn boundary_gate(
    response: &AbsorbingStackResponse,
    admittance: f64,
    cancelled: &mut impl FnMut() -> bool,
) -> Result<()> {
    let tolerance = 8192.0 * f64::EPSILON;
    let mut left = response.plain_field_at(StackRegion::TotalIncidentSide, [0.0; 3])?;
    let incident = tangential(
        response.incident_electric_v_m,
        response.incident_magnetic_a_m,
    );
    let reflected = tangential(
        response.reflected_electric_v_m,
        response.reflected_magnetic_a_m,
    );
    let mut left_bound =
        core::array::from_fn::<_, 4, _>(|j| incident[j].abs() + reflected[j].abs());
    let continuity = |a: &OrientedCrystalField,
                      b: &OrientedCrystalField,
                      a_bound: [f64; 4],
                      b_bound: [f64; 4]|
     -> Result<()> {
        let electric_scale = [a, b].iter().fold(1.0_f64, |m, field| {
            m.max(field.electric_v_m[0].abs())
                .max(field.electric_v_m[1].abs())
        });
        let magnetic_scale = [a, b].iter().fold(admittance, |m, field| {
            m.max(field.magnetic_a_m[0].abs())
                .max(field.magnetic_a_m[1].abs())
        });
        for j in 0..2 {
            // Independently stored wave amplitudes each carry f64 rounding.
            // At grazing TM incidence, for example, O(1/cosine) incoming and
            // reflected H cancel to a finite trace. Account for those actual
            // terms separately from the coherent-field solver tolerance.
            // Eight ulps per term cover complex propagation and the four-wave
            // trace; the physical flux/passivity gates remain independent.
            let electric_rounding = 8.0 * f64::EPSILON * (a_bound[j] + b_bound[j]);
            let magnetic_rounding = 8.0 * f64::EPSILON * (a_bound[j + 2] + b_bound[j + 2]);
            if !electric_rounding.is_finite()
                || !magnetic_rounding.is_finite()
                || (a.electric_v_m[j] - b.electric_v_m[j]).abs()
                    > tolerance * electric_scale + electric_rounding
                || (a.magnetic_a_m[j] - b.magnetic_a_m[j]).abs()
                    > tolerance * magnetic_scale + magnetic_rounding
            {
                return Err(bad(
                    "absorbing stack tangential Maxwell continuity is unresolved",
                ));
            }
        }
        Ok(())
    };
    for layer in &response.layers {
        poll(cancelled)?;
        let front = layer.face(false)?;
        let back = layer.face(true)?;
        continuity(&left, &front, left_bound, face_term_bound(layer, false)?)?;
        let before = front.poynting_w_m2[2] / response.incident_normal_power_w_m2;
        let after = back.poynting_w_m2[2] / response.incident_normal_power_w_m2;
        let absorbed = (Dd::from_f64(before) - Dd::from_f64(after)).to_f64();
        if !before.is_finite()
            || !after.is_finite()
            || !absorbed.is_finite()
            || before < -tolerance
            || after < -tolerance
            || absorbed < -tolerance
            || layer.loss_f_m == [0.0; 2] && absorbed.abs() > tolerance
        {
            return Err(bad(
                "absorbing stack layer violates passive normal-power balance",
            ));
        }
        left = back;
        left_bound = face_term_bound(layer, true)?;
    }
    let exit_bound = tangential(
        response.transmitted_electric_v_m,
        response.transmitted_magnetic_a_m,
    )
    .map(|v| v.abs());
    continuity(
        &left,
        &response.plain_field_at(StackRegion::Transmitted, [0.0, 0.0, response.thickness_m])?,
        left_bound,
        exit_bound,
    )
}

// Magnitudes at the actual face, after decaying-side propagation. Raw remote
// anchor amplitudes would overstate the rounding of attenuated contributions.
fn face_term_bound(layer: &AbsorbingStackLayerField, end: bool) -> Result<[f64; 4]> {
    let mut bound = [0.0; 4];
    for (waves, at_anchor, distance) in [
        (&layer.forward, !end, layer.thickness_m),
        (&layer.backward, end, -layer.thickness_m),
    ] {
        for wave in waves {
            let phase = if at_anchor {
                C64::ONE
            } else {
                propagation(wave.normal_wavenumber_rad_m, distance)
                    .ok_or_else(|| bad("absorbing stack face phase is unrepresentable"))?
            };
            for (j, value) in tangential(wave.electric_v_m, wave.magnetic_a_m)
                .into_iter()
                .enumerate()
            {
                bound[j] += (value * phase).abs();
            }
        }
    }
    if bound.iter().any(|v| !v.is_finite()) {
        return Err(bad(
            "absorbing stack face rounding bound is unrepresentable",
        ));
    }
    Ok(bound)
}

#[cfg(test)]
mod tests;
