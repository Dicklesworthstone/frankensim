//! Shear-deformable (Timoshenko) straight 3-D beam segment between two of its own centroidal nodes.
//!
//! Bending stiffness is the two-node Timoshenko element of Przemieniecki (1968, section 5.6), exact for
//! end loads on a prismatic beam: in the x-y plane, with Phi_y = 12 E I_z / (G A_sy L^2),
//! K = E I_z / ((1 + Phi_y) L^3) [[12, 6L, -12, 6L], [6L, (4 + Phi_y) L^2, -6L, (2 - Phi_y) L^2], ...],
//! and the same in the x-z plane with Phi_z = 12 E I_y / (G A_sz L^2). Axial and Saint-Venant torsion
//! terms are those of the Euler-Bernoulli `ShellBeam`. The shear areas A_sy, A_sz are cross-section
//! data (the area weighted by the section's transverse shear-stress distribution), never a fitted factor.
//!
//! Inertia integrates both translational and section-rotation kinetic energy with the exact static
//! Timoshenko fields. G -> infinity recovers the cubic Hermite/Rayleigh beam in both K and M. Rigid
//! translation and rotation reproduce the continuum mass and full inertia at any segment length.
//! Explicit `BeamMass::Lumped` retains the old endpoint quadrature for legacy pencil comparisons.
//!
//! DOF order at each node: [u, v, w, theta_x, theta_y, theta_z] in global axes. Rotations about the
//! section axes follow `ShellBeam`: the slope of v is theta_z, the slope of w is -theta_y.
//! No offsets: an off-node beam is attached to a shell by a tie, not by a rigid arm here.
use super::{BeamMass, BeamSection, bad, congruence, cross, dot, inertia, unit};
use crate::PlateError;
use fs_sparse::{Coo, Csr};

/// Shear areas [m^2] of the cross-section for transverse shear along the section y and z axes.
#[derive(Debug, Clone, Copy)]
pub struct ShearAreas {
    /// Shear area for shear along y, which governs bending in the x-y plane (with I_z).
    pub y_m2: f64,
    /// Shear area for shear along z, which governs bending in the x-z plane (with I_y).
    pub z_m2: f64,
}

/// Beam theory of one segment.
#[derive(Debug, Clone, Copy)]
pub enum BeamKind {
    /// Euler-Bernoulli (no shear deformation): the Timoshenko element with Phi = 0.
    EulerBernoulli,
    /// Timoshenko, with the cross-section's shear areas.
    Timoshenko(ShearAreas),
}

/// Global 12x12 stiffness and consistent inertia of one straight centroidal segment from `p0_m` to `p1_m`.
/// `section_z` is the direction of the section's +z axis; its part normal to the chord is used.
///
/// # Errors
/// Invalid section (nonpositive or non-finite E, G, rho, A, I_y, I_z, J or shear area), coincident or
/// non-finite end points, a section direction parallel to the chord, or numeric overflow.
pub fn segment_matrices(
    p0_m: [f64; 3],
    p1_m: [f64; 3],
    section_z: [f64; 3],
    section: BeamSection,
    kind: BeamKind,
) -> Result<([f64; 144], [f64; 144]), PlateError> {
    segment_matrices_with_mass(p0_m, p1_m, section_z, section, kind, BeamMass::Consistent)
}

/// Segment matrices with an explicit inertia choice; see [`segment_matrices`].
/// # Errors
/// The same section, geometry and numeric refusals as [`segment_matrices`].
pub fn segment_matrices_with_mass(
    p0_m: [f64; 3],
    p1_m: [f64; 3],
    section_z: [f64; 3],
    section: BeamSection,
    kind: BeamKind,
    mass: BeamMass,
) -> Result<([f64; 144], [f64; 144]), PlateError> {
    section.validate()?;
    if let BeamKind::Timoshenko(shear) = kind {
        if [shear.y_m2, shear.z_m2]
            .iter()
            .any(|v| !v.is_finite() || *v <= 0.0)
        {
            return Err(bad("Timoshenko beam requires finite positive shear areas"));
        }
    }
    if p0_m
        .iter()
        .chain(&p1_m)
        .chain(&section_z)
        .any(|x| !x.is_finite())
    {
        return Err(bad("invalid beam end points or orientation"));
    }
    let chord = std::array::from_fn(|i| p1_m[i] - p0_m[i]);
    let length = dot(chord, chord).sqrt();
    let x = unit(chord)?;
    let hint = unit(section_z)?;
    let projection = dot(hint, x);
    let z = unit(std::array::from_fn(|i| hint[i] - projection * x[i]))?;
    let axes = [x, cross(z, x), z];
    let (k, m) = local_matrices(section, kind, length, mass)?;
    let mut transform = [0.0; 144];
    for node in 0..2 {
        let base = 6 * node;
        for (a, axis) in axes.iter().enumerate() {
            for (b, &component) in axis.iter().enumerate() {
                transform[12 * (base + a) + base + b] = component;
                transform[12 * (base + 3 + a) + base + 3 + b] = component;
                // The rigid-arm column of `ShellBeam` with a zero arm: kept so that the arithmetic,
                // and so the Euler-Bernoulli limit, is the same bit for bit.
                let mut basis = [0.0; 3];
                basis[b] = 1.0;
                transform[12 * (base + a) + base + 3 + b] = dot(*axis, cross(basis, [0.0; 3]));
            }
        }
    }
    let k = congruence(&k, &transform);
    let m = congruence(&m, &transform);
    if k.iter().chain(m.iter()).any(|x| !x.is_finite()) {
        return Err(bad("beam stiffness/inertia overflow"));
    }
    Ok((k, m))
}

/// Local (chord axes) stiffness and inertia of a segment of length `l`.
fn local_matrices(
    s: BeamSection,
    kind: BeamKind,
    l: f64,
    mass: BeamMass,
) -> Result<([f64; 144], [f64; 144]), PlateError> {
    let mut k = [0.0; 144];
    for (a, b, coefficient) in [
        (0, 6, s.young_pa * s.area_m2 / l),
        (3, 9, s.shear_pa * s.torsion_m4 / l),
    ] {
        k[12 * a + a] += coefficient;
        k[12 * b + b] += coefficient;
        k[12 * a + b] -= coefficient;
        k[12 * b + a] -= coefficient;
    }
    let (shear_y, shear_z) = match kind {
        BeamKind::EulerBernoulli => (None, None),
        BeamKind::Timoshenko(areas) => (Some(areas.y_m2), Some(areas.z_m2)),
    };
    let mut phi_bending = [0.0; 2];
    for (plane, (slots, signs, inertia, shear_area)) in [
        ([1, 5, 7, 11], [1.0; 4], s.iz_m4, shear_y),
        ([2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0], s.iy_m4, shear_z),
    ]
    .into_iter()
    .enumerate()
    {
        let phi = shear_area.map_or(0.0, |area| {
            12.0 * s.young_pa * inertia / (s.shear_pa * area * l * l)
        });
        phi_bending[plane] = phi;
        let block = [
            [12.0, 6.0 * l, -12.0, 6.0 * l],
            [6.0 * l, (4.0 + phi) * l * l, -6.0 * l, (2.0 - phi) * l * l],
            [-12.0, -6.0 * l, 12.0, -6.0 * l],
            [6.0 * l, (2.0 - phi) * l * l, -6.0 * l, (4.0 + phi) * l * l],
        ];
        let coefficient = s.young_pa * inertia / ((1.0 + phi) * l * l * l);
        for i in 0..4 {
            for j in 0..4 {
                k[12 * slots[i] + slots[j]] += coefficient * block[i][j] * signs[i] * signs[j];
            }
        }
    }
    Ok((k, inertia::local_mass(s, l, phi_bending, mass)?))
}

/// Stiffness and inertia of a chain of segments through `nodes_m` (6 DOF per node, global axes, no
/// supports), each segment with the same section, section direction and beam theory.
///
/// # Errors
/// Fewer than two nodes, and every refusal of [`segment_matrices`].
pub fn assemble_beam_chain(
    nodes_m: &[[f64; 3]],
    section_z: [f64; 3],
    section: BeamSection,
    kind: BeamKind,
) -> Result<(Csr, Csr), PlateError> {
    assemble_beam_chain_with_mass(nodes_m, section_z, section, kind, BeamMass::Consistent)
}

/// Beam-chain form of [`segment_matrices_with_mass`].
/// # Errors
/// The same node, section and numeric refusals as [`assemble_beam_chain`].
pub fn assemble_beam_chain_with_mass(
    nodes_m: &[[f64; 3]],
    section_z: [f64; 3],
    section: BeamSection,
    kind: BeamKind,
    mass: BeamMass,
) -> Result<(Csr, Csr), PlateError> {
    if nodes_m.len() < 2 {
        return Err(bad("a beam chain needs at least two nodes"));
    }
    let n = 6 * nodes_m.len();
    let mut k = Coo::new(n, n);
    let mut m = Coo::new(n, n);
    for (e, pair) in nodes_m.windows(2).enumerate() {
        let (ke, me) =
            segment_matrices_with_mass(pair[0], pair[1], section_z, section, kind, mass)?;
        for i in 0..12 {
            for j in 0..12 {
                let (row, col) = (6 * e + i, 6 * e + j);
                if ke[12 * i + j] != 0.0 {
                    k.push(row, col, ke[12 * i + j]);
                }
                if me[12 * i + j] != 0.0 {
                    m.push(row, col, me[12 * i + j]);
                }
            }
        }
    }
    Ok((k.assemble(), m.assemble()))
}

#[cfg(test)]
mod tests {
    use super::super::{ShellBeam, element};
    use super::*;
    use crate::shell::ShellMesh;

    fn section() -> BeamSection {
        BeamSection {
            young_pa: 12e9,
            shear_pa: 0.8e9,
            density_kg_m3: 450.0,
            area_m2: 0.001,
            iy_m4: 1.3e-7,
            iz_m4: 4.5e-8,
            torsion_m4: 8e-8,
        }
    }
    fn shear() -> ShearAreas {
        ShearAreas {
            y_m2: 0.001 * 5.0 / 6.0,
            z_m2: 0.001 * 5.0 / 6.0,
        }
    }
    fn quadratic(a: &[f64; 144], u: &[f64; 12]) -> f64 {
        (0..12)
            .map(|i| u[i] * (0..12).map(|j| a[12 * i + j] * u[j]).sum::<f64>())
            .sum()
    }

    #[test]
    fn consistent_mass_integrates_independent_static_velocity_fields() {
        // Construct polynomial fields directly from EI theta''=-V and
        // v'=theta+V/(GA_s), then integrate their squares analytically. This
        // oracle does not use the element's nodal interpolation or quadrature.
        let l: f64 = 0.31;
        let square_integral = |c: &[f64]| -> f64 {
            c.iter()
                .enumerate()
                .map(|(i, a)| {
                    c.iter()
                        .enumerate()
                        .map(|(j, b)| a * b * l.powi((i + j + 1) as i32) / (i + j + 1) as f64)
                        .sum::<f64>()
                })
                .sum()
        };
        let at_tip = |c: &[f64]| c.iter().rev().fold(0.0, |v, x| v * l + x);
        for kind in [BeamKind::EulerBernoulli, BeamKind::Timoshenko(shear())] {
            for shear_pa in [8e6, 8e8, 8e10] {
                let mut s = section();
                s.shear_pa = shear_pa;
                let mut q = [0.0; 12];
                let mut want = 0.0;
                for (slots, sign, inertia, area, theta, v0) in [
                    (
                        [1, 5, 7, 11],
                        1.0,
                        s.iz_m4,
                        shear().y_m2,
                        [0.7, -0.8, 0.6],
                        -0.2,
                    ),
                    (
                        [2, 4, 8, 10],
                        -1.0,
                        s.iy_m4,
                        shear().z_m2,
                        [-0.4, 0.5, -0.3],
                        0.1,
                    ),
                ] {
                    let gamma = match kind {
                        BeamKind::EulerBernoulli => 0.0,
                        BeamKind::Timoshenko(_) => {
                            -2.0 * s.young_pa * inertia * theta[2] / (s.shear_pa * area)
                        }
                    };
                    let v = [v0, theta[0] + gamma, 0.5 * theta[1], theta[2] / 3.0];
                    for (slot, value) in slots.into_iter().zip([
                        v0,
                        sign * theta[0],
                        at_tip(&v),
                        sign * at_tip(&theta),
                    ]) {
                        q[slot] = value;
                    }
                    want += s.density_kg_m3
                        * (s.area_m2 * square_integral(&v) + inertia * square_integral(&theta));
                }
                for (a, b, inertia, field) in [
                    (0, 6, s.area_m2, [0.4, -0.7]),
                    (3, 9, s.iy_m4 + s.iz_m4, [-0.6, 0.8]),
                ] {
                    q[a] = field[0];
                    q[b] = at_tip(&field);
                    want += s.density_kg_m3 * inertia * square_integral(&field);
                }
                let (_, m) =
                    segment_matrices([0.0; 3], [l, 0.0, 0.0], [0.0, 0.0, 1.0], s, kind).unwrap();
                assert!((quadratic(&m, &q) / want - 1.0).abs() < 3e-13);
                // J controls torsional stiffness, not rotational mass.
                s.torsion_m4 *= 17.0;
                let (_, changed) =
                    segment_matrices([0.0; 3], [l, 0.0, 0.0], [0.0, 0.0, 1.0], s, kind).unwrap();
                assert_eq!(m, changed);
            }
        }
    }

    #[test]
    fn consistent_mass_recovers_rayleigh_limit_and_exact_rigid_spin() {
        let l: f64 = 0.31;
        let ends = [[-0.5 * l, 0.0, 0.0], [0.5 * l, 0.0, 0.0]];
        let mut spin = [0.0; 12];
        spin[1] = -0.5 * l;
        spin[7] = 0.5 * l;
        spin[5] = 1.0;
        spin[11] = 1.0;
        let mut errors = Vec::new();
        for g in [1e10, 1e11, 1e12, 1e13] {
            let mut s = section();
            s.shear_pa = g;
            let (_, m) = segment_matrices(
                ends[0],
                ends[1],
                [0.0, 0.0, 1.0],
                s,
                BeamKind::Timoshenko(shear()),
            )
            .unwrap();
            let (_, eb) = segment_matrices(
                ends[0],
                ends[1],
                [0.0, 0.0, 1.0],
                s,
                BeamKind::EulerBernoulli,
            )
            .unwrap();
            let orbital = s.density_kg_m3 * s.area_m2 * l.powi(3) / 12.0;
            let sectional = s.density_kg_m3 * s.iz_m4 * l;
            assert!((quadratic(&m, &spin) / (orbital + sectional) - 1.0).abs() < 5e-14);
            let (_, lumped) = segment_matrices_with_mass(
                ends[0],
                ends[1],
                [0.0, 0.0, 1.0],
                s,
                BeamKind::EulerBernoulli,
                BeamMass::Lumped,
            )
            .unwrap();
            assert!(((quadratic(&lumped, &spin) - sectional) / orbital - 3.0).abs() < 5e-14);
            errors.push(
                m.iter()
                    .zip(&eb)
                    .map(|(a, b)| (a - b).abs())
                    .fold(0.0, f64::max),
            );
        }
        for pair in errors.windows(2) {
            assert!(
                ((pair[1] / pair[0]).log10() + 1.0).abs() < 0.02,
                "consistent inertia must converge to Rayleigh as 1/G: {errors:?}"
            );
        }
    }

    #[test]
    fn g0_rigid_modes_and_shear_patch() {
        let (p0, p1) = ([0.1, -0.2, 0.3], [0.9, 0.1, 0.5]);
        let (k, _) = segment_matrices(
            p0,
            p1,
            [0.0, 0.0, 1.0],
            section(),
            BeamKind::Timoshenko(shear()),
        )
        .unwrap();
        let scale = k.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
        let ends = [p0, p1];
        for axis in 0..3 {
            for rotational in [false, true] {
                let mut direction = [0.0; 3];
                direction[axis] = 1.0;
                let mut u = [0.0; 12];
                for node in 0..2 {
                    let d = if rotational {
                        cross(direction, ends[node])
                    } else {
                        direction
                    };
                    u[6 * node..6 * node + 3].copy_from_slice(&d);
                    if rotational {
                        u[6 * node + 3..6 * node + 6].copy_from_slice(&direction);
                    }
                }
                for i in 0..12 {
                    let residual = (0..12).map(|j| k[12 * i + j] * u[j]).sum::<f64>();
                    assert!(residual.abs() < 1e-12 * scale, "rigid residual {residual}");
                }
            }
        }
        // Shear patch: a cantilever (node 0 clamped) along x under a tip force P along y deflects by
        // P L^3 / (3 E I_z) + P L / (G A_sy) (exact for this element), with tip rotation P L^2 / (2 E I_z).
        let l = 0.25;
        let s = section();
        let a = shear();
        let (k, _) = segment_matrices(
            [0.0; 3],
            [l, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            s,
            BeamKind::Timoshenko(a),
        )
        .unwrap();
        // Solve the 2x2 (v2, theta_z2) block.
        let (k11, k12, k22) = (k[12 * 7 + 7], k[12 * 7 + 11], k[12 * 11 + 11]);
        let p = 1.0;
        let det = k11 * k22 - k12 * k12;
        let (v, theta) = (k22 * p / det, -k12 * p / det);
        let want_v = p * l.powi(3) / (3.0 * s.young_pa * s.iz_m4) + p * l / (s.shear_pa * a.y_m2);
        let want_theta = p * l * l / (2.0 * s.young_pa * s.iz_m4);
        assert!(
            (v / want_v - 1.0).abs() < 1e-13,
            "tip deflection {v} vs {want_v}"
        );
        assert!(
            (theta / want_theta - 1.0).abs() < 1e-13,
            "tip rotation {theta} vs {want_theta}"
        );
    }

    #[test]
    fn g0_euler_bernoulli_limit_bitwise_structure() {
        // Explicit endpoint quadrature reproduces the legacy on-node `ShellBeam` bit for bit.
        let nodes = vec![[0.1, -0.2, 0.3], [0.9, 0.1, 0.5], [0.0, 1.0, 0.0]];
        let mesh = ShellMesh::new(nodes.clone(), vec![[0, 1, 2]]).unwrap();
        let beam = ShellBeam {
            nodes: [0, 1],
            offsets_m: [[0.0; 3]; 2],
            section_z: [0.0, 0.0, 1.0],
            section: section(),
        };
        let (k_eb, m_eb) = element(&mesh, &beam).unwrap();
        let (k, m) = segment_matrices_with_mass(
            nodes[0],
            nodes[1],
            [0.0, 0.0, 1.0],
            section(),
            BeamKind::EulerBernoulli,
            BeamMass::Lumped,
        )
        .unwrap();
        assert!(k.iter().zip(&k_eb).all(|(a, b)| a.to_bits() == b.to_bits()));
        assert!(m.iter().zip(&m_eb).all(|(a, b)| a.to_bits() == b.to_bits()));
        // G -> infinity: M does not change at all, and K approaches the Euler-Bernoulli K as 1/G.
        let mut errors = Vec::new();
        for g in [1e9, 1e10, 1e11, 1e12] {
            let mut s = section();
            s.shear_pa = g;
            let mut eb = s;
            eb.shear_pa = g;
            let (kt, mt) = segment_matrices_with_mass(
                nodes[0],
                nodes[1],
                [0.0, 0.0, 1.0],
                s,
                BeamKind::Timoshenko(shear()),
                BeamMass::Lumped,
            )
            .unwrap();
            let (ke, me) = segment_matrices_with_mass(
                nodes[0],
                nodes[1],
                [0.0, 0.0, 1.0],
                eb,
                BeamKind::EulerBernoulli,
                BeamMass::Lumped,
            )
            .unwrap();
            assert!(mt.iter().zip(&me).all(|(a, b)| a.to_bits() == b.to_bits()));
            let scale = ke.iter().fold(0.0_f64, |a, x| a.max(x.abs()));
            errors.push(
                kt.iter()
                    .zip(&ke)
                    .fold(0.0_f64, |a, (x, y)| a.max((x - y).abs()))
                    / scale,
            );
        }
        for pair in errors.windows(2) {
            let slope = (pair[1] / pair[0]).log10();
            assert!(
                (slope + 1.0).abs() < 0.05,
                "K difference slope {slope} per decade of G"
            );
        }
        // A bending state with no shear (pure curvature) has the same energy in both theories.
        let mut bending = [0.0; 12];
        bending[5] = -1.0;
        bending[11] = 1.0;
        let (kt, _) = segment_matrices(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            section(),
            BeamKind::Timoshenko(shear()),
        )
        .unwrap();
        let (ke, _) = segment_matrices(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            section(),
            BeamKind::EulerBernoulli,
        )
        .unwrap();
        assert!((quadratic(&kt, &bending) / quadratic(&ke, &bending) - 1.0).abs() < 1e-13);
    }

    #[test]
    fn g0_section_refusals() {
        let ok = segment_matrices(
            [0.0; 3],
            [1.0, 0.0, 0.0],
            [0.0, 0.0, 1.0],
            section(),
            BeamKind::Timoshenko(shear()),
        );
        assert!(ok.is_ok());
        for areas in [
            ShearAreas {
                y_m2: 0.0,
                z_m2: 1e-3,
            },
            ShearAreas {
                y_m2: 1e-3,
                z_m2: -1e-3,
            },
            ShearAreas {
                y_m2: f64::NAN,
                z_m2: 1e-3,
            },
            ShearAreas {
                y_m2: 1e-3,
                z_m2: f64::INFINITY,
            },
        ] {
            assert!(
                segment_matrices(
                    [0.0; 3],
                    [1.0, 0.0, 0.0],
                    [0.0, 0.0, 1.0],
                    section(),
                    BeamKind::Timoshenko(areas)
                )
                .is_err()
            );
        }
        let mut s = section();
        s.shear_pa = 0.0;
        assert!(
            segment_matrices(
                [0.0; 3],
                [1.0, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                s,
                BeamKind::EulerBernoulli
            )
            .is_err()
        );
        assert!(
            segment_matrices(
                [0.0; 3],
                [0.0; 3],
                [0.0, 0.0, 1.0],
                section(),
                BeamKind::EulerBernoulli
            )
            .is_err()
        );
        assert!(
            segment_matrices(
                [0.0; 3],
                [1.0, 0.0, 0.0],
                [1.0, 0.0, 0.0],
                section(),
                BeamKind::EulerBernoulli
            )
            .is_err()
        );
        assert!(
            segment_matrices(
                [0.0; 3],
                [f64::NAN, 0.0, 0.0],
                [0.0, 0.0, 1.0],
                section(),
                BeamKind::EulerBernoulli
            )
            .is_err()
        );
        assert!(
            assemble_beam_chain(
                &[[0.0; 3]],
                [0.0, 0.0, 1.0],
                section(),
                BeamKind::EulerBernoulli
            )
            .is_err()
        );
    }
}
