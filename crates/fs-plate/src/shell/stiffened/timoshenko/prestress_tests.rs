//! G0/G1/G3: original beam-column energy, stability and covariance.
use super::*;
use fs_sparse::{DirectOrdering, LdltOptions, SymbolicLdlt};

fn section() -> BeamSection {
    BeamSection {
        young_pa: 12e9,
        shear_pa: 0.75e9,
        density_kg_m3: 450.0,
        area_m2: 6e-4,
        iy_m4: 4.5e-8,
        iz_m4: 2e-8,
        torsion_m4: 3e-8,
    }
}

fn quadratic(k: &[f64; 144], q: &[f64; 12]) -> f64 {
    (0..12)
        .map(|i| q[i] * (0..12).map(|j| k[12 * i + j] * q[j]).sum::<f64>())
        .sum()
}

fn close(got: f64, expected: f64, rel: f64) {
    assert!(
        (got - expected).abs() <= rel * expected.abs().max(1e-12),
        "{got:e} != {expected:e}"
    );
}

#[test]
fn geometric_energy_matches_independent_static_fields_and_euler_matrix() {
    let l: f64 = 0.43;
    let force = 8_000.0;
    let s = section();
    let areas = ShearAreas {
        y_m2: 4e-4,
        z_m2: 5e-4,
    };
    for kind in [BeamKind::EulerBernoulli, BeamKind::Timoshenko(areas)] {
        let mut kg = [0.0; 144];
        add_initial_stress(&mut kg, s, kind, l, force).unwrap();
        let mut q = [0.0; 12];
        let mut integral = 0.0;
        for (slots, sign, inertia, area, theta, v0) in [
            (
                [1, 5, 7, 11],
                1.0,
                s.iz_m4,
                areas.y_m2,
                [0.07, -0.08, 0.06],
                -0.02,
            ),
            (
                [2, 4, 8, 10],
                -1.0,
                s.iy_m4,
                areas.z_m2,
                [-0.04, 0.05, -0.03],
                0.01,
            ),
        ] {
            // Independently prescribe theta=a+b*x+c*x², then use the
            // governing equations V=-EI*theta'' and v'=theta+V/(GA_s).
            let shear_slope = match kind {
                BeamKind::EulerBernoulli => 0.0,
                BeamKind::Timoshenko(_) => {
                    -2.0 * s.young_pa * inertia * theta[2] / (s.shear_pa * area)
                }
            };
            let slope = [theta[0] + shear_slope, theta[1], theta[2]];
            let v1 = v0 + slope[0] * l + 0.5 * slope[1] * l * l + slope[2] * l * l * l / 3.0;
            q[slots[0]] = v0;
            q[slots[1]] = sign * theta[0];
            q[slots[2]] = v1;
            q[slots[3]] = sign * (theta[0] + theta[1] * l + theta[2] * l * l);
            for i in 0..3 {
                for j in 0..3 {
                    integral +=
                        slope[i] * slope[j] * l.powi((i + j + 1) as i32) / (i + j + 1) as f64;
                }
            }
        }
        close(quadratic(&kg, &q), force * integral, 2e-13);
        for i in 0..12 {
            for j in 0..12 {
                assert_eq!(kg[12 * i + j].to_bits(), kg[12 * j + i].to_bits());
            }
        }
    }
    // Independent closed-form Euler-Bernoulli initial-stress matrix.
    let mut kg = [0.0; 144];
    add_initial_stress(&mut kg, s, BeamKind::EulerBernoulli, l, force).unwrap();
    let block = [
        [36.0, 3.0 * l, -36.0, 3.0 * l],
        [3.0 * l, 4.0 * l * l, -3.0 * l, -l * l],
        [-36.0, -3.0 * l, 36.0, -3.0 * l],
        [3.0 * l, -l * l, -3.0 * l, 4.0 * l * l],
    ];
    for (slots, signs) in [
        ([1, 5, 7, 11], [1.0; 4]),
        ([2, 4, 8, 10], [1.0, -1.0, 1.0, -1.0]),
    ] {
        for i in 0..4 {
            for j in 0..4 {
                close(
                    kg[12 * slots[i] + slots[j]],
                    force / (30.0 * l) * block[i][j] * signs[i] * signs[j],
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn signed_prestress_preserves_mass_zero_bits_and_rotated_work() {
    let s = section();
    let kind = BeamKind::Timoshenko(ShearAreas {
        y_m2: 4e-4,
        z_m2: 5e-4,
    });
    let l = 0.43;
    let (base, mass) = segment_matrices([0.0; 3], [l, 0.0, 0.0], [0.0, 0.0, 1.0], s, kind).unwrap();
    for zero in [0.0, -0.0] {
        let (k, m) =
            prestressed_segment_matrices([0.0; 3], [l, 0.0, 0.0], [0.0, 0.0, 1.0], s, kind, zero)
                .unwrap();
        assert_eq!(k.map(f64::to_bits), base.map(f64::to_bits));
        assert_eq!(m.map(f64::to_bits), mass.map(f64::to_bits));
    }
    // A proper cyclic rotation (x,y,z)->(z,x,y), plus a translation.
    let rotate = |v: [f64; 3]| [v[2], v[0], v[1]];
    let p0 = [0.4, -0.3, 0.2];
    let p1 = [p0[0], p0[1] + l, p0[2]];
    let q: [f64; 12] = std::array::from_fn(|i| ((i * 17 + 3) as f64).sin() * 0.003);
    let mut rotated = [0.0; 12];
    for block in 0..4 {
        rotated[3 * block..3 * block + 3].copy_from_slice(&rotate([
            q[3 * block],
            q[3 * block + 1],
            q[3 * block + 2],
        ]));
    }
    for force in [-8_000.0, 8_000.0] {
        let (k, m) =
            prestressed_segment_matrices([0.0; 3], [l, 0.0, 0.0], [0.0, 0.0, 1.0], s, kind, force)
                .unwrap();
        let (kr, mr) =
            prestressed_segment_matrices(p0, p1, rotate([0.0, 0.0, 1.0]), s, kind, force).unwrap();
        assert_eq!(m.map(f64::to_bits), mass.map(f64::to_bits));
        close(quadratic(&kr, &rotated), quadratic(&k, &q), 2e-13);
        close(quadratic(&mr, &rotated), quadratic(&m, &q), 2e-13);
        let delta: [f64; 144] = std::array::from_fn(|i| k[i] - base[i]);
        assert_eq!(
            quadratic(&delta, &q).is_sign_positive(),
            force.is_sign_positive()
        );
    }
}

fn pinned_matrix(n: usize, kind: BeamKind, compression: f64) -> Csr {
    let l = 0.6;
    let nodes: Vec<_> = (0..=n)
        .map(|i| [l * i as f64 / n as f64, 0.0, 0.0])
        .collect();
    let (k, _) = assemble_prestressed_beam_chain(
        &nodes,
        [0.0, 0.0, 1.0],
        section(),
        kind,
        &vec![-compression; n],
    )
    .unwrap();
    let dofs: Vec<_> = (0..=n)
        .flat_map(|i| {
            if i == 0 || i == n {
                vec![6 * i + 5]
            } else {
                vec![6 * i + 1, 6 * i + 5]
            }
        })
        .collect();
    let mut out = Coo::new(dofs.len(), dofs.len());
    for (r, &row) in dofs.iter().enumerate() {
        let (cols, vals) = k.row(row);
        for (&col, &v) in cols.iter().zip(vals) {
            if let Some(c) = dofs.iter().position(|&x| x == col) {
                out.push(r, c, v);
            }
        }
    }
    out.assemble()
}

#[test]
fn compression_buckling_converges_to_independent_continuum_load() {
    let s = section();
    let area = 5e-4;
    let euler = s.young_pa * s.iz_m4 * (std::f64::consts::PI / 0.6).powi(2);
    for kind in [
        BeamKind::EulerBernoulli,
        BeamKind::Timoshenko(ShearAreas {
            y_m2: area,
            z_m2: area,
        }),
    ] {
        // Derived from det([[ (GA-P)q²,-GAq],[-GAq,EIq²+GA]])=0;
        // independent of the finite-element interpolation and quadrature.
        let exact = match kind {
            BeamKind::EulerBernoulli => euler,
            BeamKind::Timoshenko(_) => euler / (1.0 + euler / (s.shear_pa * area)),
        };
        let mut errors = Vec::new();
        for n in [8, 16, 32] {
            let positive = |p| {
                let a = pinned_matrix(n, kind, p);
                SymbolicLdlt::analyze(&a, DirectOrdering::Natural)
                    .unwrap()
                    .factor(
                        &a,
                        &LdltOptions {
                            zero_pivot_rel: f64::EPSILON,
                            ..LdltOptions::default()
                        },
                    )
                    .map(|factor| factor.inertia().negative == 0)
            };
            let (mut lo, mut hi) = (0.9 * exact, 1.1 * exact);
            assert!(positive(lo).unwrap());
            assert!(!positive(hi).unwrap());
            for _ in 0..28 {
                let mid = 0.5 * (lo + hi);
                match positive(mid) {
                    Ok(true) => lo = mid,
                    Ok(false) => hi = mid,
                    Err(fs_sparse::LdltError::PivotBreakdown { .. }) if hi - lo < 2e-8 * exact => {
                        break;
                    }
                    Err(e) => panic!("unresolved buckling bracket [{lo},{hi}]: {e:?}"),
                }
            }
            errors.push((0.5 * (lo + hi) / exact - 1.0).abs());
        }
        assert!(
            errors.windows(2).all(|w| w[1] < w[0] / 3.5),
            "buckling h-ladder: {errors:?}"
        );
        assert!(errors[2] < 3e-5, "buckling accuracy: {errors:?}");
    }
}

#[test]
fn prestressed_frequency_energy_converges_with_original_rotary_mass() {
    let s = section();
    let l = 0.6;
    let q = std::f64::consts::PI / l;
    let shear = s.shear_pa * 5e-4;
    for kind in [
        BeamKind::EulerBernoulli,
        BeamKind::Timoshenko(ShearAreas {
            y_m2: 5e-4,
            z_m2: 5e-4,
        }),
    ] {
        for force in [-2_000.0, 0.0, 8_000.0] {
            let (lambda, theta) = match kind {
                BeamKind::EulerBernoulli => (
                    (s.young_pa * s.iz_m4 * q.powi(4) + force * q * q)
                        / (s.density_kg_m3 * (s.area_m2 + s.iz_m4 * q * q)),
                    q,
                ),
                BeamKind::Timoshenko(_) => {
                    let a = (shear + force) * q * q / (s.density_kg_m3 * s.area_m2);
                    let b = (s.young_pa * s.iz_m4 * q * q + shear) / (s.density_kg_m3 * s.iz_m4);
                    let c = shear * q / (s.density_kg_m3 * (s.area_m2 * s.iz_m4).sqrt());
                    let high = 0.5 * (a + b + (a - b).hypot(2.0 * c));
                    let det = q
                        * q
                        * (shear * s.young_pa * s.iz_m4 * q * q
                            + force * (s.young_pa * s.iz_m4 * q * q + shear))
                        / (s.density_kg_m3 * s.density_kg_m3 * s.area_m2 * s.iz_m4);
                    let lambda = det / high;
                    (
                        lambda,
                        shear * q
                            / (s.young_pa * s.iz_m4 * q * q + shear
                                - s.density_kg_m3 * s.iz_m4 * lambda),
                    )
                }
            };
            let mut errors = Vec::new();
            for n in [8, 16, 32, 64] {
                let nodes: Vec<_> = (0..=n)
                    .map(|i| [l * i as f64 / n as f64, 0.0, 0.0])
                    .collect();
                let (k, m) = assemble_prestressed_beam_chain(
                    &nodes,
                    [0.0, 0.0, 1.0],
                    s,
                    kind,
                    &vec![force; n],
                )
                .unwrap();
                let mut v = vec![0.0; 6 * (n + 1)];
                for (i, p) in nodes.iter().enumerate() {
                    v[6 * i + 1] = (q * p[0]).sin();
                    v[6 * i + 5] = theta * (q * p[0]).cos();
                }
                let energy = |a: &Csr| {
                    let mut av = vec![0.0; v.len()];
                    a.spmv(&v, &mut av);
                    v.iter().zip(av).map(|(x, y)| x * y).sum::<f64>()
                };
                errors.push((energy(&k) / energy(&m) / lambda - 1.0).abs());
            }
            assert!(
                errors.windows(2).all(|w| w[1] < w[0] / 3.5),
                "frequency energy h-ladder: {errors:?}"
            );
            assert!(errors[3] < 1e-5, "frequency energy accuracy: {errors:?}");
        }
    }
}

#[test]
fn force_fields_are_complete_finite_and_assembled_overflow_refuses() {
    let nodes = [[0.0; 3], [1.0, 0.0, 0.0], [2.0, 0.0, 0.0]];
    let s = section();
    let original =
        assemble_beam_chain(&nodes, [0.0, 0.0, 1.0], s, BeamKind::EulerBernoulli).unwrap();
    let zero = assemble_prestressed_beam_chain(
        &nodes,
        [0.0, 0.0, 1.0],
        s,
        BeamKind::EulerBernoulli,
        &[0.0, -0.0],
    )
    .unwrap();
    for row in 0..18 {
        assert_eq!(original.0.row(row), zero.0.row(row));
        assert_eq!(original.1.row(row), zero.1.row(row));
    }
    for forces in [
        &[][..],
        &[1.0][..],
        &[1.0, 2.0, 3.0][..],
        &[f64::NAN, 0.0][..],
        &[0.0, f64::INFINITY][..],
    ] {
        assert!(
            assemble_prestressed_beam_chain(
                &nodes,
                [0.0, 0.0, 1.0],
                s,
                BeamKind::EulerBernoulli,
                forces
            )
            .is_err()
        );
    }
    assert!(
        prestressed_segment_matrices(
            nodes[0],
            nodes[1],
            [0.0, 0.0, 1.0],
            s,
            BeamKind::EulerBernoulli,
            f64::NAN
        )
        .is_err()
    );
    // Each element is representable; its shared-node sum is not.
    assert!(
        prestressed_segment_matrices(
            nodes[0],
            nodes[1],
            [0.0, 0.0, 1.0],
            s,
            BeamKind::EulerBernoulli,
            8e307
        )
        .is_ok()
    );
    assert!(
        assemble_prestressed_beam_chain(
            &nodes,
            [0.0, 0.0, 1.0],
            s,
            BeamKind::EulerBernoulli,
            &[8e307, 8e307]
        )
        .is_err()
    );
}
