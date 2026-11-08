use super::*;
use crate::shell::stiffened::beam::{RoundBeamModes, RoundBeamSpec};

fn quadratic(a: &[f64; 16], x: [f64; 4]) -> f64 {
    (0..4)
        .map(|i| x[i] * (0..4).map(|j| a[4 * i + j] * x[j]).sum::<f64>())
        .sum()
}

fn free_matrix(a: &Csr, fixed: &[usize]) -> Vec<f64> {
    let ids: Vec<_> = (0..a.nrows()).filter(|i| !fixed.contains(i)).collect();
    let mut out = vec![0.0; ids.len() * ids.len()];
    for (i, &row) in ids.iter().enumerate() {
        let (columns, values) = a.row(row);
        for (&column, &value) in columns.iter().zip(values) {
            if let Ok(j) = ids.binary_search(&column) {
                out[i * ids.len() + j] = value;
            }
        }
    }
    out
}

fn modes(pencil: &RoundBeamPencil, fixed: &[usize]) -> Vec<fs_modal::ModePair> {
    fs_modal::eigh_gen_dense(
        &free_matrix(&pencil.k, fixed),
        &free_matrix(&pencil.m, fixed),
        pencil.k.nrows() - fixed.len(),
    )
    .unwrap()
}

#[test]
fn g0_signed_graded_force_has_exact_geometric_work_and_original_mass() {
    let l = 0.7;
    let radii = [0.003, 0.009];
    let e = [5e9, 21e9];
    let rho = [430.0, 1800.0];
    let original = round_element_with_reference_material(l, radii, e, rho).unwrap();
    assert_eq!(
        original,
        round_element_with_axial_force(l, radii, e, rho, [0.0, -0.0]).unwrap()
    );
    // Independent polynomial derivatives in t=x/l, not the production Gauss rule.
    let d = [
        [0.0, -6.0 / l, 6.0 / l],
        [1.0, -4.0, 3.0],
        [0.0, 6.0 / l, -6.0 / l],
        [0.0, -2.0, 3.0],
    ];
    for force in [[180.0, 220.0], [-180.0, -220.0], [220.0, -180.0]] {
        let (k, m) = round_element_with_axial_force(l, radii, e, rho, force).unwrap();
        assert_eq!(m, original.1);
        let geometric = std::array::from_fn(|i| k[i] - original.0[i]);
        for i in 0..4 {
            for j in 0..4 {
                let mut expected = 0.0;
                for a in 0..3 {
                    for b in 0..3 {
                        expected += l
                            * d[i][a]
                            * d[j][b]
                            * (force[0] / (a + b + 1) as f64
                                + (force[1] - force[0]) / (a + b + 2) as f64);
                    }
                }
                assert!((geometric[4 * i + j] - expected).abs() < 5e-12 * expected.abs().max(1.0));
            }
        }
        assert!(quadratic(&geometric, [1.0, 0.0, 1.0, 0.0]).abs() < 1e-10);
        // A dead axial load is direction-fixed: rotation is NOT a null mode.
        assert!(
            (quadratic(&geometric, [0.0, 1.0, l, 1.0]) - l * (force[0] + force[1]) / 2.0).abs()
                < 1e-10
        );
    }
}

#[test]
fn g0_force_interfaces_and_attached_mass_survive_subdivision() {
    let stations = [
        RoundStation {
            x_m: 0.0,
            radius_m: 0.003,
        },
        RoundStation {
            x_m: 0.37,
            radius_m: 0.008,
        },
        RoundStation {
            x_m: 0.8,
            radius_m: 0.004,
        },
    ];
    let e = [[12e9, 6e9], [25e9, 18e9]];
    let rho = [[400.0, 800.0], [1600.0, 900.0]];
    let force = [[180.0, 220.0], [-150.0, -80.0]];
    let attached = [RoundInertia {
        x_m: 0.61,
        mass_kg: 0.01,
        rotary_kg_m2: 2e-7,
    }];
    // Independent integral N(x)*(2x)^2 dx for the exact trial field w=x^2.
    let mut expected = 0.0;
    for (pair, n) in stations.windows(2).zip(force) {
        let x = pair[0].x_m;
        let l = pair[1].x_m - x;
        expected += 4.0
            * l
            * (n[0] * (x * x + x * l + l * l / 3.0)
                + (n[1] - n[0]) * (x * x / 2.0 + 2.0 * x * l / 3.0 + l * l / 4.0));
    }
    for subdivisions in [2, 4, 8] {
        let unloaded = assemble_prestressed_round_beam(
            &stations,
            subdivisions,
            &e,
            &rho,
            &[[0.0; 2]; 2],
            &attached,
        )
        .unwrap();
        let loaded =
            assemble_prestressed_round_beam(&stations, subdivisions, &e, &rho, &force, &attached)
                .unwrap();
        assert_eq!(loaded.m, unloaded.m);
        assert_eq!(loaded.nodes_m, unloaded.nodes_m);
        assert!(loaded.nodes_m.contains(&0.37));
        let x: Vec<_> = loaded
            .nodes_m
            .iter()
            .flat_map(|x| [x * x, 2.0 * x])
            .collect();
        let a = free_matrix(&loaded.k, &[]);
        let b = free_matrix(&unloaded.k, &[]);
        let n = x.len();
        let work: f64 = (0..n)
            .map(|i| {
                x[i] * (0..n)
                    .map(|j| (a[n * i + j] - b[n * i + j]) * x[j])
                    .sum::<f64>()
            })
            .sum();
        assert!((work - expected).abs() < 2e-9 * expected.abs().max(1.0));
    }
}

#[test]
fn g3_unloaded_full_pencil_reproduces_original_pin_supported_spectrum() {
    let stations = [
        RoundStation {
            x_m: 0.0,
            radius_m: 0.005,
        },
        RoundStation {
            x_m: 0.12,
            radius_m: 0.005,
        },
        RoundStation {
            x_m: 0.4,
            radius_m: 0.005,
        },
    ];
    let spec = RoundBeamSpec {
        young_pa: 12e9,
        density_kg_m3: 800.0,
        pivot_m: 0.12,
        contact_m: 0.39,
        hand_m: 0.2,
        subdivisions: 4,
        maximum_hz: 1500.0,
        maximum_modes: 12,
    };
    let e = [[spec.young_pa; 2]; 2];
    let rho = [[spec.density_kg_m3; 2]; 2];
    let attached = [RoundInertia {
        x_m: 0.35,
        mass_kg: 0.01,
        rotary_kg_m2: 2e-7,
    }];
    let original =
        RoundBeamModes::with_reference_material(&stations, spec, &e, &rho, &attached).unwrap();
    let pencil =
        assemble_prestressed_round_beam(&stations, 4, &e, &rho, &[[0.0; 2]; 2], &attached).unwrap();
    assert_eq!(pencil.nodes_m, original.nodes_m);
    let pivot = 2 * pencil
        .nodes_m
        .iter()
        .position(|x| *x == spec.pivot_m)
        .unwrap();
    let full = modes(&pencil, &[pivot]);
    assert!(full[0].lambda.abs() < 1e-8 * full[1].lambda);
    for (mode, omega) in full.iter().zip(&original.omega).skip(1) {
        assert!((mode.lambda.sqrt() / omega - 1.0).abs() < 2e-7);
    }
}

fn uniform(subdivisions: usize, force: f64) -> RoundBeamPencil {
    assemble_prestressed_round_beam(
        &[
            RoundStation {
                x_m: 0.0,
                radius_m: 0.005,
            },
            RoundStation {
                x_m: 0.6,
                radius_m: 0.005,
            },
        ],
        subdivisions,
        &[[12e9; 2]],
        &[[800.0; 2]],
        &[[force; 2]],
        &[],
    )
    .unwrap()
}

#[test]
fn g1_tension_and_compression_modes_converge_with_original_rayleigh_mass() {
    let area = std::f64::consts::PI * 0.005_f64.powi(2);
    let inertia = area * 0.005_f64.powi(2) / 4.0;
    let q = std::f64::consts::PI / 0.6;
    let critical = 12e9 * inertia * q * q;
    for fraction in [-0.6, 0.0, 1.3] {
        let force = fraction * critical;
        // Exact pinned-pinned Rayleigh continuum: w=sin(q*x), theta=w'.
        let expected =
            (12e9 * inertia * q.powi(4) + force * q * q) / (800.0 * (area + inertia * q * q));
        let mut errors = Vec::new();
        for subdivisions in [4, 8, 16] {
            let pencil = uniform(subdivisions, force);
            let fixed = [0, pencil.k.nrows() - 2];
            let first = modes(&pencil, &fixed)[0].lambda;
            assert!(first > 0.0);
            errors.push((first / expected - 1.0).abs());
        }
        assert!(
            errors.windows(2).all(|e| e[1] < 0.2 * e[0]),
            "frequency ladder: {errors:?}"
        );
        assert!(errors[2] < 6e-6, "frequency accuracy: {errors:?}");
    }
}

#[test]
fn g2_euler_buckling_converges_and_unstable_tangents_remain_signed() {
    let exact = 12e9 * std::f64::consts::PI * 0.005_f64.powi(4) / 4.0
        * (std::f64::consts::PI / 0.6).powi(2);
    let mut errors = Vec::new();
    for subdivisions in [4, 8, 16] {
        let unloaded = uniform(subdivisions, 0.0);
        let tension = uniform(subdivisions, 1.0);
        let n = unloaded.k.nrows() - 2;
        let fixed = [0, unloaded.k.nrows() - 2];
        let k = free_matrix(&unloaded.k, &fixed);
        let kg: Vec<_> = free_matrix(&tension.k, &fixed)
            .iter()
            .zip(&k)
            .map(|(a, b)| a - b)
            .collect();
        // Reuse the native symmetric-definite service on K phi = P Kg phi;
        // positive Kg is the geometric work of ONE newton in tension.
        let critical = fs_modal::eigh_gen_dense(&k, &kg, n).unwrap()[0].lambda;
        errors.push((critical / exact - 1.0).abs());
    }
    assert!(
        errors.windows(2).all(|e| e[1] < 0.2 * e[0]),
        "buckling ladder: {errors:?}"
    );
    assert!(errors[2] < 4e-6, "buckling accuracy: {errors:?}");
    for (factor, stable) in [(0.99, true), (1.01, false)] {
        let pencil = uniform(16, -factor * exact);
        let first = &modes(&pencil, &[0, pencil.k.nrows() - 2])[0];
        if stable {
            assert!(first.interval.0 > 0.0);
        } else {
            assert!(
                first.interval.1 < 0.0,
                "compression must retain its negative eigenvalue"
            );
        }
    }
}

#[test]
fn g0_prestress_refuses_invalid_fields_and_accumulated_overflow() {
    let stations = [
        RoundStation {
            x_m: 0.0,
            radius_m: 0.005,
        },
        RoundStation {
            x_m: 1.0,
            radius_m: 0.005,
        },
        RoundStation {
            x_m: 2.0,
            radius_m: 0.005,
        },
    ];
    for force in [
        vec![],
        vec![[0.0; 2]],
        vec![[0.0; 2]; 3],
        vec![[f64::NAN, 0.0]; 2],
        vec![[0.0, f64::INFINITY]; 2],
    ] {
        assert!(
            assemble_prestressed_round_beam(
                &stations,
                1,
                &[[12e9; 2]; 2],
                &[[800.0; 2]; 2],
                &force,
                &[]
            )
            .is_err()
        );
    }
    assert!(
        assemble_prestressed_round_beam(
            &stations,
            0,
            &[[12e9; 2]; 2],
            &[[800.0; 2]; 2],
            &[[0.0; 2]; 2],
            &[]
        )
        .is_err()
    );
    assert!(
        assemble_prestressed_round_beam(
            &stations,
            1,
            &[[0.0; 2]; 2],
            &[[800.0; 2]; 2],
            &[[0.0; 2]; 2],
            &[]
        )
        .is_err()
    );
    assert!(
        round_element_with_axial_force(1.0, [0.005; 2], [12e9; 2], [800.0; 2], [0.0, f64::NAN])
            .is_err()
    );
    // Each tangent is finite; their shared-node diagonal overflows.
    assert!(
        round_element_with_axial_force(1.0, [0.005; 2], [12e9; 2], [800.0; 2], [8e307; 2]).is_ok()
    );
    assert!(
        assemble_prestressed_round_beam(
            &stations,
            1,
            &[[12e9; 2]; 2],
            &[[800.0; 2]; 2],
            &[[8e307; 2]; 2],
            &[]
        )
        .is_err()
    );
}
