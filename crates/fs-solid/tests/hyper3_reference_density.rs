//! G0/G1/G3: exact graded reference inventory in native finite-strain inertia.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_solid::hyper3::pressure::ClosedTetPressureSurface;
use fs_solid::hyper3::vibration::{VibrationError, VibrationSettings};
use fs_solid::{HyperTetError, HyperTetProblem, HyperTetSettings, HyperelasticLawField};
use fs_sparse::Csr;

const NODES: [[f64; 3]; 4] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
// rho(X) = 3 + 12 X, with mass 1 kg on the unit tetrahedron.
const DENSITY: [f64; 4] = [3.0, 15.0, 3.0, 3.0];

fn with_cx(run: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 31,
                    kernel_id: 17,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        );
    });
}

fn material() -> Hyperelastic {
    Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 2.0,
            lambda: 3.0,
        },
        5.0,
    )
    .unwrap()
}

fn problem<'a>(
    nodes: &'a [[f64; 3]],
    tets: &'a [[usize; 4]],
    material: &'a Hyperelastic,
) -> HyperTetProblem<'a> {
    HyperTetProblem {
        nodes_m: nodes,
        tetrahedra: tets,
        material: HyperelasticLawField::Uniform(material),
        prescribed_m: &[],
        nodal_forces_n: &[],
        body_force_n_m3: [0.0; 3],
        budget: Default::default(),
        settings: HyperTetSettings::default(),
    }
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-13 * expected.abs().max(1.0),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn kinetic(mass: &Csr, velocity: &[f64]) -> f64 {
    let mut momentum = vec![0.0; velocity.len()];
    mass.spmv(velocity, &mut momentum);
    0.5 * velocity
        .iter()
        .zip(momentum)
        .map(|(v, p)| v * p)
        .sum::<f64>()
}

#[test]
fn exact_mass_centroid_and_cubic_moments_retain_supported_coordinates() {
    with_cx(|cx, _| {
        let law = material();
        let mut p = problem(&NODES, &TETS, &law);
        p.prescribed_m = &[(0, 0.0), (1, 0.0), (2, 0.0)];
        let mass = p
            .reference_mass_matrix_with_nodal_density(&DENSITY, cx)
            .unwrap();
        assert_eq!(mass.nrows(), 12);
        // Independent barycentric monomial integrals: product(alpha_i!)/6!
        // for a cubic monomial on this unit tetrahedron.
        let factorial = [1.0, 1.0, 2.0, 6.0];
        for a in 0..4 {
            for b in 0..4 {
                let exact: f64 = DENSITY
                    .iter()
                    .enumerate()
                    .map(|(i, rho)| {
                        let mut powers = [0; 4];
                        powers[i] += 1;
                        powers[a] += 1;
                        powers[b] += 1;
                        rho * powers
                            .into_iter()
                            .map(|power| factorial[power])
                            .product::<f64>()
                            / 720.0
                    })
                    .sum();
                for component in 0..3 {
                    close(mass.get(3 * a + component, 3 * b + component), exact);
                    close(
                        mass.get(3 * a + component, 3 * b + (component + 1) % 3),
                        0.0,
                    );
                }
            }
        }
        let nodal_mass: Vec<_> = (0..4)
            .map(|node| mass.row(3 * node).1.iter().sum::<f64>())
            .collect();
        close(nodal_mass.iter().sum(), 1.0);
        for (actual, expected) in nodal_mass.iter().zip([0.225, 0.325, 0.225, 0.225]) {
            close(*actual, expected);
        }
        for axis in 0..3 {
            close(
                (0..4)
                    .map(|node| nodal_mass[node] * NODES[node][axis])
                    .sum(),
                [0.325, 0.225, 0.225][axis],
            );
        }
        // Equal total mass does not imply the same mass distribution.
        let averaged = p.reference_mass_matrix(&[6.0], cx).unwrap();
        close(averaged.row(3).1.iter().sum(), 0.25);
        assert!((mass.get(3, 3) - averaged.get(3, 3)).abs() > 0.04);
    });
}

#[test]
fn positive_inertia_matches_independent_physical_kinetic_energy() {
    with_cx(|cx, _| {
        let law = material();
        let p = problem(&NODES, &TETS, &law);
        let mass = p
            .reference_mass_matrix_with_nodal_density(&DENSITY, cx)
            .unwrap();
        // Positive Cholesky pivots check the whole coupled shape mass, beyond
        // positive diagonal entries or one chosen velocity.
        let mut lower = [[0.0; 12]; 12];
        for i in 0..12 {
            for j in 0..=i {
                assert_eq!(mass.get(i, j).to_bits(), mass.get(j, i).to_bits());
                let remainder =
                    mass.get(i, j) - (0..j).map(|k| lower[i][k] * lower[j][k]).sum::<f64>();
                lower[i][j] = if i == j {
                    assert!(remainder > 0.0);
                    remainder.sqrt()
                } else {
                    remainder / lower[j][j]
                };
            }
        }
        let velocity = |x: [f64; 3]| {
            [
                1.0 + 2.0 * x[0] - 0.5 * x[1],
                -0.7 + 1.2 * x[2] + 0.3 * x[0],
                0.2 + x[1] - x[2],
            ]
        };
        let nodal_velocity: Vec<_> = NODES.into_iter().flat_map(velocity).collect();
        // Independent tensor Gauss integration after the Duffy map. Three
        // points integrate the transformed degree-five polynomial exactly.
        let offset = (0.6_f64).sqrt() / 2.0;
        let gauss = [
            (0.5 - offset, 5.0 / 18.0),
            (0.5, 4.0 / 9.0),
            (0.5 + offset, 5.0 / 18.0),
        ];
        let mut exact = 0.0;
        for (u, wu) in gauss {
            for (v, wv) in gauss {
                for (w, ww) in gauss {
                    let x = [u, (1.0 - u) * v, (1.0 - u) * (1.0 - v) * w];
                    let jacobian = (1.0 - u).powi(2) * (1.0 - v);
                    let speed2 = velocity(x).iter().map(|value| value * value).sum::<f64>();
                    exact += 0.5 * (3.0 + 12.0 * x[0]) * speed2 * jacobian * wu * wv * ww;
                }
            }
        }
        close(kinetic(&mass, &nodal_velocity), exact);
        assert!(exact > 0.0);
    });
}

#[test]
fn uniform_bits_and_shared_vertex_scatter_survive_element_permutation() {
    with_cx(|cx, _| {
        let law = material();
        let nodes = [NODES[0], NODES[1], NODES[2], NODES[3], [0.0, 0.0, -1.0]];
        let tets = [[0, 1, 2, 3], [0, 2, 1, 4]];
        let mut p = problem(&nodes, &tets, &law);
        let old = p.reference_mass_matrix(&[7.3; 2], cx).unwrap();
        let uniform = p
            .reference_mass_matrix_with_nodal_density(&[7.3; 5], cx)
            .unwrap();
        for (a, b) in old.to_dense().iter().zip(uniform.to_dense()) {
            assert_eq!(a.to_bits(), b.to_bits());
        }
        let rho = [3.0, 15.0, 3.0, 3.0, 9.0];
        let original = p
            .reference_mass_matrix_with_nodal_density(&rho, cx)
            .unwrap();
        let permuted = [[2, 1, 0, 4], [1, 2, 0, 3]];
        p.tetrahedra = &permuted;
        let reordered = p
            .reference_mass_matrix_with_nodal_density(&rho, cx)
            .unwrap();
        for (a, b) in original.to_dense().iter().zip(reordered.to_dense()) {
            close(*a, b);
        }
        close(
            (0..5)
                .map(|node| original.row(3 * node).1.iter().sum::<f64>())
                .sum(),
            2.25,
        );
        // The two tetrahedra contribute opposite z-centroid moments, with
        // greater mass below z=0; shared vertices do not duplicate inventory.
        close(
            (0..5)
                .map(|node| nodes[node][2] * original.row(3 * node).1.iter().sum::<f64>())
                .sum(),
            -0.1,
        );
    });
}

#[test]
fn invalid_inventory_geometry_budget_and_cancellation_return_no_matrix() {
    with_cx(|cx, gate| {
        let law = material();
        let mut p = problem(&NODES, &TETS, &law);
        for density in [
            vec![1.0; 3],
            vec![1.0, 2.0, 3.0, 0.0],
            vec![1.0, 2.0, 3.0, -1.0],
            vec![1.0, 2.0, 3.0, f64::NAN],
            vec![1.0, 2.0, 3.0, f64::INFINITY],
            vec![f64::from_bits(1); 4],
        ] {
            assert!(matches!(
                p.reference_mass_matrix_with_nodal_density(&density, cx),
                Err(HyperTetError::InvalidInput { .. })
            ));
        }
        let orphan = [NODES[0], NODES[1], NODES[2], NODES[3], [1.0; 3]];
        p.nodes_m = &orphan;
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&[1.0; 5], cx),
            Err(HyperTetError::InvalidInput { .. })
        ));
        let large = NODES.map(|point| point.map(|x| 12.0 * x));
        p.nodes_m = &large;
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&[f64::MAX; 4], cx),
            Err(HyperTetError::InvalidInput { .. })
        ));
        p.nodes_m = &NODES;
        p.budget.maximum_nodes = 3;
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&DENSITY, cx),
            Err(HyperTetError::InvalidInput { .. })
        ));
        gate.request();
        assert!(matches!(
            p.reference_mass_matrix_with_nodal_density(&DENSITY, cx),
            Err(HyperTetError::Cancelled)
        ));
    });
}

fn vibration_settings() -> VibrationSettings {
    VibrationSettings {
        maximum_free_dofs: 128,
        relative_residual_tolerance: 1e-9,
    }
}

fn modal_close(actual: f64, expected: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= 2e-8 * expected.abs(),
        "{actual:.17e} != {expected:.17e}"
    );
}

fn window_plan(lower: f64, upper: f64) -> fs_modal::WindowPlan {
    fs_modal::WindowPlan {
        windows: vec![(lower, upper)],
        threads: 1,
        maximum_windows: 1,
        maximum_bytes: 32 * 1024 * 1024,
        maximum_relative_halfwidth: 1e-8,
        coordinates: None,
    }
}

fn window_options() -> fs_modal::SliceOptions {
    fs_modal::SliceOptions {
        max_lanczos: 256,
        ritz_tol: 1e-12,
        ldlt: fs_sparse::LdltOptions {
            max_factor_bytes: Some(16 * 1024 * 1024),
            ..Default::default()
        },
        ..Default::default()
    }
}

#[test]
fn graded_vibration_retains_exact_modal_mass_through_finite_stretch() {
    with_cx(|cx, _| {
        let law = material();
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        for stretch in [1.0_f64, 1.3, 1.6] {
            // F=diag(s,1,1); direct differentiation of the neo-Hookean
            // energy gives P11=2s+(3ln(s)-2)/s and dP11/ds below.
            let piola = 2.0 * stretch + (3.0 * stretch.ln() - 2.0) / stretch;
            let tangent = 2.0 + (2.0 + 3.0 * (1.0 - stretch.ln())) / stretch.powi(2);
            let forces = [(3, piola / 6.0)];
            let mut p = problem(&NODES, &TETS, &law);
            p.prescribed_m = &fixed;
            p.nodal_forces_n = &forces;
            p.settings.load_steps = 4;
            let spectrum = p
                .solve_vibration_with_nodal_density(&DENSITY, vibration_settings(), cx)
                .unwrap();
            assert_eq!(spectrum.free_dofs(), &[3]);
            assert_eq!(spectrum.reference_nodal_density_kg_m3(), Some(&DENSITY[..]));
            close(spectrum.reference_density_kg_m3()[0], 6.0);
            modal_close(1.0 + spectrum.equilibrium().displacement_m[1][0], stretch);
            // Integral (3+12X) X^2 dV = 3/60 + 12/120 = 0.15 kg.
            // This mass is unchanged even when the current volume grows by s.
            close(spectrum.mass_kg()[0], 0.15);
            modal_close(spectrum.stiffness_n_m()[0], tangent / 6.0);
            modal_close(spectrum.modes()[0].lambda, tangent / 0.9);
            close(
                spectrum.modes()[0].phi[0].powi(2) * spectrum.mass_kg()[0],
                1.0,
            );
        }
    });
}

#[test]
fn equal_total_mass_keeps_density_distribution_and_uniform_vibration_bits() {
    with_cx(|cx, _| {
        let law = material();
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let mut p = problem(&NODES, &TETS, &law);
        p.prescribed_m = &fixed;
        let moved = [15.0, 3.0, 3.0, 3.0];
        let dense = p
            .solve_vibration_with_nodal_density(&DENSITY, vibration_settings(), cx)
            .unwrap();
        let relocated = p
            .solve_vibration_with_nodal_density(&moved, vibration_settings(), cx)
            .unwrap();
        assert_eq!(dense.equilibrium(), relocated.equilibrium());
        assert_eq!(dense.stiffness_n_m(), relocated.stiffness_n_m());
        close(
            dense.reference_density_kg_m3()[0],
            relocated.reference_density_kg_m3()[0],
        );
        close(relocated.mass_kg()[0], 1.0 / 12.0);
        close(dense.modes()[0].lambda, 70.0 / 9.0);
        close(relocated.modes()[0].lambda, 14.0);
        let scaled = p
            .solve_vibration_with_nodal_density(
                &DENSITY.map(|rho| 4.0 * rho),
                vibration_settings(),
                cx,
            )
            .unwrap();
        assert_eq!(dense.equilibrium(), scaled.equilibrium());
        assert_eq!(dense.stiffness_n_m(), scaled.stiffness_n_m());
        close(scaled.mass_kg()[0], 4.0 * dense.mass_kg()[0]);
        close(scaled.modes()[0].lambda, 0.25 * dense.modes()[0].lambda);

        let cellwise = p.solve_vibration(&[7.3], vibration_settings(), cx).unwrap();
        let uniform = p
            .solve_vibration_with_nodal_density(&[7.3; 4], vibration_settings(), cx)
            .unwrap();
        assert_eq!(cellwise.reference_nodal_density_kg_m3(), None);
        assert_eq!(
            cellwise.reference_density_kg_m3(),
            uniform.reference_density_kg_m3()
        );
        assert_eq!(cellwise.equilibrium(), uniform.equilibrium());
        assert_eq!(cellwise.mass_kg(), uniform.mass_kg());
        assert_eq!(cellwise.stiffness_n_m(), uniform.stiffness_n_m());
        assert_eq!(
            cellwise.modes()[0].lambda.to_bits(),
            uniform.modes()[0].lambda.to_bits()
        );
        assert_eq!(cellwise.modes()[0].phi, uniform.modes()[0].phi);
    });
}

#[test]
fn graded_initial_iterates_and_sparse_windows_reuse_the_same_physical_pencil() {
    with_cx(|cx, _| {
        let law = material();
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let stretch = 1.2_f64;
        let forces = [(
            3,
            (2.0 * stretch + (3.0 * stretch.ln() - 2.0) / stretch) / 6.0,
        )];
        let mut p = problem(&NODES, &TETS, &law);
        p.prescribed_m = &fixed;
        p.nodal_forces_n = &forces;
        p.settings.load_steps = 1;
        let mut initial = [0.0; 12];
        initial[3] = 0.18; // Deliberately not the accepted endpoint, 0.20.
        let dense = p
            .solve_vibration_from_with_nodal_density(&initial, &DENSITY, vibration_settings(), cx)
            .unwrap();
        modal_close(1.0 + dense.equilibrium().displacement_m[1][0], stretch);
        let lambda = dense.modes()[0].lambda;
        let plan = window_plan(0.8 * lambda, 1.2 * lambda);
        let options = window_options();
        let sparse = p
            .solve_vibration_windows_from_with_nodal_density(
                &initial, &DENSITY, &plan, &options, cx,
            )
            .unwrap();
        assert_eq!(dense.equilibrium(), sparse.equilibrium());
        assert_eq!(sparse.reference_nodal_density_kg_m3(), Some(&DENSITY[..]));
        assert_eq!(
            sparse.reference_density_kg_m3(),
            dense.reference_density_kg_m3()
        );
        assert_eq!(
            sparse.mass_kg().get(0, 0).to_bits(),
            dense.mass_kg()[0].to_bits()
        );
        assert_eq!(
            sparse.stiffness_n_m().get(0, 0).to_bits(),
            dense.stiffness_n_m()[0].to_bits()
        );
        assert_eq!(sparse.spectrum().reports[0].expected, 1);
        modal_close(sparse.spectrum().reports[0].modes[0].lambda, lambda);
        let direct = p
            .solve_vibration_windows_with_nodal_density(&DENSITY, &plan, &options, cx)
            .unwrap();
        modal_close(direct.spectrum().reports[0].modes[0].lambda, lambda);
        close(direct.mass_kg().get(0, 0), 0.15);
    });
}

#[test]
fn graded_pressure_spectrum_keeps_follower_tangent_and_original_mass() {
    with_cx(|cx, _| {
        let law = material();
        let fixed: Vec<_> = (0..12)
            .filter(|dof| ![3, 7, 11].contains(dof))
            .map(|dof| (dof, 0.0))
            .collect();
        let mut p = problem(&NODES, &TETS, &law);
        p.prescribed_m = &fixed;
        p.settings.load_steps = 6;
        let surface = ClosedTetPressureSurface::new(&NODES, &TETS, p.budget, 16, cx).unwrap();
        let stretch = 0.9_f64;
        let log_j = 3.0 * stretch.ln();
        let pressure = -(2.0 * stretch + (3.0 * log_j - 2.0) / stretch) / stretch.powi(2);
        let diagonal = 2.0 + (2.0 + 3.0 * (1.0 - log_j)) / stretch.powi(2);
        let off_diagonal = 3.0 / stretch.powi(2) + pressure * stretch;
        let dense = p
            .solve_vibration_with_nodal_density_and_pressure(
                &surface,
                pressure,
                &DENSITY,
                vibration_settings(),
                cx,
            )
            .unwrap();
        let (retained_surface, retained_pressure) = dense.held_pressure().unwrap();
        assert!(core::ptr::eq(retained_surface, &surface));
        assert_eq!(retained_pressure.to_bits(), pressure.to_bits());
        let masses = [0.15, 1.0 / 12.0, 1.0 / 12.0];
        for i in 0..3 {
            modal_close(1.0 + dense.equilibrium().displacement_m[i + 1][i], stretch);
            for j in 0..3 {
                close(
                    dense.mass_kg()[3 * i + j],
                    if i == j { masses[i] } else { 0.0 },
                );
                modal_close(
                    dense.stiffness_n_m()[3 * i + j],
                    (if i == j { diagonal } else { off_diagonal }) / 6.0,
                );
            }
        }
        // The y-z antisymmetric mode decouples. The other two eigenvalues
        // are the roots of the independent two-dimensional whitened pencil.
        let a = diagonal / (6.0 * masses[0]);
        let d = (diagonal + off_diagonal) / (6.0 * masses[1]);
        let b = 2.0_f64.sqrt() * off_diagonal / (6.0 * (masses[0] * masses[1]).sqrt());
        let gap = ((a - d).powi(2) + 4.0 * b * b).sqrt();
        let mut expected = [
            (a + d - gap) / 2.0,
            (diagonal - off_diagonal) / (6.0 * masses[1]),
            (a + d + gap) / 2.0,
        ];
        expected.sort_by(f64::total_cmp);
        for (mode, target) in dense.modes().iter().zip(expected) {
            modal_close(mode.lambda, target);
        }
        let plan = window_plan(0.5 * expected[0], 1.5 * expected[2]);
        let sparse = p
            .solve_vibration_windows_with_nodal_density_and_pressure(
                &surface,
                pressure,
                &DENSITY,
                &plan,
                &window_options(),
                cx,
            )
            .unwrap();
        assert_eq!(dense.equilibrium(), sparse.equilibrium());
        assert_eq!(sparse.reference_nodal_density_kg_m3(), Some(&DENSITY[..]));
        assert_eq!(sparse.spectrum().reports[0].expected, 3);
        for (mode, target) in sparse.spectrum().reports[0].modes.iter().zip(expected) {
            modal_close(mode.lambda, target);
        }
        assert_eq!(sparse.mass_kg().to_dense(), dense.mass_kg());
        assert_eq!(sparse.stiffness_n_m().to_dense(), dense.stiffness_n_m());
    });
}

#[test]
fn graded_vibration_refuses_invalid_density_before_equilibrium_and_cancels() {
    with_cx(|cx, gate| {
        let law = material();
        let fixed: Vec<_> = (0..12)
            .filter(|&dof| dof != 3)
            .map(|dof| (dof, 0.0))
            .collect();
        let mut p = problem(&NODES, &TETS, &law);
        p.prescribed_m = &fixed;
        p.settings.load_steps = 1;
        let mut invalid_initial = [0.0; 12];
        invalid_initial[3] = 5.0; // Outside the original material's domain.
        let plan = window_plan(1.0, 100.0);
        let options = window_options();
        for density in [
            &DENSITY[..3],
            &[3.0, 0.0, 3.0, 3.0],
            &[3.0, f64::NAN, 3.0, 3.0],
        ] {
            assert!(matches!(
                p.solve_vibration_from_with_nodal_density(
                    &invalid_initial,
                    density,
                    vibration_settings(),
                    cx
                ),
                Err(VibrationError::Solid(HyperTetError::InvalidInput { .. }))
            ));
            assert!(matches!(
                p.solve_vibration_windows_from_with_nodal_density(
                    &invalid_initial,
                    density,
                    &plan,
                    &options,
                    cx
                ),
                Err(VibrationError::Solid(HyperTetError::InvalidInput { .. }))
            ));
        }
        gate.request();
        assert!(matches!(
            p.solve_vibration_with_nodal_density(&DENSITY, vibration_settings(), cx),
            Err(VibrationError::Solid(HyperTetError::Cancelled))
        ));
        assert!(matches!(
            p.solve_vibration_windows_with_nodal_density(&DENSITY, &plan, &options, cx),
            Err(VibrationError::Solid(HyperTetError::Cancelled))
        ));
    });
}
