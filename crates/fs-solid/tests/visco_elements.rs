//! IRP M1.3-W2a: synthetic constitutive fixtures, not specimen validation.
use fs_exec::Cx;
use fs_material::visco::ViscoError;
use fs_material::visco::tensor::{
    KelvinFractional, TensorComplexModulus, TensorLoweringTolerance, lower_tensor_to_prony,
};
use fs_material::visco::tensor::{Mandel6, PronyBranch, TensorProny};
use fs_material::{OrthotropicElastic, visco::FractionalZener};
use fs_math::det;
use fs_solid::linear3::hex::{
    HEX20_REFERENCE_NODES, HEX27_REFERENCE_NODES, HexConnectivity, HexElasticProblem, HexFamily,
    HexQuadrature,
};
use fs_solid::linear3::tet10::{QuadraticTetProblem, TetQuadrature, promote_p1_to_p2};
use fs_solid::linear3::visco::*;
use fs_solid::linear3::{TetAssemblyBudget, TetElasticMaterial, TetMaterialField};

const Q: TetQuadrature = TetQuadrature {
    stiffness_points_per_dir: 3,
    mass_points_per_dir: 4,
};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>, &fs_exec::CancelGate) -> R) -> R {
    let gate = fs_exec::CancelGate::new_clock_free();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey {
                seed: 60,
                kernel_id: 60,
                tile: 0,
                iteration: 0,
            },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        f(&cx, &gate)
    })
}

fn diagonal(v: [f64; 6]) -> Mandel6 {
    core::array::from_fn(|i| core::array::from_fn(|j| if i == j { v[i] } else { 0.0 }))
}
fn fixture(n: usize) -> TetViscoMaterial {
    let eq = diagonal([11e9, 0.9e9, 0.5e9, 1.4e9, 0.12e9, 1.3e9]);
    let branch = PronyBranch::new(diagonal([2e8, 3e7, 2e7, 4e7, 1e7, 3e7]), 1e-3).unwrap();
    TetViscoMaterial {
        law: TensorProny::new(eq, None, vec![branch; n], (0.0, 1e12)).unwrap(),
        density_kg_m3: 450.0,
        material_state_identity: [60; 32],
    }
}
fn nodes() -> Vec<[f64; 3]> {
    vec![[0.0; 3], [0.1, 0.0, 0.0], [0.0, 0.1, 0.0], [0.0, 0.0, 0.1]]
}
fn displacement(nodes: &[[f64; 3]], strain: [f64; 6]) -> Vec<f64> {
    let s = core::f64::consts::FRAC_1_SQRT_2;
    nodes
        .iter()
        .flat_map(|p| {
            [
                strain[0] * p[0] + s * (strain[3] * p[1] + strain[5] * p[2]),
                strain[1] * p[1] + s * (strain[3] * p[0] + strain[4] * p[2]),
                strain[2] * p[2] + s * (strain[5] * p[0] + strain[4] * p[1]),
            ]
        })
        .collect()
}
fn quadratic(a: &fs_sparse::Csr, u: &[f64]) -> f64 {
    let mut y = vec![0.0; u.len()];
    a.spmv(u, &mut y);
    u.iter().zip(y).map(|(x, y)| x * y).sum()
}
fn q6(a: &Mandel6, u: &[f64; 6]) -> f64 {
    (0..6)
        .map(|i| u[i] * (0..6).map(|j| a[i][j] * u[j]).sum::<f64>())
        .sum()
}
fn close(a: f64, b: f64, rel: f64) {
    assert!(
        (a - b).abs() <= rel * a.abs().max(b.abs()).max(1e-25),
        "{a:e} != {b:e}, relative={:e}",
        (a - b).abs() / a.abs().max(b.abs())
    );
}
fn system(mesh: SolidMeshRef<'_>, m: &TetViscoMaterial, cx: &Cx<'_>) -> ViscoElements {
    ViscoElements::new(
        mesh,
        TetViscoField::Uniform(m),
        &[],
        TetAssemblyBudget::standard(),
        ViscoStateBudget {
            max_internal_variables: 100_000,
        },
        cx,
    )
    .unwrap()
}

#[test]
fn g0_hex_tensor_spectra_and_history_preserve_material_work() {
    with_cx(|cx, _| {
        let physical = |p: [i8; 3]| {
            [
                0.1 * f64::from(p[0]),
                0.05 * f64::from(p[1]),
                0.025 * f64::from(p[2]),
            ]
        };
        let n20 = HEX20_REFERENCE_NODES.map(physical);
        let n27 = HEX27_REFERENCE_NODES.map(physical);
        let h20 = [core::array::from_fn(|i| i)];
        let h27 = [core::array::from_fn(|i| i)];
        let mut material = fixture(2);
        material.law = material
            .law
            .rotated(fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37))
            .unwrap();
        let strain = [2e-5, -1e-5, 3e-5, 4e-5, -2e-5, 1e-5];
        let volume = 0.001;
        for (nodes, elements) in [
            (n20.as_slice(), HexConnectivity::Hex20(&h20)),
            (n27.as_slice(), HexConnectivity::Hex27(&h27)),
        ] {
            let solid = system(
                SolidMeshRef::Hex {
                    nodes_m: nodes,
                    elements,
                    quadrature: HexQuadrature { points_per_dir: 3 },
                },
                &material,
                cx,
            );
            assert_eq!(solid.internal_variables(), 27 * 2 * 6);
            let u = displacement(nodes, strain);
            for omega in [0.0, 1000.0, 1e8] {
                let k = solid.complex_stiffness(omega, cx).unwrap();
                let (storage, loss) = material.law.complex_modulus(omega).unwrap();
                close(
                    quadratic(&k.real, &u),
                    volume * q6(&storage, &strain),
                    2e-12,
                );
                if omega > 0.0 {
                    close(quadratic(&k.imag, &u), volume * q6(&loss, &strain), 2e-12);
                } else {
                    assert_eq!(k.real, solid.equilibrium().stiffness);
                }
            }
            let mut history = solid.state(cx).unwrap();
            let mut point = material.law.state();
            for (amplitude, dt) in [(1.0, 0.0002), (-0.6, 0.001), (0.3, 0.02), (0.3, 0.005)] {
                let epsilon = strain.map(|v| amplitude * v);
                let stress = material.law.step(&mut point, epsilon, dt).unwrap();
                let report = history
                    .step_internal_force(&displacement(nodes, epsilon), dt, cx)
                    .unwrap();
                assert!(report.step_dissipated_j >= 0.0);
                close(report.energy.stored_j, volume * point.stored, 2e-12);
                close(report.energy.work_j, volume * point.work, 2e-12);
                close(report.energy.dissipated_j, volume * point.dissipated, 2e-12);
                close(
                    report.energy.work_j,
                    report.energy.stored_j + report.energy.dissipated_j,
                    2e-12,
                );
                // Virtual work of the full assembled force, including all
                // boundary rows, is the original material stress contraction.
                let work: f64 = report.full_force_n.iter().zip(&u).map(|(f, v)| f * v).sum();
                close(
                    work,
                    volume * stress.iter().zip(strain).map(|(s, e)| s * e).sum::<f64>(),
                    2e-12,
                );
            }
        }
    });
}

#[test]
fn g0_hex_visco_reuses_elastic_mass_and_admission() {
    with_cx(|cx, _| {
        let physical = |p: [i8; 3]| p.map(|v| 0.1 * f64::from(v));
        let n20 = HEX20_REFERENCE_NODES.map(physical);
        let n27 = HEX27_REFERENCE_NODES.map(physical);
        let h20 = [core::array::from_fn(|i| i)];
        let h27 = [core::array::from_fn(|i| i)];
        for (nodes, elements, family) in [
            (
                n20.as_slice(),
                HexConnectivity::Hex20(&h20),
                HexFamily::Serendipity20,
            ),
            (
                n27.as_slice(),
                HexConnectivity::Hex27(&h27),
                HexFamily::Lagrange27,
            ),
        ] {
            let material = fixture(0);
            let elastic = TetElasticMaterial::try_new_mandel(
                material.density_kg_m3,
                *material.law.c_eq(),
                material.material_state_identity,
            )
            .unwrap();
            let fixed = [0, 1, 2];
            let quadrature = HexQuadrature { points_per_dir: 3 };
            let expected = HexElasticProblem {
                nodes_m: nodes,
                family,
                elements,
                materials: TetMaterialField::Uniform(&elastic),
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature,
            }
            .assemble(cx)
            .unwrap();
            let make = |material: &TetViscoMaterial, points_per_dir, cap| {
                ViscoElements::new(
                    SolidMeshRef::Hex {
                        nodes_m: nodes,
                        elements,
                        quadrature: HexQuadrature { points_per_dir },
                    },
                    TetViscoField::Uniform(material),
                    &fixed,
                    TetAssemblyBudget::standard(),
                    ViscoStateBudget {
                        max_internal_variables: cap,
                    },
                    cx,
                )
            };
            let solid = make(&material, 3, 0).unwrap();
            assert_eq!(solid.equilibrium().stiffness, expected.stiffness);
            assert_eq!(solid.equilibrium().mass, expected.mass);
            assert_eq!(solid.equilibrium().free_dofs, expected.free_dofs);
            let damped = fixture(1);
            assert_eq!(make(&damped, 3, 162).unwrap().internal_variables(), 162);
            assert!(matches!(
                make(&damped, 3, 161),
                Err(ViscoElementError::InternalStateBudget {
                    needed: 162,
                    cap: 161
                })
            ));
            assert!(matches!(
                make(&damped, 2, 10000),
                Err(ViscoElementError::Elastic(
                    fs_solid::linear3::TetElasticError::UnderIntegrated { .. }
                ))
            ));
        }
    });
}

#[test]
fn g1_curved_hex_loss_quadrature_converges_to_analytic_integral() {
    with_cx(|cx, _| {
        // x=L*xi*(1+a*eta^2), y=L*eta, z=L*zeta is an exactly represented
        // Q2 geometry. u_x=b*xi^2 gives rational physical strain, so refining
        // quadrature changes the answer even though the interpolation is fixed.
        let (length, warp, amplitude) = (0.1, 0.2_f64, 1e-5);
        let nodes = HEX27_REFERENCE_NODES.map(|p| {
            let [x, y, z] = p.map(f64::from);
            [length * x * (1.0 + warp * y * y), length * y, length * z]
        });
        let cells = [core::array::from_fn(|i| i)];
        let displacement: Vec<f64> = HEX27_REFERENCE_NODES
            .iter()
            .flat_map(|p| [amplitude * f64::from(p[0]).powi(2), 0.0, 0.0])
            .collect();
        let material = fixture(1);
        let (storage, loss) = material.law.complex_modulus(1000.0).unwrap();
        // Integrate xi and zeta analytically; the two remaining eta integrals
        // are I0=int 1/(1+a*eta^2), I2=int eta^2/(1+a*eta^2), over [-1,1].
        let i0 = 2.0 * warp.sqrt().atan() / warp.sqrt();
        let i2 = (2.0 - i0) / warp;
        let axial = (16.0 / 3.0) * amplitude * amplitude * length * i0;
        let shear = (32.0 / 5.0) * amplitude * amplitude * length * warp * warp * i2;
        let expected = [
            storage[0][0] * axial + storage[3][3] * shear,
            loss[0][0] * axial + loss[3][3] * shear,
        ];
        let mut previous = [f64::INFINITY; 2];
        for order in [3, 5, 7] {
            let solid = system(
                SolidMeshRef::Hex {
                    nodes_m: &nodes,
                    elements: HexConnectivity::Hex27(&cells),
                    quadrature: HexQuadrature {
                        points_per_dir: order,
                    },
                },
                &material,
                cx,
            );
            let k = solid.complex_stiffness(1000.0, cx).unwrap();
            let values = [
                quadratic(&k.real, &displacement),
                quadratic(&k.imag, &displacement),
            ];
            for channel in 0..2 {
                let error = (values[channel] / expected[channel] - 1.0).abs();
                assert!(
                    error < previous[channel],
                    "order={order}, channel={channel}, error={error:e}, previous={:e}",
                    previous[channel]
                );
                if order == 7 {
                    assert!(error < 1e-8, "channel={channel}, error={error:e}");
                }
                previous[channel] = error;
            }
        }
    });
}

/// Freudenthal cells, with x=L, y=R, z=T. This fixed spatial mesh isolates
/// material lowering error; it makes no claim of spatial convergence for wood.
fn bar(cells: usize, extent: [f64; 3], cx: &Cx<'_>) -> (Vec<[f64; 3]>, Vec<[usize; 10]>) {
    let nodes: Vec<_> = (0..=cells)
        .flat_map(|i| {
            (0..4).map(move |j| {
                [
                    extent[0] * i as f64 / cells as f64,
                    extent[1] * (j / 2) as f64,
                    extent[2] * (j % 2) as f64,
                ]
            })
        })
        .collect();
    let mut tets = Vec::new();
    for i in 0..cells {
        let v = [
            4 * i,
            4 * i + 4,
            4 * i + 2,
            4 * i + 6,
            4 * i + 1,
            4 * i + 5,
            4 * i + 3,
            4 * i + 7,
        ];
        tets.extend([
            [v[0], v[1], v[3], v[7]],
            [v[0], v[3], v[2], v[7]],
            [v[0], v[2], v[6], v[7]],
            [v[0], v[6], v[4], v[7]],
            [v[0], v[4], v[5], v[7]],
            [v[0], v[5], v[1], v[7]],
        ]);
    }
    promote_p1_to_p2(&nodes, &tets, None, cx).unwrap()
}

fn dense(a: &fs_sparse::Csr) -> Vec<f64> {
    let mut result = vec![0.0; a.nrows() * a.ncols()];
    for i in 0..a.nrows() {
        let (cols, values) = a.row(i);
        for (&j, &v) in cols.iter().zip(values) {
            result[i * a.ncols() + j] = v;
        }
    }
    result
}

/// M0.2-W2b's Kelvin wood fixture, expressed in the bar's L/R/T frame.
/// Constants and Kelvin-mode ordering match prony_tensor_lowering.rs. They are
/// synthetic inputs, not fitted specimen parameters or measured uncertainty.
fn fractional_wood() -> KelvinFractional {
    let ortho = OrthotropicElastic::new(
        [11e9, 0.90e9, 0.50e9],
        [0.37, 0.47, 0.43],
        [0.70e9, 0.06e9, 0.65e9],
        0.01,
    )
    .unwrap();
    let (evals, _) = fs_la::eigen::jacobi_eigh(&ortho.stiffness().concat(), 6);
    let alpha = [0.25, 0.30, 0.35, 0.40, 0.45, 0.50];
    let eta = [0.005, 0.008, 0.011, 0.014, 0.017, 0.020];
    let tau = [1.2e-3, 9e-4, 7e-4, 5e-4, 3.5e-4, 2.5e-4];
    let modes = core::array::from_fn(|i| {
        let angle = alpha[i] * core::f64::consts::FRAC_PI_2;
        let peak_im = det::sin(angle) / (2.0 * (1.0 + det::cos(angle)));
        let delta = eta[i] * evals[i] / (peak_im - 0.5 * eta[i]);
        FractionalZener::new(evals[i], evals[i] + delta, alpha[i], tau[i]).unwrap()
    });
    KelvinFractional::from_eigenbasis_of(&ortho.stiffness(), modes).unwrap()
}

#[test]
fn g0_isotropic_zener_bar_reduction() {
    with_cx(|cx, _| {
        // nu=0 and transverse rollers isolate the longitudinal rod mode.
        // TET10 retains its full consistent mass. The analytic real-frequency
        // complex eigenvalue is lambda*(w) = (pi/2L)^2 E*(w)/rho.
        // This is not a claim about the damped pole at complex frequency.
        let (nodes, elements) = bar(16, [1.0, 0.1, 0.1], cx);
        let fixed: Vec<_> = nodes
            .iter()
            .enumerate()
            .flat_map(|(i, p)| {
                (0..3)
                    .filter(move |&axis| axis != 0 || p[0] == 0.0)
                    .map(move |axis| 3 * i + axis)
            })
            .collect();
        let (e0, rho, alpha, tau) = (1e7, 1000.0, 0.8, 0.01);
        let c0 = TetElasticMaterial::try_new(rho, e0, 0.0, [60; 32]).unwrap();
        let stiffness = *c0.stiffness_mandel_pa();
        let material = TetViscoMaterial {
            law: TensorProny::new(
                stiffness,
                None,
                vec![PronyBranch::new(stiffness.map(|row| row.map(|x| alpha * x)), tau).unwrap()],
                (0.0, 1e9),
            )
            .unwrap(),
            density_kg_m3: rho,
            material_state_identity: [60; 32],
        };
        let system = ViscoElements::new(
            SolidMeshRef::Tet10 {
                nodes_m: &nodes,
                elements: &elements,
                quadrature: Q,
            },
            TetViscoField::Uniform(&material),
            &fixed,
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 16 * 6 * 27 * 6,
            },
            cx,
        )
        .unwrap();
        let eq = system.equilibrium();
        let first =
            fs_modal::eigh_gen_dense(&dense(&eq.stiffness), &dense(&eq.mass), eq.free_dofs.len())
                .unwrap()
                .remove(0);
        let analytic_lambda = core::f64::consts::FRAC_PI_2.powi(2) * e0 / rho;
        close(first.lambda, analytic_lambda, 1e-6);
        for omega in [20.0, 100.0, 1000.0] {
            let x = omega * tau;
            let storage = 1.0 + alpha * x * x / (1.0 + x * x);
            let loss = alpha * x / (1.0 + x * x);
            let k = system.complex_stiffness(omega, cx).unwrap();
            let modal_mass = quadratic(&eq.mass, &first.phi);
            let real = quadratic(&k.real, &first.phi) / modal_mass;
            let imag = quadratic(&k.imag, &first.phi) / modal_mass;
            close(real, analytic_lambda * storage, 1e-6);
            close(imag, analytic_lambda * loss, 1e-6);
            let complex_error = det::sqrt(
                (real - analytic_lambda * storage).powi(2)
                    + (imag - analytic_lambda * loss).powi(2),
            ) / (analytic_lambda * det::sqrt(storage * storage + loss * loss));
            println!(
                "{{\"suite\":\"irp-s1\",\"test\":\"g0_isotropic_zener_bar_reduction\",\"tier\":\"G0\",\"omega_rad_s\":{omega},\"relative_complex_eigenvalue_error\":{:.9e},\"tolerance\":1e-6}}",
                complex_error
            );
        }
    });
}

#[test]
fn g1_prony_ladder_modal_loss_certified() {
    with_cx(|cx, _| {
        let started = std::time::Instant::now();
        let target = fractional_wood();
        let (nodes, elements) = bar(1, [0.24, 0.02, 0.01], cx);
        let fixed: Vec<_> = nodes
            .iter()
            .enumerate()
            .filter(|(_, p)| p[0] == 0.0)
            .flat_map(|(i, _)| [3 * i, 3 * i + 1, 3 * i + 2])
            .collect();
        let assemble = |tensor| {
            let material = TetElasticMaterial::try_new_mandel(450.0, tensor, [60; 32]).unwrap();
            QuadraticTetProblem {
                nodes_m: &nodes,
                elements: &elements,
                materials: TetMaterialField::Uniform(&material),
                fixed_dofs: &fixed,
                budget: TetAssemblyBudget::standard(),
                quadrature: Q,
            }
            .assemble(cx)
            .unwrap()
        };
        let eq = assemble(target.equilibrium_modulus(1.0).unwrap());
        let mass = dense(&eq.mass);
        let n = eq.free_dofs.len();
        let initial = fs_modal::eigh_gen_dense(&dense(&eq.stiffness), &mass, n).unwrap();
        // Independent direct-fractional reference: evaluate C*(omega), then
        // use contraction iteration (not the production Prony bisection).
        // The analytic storage log-slope below is < 2, proving contraction.
        let mut reference = Vec::new();
        for (index, initial_mode) in initial.iter().take(6).enumerate() {
            let mut omega = det::sqrt(initial_mode.lambda);
            let mut solved = None;
            for _ in 0..16 {
                let (storage, loss) = target.complex_modulus(omega).unwrap();
                let real = assemble(storage).stiffness;
                let pairs = fs_modal::eigh_gen_dense(&dense(&real), &mass, n).unwrap();
                let pair = &pairs[index];
                assert!(pair.residual / pair.lambda < 1e-9);
                let next = det::sqrt(pair.lambda);
                if (next / omega - 1.0).abs() < 1e-10 {
                    let k = ComplexStiffness {
                        real,
                        imag: assemble(loss).stiffness,
                        omega_rad_s: omega,
                    };
                    let eta = modal_loss_energy_quotient(&k, &pair.phi).unwrap();
                    solved = Some((omega, eta, pairs));
                    break;
                }
                omega = next;
            }
            let mode = solved.expect("direct fractional storage iteration converges");
            assert!((20.0..20_000.0).contains(&(mode.0 / core::f64::consts::TAU)));
            reference.push(mode);
        }
        // For each scalar fractional Zener with r=(omega*tau)^alpha:
        // |dE*/dlnw|=alpha*DeltaE*r/(1+2*r*cos(alpha*pi/2)+r²)
        // <= alpha*DeltaE/[2*(1+cos(alpha*pi/2))]. C' >= E0.
        // Also |dln(E'')/dlnw| <= alpha. These Loewner bounds survive
        // positive-weight assembly in this fixed Kelvin eigenbasis.
        let storage_slope = target
            .modes
            .iter()
            .map(|m| {
                m.alpha * (m.e_inf - m.e0)
                    / (2.0 * (1.0 + det::cos(m.alpha * core::f64::consts::FRAC_PI_2)) * m.e0)
            })
            .fold(0.0f64, f64::max);
        let loss_slope = target.modes.iter().map(|m| m.alpha).fold(0.0f64, f64::max);
        assert!(storage_slope < 2.0);
        let mut previous = (f64::INFINITY, f64::INFINITY);
        let mut previous_per_mode = [(f64::INFINITY, f64::INFINITY); 6];
        for terms in [2, 4, 8, 16] {
            let lowered = lower_tensor_to_prony(
                &target,
                20.0,
                20_000.0,
                terms,
                TensorLoweringTolerance {
                    storage_rel: 0.5,
                    loss_rel: 0.99,
                },
            )
            .unwrap();
            let eps_s = lowered.certificate.total_storage_bound();
            let eps_l = lowered.certificate.total_loss_bound();
            let material = TetViscoMaterial {
                law: lowered.law,
                density_kg_m3: 450.0,
                material_state_identity: [60; 32],
            };
            assert_eq!(material.law.branches().len(), terms);
            let system = ViscoElements::new(
                SolidMeshRef::Tet10 {
                    nodes_m: &nodes,
                    elements: &elements,
                    quadrature: Q,
                },
                TetViscoField::Uniform(&material),
                &fixed,
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: 6 * 27 * 6 * terms,
                },
                cx,
            )
            .unwrap();
            assert_eq!(system.equilibrium().mass, eq.mass);
            let modes = system
                .modes(
                    ViscoModeBudget {
                        modes: 6,
                        max_free_dofs: n,
                        max_iterations: 64,
                        frequency_bracket_rad_s: (
                            core::f64::consts::TAU * 20.0,
                            core::f64::consts::TAU * 20_000.0,
                        ),
                        relative_frequency_tolerance: 1e-8,
                    },
                    cx,
                )
                .unwrap();
            // Min-max plus the storage log-slope bounds the movement of each
            // self-consistent root: |ln(w_fit/w_true)| <= -ln(1-eps_s)/(2-Ls).
            // Numerical solve margins are explicit; these are not outward
            // roundoff certificates from fs-modal.
            let delta = -det::ln(1.0 - eps_s) / (2.0 - storage_slope) + 1e-7;
            let frequency_bound = det::exp(delta) - 1.0;
            let tensor_delta = ((1.0 + eps_s) * det::exp(storage_slope * delta) - 1.0)
                .max(1.0 - (1.0 - eps_s) * det::exp(-storage_slope * delta))
                + 1e-7;
            let quotient_bound = ((1.0 + eps_l) / (1.0 - eps_s)
                * det::exp((storage_slope + loss_slope) * delta)
                - 1.0)
                .max(
                    1.0 - (1.0 - eps_l) / (1.0 + eps_s)
                        * det::exp(-(storage_slope + loss_slope) * delta),
                );
            let mut worst = (0.0f64, 0.0f64);
            for (index, (mode, (omega, eta, spectrum))) in modes.iter().zip(&reference).enumerate()
            {
                let frequency_error = (mode.omega_rad_s / omega - 1.0).abs();
                let loss_error = (mode.loss_factor / eta - 1.0).abs();
                // In target energy coordinates, the off-mode component obeys
                // sin(theta) <= delta_K * max_{j!=i} lambda_j/|lambda_j-lambda_fit|.
                // This follows by expanding (A-lambda_fit)y=-Ey in the FULL
                // target eigenbasis and using ||A^-1/2 E A^-1/2|| <= delta_K.
                let lambda_fit = mode.storage_eigenvalue_per_s2;
                let gap_factor = spectrum
                    .iter()
                    .enumerate()
                    .filter(|(j, _)| *j != index)
                    .map(|(_, p)| {
                        let gap = (p.lambda - lambda_fit).abs() - p.residual;
                        if gap > 0.0 {
                            p.interval.1 / gap
                        } else {
                            f64::INFINITY
                        }
                    })
                    .fold(0.0f64, f64::max);
                // Charge the approximate eigenvector's M^-1 residual too:
                // ||A^-1/2 r|| / ||A^1/2 y|| <= ||r|| / lambda_min(A)
                // for this mass-normalized vector. Root and eigen residuals
                // are distinct; omega² cannot replace the Ritz value here.
                assert!(spectrum[0].interval.0 > 0.0);
                let solver_delta = mode.eigen_residual_per_s2 / spectrum[0].interval.0;
                let angle = ((tensor_delta + solver_delta) * gap_factor).min(1.0);
                let etas: Vec<_> = target
                    .modes
                    .iter()
                    .map(|m| {
                        let (re, im) = m.modulus(*omega);
                        im / re
                    })
                    .collect();
                let eta_min = etas.iter().copied().fold(f64::INFINITY, f64::min);
                let eta_max = etas.iter().copied().fold(0.0f64, f64::max);
                // For energy-unit vectors u,v, P_v-P_u has eigenvalues
                // +/-sin(theta). After subtracting eta_min*I, the trace
                // bound is (eta_max-eta_min)*sin(theta), including the
                // quadratic-in-sin term; it is not a first-order expansion.
                let shape_bound = (eta_max - eta_min) * angle / eta;
                let loss_bound = quotient_bound * (1.0 + shape_bound) + shape_bound;
                assert!(
                    frequency_error <= frequency_bound,
                    "N={terms} mode={index}: frequency {frequency_error:e} > {frequency_bound:e}"
                );
                assert!(
                    loss_error <= loss_bound,
                    "N={terms} mode={index}: loss {loss_error:e} > {loss_bound:e}"
                );
                assert!(
                    frequency_error < previous_per_mode[index].0
                        && loss_error < previous_per_mode[index].1,
                    "N={terms} mode={index}: frequency/loss errors ({frequency_error:e}, {loss_error:e}) do not decrease from {:?}",
                    previous_per_mode[index]
                );
                previous_per_mode[index] = (frequency_error, loss_error);
                worst.0 = worst.0.max(frequency_error);
                worst.1 = worst.1.max(loss_error);
                println!(
                    "{{\"suite\":\"irp-s1\",\"test\":\"g1_prony_ladder_modal_loss_certified\",\"tier\":\"G1\",\"terms\":{terms},\"mode\":{index},\"frequency_hz\":{:.9e},\"loss_factor\":{:.9e},\"frequency_error\":{frequency_error:.9e},\"loss_error\":{loss_error:.9e},\"frequency_bound\":{frequency_bound:.9e},\"loss_bound\":{loss_bound:.9e},\"storage_certificate\":{eps_s:.9e},\"loss_certificate\":{eps_l:.9e},\"elapsed_s\":{:.6}}}",
                    mode.omega_rad_s / core::f64::consts::TAU,
                    mode.loss_factor,
                    started.elapsed().as_secs_f64()
                );
            }
            assert!(
                worst.0 < previous.0 && worst.1 < previous.1,
                "N={terms}: six-mode worst frequency/loss errors {worst:?} do not decrease from {previous:?}"
            );
            if terms == 16 {
                assert!(worst.0 <= 5e-4, "frequency error {worst:?}");
                assert!(worst.1 <= 1e-2, "loss error {worst:?}");
            }
            previous = worst;
        }
    });
}

#[test]
fn g0_uniform_strain_complex_energy_exact() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(1);
        let strain = [2e-5, -1e-5, 3e-5, 4e-5, -2e-5, 1e-5];
        for (mesh, xyz) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                &v,
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                &p,
            ),
        ] {
            let s = system(mesh, &m, cx);
            let u = displacement(xyz, strain);
            for w in [0.0, 1000.0, 1e12] {
                let k = s.complex_stiffness(w, cx).unwrap();
                let (r, i) = m.law.complex_modulus(w).unwrap();
                close(quadratic(&k.real, &u), q6(&r, &strain) * 0.001 / 6.0, 1e-13);
                if w > 0.0 {
                    close(quadratic(&k.imag, &u), q6(&i, &strain) * 0.001 / 6.0, 1e-13);
                }
                if w == 0.0 {
                    assert_eq!(k.real, s.equilibrium().stiffness);
                }
                if w == 1e12 {
                    let instant = core::array::from_fn(|i| {
                        core::array::from_fn(|j| {
                            m.law.c_eq()[i][j] + m.law.branches()[0].c_pa[i][j]
                        })
                    });
                    close(
                        quadratic(&k.real, &u),
                        q6(&instant, &strain) * 0.001 / 6.0,
                        1e-13,
                    );
                }
            }
        }
    });
}

#[test]
fn g0_empty_branches_bitwise_equal_elastic() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(0);
        for mesh in [
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            SolidMeshRef::Tet10 {
                nodes_m: &p,
                elements: &t10,
                quadrature: Q,
            },
        ] {
            let s = system(mesh, &m, cx);
            let k = s.complex_stiffness(1000.0, cx).unwrap();
            assert_eq!(k.real, s.equilibrium().stiffness);
            assert!(k.imag.row(0).1.iter().all(|v| *v == 0.0));
        }
    });
}

#[test]
fn g0_passivity_under_random_multiaxial_histories() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(3);
        for (mesh, n) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                v.len(),
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                p.len(),
            ),
        ] {
            let s = system(mesh, &m, cx);
            let mut state = s.state(cx).unwrap();
            let mut rng = fs_rand::StreamKey {
                seed: 0x4d31_3357_3261,
                kernel: 60,
                tile: n as u32,
            }
            .stream();
            for _ in 0..4096 {
                let u: Vec<f64> = (0..3 * n).map(|_| 2e-7 * (rng.next_f64() - 0.5)).collect();
                let out = state
                    .step_internal_force(&u, 1e-5 + 1e-3 * rng.next_f64(), cx)
                    .unwrap();
                assert!(out.step_dissipated_j >= 0.0);
                close(
                    out.energy.work_j,
                    out.energy.stored_j + out.energy.dissipated_j,
                    1e-12,
                );
            }
            println!(
                "{{\"suite\":\"irp-s1\",\"test\":\"g0_passivity_under_random_multiaxial_histories\",\"tier\":\"G0\",\"nodes\":{n},\"steps\":4096,\"seed\":\"0x4d3133573261\",\"work_balance_tolerance\":1e-12,\"verdict\":\"pass\"}}"
            );
        }
    });
}

#[test]
fn g1_time_domain_element_force_matches_complex_stiffness() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let m = fixture(1);
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let w = 1000.0;
        let u = displacement(&v, [1e-5, -2e-5, 3e-5, 1e-5, 2e-5, -1e-5]);
        let k = s.complex_stiffness(w, cx).unwrap();
        let mut kr = vec![0.0; 12];
        let mut ki = kr.clone();
        k.real.spmv(&u, &mut kr);
        k.imag.spmv(&u, &mut ki);
        let mut errors = Vec::new();
        for n in [32, 64, 128] {
            let mut state = s.state(cx).unwrap();
            let dt = core::f64::consts::TAU / (w * f64::from(n));
            let mut err = 0.0;
            let mut norm = 0.0;
            for j in 1..=20 * n {
                let phase = w * dt * f64::from(j);
                let (a, b) = (det::sin(phase), det::cos(phase));
                let out = state
                    .step_internal_force(&u.iter().map(|v| v * a).collect::<Vec<_>>(), dt, cx)
                    .unwrap();
                if j > 19 * n {
                    for i in 0..12 {
                        let exact = kr[i] * a + ki[i] * b;
                        let difference = out.full_force_n[i] - exact;
                        err += difference * difference;
                        norm += exact * exact;
                    }
                }
            }
            errors.push(det::sqrt(err / norm));
        }
        for pair in errors.windows(2) {
            let order = det::ln(pair[0] / pair[1]) / det::ln(2.0);
            assert!(
                (order - 2.0).abs() < 0.2,
                "errors={errors:?}, order={order}"
            );
        }
        println!(
            "{{\"suite\":\"irp-s1\",\"test\":\"g1_time_domain_element_force_matches_complex_stiffness\",\"tier\":\"G1\",\"samples_per_cycle\":[32,64,128],\"relative_force_errors\":{errors:?},\"order_band\":[1.8,2.2],\"verdict\":\"pass\"}}"
        );
    });
}

#[test]
fn g3_rotation_metamorphic() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let m = fixture(2);
        let r = fs_material::tensor::rotation([0.0, 0.0, 1.0], 0.37);
        let rot = |x: [f64; 3]| core::array::from_fn(|i| (0..3).map(|j| r[i][j] * x[j]).sum());
        let vr: Vec<_> = v.iter().copied().map(rot).collect();
        let mr = TetViscoMaterial {
            law: m.law.rotated(r).unwrap(),
            ..m.clone()
        };
        let a = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let b = system(
            SolidMeshRef::Tet4 {
                nodes_m: &vr,
                elements: &t,
            },
            &mr,
            cx,
        );
        let u = displacement(&v, [1e-5, 2e-5, -1e-5, 3e-5, -2e-5, 4e-5]);
        let ur: Vec<f64> = u
            .chunks_exact(3)
            .flat_map(|p| rot([p[0], p[1], p[2]]))
            .collect();
        let ka = a.complex_stiffness(1000.0, cx).unwrap();
        let kb = b.complex_stiffness(1000.0, cx).unwrap();
        close(quadratic(&ka.real, &u), quadratic(&kb.real, &ur), 1e-12);
        close(
            modal_loss_energy_quotient(&ka, &u).unwrap(),
            modal_loss_energy_quotient(&kb, &ur).unwrap(),
            1e-12,
        );
    });
}

#[test]
fn g0_refuses_frequency_outside_law_band() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let mut m = fixture(1);
        m.law = m.law.with_band((100.0, 2000.0)).unwrap();
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        for w in [0.0, 99.0, 2001.0] {
            assert!(matches!(
                s.complex_stiffness(w, cx),
                Err(ViscoElementError::Material(ViscoError::OutOfBand { .. }))
            ));
        }
        assert!(s.complex_stiffness(f64::NAN, cx).is_err());
        assert!(
            modal_loss_energy_quotient(&s.complex_stiffness(1000.0, cx).unwrap(), &[0.0; 12])
                .is_err()
        );
    });
}

#[test]
fn g0_internal_state_budget_cap_and_cap_plus_one() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        let m = fixture(64);
        for (mesh, needed) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                384,
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                27 * 384,
            ),
        ] {
            let make = |cap| {
                ViscoElements::new(
                    mesh,
                    TetViscoField::Uniform(&m),
                    &[],
                    TetAssemblyBudget::standard(),
                    ViscoStateBudget {
                        max_internal_variables: cap,
                    },
                    cx,
                )
            };
            assert_eq!(make(needed).unwrap().internal_variables(), needed);
            assert!(
                matches!(make(needed-1),Err(ViscoElementError::InternalStateBudget {needed:n,cap}) if n==needed && cap==needed-1)
            );
        }
        assert!(
            TensorProny::new(
                *m.law.c_eq(),
                None,
                vec![m.law.branches()[0]; 65],
                (0.0, 1e12)
            )
            .is_err()
        );
    });
}

#[test]
fn g0_invalid_history_and_cancellation_preserve_state() {
    with_cx(|cx, gate| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let m = fixture(1);
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let mut state = s.state(cx).unwrap();
        state.step_internal_force(&[1e-7; 12], 1e-4, cx).unwrap();
        let energy = state.energy();
        for (u, dt) in [
            (vec![0.0; 11], 1e-3),
            (vec![f64::NAN; 12], 1e-3),
            (vec![0.0; 12], 0.0),
        ] {
            assert!(state.step_internal_force(&u, dt, cx).is_err());
            assert_eq!(state.energy(), energy);
        }
        // Overflows partway through a candidate; the original history survives.
        assert!(
            state
                .step_internal_force(&[f64::MAX; 12], 1e-6, cx)
                .is_err()
        );
        assert_eq!(state.energy(), energy);
        gate.request();
        assert!(state.step_internal_force(&[0.0; 12], 1e-3, cx).is_err());
        assert_eq!(state.energy(), energy);
    });
}

#[test]
fn g0_isotropic_zener_uniform_strain_loss() {
    with_cx(|cx, _| {
        // Uniaxial constrained strain: the independent P-wave Zener modulus is
        // lambda+2mu; this is not a free-side Poisson-relaxation beam oracle.
        let elastic = TetElasticMaterial::try_new(1200.0, 3e9, 0.25, [60; 32]).unwrap();
        let c = *elastic.stiffness_mandel_pa();
        let fraction = 0.2;
        let tau = 0.002;
        let branch = PronyBranch::new(c.map(|row| row.map(|v| fraction * v)), tau).unwrap();
        let m = TetViscoMaterial {
            law: TensorProny::new(c, None, vec![branch], (0.0, 1e12)).unwrap(),
            density_kg_m3: 1200.0,
            material_state_identity: [60; 32],
        };
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let s = system(
            SolidMeshRef::Tet4 {
                nodes_m: &v,
                elements: &t,
            },
            &m,
            cx,
        );
        let u = displacement(&v, [1e-5, 0.0, 0.0, 0.0, 0.0, 0.0]);
        for w in [100.0, 500.0, 2000.0] {
            let x = w * tau;
            let k = s.complex_stiffness(w, cx).unwrap();
            close(
                modal_loss_energy_quotient(&k, &u).unwrap(),
                fraction * x / (1.0 + (1.0 + fraction) * x * x),
                1e-13,
            );
        }
    });
}

#[test]
fn g0_heterogeneous_material_force_energy_and_support_reactions() {
    with_cx(|cx, _| {
        let mut xyz = nodes();
        xyz.extend(nodes().into_iter().map(|p| [p[0] + 0.2, p[1], p[2]]));
        let tets = [[0, 1, 2, 3], [4, 5, 6, 7]];
        let materials = [fixture(1), fixture(3)];
        let s = ViscoElements::new(
            SolidMeshRef::Tet4 {
                nodes_m: &xyz,
                elements: &tets,
            },
            TetViscoField::PerElement(&materials),
            &[0, 1, 2],
            TetAssemblyBudget::standard(),
            ViscoStateBudget {
                max_internal_variables: 24,
            },
            cx,
        )
        .unwrap();
        let strains = [
            [2e-5, -1e-5, 3e-5, 4e-5, -2e-5, 1e-5],
            [-1e-5, 3e-5, -2e-5, 2e-5, 1e-5, -3e-5],
        ];
        let mut u = displacement(&xyz[..4], strains[0]);
        u.extend(displacement(&xyz[4..], strains[1]));
        let mut state = s.state(cx).unwrap();
        let out = state.step_internal_force(&u, 0.0003, cx).unwrap();
        let mut expected_storage = 0.0;
        for e in 0..2 {
            let mut point = materials[e].law.state();
            let stress = materials[e]
                .law
                .step(&mut point, strains[e], 0.0003)
                .unwrap();
            expected_storage += point.stored * 0.001 / 6.0;
            // For this axis-aligned tetrahedron, the positive x/y/z corner
            // forces are its stress-tensor columns times V/h.
            let q = core::f64::consts::FRAC_1_SQRT_2;
            let columns = [
                [stress[0], q * stress[3], q * stress[5]],
                [q * stress[3], stress[1], q * stress[4]],
                [q * stress[5], q * stress[4], stress[2]],
            ];
            for a in 0..3 {
                for c in 0..3 {
                    close(
                        out.full_force_n[12 * e + 3 * (a + 1) + c],
                        columns[a][c] * 0.001 / 0.6,
                        1e-13,
                    );
                }
            }
            for c in 0..3 {
                close(
                    out.full_force_n[12 * e + c],
                    -columns.iter().map(|v| v[c]).sum::<f64>() * 0.001 / 0.6,
                    1e-13,
                );
            }
        }
        close(out.energy.stored_j, expected_storage, 1e-13);
        assert_eq!(out.reduced_force_n, out.full_force_n[3..]);
        // The first element can finish before a later element refuses. Verify
        // internal strains as well as the published ledger survive rollback.
        let mut control = state.clone();
        let mut invalid = u.clone();
        invalid[12..].fill(f64::MAX);
        assert!(state.step_internal_force(&invalid, 0.0002, cx).is_err());
        let next: Vec<_> = u.iter().map(|v| 0.7 * v).collect();
        let actual = state.step_internal_force(&next, 0.0002, cx).unwrap();
        let expected = control.step_internal_force(&next, 0.0002, cx).unwrap();
        assert_eq!(actual.full_force_n, expected.full_force_n);
        assert_eq!(actual.energy, expected.energy);
        let k = s.complex_stiffness(1000.0, cx).unwrap();
        assert_eq!(k.real.nrows(), 21);
        assert_eq!(s.equilibrium().free_dofs, (3..24).collect::<Vec<_>>());
        assert!(
            ViscoElements::new(
                SolidMeshRef::Tet4 {
                    nodes_m: &xyz,
                    elements: &tets
                },
                TetViscoField::PerElement(&materials[..1]),
                &[],
                TetAssemblyBudget::standard(),
                ViscoStateBudget {
                    max_internal_variables: 24
                },
                cx,
            )
            .is_err()
        );
    });
}

#[test]
fn g1_dispersive_zener_mode_matches_closed_form() {
    with_cx(|cx, _| {
        let v = nodes();
        let t = [[0, 1, 2, 3]];
        let (p, t10) = promote_p1_to_p2(&v, &t, None, cx).unwrap();
        for (mesh, n) in [
            (
                SolidMeshRef::Tet4 {
                    nodes_m: &v,
                    elements: &t,
                },
                12,
            ),
            (
                SolidMeshRef::Tet10 {
                    nodes_m: &p,
                    elements: &t10,
                    quadrature: Q,
                },
                30,
            ),
        ] {
            // One retained displacement coordinate yields an independent
            // analytic quadratic for s=omega². Both element orders use their
            // own consistent mass. This is a discrete Zener fixture, not a
            // continuum bar or an empirical wood validation.
            let fixed: Vec<_> = (0..n).filter(|&i| i != 3).collect();
            let elastic = fixture(0);
            let make = |m: &TetViscoMaterial| {
                ViscoElements::new(
                    mesh,
                    TetViscoField::Uniform(m),
                    &fixed,
                    TetAssemblyBudget::standard(),
                    ViscoStateBudget {
                        max_internal_variables: 1000,
                    },
                    cx,
                )
                .unwrap()
            };
            let base = make(&elastic);
            let lambda0 =
                base.equilibrium().stiffness.row(0).1[0] / base.equilibrium().mass.row(0).1[0];
            let omega0 = det::sqrt(lambda0);
            let alpha = 4.0;
            let tau = 1.0 / omega0;
            let c = *elastic.law.c_eq();
            let branch = PronyBranch::new(c.map(|r| r.map(|v| alpha * v)), tau).unwrap();
            let m = TetViscoMaterial {
                law: TensorProny::new(c, None, vec![branch], (0.0, 1e12)).unwrap(),
                ..elastic
            };
            let s = make(&m);
            let b = ViscoModeBudget {
                modes: 1,
                max_free_dofs: 1,
                max_iterations: 64,
                frequency_bracket_rad_s: (omega0 / 10.0, omega0 * 10.0),
                relative_frequency_tolerance: 1e-10,
            };
            let modes = s.modes(b, cx).unwrap();
            let mode = &modes[0];
            let expected = omega0 * det::sqrt(0.5 * (alpha + det::sqrt(alpha * alpha + 4.0)));
            close(mode.omega_rad_s, expected, 2e-10);
            let x = expected * tau;
            close(
                mode.loss_factor,
                alpha * x / (1.0 + (1.0 + alpha) * x * x),
                3e-10,
            );
            assert!(mode.relative_frequency_residual <= 4e-10);
            assert!(
                mode.omega_rad_s > 2.0 * omega0,
                "equilibrium-only modes must fail this fixture"
            );
            assert!(matches!(
                s.modes(
                    ViscoModeBudget {
                        max_iterations: 1,
                        ..b
                    },
                    cx
                ),
                Err(ViscoElementError::ModalConvergence { mode: 0 })
            ));
            assert!(matches!(
                s.modes(
                    ViscoModeBudget {
                        frequency_bracket_rad_s: (omega0 / 10.0, omega0 / 2.0),
                        ..b
                    },
                    cx
                ),
                Err(ViscoElementError::UnbracketedMode { mode: 0 })
            ));
            assert!(matches!(
                s.modes(
                    ViscoModeBudget {
                        max_free_dofs: 0,
                        ..b
                    },
                    cx
                ),
                Err(ViscoElementError::InvalidModalBudget)
            ));
        }
    });
}
