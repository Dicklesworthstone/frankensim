//! G0/G1/G3 finite Hex8 checks with independent mechanics and energy references.
//! Synthetic coefficients qualify the implementation, not a measured material.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::thermal::{
    AffineThermalFiniteJ2, FiniteJ2TemperatureSlopes, FiniteJ2ThermalOverstress,
    FiniteJ2ThermalParameters,
};
use fs_material::plastic::finite::{FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::hex::{
    FiniteHexFormulation, FiniteThermoplasticHexError, FiniteThermoplasticHexSolid,
};
use fs_solid::plastic3::hex::HexAssemblyBudget;
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 7,
                kernel_id: 177,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ))
    })
}

fn family(mu: f64, bulk: f64, yield_pa: f64, hardening: f64) -> AffineThermalFiniteJ2 {
    AffineThermalFiniteJ2::new(
        FiniteJ2ThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [280.0, 500.0],
            reference: FiniteJ2Parameters {
                shear_modulus_pa: mu,
                bulk_modulus_pa: bulk,
                yield_stress_pa: yield_pa,
                hardening_modulus_pa: hardening,
            },
            slopes: FiniteJ2TemperatureSlopes {
                shear_modulus_pa_k: 0.0,
                bulk_modulus_pa_k: 0.0,
                yield_stress_pa_k: 0.0,
                hardening_modulus_pa_k: 0.0,
            },
            logarithmic_expansion_per_k: 0.0,
            logarithmic_expansion_slope_per_k2: 0.0,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.8,
            maximum_plastic_log_strain: 3.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 2e-6,
        maximum_endpoint_work_defect_j: 1e8,
        yield_tolerance_pa: 2e-6,
    }
}
fn settings() -> IncrementSettings {
    IncrementSettings {
        relative_tolerance: 1e-10,
        ..IncrementSettings::default()
    }
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}
fn mesh(
    n: [usize; 3],
    map: impl Fn(f64, f64, f64) -> [f64; 3],
) -> (Vec<[f64; 3]>, Vec<[usize; 8]>) {
    let id = |i, j, k| (k * (n[1] + 1) + j) * (n[0] + 1) + i;
    let mut nodes = Vec::new();
    let mut cells = Vec::new();
    for k in 0..=n[2] {
        for j in 0..=n[1] {
            for i in 0..=n[0] {
                nodes.push(map(
                    i as f64 / n[0] as f64,
                    j as f64 / n[1] as f64,
                    k as f64 / n[2] as f64,
                ));
            }
        }
    }
    for k in 0..n[2] {
        for j in 0..n[1] {
            for i in 0..n[0] {
                cells.push([
                    id(i, j, k),
                    id(i + 1, j, k),
                    id(i + 1, j + 1, k),
                    id(i, j + 1, k),
                    id(i, j, k + 1),
                    id(i + 1, j, k + 1),
                    id(i + 1, j + 1, k + 1),
                    id(i, j + 1, k + 1),
                ]);
            }
        }
    }
    (nodes, cells)
}
fn prescribed(u: &[[f64; 3]]) -> Vec<(usize, f64)> {
    u.iter().flatten().copied().enumerate().collect()
}
fn distorted() -> (Vec<[f64; 3]>, Vec<[usize; 8]>) {
    mesh([1, 1, 1], |s, t, v| {
        [
            s + 0.12 * s * t + 0.05 * t * v,
            t + 0.09 * s * v,
            v + 0.07 * s * t,
        ]
    })
}
fn motion(nodes: &[[f64; 3]], scale: f64) -> Vec<[f64; 3]> {
    nodes
        .iter()
        .map(|&[x, y, z]| {
            [
                scale * (0.18 * x + 0.07 * y + 0.04 * x * y),
                scale * (-0.03 * y + 0.02 * x * z),
                scale * (0.02 * z + 0.04 * x * y),
            ]
        })
        .collect()
}
fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}
fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    (0..9).map(|i| a[i] * b[i]).sum()
}
fn kirchhoff(p: &[f64; 9], f: &[f64; 9]) -> [f64; 9] {
    std::array::from_fn(|ij| {
        (0..3)
            .map(|k| p[3 * (ij / 3) + k] * f[3 * (ij % 3) + k])
            .sum()
    })
}
fn equivalent(tau: &[f64; 9]) -> f64 {
    let mean = (tau[0] + tau[4] + tau[8]) / 3.0;
    (1.5 * (0..9)
        .map(|i| (tau[i] - if i % 4 == 0 { mean } else { 0.0 }).powi(2))
        .sum::<f64>())
    .sqrt()
}

#[test]
fn nonuniform_yielded_hex_has_the_complete_projected_tangent_and_objective_force() {
    with_cx(|cx| {
        let (nodes, cells) = distorted();
        let initial = FiniteThermoplasticHexSolid::try_new(
            &nodes,
            &cells,
            vec![family(8000.0, 800_000.0, 200.0, 4000.0)],
            &[300.0],
            FiniteHexFormulation::MeanDilatation,
            HexAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let old_u = motion(&nodes, 0.45);
        let old = initial
            .advance(
                &[300.0],
                &prescribed(&old_u),
                &[],
                0.2,
                settings(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let alphas: Vec<_> = old.quadrature_states().iter().map(|s| s.alpha).collect();
        assert!(alphas.iter().all(|a| *a > 0.0));
        assert!(
            alphas.iter().copied().fold(f64::NEG_INFINITY, f64::max)
                - alphas.iter().copied().fold(f64::INFINITY, f64::min)
                > 1e-3
        );
        let u = motion(&nodes, 1.0);
        let result = old.evaluate_trial(&[300.0], &u, 0.2, true, cx).unwrap();
        let matrix = result.tangent_n_m.as_ref().unwrap();
        let cheap = old.evaluate_trial(&[300.0], &u, 0.2, false, cx).unwrap();
        near(
            cheap.incremental_potential_j,
            result.incremental_potential_j,
            1e-9,
            1e-13,
        );
        for column in 0..24 {
            let h = 2e-6;
            let mut plus = u.clone();
            let mut minus = u.clone();
            plus[column / 3][column % 3] += h;
            minus[column / 3][column % 3] -= h;
            let p = old.evaluate_trial(&[300.0], &plus, 0.2, false, cx).unwrap();
            let m = old
                .evaluate_trial(&[300.0], &minus, 0.2, false, cx)
                .unwrap();
            near(
                (p.incremental_potential_j - m.incremental_potential_j) / (2.0 * h),
                result.internal_force_n[column / 3][column % 3],
                2e-4,
                3e-7,
            );
            for row in 0..24 {
                let derivative = (p.internal_force_n[row / 3][row % 3]
                    - m.internal_force_n[row / 3][row % 3])
                    / (2.0 * h);
                near(matrix.get(row, column), derivative, 2e-3, 4e-6);
                near(
                    matrix.get(row, column),
                    matrix.get(column, row),
                    1e-7,
                    1e-11,
                );
            }
        }
        let (s, c) = 0.91_f64.sin_cos();
        let rotate = |v: [f64; 3]| [c * v[0] - s * v[1], s * v[0] + c * v[1], v[2]];
        let rotated: Vec<_> = nodes
            .iter()
            .zip(&u)
            .map(|(x, u)| {
                let p = rotate(std::array::from_fn(|i| x[i] + u[i]));
                std::array::from_fn(|i| p[i] - x[i] + [0.17, -0.08, 0.23][i])
            })
            .collect();
        let objective = old
            .evaluate_trial(&[300.0], &rotated, 0.2, false, cx)
            .unwrap();
        near(
            objective.incremental_potential_j,
            result.incremental_potential_j,
            1e-7,
            2e-11,
        );
        for (actual, original) in objective
            .internal_force_n
            .iter()
            .zip(&result.internal_force_n)
        {
            let expected = rotate(*original);
            for i in 0..3 {
                near(actual[i], expected[i], 2e-6, 2e-10);
            }
        }
        assert_eq!(
            old.quadrature_states()
                .iter()
                .map(|s| s.alpha)
                .collect::<Vec<_>>(),
            alphas
        );
        assert_eq!(old.displacement_m(), old_u);

        let (cube, hexes) = mesh([1, 1, 1], |x, y, z| [x, y, z]);
        let virgin = FiniteThermoplasticHexSolid::try_new(
            &cube,
            &hexes,
            vec![family(8000.0, 16000.0, 200.0, 4000.0)],
            &[300.0],
            FiniteHexFormulation::MeanDilatation,
            HexAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let mut folded = vec![[0.0; 3]; 8];
        folded[7][0] = -1.3;
        // detF = 1 - 1.3*y*z: mean0.675 and Gauss determinants>0.19,
        // but the top corner has detF=-0.3. A valid mean cannot hide it.
        let error = virgin
            .evaluate_trial(&[300.0], &folded, 0.2, false, cx)
            .unwrap_err();
        assert!(matches!(error, FiniteThermoplasticHexError::Solid(
            SmallStrainTetError::InvalidInput { what }) if what.contains("physical Jacobian")));
        assert!(virgin.quadrature_states().iter().all(|s| s.alpha == 0.0));
    });
}

fn cylinder(nr: usize, nt: usize, formulation: FiniteHexFormulation, cx: &Cx<'_>) -> f64 {
    let height = 0.2;
    let (nodes, cells) = mesh([nr, nt, 1], |s, t, v| {
        let angle = 0.5 * core::f64::consts::PI * t;
        [(1.0 + s) * angle.cos(), (1.0 + s) * angle.sin(), height * v]
    });
    let id = |i, j, k| (k * (nt + 1) + j) * (nr + 1) + i;
    let mut fixed = Vec::new();
    for k in 0..=1 {
        for j in 0..=nt {
            for i in 0..=nr {
                let node = id(i, j, k);
                fixed.push((3 * node + 2, 0.0));
                if j == 0 {
                    fixed.push((3 * node + 1, 0.0));
                }
                if j == nt {
                    fixed.push((3 * node, 0.0));
                }
            }
        }
    }
    let mut forces = Vec::new();
    for j in 0..nt {
        let a = id(0, j, 0);
        let b = id(0, j + 1, 0);
        let dx = nodes[b][0] - nodes[a][0];
        let dy = nodes[b][1] - nodes[a][1];
        for node in [a, b, id(0, j, 1), id(0, j + 1, 1)] {
            // Reference-area pressure work; finite geometry effects are O(p/E).
            forces.push((3 * node, 0.01 * height * dy / 4.0));
            forces.push((3 * node + 1, -0.01 * height * dx / 4.0));
        }
    }
    let (young, nu) = (1000.0, 0.4999);
    let law = family(
        young / (2.0 * (1.0 + nu)),
        young / (3.0 * (1.0 - 2.0 * nu)),
        1e6,
        0.0,
    );
    let solid = FiniteThermoplasticHexSolid::try_new(
        &nodes,
        &cells,
        vec![law; cells.len()],
        &vec![300.0; cells.len()],
        formulation,
        HexAssemblyBudget::default(),
        cx,
    )
    .unwrap();
    let settings = IncrementSettings {
        relative_tolerance: 1e-6,
        linear_tolerance: 1e-10,
        ..IncrementSettings::default()
    };
    let step = solid
        .advance(
            &vec![300.0; cells.len()],
            &fixed,
            &forces,
            1.0,
            settings,
            controls(),
            cx,
        )
        .unwrap_or_else(|e| panic!("finite cylinder {nr}x{nt} {formulation:?}: {e}"));
    assert!(
        step.solid
            .quadrature_states()
            .iter()
            .all(|s| s.alpha == 0.0)
    );
    let mut radial = 0.0;
    for k in 0..=1 {
        for j in 0..=nt {
            let n = id(0, j, k);
            let x = nodes[n];
            let u = step.solid.displacement_m()[n];
            radial += (x[0] * u[0] + x[1] * u[1]) / x[0].hypot(x[1]);
        }
    }
    radial / (2 * (nt + 1)) as f64
}

#[test]
fn finite_mean_dilatation_recovers_the_nearly_incompressible_lame_limit_and_refines() {
    with_cx(|cx| {
        let nu = 0.4999;
        let exact = (1.0 + nu) / 1000.0 * ((1.0 - 2.0 * nu) * 0.01 / 3.0 + 0.04 / 3.0);
        let standard = (cylinder(2, 8, FiniteHexFormulation::Standard, cx) / exact - 1.0).abs();
        let coarse = (cylinder(2, 8, FiniteHexFormulation::MeanDilatation, cx) / exact - 1.0).abs();
        let fine = (cylinder(4, 16, FiniteHexFormulation::MeanDilatation, cx) / exact - 1.0).abs();
        eprintln!(
            "finite Lamé displacement relative errors: standard {standard:.9e}, mean coarse {coarse:.9e}, mean fine {fine:.9e}"
        );
        assert!(
            standard > 0.2,
            "fixture did not excite volumetric locking: {standard}"
        );
        assert!(coarse < 0.04, "coarse mean-dilatation error {coarse}");
        assert!(
            fine < 0.007 && fine < 0.35 * coarse,
            "refinement: coarse {coarse}, fine {fine}"
        );
    });
}

fn axial_supports(nodes: &[[f64; 3]], log_axial: f64) -> Vec<(usize, f64)> {
    let mut result = Vec::new();
    for (node, x) in nodes.iter().enumerate() {
        if x[0] == 0.0 {
            result.push((3 * node, 0.0));
        }
        if x[0] == 1.0 {
            result.push((3 * node, log_axial.exp_m1()));
        }
        if x[1] == 0.0 {
            result.push((3 * node + 1, 0.0));
        }
        if x[2] == 0.0 {
            result.push((3 * node + 2, 0.0));
        }
    }
    result
}

#[test]
fn finite_yield_unload_and_reversal_retain_eight_histories_with_free_lateral_motion() {
    with_cx(|cx| {
        let (nodes, cells) = mesh([1, 1, 1], |x, y, z| [x, y, z]);
        let (mu, bulk, yield_pa, hardening) = (8000.0, 16000.0, 200.0, 4000.0);
        let compliance = 1.0 / (3.0 * mu) + 1.0 / (9.0 * bulk);
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            let mut solid = FiniteThermoplasticHexSolid::try_new(
                &nodes,
                &cells,
                vec![family(mu, bulk, yield_pa, hardening)],
                &[300.0],
                formulation,
                HexAssemblyBudget::default(),
                cx,
            )
            .unwrap();
            let (mut signed_p, mut alpha, mut old_a, mut storage) = (0.0, 0.0, 0.0_f64, 0.0);
            for a in [0.004, 0.14, 0.13, -0.10, 0.03] {
                let trial = (a - signed_p) / compliance;
                let dp = (trial.abs() - yield_pa - hardening * alpha).max(0.0)
                    / (1.0 / compliance + hardening);
                signed_p += trial.signum() * dp;
                alpha += dp;
                let tau = (a - signed_p) / compliance;
                let b = tau / (9.0 * bulk) - tau / (6.0 * mu) - 0.5 * signed_p;
                let want_storage = 0.5 * compliance * tau * tau + 0.5 * hardening * alpha * alpha;
                let saved = solid.quadrature_states().to_vec();
                let step = solid
                    .advance(
                        &[300.0],
                        &axial_supports(&nodes, a),
                        &[],
                        0.2,
                        settings(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                for (x, u) in nodes.iter().zip(step.solid.displacement_m()) {
                    for axis in 0..3 {
                        near(
                            u[axis],
                            x[axis] * (if axis == 0 { a } else { b }).exp_m1(),
                            2e-9,
                            2e-9,
                        );
                    }
                }
                for (history, p) in step
                    .solid
                    .quadrature_states()
                    .iter()
                    .zip(step.solid.quadrature_first_piola_pa())
                {
                    near(history.alpha, alpha, 1e-10, 2e-9);
                    assert_eq!(history.creep_alpha, 0.0);
                    near(
                        history.inverse_plastic_cauchy_green[0],
                        (-2.0 * signed_p).exp(),
                        2e-10,
                        2e-9,
                    );
                    for i in [1, 2] {
                        near(
                            history.inverse_plastic_cauchy_green[i],
                            signed_p.exp(),
                            2e-10,
                            2e-9,
                        );
                    }
                    near(p[0], tau / a.exp(), 2e-6, 2e-9);
                    near(p[4], 0.0, 2e-6, 0.0);
                    near(p[8], 0.0, 2e-6, 0.0);
                }
                near(step.solid.free_energy_j(), want_storage, 2e-6, 2e-9);
                near(step.energy.plastic_dissipation_j, yield_pa * dp, 2e-7, 2e-9);
                near(step.energy.creep_dissipation_j, 0.0, 0.0, 0.0);
                near(step.energy.caloric_heat_j, yield_pa * dp, 2e-7, 2e-9);
                let work = tau / a.exp() * (a.exp() - old_a.exp());
                near(step.energy.support_endpoint_work_j, work, 2e-6, 2e-9);
                near(
                    step.energy.endpoint_work_defect_j,
                    work - (want_storage - storage) - yield_pa * dp,
                    3e-6,
                    2e-9,
                );
                near(step.energy.projection_work_defect_j, 0.0, 2e-7, 0.0);
                assert_eq!(solid.quadrature_states(), saved);
                solid = step.solid;
                old_a = a;
                storage = want_storage;
            }
            assert!(alpha > 0.25 && signed_p < 0.0);
        }
    });
}

#[test]
fn nonuniform_projected_creep_yield_heating_keeps_physical_loss_and_projection_work_separate() {
    with_cx(|cx| {
        let (nodes, cells) = distorted();
        let (mu, bulk, yield_pa, hardening, eta, activation) =
            (8000.0, 16000.0, 200.0, 4000.0, 80000.0, 3000.0);
        let law = family(mu, bulk, yield_pa, hardening)
            .with_norton_creep(
                FiniteJ2ThermalOverstress {
                    reference_viscosity_pa_s: eta,
                    stress_scale_pa: 400.0,
                    exponent: 1.0,
                    activation_temperature_k: activation,
                },
                3.0,
            )
            .unwrap();
        let mut solid = FiniteThermoplasticHexSolid::try_new(
            &nodes,
            &cells,
            vec![law.clone()],
            &[300.0],
            FiniteHexFormulation::MeanDilatation,
            HexAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let mut saw_projection_defect = false;
        for (index, (scale, temperature, dt)) in
            [(0.7, 300.0, 0.1), (1.0, 330.0, 0.2), (1.0, 350.0, 0.4)]
                .into_iter()
                .enumerate()
        {
            let u = motion(&nodes, scale);
            let old = solid.quadrature_states().to_vec();
            let old_u = solid.displacement_m().to_vec();
            let old_f = solid.quadrature_assumed_deformation().to_vec();
            let old_free = solid.free_energy_j();
            let probe = solid
                .evaluate_trial(&[temperature], &u, dt, false, cx)
                .unwrap();
            let nodal_work: f64 = probe
                .internal_force_n
                .iter()
                .zip(u.iter().zip(&old_u))
                .map(|(force, (new, old))| {
                    (0..3).map(|i| force[i] * (new[i] - old[i])).sum::<f64>()
                })
                .sum();
            if index == 1 {
                let mut rejected = controls();
                rejected.maximum_endpoint_work_defect_j = 0.0;
                assert!(
                    solid
                        .advance(
                            &[temperature],
                            &prescribed(&u),
                            &[],
                            dt,
                            settings(),
                            rejected,
                            cx
                        )
                        .is_err()
                );
                assert_eq!(solid.quadrature_states(), old);
                assert_eq!(solid.displacement_m(), old_u);
                assert_eq!(solid.temperatures_k(), &[300.0]);
            }
            let step = solid
                .advance(
                    &[temperature],
                    &prescribed(&u),
                    &[],
                    dt,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            let viscosity = eta * (activation * (1.0 / temperature - 1.0 / 300.0)).exp();
            let volume = step.solid.reference_volume_m3(0).unwrap();
            let mut mean_j = 0.0;
            for q in 0..8 {
                mean_j += step.solid.quadrature_volume_m3(0, q).unwrap()
                    * determinant(&step.solid.quadrature_deformation()[q])
                    / volume;
            }
            let (mut free, mut plastic_loss, mut creep_loss, mut coordinate_work) =
                (0.0, 0.0, 0.0, 0.0);
            for q in 0..8 {
                let state = &step.solid.quadrature_states()[q];
                let f = &step.solid.quadrature_assumed_deformation()[q];
                let p = &step.solid.quadrature_first_piola_pa()[q];
                let weight = step.solid.quadrature_volume_m3(0, q).unwrap();
                let tau = kirchhoff(p, f);
                let stress = equivalent(&tau);
                let dp = state.alpha - old[q].alpha;
                let dc = state.creep_alpha - old[q].creep_alpha;
                near(determinant(f), mean_j, 2e-13, 2e-13);
                near(dc, dt * stress / viscosity, 1e-11, 2e-9);
                assert!(dc > 0.0);
                if index == 2 {
                    assert_eq!(state.alpha, old[q].alpha);
                } else {
                    assert!(dp > 0.0);
                    near(stress, yield_pa + hardening * state.alpha, 2e-6, 2e-9);
                }
                let local_free = stress * stress / (6.0 * mu)
                    + 0.5 * bulk * mean_j.ln().powi(2)
                    + 0.5 * hardening * state.alpha.powi(2);
                free += weight * local_free;
                plastic_loss += weight * yield_pa * dp;
                creep_loss += weight * stress * dc;
                coordinate_work += weight * dot(p, &std::array::from_fn(|i| f[i] - old_f[q][i]));
                let account = &step.quadrature_energy[q];
                near(account.free_after_j_m3, local_free, 2e-6, 2e-9);
                near(account.plastic_dissipation_j_m3, yield_pa * dp, 2e-7, 2e-9);
                near(account.creep_dissipation_j_m3, stress * dc, 2e-7, 2e-9);
            }
            near(step.solid.free_energy_j(), free, 3e-6, 2e-9);
            near(step.solid.internal_energy_j(), free, 3e-6, 2e-9);
            near(step.energy.plastic_dissipation_j, plastic_loss, 2e-7, 2e-9);
            near(step.energy.creep_dissipation_j, creep_loss, 2e-7, 2e-9);
            near(
                step.energy.caloric_heat_j,
                plastic_loss + creep_loss,
                3e-7,
                2e-9,
            );
            near(step.energy.reversible_heat_j, 0.0, 1e-9, 0.0);
            near(
                step.energy.incremental_potential_j,
                free + plastic_loss + 0.5 * creep_loss,
                3e-6,
                2e-9,
            );
            near(step.energy.internal_endpoint_work_j, nodal_work, 2e-6, 2e-9);
            near(
                step.energy.projection_work_defect_j,
                nodal_work - coordinate_work,
                2e-6,
                2e-9,
            );
            near(
                step.energy.material_endpoint_work_defect_j,
                coordinate_work - (free - old_free) - plastic_loss - creep_loss,
                3e-6,
                2e-9,
            );
            near(
                step.energy.endpoint_work_defect_j,
                nodal_work - (free - old_free) - plastic_loss - creep_loss,
                4e-6,
                2e-9,
            );
            near(step.energy.balance_error_j, 0.0, 3e-6, 0.0);
            saw_projection_defect |= step.energy.projection_work_defect_j.abs() > 1e-4;
            assert_eq!(step.solid.families()[0], law);
            assert_eq!(solid.quadrature_states(), old);
            solid = step.solid;
        }
        assert!(saw_projection_defect);
        assert!(
            solid
                .quadrature_states()
                .iter()
                .all(|s| s.alpha > 0.03 && s.creep_alpha > 0.0)
        );
    });
}
