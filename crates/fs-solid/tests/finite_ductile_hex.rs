//! G1/G3: native progressive ductility at all eight original H8 material points.
//! Synthetic coefficients test mechanics; the explicit length is not a claim
//! of mesh-objective fracture or an experimental alloy calibration.
#![deny(unsafe_code)]

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::ductile::evolution::DuctileJ2Settings;
use fs_material::ductile::{DuctileDomain, TabulatedDuctileLocus};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::{FiniteDuctileError, FiniteDuctileJ2};
use fs_material::plastic::finite::rate::FiniteOverstressParameters;
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::ductile::{DuctileStepConditions, DuctileTemperatures};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::hex::FiniteHexFormulation;
use fs_solid::hyper3::plastic::hex::ductile::{FiniteDuctileHexError, FiniteDuctileHexSolid};
use fs_solid::plastic3::hex::HexAssemblyBudget;
use fs_solid::plastic3::{IncrementSettings, SmallStrainTetError};

const NODES: [[f64; 3]; 8] = [
    [0.0, 0.0, 0.0],
    [1.0, 0.0, 0.0],
    [1.0, 1.0, 0.0],
    [0.0, 1.0, 0.0],
    [0.0, 0.0, 1.0],
    [1.0, 0.0, 1.0],
    [1.0, 1.0, 1.0],
    [0.0, 1.0, 1.0],
];
const HEXES: [[usize; 8]; 1] = [[0, 1, 2, 3, 4, 5, 6, 7]];
const MU: f64 = 500.0;
const K: f64 = 1000.0;
const Y: f64 = 3.0;
const H: f64 = 40.0;
const LENGTH: f64 = 0.02;
const GP: f64 = 0.004;
const ONSET: f64 = 0.002;
const DT: f64 = 0.2;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 67,
                kernel_id: 199,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&gate, &cx)
    })
}

fn model(onset: f64, damage_limit: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    let native = FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: K,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.6,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap();
    let locus = TabulatedDuctileLocus::new(
        vec![(-50.0, onset), (50.0, onset)],
        DuctileDomain {
            triaxiality: [-50.0, 50.0],
            strain_rate_s_inv: [1e-12, 100.0],
            temperature_k: [300.0, 300.0],
        },
    )
    .unwrap();
    FiniteDuctileJ2::new(
        native,
        locus,
        DuctileJ2Settings {
            characteristic_length_m: LENGTH,
            post_initiation_plastic_work_j_m2: GP,
            maximum_damage_increment: damage_limit,
            energy_tolerances: EnergyTolerances {
                energy_j_m3: 1e-9,
                yield_pa: 1e-8,
            },
        },
    )
    .unwrap()
}

fn conditions() -> DuctileStepConditions<'static> {
    DuctileStepConditions {
        duration_s: DT,
        temperatures_k: DuctileTemperatures::Uniform(300.0),
    }
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1e-8,
        maximum_endpoint_work_defect_j: 1000.0,
        yield_tolerance_pa: 1e-7,
    }
}

fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 + 2e-8 * actual.abs().max(expected.abs()),
        "actual={actual:.16e}, expected={expected:.16e}"
    );
}

fn integrity(alpha: f64, onset: f64) -> f64 {
    let post = (alpha - onset).max(0.0);
    (-LENGTH * post * (Y + H * (onset + 0.5 * post)) / GP).exp()
}

fn coaxial(e: f64, dilation: f64, angle: f64) -> [f64; 9] {
    let axial = (dilation + e).exp();
    let lateral = (dilation - 0.5 * e).exp();
    let (s, c) = angle.sin_cos();
    [
        c * axial,
        -s * lateral,
        0.0,
        s * axial,
        c * lateral,
        0.0,
        0.0,
        0.0,
        lateral,
    ]
}

fn affine(f: &[f64; 9]) -> Vec<[f64; 3]> {
    NODES
        .iter()
        .map(|x| core::array::from_fn(|i| (0..3).map(|j| f[3 * i + j] * x[j]).sum::<f64>() - x[i]))
        .collect()
}

fn motion(scale: f64) -> Vec<[f64; 3]> {
    NODES
        .iter()
        .map(|&[x, y, z]| {
            [
                scale * (0.06 * x + 0.02 * x * y + 0.01 * y * z),
                scale * (-0.015 * y + 0.008 * x * z),
                scale * (-0.01 * z + 0.012 * x * y),
            ]
        })
        .collect()
}

fn prescribed(u: &[[f64; 3]]) -> Vec<(usize, f64)> {
    u.as_flattened().iter().copied().enumerate().collect()
}

fn determinant(f: &[f64; 9]) -> f64 {
    f[0] * (f[4] * f[8] - f[5] * f[7]) - f[1] * (f[3] * f[8] - f[5] * f[6])
        + f[2] * (f[3] * f[7] - f[4] * f[6])
}

fn dot(a: &[f64; 9], b: &[f64; 9]) -> f64 {
    (0..9).map(|i| a[i] * b[i]).sum()
}

// Independent positive-volume storage reconstructed from physical Kirchhoff
// stress. The native resistance is energetic linear hardening; damage scales
// this whole storage on the tensile volumetric branch used by these fixtures.
fn storage(p: &[f64; 9], f: &[f64; 9], alpha: f64, g: f64) -> f64 {
    let tau: [f64; 9] = core::array::from_fn(|ij| {
        (0..3)
            .map(|k| p[3 * (ij / 3) + k] * f[3 * (ij % 3) + k])
            .sum::<f64>()
            / g
    });
    let pressure = (tau[0] + tau[4] + tau[8]) / 3.0;
    let dev_squared: f64 = (0..9)
        .map(|i| (tau[i] - if i % 4 == 0 { pressure } else { 0.0 }).powi(2))
        .sum();
    g * (dev_squared / (4.0 * MU) + pressure * pressure / (2.0 * K) + 0.5 * H * alpha * alpha)
}

#[test]
fn affine_damage_and_physical_rate_match_independent_oracles_in_both_formulations() {
    with_cx(|_, cx| {
        let e = 0.04;
        let dilation = 0.001;
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            for viscosity in [0.0, 400.0] {
                let material = model(ONSET, 1.0);
                let material = if viscosity == 0.0 {
                    material
                } else {
                    material
                        .with_overstress(FiniteOverstressParameters {
                            viscosity_pa_s: viscosity,
                            stress_scale_pa: Y,
                            exponent: 1.0,
                        })
                        .unwrap()
                };
                let materials = [&material; 8];
                let initial = FiniteDuctileHexSolid::try_new(
                    &NODES,
                    &HEXES,
                    &materials,
                    formulation,
                    HexAssemblyBudget::default(),
                    cx,
                )
                .unwrap();
                let old = initial.quadrature_states().to_vec();
                let alpha = (3.0 * MU * e - Y) / (3.0 * MU + H + viscosity / DT);
                let g = integrity(alpha, ONSET);
                let tau = [
                    3.0 * K * dilation + 2.0 * MU * (e - alpha),
                    3.0 * K * dilation - MU * (e - alpha),
                    3.0 * K * dilation - MU * (e - alpha),
                ];
                let native_storage = 0.5 * K * (3.0 * dilation).powi(2)
                    + 1.5 * MU * (e - alpha).powi(2)
                    + 0.5 * H * alpha * alpha;
                for angle in [0.0_f64, 0.47] {
                    let f = coaxial(e, dilation, angle);
                    let supports = prescribed(&affine(&f));
                    let step = initial
                        .advance(
                            &supports,
                            &[],
                            conditions(),
                            IncrementSettings::default(),
                            controls(),
                            cx,
                        )
                        .unwrap();
                    let (s, c) = angle.sin_cos();
                    let px = g * tau[0] / (dilation + e).exp();
                    let py = g * tau[1] / (dilation - 0.5 * e).exp();
                    let expected_p = [c * px, -s * py, 0.0, s * px, c * py, 0.0, 0.0, 0.0, py];
                    for q in 0..8 {
                        let state = &step.solid.quadrature_states()[q];
                        near(state.native_state().alpha, alpha);
                        near(state.integrity(), g);
                        assert!(state.damage() > 0.1);
                        for i in 0..9 {
                            near(step.solid.quadrature_deformation()[q][i], f[i]);
                            near(step.solid.quadrature_assumed_deformation()[q][i], f[i]);
                            near(step.solid.quadrature_first_piola_pa()[q][i], expected_p[i]);
                        }
                        near(step.solid.quadrature_volume_m3(0, q).unwrap(), 0.125);
                    }
                    near(step.solid.reference_volume_m3(0).unwrap(), 1.0);
                    assert_eq!(step.solid.formulation(), formulation);
                    near(step.solid.stored_energy_j(), g * native_storage);
                    near(
                        step.energy.plastic_dissipation_j,
                        g * (Y + viscosity * alpha / DT) * alpha,
                    );
                    near(step.energy.damage_release_j, 0.0);
                    let delta = core::array::from_fn(|i| f[i] - if i % 4 == 0 { 1.0 } else { 0.0 });
                    let work = dot(&expected_p, &delta);
                    near(step.energy.material_endpoint_work_j, work);
                    near(step.energy.internal_endpoint_work_j, work);
                    near(step.energy.support_endpoint_work_j, work);
                    near(step.energy.projection_work_defect_j, 0.0);
                    near(step.energy.balance_error_j, 0.0);
                    if angle == 0.0 {
                        let retained = step.solid.quadrature_states().to_vec();
                        let held = step
                            .solid
                            .advance(
                                &supports,
                                &[],
                                conditions(),
                                IncrementSettings::default(),
                                controls(),
                                cx,
                            )
                            .unwrap();
                        if viscosity == 0.0 {
                            assert_eq!(held.solid.quadrature_states(), retained);
                            near(held.energy.plastic_dissipation_j, 0.0);
                            near(held.energy.damage_release_j, 0.0);
                        } else {
                            let relaxed = (3.0 * MU * e - Y + viscosity * alpha / DT)
                                / (3.0 * MU + H + viscosity / DT);
                            for state in held.solid.quadrature_states() {
                                near(state.native_state().alpha, relaxed);
                                near(state.integrity(), integrity(relaxed, ONSET));
                                assert!(state.damage() > retained[0].damage());
                            }
                            assert!(held.energy.plastic_dissipation_j > 0.0);
                            assert!(held.energy.damage_release_j > 0.0);
                        }
                        near(held.energy.internal_endpoint_work_j, 0.0);
                        near(held.energy.balance_error_j, 0.0);
                        assert_eq!(step.solid.quadrature_states(), retained);
                    }
                }
                assert_eq!(initial.quadrature_states(), old);
                assert_eq!(initial.displacement_m(), &[[0.0; 3]; 8]);
            }
        }
    });
}

#[test]
fn nonaffine_damaged_tangents_retain_every_material_and_projection_derivative() {
    with_cx(|_, cx| {
        let material = model(ONSET, 1.0);
        let materials = [&material; 8];
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            let initial = FiniteDuctileHexSolid::try_new(
                &NODES,
                &HEXES,
                &materials,
                formulation,
                HexAssemblyBudget::default(),
                cx,
            )
            .unwrap();
            let solid = initial
                .advance(
                    &prescribed(&motion(0.4)),
                    &[],
                    conditions(),
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap()
                .solid;
            let old_states = solid.quadrature_states().to_vec();
            let old_u = solid.displacement_m().to_vec();
            assert!(old_states.iter().all(|s| s.damage() > 0.0));
            let trial = motion(1.0);
            let response = solid
                .evaluate_trial(&trial, conditions(), true, cx)
                .unwrap();
            let replay = solid
                .evaluate_trial(&trial, conditions(), false, cx)
                .unwrap();
            assert_eq!(response.internal_force_n, replay.internal_force_n);
            let tangent = response.tangent_n_m.unwrap();
            assert_eq!((tangent.nrows(), tangent.ncols()), (24, 24));
            let h = 1e-6;
            let mut error = 0.0_f64;
            let mut norm = 0.0_f64;
            let mut asymmetry = 0.0_f64;
            for column in 0..24 {
                let mut plus = trial.clone();
                let mut minus = trial.clone();
                plus[column / 3][column % 3] += h;
                minus[column / 3][column % 3] -= h;
                let plus = solid
                    .evaluate_trial(&plus, conditions(), false, cx)
                    .unwrap();
                let minus = solid
                    .evaluate_trial(&minus, conditions(), false, cx)
                    .unwrap();
                for row in 0..24 {
                    let numerical = (plus.internal_force_n[row / 3][row % 3]
                        - minus.internal_force_n[row / 3][row % 3])
                        / (2.0 * h);
                    let analytical = tangent.get(row, column);
                    error = error.hypot(numerical - analytical);
                    norm = norm.hypot(analytical);
                    asymmetry = asymmetry.max((analytical - tangent.get(column, row)).abs());
                }
            }
            assert!(
                norm > 0.0 && error / norm < 2e-5,
                "relative full tangent error {}",
                error / norm
            );
            assert!(
                asymmetry > 1e-6,
                "damage force must retain its nonsymmetric derivative"
            );
            let shifted: Vec<_> = trial
                .iter()
                .map(|u| [u[0] + 2.0, u[1] - 3.0, u[2] + 4.0])
                .collect();
            let translated = solid
                .evaluate_trial(&shifted, conditions(), false, cx)
                .unwrap();
            for (&a, &b) in response
                .internal_force_n
                .as_flattened()
                .iter()
                .zip(translated.internal_force_n.as_flattened())
            {
                near(a, b);
            }
            assert_eq!(solid.quadrature_states(), old_states);
            assert_eq!(solid.displacement_m(), old_u);
        }
    });
}

#[test]
fn original_gauss_points_keep_distinct_calibrations_and_unaveraged_damage() {
    with_cx(|_, cx| {
        let onsets: [f64; 8] = core::array::from_fn(|q| {
            if q % 2 == 0 {
                ONSET + 0.001 * q as f64
            } else {
                0.2
            }
        });
        let models: Vec<_> = onsets.iter().map(|&onset| model(onset, 1.0)).collect();
        let materials: Vec<_> = models.iter().collect();
        let e = 0.04;
        let alpha = (3.0 * MU * e - Y) / (3.0 * MU + H);
        assert!(onsets.iter().sum::<f64>() / 8.0 > alpha);
        let f = coaxial(e, 0.001, 0.0);
        let temperatures = [300.0; 8];
        let at = DuctileStepConditions {
            temperatures_k: DuctileTemperatures::PerPoint(&temperatures),
            ..conditions()
        };
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            let initial = FiniteDuctileHexSolid::try_new(
                &NODES,
                &HEXES,
                &materials,
                formulation,
                HexAssemblyBudget::default(),
                cx,
            )
            .unwrap();
            let step = initial
                .advance(
                    &prescribed(&affine(&f)),
                    &[],
                    at,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let mut damaged = 0;
            for q in 0..8 {
                let state = &step.solid.quadrature_states()[q];
                assert!(std::ptr::eq(step.solid.materials()[q], &models[q]));
                assert!(models[q].owns_state(state));
                near(state.native_state().alpha, alpha);
                near(state.integrity(), integrity(alpha, onsets[q]));
                if q % 2 == 0 {
                    assert!(state.damage() > 0.0);
                    damaged += 1;
                } else {
                    assert_eq!(state.damage(), 0.0);
                }
            }
            assert_eq!(
                damaged, 4,
                "an averaged onset would miss all four initiated points"
            );
            assert!(
                initial
                    .quadrature_states()
                    .iter()
                    .all(|s| s.damage() == 0.0)
            );
        }
    });
}

#[test]
fn nonuniform_projection_defect_stays_out_of_physical_plastic_and_damage_work() {
    with_cx(|_, cx| {
        let material = model(ONSET, 1.0);
        let materials = [&material; 8];
        let temperatures = [300.0];
        let at = DuctileStepConditions {
            temperatures_k: DuctileTemperatures::PerHexElement(&temperatures),
            ..conditions()
        };
        for formulation in [
            FiniteHexFormulation::Standard,
            FiniteHexFormulation::MeanDilatation,
        ] {
            let initial = FiniteDuctileHexSolid::try_new(
                &NODES,
                &HEXES,
                &materials,
                formulation,
                HexAssemblyBudget::default(),
                cx,
            )
            .unwrap();
            let old = initial
                .advance(
                    &prescribed(&motion(0.4)),
                    &[],
                    at,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap()
                .solid;
            let retained = old.quadrature_states().to_vec();
            let target = motion(1.0);
            let probe = old.evaluate_trial(&target, at, false, cx).unwrap();
            let nodal_work: f64 = probe
                .internal_force_n
                .iter()
                .zip(target.iter().zip(old.displacement_m()))
                .map(|(force, (new, before))| {
                    (0..3).map(|i| force[i] * (new[i] - before[i])).sum::<f64>()
                })
                .sum();
            let step = old
                .advance(
                    &prescribed(&target),
                    &[],
                    at,
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            let volume = step.solid.reference_volume_m3(0).unwrap();
            let mean_j: f64 = (0..8)
                .map(|q| {
                    step.solid.quadrature_volume_m3(0, q).unwrap()
                        * determinant(&step.solid.quadrature_deformation()[q])
                        / volume
                })
                .sum();
            let (
                mut stored,
                mut plastic,
                mut released,
                mut material_work,
                mut material_defect,
                mut absolute_material,
            ) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
            for q in 0..8 {
                let state = &step.solid.quadrature_states()[q];
                let before = &retained[q];
                let f = &step.solid.quadrature_assumed_deformation()[q];
                let physical = &step.solid.quadrature_deformation()[q];
                let p = &step.solid.quadrature_first_piola_pa()[q];
                let weight = step.solid.quadrature_volume_m3(0, q).unwrap();
                let physical_j = determinant(physical);
                assert!(physical_j > 1.0 && determinant(f) > 1.0);
                let scale = if formulation == FiniteHexFormulation::MeanDilatation {
                    near(determinant(f), mean_j);
                    (mean_j / physical_j).cbrt()
                } else {
                    1.0
                };
                for i in 0..9 {
                    near(f[i], scale * physical[i]);
                }
                assert_eq!(state.deformation(), f);
                assert_eq!(state.first_piola_pa(), p);
                let alpha = state.native_state().alpha;
                let g = integrity(alpha, ONSET);
                near(state.integrity(), g);
                let local_storage = storage(p, f, alpha, g);
                let before_storage = storage(
                    before.first_piola_pa(),
                    before.deformation(),
                    before.native_state().alpha,
                    before.integrity(),
                );
                let local_plastic = g * Y * (alpha - before.native_state().alpha);
                let local_release = (before.integrity() - g) * before_storage / before.integrity();
                let delta = core::array::from_fn(|i| f[i] - before.deformation()[i]);
                let work = dot(p, &delta);
                let defect =
                    work - (local_storage - before_storage) - local_plastic - local_release;
                near(state.stored_energy_j_m3(), local_storage);
                near(
                    step.quadrature_energy[q].plastic_dissipation_j_m3,
                    local_plastic,
                );
                near(step.quadrature_energy[q].damage_release_j_m3, local_release);
                near(step.quadrature_energy[q].endpoint_work_defect_j_m3, defect);
                stored += weight * local_storage;
                plastic += weight * local_plastic;
                released += weight * local_release;
                material_work += weight * work;
                material_defect += weight * defect;
                absolute_material += weight * defect.abs();
            }
            let projection = nodal_work - material_work;
            near(step.energy.stored_before_j, old.stored_energy_j());
            near(step.energy.stored_after_j, stored);
            near(step.solid.stored_energy_j(), stored);
            near(step.energy.plastic_dissipation_j, plastic);
            near(step.energy.damage_release_j, released);
            assert!(plastic > 0.0 && released > 0.0);
            near(step.energy.material_endpoint_work_j, material_work);
            near(step.energy.internal_endpoint_work_j, nodal_work);
            near(step.energy.support_endpoint_work_j, nodal_work);
            near(step.energy.material_endpoint_work_defect_j, material_defect);
            near(step.energy.projection_work_defect_j, projection);
            near(
                step.energy.endpoint_work_defect_j,
                material_defect + projection,
            );
            near(
                step.energy.absolute_material_endpoint_work_defect_j,
                absolute_material,
            );
            near(
                step.energy.absolute_projection_work_defect_j,
                projection.abs(),
            );
            near(
                step.energy.absolute_endpoint_work_defect_j,
                absolute_material + projection.abs(),
            );
            near(step.energy.equilibrium_work_residual_j, 0.0);
            near(step.energy.balance_error_j, 0.0);
            if formulation == FiniteHexFormulation::MeanDilatation {
                assert!(
                    projection.abs() > 1e-8,
                    "fixture must exercise actual projection work"
                );
            } else {
                near(projection, 0.0);
            }
            let strict = FinitePlasticStepControls {
                maximum_endpoint_work_defect_j: 0.5 * step.energy.absolute_endpoint_work_defect_j,
                ..controls()
            };
            assert!(matches!(
                old.advance(
                    &prescribed(&target),
                    &[],
                    at,
                    IncrementSettings::default(),
                    strict,
                    cx,
                ),
                Err(FiniteDuctileHexError::Energy { .. })
            ));
            assert_eq!(old.quadrature_states(), retained);
            assert_eq!(old.displacement_m(), motion(0.4));
        }
    });
}

#[test]
fn point_damage_limits_hidden_corner_folds_and_cancellation_are_transactional() {
    with_cx(|gate, cx| {
        let material = model(ONSET, 1.0);
        let limited = model(ONSET, 0.01);
        let incomplete = [&material; 7];
        assert!(
            FiniteDuctileHexSolid::try_new(
                &NODES,
                &HEXES,
                &incomplete,
                FiniteHexFormulation::Standard,
                HexAssemblyBudget::default(),
                cx,
            )
            .is_err()
        );
        let point_materials = [
            &material, &material, &material, &material, &material, &material, &material, &limited,
        ];
        let bounded = FiniteDuctileHexSolid::try_new(
            &NODES,
            &HEXES,
            &point_materials,
            FiniteHexFormulation::MeanDilatation,
            HexAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let virgin = bounded.quadrature_states().to_vec();
        let target = affine(&coaxial(0.04, 0.001, 0.0));
        assert!(matches!(
            bounded.evaluate_trial(&target, conditions(), true, cx),
            Err(FiniteDuctileHexError::Material {
                element: 0,
                point: 7,
                source: FiniteDuctileError::DamageIncrement { .. },
            })
        ));
        assert_eq!(bounded.quadrature_states(), virgin);
        assert_eq!(bounded.displacement_m(), &[[0.0; 3]; 8]);
        let materials = [&material; 8];
        let initial = FiniteDuctileHexSolid::try_new(
            &NODES,
            &HEXES,
            &materials,
            FiniteHexFormulation::MeanDilatation,
            HexAssemblyBudget::default(),
            cx,
        )
        .unwrap();
        let solid = initial
            .advance(
                &prescribed(&motion(0.4)),
                &[],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            )
            .unwrap()
            .solid;
        let old_states = solid.quadrature_states().to_vec();
        let old_u = solid.displacement_m().to_vec();
        let old_f = solid.quadrature_deformation().to_vec();
        let old_assumed = solid.quadrature_assumed_deformation().to_vec();
        let old_storage = solid.stored_energy_j();
        let short_temperatures = [300.0; 7];
        let incomplete = DuctileStepConditions {
            temperatures_k: DuctileTemperatures::PerPoint(&short_temperatures),
            ..conditions()
        };
        assert!(
            solid
                .evaluate_trial(&target, incomplete, false, cx)
                .is_err()
        );
        assert!(
            solid
                .advance(
                    &prescribed(&target),
                    &[],
                    DuctileStepConditions {
                        duration_s: 0.0,
                        ..conditions()
                    },
                    IncrementSettings::default(),
                    controls(),
                    cx,
                )
                .is_err()
        );
        let mut folded = vec![[0.0; 3]; 8];
        folded[6][0] = -1.3;
        // This trilinear displacement gives detF=1-1.3*y*z. All eight
        // Gauss values exceed 0.19 and the mean is 0.675, yet a corner is -0.3.
        let high_gauss = 0.5 * (1.0 + 1.0 / 3.0_f64.sqrt());
        assert!(1.0 - 1.3 * high_gauss * high_gauss > 0.1);
        assert!(
            solid
                .evaluate_trial(&folded, conditions(), false, cx)
                .is_err()
        );
        gate.request();
        assert!(matches!(
            solid.advance(
                &prescribed(&target),
                &[],
                conditions(),
                IncrementSettings::default(),
                controls(),
                cx,
            ),
            Err(FiniteDuctileHexError::Solid(SmallStrainTetError::Cancelled))
        ));
        assert_eq!(solid.quadrature_states(), old_states);
        assert_eq!(solid.displacement_m(), old_u);
        assert_eq!(solid.quadrature_deformation(), old_f);
        assert_eq!(solid.quadrature_assumed_deformation(), old_assumed);
        assert_eq!(solid.stored_energy_j(), old_storage);
    });
}
