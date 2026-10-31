//! Independent ramp/hold, laminate and work oracles; synthetic numerical laws.
use fs_material::Tangent6;
use fs_material::laminate::thermal::LaminaThermalStrain;
use fs_material::laminate::thermal::history::{LaminateHistoryError, LaminateThermalHistory};
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use fs_material::visco_tensor::{TensorProny, TensorPronyBranch};
use std::sync::Arc;

const E: f64 = 100.0;
const BRANCH: f64 = 200.0;
const NU: f64 = 0.25;
const VISCOSITY: f64 = 6.0;

fn isotropic(e: f64, nu: f64) -> Tangent6 {
    let g = e / (2.0 * (1.0 + nu));
    let lambda = e * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    core::array::from_fn(|i| {
        core::array::from_fn(|j| {
            if i < 3 && j < 3 {
                lambda + if i == j { 2.0 * g } else { 0.0 }
            } else if i == j {
                2.0 * g
            } else {
                0.0
            }
        })
    })
}
fn law(branch_nu: f64, tau: f64, kelvin: bool) -> TensorProny {
    let law = TensorProny::new_mandel(
        isotropic(E, NU),
        vec![TensorPronyBranch::from_stiffness_mandel(isotropic(BRANCH, branch_nu), tau).unwrap()],
        0.05,
    )
    .unwrap();
    if kelvin {
        law.with_kelvin_viscosity_mandel(isotropic(VISCOSITY, NU))
            .unwrap()
    } else {
        law
    }
}
fn layer(law: TensorProny, h: f64, angle: f64) -> LaminaLayer {
    LaminaLayer {
        thickness_m: h,
        orientation_rad: angle,
        density_kg_m3: 1200.0,
        material: None,
        law: LayerLaw::Prony(law),
    }
}
fn history(layers: Vec<LaminaLayer>) -> LaminateThermalHistory {
    let section = Arc::new(LaminateSection::try_new(layers, ReferenceSurface::Midsurface).unwrap());
    let state = LaminateThermalHistory::try_new(Arc::clone(&section), || false).unwrap();
    assert!(core::ptr::eq(section.as_ref(), state.section()));
    state
}
fn free(e: f64) -> LaminaThermalStrain {
    LaminaThermalStrain::uniform([e, e, 0.0]).unwrap()
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_thermal_ramp_and_hold_keep_actual_branch_storage_loss_and_separate_work_ports() {
    for kelvin in [false, true] {
        for imposed in [0.0, 0.001] {
            let (h, tau, dt, thermal) = (0.03, 2.0, 0.5, 0.004);
            let mut body = history(vec![layer(law(NU, tau, kelvin), h, 0.0)]);
            let x = imposed - thermal;
            let r: f64 = dt / tau;
            let decay = (-r).exp();
            let phi = (1.0 - decay) / r;
            let eq = E / (1.0 - NU);
            let branch = BRANCH / (1.0 - NU);
            let viscous = if kelvin { VISCOSITY / (1.0 - NU) } else { 0.0 };
            let q = phi * x;
            let mean = 0.5 * eq * x + branch * x * (1.0 - phi) / r + viscous * x / dt;
            let stress = eq * x + branch * q + viscous * x / dt;
            let storage = h * (eq * x * x + branch * q * q);
            // Direct integral of the two equal biaxial Maxwell channels, plus
            // Kelvin loss. This does not define heat by subtracting work/storage.
            let loss = 2.0
                * h
                * branch
                * (x / r).powi(2)
                * (r - 2.0 * (1.0 - decay) + 0.5 * (1.0 - decay * decay))
                + 2.0 * h * viscous * x * x / dt;
            let step = body
                .advance(
                    [imposed, imposed, 0.0],
                    [0.0; 3],
                    &[free(thermal)],
                    dt,
                    || false,
                )
                .unwrap();
            near(step.membrane_force_n_m[0], h * stress, 2e-15);
            near(step.membrane_force_n_m[1], h * stress, 2e-15);
            near(step.stored_energy_j_m2, storage, 2e-18);
            near(step.branch_stored_energy_j_m2, h * branch * q * q, 2e-18);
            near(step.mechanical_work_j_m2, 2.0 * h * mean * imposed, 2e-18);
            near(
                step.thermal_strain_work_j_m2,
                -2.0 * h * mean * thermal,
                3e-18,
            );
            near(step.dissipated_j_m2, loss, 4e-18);
            near(step.balance_residual_j_m2, 0.0, 4e-18);
            assert!(step.dissipated_j_m2 > 0.0);
            if imposed == 0.0 {
                assert_eq!(step.mechanical_work_j_m2, 0.0);
            } else {
                assert!(step.mechanical_work_j_m2 < 0.0);
            }
            let held_s: f64 = 3.0;
            let held_decay = (-held_s / tau).exp();
            let held = body
                .advance(
                    [imposed, imposed, 0.0],
                    [0.0; 3],
                    &[free(thermal)],
                    held_s,
                    || false,
                )
                .unwrap();
            near(
                held.membrane_force_n_m[0],
                h * (eq * x + branch * q * held_decay),
                2e-15,
            );
            near(
                held.branch_stored_energy_j_m2,
                h * branch * (q * held_decay).powi(2),
                1e-18,
            );
            near(
                held.dissipated_j_m2,
                h * branch * q * q * (1.0 - held_decay * held_decay),
                2e-18,
            );
            assert_eq!(held.mechanical_work_j_m2, 0.0);
            assert_eq!(held.thermal_strain_work_j_m2, 0.0);
            near(held.balance_residual_j_m2, 0.0, 4e-18);
            near(body.dissipated_j_m2(), loss + held.dissipated_j_m2, 4e-18);
            assert_eq!(body.time_s(), dt + held_s);
            near(body.section().i0_kg_m2(), 1200.0 * h, 0.0);
        }
    }
}

#[test]
fn g1_compatible_affine_expansion_remains_stress_free_with_original_layer_orientations() {
    let mut layers = vec![
        layer(law(NU, 0.2, false), 0.02, 0.37),
        layer(law(NU, 3.0, false), 0.03, -0.6),
    ];
    layers.push(LaminaLayer {
        thickness_m: 0.01,
        orientation_rad: 0.4,
        density_kg_m3: 1500.0,
        material: None,
        law: LayerLaw::PlaneStress {
            stiffness_pa: [200.0, 10.0, 0.0, 10.0, 100.0, 0.0, 0.0, 0.0, 30.0],
            transverse_shear_pa: [30.0; 2],
        },
    });
    let mut body = history(layers);
    let membrane = [0.002, -0.001, 0.001];
    let gradient = [0.01, 0.005, -0.003];
    let profiles: Vec<_> = body
        .section()
        .layers()
        .iter()
        .zip(body.section().layer_integrals())
        .map(|(layer, geometry)| {
            let (s, c) = layer.orientation_rad.sin_cos();
            let local = |z| {
                let e: [f64; 3] = core::array::from_fn(|j| membrane[j] + z * gradient[j]);
                [
                    c * c * e[0] + s * s * e[1] + c * s * e[2],
                    s * s * e[0] + c * c * e[1] - c * s * e[2],
                    -2.0 * c * s * (e[0] - e[1]) + (c * c - s * s) * e[2],
                ]
            };
            LaminaThermalStrain::from_faces(local(geometry.z_bottom_m), local(geometry.z_top_m))
                .unwrap()
        })
        .collect();
    for dt in [0.3, 1.0] {
        let step = body
            .advance(membrane, gradient, &profiles, dt, || false)
            .unwrap();
        near(step.stored_energy_j_m2, 0.0, 1e-31);
        near(step.dissipated_j_m2, 0.0, 1e-31);
        for value in step
            .membrane_force_n_m
            .into_iter()
            .chain(step.bending_moment_n)
        {
            near(value, 0.0, 3e-16);
        }
        for index in 0..3 {
            let stress = body.layer_stress(index).unwrap();
            for value in stress.bottom_pa.into_iter().chain(stress.top_pa) {
                near(value, 0.0, 2e-14);
            }
        }
    }
}

#[test]
fn g1_bent_mismatched_layers_keep_residual_memory_and_distinct_relaxation_times() {
    let (h, a, b) = (0.08, 0.001, 0.005);
    let delta = b - a;
    let mut body = history(vec![
        layer(law(NU, 1.0, false), h / 2.0, 0.0),
        layer(law(NU, 3.0, false), h / 2.0, 0.0),
    ]);
    let membrane = [(a + b) / 2.0, (a + b) / 2.0, 0.0];
    let gradient = [1.5 * delta / h, 1.5 * delta / h, 0.0];
    let profiles = [free(a), free(b)];
    let eq = E / (1.0 - NU);
    let branch = BRANCH / (1.0 - NU);
    let jump = body
        .advance(membrane, gradient, &profiles, 0.0, || false)
        .unwrap();
    for value in jump
        .membrane_force_n_m
        .into_iter()
        .chain(jump.bending_moment_n)
    {
        near(value, 0.0, 2e-16);
    }
    near(
        jump.stored_energy_j_m2,
        (eq + branch) * h * delta * delta / 16.0,
        2e-18,
    );
    near(
        jump.branch_stored_energy_j_m2,
        branch * h * delta * delta / 16.0,
        2e-18,
    );
    assert_eq!(jump.dissipated_j_m2, 0.0);
    // Hold this initially free bending shape. Unequal native relaxation times
    // develop a net restraint force without erasing either layer's self-stress.
    let time: f64 = 2.0;
    let d0 = (-time).exp();
    let d1 = (-time / 3.0).exp();
    let held = body
        .advance(membrane, gradient, &profiles, time, || false)
        .unwrap();
    near(
        held.membrane_force_n_m[0],
        h * delta * branch * (d0 - d1) / 16.0,
        2e-16,
    );
    near(held.bending_moment_n[0], 0.0, 1e-17);
    let branch_storage = branch * h * delta * delta * (d0 * d0 + d1 * d1) / 32.0;
    near(held.branch_stored_energy_j_m2, branch_storage, 2e-18);
    near(
        held.stored_energy_j_m2,
        eq * h * delta * delta / 16.0 + branch_storage,
        2e-18,
    );
    near(
        held.dissipated_j_m2,
        jump.branch_stored_energy_j_m2 - branch_storage,
        2e-18,
    );
    near(held.balance_residual_j_m2, 0.0, 3e-18);
    near(
        body.layer_stress(0).unwrap().bottom_pa[0],
        -(eq + branch * d0) * delta / 4.0,
        2e-14,
    );
    near(
        body.layer_stress(1).unwrap().bottom_pa[0],
        -(eq + branch * d1) * delta / 2.0,
        2e-14,
    );
}

#[test]
fn g3_rotated_anisotropic_thermal_memory_preserves_engineering_shear_and_energy() {
    let diagonal = [100.0, 60.0, 40.0, 30.0, 10.0, 15.0];
    let eq: Tangent6 =
        core::array::from_fn(|i| core::array::from_fn(|j| if i == j { diagonal[i] } else { 0.0 }));
    let branch = eq.map(|row| row.map(|v| 0.5 * v));
    let law = TensorProny::new_mandel(
        eq,
        vec![TensorPronyBranch::from_stiffness_mandel(branch, 0.3).unwrap()],
        0.05,
    )
    .unwrap();
    let h = 0.02;
    let mut body = history(vec![layer(law, h, core::f64::consts::FRAC_PI_4)]);
    let free = [0.002, -0.001, 0.0008];
    let profile = [LaminaThermalStrain::uniform(free).unwrap()];
    body.advance([0.0; 3], [0.0; 3], &profile, 0.0, || false)
        .unwrap();
    let time: f64 = 0.4;
    let decay = (-time / 0.3).exp();
    let step = body
        .advance([0.0; 3], [0.0; 3], &profile, time, || false)
        .unwrap();
    let local =
        [-100.0 * free[0], -60.0 * free[1], -15.0 * free[2]].map(|v| (1.0 + 0.5 * decay) * v);
    let expected = [
        0.5 * (local[0] + local[1]) - local[2],
        0.5 * (local[0] + local[1]) + local[2],
        0.5 * (local[0] - local[1]),
    ];
    for j in 0..3 {
        near(step.membrane_force_n_m[j], h * expected[j], 2e-16);
    }
    let base =
        0.5 * h * (100.0 * free[0] * free[0] + 60.0 * free[1] * free[1] + 15.0 * free[2] * free[2]);
    near(
        step.stored_energy_j_m2,
        base * (1.0 + 0.5 * decay * decay),
        2e-19,
    );
    near(
        step.dissipated_j_m2,
        0.5 * base * (1.0 - decay * decay),
        2e-19,
    );
    near(step.balance_residual_j_m2, 0.0, 2e-19);
}

#[test]
fn g0_incompatible_plane_stress_spectra_refuse_and_failed_advances_leave_all_memory_unchanged() {
    let section = Arc::new(
        LaminateSection::try_new(
            vec![layer(law(0.4, 1.0, false), 0.02, 0.0)],
            ReferenceSurface::Midsurface,
        )
        .unwrap(),
    );
    assert!(matches!(
        LaminateThermalHistory::try_new(section, || false),
        Err(LaminateHistoryError::IncompatiblePlaneStress { layer: 0 })
    ));
    let mut body = history(vec![layer(law(NU, 1.0, false), 0.02, 0.0); 2]);
    body.advance([0.0; 3], [0.0; 3], &[free(0.001); 2], 0.1, || false)
        .unwrap();
    let before = body.clone();
    let mut polls = 0;
    assert_eq!(
        body.advance([0.0; 3], [0.0; 3], &[free(0.002); 2], 0.1, || {
            polls += 1;
            polls == 6
        }),
        Err(LaminateHistoryError::Cancelled)
    );
    assert_eq!(body, before);
    assert!(
        body.advance([0.0; 3], [0.0; 3], &[free(0.002), free(0.1)], 0.1, || false)
            .is_err()
    );
    assert_eq!(body, before);
    // End faces, not merely interior Gauss points, must respect native limits.
    let spike = LaminaThermalStrain::from_faces([0.0; 3], [0.06, 0.0, 0.0]).unwrap();
    assert!(
        body.advance([0.0; 3], [0.0; 3], &[free(0.002), spike], 0.1, || false)
            .is_err()
    );
    assert_eq!(body, before);
    let mut kelvin = history(vec![layer(law(NU, 1.0, true), 0.02, 0.0)]);
    let before = kelvin.clone();
    assert!(
        kelvin
            .advance([0.0; 3], [0.0; 3], &[free(0.001)], 0.0, || false)
            .is_err()
    );
    assert_eq!(kelvin, before);
}
