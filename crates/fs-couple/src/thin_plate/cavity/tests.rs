use super::*;
use crate::cavity_modal::CavityNeckLoss;
use fs_material::gas::{GasSpec, GasState};
use fs_nlmodal::{SosModalStorage, StressChannel};
use fs_scenario::{HelmholtzCavity, NeckEnd};

fn source(resistance: f64) -> HelmholtzModalCavity {
    HelmholtzModalCavity::new(
        HelmholtzCavity {
            volume_m3: 1.0,
            neck_radius_m: 0.005,
            neck_length_m: 0.05,
            inner_end: NeckEnd::Flanged,
            outer_end: NeckEnd::Unflanged,
        },
        GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap(),
        CavityNeckLoss::Resistance(resistance),
    )
    .unwrap()
}
fn body(nonlinear: bool) -> VkBody {
    let mut body = VkBody::from_modal_storage_with_damping(
        SosModalStorage {
            omegas: vec![10.0, 16.0],
            channels: if nonlinear {
                vec![StressChannel {
                    coefficient: 3e5,
                    coupling: vec![1.0, 0.2, 0.2, 0.6],
                }]
            } else {
                vec![]
            },
        },
        &[0.4, -0.1, -0.1, 0.2],
        &[1.0, 0.3],
        vec![0.04, 0.01],
    )
    .unwrap()
    .with_relaxation_branches(
        vec![fs_phs::RelaxationBranch {
            projection: vec![1.0, 0.0, -0.5, 0.0],
            stiffness: 30.0,
            relaxation_time_s: 0.2,
        }],
        1,
    )
    .unwrap();
    body.x.copy_from_slice(&[0.02, 0.03, -0.01, -0.04, 0.003]);
    body
}
fn attach(body: VkBody, r: f64) -> (VkBody, VkCavityLoad) {
    body.with_cavity_load(
        &source(r),
        &[0.02, -0.01],
        VkCavityInitial {
            pressure_increment_pa: 10.0,
            neck_outflow_m3_s: 1e-5,
        },
        32,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

// Independent physical equations use pressure and outward neck flow, not the
// implementation's gas compression/momentum coordinates or J/R matrices.
fn rhs(x: &[f64; 9], k: f64, m: f64, r: f64, nonlinear: bool) -> [f64; 9] {
    let [q, v, t, w, z, p, flow, _, _] = *x;
    let shape = q * q + 0.4 * q * t + 0.6 * t * t;
    let beta = if nonlinear { 3e5 } else { 0.0 };
    let a = 30_f64.sqrt();
    let elastic = a * (q - 0.5 * t) - z;
    [
        v,
        0.6 - 100.0 * q - beta * shape * (q + 0.2 * t) - 0.4 * v + 0.1 * w - a * elastic + 0.02 * p,
        w,
        0.18 - 256.0 * t - beta * shape * (0.2 * q + 0.6 * t) + 0.1 * v - 0.2 * w
            + 0.5 * a * elastic
            - 0.01 * p,
        elastic / 0.2,
        k * (-0.02 * v + 0.01 * w - flow),
        (p - r * flow) / m,
        0.6 * (v + 0.3 * w),
        0.4 * v * v - 0.2 * v * w + 0.2 * w * w + elastic * elastic / 0.2 + r * flow * flow,
    ]
}
fn reference(nonlinear: bool) -> [f64; 9] {
    let s = source(80.0);
    let area2 = s.neck_area_m2().powi(2);
    let k = s.compression_stiffness_n_m() / area2;
    let m = s.neck_mass_kg() / area2;
    let r = s.neck_resistance_n_s_m() / area2;
    let mut x = [0.02, 0.03, -0.01, -0.04, 0.003, 10.0, 1e-5, 0.0, 0.0];
    let dt = 0.025 / 6000.0;
    for _ in 0..6000 {
        let a = rhs(&x, k, m, r, nonlinear);
        let b = rhs(
            &core::array::from_fn(|j| x[j] + 0.5 * dt * a[j]),
            k,
            m,
            r,
            nonlinear,
        );
        let c = rhs(
            &core::array::from_fn(|j| x[j] + 0.5 * dt * b[j]),
            k,
            m,
            r,
            nonlinear,
        );
        let d = rhs(
            &core::array::from_fn(|j| x[j] + dt * c[j]),
            k,
            m,
            r,
            nonlinear,
        );
        for j in 0..9 {
            x[j] += dt * (a[j] + 2.0 * b[j] + 2.0 * c[j] + d[j]) / 6.0;
        }
    }
    x
}
fn run(nonlinear: bool, steps: usize) -> ([f64; 9], f64) {
    let original = body(nonlinear);
    let prefix = original.state().to_vec();
    let (body, view) = attach(original, 80.0);
    assert_eq!(&body.state()[..5], prefix);
    let start = view.observe(body.state()).unwrap();
    assert!((start.pressure_increment_pa - 10.0).abs() < 1e-13);
    assert!((start.neck_outflow_m3_s - 1e-5).abs() < 1e-19);
    let mut body = body.into_prepared(0.025 / steps as f64).unwrap();
    let initial_energy = body.stored_energy_j();
    let (mut work, mut heat) = (0.0, 0.0);
    for _ in 0..steps {
        let record = body.drive_and_radiate(0.6, 1.2, 1.0).unwrap();
        work += record.work.supplied;
        heat += record.work.dissipated;
    }
    let s = body.state();
    let gas = view.observe(s).unwrap();
    assert!((work - 0.6 * ((s[0] - prefix[0]) + 0.3 * (s[2] - prefix[2]))).abs() < 2e-10);
    let defect = body.stored_energy_j() + heat - initial_energy - work;
    assert!(defect.abs() < 2e-9, "balance {defect}");
    (
        [
            s[0],
            s[1],
            s[2],
            s[3],
            s[4],
            gas.pressure_increment_pa,
            gas.neck_outflow_m3_s,
            work,
            heat,
        ],
        defect,
    )
}

#[test]
fn complete_linear_and_nonlinear_material_cavity_match_independent_physical_dynamics() {
    for nonlinear in [false, true] {
        let expected = reference(nonlinear);
        let scales = [0.1, 1.0, 0.1, 1.0, 1.0, 100.0, 0.01, 1.0, 1.0];
        let mut errors = Vec::new();
        for steps in [125, 250, 500] {
            let (actual, _) = run(nonlinear, steps);
            let error = (0..9)
                .map(|i| ((actual[i] - expected[i]) / scales[i]).abs())
                .fold(0.0_f64, f64::max);
            errors.push(error);
            assert!((actual[7] - expected[7]).abs() < 2e-7);
            assert!((actual[8] - expected[8]).abs() < 2e-7);
        }
        if nonlinear {
            assert!(
                errors[1] < 0.4 * errors[0] && errors[2] < 0.4 * errors[1],
                "{errors:?}"
            );
        }
        assert!(errors[2] < 2e-6, "nonlinear={nonlinear}: {errors:?}");
    }
}

#[test]
fn pressure_uses_signed_surface_shapes_and_a_cancelling_surface_stays_silent() {
    let make = || {
        VkBody::from_modal_storage_with_damping(
            SosModalStorage {
                omegas: vec![10.0, 10.0],
                channels: vec![],
            },
            &[0.0; 4],
            &[1.0, 0.0],
            vec![0.1, 0.1],
        )
        .unwrap()
    };
    let (body, view) = make()
        .with_cavity_load(
            &source(0.0),
            &[0.0, -0.03],
            VkCavityInitial {
                pressure_increment_pa: 10.0,
                neck_outflow_m3_s: 0.0,
            },
            16,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut prepared = body.into_prepared(0.001).unwrap();
    let energy = prepared.stored_energy_j();
    for _ in 0..100 {
        prepared.drive_and_radiate(0.0, 1.2, 1.0).unwrap();
    }
    assert!(prepared.state()[0].abs() < 1e-15);
    assert!(prepared.state()[2] < -1e-4);
    assert!((prepared.stored_energy_j() - energy).abs() < 1e-12);
    assert!(
        view.observe(prepared.state())
            .unwrap()
            .neck_kinetic_energy_j
            > 0.0
    );

    let mut original = make();
    original.x.copy_from_slice(&[0.001, 0.1, 0.001, 0.1]);
    let (loaded, view) = original
        .with_cavity_load(
            &source(80.0),
            &[0.01, -0.01],
            VkCavityInitial::default(),
            16,
            &CancelGate::new_clock_free(),
        )
        .unwrap();
    let mut loaded = loaded.into_prepared(0.001).unwrap();
    for _ in 0..100 {
        loaded.drive_and_radiate(0.0, 1.2, 1.0).unwrap();
    }
    assert!((loaded.state()[0] - loaded.state()[2]).abs() < 1e-13);
    assert!(
        view.observe(loaded.state())
            .unwrap()
            .pressure_increment_pa
            .abs()
            < 1e-11
    );
}

#[test]
fn cavity_state_survives_exact_subdivision_and_failed_prepared_steps() {
    let (body_a, view) = attach(body(false), 80.0);
    let (body_b, _) = attach(body(false), 80.0);
    let mut a = body_a.into_prepared(0.01).unwrap();
    let mut b = body_b.into_prepared(0.0025).unwrap();
    a.drive_and_radiate(0.7, 1.2, 1.0).unwrap();
    for _ in 0..4 {
        b.drive_and_radiate(0.7, 1.2, 1.0).unwrap();
    }
    for (x, y) in a.state().iter().zip(b.state()) {
        assert!((x - y).abs() < 2e-11);
    }
    let before = a.state().to_vec();
    let observed = view.observe(&before).unwrap();
    assert!(a.drive_and_radiate(1.0, -1.0, 1.0).is_err());
    let mut polls = 0;
    assert!(
        a.drive_and_radiate_controlled(1.0, 1.2, 1.0, || {
            polls += 1;
            polls == 2
        })
        .is_err()
    );
    assert_eq!(a.state(), before);
    assert_eq!(view.observe(a.state()).unwrap(), observed);
    assert!(
        body(false)
            .with_cavity_load(
                &source(0.0),
                &[0.01],
                VkCavityInitial::default(),
                32,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    assert!(
        body(false)
            .with_cavity_load(
                &source(0.0),
                &[0.01, 0.02],
                VkCavityInitial::default(),
                6,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    assert!(
        body(false)
            .with_cavity_load(
                &source(0.0),
                &[0.01, f64::NAN],
                VkCavityInitial::default(),
                32,
                &CancelGate::new_clock_free()
            )
            .is_err()
    );
    assert!(view.observe(&before[..before.len() - 1]).is_err());
}
