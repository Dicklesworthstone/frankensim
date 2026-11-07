//! G1/G3/G4: original adhesive microslip, full native response and atomic history.
use super::{MA, MB, RATE, columns, network};
#[path = "../../../fs-contact/tests/adhesive_dmt_mindlin/common.rs"]
mod source_fixture;
use fs_contact::interface_binding::mindlin::loading::path::SourceMindlinPathInput;
use fs_contact::normal_patch::dmt::mindlin::SourceDmtMindlinContact;
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::schedule::force::coupled::contact::adhesive_dmt::mindlin::{
    AdhesiveDmtMindlinConfig, AdhesiveDmtMindlinError, AdhesiveDmtMindlinModalSystem,
};
use fs_exec::CancelGate;
use fs_math::det;
use fs_tribo::TangentialSolveConfig;

#[path = "mindlin/gaps.rs"]
mod gaps;

const Z: f64 = 1e-9;
const MU: f64 = 0.6;
const L: [f64; 3] = [0.7, 0.3, -0.2];
fn close(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:e} != {b:e}"
    );
}
fn config() -> AdhesiveDmtMindlinConfig {
    AdhesiveDmtMindlinConfig {
        maximum_evaluations: 96,
        maximum_step_s: 1e-4,
        maximum_coordinate_increment_m: 3.0 * Z,
        maximum_coordinate_speed_m_s: 1.0,
        coordinate_absolute_tolerance_m: 1e-24,
        coordinate_relative_tolerance: 1e-12,
        work_absolute_tolerance_j: 1e-29,
        work_relative_tolerance: 2e-9,
        shear: TangentialSolveConfig {
            max_iterations: 96,
            maximum_force_n: 1e-3,
            force_absolute_tolerance_n: 1e-18,
            force_relative_tolerance: 1e-9,
        },
    }
}
fn make(
    source: &SourceDmtMindlinContact,
    indentation: f64,
    cross: bool,
) -> AdhesiveDmtMindlinModalSystem<'_> {
    let materials = source
        .normal_source()
        .normal_state()
        .elastic()
        .interface()
        .surface_materials();
    let [normal, shear, _] = columns();
    AdhesiveDmtMindlinModalSystem::new(
        network([0.0; 3], cross, RATE),
        source,
        [&materials[0], &materials[1]],
        indentation,
        [normal, shear],
        config(),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}
fn all_states(s: &AdhesiveDmtMindlinModalSystem<'_>) -> Vec<ModalAcousticState> {
    s.components()
        .iter()
        .flat_map(|m| m.states().iter().copied())
        .collect()
}
fn relative(s: &AdhesiveDmtMindlinModalSystem<'_>) -> ([f64; 3], [f64; 3]) {
    (
        core::array::from_fn(|i| {
            s.components()[1].states()[i].displacement_m_sqrt_kg / det::sqrt(MB)
                - s.components()[0].states()[i].displacement_m_sqrt_kg / det::sqrt(MA)
        }),
        core::array::from_fn(|i| {
            s.components()[1].states()[i].velocity_m_sqrt_kg_per_s / det::sqrt(MB)
                - s.components()[0].states()[i].velocity_m_sqrt_kg_per_s / det::sqrt(MA)
        }),
    )
}
fn momentum(s: &AdhesiveDmtMindlinModalSystem<'_>) -> [f64; 3] {
    core::array::from_fn(|i| {
        det::sqrt(MA) * s.components()[0].states()[i].velocity_m_sqrt_kg_per_s
            + det::sqrt(MB) * s.components()[1].states()[i].velocity_m_sqrt_kg_per_s
    })
}
// Independent two-mass/one-oblique-spring solution, including the unforced z mode.
fn response(cross: bool, h: f64) -> ([[f64; 3]; 3], [[f64; 3]; 3]) {
    let mass = MA * MB / (MA + MB);
    let norm2 = L.iter().map(|x| x * x).sum::<f64>();
    let (parallel_q, parallel_v) = if cross {
        let omega = det::sqrt(3e4 * norm2 / mass);
        (
            (1.0 - det::cos(omega * h)) / (mass * omega * omega),
            det::sin(omega * h) / (mass * omega),
        )
    } else {
        (h * h / (2.0 * mass), h / mass)
    };
    let isotropic_q = h * h / (2.0 * mass);
    let isotropic_v = h / mass;
    (
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (if i == j { isotropic_q } else { 0.0 })
                    + (parallel_q - isotropic_q) * L[i] * L[j] / norm2
            })
        }),
        core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                (if i == j { isotropic_v } else { 0.0 })
                    + (parallel_v - isotropic_v) * L[i] * L[j] / norm2
            })
        }),
    )
}
fn released(q: [f64; 3], v: [f64; 3], cross: bool, h: f64) -> ([f64; 3], [f64; 3]) {
    if !cross {
        return (v.map(|x| h * x), v);
    }
    let mass = MA * MB / (MA + MB);
    let norm2 = L.iter().map(|x| x * x).sum::<f64>();
    let omega = det::sqrt(3e4 * norm2 / mass);
    let ql = (0..3).map(|i| L[i] * q[i]).sum::<f64>() / norm2;
    let vl = (0..3).map(|i| L[i] * v[i]).sum::<f64>() / norm2;
    let c = det::cos(omega * h);
    let sn = det::sin(omega * h);
    (
        core::array::from_fn(|i| h * v[i] + ((c - 1.0) * ql + (sn / omega - h) * vl) * L[i]),
        core::array::from_fn(|i| v[i] + (-omega * sn * ql + (c - 1.0) * vl) * L[i]),
    )
}

#[test]
fn g1_adhesive_microslip_at_zero_and_tensile_net_load_drives_the_full_network() {
    let source = source_fixture::source(0.01, Some(MU), "initial").unwrap();
    let zero = source
        .normal_source()
        .evaluate_contact_load(0.0)
        .unwrap()
        .indentation_m;
    for factor in [1.0, 0.5] {
        for cross in [false, true] {
            let indentation = factor * zero;
            let normal = source.normal_source().evaluate(indentation).unwrap();
            let shear_modulus =
                1.0 / ((2.0 - 0.25) / (1e11 / (2.0 * 1.25)) + (2.0 - 0.25) / (2e11 / (2.0 * 1.25)));
            let stiffness = 8.0 * shear_modulus * normal.contact_radius_m;
            let capacity = MU * normal.repulsive_normal_force_n;
            let threshold = 1.5 * capacity / stiffness;
            let x = 0.25;
            let slip = x * threshold;
            // Independently integrated classical initial Mindlin loading curve.
            let mean_shear = capacity * (1.0 - (2.0 / (5.0 * x)) * (1.0 - det::pow(1.0 - x, 2.5)));
            let endpoint = capacity * (1.0 - det::pow(1.0 - x, 1.5));
            let efforts = [normal.normal_force_n, mean_shear];
            let mut system = make(&source, indentation, cross);
            let h = system.sample_period_s();
            let (s, v) = response(cross, h);
            for i in 0..2 {
                for j in 0..2 {
                    close(system.compliance_m_per_n()[i][j], s[i][j], 1e-18, 2e-11);
                }
            }
            if cross {
                assert_ne!(system.compliance_m_per_n()[0][1], 0.0);
            }
            let prescribed = [
                (s[0][0] * efforts[0] + s[0][1] * efforts[1]) / h,
                (slip + s[1][0] * efforts[0] + s[1][1] * efforts[1]) / h,
            ];
            let before = system.total_energy_j().unwrap();
            let frame = system
                .step(&[0.0; 6], prescribed, [0.0, Z])
                .unwrap()
                .clone();
            close(frame.position_m[0], indentation, 2e-23, 2e-10);
            close(frame.position_m[1], slip, 2e-23, 2e-9);
            close(frame.mean_force_n[0], efforts[0], 3e-18, 2e-9);
            close(frame.mean_force_n[1], mean_shear, 3e-18, 2e-8);
            close(frame.endpoint_shear_force_n, endpoint, 3e-18, 2e-8);
            assert!(frame.endpoint_shear_force_n.abs() < capacity);
            assert!(frame.tangential_energy_j > 0.0 && frame.microslip_dissipation_j > 0.0);
            if factor < 1.0 {
                assert!(frame.normal_response.normal_force_n < 0.0);
            }
            let (qend, vend) = relative(&system);
            for i in 0..3 {
                close(
                    qend[i],
                    -s[i][0] * efforts[0] - s[i][1] * efforts[1],
                    3e-23,
                    3e-8,
                );
                close(
                    vend[i],
                    -v[i][0] * efforts[0] - v[i][1] * efforts[1],
                    3e-16,
                    3e-8,
                );
                close(momentum(&system)[i], 0.0, 2e-25, 0.0);
            }
            if cross {
                assert_ne!(qend[2], 0.0);
            }
            close(
                system.total_energy_j().unwrap()
                    + frame.microslip_dissipation_j
                    + frame.network_dissipation_j
                    + frame.numerical_dissipation_j
                    - before,
                frame.prescribed_work_j,
                2e-26,
                3e-8,
            );
            assert!(frame.contact_work_residual_j.abs() <= frame.contact_work_tolerance_j);
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            assert!(frame.evaluations <= config().maximum_evaluations);
        }
    }
}

#[test]
fn g3_loading_unloading_and_shear_reversals_retain_original_memory_and_work() {
    let source = source_fixture::source(0.01, Some(MU), "initial").unwrap();
    let zero = source
        .normal_source()
        .evaluate_contact_load(0.0)
        .unwrap()
        .indentation_m;
    let mut system = make(&source, 0.6 * zero, true);
    let initial = system.total_energy_j().unwrap();
    let h = system.sample_period_s();
    let (s, v) = response(true, h);
    let mut supplied = 0.0;
    let mut losses = 0.0;
    let mut heat = 0.0;
    for (factor, slip) in [
        (0.8, 0.002 * Z),
        (1.1, 0.003 * Z),
        (0.6, -0.002 * Z),
        (0.4, -0.004 * Z),
        (0.9, 0.003 * Z),
        (1.0, 0.002 * Z),
    ] {
        let target = factor * zero;
        let original = source
            .trial_coupled_path(
                system.state(),
                &SourceMindlinPathInput {
                    indentation_m: target,
                    free_slip_increment_m: slip,
                    dt_s: h,
                    compliance_m_per_n: 0.0,
                    solve: config().shear,
                },
                &mut || false,
            )
            .unwrap();
        let f = [
            original.mean_normal_force_n,
            original.tangential.average_force_n,
        ];
        let (q0, v0) = relative(&system);
        let (free, free_v) = released(q0, v0, true, h);
        let desired = [target - system.state().indentation_m(), slip];
        let prescribed =
            core::array::from_fn(|i| (desired[i] - free[i] + s[i][0] * f[0] + s[i][1] * f[1]) / h);
        let old_slip = system.state().slip_m();
        let frame = system
            .step(&[0.0; 6], prescribed, [0.0, Z])
            .unwrap()
            .clone();
        close(system.state().indentation_m(), target, 3e-23, 3e-9);
        close(system.state().slip_m(), old_slip + slip, 3e-23, 3e-9);
        close(
            frame.tangential_stored_change_j,
            original.tangential.stored_change_j,
            1e-28,
            2e-8,
        );
        close(
            frame.microslip_dissipation_j,
            original.tangential.dissipated_j,
            1e-28,
            2e-8,
        );
        let (q1, v1) = relative(&system);
        for i in 0..3 {
            close(
                q1[i] - q0[i],
                free[i] - s[i][0] * f[0] - s[i][1] * f[1],
                5e-23,
                5e-8,
            );
            close(
                v1[i],
                free_v[i] - v[i][0] * f[0] - v[i][1] * f[1],
                5e-16,
                5e-8,
            );
            close(momentum(&system)[i], 0.0, 5e-25, 0.0);
        }
        supplied += frame.prescribed_work_j;
        heat += frame.microslip_dissipation_j;
        losses += frame.microslip_dissipation_j
            + frame.network_dissipation_j
            + frame.numerical_dissipation_j;
        close(
            system.total_energy_j().unwrap() + losses,
            initial + supplied,
            1e-25,
            5e-8,
        );
    }
    assert!(heat > 0.0 && system.state().history_intervals() > 0);
    assert_eq!(system.samples_rendered(), 6);
}

#[test]
fn g4_bad_brackets_receiver_refusal_and_late_cancellation_preserve_the_accepted_step() {
    let source = source_fixture::source(0.01, Some(MU), "initial").unwrap();
    let initial = source
        .normal_source()
        .evaluate_contact_load(0.0)
        .unwrap()
        .indentation_m;
    let mut owner = make(&source, initial, false);
    let state = owner.state().clone();
    let native = all_states(&owner);
    let h = owner.sample_period_s();
    let target = 1.1 * initial;
    let slip = 0.001 * Z;
    let direct = source
        .trial_coupled_path(
            &state,
            &SourceMindlinPathInput {
                indentation_m: target,
                free_slip_increment_m: slip,
                dt_s: h,
                compliance_m_per_n: 0.0,
                solve: config().shear,
            },
            &mut || false,
        )
        .unwrap();
    let c = owner.compliance_m_per_n();
    let f = [
        direct.mean_normal_force_n,
        direct.tangential.average_force_n,
    ];
    let rate = [
        (target - initial + c[0][0] * f[0] + c[0][1] * f[1]) / h,
        (slip + c[1][0] * f[0] + c[1][1] * f[1]) / h,
    ];
    for bracket in [[f64::NAN, Z], [Z, 0.0], [0.0, 0.0]] {
        assert!(owner.step(&[0.0; 6], rate, bracket).is_err());
        assert_eq!(owner.state(), &state);
        assert_eq!(all_states(&owner), native);
        assert_eq!(owner.samples_rendered(), 0);
        assert!(owner.last_frame().is_none());
    }
    #[derive(Debug)]
    enum Rejected {
        Native(AdhesiveDmtMindlinError),
        Receiver,
    }
    impl From<AdhesiveDmtMindlinError> for Rejected {
        fn from(e: AdhesiveDmtMindlinError) -> Self {
            Self::Native(e)
        }
    }
    let gate = CancelGate::new_clock_free();
    let rejected: Result<(), Rejected> =
        owner.step_with_admission(&[0.0; 6], rate, [0.0, Z], &gate, |frame, candidate| {
            assert_eq!(frame.sample, 1);
            assert_ne!(candidate[0].states(), &native[..3]);
            Err(Rejected::Receiver)
        });
    match rejected {
        Err(Rejected::Receiver) => {}
        Err(Rejected::Native(e)) => panic!("native trial failed: {e}"),
        Ok(()) => panic!("receiver accepted"),
    }
    assert_eq!(owner.state(), &state);
    assert_eq!(all_states(&owner), native);
    let cancelled = CancelGate::new_clock_free();
    let rejected: Result<(), AdhesiveDmtMindlinError> =
        owner.step_with_admission(&[0.0; 6], rate, [0.0, Z], &cancelled, |_, _| {
            cancelled.request();
            Ok(())
        });
    assert!(rejected.is_err());
    assert_eq!(owner.state(), &state);
    assert_eq!(all_states(&owner), native);
    let expected = owner.step(&[0.0; 6], rate, [0.0, Z]).unwrap().clone();
    let mut fresh = make(&source, initial, false);
    assert_eq!(fresh.step(&[0.0; 6], rate, [0.0, Z]).unwrap(), &expected);
    assert_eq!(fresh.state(), owner.state());
    assert_eq!(all_states(&fresh), all_states(&owner));
}
