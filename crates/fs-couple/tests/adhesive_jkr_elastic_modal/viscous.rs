//! G1/G2: original dashpot work, damped oscillator and adhesive impact outcomes.
use super::*;

fn make_damped<'a>(
    law: &'a SourceJkrSphereContact,
    contact: JkrState,
    setup: Setup,
    divisions: f64,
) -> JkrNormalEventModalSystem<'a> {
    let ids = law.normal_state().elastic().interface().surface_materials();
    JkrNormalEventModalSystem::new_with_viscoelastic_connection(
        network(setup),
        law,
        [&ids[0], &ids[1]],
        contact,
        config(law, divisions),
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

// Independent continuum Kelvin-Voigt oscillator, initially at rest.
fn oscillator(extension: f64, omega: f64, ratio: f64, time: f64) -> (f64, f64) {
    let decay = ratio * omega;
    if ratio < 1.0 {
        let w = omega * det::sqrt(1.0 - ratio * ratio);
        let envelope = extension * det::exp(-decay * time);
        (
            envelope * (det::cos(w * time) + decay / w * det::sin(w * time)),
            -envelope * omega * omega / w * det::sin(w * time),
        )
    } else if ratio == 1.0 {
        let envelope = extension * det::exp(-omega * time);
        (
            envelope * (1.0 + omega * time),
            -envelope * omega * omega * time,
        )
    } else {
        let split = omega * det::sqrt(ratio * ratio - 1.0);
        let a = -decay + split;
        let b = -decay - split;
        let ca = -b * extension / (a - b) * det::exp(a * time);
        let cb = a * extension / (a - b) * det::exp(b * time);
        (ca + cb, a * ca + b * cb)
    }
}

#[test]
fn original_dashpot_converges_in_all_damping_regimes_with_positive_native_loss() {
    let law = source();
    let (d, f, v) = scales(&law);
    let k = 0.4 * f / d;
    let omega = det::sqrt(k / MU);
    for ratio in [0.2, 1.0, 2.0] {
        let setup = Setup {
            damping: 2.0 * ratio * det::sqrt(k * MU),
            ..Setup::new(k, 2.0 * d, 5.0 * d)
        };
        let mut errors = Vec::new();
        for divisions in [4.0, 8.0] {
            let mut owner = make_damped(
                &law,
                law.detached_state(-20.0 * d).unwrap(),
                setup,
                divisions,
            );
            assert_eq!(owner.support_damping_n_s_per_m(), setup.damping);
            let initial_energy = owner.total_energy_j().unwrap();
            let mut loss = 0.0;
            let mut spring = 0.5 * k * setup.extension * setup.extension;
            for _ in 0..8 {
                let frame = owner.step(&[0.0; 2]).unwrap().clone();
                assert_eq!(frame.path.event_count(), 0);
                assert!(frame.damping_dissipation_j > 0.0);
                assert!(frame.dissipation_trace().is_some());
                loss += frame.damping_dissipation_j;
                close(
                    owner.total_energy_j().unwrap() + loss,
                    initial_energy,
                    f * d,
                    2e-8,
                );
                close(frame.normal_work_on_network_j, 0.0, f * d, 2e-9);
                close(
                    momentum(&owner),
                    (MA + MB) * setup.common_velocity,
                    (MA + MB) * v,
                    2e-10,
                );
                let extension = owner.state().contact_state().indentation_m()
                    - owner.support_zero_force_indentation_m().unwrap();
                let force = -k * extension - setup.damping * owner.state().velocity_m_s();
                close(frame.connection_force_n, force, f, 2e-9);
                assert!(frame.connection_force_bound_n + 2e-9 * f >= force.abs());
                close(
                    owner
                        .network()
                        .last_frame()
                        .unwrap()
                        .connection_dissipation_j,
                    frame.damping_dissipation_j,
                    f * d,
                    2e-12,
                );
                // The original connection port supplies both spring work and
                // dashpot loss; only source contact work is a boundary input.
                ports(
                    &owner,
                    frame.spring_energy_j - spring + frame.damping_dissipation_j,
                    f * d,
                );
                spring = frame.spring_energy_j;
            }
            let (exact, exact_velocity) = oscillator(setup.extension, omega, ratio, 8.0 * DT);
            let error = ((owner.state().contact_state().indentation_m() + 25.0 * d - exact) / d)
                .abs()
                + ((owner.state().velocity_m_s() - exact_velocity) / v).abs();
            errors.push(error);
            assert_eq!(owner.state().contact_state().branch(), JkrBranch::Detached);
            let continuum_loss = 0.5 * k * setup.extension.powi(2)
                - 0.5 * k * exact * exact
                - 0.5 * MU * exact_velocity * exact_velocity;
            close(loss, continuum_loss, f * d, 0.015);
        }
        assert!(
            errors[0] > 3.5 * errors[1],
            "damping ratio {ratio}: {errors:?}"
        );
        assert!(errors[1] < 0.02, "damping ratio {ratio}: {errors:?}");
    }
}

#[test]
fn original_dashpot_changes_adhesive_rebound_into_capture_without_erasing_snap_energy() {
    let law = source();
    let (d, f, v) = scales(&law);
    let setup = Setup {
        velocity: 10.0 * v,
        ..Setup::new(0.05 * f / d, 2.0 * d, 0.0)
    };
    for (damping, expect_detach) in [(0.0, true), (0.5 * MU * v / d, false)] {
        let mut owner = make_damped(
            &law,
            law.detached_state(-0.25 * d).unwrap(),
            Setup { damping, ..setup },
            20.0,
        );
        let initial = owner.total_energy_j().unwrap();
        let mut released = 0.0;
        let mut dissipated = 0.0;
        let mut attachment = false;
        let mut detachment = false;
        let mut spring = 0.0;
        for _ in 0..80 {
            let frame = owner.step(&[0.0; 2]).unwrap().clone();
            released += frame.released_energy_j;
            dissipated += frame.damping_dissipation_j;
            for event in frame.path.events() {
                attachment |= event.jump.kind == JkrJumpKind::Attach;
                detachment |= event.jump.kind == JkrJumpKind::Detach;
            }
            close(
                owner.total_energy_j().unwrap() + released + dissipated,
                initial,
                f * d,
                8e-8,
            );
            close(
                frame.normal_work_on_network_j,
                -frame.path.normal_work_j,
                f * d,
                2e-8,
            );
            close(
                momentum(&owner),
                (MA + MB) * setup.common_velocity,
                (MA + MB) * v,
                2e-10,
            );
            ports(
                &owner,
                frame.spring_energy_j - spring + frame.damping_dissipation_j,
                f * d,
            );
            spring = frame.spring_energy_j;
            if detachment {
                break;
            }
        }
        assert!(attachment);
        assert_eq!(detachment, expect_detach);
        assert!(released > 0.0);
        if damping > 0.0 {
            assert!(dissipated > 0.0);
            assert_eq!(owner.state().contact_state().branch(), JkrBranch::Attached);
        } else {
            assert_eq!(dissipated, 0.0);
            assert!(owner.last_frame().unwrap().dissipation_trace().is_none());
        }
    }
}

#[test]
fn viscous_force_limits_and_receiving_refusal_preserve_original_bodies_and_loss_history() {
    let law = source();
    let (d, f, v) = scales(&law);
    let setup = Setup {
        damping: 0.2 * MU * v / d,
        velocity: -20.0 * v,
        force_limit: 4.5 * f,
        ..Setup::new(0.5 * f / d, 2.0 * d, 0.0)
    };
    let contact = law.detached_state(-20.0 * d).unwrap();
    let mut owner = make_damped(&law, contact, setup, 20.0);
    let before = states(&owner);
    let token = owner.state();
    assert!(owner.step(&[0.0; 2]).is_err());
    assert_eq!(states(&owner), before);
    assert_eq!(owner.state(), token);
    assert_eq!(owner.samples_rendered(), 0);
    assert!(owner.last_frame().is_none());
    assert!(owner.network().last_frame().is_none());

    let setup = Setup {
        force_limit: 1.0,
        ..setup
    };
    let mut owner = make_damped(&law, contact, setup, 20.0);
    let accepted = owner.step(&[0.0; 2]).unwrap().clone();
    let before = states(&owner);
    let token = owner.state();
    let rejected: Result<(), JkrNormalEventModalError> =
        owner.step_with_admission(&[0.0; 2], &CancelGate::new_clock_free(), |frame, _| {
            assert!(frame.damping_dissipation_j > 0.0);
            assert!(frame.dissipation_trace().is_some());
            Err(JkrNormalEventModalError::Compatibility {
                what: "declared viscous heat receiver refused",
                residual: 1.0,
                tolerance: 0.0,
            })
        });
    assert!(rejected.is_err());
    assert_eq!(states(&owner), before);
    assert_eq!(owner.state(), token);
    assert_eq!(owner.last_frame(), Some(&accepted));
    assert_eq!(owner.samples_rendered(), 1);
    let retry = owner.step(&[0.0; 2]).unwrap().clone();
    let mut reference = make_damped(&law, contact, setup, 20.0);
    reference.step(&[0.0; 2]).unwrap();
    assert_eq!(&retry, reference.step(&[0.0; 2]).unwrap());
    assert_eq!(states(&owner), states(&reference));

    let setup = Setup {
        damping: 0.0,
        ..setup
    };
    let mut old = make(&law, contact, setup, 20.0);
    let mut explicit = make_damped(&law, contact, setup, 20.0);
    assert_eq!(
        old.step(&[0.0; 2]).unwrap(),
        explicit.step(&[0.0; 2]).unwrap()
    );
    assert_eq!(states(&old), states(&explicit));
}
