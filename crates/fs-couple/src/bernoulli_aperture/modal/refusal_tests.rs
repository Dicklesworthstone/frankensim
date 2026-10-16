//! G0/G4 DEF-05: an external-force work defect cannot hide in 128*tolerance.
use super::*;

fn git_sha() -> String {
    let output = std::process::Command::new("git")
        .args(["rev-parse", "HEAD"])
        .output()
        .unwrap();
    String::from_utf8(output.stdout).unwrap().trim().to_owned()
}

fn model(tolerance: f64) -> Result<ModalAperture, AcousticRealizeError> {
    ModalAperture::new(
        ModalApertureSpec {
            rest_opening_m: 4e-4,
            width_m: 0.013,
            modes: vec![ApertureMode {
                mass_kg: 1e-5,
                stiffness_n_m: 500.0,
                damping_n_s_m: 0.03,
                pressure_area_m2: 1.25e-4,
                gap_projection: 1.0,
                force_projection: [1.0, 0.0],
            }],
            density_kg_m3: 1.2,
            impedance_pa_s_m3: 1e6,
            time_step_s: 1.0 / 96000.0,
            max_steps: 4096,
            limits: PairedApertureLimits {
                maximum_abs_displacement_m: 0.01,
                maximum_abs_velocity_m_s: 100.0,
                maximum_energy_j: 1.0,
                max_evaluations: 256,
                velocity_tolerance_m_s: 1e-13,
                relative_tolerance: tolerance,
            },
        },
        vec![BladeState::default()],
        Obstacle::new(
            vec![-1.0],
            1,
            1,
            vec![0.0],
            vec![1.0],
            1e8,
            2.0,
            "synthetic refusal-certificate slit".into(),
        )
        .unwrap(),
    )
}

#[test]
fn modal_aperture_limits_refuse_tolerance_that_disables_the_gate() {
    let cap = 1.0_f64 / 128.0;
    assert!(model(cap.next_down()).is_ok());
    assert!(model(cap).is_err());
    assert!(model(cap.next_up()).is_err());
    println!(
        "{{\"suite\":\"modal-refusal\",\"case\":\"relative-tolerance-cap\",\"site\":\"Admission\",\"residual\":{cap},\"budget\":{cap},\"floor\":0,\"iterations\":0,\"verdict\":\"refuse\",\"git_sha\":\"{}\"}}",
        git_sha()
    );
}

#[test]
fn modal_aperture_refuses_injected_work_imbalance_at_loose_tolerance() {
    let a = model((1.0_f64 / 128.0).next_down()).unwrap();
    let state = a.modal_states().to_vec();
    let input = ApertureDrive {
        upstream_pressure_pa: 50.0,
        incoming_pressure_pa: 0.0,
        body_flow_m3_s: 0.0,
    };
    let mut frame = a.preview(input, [0.001, 0.0]).unwrap();
    assert!(!frame.work_certificate.floor_raised);
    assert!(frame.work_certificate.iterations > 0);
    let scale = a.stored_energy_j()
        + frame.aperture.stored_energy_j
        + frame.aperture.pressure_work_j.abs()
        + frame.aperture.mechanical_work_j.abs()
        + frame.aperture.dissipated_energy_j.abs();
    let velocity = (frame.modes[0].displacement_m - state[0].displacement_m) / a.spec.time_step_s;
    assert!(velocity != 0.0);
    // A test-only external force whose work was not included in the equations.
    // Feed its 1% defect through the SAME final gate used by real preview().
    let unaccounted_force = 0.01 * scale / (a.spec.time_step_s * velocity);
    frame.aperture.mechanical_work_j += a.spec.time_step_s * unaccounted_force * velocity;
    let scale = a.stored_energy_j()
        + frame.aperture.stored_energy_j
        + frame.aperture.pressure_work_j.abs()
        + frame.aperture.mechanical_work_j.abs()
        + frame.aperture.dissipated_energy_j.abs();
    let error = a
        .admit_balance(
            &frame.aperture,
            scale,
            frame.flow_certificate.scale,
            frame.work_certificate.iterations,
        )
        .unwrap_err();
    let AcousticRealizeError::AlgebraicSolve(fs_phs::PhsError::NewtonRefused { certificate: c }) =
        error
    else {
        panic!("expected typed work certificate: {error:?}")
    };
    assert_eq!(c.site, NewtonSite::ModalWork);
    assert!(c.residual > c.effective_budget);
    assert_eq!(a.modal_states(), state);
    assert_eq!(a.accepted_steps(), 0);
    println!(
        "{{\"suite\":\"modal-refusal\",\"case\":\"unaccounted-force-work\",\"site\":\"ModalWork\",\"residual\":{},\"budget\":{},\"floor\":{},\"iterations\":{},\"verdict\":\"refuse\",\"git_sha\":\"{}\"}}",
        c.residual,
        c.effective_budget,
        c.rounding_floor,
        c.iterations,
        git_sha()
    );
}
