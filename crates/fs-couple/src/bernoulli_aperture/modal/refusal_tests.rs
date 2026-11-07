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
            time_step_s: 1.0 / 384000.0,
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

// A near-closed aperture: the net flows are about 1e-11 m^3/s while `wave` is a difference of kPa-sized
// pressures over the bore impedance, so its arithmetic carries eps |P| / z of rounding. These are the numbers of
// the step that refused a 2.65 s D4 render (residual 4.2e-21, flow scale 6.16e-11, relative tolerance 1e-11).
const CLOSURE_RESIDUAL: f64 = 4.2085786343961645e-21;
const CLOSURE_FLOW_SCALE: f64 = 6.15680493824026e-11;
const CLOSURE_PRESSURE_PA: f64 = 6178.0;
const CLOSURE_IMPEDANCE: f64 = 1.73e8;

fn closure_frame() -> (ModalAperture, ModalApertureFrame, f64) {
    let a = model(1e-11).unwrap();
    let input = ApertureDrive {
        upstream_pressure_pa: 50.0,
        incoming_pressure_pa: 0.0,
        body_flow_m3_s: 0.0,
    };
    let frame = a.preview(input, [0.001, 0.0]).unwrap();
    let scale = a.stored_energy_j()
        + frame.aperture.stored_energy_j
        + frame.aperture.pressure_work_j.abs()
        + frame.aperture.mechanical_work_j.abs()
        + frame.aperture.dissipated_energy_j.abs();
    (a, frame, scale)
}

fn refused_flow(error: AcousticRealizeError) -> fs_phs::SolveCertificate {
    let AcousticRealizeError::AlgebraicSolve(fs_phs::PhsError::NewtonRefused { certificate: c }) = error
    else {
        panic!("expected typed flow certificate: {error:?}")
    };
    assert_eq!(c.site, NewtonSite::ModalFlow);
    c
}

#[test]
fn a_near_closure_flow_residual_at_the_rounding_of_its_pressure_terms_is_admitted_with_a_recorded_floor() {
    let (a, mut frame, scale) = closure_frame();
    frame.aperture.flow_residual_m3_s = CLOSURE_RESIDUAL;
    // The floor over the net flows alone cannot admit it: the refusal as it happened.
    let alone = refused_flow(
        a.admit_balance_with_rounding(&frame.aperture, scale, CLOSURE_FLOW_SCALE, CLOSURE_FLOW_SCALE, 40)
            .unwrap_err(),
    );
    assert_eq!(alone.requested_budget, 1e-11 * CLOSURE_FLOW_SCALE);
    assert!(alone.residual > alone.effective_budget && !alone.floor_raised);
    // The floor over the magnitude the arithmetic works at: eps |P| / z is about 8e-21 here.
    let rounding_scale = CLOSURE_FLOW_SCALE + CLOSURE_PRESSURE_PA / CLOSURE_IMPEDANCE;
    let (_, flow) = a
        .admit_balance_with_rounding(&frame.aperture, scale, CLOSURE_FLOW_SCALE, rounding_scale, 40)
        .unwrap();
    assert!(flow.floor_raised);
    assert_eq!(flow.requested_budget, alone.requested_budget);
    assert_eq!(flow.effective_budget, flow.rounding_floor);
    assert!(flow.residual <= flow.effective_budget);
    println!(
        "{{\"suite\":\"modal-refusal\",\"case\":\"near-closure-flow\",\"site\":\"ModalFlow\",\"residual\":{},\"budget\":{},\"floor\":{},\"iterations\":{},\"verdict\":\"admit\",\"git_sha\":\"{}\"}}",
        flow.residual, flow.effective_budget, flow.rounding_floor, flow.iterations, git_sha()
    );
}

#[test]
fn the_recorded_floor_does_not_hide_a_real_flow_imbalance_at_the_same_scale() {
    let (a, mut frame, scale) = closure_frame();
    let rounding_scale = CLOSURE_FLOW_SCALE + CLOSURE_PRESSURE_PA / CLOSURE_IMPEDANCE;
    // Six orders above the rounding of the pressure terms: a defect, not noise.
    for defect in [1e-15, -1e-15] {
        frame.aperture.flow_residual_m3_s = defect;
        let refused = refused_flow(
            a.admit_balance_with_rounding(&frame.aperture, scale, CLOSURE_FLOW_SCALE, rounding_scale, 40)
                .unwrap_err(),
        );
        assert!(refused.residual > refused.effective_budget && refused.floor_raised);
    }
}

#[test]
fn where_the_floor_is_under_the_declared_ceiling_the_certificate_is_the_explicit_one() {
    let (a, frame, scale) = closure_frame();
    let flow_scale = frame.flow_certificate.scale;
    let (work_a, flow_a) = a.admit_balance(&frame.aperture, scale, flow_scale, 7).unwrap();
    let (work_b, flow_b) = a
        .admit_balance_with_rounding(&frame.aperture, scale, flow_scale, flow_scale, 7)
        .unwrap();
    assert_eq!((work_a, flow_a), (work_b, flow_b));
    // The certificate the gate gave before the allowance existed (explicit declaration, floor over the net flows):
    // equal field for field at a step whose floor is under the declared ceiling.
    let explicit = admit_newton_iterate(
        SolveBudget { relative: 1e-11, absolute: 0.0, max_iterations: 0 }.certificate(
            NewtonSite::ModalFlow,
            frame.aperture.flow_residual_m3_s.abs(),
            flow_scale.max(f64::MIN_POSITIVE),
            arithmetic_rounding_floor(2 * a.state.len() + 8, flow_scale),
            7,
            0,
        ),
    )
    .unwrap();
    assert_eq!(explicit, flow_a);
    println!(
        "{{\"suite\":\"modal-refusal\",\"case\":\"explicit-certificate-unchanged\",\"site\":\"ModalFlow\",\"residual\":{:e},\"requested\":{:e},\"effective\":{:e},\"floor\":{:e},\"floor_raised\":{},\"verdict\":\"equal\"}}",
        flow_a.residual, flow_a.requested_budget, flow_a.effective_budget, flow_a.rounding_floor, flow_a.floor_raised
    );
    assert_eq!(flow_a.effective_budget, flow_a.requested_budget);
    assert!(!flow_a.floor_raised && flow_a.rounding_floor <= flow_a.requested_budget);
    // A larger rounding magnitude changes only the floor, and the ceiling only where the floor passes it.
    let (work_c, flow_c) = a
        .admit_balance_with_rounding(&frame.aperture, scale, flow_scale, flow_scale + 50.0 / 1e6, 7)
        .unwrap();
    assert_eq!(work_c, work_a);
    assert!(flow_c.rounding_floor >= flow_a.rounding_floor);
    assert_eq!(flow_c.requested_budget, flow_a.requested_budget);
    assert_eq!(flow_c.effective_budget, flow_a.requested_budget.max(flow_c.rounding_floor));
    assert_eq!(flow_c.floor_raised, flow_c.rounding_floor > flow_c.requested_budget);
    assert_eq!(flow_c.residual, flow_a.residual);
}
