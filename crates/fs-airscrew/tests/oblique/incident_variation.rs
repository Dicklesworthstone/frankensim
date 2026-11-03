use super::*;
use fs_airscrew::dynamic_inflow::first_harmonic::HarmonicInduction;

fn variation(frame: DiskFrame, arm: [f64; 3]) -> [f64; 3] {
    let y = dot(arm, frame.radial_reference());
    let z = dot(arm, frame.azimuth_reference());
    let components = [1.2 * y + 0.8 * z * z - 0.4, 0.25 * z, -0.15 * y];
    core::array::from_fn(|j| {
        components[0] * frame.thrust_axis()[j]
            + components[1] * frame.radial_reference()[j]
            + components[2] * frame.azimuth_reference()[j]
    })
}

#[test]
fn spatial_air_loads_and_power_match_an_independent_blade_quadrature() {
    let rotor = rotor();
    let frame = DiskFrame::new(AXIS, RADIAL, SpinDirection::AlongThrust).unwrap();
    let rule = AzimuthQuadrature::new(64).unwrap();
    let offsets = rule.section_offsets(&rotor, frame).unwrap();
    let mut samples = 0;
    let load = rule
        .evaluate_with_incident_variation(
            &rotor,
            RHO,
            frame,
            [7.0, 0.0, 0.0],
            ZERO,
            HarmonicInduction {
                mean_mps: 2.0,
                first_mps: [0.0; 2],
            },
            OMEGA,
            |sample, arm| {
                assert_eq!(sample, samples);
                assert_eq!(arm, offsets[sample]);
                samples += 1;
                Ok(variation(frame, arm))
            },
        )
        .unwrap();
    assert_eq!(samples, rotor.stations.len() * rule.points());

    // Independent midpoint azimuth integration of the attached native law.
    // Resolve lift/drag using normalized velocity, not production sin/cos(phi),
    // and use separately declared radial interval widths.
    let mut force = ZERO;
    let mut moment = ZERO;
    let mut profile = 0.0;
    let mut air_work = 0.0;
    let points = 8192;
    for (station, dr) in rotor.stations.iter().zip([0.15, 0.2, 0.2, 0.25]) {
        let r = station.r_over_r;
        for k in 0..points {
            let theta = std::f64::consts::TAU * (k as f64 + 0.5) / points as f64;
            let arm = [0.0, r * theta.cos(), r * theta.sin()];
            let tangent = [0.0, -theta.sin(), theta.cos()];
            let delta = variation(frame, arm);
            let ua = 9.0 + delta[0];
            let ut = OMEGA * r + dot(delta, tangent);
            let speed = ua.hypot(ut);
            let alpha = station.beta_rad - ua.atan2(ut);
            assert!(alpha.abs() < 0.25);
            let cl = std::f64::consts::TAU * (alpha + 2.0 * rotor.camber_ratio);
            let cd = 0.012 + 0.8 * alpha * alpha;
            let q = 0.5 * RHO * speed * speed * station.chord_m * 3.0 * dr / points as f64;
            let axial = q * (cl * ut - cd * ua) / speed;
            let tangential = q * (cl * ua + cd * ut) / speed;
            let f = [axial, -tangential * tangent[1], -tangential * tangent[2]];
            let m = cross(arm, f);
            for j in 0..3 {
                force[j] += f[j];
                moment[j] += m[j];
            }
            profile += q * cd * speed;
            air_work += dot(f, delta);
        }
    }
    vector_near(load.force_n, force, 4e-10);
    vector_near(load.hub_moment_nm, [0.0, moment[1], moment[2]], 4e-10);
    near(load.torque_nm, -moment[0], 4e-10);
    near(load.profile_dissipation_w, profile, 3e-9);
    near(load.incident_variation_work_w, air_work, 3e-9);
    near(load.power_residual_w, 0.0, 3e-9);
    near(
        load.kinematic_power_w + load.induced_work_w + profile,
        0.0,
        3e-9,
    );
    assert!(load.hub_moment_nm[2].abs() > 1.0);
    assert!(load.incident_variation_work_w.abs() > 1.0);
    // A hub-mean-only flow cannot create this spatial loading couple.
    let mean = evaluate(64, SpinDirection::AlongThrust, [7.0, 0.0, 0.0], ZERO);
    vector_near(mean.hub_moment_nm, ZERO, 1e-12);
}

#[test]
fn incident_variation_is_frame_covariant_and_zero_variation_keeps_the_original_result() {
    let rotor = rotor();
    let rule = AzimuthQuadrature::new(32).unwrap();
    let frame = DiskFrame::new(AXIS, RADIAL, SpinDirection::AgainstThrust).unwrap();
    let rotated =
        DiskFrame::new(rotate(AXIS), rotate(RADIAL), SpinDirection::AgainstThrust).unwrap();
    let induction = HarmonicInduction {
        mean_mps: 2.0,
        first_mps: [0.2, -0.1],
    };
    let velocity = [7.0, 1.0, -0.5];
    let rate = [0.1, 0.2, -0.3];
    let evaluate = |frame, velocity, rate| {
        rule.evaluate_with_incident_variation(
            &rotor,
            RHO,
            frame,
            velocity,
            rate,
            induction,
            OMEGA,
            |_, arm| Ok(variation(frame, arm)),
        )
        .unwrap()
    };
    let a = evaluate(frame, velocity, rate);
    let b = evaluate(rotated, rotate(velocity), rotate(rate));
    vector_near(b.force_n, rotate(a.force_n), 4e-10);
    vector_near(b.hub_moment_nm, rotate(a.hub_moment_nm), 4e-10);
    near(b.torque_nm, a.torque_nm, 4e-10);
    near(
        b.incident_variation_work_w,
        a.incident_variation_work_w,
        4e-10,
    );
    near(b.power_residual_w, 0.0, 3e-9);
    let original = rule
        .evaluate_first_harmonic(&rotor, RHO, frame, velocity, rate, induction, OMEGA)
        .unwrap();
    let zero = rule
        .evaluate_with_incident_variation(
            &rotor,
            RHO,
            frame,
            velocity,
            rate,
            induction,
            OMEGA,
            |_, _| Ok(ZERO),
        )
        .unwrap();
    assert_eq!(zero, original);
    assert_eq!(zero.incident_variation_work_w, 0.0);
}

#[test]
fn invalid_or_failed_spatial_air_samples_refuse_the_whole_load_evaluation() {
    let rotor = rotor();
    let rule = AzimuthQuadrature::new(16).unwrap();
    let frame = DiskFrame::new(AXIS, RADIAL, SpinDirection::AlongThrust).unwrap();
    let induction = HarmonicInduction {
        mean_mps: 2.0,
        first_mps: [0.0; 2],
    };
    let last = rule.points() * rotor.stations.len() - 1;
    for (delta, code) in [
        ([f64::NAN, 0.0, 0.0], "oblique-rotor-invalid"),
        ([-100.0, 0.0, 0.0], "oblique-rotor-section-domain"),
    ] {
        let result = rule.evaluate_with_incident_variation(
            &rotor,
            RHO,
            frame,
            [7.0, 0.0, 0.0],
            ZERO,
            induction,
            OMEGA,
            |sample, _| Ok(if sample == last { delta } else { ZERO }),
        );
        assert_eq!(result.unwrap_err().code, code);
    }
    let refusal = fs_airscrew::Refusal {
        code: "incident-field-unavailable",
        message: "last physical sample refused".into(),
        ranked_repairs: Vec::new(),
    };
    let result = rule.evaluate_with_incident_variation(
        &rotor,
        RHO,
        frame,
        [7.0, 0.0, 0.0],
        ZERO,
        induction,
        OMEGA,
        |sample, _| {
            if sample == last {
                Err(refusal.clone())
            } else {
                Ok(ZERO)
            }
        },
    );
    assert_eq!(result.unwrap_err(), refusal);
}
