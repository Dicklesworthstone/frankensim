use super::*;
use crate::bem_loaded_modes::molecular::tests as fixture;
use fs_bem::helmholtz::bc::{IncidentSource, Monopole, PlaneWave};

fn limits() -> MolecularIncidentLimits {
    MolecularIncidentLimits {
        minimum_source_clearance_m: 0.01,
        maximum_source_point_products: 100_000,
        maximum_live_bytes: 2_000_000,
        maximum_panel_products: 1_000_000,
    }
}
fn plane() -> PlaneWave {
    PlaneWave {
        direction: [0.0, 0.0, 1.0],
        amplitude_pa: C64::new(0.8, -0.2),
    }
}
fn close(a: C64, b: C64, scale: f64) {
    assert!((a - b).abs() <= 3e-9 * scale.max(1e-20), "{a:?} != {b:?}");
}

#[test]
fn finite_source_and_plane_drive_loaded_wall_motion_and_one_total_field() {
    let exterior = fixture::exterior(true);
    let (source, boundary) = fixture::partly_driven(&exterior);
    let receivers = fixture::receivers();
    let prepared = source
        .prepare_molecular_harmonic(
            &exterior,
            700.0,
            &receivers,
            0.01,
            Some(boundary.clone()),
            fixture::limits(),
            fixture::options(),
            || false,
        )
        .unwrap();
    let finite = [Monopole {
        position_m: [0.35, -0.12, 0.27],
        volume_velocity_m3_s: C64::new(2e-4, 5e-5),
    }];
    let sources = [
        IncidentSource::PlaneWave(plane()),
        IncidentSource::Monopole(finite[0]),
    ];
    let incidents = [
        IncidentField::PlaneWave(plane()),
        IncidentField::Monopoles(&finite),
        IncidentField::Sources(&sources),
    ];
    let first = [C64::new(0.2, 0.1), C64::new(-0.1, 0.05)];
    let second = [C64::new(-0.08, 0.04), C64::new(0.03, -0.02)];
    let combined = [first[0] + second[0], first[1] + second[1]];
    let responses = prepared
        .harmonic_incident_response_batch(
            &incidents,
            &[&first, &second, &combined],
            limits(),
            || false,
        )
        .unwrap();
    let response = &responses[2];
    for mode in 0..2 {
        close(
            response.motion.modal_velocity[mode],
            responses[0].motion.modal_velocity[mode] + responses[1].motion.modal_velocity[mode],
            response.motion.modal_velocity[mode].abs(),
        );
        close(
            response.blocked_acoustic_force[mode],
            responses[0].blocked_acoustic_force[mode] + responses[1].blocked_acoustic_force[mode],
            response.blocked_acoustic_force[mode].abs(),
        );
        close(
            response.total_fluid_reaction[mode],
            response.blocked_acoustic_force[mode] + response.motion.fluid_reaction[mode],
            response.blocked_acoustic_force[mode].abs()
                + response.motion.fluid_reaction[mode].abs(),
        );
    }
    let pscale = response.total_receivers[0].pressure_pa.abs();
    close(
        response.total_receivers[0].pressure_pa,
        responses[0].total_receivers[0].pressure_pa + responses[1].total_receivers[0].pressure_pa,
        pscale,
    );
    for axis in 0..3 {
        close(
            response.total_receivers[0].particle_velocity_m_s[axis],
            responses[0].total_receivers[0].particle_velocity_m_s[axis]
                + responses[1].total_receivers[0].particle_velocity_m_s[axis],
            pscale / 400.0,
        );
    }

    // Independent original-owner solve of the ONE combined outgoing trace:
    // scattered drive = actual source motion - incident v - Y*incident p.
    // It does not use either the blocked/motion boundary or receiver bank.
    let beta = boundary.normalized_admittance(700.0, || false).unwrap();
    let rho_c = source.medium().density * source.medium().sound_speed;
    let mut drive = Vec::new();
    let mut incident_pressure = Vec::new();
    let mut incident_velocity = Vec::new();
    for panel in 0..source.panel_count() {
        let sample = exterior
            .incident_pressure_velocity_at(
                700.0,
                incidents[2],
                source.surface().centroids()[panel],
                || false,
            )
            .unwrap();
        let velocity = (0..3).fold(C64::ZERO, |sum, axis| {
            sum + sample.velocity_m_s[axis].scale(source.surface().normals()[panel][axis])
        });
        let skin_velocity = (0..2).fold(C64::ZERO, |sum, mode| {
            sum + response.motion.modal_velocity[mode]
                .scale(source.normal_velocity_weights()[mode][panel])
        });
        drive.push(skin_velocity - velocity - beta[panel] * sample.pressure_pa.scale(1.0 / rho_c));
        incident_pressure.push(sample.pressure_pa);
        incident_velocity.push(velocity);
        close(
            response.total_surface_normal_velocity_m_s[panel],
            skin_velocity
                - beta[panel] * response.total_surface_pressure_pa[panel].scale(1.0 / rho_c),
            skin_velocity.abs()
                + (beta[panel] * response.total_surface_pressure_pa[panel]).abs() / rho_c,
        );
    }
    let direct = exterior
        .solve_boundary_batch(
            source.surface(),
            700.0,
            &[&drive],
            Admittance::Local {
                beta: &beta,
                declared_active: false,
            },
            source.settings.formulation,
            &source.settings.boundary,
            || false,
        )
        .unwrap()
        .remove(0)
        .radiation;
    for panel in 0..source.panel_count() {
        close(
            response.total_surface_pressure_pa[panel],
            incident_pressure[panel] + direct.pressure_pa()[panel],
            response.total_surface_pressure_pa[panel].abs(),
        );
        close(
            response.total_surface_normal_velocity_m_s[panel],
            incident_velocity[panel] + direct.velocity_m_s()[panel],
            incident_velocity[panel].abs() + direct.velocity_m_s()[panel].abs(),
        );
    }
    let point = receivers[0].position_m;
    let incident = exterior
        .incident_pressure_velocity_at(700.0, incidents[2], point, || false)
        .unwrap();
    let outgoing = direct
        .observe(source.surface(), &[point], 0.01, || false)
        .unwrap()
        .remove(0);
    close(
        response.total_receivers[0].pressure_pa,
        incident.pressure_pa + outgoing.pressure_pa,
        pscale,
    );
    for axis in 0..3 {
        close(
            response.total_receivers[0].particle_velocity_m_s[axis],
            incident.velocity_m_s[axis]
                + outgoing.gradient_pa_per_m[axis]
                    * C64::new(0.0, -1.0 / (700.0 * source.medium().density)),
            pscale / rho_c,
        );
    }
    let power = response.power;
    let tolerance = 3e-9 * power.apparent_power_scale_w;
    assert!(power.kelvin_w > 0.0 && power.relaxation_w > 0.0 && power.total_wall_w > 0.0);
    assert!(power.relative_balance < 1e-9);
    assert!((power.blocked_exterior_w + power.blocked_wall_w).abs() < tolerance);
    assert!(
        (power.interference_exterior_w
            + power.interference_wall_w
            + power.blocked_acoustic_input_w)
            .abs()
            < tolerance
    );
    assert!(
        (power.total_exterior_w
            - power.blocked_exterior_w
            - power.motion_exterior_w
            - power.interference_exterior_w)
            .abs()
            < tolerance
    );
    assert!(
        (power.total_wall_w
            - power.blocked_wall_w
            - power.motion_wall_w
            - power.interference_wall_w)
            .abs()
            < tolerance
    );
    assert!(power.interference_wall_w.abs() > 1e-8 * power.apparent_power_scale_w);
    assert!((power.motion_exterior_w - response.motion.power.exterior_w).abs() < tolerance);
    assert!((power.motion_wall_w - response.motion.power.wall_w).abs() < tolerance);
}

#[test]
fn exact_force_cancellation_keeps_blocked_sound_and_failed_batches_leave_cache_reusable() {
    let exterior = fixture::exterior(true);
    let (source, boundary) = fixture::partly_driven(&exterior);
    let prepared = source
        .prepare_molecular_harmonic(
            &exterior,
            700.0,
            &fixture::receivers(),
            0.01,
            Some(boundary),
            fixture::limits(),
            fixture::options(),
            || false,
        )
        .unwrap();
    let incident = [IncidentField::PlaneWave(plane())];
    let zero = [C64::ZERO; 2];
    let mut calls = 0;
    let sound = prepared
        .harmonic_incident_response_batch(&incident, &[&zero], limits(), || {
            calls += 1;
            false
        })
        .unwrap()
        .remove(0);
    assert!(sound.motion.modal_velocity.iter().any(|&v| v != C64::ZERO));
    let opposing: Vec<_> = sound.blocked_acoustic_force.iter().map(|&f| -f).collect();
    let cancelled_motion = prepared
        .harmonic_incident_response_batch(&incident, &[&opposing], limits(), || false)
        .unwrap()
        .remove(0);
    assert!(
        cancelled_motion
            .motion
            .modal_velocity
            .iter()
            .all(|&v| v == C64::ZERO)
    );
    assert!(
        cancelled_motion
            .motion
            .elastic_memory_sqrt_j
            .iter()
            .all(|&e| e == C64::ZERO)
    );
    assert!(cancelled_motion.total_receivers[0].pressure_pa.abs() > 0.1);
    assert_eq!(cancelled_motion.power.kelvin_w, 0.0);
    assert_eq!(cancelled_motion.power.relaxation_w, 0.0);
    assert_eq!(cancelled_motion.power.interference_exterior_w, 0.0);
    assert_eq!(cancelled_motion.power.interference_wall_w, 0.0);
    assert!(cancelled_motion.power.total_wall_w > 0.0);

    for rejected in [
        MolecularIncidentLimits {
            maximum_live_bytes: 0,
            ..limits()
        },
        MolecularIncidentLimits {
            maximum_panel_products: 0,
            ..limits()
        },
        MolecularIncidentLimits {
            maximum_source_point_products: 1,
            ..limits()
        },
    ] {
        let mut cold_calls = 0;
        assert!(
            prepared
                .harmonic_incident_response_batch(&incident, &[&zero], rejected, || {
                    cold_calls += 1;
                    false
                })
                .is_err()
        );
        assert_eq!(cold_calls, 1, "resource admission must precede BEM work");
    }
    let mut cancelled_calls = 0;
    assert!(
        prepared
            .harmonic_incident_response_batch(&incident, &[&zero], limits(), || {
                cancelled_calls += 1;
                cancelled_calls >= calls - 2
            })
            .is_err()
    );
    let replay = prepared
        .harmonic_incident_response_batch(&incident, &[&zero], limits(), || false)
        .unwrap()
        .remove(0);
    assert_eq!(replay.motion.modal_velocity, sound.motion.modal_velocity);
    assert_eq!(replay.total_receivers, sound.total_receivers);
    assert_eq!(
        replay.total_surface_pressure_pa,
        sound.total_surface_pressure_pa
    );
    let inside = [Monopole {
        position_m: [0.0; 3],
        volume_velocity_m3_s: C64::ONE,
    }];
    assert!(
        prepared
            .harmonic_incident_response_batch(
                &[IncidentField::Monopoles(&inside)],
                &[&zero],
                limits(),
                || false,
            )
            .is_err()
    );
}
