//! G1/G3: incoming sound drives original shell motion and transfers real work.
use super::*;
use crate::bem_loaded_modes::open_screen::molecular::tests as fixture;
use crate::bem_loaded_modes::open_screen::tests as screen_fixture;
use fs_bem::helmholtz::bc::{IncidentSource, Monopole, PlaneWave};

fn options() -> MolecularOpenScreenIncidentOptions {
    let original = fixture::options();
    MolecularOpenScreenIncidentOptions {
        boundary: MolecularScreenIlluminationLimits {
            boundary: original.boundary,
            maximum_incident_cases: 8,
            maximum_incident_jump_entries: 2048,
            maximum_source_point_products: 1_000_000,
            minimum_source_clearance_m: 0.001,
            maximum_live_bytes: 8_000_000,
        },
        maximum_load_cases: 8,
        maximum_response_bytes: 200_000,
        maximum_work: 2_000_000,
        physical_residual_tolerance: 1e-10,
        power_relative_tolerance: 1e-10,
    }
}
fn close(a: C64, b: C64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance * a.abs().max(b.abs()) + 1e-24,
        "{a:?} != {b:?}"
    );
}
fn plane(amplitude: C64) -> IncidentField<'static> {
    IncidentField::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, -1.0],
        amplitude_pa: amplitude,
    })
}

#[test]
fn sound_only_shell_closes_signed_work_and_keeps_full_loaded_mechanics() {
    let exterior = fixture::gas_exterior(true, 101325.0);
    let source = fixture::source(&exterior);
    let prepared = source
        .prepare(800.0, screen_fixture::boundary_options(), || false)
        .unwrap();
    let zero = [C64::ZERO; 2];
    let external = [C64::new(0.02, 0.01), C64::new(-0.01, 0.03)];
    let incident = plane(C64::new(2.0, 0.4));
    let batch = prepared
        .harmonic_molecular_incident_response_batch(
            &exterior,
            450.0,
            &[incident, incident],
            &[&zero, &external],
            options(),
            || false,
        )
        .unwrap();
    assert!(
        batch.responses()[0]
            .motion
            .modal_velocity
            .iter()
            .any(|v| v.abs() > 0.0)
    );
    let only = &batch.responses()[0].power;
    assert_eq!(only.external_mechanical_input_w, 0.0);
    assert!(only.structural_viscosity_w > 0.0);
    assert!(only.total_exterior_w < 0.0);
    assert!(
        (only.total_exterior_w + only.structural_viscosity_w).abs()
            < 1e-10 * only.apparent_power_scale_w
    );
    for (case, force) in [&zero[..], &external[..]].into_iter().enumerate() {
        let r = &batch.responses()[case];
        let effective: Vec<_> = force
            .iter()
            .zip(&r.blocked_acoustic_force)
            .map(|(&a, &b)| a + b)
            .collect();
        // Independently use the existing force-driven owner with the physical
        // blocked force. It builds its own original molecular boundary bank.
        let direct = prepared
            .harmonic_molecular_response_batch(
                &exterior,
                450.0,
                &[&effective],
                fixture::options(),
                || false,
            )
            .unwrap();
        for i in 0..2 {
            close(
                r.motion.modal_displacement[i],
                direct.responses()[0].modal_displacement[i],
                5e-11,
            );
            close(
                r.total_fluid_reaction[i],
                r.blocked_acoustic_force[i] + r.motion.fluid_reaction[i],
                5e-11,
            );
        }
        let direct_work = -0.5
            * r.total_fluid_reaction
                .iter()
                .zip(&r.motion.modal_velocity)
                .map(|(&f, &v)| (f * v.conj()).re)
                .sum::<f64>();
        let p = &r.power;
        assert!((p.total_exterior_w - direct_work).abs() < 1e-11 * p.apparent_power_scale_w);
        assert!(
            (p.total_exterior_w - p.motion_exterior_w - p.blocked_pressure_work_w).abs()
                < 1e-10 * p.apparent_power_scale_w
        );
        assert!(p.relative_balance < 1e-10);
    }
    let points = [[0.04, 0.04, 0.12], [0.04, 0.04, -0.12]];
    let receiver = batch
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let observed = receiver.evaluate(0, || false).unwrap();
    let same = batch
        .radiation()
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap()
        .evaluate(0, &batch.responses()[0].motion.modal_velocity, || false)
        .unwrap();
    assert_eq!(observed.scalar.pressure, same.scalar.pressure);
    assert_eq!(observed.particle_velocity_m_s, same.particle_velocity_m_s);
    assert!(
        observed
            .scalar
            .pressure
            .iter()
            .all(|z| finite(*z) && z.abs() > 0.0)
    );
    assert_eq!(receiver.points_m(), points);
}

#[test]
fn coherent_mixed_sources_preserve_motion_fields_and_refuse_partial_requests() {
    let exterior = fixture::gas_exterior(true, 101325.0);
    let source = fixture::source(&exterior);
    let prepared = source
        .prepare(800.0, screen_fixture::boundary_options(), || false)
        .unwrap();
    let first = IncidentSource::PlaneWave(PlaneWave {
        direction: [0.0, 0.0, -1.0],
        amplitude_pa: C64::new(1.5, 0.2),
    });
    let second = IncidentSource::Monopole(Monopole {
        position_m: [0.1, 0.1, -0.4],
        volume_velocity_m3_s: C64::new(2e-4, -1e-4),
    });
    let both = [first, second];
    let a = [first];
    let b = [second];
    let incidents = [
        IncidentField::Sources(&a),
        IncidentField::Sources(&b),
        IncidentField::Sources(&both),
    ];
    let zero = [C64::ZERO; 2];
    let forces: [&[C64]; 3] = [&zero, &zero, &zero];
    let batch = prepared
        .harmonic_molecular_incident_response_batch(
            &exterior,
            450.0,
            &incidents,
            &forces,
            options(),
            || false,
        )
        .unwrap();
    for i in 0..2 {
        close(
            batch.responses()[2].motion.modal_velocity[i],
            batch.responses()[0].motion.modal_velocity[i]
                + batch.responses()[1].motion.modal_velocity[i],
            5e-11,
        );
    }
    let points = [[0.04, 0.04, 0.12], [0.04, 0.04, -0.12]];
    let receiver = batch
        .prepare_receivers(
            &points,
            0.001,
            near_field::Options::default(),
            1_000_000,
            || false,
        )
        .unwrap();
    let evaluations = receiver.kernel_evaluations();
    let fields: Vec<_> = (0..3)
        .map(|i| receiver.evaluate(i, || false).unwrap())
        .collect();
    for r in 0..2 {
        close(
            fields[2].scalar.pressure[r],
            fields[0].scalar.pressure[r] + fields[1].scalar.pressure[r],
            5e-11,
        );
        for axis in 0..3 {
            close(
                fields[2].particle_velocity_m_s[r][axis],
                fields[0].particle_velocity_m_s[r][axis] + fields[1].particle_velocity_m_s[r][axis],
                5e-11,
            );
        }
    }
    assert_eq!(evaluations, receiver.kernel_evaluations());
    assert_eq!(
        fields[2].scalar.pressure,
        receiver.evaluate(2, || false).unwrap().scalar.pressure
    );
    assert!(receiver.evaluate(3, || false).is_err());
    assert!(matches!(
        receiver.evaluate(0, || true),
        Err(BemModesError::Cancelled)
    ));
    assert!(
        prepared
            .harmonic_molecular_incident_response_batch(
                &exterior,
                450.0,
                &incidents,
                &[&zero],
                options(),
                || false
            )
            .is_err()
    );
    let mut limited = options();
    limited.maximum_response_bytes = 1;
    assert!(
        prepared
            .harmonic_molecular_incident_response_batch(
                &exterior,
                450.0,
                &incidents,
                &forces,
                limited,
                || false
            )
            .is_err()
    );
    assert!(matches!(
        prepared.harmonic_molecular_incident_response_batch(
            &exterior,
            450.0,
            &incidents,
            &forces,
            options(),
            || true
        ),
        Err(BemModesError::Cancelled)
    ));
    let zero_incident = plane(C64::ZERO);
    let quiet = prepared
        .harmonic_molecular_incident_response_batch(
            &exterior,
            450.0,
            &[zero_incident],
            &[&zero],
            options(),
            || false,
        )
        .unwrap();
    assert_eq!(quiet.responses()[0].power.total_exterior_w, 0.0);
    assert_eq!(quiet.responses()[0].power.relative_balance, 0.0);
    assert!(
        quiet.responses()[0]
            .motion
            .modal_velocity
            .iter()
            .all(|&v| v == C64::ZERO)
    );
}
