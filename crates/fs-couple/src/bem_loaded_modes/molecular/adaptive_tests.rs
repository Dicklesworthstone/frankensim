//! G1/G3: observation quadrature changes neither original mechanics nor work.
use super::tests as fixture;
use super::*;
use fs_bem::helmholtz::{
    bc::{Admittance, IncidentField, Monopole, PlaneWave},
    molecular::MolecularScatteringOptions,
};

fn close(actual: C64, expected: C64, scale: f64) {
    assert!(
        (actual - expected).abs() < 3e-9 * scale.max(1e-25),
        "{actual:?} != {expected:?}"
    );
}

fn unchanged_mechanics(a: &MolecularLoadedResponse, b: &MolecularLoadedResponse) {
    assert_eq!(a.modal_displacement, b.modal_displacement);
    assert_eq!(a.modal_velocity, b.modal_velocity);
    assert_eq!(a.elastic_memory_sqrt_j, b.elastic_memory_sqrt_j);
    assert_eq!(a.fluid_reaction, b.fluid_reaction);
    assert_eq!(a.mechanical_residual, b.mechanical_residual);
    assert_eq!(a.memory_residual, b.memory_residual);
    let powers = |p: MolecularModesPower| {
        [
            p.input_w,
            p.kelvin_w,
            p.relaxation_w,
            p.exterior_w,
            p.wall_w,
            p.balance_w,
            p.apparent_power_scale_w,
            p.relative_balance,
        ]
    };
    assert_eq!(powers(a.power), powers(b.power));
}

fn drive(source: &BemLoadedModes, velocity: &[C64]) -> Vec<C64> {
    (0..source.panel_count())
        .map(|panel| {
            velocity
                .iter()
                .zip(source.normal_velocity_weights())
                .fold(C64::ZERO, |sum, (&v, weights)| {
                    sum + v.scale(weights[panel])
                })
        })
        .collect()
}

#[test]
fn adaptive_loaded_and_incident_receivers_keep_original_wall_motion_and_work() {
    let exterior = fixture::exterior(true);
    let (source, boundary) = fixture::partly_driven(&exterior);
    let omega = 700.0;
    let receivers: Vec<_> = [0, 13]
        .into_iter()
        .map(|panel| ReceiverSpec {
            id: format!("original skin receiver {panel}"),
            position_m: std::array::from_fn(|axis| {
                source.surface().centroids()[panel][axis]
                    + 0.025 * source.surface().normals()[panel][axis]
            }),
        })
        .collect();
    let points: Vec<_> = receivers
        .iter()
        .map(|receiver| receiver.position_m)
        .collect();
    let near = MolecularNearFieldOptions {
        quadrature: fs_bem::near_field::Options {
            relative_tolerance: 1e-8,
            maximum_kernel_evaluations: 4_000_000,
            ..fs_bem::near_field::Options::default()
        },
        maximum_row_bytes: 1_000_000,
    };
    let centroid = source
        .prepare_molecular_harmonic(
            &exterior,
            omega,
            &receivers,
            0.01,
            Some(boundary.clone()),
            fixture::limits(),
            fixture::options(),
            || false,
        )
        .unwrap();
    let adaptive = source
        .prepare_molecular_harmonic_adaptive(
            &exterior,
            omega,
            &receivers,
            0.01,
            Some(boundary.clone()),
            fixture::limits(),
            fixture::options(),
            near,
            || false,
        )
        .unwrap();
    assert!(centroid.near_field_options().is_none());
    assert!(adaptive.near_field_options().is_some());
    assert_eq!(
        centroid.radiation().impedance(),
        adaptive.radiation().impedance()
    );
    let quadrature = adaptive.radiation().receiver_quadrature().unwrap();
    assert_eq!(
        quadrature.source_fields.len(),
        source.mode_count() * receivers.len()
    );
    let row_evaluations = quadrature.kernel_evaluations;
    assert!(row_evaluations > 80 * source.panel_count() * receivers.len());
    let forces = [
        [C64::new(0.2, 0.1), C64::new(-0.1, 0.05)],
        [C64::new(-0.08, 0.04), C64::new(0.03, -0.02)],
    ];
    let columns = [&forces[0][..], &forces[1][..]];
    let original = centroid
        .harmonic_response_batch(&columns, || false)
        .unwrap();
    let actual = adaptive
        .harmonic_response_batch(&columns, || false)
        .unwrap();
    let beta = boundary.normalized_admittance(omega, || false).unwrap();
    for (a, b) in actual.iter().zip(&original) {
        unchanged_mechanics(a, b);
        assert!(b.receiver_quadrature_error_estimates.is_none());
        let estimates = a.receiver_quadrature_error_estimates.as_ref().unwrap();
        assert_eq!(estimates.len(), receivers.len());
        assert!(estimates.iter().all(|error| {
            error.pressure_pa.is_finite()
                && error.pressure_pa > 0.0
                && error
                    .particle_velocity_m_s
                    .iter()
                    .all(|e| e.is_finite() && *e >= 0.0)
        }));
        // Independent original-owner solve of this actual combined motion,
        // including original wall response; no unit-source receiver bank.
        let velocity = drive(&source, &a.modal_velocity);
        let direct = exterior
            .solve_boundary_batch(
                source.surface(),
                omega,
                &[&velocity],
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
        let observed = direct
            .observe_near_field(source.surface(), &points, 0.01, near, || false)
            .unwrap();
        for i in 0..receivers.len() {
            let scale = observed.scalar.pressure[i].abs();
            close(
                a.receivers[i].pressure_pa,
                observed.scalar.pressure[i],
                scale,
            );
            for axis in 0..3 {
                close(
                    a.receivers[i].particle_velocity_m_s[axis],
                    observed.particle_velocity_m_s[i][axis],
                    scale / 400.0,
                );
            }
        }
    }
    assert_eq!(
        adaptive
            .radiation()
            .receiver_quadrature()
            .unwrap()
            .kernel_evaluations,
        row_evaluations
    );

    let monopoles = [Monopole {
        position_m: [0.35, -0.12, 0.27],
        volume_velocity_m3_s: C64::new(2e-4, 5e-5),
    }];
    let incidents = [
        IncidentField::PlaneWave(PlaneWave {
            direction: [0.0, 0.0, 1.0],
            amplitude_pa: C64::new(0.8, -0.2),
        }),
        IncidentField::Monopoles(&monopoles),
    ];
    let limits = MolecularIncidentLimits {
        minimum_source_clearance_m: 0.01,
        maximum_source_point_products: 100_000,
        maximum_live_bytes: 2_000_000,
        maximum_panel_products: 1_000_000,
    };
    let original = centroid
        .harmonic_incident_response_batch(&incidents, &columns, limits, || false)
        .unwrap();
    let actual = adaptive
        .harmonic_incident_response_batch(&incidents, &columns, limits, || false)
        .unwrap();
    let blocked = exterior
        .solve_scattering_batch(
            source.surface(),
            omega,
            &incidents,
            Admittance::Local {
                beta: &beta,
                declared_active: false,
            },
            source.settings.formulation,
            &MolecularScatteringOptions {
                boundary: source.settings.boundary,
                minimum_source_clearance_m: limits.minimum_source_clearance_m,
                maximum_source_point_products: limits.maximum_source_point_products,
            },
            || false,
        )
        .unwrap();
    // The original blocked batch admits one reusable source-bound row bank.
    let blocked_rows = blocked[0]
        .solution()
        .scattered
        .prepare_near_field(source.surface(), &points, 0.01, near, || false)
        .unwrap();
    let blocked_evaluations = blocked_rows.kernel_evaluations();
    let rho_c = source.medium().density * source.medium().sound_speed;
    for case in 0..actual.len() {
        let (a, b) = (&actual[case], &original[case]);
        unchanged_mechanics(&a.motion, &b.motion);
        assert_eq!(a.blocked_acoustic_force, b.blocked_acoustic_force);
        assert_eq!(a.total_fluid_reaction, b.total_fluid_reaction);
        assert_eq!(a.total_surface_pressure_pa, b.total_surface_pressure_pa);
        assert_eq!(
            a.total_surface_normal_velocity_m_s,
            b.total_surface_normal_velocity_m_s
        );
        assert_eq!(a.power.balance_w, b.power.balance_w);
        assert_eq!(a.power.total_exterior_w, b.power.total_exterior_w);
        assert_eq!(a.power.total_wall_w, b.power.total_wall_w);
        assert!(b.total_receiver_quadrature_error_estimates.is_none());
        let total_errors = a
            .total_receiver_quadrature_error_estimates
            .as_ref()
            .unwrap();
        let motion_errors = a
            .motion
            .receiver_quadrature_error_estimates
            .as_ref()
            .unwrap();
        let blocked_fields = blocked[case]
            .total_pressure_velocity_with_near_field(&blocked_rows, || false)
            .unwrap();
        for i in 0..receivers.len() {
            assert_eq!(
                total_errors[i].pressure_pa,
                blocked_fields.scalar.quadrature_error_estimate_pa[i]
                    + motion_errors[i].pressure_pa
            );
            for axis in 0..3 {
                assert_eq!(
                    total_errors[i].particle_velocity_m_s[axis],
                    blocked_fields.quadrature_error_estimate_m_s[i][axis]
                        + motion_errors[i].particle_velocity_m_s[axis]
                );
            }
        }
        // Solve one independently combined outgoing drive. Add the analytic
        // incident at the final physical receiver exactly once.
        let incident = blocked[case].solution();
        let skin = drive(&source, &a.motion.modal_velocity);
        let outgoing: Vec<_> = (0..source.panel_count())
            .map(|panel| {
                skin[panel]
                    - incident.incident_normal_velocity_m_s[panel]
                    - beta[panel] * incident.incident_pressure_pa[panel].scale(1.0 / rho_c)
            })
            .collect();
        let direct = exterior
            .solve_boundary_batch(
                source.surface(),
                omega,
                &[&outgoing],
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
        let observed = direct
            .observe_near_field(source.surface(), &points, 0.01, near, || false)
            .unwrap();
        for (i, &point) in points.iter().enumerate() {
            let incident = exterior
                .incident_pressure_velocity_at(omega, incidents[case], point, || false)
                .unwrap();
            let pressure = observed.scalar.pressure[i] + incident.pressure_pa;
            let scale = observed.scalar.pressure[i].abs() + incident.pressure_pa.abs();
            close(a.total_receivers[i].pressure_pa, pressure, scale);
            for axis in 0..3 {
                close(
                    a.total_receivers[i].particle_velocity_m_s[axis],
                    observed.particle_velocity_m_s[i][axis] + incident.velocity_m_s[axis],
                    scale / 400.0,
                );
            }
        }
    }
    assert_eq!(blocked_rows.kernel_evaluations(), blocked_evaluations);
    assert_eq!(
        adaptive
            .radiation()
            .receiver_quadrature()
            .unwrap()
            .kernel_evaluations,
        row_evaluations
    );
}
