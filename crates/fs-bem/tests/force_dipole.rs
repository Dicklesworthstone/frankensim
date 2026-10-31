//! G0/G1/G3: physical force [N] through the original diffraction and observation owners.
use fs_bem::helmholtz::bc::{
    Admittance, ForceDipole, IncidentField, IncidentSource, IncidentTraces, Monopole, incident,
};
use fs_bem::helmholtz::{DirectivityProjectionOptions, Formulation, IterativeOptions, Medium};
use fs_bem::near_field::Options;
use fs_bem::observation_transfers::{
    BemObservation, CompactSourceSpec, ForceDipoleSpec, HarmonicField, HarmonicFieldStatistics,
    HarmonicWorkBudget, IterativeBoundaryObservation, IterativeDirectivityOptions,
    IterativeObservationOptions, LocalBoundary, LocalBoundaryObservation, LocalWall,
    LowerBandWindow, MAX_OBSERVATION_FFT, PointSourceSpec, ProducerUncertainty, ReceiverSpec,
    ScatteringObservation, SourceId, TransferDeclaration, bake_observation_transfers,
};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::ContentHash;
use fs_exec::CancelGate;
use fs_math::{c64::C64, det};
use std::f64::consts::TAU;

const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 340.0,
};
const Q_POSITION: [f64; 3] = [0.02, 0.01, 0.09];
const F_POSITION: [f64; 3] = [-0.03, -0.01, 0.07];
const DIRECTION: [f64; 3] = [0.6, 0.0, 0.8];

fn close(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-18),
        "actual {actual:?}, expected {expected:?}"
    );
}
fn specs() -> [CompactSourceSpec; 2] {
    [
        CompactSourceSpec::VolumeFlow(PointSourceSpec {
            id: "flow".into(),
            position_m: Q_POSITION,
        }),
        CompactSourceSpec::ForceDipole(ForceDipoleSpec {
            id: "force-on-fluid".into(),
            position_m: F_POSITION,
            direction: DIRECTION,
        }),
    ]
}
fn receivers() -> [ReceiverSpec; 1] {
    [ReceiverSpec {
        id: "shadow".into(),
        position_m: [-0.03, 0.02, -0.08],
    }]
}
fn illumination(amplitudes: [C64; 2]) -> [IncidentSource; 2] {
    [
        IncidentSource::Monopole(Monopole {
            position_m: Q_POSITION,
            volume_velocity_m3_s: amplitudes[0],
        }),
        IncidentSource::ForceDipole(ForceDipole {
            position_m: F_POSITION,
            force_on_fluid_n: DIRECTION.map(|d| amplitudes[1].scale(d)),
        }),
    ]
}
fn observation(surface: &SpherePanels, medium: Medium) -> ScatteringObservation<'_> {
    ScatteringObservation::new_mixed(
        surface,
        medium,
        &specs(),
        &receivers(),
        TAU * 1600.0,
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap()
}
fn rigid(surface: &SpherePanels) -> LocalBoundary {
    LocalBoundary::try_new(
        MEDIUM,
        vec![LocalWall::Rigid],
        vec![0; surface.areas().len()],
    )
    .unwrap()
}
fn declaration(
    hash: ContentHash,
    sources: &[SourceId],
    receivers: &[ReceiverSpec],
    delays: Vec<f64>,
) -> TransferDeclaration {
    TransferDeclaration {
        producer_hash: hash,
        sources: sources.to_vec(),
        receivers: receivers.to_vec(),
        fs_hz: 8000.0,
        geometric_delays_s: delays,
        pass_band_upper_hz: 800.0,
        taper_stop_hz: 1600.0,
        filter_latency_samples: 256,
        fft_lengths: vec![512, 1024, 2048, 4096],
        tail_energy_tolerance: 1e-10,
        dtft_relative_tolerance: 1e-4,
        producer_uncertainty: ProducerUncertainty::Unresolved { source_hash: hash },
        maximum_frequency_entries: 8 * MAX_OBSERVATION_FFT,
    }
}

// Two separated volume injections converge to a FORCE source only with
// Q_plus = i F_axis/(omega rho epsilon), Q_minus = -Q_plus. This independent
// source construction checks sign, dimensions, near-field terms and scattering.
fn monopole_pairs(force: ForceDipole, k: f64, epsilon: f64) -> Vec<Monopole> {
    let mut sources = Vec::new();
    for axis in 0..3 {
        let flow = force.force_on_fluid_n[axis]
            * C64::new(
                0.0,
                1.0 / (k * MEDIUM.sound_speed * MEDIUM.density * epsilon),
            );
        let mut plus = force.position_m;
        let mut minus = plus;
        plus[axis] += epsilon / 2.0;
        minus[axis] -= epsilon / 2.0;
        sources.push(Monopole {
            position_m: plus,
            volume_velocity_m3_s: flow,
        });
        sources.push(Monopole {
            position_m: minus,
            volume_velocity_m3_s: -flow,
        });
    }
    sources
}

#[test]
fn g1_force_matches_independent_dipole_limit_euler_gradient_and_pressure_null() {
    let k = 0.7;
    let forces = [ForceDipole {
        position_m: [0.1, -0.2, 0.05],
        force_on_fluid_n: [C64::new(0.8, -0.3), C64::new(-0.4, 0.2), C64::new(0.1, 0.5)],
    }];
    let field = IncidentField::ForceDipoles(&forces);
    for point in [[0.3, 0.5, 0.2], [3.0, -1.0, 0.4]] {
        let exact = field.pressure_velocity_at(k, MEDIUM, point).unwrap();
        let mut errors = Vec::new();
        for epsilon in [0.02, 0.01, 0.005] {
            let pair = monopole_pairs(forces[0], k, epsilon);
            let limit = IncidentField::Monopoles(&pair)
                .pressure_velocity_at(k, MEDIUM, point)
                .unwrap();
            let velocity_error = det::sqrt(
                (0..3)
                    .map(|i| (limit.velocity_m_s[i] - exact.velocity_m_s[i]).norm_sq())
                    .sum::<f64>()
                    / exact.velocity_m_s.iter().map(|v| v.norm_sq()).sum::<f64>(),
            );
            errors
                .push(velocity_error.max((limit.pressure_pa / exact.pressure_pa - C64::ONE).abs()));
        }
        assert!(
            errors[1] < 0.27 * errors[0] && errors[2] < 0.27 * errors[1],
            "{errors:?}"
        );
        assert!(errors[2] < 1e-4, "{errors:?}");
        for axis in 0..3 {
            let h = 1e-5;
            let mut plus = point;
            let mut minus = point;
            plus[axis] += h;
            minus[axis] -= h;
            let derivative = (field
                .pressure_velocity_at(k, MEDIUM, plus)
                .unwrap()
                .pressure_pa
                - field
                    .pressure_velocity_at(k, MEDIUM, minus)
                    .unwrap()
                    .pressure_pa)
                .scale(0.5 / h);
            close(
                derivative / C64::new(0.0, k * MEDIUM.sound_speed * MEDIUM.density),
                exact.velocity_m_s[axis],
                2e-8,
            );
        }
        let dense = field
            .pressure_velocity_at(
                k,
                Medium {
                    density: 3.6,
                    ..MEDIUM
                },
                point,
            )
            .unwrap();
        assert_eq!(
            dense.pressure_pa, exact.pressure_pa,
            "force pressure must not gain density"
        );
        for axis in 0..3 {
            close(
                dense.velocity_m_s[axis].scale(3.0),
                exact.velocity_m_s[axis],
                2e-14,
            );
        }
    }
    let axial = [ForceDipole {
        position_m: [0.0; 3],
        force_on_fluid_n: [C64::ONE, C64::ZERO, C64::ZERO],
    }];
    for radius in [0.01, 0.5, 20.0] {
        let null = IncidentField::ForceDipoles(&axial)
            .pressure_velocity_at(k, MEDIUM, [0.0, radius, 0.0])
            .unwrap();
        assert_eq!(null.pressure_pa, C64::ZERO);
        assert!(
            null.velocity_m_s[0].abs() > 0.0,
            "angular derivative was lost"
        );
        assert_eq!(null.velocity_m_s[1], C64::ZERO);
        assert_eq!(null.velocity_m_s[2], C64::ZERO);
    }
}

#[test]
fn g3_mixed_force_scattering_and_harmonic_observation_keep_si_cross_terms() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let omega = TAU * 700.0;
    let k = omega / MEDIUM.sound_speed;
    let amplitudes = [C64::new(0.002, 0.001), C64::new(0.4, -0.2)];
    let sources = illumination(amplitudes);
    let fields = [
        IncidentField::Sources(&sources[..1]),
        IncidentField::Sources(&sources[1..]),
        IncidentField::Sources(&sources),
    ];
    let solved = incident::solve_scattering_batch(
        &surface,
        k,
        MEDIUM,
        Admittance::Rigid,
        &fields,
        1e-3,
        Formulation::BurtonMiller,
    )
    .unwrap();
    let points = receivers().map(|r| r.position_m);
    let total: Vec<_> = solved
        .iter()
        .map(|s| {
            s.total_pressure_velocity_at(&points, 1e-3, Options::default())
                .unwrap()
        })
        .collect();
    close(
        total[0].scalar.pressure[0] + total[1].scalar.pressure[0],
        total[2].scalar.pressure[0],
        2e-11,
    );
    for axis in 0..3 {
        close(
            total[0].particle_velocity_m_s[0][axis] + total[1].particle_velocity_m_s[0][axis],
            total[2].particle_velocity_m_s[0][axis],
            2e-11,
        );
    }
    let IncidentSource::ForceDipole(force) = sources[1] else {
        unreachable!()
    };
    let pairs = monopole_pairs(force, k, 2e-5);
    let limit = incident::solve_scattering(
        &surface,
        k,
        MEDIUM,
        Admittance::Rigid,
        IncidentField::Monopoles(&pairs),
        1e-3,
        Formulation::BurtonMiller,
    )
    .unwrap();
    let pair_field = limit
        .total_pressure_velocity_at(&points, 1e-3, Options::default())
        .unwrap();
    close(
        pair_field.scalar.pressure[0],
        total[1].scalar.pressure[0],
        3e-6,
    );
    for axis in 0..3 {
        close(
            pair_field.particle_velocity_m_s[0][axis],
            total[1].particle_velocity_m_s[0][axis],
            3e-6,
        );
    }
    close(
        limit.total_far_field([0.0, 0.0, -1.0]).unwrap(),
        solved[1].total_far_field([0.0, 0.0, -1.0]).unwrap(),
        3e-6,
    );
    assert_eq!(
        solved[1].total_far_field([0.0, 0.0, -1.0]).unwrap(),
        solved[1].total_far_field([0.0, 0.0, -1.0 + 5e-11]).unwrap()
    );

    let owner = observation(&surface, MEDIUM);
    assert!(matches!(&owner.sources()[0], SourceId::PointVolumeFlow(_)));
    assert!(matches!(&owner.sources()[1], SourceId::PointForce(_)));
    let bank = owner.prepare_harmonic(omega, || false).unwrap();
    let mut received = [HarmonicField::default()];
    bank.evaluate_coherent(&amplitudes, &mut received, || false)
        .unwrap();
    close(received[0].pressure_pa, total[2].scalar.pressure[0], 2e-11);
    let mut statistics = [HarmonicFieldStatistics::default()];
    bank.evaluate_coherence_modes(
        &amplitudes,
        1,
        HarmonicWorkBudget {
            maximum_products: 1000,
        },
        &mut statistics,
        || false,
    )
    .unwrap();
    for axis in 0..3 {
        close(
            received[0].particle_velocity_m_s[axis],
            total[2].particle_velocity_m_s[0][axis],
            2e-11,
        );
        close(
            statistics[0].complex_intensity_w_m2[axis],
            (received[0].pressure_pa * received[0].particle_velocity_m_s[axis].conj()).scale(0.5),
            2e-11,
        );
    }
    let base = owner.frequency_sample(omega, || false).unwrap();
    let dense = observation(
        &surface,
        Medium {
            density: 3.6,
            ..MEDIUM
        },
    )
    .frequency_sample(omega, || false)
    .unwrap();
    for i in 0..4 {
        close(dense[i], base[i].scale(3.0), 2e-11);
        close(dense[i + 4], base[i + 4], 2e-11);
    }
    let direct = fields[2]
        .pressure_velocity_at(k, MEDIUM, points[0])
        .unwrap();
    assert!((received[0].pressure_pa - direct.pressure_pa).abs() > 1e-3 * direct.pressure_pa.abs());

    let gate = CancelGate::new_clock_free();
    let solver = IterativeOptions {
        tol: 1e-10,
        restart: 40,
        max_iters: 200,
        leaf_capacity: 64,
        fmm_order: 4,
        ..IterativeOptions::default()
    };
    let (iterative, report, _) = incident::solve_scattering_iterative(
        &surface,
        k,
        MEDIUM,
        Admittance::Rigid,
        fields[2],
        1e-3,
        Formulation::BurtonMiller,
        &solver,
        &gate,
    )
    .unwrap();
    assert!(report.converged);
    let actual = iterative
        .total_pressure_velocity_at(&points, 1e-3, Options::default())
        .unwrap();
    close(actual.scalar.pressure[0], received[0].pressure_pa, 2e-7);
    let iterative_owner = IterativeBoundaryObservation::for_mixed_scattering(
        &surface,
        MEDIUM,
        &specs(),
        &receivers(),
        rigid(&surface),
        IterativeObservationOptions {
            highest_omega: TAU * 1600.0,
            minimum_clearance_m: 1e-3,
            receiver_quadrature: Options::default(),
            solver,
            lower_window: None,
        },
        &gate,
    )
    .unwrap();
    let sample = iterative_owner
        .physical_frequency_sample(omega, &gate)
        .unwrap();
    for (&actual, &expected) in sample.transfer.iter().zip(&base) {
        close(actual, expected, 2e-7);
    }
    // The retained solve must preserve the same mixed units and direct field
    // through both angular output paths, without another boundary solve.
    let retained = iterative_owner
        .solve_frequency(omega, 4 * surface.areas().len(), &gate)
        .unwrap();
    assert!(matches!(
        retained.source_incident(0).unwrap(),
        Some(IncidentSource::Monopole(_))
    ));
    let Some(IncidentSource::ForceDipole(unit_force)) = retained.source_incident(1).unwrap() else {
        panic!("retained force coordinate lost its physical kind");
    };
    assert_eq!(unit_force.force_on_fluid_n, DIRECTION.map(C64::from_re));
    assert!(retained.source_monopole(1).is_err());
    assert_eq!(
        retained.observe_receivers(&gate).unwrap().transfer,
        sample.transfer
    );
    let directions = [[0.0, 0.0, -1.0], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let work = 2 * directions.len() * (16 * surface.areas().len() + 1);
    let outgoing = retained.far_field(&directions, work, &gate).unwrap();
    assert_eq!(outgoing.kernel_evaluations, work);
    assert_eq!(outgoing.sources.as_slice(), owner.sources());
    assert!(retained.far_field(&directions, work - 1, &gate).is_err());
    for source in 0..2 {
        for (index, direction) in directions.into_iter().enumerate() {
            close(
                outgoing.amplitudes[source * directions.len() + index] * amplitudes[source],
                solved[source].total_far_field(direction).unwrap(),
                2e-7,
            );
        }
    }
    let projected = retained
        .coherent_directivity(
            &amplitudes,
            IterativeDirectivityOptions {
                l_max: 8,
                angular: DirectivityProjectionOptions {
                    maximum_directions: 20_000,
                    relative_tolerance: 1e-8,
                },
                maximum_kernel_evaluations: 100_000_000,
                maximum_combined_entries: 2 * surface.areas().len(),
            },
            &gate,
        )
        .unwrap();
    for direction in directions {
        close(
            projected.table.evaluate(direction),
            solved[2].total_far_field(direction).unwrap(),
            2e-5,
        );
    }
    assert!(
        iterative_owner
            .physical_frequency_sample(0.0, &gate)
            .is_err()
    );
}

#[test]
fn g0_force_geometry_directions_dc_and_surface_trace_substitution_refuse() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let owner = observation(&surface, MEDIUM);
    assert!(owner.frequency_sample(0.0, || false).is_err());
    let d = declaration(
        owner.operator_hash(),
        owner.sources(),
        &receivers(),
        owner.geometric_delays_s(),
    );
    assert!(
        owner.bake(&d, || false).is_err(),
        "the FIR's implicit zero bin must not hide force DC"
    );
    let mut called = false;
    assert!(
        bake_observation_transfers(
            |_| {
                called = true;
                Ok(vec![C64::ZERO; 8])
            },
            &d,
            || false
        )
        .is_err()
    );
    assert!(
        !called,
        "a generic producer must not bypass force-DC admission"
    );
    assert!(
        BemObservation::new(
            &surface,
            MEDIUM,
            owner.sources(),
            &vec![vec![1.0; surface.areas().len()]; 2],
            &receivers(),
            1e-3
        )
        .is_err()
    );
    for position in [[0.0; 3], surface.centroids()[0], receivers()[0].position_m] {
        let invalid = [CompactSourceSpec::ForceDipole(ForceDipoleSpec {
            id: "force".into(),
            position_m: position,
            direction: DIRECTION,
        })];
        assert!(
            ScatteringObservation::new_mixed(
                &surface,
                MEDIUM,
                &invalid,
                &receivers(),
                TAU * 1600.0,
                1e-3,
                Options::default(),
                || false
            )
            .is_err()
        );
    }
    let force = [ForceDipole {
        position_m: [0.0; 3],
        force_on_fluid_n: [C64::ONE; 3],
    }];
    assert!(
        IncidentTraces::from_field(
            &surface,
            1.0,
            MEDIUM,
            IncidentField::ForceDipoles(&force),
            1e-3
        )
        .is_err()
    );
    let mut changed = specs();
    let CompactSourceSpec::ForceDipole(source) = &mut changed[1] else {
        unreachable!()
    };
    source.direction = [2.0, 0.0, 0.0];
    assert!(
        ScatteringObservation::new_mixed(
            &surface,
            MEDIUM,
            &changed,
            &receivers(),
            TAU * 1600.0,
            1e-3,
            Options::default(),
            || false
        )
        .is_err()
    );
    let CompactSourceSpec::ForceDipole(source) = &mut changed[1] else {
        unreachable!()
    };
    source.direction = [-0.6, 0.0, -0.8];
    let reverse = ScatteringObservation::new_mixed(
        &surface,
        MEDIUM,
        &changed,
        &receivers(),
        TAU * 1600.0,
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap();
    assert_ne!(reverse.operator_hash(), owner.operator_hash());
    let mut wrong_kind = d.clone();
    wrong_kind.sources[1] = SourceId::PointVolumeFlow("force-on-fluid".into());
    assert_ne!(wrong_kind.declaration_hash(), d.declaration_hash());
    let unwindowed =
        LocalBoundaryObservation::for_scattering(owner, rigid(&surface), None).unwrap();
    let unwindowed_d = declaration(
        unwindowed.operator_hash(),
        unwindowed.sources(),
        unwindowed.receivers(),
        unwindowed.geometric_delays_s(),
    );
    assert!(unwindowed.physical_frequency_sample(0.0, || false).is_err());
    assert!(unwindowed.bake(&unwindowed_d, || false).is_err());
    assert!(
        ScatteringObservation::new_mixed(
            &surface,
            MEDIUM,
            &specs(),
            &receivers(),
            TAU * 1600.0,
            1e-3,
            Options::default(),
            || true
        )
        .is_err()
    );
}

#[test]
fn g1_force_and_volume_histories_reach_the_existing_windowed_fir_in_si_units() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let window = LowerBandWindow::try_new(150.0, 400.0).unwrap();
    let owner = LocalBoundaryObservation::for_scattering(
        observation(&surface, MEDIUM),
        rigid(&surface),
        Some(window),
    )
    .unwrap();
    assert!(owner.physical_frequency_sample(0.0, || false).is_err());
    assert_eq!(
        owner.band_limited_frequency_sample(0.0, || false).unwrap(),
        vec![C64::ZERO; 8]
    );
    let d = declaration(
        owner.operator_hash(),
        owner.sources(),
        owner.receivers(),
        owner.geometric_delays_s(),
    );
    let baked = owner.bake(&d, || false).unwrap();
    assert_eq!(baked.lower_band_window(), Some(window));
    assert_eq!(baked.transfer().sources, d.sources);
    assert_eq!(baked.transfer().ladder.len(), 4);
    let omega = TAU * 500.0;
    let physical = owner.physical_frequency_sample(omega, || false).unwrap();
    let angle = -omega * d.filter_latency_samples as f64 / d.fs_hz;
    let latency = C64::new(det::cos(angle), det::sin(angle));
    for (kernel, value) in baked.transfer().kernels.iter().zip(physical) {
        let expected = value.conj() * latency;
        close(kernel.transfer_at(omega / d.fs_hz).unwrap(), expected, 2e-4);
    }
    let mut playback = baked.observer().unwrap();
    let mut scaled = baked.observer().unwrap();
    let mut peak = 0.0_f64;
    for time in 0..512 {
        let input = if time == 0 { [1e-6, 1e-3] } else { [0.0; 2] };
        let output = playback.step(&input).unwrap();
        let doubled = scaled.step(&input.map(|v| -2.0 * v)).unwrap();
        for (a, b) in output.iter().zip(doubled) {
            assert!((b + 2.0 * a).abs() <= 1e-12 * a.abs().max(1e-18));
        }
        peak = peak.max(output[0].abs());
    }
    assert!(peak > 1e-9, "physical source histories were discarded");
}
