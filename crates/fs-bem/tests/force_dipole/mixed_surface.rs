//! G1/G3: finite moving-wall and aperture traces interfere with actual compact
//! flow/force fields through one original passive boundary operator.
use super::*;
use fs_bem::helmholtz::{self, bc::PreparedBoundary};
use fs_bem::near_field::Geometry;

fn finite_sources(surface: &SpherePanels) -> ([SourceId; 2], Vec<Vec<f64>>) {
    let mut wall = vec![0.0; surface.areas().len()];
    wall[0] = 1.0;
    wall[1] = -surface.areas()[0] / surface.areas()[1];
    let mut aperture = vec![0.0; wall.len()];
    aperture[2] = 1.0 / surface.areas()[2];
    (
        [
            SourceId::StructuralMode("wall".into()),
            SourceId::ApertureVolumeFlow("aperture".into()),
        ],
        vec![wall, aperture],
    )
}

fn boundary(surface: &SpherePanels) -> LocalBoundary {
    LocalBoundary::try_new(
        MEDIUM,
        vec![
            LocalWall::Rigid,
            LocalWall::Mechanical(fs_phs::WallPin {
                surface_density: 0.1,
                stiffness_per_area: 4e5,
                resistance: 200.0,
            }),
        ],
        (0..surface.areas().len())
            .map(|i| usize::from(i >= 3))
            .collect(),
    )
    .unwrap()
}

fn dense_pair(
    surface: &SpherePanels,
    window: Option<LowerBandWindow>,
) -> (LocalBoundaryObservation<'_>, LocalBoundaryObservation<'_>) {
    let (ids, traces) = finite_sources(surface);
    let finite = BemObservation::new(surface, MEDIUM, &ids, &traces, &receivers(), 1e-3).unwrap();
    (
        LocalBoundaryObservation::for_surface(
            finite,
            boundary(surface),
            TAU * 1600.0,
            Options::default(),
            window,
            || false,
        )
        .unwrap(),
        LocalBoundaryObservation::for_scattering(
            observation(surface, MEDIUM),
            boundary(surface),
            window,
        )
        .unwrap(),
    )
}

fn iterative_pair<'a>(
    surface: &'a SpherePanels,
    gate: &CancelGate,
) -> (
    IterativeBoundaryObservation<'a>,
    IterativeBoundaryObservation<'a>,
) {
    let options = IterativeObservationOptions {
        highest_omega: TAU * 1600.0,
        minimum_clearance_m: 1e-3,
        receiver_quadrature: Options::default(),
        lower_window: None,
        solver: IterativeOptions {
            tol: 1e-10,
            restart: 40,
            max_iters: 200,
            leaf_capacity: 64,
            fmm_order: 4,
            ..IterativeOptions::default()
        },
    };
    let (ids, traces) = finite_sources(surface);
    (
        IterativeBoundaryObservation::for_surface(
            surface,
            MEDIUM,
            &ids,
            &traces,
            &receivers(),
            boundary(surface),
            options.clone(),
            gate,
        )
        .unwrap(),
        IterativeBoundaryObservation::for_mixed_scattering(
            surface,
            MEDIUM,
            &specs(),
            &receivers(),
            boundary(surface),
            options,
            gate,
        )
        .unwrap(),
    )
}

// Independently assemble ONE Robin right-hand side from prescribed wall
// velocity and total analytic incident p/v: v_out + Y p_out = v_drive - v_inc - Y p_inc.
fn combined_solution(
    surface: &SpherePanels,
    omega: f64,
    drive: &[C64; 4],
) -> helmholtz::RadiationSolution {
    let k = omega / MEDIUM.sound_speed;
    let incident = IncidentTraces::from_sources(
        surface,
        k,
        MEDIUM,
        &illumination([drive[2], drive[3]]),
        1e-3,
    )
    .unwrap();
    let beta = boundary(surface)
        .normalized_admittance(omega, || false)
        .unwrap();
    let (_, traces) = finite_sources(surface);
    let velocity: Vec<_> = (0..surface.areas().len())
        .map(|i| {
            drive[0].scale(traces[0][i]) + drive[1].scale(traces[1][i])
                - incident.normal_velocity_m_s()[i]
                - (beta[i] * incident.pressure_pa()[i])
                    .scale(1.0 / (MEDIUM.density * MEDIUM.sound_speed))
        })
        .collect();
    PreparedBoundary::new(
        surface,
        k,
        MEDIUM,
        Admittance::Local {
            beta: &beta,
            declared_active: false,
        },
        Formulation::BurtonMiller,
    )
    .unwrap()
    .solve(&velocity)
    .unwrap()
    .radiation
}

fn analytic_far(omega: f64, drive: &[C64; 4], direction: [f64; 3]) -> C64 {
    let k = omega / MEDIUM.sound_speed;
    let phase = |position: [f64; 3]| {
        let angle = -k * (0..3).map(|i| direction[i] * position[i]).sum::<f64>();
        C64::new(det::cos(angle), det::sin(angle))
    };
    let projection = (0..3).map(|i| direction[i] * DIRECTION[i]).sum::<f64>();
    let denominator = 4.0 * std::f64::consts::PI;
    C64::new(0.0, -omega * MEDIUM.density / denominator) * drive[2] * phase(Q_POSITION)
        + C64::new(0.0, -k * projection / denominator) * drive[3] * phase(F_POSITION)
}

#[test]
fn dense_complete_bank_matches_independent_combined_boundary_and_signed_coordinate_change() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let omega = TAU * 900.0;
    let (finite, compact) = dense_pair(&surface, None);
    let mut separate = finite.physical_frequency_sample(omega, || false).unwrap();
    separate.extend(compact.physical_frequency_sample(omega, || false).unwrap());
    let mut delays = finite.geometric_delays_s();
    delays.extend(compact.geometric_delays_s());
    let combined =
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .unwrap();
    assert_eq!(
        combined.sources(),
        [
            SourceId::StructuralMode("wall".into()),
            SourceId::ApertureVolumeFlow("aperture".into()),
            SourceId::PointVolumeFlow("flow".into()),
            SourceId::PointForce("force-on-fluid".into())
        ]
    );
    assert_eq!(combined.geometric_delays_s(), delays);
    let samples = combined.physical_frequency_sample(omega, || false).unwrap();
    for (&a, &b) in samples.iter().zip(&separate) {
        close(a, b, 2e-11);
    }
    let bank = combined.prepare_harmonic(omega, || false).unwrap();
    let drive = [
        C64::new(0.002, -0.001),
        C64::new(1e-6, 0.3e-6),
        C64::new(-0.4e-6, 0.8e-6),
        C64::new(0.003, -0.002),
    ];
    let mut observed = [HarmonicField::default()];
    bank.evaluate_coherent(&drive, &mut observed, || false)
        .unwrap();
    let outgoing = combined_solution(&surface, omega, &drive);
    let points = receivers().map(|r| r.position_m);
    let geometry = Geometry::new(&surface, &points, 1e-3).unwrap();
    let field = geometry
        .prepare_velocity(omega / MEDIUM.sound_speed, MEDIUM, Options::default())
        .unwrap()
        .evaluate_velocity(&outgoing)
        .unwrap();
    let illumination = illumination([drive[2], drive[3]]);
    let direct = IncidentField::Sources(&illumination)
        .pressure_velocity_at(omega / MEDIUM.sound_speed, MEDIUM, points[0])
        .unwrap();
    close(
        observed[0].pressure_pa,
        field.scalar.pressure[0] + direct.pressure_pa,
        2e-10,
    );
    for axis in 0..3 {
        close(
            observed[0].particle_velocity_m_s[axis],
            field.particle_velocity_m_s[0][axis] + direct.velocity_m_s[axis],
            2e-10,
        );
    }

    let (ids, mut traces) = finite_sources(&surface);
    for value in &mut traces[0] {
        *value = -*value;
    }
    let finite = LocalBoundaryObservation::for_surface(
        BemObservation::new(&surface, MEDIUM, &ids, &traces, &receivers(), 1e-3).unwrap(),
        boundary(&surface),
        TAU * 1600.0,
        Options::default(),
        None,
        || false,
    )
    .unwrap();
    let (_, compact) = dense_pair(&surface, None);
    let reversed =
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .unwrap();
    assert_ne!(reversed.operator_hash(), combined.operator_hash());
    let mut reversed_drive = drive;
    reversed_drive[0] = -reversed_drive[0];
    let mut result = [HarmonicField::default()];
    reversed
        .prepare_harmonic(omega, || false)
        .unwrap()
        .evaluate_coherent(&reversed_drive, &mut result, || false)
        .unwrap();
    close(result[0].pressure_pa, observed[0].pressure_pa, 2e-11);
    for axis in 0..3 {
        close(
            result[0].particle_velocity_m_s[axis],
            observed[0].particle_velocity_m_s[axis],
            2e-11,
        );
    }
}

#[test]
fn retained_mixed_bank_preserves_direct_terms_nulls_and_partial_coherence_to_far_field() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let omega = TAU * 900.0;
    let gate = CancelGate::new_clock_free();
    let (finite, compact) = iterative_pair(&surface, &gate);
    let owner =
        IterativeBoundaryObservation::combine_surface_and_scattering(finite, compact, &gate)
            .unwrap();
    let retained = owner
        .solve_frequency(omega, 8 * surface.areas().len(), &gate)
        .unwrap();
    assert!(
        retained
            .reports()
            .iter()
            .all(|report| report.solve.converged)
    );
    assert!(retained.source_incident(0).unwrap().is_none());
    assert!(retained.source_incident(1).unwrap().is_none());
    assert!(matches!(
        retained.source_incident(2).unwrap(),
        Some(IncidentSource::Monopole(_))
    ));
    assert!(matches!(
        retained.source_incident(3).unwrap(),
        Some(IncidentSource::ForceDipole(_))
    ));
    assert!(
        retained
            .velocity_port_impedance(16, 16 * surface.areas().len() as u64, &gate)
            .is_err()
    );
    assert!(owner.aperture_support_areas_m2(&gate).is_err());
    let (finite, compact) = dense_pair(&surface, None);
    let dense = LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
        .unwrap();
    let reference = dense.physical_frequency_sample(omega, || false).unwrap();
    let sample = retained.observe_receivers(&gate).unwrap();
    for (&a, &b) in sample.transfer.iter().zip(&reference) {
        close(a, b, 2e-7);
    }

    let bank = retained.prepare_harmonic(&gate).unwrap();
    let mut cancelled = [C64::ZERO; 4];
    cancelled[0] = C64::from_re(1e-3);
    cancelled[3] = -sample.transfer[0] * cancelled[0] / sample.transfer[12];
    let mut bright = cancelled;
    bright[3] = -bright[3];
    let mut dark_field = [HarmonicField::default()];
    let mut bright_field = [HarmonicField::default()];
    bank.evaluate_coherent(&cancelled, &mut dark_field, || false)
        .unwrap();
    bank.evaluate_coherent(&bright, &mut bright_field, || false)
        .unwrap();
    assert!(dark_field[0].pressure_pa.abs() < 1e-12 * bright_field[0].pressure_pa.abs());
    assert!(
        dark_field[0]
            .particle_velocity_m_s
            .iter()
            .any(|v| v.abs() > 1e-12)
    );
    let factors: Vec<_> = cancelled.into_iter().chain(bright).collect();
    let mut statistics = [HarmonicFieldStatistics::default()];
    bank.evaluate_coherence_modes(
        &factors,
        2,
        HarmonicWorkBudget {
            maximum_products: 10000,
        },
        &mut statistics,
        || false,
    )
    .unwrap();
    close(
        C64::from_re(statistics[0].mean_square_pressure_pa2),
        C64::from_re(
            0.5 * (dark_field[0].pressure_pa.norm_sq() + bright_field[0].pressure_pa.norm_sq()),
        ),
        2e-11,
    );

    let directions = [[1.0, 0.0, 0.0], [0.0, 0.0, -1.0], [0.0, 1.0, 0.0]];
    let angular = retained
        .prepare_radiation(
            &directions,
            PreparedRadiationOptions {
                angular: DirectivityProjectionOptions::default(),
                initial_latitudes: 8,
                refinements: 3,
                maximum_kernel_evaluations: 4_000_000,
                maximum_source_entries: 6000,
            },
            &gate,
        )
        .unwrap();
    let budget = HarmonicWorkBudget {
        maximum_products: 1_000_000,
    };
    let mut a = [C64::ZERO; 3];
    let mut b = [C64::ZERO; 3];
    let pa = angular
        .evaluate_coherent(&cancelled, budget, &mut a, || false)
        .unwrap();
    let pb = angular
        .evaluate_coherent(&bright, budget, &mut b, || false)
        .unwrap();
    let mut density = [0.0; 3];
    let partial = angular
        .evaluate_coherence_modes(&factors, 2, budget, &mut density, || false)
        .unwrap();
    assert!(pa.total_power_w > 0.0); // A microphone null is not silence everywhere.
    close(
        C64::from_re(partial.total_power_w),
        C64::from_re(pa.total_power_w + pb.total_power_w),
        2e-10,
    );
    let solution = combined_solution(&surface, omega, &cancelled);
    let far = helmholtz::far_field_integrated(&surface, &solution, MEDIUM, &directions);
    for i in 0..3 {
        close(
            a[i],
            far[i] + analytic_far(omega, &cancelled, directions[i]),
            2e-7,
        );
        close(
            C64::from_re(density[i]),
            C64::from_re(
                (a[i].norm_sq() + b[i].norm_sq()) / (2.0 * MEDIUM.density * MEDIUM.sound_speed),
            ),
            2e-10,
        );
    }
}

#[test]
fn combined_admission_preserves_identity_force_window_clearance_and_cancellation() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let other_surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let (finite, _) = dense_pair(&surface, None);
    let (_, compact) = dense_pair(&other_surface, None);
    assert!(
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .is_err()
    );
    let window = Some(LowerBandWindow::try_new(100.0, 200.0).unwrap());
    let (finite, _) = dense_pair(&surface, window);
    let (_, compact) = dense_pair(&surface, None);
    assert!(
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .is_err()
    );
    let (finite, compact) = dense_pair(&surface, None);
    assert!(
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || true).is_err()
    );
    let (finite, _) = dense_pair(&surface, None);
    let compact = LocalBoundaryObservation::for_scattering(
        observation(&surface, MEDIUM),
        rigid(&surface),
        None,
    )
    .unwrap();
    assert!(
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .is_err()
    );
    let (finite, _) = dense_pair(&surface, None);
    let mut duplicate = specs();
    let CompactSourceSpec::VolumeFlow(source) = &mut duplicate[0] else {
        unreachable!()
    };
    source.id = "wall".into();
    let compact = LocalBoundaryObservation::for_scattering(
        ScatteringObservation::new_mixed(
            &surface,
            MEDIUM,
            &duplicate,
            &receivers(),
            TAU * 1600.0,
            1e-3,
            Options::default(),
            || false,
        )
        .unwrap(),
        boundary(&surface),
        None,
    )
    .unwrap();
    assert!(
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .is_err()
    );
    let CompactSourceSpec::VolumeFlow(source) = &mut duplicate[0] else {
        unreachable!()
    };
    source.position_m = [0.0; 3];
    assert!(
        ScatteringObservation::new_mixed(
            &surface,
            MEDIUM,
            &duplicate,
            &receivers(),
            TAU * 1600.0,
            1e-3,
            Options::default(),
            || false
        )
        .is_err()
    );
    let (finite, compact) = dense_pair(&surface, None);
    let combined =
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .unwrap();
    assert!(combined.physical_frequency_sample(0.0, || false).is_err());
    let decl = declaration(
        combined.operator_hash(),
        combined.sources(),
        combined.receivers(),
        combined.geometric_delays_s(),
    );
    assert!(combined.bake(&decl, || false).is_err());
    let (finite, compact) = dense_pair(&surface, window);
    let windowed =
        LocalBoundaryObservation::combine_surface_and_scattering(finite, compact, || false)
            .unwrap();
    assert_eq!(windowed.lower_band_window(), window);
    assert!(
        windowed
            .band_limited_frequency_sample(0.0, || false)
            .unwrap()
            .iter()
            .all(|&v| v == C64::ZERO)
    );
    let gate = CancelGate::new_clock_free();
    let (finite, compact) = iterative_pair(&surface, &gate);
    gate.request();
    assert!(
        IterativeBoundaryObservation::combine_surface_and_scattering(finite, compact, &gate)
            .is_err()
    );
}
