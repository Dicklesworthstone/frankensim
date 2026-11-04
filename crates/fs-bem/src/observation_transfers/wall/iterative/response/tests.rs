use super::*;
use core::f64::consts::TAU;
use fs_phs::WallPin;

fn medium() -> Medium {
    Medium {
        density: 1.2,
        sound_speed: 340.0,
    }
}

fn receivers() -> [ReceiverSpec; 2] {
    [
        ReceiverSpec {
            id: "front".into(),
            position_m: [0.12, 0.04, 0.06],
        },
        ReceiverSpec {
            id: "back".into(),
            position_m: [-0.04, 0.02, -0.09],
        },
    ]
}

fn options(leaf_capacity: usize) -> IterativeObservationOptions {
    IterativeObservationOptions {
        highest_omega: TAU * 1600.0,
        minimum_clearance_m: 1e-3,
        receiver_quadrature: Options::default(),
        lower_window: None,
        solver: IterativeOptions {
            tol: 1e-10,
            restart: 80,
            max_iters: 600,
            fmm_order: 7,
            leaf_capacity,
            ..IterativeOptions::default()
        },
    }
}

// Two independently driven finite hemispherical patches on one body, plus
// a separate passive material body. The original triangle ordering determines
// the physical source map; no point-source or isolated-radiator substitution.
fn fixture(
    level: u32,
    obstacle: Option<[f64; 3]>,
) -> (SpherePanels, Vec<SourceId>, Vec<Vec<f64>>, LocalBoundary) {
    let radiator = SpherePanels::icosphere(0.025, level).unwrap();
    let driven_panels = radiator.areas().len();
    let mut triangles = radiator.triangles().unwrap().to_vec();
    if let Some(offset) = obstacle {
        let body = SpherePanels::icosphere(0.02, level).unwrap();
        triangles.extend(body.triangles().unwrap().iter().map(|triangle| {
            triangle.map(|point| std::array::from_fn(|axis| point[axis] + offset[axis]))
        }));
    }
    let surface = SpherePanels::from_triangles(triangles).unwrap();
    let mut traces = vec![vec![0.0; surface.areas().len()]; 2];
    for (row, trace) in traces.iter_mut().enumerate() {
        let indices: Vec<_> = (0..driven_panels)
            .filter(|&panel| usize::from(surface.normals()[panel][2] < 0.0) == row)
            .collect();
        let area: f64 = indices.iter().map(|&panel| surface.areas()[panel]).sum();
        for panel in indices {
            trace[panel] = 1.0 / area;
        }
    }
    let mut assignments = vec![0; surface.areas().len()];
    assignments[driven_panels..].fill(1);
    let boundary = LocalBoundary::try_new_iterative(
        medium(),
        vec![
            LocalWall::Rigid,
            LocalWall::Mechanical(WallPin {
                surface_density: 0.1,
                stiffness_per_area: 4e5,
                resistance: 200.0,
            }),
        ],
        assignments,
    )
    .unwrap();
    (
        surface,
        vec![
            SourceId::ApertureVolumeFlow("upper".into()),
            SourceId::ApertureVolumeFlow("lower".into()),
        ],
        traces,
        boundary,
    )
}

fn relative(actual: &[C64], expected: &[C64]) -> f64 {
    let delta: f64 = actual
        .iter()
        .zip(expected)
        .map(|(&a, &b)| (a - b).norm_sq())
        .sum();
    let norm: f64 = expected.iter().map(|z| z.norm_sq()).sum();
    det::sqrt(delta / norm.max(f64::MIN_POSITIVE))
}

#[test]
fn retained_multibody_fields_supply_loading_receivers_and_coherent_wall_work() {
    let (surface, sources, traces, boundary) = fixture(0, Some([0.06, 0.0, 0.02]));
    let gate = CancelGate::new_clock_free();
    let omega = TAU * 700.0;
    let beta = boundary.normalized_admittance(omega, || false).unwrap();
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        medium(),
        &sources,
        &traces,
        &receivers(),
        boundary,
        options(128),
        &gate,
    )
    .unwrap();
    let retained = owner
        .solve_frequency(omega, 4 * surface.areas().len(), &gate)
        .unwrap();
    let ports = retained
        .velocity_port_impedance(4, 4 * surface.areas().len() as u64, &gate)
        .unwrap();
    let observed = retained.observe_receivers(&gate).unwrap();
    let streamed = owner.physical_frequency_sample(omega, &gate).unwrap();
    assert_eq!(observed.transfer, streamed.transfer);
    let harmonic = retained.prepare_harmonic(&gate).unwrap();
    let mut coherent = vec![HarmonicField::default(); receivers().len()];
    harmonic
        .evaluate_coherent(&[C64::ONE, C64::ZERO], &mut coherent, || false)
        .unwrap();
    let euler = C64::new(0.0, -1.0 / (omega * medium().density));
    for (receiver, field) in coherent.iter().enumerate() {
        assert_eq!(field.pressure_pa, observed.transfer[4 * receiver]);
        for axis in 0..3 {
            assert_eq!(
                field.particle_velocity_m_s[axis],
                euler * observed.transfer[4 * receiver + axis + 1]
            );
        }
    }
    assert_eq!(ports.sources(), owner.sources());
    assert_eq!(ports.operator_hash(), owner.operator_hash());
    assert_eq!(ports.acoustic_impedance().unwrap(), ports.impedance());
    assert_eq!(
        ports.projection_panel_products(),
        4 * surface.areas().len() as u64
    );
    assert!(ports.impedance()[1].abs() > 1e-3 * ports.impedance()[0].abs());
    assert!(ports.impedance()[2].abs() > 1e-3 * ports.impedance()[3].abs());
    assert!(ports.reciprocity_defect().is_finite());
    assert!(ports.passivity_margin().unwrap().is_finite());
    assert!(
        retained
            .reports()
            .iter()
            .all(|report| report.solve.converged_euclidean() && report.absorbed_power_w > 0.0)
    );

    // Independent physical surface work, including coherent cross terms.
    // This is the work identity of the returned discrete fields; it does not
    // promote a collocation matrix to a certified passive rational model.
    let rates = [C64::new(0.8e-6, 0.3e-6), C64::new(-0.2e-6, 0.7e-6)];
    let mut direct_drive_work = 0.0;
    let mut exterior_flux = 0.0;
    let mut absorbed = 0.0;
    let rho_c = medium().density * medium().sound_speed;
    for (panel, &area) in surface.areas().iter().enumerate() {
        let mut pressure = C64::ZERO;
        let mut velocity = C64::ZERO;
        let mut drive = C64::ZERO;
        for (source, &rate) in rates.iter().enumerate() {
            pressure = pressure + retained.solutions()[source].pressure[panel] * rate;
            velocity = velocity + retained.solutions()[source].velocity[panel] * rate;
            drive = drive + rate.scale(traces[source][panel]);
        }
        direct_drive_work += 0.5 * (pressure * drive.conj()).re * area;
        exterior_flux += 0.5 * (pressure * velocity.conj()).re * area;
        absorbed += 0.5 * beta[panel].re * pressure.norm_sq() * area / rho_c;
    }
    let input = ports.input_power_w(&rates).unwrap();
    let scale = input.abs().max(absorbed).max(f64::MIN_POSITIVE);
    assert!(absorbed > 0.0);
    assert!((input - direct_drive_work).abs() < 2e-12 * scale);
    assert!((input - exterior_flux - absorbed).abs() < 2e-12 * scale);
    let incoherent: f64 = retained
        .reports()
        .iter()
        .zip(rates)
        .map(|(report, rate)| report.absorbed_power_w * rate.norm_sq())
        .sum();
    assert!((absorbed - incoherent).abs() > 1e-6 * absorbed);

    let work = 2 * rates.len() as u64 * surface.areas().len() as u64;
    let coherent = retained
        .coherent_boundary_power(&rates, work, &gate)
        .unwrap();
    assert_eq!(coherent.omega_rad_s, omega);
    assert_eq!(coherent.operator_hash, retained.operator_hash());
    assert_eq!(coherent.source_panel_products, work);
    assert!((coherent.input_power_w - direct_drive_work).abs() < 2e-12 * scale);
    assert!((coherent.outward_power_w - exterior_flux).abs() < 2e-12 * scale);
    assert!((coherent.absorbed_power_w - absorbed).abs() < 2e-12 * scale);
    assert!(coherent.relative_balance < 2e-12);
    assert!(coherent.apparent_power_scale_w >= coherent.input_power_w.abs());

    // The same phases and physical fields, over 580 decades of power.
    // Neither squaring the raw rates nor normalizing all fields to one SI
    // scale would preserve this range without avoidable overflow/underflow.
    for amplitude in [1e-145, 1.0, 1e145] {
        let transformed = rates.map(|rate| C64::new(-rate.im, rate.re).scale(amplitude));
        let changed = retained
            .coherent_boundary_power(&transformed, work, &gate)
            .unwrap();
        for (actual, reference) in [
            (changed.input_power_w, coherent.input_power_w),
            (changed.outward_power_w, coherent.outward_power_w),
            (changed.absorbed_power_w, coherent.absorbed_power_w),
            (
                changed.apparent_power_scale_w,
                coherent.apparent_power_scale_w,
            ),
        ] {
            let expected = (reference * amplitude) * amplitude;
            assert!(
                (actual / expected - 1.0).abs() < 2e-12,
                "{actual:e} != {expected:e}"
            );
        }
        assert!(changed.relative_balance < 2e-12);
    }
    let zero = retained
        .coherent_boundary_power(&[C64::ZERO; 2], work, &gate)
        .unwrap();
    assert_eq!(zero.input_power_w, 0.0);
    assert_eq!(zero.outward_power_w, 0.0);
    assert_eq!(zero.absorbed_power_w, 0.0);
    assert_eq!(zero.apparent_power_scale_w, 0.0);
    assert_eq!(zero.relative_balance, 0.0);

    // Removing the passive body changes the same finite radiator's physical
    // load and microphones, not just a later directional display.
    let (bare_surface, bare_sources, bare_traces, bare_boundary) = fixture(0, None);
    let bare = IterativeBoundaryObservation::for_surface(
        &bare_surface,
        medium(),
        &bare_sources,
        &bare_traces,
        &receivers(),
        bare_boundary,
        options(128),
        &gate,
    )
    .unwrap();
    let bare_response = bare
        .solve_frequency(omega, 4 * bare_surface.areas().len(), &gate)
        .unwrap();
    let bare_ports = bare_response
        .velocity_port_impedance(4, u64::MAX, &gate)
        .unwrap();
    assert!(relative(ports.impedance(), bare_ports.impedance()) > 1e-4);
    let bare_observed = bare_response.observe_receivers(&gate).unwrap();
    // Compare only pressure [Pa], avoiding a mixed p/gradient norm.
    let pressure: Vec<_> = observed.transfer.iter().step_by(4).copied().collect();
    let bare_pressure: Vec<_> = bare_observed.transfer.iter().step_by(4).copied().collect();
    assert!(relative(&pressure, &bare_pressure) > 1e-4);
}

#[test]
fn coherent_boundary_work_retains_quiet_fields_between_large_cancelling_sources() {
    let (surface, _, traces, boundary) = fixture(0, Some([0.06, 0.0, 0.02]));
    let gate = CancelGate::new_clock_free();
    // The last mode reverses the original signed shape. Its positive rate
    // cancels the first mode, while the middle mode carries the quiet field.
    let sources: Vec<_> = (0..3)
        .map(|i| SourceId::StructuralMode(format!("shape-{i}")))
        .collect();
    let modes = [
        traces[0].clone(),
        traces[0].clone(),
        traces[0].iter().map(|v| -v).collect(),
    ];
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        medium(),
        &sources,
        &modes,
        &receivers(),
        boundary,
        options(128),
        &gate,
    )
    .unwrap();
    let mut response = owner
        .solve_frequency(TAU * 700.0, 6 * surface.areas().len(), &gate)
        .unwrap();
    // Preserve exact +/- columns for this arithmetic test. Independent GMRES
    // roundoff is not source cancellation accuracy and is tested elsewhere.
    let positive = response.solutions[0].clone();
    response.solutions[1] = positive.clone();
    response.solutions[2].pressure = positive.pressure.iter().map(|&v| -v).collect();
    response.solutions[2].velocity = positive.velocity.iter().map(|&v| -v).collect();
    let work = 6 * surface.areas().len() as u64;
    let reference = response
        .coherent_boundary_power(&[C64::ZERO, C64::ONE, C64::ZERO], work, &gate)
        .unwrap();
    let quiet = response
        .coherent_boundary_power(
            &[C64::from_re(1e300), C64::ONE, C64::from_re(1e300)],
            work,
            &gate,
        )
        .unwrap();
    for (actual, expected) in [
        (quiet.input_power_w, reference.input_power_w),
        (quiet.outward_power_w, reference.outward_power_w),
        (quiet.absorbed_power_w, reference.absorbed_power_w),
    ] {
        assert!(
            (actual / expected - 1.0).abs() < 2e-12,
            "{actual:e} != {expected:e}"
        );
    }
    assert!(quiet.relative_balance < 2e-12);
    // Finite input is insufficient when its nonzero physical work cannot be
    // represented. Do not report either infinity or an artificial zero loss.
    for amplitude in [1e300, 1e-300] {
        assert!(
            response
                .coherent_boundary_power(
                    &[C64::ZERO, C64::from_re(amplitude), C64::ZERO],
                    work,
                    &gate,
                )
                .is_err()
        );
    }
    // A defect in the retained actual velocity must remain visible in both
    // signed flux and closure; rebuilding velocity from beta would hide it.
    for velocity in &mut response.solutions[1].velocity {
        *velocity = -*velocity;
    }
    let perturbed = response
        .coherent_boundary_power(&[C64::ZERO, C64::ONE, C64::ZERO], work, &gate)
        .unwrap();
    assert!(perturbed.outward_power_w < 0.0);
    assert!((perturbed.outward_power_w / reference.outward_power_w + 1.0).abs() < 2e-12);
    assert!((perturbed.balance_w / reference.outward_power_w - 2.0).abs() < 2e-12);
    assert!(perturbed.relative_balance > 1e-6);
}

#[test]
fn three_fmm_orders_converge_to_dense_material_mutual_loading() {
    let (surface, sources, traces, boundary) = fixture(1, Some([0.06, 0.0, 0.02]));
    let gate = CancelGate::new_clock_free();
    let omega = TAU * 1100.0;
    let beta = boundary.normalized_admittance(omega, || false).unwrap();
    let dense = PreparedBoundary::new(
        &surface,
        omega / medium().sound_speed,
        medium(),
        Admittance::Local {
            beta: &beta,
            declared_active: false,
        },
        helmholtz::Formulation::BurtonMiller,
    )
    .unwrap();
    let mut expected = vec![C64::ZERO; 4];
    for (column, trace) in traces.iter().enumerate() {
        let velocity: Vec<_> = trace.iter().map(|&v| C64::from_re(v)).collect();
        let solution = dense.solve(&velocity).unwrap();
        // Independently integrate against the physical real source map.
        for (row, test_trace) in traces.iter().enumerate() {
            expected[2 * row + column] = solution
                .radiation
                .pressure
                .iter()
                .zip(test_trace)
                .zip(surface.areas())
                .fold(C64::ZERO, |sum, ((&p, &shape), &area)| {
                    sum + p.scale(shape * area)
                });
        }
    }
    let mut errors = Vec::new();
    for order in [4, 6, 7] {
        let mut config = options(4);
        config.solver.fmm_order = order;
        let owner = IterativeBoundaryObservation::for_surface(
            &surface,
            medium(),
            &sources,
            &traces,
            &receivers(),
            boundary.clone(),
            config,
            &gate,
        )
        .unwrap();
        let response = owner
            .solve_frequency(omega, 4 * surface.areas().len(), &gate)
            .unwrap();
        let ports = response
            .velocity_port_impedance(4, u64::MAX, &gate)
            .unwrap();
        errors.push(relative(ports.impedance(), &expected));
        assert!(response.reports().iter().all(
            |report| report.solve.converged_euclidean() && report.resources.fmm_order == order
        ));
    }
    assert!(
        errors[2] < errors[0],
        "FMM order must reduce actual port error: {errors:?}"
    );
    assert!(
        errors[2] < 2e-5,
        "finest material self/mutual matrix agreement: {errors:?}"
    );
}

#[test]
fn signed_structural_port_keeps_work_units_and_refuses_unfunded_consumers() {
    let (surface, mut sources, mut traces, boundary) = fixture(0, None);
    sources[0] = SourceId::StructuralMode("signed-mode".into());
    for panel in 0..surface.areas().len() {
        traces[0][panel] = surface.normals()[panel][2];
    }
    let gate = CancelGate::new_clock_free();
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        medium(),
        &sources,
        &traces,
        &receivers(),
        boundary,
        options(128),
        &gate,
    )
    .unwrap();
    let omega = TAU * 700.0;
    let retained_entries = 4 * surface.areas().len();
    assert!(matches!(
        owner.solve_frequency(omega, retained_entries - 1, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Cap { .. }
        ))
    ));
    assert!(owner.solve_frequency(0.0, retained_entries, &gate).is_err());
    let response = owner
        .solve_frequency(omega, retained_entries, &gate)
        .unwrap();
    assert!(response.source_monopole(0).unwrap().is_none());
    assert!(response.source_monopole(2).is_err());
    assert!(
        response
            .velocity_port_impedance(3, u64::MAX, &gate)
            .is_err()
    );
    assert!(
        response
            .velocity_port_impedance(4, retained_entries as u64 - 1, &gate)
            .is_err()
    );
    let ports = response
        .velocity_port_impedance(4, u64::MAX, &gate)
        .unwrap();
    assert!(ports.acoustic_impedance().is_err());
    assert!(ports.resisting_efforts(&[C64::ONE]).is_err());
    let rates = [C64::new(0.4, -0.2), C64::new(-2e-6, 3e-6)];
    assert!(response.coherent_boundary_power(&rates, 0, &gate).is_err());
    assert!(
        response
            .coherent_boundary_power(&rates[..1], u64::MAX, &gate)
            .is_err()
    );
    assert!(
        response
            .coherent_boundary_power(&[C64::ONE, C64::from_re(f64::NAN)], u64::MAX, &gate)
            .is_err()
    );
    let efforts = ports.resisting_efforts(&rates).unwrap();
    for row in 0..2 {
        let expected =
            ports.impedance()[2 * row] * rates[0] + ports.impedance()[2 * row + 1] * rates[1];
        assert_eq!(efforts[row], expected);
    }
    gate.request();
    assert!(matches!(
        response.observe_receivers(&gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    assert!(matches!(
        response.velocity_port_impedance(4, u64::MAX, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
    assert!(matches!(
        response.coherent_boundary_power(&rates, u64::MAX, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
}

#[test]
fn point_response_retains_total_observation_but_has_no_finite_self_port() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let gate = CancelGate::new_clock_free();
    let boundary = LocalBoundary::try_new_iterative(
        medium(),
        vec![LocalWall::Rigid],
        vec![0; surface.areas().len()],
    )
    .unwrap();
    let point = [PointSourceSpec {
        id: "point".into(),
        position_m: [0.02, 0.01, 0.09],
    }];
    let owner = IterativeBoundaryObservation::for_scattering(
        &surface,
        medium(),
        &point,
        &receivers(),
        boundary,
        options(128),
        &gate,
    )
    .unwrap();
    let omega = TAU * 700.0;
    let response = owner
        .solve_frequency(omega, 2 * surface.areas().len(), &gate)
        .unwrap();
    assert!(
        response
            .coherent_boundary_power(&[C64::ONE], u64::MAX, &gate)
            .is_err()
    );
    assert_eq!(
        response.source_monopole(0).unwrap().unwrap().position_m,
        point[0].position_m
    );
    assert_eq!(
        response.observe_receivers(&gate).unwrap().transfer,
        owner
            .physical_frequency_sample(omega, &gate)
            .unwrap()
            .transfer
    );
    assert!(
        response
            .velocity_port_impedance(1, u64::MAX, &gate)
            .is_err()
    );
    assert!(
        PreparedBoundaryWorkExposure::new(
            &response,
            "point-clock",
            &[],
            &[],
            exposure_limits(surface.areas().len()),
            &gate,
        )
        .is_err()
    );
}

fn exposure_limits(panels: usize) -> BoundaryWorkExposureLimits {
    BoundaryWorkExposureLimits {
        maximum_lines: 2,
        maximum_factors: 1,
        maximum_source_panel_products: 12 * panels,
        maximum_panel_frequency_pairs: 4 * panels,
        maximum_workspace_bytes: 1 << 20,
    }
}

#[test]
fn g1_boundary_work_exposure_uses_retained_material_traces_and_common_source_clock() {
    use fs_fft::harmonic_window::HarmonicWindow;
    let (surface, sources, traces, boundary) = fixture(0, Some([0.06, 0.0, 0.02]));
    let gate = CancelGate::new_clock_free();
    let owner = IterativeBoundaryObservation::for_surface(
        &surface,
        medium(),
        &sources,
        &traces,
        &receivers(),
        boundary,
        options(128),
        &gate,
    )
    .unwrap();
    let panels = surface.areas().len();
    let a = owner
        .solve_frequency(TAU * 700.0, 4 * panels, &gate)
        .unwrap();
    let b = owner
        .solve_frequency(TAU * 900.0, 4 * panels, &gate)
        .unwrap();
    let rates_a = [C64::new(3e-6, -2e-6), C64::new(-1e-6, 4e-6)];
    let rates_b = [C64::new(2e-6, 1e-6), C64::new(-3e-6, -2e-6)];
    let lines = [
        BoundaryWorkExposureLine {
            response: &a,
            rates: &rates_a,
        },
        BoundaryWorkExposureLine {
            response: &b,
            rates: &rates_b,
        },
    ];
    let limits = exposure_limits(panels);
    let prepared =
        PreparedBoundaryWorkExposure::new(&a, "drive-clock", &lines, &[], limits, &gate).unwrap();
    assert_eq!(prepared.sources(), owner.sources());
    assert_eq!(prepared.operator_hash(), owner.operator_hash());
    assert_eq!(prepared.source_clock(), "drive-clock");
    assert_eq!(prepared.panel_frequency_pairs(), 4 * panels);
    assert!(prepared.workspace_bytes() <= limits.maximum_workspace_bytes);
    let pa = a
        .coherent_boundary_power(&rates_a, u64::MAX, &gate)
        .unwrap();
    let pb = b
        .coherent_boundary_power(&rates_b, u64::MAX, &gate)
        .unwrap();
    // A full common period eliminates carrier and beat terms for 700/900 Hz.
    let cycle = prepared
        .observe_window(HarmonicWindow::new(0.00013, 0.01).unwrap(), &gate)
        .unwrap();
    let expected = [
        pa.input_power_w + pb.input_power_w,
        pa.outward_power_w + pb.outward_power_w,
        pa.absorbed_power_w + pb.absorbed_power_w,
    ];
    let actual = [
        cycle.input_power_w,
        cycle.outward_power_w,
        cycle.wall_power_w,
    ];
    for (got, want) in actual.into_iter().zip(expected) {
        assert!((got - want).abs() <= 2e-12 * expected[0].abs().max(1e-30));
    }
    let short = prepared
        .observe_window(HarmonicWindow::new(0.00013, 0.00017).unwrap(), &gate)
        .unwrap();
    assert!((short.input_power_w - cycle.input_power_w).abs() > 1e-4 * cycle.input_power_w.abs());
    assert!(PreparedBoundaryWorkExposure::new(&a, "", &lines, &[], limits, &gate).is_err());
    let tight = BoundaryWorkExposureLimits {
        maximum_source_panel_products: 12 * panels - 1,
        ..limits
    };
    assert!(
        PreparedBoundaryWorkExposure::new(&a, "drive-clock", &lines, &[], tight, &gate).is_err()
    );
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert!(
        PreparedBoundaryWorkExposure::new(&a, "drive-clock", &lines, &[], limits, &cancelled)
            .is_err()
    );
}
