use super::*;
use core::f64::consts::TAU;
use fs_material::fluid_acoustics::{FluidAcousticState, thermoviscous::ThermoviscousFluidState};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_phs::WallPin;

fn medium() -> Medium {
    Medium {
        density: 1.2,
        sound_speed: 340.0,
    }
}
fn sources() -> [PointSourceSpec; 2] {
    [
        PointSourceSpec {
            id: "a".into(),
            position_m: [0.02, 0.01, 0.09],
        },
        PointSourceSpec {
            id: "b".into(),
            position_m: [-0.03, -0.01, 0.07],
        },
    ]
}
fn receivers() -> [ReceiverSpec; 2] {
    [
        ReceiverSpec {
            id: "shadow".into(),
            position_m: [-0.03, 0.02, -0.08],
        },
        ReceiverSpec {
            id: "side".into(),
            position_m: [0.1, 0.0, 0.03],
        },
    ]
}
fn options() -> IterativeObservationOptions {
    IterativeObservationOptions {
        highest_omega: TAU * 1600.0,
        minimum_clearance_m: 1e-3,
        receiver_quadrature: Options::default(),
        lower_window: None,
        solver: IterativeOptions {
            tol: 1e-9,
            restart: 60,
            max_iters: 600,
            fmm_order: 7,
            leaf_capacity: 16,
            ..IterativeOptions::default()
        },
    }
}
fn mechanical() -> LocalWall {
    LocalWall::Mechanical(WallPin {
        surface_density: 0.1,
        stiffness_per_area: 4e5,
        resistance: 200.0,
    })
}
fn uniform(surface: &SpherePanels, law: LocalWall) -> LocalBoundary {
    LocalBoundary::try_new_iterative(medium(), vec![law], vec![0; surface.areas().len()]).unwrap()
}
fn declaration(owner: &IterativeBoundaryObservation<'_>) -> TransferDeclaration {
    TransferDeclaration {
        producer_hash: owner.operator_hash(),
        sources: owner.sources().to_vec(),
        receivers: owner.receivers().to_vec(),
        fs_hz: 8000.0,
        geometric_delays_s: owner.geometric_delays_s(),
        pass_band_upper_hz: 800.0,
        taper_stop_hz: 1600.0,
        filter_latency_samples: 256,
        fft_lengths: vec![512, 1024, 2048, 4096],
        tail_energy_tolerance: 1e-10,
        dtft_relative_tolerance: 1e-4,
        producer_uncertainty: ProducerUncertainty::Unresolved {
            source_hash: owner.operator_hash(),
        },
        maximum_frequency_entries: MAX_OBSERVATION_ENTRIES,
    }
}

// Normalize pressure and the Cartesian gradient separately; never mix Pa and
// Pa/m into one relative-error denominator.
fn relative_field(actual: &[C64], reference: &[C64]) -> f64 {
    assert_eq!(actual.len(), reference.len());
    quantities(actual.len())
        .map(|range| {
            let denominator: f64 = reference[range.clone()].iter().map(|z| z.norm_sq()).sum();
            let numerator: f64 = actual[range.clone()]
                .iter()
                .zip(&reference[range])
                .map(|(&a, &b)| (a - b).norm_sq())
                .sum();
            det::sqrt(numerator / denominator.max(f64::MIN_POSITIVE))
        })
        .fold(0.0_f64, f64::max)
}

#[test]
fn three_fmm_orders_converge_to_dense_total_material_diffraction() {
    let surface = SpherePanels::icosphere(0.025, 2).unwrap();
    let gate = CancelGate::new_clock_free();
    let wall = uniform(&surface, mechanical());
    let dense = LocalBoundaryObservation::for_scattering(
        ScatteringObservation::new(
            &surface,
            medium(),
            &sources(),
            &receivers(),
            options().highest_omega,
            1e-3,
            Options::default(),
            || false,
        )
        .unwrap(),
        wall.clone(),
        None,
    )
    .unwrap();
    let omega = TAU * 700.0;
    let reference = dense.physical_frequency_sample(omega, || false).unwrap();
    let mut errors = Vec::new();
    let mut identities = Vec::new();
    for order in [4, 6, 7] {
        let mut config = options();
        config.solver.fmm_order = order;
        let owner = IterativeBoundaryObservation::for_scattering(
            &surface,
            medium(),
            &sources(),
            &receivers(),
            wall.clone(),
            config,
            &gate,
        )
        .unwrap();
        let sample = owner.physical_frequency_sample(omega, &gate).unwrap();
        assert_eq!(sample.sources.len(), 2);
        for (report, source) in sample.sources.iter().zip(owner.sources()) {
            assert_eq!(&report.source, source);
            assert!(report.solve.converged_euclidean() && report.solve.rel_residual <= 1e-9);
            assert_eq!(report.resources.fmm_order, order);
            assert!(report.resources.matvecs > 0);
            assert!(report.absorbed_power_w > 0.0);
        }
        errors.push(relative_field(&sample.transfer, &reference));
        identities.push(owner.operator_hash());
    }
    assert!(
        errors[2] < errors[0],
        "FMM order must change the field error: {errors:?}"
    );
    assert!(
        errors[2] < 2e-5,
        "finest complete p/gradient dense agreement: {errors:?}"
    );
    assert_ne!(identities[0], identities[1]);
    assert_ne!(identities[1], identities[2]);
}

#[test]
fn finite_prescribed_sources_share_iterative_boundary_and_keep_passive_work() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let gate = CancelGate::new_clock_free();
    let source_ids = [
        SourceId::ApertureVolumeFlow("a".into()),
        SourceId::ApertureVolumeFlow("b".into()),
    ];
    let mut traces = vec![vec![0.0; surface.areas().len()]; 2];
    for (index, trace) in traces.iter_mut().enumerate() {
        trace[index] = 1.0 / surface.areas()[index];
    }
    let mut assignments = vec![1; surface.areas().len()];
    assignments[..2].fill(0);
    let boundary = LocalBoundary::try_new_iterative(
        medium(),
        vec![LocalWall::Rigid, mechanical()],
        assignments,
    )
    .unwrap();
    let mut config = options();
    config.solver.leaf_capacity = 128; // Exact near pass isolates GMRES from FMM approximation.
    let iterative = IterativeBoundaryObservation::for_surface(
        &surface,
        medium(),
        &source_ids,
        &traces,
        &receivers(),
        boundary.clone(),
        config.clone(),
        &gate,
    )
    .unwrap();
    let dense = LocalBoundaryObservation::for_surface(
        BemObservation::new(&surface, medium(), &source_ids, &traces, &receivers(), 1e-3).unwrap(),
        boundary,
        config.highest_omega,
        Options::default(),
        None,
        || false,
    )
    .unwrap();
    let omega = TAU * 500.0;
    let actual = iterative.physical_frequency_sample(omega, &gate).unwrap();
    let expected = dense.physical_frequency_sample(omega, || false).unwrap();
    assert!(relative_field(&actual.transfer, &expected) < 1e-7);
    assert!(
        actual
            .sources
            .iter()
            .all(|source| source.absorbed_power_w > 0.0 && source.solve.converged_euclidean())
    );
    assert!(matches!(
        IterativeBoundaryObservation::for_surface(
            &surface,
            medium(),
            &source_ids,
            &traces,
            &receivers(),
            uniform(&surface, mechanical()),
            config,
            &gate
        ),
        Err(IterativeObservationError::Transfer(TransferError::Invalid(
            _
        )))
    ));
}

#[test]
fn failed_gmres_refuses_complete_samples_and_fir_bakes_with_native_diagnostics() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let gate = CancelGate::new_clock_free();
    let mut config = options();
    config.solver.max_iters = 0;
    config.solver.leaf_capacity = 128;
    let owner = IterativeBoundaryObservation::for_scattering(
        &surface,
        medium(),
        &sources()[..1],
        &receivers()[..1],
        uniform(&surface, mechanical()),
        config,
        &gate,
    )
    .unwrap();
    for result in [
        owner
            .physical_frequency_sample(TAU * 500.0, &gate)
            .map(|_| ()),
        owner.bake(&declaration(&owner), &gate).map(|_| ()),
    ] {
        let Err(IterativeObservationError::Solver(error)) = result else {
            panic!("expected original solve refusal: {result:?}");
        };
        let IterativeError::NotConverged {
            report, receipt, ..
        } = *error
        else {
            panic!("{error:?}");
        };
        assert!(!report.converged_euclidean());
        assert_eq!(receipt.iterations, 0);
    }
    gate.request();
    assert!(matches!(
        owner.physical_frequency_sample(0.0, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
}

#[test]
fn large_geometry_admission_preserves_exterior_work_and_cancellation_limits() {
    // 20,480 panels crosses 8,192 without constructing an FMM or dense matrix.
    let surface = SpherePanels::icosphere(0.025, 5).unwrap();
    assert!(surface.areas().len() > helmholtz::MAX_DENSE_PANELS);
    let gate = CancelGate::new_clock_free();
    let wall = uniform(&surface, LocalWall::Rigid);
    assert!(
        LocalBoundary::try_new(
            medium(),
            vec![LocalWall::Rigid],
            vec![0; surface.areas().len()]
        )
        .is_err()
    );
    let owner = IterativeBoundaryObservation::for_scattering(
        &surface,
        medium(),
        &sources()[..1],
        &receivers()[..1],
        wall,
        options(),
        &gate,
    )
    .unwrap();
    assert_eq!(owner.sources().len(), 1);
    assert!(
        ScatteringObservation::new(
            &surface,
            medium(),
            &sources()[..1],
            &receivers()[..1],
            options().highest_omega,
            1e-3,
            Options::default(),
            || false
        )
        .is_err()
    );
    assert!(Geometry::new_iterative(&surface, &[[0.0; 3]], 1e-3, &gate).is_err());
    let geometry =
        Geometry::new_iterative(&surface, &[receivers()[0].position_m], 1e-3, &gate).unwrap();
    let work_limited = Options {
        maximum_kernel_evaluations: 80,
        ..Options::default()
    };
    assert!(
        geometry
            .prepare_velocity_cancellable(
                TAU * 500.0 / medium().sound_speed,
                medium(),
                work_limited,
                &gate
            )
            .is_err()
    );
    gate.request();
    assert!(matches!(
        Geometry::new_iterative(&surface, &[receivers()[0].position_m], 1e-3, &gate),
        Err(HelmholtzError::Cancelled)
    ));
    assert!(matches!(
        geometry.prepare_velocity_cancellable(
            TAU * 500.0 / medium().sound_speed,
            medium(),
            Options::default(),
            &gate
        ),
        Err(HelmholtzError::Cancelled)
    ));
}

#[test]
fn iterative_porous_bake_keeps_explicit_band_and_replays_physical_pressure() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let gate = CancelGate::new_clock_free();
    let fluid = ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(1.2, 340.0).unwrap(),
        FluidTransportState::try_new(1.2, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap();
    let law = PorousFluidLaw::try_new(PorousCard::basotect_foam02(), 1.81e-5).unwrap();
    let boundary = uniform(
        &surface,
        LocalWall::PorousRigid {
            law,
            fluid,
            thickness_m: 0.03,
        },
    );
    let mut config = options();
    config.solver.leaf_capacity = 128;
    config.solver.fmm_order = 4;
    config.lower_window = Some(LowerBandWindow::try_new(150.0, 400.0).unwrap());
    let owner = IterativeBoundaryObservation::for_scattering(
        &surface,
        medium(),
        &sources()[..1],
        &receivers()[..1],
        boundary,
        config.clone(),
        &gate,
    )
    .unwrap();
    assert!(matches!(
        owner.physical_frequency_sample(TAU * 100.0, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Material(_)
        ))
    ));
    let stopped = owner
        .band_limited_frequency_sample(TAU * 100.0, &gate)
        .unwrap();
    assert_eq!(stopped.transfer, vec![C64::ZERO; 4]);
    assert!(stopped.sources.is_empty());
    let d = declaration(&owner);
    let bake = owner.bake(&d, &gate).unwrap();
    assert_eq!(bake.observation().lower_band_window(), config.lower_window);
    assert_eq!(bake.solver_options().fmm_order, config.solver.fmm_order);
    assert_eq!(
        bake.observation().transfer().declaration_hash,
        d.declaration_hash()
    );
    assert_eq!(bake.observation().transfer().ladder.len(), 4);
    let omega = TAU * 500.0;
    let reference = owner.physical_frequency_sample(omega, &gate).unwrap();
    let latency = phase(-omega * d.filter_latency_samples as f64 / d.fs_hz);
    for (kernel, sample) in bake
        .observation()
        .transfer()
        .kernels
        .iter()
        .zip(reference.transfer)
    {
        let expected = sample.conj() * latency;
        let measured = kernel.transfer_at(omega / d.fs_hz).unwrap();
        assert!((measured - expected).abs() < 2e-4 * expected.abs().max(1.0));
    }
    let mut runtime = bake.observer().unwrap();
    let mut peak = 0.0_f64;
    for index in 0..512 {
        let output = runtime
            .step(&[if index == 0 { 1e-6 } else { 0.0 }])
            .unwrap();
        assert!(output.iter().all(|value| value.is_finite()));
        peak = peak.max(output[0].abs());
    }
    assert!(peak > 1e-4);
}

#[test]
fn prepared_iterative_boundary_refuses_changed_nonlocal_admittance() {
    use crate::helmholtz::bc::BoundaryOperator;
    use std::sync::atomic::{AtomicU64, Ordering};
    struct MutableZero {
        panels: usize,
        stamp: AtomicU64,
    }
    impl BoundaryOperator for MutableZero {
        fn len(&self) -> usize {
            self.panels
        }
        fn fingerprint(&self) -> u64 {
            self.stamp.load(Ordering::Relaxed)
        }
        fn apply(&self, _pressure: &[C64], output: &mut [C64]) -> Result<(), HelmholtzError> {
            output.fill(C64::ZERO);
            Ok(())
        }
    }
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let gate = CancelGate::new_clock_free();
    let law = MutableZero {
        panels: surface.areas().len(),
        stamp: AtomicU64::new(1),
    };
    let prepared = PreparedIterativeBoundary::new(
        &surface,
        TAU * 500.0 / medium().sound_speed,
        medium(),
        Admittance::NonLocal {
            operator: &law,
            declared_active: false,
        },
        helmholtz::Formulation::BurtonMiller,
        &options().solver,
        &gate,
    )
    .unwrap();
    law.stamp.store(2, Ordering::Relaxed);
    assert!(matches!(
        prepared.solve(&vec![C64::ONE; surface.areas().len()], &gate),
        Err(IterativeError::Helmholtz(
            HelmholtzError::BadParameter { .. }
        ))
    ));
}
