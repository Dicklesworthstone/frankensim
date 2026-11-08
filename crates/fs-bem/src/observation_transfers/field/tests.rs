use super::*;
use crate::helmholtz::{IterativeOptions, laplace::complex::ComplexOptions};
use crate::near_field::Options;
use crate::observation_transfers::{
    HarmonicField, IterativeBoundaryObservation, IterativeObservationError,
    IterativeObservationOptions, LocalBoundary, LocalBoundaryObservation, LocalWall,
    MolecularBemObservation, MolecularObservationFormulation, PointSourceSpec,
    ScatteringObservation,
};
use crate::panel3d::SpherePanels;
use core::f64::consts::TAU;
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::{
    FluidAcousticState,
    thermoviscous::{ThermoviscousFluid, ThermoviscousFluidState},
};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::porous::{PorousCard, PorousFluidLaw};
use fs_material::relaxation::{RelaxingMedium, input::read_relaxation};

fn medium() -> Medium {
    Medium {
        density: 1.2,
        sound_speed: 340.0,
    }
}

fn window() -> LowerBandWindow {
    LowerBandWindow::try_new(150.0, 400.0).unwrap()
}

fn declaration(
    sources: &[SourceId],
    receivers: &[ReceiverSpec],
    delays: Vec<f64>,
    producer: ContentHash,
    uncertainty: ProducerUncertainty,
) -> TransferDeclaration {
    TransferDeclaration {
        producer_hash: producer,
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
        producer_uncertainty: uncertainty,
        maximum_frequency_entries: fs_vfit::observation::MAX_OBSERVATION_ENTRIES,
    }
}

fn phase(angle: f64) -> C64 {
    C64::new(angle.cos(), angle.sin())
}

fn norm(values: &[C64]) -> f64 {
    values.iter().map(|v| v.abs().powi(2)).sum::<f64>().sqrt()
}

// Independent exact outgoing monopole in the native negative-time convention.
// The p/u oracle uses the radial velocity directly, not the new conversion.
fn monopole(omega: f64, source: [f64; 3], receiver: [f64; 3]) -> [C64; 4] {
    let delta: [f64; 3] = std::array::from_fn(|a| receiver[a] - source[a]);
    let radius = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
    let k = omega / medium().sound_speed;
    let green = phase(k * radius).scale(1.0 / (2.0 * TAU * radius));
    let pressure = C64::new(0.0, -omega * medium().density) * green;
    let radial = green * C64::new(1.0 / radius, -k);
    [
        pressure,
        radial.scale(delta[0] / radius),
        radial.scale(delta[1] / radius),
        radial.scale(delta[2] / radius),
    ]
}

fn monopole_gradient(omega: f64, source: [f64; 3], receiver: [f64; 3]) -> [C64; 4] {
    let delta: [f64; 3] = std::array::from_fn(|a| receiver[a] - source[a]);
    let radius = delta.iter().map(|v| v * v).sum::<f64>().sqrt();
    let k = omega / medium().sound_speed;
    let pressure = monopole(omega, source, receiver)[0];
    let radial = pressure * C64::new(-1.0 / radius, k);
    [
        pressure,
        radial.scale(delta[0] / radius),
        radial.scale(delta[1] / radius),
        radial.scale(delta[2] / radius),
    ]
}

fn assert_transfer(bake: &PressureVelocityBake, omega: f64, field: &[C64], tolerance: f64) {
    assert_eq!(bake.kernels().len(), field.len());
    let latency =
        phase(-omega * bake.reconstruction_latency_samples() as f64 / bake.sample_rate_hz());
    for (kernels, expected) in bake.kernels().chunks_exact(4).zip(field.chunks_exact(4)) {
        let actual: Vec<_> = kernels
            .iter()
            .map(|k| k.transfer_at(omega / bake.sample_rate_hz()).unwrap())
            .collect();
        let expected: Vec<_> = expected.iter().map(|z| z.conj() * latency).collect();
        assert!(
            (actual[0] - expected[0]).abs() < tolerance * expected[0].abs().max(f64::MIN_POSITIVE)
        );
        let difference: Vec<_> = (1..4).map(|a| actual[a] - expected[a]).collect();
        assert!(norm(&difference) < tolerance * norm(&expected[1..]).max(f64::MIN_POSITIVE));
    }
}

fn harmonic_field(owner: &LocalBoundaryObservation<'_>, omega: f64) -> Vec<C64> {
    let native = owner.prepare_harmonic(omega, || false).unwrap();
    let mut fields = vec![HarmonicField::default(); owner.receivers().len()];
    native
        .evaluate_coherent(&[C64::ONE], &mut fields, || false)
        .unwrap();
    fields
        .iter()
        .flat_map(|f| {
            [
                f.pressure_pa,
                f.particle_velocity_m_s[0],
                f.particle_velocity_m_s[1],
                f.particle_velocity_m_s[2],
            ]
        })
        .collect()
}

#[test]
fn analytic_near_field_retains_velocity_at_a_pressure_null_and_changing_microphones() {
    let positions = [[-0.08, 0.0, 0.0], [0.08, 0.0, 0.0]];
    let receivers = [ReceiverSpec {
        id: "midplane".into(),
        position_m: [0.0; 3],
    }];
    let sources = [
        SourceId::PointVolumeFlow("left".into()),
        SourceId::PointVolumeFlow("right".into()),
    ];
    let producer = ContentHash([17; 32]);
    let d = declaration(
        &sources,
        &receivers,
        vec![0.08 / medium().sound_speed; 2],
        producer,
        ProducerUncertainty::ExactAnalytic,
    );
    let mut sampled_frequencies = Vec::new();
    let bake = bake_with_errors::<TransferError>(
        |omega| {
            sampled_frequencies.push(omega / TAU);
            Ok(positions
                .iter()
                .flat_map(|&source| monopole_gradient(omega, source, [0.0; 3]))
                .collect())
        },
        &d,
        medium(),
        window(),
        || false,
    )
    .unwrap();
    assert!(
        sampled_frequencies
            .iter()
            .all(|&f| f > window().stop_hz() && f < d.taper_stop_hz)
    );
    assert_ne!(bake.identity(), d.declaration_hash());
    assert_eq!(bake.producer_hash(), producer);
    assert_eq!(bake.sources(), sources);
    assert_eq!(bake.receivers(), receivers);
    assert_eq!(bake.geometric_delays_s(), d.geometric_delays_s);
    assert_eq!(bake.ladder().len(), 4);
    for frequency in [450.0, 500.0, 750.0] {
        let expected: Vec<_> = positions
            .iter()
            .flat_map(|&source| monopole(TAU * frequency, source, [0.0; 3]))
            .collect();
        assert_transfer(&bake, TAU * frequency, &expected, 3e-4);
    }
    let mut replay = bake.observer().unwrap();
    assert!(replay.current_fields().is_none());
    assert!(
        replay
            .microphone_pressure_pa(0, FirstOrder::default())
            .is_err()
    );
    let left = FirstOrder::new(0.0, [-1.0, 0.0, 0.0]).unwrap();
    let right = FirstOrder::new(0.0, [1.0, 0.0, 0.0]).unwrap();
    let prefix = (d.geometric_delays_s[0] * d.fs_hz).floor() as usize;
    let mut peak_velocity = 0.0_f64;
    for n in 0..700 {
        // Opposing co-clocked volume flows cancel pressure identically at the
        // midplane, while the two physical velocity vectors ADD there.
        let q = if n == 0 { 1e-6 } else { 0.0 };
        let frame = replay.step(&[q, -q], || false).unwrap()[0];
        assert_eq!(frame.pressure_pa, 0.0);
        assert_eq!(frame.particle_velocity_m_s[1..], [0.0; 2]);
        assert_eq!(frame.instantaneous_intensity_w_m2().unwrap(), [0.0; 3]);
        if n < prefix {
            assert_eq!(frame, PressureVelocityFrame::default());
        }
        peak_velocity = peak_velocity.max(frame.particle_velocity_m_s[0].abs());
        let expected = medium().density * medium().sound_speed * frame.particle_velocity_m_s[0];
        let a = replay.microphone_pressure_pa(0, left).unwrap();
        let b = replay.microphone_pressure_pa(0, right).unwrap();
        assert!((a - expected).abs() <= 2e-15 * expected.abs());
        assert_eq!(a, -b);
        assert_eq!(
            replay
                .microphone_pressure_pa(0, FirstOrder::default())
                .unwrap(),
            0.0
        );
        assert_eq!(replay.accepted_sample_count(), n + 1);
    }
    assert!(peak_velocity > 1e-7);
    // Failed inputs/cancellation leave the complete previously accepted frame
    // and shared source history available for an exact retry.
    let fields = replay.current_fields().unwrap().to_vec();
    let count = replay.accepted_sample_count();
    let mut reference = replay.clone();
    assert!(replay.step(&[f64::NAN, 0.0], || false).is_err());
    assert!(matches!(
        replay.step(&[0.0, 0.0], || true),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(replay.accepted_sample_count(), count);
    assert_eq!(replay.current_fields().unwrap(), fields);
    assert_eq!(
        replay.step(&[3e-7, -2e-7], || false).unwrap(),
        reference.step(&[3e-7, -2e-7], || false).unwrap()
    );
    assert!(replay.microphone_pressure_pa(1, left).is_err());

    // Single-source steady excitation retains reactive near-field motion and
    // the signed instantaneous energy exchange, not just p^2/(rho*c).
    let mut steady = bake.observer().unwrap();
    let omega = TAU * 500.0;
    let original = monopole(omega, positions[0], [0.0; 3]);
    let amplitude = 1e-6;
    let settled = bake
        .kernels()
        .iter()
        .map(|k| k.integer_delay_samples + k.taps.len())
        .max()
        .unwrap();
    let mut negative_flux = false;
    let mut mean = 0.0;
    for n in 0..settled + 64 {
        let frame = steady
            .step(
                &[amplitude * (omega * n as f64 / d.fs_hz).cos(), 0.0],
                || false,
            )
            .unwrap()[0];
        if n < settled {
            continue;
        }
        let clock = phase(-omega * (n as f64 - d.filter_latency_samples as f64) / d.fs_hz);
        let expected_p = amplitude * (original[0] * clock).re;
        let expected_u = amplitude * (original[1] * clock).re;
        assert!((frame.pressure_pa - expected_p).abs() < 3e-4 * amplitude * original[0].abs());
        assert!(
            (frame.particle_velocity_m_s[0] - expected_u).abs()
                < 3e-4 * amplitude * original[1].abs()
        );
        let intensity = frame.instantaneous_intensity_w_m2().unwrap()[0];
        assert!(
            (intensity - expected_p * expected_u).abs()
                < 7e-4 * amplitude.powi(2) * original[0].abs() * original[1].abs()
        );
        negative_flux |= intensity < 0.0;
        mean += intensity / 64.0;
    }
    assert!(
        negative_flux,
        "the reactive monopole field reverses instantaneous flux"
    );
    let expected_mean = 0.5 * amplitude.powi(2) * (original[0] * original[1].conj()).re;
    assert!((mean / expected_mean - 1.0).abs() < 1e-3);
}

fn porous(fluid: ThermoviscousFluidState) -> LocalWall {
    LocalWall::PorousRigid {
        law: PorousFluidLaw::try_new(
            PorousCard::basotect_foam02(),
            fluid.transport().dynamic_viscosity_pa_s(),
        )
        .unwrap(),
        fluid,
        thickness_m: 0.03,
    }
}

fn standard_fluid() -> ThermoviscousFluidState {
    ThermoviscousFluidState::from_heat_capacity_ratio(
        FluidAcousticState::try_new(medium().density, medium().sound_speed).unwrap(),
        FluidTransportState::try_new(medium().density, 1005.0, 0.0257, 1.81e-5).unwrap(),
        1.4,
    )
    .unwrap()
}

fn source() -> PointSourceSpec {
    PointSourceSpec {
        id: "accepted-flow".into(),
        position_m: [0.02, 0.01, 0.09],
    }
}

fn receiver() -> ReceiverSpec {
    ReceiverSpec {
        id: "material-shadow".into(),
        position_m: [-0.03, 0.02, -0.08],
    }
}

fn local(
    surface: &SpherePanels,
    lower: Option<LowerBandWindow>,
    law: LocalWall,
) -> LocalBoundaryObservation<'_> {
    let scattering = ScatteringObservation::new(
        surface,
        medium(),
        &[source()],
        &[receiver()],
        TAU * 1600.0,
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap();
    let boundary =
        LocalBoundary::try_new(medium(), vec![law], vec![0; surface.areas().len()]).unwrap();
    LocalBoundaryObservation::for_scattering(scattering, boundary, lower).unwrap()
}

fn declared_local(owner: &LocalBoundaryObservation<'_>) -> TransferDeclaration {
    declaration(
        owner.sources(),
        owner.receivers(),
        owner.geometric_delays_s(),
        owner.operator_hash(),
        ProducerUncertainty::Unresolved {
            source_hash: owner.operator_hash(),
        },
    )
}

#[test]
fn finite_porous_boundary_dense_and_iterative_fields_reach_the_same_typed_replay() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let dense = local(&surface, Some(window()), porous(standard_fluid()));
    let d = declared_local(&dense);
    let bake = dense.bake_pressure_velocity(&d, || false).unwrap();
    assert_eq!(bake.producer_uncertainty(), d.producer_uncertainty);
    assert_eq!(bake.lower_band_window(), window());
    assert_eq!(bake.pass_band_hz(), [400.0, 800.0]);
    let omega = TAU * 500.0;
    let reference = harmonic_field(&dense, omega);
    assert_transfer(&bake, omega, &reference, 3e-4);
    let direct = monopole(omega, source().position_m, receiver().position_m);
    assert!(
        (reference[0] - direct[0]).abs() > 1e-3 * direct[0].abs(),
        "finite wall diffraction must be active"
    );
    let gate = CancelGate::new_clock_free();
    let boundary = LocalBoundary::try_new_iterative(
        medium(),
        vec![porous(standard_fluid())],
        vec![0; surface.areas().len()],
    )
    .unwrap();
    let options = IterativeObservationOptions {
        highest_omega: TAU * 1600.0,
        minimum_clearance_m: 1e-3,
        receiver_quadrature: Options::default(),
        lower_window: Some(window()),
        solver: IterativeOptions {
            tol: 1e-9,
            restart: 60,
            max_iters: 600,
            fmm_order: 4,
            leaf_capacity: 128,
            ..IterativeOptions::default()
        },
    };
    let iterative = IterativeBoundaryObservation::for_scattering(
        &surface,
        medium(),
        &[source()],
        &[receiver()],
        boundary,
        options,
        &gate,
    )
    .unwrap();
    let id = declaration(
        iterative.sources(),
        iterative.receivers(),
        iterative.geometric_delays_s(),
        iterative.operator_hash(),
        ProducerUncertainty::Unresolved {
            source_hash: iterative.operator_hash(),
        },
    );
    let matrix_free = iterative.bake_pressure_velocity(&id, &gate).unwrap();
    assert_ne!(bake.identity(), matrix_free.identity());
    let (native, reports) = iterative.prepare_harmonic(omega, &gate).unwrap();
    assert_eq!(reports.len(), 1);
    assert!(reports[0].solve.converged_euclidean());
    let mut fields = [HarmonicField::default()];
    native
        .evaluate_coherent(&[C64::ONE], &mut fields, || false)
        .unwrap();
    let field = fields[0];
    let vector = [
        field.pressure_pa,
        field.particle_velocity_m_s[0],
        field.particle_velocity_m_s[1],
        field.particle_velocity_m_s[2],
    ];
    assert_transfer(&matrix_free, omega, &vector, 3e-4);
    assert!((vector[0] - reference[0]).abs() < 2e-7 * reference[0].abs());
    let mut a = bake.observer().unwrap();
    let mut b = matrix_free.observer().unwrap();
    let mut scale_p = 0.0_f64;
    let mut scale_u = 0.0_f64;
    let mut error_p = 0.0_f64;
    let mut error_u = 0.0_f64;
    let microphone = FirstOrder::new(0.5, [0.0, 0.0, 1.0]).unwrap();
    let mut microphone_peak = 0.0_f64;
    for n in 0..700 {
        let input = [if n == 0 { 1e-6 } else { 0.0 }];
        let left = a.step(&input, || false).unwrap()[0];
        let right = b.step(&input, || false).unwrap()[0];
        scale_p = scale_p.max(left.pressure_pa.abs());
        scale_u = scale_u.max(
            left.particle_velocity_m_s
                .iter()
                .map(|u| u * u)
                .sum::<f64>()
                .sqrt(),
        );
        error_p = error_p.max((left.pressure_pa - right.pressure_pa).abs());
        error_u = error_u.max(
            (0..3)
                .map(|j| (left.particle_velocity_m_s[j] - right.particle_velocity_m_s[j]).powi(2))
                .sum::<f64>()
                .sqrt(),
        );
        let expected = 0.5 * left.pressure_pa
            - 0.5 * medium().density * medium().sound_speed * left.particle_velocity_m_s[2];
        let observed = a.microphone_pressure_pa(0, microphone).unwrap();
        assert!((observed - expected).abs() < 1e-15 + 2e-14 * expected.abs());
        microphone_peak = microphone_peak.max(observed.abs());
        assert!(
            left.instantaneous_intensity_w_m2()
                .unwrap()
                .iter()
                .all(|x| x.is_finite())
        );
    }
    assert!(scale_p > 0.0 && scale_u > 0.0 && microphone_peak > 0.0);
    assert!(error_p < 2e-6 * scale_p);
    assert!(error_u < 2e-6 * scale_u);
}

#[test]
fn physical_owners_refuse_missing_window_cross_wiring_and_cancelled_realization() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let without = local(
        &surface,
        None,
        LocalWall::Mechanical(fs_phs::WallPin {
            surface_density: 0.1,
            stiffness_per_area: 4e5,
            resistance: 200.0,
        }),
    );
    assert!(
        without
            .bake_pressure_velocity(&declared_local(&without), || false)
            .is_err()
    );
    let owner = local(&surface, Some(window()), porous(standard_fluid()));
    let d = declared_local(&owner);
    assert!(matches!(
        owner.bake_pressure_velocity(&d, || true),
        Err(TransferError::Cancelled)
    ));
    let mut wrong = d.clone();
    wrong.receivers[0].position_m[0] += 0.01;
    assert!(owner.bake_pressure_velocity(&wrong, || false).is_err());
    wrong = d.clone();
    wrong.producer_hash = without.operator_hash();
    assert!(owner.bake_pressure_velocity(&wrong, || false).is_err());
    wrong = d.clone();
    wrong.geometric_delays_s[0] += 1.0 / wrong.fs_hz;
    assert!(owner.bake_pressure_velocity(&wrong, || false).is_err());
    wrong = d.clone();
    wrong.producer_uncertainty = ProducerUncertainty::ExactAnalytic;
    assert!(owner.bake_pressure_velocity(&wrong, || false).is_err());
    let mut calls = 0;
    assert!(matches!(
        bake_with_errors::<TransferError>(
            |_| {
                calls += 1;
                Ok(vec![C64::ZERO; 4])
            },
            &d,
            medium(),
            window(),
            || true
        ),
        Err(TransferError::Cancelled)
    ));
    assert_eq!(calls, 0);
    // The iterative adapter preserves the actual owner's native refusal,
    // including cancellation before any hypothetical successful field artifact.
    let gate = CancelGate::new_clock_free();
    let boundary = LocalBoundary::try_new_iterative(
        medium(),
        vec![LocalWall::Rigid],
        vec![0; surface.areas().len()],
    )
    .unwrap();
    let iterative = IterativeBoundaryObservation::for_scattering(
        &surface,
        medium(),
        &[source()],
        &[receiver()],
        boundary,
        IterativeObservationOptions {
            highest_omega: TAU * 1600.0,
            minimum_clearance_m: 1e-3,
            receiver_quadrature: Options::default(),
            solver: IterativeOptions::default(),
            lower_window: Some(window()),
        },
        &gate,
    )
    .unwrap();
    let declared = declaration(
        iterative.sources(),
        iterative.receivers(),
        iterative.geometric_delays_s(),
        iterative.operator_hash(),
        ProducerUncertainty::Unresolved {
            source_hash: iterative.operator_hash(),
        },
    );
    gate.request();
    assert!(matches!(
        iterative.bake_pressure_velocity(&declared, &gate),
        Err(IterativeObservationError::Transfer(
            TransferError::Cancelled
        ))
    ));
}

#[test]
fn molecular_gas_and_positive_band_porous_wall_keep_original_momentum_and_frozen_front() {
    let material = read_relaxation(
        &mut "declared 1 synthetic 1 vibration 1 900 1 power 101.325 293.15 0 200 500"
            .split_ascii_whitespace(),
    )
    .unwrap();
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101_325.0).unwrap();
    let state = material
        .evaluate(
            &RelaxingMedium::new(gas, vec![(material.species()[0].species.clone(), 1.0)]).unwrap(),
        )
        .unwrap();
    let gas = *state.gas();
    let medium = Medium::from(&gas);
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let panels = surface.areas().len();
    let trace: Vec<_> = (0..panels)
        .map(|i| if i < panels / 2 { 1.0 } else { 0.0 })
        .collect();
    let sources = [SourceId::StructuralMode(
        "original signed patch velocity".into(),
    )];
    let receivers = [ReceiverSpec {
        id: "relaxing receiver".into(),
        position_m: [0.09, -0.02, 0.03],
    }];
    let boundary = LocalBoundary::try_new(
        medium,
        vec![LocalWall::Rigid, porous(gas.thermoviscous_state().unwrap())],
        (0..panels).map(|i| usize::from(i >= panels / 2)).collect(),
    )
    .unwrap();
    let owner = MolecularBemObservation::new_with_formulation(
        &surface,
        state,
        &sources,
        &[trace],
        &receivers,
        1e-3,
        [TAU * 150.0, TAU * 1600.0],
        MolecularObservationFormulation::CentroidBurtonMiller,
        ComplexOptions {
            max_panels: panels,
            max_bytes: 32 * panels * panels,
            max_work: 8 * (panels * panels) as u64,
        },
    )
    .unwrap()
    .with_local_boundary(boundary)
    .unwrap();
    let d = declaration(
        &sources,
        &receivers,
        owner.geometric_delays_s(),
        owner.operator_hash(),
        ProducerUncertainty::Unresolved {
            source_hash: owner.operator_hash(),
        },
    );
    // This genuinely positive-band material wall remains unsupported by the
    // old unwindowed p/gradient bake; the new typed path supplies the missing
    // explicit window, rather than asking the material for a DC extrapolation.
    assert!(owner.bake(&d, || false).is_err());
    assert!(
        owner
            .bake_pressure_velocity(&d, LowerBandWindow::try_new(100.0, 400.0).unwrap(), || {
                false
            })
            .is_err()
    );
    let bake = owner
        .bake_pressure_velocity(&d, window(), || false)
        .unwrap();
    assert_eq!(bake.medium(), medium);
    assert_eq!(bake.geometric_delays_s(), owner.geometric_delays_s());
    assert_eq!(bake.producer_hash(), owner.operator_hash());
    assert_eq!(bake.producer_uncertainty(), d.producer_uncertainty);
    let omega = TAU * 500.0;
    let harmonic = owner.prepare_harmonic(omega, || false).unwrap();
    let mut fields = [HarmonicField::default()];
    harmonic
        .evaluate_coherent(&[C64::ONE], &mut fields, || false)
        .unwrap();
    let expected = [
        fields[0].pressure_pa,
        fields[0].particle_velocity_m_s[0],
        fields[0].particle_velocity_m_s[1],
        fields[0].particle_velocity_m_s[2],
    ];
    assert_transfer(&bake, omega, &expected, 3e-4);
    let mut replay = bake.observer().unwrap();
    let microphone = FirstOrder::new(0.0, [1.0, 0.0, 0.0]).unwrap();
    let mut peak = 0.0_f64;
    for n in 0..700 {
        let frame = replay
            .step(&[if n == 0 { 1e-4 } else { 0.0 }], || false)
            .unwrap()[0];
        let signal = replay.microphone_pressure_pa(0, microphone).unwrap();
        let original_calibration = -gas.characteristic_impedance * frame.particle_velocity_m_s[0];
        assert!(
            (signal - original_calibration).abs()
                < 2e-14 * original_calibration.abs().max(f64::MIN_POSITIVE)
        );
        peak = peak.max(signal.abs());
    }
    assert!(peak > 0.0);
}
