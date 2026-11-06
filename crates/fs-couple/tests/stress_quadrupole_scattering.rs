//! Physical volume-stress moments through body diffraction and the original
//! harmonic pressure clock. Cross-owner tests keep all six tensor components.

use fs_aeroac::curle3d::quadrupole::{PreparedQuadrupoleObserver, QuadrupolePhasor};
use fs_aeroac::observation::{
    HarmonicFieldLine, HarmonicObservationBudget, HarmonicSampleBudget, HarmonicSampleClock,
    PreparedHarmonicField,
};
use fs_bem::helmholtz::bc::{
    Admittance, IncidentField, IncidentSource, StressQuadrupole, incident,
};
use fs_bem::helmholtz::{Formulation, IterativeOptions, Medium};
use fs_bem::near_field::Options;
use fs_bem::observation_transfers::{
    BemObservation, CompactSourceSpec, HarmonicField, HarmonicSourceMapLimits,
    IterativeBoundaryObservation, IterativeObservationOptions, LocalBoundary,
    LocalBoundaryObservation, LocalWall, LowerBandWindow, MAX_OBSERVATION_FFT,
    ProducerUncertainty, ReceiverSpec, ScatteringObservation, SourceId, StressQuadrupoleSpec,
    TransferDeclaration, bake_observation_transfers,
};
use fs_bem::panel3d::SpherePanels;
use fs_blake3::ContentHash;
use fs_couple::pcm_wav::observation::{PressureRenderer, harmonic::HarmonicPressureRenderer};
use fs_exec::CancelGate;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::{c64::C64, dd::Dd};
use std::f64::consts::TAU;

const MEDIUM: Medium = Medium {
    density: 1.2,
    sound_speed: 340.0,
};
const POSITIONS: [[f64; 3]; 2] = [[0.02, 0.01, 0.09], [-0.03, -0.01, 0.07]];
const COMPONENTS: [&str; 6] = ["xx", "yy", "zz", "xy", "xz", "yz"];
const FREQUENCY: f64 = 700.0;

fn close(actual: C64, expected: C64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-18),
        "actual {actual:?}, expected {expected:?}"
    );
}
fn moments() -> [QuadrupolePhasor; 2] {
    [0.001, -0.0006].map(|scale| QuadrupolePhasor {
        stress_moment_n_m: [
            C64::new(0.7, 0.2),
            C64::new(-0.3, 0.4),
            C64::new(0.9, -0.1),
            C64::new(0.23, -0.3),
            C64::new(-0.11, 0.3),
            C64::new(0.4, 0.2),
        ]
        .map(|q| q.scale(scale)),
    })
}
fn incidents(moments: &[QuadrupolePhasor; 2]) -> [IncidentSource; 2] {
    core::array::from_fn(|i| {
        IncidentSource::StressQuadrupole(StressQuadrupole {
            position_m: POSITIONS[i],
            stress_moment_n_m: moments[i].stress_moment_n_m,
        })
    })
}
fn specs() -> Vec<CompactSourceSpec> {
    let mut sources = Vec::new();
    for (cell, position_m) in POSITIONS.into_iter().enumerate() {
        for (component, name) in COMPONENTS.into_iter().enumerate() {
            sources.push(CompactSourceSpec::StressQuadrupole(StressQuadrupoleSpec {
                id: format!("cell-{cell}:{name}"),
                position_m,
                tensor_shape: core::array::from_fn(|i| if i == component { 1.0 } else { 0.0 }),
            }));
        }
    }
    sources
}
fn receivers() -> [ReceiverSpec; 2] {
    [
        ReceiverSpec {
            id: "shadow".into(),
            position_m: [-0.03, 0.02, -0.08],
        },
        ReceiverSpec {
            id: "far".into(),
            position_m: [1.2, -0.7, 0.4],
        },
    ]
}
fn observation<'a>(
    surface: &'a SpherePanels,
    sources: &[CompactSourceSpec],
) -> ScatteringObservation<'a> {
    ScatteringObservation::new_mixed(
        surface,
        MEDIUM,
        sources,
        &receivers(),
        TAU * 1600.0,
        1e-3,
        Options::default(),
        || false,
    )
    .unwrap()
}
fn rigid(surface: &SpherePanels) -> LocalBoundary {
    LocalBoundary::try_new(MEDIUM, vec![LocalWall::Rigid], vec![0; surface.areas().len()])
        .unwrap()
}
fn declaration(
    hash: ContentHash,
    sources: &[SourceId],
    delays: Vec<f64>,
) -> TransferDeclaration {
    TransferDeclaration {
        producer_hash: hash,
        sources: sources.to_vec(),
        receivers: receivers().to_vec(),
        fs_hz: 8000.0,
        geometric_delays_s: delays,
        pass_band_upper_hz: 800.0,
        taper_stop_hz: 1600.0,
        filter_latency_samples: 256,
        fft_lengths: vec![512, 1024, 2048, 4096],
        tail_energy_tolerance: 1e-10,
        dtft_relative_tolerance: 1e-4,
        producer_uncertainty: ProducerUncertainty::Unresolved { source_hash: hash },
        maximum_frequency_entries: 8 * sources.len() * MAX_OBSERVATION_FFT,
    }
}

#[test]
fn g1_original_aeroac_stress_phasors_match_bem_free_field_without_rescaling() {
    let moments = moments();
    let sources = incidents(&moments);
    let field = IncidentField::Sources(&sources);
    for frequency_hz in [30.0, FREQUENCY] {
        for receiver in receivers() {
            let acoustic = FluidAcousticState::try_new(MEDIUM.density, MEDIUM.sound_speed).unwrap();
            let original = PreparedQuadrupoleObserver::new(
                acoustic,
                frequency_hz,
                receiver.position_m,
                &POSITIONS,
                2,
            )
            .unwrap()
            .evaluate(&moments)
            .unwrap();
            let actual = field
                .pressure_velocity_at(
                    TAU * frequency_hz / MEDIUM.sound_speed,
                    MEDIUM,
                    receiver.position_m,
                )
                .unwrap();
            close(actual.pressure_pa, original.pressure_pa, 2e-12);
            for axis in 0..3 {
                close(
                    actual.velocity_m_s[axis],
                    original.particle_velocity_m_s[axis],
                    2e-12,
                );
            }
        }
    }
}

#[test]
fn g3_tensor_components_scatter_coherently_and_reach_the_original_pressure_clock() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let omega = TAU * FREQUENCY;
    let moments = moments();
    let sources = incidents(&moments);
    let field = IncidentField::Sources(&sources);
    let solved = incident::solve_scattering(
        &surface,
        omega / MEDIUM.sound_speed,
        MEDIUM,
        Admittance::Rigid,
        field,
        1e-3,
        Formulation::BurtonMiller,
    )
    .unwrap();
    let trace = solved.solution();
    let velocity_scale = trace
        .incident_normal_velocity_m_s
        .iter()
        .map(|v| v.abs())
        .fold(0.0_f64, f64::max);
    assert!(velocity_scale > 0.0);
    for total in trace.total_normal_velocity_m_s() {
        assert!(
            total.abs() < 2e-12 * velocity_scale,
            "rigid wall must cancel normal flow"
        );
    }
    let points = receivers().map(|r| r.position_m);
    let expected = solved
        .total_pressure_velocity_at(&points, 1e-3, Options::default())
        .unwrap();
    let owner = observation(&surface, &specs());
    let bank = owner.prepare_harmonic(omega, || false).unwrap();
    let coefficients: Vec<_> = moments.iter().flat_map(|q| q.stress_moment_n_m).collect();
    assert_eq!(bank.sources().len(), coefficients.len());
    for (index, id) in bank.sources().iter().enumerate() {
        assert_eq!(
            id,
            &SourceId::PointStressQuadrupole(format!(
                "cell-{}:{}",
                index / 6,
                COMPONENTS[index % 6]
            ))
        );
    }
    let mut received = [HarmonicField::default(); 2];
    bank.evaluate_coherent(&coefficients, &mut received, || false)
        .unwrap();
    for (index, value) in received.iter().enumerate() {
        close(value.pressure_pa, expected.scalar.pressure[index], 2e-11);
        for axis in 0..3 {
            close(
                value.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[index][axis],
                2e-11,
            );
        }
    }
    let direct = field
        .pressure_velocity_at(omega / MEDIUM.sound_speed, MEDIUM, points[0])
        .unwrap();
    assert!(
        (received[0].pressure_pa - direct.pressure_pa).abs() > 1e-4 * direct.pressure_pa.abs(),
        "finite-body scattering must reach the receiver"
    );

    // A retained source map must preserve the typed tensor coordinates and
    // their authored order; applying the same permutation leaves physics fixed.
    let ids: Vec<_> = bank.sources().iter().rev().cloned().collect();
    let n = ids.len();
    let mut map = vec![C64::ZERO; n * n];
    for index in 0..n {
        map[index * n + n - 1 - index] = C64::ONE;
    }
    let remapped = bank
        .compose_sources(
            &ids,
            &map,
            HarmonicSourceMapLimits {
                maximum_sources: n,
                maximum_products: 4096,
                maximum_storage_bytes: 1 << 20,
            },
            || false,
        )
        .unwrap();
    assert_eq!(remapped.sources(), ids.as_slice());
    let reversed: Vec<_> = coefficients.iter().rev().copied().collect();
    let mut reordered = [HarmonicField::default(); 2];
    remapped
        .evaluate_coherent(&reversed, &mut reordered, || false)
        .unwrap();
    for (actual, original) in reordered.iter().zip(&received) {
        close(actual.pressure_pa, original.pressure_pa, 2e-12);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                original.particle_velocity_m_s[axis],
                2e-12,
            );
        }
    }

    // The matrix-free owner consumes the identical declared tensor basis.
    let gate = CancelGate::new_clock_free();
    let iterative = IterativeBoundaryObservation::for_mixed_scattering(
        &surface,
        MEDIUM,
        &specs(),
        &receivers(),
        rigid(&surface),
        IterativeObservationOptions {
            highest_omega: TAU * 1600.0,
            minimum_clearance_m: 1e-3,
            receiver_quadrature: Options::default(),
            solver: IterativeOptions {
                tol: 1e-10,
                restart: 40,
                max_iters: 200,
                leaf_capacity: 64,
                fmm_order: 4,
                ..IterativeOptions::default()
            },
            lower_window: None,
        },
        &gate,
    )
    .unwrap();
    assert!(iterative.physical_frequency_sample(0.0, &gate).is_err());
    let (iterative_field, _) = iterative.prepare_harmonic(omega, &gate).unwrap();
    let mut iterative_received = [HarmonicField::default(); 2];
    iterative_field
        .evaluate_coherent(&coefficients, &mut iterative_received, || false)
        .unwrap();
    for (actual, original) in iterative_received.iter().zip(&received) {
        close(actual.pressure_pa, original.pressure_pa, 2e-7);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                original.particle_velocity_m_s[axis],
                2e-7,
            );
        }
    }

    // No propagation phase, rho, i*omega or normalization is added by output.
    let pressure = PreparedHarmonicField::new(
        &[HarmonicFieldLine {
            frequency_hz: FREQUENCY,
            pressure_pa: received[0].pressure_pa,
            particle_velocity_m_s: received[0].particle_velocity_m_s,
        }],
        HarmonicObservationBudget {
            max_lines: 1,
            max_line_pairs: 1,
        },
    )
    .unwrap();
    let start_s = 0.013;
    let clock = HarmonicSampleClock::new(Dd::from_f64(start_s), 8000).unwrap();
    let mut renderer = HarmonicPressureRenderer::new(
        &pressure,
        clock,
        HarmonicSampleBudget {
            max_samples: 64,
            max_tone_samples: 64,
        },
    )
    .unwrap();
    let mut samples = [0.0; 64];
    renderer.block(&mut samples).unwrap();
    let p = received[0].pressure_pa;
    for (index, &sample) in samples.iter().enumerate() {
        let (sin, cos) = (omega * (start_s + index as f64 / 8000.0)).sin_cos();
        assert!((sample - (p.re * cos + p.im * sin)).abs() < 2e-11 * p.abs());
    }
}

#[test]
fn g0_tensor_shapes_identity_and_static_pressure_are_not_silently_reinterpreted() {
    let surface = SpherePanels::icosphere(0.025, 0).unwrap();
    let source = |shape| [CompactSourceSpec::StressQuadrupole(StressQuadrupoleSpec {
        id: "shear-xy".into(),
        position_m: POSITIONS[0],
        tensor_shape: shape,
    })];
    let unit = [0.0, 0.0, 0.0, 1.0, 0.0, 0.0];
    let owner = observation(&surface, &source(unit));
    let doubled = observation(&surface, &source(unit.map(|v| 2.0 * v)));
    assert_ne!(owner.operator_hash(), doubled.operator_hash());
    let base = owner.frequency_sample(TAU * FREQUENCY, || false).unwrap();
    let twice = doubled.frequency_sample(TAU * FREQUENCY, || false).unwrap();
    for (&a, &b) in base.iter().zip(&twice) {
        close(b, a.scale(2.0), 2e-12);
    }
    for shape in [[0.0; 6], [f64::NAN; 6]] {
        assert!(
            ScatteringObservation::new_mixed(
                &surface, MEDIUM, &source(shape), &receivers(), TAU * 1600.0,
                1e-3, Options::default(), || false,
            ).is_err()
        );
    }
    assert!(owner.frequency_sample(0.0, || false).is_err());
    let d = declaration(owner.operator_hash(), owner.sources(), owner.geometric_delays_s());
    d.admit_work().unwrap();
    assert!(owner.bake(&d, || false).is_err());
    let mut called = false;
    assert!(
        bake_observation_transfers(
            |_| {
                called = true;
                Ok(vec![C64::ZERO; 8])
            },
            &d,
            || false,
        ).is_err()
    );
    assert!(!called, "generic FIR zero bin cannot erase static stress pressure");
    let mut wrong_kind = d.clone();
    wrong_kind.sources[0] = SourceId::PointForce("shear-xy".into());
    assert_ne!(wrong_kind.declaration_hash(), d.declaration_hash());
    assert!(
        BemObservation::new(
            &surface, MEDIUM, owner.sources(),
            &[vec![1.0; surface.areas().len()]], &receivers(), 1e-3,
        ).is_err(),
        "stress moments are not surface-velocity traces"
    );
    let unwindowed = LocalBoundaryObservation::for_scattering(
        observation(&surface, &source(unit)), rigid(&surface), None,
    ).unwrap();
    let d = declaration(
        unwindowed.operator_hash(), unwindowed.sources(), unwindowed.geometric_delays_s(),
    );
    assert!(unwindowed.physical_frequency_sample(0.0, || false).is_err());
    assert!(unwindowed.bake(&d, || false).is_err());
    let windowed = LocalBoundaryObservation::for_scattering(
        owner,
        rigid(&surface),
        Some(LowerBandWindow::try_new(150.0, 400.0).unwrap()),
    ).unwrap();
    assert!(windowed.physical_frequency_sample(0.0, || false).is_err());
    for frequency in [0.0, 100.0] {
        assert_eq!(
            windowed
                .band_limited_frequency_sample(TAU * frequency, || false)
                .unwrap(),
            vec![C64::ZERO; 8],
        );
    }
    let passed = windowed
        .band_limited_frequency_sample(TAU * FREQUENCY, || false)
        .unwrap();
    for (&actual, &expected) in passed.iter().zip(&base) {
        close(actual, expected, 2e-11);
    }
}
