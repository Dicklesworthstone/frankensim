use super::*;
use crate::curle3d::{CoherenceBudget, CompactSourcePhasor, PreparedCompactSourceObserver};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_flux::compressible::spatial::{AxisBoundaries, StepControl};
use fs_gas::dynamics::CaloricDomain;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_math::c64::C64;
use fs_qty::{Temperature, Time};
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, cantera_h2o2};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const ROTATION: [[f64; 3]; 3] = [[0.48, -0.6, 0.64], [0.36, 0.8, 0.48], [-0.8, 0.0, 0.6]];

fn frozen(species: usize) -> FrozenGas {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let model = mechanism.thermodynamics().models()[species].clone();
    FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![model],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap()
}

fn config(gas: &FrozenGas) -> CompressibleVolumeConfig3 {
    let rho = 1.7;
    let state = State::from_primitive(gas, rho, [0.0; 3], rho * gas.gas_constant() * 500.0)
        .unwrap()
        .primitive(gas, || false)
        .unwrap();
    CompressibleVolumeConfig3 {
        acoustic: FluidAcousticState::try_new(rho, state.sound_speed).unwrap(),
        reference_pressure_pa: state.pressure,
        origin_m: [0.2, -0.3, 0.1],
        grid_to_acoustic: ROTATION,
    }
}

fn controls(faces: usize, species: usize) -> PermeableSurfaceControls3 {
    PermeableSurfaceControls3 {
        maximum_faces: faces,
        maximum_workspace_bytes: PreparedPermeableSurface3::required_workspace_bytes(
            faces, species,
        )
        .unwrap(),
        reference_relative_tolerance: 2e-9,
    }
}

fn uniform(gas: FrozenGas, grid: CartesianGrid) -> Euler3d {
    let state = State::from_primitive(
        &gas,
        2.4,
        [2.0, -3.0, 4.0],
        2.4 * gas.gas_constant() * 620.0,
    )
    .unwrap();
    let count = grid.shape().into_iter().product();
    Euler3d::new(gas, grid, vec![state; count], PERIODIC, || false).unwrap()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * (1.0 + expected.abs()),
        "{actual:e} != {expected:e}"
    );
}
fn magnitude(vector: [C64; 3]) -> f64 {
    vector
        .into_iter()
        .map(|z| z.abs().powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn uniform_throughflow_retains_density_convected_momentum_and_closed_surface_cancellation() {
    let gas = frozen(6);
    let settings = config(&gas);
    let grid = CartesianGrid::new(5, 6, 4, 0.02, 0.03, 0.04).unwrap();
    let control_box = PermeableBox3 {
        lower: [1, 2, 1],
        upper: [4, 5, 3],
    };
    let count = control_box.face_count(&grid).unwrap();
    assert_eq!(count, 42);
    let mut flow = uniform(gas, grid);
    let mut source = PreparedPermeableSurface3::from_frozen(
        &flow,
        control_box,
        settings,
        controls(count, 0),
        || false,
    )
    .unwrap();
    let states_before = flow.states().to_vec();
    let native = flow.states()[0].primitive(flow.gas(), || false).unwrap();
    let mut total = [Dd::from_f64(0.0); 4];
    for ((face, position), actual) in source
        .faces()
        .iter()
        .zip(source.source_positions_m())
        .zip(source.source_components())
    {
        let axis = face.side.axis();
        let normal_flux = native.density * native.velocity[axis];
        let signed_area = face.side.sign() * face.area_m2;
        close(
            actual[0],
            signed_area * normal_flux / settings.acoustic.density_kg_m3(),
            2e-13,
        );
        let local: [f64; 3] = std::array::from_fn(|d| {
            normal_flux * native.velocity[d]
                + if d == axis {
                    native.pressure - settings.reference_pressure_pa
                } else {
                    0.0
                }
        });
        for (d, row) in ROTATION.iter().enumerate() {
            let expected: f64 = row
                .iter()
                .zip(local)
                .map(|(r, f)| r * f * signed_area)
                .sum();
            close(actual[d + 1], expected, 2e-13);
        }
        let left = grid.cell_center(face.adjacent_cells[0]).unwrap();
        let right = grid.cell_center(face.adjacent_cells[1]).unwrap();
        for (d, row) in ROTATION.iter().enumerate() {
            let center: f64 = row
                .iter()
                .zip(left.into_iter().zip(right))
                .map(|(r, (l, h))| r * (0.5 * l + 0.5 * h))
                .sum();
            close(position[d], settings.origin_m[d] + center, 1e-15);
        }
        for a in 0..4 {
            total[a] = total[a] + Dd::from_f64(actual[a]);
        }
    }
    // A uniform closed-surface throughflow has no net injection or resultant
    // loading. Each face still retains its full physical source before summing.
    assert!(total.into_iter().all(|v| v.to_f64().abs() < 1e-12));
    assert_eq!(flow.states(), states_before);
    assert_eq!(source.frozen_source(), Some(flow.gas()));
    assert!(source.mixture_source().is_none());
    assert_eq!(source.config(), settings);
    assert_eq!(source.control_box(), control_box);
    let initial = source.source_components().to_vec();
    flow.step(Time::new(1e-6), StepControl::default(), || false)
        .unwrap();
    source.sample_frozen(&flow, || false).unwrap();
    assert_eq!(source.sample_time_s(), flow.time().value());
    for (actual, expected) in source.source_components().iter().zip(initial) {
        for a in 0..4 {
            close(actual[a], expected[a], 2e-11);
        }
    }
}

// Independent outgoing monopole and an incident plane wave. The former must
// be reconstructed outside an enclosing surface; the latter has no enclosed
// source and its complete surface integral must cancel there.
fn analytic_field(point: [f64; 3], acoustic: FluidAcousticState, plane: bool) -> (C64, [C64; 3]) {
    let rho = acoustic.density_kg_m3();
    let c = acoustic.sound_speed_m_s();
    let k = 4.0;
    if plane {
        let direction = [0.6, 0.0, 0.8];
        let phase = k * (point[0] * direction[0] + point[2] * direction[2]);
        let p = C64::new(phase.cos(), phase.sin());
        return (p, direction.map(|d| p.scale(d / (rho * c))));
    }
    let radius = point.into_iter().map(|x| x * x).sum::<f64>().sqrt();
    let phase = k * radius;
    let g = C64::new(phase.cos(), phase.sin()).scale(1.0 / (4.0 * core::f64::consts::PI * radius));
    let q = 1e-5;
    let p = g * C64::new(0.0, -k * c * rho * q);
    let radial = g * C64::new(1.0 / radius, -k).scale(q);
    (p, point.map(|x| radial.scale(x / radius)))
}

fn surface_harmonic_error(cells_per_side: usize, plane: bool) -> f64 {
    let gas = frozen(6);
    let mut settings = config(&gas);
    let half_width = 0.35;
    let h = 2.0 * half_width / cells_per_side as f64;
    settings.origin_m = [-half_width - h; 3];
    settings.grid_to_acoustic = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
    let side = cells_per_side + 2;
    let grid = CartesianGrid::new(side, side, side, h, h, h).unwrap();
    let control_box = PermeableBox3 {
        lower: [1; 3],
        upper: [cells_per_side + 1; 3],
    };
    let count = control_box.face_count(&grid).unwrap();
    let mut peaks = vec![[C64::ZERO; 4]; count];
    let mut positions = Vec::new();
    // Four independent, accepted manufactured states recover the fundamental
    // of the complete nonlinear surface products. The native-flow history test
    // below separately checks advancement and its actual accepted time.
    for phase_index in 0..4 {
        let phase = core::f64::consts::FRAC_PI_2 * phase_index as f64;
        let real = |z: C64| z.re * phase.cos() + z.im * phase.sin();
        let states = (0..side.pow(3))
            .map(|cell| {
                let local = grid.cell_center(cell).unwrap();
                let point = std::array::from_fn(|a| settings.origin_m[a] + local[a]);
                let (p, u) = analytic_field(point, settings.acoustic, plane);
                let pressure = settings.reference_pressure_pa + real(p);
                let density = settings.acoustic.density_kg_m3()
                    + real(p) / settings.acoustic.sound_speed_m_s().powi(2);
                State::from_primitive(&gas, density, u.map(real), pressure).unwrap()
            })
            .collect();
        let flow = Euler3d::new(gas.clone(), grid, states, PERIODIC, || false).unwrap();
        let source = PreparedPermeableSurface3::from_frozen(
            &flow,
            control_box,
            settings,
            controls(count, 0),
            || false,
        )
        .unwrap();
        if phase_index == 0 {
            positions = source.source_positions_m().to_vec();
        }
        let weight = C64::new(phase.cos(), phase.sin()).scale(0.5);
        for (peak, sample) in peaks.iter_mut().zip(source.source_components()) {
            for a in 0..4 {
                peak[a] = peak[a] + weight.scale(sample[a]);
            }
        }
    }
    let peaks: Vec<_> = peaks
        .into_iter()
        .map(|v| CompactSourcePhasor {
            volume_flow_m3_s: v[0],
            force_on_fluid_n: [v[1], v[2], v[3]],
        })
        .collect();
    let receiver = [1.1, 0.7, -0.6];
    let observer = PreparedCompactSourceObserver::new(
        settings.acoustic,
        4.0 * settings.acoustic.sound_speed_m_s() / core::f64::consts::TAU,
        receiver,
        &positions,
        count,
    )
    .unwrap();
    let actual = observer.evaluate(&peaks).unwrap();
    let (incident_p, incident_u) = analytic_field(receiver, settings.acoustic, plane);
    let (wanted_p, wanted_u) = if plane {
        (C64::ZERO, [C64::ZERO; 3])
    } else {
        (incident_p, incident_u)
    };
    let p_error = (actual.pressure_pa - wanted_p).abs() / incident_p.abs();
    let u_error = magnitude(std::array::from_fn(|a| {
        actual.particle_velocity_m_s[a] - wanted_u[a]
    })) / magnitude(incident_u);
    p_error.max(u_error)
}

#[test]
fn closed_surface_recovers_radiating_monopole_p_and_vector_u_and_cancels_an_external_wave() {
    for plane in [false, true] {
        let errors = [4, 8, 16].map(|n| surface_harmonic_error(n, plane));
        assert!(errors[1] < 0.29 * errors[0], "{plane}: {errors:?}");
        assert!(errors[2] < 0.27 * errors[1], "{plane}: {errors:?}");
        assert!(
            errors[2] < if plane { 2.3e-4 } else { 1.1e-3 },
            "{plane}: {errors:?}"
        );
    }
}

#[test]
fn mixed_species_use_their_own_calorics_with_an_explicit_checked_ambient_composition() {
    let frozen = frozen(6);
    let settings = config(&frozen);
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let gas = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap();
    let species = gas.species_count();
    let grid = CartesianGrid::new(4, 4, 4, 0.02, 0.03, 0.04).unwrap();
    let mut reference = vec![0.0; species];
    reference[6] = settings.acoustic.density_kg_m3();
    let mut values = vec![0.0; 64 * (species + 4)];
    for (cell, out) in values.chunks_exact_mut(species + 4).enumerate() {
        let x = (cell % 4) as f64;
        let y = ((cell / 4) % 4) as f64;
        let mut partial = vec![0.0; species];
        partial[0] = 0.2 + 0.05 * x;
        partial[6] = 1.1 + 0.1 * y;
        state_into(
            &gas,
            &partial,
            [2.0 + 0.1 * x, -1.0 + 0.2 * y, 0.5],
            settings.reference_pressure_pa + 100.0 * x,
            out,
        )
        .unwrap();
    }
    let workspace = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    let mut flow = MixtureEuler3d::new(gas, grid, values, PERIODIC, workspace, || false).unwrap();
    let control_box = PermeableBox3 {
        lower: [1; 3],
        upper: [3; 3],
    };
    let count = control_box.face_count(&grid).unwrap();
    let mut source = PreparedPermeableSurface3::from_mixture(
        &flow,
        &reference,
        control_box,
        settings,
        controls(count, species),
        || false,
    )
    .unwrap();
    assert_eq!(source.reference_partial_densities_kg_m3(), reference);
    assert_eq!(source.mixture_source(), Some(flow.gas()));
    assert!(source.frozen_source().is_none());
    for (face, actual) in source.faces().iter().zip(source.source_components()) {
        let mut expected = [0.0; 4];
        for &cell in &face.adjacent_cells {
            let p = flow
                .cell(cell)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap();
            let normal = p.velocity[face.side.axis()] * face.side.sign();
            expected[0] +=
                0.5 * face.area_m2 * p.density * normal / settings.acoustic.density_kg_m3();
            for (a, row) in ROTATION.iter().enumerate() {
                for (b, r) in row.iter().enumerate() {
                    let pressure = if b == face.side.axis() {
                        (p.pressure - settings.reference_pressure_pa) * face.side.sign()
                    } else {
                        0.0
                    };
                    expected[a + 1] +=
                        0.5 * face.area_m2 * r * (pressure + p.density * p.velocity[b] * normal);
                }
            }
        }
        for a in 0..4 {
            close(actual[a], expected[a], 3e-12);
        }
    }
    let mut wrong = reference.clone();
    wrong[6] *= 1.01;
    assert!(
        PreparedPermeableSurface3::from_mixture(
            &flow,
            &wrong,
            control_box,
            settings,
            controls(count, species),
            || false
        )
        .is_err()
    );
    assert!(
        PreparedPermeableSurface3::from_mixture(
            &flow,
            &reference[..species - 1],
            control_box,
            settings,
            controls(count, species),
            || false
        )
        .is_err()
    );
    flow.step(Time::new(1e-6), StepControl::default(), || false)
        .unwrap();
    source.sample_mixture(&flow, || false).unwrap();
    assert_eq!(source.sample_time_s(), flow.time().value());
    assert_eq!(source.reference_partial_densities_kg_m3(), reference);
}

#[test]
fn source_medium_closed_geometry_budgets_clock_and_late_cancellation_refuse_atomically() {
    let gas = frozen(6);
    let settings = config(&gas);
    let grid = CartesianGrid::new(4, 4, 4, 0.02, 0.03, 0.04).unwrap();
    let mut flow = uniform(gas, grid);
    let control_box = PermeableBox3 {
        lower: [1; 3],
        upper: [3; 3],
    };
    let count = control_box.face_count(&grid).unwrap();
    let limits = controls(count, 0);
    let build = |settings, limits| {
        PreparedPermeableSurface3::from_frozen(&flow, control_box, settings, limits, || false)
    };
    assert!(
        build(
            settings,
            PermeableSurfaceControls3 {
                maximum_faces: count - 1,
                ..limits
            }
        )
        .is_err()
    );
    assert!(matches!(
        build(
            settings,
            PermeableSurfaceControls3 {
                maximum_workspace_bytes: limits.maximum_workspace_bytes - 1,
                ..limits
            }
        ),
        Err(CompressibleVolumeError3::Workspace { .. })
    ));
    let mut wrong = settings;
    wrong.acoustic = FluidAcousticState::try_new(
        settings.acoustic.density_kg_m3(),
        1.02 * settings.acoustic.sound_speed_m_s(),
    )
    .unwrap();
    assert!(matches!(
        build(wrong, limits),
        Err(CompressibleVolumeError3::Invalid(_))
    ));
    wrong = settings;
    wrong.grid_to_acoustic[0] = wrong.grid_to_acoustic[0].map(|x| -x);
    assert!(build(wrong, limits).is_err());
    for bounds in [
        PermeableBox3 {
            lower: [0, 1, 1],
            ..control_box
        },
        PermeableBox3 {
            upper: [4, 3, 3],
            ..control_box
        },
        PermeableBox3 {
            upper: [1; 3],
            ..control_box
        },
    ] {
        assert!(
            PreparedPermeableSurface3::from_frozen(&flow, bounds, settings, limits, || false)
                .is_err()
        );
    }
    let mut source = build(settings, limits).unwrap();
    let mut complete = build(settings, limits).unwrap();
    let original = source.source_components().to_vec();
    let original_time = source.sample_time_s();
    assert!(matches!(
        source.sample_frozen(&flow, || false),
        Err(CompressibleVolumeError3::SampleClock { .. })
    ));
    flow.step(Time::new(1e-6), StepControl::default(), || false)
        .unwrap();
    let native_before = flow.states().to_vec();
    let mut polls = 0;
    complete
        .sample_frozen(&flow, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 2 * count);
    let mut stopped = 0;
    assert!(matches!(
        source.sample_frozen(&flow, || {
            stopped += 1;
            stopped == polls
        }),
        Err(CompressibleVolumeError3::Cancelled)
    ));
    assert_eq!(source.source_components(), original);
    assert_eq!(source.sample_time_s(), original_time);
    assert_eq!(flow.states(), native_before);
    source.sample_frozen(&flow, || false).unwrap();
    assert_eq!(source.source_components(), complete.source_components());
    assert_eq!(source.sample_time_s(), complete.sample_time_s());
    let alien = uniform(frozen(0), grid);
    assert!(matches!(
        source.sample_frozen(&alien, || false),
        Err(CompressibleVolumeError3::SourceMismatch)
    ));
    let changed = uniform(
        flow.gas().clone(),
        CartesianGrid::new(5, 4, 4, 0.02, 0.03, 0.04).unwrap(),
    );
    assert!(matches!(
        source.sample_frozen(&changed, || false),
        Err(CompressibleVolumeError3::GeometryMismatch)
    ));
}

#[test]
fn native_accepted_history_reaches_coherent_welch_pressure_and_vector_velocity() {
    let gas = frozen(6);
    let settings = config(&gas);
    let grid = CartesianGrid::new(4, 4, 4, 0.02, 0.03, 0.04).unwrap();
    let states = (0..64)
        .map(|cell| {
            let x = (cell % 4) as f64;
            let z = (cell / 16) as f64;
            State::from_primitive(
                &gas,
                1.7 + 0.01 * x,
                [0.3 * x, 0.7 - z * 0.1, 0.2],
                settings.reference_pressure_pa + 400.0 * (x - z),
            )
            .unwrap()
        })
        .collect();
    let mut flow = Euler3d::new(gas, grid, states, PERIODIC, || false).unwrap();
    let control_box = PermeableBox3 {
        lower: [1; 3],
        upper: [3; 3],
    };
    let count = control_box.face_count(&grid).unwrap();
    let mut source = PreparedPermeableSurface3::from_frozen(
        &flow,
        control_box,
        settings,
        controls(count, 0),
        || false,
    )
    .unwrap();
    use crate::curle3d::array::recorded::{
        CompactSourceRecordingBudget, CompactSourceWelchStream, SourceSampleClock,
    };
    let samples = 8;
    let dt = 1e-6;
    let clock = SourceSampleClock::new(source.sample_time_s(), dt, 1e-18).unwrap();
    let recording_config = WelchConfig {
        segment_samples: samples,
        hop_samples: samples,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    };
    let recording_budget = CompactSourceRecordingBudget {
        maximum_sources: count,
        welch: WelchBudget {
            maximum_channels: 4 * count,
            maximum_samples_per_channel: samples,
            maximum_transforms: 4 * count,
            maximum_workspace_bytes: 1 << 20,
        },
    };
    let mut recording = source
        .begin_recording(clock, samples, recording_config, recording_budget, || false)
        .unwrap();
    assert_eq!(recording.accepted_sample_count(), 1);
    assert!(source.record_sample(&mut recording, || false).is_err());
    assert_eq!(recording.accepted_sample_count(), 1);
    let wrong_acoustic = FluidAcousticState::try_new(
        1.01 * settings.acoustic.density_kg_m3(),
        settings.acoustic.sound_speed_m_s(),
    ).unwrap();
    let mut alien = CompactSourceWelchStream::new(
        wrong_acoustic, source.source_positions_m(), clock, samples,
        recording_config, recording_budget, || false,
    ).unwrap();
    assert!(source.record_sample(&mut alien, || false).is_err());
    assert_eq!(alien.accepted_sample_count(), 0);
    let mut history = Vec::new();
    for n in 0..samples {
        if n != 0 {
            flow.step(Time::new(dt), StepControl::default(), || false)
                .unwrap();
            source.sample_frozen(&flow, || false).unwrap();
            source.record_sample(&mut recording, || false).unwrap();
        }
        assert_eq!(source.sample_time_s(), flow.time().value());
        assert!((source.sample_time_s() - n as f64 * dt).abs() < 1e-18);
        history.push(source.source_components().to_vec());
    }
    let recorded = recording.finish().unwrap();
    // This checks the unaltered accepted clock and coherent transform, not
    // a spatial/bandwidth qualification of this deliberately short record.
    let measured = WelchCrossSpectrum::estimate_with(
        4 * count,
        samples,
        |channel, sample| history[sample][channel / 4][channel % 4],
        dt,
        recording_config,
        recording_budget.welch,
        || false,
    )
    .unwrap();
    for bin in 0..=samples / 2 {
        assert_eq!(
            recorded.spectrum().factors(bin).unwrap().collect::<Vec<_>>(),
            measured.factors(bin).unwrap().collect::<Vec<_>>(),
        );
    }
    let observer = recorded.prepare_observer([1.1, -0.7, 0.9], 1, || false).unwrap();
    let actual = recorded
        .evaluate_bin(
            &observer,
            1,
            CoherenceBudget {
                max_modes: 1,
                max_source_mode_pairs: count,
            },
            || false,
        )
        .unwrap();
    let mut peaks = vec![CompactSourcePhasor::default(); count];
    for (n, sample) in history.iter().enumerate() {
        let phase = core::f64::consts::TAU * n as f64 / samples as f64;
        let weight = C64::new(phase.cos(), phase.sin()).scale(2.0 / samples as f64);
        for (peak, row) in peaks.iter_mut().zip(sample) {
            peak.volume_flow_m3_s = peak.volume_flow_m3_s + weight.scale(row[0]);
            for a in 0..3 {
                peak.force_on_fluid_n[a] = peak.force_on_fluid_n[a] + weight.scale(row[a + 1]);
            }
        }
    }
    let direct = observer.evaluate(&peaks).unwrap();
    let mean_square = 0.5 * direct.pressure_pa.abs().powi(2);
    assert!(mean_square > 0.0);
    assert!((actual.mean_square_pressure_pa2 - mean_square).abs() < 2e-8 * mean_square);
    for a in 0..3 {
        let u = direct.particle_velocity_m_s[a];
        let expected = 0.5 * u.abs().powi(2);
        assert!((actual.velocity_covariance_m2_s2[a][a].re - expected).abs() < 2e-8 * expected);
        assert!(
            (actual.active_intensity_w_m2[a] - direct.active_intensity_w_m2[a]).abs()
                < 2e-8 * (direct.active_intensity_w_m2[a].abs() + 1e-30)
        );
    }
}

mod molecular;
