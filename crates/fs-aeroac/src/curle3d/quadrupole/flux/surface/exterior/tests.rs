use super::*;
use crate::curle3d::CoherenceBudget;
use crate::curle3d::array::recorded::SourceSampleClock;
use crate::curle3d::array::recorded::joint::JointSourceRecordingBudget;
use crate::curle3d::quadrupole::flux::surface::{
    PermeableBox3, PermeableSurfaceControls3, PreparedPermeableSurface3,
};
use crate::curle3d::quadrupole::flux::{
    CompressibleVolumeConfig3, CompressibleVolumeError3, PreparedCompressibleVolume3,
};
use fs_fft::welch::{WelchBudget, WelchConfig, WelchCrossSpectrum, WelchDetrend, WelchWindow};
use fs_flux::compressible::spatial::mixture::{MixtureEuler3d, MixtureGas, state_into};
use fs_flux::compressible::spatial::{
    AxisBoundaries, Boundaries3d, CartesianGrid, Euler3d, FrozenGas, State, StepControl,
};
use fs_gas::dynamics::CaloricDomain;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_qty::{Temperature, Time};
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, cantera_h2o2};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis};

const DT: f64 = 1e-6;
const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const ROTATION: [[f64; 3]; 3] = [[0.48, -0.6, 0.64], [0.36, 0.8, 0.48], [-0.8, 0.0, 0.6]];
const BOX: PermeableBox3 = PermeableBox3 {
    lower: [1; 3],
    upper: [3; 3],
};
// Native x-fast indices: (0,1,1), (3,1,1), (1,1,3). Each lies immediately
// OUTSIDE a different face. Interpolation of its state onto the control
// surface does not remove that cell's exterior volume contribution.
const EXTERIOR: [usize; 3] = [20, 23, 53];

fn gas(species: usize) -> FrozenGas {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![mechanism.thermodynamics().models()[species].clone()],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap()
}

fn settings(gas: &FrozenGas) -> CompressibleVolumeConfig3 {
    let rho = 1.7;
    let ambient = State::from_primitive(gas, rho, [0.0; 3], rho * gas.gas_constant() * 500.0)
        .unwrap()
        .primitive(gas, || false)
        .unwrap();
    CompressibleVolumeConfig3 {
        acoustic: FluidAcousticState::try_new(rho, ambient.sound_speed).unwrap(),
        reference_pressure_pa: ambient.pressure,
        origin_m: [0.2, -0.3, 0.1],
        grid_to_acoustic: ROTATION,
    }
}

fn grid() -> CartesianGrid {
    CartesianGrid::new(4, 4, 4, 0.02, 0.03, 0.04).unwrap()
}

fn frozen_flow(gas: FrozenGas, grid: CartesianGrid) -> Euler3d {
    let reference = settings(&gas);
    let [nx, ny, nz] = grid.shape();
    let states = (0..nx * ny * nz)
        .map(|cell| {
            let x = (cell % nx) as f64;
            let y = ((cell / nx) % ny) as f64;
            let z = (cell / (nx * ny)) as f64;
            State::from_primitive(
                &gas,
                1.7 + 0.01 * x,
                [0.3 * x, 0.7 - 0.1 * z, 0.2 + 0.05 * y],
                reference.reference_pressure_pa + 400.0 * (x - z),
            )
            .unwrap()
        })
        .collect();
    Euler3d::new(gas, grid, states, PERIODIC, || false).unwrap()
}

fn config() -> WelchConfig {
    WelchConfig {
        segment_samples: 4,
        hop_samples: 2,
        window: WelchWindow::Rectangular,
        detrend: WelchDetrend::None,
    }
}

fn limits(species: usize, samples: usize) -> ExteriorSurfaceRecordingBudget {
    let faces = BOX.face_count(&grid()).unwrap();
    let channels = 4 * faces + 6 * EXTERIOR.len();
    let segments = 1 + (samples - config().segment_samples) / config().hop_samples;
    ExteriorSurfaceRecordingBudget {
        surface: PermeableSurfaceControls3 {
            maximum_faces: faces,
            maximum_workspace_bytes: PreparedPermeableSurface3::required_workspace_bytes(
                faces, species,
            )
            .unwrap(),
            reference_relative_tolerance: 2e-9,
        },
        maximum_volume_workspace_bytes: PreparedCompressibleVolume3::required_workspace_bytes(
            EXTERIOR.len(),
        )
        .unwrap(),
        recording: JointSourceRecordingBudget {
            maximum_compact_sources: faces,
            maximum_quadrupole_sources: EXTERIOR.len(),
            welch: WelchBudget {
                maximum_channels: channels,
                maximum_samples_per_channel: samples,
                maximum_transforms: channels * segments,
                maximum_workspace_bytes: 1 << 20,
            },
        },
    }
}

fn recorder(flow: &Euler3d, samples: usize) -> ExteriorSurfaceWelchStream {
    ExteriorSurfaceWelchStream::from_frozen(
        flow,
        BOX,
        &EXTERIOR,
        settings(flow.gas()),
        SourceSampleClock::new(flow.time().value(), DT, 1e-17).unwrap(),
        samples,
        config(),
        limits(0, samples),
        || false,
    )
    .unwrap()
}

fn frame(surface: &PreparedPermeableSurface3, volume: &PreparedCompressibleVolume3) -> Vec<f64> {
    surface
        .source_components()
        .iter()
        .flatten()
        .chain(volume.stress_moments_n_m().iter().flatten())
        .copied()
        .collect()
}

fn assert_factors(actual: &WelchCrossSpectrum, expected: &WelchCrossSpectrum) {
    assert_eq!(actual.channel_count(), expected.channel_count());
    assert_eq!(actual.segment_count(), expected.segment_count());
    assert_eq!(actual.bin_count(), expected.bin_count());
    assert_eq!(actual.timestep_s(), expected.timestep_s());
    for bin in 0..actual.bin_count() {
        for (a, b) in actual
            .factors(bin)
            .unwrap()
            .zip(expected.factors(bin).unwrap())
        {
            // Both owners use the same original segment transform, including
            // overlap, final-record normalization, channel order and phase.
            assert_eq!(a, b);
        }
    }
}

fn batch(history: &[Vec<f64>], species: usize) -> WelchCrossSpectrum {
    WelchCrossSpectrum::estimate_with(
        history[0].len(),
        history.len(),
        |channel, sample| history[sample][channel],
        DT,
        config(),
        limits(species, history.len()).recording.welch,
        || false,
    )
    .unwrap()
}

fn point(settings: CompressibleVolumeConfig3, local: [f64; 3]) -> [f64; 3] {
    std::array::from_fn(|a| {
        settings.origin_m[a]
            + (0..3)
                .map(|b| settings.grid_to_acoustic[a][b] * local[b])
                .sum::<f64>()
    })
}

#[test]
fn accepted_frozen_history_reaches_one_joint_transform_and_complete_vector_field() {
    let mut flow = frozen_flow(gas(6), grid());
    // The physical epoch is the already accepted native time, not an invented
    // zero. This short fixture checks ownership/coherence, not resolved SPL.
    flow.step(Time::new(DT), StepControl::default(), || false)
        .unwrap();
    let settings = settings(flow.gas());
    let mut stream = recorder(&flow, 8);
    let clock = stream.clock();
    assert_eq!(clock.epoch_s(), flow.time().value());
    assert_eq!(stream.accepted_sample_count(), 1);
    let budget = limits(0, 8);
    let mut surface =
        PreparedPermeableSurface3::from_frozen(&flow, BOX, settings, budget.surface, || false)
            .unwrap();
    let mut volume = PreparedCompressibleVolume3::from_frozen(
        &flow,
        &EXTERIOR,
        settings,
        budget.maximum_volume_workspace_bytes,
        || false,
    )
    .unwrap();
    let mut history = Vec::new();
    for n in 0..8 {
        if n != 0 {
            flow.step(Time::new(DT), StepControl::default(), || false)
                .unwrap();
            surface.sample_frozen(&flow, || false).unwrap();
            volume.sample_frozen(&flow, || false).unwrap();
            let original = flow.states().to_vec();
            stream.sample_frozen(&flow, || false).unwrap();
            assert_eq!(flow.states(), original);
        }
        assert_eq!(stream.accepted_sample_count(), n + 1);
        assert_eq!(stream.remaining_sample_count(), 7 - n);
        assert_eq!(stream.surface().sample_time_s(), flow.time().value());
        assert_eq!(stream.volume().sample_time_s(), flow.time().value());
        assert_eq!(
            stream.surface().source_components(),
            surface.source_components()
        );
        assert_eq!(
            stream.volume().stress_moments_n_m(),
            volume.stress_moments_n_m()
        );
        history.push(frame(&surface, &volume));
    }
    assert_eq!(stream.completed_segment_count(), 3);
    assert!(history.windows(2).any(|samples| samples[0] != samples[1]));
    let measured = batch(&history, 0);
    let recorded = stream.finish().unwrap();
    assert_eq!(recorded.clock(), clock);
    assert_factors(recorded.spectrum(), &measured);
    // A nonzero surface/volume cross spectrum proves that treating the two
    // extracted mechanisms as independent would discard measured information.
    let volume_start = 4 * surface.faces().len();
    assert!((0..volume_start).any(|s| {
        (volume_start..measured.channel_count())
            .any(|v| {
                let cross = measured.cross_spectrum(1, s, v).unwrap();
                cross.re.hypot(cross.im) > 0.0
            })
    }));
    let (compact, quadrupole) = recorded
        .prepare_observers([1.1, -0.7, 0.9], 1, || false)
        .unwrap();
    let coherence = CoherenceBudget {
        max_modes: measured.segment_count(),
        max_source_mode_pairs: (surface.faces().len() + EXTERIOR.len()) * measured.segment_count(),
    };
    let observed = recorded
        .evaluate_bin(&compact, &quadrupole, 1, coherence, || false)
        .unwrap();
    let expected = compact
        .evaluate_welch_bin_with_quadrupoles(&quadrupole, &measured, 1, coherence, || false)
        .unwrap();
    assert!(observed.mean_square_pressure_pa2 > 0.0);
    assert_eq!(
        observed.mean_square_pressure_pa2,
        expected.mean_square_pressure_pa2
    );
    assert_eq!(
        observed.velocity_covariance_m2_s2,
        expected.velocity_covariance_m2_s2
    );
    assert_eq!(
        observed.active_intensity_w_m2,
        expected.active_intensity_w_m2
    );
    assert_eq!(
        observed.reactive_intensity_w_m2,
        expected.reactive_intensity_w_m2
    );
    // Rotation and origin belong to the original geometry. The center and
    // every exact face center are excluded, even after that frame transform.
    for local in [
        [0.04, 0.06, 0.08],
        [0.02, 0.06, 0.08],
        [0.06, 0.06, 0.08],
        [0.04, 0.03, 0.08],
        [0.04, 0.09, 0.08],
        [0.04, 0.06, 0.04],
        [0.04, 0.06, 0.12],
    ] {
        assert!(
            recorded
                .prepare_observers(point(settings, local), 1, || false)
                .is_err()
        );
    }
    // A genuinely exterior point immediately beyond a face is still usable.
    assert!(
        recorded
            .prepare_observers(point(settings, [0.019, 0.06, 0.08]), 1, || false)
            .is_ok()
    );
}

#[test]
fn enclosed_cells_refuse_but_original_face_adjacent_exterior_cells_remain() {
    let flow = frozen_flow(gas(6), grid());
    let settings = settings(flow.gas());
    let build = |cells: &[usize]| {
        ExteriorSurfaceWelchStream::from_frozen(
            &flow,
            BOX,
            cells,
            settings,
            SourceSampleClock::new(0.0, DT, 1e-17).unwrap(),
            4,
            config(),
            limits(0, 4),
            || false,
        )
    };
    // Both the low and high corner cells are inside the half-open box. A
    // valid exterior entry beside one must not make a mixed selection pass.
    for cells in [&[21][..], &[42][..], &[20, 21, 23][..], &[20, 23, 42][..]] {
        assert!(build(cells).is_err());
    }
    let stream = build(&EXTERIOR).unwrap();
    assert_eq!(stream.volume().cell_indices(), EXTERIOR);
    for &cell in &EXTERIOR {
        assert!(
            stream
                .surface()
                .faces()
                .iter()
                .any(|face| face.adjacent_cells.contains(&cell))
        );
    }
    assert!(
        stream
            .volume()
            .stress_moments_n_m()
            .iter()
            .any(|row| row.iter().any(|v| *v != 0.0))
    );
    assert_eq!(stream.surface().config(), stream.volume().config());
    assert_eq!(stream.surface().frozen_source(), Some(flow.gas()));
    assert_eq!(stream.volume().frozen_source(), Some(flow.gas()));
    assert!(build(&[20, 20]).is_err());
    assert!(build(&[23, 20]).is_err());
}

#[test]
fn failed_identity_clock_and_late_transform_leave_both_sources_retryable() {
    let mut flow = frozen_flow(gas(6), grid());
    let mut stream = recorder(&flow, 4);
    let original = frame(stream.surface(), stream.volume());
    let source_time = stream.surface().sample_time_s();
    let mut alien = frozen_flow(gas(0), grid());
    alien
        .step(Time::new(DT), StepControl::default(), || false)
        .unwrap();
    assert!(matches!(
        stream.sample_frozen(&alien, || false),
        Err(ExteriorSurfaceRecordingError::Source(
            CompressibleVolumeError3::SourceMismatch
        ))
    ));
    let mut changed = frozen_flow(
        flow.gas().clone(),
        CartesianGrid::new(5, 4, 4, 0.02, 0.03, 0.04).unwrap(),
    );
    changed
        .step(Time::new(DT), StepControl::default(), || false)
        .unwrap();
    assert!(matches!(
        stream.sample_frozen(&changed, || false),
        Err(ExteriorSurfaceRecordingError::Source(
            CompressibleVolumeError3::GeometryMismatch
        ))
    ));
    let mut off_clock = frozen_flow(flow.gas().clone(), grid());
    off_clock
        .step(Time::new(1.5 * DT), StepControl::default(), || false)
        .unwrap();
    // The native time is later but off the record cadence. Early clock
    // admission must leave both original accepted source buffers unchanged.
    assert!(stream.sample_frozen(&off_clock, || false).is_err());
    assert_eq!(stream.accepted_sample_count(), 1);
    assert_eq!(frame(stream.surface(), stream.volume()), original);
    assert_eq!(stream.surface().sample_time_s(), source_time);
    assert_eq!(stream.volume().sample_time_s(), source_time);
    let mut reference = recorder(&flow, 4);
    for _ in 1..3 {
        flow.step(Time::new(DT), StepControl::default(), || false)
            .unwrap();
        stream.sample_frozen(&flow, || false).unwrap();
        reference.sample_frozen(&flow, || false).unwrap();
    }
    let previous = frame(stream.surface(), stream.volume());
    let previous_time = stream.surface().sample_time_s();
    let budget = limits(0, 4);
    let source_config = settings(flow.gas());
    let mut surface_probe =
        PreparedPermeableSurface3::from_frozen(&flow, BOX, source_config, budget.surface, || false)
            .unwrap();
    let mut volume_probe = PreparedCompressibleVolume3::from_frozen(
        &flow,
        &EXTERIOR,
        source_config,
        budget.maximum_volume_workspace_bytes,
        || false,
    )
    .unwrap();
    flow.step(Time::new(DT), StepControl::default(), || false)
        .unwrap();
    let native_before = flow.states().to_vec();
    let mut surface_polls = 0;
    surface_probe
        .sample_frozen(&flow, || {
            surface_polls += 1;
            false
        })
        .unwrap();
    let mut volume_polls = 0;
    volume_probe
        .sample_frozen(&flow, || {
            volume_polls += 1;
            false
        })
        .unwrap();
    assert!(volume_polls > 2);
    let stop_in_volume = 1 + surface_polls + volume_polls / 2;
    let mut extraction_polls = 0;
    // One joint clock poll precedes the original surface/volume samplers.
    // The surface has committed here but volume trial data are still partial:
    // restore only the former, never publish the latter's dirty trial.
    assert!(matches!(
        stream.sample_frozen(&flow, || {
            extraction_polls += 1;
            extraction_polls == stop_in_volume
        }),
        Err(ExteriorSurfaceRecordingError::Source(
            CompressibleVolumeError3::Cancelled
        ))
    ));
    assert_eq!(extraction_polls, stop_in_volume);
    assert_eq!(stream.accepted_sample_count(), 3);
    assert_eq!(stream.completed_segment_count(), 0);
    assert_eq!(frame(stream.surface(), stream.volume()), previous);
    assert_eq!(stream.surface().sample_time_s(), previous_time);
    assert_eq!(stream.volume().sample_time_s(), previous_time);
    let mut polls = 0;
    reference
        .sample_frozen(&flow, || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(polls > 3);
    let mut stopped = 0;
    // The last cancellation poll is inside the original transform commit.
    // Both native extractions have succeeded, but no stream state is accepted.
    assert!(
        stream
            .sample_frozen(&flow, || {
                stopped += 1;
                stopped == polls
            })
            .is_err()
    );
    assert_eq!(stopped, polls);
    assert_eq!(stream.accepted_sample_count(), 3);
    assert_eq!(stream.completed_segment_count(), 0);
    assert_eq!(frame(stream.surface(), stream.volume()), previous);
    assert_eq!(stream.surface().sample_time_s(), previous_time);
    assert_eq!(stream.volume().sample_time_s(), previous_time);
    assert_eq!(flow.states(), native_before);
    stream.sample_frozen(&flow, || false).unwrap();
    assert_eq!(
        stream.surface().source_components(),
        reference.surface().source_components()
    );
    assert_eq!(
        stream.volume().stress_moments_n_m(),
        reference.volume().stress_moments_n_m()
    );
    assert_eq!(stream.accepted_sample_count(), 4);
    assert_eq!(stream.completed_segment_count(), 1);
    let actual = stream.finish().unwrap();
    let expected = reference.finish().unwrap();
    assert_factors(actual.spectrum(), expected.spectrum());
}

#[test]
fn native_mixture_species_and_ambient_reference_share_one_accepted_history() {
    let settings = settings(&gas(6));
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let mixture = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap();
    let species = mixture.species_count();
    let mut ambient = vec![0.0; species];
    ambient[6] = settings.acoustic.density_kg_m3();
    let mut values = vec![0.0; 64 * (species + 4)];
    for (cell, state) in values.chunks_exact_mut(species + 4).enumerate() {
        let x = (cell % 4) as f64;
        let z = (cell / 16) as f64;
        let mut partial = vec![0.0; species];
        partial[0] = 0.1 + 0.01 * x; // Original AR and N2 calorics, not gamma.
        partial[6] = 1.6 - 0.01 * x;
        state_into(
            &mixture,
            &partial,
            [0.2 + 0.1 * x, -0.3, 0.5 - 0.05 * z],
            settings.reference_pressure_pa + 100.0 * (x - z),
            state,
        )
        .unwrap();
    }
    let workspace = MixtureEuler3d::required_workspace_bytes(&mixture, &grid()).unwrap();
    let mut flow =
        MixtureEuler3d::new(mixture, grid(), values, PERIODIC, workspace, || false).unwrap();
    let clock = SourceSampleClock::new(flow.time().value(), DT, 1e-17).unwrap();
    let budget = limits(species, 4);
    let mut wrong_ambient = ambient.clone();
    wrong_ambient[6] *= 1.01;
    assert!(
        ExteriorSurfaceWelchStream::from_mixture(
            &flow,
            &wrong_ambient,
            BOX,
            &EXTERIOR,
            settings,
            clock,
            4,
            config(),
            budget,
            || false
        )
        .is_err()
    );
    let mut stream = ExteriorSurfaceWelchStream::from_mixture(
        &flow,
        &ambient,
        BOX,
        &EXTERIOR,
        settings,
        clock,
        4,
        config(),
        budget,
        || false,
    )
    .unwrap();
    assert_eq!(
        stream.surface().reference_partial_densities_kg_m3(),
        ambient
    );
    assert_eq!(stream.surface().mixture_source(), Some(flow.gas()));
    assert_eq!(stream.volume().mixture_source(), Some(flow.gas()));
    assert!(stream.surface().frozen_source().is_none());
    assert!(stream.volume().frozen_source().is_none());
    let mut history = Vec::new();
    for n in 0..4 {
        if n != 0 {
            flow.step(Time::new(DT), StepControl::default(), || false)
                .unwrap();
            stream.sample_mixture(&flow, || false).unwrap();
        }
        let source_before = flow.states().to_vec();
        let surface = PreparedPermeableSurface3::from_mixture(
            &flow,
            &ambient,
            BOX,
            settings,
            budget.surface,
            || false,
        )
        .unwrap();
        let volume = PreparedCompressibleVolume3::from_mixture(
            &flow,
            &EXTERIOR,
            settings,
            budget.maximum_volume_workspace_bytes,
            || false,
        )
        .unwrap();
        assert_eq!(
            stream.surface().source_components(),
            surface.source_components()
        );
        assert_eq!(
            stream.volume().stress_moments_n_m(),
            volume.stress_moments_n_m()
        );
        assert_eq!(stream.surface().sample_time_s(), flow.time().value());
        assert_eq!(stream.volume().sample_time_s(), flow.time().value());
        assert_eq!(flow.states(), source_before);
        history.push(frame(&surface, &volume));
    }
    assert!(history.windows(2).any(|samples| samples[0] != samples[1]));
    let expected = batch(&history, species);
    assert_factors(stream.finish().unwrap().spectrum(), &expected);
}

mod molecular;
