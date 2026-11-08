use super::*;
use fs_flux::compressible::spatial::mixture::{MixtureGas, state_into};
use fs_flux::compressible::spatial::{AxisBoundaries, Boundaries3d, FrozenGas, State, StepControl};
use fs_gas::dynamics::CaloricDomain;
use fs_material::fluid_acoustics::FluidAcousticState;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation,
};
use fs_qty::{Pressure, Temperature, Time};
use fs_thermochem::kinetics::hydrogen::{HydrogenCalorics, cantera_h2o2};
use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
use fs_thermochem::{Composition, CompositionBasis, SpeciesId};

const DT: f64 = 5e-7;
const BOX: PermeableBox3 = PermeableBox3 {
    lower: [1; 3],
    upper: [4; 3],
};
const EXTERIOR: [usize; 1] = [60]; // (0,2,2), immediately outside the lower x face.
const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};

fn grid() -> CartesianGrid {
    CartesianGrid::new(5, 5, 5, 0.005, 0.005, 0.005).unwrap()
}
fn gas() -> FrozenGas {
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![mechanism.thermodynamics().models()[6].clone()],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap()
}

fn fixture() -> (Euler3d, CompressibleVolumeConfig3, RelaxationState) {
    let gas = gas();
    let rho = 1.7;
    let ambient =
        State::from_primitive(&gas, rho, [0.0; 3], rho * gas.gas_constant() * 500.0).unwrap();
    let primitive = ambient.primitive(&gas, || false).unwrap();
    let source = CompressibleVolumeConfig3 {
        acoustic: FluidAcousticState::try_new(rho, primitive.sound_speed).unwrap(),
        reference_pressure_pa: primitive.pressure,
        origin_m: [0.0; 3],
        grid_to_acoustic: [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    };
    // Keep the ORIGINAL frozen source's caloric state, not a second gamma-law
    // surrogate. Positive air transport is unused by this molecular-only test.
    let properties = gas
        .source()
        .evaluate_properties(
            Temperature::new(primitive.temperature),
            Pressure::new(primitive.pressure),
        )
        .unwrap();
    let cp = properties.mass_specific_cp().value();
    let cv = properties.mass_specific_cv().value();
    let mut state = GasState::try_new(
        &GasSpec::dry_air_ussa1976(),
        primitive.temperature,
        primitive.pressure,
    )
    .unwrap();
    state.density = rho;
    state.sound_speed = primitive.sound_speed;
    state.gamma = cp / cv;
    state.specific_gas_constant = gas.gas_constant();
    state.specific_heat_cp = cp;
    state.prandtl = state.dynamic_viscosity * cp / state.thermal_conductivity;
    state.characteristic_impedance = rho * primitive.sound_speed;
    let species = SpeciesId::new("N2").unwrap();
    let molecular = RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: species.clone(),
        processes: vec![RelaxationProcess {
            reservoir: InternalReservoir::Rotational {
                degrees_of_freedom: 2,
            },
            rate: RelaxationRateLaw::PartnerIndependent {
                p_tau_ref_pa_s: primitive.pressure * 1e-5,
                t_ref_k: primitive.temperature,
                exponent: 0.0,
                t_min_k: 300.0,
                t_max_k: 900.0,
            },
        }],
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(state, vec![(species, 1.0)]).unwrap())
    .unwrap();
    let mut states = vec![ambient; 125];
    // This pulse is strictly inside the closed surface and initially absent
    // from ALL face reconstructions and the selected exterior-volume cell.
    states[62] = State::from_primitive(&gas, rho, [0.0; 3], primitive.pressure * 1.01).unwrap();
    let flow = Euler3d::new(gas, grid(), states, PERIODIC, || false).unwrap();
    (flow, source, molecular)
}

fn limits(samples: usize, species: usize) -> ExteriorSurfaceHistoryBudget {
    let faces = BOX.face_count(&grid()).unwrap();
    ExteriorSurfaceHistoryBudget {
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
        recording: AcousticSourceHistoryBudget {
            maximum_compact_sources: faces,
            maximum_quadrupole_sources: EXTERIOR.len(),
            maximum_samples_per_source: samples,
            maximum_component_samples: (4 * faces + 6 * EXTERIOR.len()) * samples,
            maximum_workspace_bytes: 2 << 20,
        },
    }
}
fn clock() -> SourceSampleClock {
    SourceSampleClock::new(0.0, DT, 1e-17).unwrap()
}
fn endpoints() -> (Vec<CompactHistoryEnds>, Vec<QuadrupoleHistoryEnds>) {
    // Explicit clamped model, not an inferred derivative or a continuation.
    (
        vec![
            CompactHistoryEnds {
                initial_volume_flow_rate_m3_s2: 0.0,
                final_volume_flow_rate_m3_s2: 0.0,
                initial_force_rate_n_s: [0.0; 3],
                final_force_rate_n_s: [0.0; 3],
                initial_force_impulse_n_s: [0.0; 3],
            };
            BOX.face_count(&grid()).unwrap()
        ],
        vec![
            QuadrupoleHistoryEnds {
                initial_rate_n_m_s: [0.0; 6],
                final_rate_n_m_s: [0.0; 6],
                initial_integral_n_m_s: [0.0; 6],
            };
            EXTERIOR.len()
        ],
    )
}

#[test]
fn accepted_native_flow_history_drives_complete_causal_molecular_fields() {
    let (mut flow, source, molecular) = fixture();
    let samples = 64;
    let budget = limits(samples, 0);
    let mut recording = ExteriorSurfaceHistoryRecorder::from_frozen(
        &flow,
        BOX,
        &EXTERIOR,
        source,
        clock(),
        samples,
        budget,
        || false,
    )
    .unwrap();
    assert!(
        recording
            .recording()
            .volume_force_samples()
            .iter()
            .flatten()
            .all(|&v| v == 0.0)
    );
    assert!(
        recording
            .recording()
            .stress_samples_n_m()
            .iter()
            .flatten()
            .all(|&v| v == 0.0)
    );
    for n in 1..samples {
        flow.step(Time::new(DT), StepControl::default(), || false)
            .unwrap();
        let before = flow.states().to_vec();
        recording.sample_frozen(&flow, || false).unwrap();
        assert_eq!(flow.states(), before);
        assert_eq!(recording.accepted_sample_count(), n + 1);
        let surface =
            PreparedPermeableSurface3::from_frozen(&flow, BOX, source, budget.surface, || false)
                .unwrap();
        let volume = PreparedCompressibleVolume3::from_frozen(
            &flow,
            &EXTERIOR,
            source,
            budget.maximum_volume_workspace_bytes,
            || false,
        )
        .unwrap();
        let faces = surface.source_components().len();
        assert_eq!(
            &recording.recording().volume_force_samples()[n * faces..(n + 1) * faces],
            surface.source_components()
        );
        assert_eq!(
            &recording.recording().stress_samples_n_m()[n..n + 1],
            volume.stress_moments_n_m()
        );
        assert_eq!(
            recording.recording().sample_times_s()[n],
            flow.time().value()
        );
    }
    assert!(
        recording
            .recording()
            .volume_force_samples()
            .iter()
            .flatten()
            .any(|&v| v != 0.0)
    );
    let (ce, se) = endpoints();
    let history = recording
        .prepare_clamped_cubic(&ce, &se, 100_000.0, || false)
        .unwrap();
    let config = MolecularTransientConfig {
        timestep_s: DT,
        samples,
        transform_len: 4 * samples,
        contour_exponent: 24.0,
    };
    let receiver = [0.0045, 0.0125, 0.0125];
    let observer = history
        .prepare_molecular_observer(
            molecular.clone(),
            receiver,
            config,
            MolecularTransientBudget::default(),
            || false,
        )
        .unwrap();
    let fields = history
        .observe_initially_resting(&observer, || false)
        .unwrap();
    let front = (0..observer.paths().path_count())
        .map(|i| observer.paths().path_delay_s(i).unwrap())
        .fold(f64::INFINITY, f64::min);
    assert!(front > DT);
    for (n, field) in fields.fields().iter().enumerate() {
        assert!(field.pressure_pa.is_finite());
        assert!(field.particle_velocity_m_s.iter().all(|v| v.is_finite()));
        if n as f64 * DT < front {
            assert_eq!(field.pressure_pa, 0.0);
            assert_eq!(field.particle_velocity_m_s, [0.0; 3]);
        }
    }
    assert!(fields.fields().iter().any(|f| f.pressure_pa != 0.0));
    assert!(
        fields
            .fields()
            .iter()
            .any(|f| f.particle_velocity_m_s.iter().any(|&u| u != 0.0))
    );
    for receiver in [[0.005, 0.0125, 0.0125], [0.0125; 3]] {
        assert!(
            history
                .prepare_molecular_observer(
                    molecular.clone(),
                    receiver,
                    config,
                    MolecularTransientBudget::default(),
                    || false
                )
                .is_err()
        );
    }
    // A caller-created native observer cannot bypass the retained FW-H domain.
    let interior = PreparedMolecularSourceObserver::new(
        molecular,
        [0.0125; 3],
        history.source_positions_m(),
        config,
        MolecularTransientBudget::default(),
        || false,
    )
    .unwrap();
    assert!(
        history
            .observe_initially_resting(&interior, || false)
            .is_err()
    );
}

#[test]
fn extraction_and_last_raw_commit_cancellation_restore_the_same_native_snapshot() {
    let (mut flow, source, _) = fixture();
    let budget = limits(2, 0);
    let make = || {
        ExteriorSurfaceHistoryRecorder::from_frozen(
            &flow,
            BOX,
            &EXTERIOR,
            source,
            clock(),
            2,
            budget,
            || false,
        )
        .unwrap()
    };
    let mut actual = make();
    let mut reference = make();
    let mut surface_probe =
        PreparedPermeableSurface3::from_frozen(&flow, BOX, source, budget.surface, || false)
            .unwrap();
    let mut volume_probe = PreparedCompressibleVolume3::from_frozen(
        &flow,
        &EXTERIOR,
        source,
        budget.maximum_volume_workspace_bytes,
        || false,
    )
    .unwrap();
    let before_q = actual.surface().source_components().to_vec();
    let before_t = actual.volume().stress_moments_n_m().to_vec();
    flow.step(Time::new(DT), StepControl::default(), || false)
        .unwrap();
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
    let stop = 1 + surface_polls + volume_polls / 2;
    let mut polls = 0;
    assert!(matches!(
        actual.sample_frozen(&flow, || {
            polls += 1;
            polls == stop
        }),
        Err(ExteriorSurfaceRecordingError::Source(
            CompressibleVolumeError3::Cancelled
        ))
    ));
    assert_eq!(actual.accepted_sample_count(), 1);
    assert_eq!(actual.surface().source_components(), before_q);
    assert_eq!(actual.volume().stress_moments_n_m(), before_t);
    assert_eq!(actual.surface().sample_time_s(), 0.0);
    assert_eq!(actual.volume().sample_time_s(), 0.0);
    let mut final_polls = 0;
    reference
        .sample_frozen(&flow, || {
            final_polls += 1;
            false
        })
        .unwrap();
    polls = 0;
    assert!(matches!(
        actual.sample_frozen(&flow, || {
            polls += 1;
            polls == final_polls
        }),
        Err(ExteriorSurfaceRecordingError::Acoustic(
            AeroacError::Cancelled { .. }
        ))
    ));
    assert_eq!(actual.accepted_sample_count(), 1);
    assert_eq!(actual.recording().volume_force_samples(), before_q);
    assert_eq!(actual.recording().stress_samples_n_m(), before_t);
    assert_eq!(actual.surface().source_components(), before_q);
    assert_eq!(actual.volume().stress_moments_n_m(), before_t);
    actual.sample_frozen(&flow, || false).unwrap();
    assert_eq!(
        actual.recording().sample_times_s(),
        reference.recording().sample_times_s()
    );
    assert_eq!(
        actual.recording().volume_force_samples(),
        reference.recording().volume_force_samples()
    );
    assert_eq!(
        actual.recording().stress_samples_n_m(),
        reference.recording().stress_samples_n_m()
    );
    let (ce, se) = endpoints();
    assert!(
        actual
            .prepare_clamped_cubic(&ce, &se, 100_000.0, || true)
            .is_err()
    );
    assert_eq!(actual.accepted_sample_count(), 2);
    assert!(
        actual
            .prepare_clamped_cubic(&ce, &se, 100_000.0, || false)
            .is_ok()
    );
}

#[test]
fn native_mixture_history_preserves_species_reference_channels_and_exterior_cells() {
    let (_, source, _) = fixture();
    let mechanism = cantera_h2o2(HydrogenCalorics::Original).unwrap();
    let mixture = MixtureGas::new(
        mechanism.thermodynamics().clone(),
        CaloricDomain::new(Temperature::new(300.0), Temperature::new(900.0)),
        || false,
    )
    .unwrap();
    let species = mixture.species_count();
    let mut ambient = vec![0.0; species];
    ambient[6] = source.acoustic.density_kg_m3();
    let mut values = vec![0.0; 125 * (species + 4)];
    for (cell, state) in values.chunks_exact_mut(species + 4).enumerate() {
        let x = (cell % 5) as f64;
        let mut partial = ambient.clone();
        partial[6] -= 0.1;
        partial[0] = 0.1;
        state_into(
            &mixture,
            &partial,
            [0.1 * x, 0.2, -0.1],
            source.reference_pressure_pa + 100.0 * x,
            state,
        )
        .unwrap();
    }
    let workspace = MixtureEuler3d::required_workspace_bytes(&mixture, &grid()).unwrap();
    let mut flow =
        MixtureEuler3d::new(mixture, grid(), values, PERIODIC, workspace, || false).unwrap();
    let budget = limits(4, species);
    assert!(
        ExteriorSurfaceHistoryRecorder::from_mixture(
            &flow,
            &ambient,
            BOX,
            &[62],
            source,
            clock(),
            4,
            budget,
            || false
        )
        .is_err()
    );
    let mut wrong = ambient.clone();
    wrong[6] *= 1.01;
    assert!(
        ExteriorSurfaceHistoryRecorder::from_mixture(
            &flow,
            &wrong,
            BOX,
            &EXTERIOR,
            source,
            clock(),
            4,
            budget,
            || false
        )
        .is_err()
    );
    let mut recording = ExteriorSurfaceHistoryRecorder::from_mixture(
        &flow,
        &ambient,
        BOX,
        &EXTERIOR,
        source,
        clock(),
        4,
        budget,
        || false,
    )
    .unwrap();
    assert_eq!(recording.surface().mixture_source(), Some(flow.gas()));
    assert_eq!(
        recording.surface().reference_partial_densities_kg_m3(),
        ambient
    );
    assert_eq!(recording.volume().cell_indices(), EXTERIOR);
    let faces = BOX.face_count(&grid()).unwrap();
    for n in 0..4 {
        if n != 0 {
            flow.step(Time::new(DT), StepControl::default(), || false)
                .unwrap();
            recording.sample_mixture(&flow, || false).unwrap();
        }
        let surface = PreparedPermeableSurface3::from_mixture(
            &flow,
            &ambient,
            BOX,
            source,
            budget.surface,
            || false,
        )
        .unwrap();
        let volume = PreparedCompressibleVolume3::from_mixture(
            &flow,
            &EXTERIOR,
            source,
            budget.maximum_volume_workspace_bytes,
            || false,
        )
        .unwrap();
        assert_eq!(
            &recording.recording().volume_force_samples()[n * faces..(n + 1) * faces],
            surface.source_components()
        );
        assert_eq!(
            &recording.recording().stress_samples_n_m()[n..n + 1],
            volume.stress_moments_n_m()
        );
    }
    let (ce, se) = endpoints();
    let history = recording
        .prepare_clamped_cubic(&ce, &se, 100_000.0, || false)
        .unwrap();
    assert_eq!(history.sample_count(), 4);
    assert_eq!(history.compact_source_count(), faces);
    assert_eq!(history.quadrupole_source_count(), 1);
}
