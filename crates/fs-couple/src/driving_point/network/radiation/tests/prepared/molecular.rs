//! G1/G3/G4: original molecular lines reach the existing loaded BEM consumers.

use super::*;
use fs_bem::observation_transfers::PreparedRadiationOptions;
use fs_duct::ConeSlices;
use fs_duct::molecular::MolecularWallResolution;
use fs_duct::network::molecular::MolecularEdge;
use fs_material::gas_species::input::GasRecord;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxationState,
    RelaxingMedium, SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

fn regional() -> RadiatingNetwork {
    let base = two_source_model();
    let media = NetworkMedia::new(
        base.driving_point().network(),
        base.driving_point()
            .network()
            .edges
            .iter()
            .map(|edge| (edge.id.clone(), GasRecord::Legacy(base.gas())))
            .collect(),
        GasRecord::Legacy(base.gas()),
    )
    .unwrap();
    RadiatingNetwork::new_regional(
        base.driving,
        base.surface,
        base.patches,
        media,
        base.receivers,
        base.clearance_m,
        base.source,
        base.budget,
        &CancelGate::new_clock_free(),
    )
    .unwrap()
}

// Manufactured reservoirs are numerical constitutive tests, not empirical air.
fn state(gas: GasState, count: usize, rate_scale: f64) -> RelaxationState {
    let id = SpeciesId::new("synthetic").unwrap();
    RelaxationAbsorption::new(vec![SpeciesRelaxation {
        species: id.clone(),
        processes: (0..count)
            .map(|j| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + 600.0 * j as f64,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 101325.0 * (2e-5 + j as f64 * 3e-4) * rate_scale,
                    t_ref_k: 293.15,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    }])
    .unwrap()
    .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
    .unwrap()
}

fn edges(model: &RadiatingNetwork, count: usize, scale: f64) -> Vec<MolecularEdge> {
    model
        .driving_point()
        .network()
        .edges
        .iter()
        .map(|edge| MolecularEdge {
            edge_id: edge.id.clone(),
            state: state(
                *model.media().unwrap().edge_gas(&edge.id).unwrap(),
                count,
                scale,
            ),
            wall: None,
            cones: ConeSlices::Fixed(16),
        })
        .collect()
}

// Independent pressure-release input of the material's continuous complex
// transmission line. No duct transfer, graph or scattering routine is used.
fn shorted_line(state: &RelaxationState, omega: f64) -> C64 {
    let k = state.relaxing_wavenumber(omega).unwrap();
    let z = C64::from_re((core::f64::consts::PI * 0.006_f64.powi(2)).recip());
    let z = z.scale(omega * state.gas().density) / k;
    let phase = k.scale(0.13);
    let (s, c) = phase.re.sin_cos();
    let sine = C64::new(s * phase.im.cosh(), c * phase.im.sinh());
    let cosine = C64::new(c * phase.im.cosh(), -s * phase.im.sinh());
    C64::new(0.0, -1.0) * z * sine / cosine
}

#[test]
fn g1_empty_reservoir_limit_preserves_original_loaded_fields_and_source_identity() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let molecular = model
        .bind_molecular(edges(&model, 0, 1.0), 64, &gate)
        .unwrap();
    assert!(std::ptr::eq(molecular.model(), &model));
    assert_eq!(
        molecular.molecular_network().media(),
        model.media().unwrap()
    );
    for frequency in [800.0, 1200.0, 1600.0] {
        let omega = core::f64::consts::TAU * frequency;
        let bank = bank(omega);
        let physical = molecular
            .prepare_sources(omega, &[], Some(&bank), &[], &gate)
            .unwrap();
        let legacy = model
            .prepare_sources(omega, LossModel::Lossless, &[], Some(&bank), &[], &gate)
            .unwrap();
        assert!(physical.loss_model().is_none());
        assert!(matches!(legacy.loss_model(), Some(LossModel::Lossless)));
        assert!(legacy.molecular_network().is_none());
        assert!(std::ptr::eq(
            physical.molecular_network().unwrap(),
            molecular.molecular_network()
        ));
        assert_eq!(physical.boundary().z, legacy.boundary().z);
        let mut actual = physical.response_buffer();
        let mut expected = legacy.response_buffer();
        physical
            .drive_into(physical.default_excitations(), &mut actual, &gate)
            .unwrap();
        legacy
            .drive_into(legacy.default_excitations(), &mut expected, &gate)
            .unwrap();
        for (a, b) in actual
            .receiver_fields()
            .unwrap()
            .iter()
            .zip(expected.receiver_fields().unwrap())
        {
            close(a.pressure_pa, b.pressure_pa);
            for axis in 0..3 {
                close(a.particle_velocity_m_s[axis], b.particle_velocity_m_s[axis]);
            }
        }
        real_close(
            actual.network().unwrap().generator_power_w,
            expected.network().unwrap().generator_power_w,
        );
    }
}

#[test]
fn g1_dark_exterior_mode_retains_analytic_molecular_absorption_and_generator_work() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let molecular = model
        .bind_molecular(edges(&model, 2, 4.0), 64, &gate)
        .unwrap();
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    let prepared = molecular
        .prepare_sources(omega, &[], Some(&bank), &[], &gate)
        .unwrap();
    let mut common = prepared.response_buffer();
    let mut opposed = prepared.response_buffer();
    prepared
        .drive_into(&[C64::ONE; 2], &mut common, &gate)
        .unwrap();
    prepared
        .drive_into(&[C64::ONE, -C64::ONE], &mut opposed, &gate)
        .unwrap();
    let zin = shorted_line(
        molecular.molecular_network().edge("stem").unwrap().state(),
        omega,
    );
    let source_z = bank.impedance.z[0] - bank.impedance.z[1];
    let flow = C64::ONE / (source_z + zin);
    let network = opposed.network().unwrap();
    let aux = prepared
        .ports()
        .iter()
        .position(|port| port == "aux")
        .unwrap();
    close(network.fields.port_inflows[aux], flow);
    real_close(network.fields.network_power_w, zin.re * flow.norm_sq());
    real_close(network.source_dissipation_w, source_z.re * flow.norm_sq());
    real_close(network.generator_power_w, flow.re);
    assert!(zin.re > 0.0 && network.fields.network_power_w > 0.0);
    for (dark, bright) in opposed
        .receiver_fields()
        .unwrap()
        .iter()
        .zip(common.receiver_fields().unwrap())
    {
        assert!(dark.pressure_pa.abs() < 1e-10 * bright.pressure_pa.abs());
    }
    assert!(
        opposed.exterior_surface_power_w().unwrap().abs()
            < 1e-18 * common.exterior_surface_power_w().unwrap()
    );
    // A separately assembled physical velocity boundary checks the observable
    // field, independently of the cached unit-opening receiver composition.
    let direct = direct_fields(&model, omega, common.outward_flows_m3_s().unwrap());
    for (actual, expected) in common.receiver_fields().unwrap().iter().zip(direct) {
        close(actual.pressure_pa, expected.pressure_pa);
        for axis in 0..3 {
            close(
                actual.particle_velocity_m_s[axis],
                expected.particle_velocity_m_s[axis],
            );
        }
    }
}

#[test]
fn g3_collision_rates_change_loaded_sound_without_changing_the_exterior() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let slow = model
        .bind_molecular(edges(&model, 2, 8.0), 64, &gate)
        .unwrap();
    let fast = model
        .bind_molecular(edges(&model, 2, 0.1), 64, &gate)
        .unwrap();
    let mut effect = 0.0_f64;
    for frequency in [800.0, 1200.0, 1600.0] {
        let omega = core::f64::consts::TAU * frequency;
        let bank = bank(omega);
        let a = slow
            .prepare_sources(omega, &[], Some(&bank), &[], &gate)
            .unwrap();
        let b = fast
            .prepare_sources(omega, &[], Some(&bank), &[], &gate)
            .unwrap();
        assert_eq!(a.boundary().z, b.boundary().z);
        assert_eq!(a.boundary().source_id, b.boundary().source_id);
        let mut x = a.response_buffer();
        let mut y = b.response_buffer();
        a.drive_into(a.default_excitations(), &mut x, &gate)
            .unwrap();
        b.drive_into(b.default_excitations(), &mut y, &gate)
            .unwrap();
        for output in [&x, &y] {
            let n = output.network().unwrap();
            assert!(n.fields.network_power_w > 0.0);
            real_close(
                n.generator_power_w,
                n.source_dissipation_w + n.fields.network_power_w + n.fields.external_load_power_w,
            );
        }
        effect = effect.max(
            (x.receiver_fields().unwrap()[0].pressure_pa
                - y.receiver_fields().unwrap()[0].pressure_pa)
                .abs()
                / x.receiver_fields().unwrap()[0].pressure_pa.abs(),
        );
    }
    assert!(
        effect > 1e-3,
        "actual receiver pressure must retain the molecular collision law"
    );
}

#[test]
fn g1_molecular_opening_coherence_reaches_directional_microphones_and_far_field() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let molecular = model
        .bind_molecular(edges(&model, 2, 4.0), 64, &gate)
        .unwrap();
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    let directions = [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]];
    let options = PreparedRadiationOptions {
        angular: helmholtz::DirectivityProjectionOptions::default(),
        initial_latitudes: 8,
        refinements: 3,
        maximum_kernel_evaluations: 4_000_000,
        maximum_source_entries: 4096,
    };
    let prepared = molecular
        .prepare_sources_with_radiation(omega, &[], Some(&bank), &[], &directions, options, &gate)
        .unwrap();
    let mut bright = prepared.response_buffer();
    prepared
        .drive_into(&[C64::ONE; 2], &mut bright, &gate)
        .unwrap();
    let mut statistics = prepared.statistics_buffer();
    prepared
        .independent_into(
            &[1.0; 2],
            NetworkCoherenceBudget::default(),
            FIELD_BUDGET,
            &mut statistics,
            &gate,
        )
        .unwrap();
    for (actual, field) in statistics
        .receiver_fields()
        .unwrap()
        .iter()
        .zip(bright.receiver_fields().unwrap())
    {
        real_close(
            actual.mean_square_pressure_pa2,
            0.25 * field.pressure_pa.norm_sq(),
        );
    }
    let patterns = [
        FirstOrder::new(0.5, [1.0, 0.0, 0.0]).unwrap(),
        FirstOrder::new(0.2, [0.0, 1.0, 0.0]).unwrap(),
    ];
    let observed: Vec<_> = patterns
        .iter()
        .zip(bright.receiver_fields().unwrap())
        .map(|(pattern, field)| {
            pattern
                .observe(
                    field.pressure_pa,
                    field.particle_velocity_m_s,
                    Medium::from(&model.gas()),
                )
                .unwrap()
        })
        .collect();
    let mut microphones = prepared.microphones(&patterns).unwrap();
    let spectra = microphones
        .independent(&[1.0; 2], FIELD_BUDGET, &gate)
        .unwrap();
    for i in 0..2 {
        for j in 0..2 {
            close(
                spectra[i * 2 + j],
                (observed[i] * observed[j].conj()).scale(0.25),
            );
        }
    }
    let far = prepared.far_field().unwrap();
    let angular_budget = HarmonicWorkBudget {
        maximum_products: 1_000_000,
    };
    let mut field = [C64::ZERO; 2];
    let power = far
        .evaluate_coherent(&[C64::ONE; 2], angular_budget, &mut field, &gate)
        .unwrap();
    assert!(power.total_power_w > 0.0);
    let coherent = field;
    far.observe_coherent(&[C64::ONE, -C64::ONE], angular_budget, &mut field, &gate)
        .unwrap();
    for (dark, bright) in field.iter().zip(coherent) {
        assert!(dark.abs() < 1e-10 * bright.abs());
    }
    let mut intensity = [0.0; 2];
    far.observe_independent_sources(&[1.0; 2], angular_budget, &mut intensity, &gate)
        .unwrap();
    for (actual, amplitude) in intensity.iter().zip(coherent) {
        real_close(
            *actual,
            amplitude.norm_sq() / (4.0 * model.gas().characteristic_impedance),
        );
    }
    let work = statistics.network().unwrap().work().unwrap();
    real_close(
        work.generator_power_w,
        work.source_dissipation_w + work.network_power_w + work.external_load_power_w,
    );
}

#[test]
fn g4_original_records_edge_resolution_frequency_support_and_cancellation_refuse() {
    let model = regional();
    let gate = CancelGate::new_clock_free();
    let bindings = edges(&model, 1, 1.0);
    assert!(
        two_source_model()
            .bind_molecular(bindings.clone(), 64, &gate)
            .is_err()
    );
    assert!(
        model
            .bind_molecular(bindings[..3].to_vec(), 64, &gate)
            .is_err()
    );
    assert!(model.bind_molecular(bindings.clone(), 3, &gate).is_err());
    let mut wrong = bindings.clone();
    wrong[0].state = state(
        GasState::try_new(
            &fs_material::gas::GasSpec::dry_air_ussa1976(),
            320.0,
            101_325.0,
        )
        .unwrap(),
        1,
        1.0,
    );
    assert!(model.bind_molecular(wrong, 64, &gate).is_err());
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(matches!(
        model.bind_molecular(bindings.clone(), 64, &stopped),
        Err(ExteriorNetworkError::Cancelled)
    ));
    let wall = MolecularWallResolution {
        poles_per_field: 16,
        band_rad_s: [100.0, 10_000.0],
        maximum_relative_error: 0.02,
    };
    let bindings: Vec<_> = bindings
        .into_iter()
        .map(|mut edge| {
            edge.wall = Some(wall);
            edge
        })
        .collect();
    let molecular = model.bind_molecular(bindings, 64, &gate).unwrap();
    let omega = core::f64::consts::TAU * 1200.0;
    let bank = bank(omega);
    assert!(matches!(
        molecular.prepare_sources(omega, &[], Some(&bank), &[], &stopped),
        Err(ExteriorNetworkError::Cancelled)
    ));
    let outside = super::bank(11_000.0);
    assert!(
        molecular
            .prepare_sources(11_000.0, &[], Some(&outside), &[], &gate)
            .is_err()
    );
    let prepared = molecular
        .prepare_sources(omega, &[], Some(&bank), &[], &gate)
        .unwrap();
    let mut response = prepared.response_buffer();
    prepared
        .drive_into(prepared.default_excitations(), &mut response, &gate)
        .unwrap();
    assert!(response.is_valid());
    assert!(
        prepared
            .drive_into(prepared.default_excitations(), &mut response, &stopped)
            .is_err()
    );
    assert!(response.receiver_fields().is_none() && response.network().is_none());
    prepared
        .drive_into(prepared.default_excitations(), &mut response, &gate)
        .unwrap();
    assert!(response.is_valid());
}
