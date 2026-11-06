//! G0/G3 native brass-to-material-receiver integration. All fixture gas and
//! lip data are synthetic; these checks make no measured-instrument claim.
use super::*;
use crate::air_path::GasAbsorptionModel;
use crate::pcm_wav::baffled::gas::GasPropagationLimits;
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw,
    SpeciesRelaxation, VibrationalMode,
};
use fs_thermochem::SpeciesId;

const CLASSICAL: &str = "observer_gas stokes-kirchhoff 4096 1000000 50000000\n";
const MOLECULAR: &str = "observer_gas molecular 4096 1000000 50000000 declared 2 N2 1 vibration 1 1000 1 power 10 300 0 250 350 O2 0\n";
const AMBIENT: &str = "mixture 300 101325 mole 2 N2 0.8 O2 0.2";

fn fixture(ambient: &str, extra: &str) -> String {
    format!(
        "frankensim-lip-performance-v1\n\
         audio 24000 512\nfull_scale_pa 30000\nambient {ambient}\n\
         lip_response 200 7 0.1 0.012 0.0005 0.0001\n\
         lip_source synthetic integration fixture\nbore_source synthetic integration fixture\n\
         cup_volume_m3 1.5e-6\ntermination flanged\nmodes 1 1\n\
         cylinder 0.006 0.30\nmicrophone 0.02 0.01 0.15 2 8 1000\n\
         {extra}schedule\nfrankensim-gesture-schedule-v1\ncontrol_rate_hz\t700\n\
         tracks\t1\ntrack\tblow\tblowing-pressure\ninitial\tpressure-pa\t2.8e3\nevents\t0\n"
    )
}

fn limits() -> GasPropagationLimits {
    GasPropagationLimits {
        max_unique_paths: 4096,
        max_filter_taps: 1_000_000,
        max_kernel_terms: 50_000_000,
    }
}

fn molecular(gas: &GasRecord) -> GasPropagation {
    // Construct constitutive data independently of the native text reader.
    let law = RelaxationAbsorption::new(vec![
        SpeciesRelaxation {
            species: SpeciesId::new("N2").unwrap(),
            processes: vec![RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1000.0,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: 10.0,
                    t_ref_k: 300.0,
                    exponent: 0.0,
                    t_min_k: 250.0,
                    t_max_k: 350.0,
                },
            }],
        },
        SpeciesRelaxation {
            species: SpeciesId::new("O2").unwrap(),
            processes: vec![],
        },
    ])
    .unwrap();
    let medium = fs_material::relaxation::input::medium_from_record(gas).unwrap();
    GasPropagation {
        law: GasAbsorptionModel::MolecularRelaxation,
        relaxation: Some(law.evaluate(&medium).unwrap()),
        limits: limits(),
    }
}

#[test]
fn observer_gas_binds_original_species_in_either_record_order_and_refuses_unsupported_inputs() {
    for row in [
        CLASSICAL,
        MOLECULAR,
        "observer_gas lossless 4096 1000000 50000000\n",
    ] {
        let text = fixture(AMBIENT, row);
        let early = text
            .replace(row, "")
            .replace("ambient ", &format!("{row}ambient "));
        for source in [&text, &early] {
            let description = LipDescription::from_bytes(source.as_bytes()).unwrap();
            assert_eq!(description.observer_gas().unwrap().limits, limits());
            assert_eq!(
                description.observer_gas().unwrap().relaxation.is_some(),
                row == MOLECULAR
            );
        }
    }
    assert!(
        LipDescription::from_bytes(fixture(AMBIENT, "").as_bytes())
            .unwrap()
            .observer_gas()
            .is_none()
    );
    for row in [
        "observer_gas\n",
        "observer_gas guessed-air 4096 1000000 50000000\n",
        "observer_gas stokes-kirchhoff 0 1000000 50000000\n",
        "observer_gas stokes-kirchhoff 4096 4000001 50000000\n",
        "observer_gas stokes-kirchhoff 4096 1000000 50000000 extra\n",
        "observer_gas molecular 4096 1000000 50000000 declared 1 N2 0\n",
    ] {
        let error = LipDescription::from_bytes(fixture(AMBIENT, row).as_bytes()).unwrap_err();
        assert!(error.contains("observer_gas"), "{error}");
    }
    let duplicate = format!("{CLASSICAL}{CLASSICAL}");
    assert!(
        LipDescription::from_bytes(fixture(AMBIENT, &duplicate).as_bytes())
            .unwrap_err()
            .contains("may only appear once")
    );
    let modal = format!("{CLASSICAL}modal_baffle 1e-3 1000000\n");
    let error = LipDescription::from_bytes(fixture(AMBIENT, &modal).as_bytes()).unwrap_err();
    assert!(
        error.contains("observer_gas") && error.contains("modal_baffle"),
        "{error}"
    );
    let error = LipDescription::from_bytes(fixture("moist_air 300 101325 0", MOLECULAR).as_bytes())
        .unwrap_err();
    assert!(
        error.contains("explicit mixture or thermochemical species"),
        "{error}"
    );
}

fn render(description: &LipDescription, block: usize) -> (Vec<f64>, Vec<f64>) {
    let mut performance = description.build().unwrap();
    performance.record_mouthpiece();
    let mut pressure = vec![0.0; description.samples as usize];
    for output in pressure.chunks_mut(block) {
        performance.block(output).unwrap();
    }
    (pressure, performance.mouthpiece_pressure().to_vec())
}

fn compare_with_original_owners(description: &LipDescription, propagation: &GasPropagation) {
    let mut native = description.build().unwrap();
    native.record_mouthpiece();
    let mut original = description.clone();
    original.observer_gas = None;
    let mut direct = original.build().unwrap();
    assert!(direct.info().observer_gas.is_none());
    let ratio = direct.ratio;
    let gas = description.ambient_gas();
    let (mut receiver, report) = BaffledPressure::circular_openings_multirate_gas_differentiating(
        &[CircularOpening {
            center_m: [0.0; 2],
            radius_m: direct.info().mouth_radius_m,
        }],
        description.rate_hz * ratio as u32,
        description.rate_hz,
        description.microphone,
        gas,
        propagation,
        &mut || false,
    )
    .unwrap();
    let actual = native.info().observer_gas.as_ref().unwrap();
    assert_eq!(actual.gas, *gas);
    assert_eq!(actual.law, propagation.law);
    assert_eq!(actual.surface_filter_indices.len(), 16);
    assert!(
        actual.filters.len() > 1,
        "off-axis original path lengths are distinct"
    );
    assert_eq!(actual.surface_filter_indices, report.surface_filter_indices);
    assert_eq!(actual.filters, report.filters);
    assert_eq!(actual.kernel_terms, report.kernel_terms);
    assert_eq!(
        actual.additional_latency_samples,
        report.additional_latency_samples
    );
    assert_eq!(
        native.filter_delay_samples(),
        receiver.filter_delay_samples()
    );

    // Advance the original source owners sample by sample. Every nonlinear
    // mechanics substep enters the shared receiver once, in its original order.
    direct
        .engine
        .apply(BrassControl::SetBlowingPressure(2800.0))
        .unwrap();
    let mut expected = Vec::new();
    let mut original_pressure = Vec::new();
    let mut mouthpiece = Vec::new();
    let mut flows = vec![0.0; ratio];
    for _ in 0..description.samples {
        let mut internal = [0.0];
        match &mut direct.engine {
            Engine::Linear { voice, mouth } => {
                let mut outgoing = [0.0];
                voice
                    .step_block_outgoing(&mut internal, &mut outgoing)
                    .unwrap();
                assert_eq!(mouth.len(), 1);
                flows[0] = mouth[0].push(outgoing[0]).unwrap();
            }
            #[cfg(feature = "nonlinear-brass")]
            Engine::Nonlinear(lip) => lip.step_block(&mut internal, &mut flows).unwrap(),
        }
        expected.push(receiver.step_trace(&flows).unwrap());
        original_pressure.push(direct.receiver.step_trace(&flows).unwrap());
        mouthpiece.push(internal[0]);
    }
    let mut pressure = vec![0.0; description.samples as usize];
    for output in pressure.chunks_mut(17) {
        native.block(output).unwrap();
    }
    assert_eq!(
        pressure, expected,
        "native material observation must be the shared physical owner"
    );
    assert_eq!(
        native.mouthpiece_pressure(),
        mouthpiece,
        "propagation cannot alter the source"
    );
    assert!(pressure.iter().all(|p| p.is_finite()));
    assert!(pressure.iter().any(|p| p.abs() > 1e-12));
    assert_ne!(
        pressure, original_pressure,
        "the authored gas law must reach audible pressure"
    );
    let partitioned = render(description, 31);
    assert_eq!(pressure, partitioned.0);
    assert_eq!(mouthpiece, partitioned.1);
}

#[test]
fn molecular_linear_performance_matches_source_flows_and_shared_material_receiver() {
    let description = LipDescription::from_bytes(fixture(AMBIENT, MOLECULAR).as_bytes()).unwrap();
    compare_with_original_owners(&description, &molecular(description.gas_record()));
}

#[cfg(feature = "nonlinear-brass")]
#[test]
fn classical_nonlinear_performance_retains_every_mechanics_substep_and_original_material() {
    use fs_material::gas_species::input::ThermochemicalPhase;
    let packs = gas_test_support::packs("specimen-gas", 0.029, 3.5, 200.0);
    let bindings = [packs.binding()];
    let sources = ThermochemicalSources {
        species: &bindings,
        water_species: None,
        phase: ThermochemicalPhase::FrozenGas,
    };
    let rows = format!("bore nonlinear 0.01 4 4 64 20 8000 1e-2\n{CLASSICAL}");
    let source = fixture("thermochemical 300 101325 mole 1 specimen-gas 1", &rows);
    let description =
        LipDescription::from_bytes_with_gas_sources(source.as_bytes(), Some(&sources)).unwrap();
    assert_eq!(
        description.nonlinear_bore().unwrap().0.mechanics_per_output,
        4
    );
    let GasRecord::Thermochemical(gas) = description.gas_record() else {
        panic!("retain the original thermochemical source");
    };
    assert_eq!(gas.source_packs()[0].nasa9_pack, packs.nasa9.content_hash());
    compare_with_original_owners(&description, &GasPropagation::classical(limits()));
}

#[test]
fn declared_observer_budgets_reach_the_actual_surface_owner() {
    for row in [
        "observer_gas stokes-kirchhoff 1 1000000 50000000\n",
        "observer_gas stokes-kirchhoff 4096 1 50000000\n",
        "observer_gas stokes-kirchhoff 4096 1000000 1\n",
    ] {
        let description = LipDescription::from_bytes(fixture(AMBIENT, row).as_bytes()).unwrap();
        // Exercise the exact shared construction used by both engines without
        // rebuilding an unrelated bore bank for each cold-work refusal.
        let error = description
            .scalar_receiver(0.006, description.rate_hz)
            .err()
            .unwrap();
        assert!(error.contains("observer_gas"), "{error}");
    }
}
