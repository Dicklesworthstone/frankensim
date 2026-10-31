//! Real-gas volume work through the original PR/NASA and parcel owners.
use super::*;
use crate::thermochemical_gas::energy::entropy::{GasEntropyError, GasEntropySolve};
use crate::thermochemical_gas::tests::model_in_domains;
use fs_thermochem::{Nasa9RegionV1, Nasa9StandardStateModelV1};

fn entropy_controls() -> GasEntropySolve {
    GasEntropySolve {
        minimum_temperature_k: 210.0,
        maximum_temperature_k: 450.0,
        entropy_tolerance_j_kg_k: 1.0e-9,
        maximum_iterations: 128,
    }
}

#[test]
fn real_gas_volume_preserves_mass_and_returns_signed_work_complete_views_and_round_trip() {
    let initial = parcel();
    let saved = initial.clone();
    // Independent 75-digit Decimal integration of PR p(T,v) and NASA Cv/T,
    // including R*ln((v-b)/(v0-b)) and the analytic a'(T) volume logarithm.
    // This oracle does not use the production residual Gibbs/entropy identity.
    for (volume, expected_temperature, expected_pressure, expected_work) in [
        (
            0.008,
            331.606_449_235_535_6,
            6_981_183.016_655_789,
            -11_791.630_070_073_763,
        ),
        (
            0.011,
            287.505_188_450_419_14,
            4_343_285.645_717_975,
            4_658.361_642_714_088,
        ),
    ] {
        let (candidate, step) = initial
            .at_volume_isentropic(volume, entropy_controls(), || false)
            .unwrap();
        near(
            candidate.state().temperature_k(),
            expected_temperature,
            2.0e-8,
        );
        near(candidate.state().pressure_pa(), expected_pressure, 1.0e-3);
        near(step.work_by_gas_j, expected_work, 2.0e-5);
        assert_eq!(initial, saved);
        assert_eq!(candidate.mass_kg(), initial.mass_kg());
        assert_eq!(candidate.volume_m3(), volume);
        assert_eq!(candidate.density_kg_m3(), initial.mass_kg() / volume);
        assert_eq!(candidate.model(), initial.model());
        assert_eq!(candidate.point().axes().get("material-lot"), Some(&12.0));
        assert_eq!(step.work_by_gas_j, -step.stored_energy_change_j);
        assert_eq!(step.initial_volume_m3, initial.volume_m3());
        assert_eq!(step.final_volume_m3, volume);
        assert_eq!(
            step.work_by_gas_j.is_sign_positive(),
            volume > initial.volume_m3()
        );
        assert!(
            step.entropy_change_j_k.abs()
                <= initial.mass_kg() * entropy_controls().entropy_tolerance_j_kg_k
        );
        let pressure = candidate.state().pressure_pa();
        let original_pressure = initial.state().pressure_pa();
        let secant = step.work_by_gas_j / (volume - initial.volume_m3());
        assert!(
            secant > pressure.min(original_pressure) && secant < pressure.max(original_pressure)
        );
        let fluid = candidate.state().transport_state().unwrap();
        let acoustic = candidate.state().acoustic_state().unwrap();
        assert_eq!(fluid.density_kg_m3(), candidate.density_kg_m3());
        assert_eq!(acoustic.density_kg_m3(), candidate.density_kg_m3());
        assert!(candidate.state().specific_heat_cv_j_kg_k() > 0.0);
        assert!(candidate.state().root().pressure_derivatives().1 < 0.0);
        let (returned, reverse) = candidate
            .at_volume_isentropic(initial.volume_m3(), entropy_controls(), || false)
            .unwrap();
        near(
            returned.state().temperature_k(),
            initial.state().temperature_k(),
            2.0e-8,
        );
        near(returned.state().pressure_pa(), original_pressure, 1.0e-3);
        near(step.work_by_gas_j + reverse.work_by_gas_j, 0.0, 2.0e-5);
        assert_eq!(returned.mass_kg(), initial.mass_kg());
    }
}

#[test]
fn real_gas_volume_work_matches_independent_pressure_quadrature() {
    let initial = parcel();
    let final_volume = 0.008;
    let (_, report) = initial
        .at_volume_isentropic(final_volume, entropy_controls(), || false)
        .unwrap();
    // Independently integrate source-resolved endpoint pressure along the same
    // isentrope; the reported work uses only the source internal-energy change.
    let intervals = 32;
    let h = (final_volume - initial.volume_m3()) / f64::from(intervals);
    let mut weighted_pressure = 0.0;
    for i in 0..=intervals {
        let volume = initial.volume_m3() + f64::from(i) * h;
        let (state, _) = initial
            .at_volume_isentropic(volume, entropy_controls(), || false)
            .unwrap();
        let weight = if i == 0 || i == intervals {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        weighted_pressure += weight * state.state().pressure_pa();
    }
    let quadrature = weighted_pressure * h / 3.0;
    assert!((report.work_by_gas_j / quadrature - 1.0).abs() < 2.0e-7);
}

#[test]
fn real_gas_volume_approaches_analytic_dilute_gamma_law() {
    // Explicit synthetic transport support includes the dilute pressure range;
    // no production material domain is extrapolated to create this oracle.
    let model = bound(model_in_domains(
        vec![0.75, 0.25],
        0.0,
        [200.0, 2000.0],
        [200.0, 2000.0],
        [0.01, 1.0e7],
    ));
    let initial =
        IsochoricRealGasParcel::new(model, point(300.0, 0.1), 0.01, PengRobinsonPhaseV1::Vapor)
            .unwrap();
    let gamma = 3.25 / 2.25;
    let cv = 2.25 * fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K / 0.022;
    for volume in [0.005, 0.02] {
        let (candidate, step) = initial
            .at_volume_isentropic(volume, entropy_controls(), || false)
            .unwrap();
        let ratio = initial.volume_m3() / volume;
        let temperature = 300.0 * ratio.powf(gamma - 1.0);
        let pressure = 0.1 * ratio.powf(gamma);
        let work = initial.mass_kg() * cv * (300.0 - temperature);
        assert!((candidate.state().temperature_k() / temperature - 1.0).abs() < 2.0e-8);
        assert!((candidate.state().pressure_pa() / pressure - 1.0).abs() < 2.0e-8);
        assert!((step.work_by_gas_j / work - 1.0).abs() < 2.0e-7);
    }
}

pub(super) fn with_references(energy_reference: f64, entropy_reference: f64) -> IsochoricRealGasParcel {
    let original = model(vec![0.75, 0.25], 0.001);
    let components = original
        .thermodynamics()
        .components()
        .iter()
        .map(|source| {
            let regions = source
                .regions()
                .iter()
                .map(|region| {
                    let mut card = region.card().clone();
                    card.parameters.get_mut("a7").unwrap().value += energy_reference;
                    card.parameters.get_mut("a8").unwrap().value += entropy_reference;
                    Nasa9RegionV1::from_card(card).unwrap()
                })
                .collect();
            Nasa9StandardStateModelV1::new(
                source.species().clone(),
                source.molar_mass(),
                source.convention().clone(),
                regions,
            )
            .unwrap()
        })
        .collect();
    let source = ThermochemicalGasModel::new(
        FrozenIdealGasMixtureModelV1::new(
            components,
            original.thermodynamics().declared_composition().clone(),
        )
        .unwrap(),
        original.transport_materials().to_vec(),
        None,
    )
    .unwrap();
    IsochoricRealGasParcel::new(
        bound(source),
        point(300.0, 5.0e6),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

#[test]
fn real_gas_volume_ignores_large_energy_and_entropy_reference_offsets() {
    let plain = with_references(0.0, 0.0);
    let shifted = with_references(1.0e20, 1.0e24);
    for volume in [0.008, 0.012] {
        let (a, first) = plain
            .at_volume_isentropic(volume, entropy_controls(), || false)
            .unwrap();
        let (b, second) = shifted
            .at_volume_isentropic(volume, entropy_controls(), || false)
            .unwrap();
        assert_eq!(first, second);
        assert_eq!(a.mass_kg(), b.mass_kg());
        assert_eq!(a.density_kg_m3(), b.density_kg_m3());
        assert_eq!(a.state().temperature_k(), b.state().temperature_k());
        assert_eq!(a.state().pressure_pa(), b.state().pressure_pa());
        assert_eq!(a.state().sound_speed_m_s(), b.state().sound_speed_m_s());
    }
}

#[test]
fn real_gas_volume_candidate_continues_original_heating_at_the_new_density() {
    let initial = parcel();
    let (mut candidate, _) = initial
        .at_volume_isentropic(0.008, entropy_controls(), || false)
        .unwrap();
    let before = candidate.clone();
    let heat = candidate.mass_kg() * 500.0;
    let report = candidate.add_heat(heat, controls(), || false).unwrap();
    assert_eq!(candidate.mass_kg(), initial.mass_kg());
    assert_eq!(candidate.volume_m3(), before.volume_m3());
    assert_eq!(candidate.density_kg_m3(), before.density_kg_m3());
    assert!(candidate.state().temperature_k() > before.state().temperature_k());
    assert!(candidate.state().pressure_pa() > before.state().pressure_pa());
    near(
        report.stored_energy_change_j,
        heat,
        candidate.mass_kg() * controls().energy_tolerance_j_kg,
    );
}

#[test]
fn real_gas_volume_refuses_domain_budget_and_cancellation_without_exposing_a_candidate() {
    let initial = parcel();
    let saved = initial.clone();
    let (same, identity) = initial
        .at_volume_isentropic(initial.volume_m3(), entropy_controls(), || false)
        .unwrap();
    assert_eq!(same, initial);
    assert_eq!(identity.work_by_gas_j, 0.0);
    assert_eq!(identity.iterations, 0);
    for volume in [0.0, -1.0, f64::NAN, f64::INFINITY, 1.0e-12, 1.0e6] {
        assert!(
            initial
                .at_volume_isentropic(volume, entropy_controls(), || false)
                .is_err()
        );
        assert_eq!(initial, saved);
    }
    for controls in [
        GasEntropySolve {
            maximum_iterations: 0,
            ..entropy_controls()
        },
        GasEntropySolve {
            minimum_temperature_k: 310.0,
            ..entropy_controls()
        },
        GasEntropySolve {
            entropy_tolerance_j_kg_k: f64::NAN,
            ..entropy_controls()
        },
    ] {
        assert!(
            initial
                .at_volume_isentropic(initial.volume_m3(), controls, || false)
                .is_err()
        );
    }
    assert!(matches!(
        initial.at_volume_isentropic(
            0.008,
            GasEntropySolve {
                maximum_iterations: 1,
                ..entropy_controls()
            },
            || false
        ),
        Err(GasEntropyError::NotConverged { .. })
    ));
    let dense = IsochoricRealGasParcel::new(
        bound(model(vec![0.75, 0.25], 0.001)),
        point(300.0, 9.0e6),
        0.01,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap();
    // Final transport lies outside the original card's 10 MPa pressure domain.
    assert!(matches!(
        dense.at_volume_isentropic(0.008, entropy_controls(), || false),
        Err(GasEntropyError::Material(_))
    ));
    let mut total_polls = 0;
    let expected = initial
        .at_volume_isentropic(0.008, entropy_controls(), || {
            total_polls += 1;
            false
        })
        .unwrap();
    for stop in [1, total_polls / 2, total_polls] {
        let mut polls = 0;
        assert!(matches!(
            initial.at_volume_isentropic(0.008, entropy_controls(), || {
                polls += 1;
                polls >= stop
            }),
            Err(GasEntropyError::Cancelled)
        ));
        assert_eq!(initial, saved);
    }
    assert_eq!(
        initial
            .at_volume_isentropic(0.008, entropy_controls(), || false)
            .unwrap(),
        expected
    );
}
