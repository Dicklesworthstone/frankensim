//! Original material/PR owners: simultaneous energy carriage and thermal windows.
use super::super::thermal::{
    RealGasThermodiffusionModel, RealGasThermodiffusionSolve,
    network::{
        self as thermal, RealGasThermodiffusionLink, RealGasThermodiffusionNetworkAdvanceControls,
    },
};
use super::*;
use crate::state_point::{THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY};
use crate::thermochemical_gas::{
    energy::parcel::transfer::species::diffusion::{
        network::adaptive::{DiffusionNetworkAdvanceControls, DiffusionNetworkAdvanceError},
        thermal::{
            GasThermodiffusionMaterial, HEAT_OF_TRANSPORT_DIMS, HEAT_OF_TRANSPORT_PROPERTY,
            SpeciesHeatOfTransportMaterial,
        },
    },
    nonideal::energy::parcel::transfer::species::network as native,
};
use fs_qty::Dims;

fn scalar(property: &str, dims: Dims, value: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(property, dims),
            value: PropertyValue::Curve {
                abscissa: "temperature".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(250., value), (450., value)],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "temperature",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::AbsoluteTemperature,
                        ValueForm::Static,
                    )),
                    250.,
                    450.,
                )
                .with_quantity(
                    "pressure",
                    QuantitySpec::semantic(SemanticType::new(
                        QuantityKind::Pressure,
                        ValueForm::Static,
                    )),
                    1000.,
                    5e6,
                ),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "authored PR thermal fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: property.into(),
            phase: "gas".into(),
            process: "thermal fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn prepared<'a>(
    gas: &IsochoricRealGasParcel,
    binary: &'a MaterialCard,
    conductivity: &'a MaterialCard,
    heats: [&'a MaterialCard; 2],
    context: &'a QueryPoint,
) -> RealGasThermodiffusionModel<'a> {
    let source = |card: &'a MaterialCard| GasThermodiffusionMaterial {
        material: card,
        context,
        selection: MaterialPropertySelection::SingleClaimOnly,
    };
    RealGasThermodiffusionModel::new(
        gas.model(),
        &[binding(binary, context)],
        source(conductivity),
        &[
            SpeciesHeatOfTransportMaterial {
                species: SpeciesId::new("A").unwrap(),
                source: source(heats[0]),
            },
            SpeciesHeatOfTransportMaterial {
                species: SpeciesId::new("B").unwrap(),
                source: source(heats[1]),
            },
        ],
        || false,
    )
    .unwrap()
}
fn thermal_controls() -> RealGasThermodiffusionSolve {
    RealGasThermodiffusionSolve {
        maximum_temperature_difference_k: 80.,
        ..controls()
    }
}
fn original(references: [f64; 2]) -> BTreeMap<u64, IsochoricRealGasParcel> {
    BTreeMap::from([
        (
            10,
            at(real_model(source(0.75, references), 0.03), 310., 4e6, 11.),
        ),
        (
            20,
            at(real_model(source(0.35, references), 0.03), 350., 4e6, 22.),
        ),
    ])
}
fn temporal() -> RealGasThermodiffusionNetworkAdvanceControls {
    RealGasThermodiffusionNetworkAdvanceControls {
        temporal: DiffusionNetworkAdvanceControls {
            maximum_step_s: 2.,
            minimum_step_s: 1e-6,
            maximum_evaluations: 1024,
            absolute_species_tolerance_kg: 1e-6,
            relative_species_tolerance: 0.,
            temperature_tolerance_k: 0.1,
        },
        maximum_face_evaluations: 1024,
        absolute_species_transfer_tolerance_kg: 1e-6,
        absolute_additional_heat_tolerance_j: 1e-9,
    }
}

#[test]
fn dense_thermal_face_enters_one_native_species_heat_reconstruction() {
    let binary = card();
    let k = scalar(
        THERMAL_CONDUCTIVITY_PROPERTY,
        THERMAL_CONDUCTIVITY_DIMS,
        0.08,
    );
    let a = scalar(HEAT_OF_TRANSPORT_PROPERTY, HEAT_OF_TRANSPORT_DIMS, 6000.);
    let b = scalar(HEAT_OF_TRANSPORT_PROPERTY, HEAT_OF_TRANSPORT_DIMS, -9000.);
    let context = point(330., 4e6);
    let mut outputs = Vec::new();
    for refs in [[0., 0.], [1e19, -2e19]] {
        let before = original(refs);
        let saved = before.clone();
        let model = prepared(&before[&10], &binary, &k, [&a, &b], &context);
        let links = [RealGasThermodiffusionLink {
            id: 1,
            parcels: [10, 20],
            model: &model,
            geometry: geometry(),
            solve: thermal_controls(),
        }];
        let mut work = model.workspace().unwrap();
        let flux = model
            .face_flux(
                &before[&10],
                &before[&20],
                geometry(),
                thermal_controls(),
                &mut work,
                || false,
            )
            .unwrap();
        assert!(flux.entropy_production_w_k > 0.);
        close(
            flux.source_entropy_production_w_k,
            flux.entropy_production_w_k,
            1e-10,
        );
        let transfers: Vec<_> = flux
            .species
            .iter()
            .filter(|s| s.mass_rate_left_to_right_kg_s != 0.)
            .map(|s| {
                let (donor, receiver) = if s.mass_rate_left_to_right_kg_s > 0. {
                    (10, 20)
                } else {
                    (20, 10)
                };
                native::SpeciesNetworkTransfer {
                    donor,
                    receiver,
                    species: s.species.clone(),
                    mass_kg: s.mass_rate_left_to_right_kg_s.abs(),
                }
            })
            .collect();
        let heat = [native::SpeciesNetworkHeatTransfer {
            parcels: [10, 20],
            heat_left_to_right_j: flux.additional_heat_rate_left_to_right_w,
        }];
        // The face's mass-frame allowance, from the owner's published residual
        // and the same one-second masses, for both of its parcels.
        let masses: Vec<_> = flux
            .species
            .iter()
            .map(|s| s.mass_rate_left_to_right_kg_s)
            .collect();
        let allowance = super::super::network::face_allowance(
            &masses,
            flux.mass_frame_residual_kg_m2_s,
            geometry().area_m2,
            1.,
        )
        .unwrap();
        // (c) On this two-species face the owner closes the frame explicitly,
        // so its residual is at rounding level (about 7u of sum|J_i| from the
        // closure's and the products' own roundings), far inside its 1e-10
        // relative tolerance.
        let gross: f64 = masses.iter().map(|m| m.abs()).sum();
        assert!(
            flux.mass_frame_residual_kg_m2_s.abs() * geometry().area_m2
                <= 8.0 * f64::EPSILON * gross
        );
        let expected = native::exchange_candidate_with_heat(
            &before,
            &transfers,
            &heat,
            &BTreeMap::from([(10, allowance), (20, allowance)]),
            &endpoints(&before),
            physical().exchange,
            || false,
        )
        .unwrap();
        let actual =
            thermal::step_candidate(&before, &links, 1., &endpoints(&before), physical(), || {
                false
            })
            .unwrap();
        assert_eq!(actual.parcels(), expected.parcels());
        assert_eq!(actual.step().exchange, *expected.step());
        assert_eq!(actual.step().links[0].flux.as_ref().unwrap(), &flux);
        assert!(actual.step().energy_change_j.abs() <= physical().exchange.energy_tolerance_j);
        assert!(actual.step().entropy_change_j_k > 0.);
        for row in &transfers {
            close(
                species_mass(&before[&row.donor], &row.species)
                    - species_mass(&actual.parcels()[&row.donor], &row.species),
                row.mass_kg,
                physical().exchange.mass_tolerance_kg,
            );
        }
        for (&id, old) in &before {
            assert_eq!(old.volume_m3(), actual.parcels()[&id].volume_m3());
            close(
                old.mass_kg(),
                actual.parcels()[&id].mass_kg(),
                physical().exchange.mass_tolerance_kg,
            );
        }
        if refs == [0., 0.] {
            let energy: f64 = before
                .iter()
                .map(|(&id, old)| {
                    let new = &actual.parcels()[&id];
                    new.mass_kg() * new.state().specific_internal_energy_j_kg()
                        - old.mass_kg() * old.state().specific_internal_energy_j_kg()
                })
                .sum();
            let entropy: f64 = before
                .iter()
                .map(|(&id, old)| {
                    let new = &actual.parcels()[&id];
                    new.mass_kg() * new.state().specific_entropy_j_kg_k().unwrap()
                        - old.mass_kg() * old.state().specific_entropy_j_kg_k().unwrap()
                })
                .sum();
            close(energy, 0., physical().exchange.energy_tolerance_j);
            close(entropy, actual.step().entropy_change_j_k, 1e-9);
        }
        assert_eq!(before, saved);
        outputs.push(actual);
    }
    for id in [10, 20] {
        close(
            outputs[0].parcels()[&id].state().temperature_k(),
            outputs[1].parcels()[&id].state().temperature_k(),
            1e-9,
        );
    }
    close(
        outputs[0].step().entropy_change_j_k,
        outputs[1].step().entropy_change_j_k,
        1e-10,
    );
}

#[test]
fn adaptive_original_heat_integrals_refine_and_spend_all_face_work() {
    let binary = card();
    let k = scalar(
        THERMAL_CONDUCTIVITY_PROPERTY,
        THERMAL_CONDUCTIVITY_DIMS,
        0.08,
    );
    let a = scalar(HEAT_OF_TRANSPORT_PROPERTY, HEAT_OF_TRANSPORT_DIMS, 6000.);
    let b = scalar(HEAT_OF_TRANSPORT_PROPERTY, HEAT_OF_TRANSPORT_DIMS, -9000.);
    let context = point(330., 4e6);
    let before = original([0., 0.]);
    let saved = before.clone();
    let model = prepared(&before[&10], &binary, &k, [&a, &b], &context);
    let links = [RealGasThermodiffusionLink {
        id: 1,
        parcels: [10, 20],
        model: &model,
        geometry: geometry(),
        solve: thermal_controls(),
    }];
    let mut budgets = physical();
    budgets.exchange.energy_tolerance_j = 1e-6;
    let candidate = thermal::advance_candidate(
        &before,
        &links,
        2.,
        &endpoints(&before),
        budgets,
        temporal(),
        || false,
    )
    .unwrap();
    let r = candidate.report();
    assert!(r.accepted_intervals > 1 && r.rejected_intervals > 0);
    assert!(r.maximum_transfer_error_ratio.is_some_and(|x| x <= 1.));
    assert_eq!(r.face_evaluations, r.evaluations);
    assert_eq!(r.heat.len(), 1);
    assert!(r.heat[0].heat_left_to_right_j < 0.);
    assert!(r.endpoint_energy_change_j.abs() <= budgets.exchange.energy_tolerance_j);
    assert!(r.entropy_change_j_k > 0.);
    let mut cap = temporal();
    cap.maximum_face_evaluations = 1;
    assert!(matches!(
        thermal::advance_candidate(
            &before,
            &links,
            2.,
            &endpoints(&before),
            budgets,
            cap,
            || false
        ),
        Err(DiffusionNetworkAdvanceError::FaceWorkLimit {
            face_evaluations: 1,
            accepted_intervals: 0,
            ..
        })
    ));
    assert_eq!(before, saved);
}

#[test]
fn closed_thermal_paths_and_late_cancellation_do_not_publish_heat_or_state() {
    let binary = card();
    let k = scalar(
        THERMAL_CONDUCTIVITY_PROPERTY,
        THERMAL_CONDUCTIVITY_DIMS,
        0.08,
    );
    let a = scalar(HEAT_OF_TRANSPORT_PROPERTY, HEAT_OF_TRANSPORT_DIMS, 0.);
    let b = scalar(HEAT_OF_TRANSPORT_PROPERTY, HEAT_OF_TRANSPORT_DIMS, 0.);
    let context = point(330., 4e6);
    let before = original([0., 0.]);
    let saved = before.clone();
    let model = prepared(&before[&10], &binary, &k, [&a, &b], &context);
    let link = RealGasThermodiffusionLink {
        id: 1,
        parcels: [10, 20],
        model: &model,
        geometry: geometry(),
        solve: thermal_controls(),
    };
    let closed = RealGasThermodiffusionLink {
        geometry: GasDiffusionGeometry {
            area_m2: 0.,
            ..geometry()
        },
        ..link
    };
    let idle = thermal::advance_candidate(
        &before,
        &[closed],
        2.,
        &endpoints(&before),
        physical(),
        temporal(),
        || false,
    )
    .unwrap();
    assert_eq!(idle.parcels(), &before);
    assert!(idle.report().heat.is_empty());
    assert_eq!(idle.report().face_evaluations, 0);
    let mut controls = temporal();
    controls.absolute_additional_heat_tolerance_j = 1e-3;
    let mut polls = 0;
    thermal::advance_candidate(
        &before,
        &[link],
        2.,
        &endpoints(&before),
        physical(),
        controls,
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut seen = 0;
    assert!(matches!(
        thermal::advance_candidate(
            &before,
            &[link],
            2.,
            &endpoints(&before),
            physical(),
            controls,
            || {
                seen += 1;
                seen == polls
            }
        ),
        Err(DiffusionNetworkAdvanceError::Cancelled)
    ));
    assert_eq!(before, saved);
}
