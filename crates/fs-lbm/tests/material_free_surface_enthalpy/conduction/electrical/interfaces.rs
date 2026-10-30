//! Pinned interface data drives the original current, heat and native liquid owners.
//! All cards below are synthetic numerical fixtures, not measured liquid data.

use super::*;
use fs_lbm::conjugate::ElectricalContact;
use fs_lbm::material::free_surface::enthalpy::conduction::electrical::{
    FreeSurfaceMagneticInput,
    interfaces::{
        FreeSurfaceElectricalConnectivity as Connectivity,
        FreeSurfaceElectricalInterfaceBinding as Binding,
        FreeSurfaceElectricalInterfaceUpdate as Update,
    },
};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as CONTACT_DIMS,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as CONTACT_PROPERTY, resolve_electrical_interface_state,
};
use std::sync::Arc;

const AREA: f64 = 1e-6;
const SUBSTEP: f64 = 0.00025;
const CAPACITY: f64 = 0.002;

fn frames() -> [String; 2] {
    ["specimen".into(), "left-surface".into()]
}

fn interface_card(
    material: &MaterialCard,
    resistance_ohm: f64,
    slope_per_k: f64,
    upper_k: f64,
) -> Arc<InterfaceSystemCard> {
    let mut claims = ClaimSet::new();
    let mut claim = PropertyClaim {
        key: PropertyKey::new(CONTACT_PROPERTY, CONTACT_DIMS),
        value: PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (400.0, AREA * resistance_ohm),
                (500.0, AREA * resistance_ohm * (1.0 + 100.0 * slope_per_k)),
            ],
            dims: CONTACT_DIMS,
        },
        validity: ValidityDomain::unconstrained()
            .with("T", 400.0, upper_k)
            .with("pressure", 100000.0, 100000.0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic ordered liquid electrical interface fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    };
    claims.insert_claim(claim.clone()).unwrap();
    if upper_k < 500.0 {
        claim.validity = ValidityDomain::unconstrained()
            .with("T", upper_k + 1e-6, 500.0)
            .with("pressure", 100000.0, 100000.0);
        claim.provenance.source = "new fallback must not replace the pinned source".into();
        claims.insert_claim(claim).unwrap();
    }
    let [first, second] = frames();
    Arc::new(
        InterfaceSystemCard::assemble(
            SurfaceSpec {
                material: material.id().clone(),
                texture_frame: first,
            },
            SurfaceSpec {
                material: material.id().clone(),
                texture_frame: second,
            },
            SystemContext {
                medium: "conducting film".into(),
                third_body: None,
                environment: "fixture atmosphere".into(),
                history: "specified initial surfaces".into(),
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}

fn binding(card: &Arc<InterfaceSystemCard>, actual_frames: [String; 2]) -> Binding {
    let reference = resolve_electrical_interface_state(
        card,
        &point(400.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    Binding::try_new(Arc::clone(card), reference, actual_frames).unwrap()
}

fn authored(cells: [usize; 2]) -> ElectricalContact {
    // These old fixed coefficients must be REPLACED by the selected source.
    ElectricalContact {
        cells,
        resistance_ohm_m2: 0.1,
        heat_fraction_to_first: 0.9,
    }
}

#[test]
fn interface_source_temperature_changes_ordered_heat_current_and_actual_liquid_motion() {
    let run = |scale, magnetic: bool| {
        let basis = ElectricalPropertyBasis::Resistivity;
        let card = electrical_card(basis, 1.0, false, false, 500.0);
        let contact_card = interface_card(&card, scale, 0.5, 500.0);
        with_unbound(&card, basis, false, |liquid, gate| {
            let mut input = current(4.0);
            // The first surface is the specimen, even though its index is higher.
            input.contacts = vec![authored([1, 0])];
            let mut liquid = liquid
                .with_electrical_core(input, gate)
                .unwrap()
                .with_material_electrical_interfaces(
                    vec![binding(&contact_card, frames())],
                    |_, _, _, sides| {
                        Ok(Update {
                            query: point(sides[0].temperature_k),
                            connectivity: Connectivity::Closed {
                                heat_fraction_to_first: 0.25,
                            },
                        })
                    },
                    gate,
                )
                .unwrap();
            let bz = 0.001;
            if magnetic {
                liquid = liquid
                    .with_imposed_magnetic_field(FreeSurfaceMagneticInput {
                        magnetic_flux_density_t: vec![[0.0, 0.0, bz]; 3],
                        maximum_velocity_change_m_s: 0.03,
                    })
                    .unwrap();
            }
            let before = liquid.total_enthalpy_j().unwrap();
            let report = liquid
                .step_with_conduction(&[], &[], controls(), gate)
                .unwrap();
            let electric = report.electrical.as_ref().unwrap();
            // Independent explicit circuit/storage recurrence: the one bulk
            // cell contributes 1 ohm. R_contact(T) is in series, not a replacement
            // for bulk resistivity. Both real liquid cells store their heat.
            let (mut specimen, mut left, mut last_t) = (400.0, 400.0, 400.0);
            let (mut heat, mut impulse) = (0.0, 0.0);
            for _ in 0..4 {
                last_t = specimen;
                let rc = scale * (1.0 + 0.5 * (specimen - 400.0));
                let amps = 4.0 / (1.0 + rc);
                let bulk = amps * amps;
                let contact = bulk * rc;
                specimen += SUBSTEP * (bulk + 0.25 * contact) / CAPACITY;
                left += SUBSTEP * 0.75 * contact / CAPACITY;
                heat += SUBSTEP * (bulk + contact);
                impulse += SUBSTEP * amps * 0.001 * bz;
            }
            assert_eq!(electric.substeps, 4);
            assert_eq!(electric.current_solves, 8);
            near(electric.last_applied.temperature_k[1], last_t, 2e-12);
            near(electric.thermal_endpoint.temperature_k[1], specimen, 2e-12);
            near(electric.thermal_endpoint.temperature_k[0], left, 2e-12);
            let state = &electric.thermal_endpoint;
            let rc = scale * (1.0 + 0.5 * (specimen - 400.0));
            let amps = 4.0 / (1.0 + rc);
            near(state.solution.electrodes[1].current_a, amps, 2e-11);
            near(state.solution.joule_w[0], 0.75 * amps * amps * rc, 2e-11);
            near(
                state.solution.joule_w[1],
                amps * amps * (1.0 + 0.25 * rc),
                2e-11,
            );
            let epoch = state.interfaces.as_ref().unwrap();
            assert_eq!(epoch.time_s, state.time_s);
            let sample = &epoch.contacts[0];
            assert_eq!(sample.authored, authored([1, 0]));
            let applied = sample.applied().unwrap();
            assert_eq!(applied.cells, [1, 0]);
            assert_eq!(applied.heat_fraction_to_first, 0.25);
            near(applied.resistance_ohm_m2, AREA * rc, 2e-11);
            let source = sample.response.material();
            assert_eq!(
                source.resolved().card_identity(),
                contact_card.content_hash()
            );
            near(
                source.resolved().state_coordinate("T").unwrap(),
                specimen,
                2e-12,
            );
            assert_eq!(source.surfaces()[0].texture_frame, frames()[0]);
            near(electric.joule_heat_j, heat, 2e-11);
            near(electric.electrode_energy_j, heat, 2e-11);
            near(report.flow.external_heat_j, heat, 2e-11);
            assert!((liquid.total_enthalpy_j().unwrap() - before - heat).abs() < 1e-9);
            assert!(electric.energy_defect_j.abs() < 1e-12);
            if magnetic {
                near(electric.lorentz_impulse_n_s[1][1], impulse, 2e-11);
                assert!(electric.sampled_lorentz_work_j.abs() < 1e-16);
                assert!(
                    electric
                        .last_applied
                        .magnetic
                        .as_ref()
                        .unwrap()
                        .velocity_m_s[1]
                        .iter()
                        .all(|u| u.abs() < 1e-15)
                );
            } else {
                assert!(electric.lorentz_impulse_n_s.is_empty());
                assert!(electric.last_applied.magnetic.is_none());
            }
            let velocity = liquid
                .flow()
                .lattice()
                .sample_at(CORE[0], CORE[1], CORE[2])
                .unwrap()
                .velocity[1];
            (velocity, liquid.flow().lattice().grid().state_bits(), heat)
        })
    };
    let stationary_current = run(1.0, false);
    let original = run(1.0, true);
    let more_resistive = run(2.0, true);
    assert!(original.0 > stationary_current.0);
    assert!(more_resistive.0 > stationary_current.0);
    assert!(original.0 > more_resistive.0);
    assert_ne!(original.1, more_resistive.1);
    assert!(original.2 > more_resistive.2);
    near(original.2, stationary_current.2, 2e-11);
}

fn loop_velocities(liquid: &MaterialFreeSurfaceEnthalpy3<'_, '_>) -> Vec<[f64; 3]> {
    (0..4)
        .map(|cell| {
            liquid
                .flow()
                .lattice()
                .sample_at(1 + cell % 2, 3 + cell / 2, 3)
                .unwrap()
                .velocity
        })
        .collect()
}

#[test]
fn source_contact_opens_the_isolated_loop_during_the_real_thermal_step_then_recloses() {
    let basis = ElectricalPropertyBasis::Conductivity;
    let card = electrical_card(basis, 1.0, false, false, 500.0);
    let contact_card = interface_card(&card, 4.0, 0.0, 500.0);
    with_initial_flow(
        &card,
        basis,
        false,
        Some([0.0, 0.02, 0.0]),
        [0.0; 3],
        |liquid, gate| {
            let bz = 50.0;
            let input = FreeSurfaceElectricalInput {
                origin: [1, 3, 3],
                dimensions: [2, 2, 1],
                setup: ElectricSetup::default(),
                contacts: vec![authored([0, 1])],
                maximum_cells: 4,
                maximum_solves: 32,
                power_tolerance_w: 1e-12,
            };
            let bindings = vec![binding(&contact_card, frames())];
            let mut liquid = liquid
                .with_isolated_electrical_core(input, 1e-12, gate)
                .unwrap()
                .with_imposed_magnetic_field(FreeSurfaceMagneticInput {
                    magnetic_flux_density_t: vec![
                        [0.0, 0.0, bz],
                        [0.0, 0.0, bz],
                        [0.0, 0.0, -bz],
                        [0.0, 0.0, -bz],
                    ],
                    maximum_velocity_change_m_s: 0.03,
                })
                .unwrap()
                .with_material_electrical_interfaces(
                    bindings.clone(),
                    |_, _, _, sides| {
                        Ok(Update {
                            query: point(sides[0].temperature_k),
                            connectivity: if sides[0].temperature_k > 400.0 + 1e-9 {
                                Connectivity::Open
                            } else {
                                Connectivity::Closed {
                                    heat_fraction_to_first: 0.5,
                                }
                            },
                        })
                    },
                    gate,
                )
                .unwrap();
            let initial_velocity = loop_velocities(&liquid);
            let before = liquid.total_enthalpy_j().unwrap();
            let first = liquid
                .step_with_conduction(&[], &[], controls(), gate)
                .unwrap();
            let electric = first.electrical.as_ref().unwrap();
            // Four 1-ohm bulk edges plus a 4-ohm contact. It opens after the FIRST
            // heat interval; a tree has no circulating induced current. The source
            // remains present, but the opened face has no current, heat or force.
            let amps = (2.0 * 0.02 * bz * 0.001) / 8.0;
            let heat = SUBSTEP * 8.0 * amps * amps;
            let impulse = -SUBSTEP * amps * 0.001 * bz / 2.0;
            near(electric.joule_heat_j, heat, 2e-7);
            near(electric.sampled_lorentz_work_j, -heat, 2e-7);
            assert_eq!(electric.electrode_energy_j, 0.0);
            assert!((liquid.total_enthalpy_j().unwrap() - before - heat).abs() < 2e-12);
            for cell in 0..4 {
                near(electric.lorentz_impulse_n_s[cell][1], impulse, 2e-7);
            }
            for state in [&electric.last_applied, &electric.thermal_endpoint] {
                assert!(state.solution.electrodes.is_empty());
                assert!(state.solution.joule_w.iter().all(|q| q.abs() < 1e-20));
                assert_eq!(
                    state.magnetic.as_ref().unwrap().velocity_m_s,
                    initial_velocity
                );
                let sample = &state.interfaces.as_ref().unwrap().contacts[0];
                assert!(sample.response.is_open());
                assert!(sample.applied().is_none());
                assert_eq!(
                    sample.response.material().resolved().card_identity(),
                    contact_card.content_hash()
                );
            }
            let accepted_velocity = loop_velocities(&liquid);
            assert!(
                accepted_velocity
                    .iter()
                    .zip(&initial_velocity)
                    .all(|(actual, initial)| actual[1] < initial[1])
            );
            let second = liquid
                .step_with_conduction(&[], &[], controls(), gate)
                .unwrap();
            let second = second.electrical.as_ref().unwrap();
            assert!(second.joule_heat_j.abs() < 1e-20);
            assert!(
                second
                    .lorentz_impulse_n_s
                    .iter()
                    .flatten()
                    .all(|f| f.abs() < 1e-15)
            );
            assert_eq!(
                second.last_applied.magnetic.as_ref().unwrap().velocity_m_s,
                accepted_velocity
            );
            // Closing restores the ORIGINAL source on the actual accepted velocity,
            // rather than resetting flow or reusing its old electrical force frame.
            let before_reclose = loop_velocities(&liquid);
            let mut liquid = liquid
                .with_material_electrical_interfaces(
                    bindings,
                    |_, _, _, sides| {
                        Ok(Update {
                            query: point(sides[0].temperature_k),
                            connectivity: Connectivity::Closed {
                                heat_fraction_to_first: 0.5,
                            },
                        })
                    },
                    gate,
                )
                .unwrap();
            let third = liquid
                .step_with_conduction(&[], &[], controls(), gate)
                .unwrap();
            let third = third.electrical.as_ref().unwrap();
            assert!(third.joule_heat_j > 0.0);
            assert!(third.sampled_lorentz_work_j < 0.0);
            assert!(
                !third.last_applied.interfaces.as_ref().unwrap().contacts[0]
                    .response
                    .is_open()
            );
            assert_eq!(
                third.last_applied.magnetic.as_ref().unwrap().velocity_m_s,
                before_reclose
            );
            assert!(third.energy_defect_j.abs() < 1e-15);
        },
    );
}

#[test]
fn ordered_surface_mismatch_and_late_pinned_source_failure_roll_back_every_participant() {
    for invalid_texture in [false, true] {
        let basis = ElectricalPropertyBasis::Resistivity;
        let card = electrical_card(basis, 1.0, false, false, 500.0);
        let contact_card = interface_card(
            &card,
            1.0,
            0.0,
            if invalid_texture { 500.0 } else { 400.01 },
        );
        if !invalid_texture {
            let original = resolve_electrical_interface_state(
                &contact_card,
                &point(400.0),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            let fallback = resolve_electrical_interface_state(
                &contact_card,
                &point(401.0),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap();
            assert_ne!(
                original.resolved().properties()[0]
                    .answer()
                    .receipt
                    .selected,
                fallback.resolved().properties()[0]
                    .answer()
                    .receipt
                    .selected,
            );
            assert!(
                original
                    .at_query_point(&contact_card, &point(401.0))
                    .is_err()
            );
        }
        let solid = body::wall_card(0.5, 500.0);
        let chart = body::wall_chart(&solid);
        let context = point(450.0);
        with_unbound(&card, basis, false, |liquid, gate| {
            let mut input = current(4.0);
            input.contacts = vec![authored([1, 0])];
            let mut actual_frames = frames();
            if invalid_texture {
                actual_frames.swap(0, 1);
            }
            let mut liquid = liquid
                .with_time_origin(7.0)
                .unwrap()
                .with_electrical_core(input, gate)
                .unwrap()
                .with_material_electrical_interfaces(
                    vec![binding(&contact_card, actual_frames)],
                    |_, _, _, sides| {
                        Ok(Update {
                            query: point(sides[0].temperature_k),
                            connectivity: Connectivity::Closed {
                                heat_fraction_to_first: 0.25,
                            },
                        })
                    },
                    gate,
                )
                .unwrap();
            let mut atmosphere = gas(100.0);
            let mut bodies = [body::wall(&solid, &chart, &context, 450.0)];
            let bits = liquid.flow().lattice().grid().state_bits();
            let energy = liquid.total_enthalpy_j().unwrap();
            let material = liquid.material_at(CORE).unwrap().clone();
            let gas_before = atmosphere.state().clone();
            let gas_point = atmosphere.point().clone();
            let wall_before = bodies[0].specific_enthalpy_j_kg();
            let wall_source = bodies[0].material().clone();
            for _ in 0..2 {
                let error = liquid
                    .step_with_finite_gas_environment(
                        &[],
                        &[],
                        controls(),
                        &mut bodies,
                        &[body::contact()],
                        &[],
                        &mut atmosphere,
                        &[gas_face()],
                        gas_controls(),
                        false,
                        gate,
                    )
                    .unwrap_err();
                if invalid_texture {
                    assert!(matches!(
                        error,
                        FreeSurfaceEnthalpyError::Invalid(
                            "electrical interface source and actual ordered texture frames differ"
                        )
                    ));
                } else {
                    assert!(matches!(error, FreeSurfaceEnthalpyError::Material(_)));
                }
                assert_eq!(liquid.flow().lattice().grid().state_bits(), bits);
                assert_eq!(liquid.total_enthalpy_j().unwrap(), energy);
                assert_eq!(liquid.material_at(CORE).unwrap(), &material);
                assert_eq!(atmosphere.state(), &gas_before);
                assert_eq!(atmosphere.point(), &gas_point);
                assert_eq!(bodies[0].specific_enthalpy_j_kg(), wall_before);
                assert_eq!(bodies[0].material(), &wall_source);
                assert_eq!(liquid.accepted_steps(), 0);
                assert_eq!(liquid.physical_time_s(), 7.0);
            }
        });
    }
}
