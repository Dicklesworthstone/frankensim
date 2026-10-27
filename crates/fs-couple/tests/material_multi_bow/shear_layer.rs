//! Manufactured materials exercise the real instrument pressure path.
//! These values identify no rosin, coating, bow or measured instrument.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::ModalFrictionRegime;
use fs_material::state_point::ResolvedShearLayer;

fn layer(g: f64, yield_stress: f64, viscosity: f64, maximum_omega: f64) -> ResolvedShearLayer {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("shear_modulus", Pressure::DIMS, g),
        ("shear_yield_stress", Pressure::DIMS, yield_stress),
        ("plastic_shear_viscosity", DynViscosity::DIMS, viscosity),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained()
                    .with("T", 280.0, 310.0)
                    .with("omega", 0.0, maximum_omega),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured thin simple-shear layer".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
        requirements.push(
            ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap(),
        );
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured-layer".into(),
            phase: "solid".into(),
            process: "fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let material = resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 0.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    // 20 um film over a 2 mm^2 patch; explicit 5% elastic-strain ceiling.
    ResolvedShearLayer::new(material, 2e-6, 20e-6, 0.05).unwrap()
}

fn coated(cx: &Cx<'_>, gate: &CancelGate, g: f64, viscosity: f64) -> MultiBowStringBank {
    reacting_surface(joint(cx, gate, 7800.0), gate)
        .with_bow_shear_layers(vec![(0, layer(g, 2000.0, viscosity, 1e6))], gate)
        .unwrap()
}

#[test]
fn layer_material_swap_changes_shared_instrument_pressure_and_keeps_physical_energy() {
    with_cx(|cx, gate| {
        let mut a = coated(cx, gate, 8e6, 10.0);
        let mut b = coated(cx, gate, 16e6, 40.0);
        assert_eq!(
            a.system().acceleration_pressure_weights(),
            b.system().acceleration_pressure_weights()
        );
        for (left, right) in a.system().components().iter().zip(b.system().components()) {
            assert_eq!(left.modes(), right.modes());
        }
        assert!(a.source().bank().material_cavity().is_some());
        assert_eq!(a.source().body_surfaces().len(), 2);
        assert_eq!(a.bow_shear_layer(0).unwrap().temperature_kelvin(), 293.15);
        assert!(a.bow_shear_layer(1).is_none());
        let zero = vec![0.0; a.system().mode_count()];
        let (mut pressure_difference, mut traction_difference, mut heat) = (0.0, 0.0, 0.0);
        for _ in 0..256 {
            let drives = [drive(0.001, 0.1), drive(0.0, 0.0)];
            let fa = a.step(&zero, &drives).unwrap().clone();
            let fb = b.step(&zero, &drives).unwrap().clone();
            for frame in [&fa, &fb] {
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
                let port = frame.ports[0];
                let layer = port.shear_layer.unwrap();
                assert!(matches!(port.regime, ModalFrictionRegime::Shearing { .. }));
                assert_eq!(port.dissipation_j, layer.dissipated_j);
                assert!(layer.dissipated_j >= 0.0 && layer.stored_energy_j >= 0.0);
                assert_eq!(frame.shear_layer_energy_j, layer.stored_energy_j);
                assert!(
                    (layer.stored_change_j + layer.dissipated_j
                        - port.traction_n * port.slip_distance_m)
                        .abs()
                        <= frame.energy_tolerance_j
                );
                assert!(frame.ports[1].shear_layer.is_none());
            }
            heat += fa.ports[0].dissipation_j;
            pressure_difference += (fa.observer_pressure_pa - fb.observer_pressure_pa).abs();
            traction_difference += (fa.ports[0].traction_n - fb.ports[0].traction_n).abs();
        }
        assert!(heat > 0.0 && traction_difference > 1e-10 && pressure_difference > 1e-12);
        assert!(
            a.system().components()[1]
                .states()
                .iter()
                .any(|s| *s != ModalAcousticState::default()),
            "the undriven material string lost its shared-body reaction"
        );
        let f = a.system().last_frame().unwrap();
        assert!(
            (a.system().total_energy_j().unwrap() - f.network_energy_j - f.shear_layer_energy_j)
                .abs()
                < 1e-12
        );
    });
}

#[test]
fn elastic_unloading_returns_storage_and_opening_refuses_without_erasing_history() {
    with_cx(|cx, gate| {
        // Yield lies beyond the admitted small elastic strain, so this loading
        // cycle is reversible in the layer even though the bodies can dissipate.
        let mut bank = reacting_surface(joint(cx, gate, 7800.0), gate)
            .with_bow_shear_layers(vec![(0, layer(8e6, 2e6, 10.0, 1e6))], gate)
            .unwrap();
        let zero = vec![0.0; bank.system().mode_count()];
        let mut returned = false;
        for sample in 0..64 {
            let speed = if sample < 16 { 0.0003 } else { -0.0003 };
            let frame = bank
                .step(&zero, &[drive(speed, 0.1), drive(0.0, 0.0)])
                .unwrap();
            let port = frame.ports[0];
            let layer = port.shear_layer.unwrap();
            assert_eq!(
                port.regime,
                ModalFrictionRegime::Shearing {
                    plastic_flow: false
                }
            );
            assert_eq!(layer.dissipated_j, 0.0);
            returned |= layer.stored_change_j < 0.0 && port.traction_n * port.slip_distance_m < 0.0;
        }
        assert!(
            returned,
            "elastic unloading was incorrectly forced to dissipate all work"
        );
        let old_states = states(&bank);
        let old_layer = bank.system().material_shear_state(0).unwrap().clone();
        let old_frame = bank.system().last_frame().unwrap().clone();
        let old_energy = bank.system().total_energy_j().unwrap();
        assert!(
            bank.step(&zero, &[drive(0.0, 0.0), drive(0.0, 0.0)])
                .is_err()
        );
        assert_eq!(states(&bank), old_states);
        assert_eq!(bank.system().material_shear_state(0), Some(&old_layer));
        assert_eq!(bank.system().last_frame(), Some(&old_frame));
        assert_eq!(bank.system().total_energy_j().unwrap(), old_energy);
        assert_eq!(bank.system().samples_rendered(), old_frame.sample);
        // Retrying the unchanged retained-contact boundary is permitted.
        bank.step(&zero, &[drive(0.0, 0.1), drive(0.0, 0.0)])
            .unwrap();
    });
}

#[test]
fn material_layer_and_series_hair_share_one_compliance_and_energy_transaction() {
    with_cx(|cx, gate| {
        let (_, mut system) = coated(cx, gate, 8e6, 10.0).into_parts();
        system = system
            .with_series_compliant_drives(vec![(0, 2000.0, 0.4)], gate)
            .unwrap();
        let zero = vec![0.0; system.mode_count()];
        let (mut work, mut losses) = (0.0, 0.0);
        for _ in 0..192 {
            let f = system
                .step(&zero, &[drive(0.001, 0.1), drive(0.0, 0.0)])
                .unwrap()
                .clone();
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
            assert!(f.ports[0].shear_layer.is_some());
            work += f.drive_work_j;
            losses += f.network_dissipation_j + f.friction_dissipation_j;
            assert!((system.total_energy_j().unwrap() + losses - work).abs() < 2e-10);
        }
        let f = system.last_frame().unwrap();
        assert!(
            system.total_energy_j().unwrap() > f.network_energy_j + f.shear_layer_energy_j,
            "series hair elastic energy disappeared when material shear was bound"
        );
    });
}

#[test]
fn material_layer_series_drive_and_string_strain_share_storage_and_transactional_history() {
    use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::NonlinearModalConfig;
    with_cx(|cx, gate| {
        let build = |tight: bool| {
            let bank = base(cx, gate, 7800.0).with_geometric_stretching(0.2,
                NonlinearModalConfig { max_iterations: if tight { 1 } else { 64 },
                    maximum_generalized_force: 1e5,
                    force_absolute_tolerance: if tight { 1e-30 } else { 1e-11 },
                    force_relative_tolerance: if tight { 1e-30 } else { 1e-10 } }, gate)
                .unwrap().with_additional_bows(vec![extra()], joint_limits(), gate).unwrap();
            let (_, system) = reacting_surface(bank, gate)
                // Zero yield in this manufactured film guarantees that the
                // small driving speed exercises plastic loss as well as storage.
                .with_bow_shear_layers(vec![(0, layer(8e6, 0.0, 10.0, 1e6))], gate)
                .unwrap().into_parts();
            system.with_series_compliant_drives(vec![(0, 2000.0, 0.4)], gate).unwrap()
        };
        let mut system = build(false);
        let mut replay = build(false);
        let zero = vec![0.0; system.mode_count()];
        let (mut work, mut losses, mut strain, mut plastic, mut hair_loss) = (0.0, 0.0, 0.0_f64, 0.0, 0.0);
        for sample in 0..192 {
            if sample == 96 {
                let frame = system.last_frame().unwrap().clone();
                let history = system.material_shear_state(0).unwrap().clone();
                let before = system.components().iter().map(|m| m.states().to_vec()).collect::<Vec<_>>();
                let energy = system.total_energy_j().unwrap();
                assert!(system.step(&zero, &[drive(0.0, 0.0); 2]).is_err());
                assert_eq!(system.last_frame(), Some(&frame));
                assert_eq!(system.material_shear_state(0), Some(&history));
                assert_eq!(system.components().iter().map(|m| m.states().to_vec()).collect::<Vec<_>>(), before);
                assert_eq!(system.total_energy_j().unwrap(), energy);
                assert_eq!(system.samples_rendered(), frame.sample);
            }
            let drives = [drive(if sample < 128 { 0.001 } else { -0.001 }, 0.1), drive(0.0, 0.0)];
            let f = system.step(&zero, &drives).unwrap().clone();
            assert_eq!(&f, replay.step(&zero, &drives).unwrap());
            let nonlinear = f.nonlinear.unwrap();
            assert!(nonlinear.reaction_residual <= nonlinear.reaction_tolerance);
            assert!(f.energy_residual_j.abs() <= f.energy_tolerance_j);
            let layer = f.ports[0].shear_layer.unwrap();
            assert!((layer.stored_change_j + layer.dissipated_j
                - f.ports[0].traction_n * f.ports[0].slip_distance_m).abs() <= f.energy_tolerance_j);
            work += f.external_work_j + f.drive_work_j;
            losses += f.network_dissipation_j + f.friction_dissipation_j;
            strain = strain.max(nonlinear.nonlinear_energy_j);
            plastic += layer.dissipated_j;
            hair_loss += f.ports[0].series_drive_dissipation_j;
            assert!((system.total_energy_j().unwrap() + losses - work).abs() < 2e-10);
        }
        let frame = system.last_frame().unwrap();
        assert!(strain > 0.0 && plastic > 0.0 && hair_loss > 0.0);
        assert!(system.total_energy_j().unwrap() > frame.network_energy_j + frame.shear_layer_energy_j);

        // Exhaust the OUTER constitutive solve after it staged actual film and
        // hair histories. None of those candidates may become accepted state.
        let mut refused = build(true);
        let mut initial = refused.components()[0].states().to_vec();
        initial[0].displacement_m_sqrt_kg = 1e-7;
        refused.restore_component_states(0, &initial).unwrap();
        let states = refused.components().iter().map(|m| m.states().to_vec()).collect::<Vec<_>>();
        let history = refused.material_shear_state(0).unwrap().clone();
        let energy = refused.total_energy_j().unwrap();
        for _ in 0..2 {
            assert!(refused.step(&zero, &[drive(0.001, 0.1), drive(0.0, 0.0)]).is_err());
            assert_eq!(refused.components().iter().map(|m| m.states().to_vec()).collect::<Vec<_>>(), states);
            assert_eq!(refused.material_shear_state(0), Some(&history));
            assert_eq!(refused.total_energy_j().unwrap(), energy);
            assert_eq!(refused.samples_rendered(), 0);
            assert!(refused.last_frame().is_none());
        }
    });
}

fn retained_schedule() -> GestureSchedule {
    GestureSchedule::try_new(
        700,
        (0..2)
            .map(|i| {
                let station = if i == 0 { 0.2 } else { 0.3 };
                GestureTrack {
                    id: format!("bow-{i}"),
                    target: GestureTarget::BowStroke { string: i },
                    initial: value(0.0004, 0.1, station),
                    events: vec![GestureEvent {
                        time_s: 0.002,
                        transition_s: 0.003,
                        value: value(-0.0002, 0.08, station),
                    }],
                }
            })
            .collect(),
    )
    .unwrap()
}

#[test]
fn layer_history_and_pressure_are_callback_invariant_and_invalid_release_is_admitted_early() {
    with_cx(|cx, gate| {
        let build = || {
            coated(cx, gate, 8e6, 10.0)
                .into_scheduled(retained_schedule(), tracks(), performance_config())
                .unwrap()
        };
        let mut a = build();
        let mut b = build();
        let mut first = vec![0.0; 513];
        let mut second = first.clone();
        for block in first.chunks_mut(1) {
            a.block(block).unwrap();
        }
        for block in second.chunks_mut(127) {
            b.block(block).unwrap();
        }
        assert_eq!(first, second);
        assert_eq!(states(a.bank()), states(b.bank()));
        assert_eq!(
            a.bank().system().material_shear_state(0),
            b.bank().system().material_shear_state(0)
        );
        assert!(first.iter().any(|p| p.abs() > 1e-12));
        assert!(
            coated(cx, gate, 8e6, 10.0)
                .into_scheduled(schedule(), tracks(), performance_config())
                .is_err()
        );
        assert!(
            joint(cx, gate, 7800.0)
                .with_bow_shear_layers(vec![(0, layer(8e6, 2000.0, 10.0, 1000.0))], gate)
                .is_err()
        );
        let (_, mut system) = coated(cx, gate, 8e6, 10.0).into_parts();
        let zero = vec![0.0; system.mode_count()];
        assert!(system.initialize_released_equilibrium(&zero, gate).is_err());
        assert_eq!(system.samples_rendered(), 0);
    });
}

#[test]
fn solved_normal_contact_reaches_layer_and_opening_retains_the_whole_instrument() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::multiple::{
        contact_load::BowContactLoad,
        performance::motion::ObstacleMotionLimits,
        stops::{MaterialStringStop, StringStopSide},
    };
    use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{
        TransverseDirection, TransverseDirections,
    };
    use fs_couple::render::schedule::force::coupled::contact::ModalContactConfig;
    use fs_couple::render::schedule::force::coupled::contact::friction::multiple::obstacles::FrictionObstacleConfig;
    use fs_dcontact::Obstacle;
    with_cx(|cx, gate| {
        let reference = base(cx, gate, 7800.0);
        let body = reference.source().material_plate().unwrap().clone();
        let cavity = reference.source().material_cavity().unwrap().clone();
        let strings = (0..2)
            .map(|_| {
                source(7800.0)
                    .with_transverse_polarizations(
                        TransverseDirections::try_new([1.0, 0.0], [0.6, 0.8]).unwrap(),
                    )
                    .unwrap()
            })
            .collect();
        let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
        let stop = MaterialStringStop {
            location: options().bow,
            side: StringStopSide::Along(TransverseDirection::try_new([0.0, -1.0]).unwrap()),
            body: None,
            law: Obstacle::new(
                vec![-1.0],
                1,
                1,
                vec![-1e-5],
                vec![1.0],
                5000.0,
                1.0,
                "manufactured retained normal contact".into(),
            )
            .unwrap(),
            limits: ModalContactConfig {
                max_iterations: 96,
                maximum_force_n: 10.0,
                maximum_penetration_m: 0.001,
                force_absolute_tolerance_n: 1e-11,
                force_relative_tolerance: 1e-9,
            },
        };
        let mut bank = ReciprocalStringBank::from_material_plate_with_cavity(
            cx,
            strings,
            gas,
            body,
            cavity,
            law(),
            options(),
            gate,
        )
        .unwrap()
        .with_additional_bows(vec![extra()], joint_limits(), gate)
        .unwrap()
        .with_string_stops(
            vec![stop],
            FrictionObstacleConfig {
                max_contacts: 4,
                max_setup_terms: 100000,
            },
            gate,
        )
        .unwrap()
        .with_bow_contact_loads(vec![BowContactLoad { bow: 0, stop: 0 }], gate)
        .unwrap()
        .with_bow_shear_layers(vec![(0, layer(8e6, 2000.0, 10.0, 1e6))], gate)
        .unwrap()
        .with_stop_motion(
            ObstacleMotionLimits {
                maximum_abs_offset_m: 0.001,
                maximum_speed_m_s: 100.0,
            },
            gate,
        )
        .unwrap();
        let zero = vec![0.0; bank.system().mode_count()];
        // Authored normal load is zero because the original normal owner
        // supplies this layer's actual simultaneous compression.
        let drives = [drive(0.0001, 0.0), drive(0.0, 0.0)];
        for _ in 0..4 {
            let frame = bank.step(&zero, &drives).unwrap().clone();
            let normal = bank.system().last_obstacle_frames().unwrap()[0].normal_force_n;
            assert!(normal > 0.0);
            assert_eq!(frame.ports[0].normal_force_n, normal);
            assert!(frame.ports[0].shear_layer.is_some());
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
        let old_states = states(&bank);
        let old_layer = bank.system().material_shear_state(0).unwrap().clone();
        let old_frame = bank.system().last_frame().unwrap().clone();
        let old_normal = bank.system().last_obstacle_frames().unwrap().to_vec();
        assert!(
            bank.step_with_stop_offsets(&zero, &drives, &[0.001])
                .is_err()
        );
        assert_eq!(states(&bank), old_states);
        assert_eq!(bank.system().material_shear_state(0), Some(&old_layer));
        assert_eq!(bank.system().last_frame(), Some(&old_frame));
        assert_eq!(
            bank.system().last_obstacle_frames(),
            Some(old_normal.as_slice())
        );
        assert_eq!(bank.system().samples_rendered(), old_frame.sample);
        bank.step(&zero, &drives).unwrap();
    });
}
