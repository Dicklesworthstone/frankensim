//! Signed layer dashpots and Maxwell histories share the original native state.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::LongitudinalMaterialLoss;

const BRIDGE: [f64; 2] = [0.6, 0.8];

fn viscous_source(cx: &Cx<'_>, gate: &CancelGate, limit: f64, eta: f64) -> NonlinearMaterialString {
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: CORE_R,
                material: material(CORE_E, 0.0, 1e12, 0.0),
            },
            ConcentricMaterialLayer {
                outer_radius_m: RADIUS,
                material: material(SLEEVE_E, DELTA, limit, eta),
            },
        ],
    )
    .unwrap()
    .with_kelvin_voigt_bending(ETA)
    .unwrap();
    polar(
        cx,
        gate,
        MaterialStringSource::Concentric {
            template: template(3),
            section,
            prestress: StringPrestress::FixedTension(20.0),
        },
        [0.0; 2],
        BRIDGE,
    )
    .with_tensile_limit(LIMIT)
    .unwrap()
}

fn lossy_image(
    source: NonlinearMaterialString,
    support: TranslatingSupport,
    modes: &[usize],
    axial: Vec<ModalAcousticState>,
    finite_strain: bool,
    hereditary: bool,
    gate: &CancelGate,
) -> SupportedMaterialStrings {
    let property = YOUNG_VISCOSITY_PROPERTY.into();
    let loss = if finite_strain {
        LongitudinalMaterialLoss::FiniteStrainKelvinVoigt {
            viscosity_property: property,
        }
    } else {
        LongitudinalMaterialLoss::KelvinVoigt {
            viscosity_property: property,
        }
    };
    let image = SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves_and_loss(
        vec![source],
        support,
        vec![attachment(modes, axial, BRIDGE)],
        loss,
        300.0,
        gate,
    )
    .unwrap();
    if hereditary {
        image
            .with_material_longitudinal_relaxation(
                selection(),
                InitialMemory::Unrelaxed,
                modes.len() + 1,
                gate,
            )
            .unwrap()
    } else {
        image
    }
}

#[test]
fn g1_signed_layer_viscosity_and_memory_use_actual_sparse_modes_and_support_velocity() {
    with_cx(|cx, gate| {
        let eta = 2e6;
        let displacement = 2e-5;
        let amplitudes = [2e-6, 1e-6];
        let modes = [1, 3];
        let wave = core::f64::consts::PI / LENGTH;
        let area = core::f64::consts::PI * (RADIUS * RADIUS - CORE_R * CORE_R);
        for hereditary in [false, true] {
            let active = if hereditary { DELTA } else { 0.0 };
            for finite_strain in [false, true] {
                for sign in [-1.0, 1.0] {
                    let speed = sign * 0.02;
                    // With any unrelaxed Maxwell memories these velocities cancel or
                    // reinforce EVERY sleeve harmonic before the spatial bound.
                    let velocities = amplitudes.map(|a| sign * (SLEEVE_E + active) * a / eta);
                    let axial = amplitudes
                        .iter()
                        .zip(velocities)
                        .map(|(&a, v)| {
                            let mut state = initial(1, a)[0];
                            state.velocity_m_sqrt_kg_per_s =
                                initial(1, v)[0].displacement_m_sqrt_kg;
                            state
                        })
                        .collect();
                    let source = viscous_source(cx, gate, 1e9, eta);
                    let ea = source.axial_rigidity_n();
                    let loaded = lossy_image(
                        source,
                        TranslatingSupport {
                            initial_velocity_m_s: speed,
                            ..rail(displacement)
                        },
                        &modes,
                        axial,
                        finite_strain,
                        hereditary,
                        gate,
                    );
                    let run = loaded.prepare(config_at(RATE), gate).unwrap();
                    let memory_count = if hereditary { modes.len() + 1 } else { 0 };
                    assert_eq!(run.system().relaxation_branch_count(), memory_count);
                    assert_eq!(
                        run.system().state().len(),
                        2 * loaded.mode_count() + memory_count
                    );
                    let assessed = run.tensile_assessment(0).unwrap().unwrap();
                    let mean =
                        0.6 * displacement / LENGTH + 0.5 * (0.8 * displacement / LENGTH).powi(2);
                    let rate = 0.6 * speed / LENGTH
                        + if finite_strain {
                            0.8_f64.powi(2) * displacement * speed / (LENGTH * LENGTH)
                        } else {
                            0.0
                        };
                    let mut stress = SLEEVE_E * (20.0 / ea + mean) + active * mean + eta * rate;
                    let mut tension = 20.0 + (ea + active * area) * mean + eta * area * rate;
                    for ((&mode, &a), &v) in modes.iter().zip(&amplitudes).zip(&velocities) {
                        let epsilon = mode as f64 * wave * a;
                        let epsilon_dot = mode as f64 * wave * v;
                        stress += ((SLEEVE_E + active) * epsilon + eta * epsilon_dot).abs();
                        tension +=
                            ((ea + active * area) * epsilon + eta * area * epsilon_dot).abs();
                    }
                    assert_eq!(assessed.governing_layer, 1);
                    near(assessed.governing_tensile_stress_bound_pa, stress);
                    near(assessed.tension_n, tension);
                    assert_eq!(assessed.curvature_bound_per_m, 0.0);
                    // Original odd kinetic coordinates are normalized velocities,
                    // even after a distinct moving body shifts the mechanical
                    // prefix and every Maxwell memory's absolute state offset.
                    let mut bodies = vec![ImpactBody::free_mass(0.3, 0.25, 3.0).unwrap().0];
                    bodies.extend_from_slice(loaded.bodies());
                    let composed =
                        ImpactSystem::new(bodies, vec![], vec![], vec![], config_at(RATE)).unwrap();
                    assert_eq!(
                        assessed,
                        composed.support_tensile_assessment(1, 0).unwrap().unwrap()
                    );
                    assert!(
                        loaded
                            .supported()
                            .tensile_assessment(0, &initial_coordinates(&loaded))
                            .is_err()
                    );
                }
            }
        }
    });
}

#[test]
fn g0_actual_viscous_tensile_crossing_rolls_back_native_and_prepared_histories() {
    with_cx(|cx, gate| {
        let eta = 1e6;
        let prototype = viscous_source(cx, gate, 1e9, eta);
        let limit = SLEEVE_E * 20.0 / prototype.axial_rigidity_n() + 1e6;
        for hereditary in [false, true] {
            for finite_strain in [false, true] {
                let reference = lossy_image(
                    prototype.clone(),
                    rail(0.0),
                    &[1],
                    vec![ModalAcousticState::default()],
                    finite_strain,
                    hereditary,
                    gate,
                );
                let guarded = lossy_image(
                    viscous_source(cx, gate, limit, eta),
                    rail(0.0),
                    &[1],
                    vec![ModalAcousticState::default()],
                    finite_strain,
                    hereditary,
                    gate,
                );
                let mut force = vec![0.0; guarded.mode_count()];
                force[guarded.supported().axial_ranges()[0]
                    .as_ref()
                    .unwrap()
                    .start] = 10_000.0;
                let mut reference = reference.prepare(config_at(RATE), gate).unwrap();
                reference.step(&force, gate).unwrap();
                let total = reference.tensile_assessment(0).unwrap().unwrap();
                assert!(total.governing_tensile_stress_bound_pa > limit);
                // Removing only the selected dashpot leaves the original model
                // under the same applied force below the SAME limit. Actual rate
                // stress drives this crossing; no derivative is inferred from a
                // displacement increment or aggregate damping force.
                let quiet = lossy_image(
                    viscous_source(cx, gate, limit, 0.0),
                    rail(0.0),
                    &[1],
                    vec![ModalAcousticState::default()],
                    finite_strain,
                    hereditary,
                    gate,
                );
                let mut quiet = quiet.prepare(config_at(RATE), gate).unwrap();
                quiet.step(&force, gate).unwrap();
                assert!(
                    quiet
                        .tensile_assessment(0)
                        .unwrap()
                        .unwrap()
                        .governing_tensile_stress_bound_pa
                        < limit
                );

                let mut native = ImpactSystem::new(
                    guarded.bodies().to_vec(),
                    vec![],
                    vec![],
                    vec![],
                    config_at(RATE),
                )
                .unwrap();
                let before = native.state().to_vec();
                let energy = native.stored_energy_j();
                let error = native.step(&force, gate).unwrap_err();
                assert!(format!("{error:?}").contains("tensile"), "{error:?}");
                assert_eq!(native.state(), before);
                assert_eq!(native.stored_energy_j(), energy);
                assert_eq!(native.samples(), 0);
                let mut prepared = guarded.prepare(config_at(RATE), gate).unwrap();
                let assessed = prepared.tensile_assessment(0).unwrap().unwrap();
                let error = prepared.step(&force, gate).unwrap_err();
                assert!(format!("{error:?}").contains("tensile"), "{error:?}");
                assert_eq!(prepared.system().state(), before);
                assert_eq!(prepared.system().stored_energy_j(), energy);
                assert_eq!(prepared.system().samples(), 0);
                assert_eq!(prepared.tensile_assessment(0).unwrap().unwrap(), assessed);
            }
        }
    });
}

#[test]
fn g0_raw_axial_viscosity_cannot_silently_keep_an_elastic_tensile_guard() {
    with_cx(|cx, gate| {
        let source = viscous_source(cx, gate, 1e9, 1e6);
        let original = image(
            source,
            attachment(&[1], vec![ModalAcousticState::default()], BRIDGE),
            rail(0.0),
            gate,
        )
        .unwrap();
        for finite_strain in [false, true] {
            let native = if finite_strain {
                original
                    .supported()
                    .clone()
                    .with_finite_strain_axial_viscosity(vec![1.0])
            } else {
                original.supported().clone().with_axial_viscosity(vec![1.0])
            }
            .unwrap();
            assert!(
                native
                    .tensile_assessment(0, &initial_coordinates(&original))
                    .is_err()
            );
            let mut bodies = original.bodies().to_vec();
            bodies[0].potential = BodyPotential::SupportedStrings(native);
            let error = ImpactSystem::new(bodies, vec![], vec![], vec![], config_at(RATE));
            assert!(matches!(error, Err(ref e) if format!("{e:?}").contains("viscous")));
        }
    });
}

#[test]
fn g0_empty_raw_relaxation_preserves_viscous_initial_and_step_admission() {
    with_cx(|cx, gate| {
        let prototype = viscous_source(cx, gate, 1e9, 1e6);
        let limit = SLEEVE_E * 20.0 / prototype.axial_rigidity_n() + 1e6;
        for finite_strain in [false, true] {
            let guarded = lossy_image(
                viscous_source(cx, gate, limit, 1e6),
                rail(0.0),
                &[1],
                vec![ModalAcousticState::default()],
                finite_strain,
                false,
                gate,
            );
            let empty = guarded
                .supported()
                .clone()
                .with_longitudinal_relaxation(vec![], InitialMemory::Relaxed, 1)
                .unwrap();
            let mut bodies = guarded.bodies().to_vec();
            bodies[0].potential = BodyPotential::SupportedStrings(empty);
            let index = guarded.supported().axial_ranges()[0]
                .as_ref()
                .unwrap()
                .start;
            let mut native =
                ImpactSystem::new(bodies.clone(), vec![], vec![], vec![], config_at(RATE)).unwrap();
            assert_eq!(native.relaxation_branch_count(), 0);
            assert!(native.relaxation_observation().is_none());
            assert_eq!(native.state().len(), 2 * guarded.mode_count());
            let before = native.state().to_vec();
            let energy = native.stored_energy_j();
            let mut force = vec![0.0; guarded.mode_count()];
            force[index] = 10_000.0;
            let error = native.step(&force, gate).unwrap_err();
            assert!(format!("{error:?}").contains("tensile"), "{error:?}");
            assert_eq!(native.state(), before);
            assert_eq!(native.stored_energy_j(), energy);
            assert_eq!(native.samples(), 0);
            let mut prepared = native.prepare_analytic().unwrap();
            let error = prepared.step(&force, gate).unwrap_err();
            assert!(format!("{error:?}").contains("tensile"), "{error:?}");
            assert_eq!(prepared.state(), before);
            assert_eq!(prepared.stored_energy_j(), energy);
            assert_eq!(prepared.samples(), 0);

            // A public empty-arm declaration contributes no actual history
            // gate. It must not suppress the original initial velocity gate.
            bodies[0].initial[index].velocity_m_sqrt_kg_per_s = 0.02;
            let error = ImpactSystem::new(bodies, vec![], vec![], vec![], config_at(RATE));
            assert!(matches!(error, Err(ref e) if format!("{e:?}").contains("tensile")));
        }
    });
}

#[test]
fn g1_initial_signed_viscous_relief_is_admitted_after_complete_source_binding() {
    with_cx(|cx, gate| {
        let eta = 1e6;
        let prototype = viscous_source(cx, gate, 1e9, eta);
        let ea = prototype.axial_rigidity_n();
        let limit = SLEEVE_E * 20.0 / ea + 1e5;
        let displacement = 1e-4;
        let speed = -0.2;
        let support = TranslatingSupport {
            initial_velocity_m_s: speed,
            ..rail(displacement)
        };
        let source = viscous_source(cx, gate, limit, eta);
        let elastic = image(
            source.clone(),
            attachment(&[1], vec![ModalAcousticState::default()], BRIDGE),
            support,
            gate,
        );
        assert!(matches!(elastic, Err(ref e) if format!("{e:?}").contains("tensile")));
        for finite_strain in [false, true] {
            let loaded = lossy_image(
                source.clone(),
                support,
                &[1],
                vec![ModalAcousticState::default()],
                finite_strain,
                false,
                gate,
            );
            let run = loaded.prepare(config_at(RATE), gate).unwrap();
            let assessed = run.tensile_assessment(0).unwrap().unwrap();
            let mean = 0.6 * displacement / LENGTH + 0.5 * (0.8 * displacement / LENGTH).powi(2);
            let rate = 0.6 * speed / LENGTH
                + if finite_strain {
                    0.8_f64.powi(2) * displacement * speed / (LENGTH * LENGTH)
                } else {
                    0.0
                };
            let elastic_stress = SLEEVE_E * (20.0 / ea + mean);
            assert!(elastic_stress > limit);
            near(
                assessed.governing_tensile_stress_bound_pa,
                elastic_stress + eta * rate,
            );
            assert!(!assessed.exceeds_limit());
            assert!(assessed.tension_n > 0.0);
            assert_eq!(run.system().samples(), 0);
        }
    });
}
