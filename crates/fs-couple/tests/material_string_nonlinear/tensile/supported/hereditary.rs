//! Actual layer histories, independent constitutive references and rollback.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::longitudinal::relaxation::LongitudinalRelaxationSelection;
use fs_couple::render::plate::impact::ImpactBody;
use fs_material::section::relaxation::LayerRelaxationProperties;
use fs_qty::Time;

const CORE_E: f64 = 100e9;
const SLEEVE_E: f64 = 2e9;
const DELTA: f64 = 20e9;
const CORE_R: f64 = 0.0002;

fn material(young: f64, delta: f64, limit: f64) -> ResolvedMaterialStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [
        ("density", Density::DIMS, 7800.0),
        ("young_modulus", Pressure::DIMS, young),
        (ETA, DynViscosity::DIMS, 0.0),
        (LIMIT, Pressure::DIMS, limit),
        ("delta", Pressure::DIMS, delta),
        ("tau", Time::DIMS, 0.001),
    ] {
        let validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: if name == LIMIT {
                    validity
                } else {
                    validity.with("omega", 0.0, 1e8)
                },
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured layered tensile-memory fixture".into(),
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
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "hereditary tensile".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_material_state_point(
        &card,
        &QueryPoint::new()
            .with("T", 293.15)
            .unwrap()
            .with("omega", 1000.0)
            .unwrap(),
        &requirements,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn source(cx: &Cx<'_>, gate: &CancelGate, limit: f64) -> NonlinearMaterialString {
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![
            ConcentricMaterialLayer {
                outer_radius_m: CORE_R,
                material: material(CORE_E, 0.0, 1e12),
            },
            ConcentricMaterialLayer {
                outer_radius_m: RADIUS,
                material: material(SLEEVE_E, DELTA, limit),
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
        [1.0, 0.0],
    )
    .with_tensile_limit(LIMIT)
    .unwrap()
}

fn selection() -> Vec<LongitudinalRelaxationSelection> {
    vec![LongitudinalRelaxationSelection {
        string: 0,
        layers: vec![
            LayerRelaxationProperties {
                equilibrium_modulus: "young_modulus".into(),
                branches: vec![],
            },
            LayerRelaxationProperties {
                equilibrium_modulus: "young_modulus".into(),
                branches: vec![("delta".into(), "tau".into())],
            },
        ],
    }]
}

#[test]
fn g1_layer_memory_reinforcement_and_relaxed_cancellation_keep_original_state_offsets() {
    with_cx(|cx, gate| {
        let displacement = 2e-5;
        let amplitudes = [2e-6, 1e-6];
        let source = source(cx, gate, 1e9);
        let ea = source.axial_rigidity_n();
        let original = image(
            source,
            attachment(
                &[1, 3],
                amplitudes.iter().map(|&a| initial(1, a)[0]).collect(),
                [1.0, 0.0],
            ),
            rail(displacement),
            gate,
        )
        .unwrap();
        let wave = core::f64::consts::PI / LENGTH;
        let dynamic_peak = 0.6 * displacement / LENGTH
            + 0.5 * (0.8 * displacement / LENGTH).powi(2)
            + wave * (amplitudes[0] + 3.0 * amplitudes[1]);
        let area = core::f64::consts::PI * (RADIUS * RADIUS - CORE_R * CORE_R);
        for initial_memory in [InitialMemory::Relaxed, InitialMemory::Unrelaxed] {
            let loaded = original
                .clone()
                .with_material_longitudinal_relaxation(selection(), initial_memory, 3, gate)
                .unwrap();
            assert!(
                loaded
                    .supported()
                    .tensile_assessment(0, &initial_coordinates(&loaded))
                    .is_err()
            );
            // A later raw viscosity attachment cannot bypass total axial stress
            // admission by discarding its not-yet-resolved layer dashpot stress.
            let lossy = loaded
                .supported()
                .clone()
                .with_axial_viscosity(vec![1e-8])
                .unwrap();
            let mut lossy_bodies = loaded.bodies().to_vec();
            lossy_bodies[0].potential = BodyPotential::SupportedStrings(lossy);
            let error = ImpactSystem::new(lossy_bodies, vec![], vec![], vec![], config_at(RATE));
            assert!(matches!(error, Err(ref e) if format!("{e:?}").contains("viscous")));
            let mut bodies = vec![ImpactBody::free_mass(0.3, 0.0, 0.0).unwrap().0];
            bodies.extend_from_slice(loaded.bodies());
            let composed =
                ImpactSystem::new(bodies, vec![], vec![], vec![], config_at(RATE)).unwrap();
            assert!(composed.support_tensile_assessment(0, 0).unwrap().is_none());
            let run = loaded.prepare(config_at(RATE), gate).unwrap();
            let assessed = run.tensile_assessment(0).unwrap().unwrap();
            assert_eq!(
                assessed,
                composed.support_tensile_assessment(1, 0).unwrap().unwrap()
            );
            let active = if initial_memory == InitialMemory::Unrelaxed {
                DELTA
            } else {
                0.0
            };
            assert_eq!(assessed.governing_layer, 1);
            near(
                assessed.governing_tensile_stress_bound_pa,
                SLEEVE_E * (20.0 / ea + dynamic_peak) + active * dynamic_peak,
            );
            near(
                assessed.tension_n,
                20.0 + (ea + active * area) * dynamic_peak,
            );
            assert_eq!(assessed.curvature_bound_per_m, 0.0);
            // An aggregate force cannot reconstruct each layer's Maxwell stress.
            if active > 0.0 {
                assert!(
                    assessed.governing_tensile_stress_bound_pa > SLEEVE_E * assessed.tension_n / ea
                );
            }
        }
        let mut wrong_layers = selection();
        wrong_layers[0].layers.pop();
        assert!(
            original
                .clone()
                .with_material_longitudinal_relaxation(
                    wrong_layers,
                    InitialMemory::Relaxed,
                    3,
                    gate
                )
                .is_err()
        );
        let mut wrong_equilibrium = selection();
        wrong_equilibrium[0].layers[1].equilibrium_modulus = "delta".into();
        assert!(
            original
                .with_material_longitudinal_relaxation(
                    wrong_equilibrium,
                    InitialMemory::Relaxed,
                    3,
                    gate
                )
                .is_err()
        );
    });
}

#[test]
fn g0_actual_hereditary_tensile_crossing_rolls_back_native_and_prepared_histories() {
    with_cx(|cx, gate| {
        let prototype = source(cx, gate, 1e9);
        let limit = SLEEVE_E * 20.0 / prototype.axial_rigidity_n() + 1e6;
        let axial = vec![ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: 0.5,
        }];
        let support = attachment(&[1], axial, [1.0, 0.0]);
        let reference = image(prototype, support.clone(), rail(0.0), gate)
            .unwrap()
            .with_material_longitudinal_relaxation(selection(), InitialMemory::Relaxed, 2, gate)
            .unwrap();
        let elastic = image(source(cx, gate, limit), support, rail(0.0), gate).unwrap();
        let guarded = elastic
            .clone()
            .with_material_longitudinal_relaxation(selection(), InitialMemory::Relaxed, 2, gate)
            .unwrap();
        let force = vec![0.0; guarded.mode_count()];
        let mut reference = reference.prepare(config_at(RATE), gate).unwrap();
        reference.step(&force, gate).unwrap();
        assert!(
            reference
                .tensile_assessment(0)
                .unwrap()
                .unwrap()
                .governing_tensile_stress_bound_pa
                > limit
        );
        let q: Vec<_> = reference.system().state()[..2 * elastic.mode_count()]
            .chunks_exact(2)
            .map(|s| s[0])
            .collect();
        assert!(
            !elastic
                .supported()
                .tensile_assessment(0, &q)
                .unwrap()
                .unwrap()
                .exceeds_limit()
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
    });
}
