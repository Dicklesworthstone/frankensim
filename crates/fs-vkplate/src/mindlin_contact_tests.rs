//! G1/G3: source-derived contact memory participates in the actual body work.
//!
//! These synthetic material cards exercise the real source admission; they do
//! not claim that this point-contact fixture is a measured instrument joint.

#[path = "../../fs-contact/tests/source_point/common.rs"]
mod source_fixture;

use super::contact_island_tests::{dense_solve, loaded_body};
use super::mindlin_contact::MindlinContact;
use super::*;
use fs_contact::interface_binding::mindlin::SourceMindlinContact;
use fs_contact::interface_binding::*;
use fs_contact::normal_patch::source::{PointContactGeometry, SourcePointContact};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    InterfaceSystemCard, InterpolationPolicy, PropertyClaim, PropertyKey, PropertyValue,
    Provenance, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    resolve_interface_state_point, resolve_isotropic_solid_state_point,
};
use fs_qty::Dims;
use fs_tribo::{InputAuthority, TangentialLaw, TangentialSolveConfig, TangentialStepInput};

fn source(young: f64) -> SourceMindlinContact {
    let a = source_fixture::solid("resonance contact A", 1000.0, young);
    let b = source_fixture::solid("resonance contact B", 2000.0, 2e7);
    let normal = source_fixture::bind(&a, &b, None, 0.0);
    let original = normal.normal_state().elastic().interface().card();
    let mut claims = original.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("local-coulomb-coefficient", Dims::NONE),
            value: PropertyValue::Scalar {
                value: 0.6,
                dims: Dims::NONE,
            },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 320.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic ordered-interface coefficient for the resonance fixture".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = InterfaceSystemCard::assemble(
        original.surface_a().clone(),
        original.surface_b().clone(),
        original.context().clone(),
        claims,
        original.models().to_vec(),
    )
    .unwrap();
    let query = source_fixture::query();
    let selection = MaterialPropertySelection::SingleClaimOnly;
    let sa = resolve_isotropic_solid_state_point(&a, &query, selection.clone()).unwrap();
    let sb = resolve_isotropic_solid_state_point(&b, &query, selection.clone()).unwrap();
    let dry = bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap();
    let elastic = bind_isotropic_elastic_interface(&dry, &sa, &sb).unwrap();
    let interface = resolve_interface_state_point(
        &card,
        &query,
        &[ScalarPropertyRequirement::try_new(
            ADHESION_ENERGY_PROPERTY,
            ADHESION_ENERGY_DIMS,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()],
        selection.clone(),
    )
    .unwrap();
    let state = bind_normal_interface_state(&elastic, &interface).unwrap();
    let normal =
        bind_normal_contact_model(&state, NormalContactModelSelection::SingleSupported).unwrap();
    let point = SourcePointContact::new(
        normal,
        PointContactGeometry::SpherePlane {
            effective_radius_m: 0.1,
        },
        source_fixture::envelope(),
    )
    .unwrap();
    SourceMindlinContact::resolve(
        point,
        2e-5,
        ScalarPropertyRequirement::try_new(
            "local-coulomb-coefficient",
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap(),
        selection,
        64,
    )
    .unwrap()
}

fn law_controls() -> TangentialSolveConfig {
    TangentialSolveConfig {
        max_iterations: 96,
        maximum_force_n: 100.0,
        force_absolute_tolerance_n: 1e-12,
        force_relative_tolerance: 1e-12,
    }
}

fn input(source: &SourceMindlinContact, slip: f64, k: f64) -> TangentialStepInput {
    TangentialStepInput {
        normal_n: source.normal_response().normal_force_n,
        slip_increment_m: [slip, 0.0],
        dt_s: k,
        contact_temperature_k: Some(source.temperature_kelvin()),
        nominal_pressure_pa: None,
        compliance_m_per_n: 0.0,
        solve: law_controls(),
    }
}

fn solve(
    bodies: &mut [Body],
    contacts: &mut [Contact],
    frictions: &mut [Friction],
    mindlin: &mut [MindlinContact],
    k: f64,
    controls: FrictionSolveControls,
) -> Result<(ContactSolveReport, FrictionSolveReport), ContactSolveError> {
    friction_response::solve_contacts_with_friction_checked(
        bodies,
        contacts,
        frictions,
        mindlin,
        &mut ContactQueue::default(),
        &mut contact_response::Workspace::default(),
        &mut friction_response::Workspace::default(),
        k,
        ContactSolveControls::default(),
        controls,
    )
}

#[test]
fn stored_shear_energy_can_return_to_a_resting_body_without_becoming_negative_heat() {
    let k = 0.001;
    let law = source(1e7);
    let loaded = law
        .advance(
            &law.initial_state(),
            &input(&law, 0.25 * law.law().gross_slip_displacement_m(), k),
        )
        .unwrap();
    let before = loaded.state.clone();
    let stored_before = law.stored_energy_j(&before);
    let mut ports = vec![
        MindlinContact::new(law.clone(), vec![(0, 0, 1.0)], loaded.state, law_controls()).unwrap(),
    ];
    let mut body = Body::point("shear energy release", 1.0);
    body.predict(k);
    let predicted = body.wn[0];
    solve(
        std::slice::from_mut(&mut body),
        &mut [],
        &mut [],
        &mut ports,
        k,
        FrictionSolveControls::default(),
    )
    .unwrap();
    let traction = body.mass[0] * (predicted - body.wn[0]) / (k * k);
    let increment = 0.5 * (body.wn[0] - body.wp[0]);
    let expected = law.advance(&before, &input(&law, increment, k)).unwrap();
    let report = ports[0].report();
    let kinetic = 0.5 * body.mass[0] * ((body.wn[0] - body.w[0]) / k).powi(2);
    assert!(traction > 0.0 && increment < 0.0);
    assert!(report.mechanical_work_j > 0.0 && report.stored_change_j < 0.0);
    assert!(report.dissipated_j >= 0.0);
    assert!((report.traction_n - traction).abs() < 2e-10);
    assert!((report.slip_increment_m - increment).abs() < 2e-15);
    assert!((report.mechanical_work_j - kinetic).abs() < 2e-16);
    assert_eq!(ports[0].state(), &expected.state);
    assert!(
        (ports[0].stored_energy_j() + kinetic + report.dissipated_j - stored_before).abs()
            < 1e-16 + 2e-9 * stored_before
    );
}

#[test]
fn admissible_body_response_survives_trial_and_checkpoint_force_caps() {
    let k = 0.0005;
    let law = source(1e7);
    let gross_slip = law.law().gross_slip_displacement_m();
    let capacity = law.law().gross_slip_capacity_n();
    let budget = TangentialSolveConfig {
        maximum_force_n: 0.2 * capacity,
        ..law_controls()
    };
    // Half of the body's centered compliance is the constitutive load line:
    // C_midpoint = k^2/(2m) = 5*d/Q. The released/free mean force exceeds the
    // budget, while the compatible body/contact reaction is below it.
    let mass = k * k * capacity / (10.0 * gross_slip);
    for loaded_checkpoint in [false, true] {
        let state = if loaded_checkpoint {
            let loaded = law
                .advance(
                    &law.initial_state(),
                    &TangentialStepInput {
                        solve: budget,
                        ..input(&law, 0.2 * gross_slip, k)
                    },
                )
                .unwrap();
            assert!(loaded.traction_n[0] < budget.maximum_force_n);
            assert!(loaded.state.endpoint_force_n() > budget.maximum_force_n);
            loaded.state
        } else {
            assert!(
                law.advance(
                    &law.initial_state(),
                    &TangentialStepInput {
                        solve: budget,
                        ..input(&law, gross_slip, k)
                    },
                )
                .is_err()
            );
            law.initial_state()
        };
        let before = state.clone();
        let mut ports =
            vec![MindlinContact::new(law.clone(), vec![(0, 0, 1.0)], state, budget).unwrap()];
        let mut body = Body::point("force-limited retained shear", mass);
        if !loaded_checkpoint {
            body.w[0] = gross_slip;
        }
        body.predict(k);
        solve(
            std::slice::from_mut(&mut body),
            &mut [],
            &mut [],
            &mut ports,
            k,
            FrictionSolveControls::default(),
        )
        .unwrap();
        let du = 0.5 * (body.wn[0] - body.wp[0]);
        let expected = law
            .advance(
                &before,
                &TangentialStepInput {
                    solve: budget,
                    ..input(&law, du, k)
                },
            )
            .unwrap();
        let report = ports[0].report();
        assert!(report.traction_n.abs() < budget.maximum_force_n);
        assert!((report.traction_n - expected.traction_n[0]).abs() < 2e-11);
        assert_eq!(ports[0].state(), &expected.state);
        if loaded_checkpoint {
            assert!(report.mechanical_work_j > 0.0);
        } else {
            assert!((0.15..0.16).contains(&(report.traction_n / capacity)));
        }
    }
}

fn oscillator(young: f64, fraction: f64, k: f64) -> Instrument {
    let law = source(young);
    let displacement = fraction * law.law().gross_slip_displacement_m();
    let loaded = law
        .advance(&law.initial_state(), &input(&law, displacement, k))
        .unwrap();
    let port = MindlinContact::new(law, vec![(0, 0, 1.0)], loaded.state, law_controls()).unwrap();
    let mut body = Body::point("material-derived shear resonance", 0.002);
    body.w[0] = displacement;
    body.wp[0] = displacement;
    let mut instrument = Instrument::new(k.recip());
    instrument.bodies.push(body);
    instrument.attach_mindlin_contact(port).unwrap();
    instrument
}

#[test]
fn retained_microslip_memory_closes_instrument_ringdown_and_clone_continuation() {
    let k = 2e-5;
    let mut instrument = oscillator(1e7, 0.25, k);
    let initial = instrument.mindlin_contacts()[0].stored_energy_j();
    let (mut worst, mut returned) = (0.0_f64, 0);
    for _ in 0..1200 {
        let old_w = instrument.bodies[0].w[0];
        let old_wp = instrument.bodies[0].wp[0];
        instrument.try_step(&mut []).unwrap();
        let body = &instrument.bodies[0];
        let port = &instrument.mindlin_contacts()[0];
        let report = port.report();
        let actual_increment = 0.5 * (body.w[0] - old_wp);
        let kinetic = 0.5 * body.mass[0] * ((body.w[0] - old_w) / k).powi(2);
        let total = kinetic + port.stored_energy_j() + port.dissipated_j();
        worst = worst.max(((total - initial) / initial).abs());
        assert!((report.slip_increment_m - actual_increment).abs() < 2e-14);
        assert!(report.dissipated_j >= 0.0);
        returned += usize::from(report.mechanical_work_j > 0.0);
    }
    assert!(returned > 10);
    assert!(instrument.mindlin_contacts()[0].dissipated_j() > 0.01 * initial);
    assert!(worst < 3e-8, "Mindlin resonance energy closure {worst:e}");
    let mut replay = instrument.clone();
    for _ in 0..12 {
        instrument.try_step(&mut []).unwrap();
        replay.try_step(&mut []).unwrap();
        assert_eq!(instrument.bodies[0].w, replay.bodies[0].w);
        assert_eq!(
            instrument.mindlin_contacts()[0].state(),
            replay.mindlin_contacts()[0].state()
        );
        assert_eq!(
            instrument.mindlin_contacts()[0].dissipated_j(),
            replay.mindlin_contacts()[0].dissipated_j()
        );
    }
    eprintln!("Mindlin instrument: 1200-step relative storage+heat closure {worst:.12e}");
}

#[test]
fn bulk_elasticity_sets_the_small_amplitude_contact_resonance() {
    let k = 2e-5;
    for young in [1e7, 3e7] {
        let mut instrument = oscillator(young, 0.001, k);
        let mut last = instrument.bodies[0].w[0];
        let mut crossings = Vec::new();
        for step in 0..800 {
            instrument.try_step(&mut []).unwrap();
            let current = instrument.bodies[0].w[0];
            if last > 0.0 && current <= 0.0 {
                crossings.push((step as f64 + last / (last - current)) * k);
            }
            last = current;
        }
        assert!(crossings.len() >= 5);
        let measured = (crossings.len() - 1) as f64 / (crossings.last().unwrap() - crossings[0]);
        // Independent circular Hertz radius and Mindlin full-stick stiffness
        // from the two declared bulk moduli, not a caller-supplied spring.
        let radius = (0.1_f64 * 2e-5).sqrt();
        let stiffness = 8.0 * radius / (1.75 / (young / 2.5) + 1.75 / (2e7 / 2.5));
        let omega = (stiffness / instrument.bodies[0].mass[0]).sqrt();
        // The midpoint shear coordinate gives this linear-limit dispersion.
        // Small but finite partial slip is allowed a 0.2% frequency difference.
        let expected = 2.0 * (0.5 * k * omega).atan() / (2.0 * PI * k);
        assert!(
            (measured / expected - 1.0).abs() < 0.002,
            "E={young:e} Pa: measured {measured} Hz, linear limit {expected} Hz"
        );
        eprintln!(
            "Mindlin resonance E={young:e} Pa: {measured:.9} Hz, linear limit {expected:.9} Hz"
        );
    }
}

fn loaded_operator(body: &Body, k: f64) -> Vec<f64> {
    let load = body.head_radiation.as_ref().unwrap();
    let modes = load.cols.len();
    let mut operator = vec![0.0; 16];
    for i in 0..4 {
        for j in 0..4 {
            let mut value = body.air_mass * body.h.powi(4) / (k * k) + body.g[i] * body.g[j] / 4.0;
            if i == j {
                value += body.cf[i][0].recip();
            }
            for r in 0..modes {
                for s in 0..modes {
                    value += load.cols[r][i]
                        * load.cols[s][j]
                        * (load.mass[r * modes + s] / (k * k)
                            + load.damping[r * modes + s] / (2.0 * k));
                }
            }
            operator[4 * i + j] = value;
        }
    }
    operator
}

#[test]
fn normal_coulomb_and_retained_shear_share_the_full_loaded_final_pose() {
    let k = 0.0005;
    let mut body = loaded_body(k);
    let law = source(1e7);
    let n = [1.0, -0.2, 0.3, 0.0];
    let tx = [0.2, 1.0, -0.3, 0.4];
    let ty = [-0.1, 0.2, 0.5, 1.0];
    let shear = [8.0, -6.0, 4.0, 5.0];
    let target = [6e-5, 2e-5, 4e-5, -1e-5];
    let delta = [2e-6, -1e-6, 1.5e-6, -0.5e-6];
    let previous: [f64; 4] = std::array::from_fn(|i| target[i] - delta[i]);
    let dot = |a: &[f64; 4], b: &[f64; 4]| a.iter().zip(b).map(|(a, b)| a * b).sum::<f64>();
    let (d, dp) = (dot(&n, &target), dot(&n, &previous));
    let stiffness = 2e7;
    let normal_force = stiffness * (d.powf(2.5) - dp.powf(2.5)) / (2.5 * (d - dp));
    let slip = [dot(&tx, &delta), dot(&ty, &delta)];
    let capacity = 0.25 * normal_force;
    let traction = [
        capacity * slip[0] / slip[0].hypot(slip[1]),
        capacity * slip[1] / slip[0].hypot(slip[1]),
    ];
    let increment = 0.5 * dot(&shear, &delta);
    let loaded = law
        .advance(
            &law.initial_state(),
            &input(&law, 0.1 * law.law().gross_slip_displacement_m(), k),
        )
        .unwrap();
    let expected = law
        .advance(&loaded.state, &input(&law, increment, k))
        .unwrap();
    let reaction = (0..4)
        .map(|i| {
            n[i] * normal_force
                + tx[i] * traction[0]
                + ty[i] * traction[1]
                + shear[i] * expected.traction_n[0]
        })
        .collect();
    let correction = dense_solve(loaded_operator(&body, k), reaction);
    body.wp = previous.to_vec();
    body.wn = (0..4).map(|i| target[i] + correction[i]).collect();
    let terms = |column: &[f64; 4]| column.iter().enumerate().map(|(i, &v)| (0, i, v)).collect();
    let mut contacts = vec![Contact::new(
        terms(&n),
        vec![],
        0.0,
        Law::Power {
            k: stiffness,
            alpha: 1.5,
            mu: 0.0,
        },
    )];
    let mut frictions = vec![Friction {
        contact: 0,
        sx: terms(&tx),
        sy: terms(&ty),
        mu: 0.25,
        force: 0.0,
        stuck: false,
    }];
    let mut ports =
        vec![MindlinContact::new(law, terms(&shear), loaded.state, law_controls()).unwrap()];
    solve(
        std::slice::from_mut(&mut body),
        &mut contacts,
        &mut frictions,
        &mut ports,
        k,
        FrictionSolveControls::default(),
    )
    .unwrap();
    for (actual, expected) in body.wn.iter().zip(target) {
        assert!((actual - expected).abs() < 3e-13);
    }
    assert!((contacts[0].force / normal_force - 1.0).abs() < 3e-9);
    assert!((frictions[0].force / capacity - 1.0).abs() < 3e-9);
    assert!(!frictions[0].stuck);
    assert!((ports[0].report().traction_n - expected.traction_n[0]).abs() < 2e-9);
    assert!((ports[0].state().slip_m() - expected.state.slip_m()).abs() < 2e-12);
    assert!(
        (ports[0].stored_energy_j() - ports[0].source().stored_energy_j(&expected.state)).abs()
            < 2e-13
    );
    assert!((ports[0].report().dissipated_j - expected.dissipated_j).abs() < 2e-13);
}

#[test]
fn a_joint_refusal_preserves_retained_history_and_every_physical_contact_update() {
    let k = 0.0005;
    for force_budget in [false, true] {
        let mut bodies: Vec<_> = [0.001, 0.0001, 0.0]
            .into_iter()
            .map(|q| {
                let mut body = Body::point("atomic retained contact", 1.0);
                body.w[0] = q;
                body.wp[0] = 0.9 * q;
                body.predict(k);
                body
            })
            .collect();
        let mut contacts = vec![Contact::new(
            vec![(0, 0, 1.0)],
            vec![],
            0.0,
            Law::Power {
                k: 1e6,
                alpha: 1.5,
                mu: 0.0,
            },
        )];
        contacts[0].crush = Some(Crush::new(0.01, 3e6).unwrap());
        contacts[0].peak_force = 9.0;
        contacts[0].d_now = -0.001;
        let mut frictions = vec![Friction {
            contact: 0,
            sx: vec![(1, 0, 1.0)],
            sy: vec![(2, 0, 1.0)],
            mu: 0.1,
            force: 7.0,
            stuck: true,
        }];
        let law = source(1e7);
        let mut budget = law_controls();
        if force_budget {
            budget.maximum_force_n = 1e-6;
        }
        let state = law.initial_state();
        let mut ports = vec![
            MindlinContact::new(law, vec![(0, 0, 0.2), (1, 0, 1.0)], state.clone(), budget)
                .unwrap(),
        ];
        let old_report = ports[0].report();
        let old_pose: Vec<_> = bodies.iter().map(|b| b.wn.clone()).collect();
        let old_contacts = contacts.clone();
        let controls = if force_budget {
            FrictionSolveControls::default()
        } else {
            FrictionSolveControls {
                max_newton_iterations: 0,
                ..FrictionSolveControls::default()
            }
        };
        let error = solve(
            &mut bodies,
            &mut contacts,
            &mut frictions,
            &mut ports,
            k,
            controls,
        )
        .unwrap_err();
        if force_budget {
            assert!(error.to_string().contains("maximum_force_n"), "{error}");
        }
        for (body, previous) in bodies.iter().zip(old_pose) {
            assert_eq!(body.wn, previous);
        }
        assert_eq!(contacts[0].force, old_contacts[0].force);
        assert_eq!(contacts[0].peak_force, old_contacts[0].peak_force);
        assert_eq!(contacts[0].d_now, old_contacts[0].d_now);
        assert_eq!(contacts[0].crush, old_contacts[0].crush);
        assert_eq!(frictions[0].force, 7.0);
        assert!(frictions[0].stuck);
        assert_eq!(ports[0].state(), &state);
        assert_eq!(ports[0].report().traction_n, old_report.traction_n);
        assert_eq!(
            ports[0].report().stored_change_j,
            old_report.stored_change_j
        );
        assert_eq!(ports[0].report().dissipated_j, old_report.dissipated_j);
        assert_eq!(
            ports[0].report().mechanical_work_j,
            old_report.mechanical_work_j
        );
        assert_eq!(ports[0].stored_energy_j(), 0.0);
        assert_eq!(ports[0].dissipated_j(), 0.0);
    }
}
