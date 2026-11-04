//! G1/G3: original interface sources in actual quasistatic FE contact.
use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::contact::ContactFrictionBinding;
use fs_tribo::InputAuthority;

fn source(axis: &str, knots: [(f64, f64); 2]) -> SourceThermalCoulomb {
    source_with_kinetic_factor(axis, knots, 1.0)
}

fn source_with_kinetic_factor(
    axis: &str,
    knots: [(f64, f64); 2],
    kinetic_factor: f64,
) -> SourceThermalCoulomb {
    let (dims, kind) = match axis {
        "T" => (Temperature::DIMS, None),
        "P" => (Pressure::DIMS, None),
        "speed" => (
            Dims([1, 0, -1, 0, 0, 0]),
            Some(FrictionQueryAxisKind::SlidingSpeed),
        ),
        "path" => (
            Dims([1, 0, 0, 0, 0, 0]),
            Some(FrictionQueryAxisKind::AccumulatedDistance),
        ),
        _ => panic!("unknown synthetic source axis"),
    };
    let quantity = if axis == "P" {
        QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
    } else {
        QuantitySpec::dimensional(dims)
    };
    let mut claims = ClaimSet::new();
    let property = PropertyClaim {
        key: PropertyKey::new("contact-mu", Dims::NONE),
        value: PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: dims,
            knots: knots.to_vec(),
            dims: Dims::NONE,
        },
        validity: ValidityDomain::unconstrained()
            .with_quantity(axis, quantity, knots[0].0, knots[1].0),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::LinearInside,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic FE contact friction fixture; not measured data".into(),
            license: "CC0-1.0".into(),
            artifact: None,
        },
    };
    claims.insert_claim(property.clone()).unwrap();
    if kinetic_factor != 1.0 {
        let mut kinetic = property;
        kinetic.key = PropertyKey::new("contact-kinetic", Dims::NONE);
        kinetic.value = PropertyValue::Curve {
            abscissa: axis.into(),
            abscissa_dims: dims,
            knots: knots
                .iter()
                .map(|&(x, mu)| (x, kinetic_factor * mu))
                .collect(),
            dims: Dims::NONE,
        };
        claims.insert_claim(kinetic).unwrap();
    }
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "synthetic fixture".into(),
            revision: 0,
        },
        texture_frame: name.into(),
    };
    let interface = InterfaceSystemCard::assemble(
        surface("solid"),
        surface("obstacle"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "synthetic fixture".into(),
            history: "virgin".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let mut point = fs_matdb::QueryPoint::new().with("T", 300.0).unwrap();
    if axis != "T" {
        point = point.with_quantity(axis, quantity, knots[0].0).unwrap();
    }
    let requirement = ScalarPropertyRequirement::try_new(
        "contact-mu",
        Dims::NONE,
        ScalarAdmissibility::NonNegative,
    )
    .unwrap();
    let kinetic_requirement = if kinetic_factor == 1.0 {
        requirement.clone()
    } else {
        ScalarPropertyRequirement::try_new(
            "contact-kinetic",
            Dims::NONE,
            ScalarAdmissibility::NonNegative,
        )
        .unwrap()
    };
    SourceThermalCoulomb::resolve_with_query_axes(
        bind_dry_interface_system_card(&interface, InputAuthority::SyntheticFixture).unwrap(),
        point,
        requirement,
        kinetic_requirement,
        MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: (axis == "P").then(|| "P".into()),
        },
        kind.map(|kind| {
            vec![FrictionQueryAxis {
                name: axis.into(),
                kind,
            }]
        })
        .unwrap_or_default(),
    )
    .unwrap()
}

fn bindings(law: &SourceThermalCoulomb, n: usize) -> Vec<Option<ContactFrictionBinding<'_>>> {
    vec![
        Some(ContactFrictionBinding {
            law,
            conditions: law.initial_source().query(),
            nominal_contact_width_m: law.axes().nominal_pressure.as_ref().map(|_| 0.5),
        });
        n
    ]
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-11 * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_constant_source_recovers_original_unpinned_equilibrium_and_work() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let material = card();
    let pull = |_: f64, _: f64| [0.02, -0.2];
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![(Patch::Top, &pull)],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [0.0, 1.0],
        offset: -0.02,
    };
    let problem = timed_problem(&hyper, &floor);
    let reference = vec![0.0; 2 * mesh.nodes.len()];
    let stationary = vec![[0.0; 2]; mesh.nodes.len()];
    let original = problem
        .solve_timed_from(1.0, &reference, &stationary, 0.25, 0.001)
        .unwrap();
    let law = source("T", [(280.0, 0.6), (320.0, 0.6)]);
    let bindings = bindings(&law, mesh.nodes.len());
    let sourced = problem
        .solve_timed_from_with_sources(1.0, &reference, &stationary, 0.25, 0.001, &bindings)
        .unwrap();
    assert_eq!(sourced.u, original.u);
    assert_eq!(sourced.friction_reactions, original.friction_reactions);
    assert_eq!(sourced.friction_work, original.friction_work);
    assert_eq!(sourced.iterations, original.iterations);
    assert!(original.friction_sources.is_empty());
    assert_eq!(sourced.friction_sources.len(), mesh.nodes.len());
    for (sample, reaction) in sourced.friction_sources.iter().zip(&sourced.reactions) {
        if let Some(sample) = sample {
            close(
                sample.normal_reaction_n_per_m,
                reaction[0].hypot(reaction[1]),
            );
            assert_eq!(sample.law_identity, law.identity());
            assert_eq!(
                sample.source.card_identity(),
                law.initial_source().card_identity()
            );
        } else {
            assert_eq!(*reaction, [0.0; 2]);
        }
    }
}

#[test]
fn g1_local_temperature_pressure_and_speed_control_actual_interface_traction() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let material = card();
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [0.0, 1.0],
        offset: -0.02,
    };
    let mut problem = timed_problem(&hyper, &floor);
    // A supported fixture with prescribed zero displacement gives an independent
    // nodal traction oracle while the moving surface performs all interface work.
    problem.pins = (0..2 * mesh.nodes.len()).map(|i| (i, 0.0)).collect();
    problem.friction.as_mut().unwrap().mu = 0.99;
    let reference = vec![0.0; 2 * mesh.nodes.len()];
    let surface = vec![[0.001, 0.0]; mesh.nodes.len()];
    for axis in ["T", "P", "speed"] {
        let law = match axis {
            "T" => source("T", [(280.0, 0.8), (320.0, 0.2)]),
            "P" => source("P", [(0.0, 0.3), (10.0, 0.9)]),
            _ => source("speed", [(0.0, 0.2), (0.01, 0.8)]),
        };
        let mut bindings = bindings(&law, mesh.nodes.len());
        if axis == "T" {
            for (binding, node) in bindings.iter_mut().zip(&mesh.nodes) {
                binding.as_mut().unwrap().conditions = law
                    .initial_source()
                    .query()
                    .with("T", 280.0 + 40.0 * node[0])
                    .unwrap();
            }
        }
        let result = problem
            .solve_timed_from_with_sources(1.0, &reference, &surface, 0.5, 0.001, &bindings)
            .unwrap();
        let mut total = 0.0;
        for (i, sample) in result.friction_sources.iter().enumerate() {
            let Some(sample) = sample else {
                continue;
            };
            let normal = -problem.barrier.d1(mesh.nodes[i][1] + 0.02);
            let coefficient = match axis {
                "T" => 0.8 - 0.6 * mesh.nodes[i][0],
                "P" => 0.3 + 0.06 * normal / 0.5,
                _ => 0.2 + 60.0 * 0.002,
            };
            close(sample.coefficient, coefficient);
            close(sample.normal_reaction_n_per_m, normal);
            close(sample.relative_speed_m_s, 0.002);
            close(result.friction_reactions[i][0], coefficient * normal);
            close(result.friction_reactions[i][1], 0.0);
            if axis == "P" {
                close(sample.source.query().axes()["P"], normal / 0.5);
            }
            if axis == "speed" {
                close(sample.source.query().axes()["speed"], 0.002);
            }
            total += coefficient * normal * 0.001;
        }
        close(result.friction_work, total);
        close(result.friction_support_work, total);
        assert_eq!(result.u, reference);
    }
}

#[test]
fn g3_speed_source_closes_on_final_mechanical_motion_not_initial_lag() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let material = card();
    let pull = |_: f64, _: f64| [0.02, -0.2];
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![(Patch::Top, &pull)],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [0.0, 1.0],
        offset: -0.02,
    };
    let mut problem = timed_problem(&hyper, &floor);
    problem.friction.as_mut().unwrap().rounds = 40;
    let law = source("speed", [(0.0, 0.4), (0.01, 0.9)]);
    let bindings = bindings(&law, mesh.nodes.len());
    let reference = vec![0.0; 2 * mesh.nodes.len()];
    let stationary = vec![[0.0; 2]; mesh.nodes.len()];
    let result = problem
        .solve_timed_from_with_sources(1.0, &reference, &stationary, 0.5, 0.001, &bindings)
        .unwrap();
    let mut live = false;
    for (i, sample) in result.friction_sources.iter().enumerate() {
        if let Some(sample) = sample {
            let speed = result.u[2 * i].abs() / 0.5;
            close(sample.source.query().axes()["speed"], speed);
            close(sample.coefficient, 0.4 + 50.0 * speed);
            live |= sample.coefficient > 0.4001;
        }
    }
    assert!(live, "resolved motion must change the source coefficient");
    assert!(result.residual < problem.tol);
    assert!((result.friction_reactions.iter().map(|f| f[0]).sum::<f64>() + 0.02).abs() < 1e-8);
}

#[test]
fn g4_source_refusals_leave_native_time_displacement_and_work_unchanged() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let material = card();
    let pull = |_: f64, _: f64| [0.02, -0.2];
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![(Patch::Top, &pull)],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [0.0, 1.0],
        offset: -0.02,
    };
    let problem = timed_problem(&hyper, &floor);
    let reference = vec![0.0; 2 * mesh.nodes.len()];
    let stationary = vec![[0.0; 2]; mesh.nodes.len()];
    let law = source("T", [(280.0, 0.6), (320.0, 0.6)]);
    let mut bindings = bindings(&law, mesh.nodes.len());
    let mut trajectory = ContactTrajectory::new(&problem, 2.0, &reference, 0.001).unwrap();
    let solution = trajectory
        .advance_with_sources(0.5, 1.0, &stationary, &bindings)
        .unwrap();
    let missing = solution
        .friction_sources
        .iter()
        .position(Option::is_some)
        .unwrap();
    let accepted = trajectory.displacement().to_vec();
    let work = trajectory.work();
    let distance = trajectory.sliding_distance().to_vec();
    let original = bindings[missing].take().unwrap();
    assert!(
        trajectory
            .advance_with_sources(0.5, 1.0, &stationary, &bindings)
            .is_err()
    );
    bindings[missing] = Some(original);
    bindings[missing].as_mut().unwrap().conditions =
        law.initial_source().query().with("T", 321.0).unwrap();
    assert!(
        trajectory
            .advance_with_sources(0.5, 1.0, &stationary, &bindings)
            .is_err()
    );
    assert_eq!(trajectory.time_s(), 2.5);
    assert_eq!(trajectory.displacement(), accepted);
    assert_eq!(trajectory.work(), work);
    assert_eq!(trajectory.sliding_distance(), distance);
    let history_source = source("path", [(0.0, 0.6), (1.0, 0.6)]);
    let history_binding = self::bindings(&history_source, mesh.nodes.len());
    assert!(
        problem
            .solve_timed_from_with_sources(
                1.0,
                &accepted,
                &stationary,
                0.5,
                0.001,
                &history_binding,
            )
            .is_err()
    );
    assert_eq!(trajectory.displacement(), accepted);
    assert_eq!(trajectory.work(), work);
    let unequal = source_with_kinetic_factor("T", [(280.0, 0.6), (320.0, 0.6)], 0.5);
    let unequal_binding = self::bindings(&unequal, mesh.nodes.len());
    assert!(
        trajectory
            .advance_with_sources(0.5, 1.0, &stationary, &unequal_binding)
            .is_err()
    );
    assert_eq!(trajectory.time_s(), 2.5);
    assert_eq!(trajectory.displacement(), accepted);
    assert_eq!(trajectory.work(), work);
    assert_eq!(trajectory.sliding_distance(), distance);
}

#[test]
fn g3_running_in_history_drives_traction_counts_reversals_and_rolls_back() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let material = card();
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [0.0, 1.0],
        offset: -0.02,
    };
    let mut problem = timed_problem(&hyper, &floor);
    problem.pins = (0..2 * mesh.nodes.len()).map(|i| (i, 0.0)).collect();
    let reference = vec![0.0; 2 * mesh.nodes.len()];
    let prior_distance = vec![0.2; mesh.nodes.len()];
    let law = source("path", [(0.0, 0.8), (0.21, 0.38)]);
    let bindings = bindings(&law, mesh.nodes.len());
    let mut trajectory = ContactTrajectory::new_with_sliding_history(
        &problem,
        2.0,
        &reference,
        0.001,
        &prior_distance,
    )
    .unwrap();
    let mut expected_distance = 0.2;
    // Each interval resolves its reversal. Normal motion contributes no sliding,
    // and detached nodes do not acquire running-in history from surface travel.
    for surface in [[0.003, 0.0], [-0.001, 0.0], [0.0, 0.002]] {
        expected_distance += f64::abs(surface[0]);
        let solution = trajectory
            .advance_with_sources(0.5, 1.0, &vec![surface; mesh.nodes.len()], &bindings)
            .unwrap();
        for (i, sample) in solution.friction_sources.iter().enumerate() {
            let Some(sample) = sample else {
                close(trajectory.sliding_distance()[i], 0.2);
                close(solution.sliding_distance_increment[i], 0.0);
                continue;
            };
            close(trajectory.sliding_distance()[i], expected_distance);
            close(solution.sliding_distance_increment[i], surface[0].abs());
            close(sample.accumulated_distance_m.unwrap(), expected_distance);
            close(sample.source.query().axes()["path"], expected_distance);
            let mu = 0.8 - 2.0 * expected_distance;
            close(sample.coefficient, mu);
            let normal = -problem.barrier.d1(mesh.nodes[i][1] + 0.02);
            let direction = if surface[0] == 0.0 {
                0.0
            } else {
                surface[0].signum()
            };
            close(solution.friction_reactions[i][0], direction * mu * normal);
        }
        close(solution.friction_work, solution.friction_support_work);
        assert_eq!(solution.u, reference);
        close(trajectory.work().balance_residual(), 0.0);
    }
    // The source range ends at 0.21 m. Rejecting a physical advance must not age
    // the interface, and a subsequent admissible advance starts at 0.204 m.
    let accepted_distance = trajectory.sliding_distance().to_vec();
    let accepted_work = trajectory.work();
    assert!(
        trajectory
            .advance_with_sources(0.5, 1.0, &vec![[0.008, 0.0]; mesh.nodes.len()], &bindings)
            .is_err()
    );
    assert_eq!(trajectory.time_s(), 3.5);
    assert_eq!(trajectory.displacement(), reference);
    assert_eq!(trajectory.sliding_distance(), accepted_distance);
    assert_eq!(trajectory.work(), accepted_work);
    trajectory
        .advance_with_sources(0.5, 1.0, &vec![[0.001, 0.0]; mesh.nodes.len()], &bindings)
        .unwrap();
    for (i, node) in mesh.nodes.iter().enumerate() {
        close(
            trajectory.sliding_distance()[i],
            if node[1] == 0.0 { 0.205 } else { 0.2 },
        );
        assert_eq!(bindings[i].as_ref().unwrap().conditions.axes()["path"], 0.0);
    }
    assert!(
        ContactTrajectory::new_with_sliding_history(&problem, 0.0, &reference, 0.001, &[-1.0],)
            .is_err()
    );
}

#[test]
fn g3_distance_source_uses_final_equilibrium_slip_once_per_physical_interval() {
    let mesh = Mesh2::quads(1.0, 1.0, 2, 2);
    let material = card();
    let pull = |_: f64, _: f64| [0.02, -0.2];
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &material,
        dirichlet: vec![],
        traction: vec![(Patch::Top, &pull)],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [0.0, 1.0],
        offset: -0.02,
    };
    let mut problem = timed_problem(&hyper, &floor);
    problem.friction.as_mut().unwrap().rounds = 40;
    let reference = vec![0.0; 2 * mesh.nodes.len()];
    let surface = vec![[0.0; 2]; mesh.nodes.len()];
    let accepted_distance = vec![0.001; mesh.nodes.len()];
    let law = source("path", [(0.0, 0.8), (0.01, 0.2)]);
    let bindings = bindings(&law, mesh.nodes.len());
    let solution = problem
        .solve_timed_from_with_source_history(
            1.0,
            &reference,
            &surface,
            0.5,
            0.001,
            &bindings,
            &accepted_distance,
        )
        .unwrap();
    assert!(solution.residual < problem.tol);
    for (i, sample) in solution.friction_sources.iter().enumerate() {
        if let Some(sample) = sample {
            let distance = accepted_distance[i] + solution.u[2 * i].abs();
            close(sample.source.query().axes()["path"], distance);
            close(sample.accumulated_distance_m.unwrap(), distance);
            close(sample.coefficient, 0.8 - 60.0 * distance);
            close(
                solution.sliding_distance_increment[i],
                solution.u[2 * i].abs(),
            );
        }
    }
    assert!(
        (solution
            .friction_reactions
            .iter()
            .map(|f| f[0])
            .sum::<f64>()
            + 0.02)
            .abs()
            < 1e-8
    );
    problem.friction.as_mut().unwrap().rounds = 80;
    let more = problem
        .solve_timed_from_with_source_history(
            1.0,
            &reference,
            &surface,
            0.5,
            0.001,
            &bindings,
            &accepted_distance,
        )
        .unwrap();
    assert_eq!(solution.u, more.u);
    assert_eq!(
        solution.sliding_distance_increment,
        more.sliding_distance_increment
    );
    assert_eq!(solution.friction_sources, more.friction_sources);
    assert_eq!(accepted_distance, vec![0.001; mesh.nodes.len()]);
}
