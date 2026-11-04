//! Native solid friction and finite-cell heat against a scalar elimination.
use super::*;
use fs_contact::interface_binding::bind_dry_interface_system_card;
use fs_contact::interface_binding::thermal_friction::{
    FrictionQueryAxis, FrictionQueryAxisKind, SourceThermalCoulomb, ThermalFrictionAxes,
};
use fs_couple::friction_heat::solid::*;
use fs_cutfem::sdf::HalfPlane;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::hyper::{Hyperelastic, HyperelasticModel};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, QuantitySpec, Temperature};
use fs_solid::contact::{
    Barrier, ContactFrictionBinding, ContactProblem, ContactTrajectory, Friction,
};
use fs_solid::hyper2d::{HyperProblem, NewtonSettings};
use fs_solid::mesh2::Mesh2;
use fs_solver::NewtonKrylovConfig;
use fs_tribo::InputAuthority;
use std::cell::Cell;

fn friction_source(extra_temperature: bool, declare_temperature: bool) -> SourceThermalCoulomb {
    let length = QuantitySpec::dimensional(Dims([1, 0, 0, 0, 0, 0]));
    let mut validity = ValidityDomain::unconstrained()
        .with("T", 299.0, 301.0)
        .with_quantity("path", length, 0.0, 1.0);
    let mut point = QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with_quantity("path", length, 0.0)
        .unwrap();
    let mut axes = vec![FrictionQueryAxis {
        name: "path".into(),
        kind: FrictionQueryAxisKind::AccumulatedDistance,
    }];
    if extra_temperature {
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        validity = validity.with_quantity("body_T", temperature, 299.0, 301.0);
        point = point.with_quantity("body_T", temperature, 300.0).unwrap();
        if declare_temperature {
            axes.push(FrictionQueryAxis {
                name: "body_T".into(),
                kind: FrictionQueryAxisKind::Temperature,
            });
        }
    }
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("mu", Dims::NONE),
            value: PropertyValue::Curve {
                abscissa: "T".into(),
                abscissa_dims: Temperature::DIMS,
                knots: vec![(299.0, 1.3), (301.0, 0.3)],
                dims: Dims::NONE,
            },
            validity,
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic linear friction/heat oracle; not measured data".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let surface = |name: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: name.into(),
            phase: "solid".into(),
            process: "analytical fixture".into(),
            revision: 0,
        },
        texture_frame: "native face".into(),
    };
    let card = InterfaceSystemCard::assemble(
        surface("body"),
        surface("belt"),
        SystemContext {
            medium: "dry".into(),
            third_body: None,
            environment: "fixture".into(),
            history: "explicit distance".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    let coefficient =
        ScalarPropertyRequirement::try_new("mu", Dims::NONE, ScalarAdmissibility::NonNegative)
            .unwrap();
    SourceThermalCoulomb::resolve_with_query_axes(
        bind_dry_interface_system_card(&card, InputAuthority::SyntheticFixture).unwrap(),
        point,
        coefficient.clone(),
        coefficient,
        MaterialPropertySelection::SingleClaimOnly,
        ThermalFrictionAxes {
            temperature: "T".into(),
            nominal_pressure: None,
        },
        axes,
    )
    .unwrap()
}

#[test]
#[allow(clippy::too_many_lines)]
fn g1_g4_native_solid_heat_feedback_history_and_atomic_refusal() {
    let mut mesh = Mesh2::quads(1.0, 1.0, 1, 1);
    // A valid trapezoid with only node zero inside the plane's active distance.
    mesh.nodes[2][0] = 0.2;
    // Rotate, preserving orientation, so the solid is thermal side A (x < 0).
    for node in &mut mesh.nodes {
        node[0] = -node[0];
        node[1] = -node[1];
    }
    let elastic = Hyperelastic::new(
        HyperelasticModel::NeoHookean {
            mu: 1.0,
            lambda: 3.0,
        },
        3.0,
    )
    .unwrap();
    let hyper = HyperProblem {
        mesh: &mesh,
        material: &elastic,
        dirichlet: vec![],
        traction: vec![],
        settings: NewtonSettings::default(),
    };
    let floor = HalfPlane {
        normal: [-1.0, 0.0],
        offset: -0.02,
    };
    let problem = ContactProblem {
        hyper: &hyper,
        sdf: &floor,
        barrier: Barrier {
            kappa: 10.0,
            dhat: 0.05,
        },
        friction: Some(Friction {
            mu: 0.1,
            eps_v: 0.001,
            rounds: 20,
        }),
        pins: (0..8).map(|i| (i, 0.0)).collect(),
        max_newton: 100,
        tol: 1e-9,
    };
    let initial = vec![0.0; 8];
    let accepted =
        ContactTrajectory::new_with_sliding_history(&problem, 0.0, &initial, 0.001, &[0.4; 4])
            .unwrap();
    let thermal = network(1, [1.0; 2], [1.0; 2], [1.0; 2], [0.0; 2], [0.0; 2]);
    let route = FrictionHeatRoute::from_link(&thermal, CONTACT, [1.0; 2]).unwrap();
    let law = friction_source(false, false);
    let contact = SolidContactHeatContact {
        node: 0,
        friction: ContactFrictionBinding {
            law: &law,
            conditions: law.initial_source().query(),
            nominal_contact_width_m: None,
        },
        route: &route,
        heat_fraction: 0.75,
        temperature_guess_k: 300.0,
    };
    let controls = SolidContactHeatControls {
        maximum_cells: 2,
        maximum_interfaces: 1,
        maximum_unknowns: 3,
        maximum_trials: 512,
        maximum_newton_iterations: 30,
        difference_step_j: 1e-6,
        temperature_scale_j_per_k: 1.0,
        heat: tolerances(),
        port_energy_tolerance_j: 2e-10,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 3,
            forcing_minimum: 1e-6,
            ..Default::default()
        },
        audit_entropy: true,
    };
    let input = |ids: &[u64]| SolidContactHeatSource {
        interfaces: ids
            .iter()
            .map(|&link_id| SolidContactHeatConductivity {
                link_id,
                conductivity_w_m_k: [1.0; 2],
            })
            .collect(),
        transport: StepCoefficients {
            sources_w: BTreeMap::new(),
            conductances_w_k: BTreeMap::new(),
            payload: (),
        },
    };
    let displacement = [[0.0, 0.1]; 4];
    let depth = 2.0;
    let gate = CancelGate::new_clock_free();
    let solved = propose_solid_contact_heat(
        &accepted,
        &thermal,
        1.0,
        &displacement,
        1.0,
        depth,
        &[contact.clone()],
        controls,
        &[0.0; 2],
        |_, ids, _| Ok::<_, &'static str>(input(ids)),
        &gate,
    )
    .unwrap();
    // Each body's C=1 J/K and half-face conductance G=2 W/K give
    // delta(T_body)=(2/3)*delta(T_interface). With mu=0.8-0.5*delta(T_i),
    // the independent scalar balance is (4/3 + 0.5*a)*delta(T_i)=0.8*a.
    let normal = -problem.barrier.d1(0.02);
    let a = contact.heat_fraction * depth * normal * 0.1;
    let rise = 0.8 * a / (4.0 / 3.0 + 0.5 * a);
    let mu = 0.8 - 0.5 * rise;
    let work_j = depth * mu * normal * 0.1;
    close(
        solved.heat.interfaces[0].interface_temperature_k,
        300.0 + rise,
        2e-9,
    );
    for id in [CellId(1), CellId(2)] {
        close(
            solved.thermal.temperature_k(id).unwrap(),
            300.0 + 2.0 * rise / 3.0,
            2e-9,
        );
    }
    let sample = solved.mechanical.friction_sources[0].as_ref().unwrap();
    close(sample.coefficient, mu, 2e-10);
    assert!(sample.coefficient < 0.8);
    close(sample.source.query().axes()["T"], 300.0 + rise, 2e-9);
    close(sample.source.query().axes()["path"], 0.5, 1e-14);
    close(solved.trajectory.sliding_distance()[0], 0.5, 1e-14);
    assert_eq!(&solved.trajectory.sliding_distance()[1..], &[0.4; 3]);
    close(
        solved.mechanical.friction_work_per_node[0] * depth,
        work_j,
        2e-10,
    );
    close(
        solved.mechanical.friction_reactions[0][1],
        mu * normal,
        2e-10,
    );
    close(solved.heat.dissipated_j, 0.75 * work_j, 2e-10);
    close(solved.exported_nonthermal_work_j, 0.25 * work_j, 2e-10);
    close(
        solved.port_energy_residual_j,
        0.0,
        controls.port_energy_tolerance_j,
    );
    assert!(solved.heat.entropy.is_some());
    assert!(solved.mechanical_solves <= solved.trials);
    assert_eq!(solved.trajectory.time_s(), 1.0);
    assert_eq!(solved.thermal.time_s(), 1.0);

    // A later hot trial fails only after mechanics has already proposed motion
    // and distance. Neither original owner may publish any part of that trial.
    let calls = Cell::new(0);
    let failed = propose_solid_contact_heat(
        &accepted,
        &thermal,
        1.0,
        &displacement,
        1.0,
        depth,
        &[contact.clone()],
        controls,
        &[0.0; 2],
        |_, ids, temperatures| {
            calls.set(calls.get() + 1);
            if temperatures[0] > 300.02 {
                Err("late thermal source refusal")
            } else {
                Ok(input(ids))
            }
        },
        &gate,
    );
    assert!(matches!(
        failed,
        Err(SolidContactHeatError::Source("late thermal source refusal"))
    ));
    assert!(calls.get() > 1);
    assert_eq!(accepted.time_s(), 0.0);
    assert_eq!(accepted.displacement(), initial);
    assert_eq!(accepted.sliding_distance(), &[0.4; 4]);
    assert_eq!(accepted.work().relative, 0.0);
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.temperature_k(CellId(1)).unwrap(), 300.0);

    // Both declared evolving coordinates and otherwise legal static external
    // temperatures need a thermal body/side owner in this coupled entry.
    for declared in [true, false] {
        let unmapped_temperature = friction_source(true, declared);
        let mut other = contact.clone();
        other.friction.law = &unmapped_temperature;
        other.friction.conditions = unmapped_temperature.initial_source().query();
        assert!(matches!(
            propose_solid_contact_heat(
                &accepted,
                &thermal,
                1.0,
                &displacement,
                1.0,
                depth,
                &[other],
                controls,
                &[0.0; 2],
                |_, ids, _| Ok::<_, &'static str>(input(ids)),
                &gate
            ),
            Err(SolidContactHeatError::Invalid(
                "additional source temperatures require explicit thermal ownership"
            ))
        ));
    }

    let mut inactive = contact.clone();
    inactive.node = 1;
    assert!(matches!(
        propose_solid_contact_heat(
            &accepted,
            &thermal,
            1.0,
            &displacement,
            1.0,
            depth,
            &[inactive],
            controls,
            &[0.0; 2],
            |_, ids, _| Ok::<_, &'static str>(input(ids)),
            &gate
        ),
        Err(SolidContactHeatError::Invalid(
            "fixed heat footprint is not initially in contact"
        ))
    ));

    // Declaring all loss outside the thermal control volume still advances the
    // mechanical history, without silently converting that loss into heat.
    let mut outside = contact;
    outside.heat_fraction = 0.0;
    let cold_sources = Cell::new(0usize);
    let cold = propose_solid_contact_heat(
        &accepted,
        &thermal,
        1.0,
        &displacement,
        1.0,
        depth,
        &[outside],
        controls,
        &[0.0; 2],
        |_, ids, _| {
            cold_sources.set(cold_sources.get() + 1);
            let data = input(ids);
            Ok::<_, &'static str>(SolidContactHeatSource {
                interfaces: data.interfaces,
                transport: StepCoefficients {
                    sources_w: data.transport.sources_w,
                    conductances_w_k: data.transport.conductances_w_k,
                    payload: cold_sources.get(),
                },
            })
        },
        &gate,
    )
    .unwrap();
    assert_eq!(cold.heat.dissipated_j, 0.0);
    assert_eq!(cold.thermal.temperature_k(CellId(1)).unwrap(), 300.0);
    close(
        cold.exported_nonthermal_work_j,
        depth * 0.8 * normal * 0.1,
        2e-10,
    );
    close(cold.trajectory.sliding_distance()[0], 0.5, 1e-14);
    // Initial residual and final entropy audit have identical interface T.
    // They reuse mechanics, but each resolves its own exact thermal payload.
    assert_eq!(cold.mechanical_solves, 1);
    assert!(cold.trials > cold.mechanical_solves);
    assert_eq!(cold_sources.get(), cold.trials);
    assert_eq!(cold.source, cold_sources.get());
}
