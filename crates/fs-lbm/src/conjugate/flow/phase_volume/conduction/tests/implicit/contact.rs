//! Original phase/pressure contact data in the finite-fluid/wall heat solve.

use super::*;
use crate::conjugate::FaceContactResistance;
use fs_conduction::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLaw, ThermalLink, TimedEvent,
};
use fs_conduction_full::transient::deposition::contact::implicit::PreparedMaterialContact;
use fs_conduction_full::transient::deposition::contact::{ContactCell, MaterialContactLink};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermal::interface::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    ThermalInterfacePoint,
};
use fs_qty::Pressure;
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};

// The contact response gate admits a fixed external coordinate only when it is typed (`is_fixed_external_coordinate`,
// fs-material state_point.rs, d21f10ea46, 2026-10-22): a pressure axis declared as a plain dimensional pascal is refused
// ("contact response requires declared thermal/phase or typed fixed external coordinates"), as the passing fs-couple
// fixtures do by declaring `typed(QuantityKind::Pressure)`.
fn typed_pressure() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static))
}

fn interface(wall: &MaterialCard, fluid: &MaterialCard) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            // R''=.5+f_wall-p/(4e5). The mechanical pressure stays fixed;
            // the original latent fraction changes in the implicit residual.
            value: PropertyValue::Surface2D {
                axes: ["fa".into(), "p".into()],
                axis_dims: [Dims::NONE, Pressure::DIMS],
                samples: [[0.0, 1e5], [1.0, 1e5], [0.0, 1.5e5], [1.0, 1.5e5]]
                    .into_iter()
                    .map(|at| (at, 0.5 + at[0] - at[1] / 4e5))
                    .collect(),
                triangles: vec![[0, 1, 2], [1, 3, 2]],
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with_quantity(
                    "TA",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    400.0,
                )
                .with_quantity(
                    "TB",
                    QuantitySpec::dimensional(Temperature::DIMS),
                    300.0,
                    400.0,
                )
                .with_quantity("fa", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
                .with_quantity("fb", QuantitySpec::dimensional(Dims::NONE), 0.0, 1.0)
                .with_quantity("p", typed_pressure(), 1e5, 1.5e5),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured finite-wall contact calorimetry".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: wall.id().clone(),
            texture_frame: "wall A".into(),
        },
        SurfaceSpec {
            material: fluid.id().clone(),
            texture_frame: "fluid B".into(),
        },
        SystemContext {
            medium: "thermal contact".into(),
            third_body: None,
            environment: "sealed calorimetry".into(),
            history: "fixed externally supplied compression".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

/// Prepare the source with an actual native thermal inventory using the SAME
/// original curves, initial states, masses and graded geometry as the bridge.
/// Its original interface order is wall A (voxel 1), fluid B (voxel 0).
fn reference<'a>(
    fluid: &'a EquilibriumEnthalpyPhaseCurve,
    wall: &'a EquilibriumEnthalpyPhaseCurve,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (index, (mass, h, start, end, volume, area)) in [
        (2.0, 51.0, 0.0, 1.0, 1.0, 6.0),
        (4.0, 25.0, 1.0, 3.0, 2.0, 10.0),
    ]
    .into_iter()
    .enumerate()
    {
        let id = index as u64;
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(id),
                segment_id: id,
                body_id: 0,
                material_id: MaterialId(id),
                mass_kg: mass,
                specific_h_j_kg: h,
                geometry: CellGeometry {
                    start_m: [start, 0.5, 0.5],
                    end_m: [end, 0.5, 0.5],
                    volume_m3: volume,
                    surface_area_m2: area,
                },
                boundary: Boundary {
                    ambient_k: [351.0, 320.0][index],
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: 0,
            a: CellId(1),
            b: CellId(0),
            area_m2: 1.0,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(0), fluid as &dyn ThermalLaw),
            (MaterialId(1), wall as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 1,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 0.0,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn binding<'a>(
    cards: [&'a MaterialCard; 2],
    curves: [&'a EquilibriumEnthalpyPhaseCurve; 2],
    interface: &'a InterfaceSystemCard,
    context: &'a QueryPoint,
) -> MaterialContactLink<'a> {
    MaterialContactLink {
        link: 0,
        sides: [1, 0].map(|side| ContactCell {
            cell: CellId(side),
            material_id: MaterialId(side),
            curve: curves[side as usize],
        }),
        point: ThermalInterfacePoint {
            interface,
            materials: [cards[1], cards[0]],
            context,
            temperature_axes: ["TA", "TB"],
            liquid_fraction_axes: [Some("fa"), Some("fb")],
            mean_temperature_axis: None,
        },
        selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
    }
}

#[test]
fn source_phase_contact_preserves_order_pressure_and_finite_wall_energy() {
    let fluid_card = card("T", &[(300.0, 2.0), (400.0, 2.0)]);
    let fluid_curve = curve(
        &fluid_card,
        &[(0.0, 300.0, 1.0, 2.0), (100.0, 400.0, 1.0, 2.0)],
    );
    let fluid = source(&fluid_card, &fluid_curve, [0.0, 100.0]);
    let wall_card = card("f", &[(0.0, 4.0), (1.0, 4.0)]);
    let wall_curve = curve(
        &wall_card,
        &[
            (-20.0, 300.0, 0.0, 2.0),
            (0.0, 320.0, 0.0, 2.0),
            (100.0, 320.0, 1.0, 2.0),
            (180.0, 400.0, 1.0, 2.0),
        ],
    );
    let wall = source(&wall_card, &wall_curve, [-20.0, 180.0]);
    let old = reference(&fluid_curve, &wall_curve);
    let interface = interface(&wall_card, &fluid_card);
    let context = QueryPoint::new()
        .with_quantity("TA", QuantitySpec::dimensional(Temperature::DIMS), 320.0)
        .unwrap()
        .with_quantity("TB", QuantitySpec::dimensional(Temperature::DIMS), 351.0)
        .unwrap()
        .with_quantity("fa", QuantitySpec::dimensional(Dims::NONE), 0.25)
        .unwrap()
        .with_quantity("fb", QuantitySpec::dimensional(Dims::NONE), 1.0)
        .unwrap()
        .with_quantity("p", typed_pressure(), 1e5)
        .unwrap();
    let binding = binding(
        [&fluid_card, &wall_card],
        [&fluid_curve, &wall_curve],
        &interface,
        &context,
    );
    let contact =
        PreparedMaterialContact::new(&old, binding, [[-20.0, 180.0], [0.0, 100.0]], || Ok(()))
            .unwrap();
    let selected = contact.reference().interface.receipt().selected;
    let sample = contact.state_at_specific_enthalpies([30.0, 41.0]).unwrap();
    close(sample.resistance_m2_k_w, 0.55);
    close(sample.resistance_enthalpy_derivatives_m2_kg_k_w_j[0], 0.01);
    assert_eq!(sample.resistance_enthalpy_derivatives_m2_kg_k_w_j[1], 0.0);
    assert!(contact.source_surface().is_some());

    let mut domain = VoxelDomain::graded([vec![1.0, 2.0], vec![1.0], vec![1.0]]).unwrap();
    domain.set(1, 0, 0, Voxel::Solid(0));
    let field = FlowField::quiescent(&domain);
    let phases = [
        Some(fluid_curve.state_at_specific_enthalpy(51.0).unwrap()),
        None,
    ];
    let mass = [2.0, 0.0];
    let heat = [0.0; 2];
    let data = input(&mass, &phases, &heat);
    let walls = [ThermalVolumeSolidCell {
        cell: 1,
        conductivity: &wall,
        mass_kg: 4.0,
        phase: wall_curve.state_at_specific_enthalpy(25.0).unwrap(),
        net_enthalpy_increment_j: 0.0,
    }];
    let contacts = [ThermalVolumeContact::Source {
        cells: [1, 0],
        law: &contact,
    }];
    let gate = CancelGate::new_clock_free();
    let mut workspace = ImplicitThermalVolumeWorkspace::new(2).unwrap();
    let run = |contacts: &[ThermalVolumeContact<'_>], workspace: &mut _, gate: &_| {
        field.advance_isotropic_implicit_conjugate_thermal_volume_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            &fluid,
            ThermalVolumeBoundaries::default(),
            &walls,
            contacts,
            data,
            controls(1.0),
            solve_controls(),
            workspace,
            gate,
        )
    };
    let result = run(&contacts, &mut workspace, &gate).unwrap();
    // Half-cell resistance=.5/2+1/4=.5. At Q*=20, h_fluid*=41,
    // h_wall*=30, f_wall*=.30, R''*=.55. Q=(341-320)/1.05=20.
    close(result.conductive_enthalpy_increment_j[0], -20.0);
    assert_eq!(result.conductive_enthalpy_increment_j[1], 0.0);
    close(
        result.thermal.volume.phase_states[0]
            .unwrap()
            .specific_enthalpy_j_kg(),
        41.0,
    );
    close(
        result.solid_endpoints[0].phase.specific_enthalpy_j_kg(),
        30.0,
    );
    close(result.solid_endpoints[0].phase.temperature_k(), 320.0);
    close(result.solid_endpoints[0].phase.liquid_mass_fraction(), 0.30);
    close(
        result.solid_endpoints[0].conductive_enthalpy_increment_j,
        20.0,
    );
    close(result.report.solid_to_fluid_heat_j, -20.0);
    close(result.report.conduction.solve.energy.external_heat_j, 0.0);
    close(
        result.thermal.volume.report.material_enthalpy_change_j
            + result.report.solid_enthalpy_change_j,
        0.0,
    );
    assert!(result.report.conduction.endpoint_heat_defect_j <= 1e-9);
    assert!(result.report.conduction.absolute_enthalpy_defect_j <= 1e-9);
    assert_eq!(result.thermal.volume.transport.mass_kg, mass);
    assert!(result.thermal.volume.phase_states[1].is_none());

    // Freezing R at its initial value gives Q=31/1.5, independently
    // distinguishing endpoint contact physics from held-coefficient coupling.
    let fixed = [ThermalVolumeContact::Resistance(FaceContactResistance {
        cells: [0, 1],
        resistance_m2_k_w: 0.5,
    })];
    let frozen = run(&fixed, &mut workspace, &gate).unwrap();
    close(
        frozen.solid_endpoints[0].conductive_enthalpy_increment_j,
        31.0 / 1.5,
    );
    assert!(frozen.solid_endpoints[0].conductive_enthalpy_increment_j > 20.5);

    // Fixed externally authored compression changes the source response; the
    // bridge never substitutes fluid projection pressure or evolves this axis.
    let compressed_context = context
        .clone()
        .with_quantity("p", typed_pressure(), 1.5e5)
        .unwrap();
    let compressed = PreparedMaterialContact::new(
        &old,
        MaterialContactLink {
            point: ThermalInterfacePoint {
                context: &compressed_context,
                ..binding.point
            },
            ..binding
        },
        [[-20.0, 180.0], [0.0, 100.0]],
        || Ok(()),
    )
    .unwrap();
    let higher_pressure = run(
        &[ThermalVolumeContact::Source {
            cells: [1, 0],
            law: &compressed,
        }],
        &mut workspace,
        &gate,
    )
    .unwrap();
    let expected = 24800.0 / (352100.0_f64.sqrt() + 550.0);
    close(
        higher_pressure.solid_endpoints[0].conductive_enthalpy_increment_j,
        expected,
    );
    assert!(expected > 20.0);

    let zero_momentum = [[0.0; 3]; 2];
    let momentum = field
        .advance_isotropic_implicit_conjugate_thermal_volume_momentum_with_baths(
            &domain,
            [FlowFace::Wall; 6],
            None,
            &fluid,
            ThermalVolumeBoundaries::default(),
            &walls,
            &contacts,
            ThermalVolumeMomentumInput {
                volume: data,
                momentum_kg_m_s: &zero_momentum,
                acceleration_m_s2: [0.0; 3],
                force_density_n_m3: &zero_momentum,
                boundary_velocity_m_s: [None; 6],
            },
            super::super::momentum::coupled_controls(1.0),
            solve_controls(),
            &mut workspace,
            &gate,
        )
        .unwrap();
    assert_eq!(momentum.evolution.thermal, result.thermal);
    assert_eq!(momentum.solid_endpoints, result.solid_endpoints);
    assert_eq!(momentum.report, result.report);
    assert_eq!(momentum.evolution.momentum.momentum_kg_m_s, zero_momentum);

    // Neither duplicate ownership nor a source-domain failure publishes a
    // partial wall. The same workspace replays the accepted original interval.
    assert!(run(&[contacts[0], fixed[0]], &mut workspace, &gate).is_err());
    let narrow =
        PreparedMaterialContact::new(&old, binding, [[-20.0, 29.0], [0.0, 100.0]], || Ok(()))
            .unwrap();
    assert!(
        run(
            &[ThermalVolumeContact::Source {
                cells: [1, 0],
                law: &narrow
            }],
            &mut workspace,
            &gate
        )
        .is_err()
    );
    gate.request();
    assert!(matches!(
        run(&contacts, &mut workspace, &gate),
        Err(ChtError::Cancelled)
    ));
    let fresh_gate = CancelGate::new_clock_free();
    assert_eq!(result, run(&contacts, &mut workspace, &fresh_gate).unwrap());
    assert_eq!(contact.reference().interface.receipt().selected, selected);
    assert_eq!(contact.reference().cells, [CellId(1), CellId(0)]);
    assert_eq!(walls[0].phase.specific_enthalpy_j_kg(), 25.0);
    assert_eq!(phases[0].unwrap().specific_enthalpy_j_kg(), 51.0);
    assert_eq!(field, FlowField::quiescent(&domain));
}
