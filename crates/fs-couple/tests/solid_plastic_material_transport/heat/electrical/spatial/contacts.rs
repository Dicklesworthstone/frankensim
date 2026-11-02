//! Matching electrical traces on prescribed native solid kinematics, not a
//! mechanical contact/tie solver. All material/interface laws are manufactured.
#[path = "contacts/pressure.rs"]
mod pressure;
#[path = "contacts/thermal.rs"]
mod thermal;
use super::terminals::{PRISM_NODES, PRISM_TETS, current_drive, prism_thermal_cells_with_links};
use super::*;
use fs_couple::circuit_enthalpy::{
    conduction::{contacts::ConductionContactThermal, spatial::SpatialConductionThermalStep},
    contacts::heat::HeatNetworkSample,
};
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;
use fs_material::electrical_interface::resolve_electrical_interface_state;

const IDS: [CellId; 6] = [
    CellId(100),
    CellId(101),
    CellId(102),
    CellId(103),
    CellId(104),
    CellId(105),
];
const BULK_RHO: f64 = 0.1;
const AREA: f64 = L * L / 2.0;
const CURRENT: f64 = 0.015;
const FRACTION: f64 = 0.25;
const TA_SLOPE: f64 = 0.003;
const TB_SLOPE: f64 = 0.007;
const CONTACT_LINK: u64 = 200;

fn contact_controls() -> SpatialElectricalSolidSolve {
    let mut controls = spatial_controls();
    controls.coupled.coupled.maximum_cells = IDS.len();
    controls
}

fn specimen_mesh() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let mut nodes = PRISM_NODES.to_vec();
    nodes.extend(PRISM_NODES.map(|[x, y, z]| [x + L, y, z]));
    let mut tets = PRISM_TETS.to_vec();
    tets.extend(PRISM_TETS.map(|tet| tet.map(|v| v + PRISM_NODES.len())));
    (nodes, tets)
}

fn interface(card: &MaterialCard, r0: f64) -> InterfaceSystemCard {
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let mut claim = scalar(ER, ER_DIMS, r0);
    claim.value = PropertyValue::Surface2D {
        axes: ["Ta".into(), "Tb".into()],
        axis_dims: [Temperature::DIMS; 2],
        samples: [
            [280.0, 280.0],
            [320.0, 280.0],
            [280.0, 320.0],
            [320.0, 320.0],
        ]
        .map(|p| {
            (
                p,
                r0 * (1.0 + TA_SLOPE * (p[0] - 300.0) + TB_SLOPE * (p[1] - 300.0)),
            )
        })
        .to_vec(),
        triangles: vec![[0, 1, 2], [1, 3, 2]],
        dims: ER_DIMS,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim.validity = ValidityDomain::unconstrained()
        .with_quantity("Ta", temperature, 280.0, 320.0)
        .with_quantity("Tb", temperature, 280.0, 320.0)
        .with_quantity("p", pressure_quantity(), 199999.0, 200001.0);
    let mut claims = ClaimSet::new();
    claims.insert_claim(claim).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: card.id().clone(),
            texture_frame: "left face".into(),
        },
        SurfaceSpec {
            material: card.id().clone(),
            texture_frame: "right face".into(),
        },
        SystemContext {
            medium: "prescribed matching interface".into(),
            third_body: None,
            environment: "authored pressure, independent of solid stress".into(),
            history: "fixed preparation".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn with_contact<R>(
    r0: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>,
        &SpatialConduction,
        &PreparedHeatNetwork<'_, '_>,
        &[ConductionContactThermal<'_>],
        &[Vec<usize>; 2],
    ) -> R,
) -> R {
    let f = fixture_with_property(
        EMISSIVITY,
        Some(scalar(
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            BULK_RHO,
        )),
    );
    // The original second card has fixed 40 MPa yield, E, nu and Cp. Both
    // prisms stay elastic with alpha explicitly zero; only interface R'' varies.
    let interface = interface(&f.cards[1], r0);
    let temperature = QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ));
    let point = QueryPoint::new()
        .with_quantity("Ta", temperature, 300.0)
        .unwrap()
        .with_quantity("Tb", temperature, 300.0)
        .unwrap()
        .with_quantity("p", pressure_quantity(), PRESSURE_PA)
        .unwrap();
    with_contact_fixture(
        ContactFixture {
            material: f,
            interface,
            point,
            fraction: FRACTION,
            mean_temperature_axis: None,
            expansion: false,
            thermal_contact_area_m2: None,
        },
        run,
    )
}

// Both fixed-pressure and solved-pressure tests use this original native
// geometry/caloric/material assembly; only their source laws/owned roles differ.
struct ContactFixture {
    material: Fixture,
    interface: InterfaceSystemCard,
    point: QueryPoint,
    fraction: f64,
    mean_temperature_axis: Option<&'static str>,
    expansion: bool,
    thermal_contact_area_m2: Option<f64>,
}

fn with_contact_fixture<R>(
    fixture: ContactFixture,
    run: impl FnOnce(
        &Cx<'_>,
        &DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>,
        &SpatialConduction,
        &PreparedHeatNetwork<'_, '_>,
        &[ConductionContactThermal<'_>],
        &[Vec<usize>; 2],
    ) -> R,
) -> R {
    let ContactFixture {
        material: f,
        interface,
        point,
        fraction,
        mean_temperature_axis,
        expansion,
        thermal_contact_area_m2,
    } = fixture;
    let (nodes, tets) = specimen_mesh();
    let material = &f.cards[1];
    let caloric = &f.caloric[1];
    let bindings = [ConductionContactThermal {
        contact: 0,
        point: ThermalInterfacePoint {
            interface: &interface,
            materials: [material, material],
            context: &point,
            temperature_axes: ["Ta", "Tb"],
            liquid_fraction_axes: [None, None],
            mean_temperature_axis,
        },
        heat_fraction_to_first: fraction,
    }];
    with_cx(|_, cx| {
        let links: Vec<_> = thermal_contact_area_m2
            .into_iter()
            .map(|area_m2| ThermalLink {
                id: CONTACT_LINK,
                a: IDS[2],
                b: IDS[3],
                area_m2,
                conductance_w_k: 0.0,
            })
            .collect();
        let thermal =
            prism_thermal_cells_with_links(caloric.curve(), MATS[1], &nodes, &tets, &IDS, &links);
        let solid = ThermoplasticSolid::try_new(
            &thermal,
            &nodes,
            &tets,
            IDS.iter()
                .map(|&cell| ThermoplasticCell {
                    cell,
                    query: SourceQuery {
                        card: material,
                        caloric,
                        context: f.point.clone(),
                        temperature_axis: "T".into(),
                        selection: MaterialPropertySelection::SingleClaimOnly,
                        hardening_basis: BASIS,
                        thermal_expansion: expansion
                            .then_some(MaterialPropertySelection::SingleClaimOnly),
                        strain_limit: 0.04,
                    },
                })
                .collect(),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let selections = vec![
            thermal_contact_area_m2.map(|_| (
                f.conductivity.clone(),
                MaterialPropertySelection::SingleClaimOnly,
            ));
            IDS.len()
        ];
        let cells = solid.thermal_material_cells(&selections).unwrap();
        let conduction: Vec<_> = links
            .iter()
            .map(|link| ConductionLink {
                link: link.id,
                distance_a_m: L / 4.0,
                distance_b_m: L / 4.0,
            })
            .collect();
        let contacts: Vec<_> = links
            .iter()
            .map(|link| MaterialContactLink {
                link: link.id,
                sides: [link.a, link.b].map(|cell| ContactCell {
                    cell,
                    material_id: MATS[1],
                    curve: caloric.curve(),
                }),
                point: bindings[0].point,
                selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
            })
            .collect();
        // Optional original reference-area Fourier/contact path. Each side
        // keeps its finite caloric cell, even with no authored heat links.
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &conduction,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let network = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &contacts,
            radiation: None,
        };
        let prepared = network.prepare(cx, &thermal).unwrap();
        let conductor = resolve_conductor_state_with_basis(
            material,
            &f.point,
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap();
        let source = resolve_electrical_interface_state(
            &interface,
            &point,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let field = SpatialConduction::new_with_contacts(
            nodes.clone(),
            tets.iter().map(|t| t.map(|v| v as u32)).collect(),
            vec![conductor; IDS.len()],
            vec![SpatialConductionContact {
                first_vertices: [3, 4, 5],
                second_vertices: [6, 7, 8],
                source,
            }],
            contact_controls().assembly,
            &mut || false,
        )
        .unwrap();
        let electrodes: [Vec<usize>; 2] = [2.0 * L, 0.0].map(|x| {
            field
                .exterior_faces()
                .filter(|&face| {
                    field.mesh().faces[face]
                        .iter()
                        .all(|&v| nodes[v as usize][0] == x)
                })
                .collect()
        });
        run(
            cx,
            &thermal,
            &solid,
            &field,
            &prepared,
            &bindings,
            &electrodes,
        )
    })
}

fn prescribed(stretch: f64, translation: f64) -> Vec<(usize, f64)> {
    specimen_mesh()
        .0
        .iter()
        .enumerate()
        .flat_map(|(node, x)| {
            (0..3).map(move |axis| {
                (
                    3 * node + axis,
                    if axis == 0 {
                        translation
                    } else {
                        (stretch - 1.0) * x[axis]
                    },
                )
            })
        })
        .collect()
}

fn endpoint(old: [f64; 6], r0: f64, stretch: f64, current: f64, dt: f64) -> [f64; 6] {
    // Two independent finite-capacity balances eliminate the affine two-sided
    // contact law. Other cells receive only their original bulk Joule power.
    let capacity = MASS * 5.0;
    let area = AREA * stretch * stretch;
    let bulk = dt * current * current * BULK_RHO * L / (3.0 * area);
    let contact = dt * current * current * r0 / area;
    let base = 1.0 + TA_SLOPE * (old[2] - 300.0) + TB_SLOPE * (old[3] - 300.0);
    let a = capacity - FRACTION * contact * TA_SLOPE;
    let b = -FRACTION * contact * TB_SLOPE;
    let c = -(1.0 - FRACTION) * contact * TA_SLOPE;
    let d = capacity - (1.0 - FRACTION) * contact * TB_SLOPE;
    let rhs_a = bulk + FRACTION * contact * base;
    let rhs_b = bulk + (1.0 - FRACTION) * contact * base;
    let mut result = old.map(|t| t + bulk / capacity);
    result[2] = old[2] + (rhs_a * d - b * rhs_b) / (a * d - b * c);
    result[3] = old[3] + (a * rhs_b - c * rhs_a) / (a * d - b * c);
    result
}

fn check_endpoint(
    thermal: &DepositionSimulation<'_>,
    problem: &SpatialConduction,
    report: &ElectricalSolidReport<HeatNetworkSample, SpatialConductionThermalStep>,
    old: [f64; 6],
    r0: f64,
    stretch: f64,
    current: f64,
    dt: f64,
) -> [f64; 6] {
    let expected = endpoint(old, r0, stretch, current, dt);
    let actual = IDS.map(|cell| thermal.temperature_k(cell).unwrap());
    for i in 0..IDS.len() {
        close(actual[i], expected[i], 3e-6);
        close(thermal.cells()[&IDS[i]].mass_kg, MASS, 0.0);
    }
    let area = AREA * stretch * stretch;
    let bulk_r = BULK_RHO * L / area;
    let interface_r = r0 * (1.0 + TA_SLOPE * (actual[2] - 300.0) + TB_SLOPE * (actual[3] - 300.0));
    let contact_r = interface_r / area;
    let step = &report.circuit_thermal.electrical;
    let field = &step.field;
    close(
        field.voltage_v[9],
        current * (2.0 * bulk_r + contact_r),
        1e-9,
    );
    for i in 0..IDS.len() {
        close(field.current_density_a_m2[i][0], -current / area, 1e-7);
        close(
            field.cell_power_w[i],
            current * current * bulk_r / 3.0,
            1e-10,
        );
    }
    let contact = &step.contacts[0];
    assert_eq!(step.contacts.len(), 1);
    assert_eq!(contact.cells, [IDS[2], IDS[3]]);
    close(contact.heat_fraction_to_first, FRACTION, 0.0);
    close(contact.response.current_a, -current, 1e-10);
    close(
        contact.response.power_w,
        current * current * contact_r,
        1e-10,
    );
    for jump in contact.response.voltage_jump_v {
        close(jump, -current * contact_r, 1e-9);
    }
    let source = &problem.contacts()[0].source;
    close(source.resistance_ohm_m2(), interface_r, 1e-12);
    close(
        source.resolved().state_coordinate("Ta").unwrap(),
        actual[2],
        1e-10,
    );
    close(
        source.resolved().state_coordinate("Tb").unwrap(),
        actual[3],
        1e-10,
    );
    close(
        source.resolved().state_coordinate("p").unwrap(),
        PRESSURE_PA,
        0.0,
    );
    // Matched coordinates are retained exactly, while their actual metric grows.
    for (a, b) in [(3, 6), (4, 7), (5, 8)] {
        assert_eq!(problem.positions_m()[a], problem.positions_m()[b]);
    }
    let p = problem.positions_m();
    close((p[4][1] - p[3][1]) * (p[5][2] - p[3][2]) / 2.0, area, 1e-15);
    let heat = dt * current * current * (2.0 * bulk_r + contact_r);
    close(step.heat_j, heat, 1e-9);
    close(step.supplied_work_j, heat, 1e-9);
    close(report.circuit_thermal.deposited_joule_j, heat, 1e-9);
    close(
        report.circuit_thermal.thermal.stored_energy_change_j,
        heat,
        1e-9,
    );
    close(report.deposited_mechanical_heat_j, 0.0, 1e-10);
    close(report.combined_residual_j, 0.0, 1e-8);
    actual
}

#[test]
fn matched_surface_area_and_interface_only_swap_change_native_two_sided_heat() {
    let mut contact_temperatures = Vec::new();
    for r0 in [0.001, 0.002] {
        with_contact(
            r0,
            |cx, thermal, solid, problem, prepared, bindings, electrodes| {
                let first = solid
                    .advance_spatial_electrical_contact_candidate(
                        thermal,
                        problem,
                        prepared,
                        current_drive(electrodes, CURRENT),
                        bindings,
                        1.0,
                        &prescribed(1.01, 0.0),
                        &[],
                        IncrementSettings::default(),
                        contact_controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                let temperatures = check_endpoint(
                    first.thermal(),
                    first.problem(),
                    first.report(),
                    [300.0; 6],
                    r0,
                    1.01,
                    CURRENT,
                    1.0,
                );
                contact_temperatures.push(temperatures[3]);
                let (thermal, solid, problem, _) = first.into_parts();
                let next = solid
                    .advance_spatial_electrical_contact_candidate(
                        &thermal,
                        &problem,
                        prepared,
                        current_drive(electrodes, -CURRENT),
                        bindings,
                        1.25,
                        &prescribed(1.02, 0.0),
                        &[],
                        IncrementSettings::default(),
                        contact_controls(),
                        &[0.0; 6],
                        cx,
                    )
                    .unwrap();
                check_endpoint(
                    next.thermal(),
                    next.problem(),
                    next.report(),
                    temperatures,
                    r0,
                    1.02,
                    -CURRENT,
                    0.25,
                );
                assert_eq!(solid.time_s(), 1.0);
                assert_eq!(thermal.time_s(), 1.0);
            },
        );
    }
    assert!(contact_temperatures[1] > contact_temperatures[0] + 4.0);
}

#[test]
fn inertial_contact_continuation_retains_heat_partition_and_independent_kinetic_work() {
    with_contact(
        0.001,
        |cx, thermal, solid, problem, prepared, bindings, electrodes| {
            let first = solid
                .advance_spatial_electrical_contact_candidate(
                    thermal,
                    problem,
                    prepared,
                    current_drive(electrodes, CURRENT),
                    bindings,
                    1.0,
                    &prescribed(1.01, 0.0),
                    &[],
                    IncrementSettings::default(),
                    contact_controls(),
                    &[0.0; 6],
                    cx,
                )
                .unwrap();
            let old = IDS.map(|cell| first.thermal().temperature_k(cell).unwrap());
            let (thermal, solid, problem, _) = first.into_parts();
            let velocity = 0.1;
            let dt = 0.1;
            let acceleration = 1.0;
            let next_velocity = velocity + dt * acceleration;
            let displacement = dt * next_velocity;
            let moving =
                ThermoplasticDynamics::try_new(solid, vec![[velocity, 0.0, 0.0]; 12], cx).unwrap();
            let next = moving
                .advance_spatial_electrical_contact_candidate(
                    &thermal,
                    &problem,
                    prepared,
                    current_drive(electrodes, -CURRENT),
                    bindings,
                    1.0 + dt,
                    &prescribed(1.01, displacement),
                    &[],
                    IncrementSettings::default(),
                    contact_controls(),
                    &[0.0; 6],
                    cx,
                )
                .unwrap();
            check_endpoint(
                next.thermal(),
                next.problem(),
                &next.report().material,
                old,
                0.001,
                1.01,
                -CURRENT,
                dt,
            );
            for v in next.dynamics().velocity_mps() {
                close(v[0], next_velocity, 1e-12);
                close(v[1], 0.0, 1e-12);
                close(v[2], 0.0, 1e-12);
            }
            let total_mass = 6.0 * MASS;
            close(
                next.report().inertia.kinetic_change_j,
                0.5 * total_mass * (next_velocity * next_velocity - velocity * velocity),
                1e-11,
            );
            close(
                next.report().inertia.numerical_loss_j,
                0.5 * total_mass * (dt * acceleration).powi(2),
                1e-11,
            );
            close(
                next.report().material.mechanical_work_j,
                total_mass * acceleration * displacement,
                1e-10,
            );
            assert!(next.report().inertia.numerical_loss_j > 1e-6);
        },
    );
}

#[test]
fn unowned_contact_and_opened_matching_trace_refuse_without_any_publication() {
    with_contact(
        0.001,
        |cx, thermal, solid, problem, prepared, bindings, electrodes| {
            let cells = thermal.cells().clone();
            let solution = solid.solution().clone();
            let sources = problem.contacts().to_vec();
            let positions = problem.positions_m().to_vec();
            let unowned = solid.advance_spatial_electrical_candidate(
                thermal,
                problem,
                prepared,
                current_drive(electrodes, CURRENT),
                1.0,
                &prescribed(1.01, 0.0),
                &[],
                IncrementSettings::default(),
                contact_controls(),
                &[0.0; 6],
                cx,
            );
            assert!(unowned.is_err());
            let mut opened = prescribed(1.01, 0.0);
            for (dof, displacement) in &mut opened {
                if *dof >= 18 && *dof % 3 == 0 {
                    *displacement += 1e-7;
                }
            }
            let refused = solid.advance_spatial_electrical_contact_candidate(
                thermal,
                problem,
                prepared,
                current_drive(electrodes, CURRENT),
                bindings,
                1.0,
                &opened,
                &[],
                IncrementSettings::default(),
                contact_controls(),
                &[0.0; 6],
                cx,
            );
            let error = refused
                .err()
                .expect("an opened trace cannot retain a zero-thickness contact");
            assert!(
                format!("{error:?}")
                    .contains("paired contact vertices must have identical physical coordinates")
            );
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(solid.solution(), &solution);
            assert_eq!(problem.contacts(), &sources);
            assert_eq!(problem.positions_m(), &positions);
            assert_eq!(solid.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
        },
    );
}
