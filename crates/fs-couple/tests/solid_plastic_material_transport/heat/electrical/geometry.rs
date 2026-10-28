//! G1/G3: actual native displacement changes the original source-bound current.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;
use fs_couple::solid_enthalpy::plastic::electrical::geometry::{
    ConductorGeometrySolve, DeformingConductorPath,
};
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
};

fn geometry_controls() -> ConductorGeometrySolve {
    ConductorGeometrySolve {
        maximum_displacement_gradient: 0.02,
        orthogonality_tolerance: 1e-10,
    }
}
fn path(direction: [f64; 3]) -> Vec<DeformingConductorPath> {
    vec![DeformingConductorPath {
        branch: 1,
        reference_direction: direction,
    }]
}
fn affine(gradient: [[f64; 3]; 3]) -> Vec<(usize, f64)> {
    NODES
        .iter()
        .enumerate()
        .flat_map(|(node, position)| {
            (0..3).map(move |axis| {
                (
                    3 * node + axis,
                    (0..3).map(|j| gradient[axis][j] * position[j]).sum(),
                )
            })
        })
        .collect()
}
fn free_z(nodes: &[usize]) -> Vec<(usize, f64)> {
    (0..3 * NODES.len())
        .filter(|dof| !(dof % 3 == 2 && nodes.contains(&(dof / 3))))
        .map(|dof| (dof, 0.0))
        .collect()
}
fn capacitor_charge(initial: f64, resistance: f64, dt: f64) -> f64 {
    let z = dt / (2.0 * resistance * CAPACITANCE);
    initial * (1.0 - z) / (1.0 + z)
}

#[test]
fn g1_affine_geometry_changes_current_and_heat_without_reauthoring_material_or_inventory() {
    with_powered(
        10.0,
        false,
        false,
        CHARGE,
        |cx, _, thermal, original, circuit, prepared| {
            let solid = original
                .clone()
                .with_deforming_conductors(path([1.0, 0.0, 0.0]), geometry_controls(), cx)
                .unwrap();
            let stretch = [1.0005, 0.9999, 1.0002];
            let prescribed = affine(std::array::from_fn(|i| {
                std::array::from_fn(|j| if i == j { stretch[i] - 1.0 } else { 0.0 })
            }));
            let next = solid
                .advance_electrical_candidate(
                    thermal,
                    circuit,
                    prepared,
                    &[],
                    1.0,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let fixed = original
                .advance_electrical_candidate(
                    thermal,
                    circuit,
                    prepared,
                    &[],
                    1.0,
                    &prescribed,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let material = prepared.resolve(cx, next.thermal()).unwrap();
            let heat = &next.report().circuit_thermal.electrical.resistors[0];
            let conductor = &heat.material.conductor;
            assert_eq!(
                conductor.material(),
                material.bulk.resistors[0].conductor.material()
            );
            close(conductor.length_m(), L * stretch[0], 1e-14);
            close(
                conductor.area_m2(),
                ELECTRICAL_AREA * stretch[1] * stretch[2],
                1e-16,
            );
            let t = next.thermal().temperature_k(CELLS[0]).unwrap();
            let r = 10.0 * (1.0 + 0.01 * (t - 300.0)) * L / ELECTRICAL_AREA * stretch[0]
                / (stretch[1] * stretch[2]);
            let q = capacitor_charge(CHARGE, r, 1.0);
            close(conductor.resistance_ohm(), r, 1e-9);
            close(next.circuit().charges_c()[0], q, 1e-12);
            close(heat.current_a, CHARGE - q, 1e-12);
            close(
                heat.heat_j,
                (CHARGE * CHARGE - q * q) / (2.0 * CAPACITANCE),
                1e-12,
            );
            assert!((next.circuit().charges_c()[0] - fixed.circuit().charges_c()[0]).abs() > 1e-10);
            close(next.report().physical_plastic_heat_j, 0.0, 1e-14);
            close(
                next.report().circuit_thermal.thermal.stored_energy_change_j,
                heat.heat_j,
                1e-9,
            );
            assert!(
                next.report().combined_residual_j.abs()
                    <= controls().coupled.total_energy_tolerance_j
            );
            for id in CELLS {
                assert_eq!(
                    next.thermal().cells()[&id].mass_kg,
                    thermal.cells()[&id].mass_kg
                );
                assert_eq!(
                    next.thermal().cells()[&id].geometry,
                    thermal.cells()[&id].geometry
                );
            }
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(circuit.charges_c(), &[CHARGE]);
        },
    );
}

#[test]
fn g1_native_inertial_stretch_changes_resistance_and_continues_from_reference_geometry() {
    with_powered(
        1e-5,
        false,
        false,
        CHARGE,
        |cx, _, thermal, original, circuit, prepared| {
            let solid = original
                .clone()
                .with_deforming_conductors(path([0.0, 0.0, 1.0]), geometry_controls(), cx)
                .unwrap();
            let dt = 2e-6;
            let mut velocity = vec![[0.0; 3]; NODES.len()];
            velocity[3][2] = 2.0;
            let moving = ThermoplasticDynamics::try_new(solid, velocity, cx).unwrap();
            let constraints = free_z(&[3]);
            let first = moving
                .advance_electrical_candidate(
                    thermal,
                    circuit,
                    prepared,
                    &[],
                    dt,
                    &constraints,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            // One consistent-mass diagonal and the exact constrained P1 axial stiffness.
            let mass = MASS / 10.0;
            let stiffness =
                (1.2e9 * (1.0 - 0.2) / ((1.0 + 0.2) * (1.0 - 2.0 * 0.2))) * VOLUME / (L * L);
            let u = mass * dt * 2.0 / (mass + stiffness * dt * dt);
            close(
                first.dynamics().solid().solution().displacement_m[3][2],
                u,
                1e-12,
            );
            let report = &first.report().material;
            let conductor = &report.circuit_thermal.electrical.resistors[0]
                .material
                .conductor;
            close(conductor.length_m(), L + u, 1e-14);
            close(conductor.area_m2(), ELECTRICAL_AREA, 1e-16);
            close(
                first.circuit().charges_c()[0],
                capacitor_charge(CHARGE, conductor.resistance_ohm(), dt),
                1e-12,
            );
            assert!(first.report().inertia.numerical_loss_j > 0.0);
            close(
                report.circuit_thermal.thermal.stored_energy_change_j,
                report.circuit_thermal.deposited_joule_j,
                1e-9,
            );
            assert!(
                report.combined_residual_j.abs() <= controls().coupled.total_energy_tolerance_j
            );
            let second = first
                .dynamics()
                .advance_electrical_candidate(
                    first.thermal(),
                    first.circuit(),
                    prepared,
                    &[],
                    2.0 * dt,
                    &constraints,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .unwrap();
            let u2 = 2.0 * mass * u / (mass + stiffness * dt * dt);
            close(
                second.dynamics().solid().solution().displacement_m[3][2],
                u2,
                1e-12,
            );
            let conductor = &second
                .report()
                .material
                .circuit_thermal
                .electrical
                .resistors[0]
                .material
                .conductor;
            close(conductor.length_m(), L + u2, 1e-14); // Not (L+u)*(1+u2/L).
            close(
                second.circuit().charges_c()[0],
                capacitor_charge(
                    first.circuit().charges_c()[0],
                    conductor.resistance_ohm(),
                    dt,
                ),
                1e-12,
            );
            assert_eq!(moving.time_s(), 0.0);
            assert_eq!(thermal.time_s(), 0.0);
            assert_eq!(circuit.charges_c(), &[CHARGE]);
        },
    );
}

fn expanding_fixture(alpha: f64) -> Fixture {
    let old = fixture_with_property(
        EMISSIVITY,
        Some(scalar(
            "electrical_resistivity",
            ELECTRICAL_RESISTIVITY_DIMS,
            1.0,
        )),
    );
    let cards = old.cards.map(|card| {
        let mut claims = card.claims().clone();
        claims
            .insert_claim(scalar(ALPHA, INVERSE_TEMPERATURE_DIMS, alpha))
            .unwrap();
        MaterialCard::assemble(card.id().clone(), claims, vec![]).unwrap()
    });
    let lower = old.point.clone().with("T", 280.0).unwrap();
    let upper = old.point.clone().with("T", 320.0).unwrap();
    let cp = PropertyKey::new(CP, CP_DIMS);
    let rho = PropertyKey::new(DENSITY_PROPERTY, Density::DIMS);
    let caloric = [0, 1].map(|i| {
        ResolvedHeatCapacityCurve::try_new(
            &cards[i],
            HeatCapacitySourceQuery {
                heat_capacity: &cp,
                density: &rho,
                temperature_axis: "T",
                lower: &lower,
                upper: &upper,
                selection: SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            -100.0,
            1e-4,
        )
        .unwrap()
    });
    Fixture {
        cards,
        caloric,
        point: old.point,
        conductivity: old.conductivity,
    }
}

#[test]
fn g3_expansion_material_data_changes_shape_current_and_heat_in_one_endpoint() {
    let mut charges = Vec::new();
    for alpha in [0.0, 2e-5] {
        with_heat_fixture(
            expanding_fixture(alpha),
            THERMAL_CONTACT,
            |cx, _, f, thermal, _, original| {
                let mut sources = f.sources();
                for source in &mut sources {
                    source.query.thermal_expansion =
                        Some(MaterialPropertySelection::SingleClaimOnly);
                }
                let solid = ThermoplasticSolid::try_new(
                    thermal,
                    &NODES,
                    &TETS,
                    sources,
                    TetAssemblyBudget::standard(),
                    cx,
                )
                .unwrap()
                .with_deforming_conductors(path([0.0, 0.0, 1.0]), geometry_controls(), cx)
                .unwrap();
                let resistors = [CellResistor {
                    destination: ResistorCell {
                        branch: 1,
                        cell: CELLS[0],
                    },
                    selection: MaterialPropertySelection::SingleClaimOnly,
                    length_m: L,
                    area_m2: ELECTRICAL_AREA,
                }];
                let bulk = MaterialNetwork {
                    cells: original.electrical.bulk.cells,
                    links: original.electrical.bulk.links,
                    resistors: &resistors,
                };
                let bases = [(1, ElectricalPropertyBasis::Resistivity)];
                let electrical = ElectricalContactNetwork {
                    bulk: &bulk,
                    contacts: &[],
                    electrical_bases: &bases,
                };
                let network = MaterialHeatNetwork {
                    electrical: &electrical,
                    thermal_contacts: original.thermal_contacts,
                    radiation: original.radiation,
                };
                let prepared = network.prepare(cx, thermal).unwrap();
                let initial = prepared.resolve(cx, thermal).unwrap();
                let circuit = MaterialCircuit::new_with_contacts(
                    CircuitGraph {
                        node_count: 2,
                        branches: vec![
                            (
                                1,
                                0,
                                Branch::Capacitor {
                                    farads: CAPACITANCE,
                                },
                            ),
                            (1, 0, Branch::Resistor { ohms: 0.0 }),
                        ],
                        transformers: vec![],
                    },
                    initial.bulk.resistors,
                    initial.contacts,
                    vec![],
                    vec![CHARGE],
                    &[],
                )
                .unwrap();
                let next = solid
                    .advance_electrical_candidate(
                        thermal,
                        &circuit,
                        &prepared,
                        &[],
                        1.0,
                        &free_z(&[3, 4]),
                        &[],
                        IncrementSettings::default(),
                        controls(),
                        &[0.0; 2],
                        cx,
                    )
                    .unwrap();
                let t = next.thermal().temperature_k(CELLS[0]).unwrap();
                // Laterally constrained thermal strain with sigma_zz=0:
                // epsilon_zz=(1+nu)/(1-nu)*alpha*dT. Independent of the FE solve.
                let extension = L * (1.0 + 0.2) / (1.0 - 0.2) * alpha * (t - 300.0);
                close(
                    next.solid().solution().displacement_m[3][2],
                    extension,
                    1e-12,
                );
                let conductor = &next.report().circuit_thermal.electrical.resistors[0]
                    .material
                    .conductor;
                close(conductor.length_m(), L + extension, 1e-14);
                close(
                    conductor.resistance_ohm(),
                    (L + extension) / ELECTRICAL_AREA,
                    1e-9,
                );
                close(
                    next.circuit().charges_c()[0],
                    capacitor_charge(CHARGE, conductor.resistance_ohm(), 1.0),
                    1e-12,
                );
                close(next.report().mechanical_work_j, 0.0, 1e-10);
                close(next.report().physical_plastic_heat_j, 0.0, 1e-12);
                assert!(
                    next.report().combined_residual_j.abs()
                        <= controls().coupled.total_energy_tolerance_j
                );
                assert_eq!(next.thermal().cells()[&CELLS[0]].mass_kg, MASS);
                charges.push(next.circuit().charges_c()[0]);
            },
        );
    }
    assert!(charges[1] > charges[0] + 1e-10);
}

#[test]
fn g0_nonprincipal_and_excessive_deformation_refuse_without_publishing_any_history() {
    with_powered(
        10.0,
        false,
        false,
        CHARGE,
        |cx, _, thermal, original, circuit, prepared| {
            let solid = original
                .clone()
                .with_deforming_conductors(path([1.0, 0.0, 0.0]), geometry_controls(), cx)
                .unwrap();
            let cells = thermal.cells().clone();
            let state = original.solution().clone();
            let charge = circuit.clone();
            let shear = affine([[0.0, 0.0, 0.0001], [0.0; 3], [0.0; 3]]);
            let error = solid
                .advance_electrical_candidate(
                    thermal,
                    circuit,
                    prepared,
                    &[],
                    1.0,
                    &shear,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .err()
                .unwrap();
            assert!(matches!(
                error,
                ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(
                    ElectricalSolidSourceError::Invalid(
                        "deformed uniform conductor requires positive volume and a principal current direction"
                    )
                ))
            ));
            let mut bound = geometry_controls();
            bound.maximum_displacement_gradient = 0.00005;
            let limited = original
                .clone()
                .with_deforming_conductors(path([1.0, 0.0, 0.0]), bound, cx)
                .unwrap();
            let stretch = affine([[0.0001, 0.0, 0.0], [0.0; 3], [0.0; 3]]);
            let error = limited
                .advance_electrical_candidate(
                    thermal,
                    circuit,
                    prepared,
                    &[],
                    1.0,
                    &stretch,
                    &[],
                    IncrementSettings::default(),
                    controls(),
                    &[0.0; 2],
                    cx,
                )
                .err()
                .unwrap();
            assert!(matches!(
                error,
                ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(
                    ElectricalSolidSourceError::Invalid(
                        "conductor displacement gradient exceeds its admitted affine geometry bound"
                    )
                ))
            ));
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(solid.solution(), &state);
            assert_eq!(circuit, &charge);
            assert!(
                original
                    .clone()
                    .with_deforming_conductors(
                        vec![path([1.0, 0.0, 0.0])[0]; 2],
                        geometry_controls(),
                        cx,
                    )
                    .is_err()
            );
            assert!(
                original
                    .clone()
                    .with_deforming_conductors(path([0.0; 3]), geometry_controls(), cx,)
                    .is_err()
            );
        },
    );
}
