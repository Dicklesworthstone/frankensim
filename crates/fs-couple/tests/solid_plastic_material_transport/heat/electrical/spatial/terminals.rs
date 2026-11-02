//! Real equipotential terminal control, sourced expansion and native motion.
//! Manufactured insulated prism: analytic uniform-field/free-expansion and
//! restrained-thermoelastic solutions, not fitted material or geometric data.
use super::*;
use fs_couple::circuit_enthalpy::conduction::spatial::SpatialConductionDrive;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
};
use std::collections::BTreeSet;

pub(super) const PRISM_NODES: [[f64; 3]; 6] = [
    [0.0, 0.0, 0.0],
    [0.0, L, 0.0],
    [0.0, 0.0, L],
    [L, 0.0, 0.0],
    [L, L, 0.0],
    [L, 0.0, L],
];
pub(super) const PRISM_TETS: [[usize; 4]; 3] = [[0, 1, 2, 5], [0, 1, 5, 4], [0, 4, 5, 3]];
const PRISM_CELLS: [CellId; 3] = [CellId(17), CellId(3), CellId(47)];
const AREA: f64 = L * L / 2.0;
const ALPHA: f64 = 2e-5;
const TOTAL_MASS: f64 = 3.0 * MASS;
const CAPACITY: f64 = TOTAL_MASS * 5.0;
const NINE_BULK: f64 = 6e9; // 9K for E=1.2 GPa, nu=0.2.

fn prism_controls() -> SpatialElectricalSolidSolve {
    let mut result = spatial_controls();
    result.coupled.coupled.maximum_cells = 3;
    result
}

fn tet_surface(nodes: &[[f64; 3]], tet: [usize; 4]) -> f64 {
    [[0, 1, 2], [0, 1, 3], [0, 2, 3], [1, 2, 3]]
        .into_iter()
        .map(|face| {
            let [a, b, c] = face.map(|i| nodes[tet[i]]);
            let u: [f64; 3] = core::array::from_fn(|i| b[i] - a[i]);
            let v: [f64; 3] = core::array::from_fn(|i| c[i] - a[i]);
            let cross = [
                u[1] * v[2] - u[2] * v[1],
                u[2] * v[0] - u[0] * v[2],
                u[0] * v[1] - u[1] * v[0],
            ];
            0.5 * cross.iter().map(|x| x * x).sum::<f64>().sqrt()
        })
        .sum()
}

fn prism_thermal(f: &Fixture) -> DepositionSimulation<'_> {
    prism_thermal_cells(
        f.caloric[0].curve(),
        MATS[0],
        &PRISM_NODES,
        &PRISM_TETS,
        &PRISM_CELLS,
    )
}

// Both terminal and matching-interface prisms retain one equal reference-mass
// caloric cell per original tet. Only their assembly and source card differ.
pub(super) fn prism_thermal_cells<'a>(
    curve: &'a dyn ThermalLaw,
    material_id: MaterialId,
    nodes: &[[f64; 3]],
    tets: &[[usize; 4]],
    cells: &[CellId],
) -> DepositionSimulation<'a> {
    assert_eq!(tets.len(), cells.len());
    let events = tets
        .iter()
        .enumerate()
        .map(|(i, &tet)| {
            let center: [f64; 3] =
                core::array::from_fn(|axis| tet.iter().map(|&v| nodes[v][axis] / 4.0).sum());
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(DepositedCell {
                    id: cells[i],
                    segment_id: i as u64,
                    body_id: 1,
                    material_id,
                    mass_kg: MASS,
                    specific_h_j_kg: 0.0,
                    geometry: CellGeometry {
                        start_m: [center[0] - L / 16.0, center[1], center[2]],
                        end_m: [center[0] + L / 16.0, center[1], center[2]],
                        volume_m3: VOLUME,
                        surface_area_m2: tet_surface(nodes, tet),
                    },
                    boundary: boundary(),
                }),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(material_id, curve)],
        DepositionBudget {
            max_cells: cells.len(),
            max_links: 0,
            max_events: cells.len(),
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(0.0, thermal_config(), || Ok(()))
    .unwrap()
    .0
}

fn with_prism<R>(
    rho: f64,
    run: impl FnOnce(
        &Cx<'_>,
        &DepositionSimulation<'_>,
        &ThermoplasticSolid<'_>,
        &SpatialConduction,
        &PreparedHeatNetwork<'_, '_>,
        &[Vec<usize>; 2],
    ) -> R,
) -> R {
    let mut resistivity = scalar("electrical_resistivity", ELECTRICAL_RESISTIVITY_DIMS, rho);
    resistivity.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![(280.0, 0.8 * rho), (320.0, 1.2 * rho)],
        dims: ELECTRICAL_RESISTIVITY_DIMS,
    };
    resistivity.interpolation = InterpolationPolicy::LinearInside;
    resistivity.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    let f = fixture_with_properties(
        EMISSIVITY,
        &[
            resistivity,
            scalar(
                LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
                INVERSE_TEMPERATURE_DIMS,
                ALPHA,
            ),
        ],
    );
    with_cx(|_, cx| {
        let thermal = prism_thermal(&f);
        let sources = PRISM_CELLS
            .iter()
            .map(|&cell| ThermoplasticCell {
                cell,
                query: SourceQuery {
                    card: &f.cards[0],
                    caloric: &f.caloric[0],
                    context: f.point.clone(),
                    temperature_axis: "T".into(),
                    selection: MaterialPropertySelection::SingleClaimOnly,
                    hardening_basis: BASIS,
                    thermal_expansion: Some(MaterialPropertySelection::SingleClaimOnly),
                    strain_limit: 0.04,
                },
            })
            .collect();
        let solid = ThermoplasticSolid::try_new(
            &thermal,
            &PRISM_NODES,
            &PRISM_TETS,
            sources,
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        // Uniform exact endpoints need no exchange between these insulated
        // finite cells. Their actual reference mass and caloric owners persist.
        let selections = vec![None; 3];
        let cells = solid.thermal_material_cells(&selections).unwrap();
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &[],
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let network = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: &[],
            radiation: None,
        };
        let prepared = network.prepare(cx, &thermal).unwrap();
        let conductor = resolve_conductor_state_with_basis(
            &f.cards[0],
            &f.point,
            MaterialPropertySelection::SingleClaimOnly,
            ElectricalPropertyBasis::Resistivity,
        )
        .unwrap();
        let problem = SpatialConduction::new(
            PRISM_NODES.to_vec(),
            PRISM_TETS.map(|tet| tet.map(|v| v as u32)).to_vec(),
            vec![conductor; 3],
            prism_controls().assembly,
            &mut || false,
        )
        .unwrap();
        let terminals: [Vec<usize>; 2] = [L, 0.0].map(|x| {
            problem
                .exterior_faces()
                .filter(|&face| {
                    problem.mesh().faces[face]
                        .iter()
                        .all(|&v| PRISM_NODES[v as usize][0] == x)
                })
                .collect()
        });
        assert_eq!(terminals.each_ref().map(Vec::len), [1, 1]);
        run(cx, &thermal, &solid, &problem, &prepared, &terminals)
    })
}

pub(super) fn current_drive(
    terminals: &[Vec<usize>; 2],
    current_a: f64,
) -> SpatialConductionDrive<'_> {
    SpatialConductionDrive::TotalCurrent {
        positive_faces: &terminals[0],
        negative_faces: &terminals[1],
        current_a,
        reference_voltage_v: 0.7,
    }
}

fn check_current(
    problem: &SpatialConduction,
    field: &SpatialConductionField,
    terminals: &[Vec<usize>; 2],
    rho: f64,
    temperature: f64,
    stretch: f64,
    current: f64,
) {
    let resistance = rho * (1.0 + RHO_SLOPE * (temperature - 300.0)) * L / (AREA * stretch);
    let voltage = current * resistance;
    for (side, faces) in terminals.iter().enumerate() {
        let nodes: BTreeSet<_> = faces
            .iter()
            .flat_map(|&face| problem.mesh().faces[face])
            .map(|v| v as usize)
            .collect();
        let mut supplied = 0.0;
        for node in nodes {
            close(
                field.voltage_v[node],
                if side == 0 { voltage } else { 0.0 },
                2e-8,
            );
            supplied += field.nodal_boundary_current_a[node].unwrap();
        }
        close(supplied, if side == 0 { current } else { -current }, 2e-10);
    }
    let density = -current / (AREA * stretch * stretch);
    for i in 0..PRISM_TETS.len() {
        close(field.current_density_a_m2[i][0], density, 2e-6);
        close(field.current_density_a_m2[i][1], 0.0, 2e-6);
        close(field.current_density_a_m2[i][2], 0.0, 2e-6);
        close(
            field.electric_field_v_m[i][0],
            -voltage / (L * stretch),
            3e-6,
        );
        close(field.cell_power_w[i], current * voltage / 3.0, 2e-9);
    }
    close(field.power_w, current * voltage, 5e-9);
    close(field.supplied_power_w, current * voltage, 5e-9);
}

fn free_temperature(old: f64, rho: f64, current: f64, dt: f64) -> f64 {
    // C*dT = dt*I^2*rho(T)*L/(A*lambda(T)); reference mass stays fixed.
    let c = dt * current * current * rho * L / (AREA * CAPACITY);
    let b = 1.0 + ALPHA * (old - 300.0) - c * RHO_SLOPE;
    let right = c * (1.0 + RHO_SLOPE * (old - 300.0));
    old + 2.0 * right / (b + (b * b + 4.0 * ALPHA * right).sqrt())
}

fn mean_temperature(thermal: &DepositionSimulation<'_>) -> f64 {
    PRISM_CELLS
        .iter()
        .map(|&cell| thermal.temperature_k(cell).unwrap() / 3.0)
        .sum()
}

fn check_free_expansion(
    candidate: &SpatialElectricalSolidCandidate<'_, '_>,
    terminals: &[Vec<usize>; 2],
    rho: f64,
    current: f64,
    old_temperature: f64,
    dt: f64,
) -> f64 {
    let temperature = free_temperature(old_temperature, rho, current, dt);
    let strain = ALPHA * (temperature - 300.0);
    for cell in PRISM_CELLS {
        close(
            candidate.thermal().temperature_k(cell).unwrap(),
            temperature,
            3e-6,
        );
        close(candidate.thermal().cells()[&cell].mass_kg, MASS, 0.0);
    }
    for (node, reference) in PRISM_NODES.iter().enumerate() {
        for axis in 0..3 {
            close(
                candidate.solid().solution().displacement_m[node][axis],
                strain * reference[axis],
                2e-10,
            );
            close(
                candidate.problem().positions_m()[node][axis],
                (1.0 + strain) * reference[axis],
                2e-10,
            );
        }
    }
    for stress in &candidate.solid().solution().element_stress {
        for &component in stress {
            close(component, 0.0, 0.03);
        }
    }
    let report = candidate.report();
    let electrical = &report.circuit_thermal.electrical;
    // Test the field at the accepted temperature as well as the independent
    // thermal root above; power tolerances need not inherit temperature error.
    let field_temperature = mean_temperature(candidate.thermal());
    check_current(
        candidate.problem(),
        &electrical.field,
        terminals,
        rho,
        field_temperature,
        1.0 + ALPHA * (field_temperature - 300.0),
        current,
    );
    let heat = CAPACITY * (temperature - old_temperature);
    close(electrical.supplied_work_j, heat, 2e-8);
    close(electrical.heat_j, heat, 2e-8);
    close(
        report.circuit_thermal.thermal.stored_energy_change_j,
        heat,
        2e-8,
    );
    close(report.deposited_mechanical_heat_j, 0.0, 2e-8);
    close(report.combined_residual_j, 0.0, 1e-8);
    temperature
}

#[test]
fn current_control_resistivity_swap_changes_solved_expansion_and_retains_terminal_faces() {
    let mut endpoints = Vec::new();
    // Exactly six rigid-body constraints; expansion is solved, not prescribed.
    let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (3, 0.0), (5, 0.0), (6, 0.0)];
    for rho in [0.1, 0.2] {
        with_prism(rho, |cx, thermal, solid, problem, prepared, terminals| {
            let current = 0.025;
            let first = solid
                .advance_spatial_electrical_candidate(
                    thermal,
                    problem,
                    prepared,
                    current_drive(terminals, current),
                    1.0,
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    prism_controls(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            let first_temperature =
                check_free_expansion(&first, terminals, rho, current, 300.0, 1.0);
            endpoints.push((
                first_temperature,
                first.solid().solution().displacement_m[3][0],
            ));
            let (thermal, solid, problem, _) = first.into_parts();
            // Reversing the held current reverses terminal voltage/flux, while
            // its I^2 heat continues on the SAME expanded, source-bound mesh.
            let next = solid
                .advance_spatial_electrical_candidate(
                    &thermal,
                    &problem,
                    prepared,
                    current_drive(terminals, -current),
                    1.25,
                    &supports,
                    &[],
                    IncrementSettings::default(),
                    prism_controls(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            check_free_expansion(&next, terminals, rho, -current, first_temperature, 0.25);
            assert_eq!(thermal.time_s(), 1.0);
            assert_eq!(solid.time_s(), 1.0);
        });
    }
    // At fixed total current higher resistivity increases heat and actual shape
    // change. No geometry, mechanical material, force or thermal source changed.
    assert!(endpoints[1].0 > endpoints[0].0 + 5.0);
    assert!(endpoints[1].1 > 2.0 * endpoints[0].1);
}

fn translation(distance: f64) -> Vec<(usize, f64)> {
    (0..3 * PRISM_NODES.len())
        .map(|dof| (dof, if dof % 3 == 0 { distance } else { 0.0 }))
        .collect()
}

fn restrained_temperature(old: f64, rho: f64, current: f64, dt: f64) -> f64 {
    // eps=0: sigma=-3K*alpha*(T-300)I. Original reversible transfer is
    // -Told*alpha*delta(tr(sigma)), so use this independent effective capacity.
    let capacity = CAPACITY - 3.0 * VOLUME * NINE_BULK * old * ALPHA * ALPHA;
    let q = dt * current * current * rho * L / AREA;
    old + q * (1.0 + RHO_SLOPE * (old - 300.0)) / (capacity - q * RHO_SLOPE)
}

#[test]
fn current_driven_inertial_constraint_heat_and_work_continue_and_refusal_is_atomic() {
    with_prism(0.1, |cx, thermal, solid, problem, prepared, terminals| {
        let rho = 0.1;
        let current = 2.0;
        let dt = 2e-4;
        let initial_speed = 0.25;
        let acceleration = 1000.0;
        let speed = initial_speed + dt * acceleration;
        let displacement = dt * speed;
        let moving = ThermoplasticDynamics::try_new(
            solid.clone(),
            vec![[initial_speed, 0.0, 0.0]; PRISM_NODES.len()],
            cx,
        )
        .unwrap();
        let first = moving
            .advance_spatial_electrical_candidate(
                thermal,
                problem,
                prepared,
                current_drive(terminals, current),
                dt,
                &translation(displacement),
                &[],
                IncrementSettings::default(),
                prism_controls(),
                &[0.0; 3],
                cx,
            )
            .unwrap();
        let temperature = restrained_temperature(300.0, rho, current, dt);
        for cell in PRISM_CELLS {
            close(
                first.thermal().temperature_k(cell).unwrap(),
                temperature,
                3e-6,
            );
        }
        for node in 0..PRISM_NODES.len() {
            close(first.dynamics().velocity_mps()[node][0], speed, 1e-11);
            close(
                first.problem().positions_m()[node][0],
                PRISM_NODES[node][0] + displacement,
                1e-14,
            );
        }
        for stress in &first.dynamics().solid().solution().element_stress {
            for &component in &stress[..3] {
                close(
                    component,
                    -NINE_BULK / 3.0 * ALPHA * (temperature - 300.0),
                    0.3,
                );
            }
        }
        let report = first.report();
        let material = &report.material;
        let electrical = &material.circuit_thermal.electrical;
        check_current(
            first.problem(),
            &electrical.field,
            terminals,
            rho,
            mean_temperature(first.thermal()),
            1.0,
            current,
        );
        close(
            report.inertia.kinetic_change_j,
            0.5 * TOTAL_MASS * (speed * speed - initial_speed * initial_speed),
            1e-12,
        );
        close(
            report.inertia.numerical_loss_j,
            0.5 * TOTAL_MASS * (speed - initial_speed).powi(2),
            1e-12,
        );
        close(
            material.mechanical_work_j,
            TOTAL_MASS * acceleration * displacement,
            1e-11,
        );
        let reversible = 3.0 * VOLUME * NINE_BULK * 300.0 * ALPHA * ALPHA * (temperature - 300.0);
        close(material.physical_plastic_heat_j, 0.0, 1e-12);
        close(material.deposited_mechanical_heat_j, reversible, 1e-8);
        close(
            material.circuit_thermal.thermal.stored_energy_change_j,
            electrical.heat_j + reversible,
            1e-8,
        );
        close(
            electrical.supplied_work_j,
            dt * current * current * rho * (1.0 + RHO_SLOPE * (temperature - 300.0)) * L / AREA,
            1e-8,
        );
        close(material.combined_residual_j, 0.0, 1e-8);
        assert!(report.inertia.numerical_loss_j > 1e-6);

        let (thermal, moving, problem, _) = first.into_parts();
        let dt2 = dt / 4.0;
        let speed2 = speed + dt2 * acceleration;
        let displacement2 = displacement + dt2 * speed2;
        let prescribed = translation(displacement2);
        let next = moving
            .advance_spatial_electrical_candidate(
                &thermal,
                &problem,
                prepared,
                current_drive(terminals, -current),
                dt + dt2,
                &prescribed,
                &[],
                IncrementSettings::default(),
                prism_controls(),
                &[0.0; 3],
                cx,
            )
            .unwrap();
        let temperature2 = restrained_temperature(temperature, rho, current, dt2);
        for cell in PRISM_CELLS {
            close(
                next.thermal().temperature_k(cell).unwrap(),
                temperature2,
                4e-6,
            );
        }
        for velocity in next.dynamics().velocity_mps() {
            close(velocity[0], speed2, 1e-11);
        }
        check_current(
            next.problem(),
            &next.report().material.circuit_thermal.electrical.field,
            terminals,
            rho,
            mean_temperature(next.thermal()),
            1.0,
            -current,
        );
        close(
            next.report().inertia.kinetic_change_j,
            0.5 * TOTAL_MASS * (speed2 * speed2 - speed * speed),
            1e-12,
        );
        close(next.report().material.combined_residual_j, 0.0, 1e-8);

        let saved_cells = thermal.cells().clone();
        let saved_solid = moving.solid().solution().clone();
        let saved_velocity = moving.velocity_mps().to_vec();
        let saved_geometry = problem.positions_m().to_vec();
        let saved_material = problem.conductors().to_vec();
        let mut bounded = prism_controls();
        assert!(next.report().material.mechanical_trials > 1);
        bounded.coupled.maximum_mechanical_trials = next.report().material.mechanical_trials - 1;
        let refused = moving.advance_spatial_electrical_candidate(
            &thermal,
            &problem,
            prepared,
            current_drive(terminals, -current),
            dt + dt2,
            &prescribed,
            &[],
            IncrementSettings::default(),
            bounded,
            &[0.0; 3],
            cx,
        );
        assert!(matches!(
            refused,
            Err(ElectricalSolidError::Coupled(CircuitEnthalpyError::Source(
                ElectricalSolidSourceError::MechanicalBudget
            )))
        ));
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(moving.solid().solution(), &saved_solid);
        assert_eq!(moving.velocity_mps(), &saved_velocity);
        assert_eq!(problem.positions_m(), &saved_geometry);
        assert_eq!(problem.conductors(), &saved_material);
        assert_eq!(moving.time_s(), dt);
    });
}
