//! G0/G3/G4: one finite material inventory receives the original rotor losses.
#[path = "../../../fs-em/tests/winding_ohmic/common.rs"]
mod material_source;

use super::*;
use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve,
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        material::{MaterialCell, MaterialNetwork},
        spatial_induction::{
            SpatialMagneticRegion, SpatialMagneticThermalLaw, SpatialMagneticThermalNetwork,
            ThermalMagneticMaterial, rotor::heat_assembly_with_boundary_candidate,
        },
    },
    material_circuit::MaterialCircuitSolve,
    scenario_sensor::{TemperatureSensorObserver, ThermalSensorSupport},
};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::SelectionPolicy;
use fs_material::{
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::SPECIFIC_HEAT_CAPACITY_DIMS,
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::Density;
use fs_scenario::{
    EntityDeclaration, EntityRef, KindExpectation, ObservationSupport, PlacementUncertainty,
    ScenarioSensor, SensorCalibration, SensorDynamics, SensorKind, SensorLocation, SensorMount,
};
use fs_solver::NewtonKrylovConfig;
use fs_thermal::enthalpy::*;

fn assembly() -> WeightedAssemblyLimits {
    WeightedAssemblyLimits {
        max_cells: 512,
        max_dofs: 4096,
        max_triplets: 100000,
    }
}

fn thermal_controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 1,
        maximum_circuit_trials: 4096,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-3,
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 2e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-8,
            power_tolerance_w: 1e-8,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}

#[test]
fn finite_temperature_rotor_retains_field_speed_and_source_bound_joule_heat() {
    finite_temperature_rotor(false);
}

#[test]
fn applied_magnetic_history_reaches_rotor_temperature_and_complete_heat_assembly() {
    finite_temperature_rotor(true);
}

fn finite_temperature_rotor(driven: bool) {
    const SECTORS: usize = 16;
    const ALPHA: f64 = 0.005;
    const SPEED: f64 = 6.0;
    let dt = 0.1;
    let basis = ElectricalPropertyBasis::Resistivity;
    // Positive rho(T) over [200,2000] K, constant mu, sourced density=1 kg/m3
    // and cp=4 J/(kg K). These are manufactured values, not a measured rotor.
    let card = material_source::absolute_affine_magnetic_card(1.0 / SIGMA, ALPHA, MU, 0.0, basis);
    let context = material_source::absolute_point(300.0);
    let bundle = material_source::bundle_at(&card, &context, basis);
    let conductor =
        ResolvedConductorState::try_from_material_with_basis(bundle.clone(), basis).unwrap();
    let (geometry, _) = problem(SECTORS);
    let problem = SpatialVectorMagnetostatics::new(
        geometry.positions_m().to_vec(),
        geometry.mesh().tets.clone(),
        vec![SpatialMagneticMaterial::Linear(bundle); geometry.mesh().tets.len()],
        assembly(),
        &mut || false,
    )
    .unwrap();
    let volume: f64 = problem
        .mesh()
        .tets
        .iter()
        .map(|tet| six_volume(tet.map(|v| problem.positions_m()[v as usize])) / 6.0)
        .sum();
    // Exact regular-polygon prism moments: Iz/M = R^2 (2+cos(2pi/N))/6.
    // The native mechanical and finite caloric masses use the same density
    // and original tet volume. No thin-disk or guessed motor inertia is used.
    let inertia = volume * (2.0 + (std::f64::consts::TAU / SECTORS as f64).cos()) / 6.0;
    let transverse = 0.5 * inertia + volume * (2.0 * HALF_LENGTH).powi(2) / 12.0;
    let model = ArticulatedModel::new(vec![Link::new(
        "same finite thermal rotor inventory",
        None,
        Se3::identity(),
        JointModel::revolute(Vec3::new(0.0, 0.0, 1.0), None).unwrap(),
        SpatialInertia::new(
            volume,
            Vec3::default(),
            Mat3 {
                m: [
                    transverse, 0.0, 0.0, 0.0, transverse, 0.0, 0.0, 0.0, inertia,
                ],
            },
        )
        .unwrap(),
    )])
    .unwrap();
    let mut rotor = SpatialInductionRotor::new(
        initial(&problem, vec![Some(conductor); problem.mesh().tets.len()]),
        model,
        shape(),
        0.0,
        SPEED,
        limits(),
        &mut || false,
    )
    .unwrap();
    let curve = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &material_source::absolute_point(200.0),
            upper: &material_source::absolute_point(2000.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-6,
    )
    .unwrap();
    let mut heat = DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(7),
                segment_id: 0,
                body_id: 0,
                material_id: MaterialId(9),
                mass_kg: volume,
                specific_h_j_kg: 400.0, // cp * (300 K - 200 K reference)
                geometry: CellGeometry {
                    start_m: [0.0, 0.0, -HALF_LENGTH],
                    end_m: [0.0, 0.0, HALF_LENGTH],
                    volume_m3: volume,
                    surface_area_m2: volume / HALF_LENGTH
                        + 4.0
                            * SECTORS as f64
                            * (std::f64::consts::PI / SECTORS as f64).sin()
                            * HALF_LENGTH,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MaterialId(9), curve.curve() as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: dt,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let regions: Vec<_> = (0..problem.mesh().tets.len())
        .map(|tetrahedron| SpatialMagneticRegion::Cell {
            tetrahedron,
            cell: CellId(7),
        })
        .collect();
    let laws: Vec<_> = problem
        .materials()
        .iter()
        .enumerate()
        .map(|(tetrahedron, source)| {
            let SpatialMagneticMaterial::Linear(source) = source else {
                unreachable!()
            };
            SpatialMagneticThermalLaw {
                tetrahedron,
                law: ThermalMagneticMaterial::new(&card, source, "T", [200.0, 2000.0]).unwrap(),
                maximum_flux_density_t: 8.0,
            }
        })
        .collect();
    let network = SpatialMagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &regions,
        volume_tolerance_m3: 1e-12,
    };
    let zero = vec![[0.0; 3]; problem.mesh().tets.len()];
    let boundary = uniform_trace(&problem, [1.0, 0.1, 0.0]);
    let original_a = rotor.field().potential_wb().to_vec();
    let mut bounded = thermal_controls();
    bounded.maximum_circuit_trials = 1;
    let refused = network.step_rotor_candidate(
        &heat,
        dt,
        &rotor,
        &zero,
        0.15,
        &laws,
        bounded,
        field_controls(),
        controls(),
        assembly(),
        &[0.0],
        || false,
    );
    assert!(matches!(refused, Err(CircuitEnthalpyError::CircuitBudget)));
    assert_eq!(rotor.field().potential_wb(), original_a);
    assert_eq!(rotor.angular_velocity_rad_s(), SPEED);
    assert_eq!(rotor.field().time_s(), 0.0);
    assert_eq!(heat.time_s(), 0.0);
    near(heat.temperature_k(CellId(7)).unwrap(), 300.0, 1e-12);

    let assembly_id = EntityDeclaration::assembly("source-driven-rotor").identity();
    let part = EntityDeclaration::part(assembly_id, "finite-rotor").identity();
    let entity = EntityDeclaration::region(part, "caloric-cell").identity();
    let instrument = ScenarioSensor::new(
        "rotor-thermocouple",
        SensorKind::Thermocouple,
        SensorLocation::new(
            EntityRef::new(entity, KindExpectation::Domain),
            [0.0; 3],
            PlacementUncertainty::declared_exact("original whole-cell support").unwrap(),
        )
        .unwrap(),
        ObservationSupport::point(1, 0).unwrap(),
        SensorMount::affine(1.0, 0.0, "fixture mount").unwrap(),
        SensorDynamics::first_order(0.25, "fixture instrument response").unwrap(),
        SensorCalibration::physical("rotor-fixture", "2026-07-01", "synthetic", 0.04).unwrap(),
        false,
    )
    .unwrap()
    .compile()
    .unwrap();
    let mut thermometer = TemperatureSensorObserver::new(
        &instrument,
        ThermalSensorSupport {
            entity,
            cells: &[CellId(7)],
        },
        &heat,
        300.0,
        || Ok(()),
    )
    .unwrap();

    for (index, shaft) in [0.15, 0.0].into_iter().enumerate() {
        let apply_boundary = driven && index == 0;
        let target = (index + 1) as f64 * dt;
        let temperature_before = heat.temperature_k(CellId(7)).unwrap();
        let speed_before = rotor.angular_velocity_rad_s();
        let angle_before = rotor.angle_rad();
        let b_before = rotor.field().flux_density_t(&mut || false).unwrap();
        let cold = if apply_boundary {
            rotor.step_with_boundary_candidate(
                target,
                &zero,
                &boundary,
                shaft,
                field_controls(),
                controls(),
                &mut || false,
            )
        } else {
            rotor.step_candidate(
                target,
                &zero,
                shaft,
                field_controls(),
                controls(),
                &mut || false,
            )
        }
        .unwrap();
        let candidate = if apply_boundary {
            network.step_rotor_with_boundary_candidate(
                &heat,
                target,
                &rotor,
                &zero,
                &boundary,
                shaft,
                &laws,
                thermal_controls(),
                field_controls(),
                controls(),
                assembly(),
                &[0.0],
                || false,
            )
        } else {
            network.step_rotor_candidate(
                &heat,
                target,
                &rotor,
                &zero,
                shaft,
                &laws,
                thermal_controls(),
                field_controls(),
                controls(),
                assembly(),
                &[0.0],
                || false,
            )
        }
        .unwrap();
        let temperature = candidate.thermal().temperature_k(CellId(7)).unwrap();
        // Observe the converged implicit thermal state; do not integrate heat
        // again or advance the instrument through rejected Newton/speed trials.
        let original_instrument = thermometer.clone();
        let sensor_candidate = thermometer
            .observe_endpoint_candidate(&heat, candidate.thermal(), dt, || Ok(()))
            .unwrap();
        let weight = -(-dt / 0.25_f64).exp_m1();
        let old_reading = thermometer.state().reading();
        let expected = old_reading
            + (temperature_before - old_reading) * weight
            + (temperature - temperature_before) * (1.0 - 0.25 * weight / dt);
        assert!(sensor_candidate.state().reading() > old_reading);
        near(sensor_candidate.state().reading(), expected, 2e-10);
        near(
            sensor_candidate.state().equilibrium_reading(),
            temperature,
            1e-12,
        );
        assert_eq!(
            sensor_candidate.state().time_s(),
            candidate.thermal().time_s()
        );
        assert_eq!(thermometer, original_instrument);
        let speed = candidate.rotor().angular_velocity_rad_s();
        let midpoint = 0.5 * (speed_before + speed);
        let report = candidate.report();
        let native = &report.coupled.electrical;
        // The independent parent quadrature integrates the original P1 J.
        // Replace its constant 1/SIGMA by the authored rho at the accepted T.
        let (_, torque, reference_heat) = integrate(&problem, &native.induction, dt);
        let rho = (1.0 + ALPHA * (temperature - 300.0)) / SIGMA;
        let joule = reference_heat * SIGMA * rho;
        let b_after = candidate
            .rotor()
            .field()
            .flux_density_t(&mut || false)
            .unwrap();
        let magnetic_change = stored(&problem, &b_after) - stored(&problem, &b_before);
        let kinetic_change = 0.5 * inertia * (speed - speed_before) * (speed + speed_before);
        let caloric_change = 4.0 * volume * (temperature - temperature_before);
        let shaft_work = dt * shaft * midpoint;
        assert!(temperature > temperature_before && joule > 0.0);
        assert!(speed > 0.0 && native.balance.electromagnetic_torque_nm < 0.0);
        assert!((speed - cold.state().angular_velocity_rad_s()).abs() > 1e-10);
        near(native.induction.heat_j, joule, 2e-9);
        near(caloric_change, joule, 2e-8);
        near(report.coupled.deposited_joule_j, joule, 2e-9);
        near(native.mechanics.kinetic_energy_change, kinetic_change, 2e-9);
        near(
            inertia * (speed - speed_before),
            dt * (torque[2] + shaft),
            2e-9,
        );
        near(
            native.induction.mechanical_work_j,
            dt * midpoint * torque[2],
            5e-9,
        );
        near(
            magnetic_change + kinetic_change + caloric_change,
            shaft_work + native.induction.magnetic_boundary_work_j,
            3e-8,
        );
        near(report.total_energy_residual_j, 0.0, 2e-8);
        near(report.external_heat_j, 0.0, 2e-9);
        near(native.induction.parameter_work_j, 0.0, 2e-9);
        near(report.coupled.source.reversible_heat_j, 0.0, 2e-9);
        near(native.balance.shaft_work_j, shaft_work, 2e-9);
        near(
            candidate.rotor().angle_rad(),
            angle_before + dt * midpoint,
            1e-10,
        );
        assert_eq!(candidate.rotor().field().time_s(), target);
        assert_eq!(candidate.thermal().time_s(), target);
        for source in candidate.rotor().field().conductors().iter().flatten() {
            near(source.resistivity_ohm_m(), rho, 2e-12);
        }
        for sample in &report.coupled.source.fields {
            near(sample.temperature_k[0], temperature_before, 1e-10);
            near(sample.temperature_k[1], temperature, 1e-10);
        }
        assert_eq!(
            report.coupled.source.fields.len(),
            problem.mesh().tets.len()
        );
        assert!(report.coupled.circuit_trials > native.balance.iterations);
        if apply_boundary {
            assert!(native.induction.magnetic_boundary_work_j.abs() > 1e-8);
            // The complete provider must use the identical original field,
            // momentum, source-temperature and caloric residual, retaining
            // external magnetic supply even with no optional thermal links.
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
            let heat_network = MaterialHeatNetwork {
                electrical: &electrical,
                thermal_contacts: &[],
                radiation: None,
            };
            let gate = CancelGate::new();
            let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
            pool.scope(|arena| {
                let cx = Cx::new(
                    &gate,
                    arena,
                    StreamKey {
                        seed: 79,
                        kernel_id: 3,
                        tile: 0,
                        iteration: 0,
                    },
                    Budget::INFINITE,
                    ExecMode::Deterministic,
                );
                let prepared = heat_network.prepare(&cx, &heat).unwrap();
                let complete = heat_assembly_with_boundary_candidate(
                    &prepared,
                    &cx,
                    &heat,
                    target,
                    &rotor,
                    &zero,
                    &boundary,
                    shaft,
                    &regions,
                    1e-12,
                    &laws,
                    thermal_controls(),
                    field_controls(),
                    controls(),
                    assembly(),
                    &[0.0],
                )
                .unwrap();
                near(
                    complete.thermal().temperature_k(CellId(7)).unwrap(),
                    temperature,
                    2e-9,
                );
                near(complete.rotor().angular_velocity_rad_s(), speed, 2e-10);
                for (&a, &b) in complete
                    .rotor()
                    .field()
                    .potential_wb()
                    .iter()
                    .zip(candidate.rotor().field().potential_wb())
                {
                    near(a, b, 2e-10);
                }
                near(
                    complete
                        .report()
                        .coupled
                        .electrical
                        .induction
                        .magnetic_boundary_work_j,
                    native.induction.magnetic_boundary_work_j,
                    2e-9,
                );
                near(complete.report().total_energy_residual_j, 0.0, 2e-8);
                near(complete.report().external_heat_j, 0.0, 2e-9);
            });
        } else {
            assert_eq!(native.induction.magnetic_boundary_work_j, 0.0);
        }
        if driven {
            // The next ordinary interval holds the newly published trace;
            // it must not restore the constructor boundary or repeat the ramp.
            for trace in &boundary.potential {
                for (edge, expected) in problem
                    .face_edges(trace.face)
                    .unwrap()
                    .into_iter()
                    .zip(trace.edge_integrals_wb)
                {
                    assert_eq!(candidate.rotor().field().potential_wb()[edge], expected);
                }
            }
        }
        assert_eq!(rotor.field().time_s(), index as f64 * dt);
        near(
            heat.temperature_k(CellId(7)).unwrap(),
            temperature_before,
            1e-12,
        );
        (heat, rotor, _) = candidate.into_parts();
        thermometer = sensor_candidate;
    }
}
