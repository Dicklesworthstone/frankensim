//! G1/G3/G4: evolving buoyancy, conservative native enthalpy and atomic refusal.
use super::*;
use fs_airflow::graph::thermal::medium::convection::deposition::source::heat::relaxation::pressure::gravity::*;
use fs_airflow::graph::thermal::medium::series::hydraulic::pressure::{
    BoussinesqDensity, SeriesGravity,
};

fn closed(gz: f64) -> SeriesGravity {
    SeriesGravity::try_new([0.0, 0.0, gz], vec![[0.0; 3], [0.0, 0.0, 2.0], [0.0; 3]]).unwrap()
}

fn gravity_path<'a>(
    sections: &'a [PressureCoolantSection],
    links: &'a [PressureConvectiveLink],
    gravity: &'a SeriesGravity,
) -> GravityCoolantLoop<'a> {
    GravityCoolantLoop {
        sections,
        links,
        gravity,
        density: BoussinesqDensity::try_new(1000.0, 0.05).unwrap(),
        driving_pressure_pa: 0.0,
        pressure_budget: path(sections, links, 0.0).pressure_budget,
        maximum_links: 2,
        maximum_omitted_conversion_w: 1.0,
        maximum_omitted_conversion_j: 0.01,
    }
}

fn buoyant_mass(s: Spec, t: [f64; 2], gz: f64) -> f64 {
    let rho = t.map(|t| 1000.0 + s.density_slope * (t - 300.0));
    let mu = t.map(|t| s.viscosity * 0.001 * (1.0 - 0.01 * (t - 300.0)));
    let head = 2.0 * gz * (rho[0] - rho[1]);
    let a = GEOMETRY_A * (mu[0] + mu[1]) / 1000.0;
    let b = 2.0 * GEOMETRY_B / 1000.0;
    (2.0 * head.abs() / (a + (a * a + 4.0 * b * head.abs()).sqrt())).copysign(head)
}

#[test]
fn evolving_source_temperatures_drive_circulation_reversal_and_one_conservative_heat_history() {
    for (gz, viscosity, coolant_initial_k) in [
        (-9.81, 1.0, [310.0, 295.0]),
        (9.81, 1.0, [310.0, 295.0]),
        (-9.81, 2.0, [310.0, 295.0]),
        (-9.81, 1.0, [295.0, 310.0]),
    ] {
        let s = Spec {
            circulating: true,
            density_slope: -1.0,
            viscosity,
            coolant_initial_k,
            ..spec()
        };
        fixture(
            s,
            |m, plain, electric, circuit, mut heat, prepared, law, initial, cx| {
                let gate = CancelGate::new_clock_free();
                let (_, material, native) = native(m);
                let mut owner = FrameRelaxationThermalSystem::new_shared_with_link_conductances(
                    vec![(material, DEST.to_vec())],
                    native,
                    &heat,
                    plain,
                    initial,
                    budget(),
                    &gate,
                )
                .unwrap();
                owner
                    .bind_material_circuit_with_link_conductances(
                        &heat,
                        electric,
                        circuit,
                        &[
                            (1, ElectricalPropertyBasis::Resistivity),
                            (2, ElectricalPropertyBasis::Resistivity),
                        ],
                        vec![],
                        initial,
                        &gate,
                    )
                    .unwrap();
                let network = SourceConvectionHeatNetwork {
                    heat: prepared,
                    cells: &[],
                    maximum_cells: 0,
                    resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t| {
                        law.at_temperature(t)
                    },
                };
                let sections = sections();
                let patches = patches();
                let gravity = closed(gz);
                let force = vec![0.0; owner.mechanical().mode_count()];
                let e0 = owner.mechanical().total_energy_j().unwrap()
                    + caloric(&heat)
                    + CHARGE * CHARGE / (2.0 * CAP);
                let mut accounted = 0.0;
                let mut first_mass = 0.0;
                let mut last_mass = 0.0;
                for step in 1..=12 {
                    let before = heat.cells().clone();
                    let report = network
                        .step_gravity_driven_relaxation_with_links(
                            cx,
                            &mut owner,
                            &mut heat,
                            &force,
                            Some(input()),
                            controls4(),
                            &[0.0; 4],
                            &gate,
                            gravity_path(&sections, &patches, &gravity),
                            |_, _, t| law.at_temperature(t),
                        )
                        .unwrap();
                    let payload = report.transport_source.as_ref().unwrap();
                    let t = [
                        heat.temperature_k(COOLANT).unwrap(),
                        heat.temperature_k(CellId(8)).unwrap(),
                    ];
                    let mass = buoyant_mass(s, t, gz);
                    near(payload.hydraulics.mass_flow_kg_s(), mass, 1e-11);
                    assert_eq!(
                        mass.is_sign_positive(),
                        gz * s.density_slope * (coolant_initial_k[0] - coolant_initial_k[1]) > 0.0
                    );
                    if step == 1 {
                        first_mass = mass;
                    }
                    last_mass = mass;
                    let q = mass / 1000.0;
                    let speed = q.abs() / FLOW_AREA;
                    let mut wall_power = 0.0;
                    for sample in payload.convection.finite.samples() {
                        let i = sample.wall_cell.0 as usize;
                        let wall = heat.temperature_k(sample.wall_cell).unwrap();
                        let film = 0.5 * (wall + t[0]);
                        let rho = 1000.0 + s.density_slope * (film - 300.0);
                        let mu = s.viscosity * 0.001 * (1.0 - 0.01 * (film - 300.0));
                        let k = s.k * (0.5 + 0.002 * (film - 300.0));
                        let h =
                            0.664 * (rho * speed * 0.1 / mu).sqrt() * (mu * s.cp / k).cbrt() * k
                                / 0.1;
                        let power = PATCHES[i] * h * (wall - t[0]);
                        near(sample.heat_into_fluid_w, power, 1e-8);
                        wall_power += power;
                    }
                    for i in 0..2 {
                        let id = sections[i].cell;
                        let cell = &heat.cells()[&id];
                        assert_eq!(cell.mass_kg, before[&id].mass_kg);
                        assert_eq!(cell.geometry, before[&id].geometry);
                        let carried = mass.abs() * s.cp * (t[1] - t[0]);
                        let expected = DT
                            * if i == 0 {
                                wall_power + carried
                            } else {
                                -carried
                            };
                        near(
                            cell.mass_kg * (cell.specific_h_j_kg - before[&id].specific_h_j_kg),
                            expected,
                            4e-8,
                        );
                        let actual = &payload.hydraulics.evaluations()[i];
                        assert_eq!(
                            payload.hydraulics.sections()[i].medium,
                            law.at_temperature(Temperature::new(t[i])).unwrap()
                        );
                        near(actual.volume_flow_m3_s, q, 1e-14);
                        assert_eq!(actual.momentum_density_kg_m3, 1000.0);
                        near(
                            payload.omitted_hydraulic_heat_w[&id],
                            actual.dissipated_power_w,
                            1e-14,
                        );
                    }
                    assert_eq!(payload.work.imposed_pressure_w, 0.0);
                    assert!(payload.work.buoyancy_w > 0.0 && payload.work.dissipation_w > 0.0);
                    near(payload.work.balance_residual_w, 0.0, 1e-12);
                    near(
                        payload.work.omitted_absolute_conversion_j,
                        DT * payload.work.omitted_absolute_conversion_w,
                        1e-15,
                    );
                    assert_eq!(
                        payload
                            .hydraulics
                            .hydrostatic()
                            .unwrap()
                            .gravitational_power_w(),
                        0.0
                    );
                    near(report.external_thermal_heat_j, 0.0, 1e-12);
                    assert_eq!(
                        report.enthalpy_flow.as_ref().unwrap().external_enthalpy_j,
                        0.0
                    );
                    let electric = report.electrical.as_ref().unwrap();
                    accounted += report.external_mechanical_loss_j
                        + report.mechanical.numerical_dissipation_j
                        - report.mechanical.external_work_j
                        - electric.native.supplied;
                    near(
                        owner.mechanical().total_energy_j().unwrap()
                            + caloric(&heat)
                            + electric.stored_after_j
                            - e0
                            + accounted,
                        0.0,
                        5e-7,
                    );
                    assert_eq!(heat.time_s(), step as f64 * DT);
                }
                assert!(
                    (last_mass - first_mass).abs() > 1e-8,
                    "buoyancy must refresh with accepted temperatures"
                );
            },
        );
    }
}

#[test]
fn uniform_temperature_or_zero_gravity_keeps_original_adiabatic_inventory_at_rest() {
    for (temperature, gz) in [([300.0; 2], -9.81), ([310.0, 295.0], 0.0)] {
        let s = Spec {
            powered: false,
            circulating: true,
            density_slope: -1.0,
            coolant_initial_k: temperature,
            adiabatic: true,
            ..spec()
        };
        fixture(s, |m, plain, _, _, mut heat, prepared, law, initial, cx| {
            let gate = CancelGate::new_clock_free();
            let (_, material, native) = native(m);
            let mut owner = FrameRelaxationThermalSystem::new_shared_with_link_conductances(
                vec![(material, DEST.to_vec())],
                native,
                &heat,
                plain,
                initial,
                budget(),
                &gate,
            )
            .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t| {
                    law.at_temperature(t)
                },
            };
            let sections = sections();
            let gravity = closed(gz);
            let path = gravity_path(&sections, &[], &gravity);
            let before = heat.cells().clone();
            let force = vec![0.0; owner.mechanical().mode_count()];
            let report = network
                .step_gravity_driven_relaxation_with_links(
                    cx,
                    &mut owner,
                    &mut heat,
                    &force,
                    None,
                    controls4(),
                    &[0.0; 4],
                    &gate,
                    path,
                    |_, _, t| law.at_temperature(t),
                )
                .unwrap();
            let payload = report.transport_source.as_ref().unwrap();
            assert_eq!(payload.hydraulics.mass_flow_kg_s(), 0.0);
            assert_eq!(payload.work.omitted_absolute_conversion_j, 0.0);
            for section in &sections {
                assert_eq!(heat.cells()[&section.cell], before[&section.cell]);
            }
            assert_eq!(heat.time_s(), DT);
        });
    }
}

#[test]
fn omitted_work_density_geometry_and_late_cancellation_preserve_all_histories_for_retry() {
    fixture(
        Spec {
            circulating: true,
            density_slope: -1.0,
            ..spec()
        },
        |m, plain, electric, circuit, mut heat, prepared, law, initial, cx| {
            let gate = CancelGate::new_clock_free();
            let (placement, material, native) = native(m);
            let mut owner = FrameRelaxationThermalSystem::new_shared_with_link_conductances(
                vec![(material, DEST.to_vec())],
                native,
                &heat,
                plain,
                initial,
                budget(),
                &gate,
            )
            .unwrap();
            owner
                .bind_material_circuit_with_link_conductances(
                    &heat,
                    electric,
                    circuit.clone(),
                    &[
                        (1, ElectricalPropertyBasis::Resistivity),
                        (2, ElectricalPropertyBasis::Resistivity),
                    ],
                    vec![],
                    initial,
                    &gate,
                )
                .unwrap();
            let network = SourceConvectionHeatNetwork {
                heat: prepared,
                cells: &[],
                maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t| {
                    law.at_temperature(t)
                },
            };
            let sections = sections();
            let patches = patches();
            let gravity = closed(-9.81);
            let open = SeriesGravity::try_new(
                [0.0, 0.0, -9.81],
                vec![[0.0; 3], [0.0, 0.0, 2.0], [1.0, 0.0, 0.0]],
            )
            .unwrap();
            let force = vec![0.0; owner.mechanical().mode_count()];
            let cells = heat.cells().clone();
            let memory = owner.mechanical().material_memory().to_vec();
            let physical = placement
                .physical_state(owner.mechanical().components(), &gate)
                .unwrap();
            for kind in 0..6 {
                let mut path = gravity_path(&sections, &patches, &gravity);
                if kind == 0 {
                    path.maximum_omitted_conversion_w = 0.0;
                }
                if kind == 1 {
                    path.maximum_omitted_conversion_j = 0.0;
                }
                if kind == 2 {
                    path.density = BoussinesqDensity::try_new(1000.0, 0.001).unwrap();
                }
                if kind == 3 {
                    path.density = BoussinesqDensity::try_new(1001.0, 0.05).unwrap();
                }
                if kind == 4 {
                    path.gravity = &open;
                }
                let cancellation = CancelGate::new_clock_free();
                let calls = std::cell::Cell::new(0);
                let result = network.step_gravity_driven_relaxation_with_links(
                    cx,
                    &mut owner,
                    &mut heat,
                    &force,
                    Some(input()),
                    controls4(),
                    &[0.0; 4],
                    &cancellation,
                    path,
                    |_, _, t| {
                        calls.set(calls.get() + 1);
                        if kind == 5 && calls.get() >= 4 {
                            cancellation.request();
                        }
                        law.at_temperature(t)
                    },
                );
                assert!(result.is_err(), "case {kind} must refuse");
                assert_eq!(heat.cells(), &cells);
                assert_eq!(heat.time_s(), 0.0);
                assert_eq!(owner.mechanical().samples_rendered(), 0);
                assert_eq!(owner.mechanical().material_memory(), memory);
                assert_eq!(owner.electrical().unwrap(), &circuit);
                assert_eq!(
                    placement
                        .physical_state(owner.mechanical().components(), &gate)
                        .unwrap(),
                    physical
                );
            }
            network
                .step_gravity_driven_relaxation_with_links(
                    cx,
                    &mut owner,
                    &mut heat,
                    &force,
                    Some(input()),
                    controls4(),
                    &[0.0; 4],
                    &gate,
                    gravity_path(&sections, &patches, &gravity),
                    |_, _, t| law.at_temperature(t),
                )
                .unwrap();
            assert_eq!(owner.mechanical().samples_rendered(), 1);
        },
    );
}
