//! Independent two-temperature resistor/field/material feedback.
use super::*;
use fs_couple::circuit_enthalpy::{
    material::CellResistor, spatial_induction::circuit::SpatialTerminalCircuitThermalNetwork,
};
use fs_material::state_point::MaterialPropertySelection;

fn two_material_bodies(curves: &[ResolvedHeatCapacityCurve; 2]) -> DepositionSimulation<'_> {
    let events = [0, 1]
        .into_iter()
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(7 + i),
                segment_id: i,
                body_id: i,
                material_id: MaterialId(9 + i),
                mass_kg: 1.0,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [2.0 * i as f64, 0.0, 0.0],
                    end_m: [2.0 * i as f64 + 1.0, 0.0, 0.0],
                    volume_m3: 1.0,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[
            (MaterialId(9), curves[0].curve() as &dyn ThermalLaw),
            (MaterialId(10), curves[1].curve() as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 2,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn material_cells<'a>(
    cards: &'a [MaterialCard; 2],
    curves: &'a [ResolvedHeatCapacityCurve; 2],
    context: &'a fs_matdb::QueryPoint,
) -> [MaterialCell<'a>; 2] {
    [0, 1].map(|i| MaterialCell {
        cell: CellId(7 + i as u64),
        material_id: MaterialId(9 + i as u64),
        curve: curves[i].curve(),
        point: ThermalMaterialPoint {
            material: &cards[i],
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    })
}

fn load() -> CellResistor {
    CellResistor {
        destination: ResistorCell {
            branch: 1,
            cell: CellId(8),
        },
        selection: MaterialPropertySelection::SingleClaimOnly,
        length_m: 1.0,
        area_m2: 1.0,
    }
}

#[test]
fn both_source_bases_close_actual_resistor_and_conductor_temperature_feedback() {
    let field_basis = ElectricalPropertyBasis::Resistivity;
    let slope = 0.0002;
    let dt = 0.1;
    for load_basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let cards = [
            source::absolute_affine_magnetic_card(0.2, 0.004, 1.0, slope, field_basis),
            source::absolute_affine_magnetic_card(0.4, 0.008, 1.0, 0.0, load_basis),
        ];
        let curves = cards.each_ref().map(curve);
        let problem = problem(&cards[0], field_basis);
        let mut circuit = connect(
            initial(&problem, field_basis, 0.0, 1.0),
            series_graph(0.4),
            vec![],
            vec![0.0],
            0.7,
        );
        let mut heat = two_material_bodies(&curves);
        let context = source::absolute_point(300.0);
        let cells = material_cells(&cards, &curves, &context);
        let regions = regions();
        let laws = laws(&cards[0], &problem, 5.0);
        let destinations = [load().destination];
        let resistors = [load()];
        let network = SpatialTerminalCircuitThermalNetwork {
            magnetic: SpatialMagneticThermalNetwork {
                cells: &cells,
                links: &[],
                magnetic_regions: &regions,
                volume_tolerance_m3: 1e-12,
            },
            graph_resistors: &destinations,
        };
        let (mut a0, mut t0, mut load0) = (0.0, 300.0, 300.0);
        for (index, drive) in [12.0, 0.0, -8.0].into_iter().enumerate() {
            // Independent exact Whitney integrals: int W_x=1/6,
            // int |W|²=1/5 and curl stiffness=4/mu. Solve the two
            // finite caloric equations without any native circuit/heat calls.
            let old_mu = 1.0 + slope * (t0 - 300.0);
            let at = |t: f64, load_t: f64| {
                let mu = 1.0 + slope * (t - 300.0);
                let conductivity = sigma(t, 0.2, 0.004, field_basis);
                let resistance = sigma(load_t, 0.4, 0.008, load_basis).recip();
                let f = conductivity / 6.0;
                let d = conductivity / 5.0 + 2.0 * dt / mu;
                let magnetic = 2.0 * (old_mu.recip() + mu.recip()) * a0;
                let y = conductivity - f * f / d;
                let i0 = f * magnetic / d;
                let voltage = (drive - resistance * i0) / (1.0 + resistance * y);
                let current = i0 + y * voltage;
                let rate = (f * voltage - magnetic) / d;
                let a = a0 + dt * rate;
                let joule = dt
                    * conductivity
                    * (voltage * voltage - voltage * rate / 3.0 + rate * rate / 5.0);
                let parameter = 2.0 * a0 * a * (mu.recip() - old_mu.recip());
                let ts = 2.0 * slope * (t * a * a / (mu * mu) - t0 * a0 * a0 / (old_mu * old_mu));
                let resistor_heat = dt * resistance * current * current;
                (
                    a,
                    voltage,
                    current,
                    joule,
                    resistor_heat,
                    -parameter - ts,
                    resistance,
                )
            };
            let response = |t: f64| {
                let load_t = bracket(load0, load0 + 20.0, |load_t| {
                    4.0 * (load_t - load0) - at(t, load_t).4
                });
                (load_t, at(t, load_t))
            };
            let t = bracket(t0, t0 + 20.0, |t| {
                let (_, (_, _, _, joule, _, reversible, _)) = response(t);
                4.0 * (t - t0) - joule - reversible
            });
            let (load_t, (a, voltage, current, joule, resistor_heat, reversible, resistance)) =
                response(t);
            let proposal = network
                .step_linear_material_candidate(
                    &heat,
                    (index + 1) as f64 * dt,
                    &circuit,
                    &[drive],
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &resistors,
                    &[(1, load_basis)],
                    &laws,
                    thermal_controls(),
                    circuit_controls(),
                    assembly(),
                    &[0.0, 0.0],
                    || false,
                )
                .unwrap();
            let step = &proposal.report().coupled.electrical;
            near(
                proposal.thermal().temperature_k(CellId(7)).unwrap(),
                t,
                8e-8,
            );
            near(
                proposal.thermal().temperature_k(CellId(8)).unwrap(),
                load_t,
                8e-8,
            );
            near(
                proposal.circuit().field().potential_wb()[diagonal(&problem)],
                a,
                8e-8,
            );
            near(step.terminal_voltage_v[0], voltage, 1e-7);
            near(step.terminal_current_a[0], current, 2e-7);
            near(step.field.heat_j, joule, 2e-7);
            near(step.resistor_heat_j, resistor_heat, 2e-7);
            near(step.graph_resistors[0].heat_j, resistor_heat, 2e-7);
            near(
                step.graph_materials[0].conductor.resistance_ohm(),
                resistance,
                1e-9,
            );
            near(
                proposal.report().coupled.source.reversible_heat_j,
                reversible,
                3e-9,
            );
            near(proposal.deposited_heat_j(), joule + resistor_heat, 2e-7);
            near(step.external_work_j, dt * drive * current, 2e-7);
            near(proposal.report().total_energy_residual_j, 0.0, 2e-8);
            assert!(resistance > 0.4);
            assert_eq!(
                proposal.circuit().graph_resistor_materials(),
                step.graph_materials
            );
            near(heat.temperature_k(CellId(7)).unwrap(), t0, 8e-8);
            near(heat.temperature_k(CellId(8)).unwrap(), load0, 8e-8);
            let parts = proposal.into_parts();
            heat = parts.0;
            circuit = parts.1;
            a0 = a;
            t0 = t;
            load0 = load_t;
        }
    }
}

#[test]
fn retained_resistor_sources_cannot_be_frozen_or_reassigned_after_heating() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let cards = [
        source::absolute_affine_magnetic_card(0.2, 0.004, 1.0, 0.0, basis),
        source::absolute_affine_magnetic_card(0.4, 0.008, 1.0, 0.0, basis),
    ];
    let curves = cards.each_ref().map(curve);
    let problem = problem(&cards[0], basis);
    let circuit = connect(
        initial(&problem, basis, 0.0, 1.0),
        series_graph(0.4),
        vec![],
        vec![0.0],
        0.7,
    );
    let heat = two_material_bodies(&curves);
    let context = source::absolute_point(300.0);
    let cells = material_cells(&cards, &curves, &context);
    let regions = regions();
    let laws = laws(&cards[0], &problem, 5.0);
    let destinations = [load().destination];
    let network = SpatialTerminalCircuitThermalNetwork {
        magnetic: SpatialMagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &regions,
            volume_tolerance_m3: 1e-12,
        },
        graph_resistors: &destinations,
    };
    let (heat, circuit, _) = network
        .step_linear_material_candidate(
            &heat,
            0.1,
            &circuit,
            &[12.0],
            &[[0.0; 3]; 6],
            &[[0.0; 3]; 8],
            &[load()],
            &[(1, basis)],
            &laws,
            thermal_controls(),
            circuit_controls(),
            assembly(),
            &[0.0, 0.0],
            || false,
        )
        .unwrap()
        .into_parts();
    let saved = (
        heat.cells().clone(),
        circuit.native_state().to_vec(),
        circuit.field().potential_wb().to_vec(),
        circuit.graph_resistor_materials().to_vec(),
    );
    // A source-owned resistor may not silently become a held numerical R.
    assert!(
        network
            .magnetic
            .step_terminal_circuit_candidate(
                &heat,
                0.2,
                &circuit,
                &[12.0],
                &[[0.0; 3]; 6],
                &[[0.0; 3]; 8],
                &destinations,
                &laws,
                thermal_controls(),
                circuit_controls(),
                assembly(),
                &[0.0, 0.0],
                || false,
            )
            .is_err()
    );
    for (length_m, area_m2) in [(0.5, 1.0), (2.0, 0.5)] {
        // First changes the real body volume; second preserves volume but
        // changes the accepted current path and therefore original resistance.
        let bad = CellResistor {
            length_m,
            area_m2,
            ..load()
        };
        assert!(
            network
                .step_linear_material_candidate(
                    &heat,
                    0.2,
                    &circuit,
                    &[12.0],
                    &[[0.0; 3]; 6],
                    &[[0.0; 3]; 8],
                    &[bad],
                    &[(1, basis)],
                    &laws,
                    thermal_controls(),
                    circuit_controls(),
                    assembly(),
                    &[0.0, 0.0],
                    || false,
                )
                .is_err()
        );
    }
    assert_eq!(heat.cells(), &saved.0);
    assert_eq!(circuit.native_state(), saved.1);
    assert_eq!(circuit.field().potential_wb(), saved.2);
    assert_eq!(circuit.graph_resistor_materials(), saved.3);
    assert_eq!(heat.time_s(), 0.1);
    assert_eq!(circuit.time_s(), 0.1);
}
