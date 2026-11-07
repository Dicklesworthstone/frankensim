//! Actual finite coolant transport reaches source-limited native motor torque.
use super::*;
use fs_couple::circuit_enthalpy::motor::{MotorWinding, heat::MotorHeatNetwork};
use fs_mbd::rotordynamics::{
    RotorBearingModel,
    transient::{RotorTransientState, drive::{RotorCurrentDrive, RotorDrivePowerControl}},
};

const FLUID_CP: f64 = 4000.0;
const INPUT_W: f64 = 3.0;
const EFFICIENCY: f64 = 0.8;

// Independent free-inertia, power-limit, two-body caloric and flat-plate equations.
fn reference_motor(viscosity_scale: f64) -> [f64; 4] {
    let evaluate = |wall: f64| {
        let resistance = 1.0 + 0.02 * (wall - 300.0);
        let current = (INPUT_W * EFFICIENCY / (resistance + DT / 2.0)).sqrt();
        let shaft = DT * DT * current * current / 4.0;
        let heat = DT * current * current * resistance / EFFICIENCY
            + (1.0 / EFFICIENCY - 1.0) * shaft;
        let fluid = 290.0 + (heat - MASS * CP * (wall - 310.0)) / (MASS * FLUID_CP);
        let film = 0.5 * (wall + fluid);
        let k = 0.5 + 0.002 * (film - 300.0);
        let mu = viscosity_scale * 0.001 * (1.0 - 0.01 * (film - 300.0));
        let h = 0.664 * (1000.0 * 0.1 / mu).sqrt()
            * (mu * FLUID_CP / k).cbrt() * k / 0.1;
        let into_fluid = SURFACE * h * (wall - fluid);
        (MASS * CP * (wall - 310.0) - heat + DT * into_fluid,
            [wall, fluid, current, into_fluid])
    };
    let (mut lower, mut upper) = (290.0, 312.0);
    assert!(evaluate(lower).0 < 0.0 && evaluate(upper).0 > 0.0);
    for _ in 0..80 {
        let wall = 0.5 * (lower + upper);
        if evaluate(wall).0 > 0.0 { upper = wall; } else { lower = wall; }
    }
    let (residual, result) = evaluate(0.5 * (lower + upper));
    assert!(residual.abs() < 1e-10);
    result
}

#[test]
fn finite_coolant_viscosity_reaches_winding_resistance_source_limit_and_motion() {
    let mut endpoints = Vec::new();
    for viscosity_scale in [1.0, 4.0] {
        let (body_card, body_curve) = body_material();
        let (fluid_card, fluid_curve) = fluid_material(viscosity_scale, FLUID_CP);
        let body_context = body_point(310.0);
        let fluid_context = fluid_point(290.0);
        let thermal = initial_pair(&body_curve, &fluid_curve, FLUID_CP, false, 0.0);
        let saved = thermal.cells().clone();
        let fluid = liquid_law(&fluid_card);
        let cells = [
            MaterialCell {
                cell: BODY, material_id: MaterialId(0), curve: body_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &body_card, context: &body_context,
                    temperature_axis: "T", liquid_fraction_axis: None,
                },
                conductivity: None,
            },
            MaterialCell {
                cell: FLUID, material_id: MaterialId(1), curve: fluid_curve.curve(),
                point: ThermalMaterialPoint {
                    material: &fluid_card, context: &fluid_context,
                    temperature_axis: "T", liquid_fraction_axis: None,
                },
                conductivity: None,
            },
        ];
        let bulk = MaterialNetwork { cells: &cells, resistors: &[], links: &[] };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk, contacts: &[], electrical_bases: &[],
        };
        let heat = MaterialHeatNetwork {
            electrical: &electrical, thermal_contacts: &[], radiation: None,
        };
        let links = [binding()];
        let initial = evaluate_source_link_convection(
            &thermal, &links, 1, |_, _, t| fluid.at_temperature(t), || false,
        ).unwrap();
        let gate = fs_exec::CancelGate::new();
        fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default()).scope(|arena| {
            let cx = fs_exec::Cx::new(
                &gate, arena,
                fs_exec::StreamKey { seed: 173, kernel_id: 29, tile: 0, iteration: 0 },
                fs_exec::Budget::INFINITE, fs_exec::ExecMode::Deterministic,
            );
            let prepared = heat.prepare_with_link_conductances(
                &cx, &thermal, initial.conductances_w_k(),
            ).unwrap();
            let convection = SourceConvectionHeatNetwork {
                heat: &prepared, cells: &[], maximum_cells: 0,
                resolve_fluid: |_: &SourceConvectiveCell, _: &DepositionSimulation<'_>, t: Temperature| {
                    fluid.at_temperature(t)
                },
            };
            let assembly = convection.with_links(&links, 1, |_, _, t| fluid.at_temperature(t));
            // Both native loss channels are explicitly homogenized in this body.
            let network = MotorHeatNetwork {
                assembly: &assembly,
                winding: MotorWinding {
                    cell: BODY, basis: ElectricalPropertyBasis::Resistivity,
                    selection: MaterialPropertySelection::SingleClaimOnly,
                    length_m: 0.1, area_m2: 1e-5,
                },
                bridge_cell: BODY,
            };
            let rotor = RotorBearingModel::new(1, &[2.0], &[0.0], &[0.0], &[0.0], || false).unwrap();
            let native = rotor.at_speed(0.0, || false).unwrap()
                .midpoint_stepper(DT, 1e-12, || false).unwrap();
            let motor = network.bind(&cx, &thermal, &native,
                RotorTransientState { displacement: vec![0.0], velocity: vec![0.0] },
                RotorCurrentDrive {
                    shaft_coordinate: 0, torque_constant_nm_per_a: 1.0,
                    winding_resistance_ohm: 999.0, maximum_current_a: 10.0,
                    maximum_terminal_voltage_v: 100.0,
                    minimum_speed_rad_s: -100.0, maximum_speed_rad_s: 100.0,
                    motoring_efficiency: EFFICIENCY, regenerative_efficiency: 0.6,
                },
            ).unwrap();
            let mut solve = controls();
            solve.maximum_cells = 2;
            solve.newton.linear_restart = 2;
            let next = network.step_candidate(
                &cx, &thermal, DT, &motor, &[0.0], 2.0,
                RotorDrivePowerControl {
                    maximum_source_input_w: INPUT_W, maximum_source_return_w: 1000.0,
                    absolute_power_tolerance_w: 1e-9, absolute_voltage_tolerance_v: 1e-10,
                },
                solve, &[0.0; 2],
            ).unwrap();
            let expected = reference_motor(viscosity_scale);
            let wall = next.thermal().temperature_k(BODY).unwrap();
            let fluid_t = next.thermal().temperature_k(FLUID).unwrap();
            let report = next.report();
            let native_step = &report.electrical.drive;
            let finite = report.source.network.finite_convection.as_ref().unwrap();
            let source = &finite.samples()[0];
            near(wall, expected[0], 1e-8);
            near(fluid_t, expected[1], 1e-8);
            near(native_step.current_a, expected[2], 1e-8);
            near(source.heat_into_fluid_w, expected[3], 1e-7);
            near(next.motor().state().velocity[0], DT * expected[2] / 2.0, 1e-9);
            near(next.motor().parameters().winding_resistance_ohm, 1.0 + 0.02 * (wall - 300.0), 1e-10);
            near(source.property_temperature_k, 0.5 * (wall + fluid_t), 1e-12);
            assert_eq!(source.material.material().card_identity(), fluid_card.content_hash());
            near(MASS * FLUID_CP * (fluid_t - 290.0), DT * source.heat_into_fluid_w, 1e-8);
            near(MASS * CP * (wall - 310.0), report.deposited_joule_j - DT * source.heat_into_fluid_w, 1e-8);
            near(report.non_electrical_heat_j, 0.0, 1e-12);
            near(next.thermal().energy_j() - thermal.energy_j()
                + native_step.motion.kinetic_energy_change_j, native_step.electrical_draw_j, 1e-8);
            near(report.total_energy_residual_j, 0.0, 1e-8);
            assert_eq!(thermal.cells(), &saved);
            assert_eq!(motor.state().velocity, vec![0.0]);
            endpoints.push([wall, native_step.current_a]);
        });
    }
    assert!(endpoints[1][0] > endpoints[0][0]);
    assert!(endpoints[1][1] < endpoints[0][1]);
}
