//! Independently integrated physical AC electrodes with finite temperature.
use super::*;
use fs_couple::circuit_enthalpy::spatial_induction::harmonic::terminals::{
    current_heat_assembly_candidate, voltage_heat_assembly_candidate,
};
use fs_em::spatial::vector_potential::induction::terminals::{
    SpatialTerminalInduction, current::SpatialCurrentInduction,
};

#[derive(Clone, Copy)]
enum Drive {
    Voltage,
    Current,
}
fn amplitude(mode: Drive) -> C64 {
    match mode {
        Drive::Voltage => C64::new(2.5, -1.5),
        Drive::Current => C64::new(8.0, -4.0),
    }
}
fn impressed() -> C64 {
    C64::new(4.0, 3.0)
}
fn square(value: C64) -> f64 {
    value.re * value.re + value.im * value.im
}
fn complex_near(actual: C64, expected: C64, tolerance: f64) {
    near(actual.re, expected.re, tolerance);
    near(actual.im, expected.im, tolerance);
}
fn preparation() -> SpatialInductionLimits {
    SpatialInductionLimits {
        maximum_unknowns: 1000,
        maximum_conducting_cells: 100,
        maximum_assembly_terms: 10000,
    }
}
fn electrodes(problem: &SpatialVectorMagnetostatics) -> Vec<Vec<usize>> {
    [0.0, 1.0]
        .into_iter()
        .map(|x| {
            problem
                .exterior_faces()
                .into_iter()
                .filter(|&face| {
                    problem.mesh().faces[face]
                        .iter()
                        .all(|&vertex| problem.positions_m()[vertex as usize][0] == x)
                })
                .collect()
        })
        .collect()
}

struct Oracle {
    voltage: C64,
    current: C64,
    heat_w: f64,
    terminal_va: C64,
    impressed_va: C64,
    magnetic_j: f64,
}
// One Whitney mode W: integral W_x=1/6, integral |W|²=1/5, K=4/mu.
// Fixed electrodes give phi=-U*x and I=sigma*(U-Adot/6). The impressed
// uniform vector supplies F=J/2. Elimination under total-current control
// changes the dynamic mass to sigma*(1/5-1/36), not a fitted AC resistance.
fn oracle(temperature: f64, mode: Drive, impressed: C64, drive: C64) -> Oracle {
    let sigma = 1.0 / (RHO * (1.0 + ALPHA * (temperature - 300.0)));
    let f = impressed.scale(0.5);
    let a = match mode {
        Drive::Voltage => (f + drive.scale(sigma / 6.0)) / C64::new(4.0, -OMEGA * sigma / 5.0),
        Drive::Current => {
            (f + drive.scale(1.0 / 6.0)) / C64::new(4.0, -OMEGA * sigma * (1.0 / 5.0 - 1.0 / 36.0))
        }
    };
    let rate = C64::new(0.0, -OMEGA) * a;
    let (voltage, current) = match mode {
        Drive::Voltage => (drive, (drive - rate.scale(1.0 / 6.0)).scale(sigma)),
        Drive::Current => (drive.scale(1.0 / sigma) + rate.scale(1.0 / 6.0), drive),
    };
    Oracle {
        voltage,
        current,
        heat_w: 0.5
            * sigma
            * (square(voltage) - (voltage * rate.conj()).re / 3.0 + square(rate) / 5.0),
        terminal_va: (voltage * current.conj()).scale(0.5),
        impressed_va: (rate * f.conj()).scale(0.5),
        magnetic_j: square(a),
    }
}
fn reference(old: [f64; 2], mode: Drive, frozen: Option<f64>) -> [f64; 2] {
    let sink = |ta| {
        root(200.0, 1000.0, |tb| {
            4.0 * (tb - old[1]) + DT * transfer(ta, tb, 0.5)
        })
    };
    let ta = root(200.0, 1000.0, |ta| {
        let heat = frozen.unwrap_or_else(|| oracle(ta, mode, impressed(), amplitude(mode)).heat_w);
        4.0 * (ta - old[0]) - DT * (heat + transfer(ta, sink(ta), 0.5))
    });
    [ta, sink(ta)]
}

#[allow(clippy::too_many_arguments)]
fn terminal_step<'t, A: HeatAssembly>(
    heat: &A,
    cx: &Cx<'_>,
    thermal: &DepositionSimulation<'t>,
    target: f64,
    mode: Drive,
    voltage: &SpatialTerminalInduction<'_>,
    current: &SpatialCurrentInduction<'_>,
    terminals: &[C64],
    controls: CircuitEnthalpySolve,
) -> Result<
    SpatialHarmonicThermalCandidate<'t, A::Sample>,
    CircuitEnthalpyError<SpatialInductionHeatError<A::Error>>,
> {
    match mode {
        Drive::Voltage => voltage_heat_assembly_candidate(
            heat,
            cx,
            thermal,
            target,
            voltage,
            terminals,
            &[[impressed(); 3]; 6],
            OMEGA,
            &regions(),
            1e-10,
            controls,
            ac_controls(),
            assembly(),
            &[0.0; 2],
        ),
        Drive::Current => current_heat_assembly_candidate(
            heat,
            cx,
            thermal,
            target,
            current,
            terminals,
            &[[impressed(); 3]; 6],
            OMEGA,
            &regions(),
            1e-10,
            controls,
            ac_controls(),
            assembly(),
            &[0.0; 2],
        ),
    }
}
fn terminals(mode: Drive) -> [C64; 2] {
    let drive = amplitude(mode);
    match mode {
        Drive::Voltage => [drive, C64::ZERO],
        Drive::Current => [drive, -drive],
    }
}

#[test]
fn physical_voltage_and_current_heating_follow_independent_coupled_ac_oracles() {
    for mode in [Drive::Voltage, Drive::Current] {
        with_case(false, 0.5, |cx, _, heat, mut thermal, card| {
            let prepared = heat.prepare(cx, &thermal).unwrap();
            let problem = problem(card, BASIS);
            let voltage = SpatialTerminalInduction::new(
                initial(&problem, BASIS, 0.7, 1.0),
                &electrodes(&problem),
                preparation(),
                &mut || false,
            )
            .unwrap();
            let current = voltage
                .clone()
                .into_current_driven(preparation(), &mut || false)
                .unwrap();
            let old_a = voltage.potential_wb().to_vec();
            let old_sources = voltage.conductors().to_vec();
            let mut temperatures = [300.0, 340.0];
            for interval in 1..=2 {
                let expected = reference(temperatures, mode, None);
                let frozen = reference(
                    temperatures,
                    mode,
                    Some(oracle(temperatures[0], mode, impressed(), amplitude(mode)).heat_w),
                );
                assert!((expected[0] - frozen[0]).abs() > 0.01);
                let expected_field = oracle(expected[0], mode, impressed(), amplitude(mode));
                let separate = oracle(expected[0], mode, C64::ZERO, amplitude(mode)).heat_w
                    + oracle(expected[0], mode, impressed(), C64::ZERO).heat_w;
                assert!(
                    (expected_field.heat_w - separate).abs() > 0.01,
                    "coherent excitation must retain cross heat"
                );
                let old_cells = thermal.cells().clone();
                let next = terminal_step(
                    &prepared,
                    cx,
                    &thermal,
                    f64::from(interval) * DT,
                    mode,
                    &voltage,
                    &current,
                    &terminals(mode),
                    thermal_controls(),
                )
                .unwrap();
                for i in 0..2 {
                    near(
                        next.thermal().temperature_k([BODY, SINK][i]).unwrap(),
                        expected[i],
                        2e-8,
                    );
                }
                let report = next.report();
                let field = &report.electrical;
                let b = field.balance;
                near(b.mean_ohmic_power_w, expected_field.heat_w, 2e-9);
                near(b.mean_terminal_power_w, expected_field.terminal_va.re, 2e-9);
                near(
                    b.reactive_terminal_power_var,
                    expected_field.terminal_va.im,
                    2e-9,
                );
                near(
                    b.mean_impressed_power_w,
                    expected_field.impressed_va.re,
                    2e-9,
                );
                near(
                    b.reactive_impressed_power_var,
                    expected_field.impressed_va.im,
                    2e-9,
                );
                near(b.mean_magnetic_energy_j, expected_field.magnetic_j, 2e-10);
                near(b.power_residual_w, 0.0, 1e-9);
                near(b.reactive_power_residual_var, 0.0, 1e-9);
                near(b.mean_mechanical_power_w, 0.0, 1e-15);
                near(b.mean_supplied_power_w, expected_field.heat_w, 2e-9);
                assert_eq!(field.terminals.len(), 2);
                for (index, terminal) in field.terminals.iter().enumerate() {
                    assert_eq!(terminal.terminal, index);
                    assert_eq!(terminal.reference_terminal, 0);
                }
                assert_eq!(field.terminals[0].component, field.terminals[1].component);
                complex_near(
                    field.terminals[0].incoming_current_a,
                    expected_field.current,
                    2e-9,
                );
                complex_near(
                    field.terminals[1].incoming_current_a,
                    -expected_field.current,
                    2e-9,
                );
                complex_near(
                    field.terminals[1].voltage_to_reference_v,
                    -expected_field.voltage,
                    2e-9,
                );
                complex_near(
                    field.terminals[0].supplied_power_va + field.terminals[1].supplied_power_va,
                    expected_field.terminal_va,
                    2e-9,
                );
                assert!(field.maximum_charge_residual_a < ac_controls().current_tolerance_a);
                assert!(field.field.converged_euclidean());
                assert_eq!(field.cells.len(), 6);
                for (index, cell) in field.cells.iter().enumerate() {
                    assert_eq!(cell.tetrahedron, index);
                    assert_eq!(cell.cell, BODY);
                }
                near(field.heat_j, DT * expected_field.heat_w, 5e-9);
                near(report.deposited_joule_j, field.heat_j, 5e-9);
                near(report.non_electrical_heat_j, 0.0, 2e-8);
                near(
                    next.thermal().energy_j() - thermal.energy_j(),
                    field.supplied_work_j,
                    2e-8,
                );
                near(report.total_energy_residual_j, 0.0, 2e-8);
                for source in &report.source.magnetic {
                    near(
                        *source.material.query_point().get("T").unwrap(),
                        expected[0],
                        2e-8,
                    );
                    near(
                        source
                            .material
                            .property(BASIS.property_name())
                            .unwrap()
                            .value_si(),
                        RHO * (1.0 + ALPHA * (expected[0] - 300.0)),
                        1e-10,
                    );
                    assert_eq!(
                        source.material.card_identity(),
                        old_sources[source.tetrahedron]
                            .as_ref()
                            .unwrap()
                            .material()
                            .card_identity()
                    );
                }
                assert_eq!(thermal.cells(), &old_cells);
                assert_eq!(voltage.potential_wb(), old_a);
                assert_eq!(current.potential_wb(), old_a);
                assert_eq!(voltage.conductors(), old_sources);
                assert_eq!(current.conductors(), old_sources);
                assert_eq!([voltage.time_s(), current.time_s()], [0.0; 2]);
                temperatures = expected;
                thermal = next.into_parts().0;
            }
        });
    }
}

#[test]
fn common_voltage_mode_and_failed_electrode_trials_do_not_change_either_history() {
    with_case(false, 0.5, |cx, gate, heat, thermal, card| {
        let prepared = heat.prepare(cx, &thermal).unwrap();
        let problem = problem(card, BASIS);
        let voltage = SpatialTerminalInduction::new(
            initial(&problem, BASIS, 0.7, 1.0),
            &electrodes(&problem),
            preparation(),
            &mut || false,
        )
        .unwrap();
        let current = voltage
            .clone()
            .into_current_driven(preparation(), &mut || false)
            .unwrap();
        let old_a = voltage.potential_wb().to_vec();
        let old_cells = thermal.cells().clone();
        let baseline = terminal_step(
            &prepared,
            cx,
            &thermal,
            DT,
            Drive::Voltage,
            &voltage,
            &current,
            &terminals(Drive::Voltage),
            thermal_controls(),
        )
        .unwrap();
        let offset = C64::new(1073741824.0, -536870912.0);
        let shifted = terminal_step(
            &prepared,
            cx,
            &thermal,
            DT,
            Drive::Voltage,
            &voltage,
            &current,
            &[amplitude(Drive::Voltage) + offset, offset],
            thermal_controls(),
        )
        .unwrap();
        assert_eq!(baseline.thermal().cells(), shifted.thermal().cells());
        assert_eq!(
            baseline.report().electrical.terminals,
            shifted.report().electrical.terminals
        );
        assert_eq!(
            baseline.report().electrical.heat_j.to_bits(),
            shifted.report().electrical.heat_j.to_bits()
        );
        let mut budget = thermal_controls();
        budget.maximum_circuit_trials = 1;
        assert!(matches!(
            terminal_step(
                &prepared,
                cx,
                &thermal,
                DT,
                Drive::Voltage,
                &voltage,
                &current,
                &terminals(Drive::Voltage),
                budget
            ),
            Err(CircuitEnthalpyError::CircuitBudget)
        ));
        let drive = amplitude(Drive::Current);
        assert!(
            terminal_step(
                &prepared,
                cx,
                &thermal,
                DT,
                Drive::Current,
                &voltage,
                &current,
                &[drive, -drive + C64::ONE],
                thermal_controls()
            )
            .is_err()
        );
        let tagged = Tagged {
            inner: &prepared,
            gate,
            calls: Cell::new(0),
            fail_at: Cell::new(Some(3)),
            cancel_at: Cell::new(None),
        };
        assert!(matches!(
            terminal_step(
                &tagged,
                cx,
                &thermal,
                DT,
                Drive::Current,
                &voltage,
                &current,
                &terminals(Drive::Current),
                thermal_controls()
            ),
            Err(CircuitEnthalpyError::Source(
                SpatialInductionHeatError::Heat(TaggedError::Refused)
            ))
        ));
        tagged.calls.set(0);
        tagged.fail_at.set(None);
        let retry = terminal_step(
            &tagged,
            cx,
            &thermal,
            DT,
            Drive::Current,
            &voltage,
            &current,
            &terminals(Drive::Current),
            thermal_controls(),
        )
        .unwrap();
        assert_eq!(retry.report().source.network.resolution, tagged.calls.get());
        tagged.calls.set(0);
        tagged.cancel_at.set(Some(4));
        assert!(matches!(
            terminal_step(
                &tagged,
                cx,
                &thermal,
                DT,
                Drive::Current,
                &voltage,
                &current,
                &terminals(Drive::Current),
                thermal_controls()
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(thermal.cells(), &old_cells);
        assert_eq!(voltage.potential_wb(), old_a);
        assert_eq!(current.potential_wb(), old_a);
        assert_eq!([voltage.time_s(), current.time_s()], [0.0; 2]);
    });
}
