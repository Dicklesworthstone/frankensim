use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatError, ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

const DX: f64 = 0.01;
const DT: f64 = 0.01;
fn transport(cp: f64) -> FluidTransportState { FluidTransportState::try_new(1000.0, cp, 100.0, 1.0).unwrap() }
fn slab3(cp: f64) -> MaterialSlab3 {
    plan_scaling_from_material(&transport(cp), DX, DT, 0.01).unwrap().thermal_slab3(MaterialThermalSlab3 {
        nx: 2, ny: 1, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&transport(1000.0), DX, DT, 0.01).unwrap().thermal_slab(MaterialThermalSlab {
        nx: 2, ny: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    }).unwrap()
}
fn electrodes(normal: usize, voltage: f64) -> ElectricSetup {
    let (lo, hi, upper_lo, upper_hi) = if normal == 2 {
        ([0,0,1], [2,1,2], [0,0,6], [2,1,7])
    } else { ([0,1,0], [2,2,1], [0,6,0], [2,7,1]) };
    ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo, hi, drive: ElectrodeDrive::Voltage(voltage) },
        Electrode { lo: upper_lo, hi: upper_hi, drive: ElectrodeDrive::Voltage(0.0) },
    ] }
}
fn limits() -> ElectricHeatLimits { ElectricHeatLimits {
    thermal: HeatSourceLimits { max_cells: 16, temperature_range_k: [290.0, 320.0], balance_tolerance_j: 1e-9 },
    coupling_tolerance_j: 1e-10,
} }
fn input(setup: &ElectricSetup) -> ElectricHeatInput<'_> { ElectricHeatInput { setup, contacts: &[], limits: limits() } }
fn close(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a} != {b}"); }

#[test]
fn actual_current_solve_drives_actual_flow_heat_and_keeps_cell_power() {
    let setup = electrodes(2, 2.0);
    let mut coupled = slab3(1000.0);
    let mut separate = slab3(1000.0);
    let force = vec![[1.0, 0.0, 0.0]; 16];
    let report = coupled.step_electrothermal(input(&setup), &force, &CancelGate::new(),
        |_, _, _| Ok(Some(1.0))).unwrap();
    // Two parallel paths: each contains 4*dx of bulk, cross section dx^2.
    let expected_power = 2.0 * 2.0_f64.powi(2) / (4.0 * DX / DX.powi(2));
    close(report.electric.total_joule_w, expected_power, 1e-12);
    close(report.electrode_energy_j, expected_power * DT, 1e-12);
    let power: Vec<_> = report.electric.joule_w.iter().map(|w| w / DX.powi(3)).collect();
    let direct = separate.step_with_power_density(&power, &force, limits().thermal).unwrap();
    assert_eq!(report.heat, direct);
    assert_eq!(coupled.lattice.grid.f, separate.lattice.grid.f);
    assert_eq!(coupled.lattice.gpop, separate.lattice.gpop);
    assert!(report.combined_defect_j.abs() < 1e-9);
    assert!(coupled.lattice.grid.moments3(8).1[0] > 0.0);
}

#[test]
fn data_only_heat_capacity_and_resistivity_substitutions_change_response() {
    let setup = electrodes(2, 10.0);
    let off = electrodes(2, 0.0);
    let mut a = slab3(1000.0);
    let mut b = slab3(2000.0);
    let mut cold_a = slab3(1000.0);
    let mut cold_b = slab3(2000.0);
    let mut resistant = slab3(1000.0);
    let force = vec![[0.0;3];16];
    let gate = CancelGate::new();
    let ra = a.step_electrothermal(input(&setup), &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let rb = b.step_electrothermal(input(&setup), &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let rr = resistant.step_electrothermal(input(&setup), &force, &gate, |_,_,_| Ok(Some(2.0))).unwrap();
    cold_a.step_electrothermal(input(&off), &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    cold_b.step_electrothermal(input(&off), &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    close(ra.joule_energy_j, rb.joule_energy_j, 1e-12);
    close(ra.joule_energy_j, 2.0 * rr.joule_energy_j, 1e-12);
    let ta = a.lattice.temperature(0,0,3) - cold_a.lattice.temperature(0,0,3);
    let tb = b.lattice.temperature(0,0,3) - cold_b.lattice.temperature(0,0,3);
    assert!(ta > 0.0);
    close(ta, 2.0*tb, 1e-11);
}

#[test]
fn extruded_2d_and_3d_share_electrical_geometry_and_depth_scales_energy_only() {
    let setup2 = electrodes(1, 2.0);
    let setup3 = electrodes(2, 2.0);
    let mut a = slab2(); let mut b = slab2(); let mut c = slab3(1000.0);
    let gate = CancelGate::new();
    let ra = a.step_electrothermal(input(&setup2), &[[0.0;2];16], DX, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let rb = b.step_electrothermal(input(&setup2), &[[0.0;2];16], 2.0*DX, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let rc = c.step_electrothermal(input(&setup3), &[[0.0;3];16], &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    close(ra.joule_energy_j, rc.joule_energy_j, 1e-12);
    close(rb.joule_energy_j, 2.0*ra.joule_energy_j, 1e-12);
    for (pa,pb) in a.lattice.gpop.iter().flatten().zip(b.lattice.gpop.iter().flatten()) { close(*pa,*pb,1e-12); }
}

#[test]
fn contact_heat_partition_is_preserved_without_double_counting() {
    let setup = electrodes(2, 2.0);
    let contacts = [
        ElectricalContact { cells: [6,8], resistance_ohm_m2: 0.02, heat_fraction_to_first: 0.75 },
        ElectricalContact { cells: [7,9], resistance_ohm_m2: 0.02, heat_fraction_to_first: 0.75 },
    ];
    let mut slab = slab3(1000.0);
    let report = slab.step_electrothermal(ElectricHeatInput { contacts: &contacts, ..input(&setup) },
        &[[0.0;3];16], &CancelGate::new(), |_,_,_| Ok(Some(1.0))).unwrap();
    let bulk_r = 4.0 * DX / DX.powi(2);
    let contact_r = 0.02 / DX.powi(2);
    let current = 2.0 / (bulk_r + contact_r);
    close(report.electric.total_joule_w, 2.0 * 2.0 * current, 1e-12);
    close(report.electric.joule_w[6] - report.electric.joule_w[8], 0.5*contact_r*current*current, 1e-12);
    close(report.heat.source_j, report.electric.total_joule_w * DT, 1e-12);
}

#[test]
fn invalid_or_cancelled_trials_leave_all_physical_state_unchanged() {
    let mut slab = slab3(1000.0);
    let before_f = slab.lattice.grid.f.clone(); let before_g = slab.lattice.gpop.clone();
    let before_force = slab.lattice.grid.fext.clone(); let before_tau = slab.lattice.grid.tau.clone();
    let force = [[0.0;3];16];
    let setup = electrodes(2, 2.0);
    let hot = electrodes(2, 1e5);
    assert!(slab.step_electrothermal(input(&hot), &force, &CancelGate::new(), |_,_,_| Ok(Some(1.0))).is_err());
    assert!(slab.step_electrothermal(input(&setup), &force, &CancelGate::new(), |_,_,_| Err(ChtError::InvalidInput {
        field: "test.material", reason: "missing source record".into(),
    })).is_err());
    let gate = CancelGate::new();
    let result = slab.step_electrothermal(input(&setup), &force, &gate, |_,_,_| {
        gate.request(); Ok(Some(1.0))
    });
    assert!(matches!(result, Err(ElectricHeatError::Cancelled)));
    assert_eq!(slab.lattice.grid.f, before_f); assert_eq!(slab.lattice.gpop, before_g);
    assert_eq!(slab.lattice.grid.fext, before_force); assert_eq!(slab.lattice.grid.tau, before_tau);
    let mut wall = setup.clone(); wall.electrodes[0].lo[2] = 0;
    let mut calls = 0;
    assert!(slab.step_electrothermal(input(&wall), &force, &CancelGate::new(), |_,_,_| {
        calls += 1; Ok(Some(1.0))
    }).is_err());
    assert_eq!(calls, 0);
    let bounded = ElectricHeatInput { limits: ElectricHeatLimits {
        thermal: HeatSourceLimits { max_cells: 15, ..limits().thermal }, ..limits()
    }, ..input(&setup) };
    assert!(slab.step_electrothermal(bounded, &force, &CancelGate::new(), |_,_,_| {
        calls += 1; Ok(Some(1.0))
    }).is_err());
    assert_eq!(calls, 0);
}

#[test]
fn next_step_uses_evolved_kelvin_and_the_retained_material_snapshot() {
    let setup = electrodes(2, 100.0);
    let mut slab = slab3(1000.0);
    let force = [[0.0;3];16];
    let gate = CancelGate::new();
    slab.step_electrothermal(input(&setup), &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let expected: Vec<_> = slab.lattice.gpop.iter().map(|p| slab.temperature_si(p.iter().sum()).unwrap()).collect();
    let mut calls = 0;
    slab.step_electrothermal(input(&setup), &force, &gate, |cell, t, material| {
        close(t, expected[cell], 0.0);
        assert_eq!(*material, transport(1000.0));
        calls += 1;
        Ok(Some(1.0 + 0.002*(t-300.0)))
    }).unwrap();
    assert_eq!(calls, 8); // Four bulk rows, two cells each; no wall/electrode calls.
}

use fs_lbm::material::electrothermal::midpoint::{MidpointError, MidpointOptions};
fn midpoint_options() -> MidpointOptions { MidpointOptions {
    max_iterations: 32, temperature_tolerance_k: 1e-9, source_tolerance_j: 1e-8, relaxation: 1.0,
} }

#[test]
fn constant_law_midpoint_is_the_original_step_not_repeated_time_advancement() {
    let setup = electrodes(2, 20.0);
    let mut a = slab3(1000.0); let mut b = slab3(1000.0);
    let force = [[1.0, 0.0, 0.0];16]; let gate = CancelGate::new();
    let direct = a.step_electrothermal(input(&setup), &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let coupled = b.step_electrothermal_midpoint(input(&setup), &force, midpoint_options(), &gate,
        |_,_,_| Ok(Some(1.0))).unwrap();
    assert_eq!(coupled.step.heat, direct.heat);
    assert_eq!(b.lattice.grid.f, a.lattice.grid.f);
    assert_eq!(b.lattice.gpop, a.lattice.gpop);
    assert!(coupled.iterations <= 2);
    assert!(coupled.electrical_solves <= 3);
    assert_eq!(coupled.joule_residual_j, 0.0);
}

#[test]
fn current_driven_temperature_feedback_closes_at_actual_midpoint() {
    let mut setup = electrodes(2, 0.0);
    setup.electrodes[0].drive = ElectrodeDrive::Current(2.0);
    let mut slab = slab3(1000.0); let mut explicit = slab3(1000.0);
    let old: Vec<_> = slab.lattice.gpop.iter().map(|p| slab.temperature_si(p.iter().sum()).unwrap()).collect();
    let force = [[0.0;3];16]; let gate = CancelGate::new();
    let law = |_: usize, t: f64, _: &FluidTransportState| Ok(Some(1.0+0.02*(t-300.0)));
    let lagged = explicit.step_electrothermal(input(&setup), &force, &gate, law).unwrap();
    let report = slab.step_electrothermal_midpoint(input(&setup), &force, midpoint_options(), &gate, law).unwrap();
    assert!(report.step.joule_energy_j > lagged.joule_energy_j);
    assert!(report.iterations > 1);
    let exact_midpoint: Vec<_> = old.iter().zip(&slab.lattice.gpop).map(|(a,p)| {
        0.5*a+0.5*slab.temperature_si(p.iter().sum()).unwrap()
    }).collect();
    for (a,b) in exact_midpoint.iter().zip(&report.evaluation_temperature_k) {
        close(*a,*b,midpoint_options().temperature_tolerance_k);
    }
    let domain = fs_lbm::conjugate::VoxelDomain::new(2,1,8,DX).unwrap();
    let verified = setup.solve_with_cell_resistivity(&domain, &exact_midpoint, |cell,_,t| {
        if !(2..14).contains(&cell) { Ok(None) } else { Ok(Some(1.0+0.02*(t-300.0))) }
    }, &gate).unwrap();
    let discrepancy: f64 = verified.joule_w.iter().zip(&report.step.electric.joule_w).map(|(a,b)| (a-b).abs()*DT).sum();
    close(discrepancy,report.joule_residual_j,1e-12);
    assert!(discrepancy <= midpoint_options().source_tolerance_j);
}

#[test]
fn midpoint_budget_and_mid_iteration_cancellation_do_not_publish_trials() {
    let setup = electrodes(2, 20.0);
    let mut slab = slab3(1000.0); let force = [[0.0;3];16];
    let old_f = slab.lattice.grid.f.clone(); let old_g = slab.lattice.gpop.clone();
    let too_short = MidpointOptions { max_iterations: 1, ..midpoint_options() };
    let result = slab.step_electrothermal_midpoint(input(&setup), &force, too_short, &CancelGate::new(), |_,_,_| Ok(Some(1.0)));
    assert!(matches!(result, Err(MidpointError::NonConvergence { iterations: 1, .. })));
    assert_eq!(slab.lattice.grid.f, old_f); assert_eq!(slab.lattice.gpop, old_g);
    // A loose temperature criterion must not bypass the independent source recheck.
    let loose_temperature = MidpointOptions { max_iterations: 1, temperature_tolerance_k: 1e3,
        source_tolerance_j: 0.0, ..midpoint_options() };
    let result = slab.step_electrothermal_midpoint(input(&setup), &force, loose_temperature,
        &CancelGate::new(), |_,t,_| Ok(Some(1.0+0.02*(t-300.0))));
    assert!(matches!(result, Err(MidpointError::NonConvergence {
        joule_residual_j: Some(discrepancy), ..
    }) if discrepancy > 0.0));
    assert_eq!(slab.lattice.grid.f, old_f); assert_eq!(slab.lattice.gpop, old_g);
    let gate = CancelGate::new(); let mut calls = 0;
    let result = slab.step_electrothermal_midpoint(input(&setup), &force, midpoint_options(), &gate, |_,_,_| {
        calls += 1; if calls == 9 { gate.request(); } Ok(Some(1.0))
    });
    assert!(matches!(result, Err(MidpointError::Step(ElectricHeatError::Cancelled))));
    assert_eq!(slab.lattice.grid.f, old_f); assert_eq!(slab.lattice.gpop, old_g);
    let invalid = MidpointOptions { relaxation: 0.0, ..midpoint_options() };
    assert!(matches!(slab.step_electrothermal_midpoint(input(&setup), &force, invalid, &CancelGate::new(), |_,_,_| {
        panic!("invalid controls must refuse before querying material")
    }), Err(MidpointError::InvalidOptions)));
}

#[test]
fn two_dimensional_midpoint_drive_switch_off_keeps_thermal_history() {
    let on = electrodes(1, 20.0); let off = electrodes(1, 0.0);
    let mut slab = slab2(); let initial = slab.lattice.gpop.clone(); let gate = CancelGate::new();
    let force = [[0.0;2];16];
    for _ in 0..5 {
        slab.step_electrothermal_midpoint(input(&on), &force, DX, midpoint_options(), &gate,
            |_,t,_| Ok(Some(1.0+0.01*(t-300.0)))).unwrap();
    }
    let result = slab.step_electrothermal_midpoint(input(&off), &force, DX, midpoint_options(), &gate,
        |_,t,_| Ok(Some(1.0+0.01*(t-300.0)))).unwrap();
    close(result.step.heat.source_j, 0.0, 0.0);
    assert_ne!(slab.lattice.gpop, initial);
    assert!(result.step.combined_defect_j.abs() < 1e-9);
}
