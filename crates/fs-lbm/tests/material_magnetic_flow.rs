//! Real electrical/thermal/LBM composition with explicitly synthetic materials.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ChtError, ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatError, ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::fluid_transport::FluidTransportState;

const DX: f64 = 0.01;
const DT: f64 = 0.01;
fn material(cp: f64) -> FluidTransportState { FluidTransportState::try_new(1000.0, cp, 100.0, 1.0).unwrap() }
fn slab3(cp: f64) -> MaterialSlab3 {
    plan_scaling_from_material(&material(cp), DX, DT, 0.1).unwrap().thermal_slab3(MaterialThermalSlab3 {
        nx: 2, ny: 1, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&material(1000.0), DX, DT, 0.1).unwrap().thermal_slab(MaterialThermalSlab {
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
fn near(a: f64, b: f64) { assert!((a-b).abs() < 1e-10, "{a} versus {b}"); }

#[test]
fn g1_lorentz_load_delivers_one_physical_impulse_in_two_and_three_dimensions() {
    // The four bulk cells between the ideal electrodes have total conducting
    // length 4*DX. Thus J=V/(rho_e*4*DX), F=J cross B, and at an initially
    // uniform interior rho_m*delta_u=F*DT. No wall or contact crosses the
    // sampled cells' one-step lattice domain of dependence.
    let expected_speed = 2.0 / (4.0 * DX) * 2.0 * DT / material(1000.0).density_kg_m3();
    let mut two = slab2();
    let mut three = slab3(1000.0);
    let setup2 = electrodes(1, 2.0);
    let setup3 = electrodes(2, 2.0);
    let gate = CancelGate::new();
    two.step_magnetothermal(input(&setup2), &[[0.0, 0.0, 2.0]; 16], &[[0.0; 2]; 16],
        DX, &gate, |_, _, _| Ok(Some(1.0))).unwrap();
    three.step_magnetothermal(input(&setup3), &[[0.0, 2.0, 0.0]; 16], &[[0.0; 3]; 16],
        &gate, |_, _, _| Ok(Some(1.0))).unwrap();
    near(two.lattice.grid.moments(8).u[0] * DX / DT, expected_speed);
    near(three.lattice.grid.moments3(8).1[0] * DX / DT, -expected_speed);
}

#[test]
fn magnetic_reaction_and_contact_heat_reach_the_actual_flow_and_thermal_kernels() {
    let setup = electrodes(2, 2.0);
    let contacts = [
        ElectricalContact { cells: [6,8], resistance_ohm_m2: 0.02, heat_fraction_to_first: 0.75 },
        ElectricalContact { cells: [7,9], resistance_ohm_m2: 0.02, heat_fraction_to_first: 0.75 },
    ];
    let mut a = slab3(1000.0); let mut b = slab3(1000.0);
    let force = [[1.0,0.0,0.0];16]; let field = [[0.0,2.0,0.0];16];
    let r = a.step_magnetothermal(ElectricHeatInput { contacts: &contacts, ..input(&setup) },
        &field, &force, &CancelGate::new(), |_,_,_| Ok(Some(1.0))).unwrap();
    let power: Vec<_> = r.motion.electric.joule_w.iter().map(|w| w / DX.powi(3)).collect();
    let combined: Vec<[f64;3]> = r.motion.lorentz_force_density_n_m3.iter()
        .zip(force).map(|(m, f)| [m[0]+f[0], m[1]+f[1], m[2]+f[2]]).collect();
    let direct = b.step_with_power_density(&power, &combined, limits().thermal).unwrap();
    assert_eq!(r.heat, direct);
    assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
    assert_eq!(a.lattice.gpop, b.lattice.gpop);
    assert!(a.lattice.grid.moments3(8).1[0] < 0.0);
    assert!(r.motion.electric.joule_w[6] > r.motion.electric.joule_w[8]);
    assert!(r.coupling_defect_j.abs() < 1e-9);
    near(r.lagged_mechanical_work_j, 0.0); // starts at rest, not a kinetic-energy claim
}

#[test]
fn short_circuited_moving_fluid_brakes_and_heats_without_electrode_work() {
    let setup = electrodes(2, 0.0); let field = [[0.0,20.0,0.0];16];
    let mut a = slab3(1000.0); let mut control = slab3(1000.0);
    for c in 2..14 {
        a.lattice.grid.f[c] = fs_lbm::equilibrium3(1.0, [0.02,0.0,0.0]);
        control.lattice.grid.f[c] = a.lattice.grid.f[c];
    }
    let r = a.step_magnetothermal(input(&setup), &field, &[[0.0;3];16],
        &CancelGate::new(), |_,_,_| Ok(Some(1.0))).unwrap();
    control.step_electrothermal(input(&setup), &[[0.0;3];16],
        &CancelGate::new(), |_,_,_| Ok(Some(1.0))).unwrap();
    assert!(r.motion.mechanical_power_w < 0.0);
    assert!(r.joule_energy_j > 0.0);
    near(r.electrode_energy_j, 0.0);
    near(r.joule_energy_j, -r.lagged_mechanical_work_j);
    assert!(a.lattice.grid.moments3(8).1[0] < control.lattice.grid.moments3(8).1[0]);
    assert!(a.lattice.temperature(0,0,3) > control.lattice.temperature(0,0,3));
}

#[test]
fn field_reversal_and_material_substitution_change_force_without_code_switches() {
    let setup = electrodes(2, 2.0); let mut a = slab3(1000.0);
    let mut reverse = slab3(1000.0); let mut resistive = slab3(1000.0);
    let mut capacity = slab3(2000.0); let gate = CancelGate::new();
    let r = a.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];16], &[[0.0;3];16], &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let v = reverse.step_magnetothermal(input(&setup), &[[0.0,-2.0,0.0];16], &[[0.0;3];16], &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let q = resistive.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];16], &[[0.0;3];16], &gate, |_,_,_| Ok(Some(2.0))).unwrap();
    let h = capacity.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];16], &[[0.0;3];16], &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    near(r.joule_energy_j, v.joule_energy_j);
    near(r.joule_energy_j, 2.0*q.joule_energy_j);
    near(r.joule_energy_j, h.joule_energy_j);
    for c in 4..12 {
        near(r.motion.lorentz_force_density_n_m3[c][0], -v.motion.lorentz_force_density_n_m3[c][0]);
        near(r.motion.lorentz_force_density_n_m3[c][0], 2.0*q.motion.lorentz_force_density_n_m3[c][0]);
    }
    near(a.lattice.grid.moments3(8).1[0], -reverse.lattice.grid.moments3(8).1[0]);
    assert!(capacity.lattice.temperature(0,0,3) < a.lattice.temperature(0,0,3));
}

#[test]
fn zero_field_preserves_old_step_and_subsequent_current_reads_live_velocity() {
    let setup = electrodes(2, 20.0); let gate = CancelGate::new();
    let mut a = slab3(1000.0); let mut b = slab3(1000.0);
    let force = [[1.0,0.0,0.0];16];
    let old = a.step_electrothermal(input(&setup), &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let new = b.step_magnetothermal(input(&setup), &[[0.0;3];16], &force, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    assert_eq!(old.heat, new.heat);
    assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
    assert_eq!(a.lattice.gpop, b.lattice.gpop);
    let mut moving = slab3(1000.0);
    let first = moving.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];16], &[[0.0;3];16], &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let second = moving.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];16], &[[0.0;3];16], &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    assert!(second.motion.electric.electrodes[0].current_a < first.motion.electric.electrodes[0].current_a);
    assert!(second.motion.mechanical_power_w > 0.0); // current drives motion
}

#[test]
fn two_dimensional_depth_scales_energy_and_out_of_plane_forces_are_not_discarded() {
    let setup = electrodes(1, 2.0); let gate = CancelGate::new();
    let mut a = slab2(); let mut b = slab2();
    let r = a.step_magnetothermal(input(&setup), &[[0.0,0.0,2.0];16], &[[0.0;2];16], DX, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    let s = b.step_magnetothermal(input(&setup), &[[0.0,0.0,2.0];16], &[[0.0;2];16], 2.0*DX, &gate, |_,_,_| Ok(Some(1.0))).unwrap();
    near(s.joule_energy_j, 2.0*r.joule_energy_j);
    assert!(a.lattice.grid.moments(8).u[0] > 0.0); // Jy x Bz = +Fx
    for (x,y) in a.lattice.grid.f.iter().flatten().zip(b.lattice.grid.f.iter().flatten()) { near(*x,*y); }
    let before = a.lattice.grid.f.clone();
    assert!(a.step_magnetothermal(input(&setup), &[[1.0,0.0,0.0];16], &[[0.0;2];16], DX, &gate,
        |_,_,_| panic!("2D unsupported field must refuse before source query")).is_err());
    assert_eq!(a.lattice.grid.f, before);
}

#[test]
fn source_failure_cancellation_and_late_mach_refusal_preserve_the_accepted_state() {
    let setup = electrodes(2, 2.0); let mut slab = slab3(1000.0);
    let before_f=slab.lattice.grid.f.clone(); let before_g=slab.lattice.gpop.clone();
    let before_force=slab.lattice.grid.fext.clone(); let before_tau=slab.lattice.grid.tau.clone();
    let gate=CancelGate::new();
    assert!(slab.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];15], &[[0.0;3];16], &gate,
        |_,_,_| panic!("field length precedes material queries")).is_err());
    assert!(slab.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];16], &[[0.0;3];16], &gate,
        |_,_,_| Err(ChtError::InvalidInput {field:"source",reason:"unsupported data".into()})).is_err());
    assert!(slab.step_magnetothermal(input(&setup), &[[0.0,1000.0,0.0];16], &[[0.0;3];16], &gate,
        |_,_,_| Ok(Some(1.0))).is_err());
    let result = slab.step_magnetothermal(input(&setup), &[[0.0,2.0,0.0];16], &[[0.0;3];16], &gate,
        |_,_,_| {gate.request(); Ok(Some(1.0))});
    assert!(matches!(result, Err(ElectricHeatError::Cancelled)));
    assert_eq!(slab.lattice.grid.f,before_f); assert_eq!(slab.lattice.gpop,before_g);
    assert_eq!(slab.lattice.grid.fext,before_force); assert_eq!(slab.lattice.grid.tau,before_tau);
}
