//! G0/G1/G3/G4: actual current/contact/local-transport composition; synthetic data.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::material::electrothermal::magnetic::midpoint::MagneticMidpointOptions;
use fs_lbm::material::electrothermal::magnetic::transport::ElectricalTransport;
use fs_lbm::material::electrothermal::magnetic::transport::midpoint::MaterialMidpointOptions;
use fs_lbm::material::electrothermal::magnetic::transport::midpoint::interfaces::*;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::{MaterialError, fluid_transport::{FluidTransport, FluidTransportState}};

const DX: f64 = 0.01;
fn material() -> FluidTransportState { FluidTransportState::try_new(1000.0, 1000.0, 100.0, 1.0).unwrap() }
fn bulk(_: usize, _: f64, _: &FluidTransportState) -> Result<ElectricalTransport<FluidTransportState>, MaterialError> {
    ElectricalTransport::try_new(material(), Some(1.0))
}
fn slab3() -> MaterialSlab3 {
    plan_scaling_from_material(&material(), DX, DX, 0.05).unwrap().thermal_slab3(MaterialThermalSlab3 {
        nx: 2, ny: 1, nz: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&material(), DX, DX, 0.05).unwrap().thermal_slab(MaterialThermalSlab {
        nx: 2, ny: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
        gravity_m_s2: 0.0, expansion_per_k: 0.0, max_cells: 16,
    }).unwrap()
}
fn setup(dim: usize, current: f64) -> ElectricSetup {
    let (lo, hi, upper_lo, upper_hi) = if dim == 3 {
        ([0,0,1], [2,1,2], [0,0,6], [2,1,7])
    } else { ([0,1,0], [2,2,1], [0,6,0], [2,7,1]) };
    ElectricSetup { conductors: vec![], electrodes: vec![
        Electrode { lo, hi, drive: ElectrodeDrive::Current(current) },
        Electrode { lo: upper_lo, hi: upper_hi, drive: ElectrodeDrive::Voltage(0.0) },
    ] }
}
fn contacts() -> [ElectricalContact; 2] { [
    ElectricalContact { cells: [6,8], resistance_ohm_m2: 0.02, heat_fraction_to_first: 0.75 },
    ElectricalContact { cells: [7,9], resistance_ohm_m2: 0.02, heat_fraction_to_first: 0.75 },
] }
fn input<'a>(s: &'a ElectricSetup, c: &'a [ElectricalContact]) -> ElectricHeatInput<'a> {
    ElectricHeatInput { setup: s, contacts: c, limits: ElectricHeatLimits {
        thermal: HeatSourceLimits { max_cells: 16, temperature_range_k: [290.0,330.0], balance_tolerance_j: 1e-8 },
        coupling_tolerance_j: 1e-9,
    } }
}
fn options() -> InterfaceMidpointOptions { InterfaceMidpointOptions {
    material: MaterialMidpointOptions { coupling: MagneticMidpointOptions {
        max_iterations: 128, temperature_tolerance_k: 1e-9, velocity_tolerance_m_s: 1e-10,
        joule_tolerance_j: 1e-9, impulse_tolerance_n_s: 1e-11, work_tolerance_j: 1e-9, relaxation: 0.5,
    }, transport_relative_tolerance: 1e-10 },
    resistance_relative_tolerance: 1e-10, partition_absolute_tolerance: 1e-10, circuit: InterfaceCircuit::Wired,
} }
fn constant(c: &ElectricalContact) -> Result<InterfaceResponse<()>, MaterialError> {
    Ok(InterfaceResponse { material: (), resistance_ohm_m2: c.resistance_ohm_m2, heat_fraction_to_first: c.heat_fraction_to_first })
}
fn near(a: f64, b: f64, tol: f64) { assert!((a-b).abs() <= tol, "{a} versus {b}"); }

#[test]
fn constant_interfaces_reproduce_the_existing_material_midpoint_populations() {
    let s = setup(3, 0.1); let c = contacts(); let gate = CancelGate::new();
    let field = [[0.0,2.0,0.0];16]; let forces = [[1.0,0.0,0.0];16];
    let mut a = slab3(); let mut b = slab3();
    let old = a.step_magnetothermal_midpoint_with_material_transport(input(&s,&c), &field, &forces,
        options().material, &gate, bulk).unwrap();
    let new = b.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &field, &forces,
        options(), &gate, bulk, |_,c,_| constant(c)).unwrap();
    assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
    assert_eq!(a.lattice.gpop, b.lattice.gpop);
    assert_eq!(old.step.heat, new.bulk.step.heat);
    assert_eq!(old.step.iterations, new.bulk.step.iterations);
    assert_eq!(old.step.motion.electric.joule_w, new.bulk.step.motion.electric.joule_w);
    assert_eq!(new.interface_residual.resistance_relative, 0.0);
    assert_eq!(new.interface_residual.partition_absolute, 0.0);
}

#[test]
fn contact_temperature_changes_current_power_and_local_heat_without_double_deposition() {
    let s = setup(3, 1.0); let c = contacts(); let mut a = slab3(); let mut fixed = slab3();
    let gate = CancelGate::new(); let zero = [[0.0;3];16];
    let old = fixed.step_magnetothermal_midpoint_with_material_transport(input(&s,&c), &zero, &zero,
        options().material, &gate, bulk).unwrap();
    let r = a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &zero, &zero, options(), &gate, bulk,
        |i, site, pair| Ok::<_, MaterialError>(InterfaceResponse {
            material: (i, [pair[0].temperature_k, pair[1].temperature_k]),
            resistance_ohm_m2: site.resistance_ohm_m2*(1.0+0.05*(0.5*pair[0].temperature_k+0.5*pair[1].temperature_k-300.0)),
            heat_fraction_to_first: 0.5+0.002*(pair[0].temperature_k-pair[1].temperature_k),
        })).unwrap();
    assert!(r.bulk.step.joule_energy_j > old.step.joule_energy_j);
    assert!(r.bulk.step.motion.electric.electrodes[0].voltage_v > old.step.motion.electric.electrodes[0].voltage_v);
    for sample in &r.applied_interfaces {
        assert_eq!(sample.material, (sample.index, sample.temperature_k));
        for a in 0..2 {
            near(sample.temperature_k[a], r.bulk.step.evaluation_temperature_k[sample.contact.cells[a]], 0.0);
        }
        let contact_w = 0.25 * sample.contact.resistance_ohm_m2 / DX.powi(2); // I_path = 0.5 A
        let [a,b] = sample.contact.cells;
        near(r.bulk.step.motion.electric.joule_w[a], 25.0 + sample.contact.heat_fraction_to_first*contact_w, 1e-8);
        near(r.bulk.step.motion.electric.joule_w[b], 25.0 + (1.0-sample.contact.heat_fraction_to_first)*contact_w, 1e-8);
    }
    near(r.bulk.step.heat.source_j, r.bulk.step.motion.electric.total_joule_w*DX, 1e-9);
    assert!(r.bulk.step.coupling_defect_j.abs() < 1e-8);
    assert!(r.interface_residual.resistance_relative <= options().resistance_relative_tolerance);
    assert!(r.interface_residual.partition_absolute <= options().partition_absolute_tolerance);
}

#[test]
fn zero_current_does_not_bypass_interface_constitutive_convergence() {
    let s = setup(3, 0.0); let c = contacts(); let mut a = slab3(); a.lattice.gpop[6].iter_mut().for_each(|v| *v *= 1.1);
    let before_f = a.lattice.grid.f.clone(); let before_g = a.lattice.gpop.clone();
    let mut opt = options(); opt.material.coupling.max_iterations = 1;
    opt.material.coupling.temperature_tolerance_k = 1e3; opt.material.coupling.velocity_tolerance_m_s = 1e3;
    opt.resistance_relative_tolerance = 0.0;
    let result = a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &[[0.0;3];16], &[[0.0;3];16],
        opt, &CancelGate::new(), bulk, |_,site,pair| Ok::<_,MaterialError>(InterfaceResponse {
            material: (), resistance_ohm_m2: site.resistance_ohm_m2*(1.0+0.1*(pair[0].temperature_k-300.0)),
            heat_fraction_to_first: 0.5,
        }));
    match result {
        Err(InterfaceMidpointError::NonConvergence { last_interface: Some(r), .. }) => assert!(r.resistance_relative > 0.0),
        _ => panic!("zero-current electrical checks must not accept unconverged interfaces"),
    }
    assert_eq!(a.lattice.grid.f, before_f); assert_eq!(a.lattice.gpop, before_g);
}

#[derive(Debug, PartialEq)]
struct Refusal(&'static str);
impl core::fmt::Display for Refusal { fn fmt(&self, f: &mut core::fmt::Formatter<'_>) -> core::fmt::Result { write!(f, "{}", self.0) } }
impl std::error::Error for Refusal {}

#[test]
fn endpoint_typed_failure_and_cancellation_leave_accepted_history_unchanged() {
    let s = setup(3,0.0); let c = contacts(); let mut a = slab3();
    let before_f = a.lattice.grid.f.clone(); let before_g = a.lattice.gpop.clone(); let before_force = a.lattice.grid.fext.clone();
    // Quiescent linear conduction: applied pair, rechecked pair, endpoint pair.
    let mut calls = 0;
    let result = a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &[[0.0;3];16], &[[0.0;3];16],
        options(), &CancelGate::new(), bulk, |_,site,_| {
            calls += 1;
            if calls == 5 { return Err(Refusal("endpoint outside interface source support")); }
            Ok(InterfaceResponse { material: (), resistance_ohm_m2: site.resistance_ohm_m2, heat_fraction_to_first: 0.75 })
        });
    match result {
        Err(InterfaceMidpointError::Interface { index:0, cells:[6,8], source, .. }) => assert_eq!(source, Refusal("endpoint outside interface source support")),
        _ => panic!("original typed interface refusal must survive"),
    }
    let gate = CancelGate::new(); let mut calls = 0;
    assert!(a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &[[0.0;3];16], &[[0.0;3];16],
        options(), &gate, bulk, |_,site,_| {
            calls += 1; if calls == 5 { gate.request(); } constant(site)
        }).is_err());
    assert!(gate.is_requested());
    assert_eq!(a.lattice.grid.f, before_f); assert_eq!(a.lattice.gpop, before_g); assert_eq!(a.lattice.grid.fext,before_force);
}

#[test]
fn malformed_faces_refuse_before_material_callbacks_and_bad_laws_do_not_publish() {
    let s = setup(3,0.0); let good = contacts(); let mut a = slab3(); let before = a.lattice.gpop.clone();
    for c in [vec![ElectricalContact { cells:[0,2], ..good[0] }], vec![good[0],good[0]],
        vec![ElectricalContact { cells:[4,10], ..good[0] }]] {
        let mut calls = 0;
        assert!(a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &[[0.0;3];16], &[[0.0;3];16],
            options(), &CancelGate::new(), |i,t,r| {calls += 1; bulk(i,t,r)}, |_,c,_| constant(c)).is_err());
        assert_eq!(calls,0);
    }
    for (r,p) in [(-1.0,0.5),(f64::NAN,0.5),(0.1,1.01),(0.1,f64::INFINITY)] {
        assert!(a.step_magnetothermal_midpoint_with_interfaces(input(&s,&good), &[[0.0;3];16], &[[0.0;3];16],
            options(), &CancelGate::new(), bulk, |_,_,_| Ok::<_,MaterialError>(InterfaceResponse {
                material: (), resistance_ohm_m2:r, heat_fraction_to_first:p,
            })).is_err());
    }
    assert_eq!(a.lattice.gpop,before);
}

#[test]
fn ordered_surface_reversal_and_2d_extrusion_preserve_the_same_response() {
    let mut s = setup(2,0.0); s.electrodes[0].drive = ElectrodeDrive::Voltage(2.0);
    let c = contacts(); let reversed = c.map(|v| ElectricalContact { cells:[v.cells[1],v.cells[0]],
        heat_fraction_to_first:1.0-v.heat_fraction_to_first, ..v });
    let mut a = slab2(); let mut b = slab2(); let mut reverse = slab2(); let gate = CancelGate::new();
    let r = a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &[[0.0,0.0,2.0];16], &[[0.0;2];16],
        DX, options(), &gate, bulk, |_,site,_| constant(site)).unwrap();
    let q = b.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &[[0.0,0.0,2.0];16], &[[0.0;2];16],
        2.0*DX, options(), &gate, bulk, |_,site,_| constant(site)).unwrap();
    reverse.step_magnetothermal_midpoint_with_interfaces(input(&s,&reversed), &[[0.0,0.0,2.0];16], &[[0.0;2];16],
        DX, options(), &gate, bulk, |_,site,_| constant(site)).unwrap();
    near(q.bulk.step.joule_energy_j, 2.0*r.bulk.step.joule_energy_j, 1e-9);
    for (x,y) in a.lattice.grid.f.iter().flatten().zip(b.lattice.grid.f.iter().flatten()) { near(*x,*y,1e-11); }
    assert_eq!(a.lattice.grid.f, reverse.lattice.grid.f); assert_eq!(a.lattice.gpop,reverse.lattice.gpop);
}

#[test]
fn isolated_conductor_gauges_and_contact_heat_use_the_existing_owner() {
    let s = ElectricSetup::default(); let c = contacts(); let gate = CancelGate::new();
    let mut a = slab3(); let mut b = slab3();
    for cell in 2..14 {
        a.lattice.grid.f[cell] = fs_lbm::equilibrium3(1.0,[0.02,0.0,0.0]);
        b.lattice.grid.f[cell] = a.lattice.grid.f[cell];
    }
    let field: Vec<_> = (0..16).map(|i| [0.0,if i%2==0 {20.0} else {60.0},0.0]).collect();
    let mut opt = options(); opt.circuit = InterfaceCircuit::Isolated { current_tolerance_a:1e-9 };
    let r = a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &field, &[[0.0;3];16],
        opt, &gate, bulk, |_,site,_| constant(site)).unwrap();
    let q = b.step_magnetothermal_midpoint_isolated_with_material_transport(input(&s,&c), &field, &[[0.0;3];16],
        1e-9, opt.material, &gate, bulk).unwrap();
    assert_eq!(a.lattice.grid.f,b.lattice.grid.f); assert_eq!(a.lattice.gpop,b.lattice.gpop);
    assert_eq!(r.bulk.step.heat,q.step.heat); assert_eq!(r.bulk.step.electrode_energy_j,0.0);
    assert!(r.bulk.step.joule_energy_j > 0.0);
    near(r.bulk.step.joule_energy_j,-r.bulk.step.midpoint_mechanical_work_j,1e-9);
}

#[test]
fn nonclone_source_objects_are_sampled_once_and_retained_for_both_surfaces() {
    struct Once(std::cell::Cell<bool>);
    impl FluidTransport for Once {
        fn transport_state(&self) -> Result<FluidTransportState,MaterialError> {
            assert!(!self.0.replace(true),"must not re-query an evaluated bulk object"); Ok(material())
        }
    }
    struct Receipt { index:usize, temperatures:[f64;2] }
    let s = setup(3,0.1); let c = contacts(); let mut a = slab3(); let mut bulk_calls = 0; let mut interface_calls = 0;
    let r = a.step_magnetothermal_midpoint_with_interfaces(input(&s,&c), &[[0.0;3];16], &[[0.0;3];16],
        options(), &CancelGate::new(), |_,_,_| {
            bulk_calls += 1; ElectricalTransport::try_new(Once(std::cell::Cell::new(false)),Some(1.0))
        }, |index,site,pair| {
            interface_calls += 1;
            assert!(pair.iter().all(|p|p.material.material().0.get()));
            Ok::<_,MaterialError>(InterfaceResponse { material:Receipt {index,temperatures:pair.map(|p|p.temperature_k)},
                resistance_ohm_m2:site.resistance_ohm_m2, heat_fraction_to_first:site.heat_fraction_to_first })
        }).unwrap();
    assert_eq!(bulk_calls,12*(r.bulk.step.electrical_solves+1));
    assert_eq!(interface_calls,2*(r.bulk.step.electrical_solves+1));
    for sample in r.applied_interfaces.iter().chain(&r.endpoint_interfaces) {
        assert_eq!(sample.material.index,sample.index); assert_eq!(sample.material.temperatures,sample.temperature_k);
    }
}
