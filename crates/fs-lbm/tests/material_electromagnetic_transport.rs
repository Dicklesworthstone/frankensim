//! Real current/heat/momentum composition, with synthetic material data only.
use fs_exec::CancelGate;
use fs_lbm::conjugate::{ElectricSetup, ElectricalContact, Electrode, ElectrodeDrive};
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::electrothermal::{ElectricHeatInput, ElectricHeatLimits};
use fs_lbm::material::electrothermal::magnetic::transport::{ElectricalTransport, MaterialMagneticError};
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::{MaterialError, fluid_transport::{FluidTransport, FluidTransportState}};
use std::{cell::Cell, rc::Rc};

const N: usize = 16;
const DX: f64 = 0.01;
fn fluid(t: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0, 1000.0, 50.0 + 8.0*(t-290.0), 0.5+0.03*(t-290.0)).unwrap()
}
fn initial() -> FluidTransportState { FluidTransportState::try_new(1000.0,1000.0,100.0,1.0).unwrap() }
fn law(t: f64) -> Result<ElectricalTransport<FluidTransportState>, MaterialError> {
    ElectricalTransport::try_new(fluid(t), Some(1.0+0.01*(t-300.0)))
}
fn slab3() -> MaterialSlab3 {
    plan_scaling_from_material(&initial(), DX, 0.01, 0.1).unwrap().thermal_slab3(MaterialThermalSlab3 {
        nx:2,ny:1,nz:6,bottom_temperature_k:310.0,top_temperature_k:300.0,
        gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:N,
    }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&initial(), DX, 0.01, 0.1).unwrap().thermal_slab(MaterialThermalSlab {
        nx:2,ny:6,bottom_temperature_k:310.0,top_temperature_k:300.0,
        gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:N,
    }).unwrap()
}
fn setup(normal: usize, voltage: f64) -> ElectricSetup {
    let (lo,hi,lo2,hi2)=if normal==2 {([0,0,1],[2,1,2],[0,0,6],[2,1,7])}
        else {([0,1,0],[2,2,1],[0,6,0],[2,7,1])};
    ElectricSetup {conductors:vec![],electrodes:vec![
        Electrode {lo,hi,drive:ElectrodeDrive::Voltage(voltage)},
        Electrode {lo:lo2,hi:hi2,drive:ElectrodeDrive::Voltage(0.0)},
    ]}
}
fn input(setup: &ElectricSetup) -> ElectricHeatInput<'_> {
    ElectricHeatInput {setup,contacts:&[],limits:ElectricHeatLimits {
        thermal:HeatSourceLimits {max_cells:N,temperature_range_k:[290.0,330.0],balance_tolerance_j:1e-9},
        coupling_tolerance_j:1e-10,
    }}
}
fn near(a:f64,b:f64) {assert!((a-b).abs()<1e-10,"{a} != {b}");}

#[test]
fn constant_material_recovers_existing_wired_dynamics_including_zero_field() {
    for b in [0.0,2.0] {
        let electric=setup(2,20.0); let gate=CancelGate::new();
        let mut a=slab3(); let mut old=slab3();
        let force=[[1.0,0.0,0.0];N]; let field=[[0.0,b,0.0];N];
        let r=a.step_magnetothermal_with_material_transport(input(&electric),&field,&force,&gate,
            |_,_,base| ElectricalTransport::try_new(*base,Some(1.0))).unwrap();
        let expected=old.step_magnetothermal(input(&electric),&field,&force,&gate,|_,_,_| Ok(Some(1.0))).unwrap();
        assert_eq!(r.step.heat,expected.heat);
        assert_eq!(r.step.motion.electric.joule_w,expected.motion.electric.joule_w);
        assert_eq!(a.lattice.grid.f,old.lattice.grid.f);
        assert_eq!(a.lattice.gpop,old.lattice.gpop);
        assert_eq!(a.lattice.grid.tau,old.lattice.grid.tau);
        assert_eq!(r.applied.len(),12); assert_eq!(r.endpoint.len(),12);
    }
}

#[test]
fn same_cell_material_drives_actual_current_contact_heat_and_local_transport() {
    let electric=setup(2,20.0); let gate=CancelGate::new();
    let contacts=[ElectricalContact {cells:[6,8],resistance_ohm_m2:0.02,heat_fraction_to_first:0.75},
        ElectricalContact {cells:[7,9],resistance_ohm_m2:0.02,heat_fraction_to_first:0.75}];
    let input=ElectricHeatInput {contacts:&contacts,..input(&electric)};
    let force=[[1.0,0.0,0.0];N]; let field=[[0.0,2.0,0.0];N];
    let mut coupled=slab3(); let mut source_owner=slab3(); let mut reference=slab3();
    let source=source_owner.step_magnetothermal(input,&field,&force,&gate,
        |_,t,_| Ok(Some(1.0+0.01*(t-300.0)))).unwrap();
    let power:Vec<_>=source.motion.electric.joule_w.iter().map(|w| w/DX.powi(3)).collect();
    let forces:Vec<_>=source.motion.lorentz_force_density_n_m3.iter()
        .zip(force).map(|(a,b)| [a[0]+b[0],a[1]+b[1],a[2]+b[2]]).collect();
    let separate=reference.step_with_material_transport(&power,&forces,input.limits.thermal,&gate,|_,t,_| Ok(fluid(t))).unwrap();
    let r=coupled.step_magnetothermal_with_material_transport(input,&field,&force,&gate,|_,t,_| law(t)).unwrap();
    assert_eq!(r.step.heat,separate.heat);
    assert_eq!(r.step.motion.electric.joule_w,source.motion.electric.joule_w);
    assert_eq!(coupled.lattice.grid.f,reference.lattice.grid.f);
    assert_eq!(coupled.lattice.gpop,reference.lattice.gpop);
    assert_ne!(coupled.lattice.gpop,source_owner.lattice.gpop);
    assert!(r.applied.windows(2).any(|p| p[0].thermal_relaxation!=p[1].thermal_relaxation));
    near(r.step.joule_energy_j,r.step.heat.source_j);
}

struct Receipt { id:usize, state:FluidTransportState, reads:Rc<Cell<usize>> }
impl FluidTransport for Receipt {
    fn transport_state(&self)->Result<FluidTransportState,MaterialError> {
        self.reads.set(self.reads.get()+1); Ok(self.state)
    }
}
#[test]
fn original_source_payload_is_retained_and_queried_once_per_cell_epoch() {
    let electric=setup(2,20.0); let mut a=slab3(); let mut calls=0;
    let reads=Rc::new(Cell::new(0));
    let r=a.step_magnetothermal_with_material_transport(input(&electric),&[[0.0,2.0,0.0];N],&[[0.0;3];N],&CancelGate::new(),
        |cell,t,base| {
            assert_eq!(cell,2+calls%12); assert_eq!(*base,initial());
            let id=calls; calls+=1;
            ElectricalTransport::try_new(Receipt {id,state:fluid(t),reads:reads.clone()},Some(1.0))
        }).unwrap();
    assert_eq!(calls,24); assert_eq!(reads.get(),24);
    for (i,s) in r.applied.iter().enumerate() {
        assert_eq!(s.material.material().id,i);
        assert_eq!(s.transport,fluid(s.temperature_k));
        assert_eq!(r.endpoint[i].material.material().id,12+i);
        let t=a.temperature_si(a.lattice.gpop[s.cell].iter().sum()).unwrap();
        near(r.endpoint[i].temperature_k,t);
    }
}

#[test]
fn electrode_free_materials_keep_loop_heating_and_force_without_supply_power() {
    let electric=ElectricSetup::default(); let gate=CancelGate::new();
    let field:Vec<_>=(0..N).map(|c| [0.0,if c%2==0 {20.0} else {40.0},0.0]).collect();
    let mut coupled=slab3(); let mut source_owner=slab3(); let mut reference=slab3();
    for c in 2..14 {
        let f=fs_lbm::equilibrium3(1.0,[0.02,0.0,0.0]);
        coupled.lattice.grid.f[c]=f; source_owner.lattice.grid.f[c]=f; reference.lattice.grid.f[c]=f;
    }
    let source=source_owner.step_magnetothermal_isolated(input(&electric),&field,&[[0.0;3];N],1e-10,&gate,
        |_,t,_| Ok(Some(1.0+0.01*(t-300.0)))).unwrap();
    let power:Vec<_>=source.motion.electric.joule_w.iter().map(|w| w/DX.powi(3)).collect();
    reference.step_with_material_transport(&power,&source.motion.lorentz_force_density_n_m3,
        input(&electric).limits.thermal,&gate,|_,t,_| Ok(fluid(t))).unwrap();
    let r=coupled.step_magnetothermal_isolated_with_material_transport(input(&electric),&field,&[[0.0;3];N],1e-10,&gate,
        |_,t,_| law(t)).unwrap();
    assert_eq!(coupled.lattice.grid.f,reference.lattice.grid.f);
    assert_eq!(coupled.lattice.gpop,reference.lattice.gpop);
    assert_eq!(r.step.electrode_energy_j,0.0);
    assert!(r.step.joule_energy_j>0.0); assert!(r.step.lagged_mechanical_work_j<0.0);
    near(r.step.joule_energy_j,-r.step.lagged_mechanical_work_j);
    assert_eq!(r.applied.len(),12); // gauge cells retain thermal properties
}

#[test]
fn two_dimensional_depth_and_isolated_zero_field_share_the_same_material_path() {
    let electric=setup(1,20.0); let gate=CancelGate::new(); let mut a=slab2(); let mut b=slab2();
    let r=a.step_magnetothermal_with_material_transport(input(&electric),&[[0.0,0.0,2.0];N],&[[0.0;2];N],DX,&gate,|_,t,_| law(t)).unwrap();
    let s=b.step_magnetothermal_with_material_transport(input(&electric),&[[0.0,0.0,2.0];N],&[[0.0;2];N],2.0*DX,&gate,|_,t,_| law(t)).unwrap();
    near(s.step.joule_energy_j,2.0*r.step.joule_energy_j);
    for (x,y) in a.lattice.gpop.iter().flatten().zip(b.lattice.gpop.iter().flatten()) {near(*x,*y);}
    let isolated=ElectricSetup::default(); let mut c=slab2(); let mut d=slab2();
    let out=c.step_magnetothermal_isolated_with_material_transport(input(&isolated),&[[0.0;3];N],&[[0.0;2];N],DX,1e-10,&gate,|_,t,_| law(t)).unwrap();
    d.step_with_material_transport(&[0.0;N],&[[0.0;2];N],DX,input(&isolated).limits.thermal,&gate,|_,t,_| Ok(fluid(t))).unwrap();
    assert_eq!(out.step.joule_energy_j,0.0); assert_eq!(c.lattice.gpop,d.lattice.gpop);
    let before=a.lattice.grid.f.clone();
    assert!(a.step_magnetothermal_with_material_transport::<FluidTransportState>(input(&electric),&[[1.0,0.0,0.0];N],&[[0.0;2];N],DX,&gate,
        |_,_,_| panic!("unsupported field must precede material queries")).is_err());
    assert_eq!(a.lattice.grid.f,before);
}

#[test]
fn endpoint_refusal_cancellation_and_storage_changes_do_not_publish() {
    let electric=setup(2,20.0); let mut a=slab3();
    let f=a.lattice.grid.f.clone(); let g=a.lattice.gpop.clone();
    let tau=a.lattice.grid.tau.clone(); let force=a.lattice.grid.fext.clone();
    let mut calls=0;
    assert!(a.step_magnetothermal_with_material_transport(input(&electric),&[[0.0,2.0,0.0];N],&[[0.0;3];N],&CancelGate::new(),
        |_,t,_| {calls+=1; if calls>12 {Err(MaterialError::State {what:"endpoint outside source support".into()})} else {law(t)}}).is_err());
    assert_eq!(calls,13);
    let gate=CancelGate::new();
    let cancelled=a.step_magnetothermal_with_material_transport(input(&electric),&[[0.0,2.0,0.0];N],&[[0.0;3];N],&gate,
        |_,t,_| {gate.request(); law(t)});
    assert!(matches!(cancelled,Err(MaterialMagneticError::Material(MaterialTransportError::Cancelled))));
    assert!(a.step_magnetothermal_with_material_transport(input(&electric),&[[0.0;3];N],&[[0.0;3];N],&CancelGate::new(),
        |_,_,_| ElectricalTransport::try_new(FluidTransportState::try_new(1000.0,2000.0,100.0,1.0).unwrap(),Some(1.0))).is_err());
    assert_eq!(a.lattice.grid.f,f); assert_eq!(a.lattice.gpop,g);
    assert_eq!(a.lattice.grid.tau,tau); assert_eq!(a.lattice.grid.fext,force);
    // Retry starts at the original state rather than a partially heated trial.
    let mut fresh=slab3();
    let out=a.step_magnetothermal_with_material_transport(input(&electric),&[[0.0,2.0,0.0];N],&[[0.0;3];N],&CancelGate::new(),|_,t,_| law(t)).unwrap();
    let check=fresh.step_magnetothermal_with_material_transport(input(&electric),&[[0.0,2.0,0.0];N],&[[0.0;3];N],&CancelGate::new(),|_,t,_| law(t)).unwrap();
    assert_eq!(out.step.heat,check.step.heat); assert_eq!(a.lattice.grid.f,fresh.lattice.grid.f);
}

fn buoyant_fluid(t: f64) -> FluidTransportState {
    FluidTransportState::try_new(
        1000.0 - 0.2 * (t - 300.0), 1000.0 + 0.2 * (t - 300.0),
        50.0 + 8.0 * (t - 290.0), 0.5 + 0.03 * (t - 290.0),
    ).unwrap()
}
fn buoyancy<const D: usize>(gravity_m_s2: [f64; D]) -> BoussinesqConfig<D> {
    BoussinesqConfig { gravity_m_s2, max_relative_density_departure: 0.01,
        max_relative_specific_heat_departure: 0.01 }
}
fn buoyant_law(t: f64) -> Result<ElectricalTransport<FluidTransportState>, MaterialError> {
    ElectricalTransport::try_new(buoyant_fluid(t), Some(1.0 + 0.01 * (t - 300.0)))
}

#[test]
fn material_buoyancy_and_lorentz_force_share_the_existing_heat_and_flow_owners() {
    let electric = setup(2, 20.0); let gate = CancelGate::new();
    let mut coupled = slab3(); let mut source_owner = slab3(); let mut reference = slab3();
    // A legacy expansion coefficient must not contribute a second buoyancy.
    coupled.lattice.gbeta = 1e-3; reference.lattice.gbeta = 1e-3;
    let external = [[1.0, 2.0, 0.0]; N]; let field = [[0.0, 2.0, 0.0]; N];
    let config = buoyancy([0.0, 0.0, -9.81]);
    let source = source_owner.step_magnetothermal(input(&electric), &field, &external, &gate,
        |_, t, _| Ok(Some(1.0 + 0.01 * (t - 300.0)))).unwrap();
    let power: Vec<_> = source.motion.electric.joule_w.iter().map(|w| w / DX.powi(3)).collect();
    let forces: Vec<_> = source.motion.lorentz_force_density_n_m3.iter().zip(external)
        .map(|(a, b)| [a[0] + b[0], a[1] + b[1], a[2] + b[2]]).collect();
    let separate = reference.step_with_boussinesq_material_transport(
        &power, &forces, input(&electric).limits.thermal, config, &gate,
        |_, t, _| Ok(buoyant_fluid(t)),
    ).unwrap();
    let out = coupled.step_magnetothermal_with_boussinesq_material_transport(
        input(&electric), &field, &external, config, &gate, |_, t, _| buoyant_law(t),
    ).unwrap();
    assert_eq!(coupled.lattice.grid.f, reference.lattice.grid.f);
    assert_eq!(coupled.lattice.gpop, reference.lattice.gpop);
    assert_eq!(out.step.heat, separate.heat);
    assert_eq!(out.step.motion.electric.joule_w, source.motion.electric.joule_w);
    assert_eq!(coupled.lattice.gbeta, 1e-3);
    near(out.step.heat.source_j, out.step.joule_energy_j);
    let base = initial();
    let scale = plan_scaling_from_material(&base, DX, 0.01, 0.1).unwrap();
    for sample in &out.applied {
        let raw = sample.material.material().resolved_transport();
        assert_eq!(*raw, buoyant_fluid(sample.temperature_k));
        assert_eq!(sample.transport.density_kg_m3(), base.density_kg_m3());
        assert_eq!(sample.transport.specific_heat_j_kg_k(), base.specific_heat_j_kg_k());
        let density_force = (raw.density_kg_m3() - base.density_kg_m3()) * -9.81;
        near(coupled.lattice.grid.fext[sample.cell][2], scale.force_density_lattice(
            forces[sample.cell][2] + density_force,
        ).unwrap());
    }
}

#[test]
fn buoyant_zero_field_recovers_original_transport_in_both_dimensions() {
    let electric = ElectricSetup::default(); let gate = CancelGate::new();
    let mut a = slab2(); let mut b = slab2();
    let config2 = buoyancy([1.0, -9.81]);
    let out = a.step_magnetothermal_isolated_with_boussinesq_material_transport(
        input(&electric), &[[0.0; 3]; N], &[[1.0, 0.0]; N], DX, 1e-10, config2, &gate,
        |_, t, _| buoyant_law(t),
    ).unwrap();
    b.step_with_boussinesq_material_transport(&[0.0; N], &[[1.0, 0.0]; N], DX,
        input(&electric).limits.thermal, config2, &gate, |_, t, _| Ok(buoyant_fluid(t))).unwrap();
    assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
    assert_eq!(a.lattice.gpop, b.lattice.gpop);
    assert_eq!(out.step.joule_energy_j, 0.0);
    let mut a = slab3(); let mut b = slab3();
    let config3 = buoyancy([1.0, 0.0, -9.81]);
    a.step_magnetothermal_isolated_with_boussinesq_material_transport(
        input(&electric), &[[0.0; 3]; N], &[[1.0, 0.0, 0.0]; N], 1e-10, config3, &gate,
        |_, t, _| buoyant_law(t),
    ).unwrap();
    b.step_with_boussinesq_material_transport(&[0.0; N], &[[1.0, 0.0, 0.0]; N],
        input(&electric).limits.thermal, config3, &gate, |_, t, _| Ok(buoyant_fluid(t))).unwrap();
    assert_eq!(a.lattice.grid.f, b.lattice.grid.f);
    assert_eq!(a.lattice.gpop, b.lattice.gpop);
}

#[test]
fn two_dimensional_buoyant_current_preserves_physical_depth_scaling() {
    let electric = setup(1, 20.0); let gate = CancelGate::new();
    let config = buoyancy([0.0, -9.81]);
    let mut a = slab2(); let mut b = slab2();
    let r = a.step_magnetothermal_with_boussinesq_material_transport(
        input(&electric), &[[0.0, 0.0, 2.0]; N], &[[0.0; 2]; N], DX, config, &gate,
        |_, t, _| buoyant_law(t),
    ).unwrap();
    let s = b.step_magnetothermal_with_boussinesq_material_transport(
        input(&electric), &[[0.0, 0.0, 2.0]; N], &[[0.0; 2]; N], 2.0 * DX, config, &gate,
        |_, t, _| buoyant_law(t),
    ).unwrap();
    near(s.step.joule_energy_j, 2.0 * r.step.joule_energy_j);
    near(s.step.heat.stored_change_j, 2.0 * r.step.heat.stored_change_j);
    for (x, y) in a.lattice.grid.f.iter().flatten().zip(b.lattice.grid.f.iter().flatten()) { near(*x, *y); }
    for (x, y) in a.lattice.gpop.iter().flatten().zip(b.lattice.gpop.iter().flatten()) { near(*x, *y); }
}

#[test]
fn raw_buoyant_sources_and_electrical_coefficients_are_sampled_once_per_epoch() {
    let electric = setup(2, 20.0); let mut slab = slab3(); let gate = CancelGate::new();
    let reads = Rc::new(Cell::new(0)); let mut calls = 0;
    let out = slab.step_magnetothermal_with_boussinesq_material_transport(
        input(&electric), &[[0.0, 2.0, 0.0]; N], &[[0.0; 3]; N],
        buoyancy([0.0, 0.0, -9.81]), &gate,
        |_, t, reference| {
            assert_eq!(*reference, initial());
            let id = calls; calls += 1;
            ElectricalTransport::try_new(Receipt { id, state: buoyant_fluid(t), reads: reads.clone() },
                Some(1.0 + 0.01 * (t - 300.0)))
        },
    ).unwrap();
    assert_eq!(calls, 24); assert_eq!(reads.get(), 24);
    for (i, sample) in out.applied.iter().enumerate() {
        assert_eq!(sample.material.material().material().id, i);
        assert_eq!(out.endpoint[i].material.material().material().id, 12 + i);
        near(sample.material.resistivity_ohm_m().unwrap(), 1.0 + 0.01 * (sample.temperature_k - 300.0));
        near(out.endpoint[i].temperature_k,
            slab.temperature_si(slab.lattice.gpop[sample.cell].iter().sum()).unwrap());
    }
}

#[test]
fn data_only_material_change_alters_electrical_thermal_and_buoyant_response() {
    let electric = setup(2, 20.0); let gate = CancelGate::new();
    let mut runs = Vec::new();
    // Two synthetic material datasets; geometry, solver and boundary data stay fixed.
    for (density_slope, conductivity, viscosity, resistivity) in [
        (0.2, 100.0, 1.0, 1.0), (0.4, 200.0, 1.5, 2.0),
    ] {
        let mut slab = slab3();
        let out = slab.step_magnetothermal_with_boussinesq_material_transport(
            input(&electric), &[[0.0, 2.0, 0.0]; N], &[[0.0; 3]; N],
            buoyancy([0.0, 0.0, -9.81]), &gate,
            |_, t, _| ElectricalTransport::try_new(FluidTransportState::try_new(
                1000.0 - density_slope * (t - 300.0), 1000.0, conductivity, viscosity,
            )?, Some(resistivity)),
        ).unwrap();
        runs.push((slab, out));
    }
    let (a, x) = &runs[0]; let (b, y) = &runs[1];
    near(x.step.joule_energy_j, 2.0 * y.step.joule_energy_j);
    near(b.lattice.grid.fext[6][2], 2.0 * a.lattice.grid.fext[6][2]);
    assert_ne!(x.applied[0].momentum_relaxation, y.applied[0].momentum_relaxation);
    assert_ne!(x.applied[0].thermal_relaxation, y.applied[0].thermal_relaxation);
    assert_ne!(a.lattice.grid.f, b.lattice.grid.f);
    assert_ne!(a.lattice.gpop, b.lattice.gpop);
}

#[test]
fn buoyant_endpoint_support_and_cancellation_refuse_atomically() {
    let electric = setup(2, 20.0); let mut slab = slab3(); slab.lattice.gbeta = 1e-3;
    let initial_f = slab.lattice.grid.f.clone(); let initial_g = slab.lattice.gpop.clone();
    let initial_loads = slab.lattice.grid.fext.clone(); let initial_tau = slab.lattice.grid.tau.clone();
    let config = buoyancy([0.0, 0.0, -9.81]);
    for change_capacity in [false, true] {
        let mut calls = 0;
        let result = slab.step_magnetothermal_with_boussinesq_material_transport(
            input(&electric), &[[0.0, 2.0, 0.0]; N], &[[0.0; 3]; N], config, &CancelGate::new(),
            |_, t, _| {
                calls += 1;
                if calls <= 12 { return buoyant_law(t); }
                let (rho, cp) = if change_capacity { (1000.0, 1100.0) } else { (900.0, 1000.0) };
                ElectricalTransport::try_new(FluidTransportState::try_new(rho, cp, 100.0, 1.0)?, Some(1.0))
            },
        );
        assert!(result.is_err()); assert_eq!(calls, 13);
    }
    let gate = CancelGate::new(); let mut calls = 0;
    let result = slab.step_magnetothermal_with_boussinesq_material_transport(
        input(&electric), &[[0.0, 2.0, 0.0]; N], &[[0.0; 3]; N], config, &gate,
        |_, t, _| { calls += 1; if calls == 13 { gate.request(); } buoyant_law(t) },
    );
    assert!(matches!(result, Err(MaterialMagneticError::Material(MaterialTransportError::Cancelled))));
    assert_eq!(slab.lattice.grid.f, initial_f); assert_eq!(slab.lattice.gpop, initial_g);
    assert_eq!(slab.lattice.grid.fext, initial_loads); assert_eq!(slab.lattice.grid.tau, initial_tau);
    assert_eq!(slab.lattice.gbeta, 1e-3);
    let mut fresh = slab3(); fresh.lattice.gbeta = 1e-3;
    slab.step_magnetothermal_with_boussinesq_material_transport(
        input(&electric), &[[0.0, 2.0, 0.0]; N], &[[0.0; 3]; N], config, &CancelGate::new(),
        |_, t, _| buoyant_law(t),
    ).unwrap();
    fresh.step_magnetothermal_with_boussinesq_material_transport(
        input(&electric), &[[0.0, 2.0, 0.0]; N], &[[0.0; 3]; N], config, &CancelGate::new(),
        |_, t, _| buoyant_law(t),
    ).unwrap();
    assert_eq!(slab.lattice.grid.f, fresh.lattice.grid.f);
    assert_eq!(slab.lattice.gpop, fresh.lattice.gpop);
}
