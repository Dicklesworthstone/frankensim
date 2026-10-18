//! Real kernels with synthetic local laws; no calibration authority is implied.
use fs_exec::CancelGate;
use fs_lbm::material::{MaterialSlab, MaterialSlab3, MaterialThermalSlab, MaterialThermalSlab3, plan_scaling_from_material};
use fs_lbm::material::transport::MaterialTransportError;
use fs_lbm::material::transport::boussinesq::BoussinesqConfig;
use fs_lbm::thermal::sources::HeatSourceLimits;
use fs_material::{MaterialError, fluid_transport::{FluidTransport, FluidTransportState}};
use fs_material::gas::{GasSpec, GasState};

fn material(k: f64, mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(1000.0,1000.0,k,mu).unwrap()
}
fn slab3() -> MaterialSlab3 {
    plan_scaling_from_material(&material(1000.0,1.0),0.01,0.01,0.1).unwrap()
        .thermal_slab3(MaterialThermalSlab3 { nx:2,ny:1,nz:6,
            bottom_temperature_k:310.0,top_temperature_k:300.0,
            gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:16 }).unwrap()
}
fn slab2() -> MaterialSlab {
    plan_scaling_from_material(&material(1000.0,1.0),0.01,0.01,0.1).unwrap()
        .thermal_slab(MaterialThermalSlab { nx:2,ny:6,
            bottom_temperature_k:310.0,top_temperature_k:300.0,
            gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:16 }).unwrap()
}
fn limits(n: usize) -> HeatSourceLimits {
    HeatSourceLimits { max_cells:n,temperature_range_k:[290.0,330.0],balance_tolerance_j:1e-8 }
}
fn near(a:f64,b:f64,tol:f64) { assert!((a-b).abs()<=tol,"{a} != {b}"); }
fn power() -> Vec<f64> { let mut p=vec![0.0;16];p[4..12].fill(1e6);p }

#[test]
fn constant_material_is_the_existing_physical_step_in_both_dimensions() {
    let mut a=slab3();let mut b=slab3();let gate=CancelGate::new();
    let force=[[1.0,0.0,0.0];16];let p=power();
    let old=a.step_with_power_density(&p,&force,limits(16)).unwrap();
    let new=b.step_with_material_transport(&p,&force,limits(16),&gate,|_,_,_|Ok(material(1000.0,1.0))).unwrap();
    assert_eq!(new.heat,old);assert_eq!(a.lattice.grid.f,b.lattice.grid.f);assert_eq!(a.lattice.gpop,b.lattice.gpop);
    assert_eq!(new.applied.len(),12);assert_eq!(new.endpoint.len(),12);
    for s in &new.applied {
        near(s.momentum_relaxation,0.8,1e-15);near(s.thermal_relaxation,0.9,1e-15);
    }
    let mut a=slab2();let mut b=slab2();let force=[[1.0,0.0];16];
    let old=a.step_with_power_density(&p,&force,0.01,limits(16)).unwrap();
    let new=b.step_with_material_transport(&p,&force,0.01,limits(16),&gate,|_,_,_|Ok(material(1000.0,1.0))).unwrap();
    assert_eq!(new.heat,old);assert_eq!(a.lattice.grid.f,b.lattice.grid.f);assert_eq!(a.lattice.gpop,b.lattice.gpop);
    near(new.applied[0].thermal_relaxation,0.8,1e-15); // not the D3Q7 value
}

#[derive(Debug)]
struct SyntheticTicket { value:FluidTransportState, query_k:f64, condition:u32 }
impl FluidTransport for SyntheticTicket {
    fn transport_state(&self)->Result<FluidTransportState,MaterialError> { Ok(self.value) }
}

#[test]
fn actual_local_temperatures_and_concrete_source_payloads_survive_each_step() {
    let mut a=slab3();let mut control=slab3();let gate=CancelGate::new();
    a.lattice.grid.f[6][0]-=0.002;a.lattice.grid.f[6][1]+=0.001;a.lattice.grid.f[6][2]+=0.001;
    control.lattice.grid.f.clone_from(&a.lattice.grid.f);
    let p=power();let force=[[0.0;3];16];
    for _ in 0..2 {
        let old:Vec<_>=a.lattice.gpop.iter().map(|p|a.temperature_si(p.iter().sum()).unwrap()).collect();
        let r=a.step_with_material_transport(&p,&force,limits(16),&gate,|_,t,retained| {
            assert_eq!(*retained,material(1000.0,1.0));
            Ok(SyntheticTicket { value:material(1000.0*(1.0+(t-300.0)/10.0),1.0+(t-300.0)/10.0),query_k:t,condition:17 })
        }).unwrap();
        for s in &r.applied { assert_eq!(s.temperature_k,old[s.cell]);assert_eq!(s.material.query_k,s.temperature_k);assert_eq!(s.material.condition,17); }
        for s in &r.endpoint { assert_eq!(s.temperature_k,a.temperature_si(a.lattice.gpop[s.cell].iter().sum()).unwrap()); }
        assert!(r.applied[0].transport.thermal_conductivity_w_m_k()>r.applied[11].transport.thermal_conductivity_w_m_k());
        assert_eq!(a.lattice.grid.tau,vec![0.8;16]); // step-local coefficients, retained history
        control.step_with_power_density(&p,&force,limits(16)).unwrap();
    }
    assert_ne!(a.lattice.grid.f,control.lattice.grid.f);assert_ne!(a.lattice.gpop,control.lattice.gpop);
}

fn nonlinear_error(n:usize)->f64 {
    let dx=1.0/n as f64;let dt=80.0*dx*dx;
    let mut a=plan_scaling_from_material(&material(1000.0,1.0),dx,dt,0.0).unwrap()
        .thermal_slab3(MaterialThermalSlab3 {nx:1,ny:1,nz:n,bottom_temperature_k:310.0,
            top_temperature_k:300.0,gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:n+2}).unwrap();
    let gate=CancelGate::new();let p=vec![0.0;n+2];let forces=vec![[0.0;3];n+2];
    for _ in 0..30*n*n {
        a.step_with_material_transport(&p,&forces,limits(n+2),&gate,
            |_,t,_|Ok(material(1000.0*(1.0+(t-300.0)/10.0),1.0))).unwrap();
    }
    let mut error=0.0_f64;
    for i in 1..=n {
        let y=(i as f64-0.5)/n as f64;
        let exact=300.0+10.0*((4.0-3.0*y).sqrt()-1.0);
        error=error.max((a.temperature_si(a.lattice.temperature(0,0,i)).unwrap()-exact).abs());
    }
    error
}
#[test]
fn nonlinear_material_conduction_matches_a_kirchhoff_transform_reference() {
    let coarse=nonlinear_error(8);let fine=nonlinear_error(16);
    assert!(coarse<0.04);assert!(fine<0.3*coarse);
}

#[test]
fn density_capacity_and_endpoint_source_failures_leave_history_and_units_intact() {
    let mut a=slab3();let f=a.lattice.grid.f.clone();let g=a.lattice.gpop.clone();
    let tau=a.lattice.grid.tau.clone();let loads=a.lattice.grid.fext.clone();let gate=CancelGate::new();
    for changed in [
        FluidTransportState::try_new(999.0,1000.0,1000.0,1.0).unwrap(),
        FluidTransportState::try_new(1000.0,1001.0,1000.0,1.0).unwrap(),
    ] {
        let err=a.step_with_material_transport(&power(),&[[0.0;3];16],limits(16),&gate,|_,_,_|Ok(changed)).unwrap_err();
        assert!(matches!(err,MaterialTransportError::Invalid(_)));
    }
    let threshold=a.temperature_si(a.lattice.gpop[6].iter().sum()).unwrap()+0.01;
    let mut p=vec![0.0;16];p[6]=1e8;
    let err=a.step_with_material_transport(&p,&[[0.0;3];16],limits(16),&gate,|cell,t,_| {
        if cell==6 && t>threshold { Err(MaterialError::State { what: "outside the local source temperature support".into() }) }
        else { Ok(material(1000.0,1.0)) }
    }).unwrap_err();
    assert!(matches!(err,MaterialTransportError::State {cell:6,..}));
    assert_eq!(a.lattice.grid.f,f);assert_eq!(a.lattice.gpop,g);assert_eq!(a.lattice.grid.tau,tau);assert_eq!(a.lattice.grid.fext,loads);
    // Retry reaches the original accepted state, with no duplicated heating.
    let mut fresh=slab3();a.step_with_power_density(&power(),&[[0.0;3];16],limits(16)).unwrap();
    fresh.step_with_power_density(&power(),&[[0.0;3];16],limits(16)).unwrap();assert_eq!(a.lattice.gpop,fresh.lattice.gpop);
}

#[test]
fn cancelled_endpoint_and_invalid_inputs_never_publish_a_partial_step() {
    let mut a=slab3();let f=a.lattice.grid.f.clone();let g=a.lattice.gpop.clone();let gate=CancelGate::new();
    let mut calls=0;
    let result=a.step_with_material_transport(&power(),&[[0.0;3];16],limits(16),&gate,|_,_,_| {
        calls+=1;if calls==13 {gate.request();}Ok(material(1000.0,1.0))
    });
    assert!(matches!(result,Err(MaterialTransportError::Cancelled)));
    assert_eq!(a.lattice.grid.f,f);assert_eq!(a.lattice.gpop,g);
    for p in [vec![0.0;15],vec![1.0;16]] {
        assert!(a.step_with_material_transport(&p,&[[0.0;3];16],limits(16),&CancelGate::new(),
            |_,_,_|->Result<FluidTransportState,MaterialError>{panic!("preflight before material query")}).is_err());
    }
}

#[test]
fn depth_and_reversed_temperature_coordinates_preserve_physical_heating() {
    let gate=CancelGate::new();let mut a=slab2();let mut b=slab2();
    let r=a.step_with_material_transport(&power(),&[[0.0;2];16],0.01,limits(16),&gate,|_,t,_|Ok(material(1000.0+t,1.0))).unwrap();
    let s=b.step_with_material_transport(&power(),&[[0.0;2];16],0.02,limits(16),&gate,|_,t,_|Ok(material(1000.0+t,1.0))).unwrap();
    assert_eq!(a.lattice.gpop,b.lattice.gpop);near(s.heat.source_j,2.0*r.heat.source_j,1e-12);
    let mut inverted=plan_scaling_from_material(&material(1000.0,1.0),0.01,0.01,0.0).unwrap()
        .thermal_slab(MaterialThermalSlab {nx:2,ny:6,bottom_temperature_k:300.0,top_temperature_k:310.0,
            gravity_m_s2:0.0,expansion_per_k:0.0,max_cells:16}).unwrap();
    let before=inverted.temperature_si(inverted.lattice.temperature(0,3)).unwrap();
    let r=inverted.step_with_material_transport(&power(),&[[0.0;2];16],0.01,limits(16),&gate,|_,_,_|Ok(material(1000.0,1.0))).unwrap();
    assert!(r.heat.source_j>0.0);assert!(inverted.temperature_si(inverted.lattice.temperature(0,3)).unwrap()>before);
}

#[derive(Debug)]
struct GasQuery {
    gas: GasState,
    reads: std::cell::Cell<usize>,
}
impl FluidTransport for GasQuery {
    fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
        self.reads.set(self.reads.get() + 1);
        assert_eq!(self.reads.get(), 1, "the original material must be sampled once");
        self.gas.transport_state()
    }
}

#[test]
fn g1_g3_gas_eos_drives_buoyancy_and_transport_with_conserved_reference_storage() {
    let spec = GasSpec::dry_air_ussa1976();
    let gate = CancelGate::new();
    let mut relaxation = Vec::new();
    for pressure in [100_000.0, 200_000.0] {
        let reference = GasState::try_new(&spec, 305.0, pressure).unwrap();
        let stored = reference.transport_state().unwrap();
        let scale = plan_scaling_from_material(&reference, 0.01, 0.001, 0.01).unwrap();
        let mut fluid = scale.thermal_slab(MaterialThermalSlab {
            nx: 2, ny: 6, bottom_temperature_k: 310.0, top_temperature_k: 300.0,
            gravity_m_s2: 9.81, expansion_per_k: 1.0 / 305.0, max_cells: 16,
        }).unwrap();
        let original_beta = fluid.lattice.gbeta;
        let mass_before: f64 = fluid.lattice.grid.f[2..14].iter().flatten().sum();
        let mut heat = vec![0.0; 16];
        heat[2..14].fill(1000.0);
        let result = fluid.step_with_boussinesq_material_transport(
            &heat, &[[0.0; 2]; 16], 0.01, limits(16),
            BoussinesqConfig { gravity_m_s2: [0.0, -9.81],
                max_relative_density_departure: 0.05, max_relative_specific_heat_departure: 0.01 },
            &gate, |_, t, retained| {
                assert_eq!(*retained, stored);
                Ok(GasQuery { gas: GasState::try_new(&spec, t, pressure)?,
                    reads: std::cell::Cell::new(0) })
            },
        ).unwrap();
        near(result.heat.source_j, 12.0 * 1000.0 * 1e-6 * 0.001, 1e-15);
        assert!(result.heat.defect_j.abs() < 1e-10);
        let mass_after: f64 = fluid.lattice.grid.f[2..14].iter().flatten().sum();
        near(mass_after, mass_before, 1e-12);
        assert_eq!(fluid.lattice.gbeta, original_beta);
        for evaluated in result.applied.iter().chain(&result.endpoint) {
            assert_eq!(evaluated.material.material().reads.get(), 1);
            assert_eq!(evaluated.transport.density_kg_m3(), stored.density_kg_m3());
            assert_eq!(evaluated.transport.specific_heat_j_kg_k(), stored.specific_heat_j_kg_k());
            assert_eq!(&evaluated.transport, evaluated.material.reference_storage_transport());
            let source = evaluated.material.resolved_transport();
            // Independent constant-pressure ideal-gas relation, not beta*dT.
            near(source.density_kg_m3(), stored.density_kg_m3() * 305.0 / evaluated.temperature_k, 1e-12);
            near(evaluated.momentum_relaxation, 0.5 + 3.0 * source.dynamic_viscosity_pa_s()
                / stored.density_kg_m3() * 0.001 / 0.01_f64.powi(2), 1e-14);
            near(evaluated.thermal_relaxation, 0.5 + 3.0 * source.thermal_conductivity_w_m_k()
                / (stored.density_kg_m3() * stored.specific_heat_j_kg_k())
                * 0.001 / 0.01_f64.powi(2), 1e-14);
        }
        for applied in &result.applied {
            let expected = -9.81 * (305.0 / applied.temperature_k - 1.0)
                * 0.001_f64.powi(2) / 0.01;
            near(fluid.lattice.grid.fext[applied.cell][1], expected, 1e-16);
            assert_eq!(fluid.lattice.grid.fext[applied.cell][0], 0.0);
        }
        assert!(fluid.lattice.grid.fext[2][1] > 0.0);
        assert!(fluid.lattice.grid.fext[13][1] < 0.0);
        relaxation.push(result.applied[0].momentum_relaxation - 0.5);
    }
    // Raising pressure changes nu through the material density; buoyancy
    // acceleration still has the same ideal-gas temperature dependence.
    near(relaxation[1], 0.5 * relaxation[0], 1e-14);
}

#[test]
fn g3_boussinesq_three_dimensional_response_replaces_beta_and_keeps_external_loads() {
    let mut with_background = slab3();
    let mut without_background = slab3();
    with_background.lattice.gbeta = 0.3;
    let config = BoussinesqConfig { gravity_m_s2: [1.0, -2.0, -9.81],
        max_relative_density_departure: 0.01, max_relative_specific_heat_departure: 0.02 };
    let law = |_: usize, t: f64, _: &FluidTransportState| {
        FluidTransportState::try_new(1000.0 * (1.0 - 0.0001 * (t - 305.0)),
            1000.0 + t - 300.0, 1000.0 + 10.0 * (t - 300.0), 1.0 + (t - 300.0) / 100.0)
    };
    let force = [[2.0, 3.0, 4.0]; 16];
    let gate = CancelGate::new();
    let a = with_background.step_with_boussinesq_material_transport(
        &power(), &force, limits(16), config, &gate, law,
    ).unwrap();
    let b = without_background.step_with_boussinesq_material_transport(
        &power(), &force, limits(16), config, &gate, law,
    ).unwrap();
    assert_eq!(with_background.lattice.gbeta, 0.3);
    assert_eq!(without_background.lattice.gbeta, 0.0);
    assert_eq!(with_background.lattice.grid.f, without_background.lattice.grid.f);
    assert_eq!(with_background.lattice.gpop, without_background.lattice.gpop);
    assert_eq!(a.heat, b.heat);
    for applied in &a.applied {
        let raw = applied.material.resolved_transport();
        assert_ne!(raw.specific_heat_j_kg_k(), 1000.0);
        assert_eq!(applied.transport.specific_heat_j_kg_k(), 1000.0);
        near(applied.thermal_relaxation,
            0.5 + 4.0 * raw.thermal_conductivity_w_m_k() / 1e6 * 100.0, 1e-14);
        for axis in 0..3 {
            let expected = (force[applied.cell][axis]
                + (raw.density_kg_m3() - 1000.0) * config.gravity_m_s2[axis]) * 1e-5;
            near(with_background.lattice.grid.fext[applied.cell][axis], expected, 1e-16);
        }
    }
}

#[test]
fn g0_boussinesq_bounds_and_g4_endpoint_failure_preserve_all_accepted_state() {
    let mut fluid = slab3();
    let before_f = fluid.lattice.grid.f.clone();
    let before_g = fluid.lattice.gpop.clone();
    let before_tau = fluid.lattice.grid.tau.clone();
    let before_force = fluid.lattice.grid.fext.clone();
    let config = BoussinesqConfig { gravity_m_s2: [0.0, 0.0, -9.81],
        max_relative_density_departure: 0.001, max_relative_specific_heat_departure: 0.0 };
    let gate = CancelGate::new();
    for (rho, cp) in [(998.0, 1000.0), (1000.0, 1001.0)] {
        let result = fluid.step_with_boussinesq_material_transport(
            &power(), &[[0.0; 3]; 16], limits(16), config, &gate,
            |_, _, _| FluidTransportState::try_new(rho, cp, 1000.0, 1.0),
        );
        assert!(matches!(result, Err(MaterialTransportError::State { .. })));
    }
    let threshold = fluid.temperature_si(fluid.lattice.gpop[6].iter().sum()).unwrap() + 0.1;
    let mut heat = vec![0.0; 16];
    heat[6] = 1e8;
    let result = fluid.step_with_boussinesq_material_transport(
        &heat, &[[0.0; 3]; 16], limits(16), config, &gate, |cell, t, _| {
            FluidTransportState::try_new(if cell == 6 && t > threshold { 998.0 } else { 1000.0 },
                1000.0, 1000.0, 1.0)
        },
    );
    assert!(matches!(result, Err(MaterialTransportError::State { cell: 6, .. })));
    let cancel = CancelGate::new();
    let mut queries = 0;
    let result = fluid.step_with_boussinesq_material_transport(
        &power(), &[[0.0; 3]; 16], limits(16), config, &cancel, |_, _, _| {
            queries += 1;
            if queries == 13 { cancel.request(); }
            Ok(material(1000.0, 1.0))
        },
    );
    assert!(matches!(result, Err(MaterialTransportError::Cancelled)));
    assert_eq!(fluid.lattice.grid.f, before_f);
    assert_eq!(fluid.lattice.gpop, before_g);
    assert_eq!(fluid.lattice.grid.tau, before_tau);
    assert_eq!(fluid.lattice.grid.fext, before_force);
    assert_eq!(fluid.lattice.gbeta, 0.0);
    let mut invalid = config;
    invalid.gravity_m_s2[0] = f64::NAN;
    assert!(fluid.step_with_boussinesq_material_transport(
        &power(), &[[0.0; 3]; 16], limits(16), invalid, &gate,
        |_, _, _| -> Result<FluidTransportState, MaterialError> { panic!("configuration preflight before source query") },
    ).is_err());
    let mut fresh = slab3();
    let accepted = fluid.step_with_boussinesq_material_transport(
        &power(), &[[0.0; 3]; 16], limits(16), config, &gate,
        |_, _, _| Ok(material(1000.0, 1.0)),
    ).unwrap();
    let reference = fresh.step_with_power_density(&power(), &[[0.0; 3]; 16], limits(16)).unwrap();
    assert_eq!(accepted.heat, reference);
    assert_eq!(fluid.lattice.grid.f, fresh.lattice.grid.f);
    assert_eq!(fluid.lattice.gpop, fresh.lattice.gpop);
}
