//! G0/G1/G3: original full pHS equations, analytic single-mode response and source swaps.
use super::*;
use crate::air_path::cylinder::{CylinderFlowResolution, frame::{FrameFluidBudget, FrameFluidRegion}};
use crate::frame_vibration::FrameVibrationBudget;
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_qty::{Density, DynViscosity, Pressure};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, FrameViscousSection, Geometry};

fn section(rho: f64, young: f64, eta: f64) -> FrameViscousSection {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, dims, value) in [("density",Density::DIMS,rho),
        ("young_modulus",Pressure::DIMS,young),("eta",DynViscosity::DIMS,eta)] {
        let mut validity = ValidityDomain::unconstrained().with("T",300.0,400.0);
        if name == "eta" { validity = validity.with("omega",0.0,1e6); }
        claims.insert_claim(PropertyClaim { key: PropertyKey::new(name,dims),
            value: PropertyValue::Scalar { value,dims }, validity,
            uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![], provenance: Provenance { source:"manufactured harmonic fixture".into(),
                license:"CC0-1.0".into(),artifact:None } }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry:"synthetic".into(),phase:"solid".into(),
        process:"test-only constants".into(),revision:0 },claims,vec![]).unwrap();
    let query = QueryPoint::new().with("T",300.0).unwrap().with("omega",1.0).unwrap();
    let material = resolve_material_state_point(&card,&query,&requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let elastic = ConcentricElasticSection::try_new(0.0,vec![ConcentricMaterialLayer {
        outer_radius_m:0.02,material }]).unwrap();
    FrameViscousSection::try_new(elastic,"eta","eta").unwrap()
}
fn frame(two: bool, rho: f64, young: f64, eta: f64) -> Frame2d {
    let mut result = Frame2d::with_mass_formulation(Geometry::Linear,FrameMassFormulation::ConsistentReference);
    result.add_node(0.0,0.0);
    result.add_node(if two {0.5} else {1.0},0.0);
    if two { result.add_node(1.0,0.0); }
    result.fix(0,[true;3]);
    result.fix(1,[true,false,!two]);
    if two { result.fix(2,[true,false,false]); }
    for element in 0..(if two {2} else {1}) {
        result.add_material_viscous_element(element,element+1,3,section(rho,young,eta)).unwrap();
    }
    result
}
fn fluid(rho: f64, mu: f64) -> FluidTransportState {
    FluidTransportState::try_new(rho,4000.0,0.5,mu).unwrap()
}
fn model(frame: &Frame2d, modes: usize, regions: Vec<FrameFluidRegion<FluidTransportState>>)
    -> MaterialFluidFrame<'_,FluidTransportState> {
    MaterialFluidFrame::new_with_regions(frame,modes,regions,
        CylinderFlowResolution {minimum_wavenumber_radius:0.01,maximum_wavenumber_radius:20.0,intervals:4},
        FrameVibrationBudget {max_frame_dofs:12,max_dense_entries:144},
        FrameFluidBudget {max_states:128,max_dense_entries:16384,max_projection_terms:4096,mass_inverse_tolerance:1e-10},
        &CancelGate::new_clock_free()).unwrap()
}
fn budget() -> FrameHarmonicBudget {
    FrameHarmonicBudget {max_modal_entries:16,max_setup_terms:4096,max_response_values:256,
        equation_relative_tolerance:1e-11,power_absolute_tolerance_w:1e-12,power_relative_tolerance:1e-8}
}
fn near(a: C64,b: C64,tolerance: f64) {
    assert!((a-b).abs() <= tolerance*(1.0+a.abs()+b.abs()),"{a:?} != {b:?}");
}

#[test]
fn condensed_response_satisfies_every_original_full_state_equation() {
    let frame = frame(true,1000.0,2e6,500.0);
    let model = model(&frame,4,vec![
        FrameFluidRegion {element:1,span:[0.0,0.37],fluid:fluid(1000.0,0.01)},
        FrameFluidRegion {element:1,span:[0.37,1.0],fluid:fluid(3.0,0.002)},
    ]);
    let n = model.system().state_dim();
    let m = 4;
    let (j,r,g) = model.system().structure();
    let q = model.system().storage().exact_quadratic_matrix().unwrap();
    let gate = CancelGate::new_clock_free();
    for omega in [0.1,1.0,10.0,100.0,1000.0] {
        let prepared = model.prepare_harmonic(omega,budget(),&gate).unwrap();
        assert!(core::ptr::eq(prepared.model(),&model));
        assert_eq!(prepared.factor_dimension(),4);
        assert_eq!(prepared.dynamic_stiffness().len(),16);
        let force = [C64::new(0.01,0.002),C64::new(-0.003,0.008),
            C64::new(0.004,-0.002),C64::new(-0.001,-0.003)];
        let response = prepared.solve_modal(&force,&gate).unwrap();
        let mut pencil = vec![C64::ZERO;n*n];
        let mut forcing = vec![C64::ZERO;n];
        for i in 0..n {
            for k in 0..n {
                let a: f64 = (0..n).map(|p| (j[i*n+p]-r[i*n+p])*q[p*n+k]).sum();
                pencil[i*n+k] = C64::new(-a,if i==k {omega} else {0.0});
            }
            for k in 0..m { forcing[i] = forcing[i]+force[k].scale(g[i*m+k]); }
        }
        // Independent full pencil, not the condensed matrix or branch formula.
        let full_factor = lu_complex_with_cancel(&pencil,n,||false).unwrap();
        let mut full = forcing.clone();
        full_factor.solve_with_cancel(&mut full,||false).unwrap();
        for (&actual,&expected) in response.state.iter().zip(&full) { near(actual,expected,2e-9); }
        assert!(equation_residual(&pencil,&response.state,&forcing,&gate).unwrap()<1e-10);
        let mut effort = vec![C64::ZERO;n];
        for i in 0..n { for k in 0..n { effort[i] = effort[i]+response.state[k].scale(q[i*n+k]); }}
        let energy: f64 = (0..n).map(|i| 0.25*(response.state[i].conj()*effort[i]).re).sum();
        let mut loss = 0.0;
        for i in 0..n { for k in 0..n { loss += 0.5*r[i*n+k]*(effort[i].conj()*effort[k]).re; }}
        near(C64::new(response.mean_stored_energy_j,0.0),C64::new(energy,0.0),1e-10);
        near(C64::new(response.supplied_power_w,0.0),C64::new(loss,0.0),1e-10);
        assert!(response.fluid_region_loss_w.iter().all(|p| *p>0.0));
        assert!(response.material_port_loss_w.iter().all(|p| *p>=0.0));
        assert!(response.power_residual_w.abs()<1e-11);
    }
}

#[test]
fn dry_single_mode_matches_physical_dynamic_stiffness_and_amplitude_scaling() {
    let frame = frame(false,1000.0,2e6,500.0);
    let model = model(&frame,1,vec![]);
    let gate = CancelGate::new_clock_free();
    let omega = 12.0;
    let prepared = model.prepare_harmonic(omega,budget(),&gate).unwrap();
    let k = model.basis().modes()[0].eigenvalue_per_s2;
    let shape = model.basis().modes()[0].displacement[4];
    let c = 12.0*500.0*core::f64::consts::PI*0.02_f64.powi(4)/4.0*shape*shape;
    let force = C64::new(0.01,0.003);
    let a = prepared.solve_modal(&[force],&gate).unwrap();
    near(a.modal_displacement[0],force/C64::new(k-omega*omega,omega*c),1e-12);
    let scale = C64::new(-2.0,3.0);
    let b = prepared.solve_modal(&[force*scale],&gate).unwrap();
    near(b.modal_displacement[0],a.modal_displacement[0]*scale,1e-12);
    near(C64::new(b.supplied_power_w,0.0),C64::new(13.0*a.supplied_power_w,0.0),1e-12);
    near(C64::new(b.mean_stored_energy_j,0.0),C64::new(13.0*a.mean_stored_energy_j,0.0),1e-12);
    assert!(a.fluid_region_loss_w.is_empty());
}

#[test]
fn partitioning_uniform_fluid_preserves_complex_motion_and_total_heat() {
    let frame = frame(false,1000.0,2e6,500.0);
    let a = model(&frame,1,vec![FrameFluidRegion {element:0,span:[0.0,1.0],fluid:fluid(400.0,0.01)}]);
    let b = model(&frame,1,vec![
        FrameFluidRegion {element:0,span:[0.0,0.37],fluid:fluid(400.0,0.01)},
        FrameFluidRegion {element:0,span:[0.37,1.0],fluid:fluid(400.0,0.01)},
    ]);
    let gate = CancelGate::new_clock_free();
    for omega in [1.0,30.0,100.0] {
        let a = a.prepare_harmonic(omega,budget(),&gate).unwrap().solve_modal(&[C64::new(0.01,0.003)],&gate).unwrap();
        let b = b.prepare_harmonic(omega,budget(),&gate).unwrap().solve_modal(&[C64::new(0.01,0.003)],&gate).unwrap();
        near(a.modal_displacement[0],b.modal_displacement[0],1e-11);
        near(C64::new(a.fluid_region_loss_w.iter().sum(),0.0),C64::new(b.fluid_region_loss_w.iter().sum(),0.0),1e-11);
        near(C64::new(a.mean_stored_energy_j,0.0),C64::new(b.mean_stored_energy_j,0.0),1e-11);
    }
}

#[test]
fn changing_only_source_data_changes_the_forced_response() {
    let frames = [frame(false,1000.0,2e6,500.0),frame(false,2000.0,2e6,500.0),
        frame(false,1000.0,4e6,500.0)];
    let mut results = Vec::new();
    let gate = CancelGate::new_clock_free();
    for (frame,rho,mu) in [(&frames[0],400.0,0.01),(&frames[1],400.0,0.01),
        (&frames[2],400.0,0.01),(&frames[0],1600.0,0.01),(&frames[0],400.0,0.1)] {
        let model = model(frame,1,vec![FrameFluidRegion {element:0,span:[0.0,1.0],fluid:fluid(rho,mu)}]);
        let shape = model.basis().modes()[0].displacement[4];
        let response = model.prepare_harmonic(30.0,budget(),&gate).unwrap()
            .solve_modal(&[C64::new(0.01*shape,0.0)],&gate).unwrap();
        results.push(response.modal_displacement[0].scale(shape));
    }
    for &changed in &results[1..] { assert!((changed-results[0]).abs()>1e-9); }
}

#[test]
fn refusals_leave_source_factor_and_callers_unchanged_for_retry() {
    let frame = frame(false,1000.0,2e6,500.0);
    let model = model(&frame,1,vec![]);
    let original = frame.displacements().to_vec();
    let gate = CancelGate::new_clock_free();
    for omega in [0.0,-1.0,f64::NAN,f64::INFINITY,1e7,f64::MIN_POSITIVE] {
        assert!(model.prepare_harmonic(omega,budget(),&gate).is_err());
    }
    for bad in [FrameHarmonicBudget {max_modal_entries:0,..budget()},
        FrameHarmonicBudget {max_setup_terms:0,..budget()},
        FrameHarmonicBudget {max_response_values:0,..budget()},
        FrameHarmonicBudget {equation_relative_tolerance:0.0,..budget()}] {
        assert!(model.prepare_harmonic(10.0,bad,&gate).is_err());
    }
    let prepared = model.prepare_harmonic(10.0,budget(),&gate).unwrap();
    let force = [C64::new(0.01,0.003)];
    let expected = prepared.solve_modal(&force,&gate).unwrap();
    assert!(prepared.solve_modal(&[],&gate).is_err());
    assert!(prepared.solve_modal(&[C64::new(f64::NAN,0.0)],&gate).is_err());
    let cancelled = CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(prepared.solve_modal(&force,&cancelled),Err(FrameHarmonicError::Cancelled)));
    let actual = prepared.solve_modal(&force,&gate).unwrap();
    assert_eq!(actual.state,expected.state);
    assert_eq!(frame.displacements(),original.as_slice());
    assert_eq!(force,[C64::new(0.01,0.003)]);
}

#[test]
fn fast_memory_slip_is_not_lost_to_subtraction_and_zero_forcing_stays_zero() {
    let (follow,slip) = memory_transfer(1e16,1.0).unwrap();
    assert_eq!(follow.re,1.0);
    assert!(slip.re>0.0 && slip.im>0.0);
    near(follow+slip,C64::new(1.0,0.0),1e-15);
    let frame = frame(false,1000.0,2e6,500.0);
    let model = model(&frame,1,vec![FrameFluidRegion {element:0,span:[0.0,1.0],fluid:fluid(400.0,0.01)}]);
    let gate = CancelGate::new_clock_free();
    let response = model.prepare_harmonic(30.0,budget(),&gate).unwrap().solve_modal(&[C64::ZERO],&gate).unwrap();
    assert!(response.state.iter().all(|x| *x==C64::ZERO));
    assert_eq!(response.supplied_power_w,0.0);
    assert_eq!(response.mean_stored_energy_j,0.0);
    assert_eq!(response.equation_relative_residual,0.0);
}

use super::super::harmonic_ports::{FrameHarmonicPortBudget,FrameHarmonicPortError};
fn port_budget() -> FrameHarmonicPortBudget {
    FrameHarmonicPortBudget {max_ports:8,max_projection_terms:4096,max_mobility_entries:64,
        power_absolute_tolerance_w:1e-12,power_relative_tolerance:1e-8}
}
fn node_columns() -> Vec<Vec<f64>> {
    [4,5,7,8].into_iter().map(|dof| {let mut c=vec![0.0;9];c[dof]=1.0;c}).collect()
}

#[test]
fn physical_force_and_moment_multiports_are_reciprocal_and_keep_coherent_power() {
    let frame=frame(true,1000.0,2e6,500.0);
    let model=model(&frame,4,vec![
        FrameFluidRegion {element:1,span:[0.0,0.37],fluid:fluid(1000.0,0.01)},
        FrameFluidRegion {element:1,span:[0.37,1.0],fluid:fluid(3.0,0.002)},
    ]);
    let gate=CancelGate::new_clock_free();
    let columns=node_columns();let ports=model.harmonic_ports(&columns,port_budget(),&gate).unwrap();
    assert!(core::ptr::eq(ports.model(),&model));
    let efforts=[C64::new(0.01,0.02),C64::new(-0.007,0.002),C64::new(0.008,-0.01),C64::new(-0.009,-0.003)];
    for omega in [1.0,30.0,100.0] {
        let prepared=model.prepare_harmonic(omega,budget(),&gate).unwrap();
        let y=ports.mobility(&prepared,&gate).unwrap();
        assert_eq!(y.port_count,4);
        for i in 0..4 {for j in 0..4 {near(y.entries[i*4+j],y.entries[j*4+i],1e-10);}}
        let response=ports.solve(&prepared,&efforts,&gate).unwrap();
        for i in 0..4 {
            let expected=(0..4).fold(C64::ZERO,|s,j|s+y.entries[i*4+j]*efforts[j]);
            near(response.port_velocity[i],expected,1e-10);
            near(response.port_displacement[i],response.nodal_displacement[[4,5,7,8][i]],1e-12);
        }
        let independent:f64=response.nodal_force.iter().zip(&response.nodal_velocity)
            .map(|(f,v)|0.5*(*f*v.conj()).re).sum();
        near(C64::new(independent,0.0),C64::new(response.mechanical.supplied_power_w,0.0),1e-12);
        assert!(response.mechanical.supplied_power_w>=0.0);
        assert!(response.mapping_power_residual_w.abs()<1e-11);
        if omega==30.0 {assert!(response.port_power_w.iter().any(|p|*p<0.0));}
    }
}

#[test]
fn equal_opposite_drives_cancel_before_fluid_and_solid_heat_are_computed() {
    let frame=frame(false,1000.0,2e6,500.0);
    let model=model(&frame,1,vec![FrameFluidRegion {element:0,span:[0.0,1.0],fluid:fluid(400.0,0.01)}]);
    let gate=CancelGate::new_clock_free();
    let ports=model.harmonic_transverse_ports(&[(0,0.7),(0,0.7)],port_budget(),&gate).unwrap();
    let prepared=model.prepare_harmonic(30.0,budget(),&gate).unwrap();
    let u=C64::new(0.01,0.005);
    let cancelled=ports.solve(&prepared,&[u,u.scale(-1.0)],&gate).unwrap();
    assert!(cancelled.mechanical.state.iter().all(|x|*x==C64::ZERO));
    assert_eq!(cancelled.mechanical.supplied_power_w,0.0);
    assert!(cancelled.mechanical.fluid_region_loss_w.iter().all(|p|*p==0.0));
    let one=ports.solve(&prepared,&[u,C64::ZERO],&gate).unwrap();
    let both=ports.solve(&prepared,&[u,u],&gate).unwrap();
    near(both.port_velocity[0],one.port_velocity[0].scale(2.0),1e-12);
    near(C64::new(both.mechanical.supplied_power_w,0.0),C64::new(4.0*one.mechanical.supplied_power_w,0.0),1e-12);
}

#[test]
fn transverse_point_ports_use_original_geometry_and_preserve_rotation_work() {
    let frame=frame(true,1000.0,2e6,500.0);
    let model=model(&frame,4,vec![]);let gate=CancelGate::new_clock_free();
    let points=[(0,0.3),(1,0.6)];
    let ports=model.harmonic_transverse_ports(&points,port_budget(),&gate).unwrap();
    for (column,&(element,xi)) in ports.nodal_columns().iter().zip(&points) {
        let (indices,weights,_)=frame.reference_transverse_trace(element,xi).unwrap();
        let mut expected=vec![0.0;9];
        for (index,weight) in indices.into_iter().zip(weights) {expected[index]+=weight;}
        assert_eq!(column,&expected);
        assert!(indices.into_iter().filter(|i|*i%3==2).any(|i|column[i]!=0.0));
    }
    let prepared=model.prepare_harmonic(30.0,budget(),&gate).unwrap();
    let result=ports.solve(&prepared,&[C64::new(0.01,0.0),C64::new(0.0,0.02)],&gate).unwrap();
    assert!(result.nodal_force.iter().enumerate().filter(|(i,_)|*i%3==2).any(|(_,f)|*f!=C64::ZERO));
    assert!(result.mapping_power_residual_w.abs()<1e-12);
}

#[test]
fn signed_port_scaling_changes_coordinates_not_forces_motion_or_heat() {
    let frame=frame(true,1000.0,2e6,500.0);
    let model=model(&frame,4,vec![FrameFluidRegion {element:1,span:[0.0,1.0],fluid:fluid(400.0,0.01)}]);
    let gate=CancelGate::new_clock_free();let columns=node_columns();
    let scales=[2.0,-3.0,0.5,-0.25];
    let scaled:Vec<_>=columns.iter().zip(scales).map(|(c,s)|c.iter().map(|v|v*s).collect::<Vec<_>>()).collect();
    let a=model.harmonic_ports(&columns,port_budget(),&gate).unwrap();
    let b=model.harmonic_ports(&scaled,port_budget(),&gate).unwrap();
    let prepared=model.prepare_harmonic(30.0,budget(),&gate).unwrap();
    let u=[C64::new(0.01,0.02),C64::new(-0.007,0.002),C64::new(0.008,-0.01),C64::new(-0.009,-0.003)];
    let v:Vec<_>=u.iter().zip(scales).map(|(u,s)|u.scale(1.0/s)).collect();
    let ra=a.solve(&prepared,&u,&gate).unwrap();let rb=b.solve(&prepared,&v,&gate).unwrap();
    for (&a,&b) in ra.mechanical.state.iter().zip(&rb.mechanical.state) {near(a,b,1e-12);}
    near(C64::new(ra.mechanical.supplied_power_w,0.0),C64::new(rb.mechanical.supplied_power_w,0.0),1e-12);
    let ya=a.mobility(&prepared,&gate).unwrap();let yb=b.mobility(&prepared,&gate).unwrap();
    for i in 0..4 {for j in 0..4 {near(yb.entries[i*4+j],ya.entries[i*4+j].scale(scales[i]*scales[j]),1e-10);}}
}

#[test]
fn source_mismatch_invalid_ports_and_cancellation_refuse_without_corrupting_factor() {
    let frame=frame(true,1000.0,2e6,500.0);
    let a=model(&frame,4,vec![]);let b=model(&frame,4,vec![]);let gate=CancelGate::new_clock_free();
    let columns=node_columns();let ports=a.harmonic_ports(&columns,port_budget(),&gate).unwrap();
    let prepared=a.prepare_harmonic(30.0,budget(),&gate).unwrap();
    let foreign=b.prepare_harmonic(30.0,budget(),&gate).unwrap();
    assert!(ports.mobility(&foreign,&gate).is_err());
    assert!(a.harmonic_ports(&[],port_budget(),&gate).is_err());
    assert!(a.harmonic_ports(&[vec![0.0;8]],port_budget(),&gate).is_err());
    assert!(a.harmonic_ports(&[vec![f64::NAN;9]],port_budget(),&gate).is_err());
    assert!(a.harmonic_transverse_ports(&[(0,1.1)],port_budget(),&gate).is_err());
    assert!(a.harmonic_ports(&columns,FrameHarmonicPortBudget {max_mobility_entries:15,..port_budget()},&gate).is_err());
    let original=ports.mobility(&prepared,&gate).unwrap();
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    assert!(matches!(ports.mobility(&prepared,&cancelled),Err(FrameHarmonicPortError::Cancelled)));
    assert_eq!(ports.mobility(&prepared,&gate).unwrap().entries,original.entries);
}
