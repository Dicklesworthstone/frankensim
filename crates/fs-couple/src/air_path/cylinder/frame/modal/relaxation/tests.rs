//! Original pHS and independent layer/strain matrices, not a fitted decay curve.
use super::*;
use crate::air_path::cylinder::{CylinderFlowResolution, frame::{FrameFluidBinding, FrameFluidBudget}};
use crate::frame_vibration::FrameVibrationBudget;
use crate::modal_acoustic_time::{ModalAcousticState, ModalAcousticTimeBudget};
use crate::render::schedule::force::coupled::ModalCouplingConfig;
use fs_evidence::ValidityDomain;
use fs_la::eigen_complex::lu_complex_with_cancel;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::fluid_transport::FluidTransportState;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::state_point::{MaterialPropertySelection, ScalarAdmissibility,
    ScalarPropertyRequirement, resolve_material_state_point};
use fs_math::c64::C64;
use fs_phs::{PortHamiltonian, QuadraticStorage};
use fs_phs::linear_exact::{LinearSplitBudget, PreparedExactZoh};
use fs_qty::{Density, Pressure, Time};
use fs_solid::frame2d::{Frame2d, FrameMassFormulation, Geometry};

const RATE: u32 = 256;
const RADII: [f64; 2] = [0.06, 0.1];
const DELTAS: [f64; 2] = [600.0, 900.0];
const TIMES: [f64; 2] = [0.025, 0.12];
fn frame(delta_scale: f64, tau_scale: f64, one_dof: bool) -> Frame2d {
    let mut layers = Vec::new();
    for layer in 0..2 {
        let e = [1000.0, 1500.0][layer];
        let values = [("density", Density::DIMS, [3.0, 5.0][layer]),
            ("young_modulus", Pressure::DIMS, e), ("equilibrium", Pressure::DIMS, e),
            ("increment", Pressure::DIMS, DELTAS[layer] * delta_scale),
            ("time", Time::DIMS, TIMES[layer] * tau_scale),
            ("zero", Pressure::DIMS, 0.0), ("zero_time", Time::DIMS, 0.2)];
        let mut claims = ClaimSet::new();
        let mut requirements = Vec::new();
        for (name, dims, value) in values {
            let mut validity = ValidityDomain::unconstrained().with("T", 290.0, 310.0);
            if dims == Time::DIMS { validity = validity.with("omega", 0.0, 10000.0); }
            claims.insert_claim(PropertyClaim { key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims }, validity,
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![], provenance: Provenance { source: "manufactured layer E(t)".into(),
                    license: "test-only".into(), artifact: None } }).unwrap();
            requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
        }
        let card = MaterialCard::assemble(MaterialStateId { chemistry: format!("synthetic layer {layer}"),
            phase: "solid".into(), process: "test only".into(), revision: 0 }, claims, vec![]).unwrap();
        let query = QueryPoint::new().with("T", 300.0).unwrap().with("omega", 1.0).unwrap();
        layers.push(ConcentricMaterialLayer { outer_radius_m: RADII[layer],
            material: resolve_material_state_point(&card, &query, &requirements,
                MaterialPropertySelection::SingleClaimOnly).unwrap() });
    }
    let section = ConcentricElasticSection::try_new(0.0, layers).unwrap();
    let mut frame = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    frame.add_node(0.0, 0.0); frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    if one_dof { frame.fix(1, [true, false, true]); }
    // No invented eta property is required just to use a Maxwell material.
    frame.add_material_element(0, 1, 3, section).unwrap();
    frame
}
fn model(frame: &Frame2d, wet: bool, modes: usize) -> MaterialFluidFrame<'_, FluidTransportState> {
    MaterialFluidFrame::new(frame, modes, if wet { vec![FrameFluidBinding { element: 0,
        fluid: FluidTransportState::try_new(2.0, 4000.0, 0.5, 0.02).unwrap() }] } else { vec![] },
        CylinderFlowResolution { minimum_wavenumber_radius: 0.01, maximum_wavenumber_radius: 2.0, intervals: 2 },
        FrameVibrationBudget { max_frame_dofs: 6, max_dense_entries: 36 },
        FrameFluidBudget { max_states: 32, max_dense_entries: 1024,
            max_projection_terms: 4096, mass_inverse_tolerance: 1e-10 }, &CancelGate::new_clock_free()).unwrap()
}
fn config() -> ModalCouplingConfig {
    ModalCouplingConfig { max_modes: 32, max_connections: 128, max_setup_terms: 2_000_000,
        nyquist_guard_fraction: 0.9, maximum_total_energy_j: 10.0,
        maximum_abs_pressure_pa: 100.0, maximum_abs_connection_force_n: 1000.0,
        solve_relative_tolerance: 1e-10, energy_absolute_tolerance_j: 1e-10, energy_relative_tolerance: 1e-9 }
}
fn budget() -> FrameRelaxationBudget {
    FrameRelaxationBudget { max_members: 2, max_source_arms: 8,
        max_memory_coordinates: 24, max_setup_terms: 100_000 }
}
fn selections() -> [FrameLongitudinalRelaxation; 1] {
    [FrameLongitudinalRelaxation { element: 0, layers: vec![LayerRelaxationProperties {
        equilibrium_modulus: "equilibrium".into(), branches: vec![("increment".into(), "time".into()),
            ("zero".into(), "zero_time".into())] }; 2] }]
}
fn initial(model: &MaterialFluidFrame<'_, FluidTransportState>) -> Vec<f64> {
    let motion: Vec<_> = (0..model.basis().modes().len()).map(|i| ModalAcousticState {
        displacement_m_sqrt_kg: 0.0003 * (i + 1) as f64,
        velocity_m_sqrt_kg_per_s: -0.001 * (i + 1) as f64 }).collect();
    model.state_from_motion(&motion, &vec![0.0001; model.memory_masses_kg().len()]).unwrap()
}
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance * (1.0 + a.abs() + b.abs()), "{a:e} != {b:e}");
}

// Independent assembly: add scalar E_j(t) to the ORIGINAL canonical pHS
// displacement/momentum matrix, not the adapter's normalized native branches.
fn reference(model: &MaterialFluidFrame<'_, FluidTransportState>, delta_scale: f64, tau_scale: f64)
    -> (PortHamiltonian, Vec<(Vec<f64>, f64, f64)>) {
    let m = model.basis().modes().len();
    let old_n = model.system().state_dim();
    let n = old_n + 6;
    let (old_j, old_r, old_g) = model.system().structure();
    let old_q = model.system().storage().exact_quadratic_matrix().unwrap();
    let (mut j, mut r, mut q, mut g) = (vec![0.0; n*n], vec![0.0; n*n], vec![0.0; n*n], vec![0.0; n*m]);
    for i in 0..old_n {
        for k in 0..old_n { j[i*n+k] = old_j[i*old_n+k]; r[i*n+k] = old_r[i*old_n+k]; q[i*n+k] = old_q[i*old_n+k]; }
        for k in 0..m { g[i*m+k] = old_g[i*m+k]; }
    }
    let maps = [[-1.0, 0.0, 0.0, 1.0, 0.0, 0.0],
        [0.0, 2.0, 1.0, 0.0, -2.0, 1.0], [0.0, 0.0, 1.0, 0.0, 0.0, -1.0]];
    let mut arms = Vec::new();
    for layer in 0..2 {
        let inner = if layer == 0 { 0.0_f64 } else { RADII[layer - 1] };
        let area = core::f64::consts::PI * (RADII[layer].powi(2) - inner.powi(2));
        let moment = core::f64::consts::PI * (RADII[layer].powi(4) - inner.powi(4)) / 4.0;
        for (weights, geometry) in maps.iter().zip([area, 3.0 * moment, moment]) {
            let b: Vec<_> = model.basis().modes().iter().map(|mode|
                weights.iter().zip(&mode.displacement).map(|(w, u)| w*u).sum::<f64>()).collect();
            let stiffness = DELTAS[layer] * delta_scale * geometry;
            let tau = TIMES[layer] * tau_scale;
            let index = old_n + arms.len();
            for a in 0..m {
                for z in 0..m { q[a*n+z] += stiffness*b[a]*b[z]; }
                q[a*n+index] = -stiffness.sqrt()*b[a];
                q[index*n+a] = q[a*n+index];
            }
            q[index*n+index] = 1.0;
            r[index*n+index] = 1.0 / tau;
            arms.push((b, stiffness, tau));
        }
    }
    (PortHamiltonian::new(n, m, j, r, g, Box::new(QuadraticStorage::new(q, n).unwrap())).unwrap(), arms)
}

#[test]
fn sourced_layer_memories_match_full_canonical_trajectory_and_independent_heat() {
    for wet in [false, true] { for relaxed in [false, true] {
        let frame = frame(1.0, 1.0, false); let model = model(&frame, wet, 3);
        let state = initial(&model); let gate = CancelGate::new_clock_free();
        let (placement, native) = model.coupled_component(&state, RATE,
            ModalAcousticTimeBudget::audible_reference(), config(), &gate).unwrap().into_system(None, &gate).unwrap();
        let (binding, mut native) = placement.with_longitudinal_relaxation(native, &selections(),
            if relaxed { None } else { Some(vec![0.0; 6]) }, budget(), &gate).unwrap();
        assert!(core::ptr::eq(binding.source(), &model));
        assert_eq!(binding.coordinates().len(), 6); // explicit zero arms remain in the source only
        assert_eq!(binding.members()[0].source.branches().len(), 4);
        let (reference, arms) = reference(&model, 1.0, 1.0);
        let mut x = state.clone();
        for (b, k, _) in &arms { x.push(if relaxed { k.sqrt()*b.iter().zip(&state).map(|(b,q)| b*q).sum::<f64>() } else { 0.0 }); }
        let mut exact = PreparedExactZoh::new(&reference, 1.0/f64::from(RATE),
            LinearSplitBudget { max_modes: x.len(), max_omega_dt: 3.0 }).unwrap();
        for step in 0..24 {
            let mut nodal = [0.0; 6]; nodal[3] = 0.001; nodal[4] = if step < 12 { 0.002 } else { -0.001 }; nodal[5] = 0.0001;
            let mut dry_force = vec![0.0; 3];
            model.basis().project_force_into(&nodal, [0.0; 2], &mut dry_force, || false).unwrap();
            let mut force = placement.project_nodal_force(&nodal, &gate).unwrap(); force.resize(native.mode_count(), 0.0);
            let report = native.step_under_gate(&force, &gate).unwrap().clone();
            let mut next = vec![0.0; x.len()]; let mut average = vec![0.0; 3];
            let expected = exact.step_into(&x, &dry_force, &mut next, &mut average).unwrap();
            x = next;
            let mut physical = placement.physical_state(native.components(), &gate).unwrap();
            physical.extend_from_slice(&native.material_memory()[binding.memory_range()]);
            for (&a,&b) in physical.iter().zip(&x) { near(a,b,2e-9); }
            near(native.total_energy_j().unwrap(), reference.hamiltonian(&x), 2e-10);
            near(report.external_work_j, expected.supplied, 2e-10);
            near(report.material_dissipation_j + report.connection_dissipation_j + report.component_dissipation_j,
                expected.dissipated, 2e-10);
            let observed = binding.states(&native, &gate).unwrap();
            near(observed.iter().map(|s|s.stored_energy_j).sum(), native.material_energy_j().unwrap(), 2e-12);
            assert!(observed.iter().all(|s| s.stored_energy_j >= 0.0 && s.dissipated_power_w >= 0.0));
        }
    }}
}

#[test]
fn native_frequency_response_preserves_full_complex_longitudinal_spectrum() {
    for (delta_scale, tau_scale) in [(1.0, 1.0), (2.0, 0.4)] {
        let frame = frame(delta_scale, tau_scale, false); let model = model(&frame, true, 3);
        let gate = CancelGate::new_clock_free(); let state = initial(&model);
        let (placement, native) = model.coupled_component(&state, RATE,
            ModalAcousticTimeBudget::audible_reference(), config(), &gate).unwrap().into_system(None,&gate).unwrap();
        let (_, native) = placement.with_longitudinal_relaxation(native,&selections(),None,budget(),&gate).unwrap();
        let (reference, _) = reference(&model, delta_scale, tau_scale);
        let (j,r,g) = reference.structure(); let q = reference.storage().exact_quadratic_matrix().unwrap(); let n=reference.state_dim();
        let dry_force=[0.002,-0.001,0.0003];
        let chart=FrameFluidCoordinates { source:&model, factor:cholesky(model.modal_mass_matrix(),3).unwrap() };
        let mut force:Vec<_>=chart.project_column(&dry_force).unwrap().into_iter().map(|f| C64::new(f,0.0)).collect();
        force.resize(native.mode_count(),C64::ZERO);
        for omega in [0.5, 5.0, 20.0, 100.0] {
            let actual=native.harmonic_response(omega,&force,1_000_000,&gate).unwrap();
            let mut pencil=vec![C64::ZERO;n*n];let mut rhs=vec![C64::ZERO;n];
            for i in 0..n { for k in 0..n {
                pencil[i*n+k]=C64::new(-(0..n).map(|a|(j[i*n+a]-r[i*n+a])*q[a*n+k]).sum::<f64>(),
                    if i==k {-omega} else {0.0});
            } for k in 0..3 { rhs[i]=rhs[i]+C64::new(g[i*3+k]*dry_force[k],0.0); } }
            let factor=lu_complex_with_cancel(&pencil,n,||false).unwrap();factor.solve_with_cancel(&mut rhs,||false).unwrap();
            for imaginary in [false,true] {
                let part=|z:C64|if imaginary {z.im} else {z.re};
                let normalized:Vec<_>=(0..3).map(|i|ModalAcousticState {
                    displacement_m_sqrt_kg:part(actual.displacement_m_sqrt_kg[i]),
                    velocity_m_sqrt_kg_per_s:part(actual.velocity_m_sqrt_kg_per_s[i]) }).collect();
                let dry=chart.original_states(&normalized).unwrap();
                for i in 0..3 {near(dry[i].displacement_m_sqrt_kg,part(rhs[i]),2e-10);}
                for k in 0..6 {near(part(actual.material_memory[k]),part(rhs[model.system().state_dim()+k]),2e-10);}
            }
            assert!(actual.material_dissipation_w>0.0);
            near(actual.total_dissipation_w(),actual.external_power_w,2e-10);
        }
    }
}

#[test]
fn constrained_positive_arms_keep_initial_energy_and_zero_sources_add_none() {
    let frame=frame(1.0,1.0,true);let model=model(&frame,false,1);let gate=CancelGate::new_clock_free();
    let state=initial(&model);
    let (placement,native)=model.coupled_component(&state,RATE,ModalAcousticTimeBudget::audible_reference(),config(),&gate)
        .unwrap().into_system(None,&gate).unwrap();
    let (binding,mut native)=placement.with_longitudinal_relaxation(native,&selections(),Some(vec![0.01;6]),budget(),&gate).unwrap();
    let old=binding.states(&native,&gate).unwrap();
    assert_eq!(old[0].extension_m,0.0);assert!(old[0].stored_energy_j>0.0);assert!(old[0].dissipated_power_w>0.0);
    let cancelled=CancelGate::new_clock_free();cancelled.request();
    let memory=native.material_memory().to_vec();
    assert!(native.step_under_gate(&[0.0],&cancelled).is_err());assert_eq!(native.material_memory(),memory.as_slice());
    native.step_under_gate(&[0.0],&gate).unwrap();
    near(native.material_memory()[0],memory[0]*(-1.0/f64::from(RATE)/TIMES[0]).exp(),1e-10);
}

#[test]
fn source_missing_roles_duplicate_members_budget_and_warm_binding_refuse() {
    let frame=frame(1.0,1.0,false);let model=model(&frame,true,3);let state=initial(&model);let gate=CancelGate::new_clock_free();
    let build=||model.coupled_component(&state,RATE,ModalAcousticTimeBudget::audible_reference(),config(),&gate)
        .unwrap().into_system(None,&gate).unwrap();
    let (placement,native)=build();
    assert!(placement.with_longitudinal_relaxation(native,&selections(),None,
        FrameRelaxationBudget {max_memory_coordinates:11,..budget()},&gate).is_err());
    let mut missing=selections();missing[0].layers[0].branches[0].1="absent_time".into();
    let (placement,native)=build();
    assert!(matches!(placement.with_longitudinal_relaxation(native,&missing,None,budget(),&gate),Err(FrameRelaxationError::Source(_))));
    let duplicate=vec![selections()[0].clone(),selections()[0].clone()];let (placement,native)=build();
    assert!(placement.with_longitudinal_relaxation(native,&duplicate,None,budget(),&gate).is_err());
    let (placement,mut native)=build();let modes=native.mode_count();native.step_under_gate(&vec![0.0;modes],&gate).unwrap();
    assert!(placement.with_longitudinal_relaxation(native,&selections(),None,budget(),&gate).is_err());
    let (placement,native)=build();
    let (_,native)=placement.with_longitudinal_relaxation(native,&selections(),None,budget(),&gate).unwrap();
    assert!(placement.with_longitudinal_relaxation(native,&selections(),None,budget(),&gate).is_err());
}

#[test]
fn two_source_frames_append_disjoint_material_histories_without_resetting_fluid_or_prefix() {
    use crate::render::schedule::force::coupled::{ModalAttachment, ModalConnection};
    let gate=CancelGate::new_clock_free();
    let frames=[frame(1.0,1.0,false),frame(2.0,0.4,false)];
    let models=[model(&frames[0],true,3),model(&frames[1],true,3)];
    let states=[initial(&models[0]),initial(&models[1])];
    let first=models[0].coupled_component(&states[0],RATE,ModalAcousticTimeBudget::audible_reference(),config(),&gate).unwrap();
    let (a,mut components,mut connections)=first.into_parts(0,0,&gate).unwrap();
    let second=models[1].coupled_component(&states[1],RATE,ModalAcousticTimeBudget::audible_reference(),config(),&gate).unwrap();
    let (b,second_components,second_connections)=second.into_parts(components.len(),connections.len(),&gate).unwrap();
    components.extend(second_components);connections.extend(second_connections);
    let nodal=[0.0,0.0,0.0,0.0,1.0,0.0];
    connections.push(ModalConnection {
        left:ModalAttachment {component:a.structural_component(),shapes:a.project_nodal_force(&nodal,&gate).unwrap()},
        right:ModalAttachment {component:b.structural_component(),shapes:b.project_nodal_force(&nodal,&gate).unwrap()},
        stiffness_n_m:2.0,damping_n_s_m:0.02,rest_extension_m:0.0,
    });
    let native=CoupledModalSystem::new(components,connections,config(),&gate).unwrap();
    let (binding_a,native)=a.with_longitudinal_relaxation(native,&selections(),Some(vec![0.001;6]),budget(),&gate).unwrap();
    let prefix=native.material_memory().to_vec();
    let fluid_before=b.physical_state(native.components(),&gate).unwrap();
    let (binding_b,mut native)=b.with_longitudinal_relaxation(native,&selections(),None,budget(),&gate).unwrap();
    assert_eq!(binding_a.memory_range(),0..6);assert_eq!(binding_b.memory_range(),6..12);
    assert_eq!(&native.material_memory()[..6],prefix.as_slice());
    assert_eq!(b.physical_state(native.components(),&gate).unwrap(),fluid_before);
    assert_eq!(binding_b.members()[0].source.spectra()[0].terms[0].1,TIMES[0]*0.4);
    let modes=native.mode_count();
    for _ in 0..8 { native.step_under_gate(&vec![0.0;modes],&gate).unwrap(); }
    let energy_a:f64=binding_a.states(&native,&gate).unwrap().iter().map(|s|s.stored_energy_j).sum();
    let energy_b:f64=binding_b.states(&native,&gate).unwrap().iter().map(|s|s.stored_energy_j).sum();
    near(energy_a+energy_b,native.material_energy_j().unwrap(),2e-12);
    assert_eq!(native.samples_rendered(),8);
}
