//! Shared material-source and original nonlinear mechanics, not measured instrument validation.
#[path = "material_string_nonlinear/support.rs"]
mod support;
use support::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::{MaterialStringSource,
    nonlinear::{NonlinearMaterialString, supported::SupportedMaterialStrings}};
use fs_couple::render::plate::impact::{ImpactSystem, supported::{SupportedStrings, TranslatingSupport},
    string::StringStretching};
use fs_couple::string_specimen::{StringPrestress, KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::gas::{GasSpec, GasState};
use fs_exec::{Cx, CancelGate};

fn rail() -> TranslatingSupport {
    TranslatingSupport { mass_kg:0.01, stiffness_n_m:15.0, damping_n_s_m:0.03,
        initial_position_m:0.0001, initial_velocity_m_s:0.01,
        maximum_travel_m:0.02, maximum_slope:0.15 }
}
fn strings(cx:&Cx<'_>, gate:&CancelGate, second_e:f64, gas:GasState) -> Vec<NonlinearMaterialString> {
    vec![source(7800.0,200e9,2).into_nonlinear(cx,gas,48_000,0.2,initial(2,0.0),gate).unwrap(),
        source(7800.0,second_e,3).into_nonlinear(cx,gas,48_000,0.2,initial(3,0.003),gate).unwrap()]
}

#[test]
fn original_materials_and_per_string_laws_match_direct_shared_support_dynamics() {
    with_cx(|cx,gate| {
        let sources=strings(cx,gate,50e9,gas());
        let wires:Vec<_>=sources.iter().map(|s|s.wire_span().clone()).collect();
        let laws:Vec<_>=sources.iter().map(|s|Some(StringStretching {
            axial_rigidity_n:s.axial_rigidity_n(), maximum_slope:0.2 })).collect();
        let relative:Vec<_>=sources.iter().flat_map(|s|s.bodies()[0].initial.iter().copied()).collect();
        let reference=SupportedStrings::with_individual_stretching(&wires,&laws,rail()).unwrap();
        let body=reference.body_with_relative_states(&relative).unwrap();
        let mut direct=ImpactSystem::new(vec![body],vec![],vec![],vec![],config()).unwrap().prepare_analytic().unwrap();
        let image=SupportedMaterialStrings::in_comoving_gas(sources.clone(),rail(),gate).unwrap();
        assert_eq!(image.bodies().len(),1);assert_eq!(image.mode_count(),6);
        assert_eq!(image.supported().total_mass_kg(),reference.total_mass_kg());
        for (a,b) in sources.iter().zip(image.sources()) {
            assert_eq!(a.axial_rigidity_n(),b.axial_rigidity_n());
            assert_eq!(a.card().linear_density_kg_m,b.card().linear_density_kg_m);
            assert_eq!(a.wire_span().damping_per_s,b.wire_span().damping_per_s);
            assert_eq!(a.bodies()[0].initial,b.bodies()[0].initial);
            let MaterialStringSource::Homogeneous(a)=a.source() else {unreachable!()};
            let MaterialStringSource::Homogeneous(b)=b.source() else {unreachable!()};
            assert_eq!(a.material(),b.material());
        }
        let row=image.surface_shape(1,0.37).unwrap();
        let mut run=image.prepare(config(),gate).unwrap();
        let before=run.system().stored_energy_j();let mut work=0.0;let mut loss=0.0;
        for tick in 0..256 {
            let forces:Vec<_>=row.iter().map(|b|b*if tick<17 {0.001} else {0.0}).collect();
            let actual=run.step(&forces,gate).unwrap();let expected=direct.step(&forces,gate).unwrap();
            assert_eq!(actual.stored_energy_j,expected.stored_energy_j);
            assert_eq!(actual.dissipated_energy_j,expected.dissipated_energy_j);
            assert_eq!(run.system().state(),direct.state());
            work+=actual.supplied_work_j;loss+=actual.dissipated_energy_j;
            assert!((actual.stored_energy_j+loss-before-work).abs()<1e-8);
        }
        assert_eq!(run.system().samples(),256);
    });
}

#[test]
fn layered_torsion_keeps_its_original_basis_and_physical_force_projection() {
    with_cx(|cx,gate| {
        let section=ConcentricElasticSection::try_new(0.0,vec![
            ConcentricMaterialLayer {outer_radius_m:0.0002,material:material(7800.0,200e9)},
            ConcentricMaterialLayer {outer_radius_m:0.0003,material:material(15600.0,50e9)},
        ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
        let source=MaterialStringSource::Concentric {template:template(3),section,
            prestress:StringPrestress::FixedTension(20.0)}.with_torsion("G","eta_s",2).unwrap();
        let twisted=source.into_nonlinear(cx,gas(),48_000,0.2,initial(5,0.0),gate).unwrap();
        let first=support::source(7800.0,200e9,2).into_nonlinear(cx,gas(),48_000,0.2,initial(2,0.0),gate).unwrap();
        let image=SupportedMaterialStrings::in_comoving_gas(vec![first,twisted.clone()],
            TranslatingSupport {initial_position_m:0.0,initial_velocity_m_s:0.0,..rail()},gate).unwrap();
        assert_eq!(image.mode_count(),8);assert_eq!(image.bodies().len(),2);
        assert_eq!(image.torsion_ranges(),&[None,Some(6..8)]);
        assert_eq!(image.bodies()[1].initial,twisted.bodies()[1].initial);
        assert_eq!(image.bodies()[1].damping_per_s,twisted.bodies()[1].damping_per_s);
        let normal=image.transverse_shape(1,0.27).unwrap();let surface=image.surface_shape(1,0.27).unwrap();
        assert_eq!(normal[..6],surface[..6]);assert_eq!(normal[6..],[0.0;2]);
        assert_eq!(surface[6..],twisted.surface_shape(0.27).unwrap()[3..]);
        for end in [0.0,1.0] {assert_eq!(image.surface_shape(1,end).unwrap(),image.support_shape());}
        // Two actual runs distinguish normal force (no torsion) from tangential traction.
        let mut plain=image.clone().prepare(config(),gate).unwrap();
        let mut tangent=image.prepare(config(),gate).unwrap();
        for _ in 0..48 {
            plain.step(&normal.iter().map(|b|b*1e-4).collect::<Vec<_>>(),gate).unwrap();
            tangent.step(&surface.iter().map(|b|b*1e-4).collect::<Vec<_>>(),gate).unwrap();
        }
        assert!(plain.system().state()[12..].iter().all(|q|*q==0.0));
        assert!(tangent.system().state()[12..].iter().any(|q|q.abs()>0.0));
    });
}

#[test]
fn material_and_comoving_gas_substitutions_reach_real_support_and_neighbor_motion() {
    with_cx(|cx,gate| {
        let mut histories=Vec::new();
        for (e,p) in [(50e9,101325.0),(200e9,101325.0),(50e9,55000.0)] {
            let g=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,p).unwrap();
            let image=SupportedMaterialStrings::in_comoving_gas(strings(cx,gate,e,g),rail(),gate).unwrap();
            assert_eq!(image.sources()[1].card().tension_n,20.0); // Never retuned to a note.
            let mut run=image.prepare(config(),gate).unwrap();
            for _ in 0..256 {run.step(&[0.0;6],gate).unwrap();}
            histories.push(run.system().state().to_vec());
        }
        assert_ne!(histories[0][0],histories[1][0]); // Material change reaches the other string.
        assert_ne!(histories[0][10],histories[1][10]); // And the common rail.
        assert_ne!(histories[0],histories[2]); // Original gas-derived modal loss is not discarded.
    });
}

#[test]
fn conflicting_frames_clocks_invalid_ports_and_refused_steps_do_not_become_partial_physics() {
    with_cx(|cx,gate| {
        assert!(SupportedMaterialStrings::in_comoving_gas(vec![],rail(),gate).is_err());
        let mut mixed=strings(cx,gate,50e9,gas());
        mixed[1]=source(7800.0,50e9,3).into_nonlinear(cx,gas(),96_000,0.2,initial(3,0.0),gate).unwrap();
        assert!(SupportedMaterialStrings::in_comoving_gas(mixed,rail(),gate).is_err());
        let mut mixed=strings(cx,gate,50e9,gas());
        let g=GasState::try_new(&GasSpec::dry_air_ussa1976(),310.0,101325.0).unwrap();
        mixed[1]=source(7800.0,50e9,3).into_nonlinear(cx,g,48_000,0.2,initial(3,0.0),gate).unwrap();
        assert!(SupportedMaterialStrings::in_comoving_gas(mixed,rail(),gate).is_err());
        let source=strings(cx,gate,50e9,gas());
        let stop=CancelGate::new_clock_free();stop.request();
        assert!(SupportedMaterialStrings::in_comoving_gas(source.clone(),rail(),&stop).is_err());
        let image=SupportedMaterialStrings::in_comoving_gas(source,rail(),gate).unwrap();
        assert!(image.surface_shape(2,0.5).is_err());assert!(image.transverse_shape(0,f64::NAN).is_err());
        let mut wrong=config();wrong.dt_s=1.0/96_000.0;
        assert!(image.clone().prepare(wrong,gate).is_err());
        let mut run=image.prepare(config(),gate).unwrap();
        let before=run.system().state().to_vec();let energy=run.system().stored_energy_j();
        assert!(run.step(&[f64::NAN;6],gate).is_err());
        assert!(run.step(&[0.0;6],&stop).is_err());
        assert_eq!(run.system().state(),before);assert_eq!(run.system().stored_energy_j(),energy);
        assert_eq!(run.system().samples(),0);
        run.step(&[0.0;6],gate).unwrap();assert_eq!(run.system().samples(),1);
    });
}
