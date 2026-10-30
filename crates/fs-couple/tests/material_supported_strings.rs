//! Shared material-source and original nonlinear mechanics, not measured instrument validation.
#[path = "material_string_nonlinear/support.rs"]
mod support;
#[path = "material_supported_strings/axial_loss.rs"]
mod axial_loss;
#[path = "material_supported_strings/axial_relaxation.rs"]
mod axial_relaxation;
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
        // Both strings share the support, so both take the twist-resolving clock.
        let twisted=source.into_nonlinear(cx,gas(),TORSION_RATE_HZ,0.2,initial(5,0.0),gate).unwrap();
        let first=support::source(7800.0,200e9,2).into_nonlinear(cx,gas(),TORSION_RATE_HZ,0.2,initial(2,0.0),gate).unwrap();
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
        let mut plain=image.clone().prepare(config_at(TORSION_RATE_HZ),gate).unwrap();
        let mut tangent=image.prepare(config_at(TORSION_RATE_HZ),gate).unwrap();
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

#[test]
fn g3_two_plane_material_swaps_reach_the_common_bridge_and_sympathetic_motion() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;
    use fs_couple::render::plate::impact::supported::StringEndMotion;
    with_cx(|cx,gate| {
        let mut histories=Vec::new();
        for (young,pressure) in [(50e9,101325.0),(200e9,101325.0),(50e9,55000.0)] {
            let gas=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,pressure).unwrap();
            let directions=TransverseDirections::try_new([1.0,0.0],[0.6,0.8]).unwrap();
            let sources=vec![
                source(7800.0,young,2).with_transverse_polarizations(directions).unwrap()
                    .into_nonlinear(cx,gas,48_000,0.2,initial(4,0.001),gate).unwrap(),
                source(7800.0,200e9,2).with_transverse_polarizations(directions).unwrap()
                    .into_nonlinear(cx,gas,48_000,0.2,initial(4,0.0),gate).unwrap(),
            ];
            let total_mass=rail().mass_kg+sources.iter()
                .map(|s|s.card().linear_density_kg_m*s.card().length_m).sum::<f64>();
            let image=SupportedMaterialStrings::in_stationary_gas(sources,
                TranslatingSupport {initial_position_m:0.0,initial_velocity_m_s:0.0,..rail()},
                vec![StringEndMotion::Second,StringEndMotion::First],300.0,gate).unwrap();
            assert_eq!(image.mode_count(),9);
            assert!((image.supported().total_mass_kg()-total_mass).abs()<1e-15);
            for law in image.supported().absolute_drag().unwrap() {
                assert_eq!(law.modal_per_s.len(),4);
                assert_eq!(law.modal_per_s[..2],law.modal_per_s[2..]);
            }
            let force_row=image.surface_shape(0,0.2).unwrap();
            let neighbor=image.transverse_shape(1,0.37).unwrap();
            let mut run=image.prepare(config(),gate).unwrap();
            let before=run.system().stored_energy_j();let mut work=0.0;let mut loss=0.0;
            for step in 0..256 {
                let force:Vec<_>=force_row.iter().map(|b|b*if step<16 {0.001} else {0.0}).collect();
                let frame=run.step(&force,gate).unwrap();
                work+=frame.supplied_work_j;loss+=frame.dissipated_energy_j;
                assert!((frame.stored_energy_j+loss-before-work).abs()<1e-8);
            }
            let state=run.system().state();
            let neighbor_motion: f64=neighbor.iter().zip(state.chunks_exact(2)).map(|(b,q)|b*q[0]).sum();
            assert!(neighbor_motion.abs()>1e-10);
            histories.push((state.to_vec(),neighbor_motion));
        }
        assert_ne!(histories[0].0[16],histories[1].0[16]); // Same bridge, changed source material.
        assert_ne!(histories[0].1,histories[1].1); // Original unplayed string remains coupled.
        assert_ne!(histories[0].0,histories[2].0); // Both planes retain actual stationary-gas loss.
    });
}

#[test]
fn g0_two_plane_moving_support_preserves_torsional_surface_work_and_source_order() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::{
        TransverseDirection,TransverseDirections,
    };
    use fs_couple::render::plate::impact::supported::StringEndMotion;
    with_cx(|cx,gate| {
        let original=source(7800.0,200e9,2).with_torsion("G","eta_s",2).unwrap()
            .with_transverse_polarizations(TransverseDirections::try_new([0.8,-0.6],[0.6,0.8]).unwrap()).unwrap()
            .into_nonlinear(cx,gas(),TORSION_RATE_HZ,0.2,initial(6,0.0),gate).unwrap();
        let image=SupportedMaterialStrings::in_stationary_gas(vec![original.clone()],
            TranslatingSupport {initial_position_m:0.0,initial_velocity_m_s:0.0,..rail()},
            vec![StringEndMotion::Second],300.0,gate).unwrap();
        assert_eq!(image.mode_count(),7);
        assert_eq!(image.torsion_ranges(),&[Some(5..7)]);
        assert_eq!(image.bodies()[1].initial,original.bodies()[1].initial);
        assert_eq!(image.bodies()[1].damping_per_s,original.bodies()[1].damping_per_s);
        let normal=image.transverse_shape(0,0.27).unwrap();
        let surface=image.surface_shape(0,0.27).unwrap();
        assert_eq!(normal[..5],surface[..5]);assert_eq!(normal[5..],[0.0;2]);
        assert_eq!(surface[5..],original.surface_shape(0.27).unwrap()[4..]);
        assert_eq!(image.surface_shape(0,1.0).unwrap(),vec![0.0;7]); // Bow is orthogonal to bridge motion.
        let bridge_direction=TransverseDirection::try_new([0.6,0.8]).unwrap();
        let endpoint=image.transverse_shape_along(0,1.0,bridge_direction).unwrap();
        assert_eq!(endpoint,image.support_shape());
        let mut normal_run=image.clone().prepare(config_at(TORSION_RATE_HZ),gate).unwrap();
        let mut surface_run=image.prepare(config_at(TORSION_RATE_HZ),gate).unwrap();
        for _ in 0..64 {
            normal_run.step(&normal.iter().map(|v|v*1e-4).collect::<Vec<_>>(),gate).unwrap();
            surface_run.step(&surface.iter().map(|v|v*1e-4).collect::<Vec<_>>(),gate).unwrap();
        }
        assert!(normal_run.system().state()[10..].iter().all(|v|*v==0.0));
        assert!(surface_run.system().state()[10..].iter().any(|v|v.abs()>0.0));
    });
}

#[test]
fn g3_material_axial_waves_exchange_work_with_an_oblique_support_and_preserve_gas_and_twist() {
    use fs_couple::bowed_string::runtime::reciprocal::bank::polarization::TransverseDirections;
    use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::{
        longitudinal::LongitudinalStringResolution,supported::longitudinal::MaterialLongitudinalSupport,
    };
    use fs_couple::render::plate::impact::supported::{StringEndMotion,longitudinal::StringSupportGeometry};
    use fs_couple::modal_acoustic_time::ModalAcousticState;
    const RATE:u32=1_536_000;
    with_cx(|cx,gate| {
        let mut histories=Vec::new();let mut frequencies=Vec::new();
        for (young,pressure) in [(200e9,101325.0),(50e9,101325.0),(200e9,55000.0)] {
            let medium=GasState::try_new(&GasSpec::dry_air_ussa1976(),293.15,pressure).unwrap();
            let directions=TransverseDirections::try_new([0.6,0.8],[1.0,0.0]).unwrap();
            let mut first_initial=initial(3,0.0002);first_initial[2].velocity_m_sqrt_kg_per_s=0.001;
            let first=source(7800.0,young,1).with_torsion("G","eta_s",1).unwrap()
                .with_transverse_polarizations(directions).unwrap()
                .into_nonlinear(cx,medium,RATE,0.2,first_initial,gate).unwrap();
            let second=source(7800.0,200e9,1).with_transverse_polarizations(directions).unwrap()
                .into_nonlinear(cx,medium,RATE,0.2,initial(2,0.0),gate).unwrap();
            let twist=first.bodies()[1].clone();
            let geometry=StringSupportGeometry::try_new([[0.0;3],[0.5,0.0,0.0]],
                [0.0,1.0,0.0],[0.6,0.8,0.0]).unwrap();
            let attachments:Vec<_>=[StringEndMotion::Second,StringEndMotion::First].into_iter().enumerate()
                .map(|(i,moving_end)|MaterialLongitudinalSupport {moving_end,geometry,
                    resolution:LongitudinalStringResolution {axial_modes:vec![1,2],maximum_slope:0.2,maximum_strain:0.01},
                    initial_axial:vec![ModalAcousticState {displacement_m_sqrt_kg:if i==0 {1e-8} else {0.0},
                        velocity_m_sqrt_kg_per_s:0.0},ModalAcousticState::default()],
                }).collect();
            let image=SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves(vec![first,second],
                TranslatingSupport {initial_position_m:0.0,initial_velocity_m_s:0.0,..rail()},attachments,300.0,gate).unwrap();
            assert_eq!(image.mode_count(),10);assert_eq!(image.torsion_ranges(),&[Some(9..10),None]);
            assert_eq!(image.bodies()[1].initial,twist.initial);assert_eq!(image.bodies()[1].damping_per_s,twist.damping_per_s);
            assert_eq!(image.longitudinal_supports().unwrap()[0].geometry,geometry);
            assert_eq!(image.supported().axial_ranges(),&[Some(2..4),Some(6..8)]);
            assert!(image.supported().absolute_drag().unwrap().iter().all(|d|d.modal_per_s.len()==2));
            frequencies.push(image.supported().axial_angular_frequencies(0).unwrap()[0]);
            let end=image.axial_shape(0,1.0).unwrap();let support=image.support_shape();
            for (a,b) in end.iter().zip(&support) {assert!((a-0.6*b).abs()<1e-14);}
            let row=image.surface_shape(0,0.27).unwrap();assert!(row[9]!=0.0);
            let support_index=image.supported().support_coordinate();
            let mut run=image.prepare(config_at(RATE),gate).unwrap();let before=run.system().stored_energy_j();
            let (mut work,mut loss,mut bridge,mut neighbor)=(0.0,0.0,0.0_f64,0.0_f64);
            for step in 0..256 {
                let force:Vec<_>=row.iter().map(|b|b*if step<16 {1e-4} else {0.0}).collect();
                let frame=run.step(&force,gate).unwrap();work+=frame.supplied_work_j;loss+=frame.dissipated_energy_j;
                assert!((frame.stored_energy_j+loss-before-work).abs()<1e-9);
                let state=run.system().state();bridge=bridge.max(state[2*support_index].abs());
                neighbor=neighbor.max(state[12].abs());
            }
            assert!(bridge>1e-14 && neighbor>1e-16,"material axial force must reach the support and neighbor");
            histories.push(run.system().state().to_vec());
        }
        assert!((frequencies[1]/frequencies[0]-0.5).abs()<1e-14);
        assert_eq!(frequencies[0],frequencies[2]); // Gas is not invented axial solid inertia.
        assert_ne!(histories[0][16],histories[1][16]);
        assert_ne!(histories[0],histories[2]);
    });
}
