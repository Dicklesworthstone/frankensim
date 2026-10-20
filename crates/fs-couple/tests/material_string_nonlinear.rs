//! Source-to-nonlinear dynamics tests, not measured instrument validation.
#[path = "material_string_nonlinear/polarization.rs"]
mod polarization;
#[path = "material_string_nonlinear/pressure.rs"]
mod pressure;
#[path = "material_string_nonlinear/support.rs"]
mod support;
#[path = "material_string_nonlinear/bank.rs"]
mod bank;
use support::*;
use fs_couple::bowed_string::BowedStringCard;
use fs_couple::bowed_string::runtime::reciprocal::bank::MaterialStringSource;
use fs_couple::render::plate::impact::{BodyPotential, ImpactSystem};
use fs_couple::modal_acoustic_time::ModalAcousticState;
use fs_couple::render::plate::impact::linear::wire::WireSpan;
use fs_couple::render::plate::impact::string::StringStretching;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_couple::string_specimen::{StringPrestress, KELVIN_VOIGT_BENDING_VISCOSITY_PROPERTY as ETA};

#[test]
fn material_source_drives_the_original_stretching_potential_and_exact_initial_inventory() {
    with_cx(|cx, gate| {
        let src = source(7800.0, 200e9, 3);
        let MaterialStringSource::Homogeneous(s) = &src else { unreachable!() };
        let original = s.clone();
        let card = BowedStringCard::from_material_specimen_with_gas(cx, s, &gas(), 48_000).unwrap();
        let state = initial(3, 0.003);
        let image = src.into_nonlinear(cx, gas(), 48_000, 0.2, state.clone(), gate).unwrap();
        assert_eq!(image.bodies().len(), 1);
        assert_eq!(image.bodies()[0].initial, state);
        assert_eq!(image.axial_rigidity_n(), original.string().axial_stiffness_n);
        let MaterialStringSource::Homogeneous(retained) = image.source() else { unreachable!() };
        assert_eq!(retained.material().identity(), original.material().identity());
        let wire = WireSpan { endpoints_m: [[0.0,0.0],[card.length_m,0.0]], linear_density_kg_m: card.linear_density_kg_m,
            tension_n: card.tension_n, bending_stiffness_n_m2: card.bending_stiffness_n_m2,
            damping_per_s: (0..3).map(|i| 2.0*card.zetas[i]*card.mode_omega_rad_s(i+1)).collect() };
        let reference = wire.stretching_body(state, StringStretching { axial_rigidity_n: original.string().axial_stiffness_n,
            maximum_slope: 0.2 }).unwrap();
        let mut a = ImpactSystem::new(image.bodies().to_vec(), vec![], vec![], vec![], config()).unwrap().prepare_analytic().unwrap();
        let mut b = ImpactSystem::new(vec![reference], vec![], vec![], vec![], config()).unwrap().prepare_analytic().unwrap();
        for _ in 0..96 {
            let fa=a.step(&[0.0;3],gate).unwrap();let fb=b.step(&[0.0;3],gate).unwrap();
            assert_eq!(fa.stored_energy_j,fb.stored_energy_j);
            assert_eq!(fa.dissipated_energy_j,fb.dissipated_energy_j);
            assert_eq!(a.state(),b.state());
        }
    });
}

#[test]
fn material_ea_sets_geometric_extension_energy_and_actual_amplitude_dependent_dynamics() {
    with_cx(|cx, gate| {
        let low = source(7800.0, 50e9, 3).into_nonlinear(cx, gas(), 48_000, 0.2, initial(3,0.003),gate).unwrap();
        let high = source(7800.0, 200e9, 3).into_nonlinear(cx, gas(), 48_000, 0.2, initial(3,0.003),gate).unwrap();
        assert_eq!(low.card().linear_density_kg_m,high.card().linear_density_kg_m);
        assert_eq!(low.card().tension_n,high.card().tension_n);
        let mut runs = Vec::new();
        for image in [&low,&high] {
            let mut run=ImpactSystem::new(image.bodies().to_vec(),vec![],vec![],vec![],config()).unwrap().prepare_analytic().unwrap();
            let observation=run.string_observation(0).unwrap();
            // Independent continuous sine-field integral, not a copied modal coefficient.
            let strain=0.25*(core::f64::consts::PI*0.003/0.5).powi(2);
            let energy=0.5*image.axial_rigidity_n()*0.5*strain*strain;
            assert!((observation.additional_strain-strain).abs()<1e-16);
            assert!((observation.tension_n-(20.0+image.axial_rigidity_n()*strain)).abs()<1e-12);
            assert!((observation.stretching_energy_j-energy).abs()<1e-15);
            let original_energy=run.stored_energy_j();let mut loss=0.0;
            for _ in 0..256 {
                let f=run.step(&[0.0;3],gate).unwrap();loss+=f.dissipated_energy_j;
                assert!((f.stored_energy_j+loss-original_energy).abs()<1e-8);
            }
            runs.push(run);
        }
        assert!((runs[0].state()[0]-runs[1].state()[0]).abs()>1e-8);
        assert!(runs.iter().all(|s|s.stored_energy_j()>0.0));
    });
}

#[test]
fn bonded_layers_and_torsion_survive_in_one_joint_body_prefix_without_false_normal_torque() {
    with_cx(|cx,gate| {
        let section=ConcentricElasticSection::try_new(0.0,vec![
            ConcentricMaterialLayer {outer_radius_m:0.0002,material:material(7800.0,200e9)},
            ConcentricMaterialLayer {outer_radius_m:0.0003,material:material(15_600.0,50e9)},
        ]).unwrap().with_kelvin_voigt_bending(ETA).unwrap();
        let expected=section.elastic().clone();
        let src=MaterialStringSource::Concentric {template:template(3),section,prestress:StringPrestress::FixedTension(20.0)}
            .with_torsion("G","eta_s",2).unwrap();
        let mut state=vec![ModalAcousticState::default();5];state[3].velocity_m_sqrt_kg_per_s=1e-5;
        let image=src.into_nonlinear(cx,gas(),TORSION_RATE_HZ,0.2,state,gate).unwrap();
        assert_eq!(image.mode_count(),5);assert_eq!(image.bodies().len(),2);
        assert_eq!(image.axial_rigidity_n(),expected.axial_stiffness_n());
        let MaterialStringSource::WithTorsion(t)=image.source() else {panic!("lost twist source")};
        assert_eq!(t.section().shear().elastic(),&expected);
        let normal=image.transverse_shape(0.27).unwrap();let surface=image.surface_shape(0.27).unwrap();
        assert_eq!(normal[..3],surface[..3]);assert_eq!(&normal[3..],&[0.0;2]);
        assert!(surface[3..].iter().all(|x|x.abs()>0.0));
        let mut run=ImpactSystem::new(image.bodies().to_vec(),vec![],vec![],vec![],config_at(TORSION_RATE_HZ)).unwrap().prepare_analytic().unwrap();
        for _ in 0..48 {run.step(&[0.0;5],gate).unwrap();}
        assert!(run.string_observation(0).is_some());assert!(run.string_observation(1).is_none());
        assert!(run.state()[6..].iter().any(|x|*x!=0.0));
        // No geometric axial coupling is invented for a small-angle twist.
        assert!(run.state()[..6].iter().all(|x|*x==0.0));
    });
}

#[test]
fn material_nonlinear_tangent_matches_a_spatial_energy_difference() {
    with_cx(|cx,gate| {
        let image=source(7800.0,200e9,3).into_nonlinear(cx,gas(),48_000,0.2,initial(3,0.0),gate).unwrap();
        let BodyPotential::String(potential)=&image.bodies()[0].potential else {panic!("lost nonlinear storage")};
        let q=[7e-5,-2e-5,1e-5];let d=[0.2,-0.1,0.3];let mut hv=[0.0;3];
        potential.hessian_vector(&q,&d,&mut hv);
        let h=1e-9;let plus=std::array::from_fn::<_,3,_>(|i|q[i]+h*d[i]);
        let minus=std::array::from_fn::<_,3,_>(|i|q[i]-h*d[i]);let mut gp=[0.0;3];let mut gm=[0.0;3];
        potential.gradient(&plus,&mut gp);potential.gradient(&minus,&mut gm);
        for i in 0..3 {assert!((hv[i]-(gp[i]-gm[i])/(2.0*h)).abs()<1e-6*hv[i].abs().max(1.0));}
    });
}

#[test]
fn initial_and_endpoint_slope_refusal_preserve_original_mechanics() {
    with_cx(|cx,gate| {
        assert!(source(7800.0,200e9,3).into_nonlinear(cx,gas(),48_000,0.001,initial(3,0.003),gate).is_err());
        let mut state=initial(3,0.0);state[0].velocity_m_sqrt_kg_per_s=0.02;
        let image=source(7800.0,200e9,3).into_nonlinear(cx,gas(),48_000,1e-7,state,gate).unwrap();
        let mut run=ImpactSystem::new(image.bodies().to_vec(),vec![],vec![],vec![],config()).unwrap().prepare_analytic().unwrap();
        let before=run.state().to_vec();let energy=run.stored_energy_j();
        assert!(run.step(&[0.0;3],gate).is_err());
        assert_eq!(run.samples(),0);assert_eq!(run.state(),before);assert_eq!(run.stored_energy_j(),energy);
    });
}

#[test]
fn invalid_dimension_clock_gas_slope_station_and_cancellation_do_not_make_partial_images() {
    with_cx(|cx,gate| {
        for n in [0,2,4] {assert!(source(7800.0,200e9,3).into_nonlinear(cx,gas(),48_000,0.2,vec![ModalAcousticState::default();n],gate).is_err());}
        for slope in [0.0,-1.0,0.31,f64::NAN] {assert!(source(7800.0,200e9,3).into_nonlinear(cx,gas(),48_000,slope,initial(3,0.0),gate).is_err());}
        let mut bad_gas=gas();bad_gas.density*=2.0;
        assert!(source(7800.0,200e9,3).into_nonlinear(cx,bad_gas,48_000,0.2,initial(3,0.0),gate).is_err());
        assert!(source(7800.0,200e9,3).with_torsion("G","eta_s",2).unwrap()
            .into_nonlinear(cx,gas(),8_000,0.2,initial(5,0.0),gate).is_err());
        let image=source(7800.0,200e9,3).into_nonlinear(cx,gas(),48_000,0.2,initial(3,0.0),gate).unwrap();
        for station in [0.0,1.0,-1.0,f64::NAN] {assert!(image.transverse_shape(station).is_err());assert!(image.surface_shape(station).is_err());}
        let stop=fs_exec::CancelGate::new_clock_free();stop.request();
        assert!(source(7800.0,200e9,3).into_nonlinear(cx,gas(),48_000,0.2,initial(3,0.0),&stop).is_err());
    });
}
