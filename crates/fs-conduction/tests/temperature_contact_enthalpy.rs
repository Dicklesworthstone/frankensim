//! Heterogeneous phase-change contact: independent eight-node P1 references.
#[path = "temperature_contact_enthalpy/radiation.rs"]
mod radiation;
#[path = "temperature_contact_enthalpy/adjoint.rs"]
mod adjoint;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::{ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    ScalarField, ThermalBoundaryBuilder, ThermalInterfaces};
use fs_conduction::interface::temperature::{TemperatureContact, TemperatureContacts, ContactIterationConfig};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig, EnthalpyError};
use fs_conduction::transient::enthalpy::heterogeneous::{
    HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial,
};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;

fn context<T>(gate: &CancelGate, f: impl FnOnce(&Cx<'_>) -> T) -> T {
    ArenaPool::new(ArenaConfig::default()).scope(|arena| f(&Cx::new(gate, arena,
        StreamKey { seed: 20, kernel_id: 820, tile: 0, iteration: 0 },
        Budget::INFINITE, ExecMode::Deterministic)))
}
fn mesh() -> ConductionMesh {
    let points = vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.],
        [0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]];
    ConductionMesh::new(TetComplex::from_tets(8,vec![[0,1,2,3],[4,5,6,7]]),points).unwrap()
}
fn card() -> InterfaceSystemCard {
    let dims=fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,dims),
        value:PropertyValue::Curve { abscissa:"T".into(),abscissa_dims:fs_conduction::TEMPERATURE_DIMS,
            knots:vec![(280.,0.1),(380.,0.7)],dims },
        validity:ValidityDomain::unconstrained().with("T",280.,380.).with("normal_pressure",1e6,1e6),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance {
            source:"synthetic temperature-contact enthalpy reference, not measurements".into(),
            license:"internal-test-use".into(),artifact:None,
        },
    }).unwrap();
    let side=|name:&str| SurfaceSpec {material:MaterialStateId {
        chemistry:name.into(),phase:"solid".into(),process:"fixture".into(),revision:0,
    },texture_frame:"fixed reference normal".into()};
    InterfaceSystemCard::assemble(side("A"),side("B"),SystemContext {
        medium:"dry".into(),third_body:None,environment:"vacuum".into(),history:"unaged".into(),
    },claims,vec![]).unwrap()
}
fn law<'c>(mesh:&ConductionMesh,card:&'c InterfaceSystemCard)->TemperatureContact<'c> {
    let point=QueryPoint::new().with("T",300.).unwrap().with("normal_pressure",1e6).unwrap();
    TemperatureContact::new("joint",ThermalInterfaces::coincident_face_pairs(mesh).unwrap(),
        card,&point,"T",SelectionPolicy::SingleClaimOnly).unwrap()
}
fn chart(hot:bool)->EquilibriumEnthalpyPhaseCurve {
    let knots=if hot { [(-200.,280.,0.),(500.,350.,0.),(600.,350.,1.),(900.,380.,1.)] }
        else { [(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)] };
    EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([if hot {21}else{20};32]),knots.into_iter()
        .map(|(h,t,f)|EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
            liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect()).unwrap()
}
fn config()->EnthalpyStepConfig {
    EnthalpyStepConfig {newton:NewtonKrylovConfig {
        relative_tolerance:1e-12,absolute_tolerance:1e-11,linear_restart:8,max_linear_cycles:16,
        forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()
    },max_newton_iterations:64,energy_tolerance_j:1e-8}
}
fn coupling()->ContactIterationConfig {
    ContactIterationConfig {max_iterations:64,relaxation:0.5,residual_tolerance_j:1e-8}
}
fn close(a:f64,b:f64,tolerance:f64){assert!((a-b).abs()<=tolerance,"{a} vs {b}, tolerance {tolerance}");}

// Independent nonlinear roots, using rho*V/4 nodal reference masses, full P1
// tetrahedral stiffness and the consistent +/- contact blocks (1+delta)/24R.
// The oracle differentiates R(mean of six contact nodes); the product iterates
// frozen contact states instead. No production solve supplies these values.
const EXPECTED_H: [[f64; 8]; 3] = [
    [129.85089908538677, 132.11480278862373, 132.11480278862376, 82.32387512215479, 430.0871357927819, 428.57296418855896, 428.57296418855896, 452.5647459377057],
    [223.46681015518823, 229.07550850956247, 229.0755085095625, 92.1356326009034, 398.4123656797969, 391.1468564655496, 391.1468564655496, 452.4171915014955],
    [288.574542070409, 306.52574214916217, 306.5257421491622, 133.78006982809976, 380.11345751355344, 366.88169717967895, 366.88169717967895, 448.42010002867187],
];

#[test]
fn g1_unequal_mass_contact_melts_spatially_and_matches_independent_enthalpy_roots() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card=card();let contacts=TemperatureContacts::new(cx,&mesh,&bc,vec![],vec![law(&mesh,&card)],1).unwrap();
        let cold=chart(false);let hot=chart(true);
        let engine=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&[
            ReferenceEnthalpyMaterial{curve:&cold,reference_density_kg_m3:1.},
            ReferenceEnthalpyMaterial{curve:&hot,reference_density_kg_m3:2.},
        ],&[0,1],EnthalpyBudget{max_vertices:8,max_elements:2}).unwrap();
        let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(
            vec![(280.,1.6),(380.,3.6)]).unwrap());
        let source=ScalarField::Nodal(vec![0.,0.,0.,0.,60.,60.,60.,60.]);
        let problem=ConductionProblem{mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let initial=vec![80.,80.,80.,80.,450.,450.,450.,450.];
        let mut old=initial.clone();let mut energy=0.;
        for (i,dt) in [0.15,0.35,0.5].into_iter().enumerate() {
            let before=old.clone();
            let next=contacts.advance_enthalpy(cx,&engine,problem,&old,dt,config(),coupling()).unwrap();
            assert_eq!(old,before);
            assert!(next.iterations>1);
            assert!(next.nonlinear_iterations>=next.conduction.newton.iterations);
            assert!(next.krylov_iterations>=next.conduction.newton.history.iter()
                .map(|r|r.linear_iterations).sum::<usize>());
            assert!(next.physical.residual_norm_j<=coupling().residual_tolerance_j);
            close(next.physical.energy_residual_j,0.,1e-8);
            close(next.physical.stored_energy_change_j,10.*dt,1e-8);
            for (&h,&reference) in next.conduction.specific_enthalpy_j_kg.iter().zip(&EXPECTED_H[i]) {
                close(h,reference,3e-6);
            }
            let state=&next.physical_contacts.states()[0];
            close(state.resistance.value_m2_k_per_w(),0.1+0.006*(state.interface_temperature_k-280.),1e-14);
            assert_eq!(state.resistance.card_identity(),card.content_hash());
            if i==0 {
                for v in 0..3 {
                    close(next.conduction.temperature[v],310.,1e-10);
                    assert!(next.conduction.liquid_mass_fraction[v]>0.2
                        && next.conduction.liquid_mass_fraction[v]<0.4);
                }
                assert_eq!(next.conduction.liquid_mass_fraction[3],0.);
            }
            if i==2 {
                assert_eq!(next.conduction.liquid_mass_fraction[0],1.);
                assert!(next.conduction.liquid_mass_fraction[3]>0.3
                    && next.conduction.liquid_mass_fraction[3]<0.4);
                // Contact heat is internal: it never adds to the 10 J input.
                let independent:f64=next.conduction.specific_enthalpy_j_kg.iter()
                    .zip(&initial).zip(engine.reference_nodal_masses_kg())
                    .map(|((&h,&h0),&m)|m*(h-h0)).sum();
                close(independent,10.,3e-8);
            }
            energy+=next.physical.stored_energy_change_j;
            old=next.conduction.specific_enthalpy_j_kg;
        }
        close(energy,10.,3e-8);
    });
}

#[test]
fn g1_isothermal_latent_histories_remain_distinct_and_stationary() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card=card();let contacts=TemperatureContacts::new(cx,&mesh,&bc,vec![],vec![law(&mesh,&card)],1).unwrap();
        let curve=chart(false);
        let engine=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&[
            ReferenceEnthalpyMaterial{curve:&curve,reference_density_kg_m3:1.},
            ReferenceEnthalpyMaterial{curve:&curve,reference_density_kg_m3:2.},
        ],&[0,1],EnthalpyBudget{max_vertices:8,max_elements:2}).unwrap();
        let material=ConductivityModel::isotropic_declared(2.).unwrap();let source=ScalarField::Uniform(0.);
        let problem=ConductionProblem{mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let old=[120.,130.,140.,150.,160.,170.,180.,190.];
        let next=contacts.advance_enthalpy(cx,&engine,problem,&old,0.5,config(),coupling()).unwrap();
        assert_eq!(next.iterations,1);
        for (i,&h) in old.iter().enumerate() {
            close(next.conduction.specific_enthalpy_j_kg[i],h,1e-8);
            close(next.conduction.temperature[i],310.,1e-10);
            close(next.conduction.liquid_mass_fraction[i],(h-100.)/100.,1e-9);
        }
        close(next.physical.stored_energy_change_j,0.,1e-8);
        // Same temperatures, different stored energy: an audit must use h.
        let wrong=old.map(|h|h+1.);
        let audit=engine.evaluate_endpoint(cx,problem,Some(next.physical_contacts.frozen_interfaces()),
            &old,&wrong,0.5,config(),None).unwrap();
        assert!(audit.residual_norm_j>0.1);
        close(audit.energy_residual_j,0.5,1e-10);
    });
}

#[test]
fn g4_contact_enthalpy_limits_domain_and_cancellation_preserve_history() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();let bc=ThermalBoundaryBuilder::new(&mesh).adiabatic_remainder().finish().unwrap();
        let card=card();let contacts=TemperatureContacts::new(cx,&mesh,&bc,vec![],vec![law(&mesh,&card)],1).unwrap();
        let cold=chart(false);let hot=chart(true);
        let engine=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&[
            ReferenceEnthalpyMaterial{curve:&cold,reference_density_kg_m3:1.},
            ReferenceEnthalpyMaterial{curve:&hot,reference_density_kg_m3:2.},
        ],&[0,1],EnthalpyBudget{max_vertices:8,max_elements:2}).unwrap();
        let material=ConductivityModel::isotropic_declared(2.).unwrap();
        let source=ScalarField::Nodal(vec![0.,0.,0.,0.,60.,60.,60.,60.]);
        let problem=ConductionProblem{mesh:&mesh,boundary:&bc,material:&material,
            element_materials:None,source:&source};
        let old=[80.,80.,80.,80.,450.,450.,450.,450.];let before=old;
        let zero=ScalarField::Uniform(0.);
        let negative=[-100.;8];
        let idle=contacts.advance_enthalpy(cx,&engine,ConductionProblem{source:&zero,..problem},
            &negative,0.15,config(),coupling()).unwrap();
        for &h in &idle.conduction.specific_enthalpy_j_kg {close(h,-100.,1e-8);}
        assert!(contacts.advance_enthalpy(cx,&engine,problem,&old,0.15,config(),
            ContactIterationConfig{max_iterations:1,..coupling()}).is_err());
        let mut outside=old;outside[0]=-201.;
        assert!(matches!(contacts.advance_enthalpy(cx,&engine,problem,&outside,0.15,config(),coupling()),
            Err(EnthalpyError::Phase{vertex:0,..})));
        let gate=CancelGate::new_clock_free();gate.request();
        context(&gate,|cancelled|assert!(matches!(contacts.advance_enthalpy(cancelled,&engine,problem,
            &old,0.15,config(),coupling()),Err(EnthalpyError::Conduction(fs_conduction::ConductionError::Cancelled{..})))));
        assert_eq!(old,before);
        let a=contacts.advance_enthalpy(cx,&engine,problem,&old,0.15,config(),coupling()).unwrap();
        let b=contacts.advance_enthalpy(cx,&engine,problem,&old,0.15,config(),coupling()).unwrap();
        assert_eq!(a.conduction.specific_enthalpy_j_kg,b.conduction.specific_enthalpy_j_kg);
        assert_eq!(a.physical_contacts.states(),b.physical_contacts.states());
        assert_eq!(a.iterations,b.iterations);
    });
}
