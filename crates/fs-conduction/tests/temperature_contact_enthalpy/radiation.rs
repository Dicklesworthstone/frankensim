//! Contact and radiation must pass one complete enthalpy endpoint audit.
use super::*;
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, SurfaceEmissivity, ThermalBc,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, STEFAN_BOLTZMANN_W_M2_K4};

fn patch(ambient:f64)->AmbientRadiationPatch {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar{value:0.2,dims:EMISSIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:Provenance{source:"synthetic contact/radiation test".into(),
            license:"internal-test-use".into(),artifact:None},
    }).unwrap();
    let card=fs_matdb::MaterialCard::assemble(MaterialStateId{chemistry:"synthetic surface".into(),
        phase:"solid".into(),process:"fixture".into(),revision:0},claims,vec![]).unwrap();
    AmbientRadiationPatch::new("exterior",SurfaceEmissivity::from_card("exterior",&card,300.,
        SelectionPolicy::SingleClaimOnly).unwrap(),ambient).unwrap()
}

// Independent P1 roots on the nonuniform trace [0,1,3], not a lumped-body
// approximation or integral of pointwise T^4. Same reference-mass h equations
// as the parent's oracle, with the full area-mean Robin radiation residual.
const FINAL_H: [[f64; 8]; 2] = [
    [258.3924181905082, 267.36664447536435, 302.2872078100072, 77.10669590831374, 372.9703500222133, 358.17586533874453, 360.3844280594463, 446.64392453772973],
    [352.93216001778364, 384.9588787913342, 318.29345193131604, 281.9445452534815, 395.76478120590525, 385.9327085068357, 381.3611586254069, 452.531794595277],
];

#[test]
fn g1_radiative_contact_retains_latent_history_and_signed_external_heat() {
    context(&CancelGate::new_clock_free(),|cx| {
        let mesh=mesh();
        let bc=ThermalBoundaryBuilder::new(&mesh).region("exterior",
            |f|f.centroid[2]<0. && f.centroid[1]==0.,ThermalBc::neumann(0.).unwrap()).unwrap()
            .adiabatic_remainder().finish().unwrap();
        let card=card();let contacts=TemperatureContacts::new(cx,&mesh,&bc,vec![],vec![law(&mesh,&card)],1).unwrap();
        let cold=chart(false);let hot=chart(true);
        let engine=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&[
            ReferenceEnthalpyMaterial{curve:&cold,reference_density_kg_m3:1.},
            ReferenceEnthalpyMaterial{curve:&hot,reference_density_kg_m3:2.},
        ],&[0,1],EnthalpyBudget{max_vertices:8,max_elements:2}).unwrap();
        let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.6),(380.,3.6)]).unwrap());
        let source=ScalarField::Nodal(vec![0.,0.,0.,0.,60.,60.,60.,60.]);
        let problem=ConductionProblem{mesh:&mesh,boundary:&bc,material:&material,element_materials:None,source:&source};
        let initial=vec![80.,80.,80.,80.,450.,450.,450.,450.];
        for (case,ambient) in [300.,340.].into_iter().enumerate() {
            let patches=[patch(ambient)];let mut old=initial.clone();let mut external=0.;
            for dt in [0.15,0.35,0.5] {
                let before=old.clone();
                let next=contacts.advance_enthalpy_with_ambient_radiation(cx,&engine,problem,&old,dt,
                    config(),coupling(),&patches,AmbientRadiationConfig::default()).unwrap();
                assert_eq!(old,before);
                assert!(next.physical.residual_norm_j<=coupling().residual_tolerance_j);
                let t=&next.conduction.temperature;
                let mean=(t[0]+t[1]+t[3])/3.;
                let q=0.2*STEFAN_BOLTZMANN_W_M2_K4*0.5*(mean.powi(4)-ambient.powi(4));
                assert_eq!(q.is_sign_positive(),ambient==300.);
                close(next.physical.radiation_out_w,q,1e-8);
                close(next.physical.convective_out_w,0.,1e-12);
                close(next.physical.stored_energy_change_j,dt*(10.-q),2e-8);
                close(next.physical.energy_residual_j,0.,1e-8);
                let report=next.radiation.as_ref().unwrap();
                assert!(next.combined_boundary.is_some());
                assert!(next.radiation_iterations>report.iterations);
                assert!(next.nonlinear_iterations>=report.solid_iterations);
                assert!(next.krylov_iterations>=report.krylov_iterations);
                let contact=&next.physical_contacts.states()[0];
                close(contact.resistance.value_m2_k_per_w(),0.1+0.006*(contact.interface_temperature_k-280.),1e-14);
                external+=dt*(10.-q);old=next.conduction.specific_enthalpy_j_kg;
            }
            for (&h,&expected) in old.iter().zip(&FINAL_H[case]) {close(h,expected,4e-6);}
            let stored:f64=old.iter().zip(&initial).zip(engine.reference_nodal_masses_kg())
                .map(|((&h,&h0),&m)|m*(h-h0)).sum();
            close(stored,external,5e-8);
            let before=initial.clone();
            assert!(contacts.advance_enthalpy_with_ambient_radiation(cx,&engine,problem,&initial,0.15,
                config(),coupling(),&patches,AmbientRadiationConfig{max_iterations:1,..Default::default()}).is_err());
            assert_eq!(initial,before);
            assert!(contacts.advance_enthalpy_with_ambient_radiation(cx,&engine,problem,&initial,0.15,
                config(),coupling(),&[],AmbientRadiationConfig::default()).is_err());
        }
    });
}
