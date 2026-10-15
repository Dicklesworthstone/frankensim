//! Contact and radiation must pass one complete enthalpy endpoint audit.
use super::*;
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, SurfaceEmissivity, ThermalBc,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, STEFAN_BOLTZMANN_W_M2_K4};

fn patch(ambient:f64)->AmbientRadiationPatch {
    patch_with_emissivity(ambient,0.2)
}
fn patch_with_emissivity(ambient:f64,emissivity:f64)->AmbientRadiationPatch {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar{value:emissivity,dims:EMISSIVITY_DIMS},
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


fn scaled_contact_card(scale:f64)->InterfaceSystemCard {
    let dims=fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS;
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(fs_conduction::AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,dims),
        value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:fs_conduction::TEMPERATURE_DIMS,
            knots:vec![(280.,0.1*scale),(380.,0.7*scale)],dims},
        validity:ValidityDomain::unconstrained().with("T",280.,380.).with("normal_pressure",1e6,1e6),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance{source:"synthetic scaled radiative contact".into(),
            license:"internal-test-use".into(),artifact:None},
    }).unwrap();
    let side=|name:&str| SurfaceSpec{material:MaterialStateId{chemistry:name.into(),phase:"solid".into(),
        process:"fixture".into(),revision:0},texture_frame:"fixed reference normal".into()};
    InterfaceSystemCard::assemble(side("A"),side("B"),SystemContext{medium:"dry".into(),
        third_body:None,environment:"vacuum".into(),history:"unaged".into()},claims,vec![]).unwrap()
}

fn radiative_final(cx:&Cx<'_>,ambient:f64,contact_scale:f64,emissivity:f64)->f64 {
    let mesh=mesh();
    let bc=ThermalBoundaryBuilder::new(&mesh).region("exterior",
        |f|f.centroid[2]<0. && f.centroid[1]==0.,ThermalBc::neumann(0.).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let card=scaled_contact_card(contact_scale);
    let contacts=TemperatureContacts::new(cx,&mesh,&bc,vec![],vec![law(&mesh,&card)],1).unwrap();
    let cold=chart(false);let hot=chart(true);
    let engine=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&[
        ReferenceEnthalpyMaterial{curve:&cold,reference_density_kg_m3:1.},
        ReferenceEnthalpyMaterial{curve:&hot,reference_density_kg_m3:2.},
    ],&[0,1],EnthalpyBudget{max_vertices:8,max_elements:2}).unwrap();
    let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(280.,1.6),(380.,3.6)]).unwrap());
    let source=ScalarField::Nodal(vec![0.,0.,0.,0.,60.,60.,60.,60.]);
    let problem=ConductionProblem{mesh:&mesh,boundary:&bc,material:&material,element_materials:None,source:&source};
    let patches=[patch_with_emissivity(ambient,emissivity)];
    let mut old=vec![80.,80.,80.,80.,450.,450.,450.,450.];
    for dt in [0.15,0.35,0.5] {
        let next=contacts.advance_enthalpy_with_ambient_radiation(cx,&engine,problem,&old,dt,
            config(),coupling(),&patches,AmbientRadiationConfig::default()).unwrap();
        old=next.conduction.specific_enthalpy_j_kg;
    }
    engine.phase_curve_at(0).unwrap().state_at_specific_enthalpy(old[0]).unwrap().temperature_k()
}


fn radiative_gradients(cx:&Cx<'_>,ambient:f64)->[f64;3] {
    let mesh=mesh();
    let bc=ThermalBoundaryBuilder::new(&mesh).region("exterior",
        |f|f.centroid[2]<0. && f.centroid[1]==0.,ThermalBc::neumann(0.).unwrap()).unwrap()
        .adiabatic_remainder().finish().unwrap();
    let card=scaled_contact_card(1.);
    let contacts=TemperatureContacts::new(cx,&mesh,&bc,vec![],vec![law(&mesh,&card)],1).unwrap();
    let cold=chart(false);let hot=chart(true);
    let engine=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&[
        ReferenceEnthalpyMaterial{curve:&cold,reference_density_kg_m3:1.},
        ReferenceEnthalpyMaterial{curve:&hot,reference_density_kg_m3:2.},
    ],&[0,1],EnthalpyBudget{max_vertices:8,max_elements:2}).unwrap();
    let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(
        vec![(280.,1.6),(380.,3.6)]).unwrap());
    let source=ScalarField::Nodal(vec![0.,0.,0.,0.,60.,60.,60.,60.]);
    let problem=ConductionProblem{mesh:&mesh,boundary:&bc,material:&material,element_materials:None,source:&source};
    let patches=[patch_with_emissivity(ambient,0.2)];
    let dts=[0.15,0.35,0.5];
    let mut histories=vec![vec![80.,80.,80.,80.,450.,450.,450.,450.]];
    let mut accepted=Vec::new();
    for &dt in &dts {
        let next=contacts.advance_enthalpy_with_ambient_radiation(cx,&engine,problem,
            histories.last().unwrap(),dt,config(),coupling(),&patches,AmbientRadiationConfig::default()).unwrap();
        histories.push(next.conduction.specific_enthalpy_j_kg.clone());
        accepted.push(next);
    }
    let linear=fs_conduction::LinearConfig{tolerance:1e-10,max_iterations:512,restart:8};
    let mut carry=None;
    let mut totals=[0.;3];
    for index in (0..dts.len()).rev() {
        let step=contacts.linearize_accepted_enthalpy_with_ambient_radiation(
            cx,&engine,problem,&histories[index],dts[index],config(),coupling().residual_tolerance_j,
            &patches,accepted.pop().unwrap(),64,1).unwrap();
        let seed=carry.take().unwrap_or_else(|| {
            let mut temperature_seed=vec![0.;8];temperature_seed[0]=1.;
            step.temperature_pullback(cx,&temperature_seed).unwrap()
        });
        let gradient=step.pullback_with_resistance_controls(cx,&seed,linear,1).unwrap();
        totals[0]+=gradient.resistance_multipliers[0].derivative;
        totals[1]+=gradient.radiation.emissivities[0];
        totals[2]+=gradient.radiation.reservoir_temperatures[0];
        carry=Some(gradient.radiation.storage.transport.previous_specific_enthalpy);
    }
    totals
}


#[test]
fn g3_radiative_contact_adjoint_matches_fresh_physical_controls() {
    context(&CancelGate::new_clock_free(),|cx| {
        let ambient=340.;
        let actual=radiative_gradients(cx,ambient);
        let e=1e-4;
        let contact=(radiative_final(cx,ambient,1.+e,0.2)
            -radiative_final(cx,ambient,1.-e,0.2))/(2.*e);
        let emissivity=(radiative_final(cx,ambient,1.,0.2+e)
            -radiative_final(cx,ambient,1.,0.2-e))/(2.*e);
        let reservoir=(radiative_final(cx,ambient+e,1.,0.2)
            -radiative_final(cx,ambient-e,1.,0.2))/(2.*e);
        for (a,b) in actual.into_iter().zip([contact,emissivity,reservoir]) {
            close(a,b,6e-5*b.abs().max(1.));
        }
        close(contact,-5.83119458,3e-5);
        close(emissivity,28.02502439,3e-4);
        close(reservoir,0.30815780,3e-5);
    });
}
