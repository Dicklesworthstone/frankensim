//! Radiation changes latent storage without inventing a temperature rise.
use super::*;
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::bc::history::march::enthalpy::radiation::RadiativeEnthalpyHistoryRun;
use fs_conduction::{AmbientRadiationConfig, SurfaceEmissivity,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, STEFAN_BOLTZMANN_W_M2_K4};
use fs_evidence::ValidityDomain;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};

fn patch(signal: TimeSignal) -> RadiationRegionHistory {
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar {value:0.8,dims:EMISSIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:Vec::new(),provenance:Provenance {source:"synthetic phase-history fixture".into(),
            license:"test-only; not measured".into(),artifact:None},
    }).unwrap();
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"gray body".into(),phase:"solid".into(),
        process:"test".into(),revision:0},claims,Vec::new()).unwrap();
    RadiationRegionHistory {region:"surface".into(),
        emissivity:SurfaceEmissivity::from_card("surface",&card,320.,SelectionPolicy::SingleClaimOnly).unwrap(),
        ambient_temperature:signal}
}
fn surface<'m>(cx: &fs_exec::Cx<'_>, mesh: &'m ConductionMesh, htc: bool) -> ThermalBoundaryHistory<'m> {
    ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(mesh).remainder("surface",
        if htc {ThermalBc::robin(5.,300.).unwrap()} else {ThermalBc::adiabatic()}).unwrap(),
        &Environment::earth_lab(),Vec::new(),0.).unwrap()
}

#[test]
fn g1_radiative_heating_and_cooling_change_latent_fraction_at_fixed_temperature() {
    let mesh=mesh(false);let curve=curve(100.,51);
    let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let table=[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}];
    with_cx(|cx| {
        let history=surface(cx,&mesh,false);
        let base=EnthalpyHistoryRun::new(cx,problem(&history,&material,
            HeatSourceHistory::Fixed(ScalarField::Uniform(0.))),
            EnthalpyHistoryMaterials {table:&table,of_element:&[0],phase_conductivity:None},
            &[150.;4],budget(&mesh,2)).unwrap();
        let original=patch(TimeSignal::Table {times:vec![0.,0.02],values:vec![320.,300.],
            dims:TEMPERATURE_DIMS,interp:Interp::Hold});
        let mut run=RadiativeEnthalpyHistoryRun::new(cx,base,vec![original.clone()],AmbientRadiationConfig::default()).unwrap();
        let mut expected=vec![150.;4];
        for (time,ambient) in [(0.01,320.0_f64),(0.02,300.0_f64)] {
            let flux=0.8*STEFAN_BOLTZMANN_W_M2_K4*(310.0_f64.powi(4)-ambient.powi(4));
            for face in mesh.boundary() {
                for &vertex in &face.vertices {
                    expected[vertex as usize]-=0.01*flux*face.area/3./(1./24.);
                }
            }
            let result=run.advance_to(cx,time,config()).unwrap();
            assert!(result.conduction.temperature.iter().all(|&t| t==310.));
            for (v,&h) in expected.iter().enumerate() {
                assert!((result.conduction.specific_enthalpy_j_kg[v]-h).abs()<1e-8);
                assert!((result.conduction.liquid_mass_fraction[v]-(h-100.)/100.).abs()<1e-9);
            }
            assert_eq!(result.convective_out_w,0.);
            assert!(result.convective_robin_fluxes.is_empty());
            assert!(result.physical_energy_residual_j.abs()<=config().energy_tolerance_j);
            assert!(result.physical_residual_norm_j<=result.physical_residual_tolerance_j);
            assert_eq!(result.radiation.patches[0].patch.emissivity(),&original.emissivity);
            assert_eq!(result.radiation.nonlinear_radiation_out_w<0.,ambient>310.);
        }
    });
}

#[test]
fn g3_radiative_phase_schedule_keeps_original_nonlinear_solver_and_heat_split() {
    let mesh=mesh(false);let curve=curve(100.,52);
    let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.6),(380.,3.6)]).unwrap());
    let table=[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}];
    let phase=LiquidMassFractionConductivity::declared(1.,3.).unwrap();
    with_cx(|cx| {
        let history=surface(cx,&mesh,true);
        let original=patch(TimeSignal::Table {times:vec![0.,0.03],values:vec![350.,290.],
            dims:TEMPERATURE_DIMS,interp:Interp::Hold});
        let source=ScalarField::Nodal(vec![40.,60.,100.,80.]);
        let mut old=vec![80.,135.,155.,240.];
        let base=EnthalpyHistoryRun::new(cx,problem(&history,&material,HeatSourceHistory::Fixed(source.clone())),
            EnthalpyHistoryMaterials {table:&table,of_element:&[0],phase_conductivity:Some(vec![phase])},
            &old,budget(&mesh,5)).unwrap();
        let mut run=RadiativeEnthalpyHistoryRun::new(cx,base,vec![original.clone()],AmbientRadiationConfig::default()).unwrap();
        let direct=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&table,&[0],budget(&mesh,5).spatial).unwrap()
            .with_phase_conductivity(cx,vec![phase]).unwrap();
        let mut previous=0.;
        for time in [0.005,0.01,0.02,0.03,0.04] {
            let boundary=history.at(cx,time).unwrap();
            let expected=direct.advance_with_ambient_radiation(cx,ConductionProblem {mesh:&mesh,boundary:&boundary,
                material:&material,element_materials:None,source:&source},None,&old,time-previous,config(),
                &[original.at(time).unwrap()],AmbientRadiationConfig::default()).unwrap();
            let actual=run.advance_to(cx,time,config()).unwrap();
            assert_eq!(actual.conduction.specific_enthalpy_j_kg,expected.conduction.specific_enthalpy_j_kg);
            assert_eq!(actual.conduction.liquid_mass_fraction,expected.conduction.liquid_mass_fraction);
            assert_eq!(actual.radiation,expected.radiation);
            assert_eq!(actual.physical_energy_residual_j,expected.physical_energy_residual_j);
            assert_eq!(actual.convective_robin_fluxes[0].mean_reference_temperature_k,300.);
            assert_eq!(actual.radiation.nonlinear_radiation_out_w<0.,time<0.03);
            old=expected.conduction.specific_enthalpy_j_kg;previous=time;
        }
        assert_eq!(run.conduction().specific_enthalpy_j_kg(),old.as_slice());
        assert_eq!(run.patches(),&[original]);
    });
}

#[test]
fn g0_radiative_refusals_leave_the_entire_phase_state_resumable() {
    let mesh=mesh(false);let curve=curve(100.,53);
    let material=ConductivityModel::isotropic_declared(2.).unwrap();
    let table=[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}];
    let history=with_cx(|cx| surface(cx,&mesh,false));
    let create=|cx: &fs_exec::Cx<'_>,max_steps| EnthalpyHistoryRun::new(cx,
        problem(&history,&material,HeatSourceHistory::Fixed(ScalarField::Uniform(0.))),
        EnthalpyHistoryMaterials {table:&table,of_element:&[0],phase_conductivity:None},
        &[150.;4],budget(&mesh,max_steps)).unwrap();
    let original=patch(TimeSignal::Table {times:vec![0.,1.],values:vec![320.,-1.],
        dims:TEMPERATURE_DIMS,interp:Interp::Hold});
    let mut run=with_cx(|cx| RadiativeEnthalpyHistoryRun::new(cx,create(cx,1),
        vec![original.clone()],AmbientRadiationConfig::default()).unwrap());
    with_cancelled_cx(|cx| assert!(run.advance_to(cx,0.01,config()).is_err()));
    assert_eq!(run.conduction().accepted_steps(),0);
    with_cx(|cx| {
        run.advance_to(cx,0.01,config()).unwrap();
        let h=run.conduction().specific_enthalpy_j_kg().to_vec();
        let t=run.conduction().temperature().to_vec();let f=run.conduction().liquid_mass_fraction().to_vec();
        assert!(run.advance_to(cx,0.02,config()).is_err());run.extend_step_budget(2).unwrap();
        assert!(run.advance_to(cx,1.,config()).is_err());
        assert_eq!(run.conduction().specific_enthalpy_j_kg(),h.as_slice());
        assert_eq!(run.conduction().temperature(),t.as_slice());assert_eq!(run.conduction().liquid_mass_fraction(),f.as_slice());
        assert_eq!(run.conduction().time_s(),0.01);assert_eq!(run.conduction().accepted_steps(),1);
        run.advance_to(cx,0.02,config()).unwrap();
        let mut twin=RadiativeEnthalpyHistoryRun::new(cx,create(cx,2),vec![original.clone()],AmbientRadiationConfig::default()).unwrap();
        twin.advance_to(cx,0.01,config()).unwrap();twin.advance_to(cx,0.02,config()).unwrap();
        assert_eq!(run.conduction().specific_enthalpy_j_kg(),twin.conduction().specific_enthalpy_j_kg());
        assert!(RadiativeEnthalpyHistoryRun::new(cx,create(cx,1),vec![original.clone(),original],AmbientRadiationConfig::default()).is_err());
    });
}
