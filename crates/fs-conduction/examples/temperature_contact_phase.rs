//! Contact-driven melting with source-backed R(T) and a radiative reservoir.
//! All material inputs are synthetic numerical examples, not measurements.
//! cargo run -p fs-conduction --example temperature_contact_phase -- 300 1
//! Arguments: reservoir kelvin and positive multiplier on the source R(T) curve.
use std::error::Error;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::{AmbientRadiationConfig, AmbientRadiationPatch, ConductionMesh,
    ConductionProblem, ConductivityModel, ConductivityTable, ScalarField, SurfaceEmissivity,
    ThermalBc, ThermalBoundaryBuilder, ThermalInterfaces, EMISSIVITY_DIMS,
    SURFACE_EMISSIVITY_PROPERTY, AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS,
    AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY, TEMPERATURE_DIMS};
use fs_conduction::interface::temperature::{ContactIterationConfig, TemperatureContact, TemperatureContacts};
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::{HeterogeneousEnthalpyBackwardEuler, ReferenceEnthalpyMaterial};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, SurfaceSpec,
    SystemContext, UncertaintyModel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_rep_mesh::TetComplex;
use fs_solver::NewtonKrylovConfig;

fn state(name:&str)->MaterialStateId {
    MaterialStateId{chemistry:name.into(),phase:"solid".into(),process:"synthetic reference".into(),revision:0}
}
fn provenance()->Provenance {
    Provenance{source:"synthetic fixed-geometry phase/contact numerical example".into(),
        license:"example-only; not experimental data".into(),artifact:None}
}
fn phase(hot:bool)->Result<EquilibriumEnthalpyPhaseCurve,Box<dyn Error>> {
    let knots=if hot {[(-200.,280.,0.),(500.,350.,0.),(600.,350.,1.),(900.,380.,1.)]}
        else {[(-200.,280.,0.),(100.,310.,0.),(200.,310.,1.),(900.,380.,1.)]};
    Ok(EquilibriumEnthalpyPhaseCurve::try_new(ContentHash([if hot {21}else{20};32]),
        knots.into_iter().map(|(h,t,f)|EnthalpyPhaseKnot{specific_enthalpy_j_kg:h,
            temperature_k:t,liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?)
}
fn run(cx:&Cx<'_>,ambient:f64,scale:f64)->Result<(),Box<dyn Error>> {
    let mesh=ConductionMesh::new(TetComplex::from_tets(8,vec![[0,1,2,3],[4,5,6,7]]),
        vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,-1.],
            [0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])?;
    let boundary=ThermalBoundaryBuilder::new(&mesh).region("exterior",
        |f|f.centroid[2]<0. && f.centroid[1]==0.,ThermalBc::neumann(0.)?)?
        .adiabatic_remainder().finish()?;
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim{
        key:PropertyKey::new(AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY,AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS),
        value:PropertyValue::Curve{abscissa:"T".into(),abscissa_dims:TEMPERATURE_DIMS,
            knots:vec![(280.,0.1*scale),(380.,0.7*scale)],dims:AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.).with("normal_pressure",1e6,1e6),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:provenance(),
    })?;
    let side=|name|SurfaceSpec{material:state(name),texture_frame:"fixed reference normal".into()};
    let card=InterfaceSystemCard::assemble(side("body A"),side("body B"),SystemContext{
        medium:"dry".into(),third_body:None,environment:"vacuum".into(),history:"unchanged".into(),
    },claims,vec![])?;
    let point=QueryPoint::new().with("T",300.)?.with("normal_pressure",1e6)?;
    let law=TemperatureContact::new("joint",ThermalInterfaces::coincident_face_pairs(&mesh)?,
        &card,&point,"T",SelectionPolicy::SingleClaimOnly)?;
    let contacts=TemperatureContacts::new(cx,&mesh,&boundary,vec![],vec![law],1)?;
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim{
        key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar{value:0.2,dims:EMISSIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
        observations:vec![],provenance:provenance(),
    })?;
    let surface=MaterialCard::assemble(state("gray surface"),claims,vec![])?;
    let emissivity=SurfaceEmissivity::from_card("exterior",&surface,300.,SelectionPolicy::SingleClaimOnly)?;
    let patches=[AmbientRadiationPatch::new("exterior",emissivity,ambient)?];
    let cold=phase(false)?;let hot=phase(true)?;
    let engine=HeterogeneousEnthalpyBackwardEuler::new(cx,&mesh,&[
        ReferenceEnthalpyMaterial{curve:&cold,reference_density_kg_m3:1.},
        ReferenceEnthalpyMaterial{curve:&hot,reference_density_kg_m3:2.},
    ],&[0,1],EnthalpyBudget{max_vertices:8,max_elements:2})?;
    let material=ConductivityModel::isotropic(ConductivityTable::declared_curve(vec![(280.,1.6),(380.,3.6)])?);
    let source=ScalarField::Nodal(vec![0.,0.,0.,0.,60.,60.,60.,60.]);
    let problem=ConductionProblem{mesh:&mesh,boundary:&boundary,material:&material,
        element_materials:None,source:&source};
    let config=EnthalpyStepConfig{newton:NewtonKrylovConfig{
        absolute_tolerance:1e-11,relative_tolerance:1e-12,linear_restart:8,max_linear_cycles:16,
        forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()
    },max_newton_iterations:64,energy_tolerance_j:1e-8};
    let coupling=ContactIterationConfig{max_iterations:64,relaxation:0.5,residual_tolerance_j:1e-8};
    let mut h=vec![80.,80.,80.,80.,450.,450.,450.,450.];let mut time=0.;
    println!("time_s,node,reservoir_k,resistance_scale,temperature_k,h_j_kg,liquid_fraction,contact_r_m2_k_w,stored_change_j,radiation_out_w,residual_j,contact_trials,radiation_trials");
    for dt in [0.15,0.35,0.5] {
        let next=contacts.advance_enthalpy_with_ambient_radiation(cx,&engine,problem,&h,dt,
            config,coupling,&patches,AmbientRadiationConfig::default())?;
        time+=dt;
        for node in 0..8 {
            println!("{time},{node},{ambient},{scale},{},{},{},{},{},{},{},{},{}",
                next.conduction.temperature[node],next.conduction.specific_enthalpy_j_kg[node],
                next.conduction.liquid_mass_fraction[node],
                next.physical_contacts.states()[0].resistance.value_m2_k_per_w(),
                next.physical.stored_energy_change_j,next.physical.radiation_out_w,
                next.physical.residual_norm_j,next.iterations,next.radiation_iterations);
        }
        h=next.conduction.specific_enthalpy_j_kg;
    }
    Ok(())
}
fn main()->Result<(),Box<dyn Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()!=2 {return Err("usage: temperature_contact_phase RESERVOIR_K RESISTANCE_SCALE".into());}
    let ambient:f64=args[0].parse()?;let scale:f64=args[1].parse()?;
    if !ambient.is_finite()||ambient<=0.||!scale.is_finite()||scale<=0. {
        return Err("reservoir kelvin and resistance scale must be finite and positive".into());
    }
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena|run(&Cx::new(&gate,arena,
        StreamKey{seed:20,kernel_id:820,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic),ambient,scale))
}
