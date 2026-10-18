//! Scheduled radiative melting/refreezing; every material number is synthetic.
//! cargo run -p fs-conduction --example scheduled_melting -- 100 350
//! Arguments are latent heat [J/kg] and hot-reservoir temperature [K].
use std::error::Error;
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{HeatSourceHistory, ThermalHistoryProblem};
use fs_conduction::bc::history::march::enthalpy::{EnthalpyHistoryBudget, EnthalpyHistoryMaterials, EnthalpyHistoryRun};
use fs_conduction::bc::history::march::enthalpy::radiation::RadiativeEnthalpyHistoryRun;
use fs_conduction::bc::history::march::radiation::RadiationRegionHistory;
use fs_conduction::transient::enthalpy::{EnthalpyBudget, EnthalpyStepConfig};
use fs_conduction::transient::enthalpy::heterogeneous::ReferenceEnthalpyMaterial;
use fs_conduction::{AmbientRadiationConfig, ConductionMesh, ConductivityModel, EMISSIVITY_DIMS,
    SURFACE_EMISSIVITY_PROPERTY, ScalarField, SurfaceEmissivity, TEMPERATURE_DIMS, ThermalBc, ThermalBoundaryBuilder};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, SelectionPolicy, UncertaintyModel};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use fs_solver::NewtonKrylovConfig;

fn run(cx: &Cx<'_>, latent: f64, hot: f64) -> Result<(),Box<dyn Error>> {
    let mesh=ConductionMesh::new(TetComplex::from_tets(4,vec![[0,1,2,3]]),
        vec![[0.,0.,0.],[1.,0.,0.],[0.,1.,0.],[0.,0.,1.]])?;
    // Synthetic identity tag includes the supplied latent-heat bits.
    let mut identity=[54;32]; identity[..8].copy_from_slice(&latent.to_le_bytes());
    let curve=EquilibriumEnthalpyPhaseCurve::try_new(ContentHash(identity),
        [(-200.,280.,0.),(100.,310.,0.),(100.+latent,310.,1.),(800.+latent,380.,1.)]
            .into_iter().map(|(h,t,f)| EnthalpyPhaseKnot {specific_enthalpy_j_kg:h,temperature_k:t,
                liquid_mass_fraction:f,bulk_density_kg_m3:1.}).collect())?;
    let boundary=ThermalBoundaryHistory::new(cx,ThermalBoundaryBuilder::new(&mesh)
        .remainder("surface",ThermalBc::adiabatic())?,&Environment::earth_lab(),Vec::new(),0.)?;
    let conductivity=ConductivityModel::isotropic_declared(2.)?;
    let materials=[ReferenceEnthalpyMaterial {curve:&curve,reference_density_kg_m3:1.}];
    let base=EnthalpyHistoryRun::new(cx,ThermalHistoryProblem {boundary:&boundary,material:&conductivity,
        element_materials:None,interfaces:None,source:HeatSourceHistory::Fixed(ScalarField::Uniform(0.))},
        EnthalpyHistoryMaterials {table:&materials,of_element:&[0],phase_conductivity:None},&[80.;4],
        EnthalpyHistoryBudget {spatial:EnthalpyBudget {max_vertices:4,max_elements:1},max_steps:100})?;
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {key:PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY,EMISSIVITY_DIMS),
        value:PropertyValue::Scalar {value:0.8,dims:EMISSIVITY_DIMS},
        validity:ValidityDomain::unconstrained().with("T",280.,380.),uncertainty:UncertaintyModel::Unstated,
        interpolation:InterpolationPolicy::ConstantWithinValidity,observations:Vec::new(),
        provenance:Provenance {source:"synthetic radiative melting example".into(),
            license:"example-only; not experimental data".into(),artifact:None}})?;
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic gray body".into(),phase:"solid".into(),
        process:"illustrative constant emissivity".into(),revision:0},claims,Vec::new())?;
    let emissivity=SurfaceEmissivity::from_card("surface",&card,320.,SelectionPolicy::SingleClaimOnly)?;
    let patch=RadiationRegionHistory {region:"surface".into(),emissivity,
        ambient_temperature:TimeSignal::Table {times:vec![0.,0.051],values:vec![hot,300.],
            dims:TEMPERATURE_DIMS,interp:Interp::Hold}};
    let mut run=RadiativeEnthalpyHistoryRun::new(cx,base,vec![patch],AmbientRadiationConfig::default())?;
    let config=EnthalpyStepConfig {newton:NewtonKrylovConfig {absolute_tolerance:1e-11,relative_tolerance:1e-12,
        linear_restart:4,max_linear_cycles:16,forcing_minimum:1e-12,forcing_maximum:1e-3,..Default::default()},
        max_newton_iterations:64,energy_tolerance_j:1e-8};
    println!("time_s,node,temperature_K,enthalpy_J_kg,liquid_mass_fraction,reservoir_K,storage_change_J,radiation_out_W,energy_defect_J");
    for step in 1..=100 {
        let time=f64::from(step)*0.001;
        let solved=run.advance_to(cx,time,config)?;
        for node in 0..4 {
            println!("{time},{node},{},{},{},{},{},{},{}",solved.conduction.temperature[node],
                solved.conduction.specific_enthalpy_j_kg[node],solved.conduction.liquid_mass_fraction[node],
                solved.radiation.patches[0].patch.ambient_temperature_k(),solved.conduction.stored_energy_change_j,
                solved.radiation.nonlinear_radiation_out_w,solved.physical_energy_residual_j);
        }
    }
    Ok(())
}
fn main() -> Result<(),Box<dyn Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    if args.len()!=2 {return Err("usage: scheduled_melting LATENT_J_KG HOT_RESERVOIR_K".into());}
    let latent:f64=args[0].parse()?;let hot:f64=args[1].parse()?;
    if !latent.is_finite() || latent<=0. || !hot.is_finite() || !(310.0..=370.0).contains(&hot) {
        return Err("latent heat must be finite/positive and the hot reservoir must be in [310,370] K".into());
    }
    let gate=CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| run(&Cx::new(&gate,arena,
        StreamKey {seed:54,kernel_id:854,tile:0,iteration:0},Budget::INFINITE,ExecMode::Deterministic),latent,hot))
}
