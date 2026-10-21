//! Equal-pressure material mixing through explicit binary diffusion data.
//! All gas fixtures and D data are synthetic; this is not an empirical fit.
#[path = "../material_gas_spring/source.rs"]
mod source;
use fs_material::thermochemical_gas::energy::{GasEnergySolve, parcel::ClosedGasParcel};
use fs_material::thermochemical_gas::energy::parcel::transfer::species::{GasSpeciesExchangeSolve,
    diffusion::{BinaryDiffusionMaterial, GasDiffusionGeometry, GasDiffusionSolve,
        BINARY_DIFFUSIVITY_DIMS, BINARY_DIFFUSIVITY_PROPERTY}};
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, UncertaintyModel};
use fs_evidence::ValidityDomain;
use fs_qty::{Temperature, QuantitySpec, semantic::{QuantityKind, SemanticType, ValueForm}};
use std::{collections::BTreeSet, io::Read};
type Result<T> = std::result::Result<T, Box<dyn std::error::Error>>;

fn typed(kind: QuantityKind) -> QuantitySpec { QuantitySpec::semantic(SemanticType::new(kind, ValueForm::Static)) }
fn diffusion_card(d: f64) -> Result<MaterialCard> {
    if !d.is_finite() || d <= 0.0 { return Err("explicit positive diffusivity [m2/s] required".into()); }
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(BINARY_DIFFUSIVITY_PROPERTY, BINARY_DIFFUSIVITY_DIMS),
        value: PropertyValue::Curve { abscissa: "temperature".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(200.0,d),(1800.0,d)], dims: BINARY_DIFFUSIVITY_DIMS },
        validity: ValidityDomain::unconstrained()
            .with_quantity("temperature", typed(QuantityKind::AbsoluteTemperature), 200.0, 1800.0)
            .with_quantity("pressure", typed(QuantityKind::Pressure), 1000.0, 1e7),
        uncertainty: UncertaintyModel::Unstated, interpolation: InterpolationPolicy::LinearInside,
        observations: vec![], provenance: Provenance { source: "explicit synthetic CLI diffusivity".into(),
            license: "CC0-1.0".into(), artifact: None },
    })?;
    Ok(MaterialCard::assemble(MaterialStateId { chemistry: "synthetic binary pair".into(),
        phase: "gas".into(), process: "example".into(), revision: 0 }, claims, vec![])?)
}
fn controls() -> GasDiffusionSolve {
    let gas = GasEnergySolve { minimum_temperature_k: 200.0, maximum_temperature_k: 1800.0,
        energy_tolerance_j_kg: 1e-7, maximum_iterations: 128 };
    GasDiffusionSolve { maximum_temperature_difference_k: 1e-5, maximum_pressure_difference_pa: 0.01,
        gradient_tolerance_per_m: 1e-9, molar_flux_tolerance_mol_m2_s: 1e-12, entropy_tolerance_j_k: 1e-10,
        exchange: GasSpeciesExchangeSolve { gas: [gas;2], maximum_species_fraction: 0.1,
            mass_tolerance_kg: 1e-14, energy_tolerance_j: 1e-8 } }
}
fn read(path: &str) -> Result<String> {
    let mut text=String::new(); std::fs::File::open(path)?.take(65537).read_to_string(&mut text)?;
    if text.len()>65536 { return Err("synthetic mixture file exceeds 64 KiB".into()); }
    Ok(text)
}
fn main() -> Result<()> {
    let args=std::env::args().skip(1).collect::<Vec<_>>();
    if !(3..=5).contains(&args.len()) {
        return Err("usage: material_gas_diffusion LEFT.tsv RIGHT.tsv D_m2_s [steps=200] [dt_s=0.05]".into());
    }
    let d=diffusion_card(args[2].parse()?)?;
    let steps:usize=args.get(3).map_or(Ok(200),|x|x.parse())?;
    let dt:f64=args.get(4).map_or(Ok(0.05),|x|x.parse())?;
    if !(1..=10000).contains(&steps) || !dt.is_finite() || dt<=0.0 || !(dt*steps as f64).is_finite() {
        return Err("bounded step count and positive finite duration required".into());
    }
    let provenance=format!("paired synthetic diffusion inputs: {} + {}",args[0],args[1]);
    let point=source::point(400.0,100000.0)?;
    let mut left=ClosedGasParcel::new(source::load(&read(&args[0])?,&provenance)?,point.clone(),0.001)?;
    let mut right=ClosedGasParcel::new(source::load(&read(&args[1])?,&provenance)?,point.clone(),0.001)?;
    let ids=[&left,&right].iter().flat_map(|p|p.model().thermodynamics().components().iter()
        .map(|s|s.species().clone())).collect::<BTreeSet<_>>().into_iter().collect::<Vec<_>>();
    if ids.len()!=2 { return Err("this demo needs exactly two species; the library supports up to 32".into()); }
    let pairs=[BinaryDiffusionMaterial { species:[ids[0].clone(),ids[1].clone()],material:&d,
        context:&point,selection:MaterialPropertySelection::SingleClaimOnly }];
    let geometry=GasDiffusionGeometry {area_m2:0.01,distance_m:0.01};
    println!("time_s,left_first_mole_fraction,right_first_mole_fraction,left_pressure_pa,right_pressure_pa,left_temperature_k,right_temperature_k,left_sound_speed_m_s,right_sound_speed_m_s,step_entropy_change_j_k,step_energy_error_j");
    for i in 0..=steps {
        let (ds,de)=if i==0 {(0.0,0.0)} else {
            let proposal=left.diffuse_with(&right,geometry,dt,&pairs,controls(),||false)?;
            let ds=proposal.entropy_change_j_k(); let de=proposal.exchange().step().energy_residual_j;
            let (exchange,_,_)=proposal.into_parts(); let (a,b,_)=exchange.into_parts();left=a;right=b;(ds,de)
        };
        let first=|p:&ClosedGasParcel| p.model().thermodynamics().components().iter()
            .zip(p.model().thermodynamics().mole_composition().fractions())
            .find(|(s,_)|s.species()==&ids[0]).map_or(0.0,|(_,x)|*x);
        println!("{:.9},{:.12},{:.12},{:.9},{:.9},{:.9},{:.9},{:.9},{:.9},{:.12e},{:.12e}",
            i as f64*dt,first(&left),first(&right),left.state().state().pressure,right.state().state().pressure,
            left.state().state().temperature,right.state().state().temperature,
            left.state().state().sound_speed,right.state().state().sound_speed,ds,de);
    }
    Ok(())
}
