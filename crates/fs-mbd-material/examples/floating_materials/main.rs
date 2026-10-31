//! Authored synthetic solid/liquid data -> shared resolvers -> floating motion.
//! These cards are demonstration inputs, not measurements of named substances.
use std::error::Error;
use fs_evidence::ValidityDomain;
use fs_ga::{Se3,So3,Twist,Vec3};
use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialCard,MaterialStateId,PropertyClaim,
    PropertyKey,PropertyValue,Provenance,QueryPoint,UncertaintyModel};
use fs_material::liquid::{ResolvedLiquidState,resolve_liquid_state};
use fs_material::state_point::{DENSITY_PROPERTY,MaterialPropertySelection,
    SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS,THERMAL_CONDUCTIVITY_PROPERTY};
use fs_mbd::Gravity;
use fs_mbd_material::{RegionBinding,ResolvedBody,VolumeMoments,resolve_body};
use fs_mbd_material::hydrostatic::{DisplacementMesh,HydrostaticLayer,LayeredHydrostaticMedium,
    HydrostaticStepControls};
use fs_qty::{Density,DynViscosity,QuantitySpec,Temperature};
use fs_qty::semantic::{QuantityKind,SemanticType,ValueForm};

type Result<T> = std::result::Result<T,Box<dyn Error>>;
fn zero() -> Vec3 { Vec3::new(0.0,0.0,0.0) }
fn sides() -> Vec3 { Vec3::new(2.0,2.0,1.0) }
fn gravity() -> Gravity { Gravity::new(fs_mbd::Vec3::new(0.0,0.0,-10.0)).expect("finite gravity") }

fn card(phase: &str, properties: Vec<(PropertyKey,f64)>) -> Result<MaterialCard> {
    let mut claims=ClaimSet::new();
    for (key,value) in properties {
        let dims=key.dims();
        claims.insert_claim(PropertyClaim {
            key,value:PropertyValue::Curve {abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
                knots:vec![(200.0,value),(400.0,value)],dims},
            validity:ValidityDomain::unconstrained().with("T",200.0,400.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
            observations:vec![],provenance:Provenance {source:"authored hydrostatic demonstration".into(),
                license:"CC0-1.0".into(),artifact:None},
        })?;
    }
    Ok(MaterialCard::assemble(MaterialStateId {chemistry:"synthetic-hydrostatic".into(),
        phase:phase.into(),process:"authored data; not measured".into(),revision:0},claims,vec![])?)
}

fn materials(solid_density: f64, liquid_density: f64) -> Result<(ResolvedBody,ResolvedLiquidState)> {
    let solid=card("solid",vec![(PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),solid_density)])?;
    let point=QueryPoint::new().with("T",300.0)?;
    let body=resolve_body(&[RegionBinding {key:"solid",geometry:VolumeMoments::cuboid(sides())?,
        material:&solid,point:&point,selection:MaterialPropertySelection::SingleClaimOnly}],1,||false)?;
    Ok((body,liquid_material(liquid_density)?))
}

fn liquid_material(density: f64) -> Result<ResolvedLiquidState> {
    // These unused heat/viscosity data belong to the demonstration card, not to
    // the hydrostatics solver. Its common snapshot also serves heat/flow APIs.
    let conductivity=QuantitySpec::semantic(SemanticType::new(QuantityKind::ThermalConductivity,ValueForm::Static));
    let liquid=card("liquid",vec![
        (PropertyKey::new(DENSITY_PROPERTY,Density::DIMS),density),
        (PropertyKey::new(SPECIFIC_HEAT_CAPACITY_PROPERTY,SPECIFIC_HEAT_CAPACITY_DIMS),4000.0),
        (PropertyKey::with_quantity(THERMAL_CONDUCTIVITY_PROPERTY,conductivity),0.5),
        (PropertyKey::new("dynamic_viscosity",DynViscosity::DIMS),0.001),
    ])?;
    let point=QueryPoint::new().with("T",300.0)?;
    let liquid=resolve_liquid_state(&liquid,&point,MaterialPropertySelection::SingleClaimOnly)?;
    Ok(liquid)
}

fn main() -> Result<()> {
    let args:Vec<String>=std::env::args().skip(1).collect();
    if args.len()!=2 && args.len()!=4 {
        return Err("usage: cargo run -p fs-mbd-material --example floating_materials -- SOLID_DENSITY_KG_M3 SURFACE_DENSITY_KG_M3 [LOWER_DENSITY_KG_M3 INTERFACE_DEPTH_M]".into());
    }
    let solid_density=args[0].parse::<f64>()?;
    let liquid_density=args[1].parse::<f64>()?;
    let (source,liquid)=materials(solid_density,liquid_density)?;
    let body=source.principal_body(1e-12,||false)?;
    let mesh=DisplacementMesh::cuboid(sides(),||false)?;
    let lower=if args.len()==4 {
        Some((args[3].parse::<f64>()?,liquid_material(args[2].parse::<f64>()?)?))
    } else { None };
    let mut layers=vec![HydrostaticLayer {depth_m:0.0,source:&liquid}];
    if let Some((depth,source))=&lower {
        layers.push(HydrostaticLayer {depth_m:*depth,source});
    }
    let medium=LayeredHydrostaticMedium::new(None,&layers,gravity(),zero(),||false)?;
    let mut state=body.state_from_kinematics(Se3::from_parts(So3::identity(),Vec3::new(0.0,0.0,0.02))?,Twist::zero())?;
    // Quasi-static pressure reservoir only: no waves, added mass, drag or
    // physical decay. This is not a water-entry/slamming model.
    eprintln!("synthetic materials; body_mass_kg={}; displaced_envelope_m3={}",body.properties().mass(),mesh.volume_m3());
    if let Some((depth,_))=&lower {
        eprintln!("lower_layer_density_kg_m3={}; interface_depth_m={depth}",
            medium.layers().nth(1).expect("lower layer was supplied").1.density_kg_m3());
    } else if solid_density<liquid_density {
        eprintln!("upright_equilibrium_draft_m={}",solid_density/liquid_density);
    } else {
        eprintln!("no strictly partially immersed upright equilibrium at this density ratio");
    }
    println!("time_s,height_m,vertical_velocity_m_s,submerged_volume_m3,displaced_fluid_mass_kg,buoyancy_n,numerical_energy_change_j");
    let dt=0.001; let mut cumulative=0.0;
    for n in 1..=2000 {
        let step=medium.step(&mesh,&body,state,HydrostaticStepControls {
            duration_s:dt,maximum_energy_defect_j:1.0,
        },||false)?;
        // Commit only after both load evaluations and the energy gate succeed.
        state=step.dynamics.state_after; cumulative+=step.energy_defect_j;
        if n%100==0 {
            println!("{:.3},{:.9},{:.9},{:.9},{:.9},{:.6},{:.9}",f64::from(n)*dt,
                body.world_from_body(state)?.translation().z,
                state.center_of_mass_velocity_world(body.properties())?.z,
                step.load_after.displaced().volume_m3(),step.load_after.displaced_mass_kg(),
                step.load_after.wrench().force_world.z,cumulative);
        }
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_material::fluid_transport::FluidTransport;

    #[test]
    fn both_material_resolvers_feed_the_same_executable_load_and_step() -> Result<()> {
        for (solid_rho,liquid_rho,expected) in [(500.0,1000.0,0.0),(750.0,1000.0,-10.0/3.0),(500.0,1200.0,2.0)] {
            let (source,liquid)=materials(solid_rho,liquid_rho)?;
            let body=source.principal_body(1e-12,||false)?;
            let mesh=DisplacementMesh::cuboid(sides(),||false)?;
            let medium=LayeredHydrostaticMedium::new(None,
                &[HydrostaticLayer {depth_m:0.0,source:&liquid}],gravity(),zero(),||false)?;
            let (layer,transport)=medium.layers().next().expect("surface layer was supplied");
            assert!(std::ptr::addr_eq(layer.source,&liquid as &dyn FluidTransport));
            assert_eq!(transport,liquid.transport_state()?);
            let state=body.state_from_kinematics(Se3::identity(),Twist::zero())?;
            let step=medium.step(&mesh,&body,state,HydrostaticStepControls {
                duration_s:0.0001,maximum_energy_defect_j:1.0,
            },||false)?;
            let acceleration=step.dynamics.state_after.center_of_mass_velocity_world(body.properties())?.z/0.0001;
            assert!((acceleration-expected).abs()<1e-9);
        }
        Ok(())
    }

    #[test]
    fn invalid_material_data_is_not_replaced_by_water_or_a_default_solid() {
        for bad in [0.0,-1.0,f64::NAN,f64::INFINITY] {
            assert!(materials(bad,1000.0).is_err());
            assert!(materials(500.0,bad).is_err());
        }
    }

    #[test]
    fn resolved_layer_materials_change_flotation_through_the_same_native_step() -> Result<()> {
        let (source,surface)=materials(375.0,500.0)?;
        let body=source.principal_body(1e-12,||false)?;
        let mesh=DisplacementMesh::cuboid(sides(),||false)?;
        let state=body.state_from_kinematics(Se3::identity(),Twist::zero())?;
        for (lower_density,expected_mass) in [(500.0,1000.0),(1000.0,1500.0),(2000.0,2500.0)] {
            let lower=liquid_material(lower_density)?;
            let medium=LayeredHydrostaticMedium::new(None,&[
                HydrostaticLayer {depth_m:0.0,source:&surface},
                HydrostaticLayer {depth_m:0.25,source:&lower},
            ],gravity(),zero(),||false)?;
            let dt=0.0001;
            let step=medium.step(&mesh,&body,state,HydrostaticStepControls {
                duration_s:dt,maximum_energy_defect_j:1.0,
            },||false)?;
            let acceleration=step.dynamics.state_after.center_of_mass_velocity_world(body.properties())?.z/dt;
            assert!((step.load_before.displaced_mass_kg()-expected_mass).abs()<1e-9);
            assert!((acceleration-10.0*(expected_mass/1500.0-1.0)).abs()<1e-9);
        }
        Ok(())
    }
}
