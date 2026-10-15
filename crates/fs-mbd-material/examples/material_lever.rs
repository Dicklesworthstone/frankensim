//! A data-only material substitution evaluated by the actual articulated solver.
//! Run: cargo run -p fs-mbd-material --example material_lever -- 1000 3000
//! Inputs are synthetic density data [kg/m3], not measured named-material claims.
use std::error::Error;
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, So3, Vec3, Wrench};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::state_point::MaterialPropertySelection;
use fs_mbd::articulated::{BaseState, JointModel, forward_dynamics, inverse_dynamics};
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_mbd_material::articulated::{MaterialLink, resolve_articulated};
use fs_qty::{Density, Temperature};

fn material(density: f64) -> Result<MaterialCard, Box<dyn Error>> {
    if !density.is_finite() || density <= 0.0 { return Err("density must be finite and positive".into()); }
    let mut claims=ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key:PropertyKey::new("density",Density::DIMS),
        value:PropertyValue::Curve { abscissa:"T".into(),abscissa_dims:Temperature::DIMS,
            knots:vec![(290.0,density),(310.0,density)],dims:Density::DIMS },
        validity:ValidityDomain::unconstrained().with("T",290.0,310.0),
        uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::LinearInside,
        observations:vec![],provenance:Provenance {source:"explicit synthetic lever density".into(),
            license:"CC0-1.0".into(),artifact:None},
    })?;
    Ok(MaterialCard::assemble(MaterialStateId {chemistry:"synthetic-material".into(),phase:"solid".into(),
        process:"illustrative uniform rigid solid".into(),revision:0},claims,vec![])? )
}

fn response(density: f64) -> Result<[f64;4],Box<dyn Error>> {
    let card=material(density)?;
    let point=QueryPoint::new().with("T",300.0)?;
    let zero=Vec3::new(0.0,0.0,0.0);
    // A two-metre lever: pivot at x=0, centroid at x=1, full width/depth 0.2m.
    let geometry=VolumeMoments::cuboid(Vec3::new(2.0,0.2,0.2))?
        .transformed(Se3::from_parts(So3::identity(),Vec3::new(1.0,0.0,0.0))?)?;
    let regions=[RegionBinding {key:"lever-solid",geometry,material:&card,point:&point,
        selection:MaterialPropertySelection::SingleClaimOnly}];
    let mechanism=resolve_articulated(&[MaterialLink {name:"lever",parent:None,
        parent_from_child_zero:Se3::identity(),joint:JointModel::revolute(Vec3::new(0.0,0.0,1.0),None)?,
        regions:&regions}],1,||false)?;
    let model=mechanism.model();let base=BaseState::stationary(Se3::identity());
    let external=[Wrench::new(zero,zero)];
    let driven=forward_dynamics(model,base,&[0.0],&[0.0],&[2.0],zero,&external)?;
    let holding=inverse_dynamics(model,base,&[0.3],&[0.0],&[0.0],Vec3::new(0.0,-9.81,0.0),&external)?;
    let inertia=mechanism.bodies()[0].inertia();
    Ok([inertia.mass(),inertia.matrix()?.m[14],driven.generalized_acceleration[0],holding.generalized_force[0]])
}

fn main() -> Result<(),Box<dyn Error>> {
    let args:Vec<_>=std::env::args().skip(1).collect();
    let densities=match args.as_slice() {
        []=>[1000.0,3000.0],
        [a,b]=>[a.parse::<f64>()?,b.parse::<f64>()?],
        _=>return Err("usage: material_lever [density_a_kg_m3 density_b_kg_m3]".into()),
    };
    println!("Synthetic rigid specimens; unchanged geometry, joint and applied torque.");
    println!("density_kg_m3,mass_kg,pivot_inertia_kg_m2,driven_acceleration_rad_s2,gravity_holding_torque_nm");
    for density in densities {
        let [mass,inertia,acceleration,holding]=response(density)?;
        println!("{density:.9},{mass:.9},{inertia:.9},{acceleration:.9},{holding:.9}");
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn executable_uses_source_density_in_actual_mechanics() {
        let a=response(1000.0).unwrap(); let b=response(3000.0).unwrap();
        for i in [0,1,3] { assert!((b[i]/a[i]-3.0).abs()<1e-11); }
        assert!((a[2]/b[2]-3.0).abs()<1e-11);
        assert!(response(0.0).is_err()); assert!(response(f64::NAN).is_err());
    }
}
