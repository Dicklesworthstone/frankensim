//! Imported source -> material law -> nonlinear spatial equilibrium.
//! This scene provides supports and loads, never textile stiffness constants.
//!
//! Usage (after compiling an existing source-qualified material pack):
//! `textile_sheet PACK --membrane 293.15 1.02 0.05 9.81 0.3 0.1`
//! The numbers are temperature [K], imposed perimeter stretch ratio, downward
//! centre-node dead force [N], downward gravity [m/s2], grain angle [rad], and
//! maximum absolute engineering Green-strain component. The last is a caller
//! operating bound, not an experimentally established source validity range.
//! The reference patch is 0.1 m square with eight P1 triangles. Output reports
//! displacements and full force residuals; constrained entries are support
//! reactions. Positive prestretch is a physical boundary condition, not a
//! fabricated numerical stiffness floor. No bending, contact, wrinkle
//! relaxation, thermal evolution or mesh-converged garment response is claimed.
use super::{Result, query};
use fs_material::sheet::{SheetCapabilities, resolve_sheet};
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::MaterialCard;
use fs_solid::membrane::{MembraneBudget, MembraneDomain, MembranePatch, MembraneSolution, MembraneSolveSettings};
use fs_solid::mesh2::{Mesh2, Patch};
use fs_solver::nonlinear::NewtonKrylovConfig;
use std::collections::BTreeSet;

// Shared by the CLI and focused tests; no CLI-shaped constitutive substitute.
fn solve(card:&MaterialCard, temperature_k:f64, stretch:f64, point_force_n:f64,
    gravity_m_s2:f64, grain_rad:f64, maximum_abs_strain:f64)
    -> Result<(MembranePatch, MembraneSolution)> {
    if ![temperature_k,stretch,point_force_n,gravity_m_s2,grain_rad,maximum_abs_strain]
        .iter().all(|v|v.is_finite()) || temperature_k<=0.0 || stretch<=1.0
        || point_force_n<0.0 || gravity_m_s2<0.0 || maximum_abs_strain<=0.0 {
        return Err("positive Kelvin temperature, stretch>1, nonnegative downward loads and positive strain bound required".into());
    }
    // Mechanical-only use must not require Cp, conductivity or bending data.
    let sheet=resolve_sheet(card,&query(temperature_k)?,
        SheetCapabilities {membrane:true,bending:false},MaterialPropertySelection::SingleClaimOnly)?;
    let mesh=Mesh2::triangles(0.1,0.1,2,2);
    let boundary:BTreeSet<_>=[Patch::Left,Patch::Right,Patch::Bottom,Patch::Top]
        .into_iter().flat_map(|p|mesh.patch_nodes(p)).collect();
    let mut prescribed=Vec::new();
    for &node in &boundary {
        let [x,y]=mesh.nodes[node];
        prescribed.extend([(3*node,(stretch-1.0)*(x-0.05)),
            (3*node+1,(stretch-1.0)*(y-0.05)),(3*node+2,0.0)]);
    }
    let patch=MembranePatch::try_new(mesh,vec![sheet;8],vec![grain_rad;8],
        MembraneDomain {minimum_area_ratio:0.01,maximum_abs_strain},
        MembraneBudget {max_nodes:9,max_triangles:8},||false)?;
    let mut loads=patch.body_load_n([0.0,0.0,-gravity_m_s2])?;
    loads[4][2]-=point_force_n;
    let solution=patch.solve(&loads,&prescribed,&[[0.0;3];9],MembraneSolveSettings {
        newton:NewtonKrylovConfig {absolute_tolerance:1e-9,relative_tolerance:1e-10,
            linear_restart:12,..NewtonKrylovConfig::default()},
        max_iterations:80,max_evaluations:2000,max_krylov_restart:12,
    },||false)?;
    Ok((patch,solution))
}

pub(super) fn run(card:&MaterialCard, temperature_k:f64, stretch:f64, point_force_n:f64,
    gravity_m_s2:f64, grain_rad:f64, maximum_abs_strain:f64) -> Result<()> {
    let (patch,result)=solve(card,temperature_k,stretch,point_force_n,gravity_m_s2,grain_rad,maximum_abs_strain)?;
    println!("# model=source-bound P1 finite-strain membrane; 2x2 integration demonstration, not mesh-converged cloth");
    println!("# material_state={:?}; reference_mass_kg={:.12e}",patch.sheets()[0].material().identity(),patch.reference_mass_kg());
    println!("# stretch={:.12e}; point_force_N={:.12e}; gravity_m_s2={:.12e}; grain_rad={:.12e}",stretch,point_force_n,gravity_m_s2,grain_rad);
    println!("# strain_energy_J={:.12e}; potential_energy_J={:.12e}; free_residual_N={:.12e}; newton_iterations={}; evaluations={}",
        result.strain_energy_j,result.potential_energy_j,result.newton.residual_norm,result.newton.iterations,result.evaluations);
    println!("node,reference_x_m,reference_y_m,ux_m,uy_m,uz_m,residual_x_N,residual_y_N,residual_z_N");
    for (i,point) in patch.mesh().nodes.iter().enumerate() {
        let u=result.displacement_m[i];let r=result.residual_n[i];
        println!("{i},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e},{:.12e}",point[0],point[1],u[0],u[1],u[2],r[0],r[1],r[2]);
    }
    Ok(())
}

#[cfg(test)]
mod tests {
    use super::*;
    use fs_evidence::ValidityDomain;
    use fs_matdb::{ClaimSet,InterpolationPolicy,MaterialStateId,PropertyClaim,PropertyKey,
        PropertyValue,Provenance,UncertaintyModel};
    use fs_qty::Dims;

    fn card(mass:f64, stiffness:f64)->MaterialCard {
        let mut claims=ClaimSet::new();
        for (key,value,dims) in [
            ("areal_density",mass,Dims([-2,1,0,0,0,0])),("thickness",0.001,Dims([1,0,0,0,0,0])),
            ("warp_stiffness",stiffness,Dims([0,1,-2,0,0,0])),("weft_stiffness",stiffness,Dims([0,1,-2,0,0,0])),
            ("poisson_warp_weft",0.0,Dims::NONE),("shear_stiffness",0.5*stiffness,Dims([0,1,-2,0,0,0])),
        ] {
            claims.insert_claim(PropertyClaim {key:PropertyKey::new(key,dims),value:PropertyValue::Scalar {value,dims},
                interpolation:InterpolationPolicy::ConstantWithinValidity,
                validity:ValidityDomain::unconstrained().with("T",290.0,310.0).with("ftx_characterization",0.0,0.0),
                uncertainty:UncertaintyModel::Unstated,observations:vec![],provenance:Provenance {
                    source:"synthetic CLI membrane fixture, not measured fabric".into(),license:"test".into(),artifact:None}}).unwrap();
        }
        MaterialCard::assemble(MaterialStateId {chemistry:"synthetic membrane".into(),phase:"solid".into(),
            process:"fixed test construction".into(),revision:1},claims,vec![]).unwrap()
    }

    #[test]
    fn mechanical_only_sources_drive_the_actual_scene_and_independent_gravity() {
        let (base,a)=solve(&card(0.2,8000.0),300.0,1.02,0.0,9.81,0.3,0.1).unwrap();
        let (heavy,b)=solve(&card(0.6,8000.0),300.0,1.02,0.0,9.81,0.3,0.1).unwrap();
        let (_,stiff)=solve(&card(0.2,16000.0),300.0,1.02,0.0,9.81,0.3,0.1).unwrap();
        assert!(a.displacement_m[4][2]<0.0);
        assert!(b.displacement_m[4][2]<a.displacement_m[4][2]);
        assert!(stiff.displacement_m[4][2]>a.displacement_m[4][2]);
        assert!(base.sheets()[0].bending().is_err());
        assert!(base.sheets()[0].material().property("specific_heat_capacity").is_none());
        assert!((heavy.reference_mass_kg()-3.0*base.reference_mass_kg()).abs()<1e-15);
        let reaction:f64=a.residual_n.iter().map(|r|r[2]).sum();
        assert!((reaction-base.reference_mass_kg()*9.81).abs()<1e-9);
    }

    #[test]
    fn unsupported_temperature_and_unprestretched_cli_inputs_refuse() {
        let c=card(0.2,8000.0);
        assert!(solve(&c,350.0,1.02,0.05,9.81,0.3,0.1).is_err());
        assert!(solve(&c,300.0,1.0,0.05,9.81,0.3,0.1).is_err());
        assert!(solve(&c,300.0,1.02,0.05,9.81,f64::NAN,0.1).is_err());
    }
}