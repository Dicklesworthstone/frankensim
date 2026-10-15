//! Source-bound finite-deformation membrane regressions, not empirical fabrics.
use fs_evidence::ValidityDomain;
use fs_material::sheet::{AREAL_DENSITY, MEMBRANE_DIMS, ResolvedSheet, SheetCapabilities, THICKNESS, resolve_sheet};
use fs_material::state_point::MaterialPropertySelection;
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_qty::Dims;
use fs_solid::membrane::{MembraneBudget, MembraneDomain, MembraneError, MembranePatch, MembraneSolveSettings};
use fs_solid::mesh2::Mesh2;
use fs_solver::nonlinear::NewtonKrylovConfig;

fn sheet(warp: f64, weft: f64, shear: f64, mass: f64) -> ResolvedSheet {
    let mut claims = ClaimSet::new();
    for (key, value, dims) in [
        (AREAL_DENSITY,mass,Dims([-2,1,0,0,0,0])),
        (THICKNESS,0.001,Dims([1,0,0,0,0,0])),
        ("warp_stiffness",warp,MEMBRANE_DIMS), ("weft_stiffness",weft,MEMBRANE_DIMS),
        ("poisson_warp_weft",0.0,Dims::NONE), ("shear_stiffness",shear,MEMBRANE_DIMS),
    ] {
        claims.insert_claim(PropertyClaim {
            key:PropertyKey::new(key,dims), value:PropertyValue::Scalar {value,dims},
            interpolation:InterpolationPolicy::ConstantWithinValidity,
            validity:ValidityDomain::unconstrained().with("T",290.0,310.0),
            uncertainty:UncertaintyModel::Unstated, observations:vec![],
            provenance:Provenance {source:"synthetic membrane test; not a measured textile".into(),license:"test".into(),artifact:None},
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId {chemistry:"synthetic sheet".into(),
        phase:"solid".into(),process:"test".into(),revision:1},claims,vec![]).unwrap();
    resolve_sheet(&card,&QueryPoint::new().with("T",300.0).unwrap(),
        SheetCapabilities {membrane:true,bending:false},MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn mesh() -> Mesh2 {
    Mesh2 { nodes:vec![[0.0,0.0],[0.1,0.0],[0.0,0.1],[-0.1,0.0],[0.0,-0.1]],
        elems:vec![vec![0,1,2],vec![0,2,3],vec![0,3,4],vec![0,4,1]],patches:vec![] }
}
fn domain() -> MembraneDomain { MembraneDomain {minimum_area_ratio:0.01,maximum_abs_strain:0.2} }
fn budget() -> MembraneBudget { MembraneBudget {max_nodes:32,max_triangles:32} }
fn patch(warp:f64,weft:f64,shear:f64,mass:f64,grain:f64) -> MembranePatch {
    MembranePatch::try_new(mesh(),vec![sheet(warp,weft,shear,mass);4],vec![grain;4],domain(),budget(),||false).unwrap()
}
fn settings() -> MembraneSolveSettings {
    MembraneSolveSettings { newton:NewtonKrylovConfig {absolute_tolerance:1e-9,
        relative_tolerance:1e-11,linear_restart:8,..NewtonKrylovConfig::default()},
        max_iterations:30,max_evaluations:1000,max_krylov_restart:8 }
}
fn close(a:f64,b:f64,atol:f64,rtol:f64) {
    assert!((a-b).abs() <= atol+rtol*a.abs().max(b.abs()),"{a:e} != {b:e}");
}
fn prescribed(stretch:f64) -> Vec<(usize,f64)> {
    let mut result = vec![(0,0.0),(1,0.0)];
    for (node,p) in mesh().nodes.iter().enumerate().skip(1) {
        result.extend([(3*node,(stretch-1.0)*p[0]),(3*node+1,(stretch-1.0)*p[1]),(3*node+2,0.0)]);
    }
    result
}
fn displaced() -> Vec<[f64;3]> {
    vec![[0.0001,-0.0002,-0.006],[0.003,0.0001,0.0],[0.0,0.001,0.0003],
        [-0.001,0.0002,-0.0001],[0.0001,-0.002,0.0002]]
}

#[test]
fn finite_rigid_motion_has_no_membrane_strain_or_force() {
    let p=patch(8000.0,2000.0,500.0,0.2,0.37);
    let (s,c)=0.61_f64.sin_cos();
    // An out-of-plane rotation and arbitrary rigid translation, not a linear
    // infinitesimal rotation test (which would not expose wrong Green strain).
    let u:Vec<_>=p.mesh().nodes.iter().map(|x|[0.27,c*x[1]-x[1]-0.13,s*x[1]+0.18]).collect();
    let e=p.evaluate(&u,||false).unwrap();
    close(e.strain_energy_j,0.0,1e-22,0.0);
    for f in e.internal_gradient_n.iter().flatten() { close(*f,0.0,1e-9,0.0); }
    let zero=p.tangent_apply(&u,&vec![[0.7,-0.3,0.4];5],||false).unwrap();
    assert!(zero.iter().flatten().all(|f|*f==0.0));
}

#[test]
fn energy_force_and_matrix_free_tangent_match_finite_differences() {
    let p=patch(8000.0,2000.0,500.0,0.2,0.37);
    let u=displaced(); let base=p.evaluate(&u,||false).unwrap(); let h=1e-7;
    for dof in 0..15 {
        let mut lo=u.clone();let mut hi=u.clone();lo[dof/3][dof%3]-=h;hi[dof/3][dof%3]+=h;
        let a=p.evaluate(&lo,||false).unwrap();let b=p.evaluate(&hi,||false).unwrap();
        close((b.strain_energy_j-a.strain_energy_j)/(2.0*h),base.internal_gradient_n[dof/3][dof%3],1e-6,2e-6);
        let mut v=vec![[0.0;3];5];v[dof/3][dof%3]=1.0;
        let jv=p.tangent_apply(&u,&v,||false).unwrap();
        for node in 0..5 { for axis in 0..3 {
            close((b.internal_gradient_n[node][axis]-a.internal_gradient_n[node][axis])/(2.0*h),jv[node][axis],2e-5,3e-6);
        } }
    }
    // Hessian reciprocity at nonzero stress, including the geometric term.
    let v=vec![[0.4,-0.2,0.7],[0.0,0.1,-0.3],[0.2,0.5,0.0],[-0.1,0.2,0.1],[0.3,0.0,-0.1]];
    let w=vec![[0.1,0.3,-0.1],[-0.2,0.1,0.4],[0.3,0.0,0.1],[0.0,0.4,0.5],[0.2,-0.1,0.0]];
    let jv=p.tangent_apply(&u,&v,||false).unwrap();let jw=p.tangent_apply(&u,&w,||false).unwrap();
    let inner=|a:&Vec<[f64;3]>,b:&Vec<[f64;3]>|a.iter().flatten().zip(b.iter().flatten()).map(|(a,b)|a*b).sum::<f64>();
    close(inner(&v,&jw),inner(&w,&jv),1e-9,1e-12);
}

#[test]
fn grain_and_source_substitutions_change_the_right_physics_only() {
    let a=patch(8000.0,2000.0,500.0,0.2,0.0);
    let b=patch(2000.0,8000.0,500.0,0.2,std::f64::consts::FRAC_PI_2);
    let heavy=patch(8000.0,2000.0,500.0,0.6,0.0);
    let u=displaced();let ra=a.evaluate(&u,||false).unwrap();let rb=b.evaluate(&u,||false).unwrap();
    close(ra.strain_energy_j,rb.strain_energy_j,1e-12,1e-12);
    close(ra.strain_energy_j,heavy.evaluate(&u,||false).unwrap().strain_energy_j,0.0,0.0);
    for (a,b) in ra.internal_gradient_n.iter().flatten().zip(rb.internal_gradient_n.iter().flatten()) {
        close(*a,*b,1e-10,1e-11);
    }
    close(a.reference_mass_kg(),0.004,1e-16,1e-14);
    close(heavy.reference_mass_kg(),3.0*a.reference_mass_kg(),1e-16,1e-14);
    let f=a.body_load_n([0.0,0.0,-9.81]).unwrap();
    close(f.iter().map(|v|v[2]).sum::<f64>(),-9.81*a.reference_mass_kg(),1e-15,1e-13);
    assert_eq!(a.sheets()[0],sheet(8000.0,2000.0,500.0,0.2));
    assert_ne!(a.sheets()[0].material().identity(),heavy.sheets()[0].material().identity());
    assert!(a.sheets()[0].bending().is_err()); // no manufactured bending capability
}

#[test]
fn geometric_stiffness_carries_physical_prestretch() {
    let p=patch(8000.0,8000.0,4000.0,0.2,0.0);
    let zero=vec![[0.0;3];5];let mut v=zero.clone();v[0][2]=1.0;
    close(p.tangent_apply(&zero,&v,||false).unwrap()[0][2],0.0,1e-9,0.0);
    let mut u=zero;
    for (dof,value) in prescribed(1.02) {u[dof/3][dof%3]=value;}
    let expected=2.0*8000.0*(1.02_f64.powi(2)-1.0);
    close(p.tangent_apply(&u,&v,||false).unwrap()[0][2],expected,1e-8,1e-12);
}

#[test]
fn loaded_sheet_equilibrium_matches_the_exact_pyramid_energy() {
    let p=patch(8000.0,8000.0,4000.0,0.2,0.0);
    let h=0.005_f64;let stretch=1.02_f64;let side=0.1_f64;
    // Four P1 triangles: dU/dh = 2 E (lambda²-1) h + 4 E h³ / L².
    let load=2.0*8000.0*(stretch.powi(2)-1.0)*h+4.0*8000.0*h.powi(3)/side.powi(2);
    let mut loads=vec![[0.0;3];5];loads[0][2]=-load;
    let fixed=prescribed(stretch);
    let result=p.solve(&loads,&fixed,&vec![[0.0;3];5],settings(),||false).unwrap();
    assert!(result.newton.converged);
    close(result.displacement_m[0][2],-h,1e-10,1e-10);
    close(result.residual_n[0][2],0.0,1e-9,0.0);
    close(result.residual_n.iter().skip(1).map(|r|r[2]).sum::<f64>(),load,1e-9,1e-10);
    let s=0.5*(stretch.powi(2)-1.0);
    let energy=2.0*8000.0*side.powi(2)*s*s+2.0*8000.0*s*h*h+8000.0*h.powi(4)/side.powi(2);
    close(result.strain_energy_j,energy,1e-10,1e-10);
    for (dof,value) in fixed { assert_eq!(result.displacement_m[dof/3][dof%3],value); }
    let repeat=p.solve(&loads,&prescribed(stretch),&vec![[0.0;3];5],settings(),||false).unwrap();
    assert_eq!(result.displacement_m,repeat.displacement_m);
    assert_eq!(result.newton,repeat.newton);
    let stiff=patch(16000.0,16000.0,8000.0,0.2,0.0);
    let changed=stiff.solve(&loads,&prescribed(stretch),&vec![[0.0;3];5],settings(),||false).unwrap();
    assert!(changed.displacement_m[0][2].abs()<h);
    assert_eq!(p.reference_mass_kg(),stiff.reference_mass_kg());
}

#[test]
fn topology_domain_cancellation_and_work_refusals_publish_no_solution() {
    let p=patch(8000.0,8000.0,4000.0,0.2,0.0);let zero=vec![[0.0;3];5];
    let mut collapsed=zero.clone();for (i,x) in p.mesh().nodes.iter().enumerate(){collapsed[i]=[-x[0],-x[1],0.0];}
    assert!(matches!(p.evaluate(&collapsed,||false),Err(MembraneError::State {..})));
    assert!(matches!(p.evaluate(&zero,||true),Err(MembraneError::Cancelled)));
    let mut polls=0;
    assert!(matches!(p.tangent_apply(&zero,&zero,||{polls+=1;polls==3}),Err(MembraneError::Cancelled)));
    let mut duplicate=mesh();duplicate.elems.push(duplicate.elems[0].clone());
    assert!(MembranePatch::try_new(duplicate,vec![sheet(8000.0,8000.0,4000.0,0.2);5],vec![0.0;5],domain(),budget(),||false).is_err());
    assert!(matches!(MembranePatch::try_new(mesh(),p.sheets().to_vec(),vec![0.0;4],domain(),
        MembraneBudget {max_nodes:4,max_triangles:4},||false),Err(MembraneError::Budget(_))));
    let mut loads=zero.clone();loads[0][2]=-3.632;
    let mut unfunded=settings();unfunded.max_evaluations=1;
    assert!(matches!(p.solve(&loads,&prescribed(1.02),&zero,unfunded,||false),Err(MembraneError::Budget(_))));
    let mut polls=0;
    assert!(matches!(p.solve(&loads,&prescribed(1.02),&zero,settings(),||{polls+=1;polls==30}),Err(MembraneError::Cancelled)));
    assert!(p.solve(&loads,&[(0,0.0),(0,0.0)],&zero,settings(),||false).is_err());
    // No artificial bending floor turns an unstressed flat transverse tangent
    // into a successful solve. Prestretch above is supplied by actual supports.
    assert!(matches!(p.solve(&loads,&prescribed(1.0),&zero,settings(),||false),Err(MembraneError::Unconverged(_))));
    assert!(p.solve(&loads,&prescribed(1.02),&zero,settings(),||false).is_ok());
}

#[path="membrane/vibration.rs"]
mod vibration;
