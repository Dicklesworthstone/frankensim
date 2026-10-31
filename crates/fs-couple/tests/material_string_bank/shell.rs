//! Source-backed synthetic shell/string coupling; no measured instrument claim.
use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::{
    plate::SharedPlateBridge, shell::SharedMaterialShellBody,
};
use fs_couple::shear_shell::{compile_material_shear_shell, ShellProjectionBudget,
    ShellRegionMaterial, TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES};
use fs_couple::thin_plate::{PlateMaterialModel, PlateRegionMaterial, PlateThicknessConstraint};
use fs_material::plane_stress_viscosity::PLANE_STRESS_VISCOSITY_PROPERTIES;
use fs_plate::{PlateMesh, PlateRegion};
use fs_plate::shell::{ShellMesh, mitc::{ConstraintFrame, MitcBudget, NodalDirectors,
    ShellConstraint, ShellMass}};

fn shell_material(rho:f64, eta:f64) -> ResolvedMaterialStatePoint {
    shell_material_with_loss_anchor(rho, eta, None)
}
fn shell_material_with_loss_anchor(rho:f64, eta:f64, anchor:Option<f64>) -> ResolvedMaterialStatePoint {
    let mut entries=vec![("density",Density::DIMS,rho),
        ("young_modulus_1",Pressure::DIMS,1.0e9),("young_modulus_2",Pressure::DIMS,4.0e8),
        ("poisson_ratio_12",Dims::NONE,0.25),("shear_modulus_12",Pressure::DIMS,1.8e8),
        ("shear_modulus_13",Pressure::DIMS,8.0e7),("shear_modulus_23",Pressure::DIMS,3.0e7)];
    entries.extend(PLANE_STRESS_VISCOSITY_PROPERTIES.into_iter()
        .zip([eta,0.7*eta,0.3*eta,0.2*eta,0.04*eta,-0.03*eta])
        .map(|(key,value)|(key,DynViscosity::DIMS,value)));
    entries.extend(TRANSVERSE_SHEAR_VISCOSITY_PROPERTIES.into_iter()
        .zip([0.4*eta,0.8*eta,0.05*eta]).map(|(key,value)|(key,DynViscosity::DIMS,value)));
    if let Some(anchor) = anchor {
        entries.push((fs_couple::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA,
            fs_couple::shear_shell::CONSTANT_LOSS_REFERENCE_OMEGA_DIMS, anchor));
    }
    let mut claims=ClaimSet::new(); let mut requirements=Vec::new();
    for (key,dims,value) in entries {
        claims.insert_claim(PropertyClaim { key:PropertyKey::new(key,dims),
            value:PropertyValue::Scalar {value,dims},
            validity:ValidityDomain::unconstrained().with("T",280.0,310.0).with("omega",0.0,100_000.0),
            uncertainty:UncertaintyModel::Unstated,interpolation:InterpolationPolicy::ConstantWithinValidity,
            provenance:Provenance {source:"synthetic curved material shell".into(),license:"CC0-1.0".into(),artifact:None},
            observations:vec![] }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(key,dims,ScalarAdmissibility::Finite).unwrap());
    }
    let card=MaterialCard::assemble(MaterialStateId {chemistry:"synthetic shell".into(),phase:"solid".into(),
        process:"fixture".into(),revision:0},claims,vec![]).unwrap();
    resolve_material_state_point(&card,&QueryPoint::new().with("T",293.15).unwrap().with("omega",1000.0).unwrap(),
        &requirements,MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn geometry(rise:f64) -> (ShellMesh,NodalDirectors,Vec<ShellConstraint>) {
    let flat=PlateMesh::rectangle(0.4,0.3,3,3);
    let mut normals=Vec::new();
    let nodes=flat.nodes.iter().map(|&(x,y)| {
        let u=x/0.4; let v=y/0.3;
        let z=16.0*rise*u*(1.0-u)*v*(1.0-v);
        let nx=-16.0*rise*(1.0-2.0*u)*v*(1.0-v)/0.4;
        let ny=-16.0*rise*u*(1.0-u)*(1.0-2.0*v)/0.3;
        let norm=(1.0+nx*nx+ny*ny).sqrt();
        normals.push([nx/norm,ny/norm,1.0/norm]); [x,y,z]
    }).collect();
    let mesh=ShellMesh::new(nodes,flat.tris).unwrap();
    let directors=NodalDirectors::new(&mesh,normals,0.8).unwrap();
    let constraints=PlateMesh::rectangle_boundary(3,3).into_iter().map(|node|ShellConstraint {
        node,mask:[true;6],frame:ConstraintFrame::Global,
    }).collect();
    (mesh,directors,constraints)
}
fn shell(rho:f64,rise:f64,angle:f64,gate:&CancelGate) -> SharedMaterialShellBody {
    let (mesh,directors,constraints)=geometry(rise);
    let nodes=mesh.nodes.len();
    let budget=MitcBudget {max_nodes:nodes,max_triangles:mesh.tris.len()};
    let region=ShellRegionMaterial {region:PlateRegionMaterial {
        region:PlateRegion {name:"wood".into(),triangle_indices:(0..mesh.tris.len()).collect()},
        material:shell_material(rho,1000.0),model:PlateMaterialModel::Orthotropic12,
        thickness_constraint:PlateThicknessConstraint::FixedThickness(0.003),material_angle_rad:angle,
    },shear_correction:5.0/6.0,axis1_global:[1.0,0.0,0.0]};
    let specimen=compile_material_shear_shell(mesh,vec![region],budget,gate).unwrap();
    let force_ports=[5usize,10].map(|node| {
        let mut p=vec![0.0;6*nodes];p[6*node+2]=1.0;p
    }).to_vec();
    SharedMaterialShellBody {specimen,directors,constraints,mass:ShellMass::Consistent,
        drilling_alpha:1e-3,eigenvalue_window:(0.0,1e9),n_modes:3,force_ports,
        bridges:vec![SharedPlateBridge {station_fraction:0.8,stiffness_n_m:3000.0,
            damping_n_s_m:0.1,rest_extension_m:0.0};2],
        assembly_budget:budget,projection_budget:ShellProjectionBudget {max_modes:8,
            max_projection_terms:32_000_000,reference_resistance_n_s_m:1.0},
        modal_budget:ModalAcousticTimeBudget::audible_reference(),listener_m:1.0}
}
fn build_shell(cx:&Cx<'_>,gate:&CancelGate,rho:f64,rise:f64,angle:f64) -> ReciprocalStringBank {
    let mut o=options(0);o.network.max_connections=16;o.network.max_setup_terms=131_072;
    ReciprocalStringBank::from_material_shear_shell(cx,vec![member(7800.0).source,layered().source],
        gas(101_325.0),shell(rho,rise,angle,gate),None,law(),o,gate).unwrap()
}

#[test]
fn g3_curvature_density_and_global_grain_change_the_shared_physical_body() {
    with_cx(|cx,gate| {
        let mut curved=build_shell(cx,gate,700.0,0.02,0.17);
        assert_eq!(curved.source().bridge_surface_materials(1).unwrap(),
            [material(15_600.0,100e9).material().clone(), shell_material(700.0,1000.0).material().clone()]);
        let dense=build_shell(cx,gate,1400.0,0.02,0.17);
        let flat=build_shell(cx,gate,700.0,0.0,0.17);
        let cross=build_shell(cx,gate,700.0,0.02,0.91);
        let frequency=|bank:&ReciprocalStringBank|bank.system().components()[2].modes()[0].angular_frequency_rad_s;
        let w=frequency(&curved);
        assert!((frequency(&dense)/w-0.5_f64.sqrt()).abs()<2e-6);
        assert!((frequency(&flat)/w-1.0).abs()>1e-3);
        assert!((frequency(&cross)/w-1.0).abs()>1e-4);
        assert_eq!(curved.cards()[0].fundamental_hz(),dense.cards()[0].fundamental_hz());
        let image=curved.source().material_shell_image().unwrap();
        assert_eq!(image.source().regions()[0].input().axis1_global,[1.0,0.0,0.0]);
        assert_eq!(image.source().regions()[0].input().region.material.identity(),shell_material(700.0,1000.0).identity());
        let matrix=image.damping().matrix_per_s();
        assert!((0..3).any(|i|(0..3).any(|j|i!=j && matrix[3*i+j].abs()>1e-8)));
        assert_eq!(image.damping().channel_count(),3);
        assert_eq!(curved.source().body_connections().len(),3);
        let weights=curved.system().acceleration_pressure_weights().unwrap();
        assert!(weights[..8].iter().all(|w|*w==0.0));
        for (i,&volume) in image.volume_shapes_m2_per_sqrt_kg().iter().enumerate() {
            let expected=gas(101_325.0).density*volume/core::f64::consts::TAU;
            assert!((weights[8+i]-expected).abs()<1e-12*expected.abs().max(1.0));
            for (port,member) in image.force_ports().iter().zip(curved.source().strings()) {
                assert_eq!(port[i],member.connections[0].body.shapes[i]);
            }
        }
        let zero=vec![0.0;curved.system().mode_count()];
        let mut peak=0.0_f64;
        for _ in 0..96 {
            let f=curved.step(&zero,drive()).unwrap();
            assert!(f.energy_residual_j.abs()<=f.energy_tolerance_j);
            peak=peak.max(f.observer_pressure_pa.abs());
        }
        assert!(peak>1e-12);
        assert!(curved.system().components()[1].states().iter().any(|s|s.velocity_m_sqrt_kg_per_s.abs()>1e-12));
    });
}
#[test]
fn g0_source_rebind_and_geometry_admission_preserve_data_or_refuse() {
    with_cx(|_,gate| {
        let body=shell(700.0,0.02,0.17,gate);
        let original=&body.specimen;
        let changed=original.with_region_material("wood",&shell_material(1400.0,2000.0),gate).unwrap();
        assert!((changed.mass_kg()/original.mass_kg()-2.0).abs()<1e-12);
        assert_eq!(changed.mesh().nodes,original.mesh().nodes);
        assert_eq!(changed.mesh().tris,original.mesh().tris);
        assert_eq!(changed.regions()[0].input().axis1_global,[1.0,0.0,0.0]);
        assert!(original.with_region_material("wood",&material(700.0,1e9),gate).is_err());
        let mut inputs=vec![original.regions()[0].input().clone()];
        inputs[0].axis1_global=[0.0;3];
        assert!(compile_material_shear_shell(original.mesh().clone(),inputs,body.assembly_budget,gate).is_err());
        let mut mesh=original.mesh().clone();mesh.tris[0].swap(0,1);
        assert!(compile_material_shear_shell(mesh,vec![original.regions()[0].input().clone()],
            body.assembly_budget,gate).is_err());
        let mut inputs=vec![original.regions()[0].input().clone()];
        inputs[0].region.region.triangle_indices.pop();
        assert!(compile_material_shear_shell(original.mesh().clone(),inputs,body.assembly_budget,gate).is_err());
    });
}
#[test]
fn g0_shell_budget_detached_force_and_cancellation_refuse_before_playback() {
    with_cx(|cx,gate| {
        for case in 0..3 {
            let mut body=shell(700.0,0.02,0.17,gate);
            match case {
                0=>body.projection_budget.max_projection_terms=1,
                1=>body.force_ports[0].fill(0.0),
                _=>body.force_ports.pop().map(|_|()).unwrap(),
            }
            assert!(ReciprocalStringBank::from_material_shear_shell(cx,vec![member(7800.0).source,layered().source],
                gas(101_325.0),body,None,law(),options(0),gate).is_err());
        }
        let body=shell(700.0,0.02,0.17,gate);gate.request();
        assert!(ReciprocalStringBank::from_material_shear_shell(cx,vec![member(7800.0).source,layered().source],
            gas(101_325.0),body,None,law(),options(0),gate).is_err());
    });
}

#[path = "shell/connections.rs"]
mod connections;

#[path = "shell/pressure.rs"]
mod pressure;

#[path = "shell/assembly.rs"]
mod assembly;

#[path = "shell/bem.rs"]
mod bem;

#[path = "shell/cavity.rs"]
mod cavity;

#[path = "shell/coating.rs"]
mod coating;

#[path = "shell/rotation.rs"]
mod rotation;
