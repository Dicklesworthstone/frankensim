//! G0/G1 laminate-to-shell checks with independently integrated section work.
//! These are numerical mechanics checks, not varnish specimen validation.

use fs_material::elastic::OrthotropicElastic;
use fs_material::laminate::{LaminaLayer, LaminateSection, LayerLaw, ReferenceSurface};
use fs_material::sheet::SheetMatrix;
use fs_plate::shell::mitc::{
    ConstraintFrame, MitcBudget, MitcShellModel, NodalDirectors, ShearShellSection,
    ShellConstraint, ShellMass, assemble_mitc3p, solve_static,
};
use fs_plate::shell::ShellMesh;

fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a-b).abs() <= tolerance*a.abs().max(b.abs()).max(1e-30), "{a:e} != {b:e}");
}

fn layer(e: f64, h: f64, rho: f64) -> LaminaLayer {
    LaminaLayer {
        thickness_m: h, orientation_rad: 0.0, density_kg_m3: rho, material: None,
        law: LayerLaw::Elastic(OrthotropicElastic::new([e;3], [0.0;3], [0.5*e;3], 0.01).unwrap()),
    }
}

fn laminate(reverse: bool, reference: ReferenceSurface) -> LaminateSection {
    let mut layers = vec![layer(2.0e8, 0.01, 400.0), layer(8.0e8, 0.004, 1200.0)];
    if reverse { layers.reverse(); }
    LaminateSection::try_new(layers, reference).unwrap()
}

fn mesh(n: usize, width: f64) -> (ShellMesh, NodalDirectors) {
    let mut points = Vec::new();
    for j in 0..=n { for i in 0..=n {
        points.push([i as f64/n as f64, width*j as f64/n as f64, 0.0]);
    } }
    let mut triangles = Vec::new();
    for j in 0..n { for i in 0..n {
        let a = j*(n+1)+i; let b = a+1; let c = a+n+1; let d = c+1;
        triangles.extend([[a,b,d], [a,d,c]]);
    } }
    let mesh = ShellMesh::new(points, triangles).unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0,0.0,1.0];mesh.nodes.len()], 0.9).unwrap();
    (mesh, directors)
}

fn assemble(mesh: &ShellMesh, directors: &NodalDirectors, section: ShearShellSection,
    constraints: &[ShellConstraint], mass: ShellMass) -> MitcShellModel
{
    let sections: Vec<_> = (0..mesh.tris.len()).map(|e| {
        let g = mesh.facet(e).unwrap();
        section.with_material_angle_rad(g.frame[1][0].atan2(g.frame[0][0])).unwrap()
    }).collect();
    assemble_mitc3p(mesh, directors, &sections, constraints, mass,
        MitcBudget { max_nodes:mesh.nodes.len(), max_triangles:mesh.tris.len() },
        1e-3, ||false).unwrap()
}

fn mass_work(model: &MitcShellModel, velocity: &[f64]) -> f64 {
    let q = model.project_load(velocity).unwrap();
    let mut action = vec![0.0;q.len()];
    model.pencil().m.spmv(&q, &mut action);
    q.iter().zip(action).map(|(v,mv)| v*mv).sum()
}

#[test]
fn g0_single_layer_flat_shell_matches_homogeneous_owner() {
    let (e,h,rho) = (2.0e8,0.01,400.0);
    let source = LaminateSection::try_new(vec![layer(e,h,rho)], ReferenceSurface::Midsurface).unwrap();
    let laminate = ShearShellSection::from_laminate(&source,0.0).unwrap();
    let homogeneous = ShearShellSection::from_orthotropic_solid(e,e,0.0,0.5*e,0.5*e,0.5*e,
        h,rho,0.0,5.0/6.0).unwrap();
    let (mesh,directors) = mesh(1,0.4);
    for mass in [ShellMass::Consistent, ShellMass::Lumped] {
        let a = assemble(&mesh,&directors,laminate,&[],mass);
        let b = assemble(&mesh,&directors,homogeneous,&[],mass);
        for (left,right) in [(&a.pencil().k,&b.pencil().k),(&a.pencil().m,&b.pencil().m)] {
            let scale = (0..a.pencil().free).map(|i| left.get(i,i).abs()).fold(0.0_f64,f64::max);
            for i in 0..a.pencil().free { for j in 0..a.pencil().free {
                assert!((left.get(i,j)-right.get(i,j)).abs() < 2e-12*scale);
            } }
        }
    }
}

#[test]
fn g1_constant_strain_curvature_keeps_signed_b_work() {
    let (mesh,directors) = mesh(3,0.4);
    let strain = [0.001,-0.0002,0.0003];
    let curvature = [0.03,-0.01,0.02];
    let field: Vec<_> = mesh.nodes.iter().flat_map(|p| {
        let [x,y,_] = *p;
        [strain[0]*x+0.5*strain[2]*y, strain[1]*y+0.5*strain[2]*x,
            -0.5*curvature[0]*x*x-0.5*curvature[1]*y*y-0.5*curvature[2]*x*y,
            -curvature[1]*y-0.5*curvature[2]*x,
            curvature[0]*x+0.5*curvature[2]*y, 0.0]
    }).collect();
    let mut energies = Vec::new();
    for reverse in [false,true] {
        let source = laminate(reverse,ReferenceSurface::Midsurface);
        let model = assemble(&mesh,&directors,ShearShellSection::from_laminate(&source,0.0).unwrap(),
            &[],ShellMass::Consistent);
        let mut expected = 0.0;
        for i in 0..3 { for j in 0..3 {
            expected += strain[i]*source.a_pa_m()[3*i+j]*strain[j]
                + 2.0*strain[i]*source.b_pa_m2()[3*i+j]*curvature[j]
                + curvature[i]*source.d_pa_m3()[3*i+j]*curvature[j];
        } }
        let observed = model.strain_energy(&field).unwrap();
        near(observed,0.5*0.4*expected,3e-11);
        energies.push(observed);
    }
    // Reversing the laminate leaves A and D fixed but reverses the cross work.
    assert!((energies[0]-energies[1]).abs() > 1e-3*energies[0]);
}

#[test]
fn g1_first_mass_moment_matches_continuum_rigid_motion() {
    let source = laminate(false,ReferenceSurface::Offset(0.005));
    let section = ShearShellSection::from_laminate(&source,0.0).unwrap();
    assert!(section.mass_first_moment_kg_m() > 0.0);
    let mesh = ShellMesh::new(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]],vec![[0,1,2]]).unwrap();
    let directors = NodalDirectors::new(&mesh,vec![[0.0,0.0,1.0];3],0.9).unwrap();
    let (speed,omega) = (0.13,0.7);
    let velocity: Vec<_> = mesh.nodes.iter().flat_map(|p|
        [speed,0.0,-omega*p[0],0.0,omega,0.0]).collect();
    for mass in [ShellMass::Consistent, ShellMass::Lumped] {
        let model = assemble(&mesh,&directors,section,&[],mass);
        let x2_integral = if mass == ShellMass::Consistent {1.0/12.0} else {1.0/6.0};
        let expected = source.i0_kg_m2()*(0.5*speed*speed+x2_integral*omega*omega)
            + source.i1_kg_m()*speed*omega + 0.5*source.i2_kg()*omega*omega;
        near(mass_work(&model,&velocity),expected,2e-12);
    }
}

#[test]
fn g1_reference_shift_preserves_physical_kinetic_work() {
    let (mesh,directors) = mesh(2,0.4);
    let (speed,omega) = (0.11,0.83);
    let mut values = Vec::new();
    for offset in [0.0,0.007,0.014] {
        let source = laminate(false,ReferenceSurface::Offset(offset));
        let section = ShearShellSection::from_laminate(&source,0.0).unwrap();
        let model = assemble(&mesh,&directors,section,&[],ShellMass::Consistent);
        // The same material velocity, written at three reference surfaces.
        let velocity: Vec<_> = mesh.nodes.iter().flat_map(|p|
            [speed+offset*omega,0.0,-omega*p[0],0.0,omega,0.0]).collect();
        values.push(mass_work(&model,&velocity));
    }
    near(values[0],values[1],2e-12);
    near(values[0],values[2],2e-12);
}

#[test]
fn g1_varying_directors_use_material_velocity_for_rotary_inertia() {
    let source = laminate(false,ReferenceSurface::Offset(0.005));
    let section = ShearShellSection::from_laminate(&source,0.0).unwrap();
    let mesh = ShellMesh::new(vec![[0.0,0.0,0.0],[1.0,0.0,0.0],[0.0,1.0,0.0]],vec![[0,1,2]]).unwrap();
    let tilt = 0.2_f64; let normal_z = (1.0-tilt*tilt).sqrt();
    let directors = NodalDirectors::new(&mesh,
        vec![[-tilt,0.0,normal_z],[tilt,0.0,normal_z],[-tilt,0.0,normal_z]],0.9).unwrap();
    let model = assemble(&mesh,&directors,section,&[],ShellMass::Consistent);
    let (speed,omega) = (0.13,0.7);
    let velocity: Vec<_> = mesh.nodes.iter().flat_map(|p|
        [speed,0.0,-omega*p[0],0.0,omega,0.0]).collect();
    // For rigid rotation about y, d(x) = [tilt*(2*x-1),0,normal_z].
    // Integrate |u_dot + z*(omega cross d)|^2 independently over the unit
    // right triangle. omega is perpendicular to every director, so the
    // numerical drilling inertia contributes exactly zero.
    let expected = source.i0_kg_m2()*(0.5*speed*speed+omega*omega/12.0)
        + source.i1_kg_m()*speed*omega*normal_z
        + source.i2_kg()*omega*omega*(0.5*normal_z*normal_z+tilt*tilt/6.0);
    near(mass_work(&model,&velocity),expected,2e-11);
}

#[test]
fn g1_axial_load_bends_unsymmetric_laminate_on_three_meshes() {
    let source = laminate(false,ReferenceSurface::Midsurface);
    let (a,b,d) = (source.a_pa_m()[0],source.b_pa_m2()[0],source.d_pa_m3()[0]);
    let resultant = 10.0;
    let curvature = -b*resultant/(a*d-b*b);
    assert!(curvature.abs() > 1e-5);
    for n in [2,4,8] {
        let (mesh,directors) = mesh(n,0.2);
        let constraints: Vec<_> = mesh.nodes.iter().enumerate().filter(|(_,p)| p[0]==0.0)
            .map(|(node,_)| ShellConstraint{node, mask:[true;6],frame:ConstraintFrame::Global}).collect();
        let section = ShearShellSection::from_laminate(&source,0.0).unwrap();
        let model = assemble(&mesh,&directors,section,&constraints,ShellMass::Consistent);
        let mut loads = vec![0.0;6*mesh.nodes.len()];
        for j in 0..=n {
            loads[6*(j*(n+1)+n)] = resultant*0.2/n as f64
                * if j==0 || j==n {0.5} else {1.0};
        }
        let u = solve_static(&model,&loads).unwrap();
        let tip_rotation: f64 = (0..=n).map(|j| u[6*(j*(n+1)+n)+4]).sum::<f64>()/(n+1) as f64;
        let tip_height: f64 = (0..=n).map(|j| u[6*(j*(n+1)+n)+2]).sum::<f64>()/(n+1) as f64;
        near(tip_rotation,curvature,0.005);
        near(tip_height,-0.5*curvature,0.005);
    }
}

#[test]
fn g0_coupling_admission_checks_full_stiffness_and_inertia() {
    let one = SheetMatrix::orthotropic(1.0,1.0,0.0,1.0).unwrap();
    let make = |b,i1| ShearShellSection::new_coupled(one,b,one,[1.0,0.0,1.0],1.0,1.0,i1,1.0,0.0);
    let diagonal = |v| [v,0.0,0.0,0.0,v,0.0,0.0,0.0,v];
    assert!(make(diagonal(0.9),0.9).is_ok());
    assert!(make(diagonal(-0.9),-0.9).is_ok());
    assert!(make(diagonal(1.0),0.0).is_err());
    assert!(make(diagonal(1.01),0.0).is_err());
    assert!(make(diagonal(0.0),1.0).is_err());
    assert!(make(diagonal(f64::NAN),0.0).is_err());
    let mut asymmetric = diagonal(0.1); asymmetric[1] = 0.1;
    assert!(make(asymmetric,0.0).is_err());
    assert!(SheetMatrix::from_coefficients(asymmetric).is_err());
}

#[test]
fn g1_off_axis_transverse_shear_preserves_full_tensor_and_material_axis_order() {
    let elastic = OrthotropicElastic::new([2.0e8;3], [0.0;3],
        [1.0e8,3.0e7,8.0e7],0.01).unwrap();
    let angle = 0.37_f64; let h = 0.01;
    let source = LaminateSection::try_new(vec![LaminaLayer {
        thickness_m:h, orientation_rad:angle, density_kg_m3:400.0, material:None,
        law:LayerLaw::Elastic(elastic),
    }],ReferenceSurface::Midsurface).unwrap();
    let [gx,gy] = [8.0e7,3.0e7]; // G31, G23 in the material owner's order.
    let (s,c) = angle.sin_cos(); let factor = h*5.0/6.0;
    let expected = [factor*(gx*c*c+gy*s*s),factor*(gx-gy)*c*s,
        factor*(gx-gy)*c*s,factor*(gx*s*s+gy*c*c)];
    for (actual,expected) in source.shear_pa_m().iter().zip(expected) {
        near(*actual,expected,3e-13);
    }
    assert!(source.shear_pa_m()[1].abs() > 1e4);
    assert!(ShearShellSection::from_laminate(&source,0.0).is_ok());
}
