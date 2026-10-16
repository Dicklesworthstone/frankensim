use super::*;
use super::super::{edge_cubic_mass, edge_cubic_transverse_gradient_shape};
use crate::{EdgeSupport, PlateMesh, PlateSection};

fn dot(a: &[f64],b: &[f64]) -> f64 { a.iter().zip(b).map(|(a,b)| a*b).sum() }
fn product(m: &[f64;81],a: &[f64;9],b: &[f64;9]) -> f64 {
    (0..9).map(|r| a[r]*dot(&m[9*r..9*r+9],b)).sum()
}
fn close(a: f64,b: f64,scale: f64) {
    assert!((a-b).abs() < 3e-12*scale.max(a.abs()).max(b.abs()), "{a} != {b}");
}
fn integral(x: usize,y: usize) -> f64 {
    let factorial = |n: usize| (1..=n).fold(1.0,|v,i| v*i as f64);
    factorial(x)*factorial(y)/factorial(x+y+2)
}
fn monomial(p: [usize;2],x: &[f64;3],y: &[f64;3]) -> [f64;9] {
    let at = |a: f64,b: f64,i: usize,j: usize| a.powi(i as i32)*b.powi(j as i32);
    std::array::from_fn(|i| {let n=i/3;match i%3 {
        0=>at(x[n],y[n],p[0],p[1]),
        1=>if p[0]==0 {0.0} else {p[0] as f64*at(x[n],y[n],p[0]-1,p[1])},
        _=>if p[1]==0 {0.0} else {p[1] as f64*at(x[n],y[n],p[0],p[1]-1)},
    }})
}

#[test]
fn complete_kinetic_energy_matches_all_quadratic_continuum_products() {
    let x=[0.0,1.0,0.0]; let y=[0.0,0.0,1.0];
    let powers=[[0,0],[1,0],[0,1],[2,0],[1,1],[0,2]];
    for h in [0.001,0.009,0.07] {
        let density=460.0;
        let rotary=density*h*h*h/12.0;
        let mr=rotary_mass(&x,&y,0.5*rotary).unwrap();
        let mt=edge_cubic_mass(&x,&y,0.5*density*h);
        for a in powers {for b in powers {
            let qa=monomial(a,&x,&y);let qb=monomial(b,&x,&y);
            let mut expected_rotation=0.0;
            if a[0]*b[0]!=0 {expected_rotation+=(a[0]*b[0]) as f64*integral(a[0]+b[0]-2,a[1]+b[1]);}
            if a[1]*b[1]!=0 {expected_rotation+=(a[1]*b[1]) as f64*integral(a[0]+b[0],a[1]+b[1]-2);}
            close(product(&mr,&qa,&qb),rotary*expected_rotation,rotary);
            close(product(&mt,&qa,&qb)+product(&mr,&qa,&qb),
                density*h*integral(a[0]+b[0],a[1]+b[1])+rotary*expected_rotation,density*h);
        }}
    }
}

#[test]
fn exact_rotary_matrix_matches_independent_positive_gauss_integration_of_the_motion_field() {
    for (x,y) in [([0.,1.,0.],[0.,0.,1.]),([0.2,1.4,-0.3],[-0.1,0.8,1.2])] {
        let mass=rotary_mass(&x,&y,0.004).unwrap();
        let root=(3.0_f64/5.0).sqrt();
        let points=[0.5*(1.0-root),0.5,0.5*(1.0+root)];
        let weights=[5.0/18.0,4.0/9.0,5.0/18.0];
        let mut oracle=[0.0;81];
        for i in 0..3 {for j in 0..3 {
            let u=points[i];let v=points[j];
            let shape=edge_cubic_transverse_gradient_shape(&x,&y,[1.0-u-(1.0-u)*v,u,(1.0-u)*v]);
            let weight=0.004*2.0*(1.0-u)*weights[i]*weights[j];
            for r in 0..9 {for c in 0..9 {
                oracle[9*r+c]+=weight*(shape[0][r]*shape[0][c]+shape[1][r]*shape[1][c]);
            }}
        }}
        for r in 0..9 {for c in 0..9 {
            assert_eq!(mass[9*r+c],mass[9*c+r]);
            close(mass[9*r+c],oracle[9*r+c],0.004);
        }}
        let rigid=[1.,0.,0.,1.,0.,0.,1.,0.,0.];
        close(product(&mass,&rigid,&rigid),0.0,0.004);
        // Independent scalar Cholesky of the total element kinetic form.
        let mt=edge_cubic_mass(&x,&y,1.0);
        let mut factor=[0.0;81];
        for i in 0..9 {for j in 0..=i {
            let value=mass[9*i+j]+mt[9*i+j]-(0..j).map(|k|factor[9*i+k]*factor[9*j+k]).sum::<f64>();
            factor[9*i+j]=if i==j {assert!(value>0.0);value.sqrt()} else {value/factor[9*j+j]};
        }}
    }
}

fn chart() -> PlateChart {
    let mesh=PlateMesh::from_unstructured(vec![(0.,0.),(1.0,0.2),(1.1,0.9),(-0.1,0.7)],vec![[0,1,2],[0,2,3]]).unwrap();
    let a=PlateSection::isotropic(12e9,0.25,0.009,460.0).unwrap();
    let b=PlateSection::isotropic(8e9,0.2,0.006,530.0).unwrap();
    PlateChart::with_boundary_and_regions(mesh,a,vec![],vec![]).unwrap()
        .with_element_sections(vec![a,b]).unwrap()
}
fn opts(support: EdgeSupport) -> AssemblyOptions { AssemblyOptions {support,pretension:200.0} }

#[test]
fn assembly_preserves_stiffness_beam_inertia_and_exact_dirichlet_restriction() {
    let chart=chart();let options=opts(EdgeSupport::Clamped);
    let beam=Stiffener {nodes:vec![0,2],e:12e9,g:0.8e9,area:0.0004,inertia:3e-8,
        torsion:2e-8,eccentricity:0.022,density:500.0};
    for beam_mass in [StiffenerMass::Lumped,StiffenerMass::ConsistentHermite,StiffenerMass::ConsistentEccentric] {
        let old=chart.assemble_with_mass(std::slice::from_ref(&beam),&options,TransverseMass::EdgeCubic,beam_mass).unwrap();
        let full=chart.assemble_edge_cubic_kinetic_mass(std::slice::from_ref(&beam),&options,beam_mass).unwrap();
        let old_bare=chart.assemble_with_mass(&[],&options,TransverseMass::EdgeCubic,beam_mass).unwrap();
        let bare=chart.assemble_edge_cubic_kinetic_mass(&[],&options,beam_mass).unwrap();
        assert_eq!(old.k,full.k);assert_eq!(old.dof_map,full.dof_map);
        assert_ne!(old.m,full.m);
        assert_eq!(full.m,chart.assemble_edge_cubic_kinetic_mass(std::slice::from_ref(&beam),&options,beam_mass).unwrap().m);
        for r in 0..full.free {for c in 0..full.free {
            close(full.m.get(r,c)-bare.m.get(r,c),old.m.get(r,c)-old_bare.m.get(r,c),0.2);
        }}
        for support in [EdgeSupport::Clamped,EdgeSupport::SimplySupported] {
            let mut fixed=chart.clone();fixed.boundary_nodes=vec![0,3];
            let reduced=fixed.assemble_edge_cubic_kinetic_mass(std::slice::from_ref(&beam),&opts(support),beam_mass).unwrap();
            for (r,rr) in reduced.dof_map.iter().enumerate() {for (c,cc) in reduced.dof_map.iter().enumerate() {
                if let (Some(rr),Some(cc))=(rr,cc) {
                    assert_eq!(reduced.m.get(*rr,*cc),full.m.get(r,c));
                    assert_eq!(reduced.k.get(*rr,*cc),full.k.get(r,c));
                }
            }}
        }
    }
}

#[test]
fn rigid_transverse_mass_and_tilt_include_each_supplied_piecewise_section_once() {
    let chart=chart();let model=chart.assemble_edge_cubic_kinetic_mass(&[],&opts(EdgeSupport::Clamped),StiffenerMass::Lumped).unwrap();
    let q:Vec<_>=chart.mesh.nodes.iter().flat_map(|&(x,y)|[0.7+0.2*x-0.3*y,0.2,-0.3]).collect();
    let mut mq=vec![0.;model.free];model.m.spmv(&q,&mut mq);
    let mut expected=0.0;
    for (i,tri) in chart.mesh.tris.iter().enumerate() {
        let s=chart.section_field().at(i);
        let p=tri.map(|n|chart.mesh.nodes[n]);
        let area=0.5*((p[1].0-p[0].0)*(p[2].1-p[0].1)-(p[2].0-p[0].0)*(p[1].1-p[0].1));
        let w=p.map(|(x,y)|0.7+0.2*x-0.3*y);
        let mean_square=(w.iter().map(|x|x*x).sum::<f64>()+w[0]*w[1]+w[1]*w[2]+w[2]*w[0])/6.0;
        expected+=s.density*s.thickness*area*(mean_square+s.thickness*s.thickness/12.0*(0.2*0.2+0.3*0.3));
    }
    close(dot(&q,&mq),expected,expected);
}

#[test]
fn malformed_or_unrepresentable_kinetic_input_refuses_instead_of_adding_a_mass_floor() {
    for scale in [0.0,-1.0,f64::NAN,f64::INFINITY] {
        assert!(rotary_mass(&[0.,1.,0.],&[0.,0.,1.],scale).is_err());
    }
    assert!(rotary_mass(&[0.,0.,0.],&[0.,0.,0.],1.0).is_err());
    let mut chart=chart();chart.element_sections.as_mut().unwrap()[0].thickness=f64::MIN_POSITIVE;
    assert!(chart.assemble_edge_cubic_kinetic_mass(&[],&opts(EdgeSupport::Clamped),StiffenerMass::Lumped).is_err());
}
