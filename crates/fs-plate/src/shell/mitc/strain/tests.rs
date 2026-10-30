use super::*;

fn section(kappa: f64) -> ShearShellSection {
    ShearShellSection::from_orthotropic_solid(
        10e9, 0.8e9, 0.25, 0.6e9, 0.5e9, 0.3e9,
        0.003, 450.0, 0.2, kappa,
    ).unwrap()
}

fn flat() -> (ShellMesh, NodalDirectors) {
    let mesh = ShellMesh::new(
        vec![[0.0; 3], [0.8, 0.0, 0.0], [0.1, 0.6, 0.0]], vec![[0, 1, 2]],
    ).unwrap();
    let directors = NodalDirectors::new(&mesh, vec![[0.0, 0.0, 1.0]; 3], 0.9).unwrap();
    (mesh, directors)
}

fn close(actual: f64, expected: f64) {
    assert!((actual - expected).abs() <= 2e-11 * expected.abs().max(1e-15),
        "actual={actual:.17e}, expected={expected:.17e}");
}

#[test]
fn harmonic_strain_affine_extension_and_engineering_shear_have_physical_norms() {
    let (mesh, directors) = flat();
    let prepared = PreparedElementStrain::new(&mesh, 0, &directors, &section(5.0/6.0),
        &mut || false).unwrap();
    let mut extension = [0.0; 18];
    let mut shear = [0.0; 18];
    let amplitude = 1e-4;
    for (node, p) in mesh.nodes.iter().enumerate() {
        extension[6*node] = amplitude * p[0];
        shear[6*node] = amplitude * p[1];
    }
    let extension = prepared.harmonic_layer_bound(&extension, &[0.0; 18],
        [-0.002, 0.003], 0.25, &mut || false).unwrap();
    close(extension.component_amplitudes[0], amplitude);
    close(extension.component_amplitudes[2], amplitude/3.0);
    close(extension.frobenius_upper_bound, amplitude * (10.0_f64/9.0).sqrt());
    let shear = prepared.harmonic_layer_bound(&shear, &[0.0; 18],
        [-0.002, 0.003], 0.25, &mut || false).unwrap();
    close(shear.component_amplitudes[3], amplitude/2.0);
    close(shear.frobenius_upper_bound, amplitude/core::f64::consts::SQRT_2);
    assert!(shear.component_amplitudes[0..3].iter().all(|x| *x < 1e-16));
}

#[test]
fn harmonic_strain_cycle_bound_preserves_correlated_component_phase() {
    let (mesh, directors) = flat();
    let prepared = PreparedElementStrain::new(&mesh, 0, &directors, &section(1.0),
        &mut || false).unwrap();
    let mut real = [0.0; 18];
    let mut imaginary = [0.0; 18];
    for (node, p) in mesh.nodes.iter().enumerate() {
        real[6*node] = 3e-4 * p[0];
        imaginary[6*node+1] = 4e-4 * p[1];
    }
    let bound = prepared.harmonic_layer_bound(&real, &imaginary, [-0.001, 0.002],
        0.0, &mut || false).unwrap();
    // 3*cos(t) in xx and 4*sin(t) in yy have peak norm 4, not 5.
    close(bound.frobenius_upper_bound, 4e-4);
    close(bound.component_amplitudes[0], 3e-4);
    close(bound.component_amplitudes[1], 4e-4);
    let phase = 0.73_f64;
    let rotated_real = core::array::from_fn(|i| real[i]*phase.cos()-imaginary[i]*phase.sin());
    let rotated_imaginary = core::array::from_fn(|i| real[i]*phase.sin()+imaginary[i]*phase.cos());
    let rotated = prepared.harmonic_layer_bound(&rotated_real, &rotated_imaginary,
        [-0.001, 0.002], 0.0, &mut || false).unwrap();
    close(rotated.frobenius_upper_bound, bound.frobenius_upper_bound);
    for i in 0..6 { close(rotated.component_amplitudes[i], bound.component_amplitudes[i]); }
}

#[test]
fn harmonic_strain_bernstein_field_encloses_original_condensed_curved_element() {
    let mesh = ShellMesh::new(
        vec![[0.0; 3], [0.8, 0.1, 0.04], [0.05, 0.6, -0.02]], vec![[0, 1, 2]],
    ).unwrap();
    let g = mesh.facet(0).unwrap();
    let vectors = (0..3).map(|i| {
        let v = [g.frame[2][0] + 0.02*(i as f64 - 1.0),
            g.frame[2][1] + 0.01*(i as f64 - 0.5), g.frame[2][2]];
        let norm = dot(v, v).sqrt();
        v.map(|x| x/norm)
    }).collect();
    let directors = NodalDirectors::new(&mesh, vectors, 0.9).unwrap();
    // Offset the reference surface so the original equilibrium condensation
    // includes nonzero laminate B, as it does for one-sided coating stacks.
    let base = section(5.0/6.0);
    let offset = 0.0004;
    let (a, d) = base.principal_matrices();
    let bending = fs_material::sheet::SheetMatrix::from_coefficients(
        core::array::from_fn(|i| d.coefficients()[i] + offset*offset*a.coefficients()[i]),
    ).unwrap();
    let section = ShearShellSection::new_coupled(*a, a.coefficients().map(|v| offset*v),
        bending, [0.5e9*0.003,0.0,0.3e9*0.003], 5.0/6.0,
        base.areal_mass_kg_m2(), offset*base.areal_mass_kg_m2(),
        base.rotary_inertia_kg()+offset*offset*base.areal_mass_kg_m2(), 0.2,
    ).unwrap();
    let prepared = PreparedElementStrain::new(&mesh, 0, &directors, &section,
        &mut || false).unwrap();
    let real = core::array::from_fn(|i| ((i*7 % 11) as f64 - 5.0) * 2e-5);
    let imaginary = core::array::from_fn(|i| ((i*3 % 13) as f64 - 6.0) * 1e-5);
    let field = prepared.harmonic_field(&real, &imaginary, &mut || false).unwrap();
    let faces = [-0.002, 0.0035];
    let nu = 0.31;
    let bound = field.isotropic_layer_bound(faces, nu, &mut || false).unwrap();
    let normals = mesh.tris[0].map(|i| directors.vectors()[i]);
    let original_real = element::ElementStrainField::new(&mesh, &mesh.tris[0], g,
        normals, &section, &real).unwrap();
    let original_imaginary = element::ElementStrainField::new(&mesh, &mesh.tris[0], g,
        normals, &section, &imaginary).unwrap();
    let mut sampled_max = 0.0_f64;
    for i in 0..=12 { for j in 0..=12-i {
        let r = i as f64/12.0;
        let s = j as f64/12.0;
        let l = 1.0-r-s;
        let bernstein = [l*l, r*r, s*s, 2.0*l*r, 2.0*r*s, 2.0*s*l];
        let (mr, kr, sr) = original_real.eval(r, s);
        let (mi, ki, si) = original_imaginary.eval(r, s);
        // Establish exact-polynomial reproduction before checking its bound.
        for (controls, m, k, shear) in [(&field.real, mr, kr, sr),
            (&field.imaginary, mi, ki, si)] {
            let original = [m[0],m[1],m[2],k[0],k[1],k[2],shear[0],shear[1]];
            for component in 0..8 {
                let reconstructed = (0..6).map(|c| bernstein[c]*controls[c][component]).sum::<f64>();
                assert!((reconstructed-original[component]).abs() < 2e-15);
            }
        }
        for through in 0..=4 {
            let z = faces[0] + (faces[1]-faces[0])*(through as f64/4.0);
            for phase in 0..48 {
                let theta = core::f64::consts::TAU * phase as f64/48.0;
                let e = core::array::from_fn::<_,3,_>(|c|
                    (mr[c]+z*kr[c])*theta.cos() + (mi[c]+z*ki[c])*theta.sin());
                let transverse = core::array::from_fn::<_,2,_>(|c|
                    sr[c]*theta.cos() + si[c]*theta.sin());
                let zz = -nu/(1.0-nu)*(e[0]+e[1]);
                let norm = (e[0]*e[0]+e[1]*e[1]+zz*zz
                    +0.5*(e[2]*e[2]+transverse[0]*transverse[0]+transverse[1]*transverse[1])).sqrt();
                sampled_max = sampled_max.max(norm);
                assert!(norm <= bound.frobenius_upper_bound*(1.0+2e-12));
            }
        }
    } }
    assert!(sampled_max > 0.0);
    assert!(bound.ruling_control_point < 6 && bound.ruling_z_face < 2);
    // Every film face enters separately, and the same prepared field is reusable.
    let lower = field.isotropic_layer_bound([faces[0];2], nu, &mut || false).unwrap();
    let upper = field.isotropic_layer_bound([faces[1];2], nu, &mut || false).unwrap();
    close(bound.frobenius_upper_bound, lower.frobenius_upper_bound.max(upper.frobenius_upper_bound));
}

#[test]
fn harmonic_strain_layer_faces_and_tensor_shear_order_are_explicit() {
    // An affine bending/shear polynomial: engineering shear is in entries 2/6/7.
    // This directly checks thickness closure and Voigt order independently of
    // the element's equilibrium choice of its two internal rotations.
    let mut real = [[0.0; 8]; 6];
    for control in &mut real {
        control[3] = 0.02;
        control[6] = 0.006;
        control[7] = 0.008;
    }
    let field = HarmonicElementStrain { real, imaginary: [[0.0;8];6] };
    let bound = field.isotropic_layer_bound([-0.3, 0.1], 0.0, &mut || false).unwrap();
    assert_eq!(bound.ruling_z_face, 0);
    close(bound.component_amplitudes[0], 0.006);
    close(bound.component_amplitudes[4], 0.004);
    close(bound.component_amplitudes[5], 0.003);
    close(bound.frobenius_upper_bound, (0.006_f64.powi(2)+2.0*(0.004_f64.powi(2)+0.003_f64.powi(2))).sqrt());
    // An edge coefficient is not its midpoint sample: p=2*lambda0*lambda1.
    let mut edge = [[0.0;8];6];
    edge[3][3] = 1.0;
    let field = HarmonicElementStrain { real:edge, imaginary:[[0.0;8];6] };
    let bound = field.isotropic_layer_bound([0.0, 1.0], 0.0, &mut || false).unwrap();
    assert_eq!(bound.ruling_control_point, 3);
    assert_eq!(bound.ruling_z_face, 1);
    close(bound.frobenius_upper_bound, 1.0);
    assert!(bound.frobenius_upper_bound > 0.5); // Actual edge midpoint value.
}

#[test]
fn harmonic_strain_refusals_leave_preparation_reusable_and_scaled_norms_finite() {
    let (mesh, directors) = flat();
    let section = section(1.0);
    assert!(matches!(PreparedElementStrain::new(&mesh, 0, &directors, &section,
        &mut || true), Err(PlateError::Cancelled)));
    assert!(PreparedElementStrain::new(&mesh, 1, &directors, &section, &mut || false).is_err());
    let prepared = PreparedElementStrain::new(&mesh, 0, &directors, &section,
        &mut || false).unwrap();
    let mut nonfinite = [0.0;18];
    nonfinite[7] = f64::NAN;
    assert!(prepared.harmonic_field(&nonfinite, &[0.0;18], &mut || false).is_err());
    assert!(matches!(prepared.harmonic_field(&[0.0;18], &[0.0;18], &mut || true),
        Err(PlateError::Cancelled)));
    let field = prepared.harmonic_field(&[0.0;18], &[0.0;18], &mut || false).unwrap();
    for (faces, nu) in [([1.0,0.0],0.3), ([0.0,f64::INFINITY],0.3),
        ([0.0,1.0],0.5), ([0.0,1.0],-1.0), ([0.0,1.0],f64::NAN)] {
        assert!(field.isotropic_layer_bound(faces, nu, &mut || false).is_err());
    }
    assert!(matches!(field.isotropic_layer_bound([0.0,1.0], 0.3, &mut || true),
        Err(PlateError::Cancelled)));
    assert_eq!(field.isotropic_layer_bound([0.0,1.0], 0.3, &mut || false)
        .unwrap().frobenius_upper_bound, 0.0);
    for scale in [1e-280, 1e280] {
        let norm = tensor_cycle_norm(&[scale,0.0,0.0,0.0,0.0,0.0],
            &[0.0,scale,0.0,0.0,0.0,0.0]).unwrap();
        assert!(norm > 0.0 && (norm/scale - 1.0).abs() < 2e-14);
    }
    assert!(tensor_cycle_norm(&[f64::MAX;6], &[0.0;6]).is_err());
}
