use super::*;
use fs_solid::linear3::visco::SolidProjectionBudget;
use fs_solid::reduce::visco::{
    HereditaryModalValve, HereditaryValveControls, reduce_hereditary_pressure_loaded_valve,
};

fn hereditary_controls() -> HereditaryValveControls {
    HereditaryValveControls {
        assembly: TetAssemblyBudget::standard(),
        internal_variables: ViscoStateBudget {
            max_internal_variables: 1000,
        },
        projection: SolidProjectionBudget {
            max_coordinates: 8,
            max_basis_scalars: 1000,
            max_material_points: 100,
            max_scalar_terms: 1_000_000,
        },
        maximum_factor_scalar_terms: 1_000_000,
        maximum_basis_error: 1e-10,
    }
}
fn coupled_material() -> TetViscoMaterial {
    let stiffness = diagonal([80., 100.1, 100., 40., 45., 50.]);
    let mut viscosity = stiffness.map(|row| row.map(|x| x * 0.001));
    viscosity[1][2] = 0.09;
    viscosity[2][1] = 0.09;
    let mut factor = [[0.0; 6]; 6];
    factor[0] = [0., 7., -9., 0., 0., 0.];
    factor[1] = [0., 0.2, 0.1, 0., 0., 0.];
    TetViscoMaterial {
        law: TensorProny::new(
            stiffness,
            Some(viscosity),
            vec![
                PronyBranch::from_gram_rows(factor, 0.02).unwrap(),
                PronyBranch::new(stiffness.map(|r| r.map(|x| x * 0.3)), 0.08).unwrap(),
            ],
            (0., 10_000.),
        )
        .unwrap(),
        density_kg_m3: 1.,
        material_state_identity: [93; 32],
    }
}
fn projected_source(image: &HereditaryModalValve, omega: f64, cx: &Cx<'_>) -> (Vec<f64>, Vec<f64>) {
    let source = image.source().complex_stiffness(omega, cx).unwrap();
    let n = image.modes().len();
    let mass = image.total_mass_kg();
    let project = |a: &fs_sparse::Csr| {
        let mut out = vec![0.; n * n];
        let mut product = vec![0.; a.nrows()];
        for j in 0..n {
            a.spmv(&image.modes()[j].phi, &mut product);
            for i in 0..n {
                out[i * n + j] = mass
                    * image.modes()[i]
                        .phi
                        .iter()
                        .zip(&product)
                        .map(|(x, y)| x * y)
                        .sum::<f64>();
            }
        }
        out
    };
    (project(&source.real), project(&source.imag))
}

#[test]
fn full_hereditary_image_matches_original_sparse_material_over_frequency() {
    with_cx(|cx| {
        let law = coupled_material();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 7 && d != 11).collect();
        let mut request = request(TetViscoField::Uniform(&law));
        request.fixed_dofs = &fixed;
        let load = [[0.; 3], [0.; 3], [0., 0.002, 0.], [0., 0., 0.001]];
        let image =
            reduce_hereditary_pressure_loaded_valve(&request, &load, &hereditary_controls(), cx)
                .unwrap();
        assert_eq!(image.modes().len(), 2);
        assert!(
            image
                .opening_projection_per_sqrt_kg()
                .iter()
                .any(|x| *x == 0.)
        );
        assert!(image.viscosity_n_s_m()[1].abs() > 1e-5);
        assert!(
            image
                .branches()
                .iter()
                .any(|b| b.rows().iter().any(|r| (r[0] * r[1]).abs() > 1e-5))
        );
        assert!(image.memory_channel_count() <= 4);
        assert!(image.basis_error() < 1e-10);
        for omega in [0., 1., 10., 30., 100., 1_000.] {
            let (real, imag) = projected_source(&image, omega, cx);
            let n = image.modes().len();
            for i in 0..n {
                for j in 0..n {
                    let mut storage = if i == j {
                        image.total_mass_kg() * image.modes()[i].lambda
                    } else {
                        0.
                    };
                    let mut loss = omega * image.viscosity_n_s_m()[i * n + j];
                    for branch in image.branches() {
                        let z = omega * branch.relaxation_time_s();
                        let gram = branch.rows().iter().map(|r| r[i] * r[j]).sum::<f64>();
                        storage += gram * z * z / (1. + z * z);
                        loss += gram * z / (1. + z * z);
                    }
                    let scale = (real[i * n + i] * real[j * n + j]).sqrt();
                    assert!((storage - real[i * n + j]).abs() < 1e-11 * scale);
                    assert!((loss - imag[i * n + j]).abs() < 1e-11 * scale);
                }
            }
        }
        let c = image.compliance();
        close(c.retained_volume_m3_per_pa, c.exact_volume_m3_per_pa, 1e-10);
        close(c.retained_opening_m_per_pa, c.exact_opening_m_per_pa, 1e-10);
    });
}

#[test]
fn strong_source_dispersion_uses_equilibrium_modes_without_a_frozen_loss_fit() {
    with_cx(|cx| {
        let law = material(0.001, 0.8, 0.02);
        let request = request(TetViscoField::Uniform(&law));
        assert!(reduce_pressure_loaded_valve(&request, &PRESSURE, &controls(), cx).is_err());
        let image = reduce_hereditary_pressure_loaded_valve(
            &request,
            &PRESSURE,
            &hereditary_controls(),
            cx,
        )
        .unwrap();
        close(image.modes()[0].lambda, 1_000., 1e-12);
        close(
            image.compliance().retained_opening_force_m_per_n,
            image.compliance().exact_opening_force_m_per_n,
            1e-12,
        );
        assert_eq!(image.branches()[0].relaxation_time_s(), 0.02);
        assert!(image.instantaneous_frequency_bound_rad_s() > 40.);
    });
}

#[test]
fn original_material_domain_complete_basis_and_factor_budgets_are_not_relaxed() {
    with_cx(|cx| {
        let law = coupled_material();
        let fixed: Vec<_> = (0..12).filter(|&d| d != 7 && d != 11).collect();
        let mut request = request(TetViscoField::Uniform(&law));
        request.fixed_dofs = &fixed;
        let mut control = hereditary_controls();
        control.projection.max_coordinates = 1;
        assert!(
            reduce_hereditary_pressure_loaded_valve(&request, &PRESSURE, &control, cx).is_err()
        );
        control = hereditary_controls();
        control.maximum_factor_scalar_terms = 1;
        assert!(
            reduce_hereditary_pressure_loaded_valve(&request, &PRESSURE, &control, cx).is_err()
        );
        let mut bounded = material(0.001, 0.8, 0.02);
        bounded.law = bounded.law.with_band((0., 40.)).unwrap();
        assert!(
            reduce_hereditary_pressure_loaded_valve(
                &super::request(TetViscoField::Uniform(&bounded)),
                &PRESSURE,
                &hereditary_controls(),
                cx
            )
            .is_err()
        );
        let stiffness = diagonal([80., 100.1, 100., 40., 45., 50.]);
        let limit = TetViscoMaterial {
            law: TensorProny::new_mandel(stiffness, vec![], 0.01).unwrap(),
            density_kg_m3: 1.,
            material_state_identity: [94; 32],
        };
        let image = reduce_hereditary_pressure_loaded_valve(
            &super::request(TetViscoField::Uniform(&limit)),
            &PRESSURE,
            &hereditary_controls(),
            cx,
        )
        .unwrap();
        assert!(image.validate_displacement_m(&[1e-4]).is_ok());
        assert!(image.validate_displacement_m(&[0.01]).is_err());
    });
}
