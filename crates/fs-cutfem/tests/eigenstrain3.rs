//! G0/G1/G3 thermal compatibility, physical cut projection and original tensors.
use std::cell::Cell;
use std::ops::ControlFlow;

use fs_blake3::ContentHash;
use fs_cutfem::elastic3::dirichlet::EmbeddedDirichletOptions3;
use fs_cutfem::elastic3::{
    CutElasticity3, ElasticMaterialField3, ElasticityError3, ElasticityFormulation3,
    ElasticityOptions3,
};
use fs_cutfem::quad3::surface::SurfaceOptions3;
use fs_cutfem::quad3::{QuadratureControl3, QuadratureOptions3};
use fs_cutfem::{CutSdf3, HeightAxis, HexCell};
use fs_ivl::Interval;
use fs_material::IsotropicElastic;
use fs_material::elastic::acoustic::ElasticAcousticState;
use fs_material::tensor::{contract, from_mandel, to_mandel};
use fs_solver::op::LinearOp;

struct Slab(f64);
impl CutSdf3 for Slab {
    fn value(&self, p: [f64; 3]) -> f64 {
        p[2] - self.0
    }
    fn enclose(&self, lo: [f64; 3], hi: [f64; 3]) -> Interval {
        Interval::new(lo[2], hi[2]) - Interval::new(self.0, self.0)
    }
    fn derivative_enclose(&self, _: [f64; 3], _: [f64; 3], axis: HeightAxis) -> Interval {
        let value = if axis == HeightAxis::Z { 1.0 } else { 0.0 };
        Interval::new(value, value)
    }
}
fn poll() -> ControlFlow<()> {
    ControlFlow::Continue(())
}
fn domain() -> HexCell {
    HexCell::try_new([0.0; 3], [1.0; 3]).unwrap()
}
fn control<T>(f: impl FnOnce(&mut QuadratureControl3<'_>) -> T) -> T {
    let mut checkpoint = |_| poll();
    let mut quadrature = QuadratureControl3::new(
        QuadratureOptions3 {
            depth: 0,
            ..Default::default()
        },
        &mut checkpoint,
    )
    .unwrap();
    f(&mut quadrature)
}
fn options(formulation: ElasticityFormulation3) -> ElasticityOptions3 {
    ElasticityOptions3 {
        formulation,
        ..Default::default()
    }
}
fn operator(
    formulation: ElasticityFormulation3,
    all_fixed: bool,
    law: &IsotropicElastic,
) -> CutElasticity3 {
    control(|control| {
        CutElasticity3::build(
            domain(),
            [1; 3],
            &Slab(0.73),
            law,
            &|p| all_fixed || p == [0.0; 3],
            options(formulation),
            control,
        )
        .unwrap()
    })
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(actual.is_finite() && expected.is_finite());
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}
fn action(op: &impl LinearOp, u: &[f64]) -> Vec<f64> {
    let mut applied = vec![0.0; op.n()];
    op.apply(u, &mut applied);
    applied
}

#[test]
fn compatible_nonuniform_free_strain_drives_the_original_operator_without_stress() {
    let amplitude = 0.004;
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let poisson = if formulation == ElasticityFormulation3::BBar {
            0.4999
        } else {
            0.25
        };
        let law = IsotropicElastic::new(1200.0, poisson, 0.1).unwrap();
        let mut op = operator(formulation, false, &law);
        let u: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| [amplitude * p[0] * p[2], 0.0, 0.0])
            .collect();
        let prepared = op
            .prepare_eigenstrain(
                &|_, p| Ok([amplitude * p[2], 0.0, 0.0, 0.0, 0.0, 0.5 * amplitude * p[0]]),
                100_000,
                poll,
            )
            .unwrap();
        // Repeat on the same geometry after a material-scale change.
        for scale in [1.0, 0.37] {
            op.set_scales(&[scale]).unwrap();
            let rhs = op.eigenstrain_load(&prepared, poll).unwrap();
            for (ku, b) in action(&op, &u).iter().zip(rhs) {
                close(*ku, b, 5e-10);
            }
            let response = op
                .eigenstrain_response(&prepared, &u, 0.1, 100_000, poll)
                .unwrap();
            close(response.iter().map(|p| p.weight).sum(), 0.73, 1e-12);
            for point in response {
                for stress in point.stress {
                    close(stress, 0.0, 5e-10);
                }
                close(point.energy_density, 0.0, 1e-18);
            }
        }
    }
}

#[test]
fn restrained_thermal_expansion_is_compressive_and_stores_elastic_energy() {
    let alpha_delta_t: f64 = 0.002;
    let law = IsotropicElastic::new(1200.0, 0.25, 0.1).unwrap();
    let bulk = law.youngs / (3.0 * (1.0 - 2.0 * law.poisson));
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = operator(formulation, true, &law);
        let free = [alpha_delta_t, alpha_delta_t, alpha_delta_t, 0.0, 0.0, 0.0];
        let prepared = op
            .prepare_eigenstrain(&|_, _| Ok(free), 100_000, poll)
            .unwrap();
        let force = op.eigenstrain_load(&prepared, poll).unwrap();
        assert!(force.iter().all(|value| *value == 0.0));
        let solution = op
            .solve_controlled(&force, 1e-12, 100, 1, |_| poll())
            .unwrap();
        let response = op
            .eigenstrain_response(&prepared, solution.coefficients(), 0.1, 100_000, poll)
            .unwrap();
        for point in response {
            for value in &point.stress[..3] {
                close(*value, -3.0 * bulk * alpha_delta_t, 1e-12);
            }
            for value in &point.stress[3..] {
                assert_eq!(*value, 0.0);
            }
            close(
                point.energy_density,
                4.5 * bulk * alpha_delta_t.powi(2),
                1e-12,
            );
        }
    }
}

#[test]
fn thermal_bbar_uses_mean_free_dilatation_over_the_physical_cut() {
    let law = IsotropicElastic::new(1200.0, 0.4999, 0.1).unwrap();
    let bulk = law.youngs / (3.0 * (1.0 - 2.0 * law.poisson));
    let op = operator(ElasticityFormulation3::BBar, true, &law);
    let prepared = op
        .prepare_eigenstrain(
            &|_, p| {
                let expansion = 0.001 + 0.002 * p[2];
                Ok([expansion, expansion, expansion, 0.0, 0.0, 0.0])
            },
            100_000,
            poll,
        )
        .unwrap();
    let response = op
        .eigenstrain_response(&prepared, &vec![0.0; op.n()], 0.1, 100_000, poll)
        .unwrap();
    let mean_expansion = 0.001 + 0.002 * 0.73 / 2.0;
    for point in response {
        for value in &point.stress[..3] {
            close(*value, -3.0 * bulk * mean_expansion, 1e-11);
        }
        close(
            point.energy_density,
            4.5 * bulk * mean_expansion.powi(2),
            1e-11,
        );
        close(
            point.free_strain[0],
            0.001 + 0.002 * point.position[2],
            1e-14,
        );
    }
}

#[test]
fn embedded_nitsche_load_has_the_eigenstress_traction_sign() {
    let amplitude = 0.002;
    let law = IsotropicElastic::new(1200.0, 0.25, 0.1).unwrap();
    for formulation in [
        ElasticityFormulation3::Displacement,
        ElasticityFormulation3::BBar,
    ] {
        let op = control(|control| {
            CutElasticity3::build_with_embedded_dirichlet(
                domain(),
                [2; 3],
                &Slab(0.73),
                &law,
                &|_| false,
                &|_, _| true,
                options(formulation),
                EmbeddedDirichletOptions3::default(),
                SurfaceOptions3::default(),
                control,
            )
            .unwrap()
        });
        let free = [amplitude, amplitude, amplitude, 0.0, 0.0, 0.0];
        let prepared = op
            .prepare_eigenstrain(&|_, _| Ok(free), 100_000, poll)
            .unwrap();
        let mut rhs = op.eigenstrain_load(&prepared, poll).unwrap();
        let lifting = op
            .prescribed_displacement_load(&|p, _| p.map(|x| amplitude * x), poll)
            .unwrap();
        for (rhs, lift) in rhs.iter_mut().zip(lifting) {
            *rhs += lift;
        }
        let exact: Vec<_> = op
            .nodes()
            .iter()
            .flat_map(|p| p.map(|x| amplitude * x))
            .collect();
        for (ku, rhs) in action(&op, &exact).iter().zip(&rhs) {
            close(*ku, *rhs, 3e-11);
        }
        let solution = op
            .solve_controlled(&rhs, 1e-12, 10_000, 8, |_| poll())
            .unwrap();
        for (actual, exact) in solution.coefficients().iter().zip(exact) {
            close(*actual, exact, 1e-10);
        }
        for point in op
            .eigenstrain_response(&prepared, solution.coefficients(), 0.1, 100_000, poll)
            .unwrap()
        {
            for stress in point.stress {
                close(stress, 0.0, 2e-7);
            }
        }
    }
}

#[test]
fn full_tensor_and_original_material_inventory_drive_restrained_response() {
    let diagonal = [1200.0, 900.0, 1500.0, 800.0, 650.0, 700.0];
    let mut stiffness = [[0.0; 6]; 6];
    for i in 0..6 {
        stiffness[i][i] = diagonal[i];
    }
    stiffness[0][3] = 80.0;
    stiffness[3][0] = 80.0;
    stiffness[1][2] = 150.0;
    stiffness[2][1] = 150.0;
    let identity = ContentHash([77; 32]);
    let material = ElasticAcousticState::try_new_mandel(stiffness, 2700.0, identity, 0.1).unwrap();
    let mut op = control(|control| {
        CutElasticity3::build_with_material_field(
            domain(),
            [1; 3],
            &Slab(0.73),
            ElasticMaterialField3::Uniform(&material),
            &|_| true,
            options(ElasticityFormulation3::Displacement),
            control,
        )
        .unwrap()
    });
    let free = [0.001, -0.0004, 0.0007, 0.0003, -0.0002, 0.0005];
    let prepared = op
        .prepare_eigenstrain(&|_, _| Ok(free), 100_000, poll)
        .unwrap();
    let mandel = to_mandel(&free);
    let expected_mandel =
        std::array::from_fn(|i| -(0..6).map(|j| stiffness[i][j] * mandel[j]).sum::<f64>());
    let expected = from_mandel(&expected_mandel);
    let mass = op.reference_cell_masses_kg().unwrap()[0];
    op.set_scales(&[0.31]).unwrap();
    assert_eq!(
        op.cell_material_state(0).unwrap().material_state_identity(),
        identity
    );
    assert_eq!(op.reference_cell_masses_kg().unwrap()[0], mass);
    for point in op
        .eigenstrain_response(&prepared, &vec![0.0; op.n()], 0.1, 100_000, poll)
        .unwrap()
    {
        for i in 0..6 {
            close(point.stress[i], 0.31 * expected[i], 1e-12);
        }
        close(
            point.energy_density,
            -0.5 * 0.31 * contract(&free, &expected),
            1e-12,
        );
    }
}

#[test]
fn expansion_and_elastic_domains_are_separate_and_preparation_is_transactional() {
    let law = IsotropicElastic::new(1200.0, 0.25, 1e-5).unwrap();
    let op = operator(ElasticityFormulation3::Displacement, false, &law);
    let amplitude = 0.001;
    let free = [amplitude, amplitude, amplitude, 0.0, 0.0, 0.0];
    let prepared = op
        .prepare_eigenstrain(&|_, _| Ok(free), 100_000, poll)
        .unwrap();
    let u: Vec<_> = op
        .nodes()
        .iter()
        .flat_map(|p| p.map(|x| amplitude * x))
        .collect();
    assert!(
        op.eigenstrain_response(&prepared, &u, 0.01, 100_000, poll)
            .is_ok()
    );
    assert!(
        op.eigenstrain_response(&prepared, &u, 1.0, 100_000, poll)
            .is_err()
    );
    assert!(
        op.eigenstrain_response(&prepared, &u, 0.0001, 100_000, poll)
            .is_err()
    );
    assert!(
        op.eigenstrain_response(&prepared, &vec![0.0; op.n()], 0.01, 100_000, poll)
            .is_err()
    );
    assert!(
        op.eigenstrain_response(&prepared, &u, 0.01, prepared.bulk_points() - 1, poll)
            .is_err()
    );
    let another = operator(ElasticityFormulation3::Displacement, false, &law);
    assert!(another.eigenstrain_load(&prepared, poll).is_err());
    let calls = Cell::new(0);
    assert!(
        op.prepare_eigenstrain(
            &|_, _| {
                calls.set(calls.get() + 1);
                Ok(free)
            },
            0,
            poll
        )
        .is_err()
    );
    assert_eq!(calls.get(), 0);
    assert!(matches!(
        op.prepare_eigenstrain(&|_, _| Ok(free), 100_000, || ControlFlow::Break(())),
        Err(ElasticityError3::Cancelled)
    ));
    assert!(
        op.prepare_eigenstrain(&|_, _| Ok([f64::NAN; 6]), 100_000, poll)
            .is_err()
    );
    assert_eq!(op.scales(), &[1.0]);
    assert!(op.eigenstrain_load(&prepared, poll).is_ok());
}
