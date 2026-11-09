use super::*;
use crate::adjoint::RobinResponse;
use crate::assemble::{
    assemble_enthalpy_jacobian_with_optional_interfaces, assemble_jacobian, assemble_operator,
    assemble_operator_scaled_with_interfaces, full_residual,
};
use crate::{
    ConductionMesh, ConductionProblem, ElementMaterials, MaterialId, MaterialTable, ScalarField,
    ThermalBc, ThermalBoundaryBuilder, element_heat_flux, element_heat_flux_assigned,
};
use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, ExecMode, StreamKey};
use fs_rep_mesh::TetComplex;

fn with_cx<R>(f: impl FnOnce(&CancelGate, &Cx<'_>) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &gate,
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 41,
                    kernel_id: 519,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
        )
    })
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1.0),
        "{actual:e} != {expected:e}"
    );
}

fn curve() -> ConductivityModel {
    ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 10.0), (350.0, 10.0), (450.0, 210.0)])
            .unwrap(),
    )
}

fn anisotropic() -> ConductivityModel {
    let s = std::f64::consts::FRAC_1_SQRT_2;
    ConductivityModel::orthotropic(
        [[s, s, 0.0], [-s, s, 0.0], [0.0, 0.0, 1.0]],
        [
            ConductivityTable::declared_curve(vec![(250.0, 10.0), (350.0, 10.0), (450.0, 210.0)])
                .unwrap(),
            ConductivityTable::declared_curve(vec![(250.0, 30.0), (350.0, 30.0), (450.0, 10.0)])
                .unwrap(),
            ConductivityTable::declared(7.0).unwrap(),
        ],
    )
    .unwrap()
}

fn mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [1.0, 0.0, 0.0],
            [0.0, 1.0, 0.0],
            [0.0, 0.0, 1.0],
        ],
    )
    .unwrap()
}

#[test]
fn beta_moments_distinguish_integrated_material_from_mean_sampling() {
    // Uniform tetrahedral barycentric coordinates give Beta(1,3) and
    // Beta(2,2) temperature distributions. These closed-form moments are
    // independent of the production spline integration algorithm.
    for (nodal, average, derivative) in [
        (
            [300.0, 300.0, 300.0, 400.0],
            13.125,
            [0.03125, 0.03125, 0.03125, 0.15625],
        ),
        (
            [300.0, 300.0, 400.0, 400.0],
            28.75,
            [0.15625, 0.15625, 0.34375, 0.34375],
        ),
    ] {
        let model = curve();
        let moments = model.p1_conductivity(None, nodal, true).unwrap();
        close(moments.tensor[0][0], average, 2e-14);
        for (actual, expected) in moments.derivatives.iter().zip(derivative) {
            close(actual[0][0], expected, 2e-14);
        }
        assert_eq!(
            moments.tensor,
            model.p1_conductivity(None, nodal, false).unwrap().tensor
        );
        assert!(model.p1_is_smooth(nodal).unwrap());
    }
}

#[test]
fn constant_affine_and_uniform_kink_paths_keep_their_declared_semantics() {
    let nodal = [300.0, 320.0, 360.0, 400.0];
    let constant =
        ConductivityModel::constant_tensor([[5.0, 1.0, 0.0], [1.0, 4.0, 0.5], [0.0, 0.5, 3.0]])
            .unwrap();
    let values = constant.p1_conductivity(None, nodal, true).unwrap();
    assert_eq!(values.tensor, constant.tensor_at(345.0).unwrap());
    assert_eq!(values.derivatives, [[[0.0; 3]; 3]; 4]);
    let affine = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 10.0), (450.0, 30.0)]).unwrap(),
    );
    let values = affine.p1_conductivity(None, nodal, true).unwrap();
    assert_eq!(values.tensor, affine.tensor_at(345.0).unwrap());
    for derivative in values.derivatives {
        assert_eq!(derivative[0][0], 0.025);
    }

    let kink = curve();
    let values = kink.p1_conductivity(None, [350.0; 4], true).unwrap();
    assert_eq!(values.tensor[0][0], 10.0);
    for derivative in values.derivatives {
        assert_eq!(derivative[0][0], 0.5);
    }
    assert!(!kink.p1_is_smooth([350.0; 4]).unwrap());
    assert!(!kink.p1_is_smooth([250.0, 300.0, 400.0, 400.0]).unwrap());
    assert!(matches!(
        kink.p1_conductivity(None, [200.0, 300.0, 300.0, 400.0], true),
        Err(ConductionError::OutsideTemperatureSpan { .. })
    ));
}

#[test]
fn adjacent_float_temperatures_retain_the_same_rescaled_beta_moments() {
    let low = 300.0_f64;
    let mut middle = low;
    for _ in 0..8 {
        middle = middle.next_up();
    }
    let mut high = middle;
    for _ in 0..8 {
        high = high.next_up();
    }
    let model = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(low, 10.0), (middle, 10.0), (high, 110.0)])
            .unwrap(),
    );
    let values = model
        .p1_conductivity(None, [low, low, low, high], true)
        .unwrap();
    close(values.tensor[0][0], 13.125, 2e-14);
    for i in 0..4 {
        close(
            values.derivatives[i][0][0] * (high - low),
            if i == 3 { 15.625 } else { 3.125 },
            2e-14,
        );
    }
}

#[test]
fn orthotropic_rotation_preserves_integrated_components_and_nodal_derivatives() {
    let values = anisotropic()
        .p1_conductivity(None, [300.0, 300.0, 300.0, 400.0], true)
        .unwrap();
    for i in 0..2 {
        close(values.tensor[i][i], (13.125 + 29.6875) / 2.0, 2e-14);
        close(values.tensor[i][1 - i], (13.125 - 29.6875) / 2.0, 2e-14);
    }
    assert_eq!(values.tensor[2][2], 7.0);
    for c in 0..4 {
        let d = if c == 3 { 0.15625 } else { 0.03125 };
        close(values.derivatives[c][0][0], 0.45 * d, 2e-14);
        close(values.derivatives[c][0][1], 0.55 * d, 2e-14);
        assert_eq!(values.derivatives[c][2][2], 0.0);
    }
}

#[test]
fn operator_flux_and_material_pullback_use_the_same_physical_average() {
    with_cx(|_, cx| {
        let mesh = mesh();
        let material = curve();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let source = ScalarField::Uniform(0.0);
        let temperature = [300.0, 300.0, 300.0, 400.0];
        let system =
            assemble_operator(cx, &mesh, &boundary, &material, &source, &temperature).unwrap();
        let residual = full_residual(&system, &temperature);
        for (actual, expected) in residual.iter().zip([-218.75, 0.0, 0.0, 218.75]) {
            close(*actual, expected, 5e-13);
        }
        let flux = element_heat_flux(&mesh, &material, &temperature).unwrap();
        close(flux[0][2], -1312.5, 2e-14);
        let assigned = ElementMaterials::new(
            MaterialTable::new([(MaterialId(3), material.clone())]).unwrap(),
            vec![MaterialId(3)],
        )
        .unwrap();
        assert_eq!(
            flux,
            element_heat_flux_assigned(&mesh, &assigned, &temperature).unwrap()
        );
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: Some(&assigned),
            source: &source,
        };
        let gradient = RobinResponse::conductivity_scale_pullback(
            cx,
            problem,
            &temperature,
            &[0.3, -0.2, 0.4, 0.1],
            1,
        )
        .unwrap();
        close(gradient[0], 43.75, 2e-13);
        // The mean now lands exactly on the source kink, but the isotherm
        // has zero volume: the integrated law has a unique smooth tangent.
        let gradient = RobinResponse::conductivity_scale_pullback(
            cx,
            problem,
            &[300.0, 300.0, 400.0, 400.0],
            &[0.3, -0.2, 0.4, 0.1],
            1,
        )
        .unwrap();
        close(gradient[0], 575.0 / 12.0, 3e-13);
    });
}

#[test]
fn integrated_temperature_and_enthalpy_jacobians_match_perturbed_residuals() {
    with_cx(|_, cx| {
        let mesh = mesh();
        let material = anisotropic();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .remainder("film", ThermalBc::robin(1.3, 295.0).unwrap())
            .unwrap()
            .finish()
            .unwrap();
        let temperature = [310.0, 335.0, 370.0, 410.0];
        let chart = [0.002, 0.0, 0.004, 0.001];
        let scale = 1.7;
        let scale_derivative = [0.03, -0.02, 0.01, 0.04];
        let source = ScalarField::Uniform(0.0);
        let temp_jacobian =
            assemble_jacobian(cx, &mesh, &boundary, &material, &temperature).unwrap();
        let enthalpy_jacobian = assemble_enthalpy_jacobian_with_optional_interfaces(
            cx,
            &mesh,
            &boundary,
            &material,
            &temperature,
            &chart,
            Some(&[scale]),
            Some(&[scale_derivative]),
            None,
            None,
            None,
        )
        .unwrap();
        for enthalpy in [false, true] {
            let jacobian = if enthalpy {
                &enthalpy_jacobian
            } else {
                &temp_jacobian
            };
            for c in 0..4 {
                let epsilon = if enthalpy { 0.002 } else { 0.0002 };
                let mut columns = [0.0; 4];
                columns[c] = 1.0;
                let mut actual = [0.0; 4];
                jacobian.spmv(&columns, &mut actual);
                let evaluate = |delta: f64| {
                    let mut perturbed = temperature;
                    perturbed[c] += delta * if enthalpy { chart[c] } else { 1.0 };
                    let scales = [scale + delta * scale_derivative[c]];
                    let system = assemble_operator_scaled_with_interfaces(
                        cx,
                        &mesh,
                        &boundary,
                        &material,
                        &source,
                        &perturbed,
                        enthalpy.then_some(&scales[..]),
                        None,
                        None,
                    )
                    .unwrap();
                    full_residual(&system, &perturbed)
                };
                let plus = evaluate(epsilon);
                let minus = evaluate(-epsilon);
                for a in 0..4 {
                    close(actual[a], (plus[a] - minus[a]) / (2.0 * epsilon), 3e-7);
                }
            }
        }
    });
}

#[test]
fn crossing_material_knots_observes_cancellation() {
    with_cx(|gate, cx| {
        gate.request();
        assert!(matches!(
            curve().p1_conductivity(Some(cx), [300.0, 300.0, 300.0, 400.0], true),
            Err(ConductionError::Cancelled {
                stage: "integrate-conductivity-knots",
                ..
            })
        ));
    });
}

fn skew_mesh() -> ConductionMesh {
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [0.0, 0.0, 0.0],
            [3.0, 0.0, 0.0],
            [1.0, 7.0, 0.0],
            [2.0, 4.0, 11.0],
        ],
    )
    .unwrap()
}

#[test]
fn uniform_skew_tetrahedron_has_no_conductive_flux_or_material_sensitivity() {
    with_cx(|_, cx| {
        let mesh = skew_mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let source = ScalarField::Uniform(0.0);
        for material in [curve(), anisotropic()] {
            let temperature = [350.1; 4];
            assert_eq!(
                element_heat_flux(&mesh, &material, &temperature).unwrap(),
                [[0.0; 3]]
            );
            let assigned = ElementMaterials::new(
                MaterialTable::new([(MaterialId(7), material.clone())]).unwrap(),
                vec![MaterialId(7)],
            )
            .unwrap();
            assert_eq!(
                element_heat_flux_assigned(&mesh, &assigned, &temperature).unwrap(),
                [[0.0; 3]]
            );
            let problem = ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                element_materials: Some(&assigned),
                source: &source,
            };
            let dual = [0.3, -0.2, 0.4, 0.1];
            assert_eq!(
                RobinResponse::conductivity_scale_pullback(cx, problem, &temperature, &dual, 1)
                    .unwrap(),
                [0.0]
            );
            assert_eq!(
                RobinResponse::conductivity_rotation_pullback(cx, problem, &temperature, &dual, 1)
                    .unwrap(),
                [[0.0; 3]]
            );
            let operator =
                assemble_operator(cx, &mesh, &boundary, &material, &source, &temperature).unwrap();
            let jacobian =
                assemble_jacobian(cx, &mesh, &boundary, &material, &temperature).unwrap();
            let latent = assemble_enthalpy_jacobian_with_optional_interfaces(
                cx, &mesh, &boundary, &material, &temperature, &[0.0; 4], Some(&[1.7]),
                Some(&[[0.3, -0.2, 0.4, 0.1]]), None, None, None,
            )
            .unwrap();
            for row in 0..4 {
                assert_eq!(
                    operator.operator.row(row), jacobian.row(row),
                    "uniform state has no K'(T)*grad(T) correction"
                );
                assert!(
                    latent.row(row).1.iter().all(|&entry| entry == 0.0),
                    "phase conductivity cannot transport heat along an isothermal latent plateau"
                );
            }
        }
    });
}

#[test]
fn affine_skew_flux_keeps_absolute_constitutive_temperature() {
    let mesh = skew_mesh();
    let material = ConductivityModel::isotropic(
        ConductivityTable::declared_curve(vec![(250.0, 10.0), (450.0, 210.0)]).unwrap(),
    );
    // T=300+2x-3y+4z. Its mean is 305.75 K and k(T)=T-240.
    // Translating T changes k, even though the physical gradient is unchanged.
    for offset in [0.0, 32.0] {
        let temperature = [300.0, 306.0, 281.0, 336.0].map(|t| t + offset);
        let flux = element_heat_flux(&mesh, &material, &temperature).unwrap();
        for (actual, slope) in flux[0].into_iter().zip([2.0, -3.0, 4.0]) {
            close(actual, -(65.75 + offset) * slope, 2e-14);
        }
    }
}
