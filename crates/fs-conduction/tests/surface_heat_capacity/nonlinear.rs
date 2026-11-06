//! G1/G3: original film k(T), exact face moments, shared storage and tangents.

use super::*;
use fs_blake3::ContentHash;
use fs_conduction::adjoint::{ConductivityDesign, robin::RobinLinearization};
use fs_conduction::transient::enthalpy::{
    EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
};
use fs_conduction::{InitialGuess, ThermalBoundary};
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
use fs_solver::NewtonKrylovConfig;

fn right_tet() -> (ConductionMesh, usize) {
    let mesh = ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]],
    )
    .unwrap();
    let face = mesh
        .boundary()
        .iter()
        .position(|f| !f.vertices.contains(&3))
        .unwrap();
    (mesh, face)
}

fn problem<'a>(
    mesh: &'a ConductionMesh,
    boundary: &'a ThermalBoundary,
    material: &'a ConductivityModel,
    source: &'a ScalarField,
) -> ConductionProblem<'a> {
    ConductionProblem {
        mesh,
        boundary,
        material,
        source,
        element_materials: None,
    }
}

fn kink_table() -> ConductivityTable {
    ConductivityTable::declared_curve(vec![(300.0, 1.0), (350.0, 2.0), (400.0, 4.0)]).unwrap()
}

#[test]
fn crossed_knots_retain_directional_flux_and_the_complete_face_tangent() {
    with_cx(|cx| {
        let (mesh, face) = right_tet();
        let film = SurfaceHeatConduction::from_orthotropic_conductivity(
            face,
            [[1.0, 0.0, 0.0], [0.0, 1.0, 0.0]],
            [kink_table(), ConductivityTable::declared(5.0).unwrap()],
            0.2,
        )
        .unwrap();
        assert!(film.is_temperature_dependent());
        assert_eq!(film.conductivity_tables().unwrap()[0], kink_table());
        let base = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let coated = base
            .clone()
            .with_surface_heat_conduction(cx, &mesh, std::slice::from_ref(&film))
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let temperature = [300.0, 350.0, 400.0, 320.0];
        // Splitting on T=350 gives <kx>=13/6, not k(<T>)=2. The original
        // xy face has area 1/2, grad(T)=[50,100,0], thickness=1/5.
        let flux = film.mean_heat_flux_w_per_m(&mesh, &temperature).unwrap();
        for (actual, expected) in flux.into_iter().zip([-65.0 / 3.0, -100.0, 0.0]) {
            close(actual, expected, 2e-12);
        }
        let residual = |boundary: &ThermalBoundary| {
            full_residual(
                &assemble_operator(cx, &mesh, boundary, &material, &source, &temperature).unwrap(),
                &temperature,
            )
        };
        let bare = residual(&base);
        let total = residual(&coated);
        for i in 0..4 {
            close(
                total[i] - bare[i],
                [-365.0 / 6.0, 65.0 / 6.0, 50.0, 0.0][i],
                2e-10,
            );
        }
        close(
            total.iter().sum::<f64>() - bare.iter().sum::<f64>(),
            0.0,
            2e-10,
        );

        // Independent differentiated subtriangle integrals give
        // d<kx>/dT=[1/120,1/100,7/600]. For this direction, the complete
        // coating tangent action is [107/80,-47/80,-3/4,0].
        let direction = [2.0, -1.0, 0.5, 4.0];
        let mut full_action = [0.0; 4];
        let mut bare_action = [0.0; 4];
        assemble_jacobian(cx, &mesh, &coated, &material, &temperature)
            .unwrap()
            .spmv(&direction, &mut full_action);
        assemble_jacobian(cx, &mesh, &base, &material, &temperature)
            .unwrap()
            .spmv(&direction, &mut bare_action);
        for i in 0..4 {
            close(
                full_action[i] - bare_action[i],
                [107.0 / 80.0, -47.0 / 80.0, -0.75, 0.0][i],
                3e-12,
            );
        }
        assert!(matches!(
            assemble_operator(
                cx,
                &mesh,
                &coated,
                &material,
                &source,
                &[300.0, 350.0, 401.0, 320.0]
            ),
            Err(ConductionError::OutsideTemperatureSpan { .. })
        ));
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                assemble_jacobian(cancelled, &mesh, &coated, &material, &temperature),
                Err(ConductionError::Cancelled { .. })
            ));
        });
    });
}

#[test]
fn film_endpoint_shares_nonlinear_be_enthalpy_and_original_surface_storage() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        // k*t = 3 + 0.4*(T-300). The three hot faces at the target endpoint
        // have mean T=300+1/3; the fourth face is uniform and transports none.
        let table = ConductivityTable::declared_curve(vec![(295.0, 2.0), (305.0, 10.0)]).unwrap();
        let template = SurfaceHeatConduction::from_isotropic_conductivity(0, table, 0.5).unwrap();
        let films: Vec<_> = (0..4).map(|face| template.on_boundary_face(face)).collect();
        assert!(std::ptr::eq(
            template.conductivity_tables().unwrap().as_ptr(),
            films[3].conductivity_tables().unwrap().as_ptr()
        ));
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap()
            .with_surface_heat_conduction(cx, &mesh, &films)
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let p = problem(&mesh, &boundary, &material, &source);
        let surface_storage: Vec<_> = (0..4)
            .map(|face| SurfaceHeatCapacity::uniform(face, areal(1.0)))
            .collect();
        let engine =
            BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(4.0).unwrap())
                .unwrap()
                .with_surface_heat_capacity(cx, &surface_storage)
                .unwrap();
        let dt = 0.2;
        let target = [303.0, 299.0, 299.0, 299.0];
        // Exact regular-tet modal balance, including the old and new sensible
        // film energy. Manufacture history for target amplitude a=1.
        let rate = (2.0 + 2.0 * 3.0_f64.sqrt() * (3.0 + 0.4 / 3.0)) / (4.0 + 3.0 * 3.0_f64.sqrt());
        let a_old = 1.0 + dt * rate;
        let old = [
            300.0 + 3.0 * a_old,
            300.0 - a_old,
            300.0 - a_old,
            300.0 - a_old,
        ];
        assert!(engine.advance(cx, p, None, &old, dt, config()).is_err());
        assert!(matches!(
            ConductivityDesign::new(p, config().linear),
            Err(ConductionError::Conductivity { .. })
        ));
        let step = engine
            .advance_nonlinear(
                cx,
                p,
                None,
                &old,
                dt,
                config(),
                NonlinearStepConfig::default(),
            )
            .unwrap();
        for i in 0..4 {
            close(step.step.temperature[i], target[i], 2e-8);
        }
        close(step.step.stored_energy_change_j, 0.0, 1e-8);
        close(step.step.energy_residual_j, 0.0, 1e-8);

        // The same independently specified endpoint with h=2*(T-290), rho=2,
        // dT/dh=1/2 exercises the complete coating tangent in enthalpy space.
        let chart = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x71; 32]),
            SolidLiquidPhase::Solid,
            [(0.0, 290.0), (40.0, 310.0)]
                .map(|(h, t)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: h,
                    temperature_k: t,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 2.0,
                })
                .to_vec(),
        )
        .unwrap();
        let enthalpy = EnthalpyBackwardEuler::uniform(
            cx,
            &mesh,
            &chart,
            2.0,
            EnthalpyBudget {
                max_vertices: 4,
                max_elements: 1,
            },
        )
        .unwrap()
        .with_surface_heat_capacity(cx, &surface_storage)
        .unwrap();
        let old_h = old.map(|t| 2.0 * (t - 290.0));
        let enthalpy_step = enthalpy
            .advance(
                cx,
                p,
                None,
                &old_h,
                dt,
                EnthalpyStepConfig {
                    newton: NewtonKrylovConfig {
                        absolute_tolerance: 1e-10,
                        relative_tolerance: 1e-12,
                        linear_restart: 4,
                        max_linear_cycles: 8,
                        forcing_minimum: 1e-12,
                        forcing_maximum: 1e-3,
                        ..NewtonKrylovConfig::default()
                    },
                    max_newton_iterations: 32,
                    energy_tolerance_j: 1e-8,
                },
            )
            .unwrap();
        for i in 0..4 {
            close(enthalpy_step.temperature[i], target[i], 2e-8);
            close(
                enthalpy_step.specific_enthalpy_j_kg[i],
                2.0 * (target[i] - 290.0),
                4e-8,
            );
        }
        close(enthalpy_step.stored_energy_change_j, 0.0, 1e-8);
        close(enthalpy_step.energy_residual_j, 0.0, 1e-8);
        assert_eq!(old_h, old.map(|t| 2.0 * (t - 290.0)));
    });
}

#[test]
fn flat_table_limit_and_nonlinear_adjoint_source_admission_are_explicit() {
    with_cx(|cx| {
        let (mesh, face) = right_tet();
        let base = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let flat = SurfaceHeatConduction::from_isotropic_conductivity(
            face,
            ConductivityTable::declared_curve(vec![(300.0, 15.0), (400.0, 15.0)]).unwrap(),
            0.2,
        )
        .unwrap();
        assert!(!flat.is_temperature_dependent());
        let tabulated = base
            .clone()
            .with_surface_heat_conduction(cx, &mesh, &[flat])
            .unwrap();
        let constant = base
            .with_surface_heat_conduction(
                cx,
                &mesh,
                &[SurfaceHeatConduction::uniform(face, 3.0).unwrap()],
            )
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let temperature = [320.0, 350.0, 380.0, 330.0];
        let a = assemble_operator(cx, &mesh, &tabulated, &material, &source, &temperature).unwrap();
        let b = assemble_operator(cx, &mesh, &constant, &material, &source, &temperature).unwrap();
        for row in 0..4 {
            assert_eq!(a.operator.row(row), b.operator.row(row));
        }

        // A constant bulk still needs the nonlinear film Jacobian. Original
        // source endpoints and a uniform true kink refuse classical gradients.
        for (t, smooth) in [(325.0, true), (350.0, false), (300.0, false)] {
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region("bath", |_| true, ThermalBc::robin(1.0, t).unwrap())
                .unwrap()
                .finish()
                .unwrap()
                .with_surface_heat_conduction(
                    cx,
                    &mesh,
                    &[
                        SurfaceHeatConduction::from_isotropic_conductivity(face, kink_table(), 0.2)
                            .unwrap(),
                    ],
                )
                .unwrap();
            let response = RobinLinearization::new(
                cx,
                problem(&mesh, &boundary, &material, &source),
                SolveConfig {
                    initial: InitialGuess::Uniform(t),
                    ..SolveConfig::default()
                },
                &["bath"],
            )
            .unwrap();
            assert!(response.uses_nonlinear_jacobian());
            assert_eq!(response.has_smooth_material_tangent(), smooth);
            let mut direction = response.zero_direction();
            direction.references_k[0] = 1.0;
            if smooth {
                for change in response.apply(cx, &direction).unwrap().temperature_k {
                    close(change, 1.0, 1e-8);
                }
            } else {
                assert!(response.apply(cx, &direction).is_err());
            }
        }
    });
}
