//! Thin bonded films preserve storage, tangential transport and source domains.
mod support;

#[path = "surface_heat_capacity/enthalpy.rs"]
mod enthalpy;

#[path = "surface_heat_capacity/tensor.rs"]
mod tensor;

#[path = "surface_heat_capacity/nonlinear.rs"]
mod nonlinear;

use fs_conduction::assemble::{assemble_jacobian, assemble_operator, full_residual};
use fs_conduction::bc::history::ThermalBoundaryHistory;
use fs_conduction::bc::history::march::{
    HeatSourceHistory, ThermalHistoryProblem, ThermalHistoryRun,
};
use fs_conduction::fixtures::{box_grid, on_box_face};
use fs_conduction::transient::backward_euler::{
    ArealHeatCapacity, BackwardEuler, NonlinearStepConfig, StepConfig, SurfaceHeatCapacity,
};
use fs_conduction::transient::{MaterialHeatCapacity, VolumetricHeatCapacity};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, SolveConfig, SurfaceHeatConduction, ThermalBc,
    ThermalBoundaryBuilder, VOLUMETRIC_SOURCE_DIMS,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS,
    SPECIFIC_HEAT_CAPACITY_PROPERTY, THERMAL_CONDUCTIVITY_DIMS, THERMAL_CONDUCTIVITY_PROPERTY,
};
use fs_qty::Density;
use fs_rep_mesh::TetComplex;
use fs_scenario::scenario::Environment;
use fs_scenario::signal::{Interp, TimeSignal};
use support::{with_cancelled_cx, with_cx};

fn mesh() -> ConductionMesh {
    let (complex, positions) = box_grid([2, 2, 2], [1.0; 3]);
    ConductionMesh::new(complex, positions).unwrap()
}

fn config() -> StepConfig {
    StepConfig {
        linear: LinearConfig {
            tolerance: 1e-11,
            max_iterations: 512,
            restart: 16,
        },
        energy_tolerance_j: 1e-8,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} differs from {expected:e} by more than {tolerance:e}"
    );
}

fn areal(value: f64) -> ArealHeatCapacity {
    ArealHeatCapacity::declared(value).unwrap()
}

#[test]
fn uniform_heating_resolves_uniform_and_p1_surface_storage() {
    with_cx(|cx| {
        let mesh = mesh();
        let c_volume = 3.0;
        let rate = 2.0;
        let dt = 0.25;
        let old = vec![300.0; mesh.vertex_count()];
        for varying in [false, true] {
            let c_surface: Vec<f64> = mesh
                .positions()
                .iter()
                .map(|point| {
                    if varying {
                        2.0 + point[0] + 2.0 * point[1]
                    } else {
                        2.0
                    }
                })
                .collect();
            let surfaces: Vec<_> = mesh
                .boundary()
                .iter()
                .enumerate()
                .map(|(face, geometry)| {
                    if varying {
                        SurfaceHeatCapacity::nodal(
                            face,
                            geometry
                                .vertices
                                .map(|vertex| areal(c_surface[vertex as usize])),
                        )
                    } else {
                        SurfaceHeatCapacity::uniform(face, areal(2.0))
                    }
                })
                .collect();
            let boundary = ThermalBoundaryBuilder::new(&mesh)
                .region(
                    "film power",
                    |_| true,
                    ThermalBc::Neumann {
                        outward_flux: ScalarField::Nodal(
                            c_surface.iter().map(|c| -rate * c).collect(),
                        ),
                    },
                )
                .unwrap()
                .finish()
                .unwrap();
            let source = ScalarField::Uniform(c_volume * rate);
            let engine = BackwardEuler::uniform(
                cx,
                &mesh,
                VolumetricHeatCapacity::declared(c_volume).unwrap(),
            )
            .unwrap()
            .with_surface_heat_capacity(cx, &surfaces)
            .unwrap();
            // Integral of a linear function over the six faces of the unit
            // cube: the area-weighted means of x and y are both one half.
            let surface_capacity = if varying { 6.0 * 3.5 } else { 6.0 * 2.0 };
            for nonlinear in [false, true] {
                let material = if nonlinear {
                    ConductivityModel::isotropic(
                        ConductivityTable::declared_curve(vec![(280.0, 1.0), (340.0, 7.0)])
                            .unwrap(),
                    )
                } else {
                    ConductivityModel::isotropic_declared(4.0).unwrap()
                };
                let problem = ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                };
                let step = if nonlinear {
                    engine
                        .advance_nonlinear(
                            cx,
                            problem,
                            None,
                            &old,
                            dt,
                            config(),
                            NonlinearStepConfig::default(),
                        )
                        .unwrap()
                        .step
                } else {
                    engine
                        .advance(cx, problem, None, &old, dt, config())
                        .unwrap()
                };
                for &temperature in &step.temperature {
                    close(temperature, 300.0 + rate * dt, 2e-9);
                }
                close(
                    step.stored_energy_change_j,
                    (c_volume + surface_capacity) * rate * dt,
                    1e-8,
                );
                close(step.source_w, c_volume * rate, 1e-11);
                close(step.neumann_out_w, -surface_capacity * rate, 1e-10);
                close(step.dirichlet_in_w, 0.0, 1e-12);
                close(step.energy_residual_j, 0.0, config().energy_tolerance_j);
            }
        }
        assert_eq!(old, vec![300.0; mesh.vertex_count()]);
    });
}

#[test]
fn prescribed_surface_jump_adds_its_storage_to_the_reaction() {
    with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "prescribed coating",
                |_| true,
                ThermalBc::dirichlet(310.0).unwrap(),
            )
            .unwrap()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let old = vec![300.0; mesh.vertex_count()];
        let volume = VolumetricHeatCapacity::declared(5.0).unwrap();
        let base = BackwardEuler::uniform(cx, &mesh, volume).unwrap();
        let surfaces: Vec<_> = (0..mesh.boundary().len())
            .map(|face| SurfaceHeatCapacity::uniform(face, areal(7.0)))
            .collect();
        let coated = BackwardEuler::uniform(cx, &mesh, volume)
            .unwrap()
            .with_surface_heat_capacity(cx, &surfaces)
            .unwrap();
        let dt = 0.5;
        let a = base
            .advance_prescribed(cx, problem, None, &old, dt, config())
            .unwrap();
        let b = coated
            .advance_prescribed(cx, problem, None, &old, dt, config())
            .unwrap();
        // The added capacity is only on prescribed nodes. It changes their
        // heat reaction, not the free-node endpoint or the conduction flux.
        assert_eq!(a.temperature, b.temperature);
        let extra_storage = 6.0 * 7.0 * (310.0 - 300.0);
        close(
            b.stored_energy_change_j - a.stored_energy_change_j,
            extra_storage,
            1e-9,
        );
        close(
            b.dirichlet_in_w - a.dirichlet_in_w,
            extra_storage / dt,
            1e-9,
        );
        close(b.energy_residual_j, 0.0, config().energy_tolerance_j);
        assert!(
            coated
                .advance(cx, problem, None, &old, dt, config())
                .is_err()
        );
        assert_eq!(old, vec![300.0; mesh.vertex_count()]);
    });
}

#[test]
fn zero_surface_storage_preserves_the_original_solution_exactly() {
    with_cx(|cx| {
        let mesh = mesh();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(6.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &material,
            element_materials: None,
            source: &source,
        };
        let old: Vec<_> = mesh.positions().iter().map(|p| 300.0 + p[0]).collect();
        let volume = VolumetricHeatCapacity::declared(5.0).unwrap();
        let base = BackwardEuler::uniform(cx, &mesh, volume)
            .unwrap()
            .advance(cx, problem, None, &old, 0.5, config())
            .unwrap();
        let zero = [SurfaceHeatCapacity::uniform(0, areal(0.0))];
        for surfaces in [&[][..], &zero[..]] {
            let step = BackwardEuler::uniform(cx, &mesh, volume)
                .unwrap()
                .with_surface_heat_capacity(cx, surfaces)
                .unwrap()
                .advance(cx, problem, None, &old, 0.5, config())
                .unwrap();
            assert_eq!(step.temperature, base.temperature);
            assert_eq!(
                step.stored_energy_change_j.to_bits(),
                base.stored_energy_change_j.to_bits()
            );
            assert_eq!(
                step.energy_residual_j.to_bits(),
                base.energy_residual_j.to_bits()
            );
            assert_eq!(
                step.relative_residual.to_bits(),
                base.relative_residual.to_bits()
            );
            assert_eq!(step.krylov_iterations, base.krylov_iterations);
        }
    });
}

fn source_capacity(low: f64, high: f64) -> MaterialHeatCapacity {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (DENSITY_PROPERTY, 4.0, Density::DIMS),
        (
            SPECIFIC_HEAT_CAPACITY_PROPERTY,
            5.0,
            SPECIFIC_HEAT_CAPACITY_DIMS,
        ),
        (
            THERMAL_CONDUCTIVITY_PROPERTY,
            2.0,
            THERMAL_CONDUCTIVITY_DIMS,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", low, high),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: Vec::new(),
                provenance: Provenance {
                    source: "manufactured coating; not measured material data".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "solid".into(),
            process: "surface storage regression".into(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .unwrap();
    MaterialHeatCapacity::from_card(
        &card,
        &QueryPoint::new()
            .with("T", f64::midpoint(low, high))
            .unwrap(),
        low,
        high,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

#[test]
fn material_coatings_retain_sources_and_enforce_their_temperature_domains() {
    with_cx(|cx| {
        let mesh = mesh();
        let material_capacity = source_capacity(290.0, 310.0);
        let surface =
            SurfaceHeatCapacity::from_material(0, &material_capacity, [0.0, 0.1, 0.2]).unwrap();
        assert_eq!(surface.boundary_face(), 0);
        assert_eq!(
            surface
                .nodal_capacity()
                .map(ArealHeatCapacity::value_j_per_m2_k),
            [0.0, 2.0, 4.0]
        );
        assert_eq!(surface.material(), Some(&material_capacity));
        let volume = VolumetricHeatCapacity::declared(3.0).unwrap();
        let engine = BackwardEuler::uniform(cx, &mesh, volume)
            .unwrap()
            .with_surface_heat_capacity(cx, &[surface])
            .unwrap();
        assert_eq!(engine.material_capacities(), &[material_capacity.clone()]);
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let problem = ConductionProblem {
            mesh: &mesh,
            boundary: &boundary,
            material: &conductivity,
            element_materials: None,
            source: &source,
        };
        let valid = vec![300.0; mesh.vertex_count()];
        assert!(
            engine
                .advance(cx, problem, None, &valid, 0.1, config())
                .is_ok()
        );
        let mut invalid = valid.clone();
        // Even the corner with zero local thickness receives a positive P1
        // storage moment from the other corners of the coated triangle.
        invalid[mesh.boundary()[0].vertices[0] as usize] = 320.0;
        assert!(
            engine
                .advance(cx, problem, None, &invalid, 0.1, config())
                .is_err()
        );
        let disjoint = source_capacity(320.0, 340.0);
        let other = SurfaceHeatCapacity::from_material(0, &disjoint, [0.1; 3]).unwrap();
        assert!(engine.with_surface_heat_capacity(cx, &[other]).is_err());
        // An absent layer contributes neither storage nor a source restriction.
        let absent = SurfaceHeatCapacity::from_material(0, &material_capacity, [0.0; 3]).unwrap();
        let no_layer = BackwardEuler::uniform(cx, &mesh, volume)
            .unwrap()
            .with_surface_heat_capacity(cx, &[absent])
            .unwrap();
        assert!(no_layer.material_capacities().is_empty());
        assert!(
            no_layer
                .advance(
                    cx,
                    problem,
                    None,
                    &vec![320.0; mesh.vertex_count()],
                    0.1,
                    config()
                )
                .is_ok()
        );
    });
}

#[test]
fn invalid_faces_coefficients_and_cancellation_refuse_preparation() {
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(ArealHeatCapacity::declared(value).is_err());
    }
    with_cx(|cx| {
        let mesh = mesh();
        let volume = VolumetricHeatCapacity::declared(3.0).unwrap();
        let repeated = SurfaceHeatCapacity::uniform(0, areal(1.0));
        assert!(
            BackwardEuler::uniform(cx, &mesh, volume)
                .unwrap()
                .with_surface_heat_capacity(cx, &[repeated, repeated])
                .is_err()
        );
        let absent = SurfaceHeatCapacity::uniform(mesh.boundary().len(), areal(1.0));
        assert!(
            BackwardEuler::uniform(cx, &mesh, volume)
                .unwrap()
                .with_surface_heat_capacity(cx, &[absent])
                .is_err()
        );
        let material_capacity = source_capacity(290.0, 310.0);
        for thickness in [-1.0, f64::NAN, f64::INFINITY] {
            assert!(
                SurfaceHeatCapacity::from_material(0, &material_capacity, [thickness; 3]).is_err()
            );
        }
        let engine = BackwardEuler::uniform(cx, &mesh, volume).unwrap();
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                engine.with_surface_heat_capacity(cancelled, &[repeated]),
                Err(ConductionError::Cancelled { .. })
            ));
        });
    });
}

#[test]
fn coated_slab_retains_exact_bulk_plus_sheet_heat_and_tangential_flux() {
    with_cx(|cx| {
        let (complex, positions) = box_grid([3, 2, 2], [1.0; 3]);
        let mesh = ConductionMesh::new(complex, positions).unwrap();
        let base = ThermalBoundaryBuilder::new(&mesh)
            .region(
                "cold",
                |f| on_box_face(f.centroid[0], 0.0),
                ThermalBc::dirichlet(300.0).unwrap(),
            )
            .unwrap()
            .region(
                "hot",
                |f| on_box_face(f.centroid[0], 1.0),
                ThermalBc::dirichlet(310.0).unwrap(),
            )
            .unwrap()
            .adiabatic_remainder()
            .finish()
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        for varying in [false, true] {
            let surfaces: Vec<_> = mesh
                .boundary()
                .iter()
                .enumerate()
                .filter(|(_, f)| {
                    !on_box_face(f.centroid[0], 0.0) && !on_box_face(f.centroid[0], 1.0)
                })
                .map(|(slot, face)| {
                    let values = face.vertices.map(|v| {
                        let p = mesh.positions()[v as usize];
                        if varying {
                            2.0 + p[1] + 2.0 * p[2]
                        } else {
                            3.0
                        }
                    });
                    SurfaceHeatConduction::nodal(slot, values).unwrap()
                })
                .collect();
            let boundary = base
                .clone()
                .with_surface_heat_conduction(cx, &mesh, &surfaces)
                .unwrap();
            let solved = fs_conduction::solve(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    source: &source,
                    element_materials: None,
                },
                SolveConfig::default(),
            )
            .unwrap();
            for (&t, p) in solved.temperature.iter().zip(mesh.positions()) {
                close(t, 300.0 + 10.0 * p[0], 2e-8);
            }
            let system = assemble_operator(
                cx,
                &mesh,
                &boundary,
                &material,
                &source,
                &solved.temperature,
            )
            .unwrap();
            let reaction = full_residual(&system, &solved.temperature);
            let hot: f64 = mesh
                .positions()
                .iter()
                .zip(&reaction)
                .filter(|(p, _)| on_box_face(p[0], 1.0))
                .map(|(_, r)| r)
                .sum();
            let cold: f64 = mesh
                .positions()
                .iter()
                .zip(&reaction)
                .filter(|(p, _)| on_box_face(p[0], 0.0))
                .map(|(_, r)| r)
                .sum();
            // G = k*A + integral_perimeter(k_f*t). The varying film has
            // four side means 3, 4, 2.5, 4.5 W/K; it remains independent of x.
            let expected = if varying { 160.0 } else { 140.0 };
            close(hot, expected, 2e-7);
            close(cold, -expected, 2e-7);
            close(reaction.iter().sum(), 0.0, 1e-9);
            for surface in boundary.surface_heat_conduction() {
                let mean = surface.nodal_conductance_w_per_k().unwrap().iter().sum::<f64>() / 3.0;
                let flux = surface
                    .mean_heat_flux_w_per_m(&mesh, &solved.temperature)
                    .unwrap();
                close(flux[0], -10.0 * mean, 2e-7);
                close(flux[1], 0.0, 2e-7);
                close(flux[2], 0.0, 2e-7);
            }
        }
    });
}

fn regular_tetrahedron() -> ConductionMesh {
    // Volume 8/3, four equal face areas 2*sqrt(3), outward orientation.
    ConductionMesh::new(
        TetComplex::from_tets(4, vec![[0, 1, 2, 3]]),
        vec![
            [1.0, 1.0, 1.0],
            [-1.0, -1.0, 1.0],
            [-1.0, 1.0, -1.0],
            [1.0, -1.0, -1.0],
        ],
    )
    .unwrap()
}

fn coated_tetrahedron_boundary(
    cx: &fs_exec::Cx<'_>,
    mesh: &ConductionMesh,
) -> fs_conduction::ThermalBoundary {
    let films: Vec<_> = (0..mesh.boundary().len())
        .map(|face| SurfaceHeatConduction::uniform(face, 3.0).unwrap())
        .collect();
    ThermalBoundaryBuilder::new(mesh)
        .adiabatic_remainder()
        .finish()
        .unwrap()
        .with_surface_heat_conduction(cx, mesh, &films)
        .unwrap()
}

#[test]
fn coating_diffusion_with_surface_storage_matches_a_tetrahedron_eigenmode() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let old = [301.0, 299.0, 300.0, 300.0];
        let surfaces: Vec<_> = (0..mesh.boundary().len())
            .map(|face| SurfaceHeatCapacity::uniform(face, areal(6.0)))
            .collect();
        let stepper =
            BackwardEuler::uniform(cx, &mesh, VolumetricHeatCapacity::declared(4.0).unwrap())
                .unwrap()
                .with_surface_heat_capacity(cx, &surfaces)
                .unwrap();
        let step = stepper
            .advance(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    element_materials: None,
                    source: &source,
                },
                None,
                &old,
                0.1,
                config(),
            )
            .unwrap();
        // Exact P1 zero-sum mode: lambda=(k+2*sqrt(3)*kappa)/(c_v+3*sqrt(3)*C_s).
        // Both bulk and sheet conduction dissipate; sheet capacitance slows decay.
        let rate = (2.0 + 2.0 * 3.0_f64.sqrt() * 3.0) / (4.0 + 3.0 * 3.0_f64.sqrt() * 6.0);
        let amplitude = 1.0 / (1.0 + 0.1 * rate);
        for (&actual, expected) in
            step.temperature
                .iter()
                .zip([300.0 + amplitude, 300.0 - amplitude, 300.0, 300.0])
        {
            close(actual, expected, 2e-9);
        }
        close(step.temperature.iter().sum::<f64>() / 4.0, 300.0, 1e-11);
        close(step.stored_energy_change_j, 0.0, 1e-9);
        close(step.energy_residual_j, 0.0, config().energy_tolerance_j);
        close(step.dirichlet_in_w, 0.0, 1e-12);
        assert_eq!(old, [301.0, 299.0, 300.0, 300.0]);
    });
}

#[test]
fn enthalpy_transport_scales_the_same_coating_tangent_by_temperature_slopes() {
    use fs_blake3::ContentHash;
    use fs_conduction::transient::enthalpy::{
        EnthalpyBackwardEuler, EnthalpyBudget, EnthalpyStepConfig,
    };
    use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase};
    use fs_solver::NewtonKrylovConfig;
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        // h=2*(T-280), rho0=2 gives c_v=4 and dT/dh=1/2, not one.
        let curve = EquilibriumEnthalpyPhaseCurve::try_single_phase(
            ContentHash([0x62; 32]),
            SolidLiquidPhase::Solid,
            [(0.0, 280.0), (160.0, 360.0)]
                .into_iter()
                .map(|(h, t)| EnthalpyPhaseKnot {
                    specific_enthalpy_j_kg: h,
                    temperature_k: t,
                    liquid_mass_fraction: 0.0,
                    bulk_density_kg_m3: 2.0,
                })
                .collect(),
        )
        .unwrap();
        let engine = EnthalpyBackwardEuler::uniform(
            cx,
            &mesh,
            &curve,
            2.0,
            EnthalpyBudget {
                max_vertices: 4,
                max_elements: 1,
            },
        )
        .unwrap();
        let old = [42.0, 38.0, 40.0, 40.0];
        let step = engine
            .advance(
                cx,
                ConductionProblem {
                    mesh: &mesh,
                    boundary: &boundary,
                    material: &material,
                    source: &source,
                    element_materials: None,
                },
                None,
                &old,
                0.1,
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
        let rate = (2.0 + 2.0 * 3.0_f64.sqrt() * 3.0) / 4.0;
        let amplitude = 1.0 / (1.0 + 0.1 * rate);
        for (&t, expected) in
            step.temperature
                .iter()
                .zip([300.0 + amplitude, 300.0 - amplitude, 300.0, 300.0])
        {
            close(t, expected, 2e-8);
        }
        close(step.stored_energy_change_j, 0.0, 1e-8);
        close(step.energy_residual_j, 0.0, 1e-8);
        assert_eq!(old, [42.0, 38.0, 40.0, 40.0]);
    });
}

#[test]
fn endpoint_and_event_histories_preserve_original_coating_transport() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let history = ThermalBoundaryHistory::from_boundary(
            cx,
            &mesh,
            &boundary,
            &Environment::earth_lab(),
            vec![],
            0.0,
        )
        .unwrap();
        assert_eq!(
            history.at(cx, 0.1).unwrap().surface_heat_conduction(),
            boundary.surface_heat_conduction()
        );
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = HeatSourceHistory::Uniform(TimeSignal::Table {
            times: vec![0.0, 0.05, 0.1],
            values: vec![0.0; 3],
            dims: VOLUMETRIC_SOURCE_DIMS,
            interp: Interp::Hold,
        });
        let mut run = ThermalHistoryRun::new(
            cx,
            ThermalHistoryProblem {
                boundary: &history,
                material: &material,
                element_materials: None,
                interfaces: None,
                source,
            },
            VolumetricHeatCapacity::declared(4.0).unwrap(),
            &[301.0, 299.0, 300.0, 300.0],
            2,
        )
        .unwrap();
        let rate = (2.0 + 2.0 * 3.0_f64.sqrt() * 3.0) / 4.0;
        let mut amplitude = 1.0;
        for end in [0.05, 0.1] {
            let step = run.advance_event_step(cx, 0.1, 0.1, config()).unwrap();
            assert_eq!(step.end_time_s, end);
            amplitude /= 1.0 + 0.05 * rate;
            for (&t, expected) in step.report.solution().temperature.iter().zip([
                300.0 + amplitude,
                300.0 - amplitude,
                300.0,
                300.0,
            ]) {
                close(t, expected, 2e-9);
            }
            close(step.report.solution().energy_residual_j, 0.0, 1e-8);
        }
    });
}

#[test]
fn nonlinear_bulk_jacobian_retains_sheet_transport_and_its_constant_null_mode() {
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = coated_tetrahedron_boundary(cx, &mesh);
        let material = ConductivityModel::isotropic(
            ConductivityTable::declared_curve(vec![(280.0, 1.0), (320.0, 5.0)]).unwrap(),
        );
        let source = ScalarField::Uniform(0.0);
        let temperature = [295.0, 303.0, 300.0, 305.0];
        let direction = [1.0, -2.0, 0.5, -1.0];
        let jacobian = assemble_jacobian(cx, &mesh, &boundary, &material, &temperature).unwrap();
        let epsilon = 1e-3;
        let plus: Vec<_> = temperature
            .iter()
            .zip(direction)
            .map(|(&t, d)| t + epsilon * d)
            .collect();
        let minus: Vec<_> = temperature
            .iter()
            .zip(direction)
            .map(|(&t, d)| t - epsilon * d)
            .collect();
        let residual = |t: &[f64]| {
            full_residual(
                &assemble_operator(cx, &mesh, &boundary, &material, &source, t).unwrap(),
                t,
            )
        };
        let plus = residual(&plus);
        let minus = residual(&minus);
        for row in 0..4 {
            let (columns, values) = jacobian.row(row);
            let analytic: f64 = columns
                .iter()
                .zip(values)
                .map(|(&col, &v)| v * direction[col])
                .sum();
            close(analytic, (plus[row] - minus[row]) / (2.0 * epsilon), 2e-7);
        }
        for surface in boundary.surface_heat_conduction() {
            assert_eq!(
                surface.mean_heat_flux_w_per_m(&mesh, &[300.0; 4]).unwrap(),
                [0.0; 3]
            );
        }
    });
}

#[test]
fn coating_zero_duplicate_invalid_and_cancelled_inputs_preserve_admission() {
    for value in [-1.0, f64::NAN, f64::INFINITY] {
        assert!(SurfaceHeatConduction::uniform(0, value).is_err());
        assert!(SurfaceHeatConduction::nodal(0, [1.0, value, 1.0]).is_err());
    }
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let base = ThermalBoundaryBuilder::new(&mesh)
            .adiabatic_remainder()
            .finish()
            .unwrap();
        for films in [
            &[][..],
            &[SurfaceHeatConduction::uniform(0, 0.0).unwrap()][..],
        ] {
            assert_eq!(
                base.clone()
                    .with_surface_heat_conduction(cx, &mesh, films)
                    .unwrap(),
                base
            );
        }
        let film = SurfaceHeatConduction::uniform(0, 3.0).unwrap();
        assert!(
            base.clone()
                .with_surface_heat_conduction(cx, &mesh, &[film.clone(), film.clone()])
                .is_err()
        );
        let coated = base
            .clone()
            .with_surface_heat_conduction(cx, &mesh, std::slice::from_ref(&film))
            .unwrap();
        assert!(
            coated
                .with_surface_heat_conduction(cx, &mesh, std::slice::from_ref(&film))
                .is_err()
        );
        let absent = SurfaceHeatConduction::uniform(mesh.boundary().len(), 1.0).unwrap();
        assert!(
            base.clone()
                .with_surface_heat_conduction(cx, &mesh, &[absent])
                .is_err()
        );
        assert!(film.mean_heat_flux_w_per_m(&mesh, &[300.0; 3]).is_err());
        with_cancelled_cx(|cancelled| {
            assert!(matches!(
                base.with_surface_heat_conduction(cancelled, &mesh, &[film]),
                Err(ConductionError::Cancelled { .. })
            ))
        });
    });
}

#[cfg(feature = "thermal-verification")]
#[test]
fn continuum_mean_certificate_refuses_unrepresented_coating_physics() {
    use fs_conduction::verification::{
        ConductionBoundError, FluxBudget, TetError, bound_temperature_mean,
    };
    with_cx(|cx| {
        let mesh = regular_tetrahedron();
        let boundary = ThermalBoundaryBuilder::new(&mesh)
            .region("bath", |_| true, ThermalBc::robin(1.0, 300.0).unwrap())
            .unwrap()
            .finish()
            .unwrap()
            .with_surface_heat_conduction(
                cx,
                &mesh,
                &[SurfaceHeatConduction::uniform(0, 1.0).unwrap()],
            )
            .unwrap();
        let material = ConductivityModel::isotropic_declared(2.0).unwrap();
        let source = ScalarField::Uniform(0.0);
        let result = bound_temperature_mean(
            cx,
            ConductionProblem {
                mesh: &mesh,
                boundary: &boundary,
                material: &material,
                source: &source,
                element_materials: None,
            },
            &[300.0; 4],
            SolveConfig::default(),
            FluxBudget::default(),
        );
        assert!(matches!(result, Err(ConductionBoundError::Verification(
            TetError::Unsupported(reason))) if reason.contains("coating")));
    });
}
