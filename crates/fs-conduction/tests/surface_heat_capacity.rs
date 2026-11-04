//! Thin boundary storage preserves heating, reactions and source domains.
mod support;

use fs_conduction::fixtures::box_grid;
use fs_conduction::transient::backward_euler::{
    ArealHeatCapacity, BackwardEuler, NonlinearStepConfig, StepConfig, SurfaceHeatCapacity,
};
use fs_conduction::transient::{MaterialHeatCapacity, VolumetricHeatCapacity};
use fs_conduction::{
    ConductionError, ConductionMesh, ConductionProblem, ConductivityModel, ConductivityTable,
    LinearConfig, ScalarField, ThermalBc, ThermalBoundaryBuilder,
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
