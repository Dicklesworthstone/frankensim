//! G1/G3: data-only bond substitution drives the existing structural solver.
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::cohesive::source::*;
use fs_material::state_point::{MaterialPropertySelection, interface_stiffness_per_area_dims};
use fs_qty::Pressure;
use fs_solid::cohesive::p1::nonmatching::{MaterialNonmatchingSurface, NonmatchingBudget};
use fs_solid::cohesive::p1::{MaterialFacePair, P1Interfaces, TraceBudget};
use fs_solid::cohesive::{Budget, CohesivePatch, CohesiveSystem, Config, Error, JumpTerm};
use fs_solver::{LinearOp, NewtonKrylovConfig};

fn source(strength: f64, fracture_energy: f64) -> ResolvedBilinearCohesive {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (
            COHESIVE_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            1e9,
        ),
        (COHESIVE_STRENGTH_PROPERTY, Pressure::DIMS, strength),
        (
            COHESIVE_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            fracture_energy,
        ),
        (
            COHESIVE_COMPRESSION_STIFFNESS_PROPERTY,
            interface_stiffness_per_area_dims(),
            2e9,
        ),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained().with("T", 290.0, 310.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic bond substitution; no measured claim".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let surface = |chemistry: &str| SurfaceSpec {
        material: MaterialStateId {
            chemistry: chemistry.into(),
            phase: "solid".into(),
            process: "manufactured bond fixture".into(),
            revision: 1,
        },
        texture_frame: format!("{chemistry} fixed reference frame"),
    };
    let card = InterfaceSystemCard::assemble(
        surface("A"),
        surface("B"),
        SystemContext {
            medium: "explicit cohesive bond".into(),
            third_body: None,
            environment: "fixed 300 K".into(),
            history: "initially intact".into(),
        },
        claims,
        vec![],
    )
    .unwrap();
    resolve_bilinear_cohesive(
        &card,
        &QueryPoint::new().with("T", 300.0).unwrap(),
        MaterialPropertySelection::SingleClaimOnly,
        0.02,
    )
    .unwrap()
}

struct Diagonal {
    n: usize,
    stiffness: f64,
}
impl LinearOp for Diagonal {
    fn n(&self) -> usize {
        self.n
    }
    fn apply(&self, x: &[f64], y: &mut [f64]) {
        for (out, value) in y.iter_mut().zip(x) {
            *out = self.stiffness * value;
        }
    }
}
fn budget() -> Budget {
    Budget {
        max_dofs: 24,
        max_patches: 6,
        max_jump_terms: 108,
        max_krylov_basis_entries: 50_000,
    }
}
fn config() -> Config {
    Config {
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-9,
            relative_tolerance: 1e-13,
            linear_restart: 8,
            forcing_maximum: 0.05,
            ..NewtonKrylovConfig::default()
        },
        max_newton_iterations: 80,
        max_operator_evaluations: 2000,
        force_tolerance_n: 1e-7,
        local_energy_tolerance_j_m2: 1e-7,
        total_cohesive_energy_tolerance_j: 1e-8,
        maximum_displacement_m: 0.02,
    }
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn changing_only_interface_data_changes_load_path_and_fracture_energy() {
    let elastic = Diagonal {
        n: 1,
        stiffness: 5e6,
    };
    let area = 1e-4;
    let mut reactions = Vec::new();
    for (strength, fracture_energy) in [(1e6, 2000.0), (1.5e6, 4000.0)] {
        let source = source(strength, fracture_energy);
        let surfaces = [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ];
        let patch = CohesivePatch::from_interface(
            7,
            area,
            &source,
            surfaces,
            vec![JumpTerm {
                dof: 0,
                local: [1.0, 0.0, 0.0],
            }],
        )
        .unwrap();
        let system = CohesiveSystem::new(&elastic, vec![patch], budget()).unwrap();
        assert_eq!(
            system.patches()[0].source().unwrap().resolved(),
            source.resolved()
        );
        let mut state = system.bonded_state();
        let mut work = 0.0;
        for actuator in [0.0004, 0.0012, 0.0022, 0.0032, 0.006] {
            let solved = state
                .solve_controlled(
                    &[elastic.stiffness * actuator],
                    &[[0.0; 3]],
                    config(),
                    &mut || false,
                )
                .unwrap();
            // Independent equilibrium: support force K(D-g) equals A*t(g).
            let onset = strength / 1e9;
            let failure = 2.0 * fracture_energy / strength;
            let elastic_jump = elastic.stiffness * actuator / (elastic.stiffness + area * 1e9);
            let softening = strength / (failure - onset);
            let expected = if elastic_jump <= onset {
                elastic_jump
            } else if actuator >= failure {
                actuator
            } else {
                (elastic.stiffness * actuator - area * softening * failure)
                    / (elastic.stiffness - area * softening)
            };
            close(solved.state.displacements_m()[0], expected, 1e-12);
            if actuator == 0.0022 {
                reactions.push(elastic.stiffness * (actuator - solved.state.displacements_m()[0]));
            }
            work += solved.cohesive_energy_j.work;
            state = solved.state;
        }
        assert!(state.points()[0].failed());
        close(work, area * fracture_energy, 1e-9);
        let before = state.points().to_vec();
        assert!(matches!(
            state.solve_controlled(&[0.0], &[[0.0; 3]], config(), &mut || true),
            Err(Error::Cancelled)
        ));
        assert_eq!(state.points(), before);
        assert_eq!(
            system.patches()[0].source().unwrap().resolved(),
            source.resolved()
        );
        assert!(
            CohesivePatch::from_interface(8, area, &source, [surfaces[1], surfaces[0]], vec![])
                .is_err()
        );
    }
    assert!(reactions[1] > reactions[0] * 1.5);
}

#[test]
fn material_faces_use_native_quadrature_and_preserve_source_order() {
    let source = source(1e6, 2000.0);
    let nodes = [
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., -1.],
        [0., 0., 0.],
        [1., 0., 0.],
        [0., 1., 0.],
        [0., 0., 1.],
    ];
    let tetrahedra = [[0, 1, 2, 3], [4, 5, 6, 7]];
    let free = [9, 10, 11, 21, 22, 23];
    let pair = MaterialFacePair {
        id: 41,
        side_a: [2, 0, 1],
        side_b: [5, 6, 4],
        source: &source,
        surface_materials: [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
    };
    let trace_budget = TraceBudget {
        maximum_nodes: 8,
        maximum_tetrahedra: 2,
        maximum_interfaces: 1,
        minimum_face_quality: 1e-10,
    };
    let faces = P1Interfaces::from_interfaces_controlled(
        &nodes,
        &tetrahedra,
        &free,
        &[pair],
        trace_budget,
        &mut || false,
    )
    .unwrap();
    for patch in faces.patches() {
        assert_eq!(patch.source().unwrap().resolved(), source.resolved());
    }
    close(
        faces.points().iter().map(|p| p.area_weight_m2).sum(),
        0.5,
        1e-15,
    );
    let elastic = Diagonal {
        n: 6,
        stiffness: 1.0,
    };
    let system = CohesiveSystem::new(&elastic, faces.patches().to_vec(), budget()).unwrap();
    let mut prescribed = [0.0; 24];
    for node in [4, 5, 6] {
        prescribed[3 * node + 2] = 0.006;
    }
    let jumps = faces.prescribed_jumps(&prescribed, &mut || false).unwrap();
    let state = system.bonded_state();
    let solved = state
        .solve_controlled(&[0.0; 6], &jumps, config(), &mut || false)
        .unwrap();
    assert!(solved.state.points().iter().all(|p| p.failed()));
    close(solved.cohesive_energy_j.dissipated, 0.5 * 2000.0, 1e-9);
    let reversed = MaterialFacePair {
        surface_materials: [pair.surface_materials[1], pair.surface_materials[0]],
        ..pair
    };
    assert!(matches!(
        P1Interfaces::from_interfaces_controlled(
            &nodes,
            &tetrahedra,
            &free,
            &[reversed],
            trace_budget,
            &mut || false
        ),
        Err(Error::Source { .. })
    ));
    assert!(matches!(
        P1Interfaces::from_interfaces_controlled(
            &nodes,
            &tetrahedra,
            &free,
            &[pair],
            trace_budget,
            &mut || true
        ),
        Err(Error::Cancelled)
    ));
}

#[test]
fn nonmatching_overlap_points_retain_original_source_and_require_surface_order() {
    let source = source(1e6, 2000.0);
    let nodes = [
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [0.5, 0.5, -1.],
        [0., 0., 0.],
        [1., 0., 0.],
        [1., 1., 0.],
        [0., 1., 0.],
        [0.5, 0.5, 1.],
    ];
    let tetrahedra = [[0, 1, 2, 4], [0, 2, 3, 4], [5, 6, 8, 9], [6, 7, 8, 9]];
    let side_a = [[0, 1, 2], [0, 2, 3]];
    let side_b = [[5, 6, 8], [6, 7, 8]];
    let surface = MaterialNonmatchingSurface {
        id: 41,
        side_a: &side_a,
        side_b: &side_b,
        source: &source,
        surface_materials: [
            &source.surfaces()[0].material,
            &source.surfaces()[1].material,
        ],
    };
    let controls = NonmatchingBudget {
        mesh: TraceBudget {
            maximum_nodes: 10,
            maximum_tetrahedra: 4,
            maximum_interfaces: 1,
            minimum_face_quality: 1e-10,
        },
        maximum_surface_faces: 4,
        maximum_pair_checks: 6,
        maximum_quadrature_points: 12,
        relative_geometry_tolerance: 1e-11,
    };
    let faces = P1Interfaces::from_nonmatching_interfaces_controlled(
        &nodes,
        &tetrahedra,
        &[],
        &[surface],
        controls,
        &mut || false,
    )
    .unwrap();
    assert_eq!(faces.points().len(), 12);
    for patch in faces.patches() {
        assert!(core::ptr::eq(patch.source().unwrap(), &source));
    }
    close(
        faces.points().iter().map(|p| p.area_weight_m2).sum(),
        1.0,
        1e-15,
    );
    let reversed = MaterialNonmatchingSurface {
        surface_materials: [surface.surface_materials[1], surface.surface_materials[0]],
        ..surface
    };
    assert!(matches!(
        P1Interfaces::from_nonmatching_interfaces_controlled(
            &nodes,
            &tetrahedra,
            &[],
            &[reversed],
            controls,
            &mut || false
        ),
        Err(Error::Source { .. })
    ));
}
