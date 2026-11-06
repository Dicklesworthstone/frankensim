//! G0/G1/G3 frozen brittle spectra through the original solid and modal owners.
//! Manufactured coefficients verify mechanics, not measured material accuracy.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::brittle::BrittleError;
use fs_material::brittle::source::dynamics::resolve_brittle_dynamic_material;
use fs_material::brittle::source::*;
use fs_material::cohesive::CohesiveError;
use fs_material::cohesive::source::FRACTURE_ENERGY_PER_AREA_DIMS;
use fs_material::state_point::{
    DENSITY_PROPERTY, MaterialPropertySelection, POISSON_RATIO_PROPERTY, YOUNG_MODULUS_PROPERTY,
};
use fs_modal::eigh_gen_dense;
use fs_qty::{Density, Dims, Pressure};
use fs_solid::brittle::TetCrackBandGeometry;
use fs_solid::brittle::solid::dynamics::BrittleTetMotion;
use fs_solid::brittle::solid::dynamics::modes::{FrozenBrittleModesControls, FrozenBrittlePencil};
use fs_solid::brittle::solid::{BrittleTetControls, BrittleTetError, BrittleTetSolid};
use fs_solid::linear3::TetAssemblyBudget;
use fs_solid::loads3::TetLoadGeometry3;
use fs_solid::plastic3::IncrementSettings;

const NODES: [[f64; 3]; 4] = [[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const YOUNG: f64 = 1200.0;
const DENSITY: f64 = 6.0;
const VOLUME: f64 = 1.0 / 6.0;

fn with_cx<R>(run: impl FnOnce(&Cx<'_>, &CancelGate) -> R) -> R {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        run(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 584,
                    kernel_id: 584,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    })
}

fn budget() -> TetAssemblyBudget {
    TetAssemblyBudget {
        maximum_nodes: 4,
        maximum_tetrahedra: 1,
        maximum_free_dofs: 12,
        minimum_scaled_jacobian: 1e-10,
    }
}

fn controls() -> FrozenBrittleModesControls {
    FrozenBrittleModesControls {
        maximum_dofs: 12,
        maximum_staged_entries: 144,
    }
}

fn source(geometry: &TetLoadGeometry3<'_>) -> ResolvedBrittleCrackBand {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, YOUNG),
        (POISSON_RATIO_PROPERTY, Dims::NONE, 0.0),
        (DENSITY_PROPERTY, Density::DIMS, DENSITY),
        (BRITTLE_TENSILE_STRENGTH_PROPERTY, Pressure::DIMS, 1.2),
        (BRITTLE_SHEAR_STRENGTH_PROPERTY, Pressure::DIMS, 1.5),
        (
            BRITTLE_MODE_I_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            0.0012,
        ),
        (
            BRITTLE_MODE_II_FRACTURE_ENERGY_PROPERTY,
            FRACTURE_ENERGY_PER_AREA_DIMS,
            0.0024,
        ),
        (BRITTLE_BK_EXPONENT_PROPERTY, Dims::NONE, 1.5),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured frozen-brittle vibration fixture".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "manufactured brittle solid".into(),
            phase: "solid".into(),
            process: "numerical verification".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let band = TetCrackBandGeometry::from_element(geometry, 0, [1., 0., 0.]).unwrap();
    resolve_brittle_dynamic_material(
        &card,
        &QueryPoint::new(),
        MaterialPropertySelection::SingleClaimOnly,
        band.material_settings(0.02, 1e-9),
    )
    .unwrap()
}

fn deform<'a>(solid: &BrittleTetSolid<'a>, extension: f64, cx: &Cx<'_>) -> BrittleTetSolid<'a> {
    let prescribed: Vec<_> = (0..12)
        .map(|dof| (dof, if dof == 3 { extension } else { 0.0 }))
        .collect();
    solid
        .advance(
            &prescribed,
            &[],
            IncrementSettings::default(),
            BrittleTetControls {
                material_energy_tolerance_j: 1e-9,
                equilibrium_work_tolerance_j: 1e-9,
                maximum_endpoint_work_defect_j: 1.0,
            },
            cx,
        )
        .unwrap()
        .solid
}

fn pencil(solid: &BrittleTetSolid<'_>, fixed: &[usize], cx: &Cx<'_>) -> FrozenBrittlePencil {
    BrittleTetMotion::try_from_solid(solid, &[[0.; 3]; 4], 0., 10, cx)
        .unwrap()
        .frozen_vibration_pencil(fixed, controls(), cx)
        .unwrap()
}

fn eigenvalue(pencil: &FrozenBrittlePencil) -> f64 {
    assert_eq!(pencil.free_dofs(), &[3]);
    let modes = eigh_gen_dense(
        &pencil.stiffness_n_m().to_dense(),
        &pencil.mass_kg().to_dense(),
        1,
    )
    .unwrap();
    assert!(modes[0].residual < 1e-10);
    modes[0].lambda
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn g1_damage_lowers_native_modal_frequency_without_using_the_softening_tangent() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let materials = [source(&geometry)];
        let virgin = BrittleTetSolid::try_new(&geometry, &materials, budget(), cx).unwrap();
        let damaged = deform(&virgin, 0.003, cx);
        let fixed: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
        let intact = pencil(&virgin, &fixed, cx);
        let reduced = pencil(&damaged, &fixed, cx);
        let state = &damaged.solution().element_states[0];
        let d = state.damage();
        assert!(d > 0.0 && d < 1.0);
        assert!(
            materials[0]
                .law()
                .trial(state.strain(), &materials[0].law().initial_state())
                .unwrap()
                .branch_tangent_pa[0][0]
                < 0.0
        );
        // P1: K_33 = E V, M_33 = rho V / 10. Exact degree-two mass,
        // not row-sum lumping. The scalar generalized eigenvalue is 10 E/rho.
        near(intact.stiffness_n_m().get(0, 0), YOUNG * VOLUME, 1e-12);
        near(intact.mass_kg().get(0, 0), DENSITY * VOLUME / 10.0, 1e-16);
        near(eigenvalue(&intact), 10.0 * YOUNG / DENSITY, 1e-10);
        near(eigenvalue(&reduced), (1.0 - d) * eigenvalue(&intact), 1e-10);
        assert!(eigenvalue(&reduced) > 0.0);
        near(
            (eigenvalue(&reduced) / eigenvalue(&intact)).sqrt(),
            (1.0 - d).sqrt(),
            1e-14,
        );
        assert_eq!(reduced.element_damage(), &[d]);
        assert_eq!(intact.mass_kg().to_dense(), reduced.mass_kg().to_dense());
        assert_eq!(intact.reference_mass_kg(), reduced.reference_mass_kg());
    });
}

#[test]
fn g0_complete_failure_keeps_zero_transfer_modes_and_original_mass() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let materials = [source(&geometry)];
        let virgin = BrittleTetSolid::try_new(&geometry, &materials, budget(), cx).unwrap();
        let failed = deform(&virgin, 0.007, cx);
        assert!(failed.solution().element_states[0].failed());
        let fixed: Vec<_> = (0..12).filter(|&dof| dof != 3).collect();
        let lost = pencil(&failed, &fixed, cx);
        assert_eq!(lost.stiffness_n_m().get(0, 0), 0.0);
        assert_eq!(eigenvalue(&lost), 0.0);
        near(lost.reference_mass_kg(), DENSITY * VOLUME, 1e-15);
        near(lost.mass_kg().get(0, 0), DENSITY * VOLUME / 10.0, 1e-16);
        let compressed = deform(&failed, -0.001, cx);
        let compression = pencil(&compressed, &fixed, cx);
        near(compression.stiffness_n_m().get(0, 0), YOUNG * VOLUME, 1e-12);
        assert_eq!(compression.element_damage(), lost.element_damage());
        assert_eq!(
            compressed.fracture_dissipation_j(),
            failed.fracture_dissipation_j()
        );
    });
}

#[test]
fn g1_full_pencil_preserves_symmetry_rigid_motion_and_consistent_mass() {
    with_cx(|cx, _| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let materials = [source(&geometry)];
        let solid = BrittleTetSolid::try_new(&geometry, &materials, budget(), cx).unwrap();
        let p = pencil(&solid, &[], cx);
        assert_eq!(p.full_dofs(), 12);
        assert_eq!(p.free_dofs(), &(0..12).collect::<Vec<_>>());
        for i in 0..12 {
            for j in 0..12 {
                near(
                    p.stiffness_n_m().get(i, j),
                    p.stiffness_n_m().get(j, i),
                    1e-13,
                );
                let expected = if i % 3 != j % 3 {
                    0.0
                } else {
                    DENSITY * VOLUME / 20.0 * if i == j { 2.0 } else { 1.0 }
                };
                assert_eq!(p.mass_kg().get(i, j), expected);
            }
            for axis in 0..3 {
                let translation_force: f64 = (0..4)
                    .map(|node| p.stiffness_n_m().get(i, 3 * node + axis))
                    .sum();
                near(translation_force, 0.0, 1e-13);
            }
        }
        // Mixed affine strain, including tensor shear, independently contracted
        // through nu=0 isotropic storage. This exercises B^T W C B shear weights.
        let u = [
            0., 0., 0., 0.001, 0.0004, 0., 0.0002, -0.0003, 0.0006, 0., 0.0002, 0.0005,
        ];
        let work: f64 = (0..12)
            .map(|i| {
                (0..12)
                    .map(|j| 0.5 * u[i] * p.stiffness_n_m().get(i, j) * u[j])
                    .sum::<f64>()
            })
            .sum();
        let strain = [0.001_f64, -0.0003, 0.0005, 0.0003, 0.0004, 0.0];
        let exact = 0.5
            * YOUNG
            * VOLUME
            * (strain[..3].iter().map(|e| e * e).sum::<f64>()
                + 2.0 * strain[3..].iter().map(|e| e * e).sum::<f64>());
        near(work, exact, 1e-18);
        let fixed = [9, 1, 7];
        let restricted = pencil(&solid, &fixed, cx);
        for (i, &a) in restricted.free_dofs().iter().enumerate() {
            for (j, &b) in restricted.free_dofs().iter().enumerate() {
                assert_eq!(restricted.mass_kg().get(i, j), p.mass_kg().get(a, b));
                assert_eq!(
                    restricted.stiffness_n_m().get(i, j),
                    p.stiffness_n_m().get(a, b)
                );
            }
        }
    });
}

#[test]
fn g0_observation_and_refusals_preserve_accepted_history_and_motion() {
    with_cx(|cx, gate| {
        let geometry = TetLoadGeometry3::new(&NODES, &TETS, budget(), 4, cx).unwrap();
        let materials = [source(&geometry)];
        let virgin = BrittleTetSolid::try_new(&geometry, &materials, budget(), cx).unwrap();
        let damaged = deform(&virgin, 0.003, cx);
        let motion =
            BrittleTetMotion::try_from_solid(&damaged, &[[0.2, -0.1, 0.3]; 4], 2., 10, cx).unwrap();
        let old = motion.solid().solution().clone();
        let loss = motion.solid().fracture_dissipation_j();
        let kinetic = motion.kinetic_energy_j();
        let first = motion
            .frozen_vibration_pencil(&[9, 1, 7], controls(), cx)
            .unwrap();
        let repeat = motion
            .frozen_vibration_pencil(&[7, 9, 1], controls(), cx)
            .unwrap();
        assert_eq!(first.free_dofs(), repeat.free_dofs());
        assert_eq!(
            first.stiffness_n_m().to_dense(),
            repeat.stiffness_n_m().to_dense()
        );
        for fixed in [vec![1, 1], vec![12], (0..12).collect()] {
            assert!(
                motion
                    .frozen_vibration_pencil(&fixed, controls(), cx)
                    .is_err()
            );
        }
        for cap in [
            FrozenBrittleModesControls {
                maximum_dofs: 11,
                ..controls()
            },
            FrozenBrittleModesControls {
                maximum_staged_entries: 143,
                ..controls()
            },
        ] {
            assert!(motion.frozen_vibration_pencil(&[], cap, cx).is_err());
        }
        let closure = deform(&damaged, 0.0, cx);
        let at_closure =
            BrittleTetMotion::try_from_solid(&closure, &[[0.; 3]; 4], 2., 10, cx).unwrap();
        assert!(matches!(
            at_closure.frozen_vibration_pencil(&[], controls(), cx),
            Err(BrittleTetError::Material {
                element: 0,
                source: BrittleError::Cohesive(CohesiveError::NormalClosureSwitch { .. })
            })
        ));
        gate.request();
        assert_eq!(
            motion
                .frozen_vibration_pencil(&[], controls(), cx)
                .unwrap_err(),
            BrittleTetError::Cancelled
        );
        assert_eq!(motion.solid().solution(), &old);
        assert_eq!(motion.solid().fracture_dissipation_j(), loss);
        assert_eq!(motion.kinetic_energy_j(), kinetic);
        assert_eq!(motion.velocity_mps(), &[[0.2, -0.1, 0.3]; 4]);
        assert_eq!(motion.time_s(), 2.0);
        assert_eq!(motion.accepted_steps(), 0);
    });
}
