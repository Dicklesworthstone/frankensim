//! Cooling-induced crystal growth reaches the existing native continuum solver.
//! Authored caloric/kinetic/elastic data; no experimental material claim.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_blake3::ContentHash;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::phase::crystallization::CrystallizationCaloricPair;
use fs_material::phase::crystallization::mechanics::{
    CrystallizationLock, CrystallizationMechanicsError, CrystallizationStrainLimits,
};
use fs_material::phase::{EquilibriumEnthalpyPhaseCurve, HeatCapacityKnot, SolidLiquidPhase};
use fs_material::plastic::reduction::ReductionConfig;
use fs_material::{IsotropicElastic, SmallStrainLaw};
use fs_solid::{
    Formulation, IncrementSettings, Mesh2, PlaneAssemblyBudget, PlaneKind, SmallStrainLawField,
    SmallStrainPlaneProblem,
};
use fs_thermal::crystallization::CrystallizationSolve;

fn branch(phase: SolidLiquidPhase, cp: f64, datum: f64, rho: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_from_heat_capacity(
        ContentHash([0x71; 32]),
        phase,
        datum,
        &[280.0, 380.0].map(|temperature_k| HeatCapacityKnot {
            temperature_k,
            specific_heat_capacity_j_kg_k: cp,
            bulk_density_kg_m3: rho,
        }),
        1e-6,
    )
    .unwrap()
}
fn limits() -> CrystallizationStrainLimits {
    CrystallizationStrainLimits {
        maximum_absolute_total_mandel_strain: 0.02,
        maximum_absolute_elastic_mandel_strain: 0.02,
        maximum_absolute_free_isotropic_strain: 0.01,
    }
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 71,
                kernel_id: 71,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        ));
    });
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual} != {expected}"
    );
}

#[test]
fn g1_cooling_crystallization_drives_native_free_shrinkage_and_constrained_tension() {
    let parent = branch(SolidLiquidPhase::Liquid, 2000.0, 100_000.0, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1500.0, 0.0, 1030.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let rate = |_: f64| Ok(0.2);
    let kinetic = pair
        .law(&rate, ContentHash([0x72; 32]), 2.0, 0.6, 400.0)
        .unwrap();
    let initial = kinetic.state_at_temperature(350.0, 0.5).unwrap();
    let lock = CrystallizationLock::new(
        &pair,
        initial.temperature_k(),
        initial.crystalline_mass_fraction(),
        [0.0; 6],
    )
    .unwrap();
    let cooled = kinetic
        .advance_specific_enthalpy(
            initial,
            initial.specific_enthalpy_j_kg() - 50_000.0,
            2.0,
            CrystallizationSolve::default(),
            || false,
        )
        .unwrap()
        .state;
    assert!(cooled.temperature_k() < initial.temperature_k());
    assert!(cooled.crystalline_mass_fraction() > initial.crystalline_mass_fraction());
    let elastic = IsotropicElastic::new(2e9, 0.3, 0.02).unwrap();
    let law = lock
        .freeze(
            &pair,
            cooled.temperature_k(),
            cooled.crystalline_mass_fraction(),
            &elastic,
            limits(),
        )
        .unwrap();
    assert!(std::ptr::eq(law.source(), &elastic));
    let volume = |x: f64| (1.0 - x) / 1000.0 + x / 1030.0;
    let expected = (volume(cooled.crystalline_mass_fraction())
        / volume(initial.crystalline_mass_fraction()))
    .cbrt()
        - 1.0;
    close(law.free_isotropic_strain(), expected, 2e-16);
    assert!(expected < 0.0);
    with_cx(|cx| {
        for mesh in [
            Mesh2::triangles(2.0, 1.0, 2, 2),
            Mesh2::quads(2.0, 1.0, 2, 2),
        ] {
            let right = mesh.nodes.iter().position(|&p| p == [2.0, 0.0]).unwrap();
            for constrained in [false, true] {
                let prescribed: Vec<_> = if constrained {
                    (0..2 * mesh.nodes.len()).map(|d| (d, 0.0)).collect()
                } else {
                    vec![(0, 0.0), (1, 0.0), (2 * right + 1, 0.0)]
                };
                let problem = SmallStrainPlaneProblem {
                    mesh: &mesh,
                    laws: SmallStrainLawField::Uniform(&law),
                    plane: if constrained {
                        PlaneKind::Strain
                    } else {
                        PlaneKind::Stress
                    },
                    formulation: Formulation::Standard,
                    thickness_m: 0.2,
                    reduction: ReductionConfig {
                        max_iterations: 8,
                        stress_tolerance_pa: 1e-6,
                        strain_limit: 0.02,
                    },
                    rotation_limit_rad: 0.02,
                    prescribed_m: &prescribed,
                    nodal_forces_n: &[],
                    load_path: &[1.0],
                    budget: PlaneAssemblyBudget::default(),
                    settings: IncrementSettings::default(),
                };
                let solved = problem.solve(cx).unwrap();
                for (point, displacement) in mesh.nodes.iter().zip(&solved.displacement_m) {
                    for i in 0..2 {
                        close(
                            displacement[i],
                            if constrained {
                                0.0
                            } else {
                                expected * point[i]
                            },
                            2e-10,
                        );
                    }
                }
                let tension = if constrained {
                    -elastic.youngs * expected / (1.0 - 2.0 * elastic.poisson)
                } else {
                    0.0
                };
                for stress in &solved.quadrature_stress_pa {
                    for (i, component) in stress.iter().enumerate() {
                        close(*component, if i < 3 { tension } else { 0.0 }, 0.05);
                    }
                }
            }
        }
    });
}

#[test]
fn g0_lock_retains_reference_shear_and_refuses_unmodelled_reset_and_strain() {
    let parent = branch(SolidLiquidPhase::Solid, 2000.0, 100_000.0, 1000.0);
    let crystal = branch(SolidLiquidPhase::Solid, 1500.0, 0.0, 1030.0);
    let pair = CrystallizationCaloricPair::new(&parent, &crystal, [280.0, 380.0]).unwrap();
    let total = [0.002, -0.001, 0.0005, 0.001, -0.0002, 0.0003];
    let lock = CrystallizationLock::new(&pair, 350.0, 0.2, total).unwrap();
    let elastic = IsotropicElastic::new(2e9, 0.3, 0.02).unwrap();
    let law = lock.freeze(&pair, 350.0, 0.2, &elastic, limits()).unwrap();
    assert_eq!(*law.free_strain(), total);
    assert_eq!(law.stress(&total, &()), [0.0; 6]);
    assert_eq!(law.tangent(&total, &()), elastic.tangent(&[0.0; 6], &()));
    assert!(lock.freeze(&pair, 350.0, 0.1, &elastic, limits()).is_err());
    assert!(lock.freeze(&pair, 500.0, 0.3, &elastic, limits()).is_err());
    assert!(law.stress(&[0.03; 6], &()).iter().all(|s| s.is_nan()));
    let tight = lock
        .freeze(
            &pair,
            350.0,
            0.2,
            &elastic,
            CrystallizationStrainLimits {
                maximum_absolute_elastic_mandel_strain: 0.001,
                ..limits()
            },
        )
        .unwrap();
    assert!(tight.elastic_strain(&[0.0; 6]).is_err());
    let other_crystal = branch(SolidLiquidPhase::Solid, 1500.0, 0.0, 1040.0);
    let other = CrystallizationCaloricPair::new(&parent, &other_crystal, [280.0, 380.0]).unwrap();
    assert!(matches!(
        lock.freeze(&other, 350.0, 0.3, &elastic, limits()),
        Err(CrystallizationMechanicsError::PairMismatch)
    ));
}
