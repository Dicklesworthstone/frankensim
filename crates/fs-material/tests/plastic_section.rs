//! G0/G1/G3 native through-thickness histories, work conjugacy and springback.
use fs_material::SmallStrainLaw;
use fs_material::elastic::IsotropicElastic;
use fs_material::plastic::J2Plasticity;
use fs_material::plastic::chaboche::{BackstressBranch, ChabocheJ2};
use fs_material::plastic::rate::LinearOverstressJ2;
use fs_material::plastic::reduction::{PlaneStressWorkspace, ReductionConfig};
use fs_material::plastic::section::{
    LayeredPlaneStressSection, PlaneStressLayer, SectionBudget, SectionError, ThicknessQuadrature,
};
use fs_material::plastic::voce::VoceJ2;

fn controls() -> ReductionConfig {
    ReductionConfig {
        max_iterations: 48,
        stress_tolerance_pa: 1e-11,
        strain_limit: 0.05,
    }
}

fn budget() -> SectionBudget {
    SectionBudget {
        maximum_layers: 8,
        maximum_quadrature_points: 128,
    }
}

fn near(actual: f64, expected: f64, allowance: f64) {
    assert!(
        (actual - expected).abs() <= allowance,
        "{actual:e} != {expected:e}; allowance {allowance:e}"
    );
}

fn elastic(youngs: f64, poisson: f64) -> IsotropicElastic {
    IsotropicElastic::new(youngs, poisson, 0.05).unwrap()
}

fn plastic() -> J2Plasticity {
    J2Plasticity::new_mixed(elastic(1200.0, 0.2), 3.0, 40.0, 160.0).unwrap()
}

fn section<L: SmallStrainLaw>(law: L, subdivisions: usize) -> LayeredPlaneStressSection<L> {
    LayeredPlaneStressSection::new(
        vec![PlaneStressLayer {
            bottom_z_m: -0.1,
            top_z_m: 0.1,
            subdivisions,
            law,
        }],
        ThicknessQuadrature::Gauss3,
        controls(),
        budget(),
    )
    .unwrap()
}

fn plane_stiffness(youngs: f64, poisson: f64) -> [[f64; 3]; 3] {
    let diagonal = youngs / (1.0 - poisson * poisson);
    [
        [diagonal, poisson * diagonal, 0.0],
        [poisson * diagonal, diagonal, 0.0],
        [0.0, 0.0, youngs / (2.0 * (1.0 + poisson))],
    ]
}

#[test]
fn g1_elastic_layers_recover_exact_abd_and_engineering_shear_work() {
    for rule in [ThicknessQuadrature::Gauss2, ThicknessQuadrature::Gauss3] {
        let layers = [(-0.1_f64, 0.02_f64, 1200.0, 0.2), (0.02, 0.12, 2300.0, 0.3)];
        let section = LayeredPlaneStressSection::new(
            layers
                .iter()
                .enumerate()
                .map(|(index, &(bottom, top, e, nu))| PlaneStressLayer {
                    bottom_z_m: bottom,
                    top_z_m: top,
                    subdivisions: 1 + 4 * index,
                    law: elastic(e, nu),
                })
                .collect(),
            rule,
            controls(),
            budget(),
        )
        .unwrap();
        let mut exact = [[0.0; 6]; 6];
        for &(bottom, top, e, nu) in &layers {
            let moments = [
                top - bottom,
                0.5 * (top * top - bottom * bottom),
                (top * top * top - bottom * bottom * bottom) / 3.0,
            ];
            let q = plane_stiffness(e, nu);
            for i in 0..6 {
                for j in 0..6 {
                    exact[i][j] += moments[i / 3 + j / 3] * q[i % 3][j % 3];
                }
            }
        }
        let strain = [0.0002, -0.0001, 0.0003, 0.005, -0.004, 0.003];
        let old = section.initial_state();
        let response = section.evaluate(strain, &old).unwrap();
        assert_eq!(old, section.initial_state());
        assert_eq!(section.point_count(), 6 * rule.points_per_subdivision());
        assert!(section.tangent_symmetric());
        assert!(
            response.tangent[0][3] > 0.0,
            "eccentric unequal layers have B coupling"
        );
        let mut twice_energy = 0.0;
        for i in 0..6 {
            let force: f64 = (0..6).map(|j| exact[i][j] * strain[j]).sum();
            near(response.resultants[i], force, 3e-14);
            twice_energy += force * strain[i];
            for j in 0..6 {
                near(response.tangent[i][j], exact[i][j], 2e-12);
                near(response.tangent[i][j], response.tangent[j][i], 2e-12);
            }
        }
        // Independent 3-D elastic energy with tensor contraction, including
        // transverse contraction and both tensor shear entries.
        let mut integrated_energy = 0.0;
        for (point, state) in section.points().iter().zip(&response.state.points) {
            let law = &section.layers()[point.layer].law;
            let (lambda, mu) = law.lame();
            let trace = state.strain[0] + state.strain[1] + state.strain[2];
            let norm_squared: f64 = state
                .strain
                .iter()
                .enumerate()
                .map(|(i, value)| value * value * if i < 3 { 1.0 } else { 2.0 })
                .sum();
            integrated_energy +=
                point.weight_m * (mu * norm_squared + 0.5 * lambda * trace * trace);
        }
        near(integrated_energy, 0.5 * twice_energy, 2e-17);
        section.validate_state(&response.state).unwrap();
    }
}

#[test]
fn g3_reference_surface_shift_preserves_strains_forces_and_virtual_work() {
    let original = section(elastic(1200.0, 0.2), 3);
    let shift = 0.037;
    let moved = LayeredPlaneStressSection::new(
        vec![PlaneStressLayer {
            bottom_z_m: -0.1 - shift,
            top_z_m: 0.1 - shift,
            subdivisions: 3,
            law: elastic(1200.0, 0.2),
        }],
        ThicknessQuadrature::Gauss3,
        controls(),
        budget(),
    )
    .unwrap();
    let a = [0.001, -0.0007, 0.0018, 0.012, -0.005, 0.014];
    let mut b = a;
    for i in 0..3 {
        b[i] += shift * a[i + 3];
    }
    let first = original.evaluate(a, &original.initial_state()).unwrap();
    let second = moved.evaluate(b, &moved.initial_state()).unwrap();
    for i in 0..3 {
        near(second.resultants[i], first.resultants[i], 1e-14);
        near(
            second.resultants[i + 3],
            first.resultants[i + 3] - shift * first.resultants[i],
            2e-15,
        );
    }
    for (left, right) in first.state.points.iter().zip(&second.state.points) {
        for i in 0..6 {
            near(left.strain[i], right.strain[i], 2e-18);
        }
    }
    let work_a: f64 = first
        .resultants
        .iter()
        .zip(a)
        .map(|(force, strain)| force * strain)
        .sum();
    let work_b: f64 = second
        .resultants
        .iter()
        .zip(b)
        .map(|(force, strain)| force * strain)
        .sum();
    near(work_a, work_b, 1e-17);
}

#[test]
fn g1_partial_thickness_yield_and_elastic_unloading_leave_permanent_curvature() {
    let section = section(plastic(), 6);
    let virgin = section.initial_state();
    let elastic = section.evaluate([0.0; 6], &virgin).unwrap();
    let peak = [0.0, 0.0, 0.0, 0.048, 0.0, 0.0];
    let loaded = section.evaluate(peak, &virgin).unwrap();
    assert_eq!(virgin, section.initial_state());
    assert!(loaded.state.points.iter().any(|p| p.history.alpha == 0.0));
    assert!(loaded.state.points.iter().any(|p| p.history.alpha > 1e-4));
    for i in 0..3 {
        near(loaded.resultants[i], 0.0, 2e-13);
    }
    assert!(loaded.tangent[3][3] < elastic.tangent[3][3]);
    for (low, high) in loaded
        .state
        .points
        .iter()
        .zip(loaded.state.points.iter().rev())
    {
        near(low.history.alpha, high.history.alpha, 2e-16);
        for i in 0..6 {
            near(
                low.history.plastic_strain[i],
                -high.history.plastic_strain[i],
                2e-16,
            );
        }
    }
    // Independent elastic springback: invert isotropic plane-stress D, with
    // no further flow. The remaining plastic strain varies through thickness.
    let ei = 1200.0 * 0.2_f64.powi(3) / 12.0;
    let mut residual = peak;
    residual[3] -= (loaded.resultants[3] - 0.2 * loaded.resultants[4]) / ei;
    residual[4] -= (loaded.resultants[4] - 0.2 * loaded.resultants[3]) / ei;
    residual[5] -= loaded.resultants[5] / (ei / 2.4);
    assert!(residual[3] > 0.0 && residual[3] < peak[3]);
    let unloaded = section.evaluate(residual, &loaded.state).unwrap();
    for force in unloaded.resultants {
        near(force, 0.0, 2e-12);
    }
    for (before, after) in loaded.state.points.iter().zip(&unloaded.state.points) {
        assert_eq!(
            before.history, after.history,
            "springback must not reset plastic histories"
        );
    }
    assert_eq!(loaded.state.generalized_strain, peak);
}

#[test]
fn g1_nonproportional_chaboche_section_keeps_nonsymmetric_consistent_tangent() {
    let law = ChabocheJ2::new(
        VoceJ2::new(plastic(), 5.0, 180.0).unwrap(),
        [
            BackstressBranch::new(800.0, 150.0).unwrap(),
            BackstressBranch::new(250.0, 20.0).unwrap(),
        ],
    )
    .unwrap();
    let section = section(law, 2);
    assert!(!section.tangent_symmetric());
    let loaded = section
        .evaluate(
            [0.014, -0.002, 0.005, 0.015, -0.01, 0.01],
            &section.initial_state(),
        )
        .unwrap();
    let frozen = loaded.state.clone();
    let strain = [-0.004, 0.003, -0.01, -0.02, 0.04, 0.02];
    let mut workspace = PlaneStressWorkspace::new().unwrap();
    let next = section
        .evaluate_with_workspace(strain, &frozen, &mut workspace)
        .unwrap();
    assert!(
        next.state.points.iter().all(|p| p
            .history
            .backstress_pa
            .iter()
            .flatten()
            .any(|v| *v != 0.0))
    );
    let mut asymmetry = 0.0_f64;
    for i in 0..6 {
        for j in 0..6 {
            asymmetry = asymmetry.max((next.tangent[i][j] - next.tangent[j][i]).abs());
        }
    }
    assert!(
        asymmetry > 1e-4,
        "nonproportional recovery is not spuriously symmetrized"
    );
    for column in 0..6 {
        let h = if column < 3 { 1e-7 } else { 1e-6 };
        let mut plus = strain;
        let mut minus = strain;
        plus[column] += h;
        minus[column] -= h;
        let right = section
            .evaluate_with_workspace(plus, &frozen, &mut workspace)
            .unwrap();
        let left = section
            .evaluate_with_workspace(minus, &frozen, &mut workspace)
            .unwrap();
        for row in 0..6 {
            near(
                (right.resultants[row] - left.resultants[row]) / (2.0 * h),
                next.tangent[row][column],
                4e-5,
            );
        }
    }
    assert_eq!(
        loaded.state, frozen,
        "all numerical trials share frozen accepted history"
    );
}

#[test]
fn g0_rate_dependent_section_advances_one_interval_only_when_candidate_is_accepted() {
    let section = section(LinearOverstressJ2::new(plastic(), 40.0, 0.1).unwrap(), 2);
    let old = section.initial_state();
    let strain = [0.012, -0.002, 0.004, 0.02, -0.01, 0.01];
    let first = section.evaluate(strain, &old).unwrap();
    assert_eq!(first, section.evaluate(strain, &old).unwrap());
    assert_eq!(old, section.initial_state());
    let relaxed = section.evaluate(strain, &first.state).unwrap();
    assert!(relaxed.resultants[0] < first.resultants[0]);
    for (before, after) in first.state.points.iter().zip(&relaxed.state.points) {
        assert!(after.history.alpha > before.history.alpha);
    }
}

#[test]
fn g0_invalid_geometry_budget_and_stale_history_refuse_atomically() {
    let layer = PlaneStressLayer {
        bottom_z_m: -0.1,
        top_z_m: 0.1,
        subdivisions: 2,
        law: plastic(),
    };
    let tiny_budget = SectionBudget {
        maximum_layers: 1,
        maximum_quadrature_points: 3,
    };
    assert!(matches!(
        LayeredPlaneStressSection::new(
            vec![layer.clone()],
            ThicknessQuadrature::Gauss2,
            controls(),
            tiny_budget
        ),
        Err(SectionError::Budget { .. })
    ));
    let mut gap = layer.clone();
    gap.bottom_z_m = 0.1001;
    gap.top_z_m = 0.2;
    assert!(
        LayeredPlaneStressSection::new(
            vec![layer.clone(), gap],
            ThicknessQuadrature::Gauss2,
            controls(),
            budget()
        )
        .is_err()
    );
    for subdivisions in [0, usize::MAX] {
        let mut invalid = layer.clone();
        invalid.subdivisions = subdivisions;
        assert!(
            LayeredPlaneStressSection::new(
                vec![invalid],
                ThicknessQuadrature::Gauss3,
                controls(),
                budget()
            )
            .is_err()
        );
    }
    let section = section(plastic(), 2);
    let accepted = section
        .evaluate(
            [0.008, -0.002, 0.003, 0.01, 0.0, 0.0],
            &section.initial_state(),
        )
        .unwrap()
        .state;
    let saved = accepted.clone();
    for strain in [[f64::NAN; 6], [f64::MAX; 6], [0.0, 0.0, 0.0, 1.0, 0.0, 0.0]] {
        assert!(section.evaluate(strain, &accepted).is_err());
        assert_eq!(accepted, saved);
    }
    let mut stale = accepted.clone();
    stale.generalized_strain[0] += 1e-4;
    assert!(section.evaluate([0.0; 6], &stale).is_err());
    let mut truncated = accepted.clone();
    truncated.points.pop();
    assert!(section.evaluate([0.0; 6], &truncated).is_err());
}

#[path = "plastic_section/energy.rs"]
mod energy;
#[path = "plastic_section/thermal.rs"]
mod thermal;
#[path = "plastic_section/creep.rs"]
mod creep;
#[path = "plastic_section/hill.rs"]
mod hill;
