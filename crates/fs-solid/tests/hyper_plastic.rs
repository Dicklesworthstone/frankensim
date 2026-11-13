//! Native finite tetrahedra: equilibrium, retained plastic history and work.
//! Synthetic coefficients and homogeneous patches verify mechanics, not a
//! material calibration or bending accuracy of a locking-prone P1 element.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Limits, FiniteJ2Parameters};
use fs_solid::hyper3::HyperTetSettings;
use fs_solid::hyper3::plastic::{
    FinitePlasticLawField, FinitePlasticStepControls, FinitePlasticTetError, FinitePlasticTetSolid,
};
use fs_solid::linear3::TetAssemblyBudget;

const MU: f64 = 80.0e6;
const KAPPA: f64 = 160.0e6;
const YIELD: f64 = 2.0e6;
const HARDENING: f64 = 8.0e6;
const E: f64 = 9.0 * KAPPA * MU / (3.0 * KAPPA + MU);
const NU: f64 = (3.0 * KAPPA - 2.0 * MU) / (2.0 * (3.0 * KAPPA + MU));

fn with_cx<R>(f: impl FnOnce(&Cx<'_>) -> R) -> R {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            StreamKey {
                seed: 3,
                kernel_id: 17,
                tile: 0,
                iteration: 0,
            },
            Budget::INFINITE,
            ExecMode::Deterministic,
        );
        f(&cx)
    })
}

fn law() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: KAPPA,
            yield_stress_pa: YIELD,
            hardening_modulus_pa: HARDENING,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.25,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}

fn controls() -> FinitePlasticStepControls {
    FinitePlasticStepControls {
        equilibrium_work_tolerance_j: 1.0e-4,
        maximum_endpoint_work_defect_j: 20_000.0,
        yield_tolerance_pa: 1.0e-4,
    }
}

// Twelve conforming tetrahedra connect a triangulated unit-cube boundary to
// an interior vertex. Its three free coordinates participate in every patch.
fn cube() -> (Vec<[f64; 3]>, Vec<[usize; 4]>) {
    let mut nodes: Vec<_> = (0..8)
        .map(|node| {
            [
                (node & 1) as f64,
                ((node >> 1) & 1) as f64,
                ((node >> 2) & 1) as f64,
            ]
        })
        .collect();
    nodes.push([0.5; 3]);
    let mut tets = Vec::new();
    for [a, b, c, d] in [
        [0, 2, 3, 1],
        [4, 5, 7, 6],
        [0, 1, 5, 4],
        [2, 6, 7, 3],
        [0, 4, 6, 2],
        [1, 3, 7, 5],
    ] {
        for triangle in [[a, b, c], [a, c, d]] {
            let mut tet = [8, triangle[0], triangle[1], triangle[2]];
            let edge = |node: usize| [0, 1, 2].map(|i| nodes[node][i] - 0.5);
            let (u, v, w) = (edge(tet[1]), edge(tet[2]), edge(tet[3]));
            let det = u[0] * (v[1] * w[2] - v[2] * w[1]) - u[1] * (v[0] * w[2] - v[2] * w[0])
                + u[2] * (v[0] * w[1] - v[1] * w[0]);
            if det < 0.0 {
                tet.swap(2, 3);
            }
            tets.push(tet);
        }
    }
    (nodes, tets)
}

fn rollers(nodes: &[[f64; 3]], axial_stretch: Option<f64>) -> Vec<(usize, f64)> {
    let mut prescribed = Vec::new();
    for (node, position) in nodes.iter().enumerate() {
        if position[0] == 0.0 {
            prescribed.push((3 * node, 0.0));
        }
        if position[0] == 1.0
            && let Some(stretch) = axial_stretch
        {
            prescribed.push((3 * node, stretch - 1.0));
        }
        if position[1] == 0.0 {
            prescribed.push((3 * node + 1, 0.0));
        }
        if position[2] == 0.0 {
            prescribed.push((3 * node + 2, 0.0));
        }
    }
    prescribed
}

fn affine(nodes: &[[f64; 3]], f: [f64; 9], interior_fixed: bool) -> Vec<(usize, f64)> {
    nodes
        .iter()
        .enumerate()
        .filter(|(node, _)| interior_fixed || *node != 8)
        .flat_map(|(node, position)| {
            (0..3).map(move |i| {
                let current = (0..3).map(|j| f[3 * i + j] * position[j]).sum::<f64>();
                (3 * node + i, current - position[i])
            })
        })
        .collect()
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:e}"
    );
}

// Independent coaxial finite uniaxial result. Plastic axial logarithmic stretch
// p need not equal accumulated alpha after reversal. Kirchhoff tau remains the
// energetic elastic stress, while nominal Pxx=tau/lambda.
fn check_uniaxial(solid: &FinitePlasticTetSolid<'_>, a: f64, p: f64, alpha: f64) {
    check_uniaxial_with_hardening(solid, a, p, alpha, 0.5 * HARDENING * alpha * alpha);
}

fn check_uniaxial_with_hardening(
    solid: &FinitePlasticTetSolid<'_>,
    a: f64,
    p: f64,
    alpha: f64,
    hardening_storage: f64,
) {
    let tau = E * (a - p);
    let axial = a.exp();
    let lateral = (-NU * tau / E - 0.5 * p).exp();
    let expected_f = [axial, 0.0, 0.0, 0.0, lateral, 0.0, 0.0, 0.0, lateral];
    for (f, state) in solid
        .element_deformation()
        .iter()
        .zip(solid.element_states())
    {
        for i in 0..9 {
            close(f[i], expected_f[i], 2.0e-8);
        }
        close(state.alpha, alpha, 3.0e-8);
        for (i, expected) in [(-2.0 * p).exp(), p.exp(), p.exp(), 0.0, 0.0, 0.0]
            .iter()
            .enumerate()
        {
            close(state.inverse_plastic_cauchy_green[i], *expected, 1.0e-7);
        }
    }
    for stress in solid.element_first_piola_pa() {
        close(stress[0], tau / axial, 0.15);
        for &component in &[1, 2, 3, 4, 5, 6, 7, 8] {
            close(stress[component], 0.0, 0.15);
        }
    }
    // Unit reference volume; finite elastic Hencky energy plus stored hardening.
    close(
        solid.stored_energy_j(),
        0.5 * tau * tau / E + hardening_storage,
        0.01,
    );
}

#[test]
fn tabulated_compression_solves_free_sides_and_retains_permanent_set_after_release() {
    use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};
    with_cx(|cx| {
        let (nodes, tets) = cube();
        let mut parameters = *law().parameters();
        parameters.hardening_modulus_pa = 0.0;
        let mut limits = *law().limits();
        limits.maximum_accumulated_plastic_strain = 0.5;
        let table = TabulatedHardening::new(vec![
            HardeningPoint { plastic_strain: 0.0, yield_stress_pa: YIELD },
            HardeningPoint { plastic_strain: 0.04, yield_stress_pa: 5e6 },
            HardeningPoint { plastic_strain: 0.1, yield_stress_pa: 6e6 },
            HardeningPoint { plastic_strain: 0.5, yield_stress_pa: 6e6 },
        ]).unwrap();
        let material = FiniteJ2::new(parameters, limits).unwrap()
            .with_tabulated_hardening(table).unwrap();
        let mut solid = FinitePlasticTetSolid::try_new(
            &nodes, &tets, FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(), cx,
        ).unwrap();
        let mut physical_heat = 0.0;
        // controls() bounds each step's endpoint-work defect by 20 kJ, about half
        // the step's stress rise times its strain step. This table rises at
        // 75 MPa per unit plastic strain, so a 0.03 step from rest leaves
        // 33 kJ; 0.015 steps reach the same stretch inside that gate.
        for index in 1..=12 {
            let a = -0.015 * f64::from(index);
            let step = solid.advance(
                &rollers(&nodes, Some(a.exp())), &[], [0.0;3],
                HyperTetSettings::default(), controls(), cx,
            ).unwrap();
            assert!(step.equilibrium.steps[0].linear_iterations > 0);
            physical_heat += step.energy.plastic_dissipation_j;
            solid = step.solid;
        }
        // Final state lies on the measured plateau: exact free-sided Hencky
        // reduction gives |log(lambda)| = alpha + Y_plateau/E.
        let alpha = 0.18 - 6e6 / E;
        let hardening_storage = 0.5 * 3e6 * 0.04
            + 0.5 * (3e6 + 4e6) * 0.06 + 4e6 * (alpha - 0.1);
        check_uniaxial_with_hardening(&solid, -0.18, -alpha, alpha, hardening_storage);
        close(physical_heat, YIELD * alpha, 0.005);
        // Release the right face and unload its actual dead traction in four
        // parts, as the Voce coupon does: one free release from 6 MPa would put
        // the whole elastic recovery, about 87.5 kJ, into one endpoint step.
        let loaded_nominal = E * (-0.18 + alpha) / (-0.18_f64).exp();
        let mut release_heat = 0.0;
        for fraction in [0.75, 0.5, 0.25, 0.0] {
            let nominal = fraction * loaded_nominal;
            let forces = [
                (3, nominal / 3.0),
                (9, nominal / 6.0),
                (21, nominal / 3.0),
                (15, nominal / 6.0),
            ];
            let step = solid
                .advance(
                    &rollers(&nodes, None),
                    &forces,
                    [0.0; 3],
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            release_heat += step.energy.plastic_dissipation_j;
            solid = step.solid;
        }
        check_uniaxial_with_hardening(&solid, -alpha, -alpha, alpha, hardening_storage);
        close(release_heat, 0.0, 0.005);
        for (node, x) in nodes.iter().enumerate() {
            close(solid.displacement_m()[node][0], ((-alpha).exp()-1.0)*x[0], 3e-8);
            close(solid.displacement_m()[node][1], ((0.5*alpha).exp()-1.0)*x[1], 3e-8);
        }
    });
}

#[test]
fn voce_coupon_solves_lateral_contraction_unloading_and_reverse_flow_with_stored_saturation() {
    with_cx(|cx| {
        let (nodes, tets) = cube();
        let saturation = 6.0e6;
        let rate = 20.0;
        let material = law().with_voce_hardening(saturation, rate).unwrap();
        let strength = |alpha: f64| YIELD + HARDENING * alpha - saturation * (-rate * alpha).exp_m1();
        let hardening_storage = |alpha: f64| {
            0.5 * HARDENING * alpha * alpha
                + saturation * (alpha + (-rate * alpha).exp_m1() / rate)
        };
        // Independent scalar reduction of a free-sided Hencky coupon:
        // |E(a-p_old)|-E*dp = Y+H(alpha_old+dp)+Q(1-exp(-b(alpha_old+dp))).
        // Bisection of this monotone equation does not use the native return
        // map, its 3mu predictor, its yield diagnostic or its hardening helper.
        let reference = |a: f64, p: f64, alpha: f64| {
            let trial = E * (a - p);
            if trial.abs() <= strength(alpha) {
                return (p, alpha);
            }
            let (mut lo, mut hi) = (0.0, trial.abs() / E);
            for _ in 0..80 {
                let increment = f64::midpoint(lo, hi);
                if trial.abs() - E * increment > strength(alpha + increment) {
                    lo = increment;
                } else {
                    hi = increment;
                }
            }
            let increment = f64::midpoint(lo, hi);
            (p + trial.signum() * increment, alpha + increment)
        };
        let mut solid = FinitePlasticTetSolid::try_new(
            &nodes,
            &tets,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let (mut old_a, mut p, mut alpha, mut storage, mut heat) = (0.0_f64, 0.0, 0.0, 0.0, 0.0);
        for reverse in [false, true] {
            let start = if reverse { p } else { 0.0 };
            let count = if reverse { 12 } else { 10 };
            for index in 1..=count {
                let a = start + if reverse { -0.012 } else { 0.01 } * f64::from(index);
                let (next_p, next_alpha) = reference(a, p, alpha);
                let tau = E * (a - next_p);
                let next_storage = 0.5 * tau * tau / E + hardening_storage(next_alpha);
                let dissipation = YIELD * (next_alpha - alpha);
                let work = (tau / a.exp()) * (a.exp() - old_a.exp());
                let step = solid
                    .advance(
                        &rollers(&nodes, Some(a.exp())),
                        &[],
                        [0.0; 3],
                        HyperTetSettings::default(),
                        controls(),
                        cx,
                    )
                    .unwrap();
                // Lateral-face coordinates and the interior vertex remain
                // unknown: this checks the assembled nonlinear equilibrium.
                assert!(step.equilibrium.steps[0].linear_iterations > 0);
                check_uniaxial_with_hardening(
                    &step.solid,
                    a,
                    next_p,
                    next_alpha,
                    hardening_storage(next_alpha),
                );
                close(step.energy.stored_before_j, storage, 0.01);
                close(step.energy.stored_after_j, next_storage, 0.01);
                close(step.energy.plastic_dissipation_j, dissipation, 0.005);
                close(step.energy.support_endpoint_work_j, work, 0.005);
                close(
                    step.energy.endpoint_work_defect_j,
                    work - (next_storage - storage) - dissipation,
                    0.02,
                );
                close(
                    step.equilibrium.strain_energy_j,
                    next_storage + dissipation,
                    0.01,
                );
                (old_a, p, alpha, storage) = (a, next_p, next_alpha, next_storage);
                heat += step.energy.plastic_dissipation_j;
                solid = step.solid;
            }
            if !reverse {
                assert!(p > 0.06, "the coupon must retain substantial finite plastic stretch");
                assert!(hardening_storage(alpha) - 0.5 * HARDENING * alpha * alpha > 150_000.0);
                // Remove the right displacement restraint and unload actual
                // dead face traction. All right-face axial DOFs are now free.
                let loaded_a = old_a;
                let loaded_nominal = E * (loaded_a - p) / loaded_a.exp();
                for fraction in [0.75, 0.5, 0.25, 0.0] {
                    let nominal = fraction * loaded_nominal;
                    let (mut lo, mut hi) = (p, loaded_a);
                    for _ in 0..80 {
                        let a = f64::midpoint(lo, hi);
                        if E * (a - p) / a.exp() < nominal {
                            lo = a;
                        } else {
                            hi = a;
                        }
                    }
                    let a = f64::midpoint(lo, hi);
                    let forces = [
                        (3, nominal / 3.0),
                        (9, nominal / 6.0),
                        (21, nominal / 3.0),
                        (15, nominal / 6.0),
                    ];
                    let step = solid
                        .advance(
                            &rollers(&nodes, None),
                            &forces,
                            [0.0; 3],
                            HyperTetSettings::default(),
                            controls(),
                            cx,
                        )
                        .unwrap();
                    check_uniaxial_with_hardening(&step.solid, a, p, alpha, hardening_storage(alpha));
                    close(step.energy.plastic_dissipation_j, 0.0, 0.005);
                    close(step.energy.support_endpoint_work_j, 0.0, 0.0);
                    close(
                        step.energy.dead_load_endpoint_work_j,
                        nominal * (a.exp() - old_a.exp()),
                        0.005,
                    );
                    old_a = a;
                    storage = 0.5 * E * (a - p).powi(2) + hardening_storage(alpha);
                    solid = step.solid;
                }
                close(old_a, p, 1e-14);
                close(solid.stored_energy_j(), hardening_storage(alpha), 0.01);
            }
        }
        assert!(p < 0.0 && alpha > 0.15, "reversal changes the metric while accumulation grows");
        close(heat, YIELD * alpha, 0.02);
    });
}

#[test]
fn finite_loading_dead_load_unloading_permanent_set_rotation_and_reverse_flow() {
    with_cx(|cx| {
        let (nodes, tets) = cube();
        let material = law();
        let mut solid = FinitePlasticTetSolid::try_new(
            &nodes,
            &tets,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let mut old_alpha = 0.0;
        let mut old_stretch = 1.0;
        let mut old_storage = 0.0;
        let mut dissipation = 0.0;
        for index in 1..=12 {
            let a = 0.01 * index as f64;
            let alpha = ((E * a - YIELD) / (E + HARDENING)).max(0.0);
            let stretch = a.exp();
            let supports = rollers(&nodes, Some(stretch));
            let step = solid
                .advance(
                    &supports,
                    &[],
                    [0.0; 3],
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_uniaxial(&step.solid, a, alpha, alpha);
            let nominal = E * (a - alpha) / stretch;
            let expected_work = nominal * (stretch - old_stretch);
            let expected_dissipation = YIELD * (alpha - old_alpha);
            let expected_storage = 0.5 * E * (a - alpha).powi(2) + 0.5 * HARDENING * alpha * alpha;
            close(step.energy.support_endpoint_work_j, expected_work, 0.005);
            close(step.energy.dead_load_endpoint_work_j, 0.0, 0.0);
            close(
                step.energy.plastic_dissipation_j,
                expected_dissipation,
                0.005,
            );
            close(
                step.energy.endpoint_work_defect_j,
                expected_work - (expected_storage - old_storage) - expected_dissipation,
                0.01,
            );
            close(
                step.equilibrium.strain_energy_j,
                step.energy.stored_after_j + step.energy.plastic_dissipation_j,
                1.0e-4,
            );
            assert!(step.equilibrium.steps[0].linear_iterations > 0);
            dissipation += step.energy.plastic_dissipation_j;
            solid = step.solid;
            old_alpha = alpha;
            old_stretch = stretch;
            old_storage = expected_storage;
        }
        close(dissipation, YIELD * old_alpha, 0.02);
        assert!(
            old_alpha > 0.10,
            "fixture must leave the small-strain regime"
        );

        // Release the right support and remove its actual nodal traction in four
        // force-controlled endpoints. Face areas are 1/2, hence corner weights
        // 1/3,1/6,1/3,1/6 on the two original triangles of the right face.
        let nominal_loaded = (YIELD + HARDENING * old_alpha) / old_stretch;
        let supports = rollers(&nodes, None);
        for fraction in [0.75, 0.5, 0.25, 0.0] {
            let nominal = fraction * nominal_loaded;
            let forces = [
                (3, nominal / 3.0),
                (9, nominal / 6.0),
                (21, nominal / 3.0),
                (15, nominal / 6.0),
            ];
            let step = solid
                .advance(
                    &supports,
                    &forces,
                    [0.0; 3],
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            for state in step.solid.element_states() {
                close(state.alpha, old_alpha, 3.0e-8);
            }
            close(step.energy.plastic_dissipation_j, 0.0, 0.005);
            close(step.energy.support_endpoint_work_j, 0.0, 0.0);
            assert!(step.energy.dead_load_endpoint_work_j <= 1.0e-5);
            for stress in step.solid.element_first_piola_pa() {
                close(stress[0], nominal, 0.15);
            }
            solid = step.solid;
        }
        check_uniaxial(&solid, old_alpha, old_alpha, old_alpha);
        assert!(
            solid.displacement_m()[1][0] > 0.10,
            "zero force must retain a permanent set"
        );

        // A finite rotation of the stress-free permanent set preserves the
        // material metric, zero stress and the retained hardening energy.
        let (sin, cos) = 0.47_f64.sin_cos();
        let (axial, lateral) = (old_alpha.exp(), (-0.5 * old_alpha).exp());
        let rotated_f = [
            cos * axial,
            -sin * lateral,
            0.0,
            sin * axial,
            cos * lateral,
            0.0,
            0.0,
            0.0,
            lateral,
        ];
        let rotation = solid
            .advance(
                &affine(&nodes, rotated_f, false),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                controls(),
                cx,
            )
            .unwrap();
        for f in rotation.solid.element_deformation() {
            for i in 0..9 {
                close(f[i], rotated_f[i], 2.0e-8);
            }
        }
        for stress in rotation.solid.element_first_piola_pa() {
            for value in stress {
                close(*value, 0.0, 0.15);
            }
        }
        close(
            rotation.energy.stored_after_j,
            solid.stored_energy_j(),
            0.005,
        );
        close(rotation.energy.plastic_dissipation_j, 0.0, 0.005);
        close(rotation.energy.endpoint_work_defect_j, 0.0, 0.01);
        for (before, after) in solid
            .element_states()
            .iter()
            .zip(rotation.solid.element_states())
        {
            close(before.alpha, after.alpha, 3.0e-8);
            for i in 0..6 {
                close(
                    before.inverse_plastic_cauchy_green[i],
                    after.inverse_plastic_cauchy_green[i],
                    1.0e-7,
                );
            }
        }

        // Reverse the unrotated accepted path. Accumulated flow increases even
        // when the net plastic metric returns toward its virgin value.
        let start_alpha = old_alpha;
        let mut reverse_dissipation = 0.0;
        for index in 1..=12 {
            let a = start_alpha - 0.012 * index as f64;
            let reverse = ((E * (start_alpha - a) - YIELD - HARDENING * start_alpha)
                / (E + HARDENING))
                .max(0.0);
            let alpha = start_alpha + reverse;
            let plastic_axial = start_alpha - reverse;
            let step = solid
                .advance(
                    &rollers(&nodes, Some(a.exp())),
                    &[],
                    [0.0; 3],
                    HyperTetSettings::default(),
                    controls(),
                    cx,
                )
                .unwrap();
            check_uniaxial(&step.solid, a, plastic_axial, alpha);
            reverse_dissipation += step.energy.plastic_dissipation_j;
            solid = step.solid;
            old_alpha = alpha;
        }
        assert!(old_alpha > 2.0 * start_alpha);
        close(reverse_dissipation, YIELD * (old_alpha - start_alpha), 0.02);

        // A rejected physical-work budget after equilibrium cannot commit any
        // part of the proposed history or the proposed displacement.
        let states_before = solid.element_states().to_vec();
        let displacement_before = solid.displacement_m().to_vec();
        let mut tight = controls();
        tight.maximum_endpoint_work_defect_j = 0.0;
        let refused = solid
            .advance(
                &rollers(&nodes, Some((-0.05_f64).exp())),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                tight,
                cx,
            )
            .unwrap_err();
        assert!(matches!(
            refused,
            FinitePlasticTetError::EnergyGate {
                quantity: "signed endpoint integration defect",
                ..
            }
        ));
        assert_eq!(solid.element_states(), states_before);
        assert_eq!(solid.displacement_m(), displacement_before);
    });
}

#[test]
fn prestressed_rotation_reports_negative_endpoint_defect_without_plastic_heat() {
    with_cx(|cx| {
        let (nodes, tets) = cube();
        let material = law();
        let solid = FinitePlasticTetSolid::try_new(
            &nodes,
            &tets,
            FinitePlasticLawField::Uniform(&material),
            TetAssemblyBudget::standard(),
            cx,
        )
        .unwrap();
        let mut generous = controls();
        generous.maximum_endpoint_work_defect_j = 1.0e8;
        let stretch = 0.95;
        let compressed = solid
            .advance(
                &affine(
                    &nodes,
                    [stretch, 0.0, 0.0, 0.0, stretch, 0.0, 0.0, 0.0, stretch],
                    true,
                ),
                &[],
                [0.0; 3],
                HyperTetSettings::default(),
                generous,
                cx,
            )
            .unwrap()
            .solid;
        let (sin, cos) = 0.35_f64.sin_cos();
        let f = [
            cos * stretch,
            -sin * stretch,
            0.0,
            sin * stretch,
            cos * stretch,
            0.0,
            0.0,
            0.0,
            stretch,
        ];
        let step = compressed
            .advance(
                &affine(&nodes, f, true),
                &[(0, 50.0), (4, -80.0)],
                [100.0, 200.0, -150.0],
                HyperTetSettings::default(),
                generous,
                cx,
            )
            .unwrap();
        let kirchhoff_pressure = 3.0 * KAPPA * stretch.ln();
        let signed_defect = 2.0 * kirchhoff_pressure * (1.0 - cos);
        assert!(signed_defect < -1.0e6);
        close(step.energy.endpoint_work_defect_j, signed_defect, 1.0e-5);
        close(step.energy.plastic_dissipation_j, 0.0, 0.0);
        close(
            step.energy.stored_after_j,
            compressed.stored_energy_j(),
            1.0e-6,
        );
        close(
            step.energy.dead_load_endpoint_work_j + step.energy.support_endpoint_work_j,
            step.energy.internal_endpoint_work_j,
            1.0e-5,
        );
        for state in step.solid.element_states() {
            assert_eq!(state.alpha, 0.0);
        }
    });
}
