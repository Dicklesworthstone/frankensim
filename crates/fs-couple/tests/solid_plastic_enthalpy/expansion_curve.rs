//! A material alpha(T) curve drives native 3D displacement, stress and heat.
use super::*;
use fs_couple::solid_enthalpy::plastic::ThermoplasticError;
use fs_material::plastic::thermal::FeedbackError;
use fs_material::state_point::{
    INVERSE_TEMPERATURE_DIMS, LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA,
};

fn expansion(knots: Vec<(f64, f64)>) -> PropertyClaim {
    let mut claim = scalar(ALPHA, INVERSE_TEMPERATURE_DIMS, 0.0);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots,
        dims: INVERSE_TEMPERATURE_DIMS,
    };
    claim.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim
}
fn fixture(cp: f64, alpha: f64) -> Fixture {
    Fixture::with_expansion_claim(
        cp,
        Some(expansion(vec![
            (280.0, alpha),
            (300.0, alpha),
            (304.0, 2.0 * alpha),
            (310.0, 3.0 * alpha),
            (320.0, 3.0 * alpha),
        ])),
    )
}
// Independent analytical primitive and coefficient for these three segments.
fn free_strain(t: f64, alpha: f64) -> f64 {
    if t <= 304.0 {
        let d = t - 300.0;
        alpha * (d + d * d / 8.0)
    } else if t <= 310.0 {
        let d = t - 304.0;
        alpha * (6.0 + 2.0 * d + d * d / 12.0)
    } else {
        alpha * (21.0 + 3.0 * (t - 310.0))
    }
}
fn coefficient(t: f64, alpha: f64) -> f64 {
    if t <= 304.0 {
        alpha * (1.0 + (t - 300.0) / 4.0)
    } else if t <= 310.0 {
        alpha * (2.0 + (t - 304.0) / 6.0)
    } else {
        3.0 * alpha
    }
}
fn root(mut residual: impl FnMut(f64) -> f64) -> f64 {
    let (mut lower, mut upper) = (300.0, 319.0);
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..80 {
        let middle = 0.5 * (lower + upper);
        if residual(middle) > 0.0 {
            upper = middle;
        } else {
            lower = middle;
        }
    }
    0.5 * (lower + upper)
}

#[test]
fn material_curve_swap_changes_free_shape_constrained_stress_and_reciprocal_heat() {
    context(|cx| {
        let cp = 500.0;
        let duration = 0.001;
        let heat = RHO * V * cp * 12.0;
        let fixed: Vec<_> = (0..12).map(|dof| (dof, 0.0)).collect();
        let nine_bulk = 3.0 * E / (1.0 - 0.4);
        let mut responses = Vec::new();
        for alpha in [1e-5, 2e-5] {
            let fixture = fixture(cp, alpha);
            let thermal = initial(&[&fixture], &[], &[heat / duration]);
            let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
            let original = solid.solution().clone();
            let cells = thermal.cells().clone();
            let free = solid
                .advance(
                    &thermal,
                    &(),
                    duration,
                    &ROLLERS,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            close(
                free.thermal().temperature_k(CellId(100)).unwrap(),
                312.0,
                1e-8,
            );
            for (node, axis) in [(1, 0), (2, 1), (3, 2)] {
                close(
                    free.solid().solution().displacement_m[node][axis],
                    L * free_strain(312.0, alpha),
                    1e-12,
                );
            }
            close(free.report().reversible_heat_j, 0.0, 1e-9);

            let clamped = solid
                .advance(
                    &thermal,
                    &(),
                    duration,
                    &fixed,
                    &[],
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let oracle = root(|t| {
                let f = free_strain(t, alpha);
                let reversible = nine_bulk * f * (t * coefficient(t, alpha) - f);
                RHO * cp * (t - 300.0) - heat / V - reversible
            });
            let t = clamped.thermal().temperature_k(CellId(100)).unwrap();
            close(t, oracle, 3e-8);
            let f = free_strain(t, alpha);
            let trace_sigma = -nine_bulk * f;
            for stress in &clamped.solid().solution().element_stress[0][..3] {
                close(*stress, trace_sigma / 3.0, 2e-3);
            }
            let storage = V * t * coefficient(t, alpha) * trace_sigma;
            let reversible = V * trace_sigma * f - storage;
            close(
                clamped.report().thermoelastic_storage_change_j,
                storage,
                1e-8,
            );
            close(clamped.report().reversible_heat_j, reversible, 1e-8);
            close(clamped.report().physical_heat_j, 0.0, 1e-12);
            close(clamped.report().combined_residual_j, 0.0, 1e-8);
            close(clamped.report().mechanical_work_j, 0.0, 1e-12);
            assert_eq!(clamped.thermal().cells()[&CellId(100)].mass_kg, RHO * V);
            assert_eq!(solid.solution(), &original);
            assert_eq!(thermal.cells(), &cells);
            assert_eq!(solid.time_s(), 0.0);
            responses.push((
                t,
                trace_sigma.abs(),
                reversible,
                clamped.solid().materials()[0].material().identity(),
            ));
        }
        assert!(responses[1].0 > responses[0].0);
        assert!(responses[1].1 > responses[0].1);
        assert!(responses[1].2 > responses[0].2);
        assert_ne!(responses[0].3, responses[1].3);
    });
}

#[test]
fn variable_expansion_drives_actual_yielding_and_the_existing_plastic_heat_account() {
    context(|cx| {
        let cp = 10000.0;
        let alpha = 1e-4;
        let duration = 0.001;
        let heat = RHO * V * cp * 12.0;
        let fixture = fixture(cp, alpha);
        let thermal = initial(&[&fixture], &[], &[heat / duration]);
        let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
        let mut prescribed = ROLLERS.to_vec();
        prescribed.push((3, 0.0));
        prescribed.sort_by_key(|entry| entry.0);
        let candidate = solid
            .advance(
                &thermal,
                &(),
                duration,
                &prescribed,
                &[],
                IncrementSettings::default(),
                thermal_config(),
                config(),
                cx,
            )
            .unwrap();
        // Uniform constrained axial strain with transverse traction zero:
        // p=(E*f-yield)/(E+H), sigma_x=-(yield+H*p).
        let oracle = root(|t| {
            let f = free_strain(t, alpha);
            let yield_pa = YIELD - SOFTENING * (t - 300.0);
            let p = ((E * f - yield_pa) / (E + H)).max(0.0);
            let sigma = E * (p - f);
            let heat_density = yield_pa * p + sigma * (f - t * coefficient(t, alpha));
            RHO * cp * (t - 300.0) - heat / V - heat_density
        });
        let t = candidate.thermal().temperature_k(CellId(100)).unwrap();
        close(t, oracle, 3e-8);
        let f = free_strain(t, alpha);
        let yield_pa = YIELD - SOFTENING * (t - 300.0);
        let p = (E * f - yield_pa) / (E + H);
        assert!(p > 0.0);
        let solution = candidate.solid().solution();
        close(solution.element_states[0].alpha, p, 2e-10);
        close(solution.element_stress[0][0], -(yield_pa + H * p), 2e-3);
        close(solution.element_stress[0][1], 0.0, 2e-3);
        close(solution.element_stress[0][2], 0.0, 2e-3);
        close(candidate.report().physical_heat_j, V * yield_pa * p, 1e-8);
        close(candidate.report().mechanical_work_j, 0.0, 1e-12);
        close(candidate.report().combined_residual_j, 0.0, 1e-8);
        assert!(candidate.report().numerical_loss_j > 0.0);
        assert_eq!(solid.solution().element_states[0].alpha, 0.0);
    });
}

#[test]
fn failed_capacity_admission_leaves_the_accepted_thermal_inventory_untouched() {
    context(|cx| {
        let fixture = Fixture::with_expansion_claim(
            500.0,
            Some(expansion(vec![
                (280.0, 0.0),
                (300.0, 0.0),
                (301.0, 1e-4),
                (302.0, 0.0),
                (320.0, 0.0),
            ])),
        );
        let thermal = initial(&[&fixture], &[], &[0.0]);
        let cells = thermal.cells().clone();
        let result = ThermoplasticSolid::try_new(
            &thermal,
            &NODES,
            &TETS,
            vec![ThermoplasticCell {
                cell: CellId(100),
                query: fixture.query(),
            }],
            TetAssemblyBudget::standard(),
            cx,
        );
        assert!(matches!(result,
            Err(ThermoplasticError::Material(FeedbackError::Invalid(message)))
                if message.contains("effective heat capacity")));
        assert_eq!(thermal.cells(), &cells);
        assert_eq!(thermal.time_s(), 0.0);
    });
}
