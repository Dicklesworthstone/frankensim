//! Temperature-dependent sourced E/H drive loaded native 3D shape and heat.
use super::*;

const SIGMA: f64 = 3.4e6;
const NU: f64 = 0.2;

fn affine_property(name: &str, reference: f64, slope: f64) -> PropertyClaim {
    let mut claim = scalar(name, Pressure::DIMS, reference);
    claim.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: [280.0, 300.0, 310.0, 320.0]
            .map(|t| (t, reference + slope * (t - 300.0)))
            .to_vec(),
        dims: Pressure::DIMS,
    };
    claim.interpolation = InterpolationPolicy::LinearInside;
    claim.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    claim
}

#[derive(Clone, Copy)]
struct Oracle {
    t: f64,
    youngs: f64,
    hardening: f64,
    yield_pa: f64,
    p: f64,
    elastic: f64,
    strain: f64,
    free_energy: f64,
    entropy: f64,
}
fn endpoint(t: f64, de: f64, dh: f64) -> Oracle {
    let youngs = E + de * (t - 300.0);
    let hardening = H + dh * (t - 300.0);
    let yield_pa = YIELD - SOFTENING * (t - 300.0);
    let p = (SIGMA - yield_pa) / hardening;
    let elastic = SIGMA / youngs;
    Oracle {
        t,
        youngs,
        hardening,
        yield_pa,
        p,
        elastic,
        strain: elastic + p,
        free_energy: 0.5 * youngs * elastic * elastic + 0.5 * hardening * p * p,
        entropy: -0.5 * de * elastic * elastic - 0.5 * dh * p * p,
    }
}
fn initial_oracle() -> Oracle {
    Oracle {
        t: 300.0,
        youngs: E,
        hardening: H,
        yield_pa: YIELD,
        p: 0.0,
        elastic: 0.0,
        strain: 0.0,
        free_energy: 0.0,
        entropy: 0.0,
    }
}
fn coefficient_change(old: Oracle, new: Oracle) -> f64 {
    0.5 * (new.youngs - old.youngs) * old.elastic * old.elastic
        + 0.5 * (new.hardening - old.hardening) * old.p * old.p
}
fn reversible_heat(old: Oracle, new: Oracle) -> f64 {
    -coefficient_change(old, new) - (new.t * new.entropy - old.t * old.entropy)
}
fn root(mut lower: f64, mut upper: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lower) < 0.0 && residual(upper) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lower + upper);
        if residual(mid) > 0.0 {
            upper = mid;
        } else {
            lower = mid;
        }
    }
    0.5 * (lower + upper)
}

#[test]
fn g1_source_modulus_swap_changes_loaded_tet_shape_plastic_history_and_reciprocal_heat() {
    context(|cx| {
        let cp = 100.0;
        let duration = 0.001;
        let heat = RHO * V * cp * 6.0;
        let force = [(3, SIGMA * V / L)];
        let mut responses = Vec::new();
        for (de, dh) in [(0.0, 0.0), (-1e7, -1e6)] {
            let fixture = Fixture::with_constitutive_claims(
                cp,
                None,
                affine_property(YOUNG_MODULUS_PROPERTY, E, de),
                affine_property(BASIS.property_name(), H, dh),
            );
            let thermal = initial(&[&fixture], &[], &[heat / duration]);
            let solid = bind(&thermal, &[&fixture], &NODES, &TETS, cx);
            let original_solution = solid.solution().clone();
            let original_cells = thermal.cells().clone();
            let first = solid
                .advance(
                    &thermal,
                    &(),
                    duration,
                    &ROLLERS,
                    &force,
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let second = first
                .solid()
                .advance(
                    first.thermal(),
                    &(),
                    2.0 * duration,
                    &ROLLERS,
                    &force,
                    IncrementSettings::default(),
                    thermal_config(),
                    config(),
                    cx,
                )
                .unwrap();
            let mut old = initial_oracle();
            for candidate in [&first, &second] {
                let expected_t = root(old.t, 319.0, |t| {
                    let new = endpoint(t, de, dh);
                    RHO * cp * (t - old.t)
                        - heat / V
                        - new.yield_pa * (new.p - old.p)
                        - reversible_heat(old, new)
                });
                let expected = endpoint(expected_t, de, dh);
                let t = candidate.thermal().temperature_k(CellId(100)).unwrap();
                close(t, expected_t, 3e-8);
                let solution = candidate.solid().solution();
                close(solution.element_states[0].alpha, expected.p, 2e-9);
                close(solution.element_stress[0][0], SIGMA, 2e-3);
                close(solution.element_stress[0][1], 0.0, 2e-3);
                close(solution.element_stress[0][2], 0.0, 2e-3);
                close(solution.displacement_m[1][0], L * expected.strain, 2e-11);
                for (node, axis) in [(2, 1), (3, 2)] {
                    close(
                        solution.displacement_m[node][axis],
                        L * (-NU * expected.elastic - 0.5 * expected.p),
                        2e-11,
                    );
                }
                let resolved = &candidate.solid().materials()[0];
                close(resolved.law().elastic.youngs, expected.youngs, 1.0);
                close(resolved.law().hardening, expected.hardening, 0.1);
                let report = candidate.report();
                close(
                    report.coefficient_storage_change_j,
                    V * coefficient_change(old, expected),
                    1e-8,
                );
                close(
                    report.stored_mechanical_change_j,
                    V * (expected.free_energy - old.free_energy),
                    1e-8,
                );
                close(
                    report.thermoelastic_storage_change_j,
                    V * (expected.t * expected.entropy - old.t * old.entropy),
                    1e-8,
                );
                close(
                    report.reversible_heat_j,
                    V * reversible_heat(old, expected),
                    1e-8,
                );
                close(
                    report.physical_heat_j,
                    V * expected.yield_pa * (expected.p - old.p),
                    1e-8,
                );
                close(
                    report.mechanical_work_j,
                    V * SIGMA * (expected.strain - old.strain),
                    1e-8,
                );
                let delta_e = expected.elastic - old.elastic;
                let delta_p = expected.p - old.p;
                close(
                    report.numerical_loss_j,
                    V * (0.5 * expected.youngs * delta_e * delta_e
                        + 0.5 * expected.hardening * delta_p * delta_p),
                    1e-8,
                );
                close(report.combined_residual_j, 0.0, 1e-8);
                assert!(report.numerical_loss_j > 0.0);
                assert_eq!(candidate.thermal().cells()[&CellId(100)].mass_kg, RHO * V);
                old = expected;
            }
            if de < 0.0 {
                assert!(second.report().coefficient_storage_change_j < 0.0);
                assert!(second.report().reversible_heat_j < 0.0);
            } else {
                assert_eq!(second.report().coefficient_storage_change_j, 0.0);
                assert_eq!(second.report().reversible_heat_j, 0.0);
            }
            assert_eq!(solid.solution(), &original_solution);
            assert_eq!(thermal.cells(), &original_cells);
            assert_eq!(solid.time_s(), 0.0);
            responses.push((old.t, old.strain, old.p));
        }
        assert!(responses[1].1 > responses[0].1);
        assert!(responses[1].2 > responses[0].2);
        assert!((responses[1].0 - responses[0].0).abs() > 1e-3);
    });
}
