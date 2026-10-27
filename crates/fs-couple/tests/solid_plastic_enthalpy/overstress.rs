//! Actual loaded-tet creep, temperature feedback and material substitution.
use super::*;
use fs_couple::solid_enthalpy::plastic::dynamics::ThermoplasticDynamics;
use fs_material::SmallStrainLaw;
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::rate::{OVERSTRESS_VISCOSITY_DIMS, OVERSTRESS_VISCOSITY_PROPERTY};
use fs_material::plastic::thermal::ThermoplasticSource;

fn fixture(cp: f64, eta: f64, temperature_slope: f64) -> Fixture {
    let original = Fixture::new(cp);
    let mut claims = ClaimSet::new();
    for (_, claim) in original.card.claims().claims_ordered() {
        claims.insert_claim(claim.clone()).unwrap();
    }
    let mut viscosity = scalar(
        OVERSTRESS_VISCOSITY_PROPERTY,
        OVERSTRESS_VISCOSITY_DIMS,
        eta,
    );
    viscosity.value = PropertyValue::Curve {
        abscissa: "T".into(),
        abscissa_dims: Temperature::DIMS,
        knots: vec![
            (280.0, eta * (1.0 + 20.0 * temperature_slope)),
            (320.0, eta * (1.0 - 20.0 * temperature_slope)),
        ],
        dims: OVERSTRESS_VISCOSITY_DIMS,
    };
    viscosity.interpolation = InterpolationPolicy::LinearInside;
    viscosity.validity = ValidityDomain::unconstrained().with("T", 280.0, 320.0);
    claims.insert_claim(viscosity).unwrap();
    let card = MaterialCard::assemble(original.card.id().clone(), claims, vec![]).unwrap();
    let lower = QueryPoint::new().with("T", 280.0).unwrap();
    let upper = QueryPoint::new().with("T", 320.0).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(
        &card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new(CP, CP_DIMS),
            density: &PropertyKey::new(DENSITY_PROPERTY, Density::DIMS),
            temperature_axis: "T",
            lower: &lower,
            upper: &upper,
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        -20.0 * cp,
        1e-5,
    )
    .unwrap();
    Fixture {
        card,
        caloric,
        expansion: false,
    }
}

fn rate_solid<'a>(
    f: &'a Fixture,
    thermal: &DepositionSimulation<'_>,
    cx: &Cx<'_>,
) -> ThermoplasticSolid<'a> {
    bind(thermal, &[f], &NODES, &TETS, cx)
        .with_linear_overstress(&[Some(MaterialPropertySelection::SingleClaimOnly)], cx)
        .unwrap()
}
fn heat_steps(dt: f64) -> StepConfig {
    StepConfig {
        max_dt_s: dt,
        ..thermal_config()
    }
}

// Independent uniaxial endpoint balance. Free transverse faces enforce sigma_yy
// = sigma_zz = 0. The true plastic strain increment is delta_alpha. Fixed E/H,
// no expansion: rho Cp dT=(sigma-H alpha_end) delta_alpha. Solve the scalar
// Perzyna equation using these formulae only, not a production return map.
fn reference(
    alpha0: f64,
    t0: f64,
    stress: f64,
    eta: f64,
    slope: f64,
    cp: f64,
    dt: f64,
) -> (f64, f64) {
    let mut lo = 0.0;
    let mut hi = 0.03 - alpha0;
    for _ in 0..100 {
        let delta = 0.5 * (lo + hi);
        let t = t0 + (stress - H * (alpha0 + delta)) * delta / (RHO * cp);
        let viscosity = eta * (1.0 - slope * (t - 300.0));
        let yield_pa = YIELD - SOFTENING * (t - 300.0);
        let residual = (H + viscosity / dt) * delta - (stress - yield_pa - H * alpha0);
        if residual > 0.0 {
            hi = delta;
        } else {
            lo = delta;
        }
    }
    let delta = 0.5 * (lo + hi);
    (
        alpha0 + delta,
        t0 + (stress - H * (alpha0 + delta)) * delta / (RHO * cp),
    )
}

#[test]
fn g1_held_load_creeps_and_heats_the_actual_tetrahedron_from_material_viscosity() {
    context(|cx| {
        let cp = 20.0;
        let dt = 0.1;
        let stress = 5e6;
        let slope = 0.02;
        let mut endpoints = Vec::new();
        for eta in [1e8, 2e8] {
            for sign in [1.0, -1.0] {
                let f = fixture(cp, eta, slope);
                let mut thermal = initial(&[&f], &[], &[0.0]);
                let mut solid = rate_solid(&f, &thermal, cx);
                let initial_mass = thermal.cells()[&CellId(100)].mass_kg;
                let mut alpha = 0.0;
                let mut temperature = 300.0;
                let mut previous_u = 0.0;
                for n in 1..=5 {
                    let expected = reference(alpha, temperature, stress, eta, slope, cp, dt);
                    let next = solid
                        .advance(
                            &thermal,
                            &(),
                            n as f64 * dt,
                            &ROLLERS,
                            &[(3, sign * stress * V / L)],
                            IncrementSettings::default(),
                            heat_steps(dt),
                            config(),
                            cx,
                        )
                        .unwrap();
                    let solution = next.solid().solution();
                    let report = next.report();
                    close(solution.element_states[0].alpha, expected.0, 2e-10);
                    close(
                        next.thermal().temperature_k(CellId(100)).unwrap(),
                        expected.1,
                        2e-7,
                    );
                    close(solution.element_stress[0][0], sign * stress, 1e-2);
                    close(
                        solution.element_strain[0][0],
                        sign * (stress / E + expected.0),
                        2e-10,
                    );
                    for value in &solution.element_stress[0][1..] {
                        close(*value, 0.0, 1e-2);
                    }
                    close(
                        report.physical_heat_j,
                        RHO * V * cp * (expected.1 - temperature),
                        1e-9,
                    );
                    close(
                        report.thermal_audit.stored_energy_change_j,
                        report.physical_heat_j,
                        1e-9,
                    );
                    close(report.combined_residual_j, 0.0, config().energy_tolerance_j);
                    let viscosity = next.solid().materials()[0]
                        .material()
                        .property(OVERSTRESS_VISCOSITY_PROPERTY)
                        .unwrap();
                    close(
                        viscosity.value_si(),
                        eta * (1.0 - slope * (expected.1 - 300.0)),
                        3.0,
                    );
                    assert!(solution.displacement_m[1][0].abs() > previous_u);
                    assert!(report.physical_heat_j > 0.0);
                    assert!(report.numerical_loss_j >= 0.0);
                    assert_eq!(next.thermal().cells()[&CellId(100)].mass_kg, initial_mass);
                    previous_u = solution.displacement_m[1][0].abs();
                    alpha = expected.0;
                    temperature = expected.1;
                    let (t, s, _, _) = next.into_parts();
                    thermal = t;
                    solid = s;
                }
                if sign > 0.0 {
                    endpoints.push((alpha, temperature));
                }
            }
        }
        assert!(endpoints[0].0 > endpoints[1].0);
        assert!(endpoints[0].1 > endpoints[1].1);
    });
}

#[test]
fn g3_material_and_physical_time_scaling_preserve_creep_but_heat_substeps_do_not_drive_it() {
    context(|cx| {
        let mut endpoints = Vec::new();
        for (eta, dt, heat_dt) in [(1e8, 0.1, 0.1), (2e8, 0.2, 0.2), (1e8, 0.1, 0.01)] {
            let f = fixture(20.0, eta, 0.02);
            let thermal = initial(&[&f], &[], &[0.0]);
            let solid = rate_solid(&f, &thermal, cx);
            let next = solid
                .advance(
                    &thermal,
                    &(),
                    dt,
                    &ROLLERS,
                    &[(3, 5e6 * V / L)],
                    IncrementSettings::default(),
                    heat_steps(heat_dt),
                    config(),
                    cx,
                )
                .unwrap();
            endpoints.push((
                next.solid().solution().element_states[0].alpha,
                next.thermal().temperature_k(CellId(100)).unwrap(),
            ));
        }
        for next in &endpoints[1..] {
            close(next.0, endpoints[0].0, 1e-10);
            close(next.1, endpoints[0].1, 1e-7);
        }
    });
}

#[test]
fn g1_inertial_owner_uses_same_rate_heat_and_keeps_inertial_loss_separate() {
    context(|cx| {
        let f = fixture(20.0, 1e8, 0.02);
        let thermal = initial(&[&f], &[], &[0.0]);
        let solid = rate_solid(&f, &thermal, cx);
        let dt = 0.1;
        let prescribed: Vec<_> = (0..12)
            .map(|i| {
                (
                    i,
                    match i {
                        3 => L * 0.01,
                        7 | 11 => -L * 0.005,
                        _ => 0.0,
                    },
                )
            })
            .collect();
        let static_step = solid
            .advance(
                &thermal,
                &(),
                dt,
                &prescribed,
                &[],
                IncrementSettings::default(),
                heat_steps(dt),
                config(),
                cx,
            )
            .unwrap();
        let dynamics = ThermoplasticDynamics::try_new(solid, vec![[0.0; 3]; 4], cx).unwrap();
        let dynamic_step = dynamics
            .advance(
                &thermal,
                &(),
                dt,
                &prescribed,
                &[],
                IncrementSettings::default(),
                heat_steps(dt),
                config(),
                cx,
            )
            .unwrap();
        close(
            dynamic_step.report().material.physical_heat_j,
            static_step.report().physical_heat_j,
            1e-9,
        );
        close(
            dynamic_step.thermal().temperature_k(CellId(100)).unwrap(),
            static_step.thermal().temperature_k(CellId(100)).unwrap(),
            1e-8,
        );
        assert!(dynamic_step.report().inertia.kinetic_change_j > 0.0);
        assert!(dynamic_step.report().inertia.numerical_loss_j > 0.0);
        close(
            dynamic_step.report().material.combined_residual_j,
            0.0,
            config().energy_tolerance_j,
        );
    });
}

#[test]
fn g4_duration_free_api_and_late_thermal_failure_do_not_advance_rate_history() {
    context(|cx| {
        // No fallback for missing viscosity or a positive reference value whose
        // pinned source path reaches zero at the caloric interval boundary.
        for refused in [Fixture::new(20.0), fixture(20.0, 1e8, 0.05)] {
            assert!(
                ThermoplasticSource::try_new(refused.query(), 0.0)
                    .unwrap()
                    .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
                    .is_err()
            );
        }
        let f = fixture(20.0, 1e8, 0.02);
        let source = ThermoplasticSource::try_new(f.query(), 0.0)
            .unwrap()
            .with_linear_overstress(MaterialPropertySelection::SingleClaimOnly)
            .unwrap();
        assert!(source.resolve_expanded_enthalpy(0.0).is_err());
        assert!(
            source
                .advance_with_energy(
                    0.0,
                    0.0,
                    &[0.0; 6],
                    &source.reference().initial_state(),
                    &[0.0; 6],
                    EnergyTolerances {
                        energy_j_m3: 1e-6,
                        yield_pa: 1e-3
                    }
                )
                .is_err()
        );
        assert!(source.resolve_expanded_enthalpy_for_step(0.0, 0.0).is_err());
        let thermal = initial(&[&f], &[], &[0.0]);
        let solid = rate_solid(&f, &thermal, cx);
        let saved = solid.solution().clone();
        let accepted = solid
            .advance(
                &thermal,
                &(),
                0.1,
                &ROLLERS,
                &[(3, 5e6 * V / L)],
                IncrementSettings::default(),
                heat_steps(0.1),
                config(),
                cx,
            )
            .unwrap();
        let mut short = config();
        short.maximum_iterations = 1;
        assert!(
            solid
                .advance(
                    &thermal,
                    &(),
                    0.1,
                    &ROLLERS,
                    &[(3, 5e6 * V / L)],
                    IncrementSettings::default(),
                    heat_steps(0.1),
                    short,
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.solution(), &saved);
        assert_eq!(solid.time_s(), 0.0);
        assert_eq!(thermal.time_s(), 0.0);
        let retry = solid
            .advance(
                &thermal,
                &(),
                0.1,
                &ROLLERS,
                &[(3, 5e6 * V / L)],
                IncrementSettings::default(),
                heat_steps(0.1),
                config(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.solid().solution(), accepted.solid().solution());
        assert_eq!(retry.thermal().cells(), accepted.thermal().cells());
    });
}
