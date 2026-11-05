use super::*;
use crate::resistive_film::power_law::{
    PowerLawControls, PowerLawFilm, PressureDrivenPowerLawFilm, RegularizedPowerLaw,
    ShearedPowerLawChannelReport, ShearedPowerLawFilm,
};
use crate::resistive_film::shear::WallVelocityPorts;

fn controls() -> PowerLawControls {
    PowerLawControls {
        maximum_iterations: 96,
        maximum_line_search_steps: 48,
        absolute_flow_tolerance_m3_s: 1e-24,
        relative_flow_tolerance: 2e-12,
        absolute_power_tolerance_w: 1e-22,
        relative_power_tolerance: 3e-12,
    }
}
fn mass_controls() -> CavitationSolveConfig {
    CavitationSolveConfig {
        maximum_iterations: 128,
        absolute_mass_tolerance_kg: 1e-22,
        relative_mass_tolerance: 2e-13,
    }
}
fn law(n: f64) -> RegularizedPowerLaw {
    RegularizedPowerLaw {
        consistency_pa_s_n: 0.1 / 5.0_f64.powf(n),
        flow_index: n,
        minimum_shear_rate_per_s: 5.0,
    }
}
fn medium() -> CavitationLiquid {
    CavitationLiquid {
        ambient_pressure_pa: 2100.0,
        ..liquid()
    }
}
fn nonlinear(
    base: &ResistiveFilm,
    laws: Vec<RegularizedPowerLaw>,
    liquid: CavitationLiquid,
) -> CavitatingPowerLawFilm {
    CavitatingPowerLawFilm::new(
        PowerLawFilm::new(base, laws, controls()).unwrap(),
        liquid,
        mass_controls(),
    )
    .unwrap()
}
#[track_caller]
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

// Independent no-slip momentum integral in the physical half gap. This uses
// the authored stress law and composite Simpson quadrature, not production
// flux formulas, effective viscosity or the nonlinear complementarity solver.
fn pressure_flow(material: RegularizedPowerLaw, edge: &FilmChannel, h: f64, dp: f64) -> f64 {
    if h <= 0.0 || dp == 0.0 {
        return 0.0;
    }
    let mu = material.consistency_pa_s_n
        * material
            .minimum_shear_rate_per_s
            .powf(material.flow_index - 1.0);
    let core = mu * material.minimum_shear_rate_per_s;
    let gradient = dp.abs() / edge.length_m;
    let half = 0.5 * h;
    let transition = half.min(core / gradient);
    let integral = |a: f64, b: f64| {
        let count = 2048;
        let step = (b - a) / f64::from(count);
        let integrand = |y: f64| {
            let stress = gradient * y;
            let rate = if stress <= core {
                stress / mu
            } else {
                (stress / material.consistency_pa_s_n).powf(1.0 / material.flow_index)
            };
            y * rate
        };
        let mut sum = integrand(a) + integrand(b);
        for i in 1..count {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(a + f64::from(i) * step);
        }
        sum * step / 3.0
    };
    (2.0 * edge.width_m * (integral(0.0, transition) + integral(transition, half))).copysign(dp)
}

#[test]
fn g3_newtonian_limit_retains_native_cavitation_trajectory_and_inventory() {
    let base = geometry(false);
    let mu = 0.03;
    let material = RegularizedPowerLaw {
        consistency_pa_s_n: mu,
        flow_index: 1.0,
        minimum_shear_rate_per_s: 1.0,
    };
    let supplied = [Some(CavitationReservoir {
        pressure_pa: 2050.0,
        liquid_fraction: 1.0,
    })];
    let native = CavitatingFilm::new(base.with_viscosity(mu).unwrap(), medium(), mass_controls())
        .unwrap()
        .with_reservoirs(&supplied)
        .unwrap();
    let film = nonlinear(&base, vec![material], medium())
        .with_reservoirs(&supplied)
        .unwrap();
    let (mut a, mut b) = (
        native.initial_state(&[0.0]).unwrap(),
        film.initial_state(&[0.0]).unwrap(),
    );
    for (q, entrainment, relative) in [
        (-0.001, 1e-8, 0.003),
        (-0.0005, -1e-8, -0.003),
        (0.0, 0.0, 0.0),
        (0.0001, 1e-8, 0.01),
        (-0.0008, -1e-8, -0.01),
    ] {
        let expected = native
            .advance(&a, &[q], &[entrainment], 0.1, || false)
            .unwrap();
        let actual = film
            .advance(&b, &[q], &[entrainment], &[relative], 0.1, || false)
            .unwrap();
        assert_eq!(actual.state.configuration(), expected.state.configuration());
        assert_eq!(actual.state.pressure_pa(), expected.state.pressure_pa());
        assert_eq!(
            actual.state.liquid_fraction(),
            expected.state.liquid_fraction()
        );
        assert_eq!(
            actual.channel_volume_flow_m3_s,
            expected.channel_volume_flow_m3_s
        );
        assert_eq!(actual.generalized_force, expected.generalized_force);
        assert_eq!(actual.report, expected.report);
        a = expected.state;
        b = actual.state;
    }
}

#[test]
fn g1_full_liquid_recovers_the_combined_power_law_pressure_and_sliding_profile() {
    let cell_gap = GapPort {
        reference_m: 2e-4,
        closure: vec![1.0, 0.0, 0.0],
    };
    let base = ResistiveFilm::new(
        vec![FilmCell {
            area_m2: 1e-4,
            gap: cell_gap,
        }],
        vec![FilmChannel {
            from: 0,
            to: None,
            width_m: 0.01,
            length_m: 0.02,
            gap: GapPort {
                reference_m: 2e-4,
                closure: vec![0.0; 3],
            },
        }],
        3,
        0.02,
        FilmLimits {
            minimum_cell_gap_m: 1e-7,
            maximum_gap_m: 1e-2,
            maximum_pressure_pa: 1e7,
        },
    )
    .unwrap();
    let dt = 0.01;
    let q = [1e-7, 0.0, 0.0];
    let velocity = [q[0] / dt, -0.0004, 0.0006];
    let entrainment = 0.01 * 2e-4 * 0.5 * (velocity[1] + velocity[2]);
    for n in [0.5, 1.7, 3.0] {
        let rheology = PowerLawFilm::new(&base, vec![law(n)], controls()).unwrap();
        let full = ShearedPowerLawFilm::new(
            PressureDrivenPowerLawFilm::new(rheology.clone(), &[Some(0.0)]).unwrap(),
            vec![WallVelocityPorts {
                lower: vec![0.0, 1.0, 0.0],
                upper: vec![0.0, 0.0, 1.0],
            }],
        )
        .unwrap();
        let film = CavitatingPowerLawFilm::new(rheology, liquid(), mass_controls()).unwrap();
        let old = film.initial_state(&[0.0; 3]).unwrap();
        let step = film
            .advance(
                &old,
                &q,
                &[entrainment],
                &[velocity[2] - velocity[1]],
                dt,
                || false,
            )
            .unwrap();
        let (mut force, mut pressure, mut channels) =
            ([0.0; 3], [0.0], [ShearedPowerLawChannelReport::default()]);
        full.evaluate_into(&q, &velocity, &mut force, &mut pressure, &mut channels)
            .unwrap();
        assert_eq!(step.state.liquid_fraction(), &[1.0]);
        near(
            step.state.pressure_pa()[0] - liquid().ambient_pressure_pa,
            pressure[0],
            2e-5,
            2e-8,
        );
        near(step.generalized_force[0], force[0], 2e-9, 2e-8);
        assert_eq!(&step.generalized_force[1..], &[0.0, 0.0]);
        near(
            step.channel_volume_flow_m3_s[0],
            channels[0].volume_flow_m3_s,
            2e-18,
            2e-8,
        );
        near(
            step.report.liquid_mass_before_kg - step.report.liquid_mass_after_kg,
            step.report.boundary_mass_out_kg,
            3e-18,
            2e-10,
        );
    }
}

#[test]
fn g1_thinning_and_thickening_rupture_then_refill_match_independent_momentum_integrals() {
    let base = geometry(false);
    for n in [0.5, 1.7, 3.0] {
        let material = law(n);
        let film = nonlinear(&base, vec![material], medium());
        let old = film.initial_state(&[0.0]).unwrap();
        let inflow = pressure_flow(material, &base.channels[0], 0.001, 100.0);
        let dt = 0.1;
        let opened = film
            .advance(&old, &[-0.001], &[0.0], &[0.0], dt, || false)
            .unwrap();
        assert_eq!(opened.state.pressure_pa(), &[2000.0]);
        near(
            opened.state.liquid_fraction()[0],
            (1e-6 + dt * inflow) / 2e-6,
            2e-13,
            2e-11,
        );
        near(opened.channel_volume_flow_m3_s[0], -inflow, 1e-20, 3e-10);
        near(
            opened.report.boundary_mass_out_kg,
            -1000.0 * dt * inflow,
            1e-18,
            3e-10,
        );
        let refilled = film
            .advance(&opened.state, &[0.0], &[0.0], &[0.0], dt, || false)
            .unwrap();
        assert_eq!(refilled.state.liquid_fraction(), &[1.0]);
        near(refilled.state.pressure_pa()[0], 2200.0, 2e-5, 2e-10);
        near(refilled.channel_volume_flow_m3_s[0], inflow, 3e-18, 3e-8);
        assert_eq!(
            refilled.report.liquid_mass_before_kg,
            opened.report.liquid_mass_after_kg
        );
        near(
            refilled.report.liquid_mass_after_kg,
            old.liquid_fraction()[0] * 1e-3,
            1e-18,
            2e-13,
        );
    }
}

#[test]
fn g1_nonlinear_reformation_front_conserves_shared_face_mass_and_reverses() {
    let native = two_cells();
    let base = native.geometry();
    let laws = vec![law(0.5), law(1.7), law(0.8)];
    let supplies = [
        Some(CavitationReservoir {
            pressure_pa: 2100.0,
            liquid_fraction: 1.0,
        }),
        None,
        Some(CavitationReservoir {
            pressure_pa: 2000.0,
            liquid_fraction: 0.4,
        }),
    ];
    let forward = nonlinear(base, laws.clone(), native.liquid())
        .with_reservoirs(&supplies)
        .unwrap();
    let reversed = nonlinear(base, laws.iter().rev().copied().collect(), native.liquid())
        .with_reservoirs(&[supplies[2], None, supplies[0]])
        .unwrap();
    let (mut a, mut b) = (
        forward.initialize(&[0.0], &[0.2, 0.1]).unwrap(),
        reversed.initialize(&[0.0], &[0.1, 0.2]).unwrap(),
    );
    let e = [-1e-8, 1e-8, 1e-8];
    let dt = 30.0;
    let mut mixed_front = false;
    for _ in 0..6 {
        let first = forward
            .advance(&a, &[0.0], &e, &[0.0; 3], dt, || false)
            .unwrap();
        let second = reversed
            .advance(&b, &[0.0], &e.map(|v| -v), &[0.0; 3], dt, || false)
            .unwrap();
        for i in 0..2 {
            near(
                first.state.pressure_pa()[i],
                second.state.pressure_pa()[1 - i],
                2e-8,
                3e-11,
            );
            near(
                first.state.liquid_fraction()[i],
                second.state.liquid_fraction()[1 - i],
                2e-12,
                3e-11,
            );
        }
        let p = first.state.pressure_pa();
        let theta = first.state.liquid_fraction();
        let mut balance = [0.0; 2];
        for i in 0..2 {
            balance[i] = 1e-6 * (theta[i] - a.liquid_fraction()[i]);
        }
        for (k, edge) in base.channels.iter().enumerate() {
            let endpoint = edge
                .to
                .map_or(supplies[k].map_or(0.0, |r| r.pressure_pa), |j| p[j]);
            let donor = if e[k] >= 0.0 {
                theta[edge.from]
            } else {
                edge.to
                    .map_or(supplies[k].map_or(0.0, |r| r.liquid_fraction), |j| theta[j])
            };
            let expected =
                pressure_flow(laws[k], edge, 0.001, p[edge.from] - endpoint) + e[k] * donor;
            near(first.channel_volume_flow_m3_s[k], expected, 2e-18, 3e-9);
            balance[edge.from] += dt * expected;
            if let Some(j) = edge.to {
                balance[j] -= dt * expected;
            }
        }
        for residual in balance {
            assert!(
                residual.abs() < 2e-15,
                "independent cell mass defect {residual}"
            );
        }
        near(
            first.report.liquid_mass_after_kg - first.report.liquid_mass_before_kg,
            -first.report.boundary_mass_out_kg,
            3e-17,
            2e-10,
        );
        mixed_front |= theta[0] == 1.0 && theta[1] < 1.0 && p[0] > 2000.0 && p[1] == 2000.0;
        a = first.state;
        b = second.state;
    }
    assert!(
        mixed_front,
        "exercise a full pressurized cell adjoining retained cavitation"
    );
}

#[test]
fn g1_starved_refill_time_refines_to_the_nonlinear_pressure_feed_solution() {
    let base = geometry(false);
    let material = law(0.5);
    let film = nonlinear(&base, vec![material], medium());
    let inflow = pressure_flow(material, &base.channels[0], 0.001, 100.0);
    let entrainment = 0.5 * inflow;
    let volume = 1e-6;
    let exact_time = (volume / entrainment) * 1.8_f64.ln();
    for count in [8, 16, 32] {
        let dt = 2.0 * exact_time / f64::from(count);
        let mut state = film.initialize(&[0.0], &[0.2]).unwrap();
        let mut reformation = None;
        for index in 1..=count {
            let candidate = film
                .advance(&state, &[0.0], &[entrainment], &[0.0], dt, || false)
                .unwrap();
            if candidate.state.liquid_fraction()[0] == 1.0 {
                if reformation.is_none() {
                    reformation = Some(f64::from(index) * dt);
                }
            } else {
                let expected = (state.liquid_fraction()[0] + dt * inflow / volume)
                    / (1.0 + dt * entrainment / volume);
                near(candidate.state.liquid_fraction()[0], expected, 2e-12, 2e-11);
                assert_eq!(candidate.state.pressure_pa(), &[2000.0]);
            }
            state = candidate.state;
        }
        let numerical = reformation.expect("the supplied film must reform");
        assert!(numerical >= exact_time && numerical - exact_time < 2.0 * dt);
    }
}

#[test]
fn g3_viscosity_refresh_keeps_original_inventory_and_starved_transport() {
    let base = geometry(false);
    let original = nonlinear(&base, vec![law(0.5)], medium());
    let state = original.initialize(&[0.0], &[0.2]).unwrap();
    let mut changed = law(0.5);
    changed.consistency_pa_s_n *= 3.0;
    let next = original.with_laws(vec![changed]).unwrap();
    let step = next
        .advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    let expected = pressure_flow(changed, &base.channels[0], 0.001, 100.0);
    near(step.channel_volume_flow_m3_s[0], -expected, 1e-20, 2e-10);
    assert_eq!(
        step.report.liquid_mass_before_kg,
        medium().density_kg_m3 * (1e-6 * state.liquid_fraction()[0])
    );
    near(
        step.report.liquid_mass_after_kg,
        0.0002 + 1000.0 * 0.1 * expected,
        1e-17,
        2e-11,
    );
    let same = original.with_laws(original.laws().to_vec()).unwrap();
    assert_eq!(
        same.advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
            .unwrap(),
        original
            .advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
            .unwrap()
    );
    let foreign = nonlinear(&base, vec![law(0.5)], medium());
    assert!(
        foreign
            .advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
            .is_err()
    );
}

#[test]
fn g4_sealed_liquid_limits_solver_budgets_and_late_cancellation_are_atomic() {
    let sealed = nonlinear(&geometry(true), vec![law(1.7)], medium());
    let old = sealed.initial_state(&[0.0]).unwrap();
    let cavity = sealed
        .advance(&old, &[-0.001], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    assert_eq!(cavity.state.liquid_fraction(), &[0.5]);
    assert_eq!(
        cavity.report.liquid_mass_before_kg,
        cavity.report.liquid_mass_after_kg
    );
    let snapshot = cavity.state.clone();
    assert!(
        sealed
            .advance(&cavity.state, &[0.0001], &[0.0], &[0.0], 0.1, || false)
            .is_err()
    );
    assert!(
        sealed
            .advance(&cavity.state, &[-0.001], &[0.0], &[0.01], 0.1, || false)
            .is_err()
    );
    assert_eq!(cavity.state, snapshot);
    let base = geometry(false);
    let film = nonlinear(&base, vec![law(0.5)], medium());
    let state = film.initial_state(&[0.0]).unwrap();
    let before = state.clone();
    let mut polls = 0;
    film.advance(&state, &[1e-5], &[1e-8], &[0.002], 0.1, || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 5);
    for stop in [5, polls - 1, polls] {
        let mut count = 0;
        assert!(
            film.advance(&state, &[1e-5], &[1e-8], &[0.002], 0.1, || {
                count += 1;
                count == stop
            })
            .is_err()
        );
        assert_eq!(state, before);
    }
    let limited = CavitatingPowerLawFilm::new(
        PowerLawFilm::new(&base, vec![law(0.5)], controls()).unwrap(),
        medium(),
        CavitationSolveConfig {
            maximum_iterations: 1,
            ..mass_controls()
        },
    )
    .unwrap();
    let limited_state = limited.initial_state(&[0.0]).unwrap();
    // One accepted correction is within a one-correction budget. The local
    // stress remains in its physical Newtonian core, despite n != 1.
    let one_correction = limited
        .advance(&limited_state, &[1e-8], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    near(one_correction.state.pressure_pa()[0], 2102.0, 2e-8, 2e-12);
    assert!(
        limited
            .advance(&limited_state, &[1e-5], &[0.0], &[0.002], 0.1, || false)
            .is_err()
    );
    assert_eq!(limited_state.configuration(), &[0.0]);
    assert!(
        film.advance(&state, &[0.0], &[0.0], &[f64::NAN], 0.1, || false)
            .is_err()
    );
    assert_eq!(state, before);
}

#[test]
fn g3_divergence_free_large_entrainment_retains_small_squeeze_and_mass() {
    let original = geometry(false);
    let base = ResistiveFilm::new(
        original.cells.to_vec(),
        vec![original.channels[0].clone(); 2],
        1,
        0.02,
        original.limits,
    )
    .unwrap();
    let film = nonlinear(&base, vec![law(0.5); 2], medium());
    let old = film.initial_state(&[0.0]).unwrap();
    let reference = film
        .advance(&old, &[1e-8], &[0.0; 2], &[0.0; 2], 0.1, || false)
        .unwrap();
    // Conditioning case: the opposing, fully filled entrainment transfers
    // cancel in incidence before they meet the small pressure/storage defect.
    // The rounded total face flows cannot represent the pressure component.
    let translated = film
        .advance(&old, &[1e-8], &[1e8, -1e8], &[0.0; 2], 0.1, || false)
        .unwrap();
    near(
        translated.state.pressure_pa()[0],
        reference.state.pressure_pa()[0],
        2e-9,
        2e-12,
    );
    near(
        translated.generalized_force[0],
        reference.generalized_force[0],
        2e-12,
        2e-12,
    );
    assert_eq!(translated.state.liquid_fraction(), &[1.0]);
    near(
        translated.report.boundary_mass_out_kg,
        reference.report.boundary_mass_out_kg,
        1e-18,
        2e-11,
    );
    near(
        translated.report.liquid_mass_before_kg - translated.report.liquid_mass_after_kg,
        translated.report.boundary_mass_out_kg,
        1e-18,
        2e-11,
    );
}

#[test]
fn g1_branching_full_at_rupture_reforms_before_merit_globalization() {
    let template = geometry(false);
    let mut channels = Vec::new();
    let mut reservoirs = Vec::new();
    for to in 1..4 {
        let mut passage = template.channels[0].clone();
        passage.to = Some(to);
        channels.push(passage);
        reservoirs.push(None);
    }
    for from in 0..4 {
        let mut passage = template.channels[0].clone();
        passage.from = from;
        channels.push(passage);
        reservoirs.push(Some(CavitationReservoir {
            pressure_pa: 2100.0,
            liquid_fraction: 1.0,
        }));
    }
    let base = ResistiveFilm::new(
        vec![template.cells[0].clone(); 4],
        channels,
        1,
        0.02,
        template.limits,
    )
    .unwrap();
    let material = RegularizedPowerLaw {
        consistency_pa_s_n: 20.0,
        flow_index: 0.5,
        minimum_shear_rate_per_s: 1e6,
    };
    let film = nonlinear(
        &base,
        vec![material; 7],
        CavitationLiquid {
            ambient_pressure_pa: 2000.0,
            ..medium()
        },
    )
    .with_reservoirs(&reservoirs)
    .unwrap();
    let original = film.initial_state(&[0.0]).unwrap();
    assert_eq!(original.pressure_pa(), &[2000.0; 4]);
    assert_eq!(original.liquid_fraction(), &[1.0; 4]);
    // Every face is within the physical Newtonian core, but n != 1 selects
    // the nonlinear owner. At the initial full/rupture corner, a content
    // direction produces unequal pressure increments because the center has
    // degree four and each leaf degree two. Without phase relinearization,
    // the center mass defect increases at every positive line-search fraction.
    let supplied = film
        .advance(&original, &[0.0], &[0.0; 7], &[0.0; 7], 1.0, || false)
        .unwrap();
    assert_eq!(supplied.state.liquid_fraction(), &[1.0; 4]);
    for pressure in supplied.state.pressure_pa() {
        near(*pressure, 2100.0, 2e-9, 2e-12);
    }
    for flow in &supplied.channel_volume_flow_m3_s {
        assert!(flow.abs() < 1e-18);
    }
    near(supplied.report.boundary_mass_out_kg, 0.0, 1e-15, 0.0);
    for reservoir in reservoirs.iter_mut().flatten() {
        reservoir.pressure_pa = 2000.0;
    }
    let withdrawn = film.with_reservoirs(&reservoirs).unwrap();
    let released = withdrawn
        .advance(&supplied.state, &[0.0], &[0.0; 7], &[0.0; 7], 1.0, || false)
        .unwrap();
    for pressure in released.state.pressure_pa() {
        near(*pressure, 2000.0, 2e-9, 2e-12);
    }
    let separated = withdrawn
        .advance(&released.state, &[-1e-5], &[0.0; 7], &[0.0; 7], 1.0, || {
            false
        })
        .unwrap();
    assert_eq!(separated.state.pressure_pa(), &[2000.0; 4]);
    for theta in separated.state.liquid_fraction() {
        near(*theta, 0.001 / 0.00101, 2e-14, 2e-13);
    }
    assert_eq!(original.pressure_pa(), &[2000.0; 4]);
    assert_eq!(original.liquid_fraction(), &[1.0; 4]);
}

#[test]
fn g1_endpoint_density_and_nonlinear_viscosity_preserve_mass_and_physical_pressure() {
    let base = geometry(false);
    let inventory = CavitatingFilm::new(base.clone(), liquid(), mass_controls()).unwrap();
    let original = inventory.initial_state(&[0.0]).unwrap();
    let endpoint = CavitationLiquid {
        density_kg_m3: 800.0,
        cavitation_pressure_pa: 6000.0,
        ..liquid()
    };
    let dt = 0.1;
    let volume = 1e-6;
    let required_flow = (liquid().density_kg_m3 / endpoint.density_kg_m3 - 1.0) * volume / dt;
    for exponent in [0.5, 1.7] {
        let film = inventory
            .with_power_law(vec![law(exponent)], controls())
            .unwrap();
        let material = RegularizedPowerLaw {
            consistency_pa_s_n: 0.6 * law(exponent).consistency_pa_s_n,
            ..law(exponent)
        };
        let next = film
            .advance_with_liquid(
                &original,
                &[0.0],
                &[0.0],
                &[0.0],
                dt,
                endpoint,
                vec![material],
                || false,
            )
            .unwrap();
        // Independent no-slip momentum quadrature, followed by a scalar
        // bracket: this does not use the production flux or pressure root.
        let edge = &base.channels[0];
        let mut lo = 0.0;
        let mut hi = 1.0;
        while pressure_flow(material, edge, 0.001, hi) < required_flow {
            hi *= 2.0;
            assert!(hi < base.limits.maximum_pressure_pa);
        }
        for _ in 0..64 {
            let mid = 0.5 * (lo + hi);
            if pressure_flow(material, edge, 0.001, mid) < required_flow {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        let gauge = 0.5 * (lo + hi);
        near(
            next.step.state.pressure_pa()[0],
            endpoint.ambient_pressure_pa + gauge,
            1e-7,
            3e-10,
        );
        near(next.step.generalized_force[0], 0.001 * gauge, 1e-10, 3e-10);
        near(
            next.step.channel_volume_flow_m3_s[0],
            required_flow,
            1e-18,
            2e-12,
        );
        assert_eq!(next.step.state.liquid_fraction(), &[1.0]);
        assert_eq!(next.step.report.liquid_mass_before_kg, 1000.0 * volume);
        near(
            next.step.report.liquid_mass_after_kg,
            800.0 * volume,
            1e-18,
            1e-13,
        );
        near(
            next.step.report.boundary_mass_out_kg,
            200.0 * volume,
            1e-18,
            1e-12,
        );
        assert_eq!(next.film.laws(), &[material]);
        assert_eq!(next.film.liquid(), endpoint);
        assert_eq!(original.pressure_pa(), &[liquid().ambient_pressure_pa]);
    }
}

#[test]
fn g3_thermal_nonlinear_inventory_conversion_and_sealed_cycles_keep_original_mass() {
    let native = CavitatingFilm::new(geometry(true), liquid(), mass_controls()).unwrap();
    let original = native.initialize(&[0.0], &[0.25]).unwrap();
    let mut film = native.with_power_law(vec![law(0.7)], controls()).unwrap();
    // The constitutive conversion accepts the ORIGINAL state without reset.
    let same = film
        .advance(&original, &[0.0], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    assert_eq!(same.state, original);
    let mut state = original.clone();
    let mut accepted_mass = 1000.0 * (1e-6 * 0.25);
    for (density, rupture, scale) in [
        (800.0, 6000.0, 0.6),
        (1200.0, 1000.0, 1.4),
        (1000.0, 2000.0, 1.0),
    ] {
        let material = RegularizedPowerLaw {
            consistency_pa_s_n: scale * law(0.7).consistency_pa_s_n,
            ..law(0.7)
        };
        let next = film
            .advance_with_liquid(
                &state,
                &[0.0],
                &[0.0],
                &[0.0],
                0.1,
                CavitationLiquid {
                    density_kg_m3: density,
                    cavitation_pressure_pa: rupture,
                    ..liquid()
                },
                vec![material],
                || false,
            )
            .unwrap();
        assert_eq!(next.step.report.liquid_mass_before_kg, accepted_mass);
        near(
            next.step.report.liquid_mass_after_kg,
            accepted_mass,
            1e-18,
            2e-13,
        );
        assert_eq!(next.step.report.boundary_mass_out_kg, 0.0);
        near(
            next.step.state.liquid_fraction()[0],
            250.0 / density,
            1e-14,
            1e-13,
        );
        assert_eq!(next.step.state.pressure_pa(), &[rupture]);
        assert!(
            film.advance(&next.step.state, &[0.0], &[0.0], &[0.0], 0.1, || false)
                .is_err()
        );
        assert!(
            next.film
                .advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
                .is_err()
        );
        accepted_mass = next.step.report.liquid_mass_after_kg;
        state = next.step.state;
        film = next.film;
    }
    assert_eq!(original.liquid_fraction(), &[0.25]);
    let refreshed = film
        .advance_with_liquid(
            &state,
            &[0.0],
            &[0.0],
            &[0.0],
            0.1,
            film.liquid(),
            vec![law(1.4)],
            || false,
        )
        .unwrap();
    // A law-only update changes transport, so both material epochs retain the
    // original liquid identity. It does not erase or recreate filling.
    assert!(
        film.advance(&refreshed.step.state, &[0.0], &[0.0], &[0.0], 0.1, || false)
            .is_ok()
    );
    assert!(
        refreshed
            .film
            .advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
            .is_ok()
    );
}

#[test]
fn g3_newtonian_endpoint_continuation_is_the_exact_original_mass_transfer() {
    let native = CavitatingFilm::new(geometry(false), medium(), mass_controls()).unwrap();
    let state = native.initialize(&[0.0], &[0.3]).unwrap();
    // Switching from n!=1 to n=1 must also delegate by the ENDPOINT law.
    let film = native.with_power_law(vec![law(0.5)], controls()).unwrap();
    let endpoint = CavitationLiquid {
        density_kg_m3: 850.0,
        cavitation_pressure_pa: 1800.0,
        ..medium()
    };
    let material = RegularizedPowerLaw {
        consistency_pa_s_n: 0.015,
        flow_index: 1.0,
        minimum_shear_rate_per_s: 3.0,
    };
    let expected = native
        .advance_with_liquid(&state, &[-0.0002], &[1e-8], 0.1, endpoint, &[0.015], || {
            false
        })
        .unwrap();
    let actual = film
        .advance_with_liquid(
            &state,
            &[-0.0002],
            &[1e-8],
            &[0.004],
            0.1,
            endpoint,
            vec![material],
            || false,
        )
        .unwrap();
    assert_eq!(
        actual.step.state.configuration(),
        expected.step.state.configuration()
    );
    assert_eq!(
        actual.step.state.pressure_pa(),
        expected.step.state.pressure_pa()
    );
    assert_eq!(
        actual.step.state.liquid_fraction(),
        expected.step.state.liquid_fraction()
    );
    assert_eq!(
        actual.step.channel_volume_flow_m3_s,
        expected.step.channel_volume_flow_m3_s
    );
    assert_eq!(
        actual.step.generalized_force,
        expected.step.generalized_force
    );
    assert_eq!(actual.step.report, expected.step.report);
}

#[test]
fn g4_thermal_nonlinear_transition_refuses_sources_domains_and_late_cancellation_atomically() {
    let barus = CavitatingFilm::new(
        geometry(false).with_barus_viscosity(1e-6).unwrap(),
        medium(),
        mass_controls(),
    )
    .unwrap();
    assert!(barus.with_power_law(vec![law(0.5)], controls()).is_err());
    let native = CavitatingFilm::new(geometry(false), medium(), mass_controls()).unwrap();
    let film = native.with_power_law(vec![law(0.5)], controls()).unwrap();
    let state = native.initial_state(&[0.0]).unwrap();
    let saved = state.clone();
    let endpoint = CavitationLiquid {
        density_kg_m3: 800.0,
        ..medium()
    };
    let trial = |cancelled: &mut dyn FnMut() -> bool| {
        film.advance_with_liquid(
            &state,
            &[0.0],
            &[1e-8],
            &[0.002],
            0.1,
            endpoint,
            vec![law(0.7)],
            cancelled,
        )
    };
    let mut polls = 0;
    trial(&mut || {
        polls += 1;
        false
    })
    .unwrap();
    assert!(polls > 5);
    for stop in [1, 5, polls - 1, polls] {
        let mut count = 0;
        assert!(
            trial(&mut || {
                count += 1;
                count == stop
            })
            .is_err()
        );
        assert_eq!(state, saved);
        assert_eq!(film.liquid(), medium());
        assert_eq!(film.laws(), &[law(0.5)]);
    }
    let supplied = film
        .with_reservoirs(&[Some(CavitationReservoir {
            pressure_pa: 2050.0,
            liquid_fraction: 1.0,
        })])
        .unwrap();
    assert!(
        supplied
            .advance_with_liquid(
                &state,
                &[0.0],
                &[0.0],
                &[0.0],
                0.1,
                CavitationLiquid {
                    cavitation_pressure_pa: 2060.0,
                    ..endpoint
                },
                vec![law(0.7)],
                || false,
            )
            .is_err()
    );
    assert!(
        film.advance_with_liquid(
            &state,
            &[0.0],
            &[0.0],
            &[f64::NAN],
            0.1,
            endpoint,
            vec![law(0.7)],
            || false,
        )
        .is_err()
    );
    assert!(
        film.advance_with_liquid(
            &state,
            &[0.0],
            &[0.0],
            &[0.0],
            0.1,
            endpoint,
            vec![],
            || false,
        )
        .is_err()
    );
    assert_eq!(state, saved);
    // A rejected proposal leaves the accepted owner immediately retryable.
    film.advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
}
