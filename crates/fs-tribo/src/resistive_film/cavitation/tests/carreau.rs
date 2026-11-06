use super::*;
use fs_material::fluid_transport::rheology::slit::LaminarCarreauSlit;
use fs_material::fluid_transport::rheology::{
    CarreauYasuda, ThermalCarreauYasuda, ViscosityShiftMode,
};
use fs_material::visco::temperature::{TemperatureShift, TemperatureShiftModel};

const MU: f64 = 0.04;
const TIME: f64 = 0.8;

fn medium() -> CavitationLiquid {
    CavitationLiquid {
        ambient_pressure_pa: 2100.0,
        ..liquid()
    }
}

fn source_controls() -> CarreauCavitationControls {
    let mut c = CarreauCavitationControls::default();
    c.passage.local.relative_constitutive_tolerance = 1e-13;
    c.passage.local.relative_quadrature_tolerance = 2e-11;
    c.passage.relative_wall_velocity_tolerance = 1e-11;
    c
}

fn mass_controls() -> CavitationSolveConfig {
    CavitationSolveConfig {
        maximum_iterations: 128,
        absolute_mass_tolerance_kg: 1e-20,
        relative_mass_tolerance: 2e-11,
    }
}

fn passage(
    edge: &FilmChannel,
    temperature: f64,
    constant: bool,
    ceiling: f64,
) -> LaminarCarreauSlit {
    let law = ThermalCarreauYasuda::new(
        CarreauYasuda::new(
            MU,
            if constant { MU } else { 0.0 },
            TIME,
            0.5,
            1.0,
            [0.0, 1e8],
        )
        .unwrap(),
        TemperatureShift::new(
            TemperatureShiftModel::Arrhenius {
                activation_energy_j_mol: 12_000.0,
            },
            300.0,
            [280.0, 360.0],
        )
        .unwrap(),
        ViscosityShiftMode::Vertical,
    );
    LaminarCarreauSlit::new(
        law,
        temperature,
        edge.width_m,
        edge.gap.reference_m.max(1e-3),
        edge.length_m,
        ceiling,
    )
    .unwrap()
}

fn passages(base: &ResistiveFilm, temperature: f64, constant: bool) -> Vec<LaminarCarreauSlit> {
    base.channels
        .iter()
        .map(|edge| passage(edge, temperature, constant, 1e7))
        .collect()
}

fn film(base: ResistiveFilm, constant: bool) -> CavitatingCarreauFilm {
    let sources = passages(&base, 300.0, constant);
    CavitatingCarreauFilm::new(base, sources, medium(), mass_controls(), source_controls()).unwrap()
}

#[track_caller]
fn near(a: f64, b: f64, absolute: f64, relative: f64) {
    assert!(a.is_finite() && b.is_finite());
    assert!(
        (a - b).abs() <= absolute + relative * a.abs().max(b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}

// Independent physical-gap oracle: for a=1, n=1/2, mu_infinity=0,
// tau=mu_zero*gamma/sqrt(1+lambda*abs(gamma)) has an analytic stress inverse.
// Original-source routines, production quadrature and graph derivatives are
// deliberately absent. Integrate momentum directly after bisection of the
// wall-speed condition; this remains a mathematical model, not calibration.
fn rate(stress: f64, viscosity: f64) -> f64 {
    let b = stress.abs() / viscosity;
    let x = TIME * b;
    (0.5 * b * (x + (x * x + 4.0).sqrt())).copysign(stress)
}

fn integral(h: f64, f: impl Fn(f64) -> f64) -> f64 {
    let n = 1024;
    let dy = h / f64::from(n);
    let mut value = f(0.0) + f(h);
    for i in 1..n {
        value += (if i % 2 == 0 { 2.0 } else { 4.0 }) * f(f64::from(i) * dy);
    }
    value * dy / 3.0
}

fn pressure_flow(edge: &FilmChannel, h: f64, dp: f64, relative: f64, mu: f64) -> f64 {
    if h <= 0.0 {
        return 0.0;
    }
    let bound = mu * (1.0 + relative.abs() / h) + dp.abs() * h / edge.length_m;
    let (mut low, mut high) = (-bound, bound);
    for _ in 0..70 {
        let c = 0.5 * (low + high);
        let slip = integral(h, |y| rate(c + dp * (0.5 * h - y) / edge.length_m, mu));
        if slip < relative {
            low = c;
        } else {
            high = c;
        }
    }
    let c = 0.5 * (low + high);
    edge.width_m
        * integral(h, |y| {
            (0.5 * h - y) * rate(c + dp * (0.5 * h - y) / edge.length_m, mu)
        })
}

#[test]
fn g3_constant_source_recovers_native_rupture_refill_supply_and_reversal() {
    let base = geometry(false);
    let supply = [Some(CavitationReservoir {
        pressure_pa: 2050.0,
        liquid_fraction: 1.0,
    })];
    let native = CavitatingFilm::new(base.with_viscosity(MU).unwrap(), medium(), mass_controls())
        .unwrap()
        .with_reservoirs(&supply)
        .unwrap();
    let carreau = film(base, true).with_reservoirs(&supply).unwrap();
    let mut a = native.initial_state(&[0.0]).unwrap();
    let mut b = carreau.initial_state(&[0.0]).unwrap();
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
        let actual = carreau
            .advance(&b, &[q], &[entrainment], &[relative], 0.1, || false)
            .unwrap();
        near(
            actual.state.pressure_pa()[0],
            expected.state.pressure_pa()[0],
            1e-8,
            2e-10,
        );
        near(
            actual.state.liquid_fraction()[0],
            expected.state.liquid_fraction()[0],
            2e-11,
            2e-10,
        );
        near(
            actual.channel_volume_flow_m3_s[0],
            expected.channel_volume_flow_m3_s[0],
            1e-18,
            2e-9,
        );
        near(
            actual.generalized_force[0],
            expected.generalized_force[0],
            1e-11,
            2e-9,
        );
        near(
            actual.report.liquid_mass_after_kg,
            expected.report.liquid_mass_after_kg,
            1e-16,
            2e-10,
        );
        assert_eq!(actual.state.configuration(), expected.state.configuration());
        a = expected.state;
        b = actual.state;
    }
}

#[test]
fn g1_original_carreau_rupture_and_refill_match_independent_momentum_and_mass() {
    let base = geometry(false);
    let original = base.clone();
    let model = film(base, false);
    let state = model.initial_state(&[0.0]).unwrap();
    let relative = 0.003;
    let rupture = model
        .advance(&state, &[-0.001], &[0.0], &[relative], 0.1, || false)
        .unwrap();
    let q = pressure_flow(&original.channels[0], 0.001, -100.0, relative, MU);
    near(rupture.channel_volume_flow_m3_s[0], q, 1e-20, 2e-8);
    near(
        rupture.state.liquid_fraction()[0],
        (1e-6 - 0.1 * q) / 2e-6,
        1e-12,
        1e-10,
    );
    assert_eq!(rupture.state.pressure_pa(), &[2000.0]);
    assert!(rupture.state.liquid_fraction()[0] < 1.0);
    let refill = model
        .advance(&rupture.state, &[0.0001], &[0.0], &[relative], 0.1, || {
            false
        })
        .unwrap();
    assert_eq!(refill.state.liquid_fraction(), &[1.0]);
    assert!(refill.state.pressure_pa()[0] > medium().ambient_pressure_pa);
    let q = pressure_flow(
        &original.channels[0],
        0.001,
        refill.state.pressure_pa()[0] - 2100.0,
        relative,
        MU,
    );
    near(refill.channel_volume_flow_m3_s[0], q, 1e-20, 2e-8);
    near(
        refill.report.liquid_mass_before_kg,
        rupture.report.liquid_mass_after_kg,
        0.0,
        0.0,
    );
    near(
        refill.report.liquid_mass_after_kg + refill.report.boundary_mass_out_kg,
        refill.report.liquid_mass_before_kg,
        1e-16,
        1e-11,
    );
}

#[test]
fn g1_manufactured_full_film_retains_combined_sliding_pressure_and_actual_aperture() {
    let mut base = geometry(false);
    let other = base.channels[0].clone();
    // Use the public constructor so topology, coefficients and identities all
    // correspond to this two-supply graph, rather than editing prepared fields.
    base = ResistiveFilm::new(
        base.cells.to_vec(),
        vec![
            FilmChannel {
                gap: gap(0.001, 0.2),
                ..other.clone()
            },
            FilmChannel {
                gap: gap(0.001, 0.4),
                ..other
            },
        ],
        1,
        MU,
        base.limits,
    )
    .unwrap();
    let supply = [2050.0, 2300.0];
    let target_pressure = 2500.0;
    let relative = [0.05, -0.018];
    let entrainment = [-1e-8, 1.5e-8];
    let q_next = 1e-5;
    let expected: Vec<_> = base
        .channels
        .iter()
        .enumerate()
        .map(|(k, edge)| {
            let h = edge.gap.reference_m - edge.gap.closure[0] * q_next;
            pressure_flow(edge, h, target_pressure - supply[k], relative[k], MU) + entrainment[k]
        })
        .collect();
    let dt = base.cells[0].area_m2 * q_next / expected.iter().sum::<f64>();
    assert!(dt > 0.0);
    let model = film(base, false)
        .with_reservoirs(&supply.map(|pressure_pa| {
            Some(CavitationReservoir {
                pressure_pa,
                liquid_fraction: 1.0,
            })
        }))
        .unwrap();
    let old = model.initial_state(&[0.0]).unwrap();
    let step = model
        .advance(&old, &[q_next], &entrainment, &relative, dt, || false)
        .unwrap();
    assert_eq!(step.state.liquid_fraction(), &[1.0]);
    near(step.state.pressure_pa()[0], target_pressure, 1e-7, 2e-8);
    for (actual, expected) in step.channel_volume_flow_m3_s.iter().zip(expected) {
        near(*actual, expected, 1e-18, 3e-8);
    }
    near(
        step.generalized_force[0],
        0.001 * (target_pressure - 2100.0),
        1e-10,
        2e-8,
    );
}

#[test]
fn g1_nonlinear_multicell_reformation_and_reversal_retain_shared_face_mass() {
    let cells = vec![
        FilmCell {
            area_m2: 0.001,
            gap: gap(0.001, 0.0)
        };
        3
    ];
    let edge = |from, to| FilmChannel {
        from,
        to,
        width_m: 0.0012,
        length_m: 0.1,
        gap: gap(0.001, 0.0),
    };
    let base = ResistiveFilm::new(
        cells,
        vec![
            edge(0, None),
            edge(0, Some(1)),
            edge(1, Some(2)),
            edge(2, None),
        ],
        1,
        MU,
        geometry(false).limits,
    )
    .unwrap();
    let supply = [
        Some(CavitationReservoir {
            pressure_pa: 2000.0,
            liquid_fraction: 0.25,
        }),
        None,
        None,
        Some(CavitationReservoir {
            pressure_pa: 2300.0,
            liquid_fraction: 1.0,
        }),
    ];
    let model = film(base.clone(), false).with_reservoirs(&supply).unwrap();
    let mut old = model.initialize(&[0.0], &[0.2, 0.7, 1.0]).unwrap();
    let relative = [0.02, -0.01, 0.015, -0.005];
    let dt = 2.0;
    let mut saw_partial = false;
    let mut saw_pressure = false;
    for direction in [1.0, 1.0, -1.0, -1.0] {
        let entrainment = [-2e-7, 2e-7, 2e-7, 2e-7].map(|e| direction * e);
        let step = model
            .advance(&old, &[0.0], &entrainment, &relative, dt, || false)
            .unwrap();
        let mut outward = [0.0; 3];
        for (k, edge) in base.channels.iter().enumerate() {
            let to_p = edge.to.map_or_else(
                || supply[k].unwrap().pressure_pa,
                |j| step.state.pressure_pa()[j],
            );
            let donor = if entrainment[k] >= 0.0 {
                step.state.liquid_fraction()[edge.from]
            } else {
                edge.to.map_or_else(
                    || supply[k].unwrap().liquid_fraction,
                    |j| step.state.liquid_fraction()[j],
                )
            };
            let expected = pressure_flow(
                edge,
                0.001,
                step.state.pressure_pa()[edge.from] - to_p,
                relative[k],
                MU,
            ) + entrainment[k] * donor;
            near(step.channel_volume_flow_m3_s[k], expected, 1e-17, 2e-8);
            outward[edge.from] += expected;
            if let Some(j) = edge.to {
                outward[j] -= expected;
            }
        }
        for (i, outflow) in outward.iter().enumerate() {
            let theta = step.state.liquid_fraction()[i];
            saw_partial |= theta < 1.0;
            saw_pressure |= step.state.pressure_pa()[i] > 2000.0;
            near(
                1e-6 * (theta - old.liquid_fraction()[i]) + dt * outflow,
                0.0,
                2e-14,
                0.0,
            );
            assert!(theta >= 0.0 && theta <= 1.0);
            if theta < 1.0 {
                assert_eq!(step.state.pressure_pa()[i], 2000.0);
            }
        }
        near(
            step.report.liquid_mass_after_kg + step.report.boundary_mass_out_kg,
            step.report.liquid_mass_before_kg,
            1e-13,
            1e-10,
        );
        old = step.state;
    }
    assert!(saw_partial && saw_pressure);
}

#[test]
fn g3_source_refresh_and_endpoint_density_keep_original_inventory_and_identity() {
    let base = geometry(true);
    let native = CavitatingFilm::new(base.clone(), medium(), mass_controls()).unwrap();
    let state = native.initialize(&[0.0], &[0.4]).unwrap();
    let model = native
        .with_carreau(passages(&base, 300.0, false), source_controls())
        .unwrap();
    let warm = model.with_passages(passages(&base, 330.0, false)).unwrap();
    let same = warm
        .advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    assert_eq!(same.state, state);
    let next_liquid = CavitationLiquid {
        density_kg_m3: 800.0,
        cavitation_pressure_pa: 2050.0,
        ..medium()
    };
    let continued = warm
        .advance_with_liquid(
            &state,
            &[0.0],
            &[0.0],
            &[0.0],
            0.1,
            next_liquid,
            passages(&base, 340.0, false),
            || false,
        )
        .unwrap();
    near(continued.step.state.liquid_fraction()[0], 0.5, 1e-15, 1e-15);
    assert_eq!(continued.step.state.pressure_pa(), &[2050.0]);
    assert_eq!(
        continued.step.report.liquid_mass_before_kg,
        continued.step.report.liquid_mass_after_kg
    );
    assert!(
        warm.advance(&continued.step.state, &[0.0], &[0.0], &[0.0], 0.1, || false)
            .is_err()
    );
    assert!(
        continued
            .film
            .advance(&state, &[0.0], &[0.0], &[0.0], 0.1, || false)
            .is_err()
    );
    continued
        .film
        .advance(&continued.step.state, &[0.0], &[0.0], &[0.0], 0.1, || false)
        .unwrap();

    let open = geometry(false);
    let cold = film(open.clone(), false);
    let warm = cold.with_passages(passages(&open, 340.0, false)).unwrap();
    let old = cold.initial_state(&[0.0]).unwrap();
    let a = cold
        .advance(&old, &[-0.001], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    let b = warm
        .advance(&old, &[-0.001], &[0.0], &[0.0], 0.1, || false)
        .unwrap();
    assert!(b.channel_volume_flow_m3_s[0].abs() > a.channel_volume_flow_m3_s[0].abs());
    assert!(b.state.liquid_fraction()[0] > a.state.liquid_fraction()[0]);
    assert_eq!(
        b.report.liquid_mass_before_kg,
        a.report.liquid_mass_before_kg
    );
}

#[test]
fn g4_constant_and_nonlinear_source_support_cannot_be_bypassed() {
    let base = geometry(false);
    for constant in [true, false] {
        let model = CavitatingCarreauFilm::new(
            base.clone(),
            vec![passage(&base.channels[0], 300.0, constant, 10.0)],
            medium(),
            mass_controls(),
            source_controls(),
        )
        .unwrap();
        let state = model.initial_state(&[0.0]).unwrap();
        assert!(
            model
                .advance(&state, &[0.0], &[0.0], &[0.05], 0.1, || false)
                .is_err()
        );
        assert!(
            model
                .advance(&state, &[0.0005], &[0.0], &[0.0], 0.001, || false)
                .is_err()
        );
        assert_eq!(state.configuration(), &[0.0]);
        assert_eq!(state.liquid_fraction(), &[1.0]);
        assert_eq!(state.pressure_pa(), &[2100.0]);
    }
    let closed = film(geometry(true), false);
    let state = closed.initialize(&[0.0], &[0.5]).unwrap();
    assert!(
        closed
            .advance(&state, &[0.0], &[0.0], &[1e-3], 0.1, || false)
            .is_err()
    );
    assert!(
        closed
            .advance(&state, &[0.0], &[1e-10], &[0.0], 0.1, || false)
            .is_err()
    );
    let mut controls = source_controls();
    controls.maximum_line_search_steps = 0;
    assert!(
        CavitatingCarreauFilm::new(
            base.clone(),
            passages(&base, 300.0, false),
            medium(),
            mass_controls(),
            controls
        )
        .is_err()
    );
    let barus = base.with_barus_viscosity(1e-8).unwrap();
    assert!(
        CavitatingCarreauFilm::new(
            barus,
            passages(&base, 300.0, false),
            medium(),
            mass_controls(),
            source_controls()
        )
        .is_err()
    );
}

#[test]
fn g4_late_source_cancellation_and_budget_refusal_preserve_exact_retry() {
    let base = geometry(false);
    let model = film(base.clone(), false);
    let old = model.initial_state(&[0.0]).unwrap();
    let mut checks = 0;
    let good = model
        .advance(&old, &[0.0001], &[0.0], &[0.04], 0.1, || {
            checks += 1;
            false
        })
        .unwrap();
    assert!(checks > 10);
    let mut count = 0;
    assert!(
        model
            .advance(&old, &[0.0001], &[0.0], &[0.04], 0.1, || {
                count += 1;
                count == checks
            })
            .is_err()
    );
    let retry = model
        .advance(&old, &[0.0001], &[0.0], &[0.04], 0.1, || false)
        .unwrap();
    assert_eq!(good.state, retry.state);
    assert_eq!(good.report, retry.report);
    assert_eq!(
        good.channel_volume_flow_m3_s,
        retry.channel_volume_flow_m3_s
    );
    assert_eq!(good.generalized_force, retry.generalized_force);
    assert_eq!(old.configuration(), &[0.0]);
    let short = CavitatingCarreauFilm::new(
        base.clone(),
        passages(&base, 300.0, false),
        medium(),
        CavitationSolveConfig {
            maximum_iterations: 1,
            ..mass_controls()
        },
        source_controls(),
    )
    .unwrap();
    let initial = short.initial_state(&[0.0]).unwrap();
    assert!(
        short
            .advance(&initial, &[0.0001], &[0.0], &[0.04], 0.1, || false)
            .is_err()
    );
}
