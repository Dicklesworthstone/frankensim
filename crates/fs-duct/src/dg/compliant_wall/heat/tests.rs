//! Physical references for solid-only heat on a continuing native trajectory.
use super::*;
use fs_material::gas::{GasSpec, GasState};
use fs_material::visco::GeneralizedMaxwell;

const RADIUS: f64 = 0.01;
const LENGTH: f64 = 0.3;
const DT: f64 = 2e-6;

fn budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn duct(inlet: DgBoundary, outlet: DgBoundary) -> DgDuct {
    let gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap();
    DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: RADIUS,
                length: LENGTH,
            }],
        },
        &gas,
        DgMesh {
            elements: 2,
            order: 1,
        },
        inlet,
        outlet,
    )
    .unwrap()
}

fn pin(resistance: f64) -> WallPin {
    WallPin {
        surface_density: 0.2,
        stiffness_per_area: 1e7,
        resistance,
    }
}

fn close(a: f64, b: f64, relative: f64) {
    assert!(
        (a - b).abs() <= relative * a.abs().max(b.abs()).max(1e-30),
        "{a:e} != {b:e}"
    );
}

#[test]
fn g1_warm_kelvin_heat_matches_the_dimensional_sealed_wall_oscillator() {
    let material = pin(300.0);
    let dg = duct(DgBoundary::Closed, DgBoundary::Closed)
        .with_compliant_wall(material)
        .unwrap();
    let propagator = dg.propagator(DT, 1e-13).unwrap();
    let mut state = dg.project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..24 {
        propagator.step(&mut state, &[]).unwrap();
    }
    let before = state.clone();
    let channel = dg
        .prepare_compliant_wall_heat(&propagator, budget(), &mut || false)
        .unwrap();
    assert_eq!(channel.operator_identity(), dg.operator_identity());
    assert_eq!(channel.duration(), DT);
    let actual = channel
        .integrated_heat_j(&state, &[], &mut || false)
        .unwrap();

    // Independent dimensional solution: sigma*w'' + r*w' + (K+2B/a)*w = p0.
    let beta = 2.0 * dg.gas().density * dg.gas().sound_speed.powi(2) / RADIUS;
    let decay = material.resistance / (2.0 * material.surface_density);
    let frequency =
        ((material.stiffness_per_area + beta) / material.surface_density - decay * decay).sqrt();
    let (points, weights) = gauss_legendre_rule(16);
    let expected = points
        .iter()
        .zip(&weights)
        .map(|(&point, &weight)| {
            let t = 24.0 * DT + 0.5 * DT * (point + 1.0);
            let velocity = 10.0 / (material.surface_density * frequency)
                * (-decay * t).exp()
                * (frequency * t).sin();
            0.5 * DT
                * weight
                * core::f64::consts::TAU
                * RADIUS
                * LENGTH
                * material.resistance
                * velocity
                * velocity
        })
        .sum::<f64>();
    assert!(actual > 0.0);
    close(actual, expected, 2e-9);
    assert_eq!(state, before);

    let half_step = dg.propagator(DT / 2.0, 1e-13).unwrap();
    let half = dg
        .prepare_compliant_wall_heat(&half_step, budget(), &mut || false)
        .unwrap();
    assert_eq!(half.duration(), DT / 2.0);
    let first = half.integrated_heat_j(&state, &[], &mut || false).unwrap();
    let mut middle = state.clone();
    half_step.step(&mut middle, &[]).unwrap();
    let second = half.integrated_heat_j(&middle, &[], &mut || false).unwrap();
    close(first + second, actual, 2e-9);

    propagator.step(&mut state, &[]).unwrap();
    close(
        before.stored_energy_j() - state.stored_energy_j(),
        actual,
        2e-8,
    );
    assert_eq!(state.step_count(), before.step_count() + 1);
}

// Independent spatial integral of r*v^2 and z_j^2/(K_j*tau_j) on each original
// element. It uses dimensional physical fields, not the source-row coordinates.
fn physical_power(
    dg: &DgDuct,
    state: &DgState,
    pins: &[WallPin; 2],
    arms: &[(usize, f64, f64)],
) -> f64 {
    let motion = dg.compliant_wall_state(state).unwrap();
    let local = dg.mesh.order + 1;
    let (points, weights) = gauss_legendre_rule(3);
    let mut power = 0.0;
    for (e, element) in dg.elements.iter().enumerate() {
        for (&point, &weight) in points.iter().zip(&weights) {
            let shape = dg.lgl.eval_basis(point);
            let field = |values: &[f64]| -> f64 {
                (0..local).map(|i| shape[i] * values[e * local + i]).sum()
            };
            let area = core::f64::consts::TAU * RADIUS * element.length * 0.5 * weight;
            let velocity = field(&motion.velocity_m_s);
            power += area * pins[e].resistance * velocity * velocity;
            for (row, &(owner, stiffness, tau)) in arms.iter().enumerate() {
                if owner == e {
                    let pressure = field(&motion.relaxation_pressure_pa[row]);
                    power += area * pressure * pressure / (stiffness * tau);
                }
            }
        }
    }
    power
}

fn physical_heat(
    dg: &DgDuct,
    initial: &DgState,
    ports: &[PortInput],
    pins: &[WallPin; 2],
    arms: &[(usize, f64, f64)],
    panels: usize,
) -> f64 {
    let step = dg.propagator(DT / panels as f64, 1e-13).unwrap();
    let mut state = initial.clone();
    let mut integral = physical_power(dg, &state, pins, arms);
    for i in 1..=panels {
        step.step(&mut state, ports).unwrap();
        let weight = if i == panels {
            1.0
        } else if i % 2 == 0 {
            2.0
        } else {
            4.0
        };
        integral += weight * physical_power(dg, &state, pins, arms);
    }
    integral * DT / (3.0 * panels as f64)
}

#[test]
fn g3_regional_kelvin_and_maxwell_heat_excludes_gas_losses_on_the_same_driven_state() {
    let pins = [
        pin(300.0),
        WallPin {
            surface_density: 0.7,
            stiffness_per_area: 4e7,
            resistance: 80.0,
        },
    ];
    let spectra = [
        GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001), (0.0, 0.3)]).unwrap(),
        GeneralizedMaxwell::new(4e6, vec![(5e5, 0.002), (9e5, 0.0004)]).unwrap(),
    ];
    let arms = [(0, 2e6, 0.001), (1, 5e6, 0.002), (1, 9e6, 0.0004)];
    let dg = duct(DgBoundary::Port(PortId(7)), DgBoundary::Closed)
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap()
        .with_compliant_wall_relaxation_regions(&pins, &[Some(&spectra[0]), Some(&spectra[1])])
        .unwrap();
    let ports = [PortInput {
        id: PortId(7),
        effort: 1e-7,
    }];
    let propagator = dg.propagator(DT, 1e-13).unwrap();
    let mut state = dg.project_initial_state(|x| 10.0 + 7.0 * x / LENGTH, |_| 1e-8);
    for _ in 0..24 {
        propagator.step(&mut state, &ports).unwrap();
    }
    let before = state.clone();
    let motion = dg.compliant_wall_state(&state).unwrap();
    assert!(motion.relaxation_stored_energy_j > 0.0);
    close(
        physical_power(&dg, &state, &pins, &arms),
        motion.dissipated_power_w,
        2e-10,
    );
    assert!(dg.wall_dissipated_power_w(&state).unwrap() > motion.dissipated_power_w);
    let channel = dg
        .prepare_compliant_wall_heat(&propagator, budget(), &mut || false)
        .unwrap();
    let heat = channel
        .integrated_heat_j(&state, &ports, &mut || false)
        .unwrap();
    let coarse = physical_heat(&dg, &state, &ports, &pins, &arms, 128);
    let fine = physical_heat(&dg, &state, &ports, &pins, &arms, 256);
    assert!(heat > 0.0);
    close(coarse, fine, 2e-7);
    close(heat, fine, 2e-7);
    assert_eq!(
        channel
            .integrated_heat_j(&state, &ports, &mut || false)
            .unwrap(),
        heat
    );
    assert_eq!(state, before);
    propagator.step(&mut state, &ports).unwrap();
    let total_loss = state.cumulative_dissipated_j() - before.cumulative_dissipated_j();
    assert!(total_loss > heat);
    assert_eq!(state.step_count(), before.step_count() + 1);
}

#[test]
fn g4_zero_solid_loss_and_source_budget_cancellation_refusals_preserve_native_state() {
    let gas = duct(DgBoundary::Port(PortId(7)), DgBoundary::Closed)
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap();
    let gas_step = gas.propagator(DT, 1e-13).unwrap();
    assert!(
        gas.prepare_compliant_wall_heat(&gas_step, budget(), &mut || false)
            .is_err()
    );
    let dg = gas.with_compliant_wall(pin(0.0)).unwrap();
    let propagator = dg.propagator(DT, 1e-13).unwrap();
    let ports = [PortInput {
        id: PortId(7),
        effort: 1e-7,
    }];
    let mut state = dg.project_initial_state(|x| 10.0 + 7.0 * x / LENGTH, |_| 1e-8);
    for _ in 0..24 {
        propagator.step(&mut state, &ports).unwrap();
    }
    let before = state.clone();
    let channel = dg
        .prepare_compliant_wall_heat(&propagator, budget(), &mut || false)
        .unwrap();
    assert!(dg.wall_dissipated_power_w(&state).unwrap() > 0.0);
    assert_eq!(
        channel
            .integrated_heat_j(&state, &ports, &mut || false)
            .unwrap(),
        0.0
    );
    let foreign = duct(DgBoundary::Port(PortId(7)), DgBoundary::Closed)
        .with_compliant_wall(pin(300.0))
        .unwrap();
    let foreign_step = foreign.propagator(DT, 1e-13).unwrap();
    assert!(
        dg.prepare_compliant_wall_heat(&foreign_step, budget(), &mut || false)
            .is_err()
    );
    assert!(
        channel
            .integrated_heat_j(&foreign.zero_state(), &ports, &mut || false)
            .is_err()
    );
    assert!(
        channel
            .integrated_heat_j(
                &state,
                &[PortInput {
                    id: PortId(9),
                    effort: 1e-7
                }],
                &mut || false
            )
            .is_err()
    );
    for allowance in [
        DenseExponentialBudget {
            max_dimension: 2 * dg.total_dofs() - 1,
            ..budget()
        },
        DenseExponentialBudget {
            max_dense_bytes: 0,
            ..budget()
        },
        DenseExponentialBudget {
            max_multiply_terms: 0,
            ..budget()
        },
    ] {
        assert!(
            dg.prepare_compliant_wall_heat(&propagator, allowance, &mut || false)
                .is_err()
        );
    }
    assert!(matches!(
        dg.prepare_compliant_wall_heat(&propagator, budget(), &mut || true),
        Err(DuctError::Cancelled { .. })
    ));

    // Exercise cancellation in a nonzero source contraction as well as prepare.
    let mut warm = foreign.project_initial_state(|_| 10.0, |_| 0.0);
    for _ in 0..24 {
        foreign_step.step(&mut warm, &ports).unwrap();
    }
    let warm_before = warm.clone();
    let active = foreign
        .prepare_compliant_wall_heat(&foreign_step, budget(), &mut || false)
        .unwrap();
    let mut polls = 0;
    let accepted = active
        .integrated_heat_j(&warm, &ports, &mut || {
            polls += 1;
            false
        })
        .unwrap();
    assert!(accepted > 0.0 && polls > 1);
    let mut visited = 0;
    assert!(matches!(
        active.integrated_heat_j(&warm, &ports, &mut || {
            visited += 1;
            visited == polls
        }),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(visited, polls);
    assert_eq!(
        active
            .integrated_heat_j(&warm, &ports, &mut || false)
            .unwrap(),
        accepted
    );
    assert_eq!(warm, warm_before);
    assert_eq!(state, before);
}
