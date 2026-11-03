//! G1 analytical transport and G0 conservation checks on the production
//! collision, streaming and source-substep functions.

use super::*;
use fs_math::det::{cos, ln, sin, sqrt};

fn wave<const D: usize>(shape: [usize; D], modes: [usize; D]) -> (Vec<[f64; 2]>, [f64; D]) {
    let k = core::array::from_fn(|axis| {
        core::f64::consts::TAU * modes[axis] as f64 / shape[axis] as f64
    });
    let basis = (0..shape.iter().product())
        .map(|mut cell| {
            let mut phase = 0.0;
            for axis in 0..D {
                phase += k[axis] * (cell % shape[axis]) as f64;
                cell /= shape[axis];
            }
            [cos(phase), sin(phase)]
        })
        .collect();
    (basis, k)
}

fn project<const Q: usize>(populations: &[[f64; Q]], basis: &[[f64; 2]]) -> [f64; 2] {
    let mut mode = [0.0; 2];
    for (pop, [c, s]) in populations.iter().zip(basis) {
        let temperature = pop.iter().sum::<f64>() - 1.0;
        mode[0] += temperature * c;
        mode[1] += temperature * s;
    }
    mode.map(|value| 2.0 * value / populations.len() as f64)
}

#[allow(clippy::too_many_arguments)]
fn transport_wave<const D: usize, const Q: usize>(
    lattice: &AxisLattice<D, Q>,
    shape: [usize; D],
    modes: [usize; D],
    tau: f64,
    mean_velocity: [f64; D],
    oscillation: [f64; D],
    frequency: f64,
    local_rates: bool,
) -> (f64, f64) {
    let (basis, k) = wave(shape, modes);
    let cells = basis.len();
    let k2: f64 = k.iter().map(|x| x * x).sum();
    let velocity = |time: usize| {
        core::array::from_fn(|axis| {
            mean_velocity[axis] + oscillation[axis] * sin(frequency * time as f64)
        })
    };
    let mut populations: Vec<_> = basis
        .iter()
        .map(|[c, _]| lattice.equilibrium(1.0 + 0.1 * c, velocity(0)))
        .collect();
    let mut history = FluxHistory::default();
    let mut post = vec![0.0; Q * cells];
    let rates = vec![tau; cells];
    let warmup = 24;
    let steps = 256;
    let mut initial = [0.0; 2];
    for time in 0..warmup + steps {
        if time == warmup {
            initial = project(&populations, &basis);
        }
        lattice.collide(
            &populations,
            &mut post,
            0..cells,
            if local_rates { 1.8 } else { tau },
            local_rates.then_some(rates.as_slice()),
            &mut history,
            |_| velocity(time),
        );
        lattice.stream(&mut populations, &post, shape, None);
    }
    let final_mode = project(&populations, &basis);
    let initial_norm = initial[0] * initial[0] + initial[1] * initial[1];
    let ratio = [
        (final_mode[0] * initial[0] + final_mode[1] * initial[1]) / initial_norm,
        (final_mode[1] * initial[0] - final_mode[0] * initial[1]) / initial_norm,
    ];
    let amplitude = sqrt(ratio[0] * ratio[0] + ratio[1] * ratio[1]);
    let measured_alpha = -ln(amplitude) / (k2 * steps as f64);
    let mut phase = 0.0;
    for axis in 0..D {
        let oscillating_displacement = if frequency == 0.0 {
            0.0
        } else {
            oscillation[axis]
                * (cos(frequency * warmup as f64) - cos(frequency * (warmup + steps) as f64))
                / frequency
        };
        phase += k[axis] * (mean_velocity[axis] * steps as f64 + oscillating_displacement);
    }
    let phase_sine_error = (ratio[1] * cos(phase) - ratio[0] * sin(phase)) / amplitude;
    assert!(ratio[0] * cos(phase) + ratio[1] * sin(phase) > 0.0);
    let mean = populations.iter().flatten().sum::<f64>() / cells as f64;
    assert!((mean - 1.0).abs() < 5e-12, "periodic heat drift: {mean}");
    (
        measured_alpha / (lattice.cs2 * (tau - 0.5)) - 1.0,
        phase_sine_error,
    )
}

#[test]
fn g1_fourier_diffusion_is_independent_of_axial_and_oblique_advection() {
    for tau in [0.65, 1.0, 1.4] {
        for (diffusion_error, phase_error) in [
            transport_wave(
                &D2Q5,
                [128, 1],
                [1, 0],
                tau,
                [0.14, 0.0],
                [0.0; 2],
                0.0,
                false,
            ),
            transport_wave(
                &D2Q5,
                [96, 96],
                [1, 1],
                tau,
                [0.10, 0.08],
                [0.0; 2],
                0.0,
                true,
            ),
            transport_wave(
                &D3Q7,
                [128, 1, 1],
                [1, 0, 0],
                tau,
                [0.14, 0.0, 0.0],
                [0.0; 3],
                0.0,
                true,
            ),
            transport_wave(
                &D3Q7,
                [96, 1, 96],
                [1, 0, 1],
                tau,
                [0.10, 0.0, 0.08],
                [0.0; 3],
                0.0,
                false,
            ),
        ] {
            assert!(
                diffusion_error.abs() < 0.007,
                "tau={tau}: diffusion error {diffusion_error}"
            );
            assert!(
                phase_error.abs() < 0.004,
                "tau={tau}: phase error {phase_error}"
            );
        }
    }
}

#[test]
fn g1_time_varying_velocity_has_integrated_phase_and_molecular_diffusion() {
    let frequency = core::f64::consts::TAU / 512.0;
    for tau in [0.7, 1.0, 1.3] {
        for (diffusion_error, phase_error) in [
            transport_wave(
                &D2Q5,
                [128, 1],
                [1, 0],
                tau,
                [0.06, 0.0],
                [0.07, 0.0],
                frequency,
                true,
            ),
            transport_wave(
                &D3Q7,
                [128, 1, 1],
                [1, 0, 0],
                tau,
                [0.06, 0.0, 0.0],
                [0.07, 0.0, 0.0],
                frequency,
                false,
            ),
        ] {
            assert!(
                diffusion_error.abs() < 0.007,
                "tau={tau}: diffusion error {diffusion_error}"
            );
            // The finite spatial wavelength and backward derivative retain
            // truncation error; bound accumulated phase to one milliradian.
            assert!(
                phase_error.abs() < 0.001,
                "tau={tau}: phase error {phase_error}"
            );
        }
    }
}

fn source_maintained_wave<const D: usize, const Q: usize>(
    lattice: &AxisLattice<D, Q>,
    shape: [usize; D],
    modes: [usize; D],
    tau: f64,
    velocity: [f64; D],
) -> f64 {
    let (basis, k) = wave(shape, modes);
    let cells = basis.len();
    let alpha = lattice.cs2 * (tau - 0.5);
    let k2: f64 = k.iter().map(|x| x * x).sum();
    let ku: f64 = k.iter().zip(velocity).map(|(k, u)| k * u).sum();
    let source_mode = [0.1 * alpha * k2, -0.1 * ku];
    let heat: Vec<_> = basis
        .iter()
        .map(|[c, s]| source_mode[0] * c + source_mode[1] * s)
        .collect();
    let mut populations: Vec<_> = basis
        .iter()
        .map(|[c, _]| lattice.equilibrium(1.0 + 0.1 * c, velocity))
        .collect();
    let mut history = FluxHistory::default();
    let mut post = vec![0.0; Q * cells];
    for _ in 0..16000 {
        lattice.source_half_kick(&mut populations, &heat, 0..cells, &mut history, |_| {
            velocity
        });
        lattice.collide(
            &populations,
            &mut post,
            0..cells,
            tau,
            None,
            &mut history,
            |_| velocity,
        );
        lattice.stream(&mut populations, &post, shape, None);
        lattice.source_half_kick(&mut populations, &heat, 0..cells, &mut history, |_| {
            velocity
        });
    }
    let temperature = project(&populations, &basis);
    let norm = temperature[0] * temperature[0] + temperature[1] * temperature[1];
    // The real part of Q/T measures diffusivity independently of the
    // dominant advective source, which would hide this error in a bulk norm.
    let measured_alpha =
        (source_mode[0] * temperature[0] + source_mode[1] * temperature[1]) / (norm * k2);
    measured_alpha / alpha - 1.0
}

#[test]
fn g1_prescribed_heating_preserves_transport_diffusivity_away_from_tau_one() {
    for tau in [0.65, 0.85, 1.4] {
        for error in [
            source_maintained_wave(&D2Q5, [64, 1], [1, 0], tau, [0.15, 0.0]),
            source_maintained_wave(&D3Q7, [64, 1, 1], [1, 0, 0], tau, [0.15, 0.0, 0.0]),
        ] {
            assert!(
                error.abs() < 0.012,
                "tau={tau}: source-maintained diffusion error {error}"
            );
        }
    }
}

fn source_rebind<const D: usize, const Q: usize>(lattice: &AxisLattice<D, Q>, u: [f64; D]) {
    let before = lattice.equilibrium(0.9, u);
    let mut populations = vec![before];
    // A retained transport derivative and a nonzero kinetic perturbation.
    populations[0][0] += 0.01;
    populations[0][lattice.positive[0]] -= 0.01;
    let previous: [f64; D] = core::array::from_fn(|axis| 0.85 * u[axis]);
    let mut history = FluxHistory {
        flux: vec![previous],
        ready: true,
    };
    let saved_history = history.clone();
    let old_pop = populations[0];
    lattice.source_half_kick(&mut populations, &[0.12], 0..1, &mut history, |_| u);
    let after_equilibrium = lattice.equilibrium(0.96, u);
    for q in 0..Q {
        let before_nonequilibrium = old_pop[q] - before[q];
        let after_nonequilibrium = populations[0][q] - after_equilibrium[q];
        assert!((before_nonequilibrium - after_nonequilibrium).abs() < 2e-16);
    }
    for axis in 0..D {
        assert!(
            ((0.96 * u[axis] - history.flux[0][axis]) - (0.9 * u[axis] - previous[axis])).abs()
                < 2e-16
        );
    }
    assert_eq!(saved_history.flux, vec![previous]);
    assert!(saved_history.ready);
    let saved = populations.clone();
    let before_history = history.flux.clone();
    lattice.source_half_kick(&mut populations, &[-0.0], 0..1, &mut history, |_| {
        panic!("zero source sampled velocity")
    });
    assert_eq!(populations, saved);
    assert_eq!(history.flux, before_history);
}

#[test]
fn g0_source_kicks_preserve_nonequilibrium_history_and_zero_source_bits() {
    source_rebind(&D2Q5, [0.11, -0.04]);
    source_rebind(&D3Q7, [0.11, -0.04, 0.03]);
}

fn wall_balance<const D: usize, const Q: usize>(lattice: &AxisLattice<D, Q>, shape: [usize; D]) {
    let cells: usize = shape.iter().product();
    let plane: usize = shape[..D - 1].iter().product();
    let active = plane..cells - plane;
    let mut populations: Vec<_> = (0..cells)
        .map(|cell| lattice.equilibrium(0.4 + 0.2 * (cell % 7) as f64 / 7.0, [0.0; D]))
        .collect();
    let heat: Vec<_> = (0..cells)
        .map(|cell| {
            if active.contains(&cell) {
                0.0003 * ((cell % 5) as f64 - 2.0)
            } else {
                0.0
            }
        })
        .collect();
    let rates: Vec<_> = (0..cells)
        .map(|cell| 0.65 + 0.1 * (cell % 7) as f64)
        .collect();
    let mut history = FluxHistory::default();
    let mut post = vec![0.0; Q * cells];
    let velocity = |_| core::array::from_fn(|axis| if axis == 0 { 0.07 } else { 0.0 });
    for step in 0..12 {
        let before = populations.clone();
        let source: f64 = heat[active.clone()].iter().sum();
        lattice.source_half_kick(
            &mut populations,
            &heat,
            active.clone(),
            &mut history,
            velocity,
        );
        lattice.collide(
            &populations,
            &mut post,
            active.clone(),
            1.7,
            Some(&rates),
            &mut history,
            velocity,
        );
        lattice.stream(&mut populations, &post, shape, Some([0.9, 0.2]));
        let positive = lattice.positive[D - 1];
        let negative = lattice.negative[D - 1];
        let bottom: f64 = populations[plane..2 * plane]
            .iter()
            .map(|p| 2.0 * (p[positive] - lattice.weights[positive] * 0.9))
            .sum();
        let top: f64 = populations[cells - 2 * plane..cells - plane]
            .iter()
            .map(|p| 2.0 * (p[negative] - lattice.weights[negative] * 0.2))
            .sum();
        lattice.source_half_kick(
            &mut populations,
            &heat,
            active.clone(),
            &mut history,
            velocity,
        );
        let stored: f64 = populations[active.clone()]
            .iter()
            .zip(&before[active.clone()])
            .map(|(a, b)| a.iter().zip(b).map(|(a, b)| a - b).sum::<f64>())
            .sum();
        assert!(
            (stored - source - bottom - top).abs() < 2e-13,
            "step {step}: heat defect {}",
            stored - source - bottom - top
        );
    }
}

#[test]
fn g0_source_and_wall_heat_balance_survives_local_relaxation() {
    wall_balance(&D2Q5, [7, 8]);
    wall_balance(&D3Q7, [5, 3, 8]);
}

#[test]
fn g1_resting_conduction_and_startup_remain_bounded() {
    // Resting Fourier diffusion is unchanged; separately check the
    // documented zero-correction startup and an explicit history restart.
    for tau in [0.65, 1.0, 1.4] {
        let (error, phase) =
            transport_wave(&D2Q5, [128, 1], [1, 0], tau, [0.0; 2], [0.0; 2], 0.0, false);
        assert!(error.abs() < 0.004);
        assert!(phase.abs() < 1e-12);
    }
    let mut history = FluxHistory::<2>::default();
    let populations = vec![D2Q5.equilibrium(0.9, [0.1, -0.02])];
    let mut post = vec![0.0; 5];
    D2Q5.collide(
        &populations,
        &mut post,
        0..1,
        0.8,
        None,
        &mut history,
        |_| [0.1, -0.02],
    );
    for (actual, expected) in post.iter().zip(populations[0]) {
        assert!((actual - expected).abs() < 2e-16);
    }
    history.reset();
    D2Q5.collide(
        &populations,
        &mut post,
        0..1,
        0.8,
        None,
        &mut history,
        |_| [-0.1, 0.02],
    );
    assert!(post.iter().all(|x| x.is_finite()));
    assert!((post.iter().sum::<f64>() - 0.9).abs() < 2e-16);
}
