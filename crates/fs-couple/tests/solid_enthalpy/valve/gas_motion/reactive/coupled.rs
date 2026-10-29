//! Joint endpoint chemistry, native elastic work and conservative continuation.
//! Equal-mass synthetic isomers give an independent scalar elimination oracle.
#[path = "coupled/thermal.rs"]
mod thermal;

use super::*;
use fs_couple::cavity_phs::modal::reactive::coupled::CoupledReactiveModalGasSolve;
use fs_couple::cavity_phs::piston::reactive::{
    boundary::GasReactionBoundary, implicit::coupled::ImplicitMechanismSolve,
};

const COUPLED_AREA: f64 = 0.004;
const OUTSIDE: f64 = 50_000.0;

fn close(actual: f64, expected: f64, absolute: f64) {
    assert!(
        actual.is_finite() && (actual - expected).abs() <= absolute,
        "{actual:e} != {expected:e}, absolute tolerance {absolute:e}"
    );
}

fn amount() -> f64 {
    EXTERIOR * VOLUME / (R * TEMPERATURE)
}

fn implicit(reaction: GasReaction) -> GasReaction {
    let bound = 0.2 * amount();
    let policy = ImplicitMechanismSolve {
        extent_brackets_mol: reaction
            .channels()
            .map(|(channel, _)| {
                (
                    channel.reverse_rate_closure().equilibrium().reaction().clone(),
                    [-bound, bound],
                )
            })
            .collect(),
        extent_tolerance_mol: 1e-14,
        interior_fraction: 1e-4,
        jacobian_relative_step: 1e-5,
        maximum_iterations: 48,
        maximum_evaluations: 4096,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-13,
            relative_tolerance: 0.0,
            linear_restart: 1,
            max_linear_cycles: 2,
            ..Default::default()
        },
    };
    reaction.with_implicit_mechanism(policy).unwrap()
}

fn solve() -> CoupledReactiveModalGasSolve {
    let mut mechanical = controls();
    mechanical.gas = chemistry().gas;
    mechanical.maximum_trials = 256;
    CoupledReactiveModalGasSolve {
        mechanical,
        pressure_interval_pa: [0.5 * EXTERIOR, 1.5 * EXTERIOR],
        pressure_tolerance_pa: 1e-6,
        pressure_difference_fraction: 1e-5,
    }
}

// For A <-> B with equal molar mass/Cv, n is conserved and V cancels from
// dt V k (p_A/P_ref - p_B/(K_eq P_ref)). This is an algebraic BE solution;
// no production reaction, caloric inversion or nonlinear solver is used.
fn endpoint_extent(temperature: f64, scale: f64) -> f64 {
    let beta = scale * R * temperature / (RATE as f64 * EXTERIOR);
    let reverse = (-RELEASE / (R * temperature)).exp();
    beta * amount() / (1.0 + beta * (1.0 + reverse))
}

fn bracketed_root(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..90 {
        let mid = 0.5 * (lo + hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    0.5 * (lo + hi)
}

struct JointReference {
    pressure: f64,
    temperature: f64,
    extent: f64,
    volume_change: f64,
    velocity: f64,
    entropy: f64,
}

// Independent six-tet cube integrals give M=0.6 kg and K=2E/3.
// From rest the original held-force response gives
// dV=A^2 (1-cos(omega dt)) (P-Pext)/K. Combining this with pV=nRT
// leaves one scalar first-law equation, Cv*n*dT-release*xi+p*dV-Q=0.
fn joint_reference(scale: f64, young: f64, heat: f64) -> JointReference {
    let mass = 0.6;
    let stiffness = 2.0 * young / 3.0;
    let omega = (stiffness / mass).sqrt();
    let dt = 1.0 / RATE as f64;
    let compliance = 2.0 * COUPLED_AREA.powi(2) * (0.5 * omega * dt).sin().powi(2)
        / stiffness;
    let sample = |pressure: f64| {
        let volume_change = compliance * (pressure - OUTSIDE);
        let temperature = pressure * (VOLUME + volume_change) / (amount() * R);
        let extent = endpoint_extent(temperature, scale);
        let residual = amount() * 2.5 * R * (temperature - TEMPERATURE)
            - RELEASE * extent
            + pressure * volume_change
            - heat;
        (residual, volume_change, temperature, extent)
    };
    let pressure = bracketed_root(0.5 * EXTERIOR, 1.5 * EXTERIOR, |p| sample(p).0);
    let (_, volume_change, temperature, extent) = sample(pressure);
    let x = extent / amount();
    let entropy = amount() * 2.5 * R * (temperature / TEMPERATURE).ln()
        + amount() * R * ((VOLUME + volume_change) / VOLUME).ln()
        - amount() * R * (x * x.ln() + (1.0 - x) * (1.0 - x).ln());
    JointReference {
        pressure,
        temperature,
        extent,
        volume_change,
        velocity: COUPLED_AREA * (pressure - OUTSIDE) * (omega * dt).sin()
            / (mass * omega),
        entropy,
    }
}

#[test]
fn explicit_boundary_uses_final_volume_heat_work_and_original_zero_exchange_path() {
    let scale = 3000.0;
    let (reaction, gas) = reacting_gas(scale, [1.0, 0.0, 0.0]);
    let reaction = implicit(reaction);
    let original = gas.clone();
    let boundary = GasReactionBoundary {
        final_volume_m3: 1.2 * VOLUME,
        heat_j: 1.0,
        work_by_gas_j: 1.0,
    };
    let next = reaction
        .boundary_implicit_candidate(
            &gas, 1.0 / RATE as f64, boundary, chemistry(), || false,
        )
        .unwrap();
    let temperature = bracketed_root(TEMPERATURE, 800.0, |t| {
        amount() * 2.5 * R * (t - TEMPERATURE) - RELEASE * endpoint_extent(t, scale)
    });
    let extent = endpoint_extent(temperature, scale);
    close(next.gas().state().state().temperature, temperature, 2e-7);
    close(
        next.gas().state().state().pressure,
        amount() * R * temperature / boundary.final_volume_m3,
        1e-5,
    );
    close(next.step().reaction.extent_mol, extent, 2e-13);
    assert_eq!(next.gas().volume_m3(), boundary.final_volume_m3);
    assert_eq!(next.gas().mass_kg(), gas.mass_kg());
    close(next.step().reaction.stored_energy_change_j, 0.0, 1e-9);
    close(next.step().energy_residual_j, 0.0, 1e-9);
    close(next.step().endpoint_heat_entropy_j_k, 1.0 / temperature, 1e-12);
    assert!(next.step().entropy_residual_j_k > 0.0);
    assert_eq!(gas, original);

    let boundary = GasReactionBoundary {
        final_volume_m3: gas.volume_m3(), heat_j: 0.0, work_by_gas_j: 0.0,
    };
    let endpoint = reaction
        .boundary_implicit_candidate(
            &gas, 1.0 / RATE as f64, boundary, chemistry(), || false,
        )
        .unwrap();
    let native = reaction
        .parcel_candidate(&gas, 1.0 / RATE as f64, chemistry(), || false)
        .unwrap();
    assert_eq!(endpoint.gas(), native.gas());
    assert_eq!(&endpoint.step().reaction, native.step());
}

#[test]
fn source_rate_and_actual_material_stiffness_share_one_pressure_chemistry_endpoint() {
    with_cx(|cx| {
        let mut extents = Vec::new();
        for (scale, young) in [(1000.0, 1000.0), (3000.0, 1000.0), (3000.0, 2000.0)] {
            let (reaction, gas) = reacting_gas(scale, [1.0, 0.0, 0.0]);
            let reaction = implicit(reaction);
            let solid = motion(young, 0.0, 0.0, COUPLED_AREA, RATE, cx);
            let identity = solid.valve().identity();
            let mut run = GasValveMotion::new(gas, solid, cx).unwrap();
            let report = run
                .step_coupled_reactive(&reaction, 1.0, OUTSIDE, chemistry(), solve(), cx)
                .unwrap();
            let expected = joint_reference(scale, young, 1.0);
            let physical = run.motion().nodal_motion(cx).unwrap();
            close(report.held_pressure_pa, expected.pressure, 1e-5);
            close(run.gas().state().state().temperature, expected.temperature, 2e-7);
            close(report.reaction.reaction.extent_mol, expected.extent, 2e-13);
            close(
                physical.displacement_m[7][0],
                expected.volume_change / COUPLED_AREA,
                1e-10,
            );
            close(physical.velocity_m_s[7][0], expected.velocity, 1e-9);
            close(
                report.reaction.boundary.work_by_gas_j,
                expected.pressure * expected.volume_change,
                1e-9,
            );
            close(report.reaction.reaction.entropy_change_j_k, expected.entropy, 1e-8);
            close(
                report.reaction.entropy_residual_j_k,
                expected.entropy - 1.0 / expected.temperature,
                1e-8,
            );
            assert!(report.reaction.reaction.rate_material.is_some());
            assert!(report.reaction.reaction.implicit_mechanism.is_some());
            assert!(report.pressure_residual_pa.abs() <= solve().pressure_tolerance_pa);
            assert!(
                report.total_energy_residual_j.abs() <= solve().mechanical.energy_tolerance_j
            );
            assert_eq!(run.motion().valve().identity(), identity);
            assert_eq!(run.motion().accepted_steps(), 1);
            extents.push(report.reaction.reaction.extent_mol);
        }
        assert!(extents[1] > 2.5 * extents[0]);
        // Stiffness enters the same chemical endpoint through actual volume/work.
        // This deliberately small, independently predicted effect is not amplified.
        assert!((extents[2] - extents[1]).abs() > 1e-13);
    });
}

#[test]
fn continuing_reaction_retains_mesh_anchor_source_history_and_full_energy_account() {
    with_cx(|cx| {
        let (reaction, gas) = reacting_gas(1000.0, [1.0, 0.0, 0.0]);
        let reaction = implicit(reaction);
        let initial_u = gas.state().specific_internal_energy_j_kg();
        let mass = gas.mass_kg();
        let solid = motion(1000.0, 0.0, 0.05, COUPLED_AREA, RATE, cx);
        let identity = solid.valve().identity();
        let mut run = GasValveMotion::new(gas, solid, cx).unwrap();
        let mut loss = 0.0;
        let mut previous_species: Option<Vec<f64>> = None;
        for tick in 1..=6 {
            let report = run
                .step_coupled_reactive(&reaction, 0.2, OUTSIDE, chemistry(), solve(), cx)
                .unwrap();
            let chemical = &report.reaction.reaction;
            if let Some(previous) = &previous_species {
                for (species, amount) in chemical.species.iter().zip(previous) {
                    close(species.initial_mol, *amount, 1e-13);
                }
            }
            previous_species = Some(chemical.species.iter().map(|s| s.final_mol).collect());
            let physical = run.motion().nodal_motion(cx).unwrap();
            close(
                run.gas().volume_m3() - VOLUME,
                COUPLED_AREA * physical.displacement_m[7][0],
                solve().mechanical.volume_tolerance_m3,
            );
            loss += report.motion.modal.viscous_dissipation_j;
            let measured = mass
                * (run.gas().state().specific_internal_energy_j_kg() - initial_u)
                + energy(run.motion())
                + loss
                + OUTSIDE * (run.gas().volume_m3() - VOLUME)
                - tick as f64 * 0.2;
            close(measured, 0.0, tick as f64 * solve().mechanical.energy_tolerance_j);
            close(
                chemical.rate_material.as_ref().unwrap()
                    .property(ACTIVITY_RATE_SCALE_PROPERTY).unwrap().value_si(),
                1000.0,
                1e-10,
            );
            assert_eq!(run.gas().mass_kg(), mass);
            assert_eq!(run.motion().valve().identity(), identity);
            assert_eq!(run.motion().accepted_steps(), tick);
            close(report.motion.time_s, tick as f64 / RATE as f64, 1e-15);
            assert!(report.reaction.entropy_residual_j_k >= -chemistry().entropy_tolerance_j_k);
        }
        assert!(loss > 0.0);
    });
}

#[test]
fn failed_source_work_budget_and_cancellation_preserve_deterministic_retry() {
    let (reaction, gas) = reacting_gas(3000.0, [1.0, 0.0, 0.0]);
    let reaction = implicit(reaction);
    let saved = gas.clone();
    let boundary = GasReactionBoundary {
        final_volume_m3: 1.1 * VOLUME, heat_j: 1.0, work_by_gas_j: 0.1,
    };
    let mut polls = 0;
    let cancelled = reaction.boundary_implicit_candidate(
        &gas, 1.0 / RATE as f64, boundary, chemistry(), || {
            polls += 1;
            polls >= 12
        },
    );
    assert!(matches!(cancelled, Err(GasReactionError::Cancelled)));
    assert_eq!(gas, saved);
    assert!(polls >= 12);
    assert!(
        reaction.boundary_implicit_candidate(
            &gas, 1.0 / RATE as f64,
            GasReactionBoundary { heat_j: 1e12, ..boundary },
            chemistry(), || false,
        ).is_err()
    );
    assert_eq!(gas, saved);

    with_cx(|cx| {
        let mut run = GasValveMotion::new(
            gas.clone(), motion(1000.0, 0.0, 0.0, COUPLED_AREA, RATE, cx), cx,
        ).unwrap();
        let original = run.clone();
        let unchanged = |run: &GasValveMotion| {
            assert_eq!(run.gas(), original.gas());
            assert_eq!(run.motion().model().states(), original.motion().model().states());
            assert_eq!(run.motion().accepted_steps(), original.motion().accepted_steps());
            assert_eq!(run.motion().valve().identity(), original.motion().valve().identity());
        };
        let mut budget = solve();
        budget.mechanical.maximum_trials = 1;
        assert!(matches!(
            run.step_coupled_reactive(&reaction, 1.0, OUTSIDE, chemistry(), budget, cx),
            Err(ReactiveModalGasError::Motion(ModalGasError::TrialBudget))
        ));
        unchanged(&run);
        let mut interval = solve();
        interval.pressure_interval_pa = [EXTERIOR - 1.0, EXTERIOR + 1.0];
        assert!(
            run.step_coupled_reactive(&reaction, 1.0, OUTSIDE, chemistry(), interval, cx)
                .is_err()
        );
        unchanged(&run);
        let mut mismatched = solve();
        mismatched.mechanical.gas.energy_tolerance_j_kg *= 2.0;
        assert!(
            run.step_coupled_reactive(&reaction, 1.0, OUTSIDE, chemistry(), mismatched, cx)
                .is_err()
        );
        unchanged(&run);
        assert!(
            run.step_coupled_reactive(&reaction, 1e12, OUTSIDE, chemistry(), solve(), cx)
                .is_err()
        );
        unchanged(&run);
        let preview = run
            .coupled_reactive_step_candidate(&reaction, 1.0, OUTSIDE, chemistry(), solve(), cx)
            .unwrap();
        unchanged(&run);
        let report = run
            .step_coupled_reactive(&reaction, 1.0, OUTSIDE, chemistry(), solve(), cx)
            .unwrap();
        assert_eq!(run.gas(), preview.motion().gas());
        assert_eq!(run.motion().model().states(), preview.motion().motion().model().states());
        assert_eq!(report.reaction, preview.step().reaction);
        assert_eq!(run.motion().accepted_steps(), 1);
    });
}
