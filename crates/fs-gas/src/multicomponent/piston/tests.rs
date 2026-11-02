use super::*;
use crate::multicomponent::tests::{gas, source};
use crate::multicomponent::{Mesh1d, MixtureGas};
use fs_qty::{Length, Temperature};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

const AREA: f64 = 0.7;
const CLOSED: Boundaries = Boundaries::Ends {
    left: Boundary::Reflecting,
    right: Boundary::Reflecting,
};

fn mixture(references: [f64; 3]) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.4, references[0]),
        source("B", 2.0 * R, 4.5, 0.2, references[1]),
        source("absent", 0.5 * R, 3.0, 0.1, references[2]),
    ])
}

fn field(count: usize, references: [f64; 3]) -> MixtureEuler1d {
    let gas = mixture(references);
    let states = (0..count)
        .map(|i| {
            let first = if i < count / 2 { 0.8 } else { 0.2 };
            gas.conserved(&[first, 1.0 - first, 0.0], 0.0, 1.0).unwrap()
        })
        .collect();
    MixtureEuler1d::new(
        gas,
        Length::new(1.0 / count as f64),
        AREA,
        states,
        CLOSED,
        || false,
    )
    .unwrap()
}

fn controls() -> PistonControl {
    PistonControl {
        midpoint_velocity_bracket_m_s: [-0.5, 0.5],
        momentum_tolerance_n_s: 1e-12,
        energy_tolerance_j: 1e-11,
        maximum_iterations: 64,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}"
    );
}

#[test]
fn g1_heterogeneous_inertial_compression_conserves_species_work_and_gauges() {
    let mut flow = field(32, [0.0; 3]);
    let mut shifted = field(32, [1e18, -7e17, 3e17]);
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    let mut other = piston;
    let before = flow.totals();
    let load = PistonLoad {
        exterior_pressure_pa: 1.8,
        external_force_n: -0.15,
    };
    let mut work = 0.0;
    let mut external_impulse = 0.0;
    let mut support = 0.0;
    let mut steps = 0;
    while flow.time().value() < 0.1 {
        let duration = Time::new(0.1 - flow.time().value());
        let report = flow
            .step_with_piston(
                &mut piston,
                duration,
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        let reference = shifted
            .step_with_piston(
                &mut other,
                duration,
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(report, reference);
        assert_eq!(piston, other);
        assert_eq!(flow.states, shifted.states);
        assert_eq!(flow.mesh(), shifted.mesh());
        assert!(report.admitted(controls()));
        for k in 0..3 {
            assert_eq!(report.gas.boundary_change.species_masses[k], 0.0);
            close(
                report.gas.conservation_residual.species_masses[k],
                0.0,
                3e-13,
            );
        }
        for q in flow.states.chunks_exact(5) {
            assert_eq!(
                q[2], 0.0,
                "pressure-driven motion must not seed an absent species"
            );
        }
        assert_eq!(report.gas.left_boundary_change.energy, 0.0);
        assert!(report.gas.right_boundary_change.energy > 0.0);
        assert!(report.displacement_m < 0.0);
        assert!(report.gas.maximum_cfl <= 0.40000000000001);
        close(
            report.gas.right_boundary_change.energy,
            -report.midpoint_velocity_m_s * report.gas_impulse_n_s,
            2e-14,
        );
        work += report.external_force_work_j - report.exterior_pressure_work_j;
        external_impulse +=
            report.gas.elapsed.value() * (load.external_force_n - AREA * load.exterior_pressure_pa);
        support += report.fixed_wall_impulse_n_s;
        steps += 1;
        assert!(steps < 500);
    }
    let after = flow.totals();
    for k in 0..3 {
        close(after.species_masses[k], before.species_masses[k], 3e-12);
    }
    let kinetic = 0.5 * piston.mass_kg() * piston.velocity_m_s().powi(2);
    close(after.energy - before.energy + kinetic, work, 3e-10);
    close(
        after.momentum - before.momentum + piston.mass_kg() * piston.velocity_m_s(),
        external_impulse - support,
        3e-10,
    );
    assert!(piston.velocity_m_s() < -0.05);
    let left = flow
        .gas()
        .primitive(flow.cell(0).unwrap(), || false)
        .unwrap();
    let right = flow
        .gas()
        .primitive(flow.cell(31).unwrap(), || false)
        .unwrap();
    assert!(right.pressure > 1.05);
    assert!(right.pressure > left.pressure + 0.03);
    assert_eq!(flow.mesh().positions()[0], 0.0);
    assert!(flow.mesh().positions()[32] < 1.0);
    close(
        flow.mesh().volumes().iter().sum(),
        AREA * flow.mesh().positions()[32],
        5e-14,
    );
}

#[test]
fn g0_single_species_inertia_reduces_to_the_shared_frozen_mechanics() {
    use crate::dynamics::{CaloricDomain, Euler1d, FrozenGas};
    use fs_thermochem::mixture::FrozenIdealGasMixtureModelV1;
    use fs_thermochem::{Composition, CompositionBasis};
    let model = source("A", R, 3.5, 0.4, 0.0);
    let mixed_gas = gas(vec![model.clone()]);
    let frozen_gas = FrozenGas::new(
        FrozenIdealGasMixtureModelV1::new(
            vec![model],
            Composition::new(CompositionBasis::MoleFraction, vec![1.0]).unwrap(),
        )
        .unwrap(),
        CaloricDomain::new(Temperature::new(0.001), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let mixed_q = mixed_gas.conserved(&[1.0], 0.0, 1.0).unwrap();
    let frozen_q = frozen_gas.conserved(1.0, 0.0, 1.0).unwrap();
    let mut mixed = MixtureEuler1d::new(
        mixed_gas,
        Length::new(0.05),
        AREA,
        vec![mixed_q; 20],
        CLOSED,
        || false,
    )
    .unwrap();
    let mut frozen = Euler1d::new(
        frozen_gas,
        Length::new(0.05),
        AREA,
        vec![frozen_q; 20],
        CLOSED,
        || false,
    )
    .unwrap();
    let mut mixed_piston = Piston::new(0.4, 0.0).unwrap();
    let mut frozen_piston = mixed_piston;
    let load = PistonLoad {
        exterior_pressure_pa: 1.6,
        external_force_n: -0.1,
    };
    for _ in 0..12 {
        let a = mixed
            .step_with_piston(
                &mut mixed_piston,
                Time::new(0.0005),
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        let b = frozen
            .step_with_piston(
                &mut frozen_piston,
                Time::new(0.0005),
                StepControl::default(),
                load,
                controls(),
                || false,
            )
            .unwrap();
        assert_eq!(a.gas.elapsed, b.gas.elapsed);
        close(
            mixed_piston.velocity_m_s(),
            frozen_piston.velocity_m_s(),
            3e-10,
        );
        close(a.gas_impulse_n_s, b.gas_impulse_n_s, 3e-11);
        close(
            a.gas.boundary_change.energy,
            b.gas.boundary_change.energy,
            3e-11,
        );
        for (i, &q) in frozen.states().iter().enumerate() {
            let m = mixed.cell(i).unwrap();
            close(m.partial_densities[0], q.density, 5e-10);
            close(m.momentum, q.momentum, 5e-10);
            close(m.energy, q.energy, 5e-10);
            close(mixed.mesh().volumes()[i], frozen.mesh().volumes()[i], 3e-13);
        }
        for (&a, &b) in mixed
            .mesh()
            .positions()
            .iter()
            .zip(frozen.mesh().positions())
        {
            close(a, b, 3e-13);
        }
    }
}

#[test]
fn g0_final_inertial_mixture_trial_replays_prescribed_species_motion() {
    let mut coupled = field(16, [0.0; 3]);
    let mut prescribed = field(16, [0.0; 3]);
    let mut piston = Piston::new(0.2, 0.0).unwrap();
    let report = coupled
        .step_with_piston(
            &mut piston,
            Time::new(0.004),
            StepControl::default(),
            PistonLoad {
                exterior_pressure_pa: 2.0,
                external_force_n: 0.1,
            },
            controls(),
            || false,
        )
        .unwrap();
    assert!(report.gas_trials > 3);
    let replay = prescribed
        .step_with_moving_walls(
            report.gas.elapsed,
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            [0.0, report.midpoint_velocity_m_s],
            || false,
        )
        .unwrap();
    assert_eq!(report.gas.elapsed, replay.elapsed);
    assert_eq!(coupled.mesh(), prescribed.mesh());
    assert_eq!(coupled.time(), prescribed.time());
    for (&q, &reference) in coupled.states.iter().zip(&prescribed.states) {
        close(q, reference, 3e-12);
    }
    for k in 0..3 {
        assert_eq!(report.gas.boundary_change.species_masses[k], 0.0);
        close(
            report.gas.after.species_masses[k],
            replay.after.species_masses[k],
            3e-13,
        );
    }
    close(
        report.gas.right_boundary_change.momentum,
        replay.right_boundary_change.momentum,
        3e-13,
    );
    close(
        report.gas.right_boundary_change.energy,
        replay.right_boundary_change.energy,
        3e-13,
    );
}

#[test]
fn g4_native_mixture_refusals_and_final_cancellation_preserve_all_owners() {
    let duration = Time::new(0.003);
    let load = PistonLoad {
        exterior_pressure_pa: 1.8,
        external_force_n: 0.0,
    };
    let mut baseline = field(12, [0.0; 3]);
    let mut wall = Piston::new(0.4, 0.0).unwrap();
    let mut polls = 0;
    baseline
        .step_with_piston(
            &mut wall,
            duration,
            StepControl::default(),
            load,
            controls(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    for cutoff in [1, 8, polls / 2, polls] {
        let mut flow = field(12, [0.0; 3]);
        let before = flow.states.clone();
        let mesh = flow.mesh().clone();
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        let original = piston;
        let mut seen = 0;
        let result = flow.step_with_piston(
            &mut piston,
            duration,
            StepControl::default(),
            load,
            controls(),
            || {
                seen += 1;
                seen == cutoff
            },
        );
        assert_eq!(result, Err(PistonError::Cancelled));
        assert_eq!(
            seen, cutoff,
            "native cancellation must never become a retry"
        );
        assert_eq!(flow.states, before);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time(), Time::new(0.0));
        assert_eq!(piston, original);
    }
    for (bracket, source_refusal) in [([-1e-9, 1e-9], false), ([-100.0, 100.0], true)] {
        let mut flow = field(12, [0.0; 3]);
        let before = flow.states.clone();
        let mesh = flow.mesh().clone();
        let mut piston = Piston::new(0.4, 0.0).unwrap();
        let original = piston;
        let result = flow.step_with_piston(
            &mut piston,
            duration,
            StepControl {
                maximum_halvings: 0,
                ..StepControl::default()
            },
            load,
            PistonControl {
                midpoint_velocity_bracket_m_s: bracket,
                ..controls()
            },
            || false,
        );
        let Err(PistonError::StepRetries { last, .. }) = result else {
            panic!("{result:?}");
        };
        if source_refusal {
            assert!(matches!(*last, PistonError::MixtureGas(_)));
        } else {
            assert!(matches!(*last, PistonError::OutsideBracket { .. }));
        }
        assert_eq!(flow.states, before);
        assert_eq!(flow.mesh(), &mesh);
        assert_eq!(flow.time(), Time::new(0.0));
        assert_eq!(piston, original);
    }
    let gas = mixture([0.0; 3]);
    let q = gas.conserved(&[0.3, 0.7, 0.0], 0.0, 1.0).unwrap();
    let mesh = Mesh1d::new(vec![0.0, 0.5, 1.0], vec![1.0; 3], vec![0.4; 2]).unwrap();
    let mut invalid =
        MixtureEuler1d::new_on_mesh(gas, mesh.clone(), vec![q; 2], CLOSED, || false).unwrap();
    let mut piston = Piston::new(0.4, 0.0).unwrap();
    assert!(matches!(
        invalid.step_with_piston(
            &mut piston,
            duration,
            StepControl::default(),
            load,
            controls(),
            || false
        ),
        Err(PistonError::MixtureGas(Error::Geometry(_)))
    ));
    assert_eq!(invalid.mesh(), &mesh);
    assert_eq!(invalid.time(), Time::new(0.0));
    assert_eq!(piston.velocity_m_s(), 0.0);
}
