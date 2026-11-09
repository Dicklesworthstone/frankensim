use super::*;
use crate::sgs3::dynamic::lagrangian::scalar::components::{
    LagrangianComponentsConfig, LagrangianDynamicComponents,
};
use crate::sgs3::dynamic::lagrangian::scalar::{LagrangianDynamicScalar, LagrangianScalarConfig};

fn controls() -> LagrangianSgsConfig {
    LagrangianSgsConfig {
        memory_time: 1000.0,
        maximum_departure_cells: 1.0,
        maximum_coefficient_squared: 100.0,
        maximum_eddy_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn heat_controls() -> LagrangianScalarConfig {
    LagrangianScalarConfig {
        memory_time: 0.7,
        maximum_departure_cells: 1.0,
        maximum_coefficient: 100.0,
        maximum_diffusivity_to_molecular: 1000.0,
        usage: SgsUse::Exploration,
    }
}

fn update(
    shape: [usize; 3],
    spacing: [f64; 3],
    normal: Option<usize>,
    elapsed: f64,
    substeps: u8,
) -> HistoryUpdate {
    let geometry = normal.map_or(FilterGeometry::Periodic, FilterGeometry::Channel);
    let flow =
        PeriodicDynamicSgs::new_cartesian_geometry(shape, spacing, geometry, usize::MAX).unwrap();
    HistoryUpdate::new(&flow, true, Some(0.0), Some(0.0), elapsed, 1000.0, 1.0)
        .unwrap()
        .with_advection(LagrangianAdvection::Midpoint { substeps }, [[0.0; 3]; 2])
}

fn spiral(shape: [usize; 3], spacing: [f64; 3], stretch: f64, omega: f64) -> Vec<[f64; 3]> {
    (0..shape.iter().product())
        .map(|i| {
            let p = coordinates(shape, i);
            let x = (p[0] as f64 - (shape[0] / 2) as f64) * spacing[0];
            let y = (p[1] as f64 - (shape[1] / 2) as f64) * spacing[1];
            [stretch * x - omega * y, omega * x - stretch * y, 0.0]
        })
        .collect()
}

fn spiral_departure(x: f64, y: f64, stretch: f64, omega: f64, elapsed: f64) -> [f64; 2] {
    // A = [[a,-w],[w,-a]], A^2 = -(w^2-a^2) I. Therefore
    // exp(-t A) = cos(beta t) I - sin(beta t) A/beta.
    let beta = (omega * omega - stretch * stretch).sqrt();
    let c = (beta * elapsed).cos();
    let s = (beta * elapsed).sin() / beta;
    [
        x - (c * x - s * (stretch * x - omega * y)),
        y - (c * y - s * (omega * x - stretch * y)),
    ]
}

fn error(a: [f64; 3], b: [f64; 3]) -> f64 {
    a.into_iter()
        .zip(b)
        .map(|(a, b)| (a - b).powi(2))
        .sum::<f64>()
        .sqrt()
}

#[test]
fn g1_midpoint_substeps_resolve_an_independent_curved_affine_characteristic() {
    let shape = [9, 9, 5];
    let spacing = [0.5, 0.75, 1.25];
    let cell = index(shape, [5, 4, 2]);
    let velocity = spiral(shape, spacing, 0.05, 0.3);
    let elapsed = 1.0;
    let physical = spiral_departure(spacing[0], 0.0, 0.05, 0.3, elapsed);
    let expected = [physical[0] / spacing[0], physical[1] / spacing[1], 0.0];
    let euler = core::array::from_fn(|a| velocity[cell][a] * elapsed / spacing[a]);
    let mut previous = error(euler, expected);
    for n in [1, 2, 4, 8] {
        let step = update(shape, spacing, None, elapsed, n);
        let actual = step
            .midpoint_departure(cell, &velocity, &mut || false)
            .unwrap();
        let observed = error(actual, expected);
        assert!(
            previous / observed > if n == 1 { 6.0 } else { 3.5 },
            "n={n}, prior={previous}, error={observed}"
        );
        previous = observed;
    }
}

#[test]
fn g1_curved_history_refinement_changes_actual_emitted_sgs_viscosity() {
    let shape = [9, 9, 5];
    let cell = index(shape, [5, 4, 2]);
    let velocity = spiral(shape, [1.0; 3], 0.05, 0.15);
    let rho = vec![1.0; velocity.len()];
    let exact = spiral_departure(1.0, 0.0, 0.05, 0.15, 1.0);
    let advected = 0.25 + 0.02 * (5.0 - exact[0]) + 0.01 * (4.0 - exact[1]);
    let old = (-1.0_f64 / controls().memory_time).exp();
    // Independent current Germano source for this incompressible, linear,
    // isotropic-grid spiral: L:M=0, M:M=2048 Delta^4 a^4, |S|=2a.
    let denominator = old + (1.0 - old) * 2048.0 * 0.05_f64.powi(4);
    let expected = (old * advected / denominator) * 0.1;
    let mut previous = f64::INFINITY;
    for policy in [
        LagrangianAdvection::Euler,
        LagrangianAdvection::Midpoint { substeps: 1 },
        LagrangianAdvection::Midpoint { substeps: 2 },
        LagrangianAdvection::Midpoint { substeps: 4 },
    ] {
        let mut owner = LagrangianDynamicSgs::new(shape, [1.0; 3], usize::MAX)
            .unwrap()
            .with_advection(policy)
            .unwrap();
        // Manufactured accepted physical histories; every subsequent filter,
        // relaxation, coefficient and viscosity comes from the real owner.
        owner.accepted_time = Some(0.0);
        owner.history_log_scale = Some(0.0);
        for (i, pair) in owner.history.iter_mut().enumerate() {
            let p = coordinates(shape, i);
            *pair = [0.25 + 0.02 * p[0] as f64 + 0.01 * p[1] as f64, 1.0];
        }
        let trial = owner
            .prepare(controls(), 1.0, &rho, &velocity, || false)
            .unwrap();
        let observed = (trial.viscosities()[cell] - expected).abs();
        assert!(observed > 0.0);
        assert!(
            previous / observed > 3.0,
            "{policy:?}: {observed} versus {previous}"
        );
        assert!(trial.coefficient_squared(cell).unwrap() > 0.0);
        previous = observed;
    }
}

fn field(shape: [usize; 3], amplitude: f64) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let n = shape.iter().product();
    let rho = (0..n).map(|i| 1.0 + (i % 7) as f64 / 64.0).collect();
    let u = (0..n)
        .map(|i| {
            let p = coordinates(shape, i);
            let [x, y, z] = core::array::from_fn(|a| {
                core::f64::consts::TAU * (p[a] as f64 + 0.5) / shape[a] as f64
            });
            [
                amplitude * (x.sin() + (y + z).sin()),
                amplitude * ((x - z).cos() + y.sin()),
                amplitude * ((x + y).sin() + z.cos()),
            ]
        })
        .collect();
    let phi = (0..n)
        .map(|i| 0.2 + ((i * i * 3 + 7 * i) % 41) as f64 / 128.0)
        .collect();
    (rho, u, phi)
}

#[test]
fn g3_joint_stress_heat_and_species_share_one_policy_without_new_payload() {
    let shape = [5; 3];
    let mc = controls();
    let hc = heat_controls();
    let cc = LagrangianComponentsConfig {
        memory_time: hc.memory_time,
        maximum_departure_cells: 1.0,
        maximum_coefficient: hc.maximum_coefficient,
        usage: SgsUse::Exploration,
    };
    let policy = LagrangianAdvection::Midpoint { substeps: 3 };
    let mut all = LagrangianDynamicComponents::new(shape, [1.0; 3], 1, usize::MAX)
        .unwrap()
        .with_advection(policy)
        .unwrap();
    let mut joint = LagrangianDynamicScalar::new(shape, [1.0; 3], usize::MAX)
        .unwrap()
        .with_advection(policy)
        .unwrap();
    let mut euler = LagrangianDynamicScalar::new(shape, [1.0; 3], usize::MAX).unwrap();
    let mut explicit = LagrangianDynamicScalar::new(shape, [1.0; 3], usize::MAX)
        .unwrap()
        .with_advection(LagrangianAdvection::Euler)
        .unwrap();
    assert_eq!(joint.workspace_bytes(), euler.workspace_bytes());
    assert_eq!(
        all.workspace_bytes(),
        LagrangianDynamicComponents::required_bytes(125).unwrap()
    );
    let mut different = false;
    for (time, amplitude) in [(0.0, 0.2), (0.5, 0.12), (1.0, 0.17)] {
        let (rho, u, phi) = field(shape, amplitude);
        let a = all
            .prepare(mc, hc, cc, time, &rho, &u, &phi, |_, i| phi[i], || false)
            .unwrap();
        let b = joint
            .prepare(mc, hc, time, &rho, &u, &phi, || false)
            .unwrap();
        let c = euler
            .prepare(mc, hc, time, &rho, &u, &phi, || false)
            .unwrap();
        let d = explicit
            .prepare(mc, hc, time, &rho, &u, &phi, || false)
            .unwrap();
        assert_eq!(a.momentum_report(), b.momentum_report());
        assert_eq!(a.scalar_report(), b.scalar_report());
        assert_eq!(a.components_report(), b.scalar_report());
        assert_eq!(a.viscosities(), b.viscosities());
        assert_eq!(a.scalar_diffusivities(), b.diffusivities());
        assert_eq!(a.component_diffusivities(), b.diffusivities());
        assert_eq!(
            c.viscosities()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            d.viscosities()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
        assert_eq!(
            c.diffusivities()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>(),
            d.diffusivities()
                .iter()
                .map(|v| v.to_bits())
                .collect::<Vec<_>>()
        );
        different |= b
            .viscosities()
            .iter()
            .zip(c.viscosities())
            .any(|(a, b)| (a - b).abs() > 1e-10);
        a.commit();
        b.commit();
        c.commit();
        d.commit();
    }
    assert!(different);
    assert_eq!(all.advection(), policy);
    assert_eq!(joint.advection(), policy);
    let mut fork = joint.fork_for_trial(joint.workspace_bytes()).unwrap();
    assert_eq!(fork.advection(), policy);
    let (rho, u, phi) = field(shape, 0.1);
    let a = joint
        .prepare(mc, hc, 1.25, &rho, &u, &phi, || false)
        .unwrap();
    let b = fork
        .prepare(mc, hc, 1.25, &rho, &u, &phi, || false)
        .unwrap();
    assert_eq!(a.viscosities(), b.viscosities());
    assert_eq!(a.diffusivities(), b.diffusivities());
    drop(a);
    drop(b);
    assert_eq!(joint.accepted_time(), Some(1.0));
    assert_eq!(fork.accepted_time(), Some(1.0));
}

#[test]
fn g1_channel_characteristics_sample_actual_half_cell_moving_walls() {
    let shape = [5; 3];
    for normal in 0..3 {
        let tangential = (normal + 1) % 3;
        let mut p = [2; 3];
        p[normal] = 0;
        let cell = index(shape, p);
        let mut step = update(shape, [1.0; 3], Some(normal), 0.5, 1);
        step.wall_velocity[0][tangential] = 0.4;
        let mut u = [0.0; 3];
        u[normal] = 0.4;
        u[tangential] = 0.1;
        let velocity = vec![u; 125];
        // Euler midpoint is one tenth of a cell into the lower half cell:
        // wall weight .2, u_normal=.32, u_tangent=.16 at that midpoint.
        let departure = step
            .midpoint_departure(cell, &velocity, &mut || false)
            .unwrap();
        assert!((departure[normal] - 0.16).abs() < 1e-15);
        assert!((departure[tangential] - 0.08).abs() < 1e-15);
        let history = vec![[0.25, 1.0]; 125];
        let sampled = step
            .advected_with_field(cell, &velocity, &history, &mut || false)
            .unwrap();
        assert!((sampled[0] - 0.25).abs() < 1e-15);
        assert!((sampled[1] - 1.0).abs() < 1e-15);
        let mut escaped = velocity.clone();
        escaped[cell][normal] = 4.0;
        assert!(
            step.midpoint_departure(cell, &escaped, &mut || false)
                .is_err()
        );
    }
}

#[test]
fn g3_tracing_preserves_uniform_advection_and_physical_unit_changes() {
    let shape = [5; 3];
    let cell = index(shape, [3, 2, 2]);
    let h = [0.5, 0.75, 1.25];
    let v = [0.125, -0.25, 0.0625];
    let elapsed = 0.5;
    let expected = core::array::from_fn(|a| v[a] * elapsed / h[a]);
    for count in [1, 2, 4, 16] {
        let step = update(shape, h, None, elapsed, count);
        assert!(
            error(
                step.midpoint_departure(cell, &vec![v; 125], &mut || false)
                    .unwrap(),
                expected
            ) < 2e-16
        );
    }
    let base = spiral(shape, h, 0.05, 0.3);
    let expected = update(shape, h, None, elapsed, 4)
        .midpoint_departure(cell, &base, &mut || false)
        .unwrap();
    for (length, time) in [(1e-60, 1e60), (1e60, 1e-60)] {
        let u: Vec<_> = base
            .iter()
            .map(|row| row.map(|v| v * (length / time)))
            .collect();
        let step = update(shape, h.map(|v| v * length), None, elapsed * time, 4);
        let actual = step.midpoint_departure(cell, &u, &mut || false).unwrap();
        assert!(error(actual, expected) < 1e-14);
    }
}

#[test]
fn g4_midpoint_limits_adaptive_memory_and_every_poll_refuse_without_history_publication() {
    let shape = [3; 3];
    let policy = LagrangianAdvection::Midpoint { substeps: 9 };
    assert!(
        LagrangianAdvection::Midpoint { substeps: 0 }
            .validate()
            .is_err()
    );
    let mut owner = LagrangianDynamicSgs::new(shape, [1.0; 3], usize::MAX)
        .unwrap()
        .with_advection(policy)
        .unwrap()
        .with_adaptive_memory(AdaptiveLagrangianMemory {
            coefficient: 0.15,
            minimum_time: 0.01,
        })
        .unwrap();
    let (rho, u, _) = field(shape, 0.2);
    owner
        .prepare(controls(), 0.0, &rho, &u, || false)
        .unwrap()
        .commit();
    let history = owner.history.clone();
    let scale = owner.history_log_scale;
    let pointers = (
        owner.history.as_ptr(),
        owner.candidate.as_ptr(),
        owner.flow.velocity.as_ptr(),
    );
    let bytes = owner.workspace_bytes();
    let mut polls = 0;
    let trial = owner
        .prepare(controls(), 0.25, &rho, &u, || {
            polls += 1;
            false
        })
        .unwrap();
    let expected = trial.viscosities().to_vec();
    let report = trial.report();
    let range = trial.memory_time_range();
    drop(trial);
    assert!(polls > 30, "substep cancellation polls must be active");
    for stop in 1..=polls {
        let mut seen = 0;
        assert!(matches!(
            owner.prepare(controls(), 0.25, &rho, &u, || {
                seen += 1;
                seen == stop
            }),
            Err(DynamicSgsError::Cancelled)
        ));
        assert_eq!(owner.history, history);
        assert_eq!(owner.history_log_scale, scale);
        assert_eq!(owner.accepted_time(), Some(0.0));
        let retry = owner.prepare(controls(), 0.25, &rho, &u, || false).unwrap();
        assert_eq!(retry.viscosities(), expected);
        assert_eq!(retry.report(), report);
        assert_eq!(retry.memory_time_range(), range);
        drop(retry);
    }
    assert_eq!(bytes, owner.workspace_bytes());
    assert_eq!(
        pointers,
        (
            owner.history.as_ptr(),
            owner.candidate.as_ptr(),
            owner.flow.velocity.as_ptr()
        )
    );
    let mut impossible = controls();
    impossible.maximum_departure_cells = 1e-8;
    assert!(owner.prepare(impossible, 0.25, &rho, &u, || false).is_err());
    assert_eq!(owner.history, history);
    assert_eq!(owner.accepted_time(), Some(0.0));
    assert!(owner.with_advection(LagrangianAdvection::Euler).is_err());
}
