//! Independent physical checks of original, element-supported Maxwell histories.
use super::*;
use super::regional::base;
use fs_material::visco::GeneralizedMaxwell;

fn pins() -> [WallPin; 2] {
    [wall(300.0), WallPin {
        surface_density: 0.7,
        stiffness_per_area: 4e7,
        resistance: 80.0,
    }]
}
fn spectra() -> [GeneralizedMaxwell; 2] {
    [
        GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001), (0.0, 0.3)]).unwrap(),
        GeneralizedMaxwell::new(4e6, vec![(5e5, 0.002), (9e5, 0.0004)]).unwrap(),
    ]
}

/// G1: recover the original dimensional wall and constitutive equations,
/// then independently integrate actual geometric spring energy and heat.
#[test]
fn each_original_arm_obeys_its_local_equation_and_physical_storage() {
    let pins = pins();
    let spectra = spectra();
    let arms = [(0, 2e6, 0.001), (1, 5e6, 0.002), (1, 9e6, 0.0004)];
    for spatial in [false, true] {
        let dg = base(spatial).with_compliant_wall_relaxation_regions(
            &pins, &[Some(&spectra[0]), Some(&spectra[1])],
        ).unwrap();
        let local = dg.mesh.order + 1;
        let m = dg.field_dofs();
        assert_eq!(dg.total_dofs(), 4 * m + 3 * local);
        let n = dg.total_dofs();
        let mut current = dg.project_initial_state(|x| 10.0 + x, |_| 0.0);
        for (i, value) in current.tilde[2 * m..].iter_mut().enumerate() {
            *value = 1e-6 * ((i * 17 % 23) as f64 - 11.0);
        }
        let mut derivative = current.clone();
        for (i, row) in dg.l_tilde.chunks_exact(n).enumerate() {
            derivative.tilde[i] = row.iter().zip(&current.tilde).map(|(a, x)| a * x).sum();
        }
        let state = dg.compliant_wall_state(&current).unwrap();
        let rate = dg.compliant_wall_state(&derivative).unwrap();
        let (pressure, _) = dg.from_energy_coords(&current.tilde);
        assert_eq!(state.relaxation_pressure_pa.len(), arms.len());
        for node in 0..m {
            let pin = pins[node / local];
            close(rate.displacement_m[node], state.velocity_m_s[node], 1.0);
            let reaction: f64 = state.relaxation_pressure_pa.iter().map(|arm| arm[node]).sum();
            let force = pressure[node] - pin.stiffness_per_area * state.displacement_m[node]
                - pin.resistance * state.velocity_m_s[node] - reaction;
            close(pin.surface_density * rate.velocity_m_s[node], force, force.abs().max(1.0));
            for (arm, &(element, stiffness, tau)) in arms.iter().enumerate() {
                if node / local == element {
                    let expected = stiffness * state.velocity_m_s[node]
                        - state.relaxation_pressure_pa[arm][node] / tau;
                    close(rate.relaxation_pressure_pa[arm][node], expected, expected.abs().max(1.0));
                } else {
                    assert_eq!(state.relaxation_pressure_pa[arm][node], 0.0);
                    assert_eq!(rate.relaxation_pressure_pa[arm][node], 0.0);
                }
            }
        }
        let (points, weights) = gauss_legendre_rule(32);
        let perimeter = core::f64::consts::TAU * RADIUS;
        let (mut mechanical, mut memory, mut loss) = (0.0, 0.0, 0.0);
        for (e, element) in dg.elements.iter().enumerate() {
            for (&point, &weight) in points.iter().zip(&weights) {
                let shape = dg.lgl.eval_basis(point);
                let interp = |values: &[f64]| -> f64 {
                    (0..local).map(|i| shape[i] * values[e * local + i]).sum()
                };
                let area = perimeter * 0.5 * element.length * weight;
                let w = interp(&state.displacement_m);
                let v = interp(&state.velocity_m_s);
                mechanical += 0.5 * area * (
                    pins[e].surface_density * v * v + pins[e].stiffness_per_area * w * w);
                loss += area * pins[e].resistance * v * v;
                for (arm, &(owner, stiffness, tau)) in arms.iter().enumerate() {
                    if owner == e {
                        let z = interp(&state.relaxation_pressure_pa[arm]);
                        memory += 0.5 * area * z * z / stiffness;
                        loss += area * z * z / (stiffness * tau);
                    }
                }
            }
        }
        close(state.relaxation_stored_energy_j, memory, memory);
        close(state.stored_energy_j, mechanical + memory, mechanical + memory);
        close(state.dissipated_power_w, loss, loss);
        let energy_rate: f64 = current.tilde.iter().zip(&derivative.tilde).map(|(x, d)| x * d).sum();
        close(energy_rate, -loss, loss);
    }
}

#[test]
fn regional_histories_retain_gas_memory_ports_and_change_actual_pressure() {
    let gas = duct(DgBoundary::Port(PortId(0)))
        .with_viscothermal(8, [50.0, 2000.0], 0.01).unwrap();
    let pins = pins();
    let base = gas.clone().with_compliant_wall_regions(&pins).unwrap();
    let spectra = spectra();
    let mut changed = spectra[1].clone();
    changed.terms[1].1 *= 5.0;
    let mut histories = Vec::new();
    let mut operators = Vec::new();
    for right in [&spectra[1], &changed] {
        let dg = gas.clone().with_compliant_wall_relaxation_regions(
            &pins, &[Some(&spectra[0]), Some(right)],
        ).unwrap();
        let original = base.total_dofs();
        let n = dg.total_dofs();
        for i in 0..original {
            assert_eq!(&dg.l_tilde[i * n..i * n + original],
                &base.l_tilde[i * original..(i + 1) * original]);
        }
        for ((id, column), (old_id, old_column)) in dg.ports.iter().zip(&base.ports) {
            assert_eq!(id, old_id);
            assert_eq!(&column[..original], old_column.as_slice());
            assert!(column[original..].iter().all(|value| *value == 0.0));
        }
        let propagator = dg.propagator(2e-6, 1e-13).unwrap();
        let mut current = dg.zero_state();
        let mut work = 0.0;
        let mut history = Vec::new();
        for _ in 0..600 {
            let frame = propagator.step(&mut current,
                &[PortInput { id: PortId(0), effort: 1e-7 }]).unwrap();
            work += frame.ledger.port_work_step_j;
            propagator.sync_nodal(&mut current).unwrap();
            history.push(*current.pressure().last().unwrap());
        }
        assert!(current.stored_energy_j() > 0.0 && current.stored_energy_j() < work);
        assert!(dg.compliant_wall_state(&current).unwrap().relaxation_stored_energy_j > 0.0);
        histories.push(history);
        operators.push(dg);
    }
    assert_ne!(operators[0].operator_identity(), operators[1].operator_identity());
    assert!(operators[1].compliant_wall_state(&operators[0].zero_state()).is_err());
    let signal: f64 = histories[0].iter().map(|p| p * p).sum();
    let difference: f64 = histories[0].iter().zip(&histories[1]).map(|(a, b)| (a-b).powi(2)).sum();
    assert!(difference > 1e-12 * signal);
}

#[test]
fn uniform_zero_and_partial_spectra_keep_their_exact_support() {
    let law = spectra()[0].clone();
    let empty = GeneralizedMaxwell::new(4e6, vec![(0.0, 0.002)]).unwrap();
    for spatial in [false, true] {
        let gas = base(spatial);
        let direct = gas.clone().with_compliant_wall_relaxation(wall(300.0), &law).unwrap();
        let regional = gas.clone().with_compliant_wall_relaxation_regions(
            &[wall(300.0); 2], &[Some(&law); 2],
        ).unwrap();
        assert_eq!(direct.operator_identity(), regional.operator_identity());
        assert_eq!(direct.l_tilde, regional.l_tilde);
        let pins = pins();
        let plain = gas.clone().with_compliant_wall_regions(&pins).unwrap();
        let zero = gas.clone().with_compliant_wall_relaxation_regions(
            &pins, &[None, Some(&empty)],
        ).unwrap();
        assert_eq!(plain.operator_identity(), zero.operator_identity());
        let local = gas.mesh.order + 1;
        let partial = gas.with_compliant_wall_relaxation_regions(
            &pins, &[None, Some(&law)],
        ).unwrap();
        assert_eq!(partial.total_dofs(), plain.total_dofs() + local);
        let mut current = partial.zero_state();
        current.tilde[plain.total_dofs()..].fill(1e-6);
        let observation = partial.compliant_wall_state(&current).unwrap();
        assert_eq!(observation.relaxation_pressure_pa.len(), 1);
        assert!(observation.relaxation_pressure_pa[0][..local].iter().all(|p| *p == 0.0));
        assert!(observation.relaxation_pressure_pa[0][local..].iter().any(|p| *p != 0.0));
    }
}

#[test]
fn invalid_regional_spectra_and_budget_or_cancellation_never_publish_history() {
    let original = base(true);
    let pins = pins();
    let spectra = spectra();
    let identity = original.operator_identity();
    assert!(original.clone().with_compliant_wall_relaxation_regions(&pins, &[None]).is_err());
    assert!(original.clone().with_compliant_wall_relaxation_regions(&pins[..1], &[None; 2]).is_err());
    for bad in [
        GeneralizedMaxwell { e_inf: 0.0, terms: vec![] },
        GeneralizedMaxwell { e_inf: 1e6, terms: vec![(1e6, 0.0)] },
        GeneralizedMaxwell { e_inf: 1e6, terms: vec![(1e6, 1e-320)] },
        GeneralizedMaxwell::new(1e6, vec![(1e6, 0.001); 400]).unwrap(),
    ] {
        assert!(original.clone().with_compliant_wall_relaxation_regions(
            &pins, &[None, Some(&bad)],
        ).is_err());
    }
    for stop in [1, 4, 20, 80] {
        let mut polls = 0;
        assert!(matches!(original.clone().with_compliant_wall_relaxation_regions_controlled(
            &pins, &[Some(&spectra[0]), Some(&spectra[1])],
            &mut || { polls += 1; polls == stop },
        ), Err(DuctError::Cancelled { .. })));
    }
    assert_eq!(original.operator_identity(), identity);
    assert!(original.compliant_wall_at_element(0).is_none());
}
