use super::*;
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::WallPin;

fn field(temperature: [f64; 2], fraction: [f64; 2], pressure: f64) -> AxialGasField {
    crate::profile::gas_field::tests::field(temperature, fraction, pressure)
}
fn resolution() -> DgWallResolution {
    DgWallResolution {
        poles_per_field: 16,
        quadrature_points: 2,
        band_rad_s: [20.0, 4000.0],
        maximum_relative_error: 0.001,
    }
}
fn build(field: &AxialGasField) -> DgDuct {
    let pins = [
        WallPin { surface_density: 0.2, stiffness_per_area: 1e7, resistance: 30.0 },
        WallPin { surface_density: 0.7, stiffness_per_area: 4e7, resistance: 80.0 },
    ];
    let laws = [
        GeneralizedMaxwell::new(1e6, vec![(2e5, 0.001)]).unwrap(),
        GeneralizedMaxwell::new(4e6, vec![(5e5, 0.002), (9e5, 0.0004)]).unwrap(),
    ];
    DgDuct::new_with_gas_field(
        &Duct { segments: vec![Segment::Cylinder { radius: 0.001, length: 0.57 }] },
        field, DgMesh { elements: 2, order: 1 },
        DgBoundary::Closed, DgBoundary::Port(PortId(1)),
    ).unwrap().with_side_ports(&[DgSidePort { id: PortId(7), downstream_element: 1 }]).unwrap()
        .with_spatial_viscothermal(resolution()).unwrap()
        .with_compliant_wall_relaxation_regions(&pins, &[Some(&laws[0]), Some(&laws[1])]).unwrap()
}
fn memory(duct: &DgDuct) -> &SpatialViscothermal {
    match duct.viscothermal.as_ref().unwrap() {
        ViscothermalMemory::Spatial(memory) => memory,
        _ => panic!("expected native spatial diffusion"),
    }
}
fn near(a: f64, b: f64, scale: f64) {
    assert!((a - b).abs() < 5e-9 * scale.max(1e-30), "{a:e} != {b:e}");
}

// Independent physical quadrature, not the native mass roots or transfer map.
// Two linear shape functions give each weak conserved variable and wave energy.
fn weak_and_energy(
    field: &AxialGasField, element: &DgElement, pressure: &[f64], flow: &[f64],
) -> ([f64; 2], [f64; 2], f64) {
    let (points, weights) = gauss_legendre_rule(32);
    let area = core::f64::consts::PI * 0.001_f64.powi(2);
    let mut compression = [0.0; 2];
    let mut momentum = [0.0; 2];
    let mut energy = 0.0;
    for (xi, weight) in points.into_iter().zip(weights) {
        let phi = [0.5 * (1.0 - xi), 0.5 * (1.0 + xi)];
        let x = element.x_left + 0.5 * element.length * (xi + 1.0);
        let gas = *field.state_at(x).unwrap().state();
        let compliance = area / (gas.density * gas.sound_speed.powi(2));
        let inertia = gas.density / area;
        let p = phi[0] * pressure[0] + phi[1] * pressure[1];
        let q = phi[0] * flow[0] + phi[1] * flow[1];
        let measure = 0.5 * element.length * weight;
        for i in 0..2 {
            compression[i] += measure * compliance * p * phi[i];
            momentum[i] += measure * inertia * q * phi[i];
        }
        energy += 0.5 * measure * (compliance * p * p + inertia * q * q);
    }
    (compression, momentum, energy)
}

#[test]
fn regional_gas_preserves_independent_weak_and_arm_storage_and_reloads_actual_wall_motion() {
    let old_field = field([290.0, 410.0], [0.9, 0.3], 100_000.0);
    let new_field = field([480.0, 330.0], [0.2, 0.75], 150_000.0);
    let mut duct = build(&old_field);
    let old_step = duct.propagator(1e-6, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|x| 10.0 + 4.0 * x, |_| 1e-6);
    for _ in 0..16 { old_step.step(&mut state, &[]).unwrap(); }
    old_step.sync_nodal(&mut state).unwrap();
    let original = duct.clone();
    let before = state.clone();
    let old_motion = duct.compliant_wall_state(&state).unwrap();
    assert!(old_motion.relaxation_stored_energy_j > 0.0);
    let prefix = 2 * duct.field_dofs + memory(&duct).admission.memory_dofs;
    assert!(energy(&before.tilde[2 * duct.field_dofs..prefix]) > 0.0);
    let report = duct.rebind_spatial_viscothermal_gas_field(
        &mut state, &new_field, 0.001, 1.0, &mut || false,
    ).unwrap();
    assert_eq!(state.step_count(), before.step_count());
    assert_eq!(state.initial_energy_j(), before.initial_energy_j());
    assert_eq!(state.cumulative_dissipated_j(), before.cumulative_dissipated_j());
    assert_eq!(&state.tilde[prefix..], &before.tilde[prefix..]);
    assert_eq!(duct.gas_field(), Some(&new_field));
    let mut expected_work = 0.0;
    let mut commutator = 0.0_f64;
    let mut matrix_scale = 0.0_f64;
    for (e, (old, new)) in original.elements.iter().zip(&duct.elements).enumerate() {
        let range = 2 * e..2 * e + 2;
        let (ca, pa, ea) = weak_and_energy(&old_field, old, &before.p[range.clone()], &before.q[range.clone()]);
        let (cb, pb, eb) = weak_and_energy(&new_field, new, &state.p[range.clone()], &state.q[range]);
        for i in 0..2 {
            near(ca[i], cb[i], ca[i].abs().max(cb[i].abs()));
            near(pa[i], pb[i], pa[i].abs().max(pb[i].abs()));
            for j in 0..2 {
                let ab: f64 = (0..2).map(|k| old.mp_sqrt[2*i+k] * new.mp_sqrt[2*k+j]).sum();
                let ba: f64 = (0..2).map(|k| new.mp_sqrt[2*i+k] * old.mp_sqrt[2*k+j]).sum();
                commutator = commutator.max((ab - ba).abs());
                matrix_scale = matrix_scale.max(ab.abs());
            }
        }
        expected_work += eb - ea;
    }
    assert!(commutator > 1e-7 * matrix_scale, "fixture must exercise noncommuting mass roots");
    let (_, quadrature_weights) = gauss_legendre_rule(resolution().quadrature_points);
    let old_memory = memory(&original);
    let new_memory = memory(&duct);
    for (index, (a, b)) in old_memory.channels.iter().zip(&new_memory.channels).enumerate() {
        let element = index / (2 * resolution().quadrature_points);
        let point = (index / 2) % resolution().quadrature_points;
        let measure = 0.5 * original.elements[element].length * quadrature_weights[point];
        let area = core::f64::consts::PI * a.radius_m.powi(2);
        let old_gas = *old_field.state_at(a.point_m).unwrap().state();
        let new_gas = *new_field.state_at(b.point_m).unwrap().state();
        let pressure = a.field_start < original.field_dofs;
        let coefficient = |gas: GasState| if pressure {
            area / (gas.density * gas.sound_speed.powi(2))
        } else { gas.density / area };
        for (arm, (ta, tb)) in a.terms.iter().zip(&b.terms).enumerate() {
            let ma = measure * coefficient(old_gas) * ta.gain_per_s / ta.rate_per_s;
            let mb = measure * coefficient(new_gas) * tb.gain_per_s / tb.rate_per_s;
            let ya = before.tilde[a.memory_start + arm];
            let yb = state.tilde[b.memory_start + arm];
            let conserved = ma.sqrt() * ya;
            near(conserved, mb.sqrt() * yb, conserved.abs());
            expected_work += 0.5 * conserved * conserved * (1.0 / mb - 1.0 / ma);
        }
    }
    near(report.parameter_work_j, expected_work, before.stored_energy_j());
    near(report.energy_after_j - report.energy_before_j, expected_work, before.stored_energy_j());
    let motion = duct.compliant_wall_state(&state).unwrap();
    assert_eq!(motion, old_motion);
    let cold = build(&new_field);
    for (&a, &b) in duct.l_tilde.iter().zip(&cold.l_tilde) { near(a, b, b.abs().max(1.0)); }
    for ((id, a), (other_id, b)) in duct.ports.iter().zip(&cold.ports) {
        assert_eq!(id, other_id);
        for (&a, &b) in a.iter().zip(b) { near(a, b, b.abs().max(1.0)); }
    }
    assert_eq!(duct.side_port_source(PortId(7)).unwrap().gas, *new_field.state_at(0.285).unwrap().state());
    let n = duct.total_dofs();
    let mut derivative = state.clone();
    for (i, row) in duct.l_tilde.chunks_exact(n).enumerate() {
        derivative.tilde[i] = row.iter().zip(&state.tilde).map(|(a, x)| a * x).sum();
    }
    let acceleration = duct.compliant_wall_state(&derivative).unwrap();
    for i in 0..duct.field_dofs {
        let pin = duct.compliant_wall_at_element(i / 2).unwrap();
        let reaction: f64 = motion.relaxation_pressure_pa.iter().map(|arm| arm[i]).sum();
        near(pin.surface_density * acceleration.velocity_m_s[i], state.p[i]
            - pin.stiffness_per_area * motion.displacement_m[i]
            - pin.resistance * motion.velocity_m_s[i] - reaction, 20.0);
    }
    let power = -state.tilde.iter().zip(&derivative.tilde).map(|(x, dx)| x * dx).sum::<f64>();
    near(power, duct.wall_dissipated_power_w(&state).unwrap(), power);
    let step = duct.propagator(1e-6, 1e-13).unwrap();
    let mut control = before.clone();
    let mut difference = 0.0;
    for _ in 0..12 {
        let frame = step.step(&mut state, &[]).unwrap();
        old_step.step(&mut control, &[]).unwrap();
        step.sync_nodal(&mut state).unwrap();
        old_step.sync_nodal(&mut control).unwrap();
        difference += state.p.iter().zip(&control.p).map(|(a, b)| (a-b).powi(2)).sum::<f64>();
        assert_eq!(frame.ledger.cumulative_parameter_work_j, report.parameter_work_j);
        near(frame.ledger.energy_balance_residual_j, 0.0, before.initial_energy_j());
    }
    assert!(difference > 1e-4);
}

#[test]
fn uniform_observations_survive_spatial_round_trip_and_every_late_refusal_is_atomic() {
    let uniform = field([300.0; 2], [0.8; 2], 100_000.0);
    let regional = field([310.0, 430.0], [0.9, 0.2], 140_000.0);
    let mut duct = build(&uniform);
    let stepping = duct.propagator(1e-6, 1e-13).unwrap();
    let mut state = duct.project_initial_state(|x| 10.0 + x, |_| 1e-6);
    for _ in 0..10 { stepping.step(&mut state, &[]).unwrap(); }
    let motion = duct.compliant_wall_state(&state).unwrap();
    for target in [&regional, &uniform] {
        duct.rebind_spatial_viscothermal_gas_field(&mut state, target, 0.001, 1.0, &mut || false).unwrap();
        assert_eq!(duct.compliant_wall_state(&state).unwrap(), motion);
    }
    let accepted = duct.clone();
    let accepted_state = state.clone();
    let accepted_step = duct.propagator(1e-6, 1e-13).unwrap();
    for (error, budget) in [(1e-16, 1.0), (0.001, 0.0)] {
        assert!(duct.rebind_spatial_viscothermal_gas_field(&mut state, &regional, error, budget, &mut || false).is_err());
        assert_eq!(state, accepted_state);
        assert_eq!(duct.operator_identity(), accepted.operator_identity());
    }
    use fs_material::{gas_field::PrimitiveGasState, gas_species::input::ThermochemicalPhase};
    use fs_thermochem::{Composition, CompositionBasis};
    let moved = AxialGasField::new(
        vec![0.0, 0.17, 0.57], 100_000.0,
        [300.0, 330.0, 360.0].into_iter().map(|temperature_k| PrimitiveGasState {
            temperature_k,
            composition: Composition::new(CompositionBasis::MoleFraction, vec![0.8, 0.2]).unwrap(),
        }).collect(),
        uniform.model().clone(), fs_matdb::QueryPoint::new(), ThermochemicalPhase::FrozenGas,
    ).unwrap();
    assert!(duct.rebind_spatial_viscothermal_gas_field(&mut state, &moved, 0.001, 1.0, &mut || false).is_err());
    assert_eq!(state, accepted_state);
    let mut polls = 0;
    accepted.clone().rebind_spatial_viscothermal_gas_field(
        &mut accepted_state.clone(), &regional, 0.001, 1.0, &mut || { polls += 1; false },
    ).unwrap();
    let mut count = 0;
    assert!(matches!(duct.rebind_spatial_viscothermal_gas_field(
        &mut state, &regional, 0.001, 1.0, &mut || { count += 1; count == polls },
    ), Err(GasFieldDuctError::Duct(DuctError::Cancelled { .. }))));
    assert_eq!(state, accepted_state);
    assert_eq!(duct.operator_identity(), accepted.operator_identity());
    assert_eq!(duct.gas_field(), accepted.gas_field());
    let mut control = accepted_state.clone();
    assert_eq!(accepted_step.step(&mut state, &[]).unwrap(), accepted_step.step(&mut control, &[]).unwrap());
}
