//! Original differential/algebraic attachments on the actual shell material owners.
use super::*;
use fs_couple::shear_shell::coating::poles::{PoleOptions, ShellResonanceDescriptorLoad};
use fs_phs::DescriptorPortHamiltonian;

// Original momentum p and zero-storage reaction lambda:
// p' = -(drag/mass)*p + lambda, 0 = -p/mass + map*v, y = map*lambda.
fn inertance(mass: f64, drag: f64, map: f64) -> DescriptorPortHamiltonian {
    DescriptorPortHamiltonian::new(
        2, 1, 1,
        vec![0.0, 1.0, -1.0, 0.0],
        vec![drag, 0.0, 0.0, 0.0],
        vec![0.0, map],
        Box::new(QuadraticStorage::new(vec![1.0 / mass], 1).unwrap()),
    ).unwrap()
}

fn inertance_pole(image: &MaterialShellImage, w: f64, added_mass: f64, drag: f64) -> C64 {
    let mut s = C64::new(-0.01 * w, -w / (1.0 + added_mass).sqrt());
    for _ in 0..64 {
        let value = body_operator(image, s, w) + (s * s).scale(added_mass) + s.scale(drag);
        let mut derivative = s.scale(2.0 * (1.0 + added_mass))
            + C64::from_re(image.damping().matrix_per_s()[0] + drag);
        for (row, tau) in branches(image) {
            let denominator = C64::ONE + s.scale(tau);
            derivative = derivative + C64::from_re(tau * row[0] * row[0])
                / (denominator * denominator);
        }
        let correction = value / derivative;
        s = s - correction;
        if correction.abs() < 1e-14 * w { return s; }
    }
    panic!("independent constrained-inertia dispersion did not converge")
}

#[test]
fn g1_descriptor_inertia_preserves_reactions_forcing_poles_and_original_material_heat() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w = bare.equilibrium_angular_frequencies_rad_s().next().unwrap();
    let beta = image.force_ports()[0][0];
    let mass = 0.75 / (beta * beta);
    let drag = 15.0 / (beta * beta);
    let input = [C64::new(0.7, -0.2), C64::ZERO, C64::ZERO, C64::ZERO];
    for map in [1.0, -1.0] {
        let prepared = image.prepare_coating_resonance_with_descriptor_load(
            ShellResonanceDescriptorLoad {
                system: inertance(mass, drag, map),
                ports: &[ShellHarmonicPort::Force(0)],
                band_rad_s: bare.source_band_rad_s(),
            },
            budget(), &gate,
        ).unwrap();
        let at = bare.dimension();
        assert_eq!(prepared.dimension(), at + 2);
        assert_eq!(prepared.load_state_count(), 2);
        assert_eq!(prepared.load_differential_range(), Some(at..at + 1));
        assert_eq!(prepared.load_algebraic_range(), Some(at + 1..at + 2));
        for omega in [800.0, 1100.0, 5000.0] {
            let s = C64::new(0.0, -omega);
            let q = input[0].scale(beta)
                / (body_operator(&image, s, w) + (s * s).scale(0.75) + s.scale(15.0));
            let velocity = s * q;
            let flow = velocity.scale(beta);
            let momentum = flow.scale(map * mass);
            let reaction = (s.scale(mass) + C64::from_re(drag)) * flow.scale(map);
            let factor = prepared.harmonic_factor(omega, forced_budget(), &gate).unwrap();
            let result = factor.response(&input, &gate).unwrap();
            source_channels(&image, &result, q, w);
            let load = result.load.as_ref().unwrap();
            near(load.state[0], momentum);
            near(load.state[1], reaction);
            near(load.flow[0], flow);
            near(load.effort[0], reaction.scale(map));
            near(load.state[0].scale(1.0 / mass), flow.scale(map));
            near_real(load.mean_stored_energy_j, 0.25 * momentum.norm_sq() / mass);
            near_real(load.dissipated_power_w, 0.5 * drag * flow.norm_sq());
            near_real(load.absorbed_port_power_w[0], load.dissipated_power_w);
            let rotated = input.map(|u| C64::new(-0.3, 0.4) * u);
            let next = factor.response(&rotated, &gate).unwrap();
            for (&a, &b) in next.physical.state.iter().zip(&result.physical.state) {
                near(a, C64::new(-0.3, 0.4) * b);
            }
        }
        let expected = inertance_pole(&image, w, 0.75, 15.0);
        let seed_rate = expected + C64::new(-0.01 * w, 0.02 * w);
        let mut seed = prepared.displacement_seed(seed_rate, &[C64::ONE], &gate).unwrap();
        let flow = seed_rate.scale(beta * map);
        seed[at] = flow.scale(mass);
        seed[at + 1] = (seed_rate.scale(mass) + C64::from_re(drag)) * flow;
        let pole = prepared.resonance_pole(seed_rate, &seed, &PoleOptions {
            spectral_scale: w,
            step_tolerance: 1e-13,
            residual_tolerance: 1e-14,
            max_iterations: 64,
            max_dimension: 128,
            max_dense_bytes: 16_000_000,
            max_work: 1_000_000_000,
        }, &gate).unwrap();
        near(pole.mode.laplace_pole_per_s, expected);
        let load = pole.load.as_ref().unwrap();
        let s = pole.mode.laplace_pole_per_s;
        let flow = pole.velocity_m_sqrt_kg_per_s[0].scale(beta);
        near(load.state[0], flow.scale(map * mass));
        near(load.state[1], (s.scale(mass) + C64::from_re(drag)) * flow.scale(map));
        near_real(load.stored_energy_j, 0.5 * load.state[0].norm_sq() / mass);
        near_real(load.dissipated_power_w, drag * flow.norm_sq());
        assert!(load.relative_equation_residual < 1e-10);
        assert!(pole.physical_relative_residual < 1e-10);
        assert!(pole.substrate_branch_dissipation_w.iter().sum::<f64>() > 0.0);
        assert!(pole.coating_branch_dissipation_w.iter().sum::<f64>() > 0.0);
    }
}

#[test]
fn g1_pure_algebraic_mutual_resistance_keeps_two_mode_memory_and_zero_reaction_storage() {
    let gate = CancelGate::new_clock_free();
    let image = image_with_modes(&gate, 2);
    let bare = image.prepare_coating_resonance(budget(), &gate).unwrap();
    let w: Vec<_> = bare.equilibrium_angular_frequencies_rad_s().collect();
    let scale: f64 = image.force_ports().iter().flatten().map(|v| v * v).sum::<f64>() / 20.0;
    let r = [2.0 * scale, 0.5 * scale, 0.5 * scale, scale];
    let z = [4.0 / (7.0 * scale), -2.0 / (7.0 * scale),
        -2.0 / (7.0 * scale), 8.0 / (7.0 * scale)];
    let system = DescriptorPortHamiltonian::new(
        2, 0, 2, vec![0.0; 4], r.to_vec(), vec![1.0, 0.0, 0.0, 1.0],
        Box::new(QuadraticStorage::new(vec![], 0).unwrap()),
    ).unwrap();
    let prepared = image.prepare_coating_resonance_with_descriptor_load(
        ShellResonanceDescriptorLoad {
            system,
            ports: &[ShellHarmonicPort::Force(0), ShellHarmonicPort::Force(1)],
            band_rad_s: bare.source_band_rad_s(),
        }, budget(), &gate,
    ).unwrap();
    let at = bare.dimension();
    assert_eq!(prepared.load_differential_range(), Some(at..at));
    assert_eq!(prepared.load_algebraic_range(), Some(at..at + 2));
    let input = [C64::new(0.7, -0.2), C64::new(-0.15, 0.3),
        C64::from_re(40.0), C64::from_re(-15.0)];
    let s = C64::new(0.0, -1500.0);
    let all = branches(&image);
    let mut matrix = [C64::ZERO; 4];
    let mut force = [C64::ZERO; 2];
    for i in 0..2 {
        force[i] = input[0].scale(image.force_ports()[0][i])
            + input[1].scale(image.force_ports()[1][i])
            + input[2].scale(image.volume_shapes_m2_per_sqrt_kg()[i])
            + input[3].scale(image.cavity_volume_shapes_m2_per_sqrt_kg()[i]);
        for j in 0..2 {
            let mut body = s.scale(image.damping().matrix_per_s()[2 * i + j]);
            if i == j { body = body + s * s + C64::from_re(w[i] * w[i]); }
            for (row, tau) in &all {
                body = body + s.scale(tau * row[i] * row[j]) / (C64::ONE + s.scale(*tau));
            }
            let mut resistance = 0.0;
            for a in 0..2 { for b in 0..2 {
                resistance += image.force_ports()[a][i] * z[2 * a + b]
                    * image.force_ports()[b][j];
            }}
            matrix[2 * i + j] = body + s.scale(resistance);
        }
    }
    let determinant = matrix[0] * matrix[3] - matrix[1] * matrix[2];
    let q = [(matrix[3] * force[0] - matrix[1] * force[1]) / determinant,
        (matrix[0] * force[1] - matrix[2] * force[0]) / determinant];
    let v = [s * q[0], s * q[1]];
    let flow: Vec<_> = image.force_ports().iter().map(|row|
        v[0].scale(row[0]) + v[1].scale(row[1])).collect();
    let reaction = [flow[0].scale(z[0]) + flow[1].scale(z[1]),
        flow[0].scale(z[2]) + flow[1].scale(z[3])];
    let result = prepared.harmonic_factor(1500.0, forced_budget(), &gate).unwrap()
        .response(&input, &gate).unwrap();
    let load = result.load.as_ref().unwrap();
    assert_eq!(load.mean_stored_energy_j, 0.0);
    let mut energy = 0.0;
    let mut heat = 0.0;
    for i in 0..2 {
        near(result.displacement_m_sqrt_kg[i], q[i]);
        near(load.state[i], reaction[i]);
        near(load.effort[i], reaction[i]);
        near(load.flow[i], flow[i]);
        near(reaction[0].scale(r[2 * i]) + reaction[1].scale(r[2 * i + 1]), flow[i]);
        energy += 0.25 * (v[i].norm_sq() + q[i].scale(w[i]).norm_sq());
        for j in 0..2 {
            heat += 0.5 * (v[i].conj() * v[j]).re * image.damping().matrix_per_s()[2 * i + j];
        }
    }
    for (i, (row, tau)) in all.iter().enumerate() {
        let extension = q[0].scale(row[0]) + q[1].scale(row[1]);
        let memory = extension / (C64::ONE + s.scale(*tau));
        near(result.physical.state[4 + i], memory);
        energy += 0.25 * (extension - memory).norm_sq();
        heat += 0.5 * (extension - memory).norm_sq() / tau;
    }
    let mut load_heat = 0.0;
    for i in 0..2 { for j in 0..2 {
        load_heat += 0.5 * (reaction[i].conj() * reaction[j]).re * r[2 * i + j];
    }}
    near_real(load.dissipated_power_w, load_heat);
    near_real(load.absorbed_port_power_w.iter().sum(), load_heat);
    near_real(result.physical.mean_stored_energy_j, energy);
    near_real(result.physical.dissipated_power_w, heat + load_heat);
    near_real(result.physical.supplied_power_w, heat + load_heat);
}

#[test]
fn g3_all_differential_descriptor_preserves_ode_response_and_original_coating_admission() {
    let gate = CancelGate::new_clock_free();
    let image = image(&gate);
    let make = || attachment(0.1, 100_000.0, 2.0, &[1.0, -0.6]);
    let original = make();
    let (j, r, g) = original.structure();
    let q = original.storage().quadratic_part().unwrap();
    let descriptor = DescriptorPortHamiltonian::new(
        original.state_dim(), original.state_dim(), original.port_dim(),
        j.to_vec(), r.to_vec(), g.to_vec(),
        Box::new(QuadraticStorage::new(q.matrix, original.state_dim()).unwrap()),
    ).unwrap();
    let ports = [ShellHarmonicPort::Force(0), ShellHarmonicPort::Force(1)];
    let ode = image.prepare_coating_resonance_with_load(
        ShellResonanceLoad { system: original, ports: &ports, band_rad_s: (900.0, 1500.0) },
        budget(), &gate,
    ).unwrap();
    let constrained = image.prepare_coating_resonance_with_descriptor_load(
        ShellResonanceDescriptorLoad { system: descriptor, ports: &ports, band_rad_s: (900.0, 1500.0) },
        budget(), &gate,
    ).unwrap();
    assert_eq!(constrained.load_algebraic_range(), Some(constrained.dimension()..constrained.dimension()));
    let a = ode.harmonic_factor(1100.0, forced_budget(), &gate).unwrap();
    let b = constrained.harmonic_factor(1100.0, forced_budget(), &gate).unwrap();
    let unit = [C64::ONE, C64::ZERO, C64::ZERO, C64::ZERO];
    let raw = b.response(&unit, &gate).unwrap();
    assert_eq!(a.response(&unit, &gate).unwrap(), raw);
    let mut observer = image.prepare_coating_strain(strain_budget(1), &gate).unwrap();
    let limit = observer.evaluate(&raw.displacement_m_sqrt_kg, &gate).unwrap()
        .sufficient_common_amplitude_scale.unwrap();
    let high = unit.map(|u| u.scale(1.2 * limit));
    let low = unit.map(|u| u.scale(0.6 * limit));
    let mut checked = constrained.harmonic_factor_with_coating_strain(
        1100.0, forced_budget(), strain_budget(2), &gate,
    ).unwrap();
    assert!(checked.response(&high, &gate).is_err());
    let accepted = checked.response(&low, &gate).unwrap();
    assert_eq!(accepted.shell, b.response(&low, &gate).unwrap());
    assert!(accepted.coating_strain.within_authored_envelope);
    assert!(std::ptr::eq(accepted.coating_strain.films[0].source, &image.source().coatings()[0]));
    assert_eq!(checked.completed_strain_evaluations(), 2);
}

#[test]
fn g0_assembled_descriptor_constraints_budgets_source_bands_and_cancellation_remain_explicit() {
    let gate = CancelGate::new_clock_free();
    let mut components = Vec::new();
    let mut rows = Vec::new();
    for i in 0..2 {
        let (shell, directors, constraints) = fixture::coupon(&gate);
        let shell = shell.with_coatings(vec![
            fixture::layer(CoatingSide::Positive, (20.0 + i as f64 * 5.0) * 1e-6),
        ], &gate).unwrap().with_coating_relaxation(CoatingRelaxationOptions {
            terms: 16, storage_relative_tolerance: 0.001, loss_relative_tolerance: 0.001,
            max_memory_coordinates: 256,
        }, &gate).unwrap().with_substrate_relaxation(SubstrateRelaxationOptions {
            storage_reference: SubstrateStorageReference::Relaxed, max_memory_coordinates: 256,
        }, &gate).unwrap();
        let mut row = vec![0.0; shell.mesh().nodes.len() * 6];
        row[6 * 5 + 2] = if i == 0 { 1.0 } else { -0.5 };
        rows.push(row);
        components.push(MaterialShellComponent {
            name: format!("original-{i}"), assembly_budget: fixture::budget(&shell),
            specimen: shell, directors, constraints, mass: ShellMass::Consistent, drilling_alpha: 1e-3,
        });
    }
    let source = MaterialShellAssembly {
        components, edge_ties: vec![], face_ties: vec![],
        tie_budget: TieBudget {
            max_segments: 32, max_slave_dofs: 24, gap_tolerance_m: 1e-10, angle_tolerance_rad: 1e-10,
        },
    };
    let image = source.reduce(&ShellAssemblyReductionOptions {
        eigenvalue_window: (0.0, 1e9), n_modes: 2,
        projection_budget: ShellProjectionBudget {
            max_modes: 2, max_projection_terms: 100_000_000, reference_resistance_n_s_m: 1.0,
        }, max_total_nodes: 32, max_total_triangles: 36,
    }, &[ShellAssemblyForcePort { components: rows }], &gate).unwrap();
    let make = || ShellResonanceDescriptorLoad {
        system: inertance(0.1, 2.0, 1.0),
        ports: &[ShellHarmonicPort::Force(0)], band_rad_s: (900.0, 1500.0),
    };
    let prepared = image.prepare_coating_resonance_with_descriptor_load(make(), budget(), &gate).unwrap();
    for limited in [
        ShellResonanceBudget { maximum_coordinates: prepared.dimension() - 1, ..budget() },
        ShellResonanceBudget { maximum_payload_bytes: 0, ..budget() },
        ShellResonanceBudget { maximum_preparation_work: 0, ..budget() },
    ] {
        assert!(image.prepare_coating_resonance_with_descriptor_load(make(), limited, &gate).is_err());
    }
    let mut invalid = make();
    invalid.ports = &[];
    assert!(image.prepare_coating_resonance_with_descriptor_load(invalid, budget(), &gate).is_err());
    for omega in [899.0, 1501.0] {
        assert!(prepared.harmonic_factor(omega, forced_budget(), &gate).is_err());
    }
    let factor = prepared.harmonic_factor(1100.0, forced_budget(), &gate).unwrap();
    let input = [C64::new(1e-4, -2e-4), C64::ZERO, C64::ZERO];
    let before = factor.response(&input, &gate).unwrap();
    let load = before.load.as_ref().unwrap();
    near(load.state[0].scale(10.0), load.flow[0]);
    near(load.effort[0], load.state[1]);
    near_real(load.mean_stored_energy_j, 2.5 * load.state[0].norm_sq());
    let mut checked = prepared.harmonic_factor_with_coating_strain(
        1100.0, forced_budget(), strain_budget(1), &gate,
    ).unwrap();
    let result = checked.response(&input, &gate).unwrap();
    assert_eq!(result.shell, before);
    assert_eq!(result.coating_strain.films.len(), 2);
    for (i, film) in result.coating_strain.films.iter().enumerate() {
        assert_eq!(film.component_index, Some(i));
        assert!(std::ptr::eq(film.source, &image.source().components[i].specimen.coatings()[0]));
    }
    let redundant = DescriptorPortHamiltonian::new(
        2, 0, 1, vec![0.0; 4], vec![0.0; 4], vec![1.0, 1.0],
        Box::new(QuadraticStorage::new(vec![], 0).unwrap()),
    ).unwrap();
    let singular = image.prepare_coating_resonance_with_descriptor_load(
        ShellResonanceDescriptorLoad {
            system: redundant, ports: &[ShellHarmonicPort::Force(0)], band_rad_s: (900.0, 1500.0),
        }, budget(), &gate,
    ).unwrap();
    assert!(singular.harmonic_factor(1100.0, forced_budget(), &gate).is_err());
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(factor.response(&input, &stopped).is_err());
    assert!(image.prepare_coating_resonance_with_descriptor_load(make(), budget(), &stopped).is_err());
    assert_eq!(factor.response(&input, &gate).unwrap(), before);
}
