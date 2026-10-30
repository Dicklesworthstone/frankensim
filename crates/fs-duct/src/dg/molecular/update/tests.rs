//! G0/G2: physical entropy/compression transfer and continuing native dynamics.
use super::*;
use crate::profile::gas_field::tests::field;
use fs_material::relaxation::{
    InternalReservoir, RelaxationProcess, RelaxationRateLaw, SpeciesRelaxation, VibrationalMode,
};
use fs_material::visco::GeneralizedMaxwell;
use fs_phs::WallPin;
use fs_thermochem::SpeciesId;

pub(super) fn declaration(active: bool) -> RelaxationAbsorption {
    RelaxationAbsorption::new(
        ["analytic-field-a", "analytic-field-b"]
            .into_iter()
            .enumerate()
            .map(|(i, name)| SpeciesRelaxation {
                species: SpeciesId::new(name).unwrap(),
                processes: if active {
                    vec![RelaxationProcess {
                        reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                            theta_k: 800.0 + 600.0 * i as f64,
                            degeneracy: 1,
                        }]),
                        rate: RelaxationRateLaw::PartnerIndependent {
                            p_tau_ref_pa_s: 2.0 + i as f64 * 8.0,
                            t_ref_k: 300.0,
                            exponent: 1.2,
                            t_min_k: 200.0,
                            t_max_k: 500.0,
                        },
                    }]
                } else {
                    vec![]
                },
            })
            .collect(),
    )
    .unwrap()
}

fn base(source: &AxialGasField, diffusion: bool) -> DgDuct {
    let duct = DgDuct::new_with_gas_field(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.0005,
                length: source.length_m(),
            }],
        },
        source,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    if diffusion {
        duct.with_viscothermal(8, [50.0, 2000.0], 0.02).unwrap()
    } else {
        duct
    }
}

fn source(temperature: f64) -> AxialGasField {
    field([temperature; 2], [0.7; 2], 100_000.0)
}

fn near(actual: f64, expected: f64, scale: f64) {
    assert!(
        (actual - expected).abs() <= 2e-10 * scale.max(1e-30),
        "actual {actual:e}, expected {expected:e}, scale {scale:e}"
    );
}

// Reconstruct temperature, entropy and compression from the independent
// linear ideal-gas equations in PHYSICAL nodal units. Radial weights come
// from the retained Bessel spectral measure, not the refresh's coefficients.
fn physical_invariants(duct: &DgDuct, state: &DgState) -> Vec<f64> {
    let material = duct.molecular_relaxation().unwrap();
    let gas = duct.gas();
    let cp = material.frozen_cv_j_kg_k() + gas.specific_gas_constant;
    let bulk = gas.density * material.frozen_sound_speed_m_s().powi(2);
    let memory = duct.molecular.as_ref().unwrap();
    let weights: Vec<_> = match duct.viscothermal_model() {
        Some(model) => model
            .shunt_terms()
            .iter()
            .map(|term| term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0)))
            .collect(),
        None => vec![1.0],
    };
    let capacities: Vec<_> = material
        .processes()
        .iter()
        .filter_map(|p| (p.heat_capacity_j_kg_k > 0.0).then_some(p.heat_capacity_j_kg_k))
        .collect();
    let m = duct.field_dofs;
    let (pressure, flow) = duct.from_energy_coords(&state.tilde);
    let mut entropy = Vec::new();
    let mut mean_temperature = vec![0.0; m];
    let mass_inverse = |values: &[f64]| {
        let local = duct.mesh.order + 1;
        let mut result = vec![0.0; m];
        for (e, element) in duct.elements.iter().enumerate() {
            for i in 0..local {
                result[e * local + i] = (0..local)
                    .map(|j| element.mp_inv_sqrt[i * local + j] * values[e * local + j])
                    .sum::<f64>()
                    / bulk.sqrt();
            }
        }
        result
    };
    for (mode, &weight) in memory.modes.iter().zip(&weights) {
        let temperature_coordinate: Vec<_> = (0..m)
            .map(|i| {
                state.tilde[mode.start + i]
                    + (weight * gas.specific_gas_constant / material.frozen_cv_j_kg_k()).sqrt()
                        * state.tilde[i]
            })
            .collect();
        let temperature = mass_inverse(&temperature_coordinate);
        for i in 0..m {
            let theta = temperature[i] * (gas.temperature / (gas.density * cp)).sqrt();
            mean_temperature[i] += weight.sqrt() * theta;
            entropy
                .push((gas.density * cp * theta - weight.sqrt() * pressure[i]) / gas.temperature);
        }
        for (j, &capacity) in capacities.iter().enumerate() {
            let offset = mode.start + (j + 1) * m;
            let temperature = mass_inverse(&state.tilde[offset..offset + m]);
            for theta in temperature {
                let theta = theta * (gas.temperature / (gas.density * capacity)).sqrt();
                entropy.push(gas.density * capacity * theta / gas.temperature);
            }
        }
    }
    let mut result: Vec<_> = pressure
        .iter()
        .zip(&mean_temperature)
        .map(|(p, theta)| p / gas.pressure - theta / gas.temperature)
        .collect();
    result.extend(entropy);
    for (element, flows) in duct.elements.iter().zip(flow.chunks(duct.mesh.order + 1)) {
        result.extend(
            flows
                .iter()
                .zip(&element.areas)
                .map(|(q, area)| gas.density * q / area),
        );
    }
    // Every series arm's independent inertial momentum in normalized geometry.
    let series = duct
        .viscothermal_model()
        .map_or(0, |model| model.series_terms().len());
    result.extend(
        state.tilde[2 * m..(2 + series) * m]
            .iter()
            .map(|x| x * gas.density.sqrt()),
    );
    result
}

pub(super) fn excited(duct: &DgDuct) -> DgState {
    let mut state = duct.project_initial_state(|x| 10.0 + 3.0 * x, |x| 1e-7 * (1.0 - x));
    let step = duct.propagator(2e-6, 1e-13).unwrap();
    for _ in 0..12 {
        step.step(&mut state, &[]).unwrap();
    }
    step.sync_nodal(&mut state).unwrap();
    assert!(energy(&state.tilde[2 * duct.field_dofs..]) > 0.0);
    state
}

#[test]
fn molecular_refresh_preserves_physical_entropy_compression_and_all_momenta() {
    for diffusion in [false, true] {
        let original_source = source(300.0);
        let replacement = source(410.0);
        let law = declaration(true);
        let mut duct = base(&original_source, diffusion)
            .with_molecular_relaxation_field(law.clone())
            .unwrap();
        assert_eq!(duct.molecular_relaxation_law(), Some(&law));
        assert!(duct.molecular_relaxation_field().is_none());
        let mut state = excited(&duct);
        let original = duct.clone();
        let before = state.clone();
        let physical = physical_invariants(&duct, &state);
        let report = duct
            .rebind_uniform_molecular_gas_field(&mut state, &replacement, 0.02, 1.0, &mut || false)
            .unwrap();
        let refreshed = physical_invariants(&duct, &state);
        for (old, new) in physical.iter().zip(&refreshed) {
            near(*new, *old, old.abs().max(1e-10));
        }
        assert_eq!(duct.gas_field(), Some(&replacement));
        assert_eq!(duct.molecular_relaxation_law(), Some(&law));
        assert_eq!(state.step_count(), before.step_count());
        assert_eq!(state.initial_energy_j(), before.initial_energy_j());
        assert_eq!(
            state.cumulative_dissipated_j(),
            before.cumulative_dissipated_j()
        );
        near(
            report.energy_after_j - report.energy_before_j,
            report.parameter_work_j,
            report.energy_before_j + report.energy_after_j,
        );
        assert_eq!(state.cumulative_parameter_work_j(), report.parameter_work_j);
        // A pressure-only p/B transfer would erase this independently live
        // caloric correction from the nonzero thermal histories.
        let naive =
            before.pressure()[0] * duct.pressure_bulk_modulus() / original.pressure_bulk_modulus();
        assert!((state.pressure()[0] - naive).abs() > 1e-5 * naive.abs());

        let update = duct
            .rebind_uniform_molecular_gas_field(
                &mut state,
                &original_source,
                0.02,
                1.0,
                &mut || false,
            )
            .unwrap();
        for (old, new) in before.tilde.iter().zip(&state.tilde) {
            near(*new, *old, before.stored_energy_j().sqrt());
        }
        near(
            report.parameter_work_j + update.parameter_work_j,
            0.0,
            before.stored_energy_j(),
        );
    }
}

#[test]
fn molecular_refresh_continues_the_original_wall_and_energy_ledger() {
    let original_source = source(300.0);
    let pin = WallPin {
        surface_density: 0.4,
        stiffness_per_area: 4e8,
        resistance: 30.0,
    };
    let wall_law = GeneralizedMaxwell::new(2e6, vec![(8e5, 0.001), (3e5, 0.0002)]).unwrap();
    let mut duct = base(&original_source, true)
        .with_molecular_relaxation_field(declaration(true))
        .unwrap()
        .with_compliant_wall_relaxation(pin, &wall_law)
        .unwrap();
    let mut state = excited(&duct);
    let old_propagator = duct.propagator(2e-6, 1e-13).unwrap();
    let before = state.clone();
    let motion = duct.compliant_wall_state(&state).unwrap();
    let last_mode = duct.molecular.as_ref().unwrap().modes.last().unwrap();
    let wall_start = last_mode.start + (last_mode.rates.len() + 1) * duct.field_dofs;
    let report = duct
        .rebind_uniform_molecular_gas_field(&mut state, &source(360.0), 0.02, 1.0, &mut || false)
        .unwrap();
    assert_eq!(&state.tilde[wall_start..], &before.tilde[wall_start..]);
    let after = duct.compliant_wall_state(&state).unwrap();
    for (old, new) in motion
        .displacement_m
        .iter()
        .zip(&after.displacement_m)
        .chain(motion.velocity_m_s.iter().zip(&after.velocity_m_s))
        .chain(
            motion
                .relaxation_pressure_pa
                .iter()
                .flatten()
                .zip(after.relaxation_pressure_pa.iter().flatten()),
        )
    {
        near(*new, *old, old.abs().max(1e-12));
    }
    assert!(
        after
            .relaxation_pressure_pa
            .iter()
            .flatten()
            .any(|p| p.abs() > 0.0)
    );
    let checkpoint = state.clone();
    assert!(old_propagator.step(&mut state, &[]).is_err());
    assert_eq!(state, checkpoint);
    let stepping = duct.propagator(2e-6, 1e-13).unwrap();
    let mut previous_energy = state.stored_energy_j();
    for _ in 0..20 {
        let output = stepping.step(&mut state, &[]).unwrap();
        assert!(state.stored_energy_j() <= previous_energy + 1e-12 * before.initial_energy_j());
        previous_energy = state.stored_energy_j();
        assert_eq!(
            output.ledger.cumulative_parameter_work_j,
            report.parameter_work_j
        );
        near(
            output.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
        assert!(duct.molecular_dissipated_power_w(&state).unwrap() >= 0.0);
    }
}

#[test]
fn molecular_refresh_noop_is_exact_and_empty_declarations_keep_the_original_transfer() {
    let original_source = source(300.0);
    for diffusion in [false, true] {
        let mut duct = base(&original_source, diffusion)
            .with_molecular_relaxation_field(declaration(true))
            .unwrap();
        let mut state = excited(&duct);
        let original = duct.clone();
        let before = state.clone();
        let report = duct
            .rebind_uniform_molecular_gas_field(
                &mut state,
                &original_source,
                0.02,
                0.0,
                &mut || false,
            )
            .unwrap();
        assert_eq!(report.parameter_work_j, 0.0);
        assert_eq!(state.tilde, before.tilde);
        assert_eq!(duct.l_tilde, original.l_tilde);
        assert_eq!(duct.operator_identity(), original.operator_identity());

        let replacement = source(380.0);
        let mut plain = base(&original_source, diffusion);
        let mut marked = plain
            .clone()
            .with_molecular_relaxation_field(declaration(false))
            .unwrap();
        let mut plain_state = plain.project_initial_state(|_| 11.0, |_| 2e-7);
        let mut marked_state = marked.project_initial_state(|_| 11.0, |_| 2e-7);
        let expected = if diffusion {
            plain
                .rebind_uniform_viscothermal_gas_field(
                    &mut plain_state,
                    &replacement,
                    0.02,
                    1.0,
                    &mut || false,
                )
                .unwrap()
        } else {
            plain
                .rebind_uniform_gas_field(&mut plain_state, &replacement, 1.0, &mut || false)
                .unwrap()
        };
        let measured = marked
            .rebind_uniform_molecular_gas_field(
                &mut marked_state,
                &replacement,
                0.02,
                1.0,
                &mut || false,
            )
            .unwrap();
        assert_eq!(measured, expected);
        assert_eq!(plain_state, marked_state);
        assert_eq!(plain.l_tilde, marked.l_tilde);
        assert!(
            marked
                .molecular_relaxation()
                .unwrap()
                .processes()
                .is_empty()
        );
        assert_eq!(marked.molecular_relaxation_law(), Some(&declaration(false)));
    }
}

#[test]
fn molecular_refresh_refusals_and_cancellation_preserve_every_original_history() {
    let original_source = source(300.0);
    let duct = base(&original_source, true)
        .with_molecular_relaxation_field(declaration(true))
        .unwrap();
    let state = excited(&duct);
    let replacement = source(400.0);
    for (field, error, budget) in [
        (field([400.0; 2], [0.4; 2], 100_000.0), 0.02, 1.0),
        (field([300.0, 400.0], [0.7; 2], 100_000.0), 0.02, 1.0),
        (replacement.clone(), 0.02, 0.0),
        (replacement.clone(), 0.0, 1.0),
    ] {
        let mut candidate = duct.clone();
        let mut next = state.clone();
        assert!(candidate.rebind_uniform_molecular_gas_field(
            &mut next, &field, error, budget, &mut || false,
        ).is_err());
        assert_eq!(candidate.operator_identity(), duct.operator_identity());
        assert_eq!(candidate.gas_field(), duct.gas_field());
        assert_eq!(
            candidate.molecular_relaxation_law(),
            duct.molecular_relaxation_law()
        );
        assert_eq!(next, state);
    }
    let mut candidate = duct.clone();
    let mut next = state.clone();
    let mut calls = 0;
    assert!(
        candidate
            .rebind_uniform_molecular_gas_field(&mut next, &replacement, 0.02, 1.0, &mut || {
                calls += 1;
                calls > 8
            },)
            .is_err()
    );
    assert!(calls > 8);
    assert_eq!(candidate.operator_identity(), duct.operator_identity());
    assert_eq!(next, state);
    let mut unbound = base(&original_source, true)
        .with_molecular_relaxation(duct.molecular_relaxation().unwrap().clone())
        .unwrap();
    let mut unbound_state = unbound.zero_state();
    assert!(
        unbound
            .rebind_uniform_molecular_gas_field(
                &mut unbound_state,
                &replacement,
                0.02,
                1.0,
                &mut || false,
            )
            .is_err()
    );
}

// Independent ideal-gas weak density and entropy equations on a conical
// section. Three Gauss points integrate A(x)*phi_i*phi_j exactly here.
// Temperatures use the declared gas capacities and physical quadrature volume,
// never the production pressure projection or refresh coordinate map.
fn spatial_physical_invariants(duct: &DgDuct, state: &DgState) -> Vec<f64> {
    let gas = duct.gas();
    let material = duct.molecular_relaxation().unwrap();
    let cp = material.frozen_cv_j_kg_k() + gas.specific_gas_constant;
    let capacities: Vec<_> = material
        .processes()
        .iter()
        .filter_map(|p| (p.heat_capacity_j_kg_k > 0.0).then_some(p.heat_capacity_j_kg_k))
        .collect();
    let resolution = duct.spatial_viscothermal_admission().unwrap().resolution;
    let (nodes, weights) = gauss_legendre_rule(resolution.quadrature_points);
    let modes = &duct.molecular.as_ref().unwrap().modes;
    let local = duct.mesh.order + 1;
    let (pressure, flow) = duct.from_energy_coords(&state.tilde);
    let mut weak_compression = vec![0.0; duct.field_dofs];
    let mut entropy = Vec::new();
    let mut at = 0;
    for (e, element) in duct.elements.iter().enumerate() {
        for (&xi, &weight) in nodes.iter().zip(&weights) {
            let x = element.x_left + 0.5 * element.length * (xi + 1.0);
            let radius = 0.0005 + 0.0003 * x / duct.total_length();
            let volume = 0.5 * element.length * weight * core::f64::consts::PI * radius * radius;
            let phi = duct.lgl.eval_basis(xi);
            let p: f64 = phi
                .iter()
                .zip(&pressure[e * local..(e + 1) * local])
                .map(|(phi, p)| phi * p)
                .sum();
            let spectral = fs_phs::ZkDiffusive::new(
                radius,
                gas.density,
                gas.dynamic_viscosity,
                gas.gamma,
                gas.prandtl,
                resolution.poles_per_field,
                resolution.band_rad_s,
                resolution.maximum_relative_error,
            )
            .unwrap();
            let mut mean_temperature = 0.0;
            for term in spectral.shunt_terms() {
                let radial = term.gain_per_s / (term.rate_per_s * (gas.gamma - 1.0));
                let mode = &modes[at];
                at += 1;
                assert_eq!(mode.projection.width(duct.field_dofs), 1);
                let temperature_coordinate = state.tilde[mode.start]
                    + (radial * volume / (gas.density * cp * gas.temperature)).sqrt() * p;
                let temperature =
                    temperature_coordinate * (gas.temperature / (volume * gas.density * cp)).sqrt();
                mean_temperature += radial.sqrt() * temperature;
                entropy
                    .push((gas.density * cp * temperature - radial.sqrt() * p) / gas.temperature);
                for (j, &capacity) in capacities.iter().enumerate() {
                    let temperature = state.tilde[mode.start + j + 1]
                        * (gas.temperature / (volume * gas.density * capacity)).sqrt();
                    entropy.push(gas.density * capacity * temperature / gas.temperature);
                }
            }
            for (i, phi) in phi.iter().enumerate() {
                weak_compression[e * local + i] +=
                    volume * phi * (p / gas.pressure - mean_temperature / gas.temperature);
            }
        }
    }
    assert_eq!(at, modes.len());
    let mut physical = weak_compression;
    physical.extend(entropy);
    for (element, values) in duct.elements.iter().zip(flow.chunks(local)) {
        physical.extend(
            values
                .iter()
                .zip(&element.areas)
                .map(|(q, area)| gas.density * q / area),
        );
    }
    physical.extend(
        state.tilde[2 * duct.field_dofs..modes[0].start]
            .iter()
            .map(|value| gas.density.sqrt() * value),
    );
    physical
}

#[test]
fn molecular_refresh_retains_spatial_wall_histories_and_independent_weak_equation_of_state() {
    let original = source(300.0);
    let shape = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.0005,
            outlet_radius: 0.0008,
            length: original.length_m(),
        }],
    };
    let mut duct = DgDuct::new_with_gas_field(
        &shape,
        &original,
        DgMesh {
            elements: 1,
            order: 1,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap()
    .with_spatial_viscothermal(DgWallResolution {
        poles_per_field: 8,
        quadrature_points: 3,
        band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
        maximum_relative_error: 0.02,
    })
    .unwrap()
    .with_molecular_relaxation_field(declaration(true))
    .unwrap();
    let mut state = excited(&duct);
    let before = state.clone();
    let operator = duct.l_tilde.clone();
    let no_op = duct
        .rebind_uniform_molecular_gas_field(&mut state, &original, 0.02, 0.0, &mut || false)
        .unwrap();
    assert_eq!(no_op.parameter_work_j, 0.0);
    assert_eq!(state.tilde, before.tilde);
    assert_eq!(duct.l_tilde, operator);
    let physical = spatial_physical_invariants(&duct, &state);
    let resolution = duct.spatial_viscothermal_admission().unwrap().resolution;
    let report = duct
        .rebind_uniform_molecular_gas_field(&mut state, &source(410.0), 0.02, 1.0, &mut || false)
        .unwrap();
    assert_eq!(
        duct.spatial_viscothermal_admission().unwrap().resolution,
        resolution
    );
    let after = spatial_physical_invariants(&duct, &state);
    for (old, new) in physical.iter().zip(after) {
        near(new, *old, old.abs().max(1e-20));
    }
    near(
        report.energy_after_j - report.energy_before_j,
        report.parameter_work_j,
        report.energy_before_j + report.energy_after_j,
    );
    let stepping = duct.propagator(2e-6, 1e-13).unwrap();
    for _ in 0..8 {
        let step = stepping.step(&mut state, &[]).unwrap();
        near(
            step.ledger.energy_balance_residual_j,
            0.0,
            before.initial_energy_j(),
        );
    }

    // An explicitly absent species stays in the original law and field. The
    // material owner omits its evaluated process; refresh preserves that exact
    // support instead of allocating an artificial zero-capacity coordinate.
    let zero_source = field([300.0; 2], [1.0; 2], 100_000.0);
    let mut zero = base(&zero_source, true)
        .with_molecular_relaxation_field(declaration(true))
        .unwrap();
    let mut zero_state = excited(&zero);
    let old_zero: Vec<_> = zero
        .molecular_relaxation()
        .unwrap()
        .processes()
        .iter()
        .map(|p| (p.species.clone(), p.index, p.heat_capacity_j_kg_k == 0.0))
        .collect();
    let absent = SpeciesId::new("analytic-field-b").unwrap();
    assert!(old_zero.iter().all(|(species, _, _)| species != &absent));
    assert!(
        zero.molecular_relaxation_law()
            .unwrap()
            .species()
            .iter()
            .any(|entry| entry.species == absent && !entry.processes.is_empty())
    );
    zero.rebind_uniform_molecular_gas_field(
        &mut zero_state,
        &field([380.0; 2], [1.0; 2], 100_000.0),
        0.02,
        1.0,
        &mut || false,
    )
    .unwrap();
    let new_zero: Vec<_> = zero
        .molecular_relaxation()
        .unwrap()
        .processes()
        .iter()
        .map(|p| (p.species.clone(), p.index, p.heat_capacity_j_kg_k == 0.0))
        .collect();
    assert_eq!(old_zero, new_zero);
    assert_eq!(zero.molecular_relaxation_law(), Some(&declaration(true)));
}
