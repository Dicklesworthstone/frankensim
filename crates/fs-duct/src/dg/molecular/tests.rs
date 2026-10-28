//! Synthetic thermodynamic fixtures test the equations, not measured gas data.
use super::*;
use fs_la::factor::lu;
use fs_material::gas::{GasSpec, GasState};
use fs_material::relaxation::{
    InternalReservoir, RelaxationAbsorption, RelaxationProcess, RelaxationRateLaw, RelaxingMedium,
    SpeciesRelaxation, VibrationalMode,
};
use fs_math::c64::C64;
use fs_thermochem::SpeciesId;

fn gas() -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 101325.0).unwrap()
}

fn law(gas: GasState, count: usize) -> RelaxationState {
    law_rate(gas, count, 1.0)
}

fn law_rate(gas: GasState, count: usize, rate_scale: f64) -> RelaxationState {
    let id = SpeciesId::new("synthetic").unwrap();
    let species = SpeciesRelaxation {
        species: id.clone(),
        processes: (0..count)
            .map(|i| RelaxationProcess {
                reservoir: InternalReservoir::Vibrational(vec![VibrationalMode {
                    theta_k: 1300.0 + i as f64 * 600.0,
                    degeneracy: 1,
                }]),
                rate: RelaxationRateLaw::PartnerIndependent {
                    p_tau_ref_pa_s: gas.pressure * (2e-5 + i as f64 * 3e-4) * rate_scale,
                    t_ref_k: gas.temperature,
                    exponent: 0.0,
                    t_min_k: 200.0,
                    t_max_k: 500.0,
                },
            })
            .collect(),
    };
    RelaxationAbsorption::new(vec![species])
        .unwrap()
        .evaluate(&RelaxingMedium::new(gas, vec![(id, 1.0)]).unwrap())
        .unwrap()
}

fn duct(walls: bool) -> DgDuct {
    let duct = DgDuct::new(
        &Duct {
            segments: vec![Segment::Cylinder {
                radius: 0.001,
                length: 0.05,
            }],
        },
        &gas(),
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::Closed,
    )
    .unwrap();
    if walls {
        duct.with_viscothermal(8, [50.0, 2000.0], 0.02).unwrap()
    } else {
        duct
    }
}

// Eliminate the ACTUAL local energy matrix using an independent dense LU.
// No production thermal-observation helper supplies the response.
fn measured_compliance(duct: &DgDuct, omega: f64) -> C64 {
    let memory = duct.molecular.as_ref().unwrap();
    let mut rows = vec![0];
    for mode in &memory.modes {
        for j in 0..=memory.rates.len() {
            rows.push(mode.start + j * duct.field_dofs);
        }
    }
    let n = rows.len();
    let mut matrix = vec![0.0; 4 * n * n];
    let stride = 2 * n;
    for i in 0..n {
        for j in 0..n {
            let real = -duct.l_tilde[rows[i] * duct.total_dofs() + rows[j]];
            matrix[i * stride + j] = real;
            matrix[(n + i) * stride + n + j] = real;
        }
        matrix[i * stride + n + i] = omega;
        matrix[(n + i) * stride + i] = -omega;
    }
    let mut response = vec![0.0; stride];
    response[0] = 1.0;
    lu(&matrix, stride).unwrap().solve(&mut response);
    (C64::new(0.0, -omega) * C64::new(response[0], response[n]).scale(duct.pressure_bulk_modulus()))
        .recip()
}

// Independent heat-equation elimination with the ORIGINAL Bessel weights,
// evaluated heat capacities and times; no augmented matrix coefficients.
fn expected_compliance(base: &DgDuct, state: &RelaxationState, omega: f64) -> C64 {
    let s = C64::new(0.0, -omega);
    let gas = state.gas();
    let cp = state.processes().iter().fold(
        C64::from_re(state.frozen_cv_j_kg_k() + gas.specific_gas_constant),
        |sum, p| {
            sum + C64::from_re(p.heat_capacity_j_kg_k) / (C64::ONE + s.scale(p.relaxation_time_s))
        },
    );
    let temperature = match base.viscothermal_model() {
        Some(model) => model.shunt_terms().iter().fold(C64::ZERO, |sum, t| {
            let weight = t.gain_per_s / (t.rate_per_s * (gas.gamma - 1.0));
            sum + s.scale(weight) / (s * cp + C64::from_re(gas.specific_heat_cp * t.rate_per_s))
        }),
        None => cp.recip(),
    };
    (C64::ONE - temperature.scale(gas.specific_gas_constant))
        .scale(1.0 / (gas.density * gas.specific_gas_constant * gas.temperature))
}

#[test]
fn molecular_pressure_memory_matches_shared_material_and_coupled_heat_equation() {
    for walls in [false, true] {
        let base = duct(walls);
        let state = law(gas(), 2);
        let resolved = base
            .clone()
            .with_molecular_relaxation(state.clone())
            .unwrap();
        assert_eq!(resolved.gas(), base.gas());
        for omega in [0.1, 10.0, 300.0, 3000.0, 30000.0, 1e6] {
            let actual = measured_compliance(&resolved, omega);
            let expected = expected_compliance(&base, &state, omega);
            assert!(
                (actual - expected).abs() < 2e-9 * expected.abs(),
                "walls={walls} omega={omega}: {actual:?} != {expected:?}"
            );
            if !walls {
                let k = state.relaxing_wavenumber(omega).unwrap();
                let material = (k * k).scale(1.0 / (gas().density * omega * omega));
                assert!((actual - material).abs() < 2e-9 * material.abs());
            }
        }
        if walls {
            let static_compliance =
                1.0 / (gas().density * gas().specific_gas_constant * gas().temperature);
            let low = measured_compliance(&resolved, 0.1);
            assert!(
                (low.re / static_compliance - 1.0).abs() < 1e-5,
                "combined thermal+molecular law must approach isothermal storage once"
            );
        }
    }
}

#[test]
fn molecular_and_wall_losses_close_the_original_operator_and_exact_step() {
    let resolved = duct(true)
        .with_molecular_relaxation(law(gas(), 2))
        .unwrap()
        .with_compliant_wall(fs_phs::WallPin {
            surface_density: 0.4,
            stiffness_per_area: 8e8,
            resistance: 100.0,
        })
        .unwrap();
    let mut arbitrary = resolved.zero_state();
    for (i, x) in arbitrary.tilde.iter_mut().enumerate() {
        *x = 1e-3 * det::sin(i as f64 + 0.3);
    }
    let n = resolved.total_dofs();
    let mut operator_loss = 0.0;
    for i in 0..n {
        for j in 0..n {
            operator_loss -= arbitrary.tilde[i] * resolved.l_tilde[i * n + j] * arbitrary.tilde[j];
        }
    }
    let molecular = resolved.molecular_dissipated_power_w(&arbitrary).unwrap();
    let wall = resolved.wall_dissipated_power_w(&arbitrary).unwrap();
    assert!(molecular > 0.0 && wall > 0.0);
    assert!((operator_loss - molecular - wall).abs() < 1e-11 * (molecular + wall));

    let mut coarse = resolved.project_initial_state(|x| 10.0 * det::cos(20.0 * x), |_| 0.0);
    let mut fine = coarse.clone();
    let dt = 2e-6;
    let large = resolved.propagator(dt, 1e-11).unwrap();
    let small = resolved.propagator(dt * 0.5, 1e-11).unwrap();
    let initial = coarse.stored_energy_j();
    let mut previous = initial;
    for _ in 0..8 {
        large.step(&mut coarse, &[]).unwrap();
        small.step(&mut fine, &[]).unwrap();
        small.step(&mut fine, &[]).unwrap();
        assert!(coarse.stored_energy_j() <= previous * (1.0 + 1e-10));
        previous = coarse.stored_energy_j();
    }
    let difference: f64 = coarse
        .tilde
        .iter()
        .zip(&fine.tilde)
        .map(|(a, b)| (a - b) * (a - b))
        .sum();
    assert!(difference < 1e-16 * initial);
    large.sync_nodal(&mut coarse).unwrap();
    assert!(coarse.pressure().iter().all(|p| p.is_finite()));
}

#[test]
fn molecular_cold_attachment_preserves_classical_limit_and_rejects_wrong_histories() {
    let base = duct(true);
    let classical = base
        .clone()
        .with_molecular_relaxation(law(gas(), 0))
        .unwrap();
    assert_eq!(classical.l_tilde, base.l_tilde);
    assert_eq!(classical.operator_identity(), base.operator_identity());
    assert_eq!(classical.total_dofs(), base.total_dofs());
    let wrong_gas = GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, 202650.0).unwrap();
    assert!(
        base.clone()
            .with_molecular_relaxation(law(wrong_gas, 1))
            .is_err()
    );
    let resolved = base
        .clone()
        .with_molecular_relaxation(law(gas(), 2))
        .unwrap();
    assert!(
        resolved
            .clone()
            .with_molecular_relaxation(law(gas(), 2))
            .is_err()
    );
    assert!(
        resolved
            .molecular_dissipated_power_w(&base.zero_state())
            .is_err()
    );
    assert!(
        resolved
            .clone()
            .with_viscothermal(8, [50.0, 2000.0], 0.02)
            .is_err()
    );
    let mut live = resolved.clone();
    let mut state = live.zero_state();
    let before = state.clone();
    assert!(
        live.rebind_uniform_viscothermal_gas(&mut state, gas(), 0.02, 1.0, &mut || false)
            .is_err()
    );
    assert_eq!(before, state);
    assert!(
        base.with_molecular_relaxation_controlled(law(gas(), 2), &mut || true)
            .is_err()
    );
}

fn complex_sin_cos(z: C64) -> (C64, C64) {
    let (s, c) = det::sin_cos(z.re);
    let positive = det::exp(z.im);
    let negative = 1.0 / positive;
    let sinh = 0.5 * (positive - negative);
    let cosh = 0.5 * (positive + negative);
    (C64::new(s * cosh, c * sinh), C64::new(c * cosh, -s * sinh))
}

#[test]
fn molecular_bore_pressure_converges_and_responds_to_material_rate_data() {
    let geometry = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.05,
        }],
    };
    let material = law(gas(), 2);
    let mut errors = Vec::new();
    for elements in [1, 2] {
        let base = DgDuct::new(
            &geometry,
            &gas(),
            DgMesh { elements, order: 2 },
            DgBoundary::Port(PortId(7)),
            DgBoundary::IdealOpen,
        )
        .unwrap()
        .with_viscothermal(8, [50.0, 2000.0], 0.02)
        .unwrap();
        let resolved = base
            .clone()
            .with_molecular_relaxation(material.clone())
            .unwrap();
        let swapped = base
            .clone()
            .with_molecular_relaxation(law_rate(gas(), 2, 0.1))
            .unwrap();
        let mut worst = 0.0_f64;
        let mut rate_effect = 0.0_f64;
        let area = core::f64::consts::PI * 0.001_f64.powi(2);
        for f in [100.0, 500.0, 1000.0] {
            let omega = core::f64::consts::TAU * f;
            let s = C64::new(0.0, -omega);
            let series = base
                .viscothermal_model()
                .unwrap()
                .response_at(omega)
                .unwrap()
                .0
                .conj()
                .scale(gas().density / area);
            let shunt = s * expected_compliance(&base, &material, omega).scale(area);
            let propagation = (series * shunt).sqrt();
            let impedance = series / propagation;
            let k = C64::new(0.0, 1.0) * propagation;
            let (sin, cos) = complex_sin_cos(k.scale(0.05));
            let expected = C64::new(0.0, -1.0) * impedance * sin / cos;
            let response = resolved
                .harmonic_response(
                    omega,
                    &[(PortId(7), C64::ONE)],
                    1_000_000_000,
                    1e-10,
                    &mut || false,
                )
                .unwrap();
            let changed = swapped
                .harmonic_response(
                    omega,
                    &[(PortId(7), C64::ONE)],
                    1_000_000_000,
                    1e-10,
                    &mut || false,
                )
                .unwrap();
            let actual = response.ports[0].pressure_pa;
            worst = worst.max((actual - expected).abs() / expected.abs());
            rate_effect =
                rate_effect.max((actual - changed.ports[0].pressure_pa).abs() / actual.abs());
            assert!(response.molecular_dissipation_w > 0.0 && response.wall_dissipation_w > 0.0);
            assert!(response.power_residual_w.abs() < 1e-7 * response.input_power_w.abs());
        }
        assert!(
            rate_effect > 1e-5,
            "material times must change the solved bore pressure"
        );
        errors.push(worst);
    }
    assert!(errors[1] < 0.6 * errors[0], "spatial errors {errors:?}");
    assert!(errors[1] < 2e-3, "resolved bore impedance {errors:?}");
}
