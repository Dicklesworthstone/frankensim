use super::*;
use crate::{LossModel, Termination, input_impedance};
use fs_la::factor::lu;
use fs_material::gas::{GasSpec, GasState};
use fs_math::c64::C64;
use fs_phs::linear_exact::DenseExponentialBudget;

fn air(pressure: f64) -> GasState {
    GasState::try_new(&GasSpec::dry_air_ussa1976(), 293.15, pressure).unwrap()
}

fn heat_budget() -> DenseExponentialBudget {
    DenseExponentialBudget {
        max_dimension: 256,
        max_dense_bytes: 64 << 20,
        max_multiply_terms: 40_000_000_000,
    }
}

fn heat_close(actual: f64, expected: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= relative * expected.abs().max(1e-30),
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn viscous_heat_matches_dimensional_poiseuille_decay_and_split_intervals() {
    let dg = duct(1, DgBoundary::IdealOpen, DgBoundary::IdealOpen)
        .with_viscothermal(8, [50.0, 2000.0], 0.02)
        .unwrap();
    let flow = 1e-6;
    let mut state = dg.project_initial_state(|_| 0.0, |_| flow);
    let m = dg.field_dofs;
    let model = dg.viscothermal_model().unwrap();
    // Initialize the steady parabolic-flow memory, not an empty boundary layer.
    // The original dimensional DC loss must then be 8 mu L Q²/(pi a⁴).
    for (arm, term) in model.series_terms().iter().enumerate() {
        for i in 0..m {
            state.tilde[2 * m + arm * m + i] =
                (term.gain_per_s / term.rate_per_s).sqrt() * state.tilde[m + i];
        }
    }
    let power = 8.0 * dg.gas().dynamic_viscosity * 0.05 * flow * flow
        / (core::f64::consts::PI * 0.001_f64.powi(4));
    heat_close(dg.wall_dissipated_power_w(&state).unwrap(), power, 3e-12);
    let initial = state.clone();
    let dt = 1e-6;
    let step = dg.propagator(dt, 1e-13).unwrap();
    let heat = dg
        .prepare_viscous_heat(&step, heat_budget(), &mut || false)
        .unwrap();
    let actual = heat.integrated_heat_j(&state, &[], &mut || false).unwrap();
    assert_eq!(state, initial);
    assert_eq!(heat.operator_identity(), dg.operator_identity());
    assert_eq!(heat.duration(), dt);
    heat_close(actual / dt, power, 0.001);
    step.step(&mut state, &[]).unwrap();
    heat_close(
        actual,
        initial.stored_energy_j() - state.stored_energy_j(),
        2e-8,
    );
    let half = dg.propagator(dt / 2.0, 1e-13).unwrap();
    let half_heat = dg
        .prepare_viscous_heat(&half, heat_budget(), &mut || false)
        .unwrap();
    let first = half_heat
        .integrated_heat_j(&initial, &[], &mut || false)
        .unwrap();
    let mut middle = initial;
    half.step(&mut middle, &[]).unwrap();
    let second = half_heat
        .integrated_heat_j(&middle, &[], &mut || false)
        .unwrap();
    heat_close(first + second, actual, 2e-9);
}

#[test]
fn viscous_heat_retains_regional_fluid_sources_and_spatial_quadrature() {
    use fs_material::fluid_acoustics::FluidAcousticState;
    use fs_material::fluid_acoustics::thermoviscous::ThermoviscousFluidState;
    use fs_material::fluid_transport::FluidTransportState;
    let fluid = |viscosity| {
        ThermoviscousFluidState::try_new(
            FluidAcousticState::try_new(2.0, 300.0).unwrap(),
            FluidTransportState::try_new(2.0, 1000.0, 0.02, viscosity).unwrap(),
            1000.0,
        )
        .unwrap()
    };
    // Cp=Cv leaves only physical viscosity: no thermal loss can be assigned to
    // it. Equal rho/c isolates the independent, fourfold regional viscosity.
    let dg = DgDuct::from_material_regions(
        &cylinder(),
        &[
            DgFluidRegion {
                end_m: 0.025,
                material: fluid(2e-4),
            },
            DgFluidRegion {
                end_m: 0.05,
                material: fluid(8e-4),
            },
        ],
        DgMesh {
            elements: 2,
            order: 1,
        },
        DgBoundary::Port(PortId(7)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_spatial_viscothermal(DgWallResolution {
        poles_per_field: 8,
        quadrature_points: 2,
        band_rad_s: [50.0, 2000.0].map(|f| core::f64::consts::TAU * f),
        maximum_relative_error: 0.02,
    })
    .unwrap();
    let mut state = dg.project_initial_state(|x| 4.0 + 20.0 * x, |x| 1e-6 * (1.0 + x / 0.05));
    let initial = state.clone();
    let ports = [PortInput {
        id: PortId(7),
        effort: 2e-6,
    }];
    let dt = 2e-6;
    let step = dg.propagator(dt, 1e-13).unwrap();
    let heat = dg
        .prepare_viscous_heat(&step, heat_budget(), &mut || false)
        .unwrap();
    let actual = heat
        .integrated_heat_j(&initial, &ports, &mut || false)
        .unwrap();
    // Independent time quadrature of the pre-existing physical-power owner,
    // evaluated on separately propagated original states under the same drive.
    let (points, weights) = gauss_legendre_rule(12);
    let mut expected = 0.0;
    for (point, weight) in points.iter().zip(&weights) {
        let mut sample = initial.clone();
        dg.propagator(0.5 * dt * (point + 1.0), 1e-13)
            .unwrap()
            .step(&mut sample, &ports)
            .unwrap();
        expected += 0.5 * dt * weight * dg.wall_dissipated_power_w(&sample).unwrap();
    }
    assert!(actual > 0.0);
    heat_close(actual, expected, 2e-9);
    let output = step.step(&mut state, &ports).unwrap();
    heat_close(
        actual,
        initial.stored_energy_j() + output.ledger.port_work_step_j - state.stored_energy_j(),
        2e-8,
    );
    assert_eq!(
        dg.thermoviscous_regions().unwrap()[0]
            .transport()
            .dynamic_viscosity_pa_s(),
        2e-4
    );
    assert_eq!(
        dg.thermoviscous_regions().unwrap()[1]
            .transport()
            .dynamic_viscosity_pa_s(),
        8e-4
    );
}

#[test]
fn viscous_heat_excludes_thermal_and_solid_losses() {
    let dg = duct(1, DgBoundary::Closed, DgBoundary::Closed)
        .with_viscothermal(8, [50.0, 2000.0], 0.02)
        .unwrap()
        .with_compliant_wall(fs_phs::WallPin {
            surface_density: 0.2,
            stiffness_per_area: 1e7,
            resistance: 300.0,
        })
        .unwrap();
    // Uniform compression and radial wall motion excite no axial velocity.
    // Thermal and solid dissipation are live, but gas viscous heat stays zero.
    let mut state = dg.project_initial_state(|_| 30.0, |_| 0.0);
    let step = dg.propagator(1e-5, 1e-13).unwrap();
    for _ in 0..4 {
        step.step(&mut state, &[]).unwrap();
    }
    let initial = state.clone();
    let heat = dg
        .prepare_viscous_heat(&step, heat_budget(), &mut || false)
        .unwrap()
        .integrated_heat_j(&state, &[], &mut || false)
        .unwrap();
    assert!(dg.wall_dissipated_power_w(&state).unwrap() > 0.0);
    step.step(&mut state, &[]).unwrap();
    let total = initial.stored_energy_j() - state.stored_energy_j();
    assert!(total > 0.0);
    assert!(
        heat < 1e-18 * total,
        "viscous={heat:e}, all losses={total:e}"
    );
}

#[test]
fn viscous_heat_zero_source_and_native_admission_preserve_state() {
    let dg = duct(1, DgBoundary::Port(PortId(7)), DgBoundary::IdealOpen);
    let state = dg.project_initial_state(|x| 1.0 + x, |_| 1e-6);
    let initial = state.clone();
    let step = dg.propagator(1e-6, 1e-13).unwrap();
    let heat = dg
        .prepare_viscous_heat(&step, heat_budget(), &mut || false)
        .unwrap();
    assert_eq!(
        heat.integrated_heat_j(&state, &[], &mut || false).unwrap(),
        0.0
    );
    assert!(
        heat.integrated_heat_j(
            &state,
            &[PortInput {
                id: PortId(8),
                effort: 1.0
            }],
            &mut || false
        )
        .is_err()
    );
    let other = duct(2, DgBoundary::Closed, DgBoundary::Closed);
    assert!(
        dg.prepare_viscous_heat(
            &other.propagator(1e-6, 1e-13).unwrap(),
            heat_budget(),
            &mut || false
        )
        .is_err()
    );
    assert!(
        heat.integrated_heat_j(
            &other.project_initial_state(|_| 0.0, |_| 0.0),
            &[],
            &mut || false
        )
        .is_err()
    );
    assert!(
        dg.prepare_viscous_heat(
            &step,
            DenseExponentialBudget {
                max_dimension: 1,
                ..heat_budget()
            },
            &mut || false
        )
        .is_err()
    );
    assert!(matches!(
        dg.prepare_viscous_heat(&step, heat_budget(), &mut || true),
        Err(DuctError::Cancelled { .. })
    ));
    assert_eq!(state, initial);
}

fn cylinder() -> Duct {
    Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: 0.05,
        }],
    }
}
fn duct(elements: usize, inlet: DgBoundary, outlet: DgBoundary) -> DgDuct {
    DgDuct::new(
        &cylinder(),
        &air(101325.0),
        DgMesh { elements, order: 2 },
        inlet,
        outlet,
    )
    .unwrap()
}

// Harmonic response of the ACTUAL assembled operator. Condense only its
// independent memory coordinates, then use fs-la's existing real LU on the
// doubled complex system. This tests the physical coupling and mass maps,
// independently of the pole producer's analytic response helper.
fn harmonic_port(duct: &DgDuct, omega: f64) -> C64 {
    let n = duct.total_dofs();
    let base = 2 * duct.field_dofs;
    let size = 2 * base;
    let mut matrix = vec![0.0; size * size];
    for i in 0..base {
        for j in 0..base {
            let mut value = C64::new(-duct.l_tilde[i * n + j], if i == j { omega } else { 0.0 });
            for k in base..n {
                value = value
                    - C64::from_re(duct.l_tilde[i * n + k] * duct.l_tilde[k * n + j])
                        / C64::new(-duct.l_tilde[k * n + k], omega);
            }
            matrix[i * size + j] = value.re;
            matrix[i * size + base + j] = -value.im;
            matrix[(base + i) * size + j] = value.im;
            matrix[(base + i) * size + base + j] = value.re;
        }
    }
    let port = &duct.ports[0].1;
    let mut solution = vec![0.0; size];
    solution[..base].copy_from_slice(&port[..base]);
    lu(&matrix, size).unwrap().solve(&mut solution);
    (0..base).fold(C64::ZERO, |sum, i| {
        sum + C64::new(solution[i], solution[base + i]).scale(port[i])
    })
}

#[test]
fn narrow_dg_viscothermal_impedance_converges_to_exact_bessel() {
    let gas = air(101325.0);
    let mut previous = f64::INFINITY;
    for elements in [1, 2, 4] {
        let dg = duct(elements, DgBoundary::Port(PortId(0)), DgBoundary::IdealOpen)
            .with_viscothermal(16, [50.0, 8000.0], 1e-4)
            .unwrap();
        let mut worst = 0.0_f64;
        for frequency in [100.0, 500.0, 1000.0] {
            let omega = core::f64::consts::TAU * frequency;
            let exact = input_impedance(
                &cylinder(),
                &gas,
                omega,
                LossModel::Bessel,
                Termination::IdealOpen,
            )
            .unwrap()
            .impedance
            .conj();
            let actual = harmonic_port(&dg, omega);
            assert!(actual.re > 0.0);
            worst = worst.max((actual - exact).abs() / exact.abs());
        }
        println!("DG elements={elements} exact-Bessel complex error={worst:e}");
        assert!(worst < previous);
        previous = worst;
    }
    assert!(previous < 2e-3);
}

#[test]
fn dg_viscothermal_memory_uses_one_passive_exponential_and_bound_state() {
    let original = duct(1, DgBoundary::Closed, DgBoundary::Closed);
    let dg = original
        .clone()
        .with_viscothermal(8, [50.0, 2000.0], 0.01)
        .unwrap();
    assert_ne!(dg.operator_identity(), original.operator_identity());
    assert_eq!(dg.total_dofs(), 2 * dg.field_dofs() * 9);
    let initial = dg.project_initial_state(|x| det::sin(core::f64::consts::PI * x / 0.05), |_| 0.0);
    assert_eq!(initial.tilde.len(), dg.total_dofs());
    assert!(
        initial.tilde[2 * dg.field_dofs()..]
            .iter()
            .all(|v| *v == 0.0)
    );
    let coarse = dg.propagator(2e-5, 1e-11).unwrap();
    let fine = dg.propagator(1e-5, 1e-11).unwrap();
    let mut state = initial.clone();
    let mut refined = initial.clone();
    let mut previous = state.stored_energy_j();
    for _ in 0..32 {
        let row = coarse.step(&mut state, &[]).unwrap();
        fine.step(&mut refined, &[]).unwrap();
        fine.step(&mut refined, &[]).unwrap();
        assert!(!row.ledger.balance_checked); // independent loss integral is not claimed
        assert!(row.ledger.cumulative_dissipated_j >= 0.0);
        assert!(state.stored_energy_j() <= previous * (1.0 + 1e-12));
        previous = state.stored_energy_j();
    }
    assert!(previous < 0.99 * initial.stored_energy_j());
    let difference = state
        .tilde
        .iter()
        .zip(&refined.tilde)
        .map(|(a, b)| (a - b) * (a - b))
        .sum::<f64>();
    assert!(difference < 1e-16 * initial.stored_energy_j());
    coarse.sync_nodal(&mut state).unwrap();
    assert!(state.pressure().iter().all(|p| p.is_finite()));
    let mut wrong = original.zero_state();
    let before = wrong.clone();
    assert!(coarse.step(&mut wrong, &[]).is_err());
    assert_eq!(wrong, before);
}

#[test]
fn material_pressure_changes_bessel_drag_and_the_complete_duct_response() {
    let a = duct(2, DgBoundary::Port(PortId(0)), DgBoundary::IdealOpen)
        .with_viscothermal(16, [50.0, 8000.0], 1e-4)
        .unwrap();
    let b = DgDuct::new(
        &cylinder(),
        &air(60000.0),
        DgMesh {
            elements: 2,
            order: 2,
        },
        DgBoundary::Port(PortId(0)),
        DgBoundary::IdealOpen,
    )
    .unwrap()
    .with_viscothermal(16, [50.0, 8000.0], 1e-4)
    .unwrap();
    assert!(
        b.viscothermal_model().unwrap().series_drag_per_s()
            > a.viscothermal_model().unwrap().series_drag_per_s()
    );
    let za = harmonic_port(&a, core::f64::consts::TAU * 500.0);
    let zb = harmonic_port(&b, core::f64::consts::TAU * 500.0);
    assert!((za - zb).abs() / za.abs() > 0.1);
}

#[test]
fn dg_wall_memory_refuses_unsupported_sources_and_oversized_states() {
    let valid = duct(1, DgBoundary::Closed, DgBoundary::IdealOpen);
    assert!(
        valid
            .clone()
            .with_viscothermal(2, [50.0, 8000.0], 1e-4)
            .is_err()
    );
    let attached = valid.with_viscothermal(16, [50.0, 8000.0], 1e-4).unwrap();
    assert!(
        attached
            .with_viscothermal(16, [50.0, 8000.0], 1e-4)
            .is_err()
    );
    let cone = Duct {
        segments: vec![Segment::Cone {
            inlet_radius: 0.001,
            outlet_radius: 0.002,
            length: 0.05,
        }],
    };
    let cone = DgDuct::new(
        &cone,
        &air(101325.0),
        DgMesh {
            elements: 1,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    )
    .unwrap();
    assert!(cone.with_viscothermal(16, [50.0, 8000.0], 1e-4).is_err());
    let large = duct(16, DgBoundary::Closed, DgBoundary::IdealOpen);
    assert!(large.with_viscothermal(32, [50.0, 8000.0], 1e-4).is_err());
    let field = crate::profile::gas_field::tests::field([290.0, 310.0], [0.8; 2], 101325.0);
    let geometry = Duct {
        segments: vec![Segment::Cylinder {
            radius: 0.001,
            length: field.length_m(),
        }],
    };
    let nonuniform = DgDuct::new_with_gas_field(
        &geometry,
        &field,
        DgMesh {
            elements: 4,
            order: 2,
        },
        DgBoundary::Closed,
        DgBoundary::IdealOpen,
    )
    .unwrap();
    assert!(
        nonuniform
            .with_viscothermal(16, [50.0, 8000.0], 1e-4)
            .is_err()
    );
}
