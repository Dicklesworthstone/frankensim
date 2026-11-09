use super::*;
use crate::multicomponent::wall_exchange::{
    RecoilLaw, WallControl, WallThermal, exchange_cell_with_recoil,
};
use crate::multicomponent::{
    Conserved,
    tests::{gas, source},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

fn model() -> MixtureGas {
    gas(vec![source("A", R, 3.5, 0.0, 0.0)])
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:e}, expected={expected:e}, tolerance={tolerance:e}"
    );
}

fn control() -> RecoilControl {
    RecoilControl {
        wall: WallControl::default(),
        momentum_tolerance_ns: 1e-12,
    }
}

fn source_control() -> SourceCoupledDragControl {
    SourceCoupledDragControl {
        absolute_slip_tolerance_m_s: 1e-11,
        relative_slip_tolerance: 1e-8,
        absolute_heating_tolerance_j: 1e-12,
        relative_heating_tolerance: 1e-8,
        ..SourceCoupledDragControl::default()
    }
}

fn thermal() -> WallThermal {
    WallThermal::Lumped {
        temperature_k: 3.0,
        heat_capacity_j_per_k: 0.96,
        conductance_w_per_k: 0.0,
    }
}

#[allow(clippy::too_many_arguments)]
fn run(
    g: &MixtureGas,
    state: &Conserved,
    support: SlidingSupport,
    fraction: f64,
    dt: f64,
    contact: WallContact,
    wall: Option<&MaterialWallCaloric<'_>>,
    cancelled: impl FnMut() -> bool,
) -> Result<SourceCoupledRecoilUpdate, WallExchangeError> {
    exchange_cell_with_source_recoil(
        g,
        state.as_cell(),
        0.2,
        Time::new(dt),
        support,
        fraction,
        contact,
        wall,
        control(),
        source_control(),
        |_, candidate, primitive, _| {
            Ok(WallDrag::MovingLinear {
                rate_per_s: 4.0 * primitive.temperature / 2.0,
                wall_velocity_m_s: candidate.velocity_m_s,
            })
        },
        cancelled,
    )
}

// Exact continuum solution when mu(T)=mu0*T/T0 and Cv=2.5. The inverse
// squared relative slip obeys a linear ODE; reduced mass includes both bodies.
fn exact(mass: f64, fraction: f64, dt: f64) -> (f64, f64) {
    let gas_weight = mass / (0.2 + mass);
    let delta = fraction * gas_weight * 0.4 * 0.4 / (2.0 * 2.5 * 2.0);
    let decay = (-(4.0 / gas_weight) * (1.0 + delta) * dt).exp();
    let slip = 0.4 * decay * ((1.0 + delta) / (1.0 + delta * decay * decay)).sqrt();
    let dissipation = 0.5 * 0.2 * gas_weight * (0.4 * 0.4 - slip * slip);
    (slip, dissipation)
}

#[test]
fn actual_temperature_and_finite_inertia_match_independent_solution_and_galilean_work() {
    let g = model();
    for mass in [0.01, 0.3, 1e8] {
        // Keep light supports resolved instead of testing only their stiff zero.
        let dt = if mass < 0.2 {
            0.3 / (1.0 + 0.2 / mass)
        } else {
            0.3
        };
        for fraction in [0.0, 0.5, 1.0] {
            for boost in [0.0, 0.7] {
                if mass > 1.0 && boost != 0.0 {
                    continue;
                }
                for sign in [-1.0, 1.0] {
                    let state = g.conserved(&[1.0], boost + sign * 0.4, 2.0).unwrap();
                    let support = SlidingSupport {
                        mass_kg: mass,
                        velocity_m_s: boost,
                    };
                    let result = run(
                        &g,
                        &state,
                        support,
                        fraction,
                        dt,
                        WallContact {
                            drag: WallDrag::None,
                            thermal: thermal(),
                        },
                        None,
                        || false,
                    )
                    .unwrap();
                    let (slip, dissipation) = exact(mass, fraction, dt);
                    let reduced_mass = 0.2 * mass / (0.2 + mass);
                    let impulse = sign * reduced_mass * (slip - 0.4);
                    let recoil = result.recoil;
                    close(
                        recoil.gas.momentum - recoil.support.velocity_m_s,
                        sign * slip,
                        3e-9,
                    );
                    close(recoil.gas.momentum, state.momentum + impulse / 0.2, 3e-9);
                    close(recoil.support.velocity_m_s, boost - impulse / mass, 3e-9);
                    close(recoil.slip_dissipation_j, dissipation, 2e-10);
                    close(recoil.gas_slip_heating_j, fraction * dissipation, 2e-10);
                    close(
                        recoil.wall_slip_heating_j,
                        (1.0 - fraction) * dissipation,
                        2e-10,
                    );
                    let endpoint = Cell {
                        momentum: recoil.gas.momentum,
                        energy: recoil.gas.energy,
                        partial_densities: &state.partial_densities,
                    };
                    close(
                        g.primitive(endpoint, || false).unwrap().temperature,
                        2.0 + fraction * dissipation / (0.2 * 2.5),
                        2e-9,
                    );
                    let actual_impulse = 0.2 * (recoil.gas.momentum - state.momentum);
                    let work = actual_impulse * (boost - actual_impulse / (2.0 * mass));
                    close(recoil.gas.report.mechanical_work_to_gas_j, work, 2e-16);
                    close(
                        0.2 * (recoil.gas.energy - state.energy),
                        work - recoil.wall_slip_heating_j,
                        3e-16,
                    );
                    assert!(recoil.momentum_residual_ns.abs() <= control().momentum_tolerance_ns);
                    assert!(recoil.energy_residual_j.abs() < 3e-15);
                    assert!(
                        result.source.estimated_slip_error_m_s
                            <= source_control().absolute_slip_tolerance_m_s
                                + source_control().relative_slip_tolerance * 0.4
                    );
                    assert_eq!(endpoint.partial_densities, [1.0]);
                }
            }
        }
    }
}

#[test]
fn constant_resistance_recovers_original_recoil_and_thermal_completion() {
    let g = model();
    let state = g.conserved(&[1.0], 0.4, 2.0).unwrap();
    let support = SlidingSupport {
        mass_kg: 0.3,
        velocity_m_s: -0.2,
    };
    for (linear, quadratic) in [(0.0, 2.0), (4.0, 0.0), (4.0, 2.0), (1e6, 2.0)] {
        let thermal = WallThermal::Lumped {
            temperature_k: 3.0,
            heat_capacity_j_per_k: 0.96,
            conductance_w_per_k: 0.2,
        };
        let expected = exchange_cell_with_recoil(
            &g,
            state.as_cell(),
            0.2,
            Time::new(0.3),
            support,
            RecoilLaw {
                linear_rate_per_s: linear,
                quadratic_per_m: quadratic,
                gas_heat_fraction: 0.5,
            },
            thermal,
            control(),
            || false,
        )
        .unwrap();
        let actual = exchange_cell_with_source_recoil(
            &g,
            state.as_cell(),
            0.2,
            Time::new(0.3),
            support,
            0.5,
            WallContact {
                drag: WallDrag::None,
                thermal,
            },
            None,
            control(),
            source_control(),
            |_, candidate, _, _| {
                Ok(WallDrag::LinearQuadratic {
                    linear_rate_per_s: linear,
                    quadratic_per_m: quadratic,
                    wall_velocity_m_s: candidate.velocity_m_s,
                })
            },
            || false,
        )
        .unwrap();
        close(actual.recoil.gas.momentum, expected.gas.momentum, 3e-15);
        close(
            actual.recoil.support.velocity_m_s,
            expected.support.velocity_m_s,
            3e-15,
        );
        close(actual.recoil.gas.energy, expected.gas.energy, 3e-12);
        close(
            actual.recoil.slip_dissipation_j,
            expected.slip_dissipation_j,
            3e-16,
        );
        close(
            actual.recoil.gas.report.wall_energy_change_j,
            expected.gas.report.wall_energy_change_j,
            3e-12,
        );
        assert_eq!(actual.source.accepted_steps, 1);
    }
}

// Original finite storage c(T)=.7+.3T. It is deliberately nonlinear, so a
// frozen tangent cannot reproduce the exact finite frictional heat deposition.
struct Capacity;
impl ThermalLaw for Capacity {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !energy.is_finite() || energy <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(2.0 * energy / (0.7 + (0.49 + 0.6 * energy).sqrt()))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 0.7
    }
}

#[test]
fn wall_share_uses_original_material_integral_and_late_cancellation_is_atomic() {
    let g = model();
    let state = g.conserved(&[1.0], 0.4, 2.0).unwrap();
    let wall = MaterialWallCaloric::from_thermal_law(
        Capacity,
        DepositedCell {
            id: CellId(1),
            segment_id: 1,
            body_id: 1,
            material_id: MaterialId(1),
            mass_kg: 0.1,
            specific_h_j_kg: 0.7 * 3.0 + 0.15 * 3.0 * 3.0,
            geometry: CellGeometry {
                start_m: [0.0; 3],
                end_m: [1.0, 0.0, 0.0],
                volume_m3: 0.1,
                surface_area_m2: 1.0,
            },
            boundary: Boundary {
                ambient_k: 2.0,
                htc_w_m2_k: 0.0,
                reservoir: None,
                source_w: 0.0,
            },
        },
        [0.0715, 22.0],
    )
    .unwrap();
    let contact = wall.initial_contact(0.0, WallDrag::None).unwrap();
    let support = SlidingSupport {
        mass_kg: 0.3,
        velocity_m_s: 0.0,
    };
    let mut polls = 0;
    let result = run(&g, &state, support, 0.5, 0.3, contact, Some(&wall), || {
        polls += 1;
        false
    })
    .unwrap();
    let (_, dissipation) = exact(0.3, 0.5, 0.3);
    let heat = 0.5 * dissipation;
    let expected_rise = 2.0 * heat / (0.16 + (0.16_f64.powi(2) + 2.0 * 0.1 * 0.3 * heat).sqrt());
    let actual = wall
        .temperature_k(result.recoil.gas.contact.thermal)
        .unwrap();
    close(actual, 3.0 + expected_rise, 3e-10);
    assert!((actual - (3.0 + heat / 0.16)).abs() > 2e-6);
    close(result.recoil.gas.report.wall_energy_change_j, heat, 2e-10);
    let mut visits = 0;
    let cancelled = run(&g, &state, support, 0.5, 0.3, contact, Some(&wall), || {
        visits += 1;
        visits == polls
    });
    assert_eq!(cancelled, Err(WallExchangeError::Cancelled));
    assert_eq!(wall.initial_contact(0.0, WallDrag::None).unwrap(), contact);
    assert_eq!(
        run(&g, &state, support, 0.5, 0.3, contact, Some(&wall), || {
            false
        })
        .unwrap(),
        result
    );
}

#[test]
fn original_source_must_admit_actual_returned_state_after_thermal_completion() {
    let g = model();
    let state = g.conserved(&[1.0], 0.4, 2.0).unwrap();
    let support = SlidingSupport {
        mass_kg: 0.3,
        velocity_m_s: 0.0,
    };
    let mut last_temperature = 0.0;
    let result = exchange_cell_with_source_recoil(
        &g,
        state.as_cell(),
        0.2,
        Time::new(0.3),
        support,
        1.0,
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Reservoir {
                temperature_k: 10.0,
                conductance_w_per_k: 1.0,
            },
        },
        None,
        control(),
        source_control(),
        |_, candidate, primitive, _| {
            last_temperature = primitive.temperature;
            if primitive.temperature > 3.0 {
                return Err(WallExchangeError::TransportDomain {
                    quantity: "manufactured final temperature",
                    value: primitive.temperature,
                    minimum: 0.001,
                    maximum: 3.0,
                });
            }
            Ok(WallDrag::MovingLinear {
                rate_per_s: 4.0,
                wall_velocity_m_s: candidate.velocity_m_s,
            })
        },
        || false,
    );
    assert!(matches!(
        result,
        Err(WallExchangeError::TransportDomain { .. })
    ));
    assert!(last_temperature > 3.0);
    assert_eq!(state.momentum, 0.4);
    assert_eq!(
        support,
        SlidingSupport {
            mass_kg: 0.3,
            velocity_m_s: 0.0
        }
    );
}
