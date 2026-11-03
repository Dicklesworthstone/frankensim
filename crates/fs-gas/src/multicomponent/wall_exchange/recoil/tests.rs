//! Manufactured caloric source through the exact production heat/recoil owners.
use super::super::WallGasState;
use super::*;
use crate::dynamics::CaloricDomain;
use fs_qty::Temperature;

struct CaloricFixture {
    density: f64,
    cv: f64,
}
impl WallGasCaloric for CaloricFixture {
    fn domain(&self) -> CaloricDomain {
        CaloricDomain::new(Temperature::new(0.1), Temperature::new(50.0))
    }
    fn primitive(
        &self,
        momentum: f64,
        energy: f64,
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<WallGasState, WallExchangeError> {
        if cancelled() {
            return Err(WallExchangeError::Cancelled);
        }
        let velocity = momentum / self.density;
        let temperature = (energy - 0.5 * momentum * velocity) / (self.density * self.cv);
        if !(0.1..=50.0).contains(&temperature) {
            return Err(WallExchangeError::Invalid(
                "manufactured caloric source domain",
            ));
        }
        Ok(WallGasState {
            density: self.density,
            velocity,
            temperature,
        })
    }
    fn caloric(&self, temperature: f64) -> Result<(f64, f64), WallExchangeError> {
        Ok((self.density * self.cv * temperature, self.density * self.cv))
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        (a - b).abs() <= tolerance,
        "{a:e} != {b:e}, tolerance {tolerance:e}"
    );
}
fn controls() -> RecoilControl {
    RecoilControl {
        wall: WallControl {
            absolute_energy_tolerance_j: 1e-12,
            relative_energy_tolerance: 1e-12,
            maximum_iterations: 80,
        },
        momentum_tolerance_ns: 1e-12,
    }
}
fn law(fraction: f64) -> RecoilLaw {
    RecoilLaw {
        linear_rate_per_s: 0.8,
        quadratic_per_m: 0.7,
        gas_heat_fraction: fraction,
    }
}
fn thermal(conductance: f64) -> WallThermal {
    WallThermal::Lumped {
        temperature_k: 3.0,
        heat_capacity_j_per_k: 0.8,
        conductance_w_per_k: conductance,
    }
}
fn wall_temperature(thermal: WallThermal) -> f64 {
    let WallThermal::Lumped { temperature_k, .. } = thermal else {
        panic!("expected finite wall")
    };
    temperature_k
}

#[test]
fn recoil_and_heat_partition_close_actual_gas_support_and_wall_storage() {
    let gas = CaloricFixture {
        density: 1.0,
        cv: 2.5,
    };
    let (volume, dt, velocity) = (0.4, 0.2, 0.9);
    let support = SlidingSupport {
        mass_kg: 0.6,
        velocity_m_s: -0.2,
    };
    let energy = 2.0 * gas.cv + 0.5 * velocity * velocity;
    for fraction in [0.0, 0.5, 1.0] {
        for conductance in [0.0, 0.3] {
            let result = exchange_caloric_with_recoil(
                &gas,
                velocity,
                energy,
                volume,
                Time::new(dt),
                support,
                law(fraction),
                thermal(conductance),
                controls(),
                || false,
            )
            .unwrap();
            let mass = gas.density * volume;
            let reduced = mass * support.mass_kg / (mass + support.mass_kg);
            let elapsed = dt * (1.0 + mass / support.mass_kg);
            let slip = velocity - support.velocity_m_s;
            let final_slip =
                slip / ((0.8 * elapsed).exp() + 0.7 * slip.abs() / 0.8 * (0.8 * elapsed).exp_m1());
            let loss = 0.5 * reduced * (slip * slip - final_slip * final_slip);
            let gas_capacity = mass * gas.cv;
            let gas_temperature = 2.0 + fraction * loss / gas_capacity;
            let support_temperature = 3.0 + (1.0 - fraction) * loss / 0.8;
            let heat = dt * conductance * (support_temperature - gas_temperature)
                / (1.0 + dt * conductance * (1.0 / gas_capacity + 1.0 / 0.8));
            let returned = gas
                .primitive(result.gas.momentum, result.gas.energy, &mut || false)
                .unwrap();
            close(
                returned.temperature,
                gas_temperature + heat / gas_capacity,
                3e-12,
            );
            let final_wall_temperature = wall_temperature(result.gas.contact.thermal);
            close(
                final_wall_temperature,
                support_temperature - heat / 0.8,
                3e-12,
            );
            close(result.slip_dissipation_j, loss, 3e-14);
            close(result.gas.report.drag_heating_j, loss, 3e-14);
            close(result.gas_slip_heating_j, fraction * loss, 3e-14);
            close(result.wall_slip_heating_j, (1.0 - fraction) * loss, 3e-14);
            close(
                result.gas.report.heat_to_gas_j,
                heat - (1.0 - fraction) * loss,
                3e-12,
            );
            close(
                gas_capacity * (returned.temperature - 2.0),
                result.gas.report.drag_heating_j + result.gas.report.heat_to_gas_j,
                3e-12,
            );
            let wall_kinetic = 0.5
                * support.mass_kg
                * (result.support.velocity_m_s - support.velocity_m_s)
                * (result.support.velocity_m_s + support.velocity_m_s);
            let actual_balance = volume * (result.gas.energy - energy)
                + wall_kinetic
                + 0.8 * (final_wall_temperature - 3.0);
            close(actual_balance, 0.0, 3e-12);
            close(result.energy_residual_j, actual_balance, 3e-14);
            close(
                volume * (result.gas.momentum - velocity)
                    + support.mass_kg * (result.support.velocity_m_s - support.velocity_m_s),
                0.0,
                3e-14,
            );
        }
    }
}

#[test]
fn recoil_complete_source_is_frame_covariant_and_zero_duration_is_exact_identity() {
    let gas = CaloricFixture {
        density: 1.0,
        cv: 2.5,
    };
    let (velocity, volume, dt) = (0.9, 0.4, 0.2);
    let support = SlidingSupport {
        mass_kg: 0.6,
        velocity_m_s: -0.2,
    };
    let energy = 2.0 * gas.cv + 0.5 * velocity * velocity;
    let base = exchange_caloric_with_recoil(
        &gas,
        velocity,
        energy,
        volume,
        Time::new(dt),
        support,
        law(0.5),
        thermal(0.3),
        controls(),
        || false,
    )
    .unwrap();
    for boost in [-1024.0, 1024.0] {
        // The laboratory total-energy ulp at this boost is ~1e-10 J/m^3;
        // the authored absolute budget must admit its represented storage.
        let mut control = controls();
        control.wall.absolute_energy_tolerance_j = 1e-9;
        let boosted = exchange_caloric_with_recoil(
            &gas,
            velocity + boost,
            2.0 * gas.cv + 0.5 * (velocity + boost).powi(2),
            volume,
            Time::new(dt),
            SlidingSupport {
                velocity_m_s: support.velocity_m_s + boost,
                ..support
            },
            law(0.5),
            thermal(0.3),
            control,
            || false,
        )
        .unwrap();
        close(boosted.gas.momentum, base.gas.momentum + boost, 3e-13);
        close(
            boosted.support.velocity_m_s,
            base.support.velocity_m_s + boost,
            3e-13,
        );
        close(boosted.slip_dissipation_j, base.slip_dissipation_j, 3e-13);
        close(
            boosted.gas.report.drag_heating_j,
            base.gas.report.drag_heating_j,
            3e-13,
        );
        close(
            wall_temperature(boosted.gas.contact.thermal),
            wall_temperature(base.gas.contact.thermal),
            1e-9,
        );
        close(boosted.energy_residual_j, 0.0, 1e-9);
    }
    let identity = exchange_caloric_with_recoil(
        &gas,
        velocity,
        energy,
        volume,
        Time::new(0.0),
        support,
        law(0.5),
        thermal(0.3),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(identity.support, support);
    assert_eq!(identity.gas.momentum, velocity);
    assert_eq!(identity.gas.energy, energy);
    assert_eq!(identity.gas.contact.thermal, thermal(0.3));
    assert_eq!(identity.slip_dissipation_j, 0.0);
    assert_eq!(identity.energy_residual_j, 0.0);
}

#[test]
fn recoil_refusal_and_late_cancellation_publish_neither_endpoint() {
    let gas = CaloricFixture {
        density: 1.0,
        cv: 2.5,
    };
    let support = SlidingSupport {
        mass_kg: 0.6,
        velocity_m_s: -0.2,
    };
    let (velocity, energy, volume) = (0.9, 5.405, 0.4);
    let mut polls = 0;
    let accepted = exchange_caloric_with_recoil(
        &gas,
        velocity,
        energy,
        volume,
        Time::new(0.2),
        support,
        law(0.5),
        thermal(0.3),
        controls(),
        || {
            polls += 1;
            false
        },
    )
    .unwrap();
    let mut visits = 0;
    let cancelled = exchange_caloric_with_recoil(
        &gas,
        velocity,
        energy,
        volume,
        Time::new(0.2),
        support,
        law(0.5),
        thermal(0.3),
        controls(),
        || {
            visits += 1;
            visits >= polls
        },
    );
    assert_eq!(cancelled, Err(WallExchangeError::Cancelled));
    let replay = exchange_caloric_with_recoil(
        &gas,
        velocity,
        energy,
        volume,
        Time::new(0.2),
        support,
        law(0.5),
        thermal(0.3),
        controls(),
        || false,
    )
    .unwrap();
    assert_eq!(replay, accepted);
    assert!(
        exchange_caloric_with_recoil(
            &gas,
            velocity,
            energy,
            volume,
            Time::new(0.2),
            support,
            law(0.5),
            WallThermal::Adiabatic,
            controls(),
            || false
        )
        .is_err()
    );
    assert!(
        exchange_caloric_with_recoil(
            &gas,
            50.0,
            1255.0,
            volume,
            Time::new(0.2),
            SlidingSupport {
                velocity_m_s: -50.0,
                ..support
            },
            law(1.0),
            WallThermal::Adiabatic,
            controls(),
            || false
        )
        .is_err()
    );
    assert_eq!(
        support,
        SlidingSupport {
            mass_kg: 0.6,
            velocity_m_s: -0.2
        }
    );
}
