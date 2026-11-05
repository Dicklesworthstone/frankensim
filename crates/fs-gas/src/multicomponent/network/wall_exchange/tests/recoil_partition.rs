//! G1 nonlinear frictional wall heating and G4 complete network publication.
use super::*;
use crate::multicomponent::wall_exchange::{MaterialWallCaloric, SlidingSupport};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};

const WALL_MASS: f64 = 0.1;
const A: f64 = 0.4;

// Manufactured u(T)=a*T²/2 J/kg. This wall's thermal mass is already part of
// the support's 0.1 kg inertia; it is not a second mechanical inventory.
struct QuadraticStorage(ThermalStorageBasis);
impl ThermalLaw for QuadraticStorage {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.0
    }
    fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
        if !u.is_finite() || u <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((2.0 * u / A).sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        100.0
    }
}

fn cell(u: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(17),
        segment_id: 3,
        body_id: 8,
        material_id: MaterialId(29),
        mass_kg: WALL_MASS,
        specific_h_j_kg: u,
        geometry: CellGeometry {
            start_m: [0.0; 3],
            end_m: [1.0, 0.0, 0.0],
            volume_m3: 0.01,
            surface_area_m2: 1.0,
        },
        boundary: Boundary {
            ambient_k: 2.0,
            htc_w_m2_k: 0.0,
            reservoir: None,
            source_w: 0.0,
        },
    }
}
fn source(basis: ThermalStorageBasis, range: [f64; 2]) -> MaterialWallCaloric<'static> {
    MaterialWallCaloric::from_thermal_law(QuadraticStorage(basis), cell(0.5 * A * 25.0), range)
        .unwrap()
}
fn support() -> SlidingSupport {
    SlidingSupport {
        mass_kg: WALL_MASS,
        velocity_m_s: 0.5,
    }
}
fn bind(
    n: &MixtureNetwork1d,
    sources: [MaterialWallCaloric<'static>; 2],
    conductance: f64,
    fraction: Option<f64>,
    material_first: bool,
) -> Result<NetworkWallOperator<'static>, WallExchangeError> {
    let contacts = sources
        .iter()
        .map(|s| s.initial_contact(conductance, WallDrag::Linear { rate_per_s: 3.0 }))
        .collect::<Result<Vec<_>, _>>()?;
    let walls = vec![sources.into_iter().map(Some).collect()];
    let supports = vec![vec![Some(support()); 2]];
    let wall = NetworkWallOperator::new(n, vec![contacts], vec![], || false)?;
    let wall = if material_first {
        wall.with_material_heat_storage(walls, vec![], || false)?
            .with_finite_supports(supports, 1e-11, || false)?
    } else {
        wall.with_finite_supports(supports, 1e-11, || false)?
            .with_material_heat_storage(walls, vec![], || false)?
    };
    match fraction {
        Some(f) => wall.with_finite_support_heat_partition(vec![vec![f; 2]], || false),
        None => Ok(wall),
    }
}
fn stored(contact: WallContact) -> f64 {
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = contact.thermal
    else {
        panic!("expected original wall energy coordinate");
    };
    mass_kg * specific_energy_j_kg
}
fn slip_loss(dt: f64) -> f64 {
    // mu=0.2*0.1/0.3; slip=1.5*exp[-3*(1+0.2/0.1)*dt].
    0.5 * (0.2 * 0.1 / 0.3) * 1.5_f64.powi(2) * -(-18.0 * dt).exp_m1()
}

#[test]
fn g1_nonlinear_wall_partition_matches_independent_recoil_and_quadratic_heat() {
    let g = gas(1, false, 0.0);
    let dt = 0.2_f64;
    let loss = slip_loss(dt);
    let slip = 1.5 * (-9.0 * dt).exp();
    let gas_velocity = 1.5 + slip / 3.0;
    let wall_velocity = 1.5 - 2.0 * slip / 3.0;
    let kinetic = 0.5 * WALL_MASS * (wall_velocity * wall_velocity - 0.25);
    for fraction in [0.0, 0.35, 1.0] {
        for conductance in [0.0, 0.3] {
            for material_first in [false, true] {
                let mut n = network(&g, &[2.0], false);
                let before = fields(&n);
                let original = source(ThermalStorageBasis::InternalEnergy, [0.01, 50.0]);
                let mut wall = bind(
                    &n,
                    [original.clone(), original],
                    conductance,
                    Some(fraction),
                    material_first,
                )
                .unwrap();
                let report = wall
                    .advance_for(&mut n, Time::new(dt), controls(), || false)
                    .unwrap();
                let wall_heat = (1.0 - fraction) * loss;
                let gas_after_friction = 2.0 + fraction * loss / 0.5;
                // After friction: Uwall=M*a*T²/2. Eliminate the gas endpoint
                // from Q=dt*G*(Tw-Tg) to solve the independent positive quadratic.
                let quadratic = 0.5 * WALL_MASS * A;
                let transfer = dt * conductance / (1.0 + dt * conductance / 0.5);
                let constant = quadratic * 25.0 + wall_heat + transfer * gas_after_friction;
                let tw = 2.0 * constant
                    / (transfer + (transfer * transfer + 4.0 * quadratic * constant).sqrt());
                let q = transfer * (tw - gas_after_friction);
                let mut gas_change = 0.0;
                for i in 0..2 {
                    let state = n.branches()[0].cell(i).unwrap();
                    let primitive = g.primitive(state, || false).unwrap();
                    let contact = wall.branch_contacts(0).unwrap()[i];
                    let caloric = wall.branch_heat_storage(0).unwrap()[i].as_ref().unwrap();
                    let exchange = wall.accepted_branch_support_exchanges(0).unwrap()[i].unwrap();
                    close(primitive.velocity, gas_velocity, 3e-12);
                    close(primitive.temperature, gas_after_friction + q / 0.5, 3e-10);
                    close(caloric.temperature_k(contact.thermal).unwrap(), tw, 3e-10);
                    close(stored(contact), quadratic * tw * tw, 3e-11);
                    close(
                        wall.branch_supports(0).unwrap()[i].unwrap().velocity_m_s,
                        wall_velocity,
                        3e-12,
                    );
                    assert_eq!(state.partial_densities, before.0[0][i].partial_densities);
                    close(exchange.wall_slip_heating_j, wall_heat, 3e-12);
                    close(exchange.gas_energy_change_j, -kinetic - wall_heat, 3e-12);
                    close(exchange.kinetic_energy_change_j, kinetic, 3e-12);
                    close(exchange.drag_heating_j, loss, 3e-12);
                    close(exchange.momentum_residual_ns, 0.0, 3e-13);
                    close(exchange.energy_residual_j, 0.0, 3e-12);
                    gas_change += 0.2 * (state.energy - before.0[0][i].energy);
                }
                close(gas_change + 2.0 * (kinetic + wall_heat - q), 0.0, 4e-11);
                close(report.wall_energy_change_j, 2.0 * (wall_heat - q), 4e-11);
                close(report.heat_to_gas_j, 2.0 * (q - wall_heat), 4e-11);
                close(report.mechanical_work_to_gas_j, -2.0 * kinetic, 4e-12);
                close(report.drag_heating_j, 2.0 * loss, 4e-12);
                close(report.energy_residual_j, 0.0, 4e-11);
                assert_eq!(report.external_heat_to_gas_j, 0.0);
                assert_eq!(n.time(), Time::new(0.0));
                if fraction == 1.0 && material_first {
                    let mut reference = network(&g, &[2.0], false);
                    let original = source(ThermalStorageBasis::InternalEnergy, [0.01, 50.0]);
                    let mut default = bind(
                        &reference,
                        [original.clone(), original],
                        conductance,
                        None,
                        false,
                    )
                    .unwrap();
                    assert_eq!(
                        default
                            .advance_for(&mut reference, Time::new(dt), controls(), || false)
                            .unwrap(),
                        report
                    );
                    assert_eq!(fields(&reference), fields(&n));
                    assert_eq!(default.branch_contacts(0), wall.branch_contacts(0));
                    assert_eq!(
                        default.accepted_branch_support_exchanges(0),
                        wall.accepted_branch_support_exchanges(0)
                    );
                }
            }
        }
    }
}

#[test]
fn g1_native_plateau_keeps_friction_energy_before_heat_exchange() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn storage_basis(&self) -> ThermalStorageBasis {
            ThermalStorageBasis::InternalEnergy
        }
        fn temperature_k(&self, u: f64) -> Result<f64, DepositionError> {
            if !(0.0..=20.0).contains(&u) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(if u < 3.0 {
                2.0 + u
            } else if u <= 8.0 {
                5.0
            } else {
                5.0 + (u - 8.0) / 2.0
            })
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
    }
    let g = gas(1, false, 0.0);
    for conductance in [0.0, 0.3] {
        let mut n = network(&g, &[2.0], false);
        let original =
            MaterialWallCaloric::from_thermal_law(Plateau, cell(6.0), [0.0, 20.0]).unwrap();
        let mut wall = bind(
            &n,
            [original.clone(), original.clone()],
            conductance,
            Some(0.0),
            false,
        )
        .unwrap();
        let report = wall
            .advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        let loss = slip_loss(0.2);
        let q = 0.2 * conductance * 3.0 / (1.0 + 0.2 * conductance / 0.5);
        for i in 0..2 {
            let contact = wall.branch_contacts(0).unwrap()[i];
            assert_eq!(original.temperature_k(contact.thermal).unwrap(), 5.0);
            close(stored(contact), WALL_MASS * 6.0 + loss - q, 3e-11);
            close(
                g.primitive(n.branches()[0].cell(i).unwrap(), || false)
                    .unwrap()
                    .temperature,
                2.0 + q / 0.5,
                3e-10,
            );
            close(
                wall.accepted_branch_support_exchanges(0).unwrap()[i]
                    .unwrap()
                    .wall_slip_heating_j,
                loss,
                3e-12,
            );
        }
        close(report.wall_energy_change_j, 2.0 * (loss - q), 4e-11);
        close(report.heat_to_gas_j, 2.0 * (q - loss), 4e-11);
        close(report.energy_residual_j, 0.0, 4e-11);
    }
}

#[test]
fn g4_partition_cancellation_and_late_source_refusal_preserve_accepted_state_and_replay() {
    let g = gas(1, false, 0.0);
    for narrow_second_cell in [false, true] {
        let make = || {
            let mut n = network(&g, &[2.0], false);
            let range = if narrow_second_cell {
                [4.99999, 5.00001]
            } else {
                [0.01, 50.0]
            };
            let mut wall = bind(
                &n,
                [
                    source(ThermalStorageBasis::InternalEnergy, [0.01, 50.0]),
                    source(ThermalStorageBasis::InternalEnergy, range),
                ],
                0.3,
                Some(0.35),
                false,
            )
            .unwrap();
            wall.advance_for(&mut n, Time::new(1e-8), controls(), || false)
                .unwrap();
            (n, wall)
        };
        let (mut n, mut wall) = make();
        let (mut reference, mut expected) = make();
        let retry_dt = if narrow_second_cell { 1e-8 } else { 0.2 };
        let mut polls = 0;
        let expected_report = expected
            .advance_for(&mut reference, Time::new(retry_dt), controls(), || {
                polls += 1;
                false
            })
            .unwrap();
        let before = fields(&n);
        let contacts = wall.branch_contacts(0).unwrap().to_vec();
        let supports = wall.branch_supports(0).unwrap().to_vec();
        let exchanges = wall.accepted_branch_support_exchanges(0).unwrap().to_vec();
        let impulses = wall.accepted_branch_impulses().unwrap().to_vec();
        let accepted_report = wall.last_report();
        if narrow_second_cell {
            assert!(
                wall.advance_for(&mut n, Time::new(0.2), controls(), || false)
                    .is_err()
            );
        } else {
            let mut visits = 0;
            assert_eq!(
                wall.advance_for(&mut n, Time::new(0.2), controls(), || {
                    visits += 1;
                    visits == polls
                }),
                Err(WallExchangeError::Cancelled)
            );
        }
        assert_eq!(fields(&n), before);
        assert_eq!(n.time(), Time::new(0.0));
        assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
        assert_eq!(wall.branch_supports(0).unwrap(), supports);
        assert_eq!(
            wall.accepted_branch_support_exchanges(0).unwrap(),
            exchanges
        );
        assert_eq!(wall.accepted_branch_impulses().unwrap(), impulses);
        assert_eq!(wall.last_report(), accepted_report);
        assert_eq!(
            wall.advance_for(&mut n, Time::new(retry_dt), controls(), || false)
                .unwrap(),
            expected_report
        );
        assert_eq!(fields(&n), fields(&reference));
        assert_eq!(wall.branch_contacts(0), expected.branch_contacts(0));
        assert_eq!(wall.branch_supports(0), expected.branch_supports(0));
        assert_eq!(
            wall.accepted_branch_support_exchanges(0),
            expected.accepted_branch_support_exchanges(0)
        );
    }
}

#[test]
fn g0_partition_requires_valid_fractions_supports_and_original_finite_energy_storage() {
    let g = gas(1, false, 0.0);
    let n = network(&g, &[2.0], false);
    for fraction in [-0.1, 1.1, f64::NAN, f64::INFINITY] {
        let original = source(ThermalStorageBasis::InternalEnergy, [0.01, 50.0]);
        assert!(matches!(
            bind(&n, [original.clone(), original], 0.3, Some(fraction), false),
            Err(WallExchangeError::Invalid(_))
        ));
    }
    for basis in [
        ThermalStorageBasis::Enthalpy,
        ThermalStorageBasis::DeclaredSensibleStorage,
    ] {
        let original = source(basis, [0.01, 50.0]);
        assert!(matches!(
            bind(&n, [original.clone(), original], 0.3, Some(0.35), true),
            Err(WallExchangeError::Invalid(_))
        ));
    }
    assert!(
        uniform_operator(&n, lumped())
            .with_finite_support_heat_partition(vec![vec![0.35; 2]], || false)
            .is_err()
    );
    let original = source(ThermalStorageBasis::InternalEnergy, [0.01, 50.0]);
    for thermal in [
        WallThermal::Adiabatic,
        WallThermal::Reservoir {
            temperature_k: 5.0,
            conductance_w_per_k: 0.3,
        },
        original
            .initial_contact(0.3, WallDrag::None)
            .unwrap()
            .thermal,
    ] {
        let wall = uniform_operator(
            &n,
            WallContact {
                drag: WallDrag::Linear { rate_per_s: 3.0 },
                thermal,
            },
        )
        .with_finite_supports(vec![vec![Some(support()); 2]], 1e-11, || false)
        .unwrap();
        assert!(
            wall.with_finite_support_heat_partition(vec![vec![0.35; 2]], || false)
                .is_err()
        );
    }
    let wall = uniform_operator(&n, lumped())
        .with_finite_supports(vec![vec![Some(support()), None]], 1e-11, || false)
        .unwrap();
    assert!(
        wall.with_finite_support_heat_partition(vec![vec![0.35; 2]], || false)
            .is_err()
    );
}
