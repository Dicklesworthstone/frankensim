use super::*;
use crate::multicomponent::Conserved;
use crate::multicomponent::tests::{gas, source};
use crate::multicomponent::wall_exchange::{
    WallControl, exchange_cell, exchange_cell_with_heat_storage,
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;

struct Capacity(ThermalStorageBasis);
impl ThermalLaw for Capacity {
    fn storage_basis(&self) -> ThermalStorageBasis {
        self.0
    }
    fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
        if !energy.is_finite() || energy <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        // Exact inverse of u(T)=0.7*T+0.15*T², c(T)=0.7+0.3*T.
        Ok(2.0 * energy / (0.7 + (0.49 + 0.6 * energy).sqrt()))
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / 0.7
    }
}

fn storage(energy: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(17),
        segment_id: 3,
        body_id: 8,
        material_id: MaterialId(29),
        mass_kg: 0.6,
        specific_h_j_kg: energy,
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
    }
}

fn wall() -> MaterialWallCaloric<'static> {
    MaterialWallCaloric::from_thermal_law(
        Capacity(ThermalStorageBasis::InternalEnergy),
        storage(3.45),
        [0.0715, 22.0],
    )
    .unwrap()
}

fn fixture() -> (MixtureGas, [Conserved; 2]) {
    let g = gas(vec![
        source("A", R, 3.5, 0.1, 1e10),
        source("B", 2.0 * R, 4.3, 0.05, -2e10),
    ]);
    let states = [
        g.conserved(&[0.6, 0.4], 0.9, 1.6).unwrap(),
        g.conserved(&[1.0, 0.5], -0.3, 5.0).unwrap(),
    ];
    (g, states)
}

fn inputs(states: &[Conserved; 2], conductance: [f64; 2]) -> [SharedRecoilHeatCell<'_>; 2] {
    core::array::from_fn(|i| SharedRecoilHeatCell {
        state: states[i].as_cell(),
        volume_m3: [0.4, 0.8][i],
        law: RecoilLaw {
            linear_rate_per_s: 0.8,
            quadratic_per_m: 0.0,
            gas_heat_fraction: [0.25, 0.75][i],
        },
        conductance_w_per_k: conductance[i],
    })
}

fn balance() -> RecoilControl {
    RecoilControl {
        wall: WallControl {
            absolute_energy_tolerance_j: 1e-11,
            relative_energy_tolerance: 1e-12,
            maximum_iterations: 96,
        },
        momentum_tolerance_ns: 1e-12,
    }
}

fn support() -> SlidingSupport {
    SlidingSupport {
        mass_kg: 0.6,
        velocity_m_s: -0.2,
    }
}

fn advance<'a>(
    workspace: &'a mut SharedRecoilHeatWorkspace<'_>,
    gas: &MixtureGas,
    cells: &[SharedRecoilHeatCell<'_>],
    thermal: WallThermal,
    cancelled: impl FnMut() -> bool,
) -> Result<SharedRecoilHeatUpdate<'a>, WallExchangeError> {
    workspace.advance(
        gas,
        cells,
        Time::new(0.2),
        support(),
        thermal,
        SharedRecoilControl::default(),
        balance(),
        SharedWallControl::default(),
        cancelled,
    )
}

fn initial(wall: &MaterialWallCaloric<'_>) -> WallThermal {
    wall.initial_contact(0.0, WallDrag::None).unwrap().thermal
}

fn close(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 2e-10,
        "{actual:e} != {expected:e}"
    );
}

#[test]
fn shared_friction_heats_one_nonlinear_inventory_by_analytic_integral() {
    let (g, mut states) = fixture();
    states[1] = g.conserved(&[1.0, 0.5], 0.9, 5.0).unwrap();
    let mut cells = inputs(&states, [0.0; 2]);
    for cell in &mut cells {
        cell.law.quadratic_per_m = 0.7;
    }
    let wall = wall();
    let thermal = initial(&wall);
    let mut workspace = SharedRecoilHeatWorkspace::new(2, wall.clone(), thermal, || false).unwrap();
    let update = advance(&mut workspace, &g, &cells, thermal, || false).unwrap();
    // Equal gas velocities and laws form one invariant group, with an exact
    // aggregate reduced-mass solution independent of the production solver.
    let mass = 1.6;
    let exponent: f64 = 0.8 * 0.2 * (1.0 + mass / 0.6);
    let slip = 1.1 / (exponent.exp() + 0.7 * 1.1 / 0.8 * exponent.exp_m1());
    let loss = 0.5 * (mass * 0.6 / (mass + 0.6)) * (1.1 * 1.1 - slip * slip);
    let wall_heat = loss * (0.4 * 0.75 + 1.2 * 0.25) / mass;
    let capacity = 0.6 * (0.7 + 0.3 * 3.0);
    let dt =
        2.0 * wall_heat / (capacity + (capacity * capacity + 2.0 * 0.6 * 0.3 * wall_heat).sqrt());
    close(update.wall_slip_heating_j.iter().sum(), wall_heat);
    close(update.wall_energy_change_j, wall_heat);
    close(wall.temperature_k(update.thermal).unwrap(), 3.0 + dt);
    assert!((dt - wall_heat / capacity).abs() > 1e-5);
    assert_eq!(update.heat_to_gas_j, [0.0; 2]);
    assert_eq!(
        update.storage_constraint,
        WallCaloricConstraint::IsochoricInternalEnergy
    );
    close(update.energy_residual_j, 0.0);
}

#[test]
fn variable_cv_cells_share_one_material_endpoint_and_conserve_all_storage() {
    let (g, states) = fixture();
    let before = states.clone();
    let cells = inputs(&states, [0.3, 0.7]);
    let wall = wall();
    let thermal = initial(&wall);
    let mut workspace = SharedRecoilHeatWorkspace::new(2, wall.clone(), thermal, || false).unwrap();
    let update = advance(&mut workspace, &g, &cells, thermal, || false).unwrap();
    let wall_t = wall.temperature_k(update.thermal).unwrap();
    let wall_change = 0.6 * (0.7 * (wall_t - 3.0) + 0.15 * (wall_t - 3.0) * (wall_t + 3.0));
    let mut transferred = 0.0;
    let mut total_gas = 0.0;
    for (i, cell) in cells.iter().enumerate() {
        let actual = Cell {
            partial_densities: cell.state.partial_densities,
            momentum: update.mechanics.cells[i].momentum,
            energy: update.energies[i],
        };
        let temperature = g.primitive(actual, || false).unwrap().temperature;
        let initial_t = [2.0, 4.0][i];
        let [a, b] = [
            cell.state.partial_densities[0],
            cell.state.partial_densities[1],
        ];
        let capacity_a = 2.5 * a + 1.65 * b;
        let capacity_b = 0.1 * a + 0.025 * b;
        let sensible = cell.volume_m3
            * (temperature - initial_t)
            * (capacity_a + 0.5 * capacity_b * (temperature + initial_t));
        let flux = 0.2 * cell.conductance_w_per_k * (wall_t - temperature);
        close(
            sensible,
            cell.law.gas_heat_fraction * update.mechanics.cells[i].slip_dissipation_j + flux,
        );
        close(update.heat_to_gas_j[i], flux);
        total_gas += (actual.energy - cell.state.energy) * cell.volume_m3;
        transferred += flux;
    }
    close(
        wall_change,
        update.wall_slip_heating_j.iter().sum::<f64>() - transferred,
    );
    close(update.wall_energy_change_j, wall_change);
    close(
        total_gas + update.mechanics.support_kinetic_energy_change_j + wall_change,
        0.0,
    );
    close(update.energy_residual_j, 0.0);
    close(update.mechanics.momentum_residual_ns, 0.0);
    assert_eq!(states, before);
    let energies = [update.energies[0], update.energies[1]];
    let reversed = [cells[1], cells[0]];
    let reverse = advance(&mut workspace, &g, &reversed, thermal, || false).unwrap();
    close(reverse.energies[1], energies[0]);
    close(reverse.energies[0], energies[1]);
    close(wall.temperature_k(reverse.thermal).unwrap(), wall_t);
}

#[test]
fn strict_heat_budget_refines_default_gas_inversion_at_stiff_contact_endpoints() {
    let (g, states) = fixture();
    let state = states[1].as_cell();
    let control = balance().wall;
    let wall = wall();
    // Keep the original gas's default 1e-10 J/kg inversion budget. The heat
    // interface needs finer temperature resolution, including after return.
    assert_eq!(g.domain().absolute_energy_tolerance, 1e-10);
    let reservoir = exchange_cell(
        &g,
        state,
        0.8,
        Time::new(0.2),
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Reservoir {
                temperature_k: 3.0,
                conductance_w_per_k: 1000.0,
            },
        },
        control,
        || false,
    )
    .unwrap();
    let finite = exchange_cell_with_heat_storage(
        &g,
        state,
        0.8,
        Time::new(0.2),
        wall.initial_contact(1000.0, WallDrag::None).unwrap(),
        &wall,
        control,
        || false,
    )
    .unwrap();
    for update in [reservoir, finite] {
        let thermal = update.energy - 0.5 * state.momentum * (state.momentum / state.density());
        // Independent inverse of this composition's exact NASA polynomial,
        // including the unchanged source's 0.001 K sensible-energy reference.
        let absolute = thermal + 3.325 * 0.001 + 0.5 * 0.1125 * 0.001 * 0.001;
        let temperature = 2.0 * absolute / (3.325 + (3.325 * 3.325 + 0.225 * absolute).sqrt());
        let wall_temperature = match update.contact.thermal {
            WallThermal::Reservoir { temperature_k, .. } => temperature_k,
            thermal => wall.temperature_k(thermal).unwrap(),
        };
        let heat = (update.energy - state.energy) * 0.8;
        close(heat, 200.0 * (wall_temperature - temperature));
        assert!(update.report.thermal_residual_j.abs() <= control.tolerance(heat.abs()));
        assert!(update.report.energy_residual_j.abs() <= control.tolerance(heat.abs()));
    }
}

#[test]
fn wall_precision_uses_its_own_iteration_budget_and_preserves_negligible_heat() {
    let (original, _) = fixture();
    let mut domain = original.domain();
    domain.absolute_energy_tolerance = 1000.0;
    domain.maximum_iterations = 1;
    let gas = MixtureGas::new(original.source().clone(), domain, || false).unwrap();
    let state = gas.conserved(&[1.0, 0.0], 0.0, 2.0).unwrap();
    assert!(
        (gas.primitive(state.as_cell(), || false)
            .unwrap()
            .temperature
            - 2.0)
            .abs()
            > 0.5
    );
    let wall = wall();
    let contact = wall.initial_contact(1e-20, WallDrag::None).unwrap();
    let update = exchange_cell_with_heat_storage(
        &gas,
        state.as_cell(),
        0.2,
        Time::new(0.2),
        contact,
        &wall,
        WallControl {
            absolute_energy_tolerance_j: 1e-12,
            ..balance().wall
        },
        || false,
    )
    .unwrap();
    assert_eq!(update.energy, state.energy);
    assert_eq!(update.momentum, state.momentum);
    assert_eq!(update.contact, contact);
    assert_eq!(update.report.heat_to_gas_j, 0.0);
    assert_eq!(update.report.wall_energy_change_j, 0.0);
    assert!(update.report.thermal_residual_j.abs() <= 1e-12);
    assert_eq!(gas.domain(), domain);
    // Exact source fixed points need no tighter inversion even when G*dt/Cv
    // makes a hypothetical finite-error allowance underflow to zero.
    let equilibrium = exchange_cell(
        &original,
        Cell {
            partial_densities: &[1.0, 0.0],
            momentum: 0.0,
            energy: 0.0,
        },
        1.0,
        Time::new(1.0),
        WallContact {
            drag: WallDrag::None,
            thermal: WallThermal::Reservoir {
                temperature_k: original.domain().minimum_temperature,
                conductance_w_per_k: 1e308,
            },
        },
        balance().wall,
        || false,
    )
    .unwrap();
    assert_eq!(equilibrium.energy, 0.0);
    assert_eq!(equilibrium.report.heat_to_gas_j, 0.0);
    assert_eq!(equilibrium.report.thermal_residual_j, 0.0);
}

#[test]
fn shared_latent_inventory_retains_friction_and_simultaneous_contact_heat() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn storage_basis(&self) -> ThermalStorageBasis {
            ThermalStorageBasis::InternalEnergy
        }
        fn temperature_k(&self, energy: f64) -> Result<f64, DepositionError> {
            if !(0.0..=10.0).contains(&energy) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(if energy < 3.0 {
                1.5 + 0.5 * energy
            } else if energy <= 8.0 {
                3.0
            } else {
                energy - 5.0
            })
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
    }
    let (g, states) = fixture();
    let cells = inputs(&states, [0.3, 0.7]);
    let wall = MaterialWallCaloric::from_thermal_law(Plateau, storage(5.0), [0.0, 10.0]).unwrap();
    let thermal = initial(&wall);
    let mut workspace = SharedRecoilHeatWorkspace::new(2, wall.clone(), thermal, || false).unwrap();
    let update = advance(&mut workspace, &g, &cells, thermal, || false).unwrap();
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = update.thermal
    else {
        panic!()
    };
    assert_eq!(wall.temperature_k(update.thermal).unwrap(), 3.0);
    assert!((specific_energy_j_kg - 5.0).abs() > 0.01);
    close(
        mass_kg * (specific_energy_j_kg - 5.0),
        update.wall_energy_change_j,
    );
    close(
        update.wall_energy_change_j,
        update.wall_slip_heating_j.iter().sum::<f64>() - update.heat_to_gas_j.iter().sum::<f64>(),
    );
    close(update.energy_residual_j, 0.0);
    assert!(update.thermal_report.maximum_contact_residual_j < 1e-10);
}

#[test]
fn all_gas_partition_recovers_original_shared_motion_and_zero_time_identity() {
    let (g, states) = fixture();
    let mut cells = inputs(&states, [0.0; 2]);
    for cell in &mut cells {
        cell.law.gas_heat_fraction = 1.0;
    }
    let wall = wall();
    let thermal = initial(&wall);
    let mut workspace = SharedRecoilHeatWorkspace::new(2, wall, thermal, || false).unwrap();
    let mut mechanical = SharedRecoilWorkspace::new(2, || false).unwrap();
    let raw = cells.map(|cell| SharedRecoilCell {
        density_kg_m3: cell.state.density(),
        momentum_kg_m2_s: cell.state.momentum,
        volume_m3: cell.volume_m3,
        linear_rate_per_s: cell.law.linear_rate_per_s,
        quadratic_per_m: cell.law.quadratic_per_m,
    });
    let expected = mechanical
        .advance(
            &raw,
            Time::new(0.2),
            support(),
            SharedRecoilControl::default(),
            balance(),
            || false,
        )
        .unwrap();
    let update = advance(&mut workspace, &g, &cells, thermal, || false).unwrap();
    assert_eq!(update.mechanics.cells, expected.cells);
    assert_eq!(update.mechanics.support, expected.support);
    assert_eq!(update.thermal, thermal);
    assert_eq!(update.wall_energy_change_j, 0.0);
    for (i, cell) in cells.iter().enumerate() {
        assert_eq!(
            update.energies[i],
            cell.state.energy + expected.cells[i].work_to_gas_j / cell.volume_m3
        );
    }
    let zero = workspace
        .advance(
            &g,
            &cells,
            Time::new(0.0),
            support(),
            thermal,
            SharedRecoilControl::default(),
            balance(),
            SharedWallControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(zero.energies, states.each_ref().map(|s| s.energy));
    assert_eq!(zero.mechanics.support, support());
    assert_eq!(zero.thermal, thermal);
}

#[test]
fn shared_storage_domain_constraint_budget_and_late_cancel_refuse_with_exact_retry() {
    let (g, states) = fixture();
    let cells = inputs(&states, [0.3, 0.7]);
    for basis in [
        ThermalStorageBasis::Enthalpy,
        ThermalStorageBasis::DeclaredSensibleStorage,
    ] {
        let source =
            MaterialWallCaloric::from_thermal_law(Capacity(basis), storage(3.45), [0.0715, 22.0])
                .unwrap();
        assert!(
            SharedRecoilHeatWorkspace::new(2, source.clone(), initial(&source), || false).is_err()
        );
    }
    let narrow = MaterialWallCaloric::from_thermal_law(
        Capacity(ThermalStorageBasis::InternalEnergy),
        storage(3.45),
        [3.44, 3.4500001],
    )
    .unwrap();
    let mut narrow_workspace =
        SharedRecoilHeatWorkspace::new(2, narrow.clone(), initial(&narrow), || false).unwrap();
    assert!(
        advance(&mut narrow_workspace, &g, &cells, initial(&narrow), || {
            false
        })
        .is_err()
    );
    let wall = wall();
    let thermal = initial(&wall);
    let mut workspace = SharedRecoilHeatWorkspace::new(2, wall, thermal, || false).unwrap();
    let exhausted = workspace.advance(
        &g,
        &cells,
        Time::new(0.2),
        support(),
        thermal,
        SharedRecoilControl::default(),
        balance(),
        SharedWallControl {
            maximum_source_evaluations: 1,
        },
        || false,
    );
    assert!(matches!(
        exhausted,
        Err(WallExchangeError::Invalid(
            "shared-wall nested source work budget"
        ))
    ));
    let mut polls = 0;
    let accepted = advance(&mut workspace, &g, &cells, thermal, || {
        polls += 1;
        false
    })
    .unwrap();
    let expected = (
        accepted.energies.to_vec(),
        accepted.thermal,
        accepted.mechanics.support,
        accepted.energy_residual_j,
        accepted.thermal_report,
    );
    let mut visits = 0;
    assert!(matches!(
        advance(&mut workspace, &g, &cells, thermal, || {
            visits += 1;
            visits == polls
        }),
        Err(WallExchangeError::Cancelled)
    ));
    let retry = advance(&mut workspace, &g, &cells, thermal, || false).unwrap();
    assert_eq!(
        (
            retry.energies.to_vec(),
            retry.thermal,
            retry.mechanics.support,
            retry.energy_residual_j,
            retry.thermal_report
        ),
        expected
    );
}
