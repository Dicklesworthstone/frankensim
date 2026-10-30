//! Original nonlinear/latent wall storage on the continuing gas transaction.
use super::*;
use crate::multicomponent::network::dissipating::{DissipatingNetwork1d, DissipatingStepControl};
use crate::multicomponent::network::dissipation::tests::law;
use crate::multicomponent::wall_exchange::{MaterialWallCaloric, exchange_cell_with_heat_storage};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionError, MaterialId, ThermalLaw,
    ThermalStorageBasis,
};

// Manufactured original source: e(T)=a*T²/2, hence heat capacity a*T.
// This is a constitutive test specimen, not measured material authority.
struct QuadraticStorage(f64);
impl ThermalLaw for QuadraticStorage {
    fn storage_basis(&self) -> ThermalStorageBasis {
        ThermalStorageBasis::InternalEnergy
    }
    fn temperature_k(&self, e: f64) -> Result<f64, DepositionError> {
        if !e.is_finite() || e <= 0.0 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok((2.0 * e / self.0).sqrt())
    }
    fn maximum_temperature_slope(&self) -> f64 {
        100.0
    }
    fn validate_cell(&self, cell: &DepositedCell) -> Result<(), DepositionError> {
        if cell.mass_kg != 1.0 || cell.geometry.volume_m3 != 0.1 {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(())
    }
}

fn cell(e: f64) -> DepositedCell {
    DepositedCell {
        id: CellId(17),
        segment_id: 3,
        body_id: 8,
        material_id: MaterialId(29),
        mass_kg: 1.0,
        specific_h_j_kg: e,
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
fn source(a: f64, range: [f64; 2]) -> MaterialWallCaloric<'static> {
    MaterialWallCaloric::from_thermal_law(QuadraticStorage(a), cell(0.5 * a * 25.0), range).unwrap()
}
fn bind(n: &MixtureNetwork1d, a: f64) -> NetworkWallOperator<'static> {
    uniform_operator(
        n,
        source(a, [0.01, 50.0])
            .initial_contact(0.3, WallDrag::None)
            .unwrap(),
    )
    .with_material_heat_storage(
        n.branches()
            .iter()
            .map(|b| {
                (0..b.cell_count())
                    .map(|_| Some(source(a, [0.01, 50.0])))
                    .collect()
            })
            .collect(),
        n.plenums()
            .iter()
            .map(|_| Some(source(a, [0.01, 50.0])))
            .collect(),
        || false,
    )
    .unwrap()
}
fn energy(contact: WallContact) -> f64 {
    let WallThermal::Caloric {
        mass_kg,
        specific_energy_j_kg,
        ..
    } = contact.thermal
    else {
        panic!("expected original caloric storage")
    };
    mass_kg * specific_energy_j_kg
}

#[test]
fn source_capacity_changes_actual_gas_and_wall_temperature_with_exact_nonlinear_storage() {
    let g = gas(1, false, 0.0);
    let mut heats = Vec::new();
    for a in [0.4, 0.8] {
        let mut n = network(&g, &[0.0], false);
        let mut wall = bind(&n, a);
        let report = wall
            .advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap();
        // Independent quadratic: Q=0.06*(Tw-2-2Q),
        // Q=a*(25-Tw²)/2. The positive root uses the stable quadratic form.
        let aa = 0.5 * a * 1.12;
        let bb = 0.06;
        let cc = -0.5 * a * 25.0 * 1.12 - 0.12;
        let tw = -2.0 * cc / (bb + (bb * bb - 4.0 * aa * cc).sqrt());
        let q = 0.06 * (tw - 2.0) / 1.12;
        let tg = g
            .primitive(n.branches()[0].cell(0).unwrap(), || false)
            .unwrap()
            .temperature;
        close(tg, 2.0 + 2.0 * q, 3e-10);
        let source = wall.branch_heat_storage(0).unwrap()[0].as_ref().unwrap();
        assert_eq!(source.initial_cell().unwrap().material_id, MaterialId(29));
        assert_eq!(source.storage_basis(), ThermalStorageBasis::InternalEnergy);
        close(
            source
                .temperature_k(wall.branch_contacts(0).unwrap()[0].thermal)
                .unwrap(),
            tw,
            3e-10,
        );
        close(
            energy(wall.branch_contacts(0).unwrap()[0]),
            0.5 * a * tw * tw,
            3e-10,
        );
        close(report.heat_to_gas_j, 2.0 * q, 3e-10);
        close(report.wall_energy_change_j, -2.0 * q, 3e-10);
        close(report.energy_residual_j, 0.0, 1e-11);
        assert_eq!(report.external_heat_to_gas_j, 0.0);
        heats.push(q);
    }
    assert!(
        heats[1] > heats[0],
        "larger finite wall capacity sustains more heat transfer"
    );
}

#[test]
fn latent_wall_keeps_its_original_energy_coordinate_without_apparent_capacity() {
    struct Plateau;
    impl ThermalLaw for Plateau {
        fn temperature_k(&self, e: f64) -> Result<f64, DepositionError> {
            if !(0.0..=20.0).contains(&e) {
                return Err(DepositionError::MaterialDomain);
            }
            Ok(if e < 3.0 {
                2.0 + e
            } else if e <= 8.0 {
                5.0
            } else {
                5.0 + (e - 8.0) / 2.0
            })
        }
        fn maximum_temperature_slope(&self) -> f64 {
            1.0
        }
    }
    let g = gas(1, false, 0.0);
    let n = network(&g, &[0.0], false);
    let source = MaterialWallCaloric::from_thermal_law(Plateau, cell(6.0), [0.0, 12.0]).unwrap();
    let update = exchange_cell_with_heat_storage(
        &g,
        n.branches()[0].cell(0).unwrap(),
        0.2,
        Time::new(0.2),
        source.initial_contact(0.3, WallDrag::None).unwrap(),
        &source,
        controls(),
        || false,
    )
    .unwrap();
    let q = 0.06 * 3.0 / 1.12;
    close(update.report.heat_to_gas_j, q, 3e-11);
    close(energy(update.contact), 6.0 - q, 3e-11);
    assert_eq!(source.temperature_k(update.contact.thermal).unwrap(), 5.0);
    assert_eq!(source.storage_basis(), ThermalStorageBasis::Enthalpy);
    close(update.report.energy_residual_j, 0.0, 1e-11);

    // Cross the entire latent interval in both directions. Independent
    // endpoint equations use the original solid/liquid branches after the
    // transition; the latent inventory is already contained in e.
    for (initial_e, gas_t, volume, expected_q) in [
        (2.0, 30.0, 0.2, -280.0 / 26.0),
        (10.0, 2.0, 2.0, 100.0 / 13.0),
    ] {
        let q = g.conserved(&[1.0], 0.0, gas_t).unwrap();
        let wall =
            MaterialWallCaloric::from_thermal_law(Plateau, cell(initial_e), [0.0, 20.0]).unwrap();
        let changed = exchange_cell_with_heat_storage(
            &g,
            q.as_cell(),
            volume,
            Time::new(10.0),
            wall.initial_contact(1.0, WallDrag::None).unwrap(),
            &wall,
            controls(),
            || false,
        )
        .unwrap();
        close(changed.report.heat_to_gas_j, expected_q, 2e-9);
        close(energy(changed.contact), initial_e - expected_q, 2e-9);
        close(changed.report.wall_energy_change_j, -expected_q, 2e-9);
        close(changed.report.energy_residual_j, 0.0, 1e-10);
        close(changed.report.thermal_residual_j, 0.0, 1e-10);
        assert!(if expected_q < 0.0 {
            energy(changed.contact) > 8.0
        } else {
            energy(changed.contact) < 3.0
        });
    }
}

#[test]
fn source_domain_and_final_cancellation_preserve_all_gas_and_wall_state_then_replay() {
    let g = gas(1, false, 0.0);
    let mut n = network(&g, &[0.0], true);
    let initial_contact = source(0.4, [0.01, 50.0])
        .initial_contact(0.3, WallDrag::None)
        .unwrap();
    let mut wall = uniform_operator(&n, initial_contact)
        .with_material_heat_storage(
            vec![vec![
                Some(source(0.4, [0.01, 50.0])),
                Some(source(0.4, [4.999, 5.001])),
            ]],
            vec![Some(source(0.4, [0.01, 50.0]))],
            || false,
        )
        .unwrap();
    wall.advance_for(&mut n, Time::new(1e-5), controls(), || false)
        .unwrap();
    let original = fields(&n);
    let contacts = wall.branch_contacts(0).unwrap().to_vec();
    let accepted_report = wall.last_report();
    assert!(
        wall.advance_for(&mut n, Time::new(0.2), controls(), || false)
            .is_err()
    );
    assert_eq!(fields(&n), original);
    assert_eq!(wall.branch_contacts(0).unwrap(), contacts);
    assert_eq!(wall.last_report(), accepted_report);

    let mut n = network(&g, &[0.0], true);
    let original = fields(&n);
    let mut good = bind(&n, 0.4);
    let mut reference = network(&g, &[0.0], true);
    let mut reference_wall = bind(&reference, 0.4);
    let mut polls = 0;
    let expected = reference_wall
        .advance_for(&mut reference, Time::new(0.2), controls(), || {
            polls += 1;
            false
        })
        .unwrap();
    let mut seen = 0;
    assert_eq!(
        good.advance_for(&mut n, Time::new(0.2), controls(), || {
            seen += 1;
            seen == polls
        }),
        Err(WallExchangeError::Cancelled)
    );
    assert_eq!(fields(&n), original);
    assert!(good.last_report().is_none());
    assert_eq!(
        good.advance_for(&mut n, Time::new(0.2), controls(), || false)
            .unwrap(),
        expected
    );
    assert_eq!(fields(&n), fields(&reference));
    assert_eq!(good.branch_contacts(0), reference_wall.branch_contacts(0));
    assert_eq!(good.plenum_contacts(), reference_wall.plenum_contacts());
}

#[test]
fn continuing_molecular_and_wave_composition_conserves_nonlinear_wall_storage() {
    let g = gas(1, false, 1000.0);
    let n = network(&g, &[0.2], false);
    let mut wall = bind(&n, 0.4);
    let mut flow =
        DissipatingNetwork1d::new(n, law(&g, 0.004, 0.0, 0.02, 0.0), vec![], || false).unwrap();
    let initial_gas: f64 = flow.branches()[0]
        .mesh()
        .volumes()
        .iter()
        .enumerate()
        .map(|(i, v)| {
            v * g
                .absolute_energy_density(flow.branches()[0].cell(i).unwrap())
                .unwrap()
        })
        .sum();
    let initial_wall: f64 = wall
        .branch_contacts(0)
        .unwrap()
        .iter()
        .copied()
        .map(energy)
        .sum();
    for _ in 0..5 {
        let report = flow
            .step_with_wall(
                &mut wall,
                Time::new(0.01),
                DissipatingStepControl::default(),
                controls(),
                || false,
            )
            .unwrap();
        let source_halves = report.wall.unwrap();
        assert_eq!(
            source_halves.first.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        assert_eq!(
            source_halves.second.elapsed.value(),
            0.5 * report.elapsed.value()
        );
        close(report.absolute_energy_residual_j, 0.0, 3e-11);
        let current_gas: f64 = flow.branches()[0]
            .mesh()
            .volumes()
            .iter()
            .enumerate()
            .map(|(i, v)| {
                v * g
                    .absolute_energy_density(flow.branches()[0].cell(i).unwrap())
                    .unwrap()
            })
            .sum();
        let current_wall: f64 = wall
            .branch_contacts(0)
            .unwrap()
            .iter()
            .copied()
            .map(energy)
            .sum();
        close(
            (current_gas - initial_gas) + (current_wall - initial_wall),
            0.0,
            5e-10,
        );
        assert_eq!(wall.last_report().unwrap().external_heat_to_gas_j, 0.0);
    }
    assert!(flow.time().value() > 0.0);
    assert!(
        g.primitive(flow.branches()[0].cell(0).unwrap(), || false)
            .unwrap()
            .temperature
            > 2.0
    );
}
