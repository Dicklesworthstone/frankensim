use super::*;

struct Sensible {
    cp: f64,
}
impl ThermalLaw for Sensible {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        let t = 300.0 + h / self.cp;
        if !(250.0..=2_000.0).contains(&t) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(t)
    }
    fn maximum_temperature_slope(&self) -> f64 {
        1.0 / self.cp
    }
}
static LAW: Sensible = Sensible { cp: 1_000.0 };
struct Latent;
impl ThermalLaw for Latent {
    fn temperature_k(&self, h: f64) -> Result<f64, DepositionError> {
        if !(0.0..=1_000.0).contains(&h) {
            return Err(DepositionError::MaterialDomain);
        }
        Ok(if h < 40.0 {
            300.0 + h / 4.0
        } else if h <= 140.0 {
            310.0
        } else {
            310.0 + (h - 140.0) / 8.0
        })
    }
    fn maximum_temperature_slope(&self) -> f64 {
        0.25
    }
}
static LATENT: Latent = Latent;
fn geometry(cells: usize) -> WallGeometry {
    WallGeometry {
        area: Area::new(2.0),
        thickness: Length::new(0.04),
        density: Density::new(500.0),
        conductivity: ThermalConductivity::new(0.2),
        cells,
    }
}
fn wall(cells: usize) -> Wall1d<'static> {
    Wall1d::new(&LAW, geometry(cells), SpecificEnthalpy::new(0.0), || false).unwrap()
}
fn contact(area: f64, temperature: f64, htc: f64) -> GasContact {
    GasContact {
        area: Area::new(area),
        temperature: Temperature::new(temperature),
        htc: HeatTransferCoefficient::new(htc),
    }
}
fn insulated() -> WallBoundary {
    WallBoundary {
        gas: [contact(0.0, 300.0, 0.0); 2],
        ambient_temperature: Temperature::new(300.0),
        ambient_htc: HeatTransferCoefficient::new(0.0),
        external_radiant_power: Power::new(0.0),
    }
}
fn conditions() -> WallBoundary {
    WallBoundary {
        gas: [contact(0.5, 400.0, 10.0), contact(1.5, 600.0, 20.0)],
        ambient_temperature: Temperature::new(280.0),
        ambient_htc: HeatTransferCoefficient::new(8.0),
        external_radiant_power: Power::new(50.0),
    }
}
fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() < tolerance,
        "actual={actual:.17e} expected={expected:.17e} tolerance={tolerance:e}"
    );
}
fn account(report: &WallReport) {
    let combined = report.audit.stored_energy_change_j
        + report
            .gas_heat_energy
            .iter()
            .map(|e| e.value())
            .sum::<f64>()
        - report.ambient_heat_energy.value()
        - report.external_radiant_energy.value()
        - report.coupled_radiant_energy.value();
    near(combined, report.audit.residual_j, 2e-8);
    assert!(report.audit.residual_j.abs() <= WallSolve::default().total_energy_tolerance_j);
    assert!(report.maximum_cell_residual_j <= WallSolve::default().cell_energy_tolerance_j);
    for i in 0..2 {
        assert_eq!(
            report.gas_heat_energy[i].value(),
            report.gas_heat_power[i].value() * report.elapsed.value()
        );
    }
}

#[test]
fn original_cells_material_binding_and_conductances_define_the_real_wall() {
    let wall = wall(4);
    let original = wall.simulation();
    assert!(original.uses_material_law(MaterialId(1), &LAW));
    assert_eq!(original.cells().len(), 4);
    assert_eq!(original.links().len(), 3);
    assert_eq!(original.mass_kg(), 40.0);
    for (i, cell) in original.cells().values().enumerate() {
        near(cell.mass_kg, 10.0, 1e-14);
        near(cell.geometry.volume_m3, 0.02, 1e-15);
        near(cell.geometry.start_m[0], 0.01 * i as f64, 1e-15);
        let expected_area = if i == 0 || i == 3 { 2.0 } else { 0.0 };
        assert_eq!(original.exposed_area_m2(cell.id).unwrap(), expected_area);
    }
    for link in original.links().values() {
        near(link.conductance_w_k, 40.0, 1e-13);
    }
    let unchanged = wall
        .candidate(
            Time::new(1e4),
            insulated(),
            WallSolve::default(),
            &mut wall.workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(unchanged.wall.simulation().cells(), original.cells());
    assert_eq!(unchanged.report.audit.stored_energy_change_j, 0.0);
    assert_eq!(unchanged.report.gas_heat_power, [Power::new(0.0); 2]);
    assert_eq!(unchanged.wall.time(), Time::new(1e4));
    assert_eq!(wall.time(), Time::new(0.0));
}

#[test]
fn signed_enclosure_radiation_heats_and_cools_without_becoming_external_fire_input() {
    let wall = wall(1);
    for power in [-80.0, 80.0] {
        let next = wall
            .candidate_with_exchange(
                Time::new(10.0),
                insulated(),
                Power::new(power),
                WallSolve::default(),
                &mut wall.workspace().unwrap(),
                || false,
            )
            .unwrap();
        near(next.report.audit.stored_energy_change_j, 10.0 * power, 1e-6);
        near(
            next.wall.cell_temperature(0).unwrap().value(),
            300.0 + power / 4000.0,
            1e-10,
        );
        near(
            next.report.surface_temperatures[0].value(),
            300.0 + power / 4000.0 + power / 20.0,
            1e-10,
        );
        assert_eq!(
            next.report.coupled_radiant_energy,
            Energy::new(10.0 * power)
        );
        assert_eq!(next.report.external_radiant_energy, Energy::new(0.0));
        assert_eq!(next.report.gas_heat_power, [Power::new(0.0); 2]);
        account(&next.report);
    }
    assert_eq!(wall.time(), Time::new(0.0));
    assert_eq!(wall.cell_temperature(0).unwrap(), Temperature::new(300.0));
}

#[test]
fn external_fire_and_internal_exchange_retain_separate_signed_accounts() {
    let wall = wall(1);
    let mut boundary = insulated();
    boundary.external_radiant_power = Power::new(50.0);
    let next = wall
        .candidate_with_exchange(
            Time::new(10.0),
            boundary,
            Power::new(-30.0),
            WallSolve::default(),
            &mut wall.workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(next.report.external_radiant_energy, Energy::new(500.0));
    assert_eq!(next.report.coupled_radiant_energy, Energy::new(-300.0));
    near(next.report.audit.stored_energy_change_j, 200.0, 1e-6);
    account(&next.report);
    for power in [f64::NAN, f64::INFINITY, f64::NEG_INFINITY] {
        assert!(matches!(
            wall.candidate_with_exchange(
                Time::new(10.0),
                boundary,
                Power::new(power),
                WallSolve::default(),
                &mut wall.workspace().unwrap(),
                || false,
            ),
            Err(WallError::Invalid(_))
        ));
    }
    assert_eq!(wall.time(), Time::new(0.0));
}

#[test]
fn one_cell_matches_independent_backward_euler_and_surface_half_resistances() {
    let wall = wall(1);
    let mut work = wall.workspace().unwrap();
    let duration = 500.0;
    let next = wall
        .candidate(
            Time::new(duration),
            conditions(),
            WallSolve::default(),
            &mut work,
            || false,
        )
        .unwrap();
    let capacitance = 40.0 * 1_000.0;
    let bulk = 2.0 * 0.2 * 2.0 / 0.04;
    let gas_h = 0.5 * 10.0 + 1.5 * 20.0;
    let outer_h = 2.0 * 8.0;
    let inner_g = bulk * gas_h / (bulk + gas_h);
    let outer_g = bulk * outer_h / (bulk + outer_h);
    let inner_drive = bulk / (bulk + gas_h) * (0.5 * 10.0 * 400.0 + 1.5 * 20.0 * 600.0 + 50.0);
    let expected = (capacitance * 300.0 + duration * (inner_drive + outer_g * 280.0))
        / (capacitance + duration * (inner_g + outer_g));
    near(
        next.wall.cell_temperature(0).unwrap().value(),
        expected,
        1e-9,
    );
    let inner_surface =
        (bulk * expected + 0.5 * 10.0 * 400.0 + 1.5 * 20.0 * 600.0 + 50.0) / (bulk + gas_h);
    let outer_surface = (bulk * expected + outer_h * 280.0) / (bulk + outer_h);
    near(
        next.report.surface_temperatures[0].value(),
        inner_surface,
        1e-9,
    );
    near(
        next.report.surface_temperatures[1].value(),
        outer_surface,
        1e-9,
    );
    near(
        next.report.gas_heat_power[0].value(),
        0.5 * 10.0 * (inner_surface - 400.0),
        1e-8,
    );
    near(
        next.report.gas_heat_power[1].value(),
        1.5 * 20.0 * (inner_surface - 600.0),
        1e-8,
    );
    near(
        next.report.ambient_heat_energy.value(),
        duration * outer_h * (280.0 - outer_surface),
        1e-5,
    );
    account(&next.report);
    assert_eq!(wall.cell_temperature(0).unwrap(), Temperature::new(300.0));
}

fn relaxation_error(steps: usize) -> f64 {
    let mut wall = wall(1);
    let mut work = wall.workspace().unwrap();
    let boundary = WallBoundary {
        gas: [contact(2.0, 600.0, 10.0), contact(0.0, 300.0, 0.0)],
        ..insulated()
    };
    let duration = 2_000.0 / steps as f64;
    for _ in 0..steps {
        let next = wall
            .candidate(
                Time::new(duration),
                boundary,
                WallSolve::default(),
                &mut work,
                || false,
            )
            .unwrap();
        account(&next.report);
        wall = next.wall;
    }
    // Half-wall conductance and film conductance are both 20 W/K, giving
    // G=10 W/K and C=40000 J/K. The continuous lumped limit is exponential.
    let exact = 600.0 - 300.0 * fs_math::det::exp(-10.0 * 2_000.0 / 40_000.0);
    (wall.cell_temperature(0).unwrap().value() - exact).abs()
}
#[test]
fn time_refinement_converges_to_the_analytic_wall_relaxation() {
    let coarse = relaxation_error(8);
    let fine = relaxation_error(16);
    assert!(coarse > 0.1);
    assert!(
        fine < 0.55 * coarse && fine > 0.45 * coarse,
        "coarse {coarse}, fine {fine}"
    );
}

#[test]
fn multicell_implicit_conduction_closes_wall_and_gas_energy_beyond_explicit_cfl() {
    let wall = wall(6);
    let boundary = WallBoundary {
        gas: [contact(2.0, 600.0, 20.0), contact(0.0, 300.0, 0.0)],
        ambient_htc: HeatTransferCoefficient::new(10.0),
        ..insulated()
    };
    let next = wall
        .candidate(
            Time::new(1_000.0),
            boundary,
            WallSolve::default(),
            &mut wall.workspace().unwrap(),
            || false,
        )
        .unwrap();
    account(&next.report);
    let temperatures: Vec<_> = (0..6)
        .map(|i| next.wall.cell_temperature(i).unwrap().value())
        .collect();
    assert!(temperatures.iter().all(|t| *t > 300.0 && *t < 600.0));
    assert!(temperatures.windows(2).all(|w| w[0] > w[1]));
    assert!(next.report.gas_heat_energy[0].value() < 0.0);
    assert!(next.report.ambient_heat_energy.value() < 0.0);
    assert_eq!(
        next.wall.simulation().mass_kg(),
        wall.simulation().mass_kg()
    );
    assert_eq!(next.wall.simulation().links(), wall.simulation().links());
    // Explicit conduction bound is about 44 s for these cell masses/links.
    assert_eq!(next.wall.time(), Time::new(1_000.0));
}

#[test]
fn prescribed_radiant_input_heats_wall_without_a_fictitious_gas_debit() {
    let wall = wall(3);
    let boundary = WallBoundary {
        external_radiant_power: Power::new(100.0),
        ..insulated()
    };
    let next = wall
        .candidate(
            Time::new(10.0),
            boundary,
            WallSolve::default(),
            &mut wall.workspace().unwrap(),
            || false,
        )
        .unwrap();
    account(&next.report);
    assert_eq!(next.report.gas_heat_power, [Power::new(0.0); 2]);
    assert_eq!(next.report.gas_heat_energy, [Energy::new(0.0); 2]);
    assert_eq!(next.report.external_radiant_energy, Energy::new(1_000.0));
    assert_eq!(next.report.audit.external_heat_j, 1_000.0);
    near(next.report.audit.stored_energy_change_j, 1_000.0, 1e-6);
    assert!(
        next.report.surface_temperatures[0].value()
            > next.wall.cell_temperature(0).unwrap().value()
    );
}

#[test]
fn latent_storage_uses_original_enthalpy_even_when_temperature_does_not_move() {
    let geometry = WallGeometry {
        area: Area::new(0.1),
        thickness: Length::new(0.1),
        density: Density::new(1_000.0),
        conductivity: ThermalConductivity::new(1.0),
        cells: 1,
    };
    let wall = Wall1d::new(&LATENT, geometry, SpecificEnthalpy::new(50.0), || false).unwrap();
    let next = wall
        .candidate(
            Time::new(20.0),
            WallBoundary {
                external_radiant_power: Power::new(10.0),
                ..insulated()
            },
            WallSolve::default(),
            &mut wall.workspace().unwrap(),
            || false,
        )
        .unwrap();
    assert_eq!(wall.cell_temperature(0).unwrap(), Temperature::new(310.0));
    assert_eq!(
        next.wall.cell_temperature(0).unwrap(),
        Temperature::new(310.0)
    );
    near(
        next.wall.simulation().cells()[&CellId(1)].specific_h_j_kg,
        70.0,
        1e-9,
    );
    near(next.report.audit.stored_energy_change_j, 200.0, 1e-8);
    account(&next.report);
}

#[test]
fn moving_layer_contact_areas_never_rescale_wall_mass_or_geometry() {
    let wall = wall(3);
    let first = WallBoundary {
        gas: [contact(1.8, 300.0, 10.0), contact(0.2, 600.0, 10.0)],
        ..insulated()
    };
    let second = WallBoundary {
        gas: [contact(0.2, 300.0, 10.0), contact(1.8, 600.0, 10.0)],
        ..insulated()
    };
    let mut work = wall.workspace().unwrap();
    let a = wall
        .candidate(
            Time::new(10.0),
            first,
            WallSolve::default(),
            &mut work,
            || false,
        )
        .unwrap();
    let b = wall
        .candidate(
            Time::new(10.0),
            second,
            WallSolve::default(),
            &mut work,
            || false,
        )
        .unwrap();
    assert!(b.report.audit.stored_energy_change_j > a.report.audit.stored_energy_change_j);
    for candidate in [&a.wall, &b.wall] {
        assert_eq!(
            candidate.simulation().mass_kg(),
            wall.simulation().mass_kg()
        );
        assert_eq!(candidate.geometry(), wall.geometry());
        for (id, old) in wall.simulation().cells() {
            assert_eq!(candidate.simulation().cells()[id].mass_kg, old.mass_kg);
            assert_eq!(candidate.simulation().cells()[id].geometry, old.geometry);
        }
    }
}

#[test]
fn cancelled_and_refused_wall_trials_leave_history_unchanged_and_reuse_is_deterministic() {
    let wall = wall(3);
    let before = wall.simulation().cells().clone();
    let mut work = wall.workspace().unwrap();
    let mut calls = 0;
    let expected = wall
        .candidate(
            Time::new(50.0),
            conditions(),
            WallSolve::default(),
            &mut work,
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    assert!(calls > 10);
    for stop in [1, calls / 2, calls] {
        let mut count = 0;
        let failed = wall.candidate(
            Time::new(50.0),
            conditions(),
            WallSolve::default(),
            &mut work,
            || {
                count += 1;
                count == stop
            },
        );
        assert!(matches!(failed, Err(WallError::Cancelled)));
        assert_eq!(wall.simulation().cells(), &before);
        assert_eq!(wall.time(), Time::new(0.0));
        let retried = wall
            .candidate(
                Time::new(50.0),
                conditions(),
                WallSolve::default(),
                &mut work,
                || false,
            )
            .unwrap();
        assert_eq!(retried.report, expected.report);
        assert_eq!(
            retried.wall.simulation().cells(),
            expected.wall.simulation().cells()
        );
    }
    let impossible = WallBoundary {
        external_radiant_power: Power::new(1e30),
        ..insulated()
    };
    assert!(
        wall.candidate(
            Time::new(1.0),
            impossible,
            WallSolve::default(),
            &mut work,
            || false
        )
        .is_err()
    );
    assert_eq!(wall.simulation().cells(), &before);
}

#[test]
fn invalid_geometry_source_and_overlapping_contacts_refuse() {
    for cells in [0, MAX_WALL_CELLS + 1] {
        assert!(Wall1d::new(&LAW, geometry(cells), SpecificEnthalpy::new(0.0), || false).is_err());
    }
    let mut bad = geometry(2);
    bad.conductivity = ThermalConductivity::new(0.0);
    assert!(Wall1d::new(&LAW, bad, SpecificEnthalpy::new(0.0), || false).is_err());
    assert!(Wall1d::new(&LAW, geometry(2), SpecificEnthalpy::new(-1e9), || false).is_err());
    let wall = wall(2);
    let mut work = wall.workspace().unwrap();
    let overlap = WallBoundary {
        gas: [contact(1.1, 400.0, 10.0); 2],
        ..insulated()
    };
    assert!(
        wall.candidate(
            Time::new(1.0),
            overlap,
            WallSolve::default(),
            &mut work,
            || false
        )
        .is_err()
    );
    let negative = WallBoundary {
        external_radiant_power: Power::new(-1.0),
        ..insulated()
    };
    assert!(
        wall.candidate(
            Time::new(1.0),
            negative,
            WallSolve::default(),
            &mut work,
            || false
        )
        .is_err()
    );
    assert!(
        wall.candidate(
            Time::new(0.0),
            insulated(),
            WallSolve::default(),
            &mut work,
            || false
        )
        .is_err()
    );
    let mut other_work = super::tests::wall(1).workspace().unwrap();
    assert!(
        wall.candidate(
            Time::new(1.0),
            insulated(),
            WallSolve::default(),
            &mut other_work,
            || false
        )
        .is_err()
    );
}
