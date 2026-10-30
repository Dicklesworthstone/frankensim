use super::*;
use fs_couple::circuit_enthalpy::eddy::thermal::MagneticThermalLaw;
use fs_material::field_storage::thermal::magnetic::ThermalMagneticMaterial;

mod heat;

fn thermal_curves(cards: &[MaterialCard; 2]) -> [ResolvedHeatCapacityCurve; 2] {
    use fs_material::phase::{HeatCapacitySourceQuery, SolidLiquidPhase};
    use fs_material::state_point::SPECIFIC_HEAT_CAPACITY_DIMS;
    use fs_qty::Density;
    cards.each_ref().map(|card| {
        ResolvedHeatCapacityCurve::try_new(
            card,
            HeatCapacitySourceQuery {
                heat_capacity: &PropertyKey::new(
                    "specific_heat_capacity",
                    SPECIFIC_HEAT_CAPACITY_DIMS,
                ),
                density: &PropertyKey::new("density", Density::DIMS),
                temperature_axis: "T",
                lower: &field::absolute_point(200.0),
                upper: &field::absolute_point(2000.0),
                selection: fs_matdb::SelectionPolicy::SingleClaimOnly,
            },
            SolidLiquidPhase::Solid,
            0.0,
            1e-6,
        )
        .unwrap()
    })
}

fn problem(card: &MaterialCard) -> PlanarMagnetostatics {
    let mesh = TriComplex2::from_indexed_triangles(
        "eddy-enthalpy-absolute-temperature-mode",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(1.0).unwrap(),
    )
    .unwrap();
    PlanarMagnetostatics::new(
        mesh,
        vec![
            field::bundle_at(
                card,
                &field::absolute_point(300.0),
                ElectricalPropertyBasis::Resistivity
            );
            4
        ],
        field::limits(),
        &mut || false,
    )
    .unwrap()
}
fn laws(card: &MaterialCard, problem: &PlanarMagnetostatics) -> Vec<MagneticThermalLaw> {
    problem
        .materials()
        .iter()
        .enumerate()
        .map(|(triangle, material)| MagneticThermalLaw {
            triangle,
            law: ThermalMagneticMaterial::new(card, material, "T", [290.0, 360.0]).unwrap(),
            maximum_flux_density_t: 6.1,
        })
        .collect()
}

#[test]
fn permeability_entropy_and_spatial_joule_heat_share_one_actual_temperature() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let dt = 0.1;
    let mut temperatures = Vec::new();
    for slope in [-0.001, 0.0, 0.001] {
        let card = field::absolute_affine_magnetic_card(1.0 / 18.0, 0.0, 1.0, slope, basis);
        let curves = thermal_curves(&[card.clone(), card.clone()]);
        let thermal = single_thermal(&curves[0], 1.0);
        let problem = problem(&card);
        let state = initial_field(&problem, 1.0);
        let context = field::absolute_point(300.0);
        let cells = [MaterialCell {
            cell: CellId(7),
            material_id: MaterialId(9),
            curve: curves[0].curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        }];
        let routes = single_routes();
        let network = MagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: &routes,
            volume_tolerance_m3: 1e-12,
        };
        let bindings = laws(&card, &problem);
        let result = network
            .step_eddy_thermal_candidate(
                &thermal,
                dt,
                &state,
                &[0.0; 4],
                &bindings,
                common::config(),
                field_controls(),
                field::limits(),
                &[0.0],
                || false,
            )
            .unwrap();
        // Independently eliminate the one P1 field coordinate, then bracket
        // total physical U (the production path solves all enthalpy residuals).
        let amplitude = |t: f64| 3.0 * (10.0 - 2.0) / (10.0 + 2.0 / (1.0 + slope * (t - 300.0)));
        let energy = |t: f64, a: f64| {
            let mu = 1.0 + slope * (t - 300.0);
            2.0 * a * a / mu * (1.0 + t * slope / mu)
        };
        let initial_u = energy(300.0, 3.0);
        let mut lo = 290.0;
        let mut hi = 320.0;
        for _ in 0..100 {
            let mid = (lo + hi) / 2.0;
            if 4.0 * (mid - 300.0) + energy(mid, amplitude(mid)) > initial_u {
                hi = mid;
            } else {
                lo = mid;
            }
        }
        let t = result.thermal().temperature_k(CellId(7)).unwrap();
        common::near(t, (lo + hi) / 2.0);
        common::near(result.field().potential_wb_m()[4], amplitude(t));
        let report = result.report();
        common::near(report.total_energy_residual_j, 0.0);
        common::near(report.external_heat_j, 0.0);
        common::near(report.coupled.electrical.energy_residual_j, 0.0);
        let heat = report.coupled.deposited_joule_j;
        let rev = report.coupled.source.reversible_heat_j;
        common::near(4.0 * (t - 300.0), heat + rev);
        if slope == 0.0 {
            common::near(t, 302.5);
            common::near(rev, 0.0);
        } else {
            assert!(rev * slope > 0.0);
        }
        for sample in &report.coupled.source.fields {
            assert_eq!(sample.cell, CellId(7));
            common::near(sample.temperature_k[1], t);
            common::near(sample.flux_density_t[1], 2.0 * amplitude(t));
            common::near(sample.reversible_heat_j, rev / 4.0);
        }
        // Same original law can continue the accepted owned field/enthalpy
        // candidate. No initial-field reset or source replacement is needed.
        let replay = network
            .step_eddy_thermal_candidate(
                &thermal,
                dt,
                &state,
                &[0.0; 4],
                &bindings,
                common::config(),
                field_controls(),
                field::limits(),
                &[0.0],
                || false,
            )
            .unwrap();
        assert_eq!(
            replay.field().potential_wb_m(),
            result.field().potential_wb_m()
        );
        assert_eq!(replay.thermal().cells(), result.thermal().cells());
        let second = network
            .step_eddy_thermal_candidate(
                result.thermal(),
                2.0 * dt,
                result.field(),
                &[0.0; 4],
                &bindings,
                common::config(),
                field_controls(),
                field::limits(),
                &[0.0],
                || false,
            )
            .unwrap();
        common::near(second.report().total_energy_residual_j, 0.0);
        assert_eq!(second.field().time_s(), 2.0 * dt);
        assert_eq!(state.time_s(), 0.0);
        common::near(thermal.temperature_k(CellId(7)).unwrap(), 300.0);
        temperatures.push(t);
    }
    assert!(temperatures[0] < temperatures[1] && temperatures[1] < temperatures[2]);
}

#[test]
fn all_field_capacity_contributions_and_exact_source_bindings_are_required() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let card = field::absolute_affine_magnetic_card(1.0 / 18.0, 0.0, 1.0, 0.001, basis);
    let curves = thermal_curves(&[card.clone(), card.clone()]);
    let thermal = single_thermal(&curves[0], 1.0);
    let problem = problem(&card);
    let state = initial_field(&problem, 1.0);
    let context = field::absolute_point(300.0);
    let cells = [MaterialCell {
        cell: CellId(7),
        material_id: MaterialId(9),
        curve: curves[0].curve(),
        point: ThermalMaterialPoint {
            material: &card,
            context: &context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }];
    let routes = single_routes();
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let bindings = laws(&card, &problem);
    let run = |b: &[MagneticThermalLaw], config| {
        network.step_eddy_thermal_candidate(
            &thermal,
            0.1,
            &state,
            &[0.0; 4],
            b,
            config,
            field_controls(),
            field::limits(),
            &[0.0],
            || false,
        )
    };
    assert!(run(&bindings[..3], common::config()).is_err());
    let mut changed = bindings.clone();
    changed[3] = changed[2].clone();
    assert!(run(&changed, common::config()).is_err());
    let mut changed = bindings.clone();
    changed[0].maximum_flux_density_t = 5.9;
    assert!(run(&changed, common::config()).is_err());
    let mut changed = bindings.clone();
    for b in &mut changed {
        b.maximum_flux_density_t = 200.0;
    }
    assert!(run(&changed, common::config()).is_err());
    let foreign = field::absolute_affine_magnetic_card(1.0 / 18.0, 0.0, 1.0, 0.002, basis);
    let foreign_problem = self::problem(&foreign);
    let foreign_laws = laws(&foreign, &foreign_problem);
    assert!(run(&foreign_laws, common::config()).is_err());
    let mut exhausted = common::config();
    exhausted.maximum_circuit_trials = 1;
    assert!(run(&bindings, exhausted).is_err());
    let mut polls = 0;
    assert!(
        network
            .step_eddy_thermal_candidate(
                &thermal,
                0.1,
                &state,
                &[0.0; 4],
                &bindings,
                common::config(),
                field_controls(),
                field::limits(),
                &[0.0],
                || {
                    polls += 1;
                    polls > 100
                }
            )
            .is_err()
    );
    common::near(state.potential_wb_m()[4], 3.0);
    assert_eq!(state.time_s(), 0.0);
    common::near(thermal.temperature_k(CellId(7)).unwrap(), 300.0);
}

#[test]
fn spatial_temperatures_rebuild_local_permeability_without_averaging_the_body() {
    let basis = ElectricalPropertyBasis::Resistivity;
    let mu = [1.0, 2.0, 3.0, 4.0];
    let slopes = [0.001, -0.001, 0.002, 0.0];
    let sigma = [18.0, 9.0, 36.0, 3.0];
    let cards: Vec<_> = (0..4)
        .map(|i| {
            field::absolute_affine_magnetic_card(
                1.0 / sigma[i],
                0.0,
                mu[i],
                mu[i] * slopes[i],
                basis,
            )
        })
        .collect();
    let curves: Vec<_> = cards
        .iter()
        .map(|c| {
            let [first, _] = thermal_curves(&[c.clone(), c.clone()]);
            first
        })
        .collect();
    let context = field::absolute_point(300.0);
    let mesh = problem(&cards[0]).mesh().clone();
    let problem = PlanarMagnetostatics::new(
        mesh,
        cards
            .iter()
            .map(|c| field::bundle_at(c, &context, basis))
            .collect(),
        field::limits(),
        &mut || false,
    )
    .unwrap();
    let bc: Vec<_> = (0..4).map(|v| (v, 0.0)).collect();
    let k = mu.iter().map(|x| 1.0 / x).sum::<f64>();
    let initial = problem
        .solve(&[9.0 * k; 4], &bc, field::linear(), &mut || false)
        .unwrap();
    let state = PlanarEddy::new(
        &initial,
        sources(&problem),
        &[],
        &bc,
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap();
    common::near(state.potential_wb_m()[4], 3.0);
    let events = (0..4)
        .map(|i| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(i as u64),
                segment_id: i as u64,
                body_id: 0,
                material_id: MaterialId(i as u64),
                mass_kg: 0.25,
                specific_h_j_kg: 400.0,
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: 0.25,
                    surface_area_m2: 1.0,
                },
                boundary: common::boundary(0.0),
            }),
        })
        .collect();
    let caloric: Vec<_> = curves
        .iter()
        .enumerate()
        .map(|(i, c)| (MaterialId(i as u64), c.curve() as &dyn ThermalLaw))
        .collect();
    let thermal = DepositionSimulation::new(
        events,
        &caloric,
        DepositionBudget {
            max_cells: 4,
            max_links: 1,
            max_events: 4,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 0.1,
            energy_tolerance_j: 1e-7,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    let cells: Vec<_> = (0..4)
        .map(|i| MaterialCell {
            cell: CellId(i as u64),
            material_id: MaterialId(i as u64),
            curve: curves[i].curve(),
            point: ThermalMaterialPoint {
                material: &cards[i],
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: None,
        })
        .collect();
    let routes: Vec<_> = (0..4)
        .map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: CellId(triangle as u64),
        })
        .collect();
    let bindings: Vec<_> = (0..4)
        .map(|i| MagneticThermalLaw {
            triangle: i,
            law: ThermalMagneticMaterial::new(
                &cards[i],
                &problem.materials()[i],
                "T",
                [290.0, 360.0],
            )
            .unwrap(),
            maximum_flux_density_t: 6.1,
        })
        .collect();
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let next = network
        .step_eddy_thermal_candidate(
            &thermal,
            0.1,
            &state,
            &[0.0; 4],
            &bindings,
            common::config(),
            field_controls(),
            field::limits(),
            &[0.0; 4],
            || false,
        )
        .unwrap();
    let a = next.field().potential_wb_m()[4];
    let mut endpoint_stiffness = 0.0;
    for i in 0..4 {
        let t = next.thermal().temperature_k(CellId(i as u64)).unwrap();
        let actual_mu = mu[i] * (1.0 + slopes[i] * (t - 300.0));
        endpoint_stiffness += 1.0 / actual_mu;
        common::near(next.field().problem().permeability_h_m()[i], actual_mu);
        let sample = &next.report().coupled.source.fields[i];
        let ohmic = sigma[i] * (a - 3.0) * (a - 3.0) / (72.0 * 0.1);
        common::near(t - 300.0, ohmic + sample.reversible_heat_j);
        common::near(next.report().coupled.electrical.cells[i].heat_j, ohmic);
    }
    let c = sigma.iter().sum::<f64>() / 72.0;
    common::near(
        a,
        3.0 * (c / 0.1 - k / 2.0) / (c / 0.1 + endpoint_stiffness / 2.0),
    );
    assert_ne!(
        next.thermal().temperature_k(CellId(0)).unwrap(),
        next.thermal().temperature_k(CellId(1)).unwrap()
    );
    common::near(next.report().total_energy_residual_j, 0.0);
}
