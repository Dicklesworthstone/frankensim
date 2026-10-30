//! Independent endpoint chemistry/heat references on the native spatial solid.
use super::*;
use fs_couple::cavity_phs::piston::reactive::implicit::coupled::ImplicitMechanismSolve;
use fs_solver::NewtonKrylovConfig;

const R: f64 = fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K;

pub(super) fn implicit(law: GasReaction) -> GasReaction {
    let policy = ImplicitMechanismSolve {
        extent_brackets_mol: law
            .channels()
            .map(|(channel, _)| {
                (
                    channel
                        .reverse_rate_closure()
                        .equilibrium()
                        .reaction()
                        .clone(),
                    [-1.0, 1.0],
                )
            })
            .collect(),
        extent_tolerance_mol: 1e-11,
        interior_fraction: 0.05,
        jacobian_relative_step: 1e-5,
        maximum_iterations: 32,
        maximum_evaluations: 1024,
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-12,
            relative_tolerance: 1e-12,
            linear_restart: 8,
            max_linear_cycles: 8,
            forcing_maximum: 1e-6,
            ..Default::default()
        },
    };
    law.with_implicit_mechanism(policy).unwrap()
}

fn amount() -> f64 {
    100000.0 * VOLUME / (R * 300.0)
}
fn capacity() -> f64 {
    2.5 * 100000.0 * VOLUME / 300.0
}

// Exact scalar BE species elimination for these authored equal-mass isomers:
// k=k0*T/300, activity_i=(n_i/n)*T/300, K=exp(release/(RT)).
// No production rate, equilibrium, caloric inversion or nonlinear solver is used.
pub(super) fn extent(
    temperature: f64,
    initial_product: f64,
    rate: f64,
    release: f64,
    dt: f64,
) -> f64 {
    let a = dt * VOLUME * rate * (temperature / 300.0).powi(2);
    let reverse = (-release / (R * temperature)).exp();
    a * (1.0 - (1.0 + reverse) * initial_product / amount())
        / (1.0 + a * (1.0 + reverse) / amount())
}

#[test]
fn heated_and_cooled_isomers_match_analytic_be_and_signed_entropy_with_zero_heat_parity() {
    let rate = 1.0;
    let law = implicit(reaction(rate, 0.0, 1.0));
    let initial = ClosedGasParcel::from_isochoric(gas(&law));
    let saved = initial.clone();
    for heat in [-2000.0, 0.0, 2000.0] {
        let candidate = law
            .heated_implicit_candidate(&initial, 1.0, heat, chemical(), || false)
            .unwrap();
        let t = 300.0 + heat / capacity();
        let xi = extent(t, 0.0, rate, 0.0, 1.0);
        let fraction = xi / amount();
        let entropy = capacity() * (t / 300.0).ln()
            - R * amount() * ((1.0 - fraction) * (1.0 - fraction).ln() + fraction * fraction.ln());
        let report = candidate.step();
        near(candidate.gas().state().state().temperature, t, 1e-8);
        near(
            candidate.gas().state().state().pressure,
            100000.0 * t / 300.0,
            1e-5,
        );
        near(report.reaction.extent_mol, xi, 1e-10);
        near(report.reaction.species[0].final_mol, amount() - xi, 1e-10);
        near(report.reaction.species[1].final_mol, xi, 1e-10);
        near(report.heat.supplied_heat_j, heat, 0.0);
        near(report.heat.stored_energy_change_j, heat, 2e-9);
        near(report.reaction.entropy_change_j_k, entropy, 1e-8);
        near(report.endpoint_heat_entropy_j_k, heat / t, 1e-10);
        near(report.entropy_residual_j_k, entropy - heat / t, 1e-8);
        assert!(report.entropy_residual_j_k >= 0.0);
        near(
            report
                .reaction
                .rate_material
                .as_ref()
                .unwrap()
                .property(ACTIVITY_RATE_SCALE_PROPERTY)
                .unwrap()
                .value_si(),
            rate * t / 300.0,
            1e-10,
        );
        if heat < 0.0 {
            assert!(
                report.reaction.entropy_change_j_k < 0.0,
                "cooling entropy must not be rejected as negative chemical production"
            );
        }
        if heat == 0.0 {
            let adiabatic = law
                .parcel_candidate(&initial, 1.0, chemical(), || false)
                .unwrap();
            assert_eq!(candidate.gas(), adiabatic.gas());
            assert_eq!(&report.reaction, adiabatic.step());
        }
    }
    assert_eq!(initial, saved);
}

// Eliminate BOTH solid temperatures by their independent 2x2 linear heat
// balance, then bracket the one remaining total-energy equation in gas T.
// The sole bisection here is a bounded test oracle for this specific specimen.
fn joint_reference(wall_k: f64, rate: f64, release: f64, product_k: f64) -> [f64; 4] {
    let sample = |t: f64| {
        let xi = extent(t, 0.0, rate, release, 1.0);
        let kg = 0.03 * (1.0 + (product_k - 1.0) * xi / amount());
        let [g1, g2] = [0.1, 0.2].map(|area| area / (0.001 / kg + 0.02 / wall_k));
        let [a, b] = [10.5 + g1, 10.5 + g2];
        let [r1, r2] = [3000.0 + g1 * t, 3000.0 + g2 * t];
        let determinant = a * b - 0.25;
        let t1 = (b * r1 + 0.5 * r2) / determinant;
        let t2 = (0.5 * r1 + a * r2) / determinant;
        let residual = capacity() * (t - 300.0) - release * xi + 10.0 * (t1 + t2 - 600.0);
        (residual, [t1, t2, t, xi])
    };
    let (mut lo, mut hi) = (300.0, 500.0);
    assert!(sample(lo).0 < 0.0 && sample(hi).0 > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if sample(mid).0 > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    sample(0.5 * (lo + hi)).1
}

#[test]
fn joint_heat_and_extent_follow_wall_data_and_independent_native_solid_reference() {
    with_cx(|cx, _| {
        let (rate, release, product_k) = (1.0, 20000.0, 4.0);
        let law = implicit(reaction(rate, release, product_k));
        let gas = gas(&law);
        let mut extents = Vec::new();
        for wall_k in [0.1, 1.0, 10.0] {
            let (card, caloric) = wall(wall_k, false);
            let thermal = initial(caloric.curve());
            let context = point(300.0);
            let model = solid(&thermal, &card, caloric.curve(), &context, cx);
            let next = model
                .step_with_coupled_reacting_isochoric_gas(
                    &thermal,
                    &gas,
                    &law,
                    1.0,
                    0.0,
                    &contacts(),
                    None,
                    chemical(),
                    solve(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            let expected = joint_reference(wall_k, rate, release, product_k);
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                expected[0],
                2e-7,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                expected[1],
                2e-7,
            );
            near(next.gas().state().state().temperature, expected[2], 2e-7);
            near(next.reaction().reaction.extent_mol, expected[3], 1e-9);
            assert_eq!(&next.report().gas, &next.reaction().heat);
            let gas_heat = capacity() * (expected[2] - 300.0) - release * expected[3];
            near(next.reaction().heat.supplied_heat_j, gas_heat, 2e-7);
            near(next.report().solid_to_gas_heat_j, gas_heat, 2e-7);
            assert!(
                next.report().combined_energy_residual_j.abs() <= solve().total_energy_tolerance_j
            );
            let loads = model
                .geometry()
                .assemble(&[], VolumeLoad3::None, &[], cx)
                .unwrap();
            let free = next.state().solve(&loads, &FREE, config(), cx).unwrap();
            let strain = [expected[0], expected[1]].map(|t| 1e-5 * (t - 300.0));
            let displacement = free.solution.displacement.displacement_m();
            near(displacement[1][0], 0.5 * (strain[0] + strain[1]), 1e-10);
            near(displacement[3][2], strain[0], 1e-10);
            near(displacement[4][2], -strain[1], 1e-10);
            let fixed: Vec<_> = (0..15).collect();
            let clamped = next.state().solve(&loads, &fixed, config(), cx).unwrap();
            for (element, expansion) in clamped.solution.stress.elements().iter().zip(strain) {
                for stress in &element.stress_mandel_pa()[..3] {
                    near(*stress, -1000.0 * expansion, 1e-8);
                }
            }
            extents.push(expected[3]);
        }
        assert!(
            extents[0] > extents[1] && extents[1] > extents[2],
            "cooling must enter the chemical equation inside this accepted interval"
        );
        assert!(extents[0] - extents[2] > 1e-4);
    });
}

#[test]
fn continuation_reuses_endpoint_species_rates_and_one_spatial_material_network() {
    with_cx(|cx, _| {
        let (rate, release, product_k) = (1.5, 20000.0, 4.0);
        let law = implicit(reaction(rate, release, product_k));
        let original_gas = gas(&law);
        let mut gas = original_gas.clone();
        let (card, caloric) = wall(1.0, true);
        let mut thermal = initial(caloric.curve());
        let initial_solid_energy = thermal.energy_j();
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = [30, 10].map(|id| MaterialCell {
            cell: CellId(id),
            material_id: MaterialId(0),
            curve: caloric.curve(),
            point: ThermalMaterialPoint {
                material: &card,
                context: &context,
                temperature_axis: "T",
                liquid_fraction_axis: None,
            },
            conductivity: Some((&key, MaterialPropertySelection::SingleClaimOnly)),
        });
        let links = [ConductionLink {
            link: 7,
            distance_a_m: 0.03,
            distance_b_m: 0.07,
        }];
        let network = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &links,
        };
        let transport = network.prepare(&thermal, || false).unwrap();
        let mut product = 0.0;
        let mut previous_rate_card = None;
        for interval in 1..=2 {
            let initial_temperature = gas.state().state().temperature;
            let next = model
                .step_with_coupled_reacting_isochoric_gas(
                    &thermal,
                    &gas,
                    &law,
                    interval as f64,
                    0.0,
                    &contacts(),
                    Some(&transport),
                    chemical(),
                    solve(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            let t = next.gas().state().state().temperature;
            let xi = extent(t, product, rate, release, 1.0);
            let chemistry = &next.reaction().reaction;
            near(chemistry.extent_mol, xi, 1e-9);
            near(chemistry.species[1].initial_mol, product, 1e-11);
            near(chemistry.species[1].final_mol, product + xi, 1e-9);
            let rate_source = chemistry.rate_material.as_ref().unwrap();
            near(
                rate_source
                    .property(ACTIVITY_RATE_SCALE_PROPERTY)
                    .unwrap()
                    .value_si(),
                rate * t / 300.0,
                1e-10,
            );
            assert!((t - initial_temperature).abs() > 1.0);
            if let Some(identity) = previous_rate_card {
                assert_eq!(rate_source.card_identity(), identity);
            }
            previous_rate_card = Some(rate_source.card_identity());
            near(
                next.gas().state().state().thermal_conductivity,
                0.03 * (1.0 + (product_k - 1.0) * (product + xi) / amount()),
                1e-10,
            );
            let temperatures =
                [CellId(30), CellId(10)].map(|id| next.thermal().temperature_k(id).unwrap());
            let sample = next.report().transport.as_ref().unwrap();
            near(
                sample.links[0].conductivity_w_m_k[0],
                temperatures[0] / 300.0,
                1e-10,
            );
            near(
                sample.links[0].conductivity_w_m_k[1],
                temperatures[1] / 300.0,
                1e-10,
            );
            near(
                sample.links[0].conductance_w_k,
                0.5 / (0.03 / (temperatures[0] / 300.0) + 0.07 / (temperatures[1] / 300.0)),
                1e-10,
            );
            for patch in &next.report().contacts {
                let internal = sample
                    .cells
                    .iter()
                    .find(|cell| cell.cell == patch.cell)
                    .unwrap()
                    .conductivity
                    .as_ref()
                    .unwrap();
                assert_eq!(internal.query(), patch.source.query());
                assert_eq!(internal.card_identity(), patch.source.card_identity());
            }
            near(
                next.thermal().energy_j() - initial_solid_energy
                    + gas.mass_kg()
                        * (next.gas().state().specific_internal_energy_j_kg()
                            - original_gas.state().specific_internal_energy_j_kg()),
                0.0,
                4e-7,
            );
            near(
                next.reaction().heat.supplied_heat_j,
                capacity() * (t - initial_temperature) - release * xi,
                2e-7,
            );
            assert_eq!(next.gas().point().axes().get("catalyst-lot"), Some(&7.0));
            product = chemistry.species[1].final_mol;
            let (coupled, _) = next.into_parts();
            let (accepted_thermal, accepted_gas, _, _) = coupled.into_parts();
            thermal = accepted_thermal;
            gas = accepted_gas;
        }
        assert_eq!(thermal.time_s(), 2.0);
        assert!(product > 0.0);
        assert_eq!(gas.mass_kg(), original_gas.mass_kg());
        assert_eq!(gas.volume_m3(), original_gas.volume_m3());
    });
}

#[test]
fn coupled_late_refusals_and_cancellation_preserve_inputs_and_deterministic_retry() {
    with_cx(|cx, gate| {
        let law = implicit(reaction(1.0, 20000.0, 4.0));
        let gas = gas(&law);
        let saved_gas = gas.clone();
        let closed = ClosedGasParcel::from_isochoric(gas.clone());
        let mut checkpoints = 0;
        let cancelled = law.heated_implicit_candidate(&closed, 1.0, -100.0, chemical(), || {
            checkpoints += 1;
            checkpoints == 12
        });
        assert!(matches!(cancelled, Err(GasReactionError::Cancelled)));
        assert!(checkpoints >= 12);
        assert_eq!(closed, ClosedGasParcel::from_isochoric(saved_gas.clone()));
        let (card, caloric) = wall(1.0, false);
        let thermal = initial(caloric.curve());
        let saved_cells = thermal.cells().clone();
        let context = point(300.0);
        let model = solid(&thermal, &card, caloric.curve(), &context, cx);
        let patches = contacts();
        let attempt = |controls: SolidGasSolve| {
            model.step_with_coupled_reacting_isochoric_gas(
                &thermal,
                &gas,
                &law,
                1.0,
                0.0,
                &patches,
                None,
                chemical(),
                controls,
                &[0.0; 3],
                cx,
            )
        };
        let chemical_trial = law
            .heated_implicit_candidate(&closed, 1.0, 0.0, chemical(), || false)
            .unwrap();
        assert!(chemical_trial.step().reaction.extent_mol > 0.0);
        assert!(matches!(
            attempt(SolidGasSolve {
                maximum_trials: 1,
                ..solve()
            }),
            Err(SolidGasError::TrialBudget)
        ));
        assert!(
            attempt(SolidGasSolve {
                maximum_gas_temperature_increment_k: 1.0,
                ..solve()
            })
            .is_err()
        );
        let mismatched = SolidGasSolve {
            gas: GasEnergySolve {
                energy_tolerance_j_kg: 1e-8,
                ..solve().gas
            },
            ..solve()
        };
        assert!(matches!(
            attempt(mismatched),
            Err(SolidGasError::Invalid(_))
        ));
        let accepted = attempt(solve()).unwrap();
        let replay = attempt(solve()).unwrap();
        assert_eq!(accepted.gas(), replay.gas());
        assert_eq!(accepted.thermal().cells(), replay.thermal().cells());
        assert_eq!(accepted.reaction(), replay.reaction());
        gate.request();
        assert!(matches!(attempt(solve()), Err(SolidGasError::Cancelled)));
        assert_eq!(thermal.cells(), &saved_cells);
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas, saved_gas);
    });
}
