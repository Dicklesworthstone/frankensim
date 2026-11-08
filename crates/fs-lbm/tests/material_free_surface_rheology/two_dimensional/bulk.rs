//! Physical planar bulk feedback retains the original selected material bundle.
use super::super::*;
use super::{E, W, same_native};
use fs_lbm::core2::CollisionModel2;
use fs_lbm::freesurface::{ContactModel, FreeSurface};
use fs_lbm::material::free_surface::rheology::MaterialFreeSurfaceRheology2;
use fs_lbm::{Cell, Grid};
use fs_material::fluid_transport::viscous::DeformationRate;

fn bulk_card(maximum_rate: f64) -> MaterialCard {
    let original = card(0.01, 100.0, 0.001);
    let mut claims = original.claims().clone();
    let dims = Dims([-1, 1, -1, 0, 0, 0]);
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new("bulk_eta", dims),
            value: PropertyValue::Curve {
                abscissa: "rate".into(),
                abscissa_dims: Frequency::DIMS,
                knots: vec![(0.0, 0.2), (maximum_rate, bulk_viscosity(maximum_rate))],
                dims,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 300.0, 500.0)
                .with("pressure", 100000.0, 100000.0)
                .with_quantity(
                    "rate",
                    QuantitySpec::dimensional(Frequency::DIMS),
                    0.0,
                    maximum_rate,
                ),
            interpolation: InterpolationPolicy::LinearInside,
            uncertainty: UncertaintyModel::Unstated,
            observations: vec![],
            provenance: Provenance {
                source: "synthetic original planar volumetric viscosity".into(),
                license: "test-only".into(),
                artifact: None,
            },
        })
        .unwrap();
    MaterialCard::supersede(&original, claims, vec![]).unwrap()
}

fn bulk_law(source: &MaterialCard) -> IsothermalRheologyLaw<'_> {
    let original = bundle(source);
    let mut required: Vec<_> = original
        .properties()
        .iter()
        .map(|property| property.requirement().clone())
        .collect();
    required.push(
        ScalarPropertyRequirement::try_new(
            "bulk_eta",
            Dims([-1, 1, -1, 0, 0, 0]),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap(),
    );
    let material = resolve_material_state_point(
        source,
        &original.query(),
        &required,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    IsothermalRheologyLaw::try_new(
        source,
        material,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Sourced,
        "T",
        "rate",
    )
    .unwrap()
}

fn bulk_viscosity(rate: f64) -> f64 {
    0.2 + 0.04 * rate
}
fn shear_viscosity(rate: f64) -> f64 {
    0.1 + 0.01 * rate
}
fn deformation(cell: usize) -> [[f64; 2]; 2] {
    match cell % 5 {
        0 => [[0.0; 2]; 2],
        1 => [[0.0, 0.4], [0.4, 0.0]],
        2 => [[0.6, 0.0], [0.0, 0.6]],
        3 => [[-0.5, 0.0], [0.0, -0.5]],
        _ => [[1.1, 0.2], [0.2, -0.4]],
    }
}
fn rate(d: [[f64; 2]; 2]) -> f64 {
    DeformationRate::from_tensor([d[0][0], d[1][1], 0.0, d[0][1], 0.0, 0.0])
        .unwrap()
        .shear_rate_per_s()
}
fn shear_tau(mu: f64) -> f64 {
    0.5 + 3.0 * (mu / 1000.0) * DT / (DX * DX)
}
fn trace_tau(mu: f64, zeta: f64) -> f64 {
    // Independent dimensional stress reduction, rather than numerical tau weights.
    0.5 + 3.0 * ((zeta + mu / 3.0) / 1000.0) * DT / (DX * DX)
}
fn seed(grid: &mut Grid, cell: usize, d: [[f64; 2]; 2], rho: f64, force: [f64; 2], zeta: f64) {
    let gamma = rate(d);
    let mu = shear_viscosity(gamma);
    let tau_s = shear_tau(mu);
    let tau_t = trace_tau(mu, zeta);
    let mean = 0.5 * (d[0][0] + d[1][1]);
    let u = [0.01, -0.015];
    let equilibrium = fs_lbm::equilibrium(rho, u[0], u[1]);
    grid.f[cell] = core::array::from_fn(|q| {
        let e = [f64::from(E[q].0), f64::from(E[q].1)];
        let first = -0.5 * (e[0] * force[0] + e[1] * force[1]);
        let mut contraction = 0.0;
        for a in 0..2 {
            for b in 0..2 {
                let trace = if a == b { mean } else { 0.0 };
                let pi = -2.0 * rho * DT / 3.0 * (tau_s * (d[a][b] - trace) + tau_t * trace)
                    - 0.5 * (u[a] * force[b] + force[a] * u[b]);
                contraction += (e[a] * e[b] - if a == b { 1.0 / 3.0 } else { 0.0 }) * pi;
            }
        }
        equilibrium[q] + W[q] * (3.0 * first + 4.5 * contraction)
    });
    grid.fext[cell] = force;
}

#[test]
fn original_bulk_and_shear_sources_share_one_planar_rate_and_native_collision() {
    let source = bulk_card(100.0);
    let binding = bulk_law(&source);
    let gate = CancelGate::new_clock_free();
    for rho in [0.8, 1.2] {
        for force in [[0.0; 2], [0.0003, -0.0002]] {
            for selection in 0..4 {
                let mut adapter =
                    MaterialFreeSurfaceRheology2::try_new(&binding, [5, 4], DX, DT, 0.05, solve())
                        .unwrap()
                        .with_bulk_viscosity_source("bulk_eta")
                        .unwrap();
                let mut grid = Grid::uniform(5, 4, 1.6);
                grid.collision = match selection {
                    0 => CollisionModel2::Bgk,
                    1 => CollisionModel2::CentralMoment,
                    _ => CollisionModel2::CentralMomentBulk { bulk_tau: 1.7 },
                };
                if selection == 3 {
                    grid.set_central_moment_bulk_field(
                        &(0..20)
                            .map(|cell| [0.65, 1.9][cell % 2])
                            .collect::<Vec<_>>(),
                    )
                    .unwrap();
                }
                let mut expected_trace = vec![0.0; 20];
                let mut expected_shear = vec![0.0; 20];
                for cell in 0..20 {
                    let d = deformation(cell);
                    let gamma = rate(d);
                    let mu = shear_viscosity(gamma);
                    let zeta = bulk_viscosity(gamma);
                    seed(&mut grid, cell, d, rho, force, zeta);
                    expected_shear[cell] = shear_tau(mu);
                    expected_trace[cell] = trace_tau(mu, zeta);
                }
                let mut actual = FreeSurface::new(
                    grid,
                    adapter.scaling().surface_tension_lattice(),
                    ContactModel::Neutral,
                );
                let mut direct = actual.clone();
                direct.grid.tau.copy_from_slice(&expected_shear);
                direct
                    .grid
                    .set_central_moment_bulk_field(&expected_trace)
                    .unwrap();
                let expected = direct.step_with_momentum(&gate).unwrap();
                let report = adapter.step(&mut actual, &gate).unwrap();
                assert_eq!(report, expected);
                assert_eq!(actual.grid.flags, direct.grid.flags);
                for cell in 0..20 {
                    let position = [cell % 5, cell / 5];
                    let expected_rate = rate(deformation(cell));
                    let sample = adapter.rheology_at(position).unwrap();
                    let material = adapter.material_at(position).unwrap();
                    near(sample.shear_rate_per_s, expected_rate, 3e-9);
                    near(sample.collision_tau, expected_shear[cell], 3e-11);
                    assert_eq!(actual.grid.tau[cell], sample.collision_tau);
                    near(
                        actual.grid.central_moment_bulk_tau_at(cell).unwrap(),
                        expected_trace[cell],
                        3e-11,
                    );
                    near(
                        material.property("dynamic_viscosity").unwrap().value_si(),
                        shear_viscosity(expected_rate),
                        3e-11,
                    );
                    near(
                        material.property("bulk_eta").unwrap().value_si(),
                        bulk_viscosity(expected_rate),
                        3e-11,
                    );
                    assert_eq!(material.card_identity(), source.content_hash());
                    assert!(
                        material
                            .query_point()
                            .contains(&("rate".into(), sample.shear_rate_per_s),)
                    );
                    for name in ["dynamic_viscosity", "bulk_eta"] {
                        assert_eq!(
                            material.property(name).unwrap().answer().receipt.selected,
                            binding
                                .material()
                                .property(name)
                                .unwrap()
                                .answer()
                                .receipt
                                .selected,
                        );
                    }
                    for q in 0..9 {
                        near(actual.grid.f[cell][q], direct.grid.f[cell][q], 2e-12);
                    }
                    near(actual.mass[cell], direct.mass[cell], 2e-12);
                }
                near(actual.ledger_mass(), direct.ledger_mass(), 2e-11);
            }
        }
    }
}

#[test]
fn physical_bulk_keeps_the_original_interface_fill_and_dry_cell_transaction() {
    let source = bulk_card(100.0);
    let binding = bulk_law(&source);
    let mut adapter =
        MaterialFreeSurfaceRheology2::try_new(&binding, [6, 6], DX, DT, 0.05, solve())
            .unwrap()
            .with_bulk_viscosity_source("bulk_eta")
            .unwrap();
    let mut grid = Grid::uniform(6, 6, 1.6);
    grid.periodic_x = false;
    grid.periodic_y = false;
    for y in 0..6 {
        for x in 0..6 {
            let cell = grid.idx(x, y);
            grid.flags[cell] = if x == 0 || x == 5 || y == 0 || y == 5 {
                Cell::Wall
            } else if y <= 2 {
                Cell::Fluid
            } else {
                Cell::Gas
            };
        }
    }
    let mut actual = FreeSurface::new(
        grid,
        adapter.scaling().surface_tension_lattice(),
        ContactModel::Neutral,
    );
    let incoming: Vec<_> = actual
        .grid
        .flags
        .iter()
        .map(|flag| matches!(flag, Cell::Fluid | Cell::Interface))
        .collect();
    assert!(actual.grid.flags.contains(&Cell::Interface));
    assert!(actual.grid.flags.contains(&Cell::Gas));
    let mut traces = vec![trace_tau(0.1, 0.2); 36];
    for cell in 0..36 {
        if incoming[cell] {
            let d = deformation(cell);
            let gamma = rate(d);
            seed(
                &mut actual.grid,
                cell,
                d,
                1.0,
                [0.0001, -0.00005],
                bulk_viscosity(gamma),
            );
            actual.mass[cell] = if actual.grid.flags[cell] == Cell::Fluid {
                1.0
            } else {
                0.4
            };
            traces[cell] = trace_tau(shear_viscosity(gamma), bulk_viscosity(gamma));
        }
    }
    let mut direct = actual.clone();
    for cell in 0..36 {
        if incoming[cell] {
            direct.grid.tau[cell] = shear_tau(shear_viscosity(rate(deformation(cell))));
        }
    }
    direct.grid.set_central_moment_bulk_field(&traces).unwrap();
    let gate = CancelGate::new_clock_free();
    let expected = direct.step_with_momentum(&gate).unwrap();
    let report = adapter.step(&mut actual, &gate).unwrap();
    near(
        report.absolute_mass_defect,
        expected.absolute_mass_defect,
        2e-11,
    );
    for axis in 0..2 {
        near(
            report.absolute_momentum_defect[axis],
            expected.absolute_momentum_defect[axis],
            2e-11,
        );
    }
    near(
        report.kinetic_energy_change,
        expected.kinetic_energy_change,
        2e-11,
    );
    assert_eq!(actual.grid.flags, direct.grid.flags);
    for cell in 0..36 {
        assert_eq!(
            adapter.material_at([cell % 6, cell / 6]).is_some(),
            incoming[cell]
        );
        near(
            actual.grid.central_moment_bulk_tau_at(cell).unwrap(),
            traces[cell],
            3e-11,
        );
        near(actual.mass[cell], direct.mass[cell], 3e-12);
        for q in 0..9 {
            near(actual.grid.f[cell][q], direct.grid.f[cell][q], 3e-12);
        }
    }
    near(actual.ledger_mass(), direct.ledger_mass(), 3e-11);
}

#[test]
fn constant_physical_bulk_keeps_local_shear_dependence_and_fresh_source_epochs() {
    let source = card(0.01, 100.0, 0.001);
    let binding = law(&source);
    let zeta = 0.25;
    let mut adapter =
        MaterialFreeSurfaceRheology2::try_new(&binding, [5, 4], DX, DT, 0.05, solve())
            .unwrap()
            .with_bulk_viscosity(zeta)
            .unwrap();
    let mut grid = Grid::uniform(5, 4, 1.6);
    for cell in 0..20 {
        seed(&mut grid, cell, deformation(cell), 1.0, [0.0; 2], zeta);
    }
    let mut liquid = FreeSurface::new(
        grid,
        adapter.scaling().surface_tension_lattice(),
        ContactModel::Neutral,
    );
    let gate = CancelGate::new_clock_free();
    let mut previous_rates = None;
    for _ in 0..3 {
        adapter.step(&mut liquid, &gate).unwrap();
        let mut rates = Vec::new();
        for cell in 0..20 {
            let position = [cell % 5, cell / 5];
            let sample = adapter.rheology_at(position).unwrap();
            let material = adapter.material_at(position).unwrap();
            let mu = material.property("dynamic_viscosity").unwrap().value_si();
            near(mu, shear_viscosity(sample.shear_rate_per_s), 2e-12);
            near(
                liquid.grid.central_moment_bulk_tau_at(cell).unwrap(),
                trace_tau(mu, zeta),
                2e-12,
            );
            assert_eq!(liquid.grid.tau[cell], sample.collision_tau);
            rates.push(sample.shear_rate_per_s.to_bits());
        }
        if let Some(previous) = previous_rates {
            assert_ne!(
                rates, previous,
                "later collision must query its own evolved stress"
            );
        }
        previous_rates = Some(rates);
    }
    assert_eq!(adapter.scaling().bulk_viscosity_pa_s(), Some(zeta));
    for invalid in [0.0, -0.1, f64::NAN, f64::INFINITY] {
        assert!(adapter.clone().with_bulk_viscosity(invalid).is_err());
    }
    assert!(
        adapter
            .clone()
            .with_bulk_viscosity_source("missing")
            .is_err()
    );
    assert!(
        adapter
            .clone()
            .with_bulk_viscosity_source("density")
            .is_err()
    );
}

#[test]
fn physical_bulk_refusals_preserve_native_fields_and_accepted_source_epochs() {
    let source = bulk_card(1.0);
    let binding = bulk_law(&source);
    let mut adapter =
        MaterialFreeSurfaceRheology2::try_new(&binding, [5, 4], DX, DT, 0.05, solve())
            .unwrap()
            .with_bulk_viscosity_source("bulk_eta")
            .unwrap();
    let mut liquid = FreeSurface::new(
        Grid::uniform(5, 4, 0.8),
        adapter.scaling().surface_tension_lattice(),
        ContactModel::Neutral,
    );
    let gate = CancelGate::new_clock_free();
    adapter.step(&mut liquid, &gate).unwrap();
    let accepted = liquid.clone();
    let receipt = adapter.material_at([0, 0]).unwrap().clone();
    let sample = *adapter.rheology_at([0, 0]).unwrap();

    let d = [[2.0, 0.0], [0.0, 2.0]];
    seed(
        &mut liquid.grid,
        19,
        d,
        1.0,
        [0.0; 2],
        bulk_viscosity(rate(d)),
    );
    let before = liquid.clone();
    // Shear viscosity supports rates through 100/s; only the separately
    // retained bulk claim leaves its original support at this last cell.
    assert!(matches!(
        adapter.step(&mut liquid, &gate),
        Err(FreeSurfaceRheologyError::Material(_))
    ));
    same_native(&liquid, &before);
    assert_eq!(adapter.material_at([0, 0]).unwrap(), &receipt);
    assert_eq!(*adapter.rheology_at([0, 0]).unwrap(), sample);

    // A valid first source proposal changes both coefficients, then the
    // original endpoint scenario-speed gate rejects the native candidate.
    let mut late = accepted.clone();
    let d = [[0.0, 0.25], [0.25, 0.0]];
    seed(&mut late.grid, 0, d, 1.0, [0.0; 2], bulk_viscosity(rate(d)));
    late.rebind_gravity([0.07, 0.0]).unwrap();
    let saved = late.clone();
    assert!(matches!(
        adapter.step(&mut late, &gate),
        Err(FreeSurfaceRheologyError::Invalid(_))
    ));
    same_native(&late, &saved);
    assert_eq!(adapter.material_at([0, 0]).unwrap(), &receipt);
    assert_eq!(*adapter.rheology_at([0, 0]).unwrap(), sample);
    let stopped = CancelGate::new_clock_free();
    stopped.request();
    assert!(adapter.step(&mut late, &stopped).is_err());
    same_native(&late, &saved);

    let unselected = law(&source);
    assert!(
        MaterialFreeSurfaceRheology2::try_new(&unselected, [5, 4], DX, DT, 0.05, solve(),)
            .unwrap()
            .with_bulk_viscosity_source("bulk_eta")
            .is_err()
    );
}

#[test]
fn coupled_planar_root_keeps_its_original_evaluation_budget() {
    let source = bulk_card(100.0);
    let binding = bulk_law(&source);
    let mut adapter = MaterialFreeSurfaceRheology2::try_new(
        &binding,
        [5, 4],
        DX,
        DT,
        0.05,
        RheologySolve {
            maximum_iterations: 1,
            ..solve()
        },
    )
    .unwrap()
    .with_bulk_viscosity_source("bulk_eta")
    .unwrap();
    let mut grid = Grid::uniform(5, 4, 1.6);
    let d = [[0.6, 0.0], [0.0, 0.6]];
    seed(&mut grid, 0, d, 1.0, [0.0; 2], bulk_viscosity(rate(d)));
    let mut liquid = FreeSurface::new(
        grid,
        adapter.scaling().surface_tension_lattice(),
        ContactModel::Neutral,
    );
    let before = liquid.clone();
    assert!(matches!(
        adapter.step(&mut liquid, &CancelGate::new_clock_free()),
        Err(FreeSurfaceRheologyError::Convergence { iterations: 1, .. }),
    ));
    same_native(&liquid, &before);
    assert!(adapter.material_at([0, 0]).is_none());
}
