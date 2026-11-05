//! Physical local bulk transport through the original compressible NSF owner.
use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::dispersion::{C64, GasPoint, linearized_modes, longitudinal_eigenvector};
use fs_lbm::nsf3::{
    NsfBulkViscosityPoint, NsfBulkViscositySource, NsfCharacteristicAxis,
    NsfCharacteristicBoundary, NsfError, NsfGrid3, NsfSpec, WallThermal,
};
use fs_material::gas::GasSpec;
use std::f64::consts::TAU;

fn spec(dims: [usize; 3], dx_m: f64) -> NsfSpec {
    NsfSpec {
        dims,
        periodic: [true; 3],
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0,
        pressure_pa: 100_000.0,
        dx_m,
        temperature_window_k: [270.0, 330.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.15,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "native ratio control; manufactured bulk-law regression",
        maximum_cells: dims.iter().product(),
        wall: WallThermal::Adiabatic,
        workers: 1,
        transport_points: 33,
    }
}

fn near(actual: f64, expected: f64, relative: f64, absolute: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * expected.abs(),
        "actual={actual:e} expected={expected:e}"
    );
}

fn variable(point: NsfBulkViscosityPoint) -> f64 {
    // Manufactured positive T/p law, deliberately not proportional to mu.
    4.0e-5 * (1.0 + (point.temperature_k - 300.0) / 120.0) * (point.pressure_pa / 100_000.0)
}

fn source() -> NsfBulkViscositySource {
    NsfBulkViscositySource::new("manufactured independent T/p bulk law", |point| {
        Ok(variable(point))
    })
    .unwrap()
}

fn varying_grid() -> NsfGrid3<D3q27> {
    NsfGrid3::new(
        spec([12, 2, 2], 2e-5),
        |_| false,
        |x| {
            let a = TAU * x[0] as f64 / 12.0;
            (
                1.0 + 0.01 * a.cos(),
                [0.003 * a.sin(), 0.001 * a.cos(), 0.0],
                300.0 + 2.0 * a.sin(),
            )
        },
    )
    .unwrap()
}

#[test]
fn g1_bulk_source_uses_actual_local_temperature_pressure_and_endpoint() {
    let mut grid = varying_grid().with_bulk_viscosity(source()).unwrap();
    assert_eq!(
        grid.bulk_viscosity_source(),
        Some("manufactured independent T/p bulk law")
    );
    let energy = grid.total_energy_j();
    let reference = GasPoint::from_law(grid.law());
    near(reference.mu_bulk, 4e-5, 5e-15, 0.0);
    let initial = grid.bulk_viscosity_at([2, 0, 0]).unwrap().unwrap();
    for step in 0..9 {
        for x in 0..12 {
            let sample = grid.bulk_viscosity_at([x, 0, 0]).unwrap().unwrap();
            let (rho, _, temperature) = grid.cell([x, 0, 0]);
            assert_eq!(sample.point.temperature_k, temperature);
            near(
                sample.point.density_kg_m3,
                rho * grid.law().rho_ref_kg_m3,
                0.0,
                0.0,
            );
            near(
                sample.point.pressure_pa,
                sample.point.density_kg_m3 * grid.law().gas_constant * temperature,
                0.0,
                0.0,
            );
            near(
                sample.bulk_viscosity_pa_s,
                variable(sample.point),
                2e-15,
                0.0,
            );
            let tau_bulk = 0.5
                + 4.5 * (sample.bulk_viscosity_pa_s / sample.point.density_kg_m3) * grid.law().dt_s
                    / grid.law().dx_m.powi(2);
            near(sample.trace_relaxation_rate, 1.0 / tau_bulk, 3e-15, 0.0);
        }
        if step != 8 {
            let ledger = grid.step().unwrap();
            assert!(ledger.residual_rel < 3e-14);
        }
    }
    assert_ne!(initial, grid.bulk_viscosity_at([2, 0, 0]).unwrap().unwrap());
    near(grid.total_energy_j(), energy, 3e-14, 0.0);
}

#[test]
fn g3_a_matching_molecular_ratio_recovers_the_original_trajectory() {
    let mut original = varying_grid();
    let mut selected = varying_grid()
        .with_bulk_viscosity(
            NsfBulkViscositySource::new("same native molecular ratio", |p| {
                Ok(0.73 * p.dynamic_viscosity_pa_s)
            })
            .unwrap(),
        )
        .unwrap();
    assert!(original.bulk_viscosity_source().is_none());
    assert!(original.bulk_viscosity_at([0, 0, 0]).unwrap().is_none());
    for _ in 0..20 {
        let a = original.step().unwrap();
        let b = selected.step().unwrap();
        near(a.internal_j, b.internal_j, 5e-14, 0.0);
        near(a.kinetic_j, b.kinetic_j, 5e-12, 1e-22);
        for x in 0..12 {
            let (ra, ua, ta) = original.cell([x, 0, 0]);
            let (rb, ub, tb) = selected.cell([x, 0, 0]);
            near(ra, rb, 5e-14, 0.0);
            near(ta, tb, 5e-14, 0.0);
            for i in 0..3 {
                near(ua[i], ub[i], 2e-12, 5e-15);
            }
        }
    }
}

#[test]
fn g3_shear_sgs_does_not_become_molecular_bulk_viscosity() {
    use fs_lbm::nsf3::NsfSgsConfig;
    use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};
    let mut grid = varying_grid()
        .with_bulk_viscosity(
            NsfBulkViscositySource::new("manufactured twice-molecular-mu bulk law", |p| {
                Ok(2.0 * p.dynamic_viscosity_pa_s)
            })
            .unwrap(),
        )
        .unwrap()
        .with_sgs(NsfSgsConfig {
            closure: SgsConfig {
                model: SgsModel::Smagorinsky { cs: 0.17 },
                filter_width_lu: 1.0,
                maximum_eddy_to_molecular: 20.0,
                usage: SgsUse::Exploration,
            },
            turbulent_prandtl: 0.8,
        })
        .unwrap();
    for _ in 0..3 {
        let mut nonzero = false;
        for x in 0..12 {
            let bulk = grid.bulk_viscosity_at([x, 0, 0]).unwrap().unwrap();
            let sgs = grid.sgs_transport([x, 0, 0]).unwrap().unwrap();
            nonzero |= sgs.eddy_viscosity_lu > 0.0;
            let theta = grid.law().theta_of_kelvin(bulk.point.temperature_k);
            let mu =
                grid.law().mu_lattice(theta) * grid.law().rho_ref_kg_m3 * grid.law().dx_m.powi(2)
                    / grid.law().dt_s;
            near(bulk.bulk_viscosity_pa_s, 2.0 * mu, 3e-15, 0.0);
        }
        assert!(nonzero);
        grid.step().unwrap();
    }
}

fn constant(value: f64) -> NsfBulkViscositySource {
    NsfBulkViscositySource::new(
        "manufactured constant molecular bulk viscosity",
        move |_| Ok(value),
    )
    .unwrap()
}

fn acoustic_decay(n: usize, bulk: f64) -> (f64, f64) {
    let length = 0.0032;
    let setup = spec([n, 1, 1], length / n as f64);
    let law = NsfGrid3::<D3q27>::new(setup.clone(), |_| false, |_| (1.0, [0.0; 3], 300.0))
        .unwrap()
        .with_bulk_viscosity(constant(bulk))
        .unwrap()
        .law()
        .clone();
    let point = GasPoint::from_law(&law);
    let mode = linearized_modes(&point, TAU / length).unwrap()[0];
    let vector = longitudinal_eigenvector(
        &point,
        TAU / length,
        C64 {
            re: mode.re,
            im: -mode.im,
        },
    )
    .unwrap();
    let speed = (law.gamma / 3.0).sqrt();
    let mut grid = NsfGrid3::<D3q27>::new(
        setup,
        |_| false,
        |x| {
            let phase = TAU * x[0] as f64 / n as f64;
            let real = |v: C64| v.re * phase.cos() - v.im * phase.sin();
            (
                1.0 + 1e-6 * real(vector[0]),
                [1e-6 * speed * real(vector[1]), 0.0, 0.0],
                300.0 * (1.0 + 1e-6 * real(vector[2])),
            )
        },
    )
    .unwrap()
    .with_bulk_viscosity(constant(bulk))
    .unwrap();
    let period = TAU / (mode.im * law.dt_s);
    let mut sums = [0.0; 5];
    for step in 0..(10.0 * period).ceil() as usize {
        if step as f64 >= 2.0 * period && step % 8 == 0 {
            let mut amplitude = [0.0_f64; 2];
            for x in 0..n {
                let rho = grid.cell([x, 0, 0]).0 - 1.0;
                let phase = TAU * x as f64 / n as f64;
                amplitude[0] += rho * phase.cos();
                amplitude[1] -= rho * phase.sin();
            }
            let t = step as f64 * law.dt_s;
            let y = amplitude[0].hypot(amplitude[1]).ln();
            sums[0] += 1.0;
            sums[1] += t;
            sums[2] += y;
            sums[3] += t * y;
            sums[4] += t * t;
        }
        let ledger = grid.step().unwrap();
        assert!(ledger.residual_rel < 3e-14);
    }
    let [count, t, y, ty, tt] = sums;
    (-(count * ty - t * y) / (count * tt - t * t), -mode.re)
}

#[test]
fn g1_added_bulk_damps_longitudinal_sound_without_adding_shear_loss() {
    let mut errors = Vec::new();
    for n in [64, 128] {
        let low = acoustic_decay(n, 1.4e-5);
        let high = acoustic_decay(n, 3.5e-4);
        let error = ((high.0 - low.0) / (high.1 - low.1) - 1.0).abs();
        println!(
            "bulk acoustic increment n={n}: measured={} exact={} relative_error={error}",
            high.0 - low.0,
            high.1 - low.1
        );
        errors.push(error);
    }
    assert!(
        errors[0] < 1.0e-3 && errors[1] < 2.0e-4 && errors[1] < errors[0],
        "bulk attenuation refinement: {errors:?}"
    );

    let make = |bulk| {
        NsfGrid3::<D3q27>::new(
            spec([32, 1, 1], 5e-5),
            |_| false,
            |x| {
                (
                    1.0,
                    [0.0, 1e-6 * (TAU * x[0] as f64 / 32.0).sin(), 0.0],
                    300.0,
                )
            },
        )
        .unwrap()
        .with_bulk_viscosity(constant(bulk))
        .unwrap()
    };
    let mut low = make(1.4e-5);
    let mut high = make(3.5e-4);
    for _ in 0..500 {
        low.step().unwrap();
        high.step().unwrap();
    }
    for x in 0..32 {
        near(
            low.cell([x, 0, 0]).1[1],
            high.cell([x, 0, 0]).1[1],
            1e-7,
            2e-15,
        );
    }
}

#[test]
fn g4_source_endpoint_refusal_keeps_the_original_energy_and_material_epoch() {
    let mut grid = NsfGrid3::<D3q27>::new_with_sources(
        spec([3, 2, 2], 5e-5),
        |_| false,
        |_| (1.0, [0.0; 3], 300.0),
        Some(Box::new(|_, _| (0.0, [0.0; 3], 0.001))),
    )
    .unwrap()
    .with_bulk_viscosity(
        NsfBulkViscositySource::new("manufactured finite temperature support", |point| {
            if point.temperature_k > 300.01 {
                Err(NsfError::Input("bulk source temperature domain"))
            } else {
                Ok(2e-5)
            }
        })
        .unwrap(),
    )
    .unwrap();
    let before = (
        grid.state_hash(),
        grid.steps(),
        grid.cumulative_ledger(),
        grid.bulk_viscosity_at([0, 0, 0]).unwrap(),
        grid.total_energy_j(),
    );
    assert!(grid.step().is_err());
    assert_eq!(
        before,
        (
            grid.state_hash(),
            grid.steps(),
            grid.cumulative_ledger(),
            grid.bulk_viscosity_at([0, 0, 0]).unwrap(),
            grid.total_energy_j()
        )
    );
    assert!(NsfBulkViscositySource::new("   ", |_| Ok(1e-5)).is_err());
    for bulk in [
        0.0,
        -1e-5,
        f64::NAN,
        f64::INFINITY,
        f64::MIN_POSITIVE,
        f64::MAX,
    ] {
        assert!(
            varying_grid().with_bulk_viscosity(constant(bulk)).is_err(),
            "bulk={bulk}"
        );
    }
}

#[test]
fn g3_characteristic_trials_retain_bulk_source_and_publish_it_once() {
    let make = || {
        let mut setup = spec([12, 2, 2], 5e-5);
        setup.periodic[0] = false;
        NsfGrid3::<D3q27>::new(setup, |_| false, |_| (1.0, [0.0; 3], 300.0))
            .unwrap()
            .with_bulk_viscosity(source())
            .unwrap()
            .with_characteristic_axis(NsfCharacteristicAxis {
                axis: 0,
                faces: [NsfCharacteristicBoundary {
                    reference_pressure_pa: 100_000.0,
                    reference_temperature_k: 300.0,
                    reference_velocity_m_s: [0.0; 3],
                    incoming_characteristic_m_s: 0.0,
                }; 2],
                maximum_relative_perturbation: 0.05,
            })
            .unwrap()
    };
    let mut direct = make();
    let mut preview = make();
    let original = (
        preview.state_hash(),
        preview.bulk_viscosity_at([0, 0, 0]).unwrap(),
    );
    drop(
        preview
            .preview_characteristic_interval([-0.1, 0.03], 3, 144, || false)
            .unwrap(),
    );
    assert_eq!(
        original,
        (
            preview.state_hash(),
            preview.bulk_viscosity_at([0, 0, 0]).unwrap()
        )
    );
    let calls = std::cell::Cell::new(0);
    assert!(
        preview
            .preview_characteristic_interval([-0.1, 0.03], 3, 144, || {
                let n = calls.get();
                calls.set(n + 1);
                n >= 2
            })
            .is_err()
    );
    assert_eq!(
        original,
        (
            preview.state_hash(),
            preview.bulk_viscosity_at([0, 0, 0]).unwrap()
        )
    );
    direct.set_incoming_characteristics([-0.1, 0.03]).unwrap();
    for _ in 0..3 {
        direct.step().unwrap();
    }
    preview
        .preview_characteristic_interval([-0.1, 0.03], 3, 144, || false)
        .unwrap()
        .accept();
    assert_eq!(direct.state_hash(), preview.state_hash());
    assert_eq!(direct.cumulative_ledger(), preview.cumulative_ledger());
    for x in 0..12 {
        assert_eq!(
            direct.bulk_viscosity_at([x, 0, 0]),
            preview.bulk_viscosity_at([x, 0, 0])
        );
    }
}
