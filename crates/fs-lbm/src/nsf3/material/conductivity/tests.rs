use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::{
    NsfBoundaryDrive, NsfCharacteristicAxis, NsfCharacteristicBoundary,
    NsfCharacteristicThermalPolicy, NsfGrid3, NsfSpec, NsfThermalWallFace, NsfWallThermal,
    ThermochemicalNsfGas, WallThermal,
};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

#[allow(dead_code)]
#[path = "../../../../tests/nsf3_thermochemical/source.rs"]
mod source;

const KNOTS: [(f64, f64); 5] = [
    (200.0, 0.02),
    (300.0, 0.03),
    (320.0, 0.09),
    (400.0, 0.10),
    (700.0, 0.13),
];

fn control(limit: usize) -> NsfSourceConductivityControl {
    NsfSourceConductivityControl {
        quadrature: GasConductivityQuadrature {
            absolute_tolerance_w_m_k: 1e-12,
            relative_tolerance: 1e-11,
            maximum_subdivisions: 64,
            maximum_depth: 16,
        },
        maximum_evaluations_per_face: limit,
    }
}

fn gas(offset: f64, limit: Option<usize>) -> ThermochemicalNsfGas {
    let gas = ThermochemicalNsfGas::new(
        source::conductivity_model(offset, 250000.0, KNOTS.to_vec()),
        source::point(300.0, 100000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 64,
        },
    )
    .unwrap();
    match limit {
        Some(limit) => gas.with_source_conductivity(control(limit)).unwrap(),
        None => gas,
    }
}

fn spec(gas: ThermochemicalNsfGas, workers: usize) -> NsfSpec<ThermochemicalNsfGas> {
    NsfSpec {
        dims: [8, 2, 2],
        periodic: [true; 3],
        gas,
        temperature_k: 300.0,
        pressure_pa: 100000.0,
        dx_m: 5e-5,
        temperature_window_k: [250.0, 600.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 1.0,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored nonlinear conductivity fixture",
        maximum_cells: 32,
        wall: WallThermal::Adiabatic,
        workers,
        transport_points: 5,
    }
}

/// Independent exact integral of the authored affine segments.
fn integral(from: f64, to: f64) -> f64 {
    if from > to {
        return -integral(to, from);
    }
    let mut sum = 0.0;
    for points in KNOTS.windows(2) {
        let [(a, ka), (b, kb)] = [points[0], points[1]];
        let low = from.max(a);
        let high = to.min(b);
        if high > low {
            let slope = (kb - ka) / (b - a);
            sum += (high - low) * (ka + slope * (0.5 * (low + high) - a));
        }
    }
    sum
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance * expected.abs().max(1e-12),
        "actual={actual:e}, expected={expected:e}"
    );
}

#[test]
fn g1_native_wide_fourier_faces_differentiate_the_original_conductivity_primitive() {
    let temperatures = [290.0, 310.0, 345.0, 300.0, 330.0, 315.0, 305.0, 295.0];
    let grid = NsfGrid3::<D3q27>::new(
        spec(gas(0.0, Some(10000)), 3),
        |_| false,
        |x| {
            let t = temperatures[x[0]];
            (300.0 / t, [0.0; 3], t)
        },
    )
    .unwrap();
    let rho: Vec<_> = (0..grid.n)
        .map(|i| 300.0 / temperatures[grid.at(i)[0]])
        .collect();
    let (divergence, wall, boundary) = grid
        .flux_divergence(
            &grid.e,
            &rho,
            &vec![[0.0; 3]; grid.n],
            &grid.theta,
            &vec![[0.0; 6]; grid.n],
            &[],
            None,
        )
        .unwrap();
    let scale =
        grid.law.transport_scale() / grid.law.gas_constant / (3.0 * grid.law.temperature_ref_k);
    let flux = |x: usize| {
        let near = integral(temperatures[x], temperatures[(x + 1) % 8]);
        let wide = integral(temperatures[(x + 7) % 8], temperatures[(x + 2) % 8]);
        -(27.0 * near - wide) * scale / 24.0
    };
    for (i, actual) in divergence.iter().enumerate() {
        let x = grid.at(i)[0];
        close(*actual, flux(x) - flux((x + 7) % 8), 2e-11);
    }
    assert_eq!(wall, 0.0);
    assert_eq!(boundary.wall_heat_out, 0.0);
    assert!(divergence.iter().sum::<f64>().abs() < 1e-15);
    // The 300..345 K interval crosses a slope change, so a midpoint
    // coefficient cannot accidentally satisfy this independent oracle.
    let midpoint_k = 0.09 + (322.5 - 320.0) * 0.01 / 80.0;
    assert!((integral(300.0, 345.0) / 45.0 - midpoint_k).abs() > 0.005);
}

#[test]
fn g1_film_and_isothermal_walls_share_the_original_half_cell_integral() {
    let mut s = spec(gas(0.0, Some(10000)), 1);
    s.periodic[0] = false;
    let mut grid = NsfGrid3::<D3q27>::new(s, |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap();
    let mut low = 300.0;
    let mut high = 350.0;
    let coefficient = 4000.0;
    for _ in 0..80 {
        let t = 0.5 * (low + high);
        if 2.0 * integral(300.0, t) / grid.law.dx_m > coefficient * (350.0 - t) {
            high = t;
        } else {
            low = t;
        }
    }
    let surface = 0.5 * (low + high);
    assert!(surface > 320.0 && surface < 350.0);
    let expected = coefficient * (350.0 - surface) / grid.law.heat_flux_unit_w_m2();
    let rho = vec![1.0; grid.n];
    for (high, x) in [(false, 0), (true, 7)] {
        for thermal in [
            NsfWallThermal::Film {
                coefficient_w_m2_k: coefficient,
                ambient_temperature_k: 350.0,
            },
            NsfWallThermal::Isothermal(surface),
        ] {
            grid.set_thermal_walls(vec![NsfThermalWallFace {
                cell: [x, 0, 0],
                axis: 0,
                high,
                thermal,
            }])
            .unwrap();
            let actual = grid.thermal_wall_fluxes(&rho, &grid.theta).unwrap();
            close(actual[0], if high { -expected } else { expected }, 2e-10);
        }
    }
    let before = grid.total_energy_j();
    let step = grid.step().unwrap();
    close(grid.total_energy_j() - before, step.wall_heat_j, 2e-8);
    assert!(step.residual_rel < 1e-11);
}

#[test]
fn g1_characteristic_inflow_carries_upstream_enthalpy_once_with_source_integrated_heat() {
    let mut s = spec(gas(-4000.0, Some(10000)), 1);
    s.periodic[0] = false;
    let mut grid = NsfGrid3::<D3q27>::new(s, |_| false, |_| (1.0, [0.005, 0.0, 0.0], 300.0))
        .unwrap()
        .with_characteristic_axis(NsfCharacteristicAxis {
            axis: 0,
            // The linear characteristic bound must lie in (0, 0.1] (open_boundary.rs, c2d1f74449, 2026-10-24); the fixture's
            // 0.2 has been refused since the test was written (5b7be2124b, 2026-11-09). 0.1 is the largest admitted value.
            maximum_relative_perturbation: 0.1,
            faces: [NsfCharacteristicBoundary {
                reference_pressure_pa: 100000.0,
                reference_temperature_k: 300.0,
                reference_velocity_m_s: [0.0; 3],
                incoming_characteristic_m_s: 0.0,
            }; 2],
        })
        .unwrap()
        .with_characteristic_thermal_policy([NsfCharacteristicThermalPolicy::EnthalpyFlux; 2])
        .unwrap();
    // The inflow temperature is admitted only within the characteristic bound of its 300 K reference: |T / 300 - 1| <= 0.1
    // (the inflow check of open_boundary.rs, 2a01cd338d, 2026-10-30, against the bound of c2d1f74449). The fixture's 340 K is
    // 13.3 % off. 330 K is not admissible either: 330 / 300 - 1 is 0.10000000000000009 in floating point, over 0.1. 325 K is
    // 8.3 %.
    let inflow_temperature_k = 325.0;
    grid.set_boundary_drives([
        NsfBoundaryDrive::PressureTemperature {
            pressure_offset_pa: 0.0,
            inflow_temperature_k,
        },
        NsfBoundaryDrive::PressureOffset(0.0),
    ])
    .unwrap();
    let (rho, _) = grid.macro_fields(&grid.f);
    let stress = vec![[0.0; 6]; grid.n];
    for trace in grid.open_boundary_traces(false).unwrap() {
        let i = grid.index(trace.cell);
        let face = grid
            .characteristic_face(grid.open_faces.as_ref(), i, 0, -1)
            .unwrap();
        assert!(face.u[0] > 0.0);
        assert!(trace.temperature_k > 300.0 && trace.temperature_k < inflow_temperature_k);
        let (total, _, advected) = grid
            .open_energy_flux(face, i, 0, -1, &rho, &grid.theta, &stress, &[])
            .unwrap();
        let source = match &grid.spec.gas {
            NsfGasMedium::Thermochemical(source) => source,
            _ => unreachable!(),
        };
        let point = source
            .context()
            .clone()
            .with("temperature", inflow_temperature_k)
            .unwrap()
            .with("pressure", trace.pressure_pa)
            .unwrap();
        let upstream = source.model().evaluate(&point).unwrap();
        let h = (upstream.specific_internal_energy_j_kg() - source.energy_reference_j_kg()
            + upstream.state().specific_gas_constant * inflow_temperature_k)
            / grid.law.specific_energy_unit_j_kg();
        let kinetic = 0.5 * face.u.iter().map(|u| u * u).sum::<f64>();
        close(total, face.rho * face.u[0] * (h + kinetic), 2e-10);
        let q = 2.0 * integral(300.0, trace.temperature_k) * grid.law.transport_scale()
            / grid.law.gas_constant
            / (3.0 * grid.law.temperature_ref_k);
        close(total - advected, q, 2e-9);
    }
}

#[test]
fn g0_source_domain_and_shared_root_budget_refuse_before_wall_state_is_published() {
    let mut s = spec(gas(0.0, Some(40)), 1);
    s.periodic[0] = false;
    let mut grid = NsfGrid3::<D3q27>::new(s, |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap();
    let before = grid.state_hash();
    let result = grid.set_thermal_walls(vec![NsfThermalWallFace {
        cell: [0, 0, 0],
        axis: 0,
        high: false,
        thermal: NsfWallThermal::Film {
            coefficient_w_m2_k: 4000.0,
            ambient_temperature_k: 350.0,
        },
    }]);
    assert!(matches!(result, Err(NsfError::Gas(ref reason)) if reason.contains("work allowance")));
    assert_eq!(before, grid.state_hash());
    assert!(grid.thermal_walls().is_empty());
    assert!(
        grid.spec
            .gas
            .prepare_face_conductivity(
                &grid.law,
                250000.0_f64.next_up(),
                [
                    grid.law.theta_of_kelvin(300.0),
                    grid.law.theta_of_kelvin(350.0)
                ],
            )
            .is_err()
    );
    assert!(gas(0.0, None).with_source_conductivity(control(0)).is_err());
    let mut invalid = control(100);
    invalid.quadrature.relative_tolerance = f64::NAN;
    assert!(gas(0.0, None).with_source_conductivity(invalid).is_err());
}

#[test]
fn g3_g5_source_integrated_native_steps_ignore_formation_datum_and_worker_count() {
    let make = |offset, workers| {
        NsfGrid3::<D3q27>::new(
            spec(gas(offset, Some(10000)), workers),
            |_| false,
            |x| {
                let t = 320.0 + 0.3 * (core::f64::consts::TAU * x[0] as f64 / 8.0).sin();
                (300.0 / t, [0.001, 0.0, 0.0], t)
            },
        )
        .unwrap()
    };
    let mut first = make(0.0, 1);
    let mut parallel = make(0.0, 3);
    let mut shifted = make(1e12, 1);
    for _ in 0..4 {
        let a = first.step().unwrap();
        let b = parallel.step().unwrap();
        let c = shifted.step().unwrap();
        assert_eq!(a, b);
        assert_eq!(a, c);
        assert_eq!(first.state_hash(), parallel.state_hash());
        assert_eq!(first.state_hash(), shifted.state_hash());
        assert!(a.residual_rel < 1e-11);
    }
}
