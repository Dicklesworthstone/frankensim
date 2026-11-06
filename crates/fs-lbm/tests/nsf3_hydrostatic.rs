//! Physical hydrostatic fixed points of the original NSF step (synthetic G0/G1).
#[path = "nsf3_thermochemical/source.rs"]
mod source;

use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{NsfError, NsfGrid3, NsfSpec, ThermochemicalNsfGas, WallThermal};
use fs_material::gas::{GasSpec, GasState};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

fn spec() -> NsfSpec {
    NsfSpec {
        dims: [4, 3, 8],
        periodic: [true, true, false],
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 300.0,
        pressure_pa: 100_000.0,
        dx_m: 1_000.0,
        temperature_window_k: [250.0, 400.0],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.3,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012 N2 ratio; hydrostatic regression fixture",
        maximum_cells: 256,
        wall: WallThermal::Adiabatic,
        workers: 1,
        transport_points: 17,
    }
}

fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:.17e} != {expected:.17e}, tolerance {tolerance:e}"
    );
}

#[test]
fn native_column_stays_at_rest_on_every_gravity_axis_and_sign() {
    for axis in 0..3 {
        for g in [-9.80665, 9.80665] {
            let mut s = spec();
            s.dims = [3; 3];
            s.dims[axis] = 8;
            s.periodic = [true; 3];
            s.periodic[axis] = false;
            let mut flow =
                NsfGrid3::<D3q27>::new_isothermal_hydrostatic(s.clone(), |_| false, axis, g, 3)
                    .unwrap();
            let mass = flow.total_mass_kg();
            let energy = flow.total_energy_j();
            let initial: Vec<_> = (0..8)
                .map(|j| {
                    let mut x = [1; 3];
                    x[axis] = j;
                    flow.cell(x).0
                })
                .collect();
            close(initial[3], 1.0, 2e-15);
            let mut gravity = [0.0; 3];
            gravity[axis] = g;
            assert_eq!(flow.gravity_m_s2(), gravity);
            // The pressure datum is at the cell center, not the half-cell wall.
            for j in 0..7 {
                let rho0 = initial[j] * flow.law().rho_ref_kg_m3;
                let rho1 = initial[j + 1] * flow.law().rho_ref_kg_m3;
                close(
                    s.pressure_pa * (initial[j + 1] - initial[j]),
                    0.5 * g * s.dx_m * (rho0 + rho1),
                    1e-9,
                );
            }
            for _ in 0..32 {
                let ledger = flow.step().unwrap();
                assert!(ledger.residual_rel < 2e-12);
            }
            for j in 0..8 {
                let mut x = [1; 3];
                x[axis] = j;
                let (rho, u, t) = flow.cell(x);
                close(rho, initial[j], 3e-11 * initial[j]);
                for v in u {
                    close(v, 0.0, 3e-11);
                }
                close(t, s.temperature_k, 3e-8);
            }
            close(flow.total_mass_kg(), mass, 3e-12 * mass);
            close(flow.total_energy_j(), energy, 3e-12 * energy);
            close(flow.kinetic_energy_j(), 0.0, 1e-19 * energy);
            close(flow.cumulative_ledger().gravity_work_j, 0.0, 3e-12 * energy);
            close(
                flow.cumulative_ledger().viscous_heating_j,
                0.0,
                3e-12 * energy,
            );
        }
    }
}

#[test]
fn closed_walls_and_embedded_solids_preserve_rest_and_worker_identity() {
    let mut s = spec();
    s.periodic = [false; 3];
    s.wall = WallThermal::Isothermal(s.temperature_k);
    let solid = |x: [usize; 3]| x == [1, 1, 3] || x == [2, 1, 4];
    let mut serial =
        NsfGrid3::<D3q27>::new_isothermal_hydrostatic(s.clone(), solid, 2, -9.80665, 0).unwrap();
    s.workers = 3;
    let mut pooled =
        NsfGrid3::<D3q27>::new_isothermal_hydrostatic(s.clone(), solid, 2, -9.80665, 0).unwrap();
    for _ in 0..16 {
        serial.step().unwrap();
        pooled.step().unwrap();
        assert_eq!(serial.state_hash(), pooled.state_hash());
    }
    for z in 0..s.dims[2] {
        for y in 0..s.dims[1] {
            for x in 0..s.dims[0] {
                if solid([x, y, z]) {
                    continue;
                }
                let (_, u, t) = serial.cell([x, y, z]);
                for v in u {
                    close(v, 0.0, 3e-11);
                }
                close(t, 300.0, 3e-8);
            }
        }
    }
}

#[test]
fn profile_converges_quadratically_to_the_physical_barometric_law() {
    let height = 8_000.0;
    let g = -9.80665;
    let mut errors = Vec::new();
    for intervals in [8, 16, 32] {
        let mut s = spec();
        s.dims = [1, 1, intervals + 1];
        s.dx_m = height / intervals as f64;
        let flow =
            NsfGrid3::<D3q27>::new_isothermal_hydrostatic(s.clone(), |_| false, 2, g, 0).unwrap();
        let exact = (g * height / (flow.law().gas_constant * s.temperature_k)).exp();
        let actual = flow.cell([0, 0, intervals]).0;
        errors.push((actual - exact).abs());
        close(
            flow.gas_state([0, 0, 0]).unwrap().pressure,
            s.pressure_pa,
            1e-8,
        );
    }
    for pair in errors.windows(2) {
        assert!((3.9..4.1).contains(&(pair[0] / pair[1])), "{errors:?}");
    }
    assert!(errors[2] < 3e-5);
}

#[test]
fn removing_gravity_releases_the_supported_pressure_gradient() {
    let mut flow =
        NsfGrid3::<D3q27>::new_isothermal_hydrostatic(spec(), |_| false, 2, -9.80665, 0).unwrap();
    let a = 9.80665 * flow.law().dt_s.powi(2) / flow.law().dx_m;
    flow.set_gravity_si([0.0; 3]).unwrap();
    close(flow.cell([1, 1, 4]).1[2], 0.0, 2e-15);
    let ledger = flow.step().unwrap();
    assert!(flow.cell([1, 1, 4]).1[2] > 0.8 * a);
    assert_eq!(ledger.gravity_work_j, 0.0);
}

#[test]
fn source_gas_retains_its_caloric_chart_and_pressure_domain() {
    let gas = |maximum_pressure| {
        ThermochemicalNsfGas::new(
            source::model(1_000_000.0, maximum_pressure, 1.0),
            source::point(300.0, 100_000.0),
            GasEnergySolve {
                minimum_temperature_k: 250.0,
                maximum_temperature_k: 400.0,
                energy_tolerance_j_kg: 1e-6,
                maximum_iterations: 64,
            },
        )
        .unwrap()
    };
    let s = spec();
    let medium_spec = |gas| NsfSpec {
        dims: s.dims,
        periodic: s.periodic,
        gas,
        temperature_k: s.temperature_k,
        pressure_pa: s.pressure_pa,
        dx_m: s.dx_m,
        temperature_window_k: s.temperature_window_k,
        maximum_mach: s.maximum_mach,
        maximum_temperature_ratio: s.maximum_temperature_ratio,
        bulk_to_shear_ratio: s.bulk_to_shear_ratio,
        bulk_ratio_source: s.bulk_ratio_source,
        maximum_cells: s.maximum_cells,
        wall: s.wall,
        workers: s.workers,
        transport_points: s.transport_points,
    };
    let mut flow = NsfGrid3::<D3q27>::new_isothermal_hydrostatic(
        medium_spec(gas(250_000.0)),
        |_| false,
        2,
        9.80665,
        0,
    )
    .unwrap();
    let energy = flow.total_energy_j();
    for _ in 0..8 {
        flow.step().unwrap();
    }
    for j in 0..8 {
        let (_, u, t) = flow.cell([1, 1, j]);
        close(t, 300.0, 2e-6);
        for v in u {
            close(v, 0.0, 2e-9);
        }
    }
    close(flow.total_energy_j(), energy, 1e-10 * energy.abs());
    assert!(matches!(
        NsfGrid3::<D3q27>::new_isothermal_hydrostatic(
            medium_spec(gas(101_000.0)),
            |_| false,
            2,
            9.80665,
            0,
        ),
        Err(NsfError::Gas(_))
    ));
}

#[test]
fn zero_gravity_is_bit_exact_and_invalid_columns_refuse() {
    let mut s = spec();
    s.periodic = [true; 3];
    let normal = NsfGrid3::<D3q27>::new(s.clone(), |_| false, |_| (1.0, [0.0; 3], 300.0)).unwrap();
    let column =
        NsfGrid3::<D3q27>::new_isothermal_hydrostatic(s.clone(), |_| false, 2, 0.0, 3).unwrap();
    assert_eq!(column.state_hash(), normal.state_hash());
    let refused = |s, axis, g, datum| {
        assert!(
            NsfGrid3::<D3q27>::new_isothermal_hydrostatic(s, |_| false, axis, g, datum,).is_err()
        );
    };
    refused(s, 2, -9.80665, 0); // no periodic gravitational potential
    refused(spec(), 3, -9.80665, 0);
    refused(spec(), 2, -9.80665, 8);
    refused(spec(), 2, f64::NAN, 0);
    let mut hot_wall = spec();
    hot_wall.wall = WallThermal::Isothermal(310.0);
    refused(hot_wall, 2, -9.80665, 0);
    let state = GasState::try_new(&GasSpec::dry_air_ussa1976(), 300.0, 100_000.0).unwrap();
    let two_rt = 2.0 * state.specific_gas_constant * 300.0;
    refused(spec(), 2, 2.0 * two_rt / spec().dx_m, 0);
    for sign in [-1.0, 1.0] {
        let mut tall = spec();
        tall.dims = [1, 1, 200];
        refused(tall, 2, sign * 0.99 * two_rt / spec().dx_m, 0);
    }
    let mut capped = spec();
    capped.maximum_cells = 1;
    assert!(matches!(
        NsfGrid3::<D3q27>::new_isothermal_hydrostatic(capped, |_| false, 2, -9.80665, 0,),
        Err(NsfError::Cap { .. })
    ));
}
