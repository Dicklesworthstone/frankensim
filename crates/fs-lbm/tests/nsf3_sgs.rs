//! SGS uses the original temperature/pressure-dependent source and caloric
//! energy owner. The reused fixture is authored NASA9 data, not a gas benchmark.
#[path = "nsf3_thermochemical/source.rs"]
mod source;

use fs_lbm::d3q27::D3q27;
use fs_lbm::nsf3::{NsfGrid3, NsfSgsConfig, NsfSpec, ThermochemicalNsfGas, WallThermal};
use fs_lbm::sgs3::{SgsConfig, SgsModel, SgsUse};
use fs_material::thermochemical_gas::energy::GasEnergySolve;

#[test]
fn local_nasa9_heat_capacity_sets_sgs_heat_flux_and_energy_stays_conservative() {
    let gas = ThermochemicalNsfGas::new(
        source::model(-4000.0, 250000.0, 1.0),
        source::point(300.0, 100000.0),
        GasEnergySolve {
            minimum_temperature_k: 250.0,
            maximum_temperature_k: 600.0,
            energy_tolerance_j_kg: 1e-6,
            maximum_iterations: 64,
        },
    )
    .unwrap();
    let mut grid = NsfGrid3::<D3q27>::new(
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
            bulk_ratio_source: "authored synthetic ratio; not empirical",
            maximum_cells: 32,
            wall: WallThermal::Adiabatic,
            workers: 1,
            transport_points: 5,
        },
        |_| false,
        |x| {
            let phase = core::f64::consts::TAU * (x[0] as f64 + 0.5) / 8.0;
            let temperature = 300.0 + 10.0 * phase.cos();
            (
                300.0 / temperature,
                [0.0, 0.02 * phase.sin(), 0.0],
                temperature,
            )
        },
    )
    .unwrap()
    .with_sgs(NsfSgsConfig {
        closure: SgsConfig {
            model: SgsModel::Smagorinsky { cs: 0.17 },
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular: 10.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.85,
    })
    .unwrap();
    let energy = grid.total_energy_j();
    for _ in 0..8 {
        for x in 0..8 {
            let location = [x, 0, 0];
            let (rho, _, temperature) = grid.cell(location);
            let material = grid.thermochemical_state(location).unwrap().unwrap();
            let transport = grid.sgs_transport(location).unwrap().unwrap();
            assert!(transport.eddy_viscosity_lu > 0.0);
            let cp_r = material.state().specific_heat_cp / grid.law().gas_constant;
            // cp/R = 3 + .001 T in the existing synthetic source; using
            // reference cp or cv would violate this actual-cell identity.
            assert!((cp_r - (3.0 + 0.001 * temperature)).abs() < 1e-12);
            let wanted = rho * cp_r * transport.eddy_viscosity_lu / 0.85;
            assert!((transport.turbulent_conductivity_lu / wanted - 1.0).abs() < 1e-12);
        }
        let ledger = grid.step().unwrap();
        assert!(ledger.residual_rel < 1e-12);
        assert!((grid.total_energy_j() / energy - 1.0).abs() < 2e-12);
        assert_eq!(ledger.source_energy_j, 0.0);
        assert_eq!(ledger.wall_heat_j, 0.0);
    }
}
