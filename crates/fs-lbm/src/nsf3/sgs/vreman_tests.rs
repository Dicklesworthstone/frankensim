//! Vreman through the existing conservative NSF stress-work/heat-flux owner.
use super::*;
use crate::d3q27::D3q27;
use crate::sgs3::SgsUse;

fn spec(workers: usize) -> NsfSpec {
    NsfSpec {
        dims: [8; 3], periodic: [true; 3], wall: WallThermal::Adiabatic,
        gas: GasSpec::dry_air_ussa1976(), temperature_k: 293.15, pressure_pa: 101_325.0,
        dx_m: 5e-5, temperature_window_k: [273.15, 333.15], maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2, bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012, Phys. Fluids 24, 066102, N2 at 293 K",
        maximum_cells: 32 * 32 * 32, workers, transport_points: 9,
    }
}

fn config(cv: f64) -> NsfSgsConfig {
    NsfSgsConfig { closure: SgsConfig { model: SgsModel::Vreman { cv },
        filter_width_lu: 1.0, maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration }, turbulent_prandtl: 0.9 }
}

fn make(workers: usize) -> NsfGrid3<D3q27> {
    NsfGrid3::<D3q27>::new(spec(workers), |_| false, |x| {
        let a = x.map(|v| core::f64::consts::TAU * (v as f64 + 0.5) / 8.0);
        (1.0, [0.02 * a[0].sin() * a[1].cos() * a[2].cos(),
            -0.02 * a[0].cos() * a[1].sin() * a[2].cos(), 0.0], 293.15 + a[2].sin())
    }).unwrap()
}

#[test]
fn vreman_zero_keeps_the_molecular_energy_and_population_trajectory() {
    let mut ordinary = make(1);
    let mut zero = make(1).with_sgs(config(0.0)).unwrap();
    assert_eq!(zero.sgs_config(), None);
    for _ in 0..8 {
        assert_eq!(ordinary.step().unwrap(), zero.step().unwrap());
        assert_eq!(ordinary.state_hash(), zero.state_hash());
    }
    let mut forbidden = config(0.07);
    forbidden.closure.usage = SgsUse::AcceptanceRender;
    assert!(matches!(make(1).with_sgs(forbidden),
        Err(NsfError::Sgs(SgsError::AcceptanceRenderForbidden))));
}

#[test]
fn vreman_uses_matching_stress_work_and_heat_without_changing_bulk_or_worker_bits() {
    let mut serial = make(1).with_sgs(config(0.07)).unwrap();
    let mut pooled = make(3).with_sgs(config(0.07)).unwrap();
    let energy = serial.total_energy_j();
    let mut active = 0;
    for i in 0..serial.n {
        let transport = serial.sgs_transport[i];
        if transport.eddy_viscosity_lu > 0.0 {
            active += 1;
            let (rho, _, _) = serial.cell(serial.at(i));
            let molecular = serial.rates(rho, serial.theta[i], None);
            let modeled = serial.sgs_rates(rho, serial.theta[i], None, Some(&transport));
            assert!(modeled.0 < molecular.0);
            assert_eq!(modeled.1.to_bits(), molecular.1.to_bits());
            let expected = rho * serial.law.gamma / (serial.law.gamma - 1.0)
                * transport.eddy_viscosity_lu / 0.9;
            assert!((transport.turbulent_conductivity_lu / expected - 1.0).abs() < 1e-14);
        }
    }
    assert!(active > 0, "the closure must reach actual transport");
    for _ in 0..16 {
        let a = serial.step().unwrap();
        let b = pooled.step().unwrap();
        assert_eq!(a, b);
        assert_eq!(serial.state_hash(), pooled.state_hash());
        assert!(a.residual_rel < 1e-12);
        assert!((serial.total_energy_j() / energy - 1.0).abs() < 2e-12);
        assert_eq!(a.wall_heat_j, 0.0);
        assert_eq!(a.source_energy_j, 0.0);
    }
}
