use super::*;
use crate::d3q27::D3q27;
use crate::sgs3::SgsUse;

fn spec(dims: [usize; 3], periodic: [bool; 3], wall: WallThermal) -> NsfSpec {
    NsfSpec {
        dims,
        periodic,
        wall,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "Cramer 2012, Phys. Fluids 24, 066102, N2 at 293 K",
        maximum_cells: 32 * 32 * 32,
        workers: 1,
        transport_points: 9,
    }
}

fn config(model: SgsModel) -> NsfSgsConfig {
    NsfSgsConfig {
        closure: SgsConfig {
            model,
            filter_width_lu: 1.0,
            maximum_eddy_to_molecular: 100.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.9,
    }
}

fn rest(_: [usize; 3]) -> (f64, [f64; 3], f64) {
    (1.0, [0.0; 3], 293.15)
}

#[test]
fn moving_halfway_gradient_and_effective_shear_keep_bulk_molecular() {
    let mut grid = NsfGrid3::<D3q27>::new(
        spec([4, 8, 4], [true, false, true], WallThermal::Adiabatic),
        |_| false,
        |x| (1.0, [0.02 * (x[1] as f64 + 0.5) / 8.0, 0.0, 0.0], 293.15),
    )
    .unwrap();
    grid.set_edge_velocity(1, true, [0.02, 0.0, 0.0]).unwrap();
    let grid = grid
        .with_sgs(config(SgsModel::Smagorinsky { cs: 0.17 }))
        .unwrap();
    let expected = 0.17_f64.powi(2) * 0.02 / 8.0;
    for y in 0..8 {
        let x = [1, y, 1];
        let i = grid.index(x);
        let transport = grid.sgs_transport(x).unwrap().unwrap();
        assert!((transport.eddy_viscosity_lu / expected - 1.0).abs() < 1e-11);
        let (rho, _, _) = grid.cell(x);
        let molecular = grid.rates(rho, grid.theta[i], None);
        let modeled = grid.sgs_rates(rho, grid.theta[i], None, Some(&transport));
        assert!(modeled.0 < molecular.0);
        assert_eq!(modeled.1.to_bits(), molecular.1.to_bits());
        let want_heat =
            rho * grid.law.gamma / (grid.law.gamma - 1.0) * transport.eddy_viscosity_lu / 0.9;
        assert!((transport.turbulent_conductivity_lu / want_heat - 1.0).abs() < 1e-14);
    }
    // A translating pure-shear wall is not a WALE production mechanism.
    let wale = grid.with_sgs(config(SgsModel::Wale { cw: 0.5 })).unwrap();
    assert!(
        wale.sgs_transport
            .iter()
            .all(|v| v.eddy_viscosity_lu < 1e-20)
    );
}

#[test]
fn heat_diffuses_conservatively_and_has_no_modeled_wall_flux() {
    let grid = NsfGrid3::<D3q27>::new(
        spec([8, 4, 4], [true; 3], WallThermal::Adiabatic),
        |_| false,
        rest,
    )
    .unwrap();
    let rho = vec![1.0; grid.n];
    let velocity = vec![[0.0; 3]; grid.n];
    let stress = vec![[0.0; 6]; grid.n];
    let theta: Vec<_> = (0..grid.n)
        .map(|i| THETA0 + 0.001 * (core::f64::consts::TAU * grid.at(i)[0] as f64 / 8.0).cos())
        .collect();
    let transport = vec![
        NsfSgsTransport {
            eddy_viscosity_lu: 0.001,
            turbulent_conductivity_lu: 0.004,
        };
        grid.n
    ];
    let molecular = grid
        .flux_divergence(&grid.e, &rho, &velocity, &theta, &stress, &[], None)
        .unwrap()
        .0;
    let modeled = grid
        .flux_divergence(&grid.e, &rho, &velocity, &theta, &stress, &transport, None)
        .unwrap()
        .0;
    assert!(modeled.iter().sum::<f64>().abs() < 1e-16);
    // dE/dt = -div: positive theta*div means temperature variance decays.
    assert!(
        theta
            .iter()
            .zip(modeled.iter().zip(&molecular))
            .map(|(t, (a, b))| t * (a - b))
            .sum::<f64>()
            > 0.0
    );

    let wall = NsfGrid3::<D3q27>::new(
        spec(
            [8, 4, 4],
            [true, false, true],
            WallThermal::Isothermal(303.15),
        ),
        |_| false,
        rest,
    )
    .unwrap();
    let molecular = wall
        .flux_divergence(&wall.e, &rho, &velocity, &wall.theta, &stress, &[], None)
        .unwrap()
        .1;
    let modeled = wall
        .flux_divergence(
            &wall.e,
            &rho,
            &velocity,
            &wall.theta,
            &stress,
            &transport,
            None,
        )
        .unwrap()
        .1;
    assert!(molecular < 0.0);
    assert_eq!(modeled.to_bits(), molecular.to_bits());
}

#[test]
fn closed_energy_worker_bits_and_exact_zero_control() {
    let initial = |x: [usize; 3]| {
        let a = x.map(|v| core::f64::consts::TAU * (v as f64 + 0.5) / 8.0);
        (
            1.0,
            [
                0.02 * a[0].sin() * a[1].cos() * a[2].cos(),
                -0.02 * a[0].cos() * a[1].sin() * a[2].cos(),
                0.0,
            ],
            293.15 + a[2].sin(),
        )
    };
    let make = |workers| {
        let mut s = spec([8; 3], [true; 3], WallThermal::Adiabatic);
        s.workers = workers;
        NsfGrid3::<D3q27>::new(s, |_| false, initial).unwrap()
    };
    let mut molecular = make(1);
    let mut zero = make(1)
        .with_sgs(config(SgsModel::Wale { cw: 0.0 }))
        .unwrap();
    for _ in 0..8 {
        assert_eq!(molecular.step().unwrap(), zero.step().unwrap());
        assert_eq!(molecular.state_hash(), zero.state_hash());
    }
    for model in [
        SgsModel::Smagorinsky { cs: 0.17 },
        SgsModel::Wale { cw: 0.5 },
    ] {
        let mut serial = make(1).with_sgs(config(model)).unwrap();
        let mut parallel = make(3).with_sgs(config(model)).unwrap();
        let energy = serial.total_energy_j();
        assert!(
            serial
                .sgs_transport
                .iter()
                .any(|v| v.eddy_viscosity_lu > 0.0)
        );
        for _ in 0..16 {
            let a = serial.step().unwrap();
            let b = parallel.step().unwrap();
            assert_eq!(a, b);
            assert_eq!(serial.state_hash(), parallel.state_hash());
            assert!(a.residual_rel < 1e-12);
            assert!((serial.total_energy_j() / energy - 1.0).abs() < 2e-12);
            assert_eq!(a.wall_heat_j, 0.0);
            assert_eq!(a.source_energy_j, 0.0);
        }
    }
}

#[test]
fn fixed_physical_shear_refinement_removes_modeled_transport() {
    // Same physical length, velocity and source gas viscosity: fixed Re.
    // This is a smooth-gradient closure rung, not developed-turbulence DNS.
    let mut mean = Vec::new();
    for n in [8, 16, 32] {
        let mut s = spec([n; 3], [true; 3], WallThermal::Adiabatic);
        s.dx_m = 8e-4 / n as f64;
        let grid = NsfGrid3::<D3q27>::new(
            s,
            |_| false,
            |x| {
                (
                    1.0,
                    [
                        0.02 * (core::f64::consts::TAU * (x[1] as f64 + 0.5) / n as f64).sin(),
                        0.0,
                        0.0,
                    ],
                    293.15,
                )
            },
        )
        .unwrap()
        .with_sgs(config(SgsModel::Smagorinsky { cs: 0.17 }))
        .unwrap();
        let physical = (0..n)
            .map(|y| {
                grid.sgs_transport([0, y, 0])
                    .unwrap()
                    .unwrap()
                    .eddy_viscosity_lu
            })
            .sum::<f64>()
            / n as f64
            * grid.law.dx_m.powi(2)
            / grid.law.dt_s;
        mean.push(physical);
    }
    for pair in mean.windows(2) {
        let ratio = pair[0] / pair[1];
        assert!((3.5..4.05).contains(&ratio), "refinement ratio {ratio}");
    }
}

#[test]
fn refusals_keep_complete_reservoirs_and_excess_transport_out_of_the_owner() {
    let make = || {
        NsfGrid3::<D3q27>::new(
            spec([6, 4, 4], [false, true, true], WallThermal::Adiabatic),
            |_| false,
            rest,
        )
        .unwrap()
    };
    let config = config(SgsModel::Smagorinsky { cs: 0.17 });
    let reservoir = NsfOpenState {
        density_kg_m3: make().law.rho_ref_kg_m3,
        velocity_m_s: [0.0; 3],
        temperature_k: 293.15,
    };
    let mut active = make().with_sgs(config).unwrap();
    let hash = active.state_hash();
    assert!(active.set_open_face(0, false, Some(reservoir)).is_err());
    assert_eq!(active.state_hash(), hash);
    let mut prescribed = make();
    prescribed.set_open_face(0, false, Some(reservoir)).unwrap();
    assert!(prescribed.with_sgs(config).is_err());
    let mut bounded = config;
    bounded.closure.maximum_eddy_to_molecular = 0.01;
    let mut active = make().with_sgs(bounded).unwrap();
    let hash = active.state_hash();
    assert!(matches!(
        active.set_edge_velocity(0, false, [0.0, 0.02, 0.0]),
        Err(NsfError::Sgs(SgsError::ViscosityRatioExceeded { .. }))
    ));
    assert_eq!(active.state_hash(), hash);
    assert_eq!(active.edge_velocity, [[[0.0; 3]; 2]; 3]);
    let mut envelope = config;
    envelope.closure.model = SgsModel::Smagorinsky { cs: 100.0 };
    envelope.closure.maximum_eddy_to_molecular = 1e20;
    let mut active = make().with_sgs(envelope).unwrap();
    assert!(matches!(
        active.set_edge_velocity(0, false, [0.0, 0.02, 0.0]),
        Err(NsfError::Viscosity { .. })
    ));
    let mut acceptance = config;
    acceptance.closure.usage = SgsUse::AcceptanceRender;
    assert!(matches!(
        make().with_sgs(acceptance),
        Err(NsfError::Sgs(SgsError::AcceptanceRenderForbidden))
    ));
    let mut bad_heat = config;
    bad_heat.turbulent_prandtl = 0.0;
    assert!(matches!(make().with_sgs(bad_heat), Err(NsfError::Input(_))));
}

mod open;
