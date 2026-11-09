use super::*;
use crate::d3q27::D3q27;
use crate::nsf3::NsfOpenState;
use crate::sgs3::SgsUse;
use crate::sgs3::dynamic::DynamicAveraging;

const SHAPE: [usize; 3] = [5; 3];
const CELLS: usize = 125;

fn settings(normal: usize) -> NsfSpec {
    NsfSpec {
        dims: SHAPE,
        periodic: core::array::from_fn(|a| a != normal),
        wall: WallThermal::Adiabatic,
        gas: GasSpec::dry_air_ussa1976(),
        temperature_k: 293.15,
        pressure_pa: 101_325.0,
        dx_m: 5e-5,
        temperature_window_k: [273.15, 333.15],
        maximum_mach: 0.2,
        maximum_temperature_ratio: 0.2,
        bulk_to_shear_ratio: 0.73,
        bulk_ratio_source: "authored reservoir dynamic SGS fixture",
        maximum_cells: CELLS,
        workers: 1,
        transport_points: 9,
    }
}

fn controls() -> NsfDynamicThermalSgsConfig {
    NsfDynamicThermalSgsConfig {
        momentum: DynamicSgsConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: DynamicScalarConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn mixed(p: [usize; 3]) -> (f64, [f64; 3], f64) {
    let i = p[0] + 5 * (p[1] + 5 * p[2]);
    (
        1.0 + ((i * 7 + 3) % 11) as f64 / 2048.0,
        core::array::from_fn(|a| {
            (((i * (17 + 8 * a) + i * i * (a + 1) + 13 * a) % 127) as f64 - 63.0) / 32768.0
        }),
        293.15 + ((i * i * 3 + 7 * i) % 41) as f64 / 128.0,
    )
}

fn exterior(grid: &NsfGrid3<D3q27>, normal: usize, high: bool) -> NsfOpenState {
    let mut velocity = [0.0; 3];
    velocity[normal] = if high { 0.0004 } else { 0.0006 };
    velocity[(normal + 1) % 3] = if high { -0.0002 } else { 0.0003 };
    NsfOpenState {
        density_kg_m3: grid.law.rho_ref_kg_m3 * if high { 1.001 } else { 1.0 },
        velocity_m_s: velocity.map(|v| v * grid.law.dx_m / grid.law.dt_s),
        temperature_k: if high { 293.3 } else { 293.0 },
    }
}

fn molecular(normal: usize) -> NsfGrid3<D3q27> {
    let mut grid = NsfGrid3::<D3q27>::new(settings(normal), |_| false, mixed).unwrap();
    for high in [false, true] {
        let state = exterior(&grid, normal, high);
        grid.set_open_face(normal, high, Some(state)).unwrap();
    }
    grid
}

fn fitted(normal: usize) -> NsfGrid3<D3q27> {
    molecular(normal)
        .with_reservoir_dynamic_thermal_sgs(controls(), 200 * CELLS)
        .unwrap()
}

fn close(a: f64, b: f64) {
    assert!(
        (a - b).abs() < 1e-11 * b.abs().max(1e-8),
        "{a:.17e} != {b:.17e}"
    );
}

fn verify_current_transport(grid: &NsfGrid3<D3q27>, normal: usize) {
    let (rho, momentum) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &momentum, &grid.force);
    let ghosts = grid.reservoir_faces[normal].map(Option::unwrap);
    let mut core = ReservoirDynamicSgs::new(SHAPE, [1.0; 3], normal, usize::MAX).unwrap();
    let mut scalar = ReservoirDynamicScalar::new(SHAPE, normal, usize::MAX).unwrap();
    core.prepare_with_scalar(
        controls().momentum,
        controls().heat,
        &rho,
        &velocity,
        ghosts.map(|g| g.u),
        &grid.theta,
        ghosts.map(|g| g.theta),
        &mut scalar,
        || false,
    )
    .unwrap();
    for i in 0..CELLS {
        let material = grid.material.get(i).and_then(Option::as_ref);
        let cp_over_r = material.map_or(grid.law.gamma / (grid.law.gamma - 1.0), |m| {
            m.state.specific_heat_cp / grid.law.gas_constant
        });
        close(
            grid.sgs_transport[i].eddy_viscosity_lu,
            core.viscosities()[i],
        );
        close(
            grid.sgs_transport[i].turbulent_conductivity_lu,
            rho[i] * cp_over_r * scalar.diffusivities()[i],
        );
        assert_eq!(
            grid.sgs_rates(
                rho[i],
                grid.theta[i],
                material,
                Some(&grid.sgs_transport[i])
            )
            .1,
            grid.rates(rho[i], grid.theta[i], material).1,
        );
    }
}

#[test]
fn g1_reservoir_dynamic_nsf_uses_current_ghosts_and_preserves_physical_force_chart() {
    for normal in 0..3 {
        let mut grid = fitted(normal);
        verify_current_transport(&grid, normal);
        assert!(grid.sgs_transport.iter().any(|t| t.eddy_viscosity_lu > 0.0));
        assert!(
            grid.sgs_transport
                .iter()
                .any(|t| t.turbulent_conductivity_lu > 0.0)
        );
        let transport = grid.sgs_transport.clone();
        let (rho, momentum) = grid.macro_fields(&grid.f);
        let velocity = grid.velocity_field(&rho, &momentum, &grid.force);
        let energy = grid.e.clone();
        let mut changed = exterior(&grid, normal, false);
        changed.temperature_k = 294.0;
        changed.velocity_m_s[normal] *= 2.0;
        grid.set_open_face(normal, false, Some(changed)).unwrap();
        verify_current_transport(&grid, normal);
        assert_ne!(grid.sgs_transport, transport);
        let (after_rho, after_momentum) = grid.macro_fields(&grid.f);
        let after_velocity = grid.velocity_field(&after_rho, &after_momentum, &grid.force);
        for i in 0..CELLS {
            close(after_rho[i], rho[i]);
            for a in 0..3 {
                // Populations are order one, including where the physical
                // velocity is exactly zero; use their roundoff scale here.
                assert!((after_velocity[i][a] - velocity[i][a]).abs() < 3e-15);
            }
        }
        assert_eq!(grid.e, energy);

        let before = grid.state_hash();
        let transport = grid.sgs_transport.clone();
        let force = grid.force.clone();
        let divergence = grid.div.clone();
        changed.temperature_k = 400.0;
        assert!(grid.set_open_face(normal, false, Some(changed)).is_err());
        assert!(grid.set_open_face(normal, false, None).is_err());
        assert_eq!(grid.state_hash(), before);
        assert_eq!(grid.sgs_transport, transport);
        assert_eq!(grid.force, force);
        assert_eq!(grid.div, divergence);
    }
}

#[test]
fn g0_g3_g4_reservoir_dynamic_nsf_admission_and_reversible_native_steps() {
    assert!(matches!(
        molecular(0).with_reservoir_dynamic_thermal_sgs(controls(), 200 * CELLS - 1),
        Err(NsfError::DynamicSgs(DynamicSgsError::WorkspaceBytes { .. }))
    ));
    let mut incomplete = NsfGrid3::<D3q27>::new(settings(0), |_| false, mixed).unwrap();
    let state = exterior(&incomplete, 0, false);
    incomplete.set_open_face(0, false, Some(state)).unwrap();
    assert!(
        incomplete
            .with_reservoir_dynamic_thermal_sgs(controls(), 200 * CELLS)
            .is_err()
    );
    let mut forbidden = controls();
    forbidden.heat.usage = SgsUse::AcceptanceRender;
    assert!(
        molecular(0)
            .with_reservoir_dynamic_thermal_sgs(forbidden, 200 * CELLS)
            .is_err()
    );

    let mut grid = fitted(0);
    let mut direct = fitted(0);
    let before = grid.state_hash();
    let transport = grid.sgs_transport.clone();
    assert!(
        grid.preview_thermal_wall_step(Vec::new(), CELLS, || true)
            .is_err()
    );
    assert_eq!(grid.state_hash(), before);
    assert_eq!(grid.sgs_transport, transport);
    let expected = direct.step().unwrap();
    assert!(expected.open_energy_j.is_finite());
    assert!(expected.residual_rel < 1e-10);
    let trial = grid
        .preview_thermal_wall_step(Vec::new(), CELLS, || false)
        .unwrap();
    assert_eq!(*trial.energy_ledger(), expected);
    assert_eq!(
        trial.grid().reservoir_dynamic_thermal_sgs_config(),
        Some(controls())
    );
    assert_eq!(
        trial.grid().dynamic_sgs_workspace_bytes(),
        Some(200 * CELLS)
    );
    assert_eq!(trial.grid().state_hash(), direct.state_hash());
    assert_eq!(trial.grid().sgs_transport, direct.sgs_transport);
    drop(trial);
    assert_eq!(grid.state_hash(), before);
    assert_eq!(grid.sgs_transport, transport);
    grid.preview_thermal_wall_step(Vec::new(), CELLS, || false)
        .unwrap()
        .accept();
    assert_eq!(grid.state_hash(), direct.state_hash());
    assert_eq!(grid.sgs_transport, direct.sgs_transport);

    // A reservoir uses the full adjacent modeled conductivity, whereas an
    // interior face averages its two cells. This field passes the latter's
    // 1000x cap but must fail the actual exterior molecular-face cap.
    let (rho, momentum) = grid.macro_fields(&grid.f);
    let velocity = grid.velocity_field(&rho, &momentum, &grid.force);
    let ghost = grid.reservoir_faces[0][0].unwrap();
    let conductivity = grid
        .spec
        .gas
        .prepare_face_conductivity(
            &grid.law,
            0.5 * (grid.law.pressure_pa(rho[0], grid.theta[0])
                + grid.law.pressure_pa(ghost.rho, ghost.theta)),
            [grid.theta[0], ghost.theta],
        )
        .unwrap();
    let molecular = conductivity
        .mean_lattice(
            grid.theta[0],
            ghost.theta,
            0.5 * (grid.theta[0] + ghost.theta),
        )
        .unwrap();
    let mut excessive = vec![
        NsfSgsTransport {
            eddy_viscosity_lu: 0.0,
            turbulent_conductivity_lu: 0.0,
        };
        CELLS
    ];
    excessive[0].turbulent_conductivity_lu = 1001.0 * molecular;
    let error = grid
        .flux_divergence(
            &grid.e,
            &rho,
            &velocity,
            &grid.theta,
            &vec![[0.0; 6]; CELLS],
            &excessive,
            None,
        )
        .err()
        .unwrap();
    let NsfError::SgsHeatRatio {
        cell,
        ratio,
        maximum,
    } = error
    else {
        panic!("expected the actual prescribed-reservoir heat-ratio gate");
    };
    assert_eq!(cell, 0);
    close(ratio, 1001.0);
    assert_eq!(maximum, 1000.0);

    // The initial scalar covariance lives only on the low plane. Its warmer
    // interior/cooler ghost gives a negative Germano heat numerator, hence a
    // clipped-zero fit. Reversing that normal scalar gradient makes the fit
    // positive and must refuse the very small admitted heat ratio atomically.
    let mut cold = NsfGrid3::<D3q27>::new(
        settings(0),
        |_| false,
        |p| {
            let mode = (core::f64::consts::TAU * (p[1] as f64 + 0.5) / 5.0).cos();
            (
                1.0,
                [0.001 * mode, 0.0, 0.0],
                293.15 + if p[0] == 0 { 0.25 * mode } else { 0.0 },
            )
        },
    )
    .unwrap();
    let mut state = NsfOpenState {
        density_kg_m3: cold.law.rho_ref_kg_m3,
        velocity_m_s: [0.0; 3],
        temperature_k: 292.15,
    };
    cold.set_open_face(0, false, Some(state)).unwrap();
    cold.set_open_face(
        0,
        true,
        Some(NsfOpenState {
            temperature_k: 293.15,
            ..state
        }),
    )
    .unwrap();
    let mut cap = controls();
    cap.heat.maximum_diffusivity_to_molecular = 1e-12;
    let mut cold = cold
        .with_reservoir_dynamic_thermal_sgs(cap, 200 * CELLS)
        .unwrap();
    let before = cold.state_hash();
    let transport = cold.sgs_transport.clone();
    let force = cold.force.clone();
    let divergence = cold.div.clone();
    state.temperature_k = 294.15;
    assert!(matches!(
        cold.set_open_face(0, false, Some(state)),
        Err(NsfError::SgsHeatRatio { .. })
    ));
    assert_eq!(cold.state_hash(), before);
    assert_eq!(cold.sgs_transport, transport);
    assert_eq!(cold.force, force);
    assert_eq!(cold.div, divergence);
}
