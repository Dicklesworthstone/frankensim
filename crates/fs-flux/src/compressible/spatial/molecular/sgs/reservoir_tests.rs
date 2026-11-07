//! Opt-in physical bath SGS fluxes, source calorics and atomic integration.

use super::*;
use crate::compressible::spatial::molecular::sgs::{
    DynamicAveraging, DynamicSgsConfig, DynamicSubgridConfig, LagrangianSgsConfig,
    LagrangianSubgridConfig, SgsModel, SgsUse, SubgridConfig,
};

fn settings(model: SgsModel) -> SubgridConfig {
    SubgridConfig {
        model,
        filter_width_factor: 1.3,
        turbulent_prandtl: 0.8,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

#[test]
fn g1_local_bath_sgs_matches_all_six_full_tensors_and_variable_source_cp() {
    let grid = CartesianGrid::new(3, 4, 5, 0.2, 0.3, 0.4).unwrap();
    let (mu, bulk, conductivity) = (0.02, 0.003, 0.05);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let profile = |x: [f64; 3]| {
            let mut velocity = [0.0; 3];
            velocity[axis] = 0.01 + 0.02 * x[a] * x[a] + 0.03 * x[b];
            velocity[a] = 0.04 * x[a] + 0.015 * x[b];
            velocity[b] = -0.05 * x[b] + 0.025 * x[a];
            let density = 1.2 + 0.1 * x[a];
            ReservoirFaceState {
                density,
                velocity,
                pressure_pa: density * (2.8 + 0.05 * x[b]),
            }
        };
        let mut flow = field(gas(0.4), grid, OPEN, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
        let bytes = flow
            .required_reservoir_profile_workspace_bytes(&[side])
            .unwrap();
        flow.prepare_reservoir_profiles(&[side], bytes, |_, _, x| Ok(profile(x)), || false)
            .unwrap();
        for model in [
            SgsModel::Smagorinsky { cs: 0.17 },
            SgsModel::Wale { cw: 0.5 },
            SgsModel::Vreman { cv: 0.07 },
        ] {
            let config = settings(model);
            let make = || {
                MolecularTransport3d::with_sgs(
                    &flow,
                    law(flow.gas(), mu, bulk, conductivity),
                    config,
                    || false,
                )
                .unwrap()
                .with_reservoir_exchange()
            };
            let mut disabled = make();
            let mut active = make().with_local_sgs_reservoir_exchange().unwrap();
            assert!(!disabled.local_sgs_reservoir_exchange_enabled());
            assert!(active.local_sgs_reservoir_exchange_enabled());
            let base_rate = prepare(&mut disabled, &flow);
            let active_rate = prepare(&mut active, &flow);
            assert!(active_rate >= base_rate);
            assert_eq!(active.work.gradient, disabled.work.gradient);
            for face in 0..grid.face_count(side) {
                let x = grid.face_center(side, face).unwrap();
                let bath = profile(x);
                let distance = 0.5 * grid.spacing[axis];
                // Independent analytic tensor: the interior tangent trace is
                // constant and the bath trace contributes half its derivative.
                let mut g = [[0.0; 3]; 3];
                for i in 0..3 {
                    g[i][axis] = side.sign() * bath.velocity[i] / distance;
                }
                g[axis][a] = 0.02 * x[a];
                g[axis][b] = 0.015;
                g[a][a] = 0.02;
                g[a][b] = 0.0075;
                g[b][a] = 0.0125;
                g[b][b] = -0.025;
                let widths = grid.spacing.map(|h| h * config.filter_width_factor);
                let delta = grid.cell_volume().cbrt() * config.filter_width_factor;
                let nu = match model {
                    SgsModel::Vreman { cv } => fs_lbm::sgs3::vreman_eddy_viscosity(cv, widths, g),
                    model => fs_lbm::sgs3::eddy_viscosity(model, delta, g),
                }
                .unwrap();
                assert!(nu > 0.0);
                let mu_t = (0.5 + 0.5 * bath.density) * nu;
                let total_mu = mu + mu_t;
                let bath_temperature = bath.pressure_pa / bath.density;
                let cp = 3.5 + 0.4 * (1.0 + 0.5 * bath_temperature);
                let total_k = conductivity + mu_t * cp / config.turbulent_prandtl;
                let normal = ((4.0 / 3.0) * total_mu + bulk) * g[axis][axis]
                    + (bulk - (2.0 / 3.0) * total_mu) * (g[a][a] + g[b][b]);
                let shear = [
                    total_mu * (g[a][axis] + g[axis][a]),
                    total_mu * (g[b][axis] + g[axis][b]),
                ];
                let energy = -normal * bath.velocity[axis]
                    - shear[0] * bath.velocity[a]
                    - shear[1] * bath.velocity[b]
                    - total_k * side.sign() * (bath_temperature - 2.0) / distance;
                let mut c = [0; 3];
                c[axis] = if side.sign() < 0.0 {
                    0
                } else {
                    grid.shape[axis]
                };
                c[a] = face % grid.shape[a];
                c[b] = face / grid.shape[a];
                let index = grid.cell_face(axis, c);
                let flux = active.work.flux[axis][index];
                assert_eq!(flux.mass, 0.0);
                close(flux.normal_momentum, -normal, 5e-13);
                for t in 0..2 {
                    close(flux.transverse_momentum[t], -shear[t], 5e-13);
                }
                close(flux.energy, energy, 8e-13);
                assert_ne!(flux, disabled.work.flux[axis][index]);
            }
        }
    }
}

#[test]
fn g0_zero_bath_sgs_is_bit_identical_to_molecular_complete_steps() {
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let flow = heated_profile_flow();
        let source = law(flow.gas(), 0.02, 0.003, 0.05);
        let mut base = MolecularEuler3d::new(flow, source, || false)
            .unwrap()
            .with_reservoir_exchange();
        let flow = heated_profile_flow();
        let source = law(flow.gas(), 0.02, 0.003, 0.05);
        let mut active = MolecularEuler3d::with_sgs(flow, source, settings(model), || false)
            .unwrap()
            .with_local_sgs_reservoir_exchange()
            .unwrap();
        assert!(active.reservoir_exchange_enabled());
        for _ in 0..2 {
            let expected = base
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            let actual = active
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            assert_eq!(actual, expected);
            for (a, b) in active.flow().states().iter().zip(base.flow().states()) {
                assert_eq!(a.array().map(f64::to_bits), b.array().map(f64::to_bits));
            }
            for side in BoundarySide::ALL {
                assert_eq!(active.boundary_impulses(side), base.boundary_impulses(side));
            }
        }
    }
}

#[test]
fn g1_g4_bath_sgs_enters_the_diffusion_bound_and_refuses_excess_transport_atomically() {
    let grid = CartesianGrid::new(1, 1, 1, 0.2, 0.3, 0.4).unwrap();
    let field = field(gas(0.0), grid, OPEN, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
    let mut flow = Euler3d::new_with_reservoirs(
        field.gas.clone(),
        grid,
        field.states.clone(),
        OPEN,
        &[ReservoirBoundary {
            side: BoundarySide::XMin,
            density: 1.0,
            velocity: [0.0, 0.2, 0.0],
            pressure_pa: 3.0,
        }],
        || false,
    )
    .unwrap();
    let (mu, bulk, conductivity) = (0.02, 0.003, 0.05);
    let config = settings(SgsModel::Smagorinsky { cs: 0.17 });
    let mut op = MolecularTransport3d::with_sgs(
        &flow,
        law(flow.gas(), mu, bulk, conductivity),
        config,
        || false,
    )
    .unwrap()
    .with_local_sgs_reservoir_exchange()
    .unwrap();
    let rate = prepare(&mut op, &flow);
    let delta = grid.cell_volume().cbrt() * config.filter_width_factor;
    let mu_t = (0.17 * delta).powi(2) * (0.2 / (0.5 * grid.spacing[0]));
    let total_mu = mu + mu_t;
    let total_k = conductivity + mu_t * 3.5 / config.turbulent_prandtl;
    let stress_scale = (4.0 / 3.0) * total_mu
        + bulk
        + 2.0 * (bulk - (2.0 / 3.0) * total_mu).abs()
        + 4.0 * total_mu;
    let expected_rate = 4.0
        * stress_scale.max(total_k / 2.5)
        * grid
            .spacing
            .into_iter()
            .map(|h| h.recip().powi(2))
            .sum::<f64>();
    close(rate, expected_rate, 5e-13);
    let states = flow.states().to_vec();
    let mut limited = MolecularTransport3d::with_sgs(
        &flow,
        law(flow.gas(), mu, bulk, conductivity),
        SubgridConfig {
            maximum_eddy_to_molecular: 0.0,
            ..config
        },
        || false,
    )
    .unwrap()
    .with_local_sgs_reservoir_exchange()
    .unwrap();
    assert!(matches!(
        limited.advance(&mut flow, Time::new(0.001), MolecularControl::default(), || false),
        Err(MolecularError::SubgridViscosityRatio { axis: 0, face: 0, maximum: 0.0, observed })
            if observed > 0.0
    ));
    assert_eq!(flow.states(), states);
    assert_eq!(flow.time().value(), 0.0);
    for side in BoundarySide::ALL {
        assert!(
            limited
                .boundary_impulses(side)
                .iter()
                .all(|q| *q == Totals::default())
        );
    }
}

#[test]
fn g0_reservoir_sgs_requires_local_authorship_and_refuses_dynamic_histories() {
    let flow = cube(3, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
    assert!(
        MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.0, 0.05), || false)
            .unwrap()
            .with_local_sgs_reservoir_exchange()
            .is_err()
    );
    let dynamic = MolecularTransport3d::with_dynamic_sgs(
        &flow,
        law(flow.gas(), 0.02, 0.0, 0.05),
        DynamicSubgridConfig {
            model: DynamicSgsConfig {
                averaging: DynamicAveraging::Domain,
                maximum_coefficient_squared: 1.0,
                maximum_eddy_to_molecular: 100.0,
                usage: SgsUse::Exploration,
            },
            turbulent_prandtl: 0.8,
        },
        usize::MAX,
        || false,
    )
    .unwrap();
    assert!(dynamic.with_local_sgs_reservoir_exchange().is_err());
    let source = law(flow.gas(), 0.02, 0.0, 0.05);
    let memory = MolecularEuler3d::with_lagrangian_dynamic_sgs(
        flow,
        source,
        LagrangianSubgridConfig {
            model: LagrangianSgsConfig {
                memory_time: 0.01,
                maximum_departure_cells: 1.0,
                maximum_coefficient_squared: 1.0,
                maximum_eddy_to_molecular: 100.0,
                usage: SgsUse::Exploration,
            },
            turbulent_prandtl: 0.8,
        },
        usize::MAX,
        || false,
    )
    .unwrap();
    assert!(memory.with_local_sgs_reservoir_exchange().is_err());
}

#[test]
fn g3_g4_active_bath_sgs_conserves_and_replays_after_cancelled_complete_steps() {
    for stabilized in [false, true] {
        let make = || {
            let flow = heated_profile_flow();
            let source = law(flow.gas(), 0.02, 0.003, 0.05);
            MolecularEuler3d::with_sgs(
                flow,
                source,
                settings(SgsModel::Smagorinsky { cs: 0.17 }),
                || false,
            )
            .unwrap()
            .with_local_sgs_reservoir_exchange()
            .unwrap()
        };
        let run = |flow: &mut MolecularEuler3d<ConstantTransport>,
                   cancelled: &mut dyn FnMut() -> bool| {
            if stabilized {
                flow.step_stabilized(
                    Time::new(0.001),
                    MolecularStepControl::default(),
                    MolecularStabilization::default(),
                    cancelled,
                )
                .map(|r| r.transport)
            } else {
                flow.step(Time::new(0.001), MolecularStepControl::default(), cancelled)
            }
        };
        let mut reference = make();
        run(&mut reference, &mut || false).unwrap();
        let previous = reference.last_step().copied();
        reference = reference.with_local_sgs_reservoir_exchange().unwrap();
        assert_eq!(reference.last_step().copied(), previous);
        let mut successful_polls = 0;
        let expected = run(&mut reference, &mut || {
            successful_polls += 1;
            false
        })
        .unwrap();
        assert!(expected.first.boundary_change.energy > 0.0);
        assert!(expected.second.boundary_change.energy > 0.0);
        for residual in expected.conservation_residual.array() {
            close(residual, 0.0, 5e-12);
        }
        for cutoff in [successful_polls / 2, successful_polls] {
            let mut replay = make();
            run(&mut replay, &mut || false).unwrap();
            let states = replay.flow().states().to_vec();
            let time = replay.flow().time();
            let previous = replay.last_step().copied();
            let loads = BoundarySide::ALL.map(|s| replay.boundary_impulses(s).to_vec());
            let traces = BoundarySide::ALL.map(|s| replay.flow().reservoir_face(s, 0));
            let mut polls = 0;
            assert!(
                run(&mut replay, &mut || {
                    polls += 1;
                    polls == cutoff
                })
                .is_err()
            );
            assert_eq!(replay.flow().states(), states);
            assert_eq!(replay.flow().time(), time);
            assert_eq!(replay.last_step().copied(), previous);
            for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
                assert_eq!(replay.boundary_impulses(side), loads[i]);
                assert_eq!(replay.flow().reservoir_face(side, 0), traces[i]);
            }
            assert_eq!(run(&mut replay, &mut || false).unwrap(), expected);
            assert_eq!(replay.flow().states(), reference.flow().states());
        }
    }
}
