use super::*;
use crate::compressible::mixture::tests::{close, gas, source};
use crate::compressible::spatial::mixture::state_into;
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;

const PERIODIC: Boundaries3d = Boundaries3d {
    x: AxisBoundaries::Periodic,
    y: AxisBoundaries::Periodic,
    z: AxisBoundaries::Periodic,
};
const WALL: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Reflecting,
    upper: BoundaryCondition::Reflecting,
};
const WALLS: Boundaries3d = Boundaries3d {
    x: WALL,
    y: WALL,
    z: WALL,
};
const A: [[f64; 3]; 3] = [[0.15, -0.1, 0.05], [0.2, 0.1, -0.3], [-0.05, 0.25, 0.12]];

fn config(model: SgsModel) -> MixtureSubgridConfig {
    MixtureSubgridConfig {
        transport: SubgridConfig {
            model,
            filter_width_factor: 1.0,
            turbulent_prandtl: 0.7,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        turbulent_schmidt: 0.8,
    }
}

fn material(formation: f64) -> MixtureGas {
    gas(vec![
        source("A", R, 3.5, 0.2, 30.0 + formation),
        source("B", 2.0 * R, 4.5, 0.3, -60.0 - 3.0 * formation),
        source("absent-C", 0.5 * R, 4.5, 0.2, 7.0 + 2.0 * formation),
    ])
}

fn field(
    gas: MixtureGas,
    grid: CartesianGrid,
    boundaries: Boundaries3d,
    initial: impl Fn([f64; 3]) -> ([f64; 3], [f64; 3], f64),
) -> MixtureEuler3d {
    let mut states = vec![0.0; grid.len() * 7];
    for (i, row) in states.chunks_exact_mut(7).enumerate() {
        let (partial, velocity, temperature) = initial(grid.cell_center(i).unwrap());
        let pressure = temperature * (partial[0] + 0.5 * partial[1] + 2.0 * partial[2]);
        state_into(&gas, &partial, velocity, pressure, row).unwrap();
    }
    let bytes = MixtureEuler3d::required_workspace_bytes(&gas, &grid).unwrap();
    MixtureEuler3d::new(gas, grid, states, boundaries, bytes, || false).unwrap()
}

fn original(flow: &MixtureEuler3d) -> SpatialTransport<ConstantTransport> {
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "manufactured-mixture-sgs",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(0.02),
            bulk_viscosity: DynViscosity::new(0.003),
            thermal_conductivity: ThermalConductivityQuantityV1::new(0.03),
        },
        vec![0.0, 0.01, 0.012, 0.01, 0.0, 0.013, 0.012, 0.013, 0.0],
    )
    .unwrap();
    SpatialTransport::new(flow.gas().clone(), law, || false).unwrap()
}

fn operator(
    flow: &MixtureEuler3d,
    config: Option<MixtureSubgridConfig>,
) -> MolecularTransport3d<ConstantTransport> {
    let source = original(flow);
    if let Some(config) = config {
        let bytes = MolecularTransport3d::required_workspace_bytes_with_sgs(flow, &source).unwrap();
        let result = MolecularTransport3d::with_sgs(flow, source, config, bytes, || false).unwrap();
        assert_eq!(result.workspace_bytes() + flow.workspace_bytes(), bytes);
        result
    } else {
        let bytes = MolecularTransport3d::required_workspace_bytes(flow, &source).unwrap();
        MolecularTransport3d::new(flow, source, bytes, || false).unwrap()
    }
}

fn faces(op: &mut MolecularTransport3d<ConstantTransport>, flow: &MixtureEuler3d) -> Stage {
    super::super::spatial(
        &mut op.source,
        &op.grid,
        op.boundaries,
        &op.walls,
        flow.states(),
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.flux,
        op.sgs.as_mut(),
        DissipationControl::default(),
        &mut || false,
    )
    .unwrap()
}

fn linear() -> MixtureEuler3d {
    field(
        material(0.0),
        CartesianGrid::new(4, 4, 4, 0.2, 0.3, 0.4).unwrap(),
        WALLS,
        |x| {
            let first = 0.4 + 0.05 * (x[0] + x[1] + x[2]);
            let velocity = std::array::from_fn(|i| {
                [0.03, -0.07, 0.04][i] + (0..3).map(|j| A[i][j] * x[j]).sum::<f64>()
            });
            (
                [first, 1.0 - first, 0.0],
                velocity,
                1.0 + 0.1 * (x[0] + x[1] + x[2]),
            )
        },
    )
}

fn periodic(formation: f64) -> MixtureEuler3d {
    field(
        material(formation),
        CartesianGrid::new(8, 2, 2, 0.125, 0.3, 0.4).unwrap(),
        PERIODIC,
        |x| {
            let angle = 2.0 * std::f64::consts::PI * x[0];
            let first = 0.5 + 0.2 * angle.sin();
            (
                [first, 1.0 - first, 0.0],
                [0.0, 0.1 * angle.sin(), 0.1 * angle.cos()],
                1.0 + 0.05 * angle.cos(),
            )
        },
    )
}

fn balance(report: &MolecularReport) {
    assert!(report.maximum_parabolic_number <= 0.40000000000001);
    assert!(report.maximum_mass_frame_residual_kg_m2_s < 1e-12);
    for r in report.conservation_residual.species_masses {
        close(r, 0.0, 2e-12);
    }
    for r in report.conservation_residual.momentum {
        close(r, 0.0, 2e-12);
    }
    close(report.conservation_residual.energy, 0.0, 2e-12);
    close(report.absolute_energy_residual_j, 0.0, 2e-10);
}

#[test]
fn sgs_mixture_faces_match_full_deviatoric_stress_heat_species_and_enthalpy() {
    let flow = linear();
    let cfg = config(SgsModel::Smagorinsky { cs: 0.18 });
    let mut modeled = operator(&flow, Some(cfg));
    let mut molecular = operator(&flow, None);
    let stage = faces(&mut modeled, &flow);
    let baseline = faces(&mut molecular, &flow);
    assert!(stage.rate > baseline.rate);
    let strain = (2.0
        * (0..3)
            .flat_map(|i| (0..3).map(move |j| (0.5 * (A[i][j] + A[j][i])).powi(2)))
            .sum::<f64>())
    .sqrt();
    let mu = (0.18 * flow.grid.cell_volume().cbrt()).powi(2) * strain;
    let left = flow.grid.index([1, 1, 1]);
    close(
        modeled.sgs.as_ref().unwrap().dynamic_viscosity[left],
        mu,
        2e-15,
    );
    for axis in 0..3 {
        let mut c = [1, 1, 1];
        c[axis] += 1;
        let right = flow.grid.index(c);
        let face = flow.grid.cell_face(axis, c);
        let delta: [f64; 7] = std::array::from_fn(|d| {
            modeled.work.flux[axis][7 * face + d] - molecular.work.flux[axis][7 * face + d]
        });
        let pl = modeled.work.primitive[left];
        let pr = modeled.work.primitive[right];
        let t = 0.5 * pl.temperature + 0.5 * pr.temperature;
        let first = 0.5 * flow.cell(left).unwrap().partial_densities[0]
            + 0.5 * flow.cell(right).unwrap().partial_densities[0];
        let cp = first * (3.5 + 0.2 * t) + (1.0 - first) * 0.5 * (4.5 + 0.3 * t);
        let j = -mu / cfg.turbulent_schmidt * 0.05;
        close(delta[0], j, 2e-15);
        close(delta[1], -j, 2e-15);
        assert_eq!(delta[2], 0.0);
        let h = |c: f64, linear: f64, inverse_mass: f64| {
            inverse_mass
                * ((c - 1.0) * (t - 0.001) + 0.5 * linear * (t * t - 0.001_f64.powi(2)) + t)
        };
        let mut energy = -mu * cp / cfg.transport.turbulent_prandtl * 0.1
            + j * (h(3.5, 0.2, 1.0) - h(4.5, 0.3, 0.5));
        for component in 0..3 {
            let stress = mu
                * (A[component][axis] + A[axis][component]
                    - if component == axis {
                        (2.0 / 3.0) * 0.37
                    } else {
                        0.0
                    });
            close(delta[3 + component], -stress, 2e-15);
            energy -= 0.5 * (pl.velocity[component] + pr.velocity[component]) * stress;
        }
        close(delta[6], energy, 3e-15);
    }
}

#[test]
fn sgs_mixture_preserves_directional_vreman_metric_and_all_shared_models() {
    let flow = linear();
    for model in [
        SgsModel::Smagorinsky { cs: 0.17 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let mut op = operator(&flow, Some(config(model)));
        faces(&mut op, &flow);
        let expected = match model {
            SgsModel::Vreman { cv } => vreman_eddy_viscosity(cv, flow.grid.spacing, A),
            model => eddy_viscosity(model, flow.grid.cell_volume().cbrt(), A),
        }
        .unwrap();
        assert!(expected > 0.0);
        close(
            op.sgs.as_ref().unwrap().dynamic_viscosity[flow.grid.index([1, 1, 1])],
            expected,
            2e-15,
        );
    }
    let transpose = std::array::from_fn(|i| std::array::from_fn(|j| A[j][i]));
    let right = vreman_eddy_viscosity(0.07, flow.grid.spacing, A).unwrap();
    let wrong = vreman_eddy_viscosity(0.07, flow.grid.spacing, transpose).unwrap();
    assert!(
        (right - wrong).abs() > 1e-5,
        "anisotropic fixture must expose transposition"
    );
}

#[test]
fn sgs_mixture_prandtl_and_schmidt_control_distinct_conservative_fluxes() {
    let flow = linear();
    let mut cfg = config(SgsModel::Smagorinsky { cs: 0.18 });
    let mut base = operator(&flow, None);
    faces(&mut base, &flow);
    let mut standard = operator(&flow, Some(cfg));
    faces(&mut standard, &flow);
    cfg.turbulent_schmidt *= 0.5;
    let mut mixing = operator(&flow, Some(cfg));
    faces(&mut mixing, &flow);
    cfg.turbulent_schmidt *= 2.0;
    cfg.transport.turbulent_prandtl *= 0.5;
    let mut heating = operator(&flow, Some(cfg));
    faces(&mut heating, &flow);
    let index = flow.grid.cell_face(0, [2, 1, 1]) * 7;
    let s = &standard.work.flux[0][index..index + 7];
    let b = &base.work.flux[0][index..index + 7];
    let m = &mixing.work.flux[0][index..index + 7];
    let h = &heating.work.flux[0][index..index + 7];
    close(m[0] - b[0], 2.0 * (s[0] - b[0]), 2e-15);
    assert_eq!(&m[3..6], &s[3..6]);
    assert_eq!(&h[..6], &s[..6]);
    assert_ne!(m[6], s[6]);
    assert_ne!(h[6], s[6]);
}

#[test]
fn sgs_mixture_zero_coefficients_reproduce_molecular_states_and_reports_exactly() {
    let mut reference = periodic(0.0);
    let mut molecular = operator(&reference, None);
    let report = molecular
        .advance(
            &mut reference,
            Time::new(0.001),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let mut flow = periodic(0.0);
        let mut op = operator(&flow, Some(config(model)));
        let result = op
            .advance(
                &mut flow,
                Time::new(0.001),
                DissipationControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(result, report);
        assert_eq!(flow.states(), reference.states());
        assert_eq!(op.boundary_impulses, molecular.boundary_impulses);
    }
}

#[test]
fn sgs_mixture_periodic_update_conserves_inventory_and_original_energy_gauge() {
    let cfg = config(SgsModel::Smagorinsky { cs: 0.25 });
    let mut flow = periodic(0.0);
    let mut shifted = periodic(1e12);
    let mut baseline = periodic(0.0);
    let mut op = operator(&flow, Some(cfg));
    let mut shifted_op = operator(&shifted, Some(cfg));
    let mut baseline_op = operator(&baseline, None);
    let before = flow.states().to_vec();
    let bytes = op.workspace_bytes() + flow.workspace_bytes();
    let report = op
        .advance(
            &mut flow,
            Time::new(0.003),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    balance(&report);
    shifted_op
        .advance(
            &mut shifted,
            Time::new(0.003),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    baseline_op
        .advance(
            &mut baseline,
            Time::new(0.003),
            DissipationControl::default(),
            || false,
        )
        .unwrap();
    assert_ne!(flow.states(), before);
    assert_ne!(flow.states(), baseline.states());
    assert_eq!(
        flow.states(),
        shifted.states(),
        "species formation offsets cancel within each original source"
    );
    for i in 0..flow.cell_count() {
        assert_eq!(flow.cell(i).unwrap().partial_densities[2], 0.0);
    }
    assert_eq!(op.workspace_bytes() + flow.workspace_bytes(), bytes);
    assert_eq!(flow.time().value(), 0.0);
}

#[test]
fn sgs_mixture_no_slip_walls_keep_molecular_heat_stress_and_impermeability() {
    let flow = field(
        material(0.0),
        CartesianGrid::new(4, 4, 4, 0.2, 0.3, 0.4).unwrap(),
        WALLS,
        |x| {
            // At the stationary z=0 wall, u_x=a*x*z has zero physical
            // tangential divergence despite nonzero cell-center d_x(u_x).
            (
                [0.6, 0.4, 0.0],
                [0.2 * x[0] * x[2], 0.0, 0.0],
                1.0 + 0.1 * x[2],
            )
        },
    );
    let mut molecular = operator(&flow, None);
    let mut modeled = operator(&flow, Some(config(SgsModel::Smagorinsky { cs: 0.2 })));
    for side in BoundarySide::ALL {
        let wall = MolecularWall3d::NoSlip {
            velocity: [0.0; 3],
            thermal: WallThermalCondition::Isothermal(Temperature::new(0.8)),
        };
        molecular = molecular.with_wall(side, wall).unwrap();
        modeled = modeled.with_wall(side, wall).unwrap();
    }
    faces(&mut molecular, &flow);
    faces(&mut modeled, &flow);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let position = if side.sign() < 0.0 {
            0
        } else {
            flow.grid.shape[axis]
        };
        for line in 0..flow.grid.line_count(axis) {
            let face = flow.grid.boundary_face(axis, position, line);
            let range = 7 * face..7 * (face + 1);
            assert_eq!(
                &modeled.work.flux[axis][range.clone()],
                &molecular.work.flux[axis][range]
            );
            assert_eq!(&modeled.work.flux[axis][7 * face..7 * face + 3], &[0.0; 3]);
        }
    }
    let c = [1, 1, 0];
    let x = flow.grid.cell_center(flow.grid.index(c)).unwrap();
    let offset = 7 * flow.grid.cell_face(2, c);
    let lower_z = &modeled.work.flux[2][offset..offset + 7];
    assert_eq!(
        lower_z[5], 0.0,
        "no fictitious normal traction from interior tangential divergence"
    );
    close(lower_z[3], -0.02 * 0.2 * x[0], 2e-15);
    close(lower_z[6], -0.03 * (1.0 + 0.1 * x[2] - 0.8) / x[2], 2e-15);
    assert_ne!(
        modeled.work.flux, molecular.work.flux,
        "interior SGS remains active"
    );
}

#[test]
fn sgs_mixture_total_diffusion_bounds_include_small_schmidt_and_atomic_refusals() {
    let mut flow = periodic(0.0);
    let mut cfg = config(SgsModel::Smagorinsky { cs: 0.25 });
    let mut ordinary = operator(&flow, Some(cfg));
    let ordinary_rate = faces(&mut ordinary, &flow).rate;
    cfg.turbulent_schmidt = 0.002;
    let mut op = operator(&flow, Some(cfg));
    assert!(faces(&mut op, &flow).rate > 5.0 * ordinary_rate);
    let mut polls = 0;
    let report = op
        .advance(
            &mut flow,
            Time::new(0.003),
            DissipationControl::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    balance(&report);
    assert!(report.accepted_steps > 1);
    let state = flow.states().to_vec();
    let loads = op.boundary_impulses.clone();
    let mut count = 0;
    assert_eq!(
        op.advance(
            &mut flow,
            Time::new(0.003),
            DissipationControl::default(),
            || {
                count += 1;
                count == polls / 2
            }
        ),
        Err(DissipationError::Cancelled)
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.boundary_impulses, loads);
    assert_eq!(op.last_step(), Some(&report));
    op.sgs
        .as_mut()
        .unwrap()
        .config
        .transport
        .maximum_eddy_to_molecular = 0.0;
    assert!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            DissipationControl::default(),
            || false
        )
        .is_err()
    );
    assert_eq!(flow.states(), state);
    assert_eq!(op.last_step(), Some(&report));
    let source = original(&flow);
    let bytes = MolecularTransport3d::required_workspace_bytes_with_sgs(&flow, &source).unwrap();
    assert!(MolecularTransport3d::with_sgs(&flow, source, cfg, bytes - 1, || false).is_err());
    for invalid in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        cfg.turbulent_schmidt = invalid;
        assert!(cfg.validate().is_err());
    }
}
