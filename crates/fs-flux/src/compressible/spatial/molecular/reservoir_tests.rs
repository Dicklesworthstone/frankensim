//! Physical bath effort, stencil orientation and native atomic integration.

use super::*;
use crate::compressible::spatial::{
    BodyForce3d, ForceControl3d, ReservoirBoundary, ReservoirFaceState,
};
use fs_gas::multicomponent::dissipation::DissipationError;

const OPEN_AXIS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const OPEN: Boundaries3d = Boundaries3d {
    x: OPEN_AXIS,
    y: OPEN_AXIS,
    z: OPEN_AXIS,
};

fn prepare<L: TransportLaw>(op: &mut MolecularTransport3d<L>, flow: &Euler3d) -> f64 {
    spatial(
        &op.gas,
        &op.grid,
        op.boundaries,
        &op.walls,
        op.reservoir_exchange.then(|| flow.reservoirs.accepted()),
        &flow.states,
        &mut op.law,
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.viscosity,
        &mut op.work.flux,
        &mut op.work.binary,
        #[cfg(feature = "experimental-sgs")]
        op.sgs.as_mut(),
        &mut || false,
    )
    .unwrap()
}

struct ObservedLaw {
    original: ConstantTransport,
    states: Vec<[f64; 2]>,
    scale_with_temperature: bool,
}
impl TransportLaw for ObservedLaw {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.original.source()
    }
    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: Pressure,
        fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        self.states.push([temperature.value(), pressure.value()]);
        let value = self
            .original
            .evaluate(temperature, pressure, fractions, binary, cancelled)?;
        if self.scale_with_temperature {
            Ok(TransportCoefficients {
                dynamic_viscosity: DynViscosity::new(
                    value.dynamic_viscosity.value() * temperature.value(),
                ),
                bulk_viscosity: DynViscosity::new(
                    value.bulk_viscosity.value() * temperature.value(),
                ),
                thermal_conductivity: ThermalConductivityQuantityV1::new(
                    value.thermal_conductivity.value() * temperature.value(),
                ),
            })
        } else {
            Ok(value)
        }
    }
}

#[test]
fn g1_six_sides_use_half_cell_full_stress_and_heat_at_prescribed_bath_work() {
    let h = [0.2, 0.3, 0.4];
    let grid = CartesianGrid::new(3, 3, 3, h[0], h[1], h[2]).unwrap();
    let interior = [0.1, -0.2, 0.05];
    let bath = ReservoirFaceState {
        density: 1.2,
        velocity: [0.3, -0.15, -0.1],
        pressure_pa: 3.6,
    };
    let (mu, bulk, conductivity) = (0.02, 0.003, 0.05);
    for side in BoundarySide::ALL {
        for profile in [false, true] {
            let mut flow = field(gas(0.0), grid, OPEN, |_| {
                [1.0, interior[0], interior[1], interior[2], 2.2]
            });
            if profile {
                let bytes = flow
                    .required_reservoir_profile_workspace_bytes(&[side])
                    .unwrap();
                flow.prepare_reservoir_profiles(&[side], bytes, |_, _, _| Ok(bath), || false)
                    .unwrap();
            } else {
                flow = Euler3d::new_with_reservoirs(
                    flow.gas.clone(),
                    grid,
                    flow.states.clone(),
                    OPEN,
                    &[ReservoirBoundary {
                        side,
                        density: bath.density,
                        velocity: bath.velocity,
                        pressure_pa: bath.pressure_pa,
                    }],
                    || false,
                )
                .unwrap();
            }
            let observed = ObservedLaw {
                original: law(flow.gas(), mu, bulk, conductivity),
                states: Vec::new(),
                scale_with_temperature: false,
            };
            let mut op = MolecularTransport3d::new(&flow, observed, || false).unwrap();
            let axis = side.axis();
            let [a, b] = tangents(axis);
            let mut c = [1; 3];
            c[axis] = if side.sign() < 0.0 {
                0
            } else {
                grid.shape[axis]
            };
            let index = grid.cell_face(axis, c);
            prepare(&mut op, &flow);
            assert!(!op.reservoir_exchange_enabled());
            assert_eq!(op.work.flux[axis][index], SpatialFaceFlux::default());
            op = op.with_reservoir_exchange();
            assert!(prepare(&mut op, &flow) > 0.0);
            let distance = 0.5 * h[axis];
            let derivative: [f64; 3] =
                std::array::from_fn(|d| side.sign() * (bath.velocity[d] - interior[d]) / distance);
            let normal = ((4.0 / 3.0) * mu + bulk) * derivative[axis];
            let shear = [mu * derivative[a], mu * derivative[b]];
            let heat = -conductivity * side.sign() * (3.0 - 2.2) / distance;
            let flux = op.work.flux[axis][index];
            close(flux.mass, 0.0, 0.0);
            close(flux.normal_momentum, -normal, 2e-14);
            for t in 0..2 {
                close(flux.transverse_momentum[t], -shear[t], 2e-14);
            }
            close(
                flux.energy,
                -normal * bath.velocity[axis]
                    - shear[0] * bath.velocity[a]
                    - shear[1] * bath.velocity[b]
                    + heat,
                3e-14,
            );
            assert!(
                op.law
                    .states
                    .iter()
                    .any(|q| { (q[0] - 2.6).abs() < 2e-13 && (q[1] - 2.9).abs() < 2e-13 }),
                "original law must see the source-admitted midpoint T,p"
            );
            c[axis] = if side.sign() < 0.0 {
                0
            } else {
                grid.shape[axis] - 1
            };
            for d in 0..3 {
                close(
                    op.work.gradient[grid.index(c)][d][axis],
                    0.5 * derivative[d],
                    2e-14,
                );
            }
        }
    }
}

#[test]
fn g1_hot_bath_coefficients_enter_the_same_conservative_half_cell_rate() {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let mut flow = field(gas(0.0), grid, OPEN, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
    let side = BoundarySide::XMin;
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&[side])
        .unwrap();
    flow.prepare_reservoir_profiles(
        &[side],
        bytes,
        |_, _, _| {
            Ok(ReservoirFaceState {
                density: 1.0,
                velocity: [0.0; 3],
                pressure_pa: 40.0,
            })
        },
        || false,
    )
    .unwrap();
    let source = ObservedLaw {
        original: law(flow.gas(), 0.02, 0.003, 0.05),
        states: Vec::new(),
        scale_with_temperature: true,
    };
    let mut op = MolecularTransport3d::new(&flow, source, || false).unwrap();
    let legacy_rate = prepare(&mut op, &flow);
    op = op.with_reservoir_exchange();
    let rate = prepare(&mut op, &flow);
    let (mu, bulk, conductivity): (f64, f64, f64) = (0.02 * 21.0, 0.003 * 21.0, 0.05 * 21.0);
    let diffusivity = ((4.0 / 3.0) * mu + bulk + 2.0 * (bulk - (2.0 / 3.0) * mu).abs() + 4.0 * mu)
        .max(conductivity / 2.5);
    let expected = 4.0
        * diffusivity
        * grid
            .spacing
            .into_iter()
            .map(|h| h.recip().powi(2))
            .sum::<f64>();
    close(rate, expected, 3e-12);
    close(rate / legacy_rate, 10.5, 3e-13);
    assert!(rate > legacy_rate);
}

#[test]
fn g1_profiles_supply_tangential_cross_derivatives_on_all_six_faces() {
    let grid = CartesianGrid::new(3, 4, 5, 0.2, 0.3, 0.4).unwrap();
    let (mu, bulk) = (0.02, 0.003);
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let profile = |x: [f64; 3]| {
            let mut velocity = [0.0; 3];
            velocity[axis] = 0.02 * x[a] * x[a] + 0.03 * x[b];
            velocity[a] = 0.04 * x[a];
            velocity[b] = -0.05 * x[b];
            ReservoirFaceState {
                density: 1.0,
                velocity,
                pressure_pa: 2.0,
            }
        };
        let mut flow = field(gas(0.0), grid, OPEN, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
        let bytes = flow
            .required_reservoir_profile_workspace_bytes(&[side])
            .unwrap();
        flow.prepare_reservoir_profiles(&[side], bytes, |_, _, x| Ok(profile(x)), || false)
            .unwrap();
        let mut op = MolecularTransport3d::new(&flow, law(flow.gas(), mu, bulk, 0.0), || false)
            .unwrap()
            .with_reservoir_exchange();
        prepare(&mut op, &flow);
        for face in 0..grid.face_count(side) {
            let x = grid.face_center(side, face).unwrap();
            let velocity = profile(x).velocity;
            let derivative: [f64; 3] =
                velocity.map(|u| side.sign() * u / (0.5 * grid.spacing[axis]));
            let normal = ((4.0 / 3.0) * mu + bulk) * derivative[axis]
                + (bulk - (2.0 / 3.0) * mu) * 0.5 * (0.04 - 0.05);
            let shear = [
                mu * (derivative[a] + 0.02 * x[a]),
                mu * (derivative[b] + 0.015),
            ];
            let mut c = [0; 3];
            c[axis] = if side.sign() < 0.0 {
                0
            } else {
                grid.shape[axis]
            };
            c[a] = face % grid.shape[a];
            c[b] = face / grid.shape[a];
            let flux = op.work.flux[axis][grid.cell_face(axis, c)];
            close(flux.normal_momentum, -normal, 3e-14);
            for t in 0..2 {
                close(flux.transverse_momentum[t], -shear[t], 3e-14);
            }
            close(
                flux.energy,
                -normal * velocity[axis] - shear[0] * velocity[a] - shear[1] * velocity[b],
                3e-14,
            );
        }
    }
}

fn heated_profile_flow() -> Euler3d {
    let grid = CartesianGrid::new(3, 2, 2, 0.2, 0.3, 0.4).unwrap();
    let mut flow = field(gas(0.0), grid, OPEN, |_| [1.0, 0.0, 0.0, 0.0, 2.0]);
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&BoundarySide::ALL)
        .unwrap();
    flow.prepare_reservoir_profiles(
        &BoundarySide::ALL,
        bytes,
        |side, face, _| {
            Ok(ReservoirFaceState {
                density: 1.0,
                velocity: [0.04, -0.02, 0.01],
                pressure_pa: 3.0 + 0.01 * face as f64 + 0.02 * side.axis() as f64,
            })
        },
        || false,
    )
    .unwrap();
    flow
}

#[test]
fn g0_equilibrium_and_unbound_sides_retain_the_original_molecular_trajectory() {
    let grid = CartesianGrid::new(2, 3, 1, 0.2, 0.3, 0.4).unwrap();
    let velocity = [0.125, -0.0625, 0.03125];
    for bound in [false, true] {
        for stabilized in [false, true] {
            let mut flow = field(gas(0.4), grid, OPEN, |_| {
                [1.0, velocity[0], velocity[1], velocity[2], 2.0]
            });
            if bound {
                let baths = BoundarySide::ALL.map(|side| ReservoirBoundary {
                    side,
                    density: 1.0,
                    velocity,
                    pressure_pa: 2.0,
                });
                flow = Euler3d::new_with_reservoirs(
                    flow.gas.clone(),
                    grid,
                    flow.states.clone(),
                    OPEN,
                    &baths,
                    || false,
                )
                .unwrap();
            }
            let before = flow.states().to_vec();
            let mut op =
                MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.05), || false)
                    .unwrap()
                    .with_reservoir_exchange();
            let report = if stabilized {
                op.advance_stabilized(
                    &mut flow,
                    Time::new(0.01),
                    MolecularControl::default(),
                    MolecularStabilization::default(),
                    || false,
                )
                .unwrap()
                .transport
            } else {
                op.advance(
                    &mut flow,
                    Time::new(0.01),
                    MolecularControl::default(),
                    || false,
                )
                .unwrap()
            };
            assert_eq!(flow.states(), before);
            assert_eq!(report.boundary_change, Totals::default());
            assert_eq!(flow.time().value(), 0.0);
        }
    }
}

#[test]
fn g3_g4_ssp_and_rkl_bath_exchange_conserve_and_cancel_without_partial_publication() {
    for stabilized in [false, true] {
        let make = || {
            let flow = heated_profile_flow();
            let op = MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.5), || false)
                .unwrap()
                .with_reservoir_exchange();
            (flow, op)
        };
        let run = |flow: &mut Euler3d,
                   op: &mut MolecularTransport3d<ConstantTransport>,
                   cancelled: &mut dyn FnMut() -> bool| {
            if stabilized {
                op.advance_stabilized(
                    flow,
                    Time::new(0.1),
                    MolecularControl::default(),
                    MolecularStabilization::default(),
                    cancelled,
                )
                .map(|r| r.transport)
            } else {
                op.advance(flow, Time::new(0.1), MolecularControl::default(), cancelled)
            }
        };
        let (mut reference, mut reference_op) = make();
        let report = run(&mut reference, &mut reference_op, &mut || false).unwrap();
        assert!(report.after.energy > report.before.energy);
        close(report.boundary_change.mass, 0.0, 0.0);
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 5e-12);
        }
        let outward: f64 = BoundarySide::ALL
            .into_iter()
            .flat_map(|s| reference_op.boundary_impulses(s))
            .map(|q| q.energy)
            .sum();
        close(outward + report.boundary_change.energy, 0.0, 3e-12);
        let mut successful_polls = 0;
        let expected = run(&mut reference, &mut reference_op, &mut || {
            successful_polls += 1;
            false
        })
        .unwrap();
        let (mut flow, mut op) = make();
        run(&mut flow, &mut op, &mut || false).unwrap();
        let states = flow.states().to_vec();
        let bytes = flow.reservoir_profile_workspace_bytes();
        let loads = BoundarySide::ALL.map(|s| op.boundary_impulses(s).to_vec());
        let traces = BoundarySide::ALL.map(|s| flow.reservoir_face(s, 0));
        let mut polls = 0;
        assert!(
            run(&mut flow, &mut op, &mut || {
                polls += 1;
                polls == successful_polls
            })
            .is_err()
        );
        assert_eq!(flow.states(), states);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(flow.reservoir_profile_workspace_bytes(), bytes);
        for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
            assert_eq!(op.boundary_impulses(side), loads[i]);
            assert_eq!(flow.reservoir_face(side, 0), traces[i]);
        }
        if stabilized {
            assert!(
                op.advance_stabilized(
                    &mut flow,
                    Time::new(0.1),
                    MolecularControl::default(),
                    MolecularStabilization {
                        maximum_operator_evaluations: 2,
                        ..MolecularStabilization::default()
                    },
                    || false,
                )
                .is_err()
            );
            assert_eq!(flow.states(), states);
        }
        let replay = run(&mut flow, &mut op, &mut || false).unwrap();
        assert_eq!(replay, expected);
        assert_eq!(flow.states(), reference.states());
    }
}

#[test]
fn g3_coupled_forced_and_stabilized_splits_keep_the_same_held_bath_account() {
    for stabilized in [false, true] {
        for forced in [false, true] {
            let flow = heated_profile_flow();
            let source = law(flow.gas(), 0.02, 0.003, 0.05);
            let mut coupled = MolecularEuler3d::new(flow, source, || false)
                .unwrap()
                .with_reservoir_exchange();
            let traces = BoundarySide::ALL.map(|s| coupled.flow().reservoir_face(s, 0));
            let mut work = coupled.flow().forcing_workspace(40 * 12).unwrap();
            let force = BodyForce3d::UniformAcceleration([0.03, -0.02, 0.01]);
            let cap = ForceControl3d { maximum_mach: 1.0 };
            let dt = Time::new(0.001);
            for _ in 0..2 {
                let control = MolecularStepControl::default();
                let report = match (stabilized, forced) {
                    (false, false) => coupled.step(dt, control, || false).unwrap(),
                    (false, true) => coupled
                        .step_with_forcing(dt, control, force, cap, &mut work, || false)
                        .unwrap(),
                    (true, false) => {
                        coupled
                            .step_stabilized(dt, control, MolecularStabilization::default(), || {
                                false
                            })
                            .unwrap()
                            .transport
                    }
                    (true, true) => {
                        coupled
                            .step_stabilized_with_forcing(
                                dt,
                                control,
                                MolecularStabilization::default(),
                                force,
                                cap,
                                &mut work,
                                || false,
                            )
                            .unwrap()
                            .transport
                    }
                };
                assert!(coupled.reservoir_exchange_enabled());
                assert!(report.first.boundary_change.energy > 0.0);
                assert!(report.second.boundary_change.energy > 0.0);
                for residual in report.conservation_residual.array() {
                    close(residual, 0.0, 5e-12);
                }
                for (i, side) in BoundarySide::ALL.into_iter().enumerate() {
                    assert_eq!(coupled.flow().reservoir_face(side, 0), traces[i]);
                }
            }
        }
    }
}

#[cfg(feature = "experimental-sgs")]
#[test]
fn g1_local_sgs_uses_bath_ghosts_but_leaves_the_bath_flux_molecular() {
    use crate::compressible::spatial::molecular::sgs::{SgsModel, SgsUse, SubgridConfig};
    let flow = heated_profile_flow();
    let mut base = MolecularTransport3d::new(&flow, law(flow.gas(), 0.02, 0.003, 0.05), || false)
        .unwrap()
        .with_reservoir_exchange();
    let mut active = MolecularTransport3d::with_sgs(
        &flow,
        law(flow.gas(), 0.02, 0.003, 0.05),
        SubgridConfig {
            model: SgsModel::Smagorinsky { cs: 0.17 },
            filter_width_factor: 1.0,
            turbulent_prandtl: 0.8,
            maximum_eddy_to_molecular: 20.0,
            usage: SgsUse::Exploration,
        },
        || false,
    )
    .unwrap()
    .with_reservoir_exchange();
    prepare(&mut base, &flow);
    prepare(&mut active, &flow);
    assert_eq!(base.work.gradient, active.work.gradient);
    assert!(
        active
            .work
            .gradient
            .iter()
            .flatten()
            .flatten()
            .any(|&g| g != 0.0)
    );
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        for face in 0..flow.grid.face_count(side) {
            let mut c = [0; 3];
            c[axis] = if side.sign() < 0.0 {
                0
            } else {
                flow.grid.shape[axis]
            };
            c[a] = face % flow.grid.shape[a];
            c[b] = face / flow.grid.shape[a];
            let index = flow.grid.cell_face(axis, c);
            assert_eq!(base.work.flux[axis][index], active.work.flux[axis][index]);
        }
    }
    assert!(
        base.work
            .flux
            .iter()
            .zip(&active.work.flux)
            .any(|(a, b)| { a.iter().zip(b).any(|(a, b)| a != b) }),
        "reservoir gradients must affect interior modeled transport"
    );
}
