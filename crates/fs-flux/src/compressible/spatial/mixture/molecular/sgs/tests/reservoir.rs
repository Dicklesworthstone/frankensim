//! Complete-bath modeled transfer retains original species, calorics and work.

use super::*;
use crate::compressible::spatial::mixture::ReservoirFaceState;
use crate::compressible::spatial::mixture::molecular::{
    MolecularEuler3d, MolecularStepControl, StabilizedDissipationControl,
};
use std::cell::Cell;
use std::rc::Rc;

const OPEN_AXIS: AxisBoundaries = AxisBoundaries::Ends {
    lower: BoundaryCondition::Transmissive,
    upper: BoundaryCondition::Transmissive,
};
const OPEN: Boundaries3d = Boundaries3d {
    x: OPEN_AXIS,
    y: OPEN_AXIS,
    z: OPEN_AXIS,
};

fn bind(
    flow: &mut MixtureEuler3d,
    mut sample: impl FnMut(BoundarySide, [f64; 3]) -> ([f64; 3], [f64; 3], f64),
) {
    let bytes = flow
        .required_reservoir_profile_workspace_bytes(&BoundarySide::ALL)
        .unwrap();
    flow.prepare_reservoir_profiles(
        &BoundarySide::ALL,
        bytes,
        |side, _, point, densities| {
            let (partial, velocity, temperature) = sample(side, point);
            densities.copy_from_slice(&partial);
            Ok(ReservoirFaceState {
                velocity,
                pressure_pa: temperature * (partial[0] + 0.5 * partial[1] + 2.0 * partial[2]),
            })
        },
        || false,
    )
    .unwrap();
}

fn bath_field(formation: f64, shape: [usize; 3]) -> MixtureEuler3d {
    let grid = CartesianGrid::new(shape[0], shape[1], shape[2], 0.2, 0.3, 0.4).unwrap();
    let mut flow = field(material(formation), grid, OPEN, |_| {
        ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
    });
    bind(&mut flow, |_, _| {
        ([0.24, 0.96, 0.0], [0.02, 0.2, -0.05], 1.8)
    });
    flow
}

fn bath_owner(
    formation: f64,
    cfg: Option<MixtureSubgridConfig>,
    exchange: bool,
) -> MolecularEuler3d<ConstantTransport> {
    let flow = bath_field(formation, [2; 3]);
    let molecular = operator(&flow, cfg);
    let molecular = if exchange {
        molecular.with_local_sgs_reservoir_exchange().unwrap()
    } else {
        molecular
    };
    let bytes = MolecularEuler3d::required_workspace_bytes(&flow, &molecular).unwrap();
    MolecularEuler3d::new(flow, molecular, bytes, || false).unwrap()
}

#[test]
fn sgs_mixture_baths_match_six_full_tensors_midpoint_cp_and_prescribed_enthalpy_work() {
    let grid = CartesianGrid::new(3, 4, 5, 0.2, 0.3, 0.4).unwrap();
    for side in BoundarySide::ALL {
        let axis = side.axis();
        let [a, b] = tangents(axis);
        let profile = |x: [f64; 3]| {
            let density = 1.2 + 0.1 * x[a];
            let fraction = 0.25 + 0.03 * x[b];
            let mut velocity = [0.0; 3];
            velocity[axis] = 0.01 + 0.02 * x[a] * x[a] + 0.03 * x[b];
            velocity[a] = 0.04 * x[a] + 0.015 * x[b];
            velocity[b] = -0.05 * x[b] + 0.025 * x[a];
            (
                [density * fraction, density * (1.0 - fraction), 0.0],
                velocity,
                2.0 + 0.1 * x[b],
            )
        };
        let mut flow = field(material(0.0), grid, OPEN, |_| {
            ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
        });
        bind(&mut flow, |which, x| {
            if which == side {
                profile(x)
            } else {
                ([0.6, 0.4, 0.0], [0.0; 3], 1.0)
            }
        });
        for model in [
            SgsModel::Smagorinsky { cs: 0.17 },
            SgsModel::Wale { cw: 0.5 },
            SgsModel::Vreman { cv: 0.07 },
        ] {
            let mut cfg = config(model);
            cfg.transport.filter_width_factor = 1.3;
            let mut molecular_bath = operator(&flow, Some(cfg));
            let mut modeled_bath = operator(&flow, Some(cfg))
                .with_local_sgs_reservoir_exchange()
                .unwrap();
            assert!(!molecular_bath.local_sgs_reservoir_exchange_enabled());
            assert!(modeled_bath.local_sgs_reservoir_exchange_enabled());
            let base = faces(&mut molecular_bath, &flow);
            let active = faces(&mut modeled_bath, &flow);
            assert!(active.rate >= base.rate);
            assert_eq!(modeled_bath.work.gradient, molecular_bath.work.gradient);
            for face in 0..grid.face_count(side) {
                let x = grid.face_center(side, face).unwrap();
                let (partial, velocity, _) = profile(x);
                let bath = flow.reservoirs.primitive_face(side, face).unwrap();
                let mut c = [0; 3];
                c[axis] = if side.sign() < 0.0 {
                    0
                } else {
                    grid.shape[axis]
                };
                c[a] = face % grid.shape[a];
                c[b] = face / grid.shape[a];
                let index = grid.cell_face(axis, c);
                c[axis] = if side.sign() < 0.0 {
                    0
                } else {
                    grid.shape[axis] - 1
                };
                let cell = grid.index(c);
                let interior = modeled_bath.work.primitive[cell];
                let distance = 0.5 * grid.spacing[axis];
                // Analytic profile tensor: the uniform cell contributes zero
                // tangent derivative, the prescribed profile contributes half.
                let mut g = [[0.0; 3]; 3];
                for i in 0..3 {
                    g[i][axis] = side.sign() * velocity[i] / distance;
                }
                g[axis][a] = 0.02 * x[a];
                g[axis][b] = 0.015;
                g[a][a] = 0.02;
                g[a][b] = 0.0075;
                g[b][a] = 0.0125;
                g[b][b] = -0.025;
                let nu = match model {
                    SgsModel::Vreman { cv } => vreman_eddy_viscosity(
                        cv,
                        grid.spacing.map(|h| h * cfg.transport.filter_width_factor),
                        g,
                    ),
                    model => eddy_viscosity(
                        model,
                        grid.cell_volume().cbrt() * cfg.transport.filter_width_factor,
                        g,
                    ),
                }
                .unwrap();
                assert!(nu > 0.0);
                let mu = (0.5 * interior.density + 0.5 * bath.density) * nu;
                let bath_fraction = partial[0] / bath.density;
                let y = 0.5 * (0.6 / interior.density) + 0.5 * bath_fraction;
                let temperature = 0.5 * interior.temperature + 0.5 * bath.temperature;
                let cp =
                    y * (3.5 + 0.2 * temperature) + (1.0 - y) * 0.5 * (4.5 + 0.3 * temperature);
                let j = -mu / cfg.turbulent_schmidt
                    * side.sign()
                    * (bath_fraction - 0.6 / interior.density)
                    / distance;
                let heat = -mu * cp / cfg.transport.turbulent_prandtl
                    * side.sign()
                    * (bath.temperature - interior.temperature)
                    / distance;
                let enthalpy = |c: f64, slope: f64, inverse_mass: f64| {
                    inverse_mass
                        * ((c - 1.0) * (bath.temperature - 0.001)
                            + 0.5 * slope * (bath.temperature.powi(2) - 0.001_f64.powi(2))
                            + bath.temperature)
                };
                let mut energy = heat + j * (enthalpy(3.5, 0.2, 1.0) - enthalpy(4.5, 0.3, 0.5));
                let difference: [f64; 7] = std::array::from_fn(|k| {
                    modeled_bath.work.flux[axis][7 * index + k]
                        - molecular_bath.work.flux[axis][7 * index + k]
                });
                close(difference[0], j, 2e-13);
                close(difference[1], -j, 2e-13);
                assert_eq!(difference[2], 0.0);
                close(sum(difference[..3].iter().copied()), 0.0, 2e-14);
                let divergence = g[0][0] + g[1][1] + g[2][2];
                for i in 0..3 {
                    let stress = mu
                        * (g[i][axis] + g[axis][i]
                            - if i == axis {
                                (2.0 / 3.0) * divergence
                            } else {
                                0.0
                            });
                    close(difference[3 + i], -stress, 2e-13);
                    energy -= bath.velocity[i] * stress;
                }
                close(difference[6], energy, 8e-13);
            }
        }
    }
}

struct CountLaw {
    original: ConstantTransport,
    calls: Rc<Cell<usize>>,
}
impl TransportLaw for CountLaw {
    fn source(&self) -> &fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1 {
        self.original.source()
    }

    fn evaluate(
        &mut self,
        temperature: Temperature,
        pressure: fs_qty::Pressure,
        mole_fractions: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        self.calls.set(self.calls.get() + 1);
        self.original
            .evaluate(temperature, pressure, mole_fractions, binary, cancelled)
    }
}

#[test]
fn sgs_mixture_baths_reuse_each_original_transport_sample_once() {
    let flow = bath_field(0.0, [2; 3]);
    let cfg = config(SgsModel::Smagorinsky { cs: 0.2 });
    let mut results = Vec::new();
    for exchange in [false, true] {
        let calls = Rc::new(Cell::new(0));
        let source = SpatialTransport::new(
            flow.gas().clone(),
            CountLaw {
                original: original(&flow).law().clone(),
                calls: calls.clone(),
            },
            || false,
        )
        .unwrap();
        let bytes =
            MolecularTransport3d::required_workspace_bytes_with_sgs(&flow, &source).unwrap();
        let mut op = MolecularTransport3d::with_sgs(&flow, source, cfg, bytes, || false).unwrap();
        if exchange {
            op = op.with_local_sgs_reservoir_exchange().unwrap();
        }
        calls.set(0);
        faces(&mut op, &flow);
        assert_eq!(
            calls.get(),
            (0..3).map(|d| flow.grid.face_len(d)).sum::<usize>()
        );
        results.push(op.work.flux);
    }
    assert_ne!(results[0], results[1]);
}

#[test]
fn sgs_mixture_zero_bath_model_preserves_complete_native_steps_and_refuses_unowned_filters() {
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let mut molecular = bath_owner(0.0, None, false);
        let mut modeled = bath_owner(0.0, Some(config(model)), true);
        for _ in 0..2 {
            let expected = molecular
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            let actual = modeled
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            assert_eq!(actual, expected);
            assert_eq!(modeled.flow().states(), molecular.flow().states());
            for side in BoundarySide::ALL {
                assert_eq!(
                    modeled.boundary_impulses(side).collect::<Vec<_>>(),
                    molecular.boundary_impulses(side).collect::<Vec<_>>()
                );
            }
        }
    }
    let flow = field(
        material(0.0),
        CartesianGrid::new(5, 5, 5, 0.2, 0.3, 0.4).unwrap(),
        PERIODIC,
        |_| ([0.6, 0.4, 0.0], [0.0; 3], 1.0),
    );
    assert!(
        operator(&flow, None)
            .with_local_sgs_reservoir_exchange()
            .is_err()
    );
    let dynamic = MixtureDynamicSubgridConfig {
        transport: DynamicSubgridConfig {
            model: DynamicSgsConfig {
                averaging: DynamicAveraging::LocalBox,
                maximum_coefficient_squared: 1.0,
                maximum_eddy_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
            turbulent_prandtl: 0.7,
        },
        turbulent_schmidt: 0.8,
    };
    let source = original(&flow);
    let bytes =
        MolecularTransport3d::required_workspace_bytes_with_dynamic_sgs(&flow, &source).unwrap();
    let op =
        MolecularTransport3d::with_dynamic_sgs(&flow, source, dynamic, bytes, || false).unwrap();
    assert!(op.with_local_sgs_reservoir_exchange().is_err());
    let mut op = operator(&flow, None);
    op.sgs = Some(
        SubgridTransport::new_lagrangian(
            &flow,
            MixtureLagrangianSubgridConfig {
                transport: LagrangianSubgridConfig {
                    model: LagrangianSgsConfig {
                        memory_time: 0.01,
                        maximum_departure_cells: 1.0,
                        maximum_coefficient_squared: 1.0,
                        maximum_eddy_to_molecular: 1000.0,
                        usage: SgsUse::Exploration,
                    },
                    turbulent_prandtl: 0.7,
                },
                turbulent_schmidt: 0.8,
            },
        )
        .unwrap(),
    );
    assert!(op.with_local_sgs_reservoir_exchange().is_err());
}

#[test]
fn sgs_mixture_bath_species_bound_controls_ssp_and_rkl_and_refuses_atomically() {
    let mut cfg = config(SgsModel::Smagorinsky { cs: 0.25 });
    let flow = bath_field(0.0, [1; 3]);
    let mut ordinary = operator(&flow, Some(cfg))
        .with_local_sgs_reservoir_exchange()
        .unwrap();
    let ordinary_rate = faces(&mut ordinary, &flow).rate;
    cfg.turbulent_schmidt = 0.002;
    for stabilized in [false, true] {
        let mut flow = bath_field(0.0, [1; 3]);
        let source = original(&flow);
        let bytes = MolecularTransport3d::required_stabilized_workspace_bytes(&flow, &source)
            .unwrap()
            + 2 * flow.cell_count() * size_of::<f64>();
        let mut op = MolecularTransport3d::with_sgs(&flow, source, cfg, bytes, || false)
            .unwrap()
            .with_local_sgs_reservoir_exchange()
            .unwrap();
        assert!(faces(&mut op, &flow).rate > 5.0 * ordinary_rate);
        let report = if stabilized {
            op.advance_stabilized(
                &mut flow,
                Time::new(0.003),
                StabilizedDissipationControl::default(),
                || false,
            )
            .unwrap()
            .transport
        } else {
            let report = op
                .advance(
                    &mut flow,
                    Time::new(0.003),
                    DissipationControl::default(),
                    || false,
                )
                .unwrap();
            assert!(report.accepted_steps > 1);
            report
        };
        // RKL has its own stage-scaled parabolic admission, so check the
        // conservative inventories directly for either native integrator.
        for residual in report
            .conservation_residual
            .species_masses
            .into_iter()
            .chain(report.conservation_residual.momentum)
            .chain([report.conservation_residual.energy])
        {
            close(residual, 0.0, 3e-12);
        }
        close(report.absolute_energy_residual_j, 0.0, 3e-10);
        assert_eq!(flow.cell(0).unwrap().partial_densities[2], 0.0);
        let state = flow.states().to_vec();
        let loads = op.boundary_impulses.clone();
        let ViscosityClosure::Local { config, .. } = &mut op.sgs.as_mut().unwrap().closure else {
            unreachable!()
        };
        config.transport.maximum_eddy_to_molecular = 1e-12;
        assert!(
            op.advance(
                &mut flow,
                Time::new(0.001),
                DissipationControl::default(),
                || false,
            )
            .is_err()
        );
        assert_eq!(flow.states(), state);
        assert_eq!(op.boundary_impulses, loads);
        assert_eq!(op.last_step(), Some(&report));
    }
}

#[test]
fn sgs_mixture_baths_conserve_coupled_exchange_and_replay_after_late_cancellation() {
    let cfg = config(SgsModel::Smagorinsky { cs: 0.2 });
    let mut reference = bath_owner(0.0, Some(cfg), true);
    let mut shifted = bath_owner(1e12, Some(cfg), true);
    let mut unmodeled_baths = bath_owner(0.0, Some(cfg), false);
    let before = reference.flow().states().to_vec();
    let bytes = reference.workspace_bytes();
    let mut calls = 0;
    let report = reference
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    balance(&report.first);
    balance(&report.second);
    for residual in report
        .conservation_residual
        .species_masses
        .into_iter()
        .chain(report.conservation_residual.momentum)
        .chain([report.conservation_residual.energy])
    {
        close(residual, 0.0, 4e-12);
    }
    close(report.absolute_energy_residual_j, 0.0, 4e-10);
    shifted
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    unmodeled_baths
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(reference.flow().states(), shifted.flow().states());
    assert_ne!(reference.flow().states(), unmodeled_baths.flow().states());
    assert_ne!(reference.flow().states(), before);
    assert_eq!(reference.workspace_bytes(), bytes);
    let mut nonzero_species_exchange = false;
    for side in BoundarySide::ALL {
        for impulse in reference.molecular().boundary_impulses(side) {
            nonzero_species_exchange |= impulse.species_masses[0] != 0.0;
            assert_eq!(impulse.species_masses[2], 0.0);
            close(sum(impulse.species_masses[..3].iter().copied()), 0.0, 2e-14);
        }
    }
    assert!(nonzero_species_exchange);
    let mut refused = bath_owner(0.0, Some(cfg), true);
    let mut count = 0;
    assert!(matches!(
        refused.step(Time::new(0.001), MolecularStepControl::default(), || {
            count += 1;
            count == calls
        }),
        Err(DissipationError::Cancelled)
    ));
    assert_eq!(refused.flow().states(), before);
    assert_eq!(refused.flow().time().value(), 0.0);
    assert!(refused.last_step().is_none());
    assert!(refused.molecular().last_step().is_none());
    let replay = refused
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(replay, report);
    assert_eq!(refused.flow().states(), reference.flow().states());
    for side in BoundarySide::ALL {
        assert_eq!(
            refused.boundary_impulses(side).collect::<Vec<_>>(),
            reference.boundary_impulses(side).collect::<Vec<_>>()
        );
    }
}
