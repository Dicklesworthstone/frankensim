use super::*;
use crate::compressible::spatial::mixture::{BodyForce3d, ForceControl3d, ForcingWorkspace3d};
use fs_lbm::sgs3::dynamic::lagrangian::LagrangianDynamicSgs;

mod species;

type Coupled = MolecularEuler3d<ConstantTransport>;

fn memory_config() -> MixtureLagrangianSubgridConfig {
    MixtureLagrangianSubgridConfig {
        transport: LagrangianSubgridConfig {
            model: LagrangianSgsConfig {
                memory_time: 0.01,
                maximum_departure_cells: 1.0,
                maximum_coefficient_squared: 100.0,
                maximum_eddy_to_molecular: 1000.0,
                usage: SgsUse::Exploration,
            },
            turbulent_prandtl: 0.7,
        },
        turbulent_schmidt: 0.8,
    }
}

fn owner(formation: f64) -> (Coupled, ForcingWorkspace3d) {
    let flow = dynamic_field(formation);
    let molecular = operator(&flow, None);
    let work = flow.forcing_workspace(usize::MAX).unwrap();
    let bytes = Coupled::required_workspace_bytes(&flow, &molecular).unwrap()
        + Coupled::required_lagrangian_dynamic_sgs_workspace_bytes(&flow).unwrap()
        + work.workspace_bytes();
    let coupled =
        Coupled::with_lagrangian_dynamic_sgs(flow, molecular, memory_config(), bytes, || false)
            .unwrap();
    assert_eq!(coupled.workspace_bytes() + work.workspace_bytes(), bytes);
    (coupled, work)
}

fn physical_fields(flow: &MixtureEuler3d) -> (Vec<f64>, Vec<[f64; 3]>) {
    let mut density = Vec::with_capacity(flow.cell_count());
    let mut velocity = Vec::with_capacity(flow.cell_count());
    for i in 0..flow.cell_count() {
        let p = flow
            .cell(i)
            .unwrap()
            .primitive(flow.gas(), || false)
            .unwrap();
        density.push(p.density);
        velocity.push(p.velocity);
    }
    (density, velocity)
}

fn reference(flow: &MixtureEuler3d) -> LagrangianDynamicSgs {
    LagrangianDynamicSgs::new(flow.grid.shape, flow.grid.spacing, usize::MAX).unwrap()
}

fn conserved(report: &MolecularStepReport3d) {
    assert_eq!(report.elapsed, report.advection.elapsed);
    assert_eq!(report.first.elapsed, report.second.elapsed);
    assert_eq!(2.0 * report.first.elapsed.value(), report.elapsed.value());
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
    assert!(report.first.maximum_mass_frame_residual_kg_m2_s < 1e-12);
    assert!(report.second.maximum_mass_frame_residual_kg_m2_s < 1e-12);
}

#[test]
fn g0_lagrangian_mixture_faces_retain_mass_frame_stress_heat_and_enthalpy() {
    let flow = dynamic_field(0.0);
    let mut active = operator(&flow, None);
    let mut base = operator(&flow, None);
    active.sgs = Some(SubgridTransport::new_lagrangian(&flow, memory_config()).unwrap());
    let bytes = active.workspace_bytes();
    let mut history = active
        .sgs
        .as_mut()
        .unwrap()
        .take_lagrangian_filter()
        .unwrap();
    assert_eq!(
        active.workspace_bytes(),
        bytes,
        "detached history remains budgeted"
    );
    let width = history.filter_width();
    let trial = active
        .sgs
        .as_mut()
        .unwrap()
        .prepare_lagrangian(&mut history, &flow, &mut || false)
        .unwrap();
    let coefficients: Vec<_> = trial.coefficients_squared().collect();
    assert!(trial.report().maximum_coefficient_squared > 0.0);
    assert!(faces(&mut active, &flow).rate > faces(&mut base, &flow).rate);
    let mut nonzero = 0;
    for axis in 0..3 {
        for left in 0..flow.cell_count() {
            let mut c = flow.grid.coordinates(left);
            c[axis] = (c[axis] + 1) % flow.grid.shape[axis];
            let right = flow.grid.index(c);
            let face = flow.grid.cell_face(axis, c);
            let mu = dynamic_face_viscosity_reference(
                &active,
                [coefficients[left], coefficients[right]],
                width,
                left,
                right,
                axis,
            );
            if mu == 0.0 {
                continue;
            }
            nonzero += 1;
            let delta: [f64; 7] = std::array::from_fn(|k| {
                active.work.flux[axis][7 * face + k] - base.work.flux[axis][7 * face + k]
            });
            let pl = active.work.primitive[left];
            let pr = active.work.primitive[right];
            let h = flow.grid.spacing[axis];
            let yl = flow.states()[7 * left] / pl.density;
            let yr = flow.states()[7 * right] / pr.density;
            let y = 0.5 * yl + 0.5 * yr;
            let j = -mu / memory_config().turbulent_schmidt * (yr - yl) / h;
            close(delta[0], j, 3e-14);
            close(delta[1], -j, 3e-14);
            assert_eq!(delta[2], 0.0);
            let g = |i: usize, d: usize| {
                if d == axis {
                    (pr.velocity[i] - pl.velocity[i]) / h
                } else {
                    0.5 * active.work.gradient[left][d][i] + 0.5 * active.work.gradient[right][d][i]
                }
            };
            let divergence = g(0, 0) + g(1, 1) + g(2, 2);
            for component in 0..3 {
                let stress = mu
                    * (g(component, axis) + g(axis, component)
                        - if component == axis {
                            (2.0 / 3.0) * divergence
                        } else {
                            0.0
                        });
                close(delta[3 + component], -stress, 3e-14);
            }
            let t = 0.5 * pl.temperature + 0.5 * pr.temperature;
            let cp = y * (3.5 + 0.2 * t) + (1.0 - y) * 0.5 * (4.5 + 0.3 * t);
            let enthalpy = |c: f64, linear: f64, inverse_mass: f64| {
                inverse_mass
                    * ((c - 1.0) * (t - 0.001) + 0.5 * linear * (t * t - 0.001_f64.powi(2)) + t)
            };
            let heat = -mu * cp / memory_config().transport.turbulent_prandtl
                * (pr.temperature - pl.temperature)
                / h;
            let work: f64 = (0..3)
                .map(|a| (0.5 * pl.velocity[a] + 0.5 * pr.velocity[a]) * delta[3 + a])
                .sum();
            close(
                delta[6],
                heat + j * (enthalpy(3.5, 0.2, 1.0) - enthalpy(4.5, 0.3, 0.5)) + work,
                8e-14,
            );
        }
    }
    assert!(nonzero > 0);
    drop(trial);
    active
        .sgs
        .as_mut()
        .unwrap()
        .restore_lagrangian_filter(history, None);
    assert_eq!(active.workspace_bytes(), bytes);
    assert_eq!(active.sgs.as_ref().unwrap().lagrangian_sample_time(), None);
}

#[test]
fn g3_lagrangian_mixture_samples_actual_start_fields_and_preserves_source_gauge() {
    let (mut active, _) = owner(0.0);
    let (mut shifted, _) = owner(1e12);
    let mut history = reference(active.flow());
    assert_eq!(
        active.lagrangian_dynamic_subgrid_config(),
        Some(memory_config())
    );
    assert_eq!(active.lagrangian_subgrid_sample_time(), None);
    assert_eq!(active.lagrangian_subgrid_report(), None);
    let bytes = active.workspace_bytes();
    for requested in [0.001, 0.0005, 0.0015] {
        let start = active.flow().time();
        let (rho, velocity) = physical_fields(active.flow());
        let trial = history
            .prepare(
                memory_config().transport.model,
                start.value(),
                &rho,
                &velocity,
                || false,
            )
            .unwrap();
        let expected = trial.report();
        let report = active
            .step(
                Time::new(requested),
                MolecularStepControl::default(),
                || false,
            )
            .unwrap();
        shifted
            .step(
                Time::new(requested),
                MolecularStepControl::default(),
                || false,
            )
            .unwrap();
        conserved(&report);
        assert_eq!(active.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(active.lagrangian_subgrid_report(), Some(expected));
        assert_eq!(
            active.flow().time(),
            Time::new(start.value() + report.elapsed.value())
        );
        assert_eq!(active.flow().states(), shifted.flow().states());
        assert!(
            active
                .flow()
                .states()
                .chunks_exact(7)
                .all(|row| row[2] == 0.0)
        );
        assert_eq!(active.workspace_bytes(), bytes);
        trial.commit();
    }
}

#[test]
fn g4_lagrangian_mixture_shortened_forced_trials_and_cancellation_replay() {
    let force = BodyForce3d::UniformAcceleration([200.0, 0.0, 0.0]);
    let force_controls = ForceControl3d { maximum_mach: 10.0 };
    let make = || {
        let (mut flow, mut work) = owner(0.0);
        flow.step_with_forcing(
            Time::new(0.0001),
            MolecularStepControl::default(),
            force,
            force_controls,
            &mut work,
            || false,
        )
        .unwrap();
        (flow, work)
    };
    let (mut probe, mut work) = make();
    let start = probe.flow().time();
    let old_sample = probe.lagrangian_subgrid_sample_time();
    let old_fit = probe.lagrangian_subgrid_report();
    let old_states = probe.flow().states().to_vec();
    let old_step = probe.last_step().copied();
    let old_molecular = probe.molecular().last_step().copied();
    let old_advection = probe.flow().last_step().copied();
    let bytes = probe.workspace_bytes();
    let mut polls = 0;
    let accepted = probe
        .step_with_forcing(
            Time::new(1.0),
            MolecularStepControl::default(),
            force,
            force_controls,
            &mut work,
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(accepted.coupling_retries > 0);
    assert!(accepted.elapsed.value() < 1.0);
    conserved(&accepted);
    assert_eq!(probe.lagrangian_subgrid_sample_time(), Some(start));
    let fitted = probe.lagrangian_subgrid_report();
    let accepted_states = probe.flow().states().to_vec();
    let next = probe
        .step(Time::new(0.0005), MolecularStepControl::default(), || false)
        .unwrap();
    let next_fit = probe.lagrangian_subgrid_report();

    let (mut direct, mut direct_work) = make();
    let short = direct
        .step_with_forcing(
            accepted.elapsed,
            MolecularStepControl::default(),
            force,
            force_controls,
            &mut direct_work,
            || false,
        )
        .unwrap();
    assert_eq!(short.coupling_retries, 0);
    assert_eq!(short.first, accepted.first);
    assert_eq!(short.advection, accepted.advection);
    assert_eq!(short.second, accepted.second);
    assert_eq!(direct.flow().states(), accepted_states);
    assert_eq!(direct.lagrangian_subgrid_report(), fitted);
    assert_eq!(
        direct
            .step(Time::new(0.0005), MolecularStepControl::default(), || false)
            .unwrap(),
        next
    );
    assert_eq!(direct.lagrangian_subgrid_report(), next_fit);
    assert_eq!(direct.flow().states(), probe.flow().states());

    for cutoff in [polls / 2, polls] {
        let (mut flow, mut work) = make();
        let mut calls = 0;
        assert_eq!(
            flow.step_with_forcing(
                Time::new(1.0),
                MolecularStepControl::default(),
                force,
                force_controls,
                &mut work,
                || {
                    calls += 1;
                    calls == cutoff
                },
            ),
            Err(DissipationError::Cancelled)
        );
        assert_eq!(flow.flow().states(), old_states);
        assert_eq!(flow.flow().time(), start);
        assert_eq!(flow.last_step().copied(), old_step);
        assert_eq!(flow.molecular().last_step().copied(), old_molecular);
        assert_eq!(flow.flow().last_step().copied(), old_advection);
        assert_eq!(flow.lagrangian_subgrid_sample_time(), old_sample);
        assert_eq!(flow.lagrangian_subgrid_report(), old_fit);
        assert_eq!(flow.workspace_bytes(), bytes);
        assert_eq!(
            flow.step_with_forcing(
                Time::new(1.0),
                MolecularStepControl::default(),
                force,
                force_controls,
                &mut work,
                || false,
            )
            .unwrap(),
            accepted
        );
        assert_eq!(flow.flow().states(), accepted_states);
        assert_eq!(flow.lagrangian_subgrid_report(), fitted);
        assert_eq!(
            flow.step(Time::new(0.0005), MolecularStepControl::default(), || false)
                .unwrap(),
            next
        );
        assert_eq!(flow.lagrangian_subgrid_report(), next_fit);
        assert_eq!(flow.flow().states(), probe.flow().states());
    }
}

fn stabilized_owner() -> Coupled {
    let flow = dynamic_field(0.0);
    let law = ConstantTransport::new(
        flow.gas().source().clone(),
        "manufactured-lagrangian-stabilized-mixture",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 0.001,
            maximum_pressure_pa: 100.0,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(2.0),
            bulk_viscosity: DynViscosity::new(0.3),
            thermal_conductivity: ThermalConductivityQuantityV1::new(3.0),
        },
        vec![0.0, 1.0, 1.2, 1.0, 0.0, 1.3, 1.2, 1.3, 0.0],
    )
    .unwrap();
    let source = SpatialTransport::new(flow.gas().clone(), law, || false).unwrap();
    let bytes = MolecularTransport3d::required_workspace_bytes(&flow, &source).unwrap();
    let molecular = MolecularTransport3d::new(&flow, source, bytes, || false).unwrap();
    let total = Coupled::required_stabilized_workspace_bytes(&flow, &molecular).unwrap()
        + Coupled::required_lagrangian_dynamic_sgs_workspace_bytes(&flow).unwrap();
    Coupled::with_lagrangian_dynamic_sgs(flow, molecular, memory_config(), total, || false).unwrap()
}

#[test]
fn g4_lagrangian_mixture_rkl_stages_hold_one_fit_and_publish_atomically() {
    let mut active = stabilized_owner();
    let mut replay = stabilized_owner();
    let mut history = reference(active.flow());
    for _ in 0..2 {
        let start = active.flow().time();
        let (rho, velocity) = physical_fields(active.flow());
        let trial = history
            .prepare(
                memory_config().transport.model,
                start.value(),
                &rho,
                &velocity,
                || false,
            )
            .unwrap();
        let expected = trial.report();
        let mut polls = 0;
        let report = active
            .step_stabilized(
                Time::new(0.002),
                StabilizedMolecularStepControl::default(),
                || {
                    polls += 1;
                    false
                },
            )
            .unwrap();
        assert!(report.maximum_stages_used > 2);
        assert_eq!(active.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(active.lagrangian_subgrid_report(), Some(expected));
        conserved(&report.transport);
        trial.commit();
        let before = replay.flow().states().to_vec();
        let prior = replay.lagrangian_subgrid_report();
        let sample = replay.lagrangian_subgrid_sample_time();
        let mut calls = 0;
        assert_eq!(
            replay.step_stabilized(
                Time::new(0.002),
                StabilizedMolecularStepControl::default(),
                || {
                    calls += 1;
                    calls == polls
                },
            ),
            Err(DissipationError::Cancelled)
        );
        assert_eq!(replay.flow().states(), before);
        assert_eq!(replay.flow().time(), start);
        assert_eq!(replay.lagrangian_subgrid_report(), prior);
        assert_eq!(replay.lagrangian_subgrid_sample_time(), sample);
        assert_eq!(
            replay
                .step_stabilized(
                    Time::new(0.002),
                    StabilizedMolecularStepControl::default(),
                    || false,
                )
                .unwrap(),
            report
        );
        assert_eq!(replay.lagrangian_subgrid_report(), Some(expected));
        assert_eq!(replay.flow().states(), active.flow().states());
        assert_eq!(replay.workspace_bytes(), active.workspace_bytes());
    }
}

#[test]
fn g0_lagrangian_mixture_budget_refusals_and_exact_molecular_limit() {
    let flow = dynamic_field(0.0);
    let molecular = operator(&flow, None);
    let additional = Coupled::required_lagrangian_dynamic_sgs_workspace_bytes(&flow).unwrap();
    assert_eq!(additional, 216 * flow.cell_count());
    let total = Coupled::required_workspace_bytes(&flow, &molecular).unwrap() + additional;
    assert!(Coupled::with_lagrangian_dynamic_sgs(
        flow, molecular, memory_config(), total - 1, || false,
    ).is_err());
    for schmidt in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            MixtureLagrangianSubgridConfig {
                turbulent_schmidt: schmidt,
                ..memory_config()
            }
            .validate()
            .is_err()
        );
    }
    let mut forbidden = memory_config();
    forbidden.transport.model.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    let flow = dynamic_field(0.0);
    let molecular = dynamic_operator(&flow, dynamic_config());
    assert!(
        Coupled::with_lagrangian_dynamic_sgs(flow, molecular, memory_config(), usize::MAX, || {
            false
        },)
        .is_err()
    );
    for (shape, boundaries) in [([2, 4, 4], PERIODIC), ([4; 3], WALLS)] {
        let flow = field(
            material(0.0),
            CartesianGrid::new(shape[0], shape[1], shape[2], 0.25, 0.25, 0.25).unwrap(),
            boundaries,
            |_| ([0.6, 0.4, 0.0], [0.0; 3], 1.0),
        );
        let molecular = operator(&flow, None);
        assert!(
            Coupled::with_lagrangian_dynamic_sgs(
                flow,
                molecular,
                memory_config(),
                usize::MAX,
                || false,
            )
            .is_err()
        );
    }
    let flow = dynamic_field(0.0);
    let molecular = operator(&flow, None);
    let before = flow.states().to_vec();
    forbidden = memory_config();
    forbidden.transport.model.maximum_coefficient_squared = 0.0;
    let mut refused =
        Coupled::with_lagrangian_dynamic_sgs(flow, molecular, forbidden, usize::MAX, || false)
            .unwrap();
    for _ in 0..2 {
        assert!(
            refused
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .is_err()
        );
        assert_eq!(refused.flow().states(), before);
        assert_eq!(refused.flow().time(), Time::new(0.0));
        assert_eq!(refused.lagrangian_subgrid_sample_time(), None);
        assert_eq!(refused.lagrangian_subgrid_report(), None);
    }
    let uniform = || {
        field(
            material(0.0),
            CartesianGrid::new(4, 4, 4, 0.25, 0.25, 0.25).unwrap(),
            PERIODIC,
            |_| ([0.6, 0.4, 0.0], [0.1, -0.05, 0.025], 1.0),
        )
    };
    let flow = uniform();
    let molecular = operator(&flow, None);
    let mut base = Coupled::new(flow, molecular, usize::MAX, || false).unwrap();
    let flow = uniform();
    let molecular = operator(&flow, None);
    let mut active =
        Coupled::with_lagrangian_dynamic_sgs(flow, molecular, memory_config(), usize::MAX, || {
            false
        })
        .unwrap();
    for _ in 0..3 {
        let expected = base
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let actual = active
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(actual, expected);
        assert_eq!(
            active
                .lagrangian_subgrid_report()
                .unwrap()
                .maximum_coefficient_squared,
            0.0
        );
        for (actual, expected) in active.flow().states().iter().zip(base.flow().states()) {
            assert_eq!(actual.to_bits(), expected.to_bits());
        }
    }
    // A detached composition cannot accidentally age history in the molecular
    // pseudo-clock, including when handed to another split/chemistry composer.
    let (mut flow, mut molecular) = active.into_parts();
    let states = flow.states().to_vec();
    assert!(
        molecular
            .advance(
                &mut flow,
                Time::new(0.001),
                DissipationControl::default(),
                || false,
            )
            .is_err()
    );
    assert_eq!(flow.states(), states);
}
