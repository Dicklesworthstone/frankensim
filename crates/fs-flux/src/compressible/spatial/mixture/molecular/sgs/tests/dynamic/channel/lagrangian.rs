use super::*;
use crate::compressible::spatial::molecular::sgs::channel_wall_values;
use fs_lbm::sgs3::dynamic::lagrangian::scalar::components::LagrangianDynamicComponents;

type Coupled = MolecularEuler3d<ConstantTransport>;

fn memory() -> MixtureLagrangianSpeciesSubgridConfig {
    MixtureLagrangianSpeciesSubgridConfig {
        momentum: LagrangianSgsConfig {
            memory_time: 0.01,
            maximum_departure_cells: 1.0,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        heat: LagrangianScalarConfig {
            memory_time: 0.017,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        species: LagrangianComponentsConfig {
            memory_time: 0.023,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            usage: SgsUse::Exploration,
        },
        maximum_species_diffusivity_m2_s: 10.0,
    }
}

fn owner(normal: usize) -> Coupled {
    let flow = channel(normal, [0.09, 0.14, 0.20], true, 0.0);
    let op = operator(&flow, None);
    let extra = Coupled::required_lagrangian_dynamic_species_sgs_workspace_bytes(&flow).unwrap();
    assert_eq!(extra, 448 * flow.cell_count());
    let cap = Coupled::required_workspace_bytes(&flow, &op).unwrap() + extra;
    let result = Coupled::with_lagrangian_dynamic_channel_species_sgs(
        flow,
        op,
        memory(),
        walls(),
        cap,
        || false,
    )
    .unwrap();
    assert_eq!(result.workspace_bytes(), cap);
    assert_eq!(result.lagrangian_channel_normal_axis(), Some(normal));
    result
}

fn physical(flow: &MixtureEuler3d) -> (Vec<f64>, Vec<[f64; 3]>, Vec<f64>) {
    let p: Vec<_> = (0..flow.cell_count())
        .map(|i| {
            flow.cell(i)
                .unwrap()
                .primitive(flow.gas(), || false)
                .unwrap()
        })
        .collect();
    (
        p.iter().map(|p| p.density).collect(),
        p.iter().map(|p| p.velocity).collect(),
        p.iter().map(|p| p.temperature).collect(),
    )
}

fn reports(
    owner: &Coupled,
) -> (
    Option<DynamicSgsReport>,
    Option<DynamicScalarReport>,
    Option<DynamicScalarReport>,
) {
    (
        owner.lagrangian_subgrid_report(),
        owner.lagrangian_subgrid_heat_report(),
        owner.lagrangian_subgrid_species_report(),
    )
}

#[test]
fn g1_channel_history_uses_current_walls_all_axes_and_original_wall_flux() {
    for normal in 0..3 {
        let (flow, mut active) = owner(normal).into_parts();
        let mut wall_velocity = [0.0; 3];
        wall_velocity[(normal + 1) % 3] = 0.015;
        let replacement = MolecularWall3d::NoSlip {
            velocity: wall_velocity,
            thermal: WallThermalCondition::Adiabatic,
        };
        active = active
            .with_wall(BoundarySide::ALL[2 * normal], replacement)
            .unwrap();
        let mut base = baseline(&flow, original(&flow), normal, [replacement, walls()[1]]);
        let bytes = active.workspace_bytes();
        let (rho, velocity, temperature) = physical(&flow);
        let (wall_velocity, scalar_walls) = channel_wall_values(normal, &active.walls).unwrap();
        let mut reference = LagrangianDynamicComponents::new_channel(
            flow.grid.shape,
            flow.grid.spacing,
            normal,
            3,
            usize::MAX,
        )
        .unwrap();
        let cfg = memory();
        let expected = reference
            .prepare_with_walls(
                cfg.momentum,
                cfg.heat,
                cfg.species,
                flow.time().value(),
                &rho,
                &velocity,
                wall_velocity,
                &temperature,
                scalar_walls,
                |component, cell| flow.states()[7 * cell + component] / rho[cell],
                || false,
            )
            .unwrap();
        let mut history = active
            .sgs
            .as_mut()
            .unwrap()
            .take_lagrangian_filter()
            .unwrap();
        assert_eq!(active.workspace_bytes(), bytes);
        let trial = active
            .sgs
            .as_mut()
            .unwrap()
            .prepare_lagrangian_with_walls(&mut history, &flow, &active.walls, &mut || false)
            .unwrap();
        assert_eq!(
            trial.reports(),
            (
                expected.momentum_report(),
                Some(expected.scalar_report()),
                Some(expected.components_report())
            )
        );
        stage(&mut active, &flow);
        stage(&mut base, &flow);
        check_wall_flux(&active, &base, &flow, normal);
        let model = active.sgs.as_ref().unwrap();
        let ViscosityClosure::Lagrangian(fields) = &model.closure else {
            unreachable!()
        };
        for cell in 0..flow.cell_count() {
            assert_eq!(
                model.dynamic_coefficient_squared[cell],
                expected.coefficient_squared(cell).unwrap()
            );
            assert_eq!(
                fields.face_heat_coefficient(cell, cell).unwrap().0,
                expected.scalar_coefficient(cell).unwrap()
            );
            assert_eq!(
                fields.face_species_coefficient(cell, cell).unwrap().0,
                expected.component_coefficient(cell).unwrap()
            );
        }
        assert!(
            active.work.flux != base.work.flux,
            "modeled interior transport must be active"
        );
        drop(trial);
        active
            .sgs
            .as_mut()
            .unwrap()
            .restore_lagrangian_filter(history, None);
        assert_eq!(active.workspace_bytes(), bytes);
    }
}

#[test]
fn g3_channel_histories_follow_accepted_start_time_and_conserve_with_physical_walls() {
    let mut active = owner(2);
    let mut reference = LagrangianDynamicComponents::new_channel(
        active.flow().grid.shape,
        active.flow().grid.spacing,
        2,
        3,
        usize::MAX,
    )
    .unwrap();
    let bytes = active.workspace_bytes();
    for duration in [0.0005, 0.0002, 0.0007] {
        let start = active.flow().time();
        let (rho, velocity, temperature) = physical(active.flow());
        let (wall_velocity, scalar_walls) =
            channel_wall_values(2, &active.molecular().walls).unwrap();
        let cfg = memory();
        let trial = reference
            .prepare_with_walls(
                cfg.momentum,
                cfg.heat,
                cfg.species,
                start.value(),
                &rho,
                &velocity,
                wall_velocity,
                &temperature,
                scalar_walls,
                |component, cell| active.flow().states()[7 * cell + component] / rho[cell],
                || false,
            )
            .unwrap();
        let expected = (
            Some(trial.momentum_report()),
            Some(trial.scalar_report()),
            Some(trial.components_report()),
        );
        let report = active
            .step(Time::new(duration), MolecularStepControl::default(), || {
                false
            })
            .unwrap();
        assert_eq!(reports(&active), expected);
        assert_eq!(active.lagrangian_subgrid_sample_time(), Some(start));
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

    // Cancelling immediately before complete publication restores the accepted
    // flow, loads, report and every history. A retry matches a clean twin.
    let mut twin = owner(2);
    for duration in [0.0005, 0.0002, 0.0007] {
        twin.step(Time::new(duration), MolecularStepControl::default(), || {
            false
        })
        .unwrap();
    }
    let mut polls = 0usize;
    let next = twin
        .step(Time::new(0.0003), MolecularStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    let before = active.flow().states().to_vec();
    let before_reports = reports(&active);
    let before_time = active.lagrangian_subgrid_sample_time();
    let before_step = active.last_step().copied();
    let before_loads = active.molecular().boundary_impulses.clone();
    let mut count = 0usize;
    assert!(
        active
            .step(Time::new(0.0003), MolecularStepControl::default(), || {
                count += 1;
                count == polls
            })
            .is_err()
    );
    assert_eq!(active.flow().states(), before);
    assert_eq!(reports(&active), before_reports);
    assert_eq!(active.lagrangian_subgrid_sample_time(), before_time);
    assert_eq!(active.last_step().copied(), before_step);
    assert_eq!(active.molecular().boundary_impulses, before_loads);
    assert_eq!(
        active
            .step(Time::new(0.0003), MolecularStepControl::default(), || false)
            .unwrap(),
        next
    );
    assert_eq!(active.flow().states(), twin.flow().states());
    assert_eq!(reports(&active), reports(&twin));
}

#[test]
fn g0_channel_history_admission_and_zero_limit_keep_original_molecular_bits() {
    let make_flow = || {
        let grid = CartesianGrid::new(3, 3, 3, 0.1, 0.2, 0.3).unwrap();
        field(material(0.0), grid, boundaries(1), |_| {
            ([0.5, 0.5, 0.0], [0.0; 3], 1.0)
        })
    };
    let wall = MolecularWall3d::NoSlip {
        velocity: [0.0; 3],
        thermal: WallThermalCondition::Adiabatic,
    };
    let flow = make_flow();
    let op = operator(&flow, None);
    let cap = Coupled::required_workspace_bytes(&flow, &op).unwrap()
        + Coupled::required_lagrangian_dynamic_species_sgs_workspace_bytes(&flow).unwrap();
    assert!(
        Coupled::with_lagrangian_dynamic_channel_species_sgs(
            flow,
            op,
            memory(),
            [wall; 2],
            cap - 1,
            || false
        )
        .is_err()
    );
    let flow = make_flow();
    let op = operator(&flow, None);
    assert!(
        Coupled::with_lagrangian_dynamic_channel_species_sgs(
            flow,
            op,
            memory(),
            [MolecularWall3d::SlipAdiabatic; 2],
            cap,
            || false
        )
        .is_err()
    );
    let flow = make_flow();
    let op = operator(&flow, None);
    let mut active = Coupled::with_lagrangian_dynamic_channel_species_sgs(
        flow,
        op,
        memory(),
        [wall; 2],
        cap,
        || false,
    )
    .unwrap();
    let flow = make_flow();
    let op = baseline(&flow, original(&flow), 1, [wall; 2]);
    let mut molecular = Coupled::new(flow, op, cap, || false).unwrap();
    let a = active
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    let b = molecular
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(a, b);
    assert_eq!(active.flow().states(), molecular.flow().states());
    assert_eq!(
        active
            .lagrangian_subgrid_report()
            .unwrap()
            .maximum_coefficient_squared,
        0.0
    );
    assert_eq!(
        active
            .lagrangian_subgrid_heat_report()
            .unwrap()
            .maximum_coefficient,
        0.0
    );
    assert_eq!(
        active
            .lagrangian_subgrid_species_report()
            .unwrap()
            .maximum_coefficient,
        0.0
    );
    let (flow, op) = active.into_parts();
    assert!(
        op.with_wall(BoundarySide::YMin, MolecularWall3d::SlipAdiabatic)
            .is_err()
    );
    assert_eq!(flow.time(), Time::new(0.001));
}
