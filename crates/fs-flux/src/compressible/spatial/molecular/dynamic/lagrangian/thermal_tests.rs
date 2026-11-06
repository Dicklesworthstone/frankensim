use super::*;
use crate::compressible::spatial::molecular::sgs::{
    LagrangianScalarConfig, LagrangianThermalSubgridConfig,
};
use fs_lbm::sgs3::dynamic::lagrangian::scalar::LagrangianDynamicScalar;

fn settings() -> LagrangianThermalSubgridConfig {
    LagrangianThermalSubgridConfig {
        momentum: memory_config().model,
        heat: LagrangianScalarConfig {
            memory_time: 0.017,
            maximum_departure_cells: 1.0,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

fn thermal_owner(mu: f64, k: f64) -> MolecularEuler3d<ConstantTransport> {
    let flow = cube(6, mixed);
    let source = law(flow.gas(), mu, 0.0, k);
    MolecularEuler3d::with_lagrangian_dynamic_thermal_sgs(
        flow,
        source,
        settings(),
        320 * 6usize.pow(3),
        || false,
    )
    .unwrap()
}

fn temperature(flow: &Euler3d) -> Vec<f64> {
    flow.states()
        .iter()
        .map(|state| state.primitive(flow.gas(), || false).unwrap().temperature)
        .collect()
}

#[test]
fn g0_lagrangian_heat_matches_both_shared_histories_at_physical_start_times() {
    let mut flow = thermal_owner(0.002, 0.004);
    let mut history =
        LagrangianDynamicScalar::new(flow.flow().grid.shape, flow.flow().grid.spacing, usize::MAX)
            .unwrap();
    assert_eq!(flow.lagrangian_thermal_subgrid_config(), Some(settings()));
    assert_eq!(flow.lagrangian_dynamic_subgrid_config(), None);
    assert_eq!(flow.dynamic_thermal_subgrid_config(), None);
    assert_eq!(flow.lagrangian_subgrid_heat_report(), None);
    for requested in [0.002, 0.001, 0.003] {
        let start = flow.flow().time();
        let (rho, velocity) = physical_fields(flow.flow());
        let scalar = temperature(flow.flow());
        let trial = history
            .prepare(
                settings().momentum,
                settings().heat,
                start.value(),
                &rho,
                &velocity,
                &scalar,
                || false,
            )
            .unwrap();
        let expected_momentum = trial.momentum_report();
        let expected_heat = trial.scalar_report();
        let report = flow
            .step(
                Time::new(requested),
                MolecularStepControl::default(),
                || false,
            )
            .unwrap();
        assert_eq!(flow.lagrangian_subgrid_sample_time(), Some(start));
        assert_eq!(flow.lagrangian_subgrid_report(), Some(expected_momentum));
        assert_eq!(flow.lagrangian_subgrid_heat_report(), Some(expected_heat));
        assert_eq!(
            flow.flow().time(),
            Time::new(start.value() + report.elapsed.value())
        );
        trial.commit();
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 4e-13);
        }
    }
}

#[test]
fn g0_lagrangian_heat_without_modeled_stress_uses_actual_conservative_faces() {
    let n: usize = 9;
    let flow = cube(n, shear_temperature);
    let mut molecular =
        MolecularTransport3d::new(&flow, law(flow.gas(), 0.0, 0.0, 0.001), || false).unwrap();
    let mut active =
        MolecularTransport3d::new(&flow, law(flow.gas(), 0.0, 0.0, 0.001), || false).unwrap();
    active.sgs =
        Some(SubgridTransport::new_lagrangian_thermal(&flow, settings(), 320 * n.pow(3)).unwrap());
    let mut filter = active
        .sgs
        .as_mut()
        .unwrap()
        .take_lagrangian_filter()
        .unwrap();
    let width = filter.filter_width();
    let trial = active
        .sgs
        .as_mut()
        .unwrap()
        .prepare_lagrangian(&mut filter, &flow, &mut || false)
        .unwrap();
    let (momentum, heat) = trial.reports();
    assert_eq!(momentum.maximum_coefficient_squared, 0.0);
    assert!(heat.unwrap().maximum_coefficient > 0.0);
    let base_rate = prepare(&mut molecular, &flow);
    let heat_rate = prepare(&mut active, &flow);
    assert!(heat_rate > base_rate);

    let (rho, velocity) = physical_fields(&flow);
    let scalar = temperature(&flow);
    let mut oracle =
        LagrangianDynamicScalar::new(flow.grid.shape, flow.grid.spacing, usize::MAX).unwrap();
    let reference = oracle
        .prepare(
            settings().momentum,
            settings().heat,
            flow.time().value(),
            &rho,
            &velocity,
            &scalar,
            || false,
        )
        .unwrap();
    let coefficients: Vec<_> = reference.scalar_coefficients().collect();
    let mut faces = 0;
    for right in 0..flow.states().len() {
        let coordinate = flow.grid.coordinates(right);
        let mut left_coordinate = coordinate;
        left_coordinate[0] = (coordinate[0] + n - 1) % n;
        let left = flow.grid.index(left_coordinate);
        // The existing independent full face-gradient oracle applies equally
        // to C_phi Delta^2 |S|, without a turbulent-Prandtl substitution.
        let rho_alpha = face_viscosity_reference(
            &active,
            [coefficients[left], coefficients[right]],
            width,
            left,
            right,
            0,
        );
        let face = flow.grid.cell_face(0, coordinate);
        let a = active.work.flux[0][face];
        let b = molecular.work.flux[0][face];
        assert_eq!(a.normal_momentum, b.normal_momentum);
        assert_eq!(a.transverse_momentum, b.transverse_momentum);
        close(
            a.energy - b.energy,
            -rho_alpha * 3.5 * (scalar[right] - scalar[left]) / flow.grid.spacing[0],
            2e-14,
        );
        faces += usize::from(rho_alpha > 0.0);
    }
    assert!(faces > 0);
    drop(trial);
    active
        .sgs
        .as_mut()
        .unwrap()
        .restore_lagrangian_filter(filter, None);
    assert_eq!(active.sgs.as_ref().unwrap().lagrangian_heat_report(), None);

    let initial = cube(n, shear_temperature);
    let source = law(initial.gas(), 0.002, 0.0, 0.001);
    let mut baseline = MolecularEuler3d::new(initial, source, || false).unwrap();
    let initial = cube(n, shear_temperature);
    let source = law(initial.gas(), 0.002, 0.0, 0.001);
    let mut coupled = MolecularEuler3d::with_lagrangian_dynamic_thermal_sgs(
        initial,
        source,
        settings(),
        320 * n.pow(3),
        || false,
    )
    .unwrap();
    baseline
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    let report = coupled
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_ne!(coupled.flow().states(), baseline.flow().states());
    assert_eq!(
        coupled
            .lagrangian_subgrid_report()
            .unwrap()
            .maximum_coefficient_squared,
        0.0
    );
    assert!(
        coupled
            .lagrangian_subgrid_heat_report()
            .unwrap()
            .maximum_coefficient
            > 0.0
    );
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 4e-13);
    }
}

#[test]
fn g4_lagrangian_heat_observations_and_cancelled_rkl_steps_preserve_both_memories() {
    let make = || {
        let mut flow = thermal_owner(2.0, 4.0);
        flow.step(Time::new(0.0001), MolecularStepControl::default(), || false)
            .unwrap();
        flow
    };
    let mut control = make();
    let mut observed = make();
    let mut work = observed.vorticity_workspace(usize::MAX).unwrap();
    let start = observed.flow().time();
    let old_fields = observed.flow().states().to_vec();
    let old_sample = observed.lagrangian_subgrid_sample_time();
    let old_momentum = observed.lagrangian_subgrid_report();
    let old_heat = observed.lagrangian_subgrid_heat_report();
    let old_step = observed.last_step().copied();
    let first = observed.vorticity(&mut work, || false).unwrap().to_vec();
    assert_eq!(observed.vorticity(&mut work, || false).unwrap(), first);
    assert_eq!(observed.lagrangian_subgrid_sample_time(), old_sample);
    assert_eq!(observed.lagrangian_subgrid_report(), old_momentum);
    assert_eq!(observed.lagrangian_subgrid_heat_report(), old_heat);

    let mut polls = 0;
    let expected = control
        .step_stabilized(
            Time::new(0.004),
            MolecularStepControl::default(),
            MolecularStabilization::default(),
            || {
                polls += 1;
                false
            },
        )
        .unwrap();
    assert!(expected.maximum_stages_used > 2);
    let mut calls = 0;
    assert_eq!(
        observed.step_stabilized(
            Time::new(0.004),
            MolecularStepControl::default(),
            MolecularStabilization::default(),
            || {
                calls += 1;
                calls == polls
            },
        ),
        Err(MolecularError::Cancelled),
    );
    assert_eq!(observed.flow().states(), old_fields);
    assert_eq!(observed.flow().time(), start);
    assert_eq!(observed.last_step().copied(), old_step);
    assert_eq!(observed.lagrangian_subgrid_sample_time(), old_sample);
    assert_eq!(observed.lagrangian_subgrid_report(), old_momentum);
    assert_eq!(observed.lagrangian_subgrid_heat_report(), old_heat);
    assert_eq!(
        observed
            .step_stabilized(
                Time::new(0.004),
                MolecularStepControl::default(),
                MolecularStabilization::default(),
                || false,
            )
            .unwrap(),
        expected
    );
    assert_eq!(observed.lagrangian_subgrid_sample_time(), Some(start));
    assert_eq!(
        observed.lagrangian_subgrid_report(),
        control.lagrangian_subgrid_report()
    );
    assert_eq!(
        observed.lagrangian_subgrid_heat_report(),
        control.lagrangian_subgrid_heat_report()
    );
    assert_eq!(observed.flow().states(), control.flow().states());

    // Accepted histories affect the next fit; comparing one later complete
    // step detects an erroneously committed diagnostic or rejected trial.
    let expected = control
        .step(Time::new(0.0001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(
        observed
            .step(Time::new(0.0001), MolecularStepControl::default(), || false,)
            .unwrap(),
        expected
    );
    assert_eq!(
        observed.lagrangian_subgrid_heat_report(),
        control.lagrangian_subgrid_heat_report()
    );
    for (a, b) in observed.flow().states().iter().zip(control.flow().states()) {
        assert_eq!(a.array().map(f64::to_bits), b.array().map(f64::to_bits));
    }
}

#[test]
fn g0_lagrangian_heat_budget_refusal_and_exact_uniform_molecular_limit() {
    let flow = cube(9, shear_temperature);
    let bytes = MolecularEuler3d::<ConstantTransport>::
        required_lagrangian_dynamic_thermal_sgs_workspace_bytes(&flow).unwrap();
    assert_eq!(bytes, 320 * 9usize.pow(3));
    let source = law(flow.gas(), 0.002, 0.0, 0.001);
    assert!(matches!(
        MolecularEuler3d::with_lagrangian_dynamic_thermal_sgs(
            flow,
            source,
            settings(),
            bytes - 1,
            || false,
        ),
        Err(MolecularError::DynamicSubgrid(
            DynamicSgsError::WorkspaceBytes { .. }
        )),
    ));
    let mut forbidden = settings();
    forbidden.heat.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    forbidden = settings();
    forbidden.heat.maximum_diffusivity_to_molecular = 0.0;
    let flow = cube(9, shear_temperature);
    let source = law(flow.gas(), 0.002, 0.0, 0.001);
    let old = flow.states().to_vec();
    let mut refused = MolecularEuler3d::with_lagrangian_dynamic_thermal_sgs(
        flow,
        source,
        forbidden,
        bytes,
        || false,
    )
    .unwrap();
    for _ in 0..2 {
        assert!(matches!(
            refused.step(Time::new(0.001), MolecularStepControl::default(), || false),
            Err(MolecularError::SubgridThermalDiffusivityRatio { .. }),
        ));
        assert_eq!(refused.flow().states(), old);
        assert_eq!(refused.flow().time(), Time::new(0.0));
        assert_eq!(refused.lagrangian_subgrid_sample_time(), None);
        assert_eq!(refused.lagrangian_subgrid_report(), None);
        assert_eq!(refused.lagrangian_subgrid_heat_report(), None);
    }

    let initial = || cube(4, |_| [1.0, 0.125, -0.0625, 0.03125, 2.0]);
    let flow = initial();
    let source = law(flow.gas(), 0.002, 0.003, 0.0);
    let mut base = MolecularEuler3d::new(flow, source, || false).unwrap();
    let flow = initial();
    let source = law(flow.gas(), 0.002, 0.003, 0.0);
    let mut active = MolecularEuler3d::with_lagrangian_dynamic_thermal_sgs(
        flow,
        source,
        settings(),
        320 * 4usize.pow(3),
        || false,
    )
    .unwrap();
    for _ in 0..3 {
        assert_eq!(
            active
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap(),
            base.step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap(),
        );
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
        for (a, b) in active.flow().states().iter().zip(base.flow().states()) {
            assert_eq!(a.array().map(f64::to_bits), b.array().map(f64::to_bits));
        }
    }
}
