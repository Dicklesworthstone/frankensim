use super::*;
use crate::compressible::spatial::molecular::sgs::{
    DynamicAveraging, DynamicScalarConfig, DynamicSgsConfig, DynamicSubgridConfig,
    DynamicThermalSubgridConfig, SgsUse,
};
use fs_lbm::sgs3::dynamic::{DynamicSgsError, PeriodicDynamicSgs};

fn config(averaging: DynamicAveraging) -> DynamicSubgridConfig {
    DynamicSubgridConfig {
        model: DynamicSgsConfig {
            averaging,
            maximum_coefficient_squared: 100.0,
            maximum_eddy_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
        turbulent_prandtl: 0.8,
    }
}

fn mixed(p: [f64; 3]) -> [f64; 5] {
    let [x, y, z] = p.map(|v| std::f64::consts::TAU * v);
    [
        1.2 + 0.2 * x.sin() + 0.1 * y.cos(),
        0.03 * x.sin() + 0.02 * (y + z).sin() + 0.01 * (2.0 * x - y).cos(),
        0.02 * y.cos() + 0.03 * (z + x).sin() + 0.01 * (x - 2.0 * z).sin(),
        0.02 * z.sin() + 0.01 * (x + y).cos() + 0.03 * (2.0 * y - z).sin(),
        2.0 + 0.1 * (x + y + z).sin(),
    ]
}

fn prepare<L: TransportLaw>(op: &mut MolecularTransport3d<L>, flow: &Euler3d) -> f64 {
    spatial(
        &op.gas,
        &op.grid,
        op.boundaries,
        &op.walls,
        &flow.states,
        &mut op.law,
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.flux,
        &mut op.work.binary,
        op.sgs.as_mut(),
        &mut || false,
    )
    .unwrap()
}

#[test]
fn g0_g5_uniform_dynamic_model_preserves_the_original_coupled_trajectory() {
    for averaging in [DynamicAveraging::Domain, DynamicAveraging::LocalBox] {
        let field = || cube(4, |_| [1.0, 0.125, -0.0625, 0.03125, 2.0]);
        let flow = field();
        let source = law(flow.gas(), 0.002, 0.003, 0.004);
        let mut base = MolecularEuler3d::new(flow, source, || false).unwrap();
        let flow = field();
        let source = law(flow.gas(), 0.002, 0.003, 0.004);
        let mut active = MolecularEuler3d::with_dynamic_sgs(
            flow,
            source,
            config(averaging),
            176 * 4usize.pow(3),
            || false,
        )
        .unwrap();
        assert_eq!(active.subgrid_config(), None);
        assert_eq!(active.dynamic_subgrid_config(), Some(config(averaging)));
        for _ in 0..3 {
            let a = base
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            let b = active
                .step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap();
            assert_eq!(a, b);
            assert_eq!(active.flow().time(), base.flow().time());
            for (qa, qb) in active.flow().states().iter().zip(base.flow().states()) {
                assert_eq!(qa.array().map(f64::to_bits), qb.array().map(f64::to_bits));
            }
        }
    }
}

#[test]
fn g0_dynamic_favre_viscosity_drives_the_original_stress_heat_and_conservative_update() {
    let n = 6;
    let mut base = cube(n, mixed);
    let mut flow = cube(n, mixed);
    let settings = config(DynamicAveraging::LocalBox);
    let mut molecular =
        MolecularTransport3d::new(&base, law(base.gas(), 0.002, 0.0, 0.004), || false).unwrap();
    let mut active = MolecularTransport3d::with_dynamic_sgs(
        &flow,
        law(flow.gas(), 0.002, 0.0, 0.004),
        settings,
        176 * n.pow(3),
        || false,
    )
    .unwrap();
    let molecular_rate = prepare(&mut molecular, &base);
    let dynamic_rate = prepare(&mut active, &flow);
    let rho: Vec<_> = active.work.primitive.iter().map(|p| p.density).collect();
    let velocity: Vec<_> = active.work.primitive.iter().map(|p| p.velocity).collect();
    let mut core = PeriodicDynamicSgs::new([n; 3], 1.0 / n as f64, usize::MAX).unwrap();
    let fit = core
        .prepare(settings.model, &rho, &velocity, || false)
        .unwrap();
    assert!(fit.maximum_eddy_viscosity > 0.0);
    assert!(dynamic_rate > molecular_rate);
    let h = 1.0 / n as f64;
    let mut active_faces = 0;
    for i in 0..flow.states.len() {
        let c = flow.grid.coordinates(i);
        let mut left = c;
        left[0] = (c[0] + n - 1) % n;
        let left = flow.grid.index(left);
        let pl = active.work.primitive[left];
        let pr = active.work.primitive[i];
        let mu_t =
            0.5 * rho[left] * core.viscosities()[left] + 0.5 * rho[i] * core.viscosities()[i];
        if mu_t == 0.0 {
            continue;
        }
        active_faces += 1;
        let face = flow.grid.cell_face(0, c);
        let a = active.work.flux[0][face];
        let b = molecular.work.flux[0][face];
        // Remove the already-conservative stress work: only Cp/Pr heat remains.
        let mean = std::array::from_fn::<_, 3, _>(|k| 0.5 * (pl.velocity[k] + pr.velocity[k]));
        let heat = |f: SpatialFaceFlux| {
            f.energy
                - mean[0] * f.normal_momentum
                - mean[1] * f.transverse_momentum[0]
                - mean[2] * f.transverse_momentum[1]
        };
        close(
            heat(a) - heat(b),
            -mu_t * 3.5 / settings.turbulent_prandtl * (pr.temperature - pl.temperature) / h,
            2e-14,
        );
        close(
            a.transverse_momentum[0] - b.transverse_momentum[0],
            -mu_t
                * ((pr.velocity[1] - pl.velocity[1]) / h
                    + 0.5 * active.work.gradient[left][0][1]
                    + 0.5 * active.work.gradient[i][0][1]),
            2e-14,
        );
    }
    assert!(active_faces > 0);
    let original = flow.totals();
    let (kinetic0, internal0) = partition(&flow);
    molecular
        .advance(
            &mut base,
            Time::new(0.002),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let report = active
        .advance(
            &mut flow,
            Time::new(0.002),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let (kinetic1, internal1) = partition(&flow);
    assert!(kinetic1 < partition(&base).0);
    close(internal1 - internal0, kinetic0 - kinetic1, 3e-13);
    close(flow.totals().energy, original.energy, 3e-13);
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 3e-13);
    }
    assert!(report.maximum_parabolic_number <= 0.40000000000001);
}

#[test]
fn g0_dynamic_grid_budget_and_coefficient_admission_refuse_without_publishing() {
    let settings = config(DynamicAveraging::LocalBox);
    for (shape, spacing, boundary) in [
        ([2, 4, 4], [0.25; 3], PERIODIC),
        ([4; 3], [0.25, 0.5, 0.25], PERIODIC),
        ([4; 3], [0.25; 3], WALLS),
    ] {
        let [nx, ny, nz] = shape;
        let [dx, dy, dz] = spacing;
        let flow = field(
            gas(0.0),
            CartesianGrid::new(nx, ny, nz, dx, dy, dz).unwrap(),
            boundary,
            mixed,
        );
        assert!(
            MolecularTransport3d::with_dynamic_sgs(
                &flow,
                law(flow.gas(), 0.002, 0.0, 0.004),
                settings,
                usize::MAX,
                || false,
            )
            .is_err()
        );
    }
    let mut flow = cube(6, mixed);
    let bytes =
        MolecularTransport3d::<ConstantTransport>::required_dynamic_sgs_workspace_bytes(&flow)
            .unwrap();
    assert_eq!(bytes, 176 * 6usize.pow(3));
    assert!(matches!(
        MolecularTransport3d::with_dynamic_sgs(
            &flow,
            law(flow.gas(), 0.002, 0.0, 0.004),
            settings,
            bytes - 1,
            || false,
        ),
        Err(MolecularError::DynamicSubgrid(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    for pr in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        assert!(
            DynamicSubgridConfig {
                turbulent_prandtl: pr,
                ..settings
            }
            .validate()
            .is_err()
        );
    }
    let mut forbidden = settings;
    forbidden.model.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
    let mut forbidden = settings;
    forbidden.model.maximum_coefficient_squared = 0.0;
    let mut op = MolecularTransport3d::with_dynamic_sgs(
        &flow,
        law(flow.gas(), 0.002, 0.0, 0.004),
        forbidden,
        bytes,
        || false,
    )
    .unwrap();
    let state = flow.states().to_vec();
    assert!(matches!(
        op.advance(
            &mut flow,
            Time::new(0.001),
            MolecularControl::default(),
            || false
        ),
        Err(MolecularError::DynamicSubgrid(
            DynamicSgsError::Coefficient { .. }
        ))
    ));
    assert_eq!(flow.states(), state);
    assert_eq!(flow.time().value(), 0.0);
    for side in BoundarySide::ALL {
        assert!(
            op.boundary_impulses(side)
                .iter()
                .all(|q| *q == Totals::default())
        );
    }
}

#[test]
fn g4_g5_dynamic_filter_and_coupled_publication_cancel_atomically_and_replay() {
    let build = || {
        let flow = cube(4, mixed);
        let source = law(flow.gas(), 0.002, 0.0, 0.004);
        MolecularEuler3d::with_dynamic_sgs(
            flow,
            source,
            config(DynamicAveraging::LocalBox),
            176 * 4usize.pow(3),
            || false,
        )
        .unwrap()
    };
    let mut probe = build();
    probe
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    let mut count = 0;
    let reference = probe
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            count += 1;
            false
        })
        .unwrap();
    for cutoff in [1, count / 2, count] {
        let mut flow = build();
        flow.step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let before = flow.flow().states().to_vec();
        let time = flow.flow().time();
        let previous = flow.last_step().copied();
        let mut polls = 0;
        assert_eq!(
            flow.step(Time::new(0.001), MolecularStepControl::default(), || {
                polls += 1;
                polls == cutoff
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), before);
        assert_eq!(flow.flow().time(), time);
        assert_eq!(flow.last_step().copied(), previous);
        let report = flow
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(report, reference);
        for (a, b) in flow.flow().states().iter().zip(probe.flow().states()) {
            assert_eq!(a.array().map(f64::to_bits), b.array().map(f64::to_bits));
        }
    }
}

fn thermal_config() -> DynamicThermalSubgridConfig {
    DynamicThermalSubgridConfig {
        momentum: config(DynamicAveraging::LocalBox).model,
        heat: DynamicScalarConfig {
            averaging: DynamicAveraging::LocalBox,
            maximum_coefficient: 100.0,
            maximum_diffusivity_to_molecular: 1000.0,
            usage: SgsUse::Exploration,
        },
    }
}

// The momentum Germano contraction vanishes for this parallel shear, while
// the scalar's cross correlation has a nonzero, spatially varying fit.
fn shear_temperature(p: [f64; 3]) -> [f64; 5] {
    let [x, y, _] = p.map(|v| core::f64::consts::TAU * v);
    [1.0, 0.03 * y.sin(), 0.0, 0.0, 2.0 + 0.2 * (x + y).sin()]
}

#[test]
fn g0_g5_dynamic_heat_zero_fits_preserve_the_original_coupled_trajectory() {
    let field = || cube(4, |_| [1.0, 0.125, -0.0625, 0.03125, 2.0]);
    let flow = field();
    let source = law(flow.gas(), 0.002, 0.003, 0.0);
    let mut base = MolecularEuler3d::new(flow, source, || false).unwrap();
    let flow = field();
    let source = law(flow.gas(), 0.002, 0.003, 0.0);
    let settings = thermal_config();
    let mut active = MolecularEuler3d::with_dynamic_thermal_sgs(
        flow,
        source,
        settings,
        256 * 4usize.pow(3),
        || false,
    )
    .unwrap();
    assert_eq!(active.dynamic_thermal_subgrid_config(), Some(settings));
    assert_eq!(active.dynamic_subgrid_config(), None);
    assert_eq!(active.subgrid_config(), None);
    for _ in 0..3 {
        let a = base
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        let b = active
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(a, b);
        assert_eq!(active.flow().states(), base.flow().states());
        assert_eq!(active.flow().time(), base.flow().time());
    }
}

#[test]
fn g1_dynamic_heat_with_zero_modeled_stress_reaches_conservative_flux_and_stability() {
    use fs_lbm::sgs3::dynamic::scalar::PeriodicDynamicScalar;
    let n: usize = 9;
    let mut flow = cube(n, shear_temperature);
    let mut base = cube(n, shear_temperature);
    let settings = thermal_config();
    let mut molecular =
        MolecularTransport3d::new(&base, law(base.gas(), 0.0, 0.0, 0.001), || false).unwrap();
    let mut active = MolecularTransport3d::with_dynamic_thermal_sgs(
        &flow,
        law(flow.gas(), 0.0, 0.0, 0.001),
        settings,
        256 * n.pow(3),
        || false,
    )
    .unwrap();
    let base_rate = prepare(&mut molecular, &base);
    let heat_rate = prepare(&mut active, &flow);
    assert!(heat_rate > base_rate);
    let rho: Vec<_> = active.work.primitive.iter().map(|p| p.density).collect();
    let velocity: Vec<_> = active.work.primitive.iter().map(|p| p.velocity).collect();
    let temperature: Vec<_> = active
        .work
        .primitive
        .iter()
        .map(|p| p.temperature)
        .collect();
    let h = 1.0 / n as f64;
    let mut momentum = PeriodicDynamicSgs::new([n; 3], h, usize::MAX).unwrap();
    let mut heat = PeriodicDynamicScalar::new([n; 3], usize::MAX).unwrap();
    let (stress_report, heat_report) = momentum
        .prepare_with_scalar(
            settings.momentum,
            settings.heat,
            &rho,
            &velocity,
            &temperature,
            &mut heat,
            || false,
        )
        .unwrap();
    assert_eq!(stress_report.maximum_eddy_viscosity, 0.0);
    assert!(heat_report.maximum_diffusivity > 0.0);
    let mut active_faces = 0;
    for i in 0..flow.states.len() {
        let c = flow.grid.coordinates(i);
        let mut left = c;
        left[0] = (c[0] + n - 1) % n;
        let left = flow.grid.index(left);
        let rho_alpha =
            0.5 * rho[left] * heat.diffusivities()[left] + 0.5 * rho[i] * heat.diffusivities()[i];
        let face = flow.grid.cell_face(0, c);
        let a = active.work.flux[0][face];
        let b = molecular.work.flux[0][face];
        assert_eq!(a.normal_momentum, b.normal_momentum);
        assert_eq!(a.transverse_momentum, b.transverse_momentum);
        close(
            a.energy - b.energy,
            -rho_alpha * 3.5 * (temperature[i] - temperature[left]) / h,
            2e-14,
        );
        active_faces += usize::from(rho_alpha > 0.0);
    }
    assert!(active_faces > 0);
    let before = flow.totals();
    let variance = |f: &Euler3d| -> f64 {
        f.states()
            .iter()
            .map(|q| (q.primitive(f.gas(), || false).unwrap().temperature - 2.0).powi(2))
            .sum::<f64>()
            / f.states().len() as f64
    };
    molecular
        .advance(
            &mut base,
            Time::new(0.001),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let report = active
        .advance(
            &mut flow,
            Time::new(0.001),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    assert!(variance(&flow) < variance(&base));
    close(flow.totals().energy, before.energy, 3e-13);
    for residual in report.conservation_residual.array() {
        close(residual, 0.0, 3e-13);
    }
    assert!(report.maximum_parabolic_number <= 0.40000000000001);
}

#[test]
fn g0_g4_dynamic_heat_budget_and_independent_thermal_ratio_refuse_atomically() {
    let mut flow = cube(9, shear_temperature);
    let settings = thermal_config();
    let bytes =
        MolecularTransport3d::<ConstantTransport>::required_dynamic_thermal_sgs_workspace_bytes(
            &flow,
        )
        .unwrap();
    assert_eq!(bytes, 256 * 9usize.pow(3));
    assert!(matches!(
        MolecularTransport3d::with_dynamic_thermal_sgs(
            &flow,
            law(flow.gas(), 0.0, 0.0, 0.001),
            settings,
            bytes - 1,
            || false,
        ),
        Err(MolecularError::DynamicSubgrid(
            DynamicSgsError::WorkspaceBytes { .. }
        ))
    ));
    let before = flow.states().to_vec();
    for molecular_k in [0.0, 0.001] {
        let mut forbidden = settings;
        forbidden.heat.maximum_diffusivity_to_molecular = 0.0;
        let mut active = MolecularTransport3d::with_dynamic_thermal_sgs(
            &flow,
            law(flow.gas(), 0.0, 0.0, molecular_k),
            forbidden,
            bytes,
            || false,
        )
        .unwrap();
        assert!(matches!(
            active.advance(
                &mut flow,
                Time::new(0.001),
                MolecularControl::default(),
                || false,
            ),
            Err(MolecularError::SubgridThermalDiffusivityRatio { .. })
        ));
        assert_eq!(flow.states(), before);
        assert_eq!(flow.time().value(), 0.0);
        for side in BoundarySide::ALL {
            assert!(
                active
                    .boundary_impulses(side)
                    .iter()
                    .all(|q| *q == Totals::default())
            );
        }
    }
    let mut forbidden = settings;
    forbidden.heat.usage = SgsUse::AcceptanceRender;
    assert!(forbidden.validate().is_err());
}

#[test]
fn g4_g5_dynamic_heat_coupled_cancellation_replays_the_original_transaction() {
    let build = || {
        let flow = cube(5, mixed);
        let source = law(flow.gas(), 0.002, 0.0, 0.004);
        MolecularEuler3d::with_dynamic_thermal_sgs(
            flow,
            source,
            thermal_config(),
            256 * 5usize.pow(3),
            || false,
        )
        .unwrap()
    };
    let mut reference = build();
    let mut polls = 0;
    let report = reference
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            polls += 1;
            false
        })
        .unwrap();
    for cutoff in [1, polls / 2, polls] {
        let mut flow = build();
        let before = flow.flow().states().to_vec();
        let mut seen = 0;
        assert_eq!(
            flow.step(Time::new(0.001), MolecularStepControl::default(), || {
                seen += 1;
                seen == cutoff
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), before);
        assert_eq!(flow.flow().time().value(), 0.0);
        assert_eq!(flow.last_step(), None);
        assert_eq!(
            flow.step(Time::new(0.001), MolecularStepControl::default(), || false)
                .unwrap(),
            report
        );
        assert_eq!(flow.flow().states(), reference.flow().states());
    }
}
