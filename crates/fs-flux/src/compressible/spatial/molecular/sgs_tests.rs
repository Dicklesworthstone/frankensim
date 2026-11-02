use super::*;
use crate::compressible::spatial::molecular::sgs::{SgsModel, SgsUse, SubgridConfig};

fn config(model: SgsModel) -> SubgridConfig {
    SubgridConfig {
        model,
        filter_width_factor: 1.0,
        turbulent_prandtl: 0.8,
        maximum_eddy_to_molecular: 20.0,
        usage: SgsUse::Exploration,
    }
}

fn mixed(x: [f64; 3]) -> [f64; 5] {
    let phase = x.map(|v| v * 2.0 * std::f64::consts::PI);
    [
        1.0,
        0.1 * phase[2].sin() + 0.06 * phase[1].cos(),
        0.1 * phase[2].cos() + 0.06 * phase[0].sin(),
        0.06 * (phase[0].cos() + phase[1].sin()),
        2.0,
    ]
}

#[test]
fn g0_g5_zero_subgrid_controls_preserve_the_complete_molecular_euler_trajectory() {
    for model in [
        SgsModel::Smagorinsky { cs: 0.0 },
        SgsModel::Wale { cw: 0.0 },
        SgsModel::Vreman { cv: 0.0 },
    ] {
        let base = cube(4, mixed);
        let active = cube(4, mixed);
        let base_law = law(base.gas(), 0.002, 0.003, 0.004);
        let active_law = law(active.gas(), 0.002, 0.003, 0.004);
        let mut base = MolecularEuler3d::new(base, base_law, || false).unwrap();
        let mut active =
            MolecularEuler3d::with_sgs(active, active_law, config(model), || false).unwrap();
        assert_eq!(active.subgrid_config(), Some(config(model)));
        for _ in 0..3 {
            let a = base
                .step(Time::new(0.002), MolecularStepControl::default(), || false)
                .unwrap();
            let b = active
                .step(Time::new(0.002), MolecularStepControl::default(), || false)
                .unwrap();
            assert_eq!(a, b);
            for (a, b) in base.flow().states().iter().zip(active.flow().states()) {
                assert_eq!(a.array().map(f64::to_bits), b.array().map(f64::to_bits));
            }
            assert_eq!(base.flow().time(), active.flow().time());
            for side in BoundarySide::ALL {
                assert_eq!(base.boundary_impulses(side), active.boundary_impulses(side));
            }
        }
    }
}

#[test]
fn g0_subgrid_stress_converts_resolved_kinetic_energy_without_losing_total_energy() {
    for model in [
        SgsModel::Smagorinsky { cs: 0.17 },
        SgsModel::Wale { cw: 0.5 },
        SgsModel::Vreman { cv: 0.07 },
    ] {
        let mut base = cube(6, mixed);
        let mut active = cube(6, mixed);
        let original = active.totals();
        let (kinetic0, internal0) = partition(&active);
        let mut base_op =
            MolecularTransport3d::new(&base, law(base.gas(), 0.002, 0.0, 0.0), || false).unwrap();
        let mut op = MolecularTransport3d::with_sgs(
            &active,
            law(active.gas(), 0.002, 0.0, 0.0),
            config(model),
            || false,
        )
        .unwrap();
        base_op
            .advance(
                &mut base,
                Time::new(0.01),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        let report = op
            .advance(
                &mut active,
                Time::new(0.01),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        let (kinetic1, internal1) = partition(&active);
        assert!(kinetic1 < partition(&base).0, "inactive {model:?}");
        assert!(kinetic1 < kinetic0 && internal1 > internal0);
        close(internal1 - internal0, kinetic0 - kinetic1, 2e-13);
        close(active.totals().energy, original.energy, 2e-13);
        for residual in report.conservation_residual.array() {
            close(residual, 0.0, 2e-13);
        }
        assert!(report.maximum_parabolic_number <= 0.40000000000001);
    }
}

fn helical(x: [f64; 3]) -> [f64; 5] {
    let w = 2.0 * std::f64::consts::PI;
    [
        1.0,
        0.1 * (w * x[2]).sin(),
        0.1 * (w * x[2]).cos(),
        0.0,
        2.0 + 0.2 * (w * x[0]).sin(),
    ]
}

#[test]
fn g1_subgrid_heat_flux_uses_cp_prandtl_and_the_physical_filter_width() {
    let n = 8;
    let flow = cube(n, helical);
    let settings = config(SgsModel::Smagorinsky { cs: 0.17 });
    let mut op = MolecularTransport3d::with_sgs(
        &flow,
        law(flow.gas(), 0.002, 0.003, 0.004),
        settings,
        || false,
    )
    .unwrap();
    let rate = spatial(
        &op.gas,
        &op.grid,
        op.boundaries,
        &flow.states,
        &mut op.law,
        &mut op.work.primitive,
        &mut op.work.gradient,
        &mut op.work.flux,
        &mut op.work.binary,
        op.sgs.as_mut(),
        &mut || false,
    )
    .unwrap();
    let h = 1.0 / n as f64;
    // Independent centered derivative of a circularly polarized shear mode:
    // |S|=A sin(kh)/h, constant even though both shear components vary.
    let strain = 0.1 * (2.0 * std::f64::consts::PI * h).sin() / h;
    let mu_t = (0.17 * h).powi(2) * strain;
    // This manufactured source has R=1 and Cp=3.5 (Cv=2.5).
    let conductivity = 0.004 + mu_t * 3.5 / 0.8;
    let c = [1, 1, 1];
    let l = op.work.primitive[op.grid.index([0, 1, 1])];
    let r = op.work.primitive[op.grid.index(c)];
    let x_flux = op.work.flux[0][op.grid.cell_face(0, c)];
    close(
        x_flux.energy,
        -conductivity * (r.temperature - l.temperature) / h,
        3e-14,
    );
    let z_flux = op.work.flux[2][op.grid.cell_face(2, c)];
    let l = op.work.primitive[op.grid.index([1, 1, 0])];
    close(
        z_flux.transverse_momentum[0],
        -(0.002 + mu_t) * (r.velocity[0] - l.velocity[0]) / h,
        2e-15,
    );
    assert!(rate > 0.0);
}

#[test]
fn g0_subgrid_heat_transport_damps_a_resolved_temperature_mode() {
    let mut base = cube(6, helical);
    let mut active = cube(6, helical);
    let mut molecular =
        MolecularTransport3d::new(&base, law(base.gas(), 0.002, 0.0, 0.004), || false).unwrap();
    let mut subgrid = MolecularTransport3d::with_sgs(
        &active,
        law(active.gas(), 0.002, 0.0, 0.004),
        config(SgsModel::Smagorinsky { cs: 0.17 }),
        || false,
    )
    .unwrap();
    molecular
        .advance(
            &mut base,
            Time::new(0.02),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    subgrid
        .advance(
            &mut active,
            Time::new(0.02),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    let amplitude = |flow: &Euler3d| -> f64 {
        flow.states()
            .iter()
            .enumerate()
            .map(|(i, q)| {
                let x = flow.grid().cell_center(i).unwrap()[0];
                q.primitive(flow.gas(), || false).unwrap().temperature
                    * (2.0 * std::f64::consts::PI * x).sin()
            })
            .sum::<f64>()
            * 2.0
            / flow.states().len() as f64
    };
    assert!(amplitude(&active) > 0.0);
    assert!(amplitude(&active) < amplitude(&base));
    close(active.totals().energy, base.totals().energy, 2e-13);
}

#[test]
fn g0_subgrid_admission_and_ratio_refusals_preserve_the_accepted_field() {
    let mut flow = cube(4, mixed);
    let states = flow.states().to_vec();
    let mut settings = config(SgsModel::Wale { cw: 0.5 });
    for pr in [0.0, -1.0, f64::NAN, f64::INFINITY] {
        settings.turbulent_prandtl = pr;
        assert!(settings.validate().is_err());
    }
    settings.turbulent_prandtl = 0.8;
    settings.usage = SgsUse::AcceptanceRender;
    assert!(settings.validate().is_err());
    settings.usage = SgsUse::Exploration;
    settings.maximum_eddy_to_molecular = 0.0;
    let mut op =
        MolecularTransport3d::with_sgs(&flow, law(flow.gas(), 0.002, 0.0, 0.004), settings, || {
            false
        })
        .unwrap();
    let result = op.advance(
        &mut flow,
        Time::new(0.001),
        MolecularControl::default(),
        || false,
    );
    assert!(
        matches!(result, Err(MolecularError::SubgridViscosityRatio { observed, maximum: 0.0, .. }) if observed > 0.0)
    );
    for (a, b) in flow.states().iter().zip(&states) {
        assert_eq!(a.array().map(f64::to_bits), b.array().map(f64::to_bits));
    }
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
fn g0_g5_active_subgrid_coupled_flow_replays_and_closes_all_five_balances() {
    let build = || {
        let flow = cube(4, mixed);
        let source = law(flow.gas(), 0.002, 0.003, 0.004);
        MolecularEuler3d::with_sgs(flow, source, config(SgsModel::Vreman { cv: 0.07 }), || {
            false
        })
        .unwrap()
    };
    let mut a = build();
    let mut b = build();
    for _ in 0..4 {
        let ra = a
            .step(Time::new(0.002), MolecularStepControl::default(), || false)
            .unwrap();
        let rb = b
            .step(Time::new(0.002), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(ra, rb);
        assert_eq!(a.last_step(), Some(&ra));
        for residual in ra.conservation_residual.array() {
            close(residual, 0.0, 3e-13);
        }
        for (qa, qb) in a.flow().states().iter().zip(b.flow().states()) {
            assert_eq!(qa.array().map(f64::to_bits), qb.array().map(f64::to_bits));
        }
    }
}
