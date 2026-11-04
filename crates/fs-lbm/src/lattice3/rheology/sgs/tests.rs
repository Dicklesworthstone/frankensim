//! Source/SGS composition on the original collision and force-history owners.

use super::*;
use super::super::tests::{
    bits, collision, fluid, grid27, storage_addresses, with_cx, with_gate,
};
use crate::d3q27::{CollisionConfig27, CollisionModel27, D3q27, Equilibrium27, equilibrium27};
use crate::sgs3::dynamic::DynamicAveraging;
use crate::sgs3::{SgsError, SgsModel, SgsUse};
use std::{cell::Cell, rc::Rc};

fn local(cs: f64) -> SgsConfig {
    SgsConfig {
        model: SgsModel::Smagorinsky { cs },
        filter_width_lu: 1.0,
        maximum_eddy_to_molecular: 100.0,
        usage: SgsUse::Exploration,
    }
}

fn viscosity(point: RheologyPoint3) -> f64 {
    0.004 + 0.01 / (1.0 + point.shear_rate_per_s)
        + 0.0004 * point.coordinates[0] as f64 + 0.02 * point.time_s
}

fn poiseuille() -> Grid3<D3q27, 27> {
    grid27([4, 6, 3], 3, [0.0; 3], [true, false, true], |p| {
        let y = p[1] as f64 + 0.5;
        equilibrium27(1.0, [0.001 * y * (6.0 - y), 0.0, 0.0], Equilibrium27::Factorized)
    })
}

fn random(cell: usize, salt: u64) -> f64 {
    let mut h = (cell as u64 + 1).wrapping_mul(0x9e37_79b9_7f4a_7c15) ^ salt;
    h ^= h >> 30;
    h = h.wrapping_mul(0xbf58_476d_1ce4_e5b9);
    h ^= h >> 27;
    h = h.wrapping_mul(0x94d0_49bb_1331_11eb);
    h ^= h >> 31;
    (h >> 11) as f64 / (1u64 << 53) as f64 * 2.0 - 1.0
}

fn fluctuating_velocity(p: [usize; 3]) -> [f64; 3] {
    let cell = (p[2] * 4 + p[1]) * 5 + p[0];
    let z = p[2] as f64 + 0.5;
    let envelope = 4.0 * z * (5.0 - z) / 25.0;
    core::array::from_fn(|a| 0.02 * envelope * random(cell, 31 * a as u64 + 7))
}

fn fluctuating(channel: bool) -> Grid3<D3q27, 27> {
    grid27([5, 4, 5], 3, [0.0; 3], [true, true, !channel], |p| {
        equilibrium27(1.0, fluctuating_velocity(p), Equilibrium27::Factorized)
    })
}

enum DynamicScratch {
    Periodic(DynamicSgsWorkspace3<CollisionConfig27>),
    Channel(ChannelDynamicSgsWorkspace3<CollisionConfig27>),
}

impl DynamicScratch {
    fn new(grid: &Grid3<D3q27, 27>, channel: bool) -> Self {
        if channel {
            Self::Channel(grid.channel_dynamic_sgs_workspace(usize::MAX).unwrap())
        } else {
            Self::Periodic(grid.dynamic_sgs_workspace(usize::MAX).unwrap())
        }
    }
    fn selection(&mut self) -> RheologySgs3<'_, CollisionConfig27> {
        match self {
            Self::Periodic(workspace) => RheologySgs3::PeriodicDynamic {
                config: DynamicSgsConfig {
                    averaging: DynamicAveraging::LocalBox,
                    maximum_coefficient_squared: 10.0,
                    maximum_eddy_to_molecular: 100.0,
                    usage: SgsUse::Exploration,
                },
                workspace,
            },
            Self::Channel(workspace) => RheologySgs3::ChannelDynamic {
                config: ChannelDynamicSgsConfig {
                    maximum_coefficient_squared: 10.0,
                    maximum_eddy_to_molecular: 100.0,
                    usage: SgsUse::Exploration,
                },
                workspace,
            },
        }
    }
    fn viscosities(&self) -> &[f64] {
        match self {
            Self::Periodic(workspace) => workspace.eddy_viscosities(),
            Self::Channel(workspace) => workspace.eddy_viscosities(),
        }
    }
}

#[test]
fn source_shear_local_sgs_and_applied_rates_match_analytic_poiseuille_then_refresh() {
    // Neither Clone nor Send: the original producer remains on the serial side
    // even when the second native step uses a worker pool.
    struct Source {
        state: FluidTransportState,
        reads: Rc<Cell<usize>>,
    }
    impl FluidTransport for Source {
        fn transport_state(&self) -> Result<FluidTransportState, MaterialError> {
            self.reads.set(self.reads.get() + 1);
            Ok(self.state)
        }
    }
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    let mut actual = poiseuille();
    let mut native = poiseuille();
    let mut sources = actual.rheology_workspace::<Source>().unwrap();
    let addresses = storage_addresses(&sources);
    let mut sgs = actual.sgs_workspace().unwrap();
    let collision_pointer = sgs.collisions().as_ptr();
    let queries = Cell::new(0);
    let reads = Rc::new(Cell::new(0));
    let resolve = |point, _: &FluidTransportState| {
        queries.set(queries.get() + 1);
        Ok(Source { state: fluid(viscosity(point)), reads: Rc::clone(&reads) })
    };
    let mut expected = Vec::new();
    let mut molecular_power = 0.0;
    let mut modeled_power = 0.0;
    for cell in 0..72 {
        let p = coordinates(cell, [4, 6, 3]);
        let shear = 0.001 * (6.0 - 2.0 * (p[1] as f64 + 0.5));
        let point = RheologyPoint3 {
            cell, coordinates: p, time_s: 0.0, shear_rate_per_s: shear.abs() / 0.01,
        };
        let nu = 10.0 * viscosity(point);
        let eddy = 0.2_f64.powi(2) * shear.abs();
        expected.push(D3q27::with_shear_relaxation(collision(3), 0.5 + 3.0 * (nu + eddy)).unwrap());
        molecular_power += nu * shear * shear;
        modeled_power += eddy * shear * shear;
    }
    with_cx(|cx| {
        let step = actual.step_serial_with_material_rheology_and_sgs(
            cx, &reference, &mut sources,
            RheologySgs3::Local { config: local(0.2), workspace: &mut sgs }, &resolve,
        ).unwrap();
        let RheologySgsReport3::Local(report) = step.sgs else { panic!("wrong closure"); };
        assert!((report.resolved_dissipation_lu - molecular_power).abs() < 1e-12 * molecular_power);
        assert!((report.sgs_dissipation_lu - modeled_power).abs() < 1e-12 * modeled_power);
        for (cell, sample) in step.rheology.samples.iter().enumerate() {
            let sample = sample.as_ref().unwrap();
            let shear = (0.001 * (5.0 - 2.0 * sample.point().coordinates[1] as f64)).abs();
            assert!((sample.point().shear_rate_per_s - shear / 0.01).abs() < 2e-13);
            let applied = step.rheology.collisions[cell];
            let molecular = D3q27::with_shear_relaxation(collision(3), sample.shear_relaxation_time()).unwrap();
            assert!((D3q27::molecular_viscosity(applied).unwrap()
                - D3q27::molecular_viscosity(molecular).unwrap() - 0.04 * shear).abs() < 2e-16);
            let (CollisionModel27::Cumulant { rates: actual }, CollisionModel27::Cumulant { rates: original })
                = (applied.model, collision(3).model) else { panic!("collision changed"); };
            assert_eq!(actual.omega[1..], original.omega[1..]);
        }
        assert_eq!(queries.get(), 72);
        assert_eq!(reads.get(), 72);
        native.step_serial_with_collisions(cx, &expected).unwrap();
        for (a, b) in actual.populations().iter().flatten().zip(native.populations().iter().flatten()) {
            assert!((a - b).abs() < 3e-15);
        }
        // A later zero model must start from new molecular responses, never
        // from the previous report's already augmented applied rates.
        let step = actual.step_pooled_with_material_rheology_and_sgs(
            &fs_exec::TilePool::for_host(2, 93), &CancelGate::new_clock_free(), &reference,
            &mut sources, RheologySgs3::Local { config: local(0.0), workspace: &mut sgs }, &resolve,
        ).unwrap();
        for (sample, &applied) in step.rheology.samples.iter().zip(&step.rheology.collisions) {
            let sample = sample.as_ref().unwrap();
            assert_eq!(sample.point().time_s, 0.01);
            assert_eq!(applied, D3q27::with_shear_relaxation(collision(3), sample.shear_relaxation_time()).unwrap());
        }
        assert_eq!(step.rheology.sampled_at_step, 1);
        assert_eq!(queries.get(), 144);
        assert_eq!(reads.get(), 144);
        assert_eq!(actual.collision, collision(3));
    });
    assert_eq!(addresses, storage_addresses(&sources));
    assert_eq!(collision_pointer, sgs.collisions().as_ptr());
}

#[test]
fn both_dynamic_geometries_use_actual_molecular_sources_and_native_transport() {
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    for channel in [false, true] {
        let mut serial = fluctuating(channel);
        let mut native = fluctuating(channel);
        let mut pooled = fluctuating(channel);
        let mut a = serial.rheology_workspace().unwrap();
        let mut b = pooled.rheology_workspace().unwrap();
        let mut sa = DynamicScratch::new(&serial, channel);
        let mut sb = DynamicScratch::new(&pooled, channel);
        with_cx(|cx| {
            let first = serial.step_serial_with_material_rheology_and_sgs(
                cx, &reference, &mut a, sa.selection(), |point, _| Ok(fluid(viscosity(point))),
            ).unwrap();
            let RheologySgsReport3::Dynamic(report) = first.sgs else { panic!("wrong closure"); };
            assert!(report.fit.maximum_eddy_viscosity > 0.0);
            for (cell, sample) in first.rheology.samples.iter().enumerate() {
                let sample = sample.as_ref().unwrap();
                let molecular = D3q27::with_shear_relaxation(collision(3), sample.shear_relaxation_time()).unwrap();
                assert!((D3q27::molecular_viscosity(first.rheology.collisions[cell]).unwrap()
                    - D3q27::molecular_viscosity(molecular).unwrap() - sa.viscosities()[cell]).abs() < 2e-16);
                assert_eq!(sample.transport().dynamic_viscosity_pa_s(), viscosity(sample.point()));
            }
            native.step_serial_with_collisions(cx, &first.rheology.collisions).unwrap();
            let second = pooled.step_pooled_with_material_rheology_and_sgs(
                &fs_exec::TilePool::for_host(2, 93), &CancelGate::new_clock_free(), &reference,
                &mut b, sb.selection(), |point, _| Ok(fluid(viscosity(point))),
            ).unwrap();
            assert_eq!(first.sgs, second.sgs);
            assert_eq!(first.rheology.collisions, second.rheology.collisions);
            assert_eq!(bits(serial.populations()), bits(native.populations()));
            assert_eq!(bits(serial.populations()), bits(pooled.populations()));
            assert_eq!(serial.completed_steps(), 1);
        });
    }
}

#[test]
fn force_endpoints_do_not_reinterpret_incoming_rheology_or_dynamic_fit() {
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    with_cx(|cx| {
        let current: Vec<_> = (0..100).map(|cell| {
            let p = coordinates(cell, [5, 4, 5]);
            [0.001 * (p[1] as f64 - 1.5), 0.0, 0.0]
        }).collect();
        let make = || {
            let grid = grid27([5, 4, 5], 3, [0.0; 3], [true, true, false], |p| {
                let cell = (p[2] * 4 + p[1]) * 5 + p[0];
                let u = fluctuating_velocity(p);
                equilibrium27(1.0, core::array::from_fn(|a| u[a] - 0.5 * current[cell][a]), Equilibrium27::Factorized)
            });
            ForcedGrid3::try_new(grid, cx, current.clone()).unwrap()
        };
        let mut serial = make();
        let mut pooled = make();
        let mut unforced = fluctuating(true);
        let mut a = serial.rheology_workspace().unwrap();
        let mut b = pooled.rheology_workspace().unwrap();
        let mut plain = unforced.rheology_workspace().unwrap();
        let mut sa = DynamicScratch::new(&serial.grid, true);
        let mut sb = DynamicScratch::new(&pooled.grid, true);
        let mut sp = DynamicScratch::new(&unforced, true);
        let end_a = vec![[0.0; 3]; 100];
        let end_b: Vec<_> = current.iter().map(|f| [-2.0 * f[0], 0.001, 0.0]).collect();
        let first = serial.advance_serial_with_material_rheology_and_sgs(
            cx, &reference, end_a.clone(), &mut a, sa.selection(),
            |point, _| Ok(fluid(viscosity(point))),
        ).unwrap();
        let second = pooled.advance_pooled_with_material_rheology_and_sgs(
            &fs_exec::TilePool::for_host(2, 93), &CancelGate::new_clock_free(),
            &reference, end_b.clone(), &mut b, sb.selection(),
            |point, _| Ok(fluid(viscosity(point))),
        ).unwrap();
        let physical = unforced.step_serial_with_material_rheology_and_sgs(
            cx, &reference, &mut plain, sp.selection(), |point, _| Ok(fluid(viscosity(point))),
        ).unwrap();
        assert_eq!(first.sgs, second.sgs);
        assert_eq!(first.rheology.collisions, second.rheology.collisions);
        for (a, b) in first.rheology.samples.iter().zip(&physical.rheology.samples) {
            let (a, b) = (a.as_ref().unwrap(), b.as_ref().unwrap());
            assert!((a.point().shear_rate_per_s - b.point().shear_rate_per_s).abs() < 2e-12);
        }
        assert_eq!(serial.forces(), end_a);
        assert_eq!(pooled.forces(), end_b);
        assert_ne!(bits(serial.populations()), bits(pooled.populations()));
        let before = bits(serial.populations());
        let accepted = a.report().unwrap().collisions.clone();
        let queries = Cell::new(0);
        assert!(matches!(serial.advance_serial_with_material_rheology_and_sgs(
            cx, &reference, vec![[1.0, 0.0, 0.0]; 100], &mut a, sa.selection(), |point, _| {
                queries.set(queries.get() + 1);
                Ok(fluid(viscosity(point)))
            },
        ), Err(RheologySgsError3::Sgs(SgsStepError::Grid(GridError3::Mach { .. })))));
        assert_eq!(queries.get(), 100);
        assert_eq!(bits(serial.populations()), before);
        assert_eq!(serial.forces(), end_a);
        assert_eq!(serial.completed_steps(), 1);
        assert_eq!(a.report().unwrap().collisions, accepted);
        assert_eq!(a.report().unwrap().sampled_at_step, 0);
    });
}

#[test]
fn molecular_cap_source_refusal_and_cancellation_keep_the_last_applied_epoch() {
    let material = fluid(0.01);
    let reference = plan_scaling_from_material(&material, 0.001, 0.01, 0.005).unwrap();
    let mut flow = poiseuille();
    let mut sources = flow.rheology_workspace().unwrap();
    let mut sgs = flow.sgs_workspace().unwrap();
    with_cx(|cx| {
        flow.step_serial_with_material_rheology_and_sgs(
            cx, &reference, &mut sources,
            RheologySgs3::Local { config: local(0.0), workspace: &mut sgs }, |_, _| Ok(fluid(0.01)),
        ).unwrap();
        let before = bits(flow.populations());
        let accepted = sources.report().unwrap().collisions.clone();
        let capped = SgsConfig { maximum_eddy_to_molecular: 0.005, ..local(0.2) };
        // Constructor viscosity would admit this eddy field. The actual thinner
        // source liquid must refuse its larger eddy/molecular ratio.
        flow.prepare_sgs(cx, capped, None, &mut sgs).unwrap();
        assert!(matches!(flow.step_serial_with_material_rheology_and_sgs(
            cx, &reference, &mut sources,
            RheologySgs3::Local { config: capped, workspace: &mut sgs }, |_, _| Ok(fluid(0.002)),
        ), Err(RheologySgsError3::Sgs(SgsStepError::Closure(SgsError::ViscosityRatioExceeded { .. })))));
        assert!(flow.step_serial_with_material_rheology_and_sgs(
            cx, &reference, &mut sources,
            RheologySgs3::Local { config: local(0.2), workspace: &mut sgs }, |point, _| {
                FluidTransportState::try_new(1000.0, 4000.0, 0.5, if point.cell == 71 { -1.0 } else { 0.01 })
            },
        ).is_err());
        let cancelled = CancelGate::new_clock_free();
        cancelled.request();
        with_gate(&cancelled, |cx| {
            assert!(flow.step_serial_with_material_rheology_and_sgs(
                cx, &reference, &mut sources,
                RheologySgs3::Local { config: local(0.2), workspace: &mut sgs }, |_, _| Ok(fluid(0.01)),
            ).is_err());
        });
        assert_eq!(bits(flow.populations()), before);
        assert_eq!(flow.completed_steps(), 1);
        assert_eq!(sources.report().unwrap().collisions, accepted);
        assert_eq!(sources.report().unwrap().sampled_at_step, 0);
        let mut native = poiseuille();
        native.f.clone_from(&flow.f);
        native.completed_steps = flow.completed_steps;
        let mut native_sources = native.rheology_workspace().unwrap();
        let retry = flow.step_serial_with_material_rheology_and_sgs(
            cx, &reference, &mut sources,
            RheologySgs3::Local { config: local(0.0), workspace: &mut sgs }, |_, _| Ok(fluid(0.011)),
        ).unwrap();
        let original = native.step_serial_with_material_rheology_workspace(
            cx, &reference, &mut native_sources, |_, _| Ok(fluid(0.011)),
        ).unwrap();
        assert_eq!(retry.rheology.collisions, original.collisions);
        assert_eq!(bits(flow.populations()), bits(native.populations()));
        assert_eq!(retry.rheology.sampled_at_step, 1);
        assert_eq!(flow.completed_steps(), 2);
    });
}
