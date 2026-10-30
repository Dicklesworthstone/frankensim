use super::*;
use crate::compressible::molecular::{
    MolecularControl, MolecularError, MolecularEuler2d, MolecularStepControl, MolecularTransport2d,
};
use fs_gas::multicomponent::dissipation::{
    ConstantTransport, KineticTransport, TransportCoefficients, TransportDomain,
};
use fs_qty::DynViscosity;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::{
    KineticGasTransportV1, KineticTransportDomainV1, LennardJonesSpeciesV1,
    ThermalConductivityQuantityV1,
};
use fs_thermochem::{ChargeVector, ElementId, ElementalMatrix};

fn bank(g: &FrozenGas) -> IdealGasEquilibriumV1 {
    let models = g.source().components().to_vec();
    let ids: Vec<_> = models.iter().map(|m| m.species().clone()).collect();
    IdealGasEquilibriumV1::new(
        ElementalMatrix::new(
            vec![ElementId::new("H").unwrap()],
            ids.clone(),
            vec![vec![1; ids.len()]],
        )
        .unwrap(),
        ChargeVector::new(ids.clone(), vec![0; ids.len()]).unwrap(),
        models,
    )
    .unwrap()
}
fn law(flow: &Euler2d, mu: f64, bulk: f64, k: f64) -> ConstantTransport {
    ConstantTransport::new(
        bank(flow.gas()),
        "manufactured-planar-transport",
        1,
        TransportDomain {
            minimum_temperature_k: 0.001,
            maximum_temperature_k: 50.0,
            minimum_pressure_pa: 1e-8,
            maximum_pressure_pa: 1e8,
        },
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(mu),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(k),
        },
        vec![0.0],
    )
    .unwrap()
}
fn planar(n: usize, initial: impl Fn([f64; 2]) -> [f64; 4]) -> Euler2d {
    field(
        gas(0.0),
        CartesianGrid::new(n, n, 1.0 / n as f64, 1.0 / n as f64, 0.7).unwrap(),
        PERIODIC,
        initial,
    )
}
fn kinetic(flow: &Euler2d, diameter: f64) -> KineticTransport {
    let species = flow.gas().source().components()[0].species().clone();
    let model = KineticGasTransportV1::new(
        bank(flow.gas()),
        "manufactured-planar-collision",
        1,
        vec![
            LennardJonesSpeciesV1::new(species, Length::new(diameter), Temperature::new(1.0))
                .unwrap(),
        ],
        KineticTransportDomainV1 {
            minimum_temperature_k: 0.5,
            maximum_temperature_k: 40.0,
            minimum_pressure_pa: 0.01,
            maximum_pressure_pa: 1e4,
        },
        || false,
    )
    .unwrap();
    KineticTransport::new(model, DynViscosity::new(0.0)).unwrap()
}
fn energy_partition(flow: &Euler2d) -> (f64, f64) {
    let kinetic: f64 = flow
        .states()
        .iter()
        .map(|q| {
            0.5 * (q.momentum[0].powi(2) + q.momentum[1].powi(2)) / q.density
                * flow.grid().cell_volume()
        })
        .sum();
    (kinetic, flow.totals().energy - kinetic)
}

#[test]
fn g1_two_dimensional_fourier_mode_refines_with_original_cv_and_conserves_energy() {
    let w = 2.0 * std::f64::consts::PI;
    let mut errors = Vec::new();
    for n in [12, 24] {
        let shape = |x: [f64; 2]| (w * x[0]).sin() * (w * x[1]).cos();
        let mut flow = planar(n, |x| [1.0, 0.0, 0.0, 2.0 + 0.2 * shape(x)]);
        let initial = flow.totals();
        let mut operator =
            MolecularTransport2d::new(&flow, law(&flow, 0.0, 0.0, 0.05), || false).unwrap();
        let end = 0.08;
        let report = operator
            .advance(
                &mut flow,
                Time::new(end),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
        // R_specific=1 and Cp/R=3.5 give Cv=2.5; the 2D mode has |k|²=2w².
        let amplitude = 0.2 * (-0.05 / 2.5 * 2.0 * w * w * end).exp();
        let error: f64 = flow
            .states()
            .iter()
            .enumerate()
            .map(|(i, q)| {
                assert_eq!(q.density, 1.0);
                assert_eq!(q.momentum, [0.0, 0.0]);
                (q.primitive(flow.gas(), || false).unwrap().temperature
                    - (2.0 + amplitude * shape(flow.grid().cell_center(i).unwrap())))
                .abs()
            })
            .sum::<f64>()
            / (n * n) as f64;
        errors.push(error);
        assert_eq!(flow.time().value(), 0.0);
        assert_eq!(report.elapsed.value(), end);
        assert!(report.maximum_parabolic_number <= 0.40000000000001);
        close(flow.totals().energy, initial.energy, 2e-13);
        for r in report.conservation_residual.array() {
            close(r, 0.0, 2e-13);
        }
    }
    assert!(errors[0] < 4e-4, "Fourier errors: {errors:?}");
    assert!(errors[1] < 0.3 * errors[0], "Fourier errors: {errors:?}");
}

#[test]
fn g1_oblique_longitudinal_and_transverse_stress_modes_include_cross_derivatives() {
    let w = 2.0 * std::f64::consts::PI;
    let end = 0.04;
    let (mu, bulk) = (0.05, 0.03);
    for longitudinal in [false, true] {
        let mut errors = Vec::new();
        for n in [12, 24] {
            let sign = if longitudinal { 1.0 } else { -1.0 };
            let mut flow = planar(n, |x| {
                let u = 0.1 * (w * (x[0] + x[1])).sin();
                [1.0, u, sign * u, 2.0]
            });
            let initial = flow.totals();
            let (kinetic0, internal0) = energy_partition(&flow);
            let mut op =
                MolecularTransport2d::new(&flow, law(&flow, mu, bulk, 0.0), || false).unwrap();
            op.advance(
                &mut flow,
                Time::new(end),
                MolecularControl::default(),
                || false,
            )
            .unwrap();
            // For a longitudinal plane wave div(tau)=(4mu/3+bulk) Laplacian(u).
            // For the perpendicular wave div(u)=0 and only mu remains.
            let diffusivity = if longitudinal {
                4.0 / 3.0 * mu + bulk
            } else {
                mu
            };
            let amplitude = 0.1 * (-2.0 * w * w * diffusivity * end).exp();
            let error: f64 = flow
                .states()
                .iter()
                .enumerate()
                .map(|(i, q)| {
                    let x = flow.grid().cell_center(i).unwrap();
                    let exact = amplitude * (w * (x[0] + x[1])).sin();
                    close(q.momentum[1], sign * q.momentum[0], 1e-12);
                    (q.momentum[0] - exact).abs()
                })
                .sum::<f64>()
                / (n * n) as f64;
            errors.push(error);
            let (kinetic1, internal1) = energy_partition(&flow);
            assert!(kinetic1 < kinetic0);
            assert!(internal1 > internal0);
            close(internal1 - internal0, kinetic0 - kinetic1, 2e-13);
            close(flow.totals().energy, initial.energy, 2e-13);
            close(flow.totals().momentum[0], initial.momentum[0], 2e-13);
        }
        assert!(errors[0] < 2e-3, "stress errors: {errors:?}");
        assert!(errors[1] < 0.3 * errors[0], "stress errors: {errors:?}");
    }
}

#[test]
fn g3_collision_data_changes_real_pressure_and_shear_without_changing_the_eos() {
    let init = |x: [f64; 2]| {
        let s = (2.0 * std::f64::consts::PI * x[1]).sin();
        [1.0, 0.2 * s, 0.0, 2.0 + 0.2 * s]
    };
    let mut original = planar(12, init);
    let mut enlarged = planar(12, init);
    let state = original.states().to_vec();
    let mut op =
        MolecularTransport2d::new(&original, kinetic(&original, 3.5e-10), || false).unwrap();
    let mut other =
        MolecularTransport2d::new(&enlarged, kinetic(&enlarged, 7.0e-10), || false).unwrap();
    op.advance(
        &mut original,
        Time::new(0.1),
        MolecularControl::default(),
        || false,
    )
    .unwrap();
    other
        .advance(
            &mut enlarged,
            Time::new(0.1),
            MolecularControl::default(),
            || false,
        )
        .unwrap();
    assert_eq!(original.gas(), enlarged.gas());
    let k = 3 * 12;
    let du = (original.states()[k].momentum[0] - state[k].momentum[0]).abs();
    let du_other = (enlarged.states()[k].momentum[0] - state[k].momentum[0]).abs();
    let p0 = state[k]
        .primitive(original.gas(), || false)
        .unwrap()
        .pressure;
    let dp = (original.states()[k]
        .primitive(original.gas(), || false)
        .unwrap()
        .pressure
        - p0)
        .abs();
    let dp_other = (enlarged.states()[k]
        .primitive(enlarged.gas(), || false)
        .unwrap()
        .pressure
        - p0)
        .abs();
    assert!(du > 1e-8 && dp > 1e-8);
    assert!((0.24..0.26).contains(&(du_other / du)));
    assert!((0.24..0.26).contains(&(dp_other / dp)));
    let bad = field(gas(0.4), *original.grid(), PERIODIC, init);
    assert!(matches!(
        MolecularTransport2d::new(&bad, law(&original, 0.01, 0.0, 0.01), || false),
        Err(MolecularError::Invalid(_))
    ));
}

fn coupled() -> MolecularEuler2d<ConstantTransport> {
    let flow = field(
        gas(0.4),
        CartesianGrid::new(4, 3, 0.25, 0.2, 0.5).unwrap(),
        Boundaries2d { x: WALLS, y: WALLS },
        |x| [1.0 + 0.1 * x[0], 0.2, -0.1, 1.0],
    );
    let transport = law(&flow, 0.02, 0.01, 0.05);
    MolecularEuler2d::new(flow, transport, || false).unwrap()
}

#[test]
fn g2_g4_coupled_wall_impulses_energy_and_late_cancellation_preserve_one_state() {
    let mut probe = coupled();
    let before = probe.flow().totals();
    let mut calls = 0;
    let report = probe
        .step(Time::new(0.001), MolecularStepControl::default(), || {
            calls += 1;
            false
        })
        .unwrap();
    assert_eq!(report.elapsed.value(), probe.flow().time().value());
    assert_eq!(report.first.elapsed.value(), 0.5 * report.elapsed.value());
    assert_eq!(report.second.elapsed.value(), report.first.elapsed.value());
    assert_eq!(report.advection.elapsed, report.elapsed);
    assert_eq!(report.boundary_change.mass, 0.0);
    assert_eq!(report.boundary_change.energy, 0.0);
    close(probe.flow().totals().energy, before.energy, 2e-13);
    let mut impulse = Totals::default();
    for side in SIDES {
        for q in probe.boundary_impulses(side) {
            impulse.momentum[0] += q.momentum[0];
            impulse.momentum[1] += q.momentum[1];
        }
    }
    for axis in 0..2 {
        close(
            probe.flow().totals().momentum[axis] - before.momentum[axis] + impulse.momentum[axis],
            0.0,
            2e-13,
        );
    }
    for cutoff in [1, calls / 2, calls] {
        let mut flow = coupled();
        let old = flow.flow().states().to_vec();
        let mut polls = 0;
        assert_eq!(
            flow.step(Time::new(0.001), MolecularStepControl::default(), || {
                polls += 1;
                polls == cutoff
            }),
            Err(MolecularError::Cancelled)
        );
        assert_eq!(flow.flow().states(), old);
        assert_eq!(flow.flow().time().value(), 0.0);
        assert_eq!(flow.last_step(), None);
        for side in SIDES {
            assert!(
                flow.boundary_impulses(side)
                    .iter()
                    .all(|q| *q == Totals::default())
            );
        }
        let next = flow
            .step(Time::new(0.001), MolecularStepControl::default(), || false)
            .unwrap();
        assert_eq!(next, report);
        assert_eq!(flow.flow().states(), probe.flow().states());
    }
}

#[test]
fn g0_zero_transport_recovers_original_euler_and_work_refusal_is_atomic() {
    let init = |x: [f64; 2]| [1.0, 0.1 * x[0], -0.1 * x[1], 2.0];
    let mut reference = planar(4, init);
    let flow = planar(4, init);
    let coefficients = law(&flow, 0.0, 0.0, 0.0);
    let mut coupled = MolecularEuler2d::new(flow, coefficients, || false).unwrap();
    let expected = reference
        .step(Time::new(0.001), StepControl::default(), || false)
        .unwrap();
    let report = coupled
        .step(Time::new(0.001), MolecularStepControl::default(), || false)
        .unwrap();
    assert_eq!(report.advection, expected);
    assert_eq!(coupled.flow().states(), reference.states());

    let mut flow = planar(4, init);
    let old = flow.states().to_vec();
    let mut operator =
        MolecularTransport2d::new(&flow, law(&flow, 0.1, 0.0, 0.1), || false).unwrap();
    let controls = MolecularControl {
        maximum_steps: 1,
        maximum_step_s: 1e-4,
        ..MolecularControl::default()
    };
    assert!(matches!(
        operator.advance(&mut flow, Time::new(0.01), controls, || false),
        Err(MolecularError::WorkLimit { .. })
    ));
    assert_eq!(flow.states(), old);
    assert_eq!(flow.time().value(), 0.0);
    for side in SIDES {
        assert!(
            operator
                .boundary_impulses(side)
                .iter()
                .all(|q| *q == Totals::default())
        );
    }
}
