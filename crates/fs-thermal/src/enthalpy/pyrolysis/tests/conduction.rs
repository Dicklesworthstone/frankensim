use super::*;
use crate::enthalpy::conduction::{ConductiveLink, ConductivityLaw};
use crate::enthalpy::pyrolysis::conduction::PyrolysisConductivity;
use crate::enthalpy::pyrolysis::implicit::ImplicitSolve;
use crate::enthalpy::pyrolysis::slab::ConductivityMixing;

struct Constant(f64);
impl ConductivityLaw for Constant {
    fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
        Ok(self.0)
    }
}
static VIRGIN_K: Constant = Constant(10.0);
static CHAR_K: Constant = Constant(0.1);
const FACE: ConductiveLink = ConductiveLink {
    link_id: 1,
    half_distances_m: [0.005; 2],
    contact_resistance_m2_k_w: 0.0,
};

fn conductivity<'a>(
    virgin: &'a dyn ConductivityLaw,
    char: Option<&'a dyn ConductivityLaw>,
    mixing: ConductivityMixing,
) -> PyrolysisConductivity<'a> {
    PyrolysisConductivity::Reactive {
        virgin,
        char,
        densities_kg_m3: [1000.0, if char.is_some() { 200.0 } else { 0.0 }],
        mixing,
        source: "manufactured constituent conductivities and declared component-volume arrangement",
    }
}
fn run<'a>(
    old: &PyrolysisSimulation<'a>,
    implicit: bool,
    target: f64,
    dt: f64,
) -> (PyrolysisSimulation<'a>, PyrolysisAudit) {
    let controls = PyrolysisControls {
        maximum_dt_s: dt,
        ..Default::default()
    };
    if implicit {
        let (state, report) = old
            .advance_implicit_to(
                target,
                controls,
                ImplicitSolve::default(),
                &mut old.implicit_workspace().unwrap(),
                || Ok(()),
            )
            .unwrap();
        (state, report.pyrolysis)
    } else {
        old.advance_to(target, controls, || Ok(())).unwrap()
    }
}
fn account(old: &PyrolysisSimulation<'_>, next: &PyrolysisSimulation<'_>, audit: &PyrolysisAudit) {
    close(
        next.mass_kg() + audit.volatile_mass_kg,
        old.mass_kg(),
        1e-15,
    );
    close(
        next.energy_j() + audit.volatile_enthalpy_j,
        old.energy_j(),
        2e-7,
    );
    assert!(audit.absolute_energy_residual_j < 2e-7);
    assert_eq!(audit.thermal.external_heat_j, 0.0);
    for (&id, cell) in next.cells() {
        assert_eq!(cell.geometry, old.cells()[&id].geometry);
    }
}

#[test]
fn original_constant_sources_preserve_transport_and_zero_conductivity_insulates() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    let old = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[800.0, 400.0], 0.0, 5.0),
        &[model(&VIRGIN, None, &gas, 0.0, 0.0)],
        || Ok(()),
    )
    .unwrap();
    let bound = old
        .clone()
        .with_conductivity(
            &[(
                MaterialId(1),
                conductivity(&VIRGIN_K, None, ConductivityMixing::Series),
            )],
            &[FACE],
            || Ok(()),
        )
        .unwrap();
    for implicit in [false, true] {
        let (reference, _) = run(&old, implicit, 0.2, 0.02);
        let (next, audit) = run(&bound, implicit, 0.2, 0.02);
        for id in [CellId(1), CellId(2)] {
            close(
                next.temperature_k(id).unwrap(),
                reference.temperature_k(id).unwrap(),
                1e-10,
            );
        }
        account(&old, &next, &audit);
    }
    let zero = Constant(0.0);
    let insulated = old
        .clone()
        .with_conductivity(
            &[(
                MaterialId(1),
                conductivity(&zero, None, ConductivityMixing::Series),
            )],
            &[FACE],
            || Ok(()),
        )
        .unwrap();
    for implicit in [false, true] {
        let (next, audit) = run(&insulated, implicit, 0.2, 0.2);
        assert_eq!(next.cells(), old.cells());
        assert_eq!(audit.volatile_mass_kg, 0.0);
        assert_eq!(next.conductivity_w_m_k(CellId(1)), Some(0.0));
    }
    let mut overridden = bound;
    assert!(
        overridden
            .set_conductances(&BTreeMap::from([(1, 2.0)]), || Ok(()))
            .is_err()
    );
    assert!(
        old.clone()
            .with_conductivity(
                &[(
                    MaterialId(1),
                    conductivity(&VIRGIN_K, None, ConductivityMixing::Series)
                )],
                &[FACE, FACE],
                || Ok(())
            )
            .is_err()
    );
    assert!(old.with_conductivity(&[], &[FACE], || Ok(())).is_err());
}

fn mixed_k(virgin: f64, char: f64, mixing: ConductivityMixing) -> f64 {
    let (v, c) = (virgin / 1000.0, char / 200.0);
    match mixing {
        ConductivityMixing::Parallel => (10.0 * v + 0.1 * c) / (v + c),
        ConductivityMixing::Series => (v + c) / (v / 10.0 + c / 0.1),
    }
}

#[test]
fn char_volume_changes_endpoint_heat_flux_and_refines_to_independent_trajectory() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    for mixing in [ConductivityMixing::Series, ConductivityMixing::Parallel] {
        let old = PyrolysisSimulation::new(
            thermal(&VIRGIN, &[800.0, 400.0], 0.0, 5.0),
            &[model(&VIRGIN, Some(&VIRGIN), &gas, 0.2, 1.0)],
            || Ok(()),
        )
        .unwrap()
        .with_conductivity(
            &[(
                MaterialId(1),
                conductivity(&VIRGIN_K, Some(&CHAR_K), mixing),
            )],
            &[FACE],
            || Ok(()),
        )
        .unwrap();
        let (one, audit) = run(&old, true, 0.1, 0.1);
        let virgin = 0.01 / 1.1;
        let char = 0.2 * (0.01 - virgin);
        let k = mixed_k(virgin, char, mixing);
        // The native reaction step heats all original mass before its actual
        // volatile release at T_end. New component volumes set endpoint G.
        let delta = 400.0 / (1.0 + 2.0 * 0.1 * (0.5 * k) / (0.01 * 2000.0));
        close(
            one.temperature_k(CellId(1)).unwrap(),
            600.0 + 0.5 * delta,
            1e-8,
        );
        close(one.conductivity_w_m_k(CellId(1)).unwrap(), k, 1e-12);
        account(&old, &one, &audit);

        // Independent continuous solution: equal histories retain their mean
        // T, while d(DeltaT)/dt = -2G(t) DeltaT/(m(t) Cp). Simpson integrates
        // the independently known exponential virgin/char inventories.
        let rate = |time: f64| {
            let v = 0.01 * (-time).exp();
            let c = 0.2 * (0.01 - v);
            mixed_k(v, c, mixing) / ((v + c) * 2000.0)
        };
        let count = 4096;
        let dt = 1.0 / count as f64;
        let mut integral = rate(0.0) + rate(1.0);
        for i in 1..count {
            integral += if i % 2 == 0 { 2.0 } else { 4.0 } * rate(i as f64 * dt);
        }
        let exact_delta = 400.0 * (-integral * dt / 3.0).exp();
        for implicit in [false, true] {
            let mut errors = [0.0; 2];
            for (rung, steps) in [100, 200].into_iter().enumerate() {
                let (next, audit) = run(&old, implicit, 1.0, 1.0 / steps as f64);
                account(&old, &next, &audit);
                let delta =
                    next.temperature_k(CellId(1)).unwrap() - next.temperature_k(CellId(2)).unwrap();
                errors[rung] = (delta - exact_delta).abs();
            }
            assert!(
                errors[1] < 0.6 * errors[0],
                "{mixing:?} implicit={implicit}: {errors:?}"
            );
            assert!(
                errors[1] < 1.0,
                "{mixing:?} implicit={implicit}: {errors:?}"
            );
        }
    }
}

struct Changing;
impl ConductivityLaw for Changing {
    fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
        let x = 1.0 + h / 500_000.0;
        Ok(0.1 + x * x)
    }
}

#[test]
fn activated_char_queries_each_original_constituent_chart_at_its_actual_temperature() {
    struct OriginalK<'a> {
        caloric: &'a dyn CondensedCaloric,
        intercept: f64,
        slope: f64,
    }
    impl ConductivityLaw for OriginalK<'_> {
        fn conductivity_w_m_k(&self, h: f64) -> Result<f64, DepositionError> {
            Ok(self.intercept + self.slope * (self.caloric.temperature_k(h)? - 300.0))
        }
    }
    let virgin = OriginalK {
        caloric: &VIRGIN,
        intercept: 10.0,
        slope: 0.001,
    };
    let char_law = OriginalK {
        caloric: &CHAR,
        intercept: 0.1,
        slope: 0.0001,
    };
    let gas = |t| Ok(2000.0 * (t - 300.0) + 1_000_000.0);
    let mut material = model(&VIRGIN, Some(&CHAR), &gas, 0.2, 2.0);
    material.kinetics.activation_energy_j_mol = 4000.0;
    let old = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[800.0, 400.0], 0.0, 5.0),
        &[material],
        || Ok(()),
    )
    .unwrap()
    .with_conductivity(
        &[(
            MaterialId(1),
            conductivity(&virgin, Some(&char_law), ConductivityMixing::Parallel),
        )],
        &[FACE],
        || Ok(()),
    )
    .unwrap();
    let (next, audit) = run(&old, true, 0.1, 0.1);
    for id in [CellId(1), CellId(2)] {
        let inventory = next.inventories()[&id];
        assert!(inventory.char_mass_kg > 0.0);
        let v = inventory.virgin_mass_kg / 1000.0;
        let c = inventory.char_mass_kg / 200.0;
        let t = inventory.temperature_k;
        let expected =
            (v * (10.0 + 0.001 * (t - 300.0)) + c * (0.1 + 0.0001 * (t - 300.0))) / (v + c);
        close(next.conductivity_w_m_k(id).unwrap(), expected, 1e-12);
    }
    account(&old, &next, &audit);
}

#[test]
fn temperature_dependent_original_conductivity_closes_the_nonlinear_endpoint() {
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    let source = Changing;
    let old = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[900.0, 400.0], 0.0, 5.0),
        &[model(&VIRGIN, None, &gas, 0.0, 0.0)],
        || Ok(()),
    )
    .unwrap()
    .with_conductivity(
        &[(
            MaterialId(1),
            conductivity(&source, None, ConductivityMixing::Series),
        )],
        &[FACE],
        || Ok(()),
    )
    .unwrap();
    let dt = 2.0;
    let (mut lo, mut hi) = (0.0, 5000.0);
    for _ in 0..100 {
        let q = 0.5 * (lo + hi);
        let (ta, tb) = (900.0 - q / 20.0, 400.0 + q / 20.0);
        let k = |t: f64| 0.1 + (1.0 + 2000.0 * (t - 300.0) / 500_000.0).powi(2);
        let g = 1.0 / (1.0 / k(ta) + 1.0 / k(tb));
        if q > dt * g * (ta - tb) {
            hi = q;
        } else {
            lo = q;
        }
    }
    let (next, audit) = run(&old, true, dt, dt);
    close(
        next.temperature_k(CellId(1)).unwrap(),
        900.0 - 0.5 * (lo + hi) / 20.0,
        1e-8,
    );
    close(
        next.temperature_k(CellId(2)).unwrap(),
        400.0 + 0.5 * (lo + hi) / 20.0,
        1e-8,
    );
    account(&old, &next, &audit);
}

#[test]
fn absent_sources_are_not_queried_and_domain_failure_cannot_publish_reaction() {
    struct Unavailable(std::cell::Cell<usize>);
    impl ConductivityLaw for Unavailable {
        fn conductivity_w_m_k(&self, _: f64) -> Result<f64, DepositionError> {
            self.0.set(self.0.get() + 1);
            Err(DepositionError::MaterialDomain)
        }
    }
    let absent = Unavailable(std::cell::Cell::new(0));
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    for rate in [0.0, 1.0] {
        let old = PyrolysisSimulation::new(
            thermal(&VIRGIN, &[800.0, 400.0], 0.0, 5.0),
            &[model(&VIRGIN, Some(&VIRGIN), &gas, 0.2, rate)],
            || Ok(()),
        )
        .unwrap()
        .with_conductivity(
            &[(
                MaterialId(1),
                conductivity(&VIRGIN_K, Some(&absent), ConductivityMixing::Series),
            )],
            &[FACE],
            || Ok(()),
        )
        .unwrap();
        let cells = old.cells().clone();
        let inventories = old.inventories().clone();
        let controls = PyrolysisControls {
            maximum_dt_s: 0.1,
            ..Default::default()
        };
        if rate == 0.0 {
            run(&old, false, 0.1, 0.1);
            run(&old, true, 0.1, 0.1);
            assert_eq!(absent.0.get(), 0);
        } else {
            assert_eq!(
                old.advance_to(0.1, controls, || Ok(())).err(),
                Some(DepositionError::MaterialDomain)
            );
            assert!(
                old.advance_implicit_to(
                    0.1,
                    controls,
                    ImplicitSolve::default(),
                    &mut old.implicit_workspace().unwrap(),
                    || Ok(())
                )
                .is_err()
            );
            assert!(absent.0.get() > 0);
            assert_eq!(old.cells(), &cells);
            assert_eq!(old.inventories(), &inventories);
            assert_eq!(old.time_s(), 0.0);
        }
    }
}

#[test]
fn present_zero_k_char_insulates_and_late_cancellation_replays_workspace() {
    let zero = Constant(0.0);
    let gas = |t| VIRGIN.specific_enthalpy_j_kg(t);
    let old = PyrolysisSimulation::new(
        thermal(&VIRGIN, &[800.0, 400.0], 0.0, 5.0),
        &[model(&VIRGIN, Some(&VIRGIN), &gas, 0.2, 1.0)],
        || Ok(()),
    )
    .unwrap()
    .with_conductivity(
        &[(
            MaterialId(1),
            conductivity(&VIRGIN_K, Some(&zero), ConductivityMixing::Series),
        )],
        &[FACE],
        || Ok(()),
    )
    .unwrap();
    let controls = PyrolysisControls {
        maximum_dt_s: 0.1,
        ..Default::default()
    };
    let mut workspace = old.implicit_workspace().unwrap();
    let mut complete = 0;
    let (expected, receipt) = old
        .advance_implicit_to(
            0.1,
            controls,
            ImplicitSolve::default(),
            &mut workspace,
            || {
                complete += 1;
                Ok(())
            },
        )
        .unwrap();
    for id in [CellId(1), CellId(2)] {
        assert!(expected.inventories()[&id].char_mass_kg > 0.0);
        assert_eq!(expected.conductivity_w_m_k(id), Some(0.0));
        close(
            expected.temperature_k(id).unwrap(),
            old.temperature_k(id).unwrap(),
            1e-9,
        );
    }
    account(&old, &expected, &receipt.pyrolysis);
    let mut calls = 0;
    assert_eq!(
        old.advance_implicit_to(
            0.1,
            controls,
            ImplicitSolve::default(),
            &mut workspace,
            || {
                calls += 1;
                if calls == complete {
                    Err(DepositionError::Cancelled)
                } else {
                    Ok(())
                }
            }
        )
        .err(),
        Some(crate::enthalpy::pyrolysis::implicit::ImplicitError::Cancelled)
    );
    assert_eq!(old.time_s(), 0.0);
    let (next, report) = old
        .advance_implicit_to(
            0.1,
            controls,
            ImplicitSolve::default(),
            &mut workspace,
            || Ok(()),
        )
        .unwrap();
    assert_eq!(next.cells(), expected.cells());
    assert_eq!(report, receipt);
}
