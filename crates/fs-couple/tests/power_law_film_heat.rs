//! G1/G3/G4: original nonlinear film work, finite enthalpy and one accepted clock.
//! Manufactured SI rheology and caloric charts test coupling, not fluid data.

#[path = "power_law_film_heat/feedback.rs"]
mod feedback;

use fs_blake3::ContentHash;
use fs_couple::render::plate::impact::material::heating::transport::CoupledThermalStep;
use fs_couple::render::plate::impact::material::heating::{
    CaloricReservoir, HeatedImpactSystem, PowerLawFilmHeatRoute, SqueezeFilmHeatRoute,
};
use fs_couple::render::plate::impact::squeeze::power_law::{
    PowerLawControls, PowerLawFilm, RegularizedPowerLaw,
};
use fs_couple::render::plate::impact::squeeze::{
    FilmCell, FilmChannel, FilmLimits, GapPort, ResistiveFilm,
};
use fs_couple::render::plate::impact::{ImpactBody, ImpactConfig, ImpactError, ImpactSystem};
use fs_exec::CancelGate;
use fs_material::phase::{
    EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve, SolidLiquidPhase, UniformEnthalpyStep,
    UniformEnthalpyStepInput,
};

const M1: f64 = 0.2;
const M2: f64 = 0.4;
const GAP: f64 = 0.001;
const WIDTH: f64 = 0.02;
const LENGTH: f64 = 0.03;
const AREA: f64 = WIDTH * LENGTH;
const DT: f64 = 0.0001;
const LAUNCH: f64 = 0.01;
const MASS: f64 = 0.0006;
const CP: f64 = 0.2;

fn law(n: f64) -> RegularizedPowerLaw {
    RegularizedPowerLaw {
        // Comparable admitted forces in three distinct manufactured laws.
        consistency_pa_s_n: if n < 1.0 { 0.05 } else if n == 1.0 { 0.002 } else { 1e-6 },
        flow_index: n,
        minimum_shear_rate_per_s: 1.0,
    }
}
fn graph(moving: bool, viscosity: f64) -> ResistiveFilm {
    let gap = GapPort {
        reference_m: GAP,
        closure: if moving { vec![1.0 / M1.sqrt(), -1.0 / M2.sqrt()] } else { vec![0.0; 2] },
    };
    ResistiveFilm::new(
        vec![FilmCell { area_m2: AREA, gap: gap.clone() }],
        [1.0, 2.0].into_iter().map(|width| FilmChannel {
            from: 0, to: None, width_m: width * WIDTH, length_m: LENGTH, gap: gap.clone(),
        }).collect(),
        2, viscosity,
        FilmLimits { minimum_cell_gap_m: 1e-6, maximum_gap_m: 0.01, maximum_pressure_pa: 1e5 },
    ).unwrap()
}
fn mechanics(n: f64, normal: bool) -> ImpactSystem {
    let (a, _) = ImpactBody::free_mass(M1, 0.0, LAUNCH).unwrap();
    let (b, _) = ImpactBody::free_mass(M2, 0.0, 0.0).unwrap();
    let mut system = ImpactSystem::new(vec![a, b], vec![], vec![], vec![], ImpactConfig {
        dt_s: DT, max_steps: 128, maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8,
        maximum_generalized_force: 1.0,
    }).unwrap();
    let law = law(n);
    let film = PowerLawFilm::new(&graph(true, law.consistency_pa_s_n), vec![law; 2],
        PowerLawControls {
            maximum_iterations: 64, maximum_line_search_steps: 32,
            absolute_flow_tolerance_m3_s: 1e-15, relative_flow_tolerance: 1e-10,
            absolute_power_tolerance_w: 1e-13, relative_power_tolerance: 1e-10,
        }).unwrap();
    system = system.with_power_law_film(film).unwrap();
    if normal { system = system.with_squeeze_film(graph(false, 0.002)).unwrap(); }
    system
}
fn curve(upper_h: f64) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_single_phase(ContentHash([92; 32]), SolidLiquidPhase::Liquid,
        vec![
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: 0.0, temperature_k: 300.0,
                liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0 },
            EnthalpyPhaseKnot { specific_enthalpy_j_kg: upper_h, temperature_k: 300.0 + upper_h / CP,
                liquid_mass_fraction: 1.0, bulk_density_kg_m3: 1000.0 },
        ]).unwrap()
}
fn stores(c: &EquilibriumEnthalpyPhaseCurve) -> Vec<CaloricReservoir<'_>> {
    ["a", "b"].into_iter().map(|key| CaloricReservoir::new(key.into(), c,
        c.state_at_specific_enthalpy(0.0).unwrap(), MASS, 6e-7, 0.0012).unwrap()).collect()
}
fn route(channel: Option<usize>, fraction: f64) -> PowerLawFilmHeatRoute {
    PowerLawFilmHeatRoute { film: 0, channel, to_a: "a".into(), to_b: "b".into(), fraction_to_a: fraction }
}
fn heated<'a>(c: &'a EquilibriumEnthalpyPhaseCurve, n: f64,
    routes: Vec<PowerLawFilmHeatRoute>, normal: bool) -> Result<HeatedImpactSystem<'a>, ImpactError> {
    let normal_routes = if normal { vec![SqueezeFilmHeatRoute {
        film: 0, to_a: "a".into(), to_b: "b".into(), fraction_to_a: 0.5,
    }] } else { vec![] };
    HeatedImpactSystem::new_with_power_law_losses(mechanics(n, normal).prepare_analytic().unwrap(),
        stores(c), vec![], normal_routes, vec![], vec![], vec![], routes, vec![], 1e-12)
}
fn relative(x: &[f64]) -> f64 { x[1] / M1.sqrt() - x[3] / M2.sqrt() }
fn closure(x: &[f64]) -> f64 { x[0] / M1.sqrt() - x[2] / M2.sqrt() }
fn momentum(x: &[f64]) -> f64 { M1.sqrt() * x[1] + M2.sqrt() * x[3] }
fn near(a: f64, b: f64, tolerance: f64) {
    assert!((a - b).abs() <= tolerance, "{a:e} != {b:e}; tolerance {tolerance:e}");
}

// Independent resolved shear-stress integration. Split at the physical core,
// then quadrature integrates y*gamma for flow and tau*gamma for true heat.
// No native channel formula, force work or energy decrement enters this oracle.
fn profile_work(pressure: f64, gap: f64, law: RegularizedPowerLaw) -> (f64, f64) {
    let half = 0.5 * gap;
    let gradient = pressure.abs() / LENGTH;
    if gradient == 0.0 { return (0.0, 0.0); }
    let k = law.consistency_pa_s_n;
    let core = (k / gradient).min(half);
    let integrand = |y: f64, heat: bool| {
        let stress = gradient * y;
        let gamma = if stress <= k { stress / k } else { (stress / k).powf(1.0 / law.flow_index) };
        (if heat { stress } else { y }) * gamma
    };
    let integrate = |lo: f64, hi: f64, heat: bool| {
        let step = (hi - lo) / 512.0;
        let mut sum = integrand(lo, heat) + integrand(hi, heat);
        for i in 1..512 {
            sum += (if i % 2 == 0 { 2.0 } else { 4.0 }) * integrand(lo + i as f64 * step, heat);
        }
        step * sum / 3.0
    };
    let total_width = 3.0 * WIDTH;
    (2.0 * total_width * (integrate(0.0, core, false) + integrate(core, half, false)),
     2.0 * total_width * LENGTH * (integrate(0.0, core, true) + integrate(core, half, true)))
}

#[test]
fn thinning_newtonian_and_thickening_profiles_heat_finite_stores_on_the_original_step() {
    let c = curve(20.0);
    for n in [0.5, 1.0, 2.0] {
        for analytic in [false, true] {
            let mut system = heated(&c, n, vec![route(None, 0.25)], false).unwrap();
            system.set_analytic_newton(analytic);
            let mut native = mechanics(n, false).prepare_analytic().unwrap();
            native.set_analytic_newton(analytic);
            let initial = native.stored_energy_j();
            let initial_momentum = momentum(native.state());
            let mut heat = 0.0;
            for tick in 1..=16 {
                let before = system.mechanics().state().to_vec();
                let frame = system.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
                let plain = native.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
                let after = system.mechanics().state();
                assert_eq!(after, native.state());
                assert_eq!(frame.mechanical.dissipated_energy_j, plain.dissipated_energy_j);
                // Pressure is reconstructed from free-body momentum, independently
                // of the film's own pressure or dissipation report.
                let pressure = (relative(&before) - relative(after)) / (DT * (1.0 / M1 + 1.0 / M2) * AREA);
                let gap = GAP - 0.5 * (closure(&before) + closure(after));
                assert!(pressure * gap / (2.0 * LENGTH) > 4.0 * law(n).consistency_pa_s_n);
                let (flow, power) = profile_work(pressure, gap, law(n));
                near(flow, AREA * 0.5 * (relative(&before) + relative(after)), 1e-12);
                near(frame.deposited_heat_j, DT * power, 2e-15 + 1e-7 * DT * power);
                near(frame.remaining_mechanical_loss_j, 0.0, 2e-15);
                near(frame.thermal_energy_change_j, frame.deposited_heat_j, 1e-12);
                near(frame.combined_energy_residual_j, 0.0, 2e-12);
                near(momentum(after), initial_momentum, 1e-12);
                assert_eq!(frame.mechanical.sample, tick);
                near(frame.mechanical.time_s, tick as f64 * DT, 1e-15);
                heat += frame.deposited_heat_j;
            }
            near(system.reservoirs()[0].received_heat_j(), 0.25 * heat, 2e-18);
            near(system.reservoirs()[1].received_heat_j(), 0.75 * heat, 2e-18);
            near(system.reservoirs()[0].state().temperature_k() - 300.0, 0.25 * heat / (MASS * CP), 2e-11);
            near(system.mechanics().stored_energy_j() + system.reservoirs().iter()
                .map(CaloricReservoir::storage_change_j).sum::<f64>(), initial, 2e-11);
            assert_eq!(system.mechanics().power_law_film(0).unwrap().laws(), &[law(n); 2]);
        }
    }
}

#[test]
fn selected_passages_keep_unrouted_loss_and_independent_film_indices_in_canonical_order() {
    let c = curve(20.0);
    let mut partial = heated(&c, 0.5, vec![route(Some(0), 1.0)], true).unwrap();
    let mut full = heated(&c, 0.5, vec![route(Some(0), 1.0), route(Some(1), 0.0)], true).unwrap();
    let mut reversed = heated(&c, 0.5, vec![route(Some(1), 0.0), route(Some(0), 1.0)], true).unwrap();
    let initial = partial.mechanics().stored_energy_j();
    let mut outside = 0.0;
    for _ in 0..16 {
        let one = partial.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        let all = full.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        let reverse = reversed.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
        assert_eq!(partial.mechanics().state(), full.mechanics().state());
        assert_eq!(full.mechanics().state(), reversed.mechanics().state());
        assert_eq!(all.deposited_heat_j, reverse.deposited_heat_j);
        assert_eq!(full.reservoirs(), reversed.reservoirs());
        // Identical laws and gap: the widths alone set the independent 1:2 split.
        near(one.deposited_heat_j, all.deposited_heat_j / 3.0, 2e-15);
        near(one.remaining_mechanical_loss_j, 2.0 * one.deposited_heat_j, 2e-15);
        assert_eq!(partial.reservoirs()[1].received_heat_j(), 0.0);
        outside += one.remaining_mechanical_loss_j;
    }
    assert!(outside > 0.0);
    assert_eq!(partial.mechanics().squeeze_film_count(), 1);
    assert_eq!(partial.mechanics().power_law_film_count(), 1);
    assert!(partial.layer_heat_profile(0, 0).is_none());
    near(partial.mechanics().stored_energy_j() + partial.reservoirs().iter()
        .map(CaloricReservoir::storage_change_j).sum::<f64>() + outside, initial, 2e-11);
}

fn adiabatic(input: UniformEnthalpyStepInput<'_>) -> UniformEnthalpyStep<()> {
    UniformEnthalpyStep {
        state: input.curve.state_at_specific_enthalpy(input.initial.specific_enthalpy_j_kg()
            + input.internal_heat_j / input.mass_kg).unwrap(),
        external_heat_j: 0.0, energy_residual_tolerance_j: 0.0, report: (),
    }
}

#[test]
fn complete_thermal_trial_sees_original_heat_and_late_cancel_preserves_every_owner() {
    let c = curve(20.0);
    let mut system = heated(&c, 0.5, vec![route(None, 0.25)], false).unwrap();
    let mut untouched = heated(&c, 0.5, vec![route(None, 0.25)], false).unwrap();
    let motion = system.mechanics().state().to_vec();
    let thermal = system.reservoirs().to_vec();
    let late = CancelGate::new_clock_free();
    let mut calls = 0;
    let refused = system.step_with_coupled_transport(&[0.0; 2], &late, |trial| {
        calls += 1;
        assert_eq!(trial.sample(), 1);
        assert_eq!(trial.duration_s(), DT);
        let inputs = trial.reservoirs();
        assert_eq!(inputs.len(), 2);
        assert!(inputs[0].thermal.internal_heat_j > 0.0);
        near(inputs[1].thermal.internal_heat_j, 3.0 * inputs[0].thermal.internal_heat_j, 2e-18);
        for (input, old) in inputs.iter().zip(&thermal) { assert_eq!(input.thermal.initial, old.state()); }
        let transfers = inputs.iter().map(|input| adiabatic(input.thermal)).collect();
        late.request();
        Ok(CoupledThermalStep { transfers, boundary_heat_j: 0.0 })
    });
    assert!(matches!(refused, Err(ImpactError::Cancelled)));
    assert_eq!(calls, 1);
    assert_eq!(system.mechanics().state(), motion);
    assert_eq!(system.mechanics().samples(), 0);
    assert_eq!(system.reservoirs(), thermal);
    let retry = system.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
    let expected = untouched.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
    assert_eq!(retry.deposited_heat_j, expected.deposited_heat_j);
    assert_eq!(system.reservoirs(), untouched.reservoirs());
    assert_eq!(system.mechanics().state(), untouched.mechanics().state());
    // The same source reaches the original generic thermal carrier on the next
    // interval; switching carriers neither repeats the previous deposit nor clock.
    let carried = system.step_with_transport(&[0.0; 2], &CancelGate::new_clock_free(),
        |_, input| Ok(adiabatic(input))).unwrap();
    let plain = untouched.step(&[0.0; 2], &CancelGate::new_clock_free()).unwrap();
    assert_eq!(carried.impact.mechanical.sample, 2);
    assert_eq!(system.mechanics().state(), untouched.mechanics().state());
    assert_eq!(carried.impact.deposited_heat_j, plain.deposited_heat_j);
    near(carried.boundary_heat_j, 0.0, 0.0);
    for (a, b) in system.reservoirs().iter().zip(untouched.reservoirs()) {
        near(a.state().specific_enthalpy_j_kg(), b.state().specific_enthalpy_j_kg(), 1e-12);
        assert_eq!(a.received_heat_j(), b.received_heat_j());
    }
}

#[test]
fn route_overlap_wrong_owner_and_thermal_domain_refuse_without_partial_publication() {
    let c = curve(20.0);
    for routes in [
        vec![route(None, 0.5), route(Some(0), 0.5)],
        vec![route(Some(0), 0.5), route(None, 0.5)],
        vec![route(Some(0), 0.5), route(Some(0), 0.5)],
        vec![route(Some(2), 0.5)],
        vec![PowerLawFilmHeatRoute { film: 1, ..route(None, 0.5) }],
        vec![PowerLawFilmHeatRoute { to_a: "absent".into(), ..route(None, 0.5) }],
        vec![route(None, f64::NAN)],
    ] {
        assert!(heated(&c, 0.5, routes, false).is_err());
    }
    let n1 = mechanics(1.0, false);
    let core = n1.power_law_film(0).unwrap().pressure_operator().clone();
    let (bare, _) = ImpactBody::free_mass(M1, 0.0, LAUNCH).unwrap();
    let (other, _) = ImpactBody::free_mass(M2, 0.0, 0.0).unwrap();
    let normal_only = ImpactSystem::new(vec![bare, other], vec![], vec![], vec![], ImpactConfig {
        dt_s: DT, max_steps: 128, maximum_energy_j: 1.0,
        energy_absolute_tolerance_j: 1e-12, energy_relative_tolerance: 1e-8,
        maximum_generalized_force: 1.0,
    }).unwrap().with_squeeze_film(core).unwrap().prepare_analytic().unwrap();
    assert!(HeatedImpactSystem::new_with_power_law_losses(normal_only, stores(&c), vec![], vec![],
        vec![], vec![], vec![], vec![route(None, 0.5)], vec![], 1e-12).is_err());

    let narrow = curve(1e-8);
    let mut system = heated(&narrow, 0.5, vec![route(None, 0.25)], false).unwrap();
    let motion = system.mechanics().state().to_vec();
    let thermal = system.reservoirs().to_vec();
    assert!(system.step(&[0.0; 2], &CancelGate::new_clock_free()).is_err());
    assert_eq!(system.mechanics().state(), motion);
    assert_eq!(system.mechanics().samples(), 0);
    assert_eq!(system.reservoirs(), thermal);
}
