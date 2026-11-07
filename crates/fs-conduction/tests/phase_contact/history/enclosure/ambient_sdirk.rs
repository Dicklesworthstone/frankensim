//! Phase-aware contact and independently scheduled material black reservoirs.
use super::*;
use fs_conduction::bc::history::march::enthalpy::radiation::enclosure::MaterialEnclosureEnthalpyHistoryRun;
use fs_conduction::bc::history::march::enthalpy::radiation::material::MaterialRadiationRegionHistory;
use fs_conduction::transient::enthalpy::sdirk::{EnthalpySdirk2Stage, EnthalpySdirk2Step, SDIRK2_GAMMA as G};

fn ambient_patches<'a>(
    cards: &'a [MaterialCard; 2], curves: &'a [EquilibriumEnthalpyPhaseCurve; 2],
    context: &'a QueryPoint, schedules: &[TimeSignal; 2],
) -> Vec<MaterialRadiationRegionHistory<'a>> {
    [0, 1].map(|body| MaterialRadiationRegionHistory {
        region: REGIONS[body].into(), curve: &curves[body],
        material: material(&cards[body], context),
        selection: MaterialPropertySelection::SingleClaimOnly,
        ambient_temperature: schedules[body].clone(),
    }).to_vec()
}

// Independent face moments. Radiation acts on [0,1,3]/[4,5,7] while contact
// acts on [0,1,2]/[4,5,6], so one mean temperature/phase cannot replace the fields.
fn loads(h: &[f64; 8], baths: [f64; 2]) -> ([f64; 8], [f64; 2], f64) {
    let fa = h[..3].iter().sum::<f64>() / 300.0;
    let fb = h[4..7].iter().sum::<f64>() / 300.0;
    let contact = AREA * (400.0 - 500.0) / (1.0 + 2.0 * fa + fb);
    let mut outward = [0.0; 8];
    for local in 0..3 {
        outward[local] = contact / 3.0;
        outward[4 + local] = -contact / 3.0;
    }
    let mut radiation = [0.0; 2];
    for body in 0..2 {
        let temperature = if body == 0 { 400.0_f64 } else { 500.0_f64 };
        let black = STEFAN_BOLTZMANN_W_M2_K4 * (temperature.powi(4) - baths[body].powi(4));
        let epsilon = RADIATING[body].map(|node| 0.2 + 0.006 * h[node]);
        for (local, node) in RADIATING[body].into_iter().enumerate() {
            let q = AREA / 12.0 * (epsilon.iter().sum::<f64>() + epsilon[local]) * black;
            outward[node] += q;
            radiation[body] += q;
        }
    }
    (outward, radiation, contact)
}

fn stage(old: &[f64; 8], guess: [f64; 8], dt: f64, baths: [f64; 2], source: f64) -> [f64; 8] {
    let mut h = guess;
    for _ in 0..128 {
        let (outward, _, _) = loads(&h, baths);
        h = std::array::from_fn(|i| old[i] + dt * (source / 24.0 - outward[i])
            * if i < 4 { 24.0 } else { 12.0 });
    }
    assert!(h.iter().all(|h| *h > 0.0 && *h < 100.0));
    h
}

fn stages(old: &[f64; 8], dt: f64, baths: [f64; 2], source: f64) -> [[f64; 8]; 2] {
    let first = stage(old, *old, G * dt, baths, source);
    let (outward, _, _) = loads(&first, baths);
    // Pure tableau algebra; the carry is never queried as a material state.
    let carry = std::array::from_fn(|i| old[i] + (1.0 - G) * dt
        * (source / 24.0 - outward[i]) * if i < 4 { 24.0 } else { 12.0 });
    [first, stage(&carry, first, G * dt, baths, source)]
}

fn admit_stage(
    stage: &EnthalpySdirk2Stage, expected: &[f64; 8], baths: [f64; 2],
    cards: &[MaterialCard; 2],
) {
    for i in 0..8 {
        near(stage.specific_enthalpy_j_kg[i], expected[i], 4e-8);
        near(stage.liquid_mass_fraction[i], expected[i] / 100.0, 4e-10);
    }
    assert_eq!(stage.temperature, [400.0, 400.0, 400.0, 400.0, 500.0, 500.0, 500.0, 500.0]);
    let (_, heat, _) = loads(expected, baths);
    let radiation = stage.pointwise_radiation.as_ref().unwrap();
    assert!(stage.pointwise_enclosure_radiation.is_none());
    near(radiation.outward_heat_w, heat.iter().sum(), 3e-7);
    for (body, patch) in radiation.patches.iter().enumerate() {
        assert_eq!(patch.patch.ambient_temperature_k(), baths[body]);
        near(patch.outward_heat_w, heat[body], 3e-7);
        for state in &patch.material_states {
            near(state.emissivity.value(), 0.2 + 0.006 * expected[state.vertex], 4e-10);
            assert_eq!(state.emissivity.card_identity(), cards[body].content_hash());
        }
    }
    for trace in stage.physical_contacts.as_ref().unwrap().material_trace_states() {
        for body in 0..2 {
            near(trace.liquid_mass_fractions[body],
                expected[4 * body..4 * body + 3].iter().sum::<f64>() / 300.0, 4e-10);
        }
    }
    assert!(stage.stage_energy_residual_j.abs() <= config().step.energy_tolerance_j);
    assert!(stage.stage_residual_norm_j <= config().residual_tolerance_j);
}

fn admit(
    step: &EnthalpySdirk2Step, old: &[f64; 8], expected: &[[f64; 8]; 2],
    dt: f64, baths: [f64; 2], source: f64, cards: &[MaterialCard; 2],
) {
    admit_stage(&step.first_stage, &expected[0], baths, cards);
    admit_stage(&step.second_stage, &expected[1], baths, cards);
    let (_, q1, _) = loads(&expected[0], baths);
    let (_, q2, _) = loads(&expected[1], baths);
    let radiated = dt * ((1.0 - G) * q1.iter().sum::<f64>() + G * q2.iter().sum::<f64>());
    let stored: f64 = expected[1].iter().zip(old).enumerate()
        .map(|(i, (new, old))| (new - old) * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 }).sum();
    near(step.source_j, dt * source / 3.0, 2e-9);
    near(step.radiation_out_j, radiated, 2e-9);
    near(step.stored_energy_change_j, stored, 2e-9);
    near(stored + radiated - dt * source / 3.0, 0.0, 2e-9);
    assert!(step.energy_residual_j.abs() <= config().step.energy_tolerance_j);
}

#[test]
fn sdirk_contact_and_independent_bath_events_share_the_actual_material_phase_fields() {
    with_cx(|cx| {
        let mesh = tetrahedra();
        let boundary = history(cx, &mesh);
        let initial = boundary.at(cx, 0.0).unwrap();
        let cards = [material_card("cold"), material_card("hot")];
        let context = context();
        let curves = [chart(&cards[0], 400.0), chart(&cards[1], 500.0)];
        let source = source(&cards, 1.0, false, false, 1.0);
        let conductivity = ConductivityModel::isotropic_declared(2.0).unwrap();
        let schedules = [
            TimeSignal::Table { times: vec![0.0, 0.0004, 0.0008],
                values: vec![450.0, 350.0, 450.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold },
            TimeSignal::Table { times: vec![0.0, 0.0006, 0.0008],
                values: vec![350.0, 550.0, 450.0], dims: TEMPERATURE_DIMS, interp: Interp::Hold },
        ];
        let create = || {
            let forcing = HeatSourceHistory::Uniform(TimeSignal::Table {
                times: vec![0.0, 0.0006], values: vec![0.0, 60.0],
                dims: VOLUMETRIC_SOURCE_DIMS, interp: Interp::Hold });
            ContactEnthalpyHistoryRun::new_newton(cx,
                base(cx, &boundary, &conductivity, &curves, &OLD, forcing, 4),
                contacts(cx, &mesh, &initial, &cards, &curves, &source, &context,
                    sides(&mesh, 4, 2, 0.0), 0, false), config()).unwrap()
                .into_sdirk2().unwrap()
                .with_material_ambient(cx, ambient_patches(&cards, &curves, &context, &schedules)).unwrap()
        };
        let mut run = create();
        assert!(run.contacts().is_some());
        assert_eq!(run.material_ambient_patches().unwrap().len(), 2);
        assert!(run.material_enclosure_patches().is_none());
        let mut old = OLD;
        let mut from = 0.0;
        let mut heat = 0.0;
        for (to, baths, q) in [
            (0.0004, [450.0, 350.0], 0.0),
            (0.0006, [350.0, 350.0], 0.0),
            (0.0008, [350.0, 550.0], 60.0),
            (0.001, [450.0, 450.0], 60.0),
        ] {
            let dt = to - from;
            let expected = stages(&old, dt, baths, q);
            let step = run.advance_event_step(cx, 0.001, 1.0, config().step).unwrap();
            assert_eq!(step.end_time_s, to);
            admit(&step.report, &old, &expected, dt, baths, q, &cards);
            heat += step.report.source_j - step.report.radiation_out_j;
            old = expected[1];
            from = to;
        }
        assert_ne!(old[0], old[2]);
        assert_ne!(old[4], old[6]);
        let stored: f64 = old.iter().zip(OLD).enumerate()
            .map(|(i, (new, initial))| (new - initial) * if i < 4 { 1.0 / 24.0 } else { 1.0 / 12.0 }).sum();
        near(heat, stored, 3e-9);
        assert_eq!(run.conduction().accepted_steps(), 4);
        assert_eq!(run.conduction().time_s(), 0.001);

        // The same combined Newton stages supply an accepted entropy window.
        // Internal contact contributes production, never an external heat input.
        let mut audited = create();
        let report = audited.advance_adaptive_step_with_entropy(cx, 0.0004, 0.0004,
            config().step, AdaptiveHistoryConfig { absolute_temperature_k: 1e-10,
                relative_temperature: 0.0, min_step_s: 1e-8, max_step_s: 0.001, max_attempts: 12 },
            PhaseErrorTolerance { absolute_enthalpy_j_kg: 10.0,
                relative_enthalpy: 0.0, absolute_liquid_fraction: 0.1 }, 1e-8).unwrap();
        assert_eq!(report.physical.end_time_s, 0.0004);
        assert_eq!(report.physical.rejected_attempts, 0);
        let mut expected = OLD;
        let mut entropy_input = 0.0;
        let mut radiation_generation = 0.0;
        let mut contact_generation = 0.0;
        let baths = [450.0, 350.0];
        for physical in [&report.physical.first_half, &report.physical.second_half] {
            let pair = stages(&expected, 0.0002, baths, 0.0);
            admit(physical, &expected, &pair, 0.0002, baths, 0.0, &cards);
            for (stage, weight) in [(&pair[0], 1.0 - G), (&pair[1], G)] {
                let (_, radiation, contact) = loads(stage, baths);
                for body in 0..2 {
                    let temperature = if body == 0 { 400.0 } else { 500.0 };
                    entropy_input -= 0.0002 * weight * radiation[body] / baths[body];
                    radiation_generation += 0.0002 * weight * radiation[body]
                        * (1.0 / baths[body] - 1.0 / temperature);
                }
                contact_generation += 0.0002 * weight * contact * (1.0 / 500.0 - 1.0 / 400.0);
            }
            expected = pair[1];
        }
        let storage_entropy: f64 = expected.iter().zip(OLD).enumerate()
            .map(|(i, (new, old))| (new - old)
                * if i < 4 { 1.0 / (24.0 * 400.0) } else { 1.0 / (12.0 * 500.0) }).sum();
        near(report.entropy.thermal_storage_change_j_k, storage_entropy, 2e-10);
        near(report.entropy.radiation_entropy_input_j_k, entropy_input, 2e-10);
        near(report.entropy.radiation_generation_j_k, radiation_generation, 2e-10);
        near(report.entropy.conduction_generation_j_k, contact_generation, 2e-10);
        near(report.entropy.temporal_defect_j_k, 0.0, 2e-10);
        assert!(!report.entropy.has_unclassified_inputs);
        report.entropy.admit(1e-8).unwrap();
        assert_eq!(audited.conduction().accepted_steps(), 2);

        // Existing reflected exchange cannot be overwritten by a second owner.
        let enclosure = closed(&mesh, &cards, &context);
        let owner = PointwiseEnclosureRadiation::new(cx, &mesh, &initial,
            &enclosure, PointwiseEnclosureTolerance::default()).unwrap();
        let reflected = MaterialEnclosureEnthalpyHistoryRun::new(cx,
            base(cx, &boundary, &conductivity, &curves, &OLD, zero(), 2),
            &owner, patches(&cards, &curves, &context), None).unwrap().into_sdirk2();
        assert!(reflected.with_material_ambient(cx,
            ambient_patches(&cards, &curves, &context, &schedules)).is_err());
    });
}
