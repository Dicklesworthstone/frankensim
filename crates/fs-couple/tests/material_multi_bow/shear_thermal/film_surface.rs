//! Fixed material film, finite sliding receiver, original coupled heat owner.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::multiple::layer::thermal::{
    dry::DryThermalCell, film_surface::RetainedFilmHeat,
};
use fs_material::thermal::source::ThermalMaterialPoint;

const LENGTH: f64 = 0.002;
const WIDTH: f64 = 0.001;

fn network(
    curve: &EquilibriumEnthalpyPhaseCurve,
    receiver_k: [f64; 2],
) -> DepositionSimulation<'_> {
    let original = thermal(curve, true);
    let events = original
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            if cell.id != FILM {
                cell.specific_h_j_kg = CP * (receiver_k[(cell.id.0 - 11) as usize] - 300.0);
            }
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    DepositionSimulation::new(
        events,
        &[(MaterialId(0), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 8,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1e-5,
            energy_tolerance_j: 1e-14,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}

fn descriptor<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    layer: &'a ResolvedShearLayer,
    context: &'a QueryPoint,
    index: usize,
) -> DryThermalCell<'a> {
    let origin = if index == 2 { LENGTH } else { 0.0 };
    DryThermalCell {
        cell: CellId(10 + index as u64),
        material: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        caloric: curve,
        transport: layer.resolved(),
        conductivity_property: "thermal_conductivity",
        density_property: "density",
        interval_m: [origin, origin + LENGTH],
        width_m: WIDTH,
        half_depth_m: if index == 0 {
            0.5 * THICKNESS
        } else {
            5.0 * THICKNESS
        },
    }
}

fn transfer<'a>(
    card: &'a MaterialCard,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    layer: &'a ResolvedShearLayer,
    context: &'a QueryPoint,
    heat: &DepositionSimulation<'_>,
) -> RetainedFilmHeat<'a> {
    RetainedFilmHeat::try_new(
        layer,
        binding(card, curve)[0],
        descriptor(card, curve, layer, context, 0),
        vec![
            descriptor(card, curve, layer, context, 1),
            descriptor(card, curve, layer, context, 2),
        ],
        heat,
        3,
    )
    .unwrap()
}

#[test]
fn retained_film_uses_series_material_resistance_and_actual_overlap_without_duplicate_links() {
    let (card, caloric) = coating(true, false);
    let layer = layer(&card);
    let context = point(300.0);
    let heat = network(caloric.curve(), [290.0, 295.0]);
    let moving = transfer(&card, caloric.curve(), &layer, &context, &heat);
    // Independent slab resistances: h_f/(2 k_f A) + h_b/(2 k_b A).
    let full_g = AREA / (0.5 * THICKNESS / 0.2 + 5.0 * THICKNESS / 0.2);
    let first = moving.heat_rates_w(&heat, [0.0, LENGTH], WIDTH).unwrap();
    assert!((first.sources_w[&CellId(11)] - full_g * 10.0).abs() < 1e-14);
    assert!(!first.sources_w.contains_key(&CellId(12)));
    assert_eq!(first.sources_w[&FILM], -first.sources_w[&CellId(11)]);
    let second = moving
        .heat_rates_w(&heat, [LENGTH, 2.0 * LENGTH], WIDTH)
        .unwrap();
    assert!((second.sources_w[&CellId(12)] - full_g * 5.0).abs() < 1e-14);
    assert!(!second.sources_w.contains_key(&CellId(11)));
    let split = moving
        .heat_rates_w(&heat, [0.5 * LENGTH, 1.5 * LENGTH], WIDTH)
        .unwrap();
    assert!((split.sources_w[&CellId(11)] - 0.5 * full_g * 10.0).abs() < 1e-14);
    assert!((split.sources_w[&CellId(12)] - 0.5 * full_g * 5.0).abs() < 1e-14);
    assert!(split.internal_transfer_residual_w.abs() < 1e-15);
    assert_eq!(
        moving.heat_rates_w(&heat, [0.0, LENGTH], WIDTH).unwrap(),
        first
    );
    assert!(moving.heat_rates_w(&heat, [-LENGTH, 0.0], WIDTH).is_err());
    assert!(
        moving
            .heat_rates_w(&heat, [0.0, LENGTH], 2.0 * WIDTH)
            .is_err()
    );
    let already_linked = thermal(caloric.curve(), true);
    assert!(
        RetainedFilmHeat::try_new(
            &layer,
            binding(&card, caloric.curve())[0],
            descriptor(&card, caloric.curve(), &layer, &context, 0),
            vec![descriptor(&card, caloric.curve(), &layer, &context, 1)],
            &already_linked,
            3
        )
        .is_err()
    );
    assert!(
        RetainedFilmHeat::try_new(
            &layer,
            binding(&card, caloric.curve())[0],
            descriptor(&card, caloric.curve(), &layer, &context, 0),
            vec![descriptor(&card, caloric.curve(), &layer, &context, 0)],
            &heat,
            3
        )
        .is_err()
    );
}

#[test]
fn actual_trial_slip_moves_film_heat_and_changes_traction_without_resetting_history() {
    with_cx(|cx, gate| {
        let (card, caloric) = coating(true, false);
        let original_layer = layer(&card);
        let context = point(300.0);
        let (_, mut moving) = instrument(cx, gate, &card).into_parts();
        let (_, mut stationary) = instrument(cx, gate, &card).into_parts();
        let mut heat = network(caloric.curve(), [300.0, 290.0]);
        let mut comparison = heat.clone();
        let transfer = transfer(&card, caloric.curve(), &original_layer, &context, &heat);
        let bound = binding(&card, caloric.curve());
        let forces = vec![0.0; moving.mode_count()];
        let drives = [drive(0.1, 0.1), drive(0.0, 0.0)];
        let initial_heat = total_heat(&heat);
        let (mut plastic, mut pressure_difference) = (0.0, 0.0);
        for _ in 0..128 {
            let f = moving
                .step_with_contact_thermal_feedback(
                    &mut heat,
                    &bound,
                    &[],
                    controls(3),
                    &[0.0; 3],
                    &forces,
                    &drives,
                    None,
                    None,
                    gate,
                    |endpoint, ports| {
                        let film = ports[0].shear_layer.expect("same staged layer");
                        let slip = film.elastic_slip_m + film.plastic_slip_m;
                        let rates =
                            transfer.heat_rates_w(endpoint, [slip, slip + LENGTH], WIDTH)?;
                        Ok::<_, LayerThermalError>(StepCoefficients {
                            sources_w: rates.sources_w,
                            conductances_w_k: BTreeMap::new(),
                            payload: (slip, rates.internal_transfer_residual_w),
                        })
                    },
                )
                .unwrap();
            let b = stationary
                .step_with_contact_thermal_feedback(
                    &mut comparison,
                    &bound,
                    &[],
                    controls(3),
                    &[0.0; 3],
                    &forces,
                    &drives,
                    None,
                    None,
                    gate,
                    |endpoint, _| {
                        let rates = transfer.heat_rates_w(endpoint, [0.0, LENGTH], WIDTH)?;
                        Ok::<_, LayerThermalError>(StepCoefficients {
                            sources_w: rates.sources_w,
                            conductances_w_k: BTreeMap::new(),
                            payload: (),
                        })
                    },
                )
                .unwrap();
            let accepted = f.mechanics.ports[0].shear_layer.unwrap();
            assert_eq!(
                f.source.0,
                accepted.elastic_slip_m + accepted.plastic_slip_m
            );
            assert!(f.source.1.abs() < 1e-14);
            assert!(f.total_energy_residual_j.abs() <= controls(3).total_energy_tolerance_j);
            assert!(f.external_heat_j.abs() < 1e-15);
            plastic += f.routed_plastic_heat_j;
            pressure_difference +=
                (f.mechanics.observer_pressure_pa - b.mechanics.observer_pressure_pa).abs();
        }
        assert!(heat.temperature_k(CellId(12)).unwrap() > 290.0);
        assert_eq!(comparison.temperature_k(CellId(12)).unwrap(), 290.0);
        assert!(heat.temperature_k(FILM).unwrap() < comparison.temperature_k(FILM).unwrap());
        assert!(pressure_difference > 1e-14);
        assert!((total_heat(&heat) - initial_heat - plastic).abs() < 1e-10);
        assert!(moving.material_shear_state(0).unwrap().plastic_slip_m()[0] > 0.0);
        let accepted = moving.material_shear_state(0).unwrap().clone();
        let accepted_h = total_heat(&heat);
        let accepted_time = heat.time_s();
        let fail = moving.step_with_contact_thermal_feedback(
            &mut heat,
            &bound,
            &[],
            controls(3),
            &[0.0; 3],
            &forces,
            &drives,
            None,
            None,
            gate,
            |endpoint, _| {
                let rates = transfer.heat_rates_w(endpoint, [-LENGTH, 0.0], WIDTH)?;
                Ok::<_, LayerThermalError>(StepCoefficients {
                    sources_w: rates.sources_w,
                    conductances_w_k: BTreeMap::new(),
                    payload: (),
                })
            },
        );
        assert!(fail.is_err());
        assert_eq!(*moving.material_shear_state(0).unwrap(), accepted);
        assert_eq!(total_heat(&heat), accepted_h);
        assert_eq!(heat.time_s(), accepted_time);
        assert_eq!(moving.samples_rendered(), 128);
    });
}
