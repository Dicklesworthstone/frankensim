//! G1/G3: both storage coordinates consume the same original P1 flux control.
//! Independent energy and one-free-node recurrences fix its sign and dt scale.

use super::*;
use fs_project::TransientRegionCapacity;

const DFLUX: &str = "dtemperature_doutward_flux_k_per_w_m2";

fn gradient(solved: &ThermalRun) -> &J {
    solved.transient().get("adjoint").unwrap()
}

fn flux<'a>(solved: &'a ThermalRun, target: &str) -> &'a J {
    let result = gradient(solved);
    let controls = result.get("prescribed_heat_flux").unwrap();
    close(n(controls, "control_steps"), n(result, "state_index"), 0.0);
    assert_eq!(controls.str_field("parameter_scope"),
        Some("independent-uniform-outward-heat-fluxes-held-throughout-selected-trajectory"));
    controls.get("boundaries").unwrap().as_array().unwrap().iter()
        .find(|row| row.str_field("target") == Some(target)).unwrap()
}

fn capacity_storage(spec: &mut ProjectSpec, initial_k: f64) {
    let policy = time(spec);
    policy.enthalpy = None;
    policy.initial_temperature = TransientInitialTemperature::Uniform(
        QtyAny::new(initial_k, dims::TEMPERATURE));
    policy.capacities = vec![TransientRegionCapacity {
        region: "air".into(),
        volumetric_heat_capacity: QtyAny::new(100.0, dims::VOLUMETRIC_HEAT_CAPACITY),
        source: "synthetic rhoCp matching the sensible enthalpy reference".into(),
    }];
}

fn set_flux(spec: &mut ProjectSpec, target: &str, value: f64) {
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries
        .iter_mut().find(|row| row.target == target).unwrap().condition =
        ThermalBoundaryCondition::HeatFlux {
            outward_flux: QtyAny::new(value, dims::HEAT_FLUX),
        };
}

#[test]
fn g1_native_flux_trajectory_controls_match_closed_energy_in_both_storage_coordinates() {
    for enthalpy in [false, true] {
        let (mut spec, bytes) = uniform_project(10.0, 0.0, 40.0, false);
        if !enthalpy { capacity_storage(&mut spec, 300.0); }
        time(&mut spec).horizon.value = 2.1;
        time(&mut spec).max_step.value = 0.6;
        set_flux(&mut spec, "air", 0.7);
        let cards = fixture_cards();
        let solved = run(&spec, Some(&bytes), &cards, Some(FINAL));
        let plain = run(&spec, Some(&bytes), &cards, None);
        assert_eq!(solved.temperature, plain.temperature);
        assert_eq!(solved.transient().get("fine"), plain.transient().get("fine"));
        // Four equal faces of the regular tetrahedron: area=8 sqrt(3),
        // volume=8/3. Uniform heating is an exact spatial mode of both owners.
        let area = 8.0 * fs_math::det::sqrt(3.0);
        let expected = -area * 2.1 / (100.0 * VOLUME);
        close(n(flux(&solved, "air"), "flux_w_m2"), 0.7, 0.0);
        close(n(flux(&solved, "air"), DFLUX), expected, 2e-9);
        close(solved.final_temperature(), 300.0 + 2.1 * (40.0 - area * 0.7)
            / (100.0 * VOLUME), 2e-8);
        let evaluate = |value| {
            let mut changed = spec.clone();
            set_flux(&mut changed, "air", value);
            run(&changed, Some(&bytes), &cards, None).final_temperature()
        };
        close(n(flux(&solved, "air"), DFLUX),
            (evaluate(0.71) - evaluate(0.69)) / 0.02, 2e-8);
    }
}

fn partially_prescribed(enthalpy: bool) -> (ProjectSpec, Vec<u8>, CardPackSet) {
    let bytes = tetra_stl();
    let mut spec = transient_product::driven_face_project(&bytes);
    let cards = transient_product::transient_conductivity_cards("constant-scalar");
    install_storage(&mut spec, vec![material("air", 10.0, 300.0, false)]);
    rebind_to(&mut spec, &cards);
    time(&mut spec).enthalpy.as_mut().unwrap().materials[0].material_card_identity =
        fs_blake3::ContentHash::from_hex(&cards.materials()[0].card().to_hex()).unwrap();
    if !enthalpy { capacity_storage(&mut spec, 330.0); }
    time(&mut spec).horizon.value = 0.7;
    spec.power.as_mut().unwrap()[0].watts.value = 4.0;
    spec.power.as_mut().unwrap()[0].duty = 1.0;
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries[0].condition =
        ThermalBoundaryCondition::FixedTemperature {
            temperature: QtyAny::new(295.0, dims::TEMPERATURE),
        };
    spec.assembly.as_mut().unwrap().push(EntityDecl::Surface {
        parent: "enclosure".into(), name: "flux-face".into(),
        display: "Independent original flux trace".into(), expect_id: None,
    });
    spec.assignments.as_mut().unwrap().push(GeometryAssignment {
        artifact: "enclosure".into(), target: "flux-face".into(), length_unit: "m".into(),
        selector: MeshSelector::HalfSpace {
            normal: [1.0, 0.0, 0.0], offset: 0.0,
            side: HalfSpaceSide::AtMost, tolerance: 0.0,
        }, allow_overlap: true,
    });
    spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries.push(
        ThermalBoundary { target: "flux-face".into(),
            condition: ThermalBoundaryCondition::HeatFlux {
                outward_flux: QtyAny::new(0.4, dims::HEAT_FLUX),
            },
        });
    (spec, bytes, cards)
}

#[test]
fn g1_g3_native_flux_control_masks_fixed_equations_and_initial_peak_future_steps() {
    for enthalpy in [false, true] {
        let (spec, bytes, cards) = partially_prescribed(enthalpy);
        let solved = run(&spec, Some(&bytes), &cards, Some(FINAL));
        close(n(gradient(&solved), "active_ties"), 1.0, 0.0);
        // Only the origin is free. Its lumped capacity is rhoCp/24;
        // the x=0 face injects -q*(1/2)/3=-q/6 at that free node.
        let capacity = 100.0 / 24.0;
        let mut expected = 0.0;
        let mut widths = Vec::new();
        for interval in gradient(&solved).get("intervals").unwrap().as_array().unwrap() {
            let dt = n(interval, "end_s") - n(interval, "start_s");
            widths.push(dt);
            expected = (capacity * expected - dt / 6.0) / (capacity + dt);
        }
        assert!(widths.iter().any(|dt| (dt - widths[0]).abs() > 0.01));
        close(n(flux(&solved, "flux-face"), DFLUX), expected, 2e-9);
        let evaluate = |value| {
            let mut changed = spec.clone();
            set_flux(&mut changed, "flux-face", value);
            run(&changed, Some(&bytes), &cards, None).final_temperature()
        };
        close(n(flux(&solved, "flux-face"), DFLUX),
            (evaluate(0.41) - evaluate(0.39)) / 0.02, 2e-8);
        let peak = run(&spec, Some(&bytes), &cards, Some(PEAK));
        close(n(gradient(&peak), "state_index"), 0.0, 0.0);
        close(n(flux(&peak, "flux-face"), DFLUX), 0.0, 0.0);
    }
}

fn surface_power(solved: &ThermalRun, source_index: usize) -> &J {
    let result = gradient(solved);
    let controls = result.get("surface_power").unwrap();
    close(n(controls, "control_steps"), n(result, "state_index"), 0.0);
    controls.get("sources").unwrap().as_array().unwrap().iter()
        .find(|row| n(row, "source_index") == source_index as f64).unwrap()
}

#[test]
fn g1_g3_native_surface_power_controls_retain_watts_duties_and_physical_trace() {
    for enthalpy in [false, true] {
        let (mut spec, bytes, cards) = partially_prescribed(enthalpy);
        // Reuse the x=0 face's original assignment, with three prescribed
        // vertices elsewhere and the origin as the only free equation.
        spec.cooling.as_mut().unwrap().conduction.as_mut().unwrap().boundaries
            .retain(|row| row.target != "flux-face");
        let first = spec.power.as_ref().unwrap().len();
        for (watts, duty) in [(7.0, 0.4), (11.0, 0.2), (0.0, 0.6), (13.0, 0.0)] {
            spec.power.as_mut().unwrap().push(fs_project::PowerDissipation {
                region: "flux-face".into(), watts: QtyAny::new(watts, dims::POWER), duty,
            });
        }
        let solved = run(&spec, Some(&bytes), &cards, Some(FINAL));
        let plain = run(&spec, Some(&bytes), &cards, None);
        assert_eq!(solved.temperature, plain.temperature);
        assert_eq!(solved.transient().get("fine"), plain.transient().get("fine"));
        close(n(gradient(&solved), "active_ties"), 1.0, 0.0);
        // Surface total watts inject P/3 at the free node, while volume
        // watts inject P/4. C=rhoCp*V/4 and K=1 for this tetrahedron.
        let capacity = 100.0 / 24.0;
        let mut expected = 0.0;
        for interval in gradient(&solved).get("intervals").unwrap().as_array().unwrap() {
            let dt = n(interval, "end_s") - n(interval, "start_s");
            expected = (capacity * expected + dt / 3.0) / (capacity + dt);
        }
        let volume = region(gradient(&solved), "regional_power", "air");
        close(n(volume, "dtemperature_duniform_added_watts_k_per_w"),
            expected * 0.75, 2e-9);
        for (offset, power) in spec.power.as_ref().unwrap()[first..].iter().enumerate() {
            let row = surface_power(&solved, first + offset);
            assert_eq!(row.str_field("target"), Some("flux-face"));
            close(n(row, "area_m2"), 0.5, 0.0);
            close(n(row, "watts"), power.watts.value, 0.0);
            close(n(row, "duty"), power.duty, 0.0);
            close(n(row, "dtemperature_ddelivered_watts_k_per_w"), expected, 2e-9);
            close(n(row, "dtemperature_dwatts_k_per_w"), expected * power.duty, 2e-9);
            close(n(row, "dtemperature_dduty_k"), expected * power.watts.value, 2e-9);
        }
        let evaluate = |source: usize, watts: Option<f64>, duty: Option<f64>| {
            let mut changed = spec.clone();
            let power = &mut changed.power.as_mut().unwrap()[source];
            if let Some(watts) = watts { power.watts.value = watts; }
            if let Some(duty) = duty { power.duty = duty; }
            run(&changed, Some(&bytes), &cards, None).final_temperature()
        };
        close(n(surface_power(&solved, first), "dtemperature_dwatts_k_per_w"),
            (evaluate(first, Some(7.01), None) - evaluate(first, Some(6.99), None)) / 0.02,
            2e-8);
        close(n(surface_power(&solved, first + 1), "dtemperature_dduty_k"),
            (evaluate(first + 1, None, Some(0.201)) - evaluate(first + 1, None, Some(0.199))) / 0.002,
            2e-8);
        // Independent feasible perturbations stay informative at zero watts
        // and at zero duty; no derivative divides by either baseline.
        close(n(surface_power(&solved, first + 2), "dtemperature_dwatts_k_per_w"),
            (evaluate(first + 2, Some(0.01), None) - plain.final_temperature()) / 0.01,
            2e-8);
        close(n(surface_power(&solved, first + 3), "dtemperature_dduty_k"),
            (evaluate(first + 3, None, Some(0.001)) - plain.final_temperature()) / 0.001,
            2e-8);
        // All sources on the trace can be off together. The forward owner
        // retains -0.0 outward flux from -delivered/area, and the reverse
        // normalization must admit that same trace without dividing by power.
        let mut zero = spec.clone();
        for power in &mut zero.power.as_mut().unwrap()[first..] {
            power.duty = 0.0;
        }
        let zero = run(&zero, Some(&bytes), &cards, Some(FINAL));
        close(n(surface_power(&zero, first), "dtemperature_ddelivered_watts_k_per_w"),
            expected, 2e-9);
        close(n(surface_power(&zero, first), "dtemperature_dwatts_k_per_w"), 0.0, 0.0);
        close(n(surface_power(&zero, first), "dtemperature_dduty_k"),
            expected * 7.0, 2e-9);
        let peak = run(&spec, Some(&bytes), &cards, Some(PEAK));
        close(n(gradient(&peak), "state_index"), 0.0, 0.0);
        for source in first..spec.power.as_ref().unwrap().len() {
            let row = surface_power(&peak, source);
            for key in ["dtemperature_ddelivered_watts_k_per_w",
                "dtemperature_dwatts_k_per_w", "dtemperature_dduty_k"]
            {
                close(n(row, key), 0.0, 0.0);
            }
        }
    }
}
