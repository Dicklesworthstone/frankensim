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
