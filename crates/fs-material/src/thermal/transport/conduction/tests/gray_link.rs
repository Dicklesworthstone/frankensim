//! G1/G3/G4: integrated gray gaps on the original finite Cp/Cv heat stores.

#[path = "gray_link/source.rs"]
mod source;

use super::*;
use crate::thermal::conduction::GrayLayerGap;
use fs_thermal::enthalpy::ContactPatchUpdate;

const SIGMA: f64 = 5.670_374_419e-8;
const EMISSIVITIES: [f64; 2] = [0.8, 0.6];
const GAP_HTC: f64 = 35.0;

// Independent flux-density bisection using the analytic quadratic primitives
// for k_a=.04T and k_b=.02T. Neither the original integral kernel nor its root
// is used here. The fourth-power law is evaluated independently at the faces.
fn exact_gap(t: [f64; 2], distances: [f64; 2]) -> (f64, [f64; 2]) {
    let sign = if t[0] >= t[1] { 1.0 } else { -1.0 };
    let effective = 1.0 / (1.0 / EMISSIVITIES[0] + 1.0 / EMISSIVITIES[1] - 1.0);
    let faces = |j: f64| {
        [
            (t[0] * t[0] - 2.0 * sign * j * distances[0] / 0.04).sqrt(),
            (t[1] * t[1] + 2.0 * sign * j * distances[1] / 0.02).sqrt(),
        ]
    };
    let (mut lo, mut hi) = (
        0.0,
        (t[0] * t[0] - t[1] * t[1]).abs() / (2.0 * (distances[0] / 0.04 + distances[1] / 0.02)),
    );
    for _ in 0..80 {
        let j = (lo + hi) / 2.0;
        let [a, b] = faces(j);
        let gap = GAP_HTC * (a - b) + effective * SIGMA * (a.powi(4) - b.powi(4));
        if sign * gap > j {
            lo = j;
        } else {
            hi = j;
        }
    }
    let j = (lo + hi) / 2.0;
    (sign * AREA * j, faces(j))
}

#[test]
fn g1_gray_gap_matches_independent_surfaces_and_original_two_sided_caloric_tangents() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let distances = [0.01, 0.008];
    let geometry = LayerGeometry::planar(distances, AREA).unwrap();
    let thermal = timeline([&cv, &cp], [650.0, 350.0], [0.4, 0.32], AREA);
    let link = IntegratedCaloricGrayLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        GrayLayerGap::equal_area_full_view(EMISSIVITIES, GAP_HTC).unwrap(),
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(link.link_id(), LINK);
    assert_eq!(link.geometry(), geometry);
    assert_eq!(link.native_geometry().contact_resistance_m2_k_w, 0.0);
    for side in 0..2 {
        assert!(core::ptr::eq(link.layers()[side], &layers[side]));
        assert!(core::ptr::addr_eq(
            link.caloric_laws()[side] as *const dyn ThermalLaw,
            layers[side].caloric_law() as *const dyn ThermalLaw
        ));
    }
    for t in [[650.0, 350.0], [350.0, 650.0], [475.0, 475.0]] {
        let e = [energy(0, t[0]), energy(1, t[1])];
        let response = link.evaluate(e, || Ok(())).unwrap();
        let native = link
            .state_at_specific_enthalpies(
                e,
                link.native_geometry(),
                thermal.contact_patch(LINK).unwrap(),
                &mut || Ok(()),
            )
            .unwrap();
        let (expected_power, surfaces) = exact_gap(t, distances);
        near(response.power_a_to_b_w, expected_power, 3e-9);
        near(response.interface_temperature_k[0], surfaces[0], 2e-9);
        near(response.interface_temperature_k[1], surfaces[1], 2e-9);
        let effective = 1.0 / (1.0 / EMISSIVITIES[0] + 1.0 / EMISSIVITIES[1] - 1.0);
        let radiation = AREA * effective * SIGMA * (surfaces[0].powi(4) - surfaces[1].powi(4));
        let conduction = AREA * GAP_HTC * (surfaces[0] - surfaces[1]);
        near(response.radiation_power_w, radiation, 3e-9);
        near(response.gap_conduction_power_w, conduction, 3e-9);
        near(response.power_a_to_b_w, radiation + conduction, 3e-9);
        near(response.surface_heat_residual_w, 0.0, 3e-9);
        near(native.conductance_w_k * (t[0] - t[1]), expected_power, 3e-9);
        let slope = surfaces.map(|s| GAP_HTC + 4.0 * effective * SIGMA * s.powi(3));
        let denominator = 1.0
            + slope[0] * distances[0] / (0.04 * surfaces[0])
            + slope[1] * distances[1] / (0.02 * surfaces[1]);
        let expected_temperature_derivatives = [
            AREA * slope[0] * t[0] / surfaces[0] / denominator,
            -AREA * slope[1] * t[1] / surfaces[1] / denominator,
        ];
        let capacities = [200.0 + 0.5 * t[0], 400.0 + t[1]];
        for side in 0..2 {
            near(
                response.endpoint_derivative_w_k[side],
                expected_temperature_derivatives[side],
                2e-9,
            );
            near(
                native.power_enthalpy_derivatives_w_kg_j[side],
                expected_temperature_derivatives[side] / capacities[side],
                2e-12,
            );
            let mut lower = e;
            let mut upper = e;
            lower[side] -= 0.25;
            upper[side] += 0.25;
            let power =
                |e: [f64; 2]| exact_gap([temperature(0, e[0]), temperature(1, e[1])], distances).0;
            near(
                native.power_enthalpy_derivatives_w_kg_j[side],
                (power(upper) - power(lower)) / 0.5,
                2e-8,
            );
        }
    }
}

#[test]
fn g1_native_gray_gap_evolves_original_u_and_h_once_with_entropy_and_late_retry() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let (t0, masses, distances, dt) = ([650.0, 330.0], [0.4, 0.2], [0.01, 0.005], 300.0);
    let thermal = timeline([&cv, &cp], t0, masses, AREA);
    let link = IntegratedCaloricGrayLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar(distances, AREA).unwrap(),
        GrayLayerGap::equal_area_full_view(EMISSIVITIES, GAP_HTC).unwrap(),
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let observer = ConductionObserver::new(
        &thermal,
        &[
            (CV_MATERIAL, &a as &dyn ConductivityLaw),
            (CP_MATERIAL, &b as &dyn ConductivityLaw),
        ],
        &[link.native_geometry()],
        &[],
        || Ok(()),
    )
    .unwrap();
    let chain = ImplicitConductionChain::new(
        &thermal,
        &observer,
        &[
            (CV_MATERIAL, &a as &dyn ImplicitConductionLaw),
            (CP_MATERIAL, &b as &dyn ImplicitConductionLaw),
        ],
        || Ok(()),
    )
    .unwrap()
    .with_link_laws(&thermal, &[(LINK, &link)], || Ok(()))
    .unwrap();
    let controls = ImplicitConductionControls {
        maximum_dt_s: dt,
        maximum_halvings: 0,
        cell_tolerance_j: 1e-6,
        total_tolerance_j: 1e-6,
        ..ImplicitConductionControls::default()
    };
    let mut workspace = ImplicitConductionWorkspace::new(2).unwrap();
    let (next, report) = chain
        .advance_to(&thermal, dt, controls, &mut workspace, || Ok(()))
        .unwrap();
    // A separate scalar total-heat balance determines the backward-Euler
    // endpoint without invoking the production thermal residual or flux root.
    let e0 = [energy(0, t0[0]), energy(1, t0[1])];
    let residual = |q: f64| {
        let t = [
            temperature(0, e0[0] - q / masses[0]),
            temperature(1, e0[1] + q / masses[1]),
        ];
        q - dt * exact_gap(t, distances).0
    };
    let (mut lo, mut hi) = (0.0, masses[0] * (e0[0] - energy(0, t0[1])));
    assert!(residual(lo) < 0.0 && residual(hi) > 0.0);
    for _ in 0..80 {
        let q = (lo + hi) / 2.0;
        if residual(q) < 0.0 {
            lo = q;
        } else {
            hi = q;
        }
    }
    let q = (lo + hi) / 2.0;
    let e = [
        next.cells()[&CellId(0)].specific_h_j_kg,
        next.cells()[&CellId(1)].specific_h_j_kg,
    ];
    near(e[0], e0[0] - q / masses[0], 5e-6);
    near(e[1], e0[1] + q / masses[1], 5e-6);
    near(
        masses[0] * (e[0] - e0[0]) + masses[1] * (e[1] - e0[1]),
        0.0,
        1e-6,
    );
    let t = [temperature(0, e[0]), temperature(1, e[1])];
    assert!(t[0] < t0[0] && t[1] > t0[1] && t[0] > t[1]);
    let entropy = masses[0] * (200.0 * (t[0] / t0[0]).ln() + 0.5 * (t[0] - t0[0]))
        + masses[1] * (400.0 * (t[1] / t0[1]).ln() + t[1] - t0[1]);
    near(report.entropy.thermal_storage_change_j_k, entropy, 2e-8);
    near(
        report.entropy.conduction_generation_j_k,
        q * (1.0 / t[1] - 1.0 / t[0]),
        2e-8,
    );
    assert!(entropy > 0.0 && report.energy.residual_j.abs() <= 1e-6);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(
        next.material_storage_basis(CV_MATERIAL).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(
        next.material_storage_basis(CP_MATERIAL).unwrap(),
        ThermalStorageBasis::Enthalpy
    );
    let response = link.evaluate(e, || Ok(())).unwrap();
    assert!(response.gap_conduction_power_w > 0.0 && response.radiation_power_w > 0.0);
    near(dt * response.power_a_to_b_w, q, 1e-6);
    // A late cancellation starts from an already accepted positive time and
    // leaves both original caloric stores intact; the same workspace replays.
    let mut calls = 0;
    let (continued, receipt) = chain
        .advance_to(&next, 2.0 * dt, controls, &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert!(calls > 8);
    let mut cancelled_calls = 0;
    let refusal = chain.advance_to(&next, 2.0 * dt, controls, &mut workspace, || {
        cancelled_calls += 1;
        if cancelled_calls == calls - 1 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        refusal,
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(next.time_s(), dt);
    assert_eq!(next.cells()[&CellId(0)].specific_h_j_kg, e[0]);
    assert_eq!(next.cells()[&CellId(1)].specific_h_j_kg, e[1]);
    let (retry, repeated) = chain
        .advance_to(&next, 2.0 * dt, controls, &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), continued.cells());
    assert_eq!(retry.links(), continued.links());
    assert_eq!(repeated, receipt);
}

#[test]
fn g4_gray_gap_admission_keeps_disjoint_sources_and_exact_insulating_conductive_limits() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let thermal = timeline([&cv, &cp], [650.0, 350.0], [0.4, 0.32], AREA);
    let geometry = LayerGeometry::planar([0.01, 0.008], AREA).unwrap();
    let gap = GrayLayerGap::equal_area_full_view(EMISSIVITIES, GAP_HTC).unwrap();
    let build = |thermal, layers, geometry, gap| {
        IntegratedCaloricGrayLink::new(
            thermal,
            LINK,
            layers,
            geometry,
            gap,
            layer_controls(),
            || Ok(()),
        )
    };
    assert!(build(&thermal, [&layers[1], &layers[0]], geometry, gap).is_err());
    assert!(
        build(
            &thermal,
            [&layers[0], &layers[1]],
            LayerGeometry::planar([0.01, 0.008], 2.0 * AREA).unwrap(),
            gap,
        )
        .is_err()
    );
    // Give a radial geometry its exact area, so refusal proves the view model
    // is unsupported instead of merely detecting an unrelated area mismatch.
    let radial = LayerGeometry::spherical([0.01, 0.03], 0.02).unwrap();
    let radial_thermal = timeline(
        [&cv, &cp],
        [650.0, 350.0],
        [0.4, 0.32],
        radial.contact_area_m2(),
    );
    assert!(build(&radial_thermal, [&layers[0], &layers[1]], radial, gap).is_err());
    let original_patch = thermal.contact_patch(LINK).unwrap();
    let mut asymmetric = original_patch;
    asymmetric.covered_area_m2[1] *= 0.5;
    let unequal = thermal
        .replace_contact_patches(
            &[ContactPatchUpdate {
                id: LINK,
                expected: Some(original_patch),
                replacement: Some(asymmetric),
            }],
            || Ok(()),
        )
        .unwrap();
    assert!(build(&unequal, [&layers[0], &layers[1]], geometry, gap).is_err());
    let link = build(&thermal, [&layers[0], &layers[1]], geometry, gap).unwrap();
    let e = [energy(0, 650.0), energy(1, 350.0)];
    let mut extra_contact = link.native_geometry();
    extra_contact.contact_resistance_m2_k_w = 0.001;
    assert!(
        link.state_at_specific_enthalpies(e, extra_contact, original_patch, &mut || Ok(()))
            .is_err()
    );
    assert!(
        link.state_at_specific_enthalpies(e, link.native_geometry(), asymmetric, &mut || Ok(()))
            .is_err()
    );
    assert!(link.evaluate([-1.0, e[1]], || Ok(())).is_err());
    let insulating = build(
        &thermal,
        [&layers[0], &layers[1]],
        geometry,
        GrayLayerGap::equal_area_full_view([0.0, 0.8], 0.0).unwrap(),
    )
    .unwrap();
    let zero = insulating.evaluate(e, || Ok(())).unwrap();
    assert_eq!(zero.power_a_to_b_w, 0.0);
    assert_eq!(zero.conductance_w_k, 0.0);
    assert_eq!(zero.endpoint_derivative_w_k, [0.0; 2]);
    assert!(insulating.evaluate([-1.0, e[1]], || Ok(())).is_err());
    let conductive = build(
        &thermal,
        [&layers[0], &layers[1]],
        geometry,
        GrayLayerGap::equal_area_full_view([0.0, 0.8], GAP_HTC).unwrap(),
    )
    .unwrap();
    let contact = IntegratedCaloricLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        1.0 / GAP_HTC,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let conduction = conductive.evaluate(e, || Ok(())).unwrap();
    let original = contact.evaluate(e, || Ok(())).unwrap();
    assert_eq!(conduction.power_a_to_b_w, original.power_a_to_b_w);
    assert_eq!(conduction.conductance_w_k, original.conductance_w_k);
    assert_eq!(
        conduction.endpoint_derivative_w_k,
        original.endpoint_derivative_w_k
    );

    // Real material support is disjoint. The finite gap keeps each material
    // face inside its own chart; an intermediate perfect-contact solve fails.
    let lower_cv = fa
        .lower
        .clone()
        .with_quantity("T", fa.lower.axis_quantities()["T"], 600.0)
        .unwrap();
    let upper_cp = fb
        .upper
        .clone()
        .with_quantity("T", fb.upper.axis_quantities()["T"], 400.0)
        .unwrap();
    let narrow_cv = IsochoricHeatCapacityLaw::try_new(
        &fa.card,
        IsochoricHeatCapacityQuery {
            specific_heat_cv: &fa.capacity,
            temperature_axis: "T",
            density_axis: "rho",
            lower: &lower_cv,
            upper: &fa.upper,
            selection: SELECTION,
        },
        650.0,
    )
    .unwrap();
    let narrow_cp = IsobaricHeatCapacityLaw::try_new(
        &fb.card,
        IsobaricHeatCapacityQuery {
            specific_heat_cp: &fb.capacity,
            temperature_axis: "T",
            pressure_axis: "p",
            lower: &fb.lower,
            upper: &upper_cp,
            selection: SELECTION,
        },
        350.0,
    )
    .unwrap();
    let narrow_a =
        IsochoricConductionLaw::try_new(&fa.card, &narrow_cv, &fa.conductivity, SELECTION).unwrap();
    let narrow_b =
        IsobaricConductionLaw::try_new(&fb.card, &narrow_cp, &fb.conductivity, SELECTION).unwrap();
    let narrow_layers = [
        CaloricConductionLayer::from_isochoric(&narrow_a).unwrap(),
        CaloricConductionLayer::from_isobaric(&narrow_b).unwrap(),
    ];
    let mut events: Vec<_> = thermal
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.specific_h_j_kg = 0.0;
            TimedEvent {
                time_s: 0.0,
                event: Event::Birth(cell),
            }
        })
        .collect();
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(original_patch.link),
    });
    let narrow = DepositionSimulation::new(
        events,
        &[
            (CV_MATERIAL, &narrow_cv as &dyn ThermalLaw),
            (CP_MATERIAL, &narrow_cp as &dyn ThermalLaw),
        ],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 8,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-8,
        },
        || Ok(()),
    )
    .unwrap()
    .0;
    assert!(
        IntegratedCaloricLink::new(
            &narrow,
            LINK,
            [&narrow_layers[0], &narrow_layers[1]],
            geometry,
            0.0,
            layer_controls(),
            || Ok(())
        )
        .is_err()
    );
    let admitted = build(
        &narrow,
        [&narrow_layers[0], &narrow_layers[1]],
        geometry,
        gap,
    )
    .unwrap();
    let response = admitted.evaluate([0.0; 2], || Ok(())).unwrap();
    near(
        response.power_a_to_b_w,
        exact_gap([650.0, 350.0], [0.01, 0.008]).0,
        3e-9,
    );
    assert!(response.interface_temperature_k[0] >= 600.0);
    assert!(response.interface_temperature_k[1] <= 400.0);
    assert!(
        admitted
            .evaluate([energy(0, 599.0) - energy(0, 650.0), 0.0], || Ok(()))
            .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.contact_patch(LINK), Some(original_patch));
}
