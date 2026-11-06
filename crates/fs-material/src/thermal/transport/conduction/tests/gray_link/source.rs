//! Original source-dependent nonradiative gap resistance, beside gray radiation.

use super::*;
use fs_thermal::enthalpy::conduction::implicit::contact::{
    ImplicitContactLaw, ImplicitContactState,
};
use std::cell::Cell;

#[derive(Clone)]
struct GapSource<'a> {
    laws: [&'a dyn ThermalLaw; 2],
    bounds: [[f64; 2]; 2],
    reference_r: f64,
    slopes_per_k: [f64; 2],
    maximum_difference_k: f64,
    queries: Cell<usize>,
}

impl GapSource<'_> {
    fn resistance(&self, t: [f64; 2]) -> f64 {
        self.reference_r
            * (self.slopes_per_k[0] * (t[0] - 650.0) + self.slopes_per_k[1] * (t[1] - 350.0)).exp()
    }
}

impl ImplicitContactLaw for GapSource<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.laws
    }
    fn enthalpy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        self.bounds
    }
    fn state_at_specific_enthalpies(
        &self,
        e: [f64; 2],
    ) -> Result<ImplicitContactState, DepositionError> {
        self.queries.set(self.queries.get() + 1);
        let t = [
            self.laws[0].temperature_k(e[0])?,
            self.laws[1].temperature_k(e[1])?,
        ];
        if (t[0] - t[1]).abs() > self.maximum_difference_k {
            return Err(DepositionError::MaterialDomain);
        }
        let r = self.resistance(t);
        Ok(ImplicitContactState {
            resistance_m2_k_w: r,
            resistance_enthalpy_derivatives_m2_kg_k_w_j: [
                r * self.slopes_per_k[0] / (200.0 + 0.5 * t[0]),
                r * self.slopes_per_k[1] / (400.0 + t[1]),
            ],
        })
    }
}

fn source<'a>(
    cv: &'a IsochoricHeatCapacityLaw,
    cp: &'a IsobaricHeatCapacityLaw,
    slopes: [f64; 2],
) -> GapSource<'a> {
    GapSource {
        laws: [cv, cp],
        bounds: [
            cv.specific_internal_energy_span_j_kg(),
            cp.specific_enthalpy_span_j_kg(),
        ],
        reference_r: 0.004,
        slopes_per_k: slopes,
        maximum_difference_k: f64::INFINITY,
        queries: Cell::new(0),
    }
}

// Exact quadratic material primitives, independently bisected in flux density.
// R remains explicit so its partial and its original caloric chains can both
// be checked without invoking the production Kirchhoff or gray-gap solver.
fn power(t: [f64; 2], d: [f64; 2], r: f64) -> f64 {
    let sign = if t[0] >= t[1] { 1.0 } else { -1.0 };
    let epsilon = 1.0 / (1.0 / EMISSIVITIES[0] + 1.0 / EMISSIVITIES[1] - 1.0);
    let (mut lo, mut hi) = (
        0.0,
        (t[0] * t[0] - t[1] * t[1]).abs() / (2.0 * (d[0] / 0.04 + d[1] / 0.02)),
    );
    for _ in 0..80 {
        let j = (lo + hi) / 2.0;
        let a = (t[0] * t[0] - 2.0 * sign * j * d[0] / 0.04).sqrt();
        let b = (t[1] * t[1] + 2.0 * sign * j * d[1] / 0.02).sqrt();
        if sign * ((a - b) / r + epsilon * SIGMA * (a.powi(4) - b.powi(4))) > j {
            lo = j;
        } else {
            hi = j;
        }
    }
    sign * AREA * (lo + hi) / 2.0
}

#[test]
fn g1_original_gap_resistance_keeps_radiative_feedback_and_signed_caloric_tangents() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let thermal = timeline([&cv, &cp], [650.0, 350.0], [0.4, 0.32], AREA);
    let d = [0.01, 0.008];
    let geometry = LayerGeometry::planar(d, AREA).unwrap();
    let law = source(&cv, &cp, [0.02, 0.005]);
    let gap = GrayLayerGap::equal_area_full_view(EMISSIVITIES, 0.0).unwrap();
    let link = IntegratedCaloricGrayLink::new_with_gap_resistance(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        gap,
        &law,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(law.queries.get(), 1);
    for t in [[650.0, 350.0], [350.0, 650.0], [475.0, 475.0]] {
        let e = [energy(0, t[0]), energy(1, t[1])];
        let before = law.queries.get();
        let response = link.evaluate(e, || Ok(())).unwrap();
        assert_eq!(law.queries.get(), before + 1);
        let native = link
            .state_at_specific_enthalpies(
                e,
                link.native_geometry(),
                thermal.contact_patch(LINK).unwrap(),
                &mut || Ok(()),
            )
            .unwrap();
        assert_eq!(law.queries.get(), before + 2);
        let r = law.resistance(t);
        near(response.power_a_to_b_w, power(t, d, r), 2e-8);
        let dr = r * 1e-5;
        let partial = (power(t, d, r + dr) - power(t, d, r - dr)) / (2.0 * dr);
        near(
            response.gap_resistance_derivative_w2_m2_k,
            partial,
            2e-7 * partial.abs().max(1.0),
        );
        for side in 0..2 {
            let mut lo = e;
            lo[side] -= 0.25;
            let mut hi = e;
            hi[side] += 0.25;
            let exact = |e: [f64; 2]| {
                let t = [temperature(0, e[0]), temperature(1, e[1])];
                power(t, d, law.resistance(t))
            };
            near(
                native.power_enthalpy_derivatives_w_kg_j[side],
                (exact(hi) - exact(lo)) / 0.5,
                3e-8,
            );
        }
        if t[0] == 650.0 {
            assert!(response.power_a_to_b_w > 0.0 && native.conductance_w_k > 0.0);
            assert!(response.endpoint_derivative_w_k[0] > 0.0);
            assert!(native.power_enthalpy_derivatives_w_kg_j[0] < 0.0);
        }
    }
    // The authored R must survive the pure-conduction path bit-for-bit: this
    // number changes by one ULP under R -> 1/R -> R.
    let mut plain_source = source(&cv, &cp, [0.0; 2]);
    plain_source.reference_r = 0.013;
    assert_ne!(
        plain_source.reference_r,
        1.0 / (1.0 / plain_source.reference_r)
    );
    let pure = IntegratedCaloricGrayLink::new_with_gap_resistance(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        GrayLayerGap::equal_area_full_view([0.0, 0.8], 0.0).unwrap(),
        &plain_source,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let contact = IntegratedCaloricLink::new_with_contact(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        &plain_source,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let e = [energy(0, 650.0), energy(1, 350.0)];
    let original = contact.evaluate(e, || Ok(())).unwrap();
    let plain = pure.evaluate(e, || Ok(())).unwrap();
    assert_eq!(plain.power_a_to_b_w, original.power_a_to_b_w);
    assert_eq!(
        plain.gap_resistance_derivative_w2_m2_k,
        original.contact_resistance_derivative_w2_m2_k
    );
    let state = |link: &dyn ImplicitLinkLaw| {
        link.state_at_specific_enthalpies(
            e,
            pure.native_geometry(),
            thermal.contact_patch(LINK).unwrap(),
            &mut || Ok(()),
        )
        .unwrap()
    };
    assert_eq!(state(&pure), state(&contact));
}

#[test]
fn g1_native_source_gap_evolves_both_calorics_with_one_heat_debit_and_transactional_retry() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let (t0, masses, d, dt) = ([650.0, 330.0], [0.4, 0.2], [0.01, 0.005], 300.0);
    let thermal = timeline([&cv, &cp], t0, masses, AREA);
    let law = source(&cv, &cp, [0.002, -0.001]);
    let link = IntegratedCaloricGrayLink::new_with_gap_resistance(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar(d, AREA).unwrap(),
        GrayLayerGap::equal_area_full_view(EMISSIVITIES, 0.0).unwrap(),
        &law,
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
    let e0 = [energy(0, t0[0]), energy(1, t0[1])];
    let residual = |q: f64| {
        let t = [
            temperature(0, e0[0] - q / masses[0]),
            temperature(1, e0[1] + q / masses[1]),
        ];
        q - dt * power(t, d, law.resistance(t))
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
    let mut calls = 0;
    let (continued, receipt) = chain
        .advance_to(&next, 2.0 * dt, controls, &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    assert!(calls > 8);
    let mut refused = 0;
    assert!(matches!(
        chain.advance_to(&next, 2.0 * dt, controls, &mut workspace, || {
            refused += 1;
            if refused == calls - 1 {
                Err(DepositionError::Cancelled)
            } else {
                Ok(())
            }
        }),
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(next.time_s(), dt);
    assert_eq!(next.cells()[&CellId(0)].specific_h_j_kg, e[0]);
    assert_eq!(next.cells()[&CellId(1)].specific_h_j_kg, e[1]);
    let (retry, repeated) = chain
        .advance_to(&next, 2.0 * dt, controls, &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), continued.cells());
    assert_eq!(repeated, receipt);
}

#[test]
fn g4_gap_source_admission_preserves_original_authority_and_initial_disjoint_support() {
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
    let gap = GrayLayerGap::equal_area_full_view(EMISSIVITIES, 0.0).unwrap();
    let law = source(&cv, &cp, [0.002, -0.001]);
    let admit = |source: &dyn ImplicitContactLaw, gap| {
        IntegratedCaloricGrayLink::new_with_gap_resistance(
            &thermal,
            LINK,
            [&layers[0], &layers[1]],
            geometry,
            gap,
            source,
            layer_controls(),
            || Ok(()),
        )
        .map(|_| ())
    };
    let cloned_cv = cv.clone();
    let mut foreign = law.clone();
    foreign.laws[0] = &cloned_cv;
    let mut bad_bounds = law.clone();
    bad_bounds.bounds[0][0] = -1.0;
    let mut missing = law.clone();
    missing.maximum_difference_k = 200.0;
    for source in [&foreign, &bad_bounds, &missing] {
        assert!(admit(source, gap).is_err());
    }
    for r in [0.0, -0.01, f64::INFINITY, f64::from_bits(1)] {
        let mut invalid = law.clone();
        invalid.reference_r = r;
        assert!(admit(&invalid, gap).is_err());
    }
    assert!(
        admit(
            &law,
            GrayLayerGap::equal_area_full_view(EMISSIVITIES, 1.0).unwrap()
        )
        .is_err()
    );
    let mut restricted = law.clone();
    restricted.bounds[0][0] = energy(0, 610.0);
    let limited = IntegratedCaloricGrayLink::new_with_gap_resistance(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        gap,
        &restricted,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let queries = restricted.queries.get();
    assert!(
        limited
            .evaluate([energy(0, 600.0), energy(1, 350.0)], || Ok(()))
            .is_err()
    );
    assert_eq!(restricted.queries.get(), queries);
    assert!(matches!(
        limited.evaluate([energy(0, 650.0), energy(1, 350.0)], || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));

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
        event: Event::Contact(*thermal.links().get(&LINK).unwrap()),
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
    let mut original = source(&narrow_cv, &narrow_cp, [0.002, -0.001]);
    original.reference_r = 0.02;
    let admitted = IntegratedCaloricGrayLink::new_with_gap_resistance(
        &narrow,
        LINK,
        [&narrow_layers[0], &narrow_layers[1]],
        geometry,
        gap,
        &original,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(original.queries.get(), 1);
    let response = admitted.evaluate([0.0; 2], || Ok(())).unwrap();
    near(
        response.power_a_to_b_w,
        power([650.0, 350.0], [0.01, 0.008], 0.02),
        3e-9,
    );
    assert!(
        response.interface_temperature_k[0] >= 600.0
            && response.interface_temperature_k[1] <= 400.0
    );
    assert_eq!(narrow.time_s(), 0.0);
    assert_eq!(narrow.cells()[&CellId(0)].specific_h_j_kg, 0.0);
}
