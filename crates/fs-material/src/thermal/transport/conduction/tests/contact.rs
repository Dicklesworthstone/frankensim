//! G1/G4: original state-dependent contact resistance on integrated Cp/Cv faces.

use super::*;
use fs_thermal::enthalpy::conduction::implicit::contact::{
    ImplicitContactLaw, ImplicitContactState,
};

#[derive(Clone, Copy)]
struct ExponentialContact<'a> {
    laws: [&'a dyn ThermalLaw; 2],
    bounds: [[f64; 2]; 2],
    reference_resistance: f64,
    slopes_per_k: [f64; 2],
    maximum_difference_k: f64,
}
impl ExponentialContact<'_> {
    fn resistance(&self, t: [f64; 2]) -> f64 {
        self.reference_resistance
            * (self.slopes_per_k[0] * (t[0] - 500.0) + self.slopes_per_k[1] * (t[1] - 500.0)).exp()
    }
}
impl ImplicitContactLaw for ExponentialContact<'_> {
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
        if e.iter()
            .zip(self.bounds)
            .any(|(&v, [lo, hi])| !v.is_finite() || v < lo || v > hi)
        {
            return Err(DepositionError::MaterialDomain);
        }
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

fn contact<'a>(
    cv: &'a IsochoricHeatCapacityLaw,
    cp: &'a IsobaricHeatCapacityLaw,
    slopes_per_k: [f64; 2],
) -> ExponentialContact<'a> {
    ExponentialContact {
        laws: [cv, cp],
        bounds: [
            cv.specific_internal_energy_span_j_kg(),
            cp.specific_enthalpy_span_j_kg(),
        ],
        reference_resistance: 0.004,
        slopes_per_k,
        maximum_difference_k: f64::INFINITY,
    }
}

// Independent scalar root in flux density. Quadratic Kirchhoff primitives give
// both face temperatures directly; the production layer solve is not invoked.
fn exact_power(t: [f64; 2], distances: [f64; 2], resistance: f64) -> f64 {
    let sign = if t[0] >= t[1] { 1.0 } else { -1.0 };
    let (mut lo, mut hi) = (
        0.0,
        (t[0] * t[0] - t[1] * t[1]).abs() / (2.0 * (distances[0] / 0.04 + distances[1] / 0.02)),
    );
    for _ in 0..80 {
        let j = (lo + hi) / 2.0;
        let a = (t[0] * t[0] - 2.0 * sign * j * distances[0] / 0.04).sqrt();
        let b = (t[1] * t[1] + 2.0 * sign * j * distances[1] / 0.02).sqrt();
        if sign * (a - b) - j * resistance > 0.0 {
            lo = j;
        } else {
            hi = j;
        }
    }
    sign * AREA * (lo + hi) / 2.0
}

#[test]
fn g1_source_contact_keeps_full_two_sided_caloric_tangents_and_negative_differential_flux() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let (t, distances) = ([650.0, 350.0], [0.01, 0.008]);
    let thermal = timeline([&cv, &cp], t, [0.4, 0.32], AREA);
    let law = contact(&cv, &cp, [0.02, 0.005]);
    let link = IntegratedCaloricLink::new_with_contact(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar(distances, AREA).unwrap(),
        &law,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    assert_eq!(link.native_geometry().contact_resistance_m2_k_w, 0.0);
    let e = [energy(0, t[0]), energy(1, t[1])];
    let state = link
        .state_at_specific_enthalpies(
            e,
            link.native_geometry(),
            thermal.contact_patch(LINK).unwrap(),
            &mut || Ok(()),
        )
        .unwrap();
    let response = link.evaluate(e, || Ok(())).unwrap();
    near(
        response.power_a_to_b_w,
        exact_power(t, distances, law.resistance(t)),
        3e-9,
    );
    assert!(response.power_a_to_b_w > 0.0 && state.conductance_w_k > 0.0);
    // Positive resistance is passive; its changing value need not produce a
    // positive differential hot-side flux. The original fixed-R partial still is.
    assert!(response.endpoint_derivative_w_k[0] > 0.0);
    assert!(state.power_enthalpy_derivatives_w_kg_j[0] < 0.0);
    for side in 0..2 {
        let mut lo = e;
        lo[side] -= 0.25;
        let mut hi = e;
        hi[side] += 0.25;
        let power = |e: [f64; 2]| {
            let t = [temperature(0, e[0]), temperature(1, e[1])];
            exact_power(t, distances, law.resistance(t))
        };
        near(
            state.power_enthalpy_derivatives_w_kg_j[side],
            (power(hi) - power(lo)) / 0.5,
            2e-8,
        );
        near(
            state.power_enthalpy_derivatives_w_kg_j[side],
            (link.evaluate(hi, || Ok(())).unwrap().power_a_to_b_w
                - link.evaluate(lo, || Ok(())).unwrap().power_a_to_b_w)
                / 0.5,
            2e-8,
        );
    }
}

#[test]
fn g1_native_contact_feedback_closes_mixed_caloric_heat_and_entropy_with_transactional_retry() {
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
    let law = contact(&cv, &cp, [0.002, -0.001]);
    let link = IntegratedCaloricLink::new_with_contact(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar(distances, AREA).unwrap(),
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
        q - dt * exact_power(t, distances, law.resistance(t))
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
    // A late refusal starts from an already accepted positive physical time.
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
    let (retry, repeated) = chain
        .advance_to(&next, 2.0 * dt, controls, &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), continued.cells());
    assert_eq!(repeated, receipt);
}

#[test]
fn g4_source_contact_admission_and_disjoint_material_support_keep_original_authority() {
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
    let law = contact(&cv, &cp, [0.002, -0.001]);
    let cloned_cv = cv.clone();
    let mut foreign = law;
    foreign.laws[0] = &cloned_cv;
    let mut bad_bounds = law;
    bad_bounds.bounds[0][0] = -1.0;
    let mut missing = law;
    missing.maximum_difference_k = 200.0;
    let mut negative = law;
    negative.reference_resistance = -0.01;
    for source in [&foreign, &bad_bounds, &missing, &negative] {
        assert!(
            IntegratedCaloricLink::new_with_contact(
                &thermal,
                LINK,
                [&layers[0], &layers[1]],
                geometry,
                source,
                layer_controls(),
                || Ok(())
            )
            .is_err()
        );
    }
    let mut restricted = law;
    restricted.bounds[0][0] = energy(0, 610.0);
    let limited = IntegratedCaloricLink::new_with_contact(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        &restricted,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    assert!(
        limited
            .evaluate([energy(0, 600.0), energy(1, 350.0)], || Ok(()))
            .is_err()
    );
    assert!(matches!(
        limited.evaluate([energy(0, 650.0), energy(1, 350.0)], || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        thermal.cells()[&CellId(0)].specific_h_j_kg,
        energy(0, 650.0)
    );

    // A real finite contact jump can connect disjoint original k/Cp/Cv spans.
    // Probing a fictitious perfect contact during construction would refuse it.
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
    let a =
        IsochoricConductionLaw::try_new(&fa.card, &narrow_cv, &fa.conductivity, SELECTION).unwrap();
    let b =
        IsobaricConductionLaw::try_new(&fb.card, &narrow_cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let mut events: Vec<_> = thermal
        .cells()
        .values()
        .cloned()
        .map(|mut cell| {
            cell.specific_h_j_kg = 0.0; // each original source's declared initial datum
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
    let mut source = contact(&narrow_cv, &narrow_cp, [0.0, 0.0]);
    source.reference_resistance = 0.01;
    assert!(
        IntegratedCaloricLink::new(
            &narrow,
            LINK,
            [&layers[0], &layers[1]],
            geometry,
            0.0,
            layer_controls(),
            || Ok(())
        )
        .is_err()
    );
    let admitted = IntegratedCaloricLink::new_with_contact(
        &narrow,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        &source,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let response = admitted.evaluate([0.0; 2], || Ok(())).unwrap();
    near(
        response.power_a_to_b_w,
        exact_power([650.0, 350.0], [0.01, 0.008], 0.01),
        3e-9,
    );
    assert!(response.interface_temperature_k[0] >= 600.0);
    assert!(response.interface_temperature_k[1] <= 400.0);
}
