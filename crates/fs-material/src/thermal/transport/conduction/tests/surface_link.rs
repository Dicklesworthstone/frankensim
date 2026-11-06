//! G1/G3/G4: live surface exchange on original mixed Cp/Cv inventories.
use super::*;
use fs_thermal::enthalpy::conduction::implicit::link::{
    ImplicitSurfaceExchangeLaw, ImplicitSurfaceExchangeState,
};
use std::cell::{Cell, RefCell};

const SIGMA: f64 = 5.670_374_419e-8;

#[derive(Clone, Copy)]
enum Model {
    Gray,
    Onset { origin: [f64; 2], slopes: [f64; 2] },
}

struct Surface<'a> {
    laws: [&'a dyn ThermalLaw; 2],
    bounds: [[f64; 2]; 2],
    face_bounds: [[f64; 2]; 2],
    model: Model,
    refuse: Cell<bool>,
    queries: RefCell<Vec<([f64; 2], [f64; 2])>>,
}

fn source<'a>(cv: &'a IsochoricHeatCapacityLaw, cp: &'a IsobaricHeatCapacityLaw) -> Surface<'a> {
    Surface {
        laws: [cv, cp],
        bounds: [
            cv.specific_internal_energy_span_j_kg(),
            cp.specific_enthalpy_span_j_kg(),
        ],
        face_bounds: [[300.0, 700.0]; 2],
        model: Model::Gray,
        refuse: Cell::new(false),
        queries: RefCell::new(Vec::new()),
    }
}

// Authored source coefficients. Each emissivity depends on its actual face
// temperature and the OTHER bulk energy; the gap HTC depends on all four.
fn coefficients(e: [f64; 2], s: [f64; 2]) -> (f64, [f64; 2]) {
    let x = [e[0] - energy(0, 500.0), e[1] - energy(1, 500.0)];
    let h = 35.0 + 0.01 * (s[0] - 500.0) + 0.02 * (s[1] - 500.0) + 2e-5 * x[0] + 1e-5 * x[1];
    let eps = [
        0.7 + 0.0002 * (s[0] - 500.0) + 1e-7 * x[1],
        0.55 - 0.0001 * (s[1] - 500.0) + 2e-7 * x[0],
    ];
    (h, eps)
}

impl ImplicitSurfaceExchangeLaw for Surface<'_> {
    fn caloric_laws(&self) -> [&dyn ThermalLaw; 2] {
        self.laws
    }
    fn energy_bounds_j_kg(&self) -> [[f64; 2]; 2] {
        self.bounds
    }
    fn surface_temperature_bounds_k(&self) -> [[f64; 2]; 2] {
        self.face_bounds
    }
    fn state_at_surfaces(
        &self,
        e: [f64; 2],
        s: [f64; 2],
    ) -> Result<ImplicitSurfaceExchangeState, DepositionError> {
        self.queries.borrow_mut().push((e, s));
        if self.refuse.get()
            || (0..2).any(|i| {
                !e[i].is_finite()
                    || e[i] < self.bounds[i][0]
                    || e[i] > self.bounds[i][1]
                    || !s[i].is_finite()
                    || s[i] < self.face_bounds[i][0]
                    || s[i] > self.face_bounds[i][1]
            })
        {
            return Err(DepositionError::MaterialDomain);
        }
        if let Model::Onset { origin, slopes } = self.model {
            let g = slopes[0] * (e[0] - origin[0]) + slopes[1] * (e[1] - origin[1]);
            return ImplicitSurfaceExchangeState::from_conductance(s, g, [0.0; 2], slopes);
        }
        let (h, eps) = coefficients(e, s);
        let effective = 1.0 / (1.0 / eps[0] + 1.0 / eps[1] - 1.0);
        let epsilon_slopes = eps.map(|value| (effective / value).powi(2));
        let sum = s[0] + s[1];
        let square_sum = s[0] * s[0] + s[1] * s[1];
        let radiation_secant = SIGMA * sum * square_sum;
        let dg_ds = [
            0.01 + epsilon_slopes[0] * 0.0002 * radiation_secant
                + effective * SIGMA * (square_sum + 2.0 * s[0] * sum),
            0.02 - epsilon_slopes[1] * 0.0001 * radiation_secant
                + effective * SIGMA * (square_sum + 2.0 * s[1] * sum),
        ];
        let dg_de = [
            2e-5 + epsilon_slopes[1] * 2e-7 * radiation_secant,
            1e-5 + epsilon_slopes[0] * 1e-7 * radiation_secant,
        ];
        ImplicitSurfaceExchangeState::from_conductance(
            s,
            h + effective * radiation_secant,
            dg_ds,
            dg_de,
        )
    }
}

// Independent fourth-power heat law and quadratic Kirchhoff primitives,
// bisected directly in flux density. No production surface/root/inverse helper.
fn reference(e: [f64; 2], d: [f64; 2], model: Model) -> (f64, [f64; 2]) {
    let t = [temperature(0, e[0]), temperature(1, e[1])];
    let sign = if t[0] >= t[1] { 1.0 } else { -1.0 };
    let faces = |j: f64| {
        [
            (t[0] * t[0] - 2.0 * sign * j * d[0] / 0.04).sqrt(),
            (t[1] * t[1] + 2.0 * sign * j * d[1] / 0.02).sqrt(),
        ]
    };
    let flux = |s: [f64; 2]| match model {
        Model::Gray => {
            let (h, eps) = coefficients(e, s);
            let effective = 1.0 / (1.0 / eps[0] + 1.0 / eps[1] - 1.0);
            h * (s[0] - s[1]) + effective * SIGMA * (s[0].powi(4) - s[1].powi(4))
        }
        Model::Onset { origin, slopes } => {
            (slopes[0] * (e[0] - origin[0]) + slopes[1] * (e[1] - origin[1])) * (s[0] - s[1])
        }
    };
    let (mut lo, mut hi) = (
        0.0,
        (t[0] * t[0] - t[1] * t[1]).abs() / (2.0 * (d[0] / 0.04 + d[1] / 0.02)),
    );
    for _ in 0..80 {
        let j = (lo + hi) / 2.0;
        if sign * flux(faces(j)) > j {
            lo = j;
        } else {
            hi = j;
        }
    }
    let j = (lo + hi) / 2.0;
    (sign * AREA * j, faces(j))
}

#[test]
fn g1_live_gray_surface_exchange_retains_crossed_source_and_original_caloric_tangents() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let d = [0.01, 0.008];
    let thermal = timeline([&cv, &cp], [650.0, 350.0], [0.4, 0.32], AREA);
    let law = source(&cv, &cp);
    let link = IntegratedCaloricSurfaceLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar(d, AREA).unwrap(),
        &law,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    for side in 0..2 {
        assert!(std::ptr::addr_eq(
            link.caloric_laws()[side] as *const dyn ThermalLaw,
            law.laws[side] as *const dyn ThermalLaw,
        ));
    }
    for t in [[650.0, 350.0], [350.0, 650.0], [475.0, 475.0]] {
        let e = [energy(0, t[0]), energy(1, t[1])];
        law.queries.borrow_mut().clear();
        let response = link.evaluate(e, || Ok(())).unwrap();
        let (expected, surfaces) = reference(e, d, law.model);
        near(response.power_a_to_b_w, expected, 4e-9);
        for side in 0..2 {
            near(response.interface_temperature_k[side], surfaces[side], 2e-9);
        }
        assert!(!law.queries.borrow().is_empty());
        for (queried_energy, queried_surfaces) in law.queries.borrow().iter() {
            assert_eq!(*queried_energy, e);
            for side in 0..2 {
                assert!((300.0..=700.0).contains(&queried_surfaces[side]));
            }
        }
        if t[0] != t[1] {
            assert!((surfaces[0] - t[0]).abs() > 0.1);
            assert!((surfaces[1] - t[1]).abs() > 0.1);
        }
        let native = link
            .state_at_specific_enthalpies(
                e,
                link.native_geometry(),
                thermal.contact_patch(LINK).unwrap(),
                &mut || Ok(()),
            )
            .unwrap();
        near(native.conductance_w_k * (t[0] - t[1]), expected, 4e-9);
        near(response.surface_heat_residual_w, 0.0, 4e-9);
        near(response.layer_heat_residual_w, 0.0, 4e-9);
        let ds = response.source.surface_derivatives_w_m2_k;
        let denominator =
            1.0 + ds[0] * d[0] / (0.04 * surfaces[0]) - ds[1] * d[1] / (0.02 * surfaces[1]);
        near(
            response.source_flux_derivative_scale_m2,
            AREA / denominator,
            2e-12,
        );
        let capacities = [200.0 + 0.5 * t[0], 400.0 + t[1]];
        for side in 0..2 {
            let frozen_energy_derivative = AREA * ds[side] * t[side] / surfaces[side] / denominator;
            near(
                response.endpoint_derivative_w_k[side],
                frozen_energy_derivative,
                2e-9,
            );
            let expected_derivative = frozen_energy_derivative / capacities[side]
                + AREA / denominator * response.source.energy_derivatives_w_kg_m2_j[side];
            near(
                native.power_enthalpy_derivatives_w_kg_j[side],
                expected_derivative,
                3e-12,
            );
            let mut lower = e;
            let mut upper = e;
            lower[side] -= 0.25;
            upper[side] += 0.25;
            near(
                native.power_enthalpy_derivatives_w_kg_j[side],
                (reference(upper, d, law.model).0 - reference(lower, d, law.model).0) / 0.5,
                3e-8,
            );
        }
    }
}

#[test]
fn g1_live_faces_drive_one_native_cp_cv_heat_step_with_entropy_and_transactional_retry() {
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
    let law = source(&cv, &cp);
    let link = IntegratedCaloricSurfaceLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar(d, AREA).unwrap(),
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
    let endpoint = |q: f64| [e0[0] - q / masses[0], e0[1] + q / masses[1]];
    let residual = |q: f64| q - dt * reference(endpoint(q), d, law.model).0;
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
    for side in 0..2 {
        near(e[side], endpoint(q)[side], 5e-6);
    }
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
    assert_eq!(report.energy.external_heat_j, 0.0);
    assert_eq!(report.energy.substeps, 1);
    assert_eq!(
        next.material_storage_basis(CV_MATERIAL).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(
        next.material_storage_basis(CP_MATERIAL).unwrap(),
        ThermalStorageBasis::Enthalpy
    );
    near(
        dt * link.evaluate(e, || Ok(())).unwrap().power_a_to_b_w,
        q,
        1e-6,
    );
    let mut calls = 0;
    let (continued, receipt) = chain
        .advance_to(&next, 2.0 * dt, controls, &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap();
    let mut cancelled = 0;
    let refusal = chain.advance_to(&next, 2.0 * dt, controls, &mut workspace, || {
        cancelled += 1;
        if cancelled == calls - 1 {
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
fn g4_supported_faces_need_neither_bulk_temperature_and_keep_original_identity_and_geometry() {
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
    let e = [energy(0, 650.0), energy(1, 350.0)];
    let expected = reference(e, d, Model::Gray);
    let mut law = source(&cv, &cp);
    law.face_bounds = expected.1.map(|face| [face - 1.0, face + 1.0]);
    assert!(law.face_bounds[0][0] > law.face_bounds[1][1]);
    assert!(650.0 > law.face_bounds[0][1] && 350.0 < law.face_bounds[1][0]);
    let link = IntegratedCaloricSurfaceLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        geometry,
        &law,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    near(
        link.evaluate(e, || Ok(())).unwrap().power_a_to_b_w,
        expected.0,
        4e-9,
    );
    for (_, surfaces) in law.queries.borrow().iter() {
        for side in 0..2 {
            assert!(
                (law.face_bounds[side][0]..=law.face_bounds[side][1]).contains(&surfaces[side])
            );
        }
    }
    let patch = thermal.contact_patch(LINK).unwrap();
    let mut wrong_geometry = link.native_geometry();
    wrong_geometry.contact_resistance_m2_k_w = 0.001;
    assert!(
        link.state_at_specific_enthalpies(e, wrong_geometry, patch, &mut || Ok(()))
            .is_err()
    );
    assert!(link.evaluate([-1.0, e[1]], || Ok(())).is_err());
    law.refuse.set(true);
    assert!(matches!(
        link.evaluate(e, || Ok(())),
        Err(DepositionError::MaterialDomain)
    ));
    law.refuse.set(false);
    near(
        link.evaluate(e, || Ok(())).unwrap().power_a_to_b_w,
        expected.0,
        4e-9,
    );
    let orphan_cv = fa.cv();
    let wrong_identity = source(&orphan_cv, &cp);
    assert!(
        IntegratedCaloricSurfaceLink::new(
            &thermal,
            LINK,
            [&layers[0], &layers[1]],
            geometry,
            &wrong_identity,
            layer_controls(),
            || Ok(()),
        )
        .is_err()
    );
    let mut empty_support = source(&cv, &cp);
    empty_support.face_bounds[0] = [500.0, 500.0];
    assert!(
        IntegratedCaloricSurfaceLink::new(
            &thermal,
            LINK,
            [&layers[0], &layers[1]],
            geometry,
            &empty_support,
            layer_controls(),
            || Ok(()),
        )
        .is_err()
    );
}

#[test]
fn g1_insulating_onset_keeps_both_inward_energy_derivatives_including_positive_cold_feedback() {
    let (fa, fb) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (fa.cv(), fb.cp());
    let a = IsochoricConductionLaw::try_new(&fa.card, &cv, &fa.conductivity, SELECTION).unwrap();
    let b = IsobaricConductionLaw::try_new(&fb.card, &cp, &fb.conductivity, SELECTION).unwrap();
    let layers = [
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    ];
    let e = [energy(0, 650.0), energy(1, 350.0)];
    let d = [0.01, 0.008];
    let thermal = timeline([&cv, &cp], [650.0, 350.0], [0.4, 0.32], AREA);
    let mut law = source(&cv, &cp);
    let slopes = [2e-4, 3e-4];
    law.model = Model::Onset { origin: e, slopes };
    for side in 0..2 {
        law.bounds[side][0] = e[side];
    }
    let link = IntegratedCaloricSurfaceLink::new(
        &thermal,
        LINK,
        [&layers[0], &layers[1]],
        LayerGeometry::planar(d, AREA).unwrap(),
        &law,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let response = link.evaluate(e, || Ok(())).unwrap();
    let native = link
        .state_at_specific_enthalpies(
            e,
            link.native_geometry(),
            thermal.contact_patch(LINK).unwrap(),
            &mut || Ok(()),
        )
        .unwrap();
    assert_eq!(response.power_a_to_b_w, 0.0);
    assert_eq!(native.conductance_w_k, 0.0);
    for side in 0..2 {
        let expected = AREA * 300.0 * slopes[side];
        near(
            native.power_enthalpy_derivatives_w_kg_j[side],
            expected,
            1e-12,
        );
        let mut next = e;
        next[side] += 0.1;
        near(
            link.evaluate(next, || Ok(())).unwrap().power_a_to_b_w / 0.1,
            expected,
            3e-8,
        );
    }
    assert!(native.power_enthalpy_derivatives_w_kg_j[1] > 0.0);
}
