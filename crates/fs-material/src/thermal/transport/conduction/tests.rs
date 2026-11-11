//! G1/G3/G4: original Cp/Cv storage and independently integrated face fluxes.
//! The material cards are synthetic numerical fixtures, not measured datasets.

#[path = "tests/boundary.rs"]
mod boundary;
#[path = "tests/contact.rs"]
mod contact;
#[path = "tests/exchange.rs"]
mod exchange;
#[path = "tests/gray_boundary.rs"]
mod gray_boundary;
#[path = "tests/gray_link.rs"]
mod gray_link;
#[path = "tests/surface_link.rs"]
mod surface_link;

#[path = "tests/current_geometry.rs"]
mod current_geometry;

#[path = "tests/radial_geometry.rs"]
mod radial_geometry;

use super::*;
use crate::thermal::transport::{
    CONDUCTIVITY, DENSITY, SPECIFIC_HEAT,
    isobaric::{
        IsobaricHeatCapacityLaw, IsobaricHeatCapacityQuery, conduction::IsobaricConductionLaw,
    },
    isochoric::{
        IsochoricHeatCapacityLaw, IsochoricHeatCapacityQuery, conduction::IsochoricConductionLaw,
    },
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSelection, ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_qty::{
    Pressure, QuantitySpec, Temperature,
    semantic::{QuantityKind, SemanticType, ValueForm},
};
use fs_thermal::enthalpy::{
    Boundary, CellGeometry, CellId, DepositedCell, DepositionBudget, DepositionSimulation, Event,
    MaterialId, StepConfig, ThermalLink, ThermalStorageBasis, TimedEvent,
    conduction::{
        ConductionObserver, ConductivityLaw,
        implicit::{
            ImplicitConductionChain, ImplicitConductionControls, ImplicitConductionError,
            ImplicitConductionLaw, ImplicitConductionWorkspace, link::ImplicitLinkLaw,
        },
    },
};

const RHO: f64 = 1000.0;
const AREA: f64 = 0.02;
const LINK: u64 = 7;
const CV_MATERIAL: MaterialId = MaterialId(1);
const CP_MATERIAL: MaterialId = MaterialId(2);
const SELECTION: ClaimSelection = ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly);

struct Fixture {
    card: MaterialCard,
    capacity: PropertyKey,
    conductivity: PropertyKey,
    lower: QueryPoint,
    upper: QueryPoint,
}

impl Fixture {
    fn new(isochoric: bool) -> Self {
        let temperature = QuantitySpec::semantic(SemanticType::new(
            QuantityKind::AbsoluteTemperature,
            ValueForm::Static,
        ));
        let (axis, quantity, fixed) = if isochoric {
            ("rho", QuantitySpec::dimensional(DENSITY), RHO)
        } else {
            ("p", QuantitySpec::dimensional(Pressure::DIMS), 100000.0)
        };
        let capacity = PropertyKey::new(
            if isochoric { "source_cv" } else { "source_cp" },
            SPECIFIC_HEAT,
        );
        let conductivity = PropertyKey::new("normal_conductivity", CONDUCTIVITY);
        let validity = ValidityDomain::unconstrained()
            .with_quantity("T", temperature, 300.0, 700.0)
            .with_quantity(axis, quantity, fixed, fixed);
        let capacity_temperatures = if isochoric {
            [300.0, 440.0, 700.0]
        } else {
            [300.0, 475.0, 700.0]
        };
        let conductivity_temperatures = if isochoric {
            [300.0, 450.0, 600.0, 700.0]
        } else {
            [300.0, 350.0, 550.0, 700.0]
        };
        let capacity_knots = capacity_temperatures
            .into_iter()
            .map(|t| {
                (
                    t,
                    if isochoric {
                        200.0 + 0.5 * t
                    } else {
                        400.0 + t
                    },
                )
            })
            .collect();
        let conductivity_knots = conductivity_temperatures
            .into_iter()
            .map(|t| (t, if isochoric { 0.04 * t } else { 0.02 * t }))
            .collect();
        let mut claims = ClaimSet::new();
        for (key, dims, knots) in [
            (capacity.clone(), SPECIFIC_HEAT, capacity_knots),
            (conductivity.clone(), CONDUCTIVITY, conductivity_knots),
        ] {
            claims
                .insert_claim(PropertyClaim {
                    key,
                    value: PropertyValue::Curve {
                        abscissa: "T".into(),
                        abscissa_dims: Temperature::DIMS,
                        knots,
                        dims,
                    },
                    validity: validity.clone(),
                    interpolation: InterpolationPolicy::LinearInside,
                    uncertainty: UncertaintyModel::Unstated,
                    observations: vec![],
                    provenance: Provenance {
                        source: "synthetic independent caloric and conductivity paths".into(),
                        license: "test-only".into(),
                        artifact: None,
                    },
                })
                .unwrap();
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: if isochoric {
                    "Cv fixture".into()
                } else {
                    "Cp fixture".into()
                },
                phase: "single phase".into(),
                process: "integrated conduction numerical fixture".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let point = |t| {
            QueryPoint::new()
                .with_quantity("T", temperature, t)
                .unwrap()
                .with_quantity(axis, quantity, fixed)
                .unwrap()
        };
        Self {
            card,
            capacity,
            conductivity,
            lower: point(300.0),
            upper: point(700.0),
        }
    }

    fn cv(&self) -> IsochoricHeatCapacityLaw {
        IsochoricHeatCapacityLaw::try_new(
            &self.card,
            IsochoricHeatCapacityQuery {
                specific_heat_cv: &self.capacity,
                temperature_axis: "T",
                density_axis: "rho",
                lower: &self.lower,
                upper: &self.upper,
                selection: SELECTION,
            },
            300.0,
        )
        .unwrap()
    }

    fn cp(&self) -> IsobaricHeatCapacityLaw {
        IsobaricHeatCapacityLaw::try_new(
            &self.card,
            IsobaricHeatCapacityQuery {
                specific_heat_cp: &self.capacity,
                temperature_axis: "T",
                pressure_axis: "p",
                lower: &self.lower,
                upper: &self.upper,
                selection: SELECTION,
            },
            300.0,
        )
        .unwrap()
    }
}

fn layer_controls() -> LayerConductionConfig {
    LayerConductionConfig {
        max_iterations: 64,
        absolute_heat_tolerance_w: 1e-10,
        relative_heat_tolerance: 1e-12,
    }
}

fn near(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.17e}, expected={expected:.17e}, tolerance={tolerance:.3e}"
    );
}

// Independent quadratic primitives and their cancellation-resistant inverses.
fn energy(side: usize, t: f64) -> f64 {
    let x = t - 300.0;
    if side == 0 {
        350.0 * x + 0.25 * x * x
    } else {
        700.0 * x + 0.5 * x * x
    }
}

fn temperature(side: usize, e: f64) -> f64 {
    if side == 0 {
        300.0 + 2.0 * e / (350.0 + (350.0 * 350.0 + e).sqrt())
    } else {
        300.0 + 2.0 * e / (700.0 + (700.0 * 700.0 + 2.0 * e).sqrt())
    }
}

fn timeline<'a>(
    laws: [&'a dyn ThermalLaw; 2],
    temperatures: [f64; 2],
    masses: [f64; 2],
    area: f64,
) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    let mut start = 0.0;
    for side in 0..2 {
        let length = masses[side] / (RHO * area);
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CellId(side as u64),
                segment_id: side as u64,
                body_id: 0,
                material_id: [CV_MATERIAL, CP_MATERIAL][side],
                mass_kg: masses[side],
                specific_h_j_kg: energy(side, temperatures[side]),
                geometry: CellGeometry {
                    start_m: [start, 0.0, 0.0],
                    end_m: [start + length, 0.0, 0.0],
                    volume_m3: masses[side] / RHO,
                    surface_area_m2: 2.0 * area + 4.0 * area.sqrt() * length,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        });
        start += length;
    }
    events.push(TimedEvent {
        time_s: 0.0,
        event: Event::Contact(ThermalLink {
            id: LINK,
            a: CellId(0),
            b: CellId(1),
            area_m2: area,
            conductance_w_k: 0.0,
        }),
    });
    DepositionSimulation::new(
        events,
        &[(CV_MATERIAL, laws[0]), (CP_MATERIAL, laws[1])],
        DepositionBudget {
            max_cells: 2,
            max_links: 1,
            max_events: 3,
            max_substeps: 128,
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
    .0
}

#[test]
fn g1_mixed_caloric_layers_preserve_native_sources_and_exact_contact_flux_tangents() {
    let (source_a, source_b) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (source_a.cv(), source_b.cp());
    let a = IsochoricConductionLaw::try_new(&source_a.card, &cv, &source_a.conductivity, SELECTION)
        .unwrap();
    let b = IsobaricConductionLaw::try_new(&source_b.card, &cp, &source_b.conductivity, SELECTION)
        .unwrap();
    let (layer_a, layer_b) = (
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    );
    assert!(core::ptr::addr_eq(
        layer_a.caloric_law(),
        &cv as &dyn ThermalLaw
    ));
    assert!(core::ptr::addr_eq(
        layer_b.caloric_law(),
        &cp as &dyn ThermalLaw
    ));
    assert!(core::ptr::eq(layer_a.source_line(), a.source_line()));
    assert!(core::ptr::eq(layer_b.source_line(), b.source_line()));
    assert_eq!(
        layer_a.caloric_law().storage_basis(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(
        layer_b.caloric_law().storage_basis(),
        ThermalStorageBasis::Enthalpy
    );
    assert_eq!(
        layer_a.conductivity_curve().knots(),
        &[[300.0, 12.0], [450.0, 18.0], [600.0, 24.0], [700.0, 28.0]]
    );
    assert_eq!(
        layer_b.conductivity_curve().knots(),
        &[[300.0, 6.0], [350.0, 7.0], [550.0, 11.0], [700.0, 14.0]]
    );

    // Manufacture an exact solution from the interface temperatures and flux.
    // k_a=.04T, k_b=.02T, so each integral is a quadratic in temperature.
    let (interface_a, interface_b, flux, resistance) = (450.0_f64, 420.0_f64, 30000.0_f64, 0.001);
    let distances = [0.01, 0.008];
    let ta = (interface_a * interface_a + 2.0 * flux * distances[0] / 0.04).sqrt();
    let tb = (interface_b * interface_b - 2.0 * flux * distances[1] / 0.02).sqrt();
    let thermal = timeline([&cv, &cp], [ta, tb], [0.4, 0.32], AREA);
    let link = IntegratedCaloricLink::new(
        &thermal,
        LINK,
        [&layer_a, &layer_b],
        LayerGeometry::planar(distances, AREA).unwrap(),
        resistance,
        layer_controls(),
        || Ok(()),
    )
    .unwrap();
    let energies = [energy(0, ta), energy(1, tb)];
    let response = link.evaluate(energies, || Ok(())).unwrap();
    near(response.power_a_to_b_w, AREA * flux, 2e-9);
    near(response.interface_temperature_k[0], interface_a, 2e-10);
    near(response.interface_temperature_k[1], interface_b, 2e-10);
    near(response.contact_temperature_residual_k, 0.0, 2e-10);
    let denominator =
        resistance + distances[0] / (0.04 * interface_a) + distances[1] / (0.02 * interface_b);
    let expected_tangent = [
        AREA * (ta / interface_a) / denominator,
        -AREA * (tb / interface_b) / denominator,
    ];
    let native = link
        .state_at_specific_enthalpies(
            energies,
            link.native_geometry(),
            thermal.contact_patch(LINK).unwrap(),
            &mut || Ok(()),
        )
        .unwrap();
    for side in 0..2 {
        near(
            response.endpoint_derivative_w_k[side],
            expected_tangent[side],
            2e-10,
        );
        let capacity = if side == 0 {
            200.0 + 0.5 * ta
        } else {
            400.0 + tb
        };
        near(
            native.power_enthalpy_derivatives_w_kg_j[side],
            expected_tangent[side] / capacity,
            2e-12,
        );
        let mut lower = energies;
        let mut upper = energies;
        lower[side] -= 0.25;
        upper[side] += 0.25;
        let finite_difference = (link.evaluate(upper, || Ok(())).unwrap().power_a_to_b_w
            - link.evaluate(lower, || Ok(())).unwrap().power_a_to_b_w)
            / 0.5;
        near(
            native.power_enthalpy_derivatives_w_kg_j[side],
            finite_difference,
            2e-8,
        );
    }
}

#[test]
fn g3_original_kirchhoff_geometry_retains_perfect_contact_and_reversed_heat_flow() {
    let (source_a, source_b) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (source_a.cv(), source_b.cp());
    let a = IsochoricConductionLaw::try_new(&source_a.card, &cv, &source_a.conductivity, SELECTION)
        .unwrap();
    let b = IsobaricConductionLaw::try_new(&source_b.card, &cp, &source_b.conductivity, SELECTION)
        .unwrap();
    let (layer_a, layer_b) = (
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    );
    let geometries = [
        LayerGeometry::planar([0.01, 0.008], AREA).unwrap(),
        LayerGeometry::cylindrical([0.01, 0.03], 0.02, 0.2).unwrap(),
        LayerGeometry::spherical([0.03, 0.01], 0.02).unwrap(),
    ];
    for geometry in geometries {
        let area = geometry.contact_area_m2();
        let thermal = timeline([&cv, &cp], [640.0, 360.0], [0.4, 0.32], area);
        let link = IntegratedCaloricLink::new(
            &thermal,
            LINK,
            [&layer_a, &layer_b],
            geometry,
            0.0,
            layer_controls(),
            || Ok(()),
        )
        .unwrap();
        let distances = geometry.equivalent_distances_m();
        let denominator = distances[0] / 0.04 + distances[1] / 0.02;
        for [ta, tb] in [[640.0_f64, 360.0_f64], [360.0, 640.0], [475.0, 475.0]] {
            let response = link
                .evaluate([energy(0, ta), energy(1, tb)], || Ok(()))
                .unwrap();
            let expected = area * (ta * ta - tb * tb) / (2.0 * denominator);
            near(response.power_a_to_b_w, expected, 2e-8);
            near(
                response.endpoint_derivative_w_k[0],
                area * ta / denominator,
                2e-9,
            );
            near(
                response.endpoint_derivative_w_k[1],
                -area * tb / denominator,
                2e-9,
            );
            let interface = (ta * ta - 2.0 * (expected / area) * distances[0] / 0.04).sqrt();
            near(response.interface_temperature_k[0], interface, 2e-9);
            near(response.interface_temperature_k[1], interface, 2e-9);
        }
    }
}

#[test]
fn g1_large_native_implicit_step_closes_original_u_h_and_entropy_without_frozen_k() {
    let (source_a, source_b) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (source_a.cv(), source_b.cp());
    let a = IsochoricConductionLaw::try_new(&source_a.card, &cv, &source_a.conductivity, SELECTION)
        .unwrap();
    let b = IsobaricConductionLaw::try_new(&source_b.card, &cp, &source_b.conductivity, SELECTION)
        .unwrap();
    let (layer_a, layer_b) = (
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    );
    let (initial_t, masses, distances, dt) = ([650.0, 330.0], [0.4, 0.2], [0.01, 0.005], 300.0);
    let thermal = timeline([&cv, &cp], initial_t, masses, AREA);
    let link = IntegratedCaloricLink::new(
        &thermal,
        LINK,
        [&layer_a, &layer_b],
        LayerGeometry::planar(distances, AREA).unwrap(),
        0.0,
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
    let mut calls = 0;
    let (next, report) = chain
        .advance_to(&thermal, dt, controls, &mut workspace, || {
            calls += 1;
            Ok(())
        })
        .unwrap();

    // Independent one-scalar backward-Euler balance with an analytic exact
    // Kirchhoff flux; this is not the production layer root or Newton solve.
    let initial_e = [energy(0, initial_t[0]), energy(1, initial_t[1])];
    let residual = |q: f64| {
        let ta = temperature(0, initial_e[0] - q / masses[0]);
        let tb = temperature(1, initial_e[1] + q / masses[1]);
        q - dt * AREA * (ta * ta - tb * tb) / (2.0 * (distances[0] / 0.04 + distances[1] / 0.02))
    };
    let (mut low, mut high) = (0.0, masses[0] * (initial_e[0] - energy(0, initial_t[1])));
    assert!(residual(low) < 0.0 && residual(high) > 0.0);
    for _ in 0..80 {
        let mid = (low + high) / 2.0;
        if residual(mid) < 0.0 {
            low = mid;
        } else {
            high = mid;
        }
    }
    let q = (low + high) / 2.0;
    let final_e = [
        next.cells()[&CellId(0)].specific_h_j_kg,
        next.cells()[&CellId(1)].specific_h_j_kg,
    ];
    near(final_e[0], initial_e[0] - q / masses[0], 5e-6);
    near(final_e[1], initial_e[1] + q / masses[1], 5e-6);
    let final_t = [temperature(0, final_e[0]), temperature(1, final_e[1])];
    assert!(final_t[1] > 475.0 && final_t[0] < 550.0);
    near(
        masses[0] * (final_e[0] - initial_e[0]) + masses[1] * (final_e[1] - initial_e[1]),
        0.0,
        1e-6,
    );
    let entropy = masses[0]
        * (200.0 * (final_t[0] / initial_t[0]).ln() + 0.5 * (final_t[0] - initial_t[0]))
        + masses[1] * (400.0 * (final_t[1] / initial_t[1]).ln() + final_t[1] - initial_t[1]);
    near(report.entropy.thermal_storage_change_j_k, entropy, 2e-8);
    assert!(entropy > 0.0 && report.entropy.conduction_generation_j_k > 0.0);
    near(
        report.entropy.conduction_generation_j_k,
        q * (1.0 / final_t[1] - 1.0 / final_t[0]),
        2e-8,
    );
    assert_eq!(report.energy.substeps, 1);
    assert!(report.energy.residual_j.abs() <= 1e-6);
    assert_eq!(
        next.material_storage_basis(CV_MATERIAL).unwrap(),
        ThermalStorageBasis::InternalEnergy
    );
    assert_eq!(
        next.material_storage_basis(CP_MATERIAL).unwrap(),
        ThermalStorageBasis::Enthalpy
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial_e[0]);
    assert_eq!(thermal.cells()[&CellId(1)].specific_h_j_kg, initial_e[1]);

    // Refuse late in the same real interval, then replay with the same scratch.
    assert!(calls > 8);
    let mut cancelled_calls = 0;
    let cancelled = chain.advance_to(&thermal, dt, controls, &mut workspace, || {
        cancelled_calls += 1;
        if cancelled_calls == calls - 1 {
            Err(DepositionError::Cancelled)
        } else {
            Ok(())
        }
    });
    assert!(matches!(
        cancelled,
        Err(ImplicitConductionError::Thermal(DepositionError::Cancelled))
    ));
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(thermal.cells()[&CellId(0)].specific_h_j_kg, initial_e[0]);
    let (retry, replay_report) = chain
        .advance_to(&thermal, dt, controls, &mut workspace, || Ok(()))
        .unwrap();
    assert_eq!(retry.cells(), next.cells());
    assert_eq!(retry.links(), next.links());
    assert_eq!(replay_report, report);
}

#[test]
fn g4_integrated_faces_refuse_foreign_inventory_geometry_and_unsupported_trials() {
    let (source_a, source_b) = (Fixture::new(true), Fixture::new(false));
    let (cv, cp) = (source_a.cv(), source_b.cp());
    let a = IsochoricConductionLaw::try_new(&source_a.card, &cv, &source_a.conductivity, SELECTION)
        .unwrap();
    let b = IsobaricConductionLaw::try_new(&source_b.card, &cp, &source_b.conductivity, SELECTION)
        .unwrap();
    let (layer_a, layer_b) = (
        CaloricConductionLayer::from_isochoric(&a).unwrap(),
        CaloricConductionLayer::from_isobaric(&b).unwrap(),
    );
    let thermal = timeline([&cv, &cp], [650.0, 350.0], [0.4, 0.32], AREA);
    let geometry = LayerGeometry::planar([0.01, 0.008], AREA).unwrap();
    let build = |layers, geometry| {
        IntegratedCaloricLink::new(
            &thermal,
            LINK,
            layers,
            geometry,
            0.0,
            layer_controls(),
            || Ok(()),
        )
    };
    assert!(build([&layer_b, &layer_a], geometry).is_err());
    assert!(
        build(
            [&layer_a, &layer_b],
            LayerGeometry::planar([0.01, 0.008], 2.0 * AREA).unwrap()
        )
        .is_err()
    );
    let cloned_cv = cv.clone();
    let cloned_a = IsochoricConductionLaw::try_new(
        &source_a.card,
        &cloned_cv,
        &source_a.conductivity,
        SELECTION,
    )
    .unwrap();
    let cloned_layer = CaloricConductionLayer::from_isochoric(&cloned_a).unwrap();
    assert!(build([&cloned_layer, &layer_b], geometry).is_err());
    assert!(matches!(
        IntegratedCaloricLink::new(
            &thermal,
            LINK,
            [&layer_a, &layer_b],
            geometry,
            0.0,
            layer_controls(),
            || Err(DepositionError::Cancelled)
        ),
        Err(DepositionError::Cancelled)
    ));
    let link = build([&layer_a, &layer_b], geometry).unwrap();
    for invalid in [-1.0, energy(0, 700.0) + 1.0, f64::NAN] {
        assert!(
            link.evaluate([invalid, energy(1, 350.0)], || Ok(()))
                .is_err()
        );
    }
    assert!(matches!(
        link.evaluate([energy(0, 650.0), energy(1, 350.0)], || Err(
            DepositionError::Cancelled
        )),
        Err(DepositionError::Cancelled)
    ));
    let mut wrong_geometry = link.native_geometry();
    wrong_geometry.half_distances_m[0] *= 1.01;
    assert!(
        link.state_at_specific_enthalpies(
            [energy(0, 650.0), energy(1, 350.0)],
            wrong_geometry,
            thermal.contact_patch(LINK).unwrap(),
            &mut || Ok(())
        )
        .is_err()
    );
    let mut wrong_patch = thermal.contact_patch(LINK).unwrap();
    wrong_patch.link.a = CellId(1);
    assert!(
        link.state_at_specific_enthalpies(
            [energy(0, 650.0), energy(1, 350.0)],
            link.native_geometry(),
            wrong_patch,
            &mut || Ok(())
        )
        .is_err()
    );
    assert_eq!(thermal.time_s(), 0.0);
    assert_eq!(
        thermal.cells()[&CellId(0)].specific_h_j_kg,
        energy(0, 650.0)
    );
    assert_eq!(thermal.contact_patch(LINK).unwrap().link.a, CellId(0));
}
