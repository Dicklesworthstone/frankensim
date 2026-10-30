//! G1/G3: one nonlinear endpoint for material interfaces, photons and finite gas.
//! All values are manufactured source data, not a measured specimen.
use super::*;
use fs_conduction::transient::deposition::{contact::*, radiation::*};
use fs_conduction::{
    AREA_SPECIFIC_THERMAL_RESISTANCE_DIMS as RD, AREA_SPECIFIC_THERMAL_RESISTANCE_PROPERTY as RP,
    EMISSIVITY_DIMS, SURFACE_EMISSIVITY_PROPERTY, SurfaceEmissivity,
};
use fs_couple::circuit_enthalpy::contacts::heat::*;
use fs_couple::circuit_enthalpy::contacts::material::ElectricalContactNetwork;
use fs_matdb::{ClaimSelection, InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_material::thermal::interface::ThermalInterfacePoint;

#[allow(dead_code)]
#[path = "../../../../../fs-conduction/tests/radiation_material/common.rs"]
mod radiation_geometry;

fn radiating_card(slope: f64, emissivity: f64) -> (MaterialCard, ResolvedHeatCapacityCurve) {
    let (original, _) = nonlinear_card(slope);
    let mut claims = original.claims().clone();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(SURFACE_EMISSIVITY_PROPERTY, EMISSIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: emissivity,
                dims: EMISSIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained()
                .with("T", 280.0, 600.0)
                .with("lot", 7.0, 7.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured gas/solid emissivity".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    let card = MaterialCard::assemble(original.id().clone(), claims, vec![]).unwrap();
    let caloric = caloric_for(&card);
    (card, caloric)
}

fn interface(cards: [&MaterialCard; 2], resistance: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(RP, RD),
            value: PropertyValue::Scalar {
                value: resistance,
                dims: RD,
            },
            validity: ValidityDomain::unconstrained()
                .with("aT", 280.0, 600.0)
                .with("bT", 280.0, 600.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: vec![],
            provenance: Provenance {
                source: "manufactured finite-cell thermal contact".into(),
                license: "CC0-1.0".into(),
                artifact: None,
            },
        })
        .unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec {
            material: cards[0].id().clone(),
            texture_frame: "a".into(),
        },
        SurfaceSpec {
            material: cards[1].id().clone(),
            texture_frame: "b".into(),
        },
        SystemContext {
            medium: "transparent manufactured gas".into(),
            third_body: None,
            environment: "fixed geometry".into(),
            history: "declared".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

// Independent endpoint balances: Fourier+contact series resistance, analytic
// two-gray-plate T^4 exchange, and the declared two half-layer gas conductances.
fn supplies() -> [f64; 3] {
    let [a, b, g] = [400.0_f64, 440.0, 410.0];
    let h = 0.5 / (0.3 / (0.005 * a) + 0.7 / (0.01 * b) + 0.2);
    let radiation = 5.670_374_419e-8 * (a.powi(4) - b.powi(4)) / (1.0 / 0.3 + 1.0 / 0.6 - 1.0);
    let internal = 2.0 * (h * (a - b) + radiation);
    let qa = 2.0 * 0.1 / (0.001 / 0.03 + 0.02 / (0.005 * a)) * (a - g);
    let qb = 2.0 * 0.2 / (0.001 / 0.03 + 0.02 / (0.01 * b)) * (b - g);
    let cg = 2.5 * 100000.0 * 0.1 / 390.0;
    [
        10.0 * (a - 340.0) + internal + qa,
        10.0 * (b - 360.0) - internal + qb,
        cg * (g - 390.0) - qa - qb,
    ]
}

fn run(emissivity: f64, resistance: f64, complete: bool) -> [f64; 4] {
    with_cx(|cx| {
        let (a, ca) = radiating_card(0.005, emissivity);
        let (b, cb) = radiating_card(0.01, 0.6);
        let cards = [&a, &b];
        let curves = [ca.curve(), cb.curve()];
        let heat = supplies();
        let thermal = thermal_pair(curves, [heat[0] / 2.0, heat[1] / 2.0]);
        let context = point(340.0);
        let source = SolidEnthalpy::try_new(
            &thermal,
            &PAIR_NODES,
            &PAIR_TETS,
            vec![
                super::super::super::source(30, &a, curves[0], &context),
                super::super::super::source(10, &b, curves[1], &context),
            ],
            TetAssemblyBudget::standard(),
            0.1,
            config().mesh_frame,
            cx,
        )
        .unwrap();
        let key = PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS);
        let cells = material_cells(cards, curves, &context, &key);
        let bulk = MaterialNetwork {
            cells: &cells,
            resistors: &[],
            links: &LINKS,
        };
        let electrical = ElectricalContactNetwork {
            bulk: &bulk,
            contacts: &[],
            electrical_bases: &[],
        };
        let card = interface(cards, resistance);
        let p = QueryPoint::new()
            .with("aT", 340.0)
            .unwrap()
            .with("bT", 360.0)
            .unwrap();
        let contact = [MaterialContactLink {
            link: 7,
            sides: [0, 1].map(|i| ContactCell {
                cell: cells[i].cell,
                material_id: cells[i].material_id,
                curve: curves[i],
            }),
            point: ThermalInterfacePoint {
                interface: &card,
                materials: cards,
                context: &p,
                temperature_axes: ["aT", "bT"],
                liquid_fraction_axes: [None; 2],
                mean_temperature_axis: None,
            },
            selection: ClaimSelection::Policy(SelectionPolicy::SingleClaimOnly),
        }];
        let surfaces = [0, 1].map(|i| RadiatingCell {
            surface: i,
            cell: cells[i].cell,
            material_id: cells[i].material_id,
            curve: curves[i],
            material: cells[i].point,
            selection: MaterialPropertySelection::SingleClaimOnly,
        });
        let emissivities = [0, 1].map(|i| {
            SurfaceEmissivity::from_thermal_state(
                "manufactured",
                cells[i].point,
                curves[i]
                    .state_at_specific_enthalpy(thermal.cells()[&cells[i].cell].specific_h_j_kg)
                    .unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
            .0
        });
        let enclosure = radiation_geometry::enclosure(emissivities);
        let radiation = MaterialRadiation {
            enclosure: &enclosure,
            cells: &surfaces,
            power_tolerance_w: 1e-8,
            radiosity_tolerance_w_m2: 1e-8,
        };
        let assembly = MaterialHeatNetwork {
            electrical: &electrical,
            thermal_contacts: if complete { &contact } else { &[] },
            radiation: complete.then_some(&radiation),
        };
        let prepared = assembly.prepare(cx, &thermal).unwrap();
        let gas = gas(390.0, 0.1);
        let patches = contacts([0.1, 0.2]);
        let next = source
            .step_with_isochoric_gas_heat_network(
                &thermal,
                &gas,
                2.0,
                heat[2],
                &patches,
                &prepared,
                solve(),
                &[0.0; 3],
                cx,
            )
            .unwrap();
        let report = next.report();
        let sample = report.transport.as_ref().unwrap();
        assert_eq!(sample.thermal_contacts.len(), usize::from(complete));
        assert_eq!(sample.radiation.is_some(), complete);
        near(report.external_thermal_heat_j, heat[0] + heat[1], 2e-7);
        near(report.combined_energy_residual_j, 0.0, 2e-7);
        near(
            next.thermal().energy_j() - thermal.energy_j()
                + gas.mass_kg()
                    * (next.gas().state().specific_internal_energy_j_kg()
                        - gas.state().specific_internal_energy_j_kg()),
            heat.iter().sum(),
            2e-7,
        );
        for patch in &report.contacts {
            let original = sample
                .material
                .cells
                .iter()
                .find(|cell| cell.cell == patch.cell)
                .unwrap()
                .conductivity
                .as_ref()
                .unwrap();
            assert_eq!(original.query(), patch.source.query());
            assert_eq!(
                original
                    .property_by_key(&key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected,
                patch
                    .source
                    .property_by_key(&key)
                    .unwrap()
                    .answer()
                    .receipt
                    .selected
            );
        }
        if complete && emissivity == 0.3 && resistance == 0.2 {
            near(
                next.thermal().temperature_k(CellId(30)).unwrap(),
                400.0,
                2e-7,
            );
            near(
                next.thermal().temperature_k(CellId(10)).unwrap(),
                440.0,
                2e-7,
            );
            near(next.gas().state().state().temperature, 410.0, 2e-7);
            near(
                next.gas().state().state().pressure,
                100000.0 * 410.0 / 390.0,
                1e-4,
            );
        }
        if !complete {
            let old = bulk.prepare(&thermal, || false).unwrap();
            let old = source
                .step_with_isochoric_gas_transport(
                    &thermal,
                    &gas,
                    2.0,
                    heat[2],
                    &patches,
                    &old,
                    solve(),
                    &[0.0; 3],
                    cx,
                )
                .unwrap();
            assert_eq!(old.thermal().cells(), next.thermal().cells());
            assert_eq!(old.gas().state(), next.gas().state());
            assert_eq!(old.report().transport.as_ref().unwrap(), &sample.material);
        }
        let mut refused = solve();
        refused.maximum_trials = 1;
        assert!(
            source
                .step_with_isochoric_gas_heat_network(
                    &thermal, &gas, 2.0, heat[2], &patches, &prepared, refused, &[0.0; 3], cx
                )
                .is_err()
        );
        assert_eq!(thermal.time_s(), 0.0);
        assert_eq!(gas.state().state().temperature, 390.0);
        [
            next.thermal().temperature_k(CellId(30)).unwrap(),
            next.thermal().temperature_k(CellId(10)).unwrap(),
            next.gas().state().state().temperature,
            next.gas().state().state().pressure,
        ]
    })
}

#[test]
fn complete_heat_network_matches_independent_cell_photon_gas_balances() {
    run(0.3, 0.2, true);
}

#[test]
fn source_emissivity_and_contact_swaps_change_the_same_finite_gas_state() {
    let original = run(0.3, 0.2, true);
    for changed in [run(0.8, 0.2, true), run(0.3, 2.0, true)] {
        for i in 0..4 {
            assert!((changed[i] - original[i]).abs() > 1e-6);
        }
    }
}

#[test]
fn empty_complete_heat_network_reproduces_original_material_transport() {
    run(0.3, 0.2, false);
}
