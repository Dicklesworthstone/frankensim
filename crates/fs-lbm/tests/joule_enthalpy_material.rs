#![cfg(feature = "enthalpy-coupling")]
//! G1/G3 material substitution through actual current, latent heat and transport.
//! The source cards below are synthetic numerical fixtures, not alloy data.

use fs_conduction::transient::deposition::*;
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::conjugate::enthalpy::transport::material::{MaterialCell, MaterialTransportState};
use fs_lbm::conjugate::enthalpy::transport::*;
use fs_lbm::conjugate::enthalpy::{JouleEnthalpyError, JouleEnthalpyInput};
use fs_lbm::conjugate::{ElectricSetup, Electrode, ElectrodeDrive, VoxelDomain};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::conductor::ElectricalPropertyBasis;
use fs_material::phase::{EnthalpyPhaseKnot, EquilibriumEnthalpyPhaseCurve};
use fs_material::state_point::{MaterialPropertySelection, THERMAL_CONDUCTIVITY_DIMS};
use fs_material::thermal::source::ThermalMaterialPoint;
use fs_qty::Dims;

const IDS: [CellId; 3] = [CellId(11), CellId(23), CellId(41)];
fn config() -> StepConfig {
    StepConfig {
        max_dt_s: 0.01,
        energy_tolerance_j: 1e-8,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-7, "{a} != {b}");
}
fn key() -> PropertyKey {
    PropertyKey::new("thermal_conductivity", THERMAL_CONDUCTIVITY_DIMS)
}

fn card(
    basis: ElectricalPropertyBasis,
    k: f64,
    varying: bool,
    upper_fraction: f64,
) -> MaterialCard {
    let mut claims = ClaimSet::new();
    let electrical = match (basis, varying) {
        (ElectricalPropertyBasis::Resistivity, false) => [2.0, 2.0],
        (ElectricalPropertyBasis::Conductivity, false) => [0.5, 0.5],
        (ElectricalPropertyBasis::Resistivity, true) => [1.0, 2.0],
        (ElectricalPropertyBasis::Conductivity, true) => [1.0, 0.5],
    };
    for (key, values) in [
        (key(), [k, 2.0 * k]),
        (
            PropertyKey::new(basis.property_name(), basis.dimensions()),
            electrical,
        ),
    ] {
        let dims = key.quantity().dims();
        claims
            .insert_claim(PropertyClaim {
                key,
                value: PropertyValue::Curve {
                    abscissa: "liquid".into(),
                    abscissa_dims: Dims([0; 6]),
                    knots: vec![(0.0, values[0]), (1.0, values[1])],
                    dims,
                },
                validity: ValidityDomain::unconstrained()
                    .with("T", 300.0, 810.0)
                    .with("liquid", 0.0, upper_fraction)
                    .with("pressure", 1.0, 2.0),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::LinearInside,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic electrothermal source".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic".into(),
            phase: "declared solid/liquid chart".into(),
            process: "test".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn law(card: &MaterialCard) -> EquilibriumEnthalpyPhaseCurve {
    EquilibriumEnthalpyPhaseCurve::try_new(
        card.content_hash(),
        [
            (0.0, 300.0, 0.0),
            (20.0, 310.0, 0.0),
            (120.0, 310.0, 1.0),
            (1120.0, 810.0, 1.0),
        ]
        .map(|(h, t, f)| EnthalpyPhaseKnot {
            specific_enthalpy_j_kg: h,
            temperature_k: t,
            liquid_mass_fraction: f,
            bulk_density_kg_m3: 1.0,
        })
        .to_vec(),
    )
    .unwrap()
}

fn initial(curve: &EquilibriumEnthalpyPhaseCurve) -> DepositionSimulation<'_> {
    let mut events: Vec<_> = IDS
        .iter()
        .enumerate()
        .map(|(c, &id)| TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id,
                segment_id: c as u64,
                body_id: 1,
                material_id: MaterialId(1),
                mass_kg: 1.0,
                specific_h_j_kg: if c == 1 { 40.0 } else { 0.0 },
                geometry: CellGeometry {
                    start_m: [c as f64, 0.5, 0.5],
                    end_m: [c as f64 + 1.0, 0.5, 0.5],
                    volume_m3: 1.0,
                    surface_area_m2: 6.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        })
        .collect();
    for c in 0..2 {
        events.push(TimedEvent {
            time_s: 0.0,
            event: Event::Contact(ThermalLink {
                id: c as u64,
                a: IDS[c],
                b: IDS[c + 1],
                area_m2: 1.0,
                conductance_w_k: 0.0,
            }),
        });
    }
    DepositionSimulation::new(
        events,
        &[(MaterialId(1), curve as &dyn ThermalLaw)],
        DepositionBudget {
            max_cells: 3,
            max_links: 2,
            max_events: 5,
            max_substeps: 10000,
        },
    )
    .unwrap()
    .advance_to(0.0, config(), || Ok(()))
    .unwrap()
    .0
}

fn context() -> QueryPoint {
    QueryPoint::new()
        .with("T", 300.0)
        .unwrap()
        .with("liquid", 0.0)
        .unwrap()
        .with("pressure", 1.5)
        .unwrap()
}
fn point<'a>(card: &'a MaterialCard, context: &'a QueryPoint) -> ThermalMaterialPoint<'a> {
    ThermalMaterialPoint {
        material: card,
        context,
        temperature_axis: "T",
        liquid_fraction_axis: Some("liquid"),
    }
}
fn bindings<'a>(
    thermal: &DepositionSimulation<'_>,
    curve: &'a EquilibriumEnthalpyPhaseCurve,
    card: &'a MaterialCard,
    context: &'a QueryPoint,
    basis: ElectricalPropertyBasis,
) -> Vec<MaterialCell<'a>> {
    IDS.iter()
        .map(|&id| {
            MaterialCell::try_new(
                thermal,
                id,
                curve,
                point(card, context),
                &key(),
                Some(basis),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .unwrap()
        })
        .collect()
}
fn electric() -> ElectricSetup {
    ElectricSetup {
        conductors: vec![],
        electrodes: vec![
            Electrode {
                lo: [0, 0, 0],
                hi: [1, 1, 1],
                drive: ElectrodeDrive::Voltage(10.0),
            },
            Electrode {
                lo: [2, 0, 0],
                hi: [3, 1, 1],
                drive: ElectrodeDrive::Voltage(0.0),
            },
        ],
    }
}
fn input<'a>(domain: &'a VoxelDomain, setup: &'a ElectricSetup) -> JouleEnthalpyInput<'a> {
    JouleEnthalpyInput {
        domain,
        setup,
        contacts: &[],
        cell_ids: &IDS,
        max_cells: 3,
        max_current_solves: 10001,
        power_tolerance_w: 1e-8,
        energy_tolerance_j: 1e-8,
    }
}
fn series(
    p: &ContactPatch,
    s: [&EnthalpyTransportSample<MaterialTransportState>; 2],
) -> Result<(f64, ()), JouleEnthalpyError> {
    // Explicit one-metre Cartesian cells, no additional thermal contact resistance.
    let [a, b] = s.map(|x| x.response.conductivity_w_m_k());
    Ok((
        if a == 0.0 || b == 0.0 {
            0.0
        } else {
            p.link.area_m2 / (0.5 / a + 0.5 / b)
        },
        (),
    ))
}

#[test]
fn g3_material_data_change_heat_phase_and_current_on_one_isothermal_plateau() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let electric = electric();
    let context = context();
    let mut outcomes = Vec::new();
    for k in [0.1, 2.0] {
        let card = card(ElectricalPropertyBasis::Resistivity, k, true, 1.0);
        let curve = law(&card);
        let old = initial(&curve);
        let bindings = bindings(
            &old,
            &curve,
            &card,
            &context,
            ElectricalPropertyBasis::Resistivity,
        );
        let (next, r) = advance_joule_enthalpy_with_material_transport(
            &old,
            0.5,
            config(),
            input(&domain, &electric),
            &[0, 1],
            &CancelGate::new(),
            |c, cell, t| bindings[c].resolve(cell, t),
            series,
        )
        .unwrap();
        let endpoint = &r.endpoint.cells[1].response;
        let state = endpoint.material();
        let f = state.phase().liquid_mass_fraction();
        near(state.phase().temperature_k(), 310.0);
        assert_eq!(state.material().card_identity(), card.content_hash());
        assert_eq!(state.material().state_coordinate("T"), Some(310.0));
        assert_eq!(state.material().state_coordinate("liquid"), Some(f));
        assert_eq!(state.material().state_coordinate("pressure"), Some(1.5));
        assert_eq!(state.material().properties().len(), 2);
        near(endpoint.conductivity_w_m_k(), k * (1.0 + f));
        near(endpoint.resistivity_ohm_m().unwrap(), 1.0 + f);
        near(
            r.joule.endpoint.electric.electrodes[0].current_a,
            10.0 / (1.0 + f),
        );
        let applied = &r.last_applied.as_ref().unwrap().cells[1].response;
        assert_eq!(
            applied.material().phase().temperature_k(),
            state.phase().temperature_k()
        );
        assert!(applied.material().phase().liquid_mass_fraction() < f);
        near(next.energy_j() - old.energy_j(), r.joule.deposited_joule_j);
        assert_eq!(next.mass_kg(), old.mass_kg());
        outcomes.push((f, r.joule.endpoint.electric.electrodes[0].current_a));
    }
    assert!(outcomes[0].0 > outcomes[1].0);
    assert!(outcomes[0].1 < outcomes[1].1);
}

#[test]
fn g1_conductivity_or_resistivity_source_drives_the_same_current_and_heat_owner() {
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let electric = electric();
    let context = context();
    let mut results = Vec::new();
    for basis in [
        ElectricalPropertyBasis::Resistivity,
        ElectricalPropertyBasis::Conductivity,
    ] {
        let card = card(basis, 0.5, false, 1.0);
        let curve = law(&card);
        let old = initial(&curve);
        let bindings = bindings(&old, &curve, &card, &context, basis);
        let (next, r) = advance_joule_enthalpy_with_material_transport(
            &old,
            0.2,
            config(),
            input(&domain, &electric),
            &[0, 1],
            &CancelGate::new(),
            |c, cell, t| bindings[c].resolve(cell, t),
            series,
        )
        .unwrap();
        let source = r.endpoint.cells[1].response.material().conductor().unwrap();
        assert_eq!(source.source_basis(), basis);
        assert_eq!(
            source.source_property().requirement().name(),
            basis.property_name()
        );
        near(r.joule.endpoint.electric.electrodes[0].current_a, 5.0);
        near(next.energy_j() - old.energy_j(), 10.0);
        results.push((next.cells().clone(), r.joule));
    }
    assert_eq!(results[0], results[1]);
}

#[test]
fn g0_g4_unsupported_phase_after_heating_refuses_without_publishing() {
    let card = card(ElectricalPropertyBasis::Resistivity, 0.1, true, 0.2001);
    let curve = law(&card);
    let old = initial(&curve);
    let context = context();
    let bindings = bindings(
        &old,
        &curve,
        &card,
        &context,
        ElectricalPropertyBasis::Resistivity,
    );
    let cells = old.cells().clone();
    let domain = VoxelDomain::new(3, 1, 1, 1.0).unwrap();
    let electric = electric();
    let err = advance_joule_enthalpy_with_material_transport(
        &old,
        0.01,
        config(),
        input(&domain, &electric),
        &[0, 1],
        &CancelGate::new(),
        |c, cell, t| bindings[c].resolve(cell, t),
        series,
    )
    .err()
    .unwrap();
    assert!(err.to_string().contains("enthalpy.material"));
    assert_eq!(old.cells(), &cells);
    assert_eq!(old.time_s(), 0.0);
    near(old.temperature_k(IDS[1]).unwrap(), 310.0);
}

#[test]
fn g0_foreign_law_card_cell_or_temperature_cannot_enter_the_material_path() {
    let card = card(ElectricalPropertyBasis::Resistivity, 0.1, true, 1.0);
    let curve = law(&card);
    let old = initial(&curve);
    let context = context();
    let source = point(&card, &context);
    let duplicate_curve = curve.clone();
    assert!(
        MaterialCell::try_new(
            &old,
            IDS[1],
            &duplicate_curve,
            source,
            &key(),
            Some(ElectricalPropertyBasis::Resistivity),
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let bindings = bindings(
        &old,
        &curve,
        &card,
        &context,
        ElectricalPropertyBasis::Resistivity,
    );
    assert!(bindings[0].resolve(&old.cells()[&IDS[1]], 310.0).is_err());
    assert!(bindings[1].resolve(&old.cells()[&IDS[1]], 300.0).is_err());
    let wrong = PropertyKey::new("optical_extinction", Dims([0; 6]));
    assert!(
        MaterialCell::try_new(
            &old,
            IDS[1],
            &curve,
            source,
            &wrong,
            None,
            MaterialPropertySelection::SingleClaimOnly
        )
        .is_err()
    );
    let insulating = MaterialCell::try_new(
        &old,
        IDS[1],
        &curve,
        source,
        &key(),
        None,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
    .resolve(&old.cells()[&IDS[1]], 310.0)
    .unwrap();
    assert_eq!(insulating.resistivity_ohm_m(), None);
    assert!(insulating.material().conductor().is_none());
    assert_eq!(insulating.material().material().properties().len(), 1);
}

#[test]
fn g0_material_binding_retains_both_initial_sources_across_phase_evolution() {
    use fs_material::state_point::{
        ScalarAdmissibility, ScalarPropertyRequirement, resolve_material_state_point,
    };

    let basis = ElectricalPropertyBasis::Conductivity;
    for split_property in ["thermal_conductivity", basis.property_name()] {
        let mut claims = ClaimSet::new();
        for key in [
            key(),
            PropertyKey::new(basis.property_name(), basis.dimensions()),
        ] {
            let ranges: &[(f64, f64, f64)] = if key.name() == split_property {
                &[(0.0, 0.25, 1.0), (0.3, 1.0, 3.0)]
            } else {
                &[(0.0, 1.0, 1.0)]
            };
            for &(low, high, value) in ranges {
                let dims = key.quantity().dims();
                claims
                    .insert_claim(PropertyClaim {
                        key: key.clone(),
                        value: PropertyValue::Scalar { value, dims },
                        validity: ValidityDomain::unconstrained()
                            .with("T", 300.0, 810.0)
                            .with("liquid", low, high)
                            .with("pressure", 1.0, 2.0),
                        uncertainty: UncertaintyModel::Unstated,
                        interpolation: InterpolationPolicy::ConstantWithinValidity,
                        observations: vec![],
                        provenance: Provenance {
                            source: "synthetic disjoint source supports".into(),
                            license: "test-only".into(),
                            artifact: None,
                        },
                    })
                    .unwrap();
            }
        }
        let card = MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic".into(),
                phase: "declared solid/liquid chart".into(),
                process: "test".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap();
        let curve = law(&card);
        let old = initial(&curve);
        let context = context();
        let source = point(&card, &context);
        let binding = MaterialCell::try_new(
            &old,
            IDS[1],
            &curve,
            source,
            &key(),
            Some(basis),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        let mut cell = old.cells()[&IDS[1]].clone();
        let first = binding.resolve(&cell, 310.0).unwrap();
        let selected = first
            .material()
            .material()
            .property(split_property)
            .unwrap()
            .answer()
            .receipt
            .selected;

        cell.specific_h_j_kg = 44.0; // f=0.24 is within the original source.
        let within = binding.resolve(&cell, 310.0).unwrap();
        assert_eq!(
            within
                .material()
                .material()
                .property(split_property)
                .unwrap()
                .answer()
                .receipt
                .selected,
            selected
        );

        cell.specific_h_j_kg = 60.0; // f=0.4 is supported only by another claim.
        let state = curve
            .state_at_specific_enthalpy(cell.specific_h_j_kg)
            .unwrap();
        let requirements = [
            ScalarPropertyRequirement::try_with_key(&key(), ScalarAdmissibility::NonNegative)
                .unwrap(),
            basis.requirement().unwrap(),
        ];
        let unbound = resolve_material_state_point(
            &card,
            &source.at(state).unwrap(),
            &requirements,
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap();
        assert_ne!(
            unbound
                .property(split_property)
                .unwrap()
                .answer()
                .receipt
                .selected,
            selected
        );
        assert!(binding.resolve(&cell, 310.0).is_err());
    }
}
