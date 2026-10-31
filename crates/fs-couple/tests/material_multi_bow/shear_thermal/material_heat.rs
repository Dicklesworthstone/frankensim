//! Manufactured distributed material/heat feedback, not measured horsehair or
//! rosin. This checks the original coupled dynamics and finite receiving mass.
use super::*;
use fs_couple::fiber_ribbon::{FiberRibbonGeometry, FiberRibbonResolution, TensionedFiberRibbon};
use fs_couple::render::schedule::force::coupled::coefficients::heat::MaterialHeatBudget;
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer,
    axial::YOUNG_VISCOSITY_PROPERTY, relaxation::LayerRelaxationProperties};
use fs_qty::Time;

const BODY: CellId = CellId(11);

fn fiber(gate: &CancelGate, temperature_dependent_time: bool) -> (MaterialCard, ResolvedHeatCapacityCurve, TensionedFiberRibbon) {
    let (constant, _) = coating(false, false);
    let mut claims = constant.claims().clone();
    for (name, dims, value) in [
        ("young_modulus", Pressure::DIMS, 3e9),
        (ETA, DynViscosity::DIMS, 0.0),
        (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS, 0.0),
        ("delta", Pressure::DIMS, 1e9), ("tau", Time::DIMS, 0.0005),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: if name == "tau" && temperature_dependent_time {
                PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(280.0, 0.0006), (340.0, 0.0003)], dims }
            } else { PropertyValue::Scalar { value, dims } },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 340.0).with("omega", 0.0, 1e6),
            uncertainty: UncertaintyModel::Unstated, interpolation: if name == "tau" && temperature_dependent_time {
                InterpolationPolicy::LinearInside
            } else { InterpolationPolicy::ConstantWithinValidity },
            observations: vec![], provenance: Provenance {
                source: "manufactured finite material-heat fixture".into(), license: "CC0-1.0".into(), artifact: None,
            },
        }).unwrap();
    }
    let card = MaterialCard::assemble(MaterialStateId { chemistry: "manufactured-hereditary-fiber".into(),
        phase: "solid".into(), process: "fixture".into(), revision: 0 }, claims, vec![]).unwrap();
    let requirements = [("density", Density::DIMS), ("young_modulus", Pressure::DIMS),
        (ETA, DynViscosity::DIMS), (YOUNG_VISCOSITY_PROPERTY, DynViscosity::DIMS),
        ("delta", Pressure::DIMS), ("tau", Time::DIMS)]
        .into_iter().map(|(name, dims)| ScalarPropertyRequirement::try_new(name, dims,
            ScalarAdmissibility::Finite).unwrap()).collect::<Vec<_>>();
    let source = resolve_material_state_point(&card, &point(300.0), &requirements,
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    let section = ConcentricElasticSection::try_new(0.0, vec![ConcentricMaterialLayer {
        outer_radius_m: 40e-6, material: source,
    }]).unwrap();
    let spectrum = section.clone().with_prony_bending(&[LayerRelaxationProperties {
        equilibrium_modulus: "young_modulus".into(), branches: vec![("delta".into(), "tau".into())],
    }]).unwrap();
    let ribbon = TensionedFiberRibbon::try_new(
        FiberRibbonGeometry { length_m: 0.6, width_m: 0.01, fiber_count: 180, tension_n: 40.0 },
        FiberRibbonResolution { cells: 1, axial_modes: 1, transverse_modes: 1 },
        section.clone().with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY).unwrap(),
        section.with_kelvin_voigt_bending(ETA).unwrap(), 48_000,
        fs_couple::modal_acoustic_time::ModalAcousticTimeBudget::audible_reference(), gate,
    ).unwrap().with_material_relaxation(spectrum, 16, gate).unwrap();
    let caloric = ResolvedHeatCapacityCurve::try_new(&card, HeatCapacitySourceQuery {
        heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
        density: &PropertyKey::new("density", Density::DIMS), temperature_axis: "T",
        lower: &point(280.0), upper: &point(340.0), selection: SelectionPolicy::SingleClaimOnly,
    }, SolidLiquidPhase::Solid, -20.0 * CP, 1e-8).unwrap();
    (card, caloric, ribbon)
}

fn finite_bodies<'a>(film: &'a EquilibriumEnthalpyPhaseCurve,
    body: &'a EquilibriumEnthalpyPhaseCurve, body_volume: f64) -> DepositionSimulation<'a> {
    let mut events = Vec::new();
    for (id, material_id, volume) in [(FILM, MaterialId(0), AREA * THICKNESS),
        (BODY, MaterialId(1), body_volume)] {
        events.push(TimedEvent { time_s: 0.0, event: Event::Birth(DepositedCell {
            id, segment_id: id.0, body_id: id.0, material_id, mass_kg: DENSITY * volume,
            specific_h_j_kg: 0.0, geometry: CellGeometry { start_m: [0.0; 3],
                end_m: [0.6, 0.0, 0.0], volume_m3: volume, surface_area_m2: 2.0 * AREA },
            boundary: Boundary { ambient_k: 300.0, htc_w_m2_k: 0.0, reservoir: None, source_w: 0.0 },
        }) });
    }
    events.push(TimedEvent { time_s: 0.0, event: Event::Contact(ThermalLink {
        id: 1, a: FILM, b: BODY, area_m2: AREA,
        conductance_w_k: series_conductance_w_k(AREA, [0.5*THICKNESS, 40e-6], [0.2, 0.2]).unwrap(),
    }) });
    DepositionSimulation::new(events, &[(MaterialId(0), film as &dyn ThermalLaw),
        (MaterialId(1), body as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 2, max_links: 1, max_events: 3, max_substeps: 4 }).unwrap()
        .advance_to(0.0, StepConfig { max_dt_s: 1e-5, energy_tolerance_j: 1e-14 }, || Ok(())).unwrap().0
}

#[test]
fn distributed_fiber_heat_enters_its_real_mass_and_changes_original_film_traction() {
    with_cx(|cx, gate| {
        let (film_card, film_caloric) = coating(true, true);
        let (body_card, body_caloric, fiber) = fiber(gate, false);
        let make = || base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(fiber.clone(), 0.4, ribbon_limits(), gate).unwrap()
            .with_bow_shear_layers(vec![(0, layer(&film_card))], gate).unwrap();
        let (mut bank, mut unrouted) = (make(), make());
        let budget = MaterialHeatBudget { max_groups: 1, max_dimension: 256,
            max_dense_bytes: 32 * 1024 * 1024, max_multiply_terms: 10_000_000_000 };
        let mut route = bank.fiber_material_heat_destination(0, 0, BODY, &body_card,
            body_caloric.curve(), "T", budget, gate).unwrap();
        let volume = core::f64::consts::PI * (40e-6_f64).powi(2) * 0.6 * 180.0;
        assert!((route.volume_m3 - volume).abs() < 1e-14 * volume);
        assert_eq!(route.channel.groups(), &[0..2]); // axial + normal, same material
        let (mut heat, mut reference_heat) = (
            finite_bodies(film_caloric.curve(), body_caloric.curve(), volume),
            finite_bodies(film_caloric.curve(), body_caloric.curve(), volume));
        let bound = binding(&film_card, film_caloric.curve());
        let first = bank.source().fiber_ribbon().unwrap().component_start();
        let first_mode = bank.system().components()[..first].iter().map(|m| m.modes().len()).sum::<usize>();
        let mut force = vec![0.0; bank.system().mode_count()];
        force[first_mode] = 10.0;
        let drives = [drive(0.001, 0.1)];
        // Source-axis aliases and a made-up thermal mass cannot advance either clock.
        route.temperature_axis = "omega";
        assert!(bank.step_with_material_contact_thermal_motion(&mut heat, &bound, &[], &[],
            std::slice::from_ref(&route), controls(2), &[0.0; 2], &force, &drives,
            None, None, None, gate, |s, _| other(s)).is_err());
        route.temperature_axis = "T";
        route.volume_m3 *= 2.0;
        assert!(bank.step_with_material_contact_thermal_motion(&mut heat, &bound, &[], &[],
            std::slice::from_ref(&route), controls(2), &[0.0; 2], &force, &drives,
            None, None, None, gate, |s, _| other(s)).is_err());
        route.volume_m3 *= 0.5;
        assert_eq!(bank.system().samples_rendered(), 0);
        assert_eq!(heat.time_s(), 0.0);
        let (mut body_heat, mut contact_heat, mut force_difference) = (0.0, 0.0, 0.0);
        for _ in 0..32 {
            let frame = bank.step_with_material_contact_thermal_motion(&mut heat, &bound, &[], &[],
                std::slice::from_ref(&route), controls(2), &[0.0; 2], &force, &drives,
                None, None, None, gate, |s, _| other(s)).unwrap();
            let reference = unrouted.step_with_layer_thermal_network(&mut reference_heat, &bound,
                controls(2), &[0.0; 2], &force, &drives, gate, other).unwrap();
            body_heat += frame.routed_material_heat_j;
            contact_heat += frame.routed_plastic_heat_j + frame.routed_relaxation_heat_j
                + frame.routed_reversible_heat_j;
            force_difference += (frame.mechanics.ports[0].traction_n - reference.mechanics.ports[0].traction_n).abs();
            assert!(frame.routed_material_heat_j >= 0.0);
            let native = bank.system().coupled_frame().unwrap().material_dissipation_j;
            assert!((frame.routed_material_heat_j - native).abs() < 2e-12);
            assert!(frame.total_energy_residual_j.abs() <= controls(2).total_energy_tolerance_j);
            assert_eq!(bank.bow_shear_layer(0).unwrap().temperature_kelvin(), heat.temperature_k(FILM).unwrap());
        }
        assert!(body_heat > 0.0);
        assert!((total_heat(&heat) - body_heat - contact_heat).abs() < 1e-10);
        assert!(heat.temperature_k(BODY).unwrap() > reference_heat.temperature_k(BODY).unwrap());
        assert!(force_difference > 0.0);
        assert_eq!(bank.system().samples_rendered(), 32);
    });
}

#[test]
fn frozen_temperature_dependent_fiber_time_cannot_masquerade_as_finite_thermal_feedback() {
    with_cx(|cx, gate| {
        let (film_card, film_caloric) = coating(true, true);
        let (body_card, body_caloric, fiber) = fiber(gate, true);
        let mut bank = base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(fiber, 0.4, ribbon_limits(), gate).unwrap()
            .with_bow_shear_layers(vec![(0, layer(&film_card))], gate).unwrap();
        let route = bank.fiber_material_heat_destination(0, 0, BODY, &body_card,
            body_caloric.curve(), "T", MaterialHeatBudget { max_groups: 1, max_dimension: 256,
                max_dense_bytes: 32 * 1024 * 1024, max_multiply_terms: 10_000_000_000 }, gate).unwrap();
        let mut heat = finite_bodies(film_caloric.curve(), body_caloric.curve(), route.volume_m3);
        let force = vec![0.0; bank.system().mode_count()];
        let original_memory = bank.system().material_memory().to_vec();
        assert!(bank.step_with_material_contact_thermal_motion(&mut heat,
            &binding(&film_card, film_caloric.curve()), &[], &[], &[route], controls(2),
            &[0.0; 2], &force, &[drive(0.001, 0.1)], None, None, None, gate,
            |s, _| other(s)).is_err());
        assert_eq!(bank.system().material_memory(), original_memory);
        assert_eq!(bank.system().samples_rendered(), 0);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(total_heat(&heat), 0.0);
    });
}
