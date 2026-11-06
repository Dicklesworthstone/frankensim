//! Independent P1/play balances through actual source-backed field/heat owners.
//! Manufactured fit coefficients establish accounting, not catalog validity.
#![cfg(feature = "magnetic-hysteresis")]

#[path = "../../fs-em/tests/nonlinear_eddy/common.rs"]
mod geometry;

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpyError, CircuitEnthalpySolve,
        magnetic::{MagneticRegion, MagneticThermalError, MagneticThermalNetwork},
        material::{MaterialCell, MaterialNetwork},
    },
    material_circuit::MaterialCircuitSolve,
};
use fs_em::planar::{
    eddy::{
        EddyConductor,
        hysteretic::{HystereticEddyControls, PlanarHystereticEddy},
    },
    nonlinear::{
        NonlinearMagneticControls,
        hysteretic::{PlanarHystereticControls, PlanarHystereticMagnetostatics},
    },
};
use fs_evidence::ValidityDomain;
use fs_feec::weighted::WeightedAssemblyLimits;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SelectionPolicy, UncertaintyModel,
};
use fs_material::{
    conductor::{ElectricalPropertyBasis, ResolvedConductorState},
    field_storage::{
        MAGNETIC_PERMEABILITY_DIMS,
        hysteresis::{
            source::{
                HysteresisBranchSources, HysteresisSources, ResolvedScalarHysteresis,
                ResolvedVectorHysteresis,
            },
            vector::VectorHysteresisControls,
        },
        permanent::{FIELD_STRENGTH_DIMS, FLUX_DENSITY_DIMS},
    },
    phase::{HeatCapacitySourceQuery, ResolvedHeatCapacityCurve, SolidLiquidPhase},
    state_point::{
        MaterialPropertySelection, SPECIFIC_HEAT_CAPACITY_DIMS, resolve_material_state_point,
    },
    thermal::source::ThermalMaterialPoint,
};
use fs_qty::{Density, QuantitySpec, Temperature};
use fs_rep_mesh::{Metric2, TriComplex2};
use fs_solver::{CheckedCgConfig, NewtonKrylovConfig};
use fs_thermal::enthalpy::*;
use std::sync::Arc;

const MU: f64 = 0.01;
const CHI: f64 = 0.2;
const PIN: f64 = 2.0;
const CP: f64 = 0.1;
const VOLUME: f64 = 0.3;
const DT: f64 = 0.1;
const CELL: CellId = CellId(7);
const MATERIAL: MaterialId = MaterialId(9);

fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t).unwrap()
}
fn domain(upper: f64) -> ValidityDomain {
    ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        200.0,
        upper,
    )
}
fn roles() -> HysteresisSources {
    HysteresisSources {
        background_permeability: PropertyKey::new("background_mu", MAGNETIC_PERMEABILITY_DIMS),
        field_limit: PropertyKey::new("field_limit", FIELD_STRENGTH_DIMS),
        branches: vec![HysteresisBranchSources {
            susceptibility: PropertyKey::new("branch_slope", MAGNETIC_PERMEABILITY_DIMS),
            pinning_field: PropertyKey::new("branch_pinning", FIELD_STRENGTH_DIMS),
            saturation_polarization: PropertyKey::new("branch_saturation", FLUX_DENSITY_DIMS),
        }],
    }
}
fn card(fit_upper: Option<f64>) -> Arc<MaterialCard> {
    let sources = roles();
    let branch = &sources.branches[0];
    let mut claims = ClaimSet::new();
    let mut insert = |key: PropertyKey, value, validity, interpolation| {
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured finite-heat magnetic play".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    for (key, value) in [
        (sources.background_permeability.clone(), MU),
        (sources.field_limit.clone(), 30.0),
        (branch.susceptibility.clone(), CHI),
        (branch.pinning_field.clone(), PIN),
        (branch.saturation_polarization.clone(), 1.0),
        (PropertyKey::new("density", Density::DIMS), 1.0),
        (PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS), CP),
    ] {
        let validity = if matches!(key.name(), "density" | "specific_heat_capacity") {
            domain(500.0)
        } else {
            fit_upper.map(domain).unwrap_or_else(ValidityDomain::unconstrained)
        };
        let dims = key.dims();
        insert(
            key,
            PropertyValue::Scalar { value, dims },
            validity,
            InterpolationPolicy::ConstantWithinValidity,
        );
    }
    let basis = ElectricalPropertyBasis::Resistivity;
    insert(
        PropertyKey::new(basis.property_name(), basis.dimensions()),
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            dims: basis.dimensions(),
            knots: vec![
                (200.0, resistivity(200.0)),
                (300.0, resistivity(300.0)),
                (500.0, resistivity(500.0)),
            ],
        },
        domain(500.0),
        InterpolationPolicy::LinearInside,
    );
    Arc::new(
        MaterialCard::assemble(
            MaterialStateId {
                chemistry: "synthetic magnetic conductor".into(),
                phase: "solid".into(),
                process: "manufactured thermal experiment".into(),
                revision: 0,
            },
            claims,
            vec![],
        )
        .unwrap(),
    )
}
fn resistivity(t: f64) -> f64 {
    (1.0 + 0.005 * (t - 300.0)) / 18.0
}
fn field_controls() -> HystereticEddyControls {
    HystereticEddyControls {
        magnetic: PlanarHystereticControls {
            nonlinear: NonlinearMagneticControls {
                absolute_residual_tolerance_a_m: 1e-9,
                relative_residual_tolerance: 1e-10,
                max_newton_iterations: 50,
                max_backtracks: 48,
                linear: CheckedCgConfig {
                    tolerance: 1e-13,
                    max_iterations: 100,
                    max_corrections: 2,
                },
            },
            material: VectorHysteresisControls {
                relative_tolerance: 1e-12,
                absolute_field_tolerance_a_m: 1e-12,
                absolute_flux_tolerance_t: 1e-13,
                work_tolerance_j_m3: 2e-10,
                maximum_branch_iterations: 100,
                maximum_inverse_iterations: 80,
                maximum_backtracks: 48,
            },
            energy_tolerance_j: 2e-8,
        },
        current_balance_tolerance_a: 1e-8,
    }
}
fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 1,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-9,
        thermal_energy_tolerance_j: 1e-9,
        total_energy_tolerance_j: 2e-8,
        circuit: MaterialCircuitSolve { energy_tolerance_j: 2e-8, power_tolerance_w: 1e-9 },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn problem(card: &Arc<MaterialCard>, mesh: TriComplex2) -> PlanarHystereticMagnetostatics {
    let law = ResolvedVectorHysteresis::from_isotropic_source(
        ResolvedScalarHysteresis::resolve(
            card.clone(),
            &point(300.0),
            roles(),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
    )
    .unwrap();
    let laws = vec![law; mesh.faces().len()];
    PlanarHystereticMagnetostatics::new(
        mesh,
        laws,
        WeightedAssemblyLimits { max_cells: 100, max_dofs: 100, max_triplets: 10000 },
        &mut || false,
    )
    .unwrap()
}
fn retained<'a>(
    card: &MaterialCard,
    owner: &'a PlanarHystereticMagnetostatics,
) -> PlanarHystereticEddy<'a> {
    let boundary = geometry::boundary(0.0);
    let zero = owner.demagnetized(0.0, &mut || false).unwrap();
    let initial = owner
        .candidate(&zero, &[72.0; 4], &boundary, field_controls().magnetic, &mut || false)
        .unwrap()
        .into_parts()
        .0;
    let basis = ElectricalPropertyBasis::Resistivity;
    let conductor = ResolvedConductorState::try_from_material_with_basis(
        resolve_material_state_point(
            card,
            &point(300.0),
            &[basis.requirement().unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        basis,
    )
    .unwrap();
    PlanarHystereticEddy::new(
        &initial,
        (0..4)
            .map(|triangle| EddyConductor { triangle, region: 7, conductor: conductor.clone() })
            .collect(),
        &[],
        &boundary,
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap()
}
fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
    ResolvedHeatCapacityCurve::try_new(
        card,
        HeatCapacitySourceQuery {
            heat_capacity: &PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            density: &PropertyKey::new("density", Density::DIMS),
            temperature_axis: "T",
            lower: &point(200.0),
            upper: &point(500.0),
            selection: SelectionPolicy::SingleClaimOnly,
        },
        SolidLiquidPhase::Solid,
        0.0,
        1e-9,
    )
    .unwrap()
}
fn thermal(curve: &ResolvedHeatCapacityCurve, volume: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 0,
                body_id: 0,
                material_id: MATERIAL,
                mass_kg: volume,
                specific_h_j_kg: CP * (300.0 - 200.0),
                geometry: CellGeometry {
                    start_m: [0.0; 3],
                    end_m: [1.0, 0.0, 0.0],
                    volume_m3: volume,
                    surface_area_m2: 1.0,
                },
                boundary: Boundary {
                    ambient_k: 300.0,
                    htc_w_m2_k: 0.0,
                    reservoir: None,
                    source_w: 0.0,
                },
            }),
        }],
        &[(MATERIAL, curve.curve() as &dyn ThermalLaw)],
        DepositionBudget { max_cells: 1, max_links: 1, max_events: 1, max_substeps: 1000 },
    )
    .unwrap()
    .advance_to(0.0, StepConfig { max_dt_s: DT, energy_tolerance_j: 1e-9 }, || Ok(()))
    .unwrap()
    .0
}
fn cells<'a>(
    card: &'a MaterialCard,
    curve: &'a ResolvedHeatCapacityCurve,
    context: &'a QueryPoint,
) -> [MaterialCell<'a>; 1] {
    [MaterialCell {
        cell: CELL,
        material_id: MATERIAL,
        curve: curve.curve(),
        point: ThermalMaterialPoint {
            material: card,
            context,
            temperature_axis: "T",
            liquid_fraction_axis: None,
        },
        conductivity: None,
    }]
}
fn routes(count: usize) -> Vec<MagneticRegion> {
    (0..count).map(|triangle| MagneticRegion::Cell { triangle, cell: CELL }).collect()
}
fn network<'a>(
    cells: &'a [MaterialCell<'a>],
    routes: &'a [MagneticRegion],
) -> MagneticThermalNetwork<'a> {
    MagneticThermalNetwork {
        cells,
        links: &[],
        magnetic_regions: routes,
        volume_tolerance_m3: 1e-12,
    }
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 5e-8 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
fn play(old_z: f64, h: f64) -> f64 {
    old_z.clamp(CHI * (h - PIN), CHI * (h + PIN)).clamp(-1.0, 1.0)
}
fn storage(h: f64, z: f64) -> f64 {
    0.5 * MU * h * h + z * z / (2.0 * CHI)
}
fn bisect(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
    assert!(residual(lo) <= 0.0 && residual(hi) >= 0.0);
    for _ in 0..80 {
        let mid = f64::midpoint(lo, hi);
        if residual(mid) > 0.0 {
            hi = mid;
        } else {
            lo = mid;
        }
    }
    f64::midpoint(lo, hi)
}
struct Reference {
    t: f64,
    h: f64,
    z: f64,
    b: f64,
    joule: f64,
    pin: f64,
    numerical: f64,
}
// Analytic P1 moments give projected conductor mass sigma*volume/18. The
// nested scalar solves do not call the production field, material or Newton code.
fn reference(t0: f64, old_h: f64, old_z: f64, old_b: f64) -> Reference {
    let at = |t| {
        let sigma = 1.0 / resistivity(t);
        let h =
            bisect(-30.0, 30.0, |h| h + sigma / (72.0 * DT) * (MU * h + play(old_z, h) - old_b));
        let z = play(old_z, h);
        let b = MU * h + z;
        (
            h,
            z,
            b,
            sigma * VOLUME / 18.0 * (0.5 * (b - old_b)).powi(2) / DT,
            VOLUME * PIN * (z - old_z).abs(),
        )
    };
    let t = bisect(t0, t0 + VOLUME * storage(old_h, old_z) / (CP * VOLUME) + 1.0, |t| {
        let (_, _, _, joule, pin) = at(t);
        CP * VOLUME * (t - t0) - joule - pin
    });
    let (h, z, b, joule, pin) = at(t);
    let numerical = VOLUME * (h * (b - old_b) - storage(h, z) + storage(old_h, old_z)) - pin;
    Reference { t, h, z, b, joule, pin, numerical }
}

#[test]
fn finite_temperature_refresh_preserves_play_history_and_deposits_only_physical_heat() {
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes(4);
    let net = network(&cells, &routes);
    let owner = problem(&card, geometry::mesh(VOLUME));
    let mut state = retained(&card, &owner);
    let mut heat = thermal(&curve, VOLUME);
    let (mut old_h, mut old_z, mut old_b) = (12.0, 1.0, 1.12);
    for step in 1..=16 {
        let expected = reference(heat.temperature_k(CELL).unwrap(), old_h, old_z, old_b);
        let before = (
            heat.cells().clone(),
            state.potential_wb_m().to_vec(),
            state.magnetic_state().cell_states().to_vec(),
        );
        let candidate = net
            .step_hysteretic_eddy_candidate(
                &heat,
                step as f64 * DT,
                &state,
                &[0.0; 4],
                controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        near(candidate.thermal().temperature_k(CELL).unwrap(), expected.t);
        near(candidate.field().potential_wb_m()[4], 0.5 * expected.b);
        near(report.electrical.joule_heat_j, expected.joule);
        near(report.electrical.hysteresis_heat_j, expected.pin);
        near(report.electrical.numerical_dissipation_j, expected.numerical);
        near(report.deposited_joule_j, expected.joule + expected.pin);
        near(report.thermal.stored_energy_change_j, expected.joule + expected.pin);
        near(report.total_energy_residual_j, 0.0);
        near(report.non_electrical_heat_j, 0.0);
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(state.potential_wb_m(), before.1);
        assert_eq!(state.magnetic_state().cell_states(), before.2);
        for source in candidate.field().conductors() {
            near(source.conductor.resistivity_ohm_m(), resistivity(expected.t));
            assert_eq!(source.conductor.source_basis(), ElectricalPropertyBasis::Resistivity);
        }
        if step == 1 {
            assert!(expected.joule > 0.01 && expected.pin > 0.01 && expected.numerical > 0.01);
            assert!(
                (candidate.thermal().energy_j()
                    - heat.energy_j()
                    - expected.joule
                    - expected.pin
                    - expected.numerical)
                    .abs()
                    > 0.01
            );
            assert!((1.0 / resistivity(expected.t) - 18.0).abs() > 0.1);
        }
        let (next_heat, next_state, _) = candidate.into_parts();
        heat = next_heat;
        state = next_state;
        (old_h, old_z, old_b) = (expected.h, expected.z, expected.b);
    }
    assert!(old_h.abs() < 1e-6);
    near(old_z, CHI * PIN); // Remanent B survives after macroscopic current has decayed.
    assert!(state.potential_wb_m()[4] > 0.19);
}

fn paired_squares() -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "insulating magnetic pair",
        vec![
            [0.0, 0.0],
            [1.0, 0.0],
            [1.0, 1.0],
            [0.0, 1.0],
            [2.0, 0.0],
            [2.0, 1.0],
            [0.5, 0.5],
            [1.5, 0.5],
        ],
        vec![
            [0, 1, 6],
            [1, 2, 6],
            [2, 3, 6],
            [3, 0, 6],
            [1, 4, 7],
            [4, 5, 7],
            [5, 2, 7],
            [2, 1, 7],
        ],
        Metric2::planar(VOLUME).unwrap(),
    )
    .unwrap()
}
#[test]
fn insulating_hysteresis_receives_pinning_heat_and_rejects_an_unowned_receiver() {
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let owner = problem(&card, paired_squares());
    let zero = owner.demagnetized(0.0, &mut || false).unwrap();
    let state = PlanarHystereticEddy::new(
        &zero,
        vec![],
        &(0..8).collect::<Vec<_>>(),
        &(0..6).map(|v| (v, 0.0)).collect::<Vec<_>>(),
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap();
    let cells = cells(&card, &curve, &context);
    let mut routes = routes(8);
    let heat = thermal(&curve, 2.0 * VOLUME);
    let drive = [72.0, 72.0, 72.0, 72.0, -72.0, -72.0, -72.0, -72.0];
    let next = network(&cells, &routes)
        .step_hysteretic_eddy_candidate(
            &heat,
            DT,
            &state,
            &drive,
            controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    let expected_pin = 2.0 * VOLUME * PIN;
    near(next.thermal().temperature_k(CELL).unwrap(), 300.0 + PIN / CP);
    near(next.report().electrical.joule_heat_j, 0.0);
    near(next.report().electrical.hysteresis_heat_j, expected_pin);
    near(next.report().deposited_joule_j, expected_pin);
    near(next.report().thermal.stored_energy_change_j, expected_pin);
    near(next.report().total_energy_residual_j, 0.0);
    assert!(next.report().electrical.numerical_dissipation_j > 0.1);
    assert!(next.report().electrical.cells.is_empty());
    assert_eq!(next.report().electrical.magnetic.cells.len(), 8);
    routes[0] = MagneticRegion::Fixed { triangle: 0 };
    // Match only mapped volume so refusal specifically reaches the lossy
    // fixed magnetic cell, rather than an unrelated volume mismatch.
    let incomplete_heat = thermal(&curve, 1.75 * VOLUME);
    assert!(matches!(
        network(&cells, &routes).step_hysteretic_eddy_candidate(
            &incomplete_heat,
            DT,
            &state,
            &drive,
            controls(),
            field_controls(),
            &[0.0],
            || false
        ),
        Err(CircuitEnthalpyError::Source(MagneticThermalError::Invalid(
            "every dissipating magnetic triangle needs an actual thermal heat destination"
        )))
    ));
    assert_eq!(state.time_s(), 0.0);
    assert_eq!(heat.time_s(), 0.0);
    assert!(state.potential_wb_m().iter().all(|a| *a == 0.0));
}

#[test]
fn unsupported_temperature_fit_card_and_late_cancellation_publish_neither_history() {
    for upper in [None, Some(300.0), Some(500.0)] {
        let card = card(upper);
        let curve = caloric(&card);
        let context = point(300.0);
        let owner = problem(&card, geometry::mesh(VOLUME));
        let state = retained(&card, &owner);
        let heat = thermal(&curve, VOLUME);
        let cells = cells(&card, &curve, &context);
        let routes = routes(4);
        let net = network(&cells, &routes);
        let before = (
            heat.cells().clone(),
            state.potential_wb_m().to_vec(),
            state.magnetic_state().cell_states().to_vec(),
        );
        let mut polls = 0;
        let run = net.step_hysteretic_eddy_candidate(
            &heat,
            DT,
            &state,
            &[0.0; 4],
            controls(),
            field_controls(),
            &[0.0],
            || {
                polls += 1;
                false
            },
        );
        if upper != Some(500.0) {
            assert!(run.is_err());
        } else {
            let accepted = run.unwrap();
            assert!(polls > 20);
            let mut seen = 0;
            let cancelled = net.step_hysteretic_eddy_candidate(
                &heat,
                DT,
                &state,
                &[0.0; 4],
                controls(),
                field_controls(),
                &[0.0],
                || {
                    seen += 1;
                    seen >= polls
                },
            );
            assert!(matches!(cancelled, Err(CircuitEnthalpyError::Cancelled)));
            let retry = net
                .step_hysteretic_eddy_candidate(
                    &heat,
                    DT,
                    &state,
                    &[0.0; 4],
                    controls(),
                    field_controls(),
                    &[0.0],
                    || false,
                )
                .unwrap();
            assert_eq!(retry.thermal().cells(), accepted.thermal().cells());
            assert_eq!(retry.field().potential_wb_m(), accepted.field().potential_wb_m());
            assert_eq!(
                retry.field().magnetic_state().cell_states(),
                accepted.field().magnetic_state().cell_states()
            );
            let wrong = MaterialCard::assemble(
                MaterialStateId { revision: 1, ..card.id().clone() },
                card.claims().clone(),
                vec![],
            )
            .unwrap();
            let wrong_cells = self::cells(&wrong, &curve, &context);
            assert!(
                network(&wrong_cells, &routes)
                    .step_hysteretic_eddy_candidate(
                        &heat,
                        DT,
                        &state,
                        &[0.0; 4],
                        controls(),
                        field_controls(),
                        &[0.0],
                        || false
                    )
                    .is_err()
            );
        }
        assert_eq!(heat.cells(), &before.0);
        assert_eq!(state.potential_wb_m(), before.1);
        assert_eq!(state.magnetic_state().cell_states(), before.2);
        assert_eq!(heat.time_s(), 0.0);
        assert_eq!(state.time_s(), 0.0);
    }
}

#[test]
fn complete_heat_provider_keeps_the_native_field_and_both_physical_losses() {
    use fs_couple::circuit_enthalpy::{
        contacts::{heat::MaterialHeatNetwork, material::ElectricalContactNetwork},
        eddy::thermal::heat::hysteretic_heat_assembly_candidate,
    };
    use fs_exec::{CancelGate, Cx};
    let card = card(Some(500.0));
    let curve = caloric(&card);
    let context = point(300.0);
    let owner = problem(&card, geometry::mesh(VOLUME));
    let state = retained(&card, &owner);
    let heat = thermal(&curve, VOLUME);
    let cells = cells(&card, &curve, &context);
    let routes = routes(4);
    let bulk = MaterialNetwork { cells: &cells, resistors: &[], links: &[] };
    let electrical = ElectricalContactNetwork { bulk: &bulk, contacts: &[], electrical_bases: &[] };
    let assembly =
        MaterialHeatNetwork { electrical: &electrical, thermal_contacts: &[], radiation: None };
    let gate = CancelGate::new();
    let pool = fs_alloc::ArenaPool::new(fs_alloc::ArenaConfig::default());
    pool.scope(|arena| {
        let cx = Cx::new(
            &gate,
            arena,
            fs_exec::StreamKey { seed: 59, kernel_id: 73, tile: 0, iteration: 0 },
            fs_exec::Budget::INFINITE,
            fs_exec::ExecMode::Deterministic,
        );
        let prepared = assembly.prepare(&cx, &heat).unwrap();
        let next = hysteretic_heat_assembly_candidate(
            &prepared,
            &cx,
            &heat,
            DT,
            &state,
            &[0.0; 4],
            &routes,
            1e-12,
            controls(),
            field_controls(),
            &[0.0],
        )
        .unwrap();
        let expected = reference(300.0, 12.0, 1.0, 1.12);
        near(next.thermal().temperature_k(CELL).unwrap(), expected.t);
        near(next.field().potential_wb_m()[4], 0.5 * expected.b);
        near(next.report().deposited_joule_j, expected.joule + expected.pin);
        near(next.report().total_energy_residual_j, 0.0);
        assert_eq!(next.report().source.network.material.cells.len(), 1);
        gate.request();
        assert!(matches!(
            hysteretic_heat_assembly_candidate(
                &prepared,
                &cx,
                &heat,
                DT,
                &state,
                &[0.0; 4],
                &routes,
                1e-12,
                controls(),
                field_controls(),
                &[0.0]
            ),
            Err(CircuitEnthalpyError::Cancelled)
        ));
        assert_eq!(state.time_s(), 0.0);
        assert_eq!(heat.time_s(), 0.0);
    });
}
