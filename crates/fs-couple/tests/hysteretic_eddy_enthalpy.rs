//! G1/G3/G4: native vector memory, spatial losses and finite enthalpy together.
//! Manufactured coefficients test the equations and accounting, not real materials.
#![cfg(feature = "hysteretic-eddy-heating")]

use fs_couple::{
    circuit_enthalpy::{
        CircuitEnthalpySolve,
        magnetic::{MagneticRegion, MagneticThermalNetwork},
        material::MaterialCell,
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

#[path = "hysteretic_eddy_enthalpy/rate.rs"]
mod rate;

const DEPTH: f64 = 0.3;
const MU: f64 = 0.01;
const CHI: f64 = 0.2;
const HC: f64 = 2.0;
const DENSITY: f64 = 2.0;
const CP: f64 = 3.0;
const CELL: CellId = CellId(7);
const MATERIAL: MaterialId = MaterialId(9);

fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
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
fn card(thermal_fit: bool) -> MaterialCard {
    let domain = ValidityDomain::unconstrained().with_quantity(
        "T",
        QuantitySpec::dimensional(Temperature::DIMS),
        200.0,
        500.0,
    );
    let mut claims = ClaimSet::new();
    let mut insert = |key: PropertyKey, value: PropertyValue, validity| {
        let interpolation = if matches!(value, PropertyValue::Scalar { .. }) {
            InterpolationPolicy::ConstantWithinValidity
        } else {
            InterpolationPolicy::LinearInside
        };
        claims
            .insert_claim(PropertyClaim {
                key,
                value,
                validity,
                interpolation,
                uncertainty: UncertaintyModel::Unstated,
                observations: vec![],
                provenance: Provenance {
                    source: "manufactured magnetic heat fixture".into(),
                    license: "CC0-1.0".into(),
                    artifact: None,
                },
            })
            .unwrap();
    };
    let r = roles();
    for (key, value) in [
        (r.background_permeability, MU),
        (r.field_limit, 30.0),
        (r.branches[0].susceptibility.clone(), CHI),
        (r.branches[0].pinning_field.clone(), HC),
        (r.branches[0].saturation_polarization.clone(), 1.0),
    ] {
        let dims = key.dims();
        insert(
            key,
            PropertyValue::Scalar { value, dims },
            if thermal_fit {
                domain.clone()
            } else {
                ValidityDomain::unconstrained()
            },
        );
    }
    for (name, dims, value) in [
        ("density", Density::DIMS, DENSITY),
        ("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS, CP),
    ] {
        insert(
            PropertyKey::new(name, dims),
            PropertyValue::Scalar { value, dims },
            domain.clone(),
        );
    }
    let basis = ElectricalPropertyBasis::Resistivity;
    insert(
        PropertyKey::new(basis.property_name(), basis.dimensions()),
        PropertyValue::Curve {
            abscissa: "T".into(),
            abscissa_dims: Temperature::DIMS,
            knots: vec![
                (200.0, 1.0 / 18.0),
                (300.0, 1.0 / 18.0),
                (500.0, 11.0 / 18.0),
            ],
            dims: basis.dimensions(),
        },
        domain,
    );
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic magnetic conductor".into(),
            phase: "solid".into(),
            process: "manufactured fit".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
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
        cell_energy_tolerance_j: 1e-8,
        thermal_energy_tolerance_j: 1e-8,
        total_energy_tolerance_j: 5e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 1e-9,
            power_tolerance_w: 1e-9,
        },
        newton: NewtonKrylovConfig {
            absolute_tolerance: 1e-10,
            relative_tolerance: 1e-14,
            forcing_minimum: 1e-4,
            linear_restart: 8,
            ..NewtonKrylovConfig::default()
        },
    }
}
fn problem(card: &MaterialCard, paired: bool) -> PlanarHystereticMagnetostatics {
    problem_with_viscosity(card, paired, None)
}
fn problem_with_viscosity(
    card: &MaterialCard,
    paired: bool,
    viscosity: Option<PropertyKey>,
) -> PlanarHystereticMagnetostatics {
    let (vertices, faces) = if paired {
        (
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
        )
    } else {
        (
            vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
            vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        )
    };
    let source = ResolvedScalarHysteresis::resolve(
        Arc::new(card.clone()),
        &point(300.0),
        roles(),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    let source = match viscosity {
        Some(key) => ResolvedVectorHysteresis::from_isotropic_source_with_viscosity(
            source,
            vec![key],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        None => ResolvedVectorHysteresis::from_isotropic_source(source).unwrap(),
    };
    let materials = vec![source; faces.len()];
    PlanarHystereticMagnetostatics::new(
        TriComplex2::from_indexed_triangles(
            "hysteresis-enthalpy",
            vertices,
            faces,
            Metric2::planar(DEPTH).unwrap(),
        )
        .unwrap(),
        materials,
        WeightedAssemblyLimits {
            max_cells: 100,
            max_dofs: 100,
            max_triplets: 10000,
        },
        &mut || false,
    )
    .unwrap()
}
fn state<'m>(
    card: &MaterialCard,
    problem: &'m PlanarHystereticMagnetostatics,
    conducting: bool,
) -> PlanarHystereticEddy<'m> {
    let count = problem.mesh().faces().len();
    let boundary: Vec<_> = (0..if count == 8 { 6 } else { 4 })
        .map(|v| (v, 0.0))
        .collect();
    let mut initial = problem.demagnetized(0.0, &mut || false).unwrap();
    let mut conductors = vec![];
    let mut insulators = vec![];
    if conducting {
        initial = problem
            .candidate(
                &initial,
                &[72.0; 4],
                &boundary,
                field_controls().magnetic,
                &mut || false,
            )
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
        conductors = (0..count)
            .map(|triangle| EddyConductor {
                triangle,
                region: 7,
                conductor: conductor.clone(),
            })
            .collect();
    } else {
        insulators = (0..count).collect();
    }
    PlanarHystereticEddy::new(
        &initial,
        conductors,
        &insulators,
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
fn thermal_state(curve: &ResolvedHeatCapacityCurve, volume: f64) -> DepositionSimulation<'_> {
    DepositionSimulation::new(
        vec![TimedEvent {
            time_s: 0.0,
            event: Event::Birth(DepositedCell {
                id: CELL,
                segment_id: 0,
                body_id: 0,
                material_id: MATERIAL,
                mass_kg: DENSITY * volume,
                specific_h_j_kg: CP * 100.0,
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
        DepositionBudget {
            max_cells: 1,
            max_links: 1,
            max_events: 1,
            max_substeps: 1000,
        },
    )
    .unwrap()
    .advance_to(
        0.0,
        StepConfig {
            max_dt_s: 1.0,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
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
    (0..count)
        .map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: CELL,
        })
        .collect()
}
fn near(a: f64, b: f64) {
    assert!(
        (a - b).abs() <= 5e-8 * (1.0 + b.abs()),
        "{a:.16e} != {b:.16e}"
    );
}
fn play(old: f64, h: f64) -> f64 {
    old.clamp(CHI * (h - HC), CHI * (h + HC)).clamp(-1.0, 1.0)
}

#[test]
fn insulating_minor_loops_heat_finite_cells_without_joule_or_numerical_heat() {
    let card = card(true);
    let curve = caloric(&card);
    let mut thermal = thermal_state(&curve, 2.0 * DEPTH);
    let problem = problem(&card, true);
    let mut field = state(&card, &problem, false);
    let context = point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes(8);
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let mut z = 0.0;
    let mut expected_temperature = 300.0;
    for h in [12.0, 0.0, 6.0, 1.0, 5.5, 1.0, 6.0, 0.0, -12.0, 0.0] {
        let next_z = play(z, h);
        let heat = 2.0 * DEPTH * HC * (next_z - z).abs();
        expected_temperature += heat / (DENSITY * CP * 2.0 * DEPTH);
        let mut currents = [6.0 * h; 8];
        currents[4..].fill(-6.0 * h);
        let candidate = network
            .step_hysteretic_eddy_candidate(
                &thermal,
                field.time_s() + 0.1,
                &field,
                &currents,
                controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        near(candidate.report().electrical.hysteresis_heat_j, heat);
        assert_eq!(candidate.report().electrical.joule_heat_j, 0.0);
        near(candidate.deposited_heat_j(), heat);
        near(candidate.report().thermal.stored_energy_change_j, heat);
        near(
            candidate.thermal().temperature_k(CELL).unwrap(),
            expected_temperature,
        );
        near(candidate.report().total_energy_residual_j, 0.0);
        if z == 0.0 && h == 12.0 {
            assert!(candidate.report().electrical.numerical_dissipation_j > heat);
        }
        near(
            candidate.field().magnetic_state().cell_states()[0].polarizations_t()[0][0],
            next_z,
        );
        for sample in &candidate.report().source.magnetic {
            near(
                sample
                    .material
                    .query_point()
                    .iter()
                    .find(|(axis, _)| axis == "T")
                    .unwrap()
                    .1,
                expected_temperature,
            );
            assert_eq!(sample.material.properties().len(), 5);
        }
        (thermal, field, _) = candidate.into_parts();
        z = next_z;
    }
}

// Independent scalar play return and one-dimensional spatial balance. The
// floating single-square conductor has C44=sigma*depth/18 and B=2*A_center.
fn reference(old_t: f64, old_z: f64, old_b: f64, dt: f64) -> [f64; 6] {
    let response = |t: f64| {
        let sigma = 18.0 / (1.0 + 0.05 * (t - 300.0));
        let (mut lo, mut hi) = (-30.0, 30.0);
        for _ in 0..90 {
            let h = 0.5 * (lo + hi);
            let b = MU * h + play(old_z, h);
            if h + sigma * (b - old_b) / (72.0 * dt) < 0.0 {
                lo = h;
            } else {
                hi = h;
            }
        }
        let h = 0.5 * (lo + hi);
        let z = play(old_z, h);
        let b = MU * h + z;
        let pinning = DEPTH * HC * (z - old_z).abs();
        let joule = sigma * DEPTH / 18.0 * (0.5 * (b - old_b)).powi(2) / dt;
        [h, z, b, pinning, joule]
    };
    let (mut lo, mut hi) = (old_t, 500.0);
    for _ in 0..90 {
        let t = 0.5 * (lo + hi);
        let x = response(t);
        if t - old_t - (x[3] + x[4]) / (DENSITY * CP * DEPTH) < 0.0 {
            lo = t;
        } else {
            hi = t;
        }
    }
    let t = 0.5 * (lo + hi);
    let x = response(t);
    [t, x[0], x[1], x[2], x[3], x[4]]
}

#[test]
fn conducting_hysteresis_and_endpoint_temperature_feedback_match_independent_balance() {
    let card = card(true);
    let curve = caloric(&card);
    let mut thermal = thermal_state(&curve, DEPTH);
    let problem = problem(&card, false);
    let mut field = state(&card, &problem, true);
    let context = point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes(4);
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let (mut old_t, mut old_z, mut old_b) = (300.0, 1.0, 1.12);
    for dt in [0.1, 1.0, 3.0] {
        let expected = reference(old_t, old_z, old_b, dt);
        let memory = field.magnetic_state().cell_states().to_vec();
        let cold = field
            .step_to_candidate(
                field.time_s() + dt,
                field.conductors().to_vec(),
                &[0.0; 4],
                field_controls(),
                &mut || false,
            )
            .unwrap();
        let candidate = network
            .step_hysteretic_eddy_candidate(
                &thermal,
                field.time_s() + dt,
                &field,
                &[0.0; 4],
                controls(),
                field_controls(),
                &[0.0],
                || false,
            )
            .unwrap();
        let report = candidate.report();
        near(
            candidate.thermal().temperature_k(CELL).unwrap(),
            expected[0],
        );
        near(candidate.field().potential_wb_m()[4], 0.5 * expected[3]);
        near(report.electrical.hysteresis_heat_j, expected[4]);
        near(report.electrical.joule_heat_j, expected[5]);
        near(
            report.thermal.stored_energy_change_j,
            expected[4] + expected[5],
        );
        near(candidate.deposited_heat_j(), expected[4] + expected[5]);
        near(report.total_energy_residual_j, 0.0);
        assert!(report.electrical.numerical_dissipation_j > 0.0);
        assert!(report.circuit_trials > 1);
        assert!(
            (cold.state().potential_wb_m()[4] - candidate.field().potential_wb_m()[4]).abs() > 1e-9
        );
        assert_eq!(field.magnetic_state().cell_states(), memory);
        let rho = candidate.field().conductors()[0]
            .conductor
            .material()
            .property(ElectricalPropertyBasis::Resistivity.property_name())
            .unwrap()
            .value_si();
        near(rho, (1.0 + 0.05 * (expected[0] - 300.0)) / 18.0);
        (thermal, field, _) = candidate.into_parts();
        (old_t, old_z, old_b) = (expected[0], expected[2], expected[3]);
    }
}

#[test]
fn physical_heat_requires_receivers_fixed_geometry_and_source_temperature_support() {
    for supported in [false, true] {
        let card = card(supported);
        let curve = caloric(&card);
        let thermal = thermal_state(&curve, 2.0 * DEPTH);
        let problem = problem(&card, true);
        let field = state(&card, &problem, false);
        let context = point(300.0);
        let cells = cells(&card, &curve, &context);
        let mapped = routes(8);
        let fixed: Vec<_> = (0..8)
            .map(|triangle| MagneticRegion::Fixed { triangle })
            .collect();
        let routes = if supported { &fixed } else { &mapped };
        let network = MagneticThermalNetwork {
            cells: &cells,
            links: &[],
            magnetic_regions: routes,
            volume_tolerance_m3: 1e-12,
        };
        let mut currents = [72.0; 8];
        currents[4..].fill(-72.0);
        assert!(
            network
                .step_hysteretic_eddy_candidate(
                    &thermal,
                    0.1,
                    &field,
                    &currents,
                    controls(),
                    field_controls(),
                    &[0.0],
                    || false,
                )
                .is_err()
        );
        near(thermal.temperature_k(CELL).unwrap(), 300.0);
        assert_eq!(field.energy_j(), 0.0);
        if supported {
            let mismatched = thermal_state(&curve, 1.01 * 2.0 * DEPTH);
            let network = MagneticThermalNetwork {
                cells: &cells,
                links: &[],
                magnetic_regions: &mapped,
                volume_tolerance_m3: 1e-12,
            };
            assert!(
                network
                    .step_hysteretic_eddy_candidate(
                        &mismatched,
                        0.1,
                        &field,
                        &currents,
                        controls(),
                        field_controls(),
                        &[0.0],
                        || false,
                    )
                    .is_err()
            );
        }
    }
}

#[test]
fn cancelled_and_exhausted_trials_preserve_both_accepted_histories() {
    let card = card(true);
    let curve = caloric(&card);
    let thermal = thermal_state(&curve, DEPTH);
    let problem = problem(&card, false);
    let field = state(&card, &problem, true);
    let context = point(300.0);
    let cells = cells(&card, &curve, &context);
    let routes = routes(4);
    let network = MagneticThermalNetwork {
        cells: &cells,
        links: &[],
        magnetic_regions: &routes,
        volume_tolerance_m3: 1e-12,
    };
    let memory = field.magnetic_state().cell_states().to_vec();
    let potential = field.potential_wb_m().to_vec();
    let enthalpy = thermal.cells()[&CELL].specific_h_j_kg.to_bits();
    let mut calls = 0;
    let accepted = network
        .step_hysteretic_eddy_candidate(
            &thermal,
            0.1,
            &field,
            &[0.0; 4],
            controls(),
            field_controls(),
            &[0.0],
            || {
                calls += 1;
                false
            },
        )
        .unwrap();
    let mut checkpoints = 0;
    assert!(
        network
            .step_hysteretic_eddy_candidate(
                &thermal,
                0.1,
                &field,
                &[0.0; 4],
                controls(),
                field_controls(),
                &[0.0],
                || {
                    checkpoints += 1;
                    checkpoints == calls / 2
                },
            )
            .is_err()
    );
    let mut limited = controls();
    limited.maximum_circuit_trials = 1;
    assert!(
        network
            .step_hysteretic_eddy_candidate(
                &thermal,
                0.1,
                &field,
                &[0.0; 4],
                limited,
                field_controls(),
                &[0.0],
                || false,
            )
            .is_err()
    );
    let mut field_limited = field_controls();
    field_limited.magnetic.nonlinear.max_newton_iterations = 0;
    assert!(
        network
            .step_hysteretic_eddy_candidate(
                &thermal,
                0.1,
                &field,
                &[0.0; 4],
                controls(),
                field_limited,
                &[0.0],
                || false,
            )
            .is_err()
    );
    assert_eq!(field.magnetic_state().cell_states(), memory);
    assert_eq!(field.potential_wb_m(), potential);
    assert_eq!(thermal.cells()[&CELL].specific_h_j_kg.to_bits(), enthalpy);
    let repeated = network
        .step_hysteretic_eddy_candidate(
            &thermal,
            0.1,
            &field,
            &[0.0; 4],
            controls(),
            field_controls(),
            &[0.0],
            || false,
        )
        .unwrap();
    assert_eq!(
        accepted.field().magnetic_state().cell_states(),
        repeated.field().magnetic_state().cell_states()
    );
    assert_eq!(
        accepted.thermal().cells()[&CELL].specific_h_j_kg.to_bits(),
        repeated.thermal().cells()[&CELL].specific_h_j_kg.to_bits()
    );
}
