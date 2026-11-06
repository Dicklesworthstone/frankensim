//! Shared real source, vector field and finite caloric owners for induction tests.
//! Coefficients are manufactured; these fixtures make no catalog-validity claim.
#![allow(dead_code)]

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

pub(super) const MU: f64 = 0.01;
pub(super) const CHI: f64 = 0.2;
pub(super) const PIN: f64 = 2.0;
pub(super) const CP: f64 = 0.1;
pub(super) const VOLUME: f64 = 0.3;
pub(super) const DT: f64 = 0.1;
pub(super) const CELL: CellId = CellId(7);
pub(super) const MATERIAL: MaterialId = MaterialId(9);
pub(super) const DIRECTIONS: [[f64; 2]; 4] = [[1.0, 0.0], [0.0, 1.0], [-1.0, 0.0], [0.0, -1.0]];

pub(super) fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", QuantitySpec::dimensional(Temperature::DIMS), t)
        .unwrap()
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
pub(super) fn card(fit_upper: Option<f64>) -> Arc<MaterialCard> {
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
        (
            PropertyKey::new("specific_heat_capacity", SPECIFIC_HEAT_CAPACITY_DIMS),
            CP,
        ),
    ] {
        let validity = if matches!(key.name(), "density" | "specific_heat_capacity") {
            domain(500.0)
        } else {
            fit_upper
                .map(domain)
                .unwrap_or_else(ValidityDomain::unconstrained)
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
pub(super) fn resistivity(t: f64) -> f64 {
    (1.0 + 0.005 * (t - 300.0)) / 18.0
}
pub(super) fn field_controls() -> HystereticEddyControls {
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
pub(super) fn controls() -> CircuitEnthalpySolve {
    CircuitEnthalpySolve {
        maximum_cells: 1,
        maximum_circuit_trials: 1000,
        maximum_newton_iterations: 24,
        difference_step_j: 1e-4,
        cell_energy_tolerance_j: 1e-9,
        thermal_energy_tolerance_j: 1e-9,
        total_energy_tolerance_j: 2e-8,
        circuit: MaterialCircuitSolve {
            energy_tolerance_j: 2e-8,
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
pub(super) fn problem(
    card: &Arc<MaterialCard>,
    mesh: TriComplex2,
) -> PlanarHystereticMagnetostatics {
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
        WeightedAssemblyLimits {
            max_cells: 100,
            max_dofs: 100,
            max_triplets: 10000,
        },
        &mut || false,
    )
    .unwrap()
}
pub(super) fn conductor(card: &MaterialCard, temperature: f64) -> ResolvedConductorState {
    let basis = ElectricalPropertyBasis::Resistivity;
    ResolvedConductorState::try_from_material_with_basis(
        resolve_material_state_point(
            card,
            &point(temperature),
            &[basis.requirement().unwrap()],
            MaterialPropertySelection::SingleClaimOnly,
        )
        .unwrap(),
        basis,
    )
    .unwrap()
}
pub(super) fn retained<'a>(
    card: &MaterialCard,
    owner: &'a PlanarHystereticMagnetostatics,
) -> PlanarHystereticEddy<'a> {
    let zero = owner.demagnetized(0.0, &mut || false).unwrap();
    let initial = owner
        .candidate(
            &zero,
            &[72.0; 4],
            &boundary(4),
            field_controls().magnetic,
            &mut || false,
        )
        .unwrap()
        .into_parts()
        .0;
    let conductor = conductor(card, 300.0);
    PlanarHystereticEddy::new(
        &initial,
        (0..4)
            .map(|triangle| EddyConductor {
                triangle,
                region: 7,
                conductor: conductor.clone(),
            })
            .collect(),
        &[],
        &boundary(4),
        0.0,
        field_controls(),
        &mut || false,
    )
    .unwrap()
}
pub(super) fn caloric(card: &MaterialCard) -> ResolvedHeatCapacityCurve {
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
pub(super) fn thermal(curve: &ResolvedHeatCapacityCurve, volume: f64) -> DepositionSimulation<'_> {
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
            max_dt_s: DT,
            energy_tolerance_j: 1e-9,
        },
        || Ok(()),
    )
    .unwrap()
    .0
}
pub(super) fn cells<'a>(
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
pub(super) fn routes(count: usize) -> Vec<MagneticRegion> {
    (0..count)
        .map(|triangle| MagneticRegion::Cell {
            triangle,
            cell: CELL,
        })
        .collect()
}
pub(super) fn network<'a>(
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
pub(super) fn square() -> TriComplex2 {
    TriComplex2::from_indexed_triangles(
        "hysteretic conductor square",
        vec![[0.0, 0.0], [1.0, 0.0], [1.0, 1.0], [0.0, 1.0], [0.5, 0.5]],
        vec![[0, 1, 4], [1, 2, 4], [2, 3, 4], [3, 0, 4]],
        Metric2::planar(VOLUME).unwrap(),
    )
    .unwrap()
}
pub(super) fn paired_squares() -> TriComplex2 {
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
pub(super) fn boundary(count: usize) -> Vec<(usize, f64)> {
    (0..count).map(|vertex| (vertex, 0.0)).collect()
}
pub(super) fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 5e-8 * (1.0 + b.abs()), "{a:e} != {b:e}");
}
pub(super) fn play(old_z: f64, h: f64) -> f64 {
    old_z
        .clamp(CHI * (h - PIN), CHI * (h + PIN))
        .clamp(-1.0, 1.0)
}
pub(super) fn storage(h: f64, z: f64) -> f64 {
    0.5 * MU * h * h + z * z / (2.0 * CHI)
}
pub(super) fn bisect(mut lo: f64, mut hi: f64, residual: impl Fn(f64) -> f64) -> f64 {
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
