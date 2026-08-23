//! G3 end-to-end PCB laminate handoff into the steady conduction solver.
//!
//! This is a model/QoI comparison, not physical validation. It proves that
//! receipt-backed anisotropy and its propagated principal bounds reach the
//! existing conduction path without being collapsed to isotropic FR4.

mod support;

use fs_conduction::fixtures::unit_cube;
use fs_conduction::material::{ConductivityModel, ProvenanceClass};
use fs_conduction::solve::element_heat_flux;
use fs_conduction::{
    ConductionMesh, ConductionProblem, InitialGuess, LinearConfig, Nonlinearity, ScalarField,
    SolveConfig, StopRule, ThermalBc, ThermalBoundaryBuilder, solve,
};
use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, CopperCoverage, InterpolationPolicy, MaterialCard, MaterialStateId,
    PCB_THERMAL_CONDUCTIVITY_DIMS, PcbConductivityDatum, PcbLayer, PcbPrincipalFrame,
    PcbScaleSeparation, PcbStackup, PropertyClaim, PropertyKey, PropertyValue, Provenance,
    QueryPoint, SelectionPolicy, UncertaintyModel,
};

use support::with_cx;

const PROPERTY: &str = "thermal_conductivity";

fn material_card(chemistry: &str, process: &str, conductivity: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    claims
        .insert_claim(PropertyClaim {
            key: PropertyKey::new(PROPERTY, PCB_THERMAL_CONDUCTIVITY_DIMS),
            value: PropertyValue::Scalar {
                value: conductivity,
                dims: PCB_THERMAL_CONDUCTIVITY_DIMS,
            },
            validity: ValidityDomain::unconstrained().with("T", 250.0, 400.0),
            uncertainty: UncertaintyModel::HalfWidth {
                half_width: 0.0,
                confidence: 0.95,
            },
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            observations: Vec::new(),
            provenance: Provenance {
                source: format!("{chemistry} PCB conduction fixture"),
                license: "test-only".to_string(),
                artifact: None,
            },
        })
        .expect("claim");
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: chemistry.to_string(),
            phase: "solid".to_string(),
            process: process.to_string(),
            revision: 0,
        },
        claims,
        Vec::new(),
    )
    .expect("card")
}

fn datum(card: &MaterialCard) -> PcbConductivityDatum {
    let point = QueryPoint::new().with("T", 300.0).expect("query");
    PcbConductivityDatum::from_card(card, PROPERTY, &point, SelectionPolicy::SingleClaimOnly)
        .expect("datum")
}

fn coverage(source_id: &str, nominal: f64, lower: f64, upper: f64) -> CopperCoverage {
    CopperCoverage::new(
        source_id,
        nominal,
        lower,
        upper,
        Provenance {
            source: format!("test design-data export {source_id}"),
            license: "test-only".to_string(),
            artifact: None,
        },
    )
    .expect("coverage")
}

fn homogenized_board() -> fs_matdb::PcbHomogenizedConductivity {
    let copper = material_card("C11000", "rolled-foil", 400.0);
    let fr4 = material_card("FR4", "cured-laminate", 0.25);
    let plane = PcbLayer::new(
        "L1-plane",
        0.2e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/L1", 0.95, 0.90, 1.0),
    )
    .expect("plane");
    let core = PcbLayer::new(
        "core",
        0.8e-3,
        datum(&copper),
        datum(&fr4),
        coverage("coverage/core", 0.05, 0.0, 0.10),
    )
    .expect("core");
    PcbStackup::new(
        "reference-pcb",
        vec![plane, core],
        PcbPrincipalFrame::default(),
        PcbScaleSeparation::new(25.0e-6, 0.05).expect("separation"),
    )
    .expect("stackup")
    .homogenize()
    .expect("homogenize")
}

fn solve_config() -> SolveConfig {
    SolveConfig {
        nonlinearity: Nonlinearity::default(),
        stop: StopRule {
            residual_rtol: 1.0e-12,
            residual_atol: 1.0e-24,
            step_atol: 0.0,
            max_iterations: 20,
        },
        linear: LinearConfig {
            tolerance: 1.0e-13,
            max_iterations: 20_000,
            restart: 40,
        },
        initial: InitialGuess::Uniform(305.0),
    }
}

fn solve_board(material: &ConductivityModel) -> fs_conduction::ConductionSolution {
    let (complex, positions) = unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let boundary = ThermalBoundaryBuilder::new(&mesh)
        .region(
            "cold-x",
            |face| face.centroid[0].abs() <= f64::EPSILON,
            ThermalBc::dirichlet(300.0).expect("cold"),
        )
        .expect("cold region")
        .region(
            "hot-x",
            |face| (face.centroid[0] - 1.0).abs() <= f64::EPSILON,
            ThermalBc::dirichlet(310.0).expect("hot"),
        )
        .expect("hot region")
        .adiabatic_remainder()
        .finish()
        .expect("boundary");
    let source = ScalarField::Uniform(0.0);
    with_cx(|cx| {
        solve(
            cx,
            ConductionProblem { element_materials: None, mesh: &mesh,
            boundary: &boundary,
            material,
            source: &source, },
            solve_config(),
        )
        .expect("board solve")
    })
}

fn mean_x_flux(material: &ConductivityModel, temperature: &[f64]) -> f64 {
    let (complex, positions) = unit_cube(2);
    let mesh = ConductionMesh::new(complex, positions).expect("mesh");
    let flux = element_heat_flux(&mesh, material, temperature).expect("flux");
    flux.iter().map(|value| value[0].abs()).sum::<f64>() / flux.len() as f64
}

fn diagonal(principal: [f64; 3]) -> [[f64; 3]; 3] {
    [
        [principal[0], 0.0, 0.0],
        [0.0, principal[1], 0.0],
        [0.0, 0.0, principal[2]],
    ]
}

#[test]
fn pcb_anisotropy_and_bounds_reach_the_board_solve_with_receipts() {
    let homogenized = homogenized_board();
    let nominal =
        ConductivityModel::from_pcb_homogenization(&homogenized).expect("receipt-backed model");
    let lower = ConductivityModel::constant_tensor(diagonal(homogenized.principal().lower_w_mk))
        .expect("lower tensor");
    let upper = ConductivityModel::constant_tensor(diagonal(homogenized.principal().upper_w_mk))
        .expect("upper tensor");
    let isotropic_fr4 = ConductivityModel::isotropic_declared(0.25).expect("isotropic");

    let nominal_solution = solve_board(&nominal);
    let lower_solution = solve_board(&lower);
    let upper_solution = solve_board(&upper);
    let isotropic_solution = solve_board(&isotropic_fr4);
    let nominal_qoi = mean_x_flux(&nominal, &nominal_solution.temperature);
    let lower_qoi = mean_x_flux(&lower, &lower_solution.temperature);
    let upper_qoi = mean_x_flux(&upper, &upper_solution.temperature);
    let isotropic_qoi = mean_x_flux(&isotropic_fr4, &isotropic_solution.temperature);

    assert_eq!(
        nominal_solution.report.material_provenance,
        ProvenanceClass::MatdbReceipts
    );
    assert_eq!(
        nominal_solution.report.material_receipts,
        homogenized.material_uses().len()
    );
    assert!(
        lower_qoi <= nominal_qoi && nominal_qoi <= upper_qoi,
        "propagated board QoI [{lower_qoi}, {upper_qoi}] must bracket {nominal_qoi}"
    );
    assert!(
        nominal_qoi / isotropic_qoi > 100.0,
        "the reference copper stack should exhibit the intended order-of-magnitude in-plane \
         effect: nominal={nominal_qoi}, isotropic={isotropic_qoi}"
    );
    let closure_relative_to_heat_flux =
        nominal_solution.report.energy.closure_w.abs() / nominal_qoi.max(f64::MIN_POSITIVE);
    assert!(
        closure_relative_to_heat_flux < 1.0e-9,
        "pure-Dirichlet energy closure {} W is too large relative to the independently recovered \
         heat-flux scale {nominal_qoi} W/m^2",
        nominal_solution.report.energy.closure_w
    );
    println!(
        "{{\"suite\":\"fs-conduction-pcb\",\"case\":\"board-qoi\",\"status\":\"pass\",\
         \"stackup_identity\":\"{}\",\"in_plane_k\":{},\"through_k\":{},\
         \"qoi_lower_w_m2\":{},\"qoi_nominal_w_m2\":{},\"qoi_upper_w_m2\":{},\
         \"isotropic_qoi_w_m2\":{},\"receipts\":{},\"via_correction\":\"not-modeled\",\
         \"claim\":\"model-comparison-not-physical-validation\"}}",
        homogenized.identity().to_hex(),
        homogenized.principal().nominal_w_mk[0],
        homogenized.principal().nominal_w_mk[2],
        lower_qoi,
        nominal_qoi,
        upper_qoi,
        isotropic_qoi,
        nominal_solution.report.material_receipts
    );
}
