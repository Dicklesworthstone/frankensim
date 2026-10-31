//! Native frame thermocreep: sourced expansion, creep shape change and reciprocal heat.
//! Manufactured coefficients are numerical references, not material qualification.

use fs_evidence::ValidityDomain;
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::creep::source::{
    CREEP_ACTIVATION_ENERGY_DIMS, CREEP_ACTIVATION_ENERGY_PROPERTY, CREEP_RATE_DIMS,
    CREEP_REFERENCE_RATE_PROPERTY, CREEP_REFERENCE_STRESS_PROPERTY,
    CREEP_REFERENCE_TEMPERATURE_PROPERTY, CREEP_STRESS_EXPONENT_PROPERTY, NortonCreepSource,
    NortonSourceQuery,
};
use fs_material::creep::thermal::expansion::ExpandingCreepSource;
use fs_material::creep::CreepTolerances;
use fs_material::state_point::{
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY, MaterialPropertySelection,
    POISSON_RATIO_PROPERTY, ScalarAdmissibility, ScalarPropertyRequirement,
    YOUNG_MODULUS_PROPERTY, resolve_material_state_point,
};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Dims, Pressure, QuantitySpec, Temperature};
use fs_solid::fiber::{Fiber, FiberLaw, Section};
use fs_solid::frame2d::creep::FrameCreepTolerances;
use fs_solid::{Frame2d, Geometry, SectionModel};

const E: f64 = 1.0e9;
const NU: f64 = 0.25;
const STRESS: f64 = 1.0e6;
const RATE: f64 = 1.0e-5;
const DT: f64 = 10.0;
const ALPHA: f64 = 2.0e-5;
const AREA: f64 = 5.0e-4;
const ARM: f64 = 0.05;
const FORCE: f64 = 2.0 * AREA * STRESS;

fn temperature_spec() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(
        QuantityKind::AbsoluteTemperature,
        ValueForm::Static,
    ))
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity("T", temperature_spec(), t)
        .unwrap()
}
fn claim(name: &str, dims: Dims, value: f64) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with_quantity("T", temperature_spec(), 300.0, 500.0),
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        uncertainty: UncertaintyModel::Unstated,
        observations: vec![],
        provenance: Provenance {
            source: "manufactured frame thermocreep fixture".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}
fn card() -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (name, dims, value) in [
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, E),
        (POISSON_RATIO_PROPERTY, Dims::NONE, NU),
        (CREEP_REFERENCE_STRESS_PROPERTY, Pressure::DIMS, STRESS),
        (CREEP_REFERENCE_RATE_PROPERTY, CREEP_RATE_DIMS, RATE),
        (CREEP_STRESS_EXPONENT_PROPERTY, Dims::NONE, 1.0),
        (
            CREEP_ACTIVATION_ENERGY_PROPERTY,
            CREEP_ACTIVATION_ENERGY_DIMS,
            0.0,
        ),
        (
            CREEP_REFERENCE_TEMPERATURE_PROPERTY,
            Temperature::DIMS,
            300.0,
        ),
        (
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ALPHA,
        ),
    ] {
        claims.insert_claim(claim(name, dims, value)).unwrap();
    }
    MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic thermal frame".into(),
            phase: "solid".into(),
            process: "manufactured".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap()
}
fn source(card: &MaterialCard) -> ExpandingCreepSource<NortonCreepSource<'_>> {
    let mechanical = NortonCreepSource::new(
        card,
        NortonSourceQuery {
            lower: &point(300.0),
            upper: &point(500.0),
            temperature_axis: "T",
            selection: MaterialPropertySelection::SingleClaimOnly,
            strain_limit: 0.03,
        },
        || false,
    )
    .unwrap();
    let expansion = resolve_material_state_point(
        card,
        &point(300.0),
        &[ScalarPropertyRequirement::try_new(
            LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY,
            Dims([0, 0, 0, -1, 0, 0]),
            ScalarAdmissibility::StrictlyPositive,
        )
        .unwrap()],
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    ExpandingCreepSource::new(card, mechanical, &expansion, || false).unwrap()
}
fn frame(
    law: &fs_material::creep::thermal::uniaxial::UniaxialThermalNorton,
    axial_only: bool,
) -> Frame2d {
    let mut frame = Frame2d::new(Geometry::Linear);
    frame.add_node(0.0, 0.0);
    frame.add_node(1.0, 0.0);
    frame.fix(0, [true; 3]);
    if axial_only {
        frame.fix(1, [false, true, true]);
    }
    frame
        .add_element(0, 1, 3, 1.0, &|| {
            SectionModel::Fiber(Section::new(
                [-ARM, ARM]
                    .into_iter()
                    .map(|y| Fiber {
                        y,
                        area: AREA,
                        law: FiberLaw::ThermalNorton(law.clone()),
                    })
                    .collect(),
            ))
        })
        .unwrap();
    frame
}
fn controls() -> FrameCreepTolerances {
    FrameCreepTolerances {
        material: CreepTolerances {
            stress_pa: 1e-3,
            energy_j_m3: 1e-5,
        },
        energy_j: 1e-7,
        force_n: 1e-5,
        moment_n_m: 1e-6,
        maximum_fiber_points: 6,
    }
}
fn axial_load(force: f64) -> [f64; 6] {
    [0.0, 0.0, 0.0, force, 0.0, 0.0]
}
fn close(actual: f64, expected: f64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "actual={actual:.16e}, expected={expected:.16e}, tolerance={tolerance:.3e}"
    );
}

#[test]
fn heated_loaded_member_expands_creeps_and_reports_reciprocal_heat() {
    let card = card();
    let source = source(&card);
    let cold = source.at_uniaxial_temperature(300.0, DT).unwrap();
    let accepted = frame(&cold, true);
    let hot = source.at_uniaxial_temperature(350.0, DT).unwrap();
    let prepared = accepted
        .with_creep_laws(&vec![hot; 6], "T", || false)
        .unwrap();

    let step = prepared
        .solve_creep_step(&axial_load(FORCE), controls(), || false)
        .unwrap();
    let free = ALPHA * 50.0;
    close(
        step.frame.displacements()[3],
        free + STRESS / E + DT * RATE,
        3e-11,
    );
    close(step.frame.element_forces(0)[0], FORCE, 1e-6);
    assert!(step.energy.physical_dissipation_j > 0.0);
    assert!(step.energy.thermoelastic_storage_change_j > 0.0);
    assert!(step.energy.reversible_heat_j < 0.0);
    close(
        step.energy.caloric_heat_j(),
        step.energy.physical_dissipation_j + step.energy.reversible_heat_j,
        0.0,
    );
    close(step.energy.balance_residual_j, 0.0, controls().energy_j);
    for section in 0..3 {
        for fiber in 0..2 {
            let (_, state) = step.frame.thermal_creep_fiber(0, section, fiber).unwrap();
            close(state.mechanical.creep_strain(), DT * RATE, 2e-11);
            assert_eq!(state.temperature_k(), 350.0);
        }
    }
}

#[test]
fn through_depth_temperature_gradient_bends_the_actual_cantilever_without_fake_stress() {
    let card = card();
    let source = source(&card);
    let cold = source.at_uniaxial_temperature(300.0, DT).unwrap();
    let accepted = frame(&cold, false);
    let bottom = source.at_uniaxial_temperature(300.0, DT).unwrap();
    let top = source.at_uniaxial_temperature(350.0, DT).unwrap();
    let laws: Vec<_> = (0..3)
        .flat_map(|_| [bottom.clone(), top.clone()])
        .collect();
    let prepared = accepted.with_creep_laws(&laws, "T", || false).unwrap();
    let step = prepared
        .solve_creep_step(&[0.0; 6], controls(), || false)
        .unwrap();

    let top_free = ALPHA * 50.0;
    let axial = 0.5 * top_free;
    let curvature = -top_free / (2.0 * ARM);
    for &(actual_axial, actual_curvature) in step.frame.section_deformations(0) {
        close(actual_axial, axial, 2e-11);
        close(actual_curvature, curvature, 2e-10);
    }
    close(step.frame.displacements()[4], 0.5 * curvature, 3e-10);
    close(step.frame.displacements()[5], curvature, 3e-10);
    close(step.frame.element_forces(0)[0], 0.0, 1e-7);
    close(step.frame.element_forces(0)[1], 0.0, 1e-7);
    close(step.frame.element_forces(0)[2], 0.0, 1e-7);
    close(step.energy.physical_dissipation_j, 0.0, 1e-10);
    close(step.energy.reversible_heat_j, 0.0, 1e-9);
    close(step.energy.balance_residual_j, 0.0, controls().energy_j);
}
