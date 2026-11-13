//! Full material-graph integration; the frozen source supplies every liquid unit.
use super::*;
use fs_evidence::ValidityDomain;
use fs_lbm::d3q19::freesurface3::ContactModel3;
use fs_lbm::material::free_surface::{
    MaterialFreeSurface3, MaterialFreeSurfaceDomain3, plan_free_surface_from_material,
};
use fs_matdb::{
    ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, UncertaintyModel,
};
use fs_material::liquid::surface::{
    LiquidCapillarity, ResolvedLiquidSurfaceState, resolve_liquid_surface_state,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection};
use fs_mbd::{Pose, UnitQuaternion};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, QuantitySpec};

fn source() -> ResolvedLiquidSurfaceState {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (DENSITY_PROPERTY, 1000.0, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, 0.1, DYNAMIC_VISCOSITY_DIMS),
    ] {
        claims
            .insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Scalar { value, dims },
                validity: ValidityDomain::unconstrained(),
                uncertainty: UncertaintyModel::Unstated,
                interpolation: InterpolationPolicy::ConstantWithinValidity,
                observations: vec![],
                provenance: Provenance {
                    source: "synthetic reciprocal momentum fixture; not measurements".into(),
                    license: "test-only".into(),
                    artifact: None,
                },
            })
            .unwrap();
    }
    let card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-liquid".into(),
            phase: "liquid".into(),
            process: "drag-source-fixture".into(),
            revision: 0,
        },
        claims,
        vec![],
    )
    .unwrap();
    let query = QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            300.0,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap();
    resolve_liquid_surface_state(
        &card,
        &query,
        LiquidViscosityBasis::Dynamic,
        LiquidCapillarity::Neglected,
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}
fn controls() -> LiquidRigidDragControls {
    LiquidRigidDragControls {
        maximum_cells: 64,
        maximum_stencil_cells: 1,
        maximum_speed_m_s: 0.1,
        maximum_momentum_defect_kg_m_s: 1e-18,
        maximum_angular_momentum_defect_kg_m2_s: 1e-21,
        maximum_work_defect_j: 1e-19,
    }
}
fn fixture<'a>(
    source: &'a ResolvedLiquidSurfaceState,
) -> (MaterialFreeSurface3<'a>, RigidBodyState, MassProperties) {
    let liquid = plan_free_surface_from_material(source, 1e-3, 2e-4, 0.1)
        .unwrap()
        .free_surface3(
            MaterialFreeSurfaceDomain3 {
                cells: [4, 4, 4],
                gravity_m_s2: [0.0; 3],
                contact: ContactModel3::Neutral,
                solid_tiles: &[],
                max_cells: 64,
            },
            |x, _, _| x < 3,
        )
        .unwrap();
    let mass = 3.0e-6;
    let properties = MassProperties::new(mass, Vec3::ZERO, Vec3::new(3e-12, 3e-12, 3e-12)).unwrap();
    let rigid = RigidBodyState::new(
        Pose::new(Vec3::new(2.5e-3, 0.5e-3, 0.5e-3), UnitQuaternion::IDENTITY).unwrap(),
        Vec3::new(mass * 0.04, 0.0, 0.0),
        Vec3::ZERO,
    )
    .unwrap();
    (liquid, rigid, properties)
}

#[test]
fn material_bound_source_uses_original_density_scale_and_partial_mass() {
    let source = source();
    let (mut liquid, mut rigid, properties) = fixture(&source);
    let original = rigid;
    let mass_unit = liquid.scaling().mass_unit_kg();
    let sample = liquid.lattice().sample_at(2, 0, 0).unwrap();
    let wet_mass = mass_unit * sample.density * sample.fill_fraction;
    let report = apply_rigid_drag(
        &mut liquid,
        &mut rigid,
        properties,
        &[2],
        Vec3::ZERO,
        3.0,
        0.2,
        controls(),
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let m = properties.mass();
    let expected = m * wet_mass / (m + wet_mass) * 0.04 * -(-0.6 * (1.0 + wet_mass / m)).exp_m1();
    assert!((report.liquid_impulse_kg_m_s.x - expected).abs() < 2e-22);
    assert_eq!(report.units.mass_unit_kg, mass_unit);
    assert_eq!(report.units.cell_size_m, 1e-3);
    assert_eq!(
        report.units.velocity_unit_m_s,
        liquid.scaling().velocity_si(1.0).unwrap()
    );
    assert!(core::ptr::eq(liquid.scaling().material(), &source));
    assert_eq!(liquid.lattice().grid().steps(), 0);
    assert_eq!(rigid.pose(), original.pose());
    assert!(report.drag_dissipation_j > 0.0);
}
#[test]
fn original_material_scenario_speed_is_admitted_before_either_owner_changes() {
    let source = source();
    let (mut liquid, mut rigid, properties) = fixture(&source);
    let initial = rigid;
    let bits = liquid.lattice().grid().state_bits();
    let mut invalid = controls();
    invalid.maximum_speed_m_s = 0.100_001;
    assert!(matches!(
        apply_rigid_drag(
            &mut liquid,
            &mut rigid,
            properties,
            &[2],
            Vec3::ZERO,
            3.0,
            0.2,
            invalid,
            &CancelGate::new_clock_free()
        ),
        Err(MaterialLiquidMomentumError::Fluid(_))
    ));
    assert_eq!(initial, rigid);
    assert_eq!(bits, liquid.lattice().grid().state_bits());
}

#[test]
fn directional_material_source_retains_original_units_and_partial_liquid_inertia() {
    let source = source();
    let (mut liquid, mut rigid, properties) = fixture(&source);
    let original = rigid;
    let mass_unit = liquid.scaling().mass_unit_kg();
    let sample = liquid.lattice().sample_at(2, 0, 0).unwrap();
    let wet_mass = mass_unit * sample.density * sample.fill_fraction;
    let axes = UnitQuaternion::from_axis_angle(Vec3::new(0.0, 0.0, 1.0), 0.4).unwrap();
    let law = heterogeneous::PrincipalLiquidDrag3 {
        rates_per_s: [3.0, 0.0, 0.0],
        axes_to_world: axes,
    };
    let limits = heterogeneous::HeterogeneousRigidDragControls {
        exchange: controls(),
        substeps: 1,
        maximum_pair_exchanges: 1,
    };
    let report = apply_heterogeneous_rigid_drag(
        &mut liquid,
        &mut rigid,
        properties,
        &[2],
        &[law],
        Vec3::ZERO,
        0.2,
        limits,
        &CancelGate::new_clock_free(),
    )
    .unwrap();
    let axis = axes.rotate_body_to_world(Vec3::new(1.0, 0.0, 0.0));
    let m = properties.mass();
    let expected = axis.scale(
        m * wet_mass / (m + wet_mass) * (0.04 * axis.x) * -(-0.6 * (1.0 + wet_mass / m)).exp_m1(),
    );
    for (actual, expected) in a(report.exchange.liquid_impulse_kg_m_s)
        .into_iter()
        .zip(a(expected))
    {
        assert!((actual - expected).abs() < 4e-22);
    }
    assert_eq!(report.exchange.units.mass_unit_kg, mass_unit);
    assert_eq!(report.exchange.units.cell_size_m, 1e-3);
    assert_eq!(
        report.exchange.units.velocity_unit_m_s,
        liquid.scaling().velocity_si(1.0).unwrap()
    );
    assert!(core::ptr::eq(liquid.scaling().material(), &source));
    assert_eq!(liquid.lattice().grid().steps(), 0);
    assert_eq!(rigid.pose(), original.pose());
    assert!(report.exchange.drag_dissipation_j > 0.0);
}

#[path = "tests/evolving.rs"]
mod evolving;
