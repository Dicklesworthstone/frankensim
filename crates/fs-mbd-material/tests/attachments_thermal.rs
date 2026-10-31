//! Finite native damping heat and original temperature-dependent interface feedback.
#[path = "heating/common.rs"]
mod thermal_source;
use thermal_source::{budget, caloric, card, close, point, region};
use fs_evidence::ValidityDomain;
use fs_ga::{Se3, Vec3, Wrench};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, PropertyClaim, PropertyKey,
    PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext, UncertaintyModel,
};
use fs_material::phase::{ResolvedHeatCapacityCurve, SolidLiquidPhase};
use fs_material::state_point::{
    MaterialPropertySelection, ScalarAdmissibility, ScalarPropertyRequirement,
    interface_stiffness_per_area_dims, interface_viscosity_per_area_dims,
    resolve_interface_state_point,
};
use fs_material::thermal::interface::ThermalInterfacePoint;
use fs_mbd::articulated::{FreeFloatingBaseState, JointModel};
use fs_mbd::articulated::motion::{FreeFloatingMotionLoads, FreeFloatingMotionState};
use fs_mbd::articulated::motion::spatial_loads::SpatialMotionLoads;
use fs_mbd::articulated::motion::spatial_loads::bushing::MaterialFrame;
use fs_mbd_material::articulated::{MaterialLink, ResolvedArticulated, resolve_articulated};
use fs_mbd_material::attachments::{MaterialBushing, RectangularPad, step_bushings_implicit};
use fs_mbd_material::attachments::thermal::{
    BushingHeatBinding, BushingHeatControls, BushingHeatError, HeatedBushings,
};
use fs_mbd_material::thermal::heating::LinkHeating;
use fs_mbd_material::{RegionBinding, VolumeMoments};
use fs_qty::{QuantitySpec, Temperature};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_thermal::enthalpy::{CellId, DepositionError, StepConfig};

const ZERO: Vec3 = Vec3 { x: 0.0, y: 0.0, z: 0.0 };
const PROPERTIES: [(&str, &str); 3] = [("kx", "cx"), ("ky", "cy"), ("kz", "cz")];
const DT: f64 = 0.02;

fn temperature() -> QuantitySpec {
    QuantitySpec::semantic(SemanticType::new(QuantityKind::AbsoluteTemperature, ValueForm::Static))
}
fn context() -> QueryPoint {
    QueryPoint::new().with("T", 300.).unwrap().with("omega", 0.).unwrap()
        .with_quantity("TA", temperature(), 300.).unwrap()
        .with_quantity("TB", temperature(), 300.).unwrap()
}
fn mechanism(a: &MaterialCard, b: &MaterialCard) -> ResolvedArticulated {
    let p = point(300., 1e5);
    let regions = [a, b].map(|material| [RegionBinding {
        key: "surface", material, point: &p,
        geometry: VolumeMoments::cuboid(Vec3::new(1., 1., 1.)).unwrap(),
        selection: MaterialPropertySelection::SingleClaimOnly,
    }]);
    resolve_articulated(&[
        MaterialLink {
            name: "a", parent: None, parent_from_child_zero: Se3::identity(),
            joint: JointModel::FIXED, regions: &regions[0],
        },
        MaterialLink {
            name: "b", parent: Some(0), parent_from_child_zero: Se3::identity(),
            joint: JointModel::prismatic(Vec3::new(1., 0., 0.), None).unwrap(),
            regions: &regions[1],
        },
    ], 2, || false).unwrap()
}
fn interface(a: &MaterialCard, b: &MaterialCard) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    for (axis, (k, c)) in PROPERTIES.into_iter().enumerate() {
        for (name, dims, value, scale) in [
            (k, interface_stiffness_per_area_dims(), [100., 200., 300.][axis], 2.),
            (c, interface_viscosity_per_area_dims(), [5., 7., 11.][axis], 3.),
        ] {
            claims.insert_claim(PropertyClaim {
                key: PropertyKey::new(name, dims),
                value: PropertyValue::Curve {
                    abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                    knots: vec![(300., value), (400., value * scale)], dims,
                },
                validity: ValidityDomain::unconstrained().with("T", 300., 400.)
                    .with("omega", 0., 1e4)
                    .with_quantity("TA", temperature(), 300., 400.)
                    .with_quantity("TB", temperature(), 300., 400.),
                interpolation: InterpolationPolicy::LinearInside,
                uncertainty: UncertaintyModel::Unstated, observations: vec![],
                provenance: Provenance {
                    source: "synthetic thermal bonded pad, not measurements".into(),
                    license: "CC0-1.0".into(), artifact: None,
                },
            }).unwrap();
        }
    }
    InterfaceSystemCard::assemble(
        SurfaceSpec { material: a.id().clone(), texture_frame: "pad axes".into() },
        SurfaceSpec { material: b.id().clone(), texture_frame: "pad axes".into() },
        SystemContext {
            medium: "fixture".into(), third_body: None,
            environment: "declared paired and mean temperatures".into(), history: "undegraded".into(),
        }, claims, vec![],
    ).unwrap()
}
fn binding<'a>(
    model: &'a ResolvedArticulated, card: &'a InterfaceSystemCard,
    materials: [&'a MaterialCard; 2], context: &'a QueryPoint,
) -> BushingHeatBinding<'a> {
    let requirements: Vec<_> = PROPERTIES.into_iter().flat_map(|(k, c)| [
        ScalarPropertyRequirement::try_new(k, interface_stiffness_per_area_dims(),
            ScalarAdmissibility::NonNegative).unwrap(),
        ScalarPropertyRequirement::try_new(c, interface_viscosity_per_area_dims(),
            ScalarAdmissibility::NonNegative).unwrap(),
    ]).collect();
    let source = resolve_interface_state_point(
        card, context, &requirements, MaterialPropertySelection::SingleClaimOnly,
    ).unwrap();
    BushingHeatBinding {
        bushing: MaterialBushing::rectangular_pad(
            model, ["surface"; 2], RectangularPad {
                frames: [0, 1].map(|link| MaterialFrame {
                    link, link_from_attachment: Se3::identity(),
                }),
                rest_translation_a: ZERO, size_m: [0.2, 0.3],
                maximum_displacement_m: 0.1, maximum_rotation_rad: 0.05,
            }, source, PROPERTIES, || false,
        ).unwrap(),
        point: ThermalInterfacePoint {
            interface: card, materials, context, temperature_axes: ["TA", "TB"],
            liquid_fraction_axes: [None; 2], mean_temperature_axis: Some("T"),
        },
        cells: [CellId(0), CellId(1)], fraction_to_a: 0.25,
    }
}
fn links<'a>(materials: [&'a MaterialCard; 2], curves: [&'a ResolvedHeatCapacityCurve; 2])
    -> Vec<LinkHeating<'a>> {
    ["a", "b"].into_iter().enumerate().map(|(i, link)| LinkHeating {
        link, regions: vec![region(materials[i], curves[i], "surface", i as u64,
            if i == 0 { 100. } else { 0. })],
    }).collect()
}
fn state() -> FreeFloatingMotionState {
    FreeFloatingMotionState {
        base: FreeFloatingBaseState::stationary(Se3::identity()),
        joint_position: vec![0.01], joint_velocity: vec![0.04],
    }
}
fn loads(wrenches: &[Wrench]) -> SpatialMotionLoads<'_> {
    SpatialMotionLoads {
        held: FreeFloatingMotionLoads {
            generalized_force: &[0.], external_body_wrench: wrenches, gravity_world: ZERO,
        }, world_points: &[], connections: &[],
    }
}
fn config() -> StepConfig { StepConfig { max_dt_s: DT, energy_tolerance_j: 1e-9 } }
fn controls() -> BushingHeatControls {
    BushingHeatControls {
        mechanics: Default::default(), maximum_parameter_work_j: 1., energy_tolerance_j: 1e-8,
    }
}

#[test]
fn actual_midpoint_loss_heats_original_cells_then_refreshes_every_interface_axis() {
    let a = card(5., 0., 500.);
    let b = card(10., 0., 500.);
    let curves = [&a, &b].map(|m| caloric(m, 1e5, SolidLiquidPhase::Solid));
    let model = mechanism(&a, &b);
    let interface = interface(&a, &b);
    let context = context();
    let binding = binding(&model, &interface, [&a, &b], &context);
    let wrenches = [Wrench::new(ZERO, ZERO); 2];
    let native = step_bushings_implicit(
        &[&binding.bushing], &state(), loads(&wrenches), DT, controls().mechanics, || false,
    ).unwrap();
    let initial = HeatedBushings::new(
        &model, links([&a, &b], [&curves[0], &curves[1]]), vec![binding],
        state(), vec![], budget(), config(), || Ok(()),
    ).unwrap();
    let (next, report) = initial.advance_to(
        DT, loads(&wrenches), config(), controls(), || Ok(()),
    ).unwrap();

    // Independent midpoint solution: each rigid body has mass 2, hence reduced
    // mass 1. Native x coefficients are area .06 times k=100 and c=5.
    let k = 6.;
    let c = 0.3;
    let h = DT / 2.;
    let g = DT * DT * k / 4.;
    let v = ((1. - h * c - g) * 0.04 - DT * k * 0.01) / (1. + h * c + g);
    let q = 0.01 + DT * (0.04 + v) / 2.;
    let heat = c * DT * ((v + 0.04) / 2.).powi(2);
    assert!((next.state().joint_velocity[0] - v).abs() < 2e-9);
    assert!((next.state().joint_position[0] - q).abs() < 2e-11);
    assert_eq!(next.state(), &native.spatial.motion.state_after);
    close(report.dissipated_heat_j, heat);
    close(report.thermal.external_heat_j, 2. + heat);
    close(report.external_heat_j, 2.);
    // Heat capacities are independently 2*5 and 2*10 J/K; neither side
    // receives parameter work, a solver defect or its loss fraction twice.
    let ta = 300. + (2. + 0.25 * heat) / 10.;
    let tb = 300. + 0.75 * heat / 20.;
    close(next.thermal().cells()[&CellId(0)].specific_h_j_kg, (2. + 0.25 * heat) / 2.);
    close(next.thermal().cells()[&CellId(1)].specific_h_j_kg, 0.75 * heat / 2.);
    let mean = (ta + tb) / 2.;
    for (i, axis) in next.bindings()[0].bushing.axes().iter().enumerate() {
        close(axis.resolved().state_coordinate("TA").unwrap(), ta);
        close(axis.resolved().state_coordinate("TB").unwrap(), tb);
        close(axis.resolved().state_coordinate("T").unwrap(), mean);
        close(axis.stiffness_n_m(), [6., 12., 18.][i] * (1. + (mean - 300.) / 100.));
        close(axis.damping_n_s_m(), [0.3, 0.42, 0.66][i] * (1. + 2. * (mean - 300.) / 100.));
    }
    let parameter = 0.5 * k * (mean - 300.) / 100. * q * q;
    assert!(parameter > 1e-8);
    close(report.elastic_parameter_work_j, parameter);
    close(report.absolute_parameter_work_j, parameter);
    assert!(report.mechanical_work_residual_j.abs() < 1e-8);
    assert!(report.total_energy_residual_j.abs() < 1e-8);
    assert_eq!(initial.state(), &state());
    assert_eq!(initial.thermal().time_s(), 0.);
    assert_eq!(initial.thermal().cells()[&CellId(0)].specific_h_j_kg, 0.);

    // Continuing motion consumes the newly sourced stiffness AND viscosity.
    let expected = step_bushings_implicit(
        &[&next.bindings()[0].bushing], next.state(), loads(&wrenches),
        DT, controls().mechanics, || false,
    ).unwrap();
    let (continued, _) = next.advance_to(
        2. * DT, loads(&wrenches), config(), controls(), || Ok(()),
    ).unwrap();
    assert_eq!(continued.state(), &expected.spatial.motion.state_after);
    let (whole, _) = initial.advance_to(
        2. * DT, loads(&wrenches), config(), controls(), || Ok(()),
    ).unwrap();
    assert_eq!(whole.state(), continued.state());
    close(whole.thermal().cells()[&CellId(0)].specific_h_j_kg,
        continued.thermal().cells()[&CellId(0)].specific_h_j_kg);
}

#[test]
fn explicit_cell_identity_work_budget_and_late_cancellation_preserve_both_histories() {
    let a = card(5., 0., 500.);
    let b = card(10., 0., 500.);
    let curves = [&a, &b].map(|m| caloric(m, 1e5, SolidLiquidPhase::Solid));
    let model = mechanism(&a, &b);
    let interface = interface(&a, &b);
    let context = context();
    let binding = binding(&model, &interface, [&a, &b], &context);
    let mut wrong = binding.clone();
    wrong.cells.reverse();
    assert!(HeatedBushings::new(
        &model, links([&a, &b], [&curves[0], &curves[1]]), vec![wrong],
        state(), vec![], budget(), config(), || Ok(()),
    ).is_err());
    let initial = HeatedBushings::new(
        &model, links([&a, &b], [&curves[0], &curves[1]]), vec![binding],
        state(), vec![], budget(), config(), || Ok(()),
    ).unwrap();
    let wrenches = [Wrench::new(ZERO, ZERO); 2];
    let zero_work = BushingHeatControls { maximum_parameter_work_j: 0., ..controls() };
    assert!(initial.advance_to(
        DT, loads(&wrenches), config(), zero_work, || Ok(()),
    ).is_err());

    let mut checkpoints = 0;
    initial.advance_to(DT, loads(&wrenches), config(), controls(), || {
        checkpoints += 1;
        Ok(())
    }).unwrap();
    let mut visited = 0;
    let result = initial.advance_to(DT, loads(&wrenches), config(), controls(), || {
        visited += 1;
        if visited == checkpoints { Err(DepositionError::Cancelled) } else { Ok(()) }
    });
    assert!(matches!(result, Err(BushingHeatError::Thermal(DepositionError::Cancelled))));
    assert_eq!(initial.state(), &state());
    assert_eq!(initial.thermal().time_s(), 0.);
    for cell in initial.thermal().cells().values() {
        assert_eq!(cell.specific_h_j_kg, 0.);
    }
    close(initial.bindings()[0].bushing.axes()[0].stiffness_n_m(), 6.);
}
