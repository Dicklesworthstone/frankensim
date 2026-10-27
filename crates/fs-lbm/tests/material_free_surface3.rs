//! Sourced equilibrium wetting reaches the actual sparse capillary operator.
//! Synthetic cards establish coupling, identity and geometry behavior, without
//! claiming calibrated sessile-drop or dynamic contact-line accuracy.

use core::f64::consts::{FRAC_PI_2, PI};
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::d3q19::freesurface3::{ContactModel3, FreeSurface3, FreeSurfaceError3, PlanarWall3};
use fs_lbm::material::MaterialLbmError;
use fs_lbm::material::free_surface::{
    MaterialFreeSurfaceDomain3, MaterialFreeSurfaceScaling, plan_free_surface_from_material,
};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::liquid::surface::{
    EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, LiquidCapillarity, ResolvedLiquidSurfaceState,
    ResolvedLiquidWettingState, SURFACE_TENSION_DIMS, SURFACE_TENSION_PROPERTY,
    resolve_liquid_surface_state, resolve_liquid_wetting_state,
};
use fs_material::liquid::{DYNAMIC_VISCOSITY_DIMS, DYNAMIC_VISCOSITY_PROPERTY};
use fs_material::liquid_temperature::LiquidViscosityBasis;
use fs_material::state_point::{DENSITY_PROPERTY, MaterialPropertySelection};
use fs_qty::semantic::{QuantityKind, SemanticType, ValueForm};
use fs_qty::{Density, Dims, QuantitySpec};

#[path = "material_free_surface3/wetting.rs"]
mod dynamic_wetting;

fn point(temperature: f64) -> QueryPoint {
    QueryPoint::new()
        .with_quantity(
            "T",
            QuantitySpec::semantic(SemanticType::new(
                QuantityKind::AbsoluteTemperature,
                ValueForm::Static,
            )),
            temperature,
        )
        .unwrap()
        .with_quantity(
            "pressure",
            QuantitySpec::semantic(SemanticType::new(QuantityKind::Pressure, ValueForm::Static)),
            100_000.0,
        )
        .unwrap()
}

fn claim(name: &str, value: f64, dims: Dims) -> PropertyClaim {
    PropertyClaim {
        key: PropertyKey::new(name, dims),
        value: PropertyValue::Scalar { value, dims },
        validity: ValidityDomain::unconstrained()
            .with_quantity(
                "T",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::AbsoluteTemperature,
                    ValueForm::Static,
                )),
                290.0,
                310.0,
            )
            .with_quantity(
                "pressure",
                QuantitySpec::semantic(SemanticType::new(
                    QuantityKind::Pressure,
                    ValueForm::Static,
                )),
                90_000.0,
                110_000.0,
            ),
        uncertainty: UncertaintyModel::Unstated,
        interpolation: InterpolationPolicy::ConstantWithinValidity,
        observations: vec![],
        provenance: Provenance {
            source: "synthetic planar wetting integration fixture; not measurements".into(),
            license: "test-only".into(),
            artifact: None,
        },
    }
}

fn liquid(revision: u32, temperature: f64, sigma: Option<f64>) -> ResolvedLiquidSurfaceState {
    liquid_at_pressure(revision, temperature, sigma, 100_000.0)
}

fn liquid_at_pressure(
    revision: u32,
    temperature: f64,
    sigma: Option<f64>,
    pressure: f64,
) -> ResolvedLiquidSurfaceState {
    let mut claims = ClaimSet::new();
    for (name, value, dims) in [
        (DENSITY_PROPERTY, 1000.0, Density::DIMS),
        (DYNAMIC_VISCOSITY_PROPERTY, 0.1, DYNAMIC_VISCOSITY_DIMS),
    ] {
        claims.insert_claim(claim(name, value, dims)).unwrap();
    }
    if let Some(value) = sigma {
        claims
            .insert_claim(claim(SURFACE_TENSION_PROPERTY, value, SURFACE_TENSION_DIMS))
            .unwrap();
    }
    let mut card = MaterialCard::assemble(
        MaterialStateId {
            chemistry: "synthetic-liquid".into(),
            phase: "liquid".into(),
            process: "numerical-fixture".into(),
            revision: 0,
        },
        claims.clone(),
        vec![],
    )
    .unwrap();
    for _ in 0..revision {
        card = MaterialCard::supersede(&card, claims.clone(), vec![]).unwrap();
    }
    resolve_liquid_surface_state(
        &card,
        &point(temperature).with("pressure", pressure).unwrap(),
        LiquidViscosityBasis::Dynamic,
        if sigma.is_some() {
            LiquidCapillarity::Sourced
        } else {
            LiquidCapillarity::Neglected
        },
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn interface_card(
    liquid: &ResolvedLiquidSurfaceState,
    angle: Option<f64>,
    angle_dims: Dims,
    reversed: bool,
    dynamic_axis: Option<&str>,
) -> InterfaceSystemCard {
    let solid = SurfaceSpec {
        material: MaterialStateId {
            chemistry: "synthetic-substrate".into(),
            phase: "solid".into(),
            process: "numerical-fixture".into(),
            revision: 0,
        },
        texture_frame: "named-planar-finish".into(),
    };
    let fluid = SurfaceSpec {
        material: liquid.material().material().clone(),
        texture_frame: "liquid-free-surface".into(),
    };
    let mut claims = ClaimSet::new();
    if let Some(angle) = angle {
        let mut source = claim(EQUILIBRIUM_CONTACT_ANGLE_PROPERTY, angle, angle_dims);
        if let Some(axis) = dynamic_axis {
            source.validity = source.validity.with(axis, 0.0, 10.0);
        }
        claims.insert_claim(source).unwrap();
    }
    let (first, second) = if reversed {
        (fluid, solid)
    } else {
        (solid, fluid)
    };
    InterfaceSystemCard::assemble(
        first,
        second,
        SystemContext {
            medium: "synthetic-liquid".into(),
            third_body: None,
            environment: "named-synthetic-gas".into(),
            history: "equilibrium-coupon".into(),
        },
        claims,
        vec![],
    )
    .unwrap()
}

fn wetting(liquid: &ResolvedLiquidSurfaceState, angle: f64) -> ResolvedLiquidWettingState {
    resolve_liquid_wetting_state(
        &interface_card(liquid, Some(angle), Dims([0; 6]), false, None),
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap()
}

fn plan(liquid: &ResolvedLiquidSurfaceState) -> MaterialFreeSurfaceScaling<'_> {
    plan_free_surface_from_material(liquid, 1e-3, 2e-4, 0.1).unwrap()
}

fn domain() -> MaterialFreeSurfaceDomain3<'static> {
    MaterialFreeSurfaceDomain3 {
        cells: [16, 16, 12],
        gravity_m_s2: [0.0; 3],
        contact: ContactModel3::Neutral,
        solid_tiles: &[],
        max_cells: 16 * 16 * 12,
    }
}

fn drop_shape(x: usize, y: usize, z: usize) -> bool {
    let x = x as f64 + 0.5 - 8.0;
    let y = y as f64 + 0.5 - 8.0;
    let z = z as f64 + 0.5;
    x * x + y * y + z * z < 25.0
}

#[test]
fn g1_ninety_degrees_reuses_the_shared_material_owner_and_is_exactly_neutral() {
    let liquid = liquid(0, 300.0, Some(0.05));
    let surface = wetting(&liquid, FRAC_PI_2);
    let scaling = plan(&liquid);
    let mut raw = FreeSurface3::new(
        16,
        16,
        12,
        scaling.flow().tau,
        [0.0; 3],
        scaling.surface_tension_lattice(),
        ContactModel3::Neutral,
        drop_shape,
    )
    .unwrap();
    let mut bound = scaling
        .with_planar_wetting(PlanarWall3::ZMin, &surface)
        .unwrap()
        .free_surface3(domain(), drop_shape)
        .unwrap();
    assert_eq!(bound.scaling().material(), &liquid);
    assert_eq!(
        bound.scaling().planar_wetting(),
        Some((PlanarWall3::ZMin, &surface))
    );
    assert_eq!(
        surface.resolved().properties().len(),
        1,
        "only theta belongs to the wetting interface"
    );
    assert_eq!(surface.context().environment, "named-synthetic-gas");
    assert_eq!(
        bound.lattice().planar_contact_angle(),
        Some((PlanarWall3::ZMin, FRAC_PI_2))
    );
    for _ in 0..12 {
        bound.step().unwrap();
        raw.step().unwrap();
    }
    assert_eq!(bound.lattice().grid().state_bits(), raw.grid().state_bits());
    assert_eq!(
        bound.lattice().ledger_mass().to_bits(),
        raw.ledger_mass().to_bits()
    );
}

#[test]
fn g3_source_angle_changes_actual_capillary_motion_and_is_inert_when_capillarity_is_neglected() {
    let tension = liquid(0, 300.0, Some(0.05));
    let neglected = liquid(0, 300.0, None);
    let surfaces = [
        wetting(&tension, PI / 3.0),
        wetting(&tension, 2.0 * PI / 3.0),
        wetting(&neglected, PI / 3.0),
        wetting(&neglected, 2.0 * PI / 3.0),
    ];
    let mut runs: Vec<_> = [&tension, &tension, &neglected, &neglected]
        .into_iter()
        .zip(&surfaces)
        .map(|(bulk, surface)| {
            plan(bulk)
                .with_planar_wetting(PlanarWall3::ZMin, surface)
                .unwrap()
                .free_surface3(domain(), drop_shape)
                .unwrap()
        })
        .collect();
    let initial = runs[0].ledger_mass_kg().unwrap();
    for _ in 0..12 {
        for run in &mut runs {
            run.step().unwrap();
        }
    }
    assert_ne!(
        runs[0].lattice().grid().state_bits(),
        runs[1].lattice().grid().state_bits(),
        "changing only sourced angle must reach actual pressure reconstruction and flow"
    );
    assert_eq!(
        runs[2].lattice().grid().state_bits(),
        runs[3].lattice().grid().state_bits(),
        "the wall angle must not inject force or mass when sigma is explicitly omitted"
    );
    for run in &runs {
        assert!((run.ledger_mass_kg().unwrap() - initial).abs() / initial < 1e-10);
    }
}

#[test]
fn g0_equilibrium_source_units_state_assignment_and_dynamic_support_refuse() {
    let bulk = liquid(0, 300.0, Some(0.05));
    for card in [
        interface_card(&bulk, None, Dims([0; 6]), false, None),
        interface_card(&bulk, Some(FRAC_PI_2), Density::DIMS, false, None),
        interface_card(&bulk, Some(0.0), Dims([0; 6]), false, None),
        interface_card(&bulk, Some(PI), Dims([0; 6]), false, None),
    ] {
        assert!(
            resolve_liquid_wetting_state(
                &card,
                &point(300.0),
                MaterialPropertySelection::SingleClaimOnly
            )
            .is_err()
        );
    }
    let surface = wetting(&bulk, FRAC_PI_2);
    assert!(
        plan(&bulk)
            .with_planar_wetting(PlanarWall3::ZMin, &surface)
            .unwrap()
            .free_surface3(
                MaterialFreeSurfaceDomain3 {
                    contact: ContactModel3::Wetting,
                    ..domain()
                },
                |_, _, _| panic!("conflicting wall laws must refuse before geometry"),
            )
            .is_err()
    );
    for other in [liquid(1, 300.0, Some(0.05)), liquid(0, 301.0, Some(0.05))] {
        assert!(
            plan(&other)
                .with_planar_wetting(PlanarWall3::ZMin, &surface)
                .is_err()
        );
    }
    let reversed = resolve_liquid_wetting_state(
        &interface_card(&bulk, Some(FRAC_PI_2), Dims([0; 6]), true, None),
        &point(300.0),
        MaterialPropertySelection::SingleClaimOnly,
    )
    .unwrap();
    assert!(
        plan(&bulk)
            .with_planar_wetting(PlanarWall3::ZMin, &reversed)
            .is_err()
    );
    assert!(
        resolve_liquid_wetting_state(
            &interface_card(&bulk, Some(FRAC_PI_2), Dims([0; 6]), false, None),
            &point(320.0),
            MaterialPropertySelection::SingleClaimOnly,
        )
        .is_err()
    );
    for axis in ["shear_rate", "contact_line_speed", "capillary_number"] {
        assert!(
            resolve_liquid_wetting_state(
                &interface_card(&bulk, Some(FRAC_PI_2), Dims([0; 6]), false, Some(axis)),
                &point(300.0).with(axis, 1.0).unwrap(),
                MaterialPropertySelection::SingleClaimOnly,
            )
            .is_err(),
            "a flow-supported observation is not an equilibrium law"
        );
    }
    let mut native = FreeSurface3::new(
        16,
        16,
        12,
        plan(&bulk).flow().tau,
        [0.0; 3],
        0.001,
        ContactModel3::Neutral,
        drop_shape,
    )
    .unwrap();
    assert_eq!(
        native.set_planar_contact_angle(PlanarWall3::ZMin, 1e-200),
        Err(FreeSurfaceError3::InvalidContactAngle)
    );
    assert_eq!(native.planar_contact_angle(), None);
}

#[test]
fn g4_cancelled_step_preserves_mass_populations_boundary_and_retry() {
    let bulk = liquid(0, 300.0, Some(0.05));
    let surface = wetting(&bulk, PI / 3.0);
    let setup = plan(&bulk)
        .with_planar_wetting(PlanarWall3::ZMin, &surface)
        .unwrap();
    let mut a = setup.free_surface3(domain(), drop_shape).unwrap();
    let mut b = setup.free_surface3(domain(), drop_shape).unwrap();
    let gate = CancelGate::new_clock_free();
    for _ in 0..4 {
        a.step_with_gate(&gate).unwrap();
        b.step().unwrap();
    }
    let bits = a.lattice().grid().state_bits();
    let mass = a.ledger_mass_kg().unwrap().to_bits();
    let conversions = a.lattice().conversions();
    let boundary = a.lattice().planar_contact_angle();
    let cancelled = CancelGate::new_clock_free();
    cancelled.request();
    assert_eq!(
        a.step_with_gate(&cancelled),
        Err(MaterialLbmError::FreeSurface(FreeSurfaceError3::Cancelled))
    );
    assert_eq!(a.lattice().grid().state_bits(), bits);
    assert_eq!(a.ledger_mass_kg().unwrap().to_bits(), mass);
    assert_eq!(a.lattice().conversions(), conversions);
    assert_eq!(a.lattice().planar_contact_angle(), boundary);
    a.step_with_gate(&gate).unwrap();
    b.step().unwrap();
    assert_eq!(
        a.lattice().grid().state_bits(),
        b.lattice().grid().state_bits()
    );
    assert_eq!(
        a.ledger_mass_kg().unwrap().to_bits(),
        b.ledger_mass_kg().unwrap().to_bits()
    );
}
