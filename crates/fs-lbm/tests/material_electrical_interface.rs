//! Synthetic source cards through the existing current and thermal solvers.
//! Values are deliberately simple equation fixtures, not measured contacts.
use fs_evidence::ValidityDomain;
use fs_exec::CancelGate;
use fs_lbm::conjugate::{
    Conductor, ElectricalContact, ElectricSetup, Electrode, ElectrodeDrive,
    EnergyConfig, FlowField, FluidProperties, SolidMaterial, ThermalFace,
    ThermalSetup, Voxel, VoxelDomain, solve_energy,
};
use fs_matdb::{
    ClaimSet, InterfaceSystemCard, InterpolationPolicy, MaterialStateId, PropertyClaim,
    PropertyKey, PropertyValue, Provenance, QueryPoint, SurfaceSpec, SystemContext,
    UncertaintyModel,
};
use fs_material::electrical_interface::{
    ELECTRICAL_CONTACT_RESISTANCE_DIMS as CONTACT_DIMS,
    ELECTRICAL_CONTACT_RESISTANCE_PROPERTY as CONTACT_KEY, resolve_electrical_interface_state,
};
use fs_material::state_point::MaterialPropertySelection;
use fs_qty::{Dims, Temperature};

fn material(name: &str) -> MaterialStateId {
    MaterialStateId { chemistry: name.into(), phase: "solid".into(),
        process: "synthetic fixture".into(), revision: 0 }
}
fn card(dims: Dims, scale: f64) -> InterfaceSystemCard {
    let mut claims = ClaimSet::new();
    claims.insert_claim(PropertyClaim {
        key: PropertyKey::new(CONTACT_KEY, dims),
        value: PropertyValue::Curve { abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
            knots: vec![(300.0, 3.0 * scale), (400.0, 6.0 * scale)], dims },
        validity: ValidityDomain::unconstrained().with("T", 300.0, 400.0).with("pressure", 1.0, 2.0),
        uncertainty: UncertaintyModel::RelativeHalfWidth { fraction: 0.01, confidence: 0.95 },
        interpolation: InterpolationPolicy::LinearInside, observations: vec![],
        provenance: Provenance { source: "synthetic electrical interface".into(),
            license: "test-only".into(), artifact: None },
    }).unwrap();
    InterfaceSystemCard::assemble(
        SurfaceSpec { material: material("a"), texture_frame: "finish-a".into() },
        SurfaceSpec { material: material("b"), texture_frame: "finish-b".into() },
        SystemContext { medium: "dry".into(), third_body: None,
            environment: "synthetic atmosphere".into(), history: "virgin".into() },
        claims, vec![],
    ).unwrap()
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap().with("pressure", 1.5).unwrap()
}
fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 1e-8 * b.abs().max(1.0), "{a} != {b}");
}
fn circuit() -> (VoxelDomain, ElectricSetup) {
    let domain = VoxelDomain::from_fn(4, 1, 1, 1.0,
        |p| Voxel::Solid(u16::from(p[0] >= 2.0))).unwrap();
    let electric = ElectricSetup {
        conductors: vec![Some(Conductor::new(2.0)), Some(Conductor::new(4.0))],
        electrodes: vec![
            Electrode { lo: [0,0,0], hi: [1,1,1], drive: ElectrodeDrive::Voltage(0.0) },
            Electrode { lo: [3,0,0], hi: [4,1,1], drive: ElectrodeDrive::Current(1.0) },
        ],
    };
    (domain, electric)
}

#[test]
fn source_only_substitutions_change_current_voltage_heat_and_real_temperature_field() {
    let (domain, mut electric) = circuit();
    let gate = CancelGate::new();
    let mut faces = [ThermalFace::Adiabatic; 6];
    faces[0] = ThermalFace::Temperature(300.0);
    faces[1] = ThermalFace::Temperature(300.0);
    let solids = [SolidMaterial::new("a", 1.0), SolidMaterial::new("b", 1.0)];
    let flow = FlowField::quiescent(&domain);
    let mut temperatures = Vec::new();
    for (t, scale) in [(300.0, 1.0), (350.0, 1.0), (300.0, 2.0)] {
        let source = card(CONTACT_DIMS, scale);
        let resolved = resolve_electrical_interface_state(&source, &point(t),
            MaterialPropertySelection::SingleClaimOnly).unwrap();
        resolved.require_surface_materials(&material("a"), &material("b")).unwrap();
        assert_eq!(resolved.surfaces()[0], *source.surface_a());
        assert_eq!(resolved.resolved().card_identity(), source.content_hash());
        assert_eq!(resolved.resolved().state_coordinate("T"), Some(t));
        let r = scale * (3.0 + 0.03 * (t - 300.0));
        near(resolved.resistance_ohm_m2(), r);
        let contacts = [ElectricalContact { cells: [1,2], resistance_ohm_m2: resolved.resistance_ohm_m2(),
            heat_fraction_to_first: 0.25 }];
        electric.electrodes[1].drive = ElectrodeDrive::Current(1.0);
        let field = electric.solve_with_contacts(&domain, &[t; 4], &contacts,
            |_, m, _| Ok(if m == 0 {2.0} else {4.0}), &gate).unwrap();
        near(field.electrodes[1].voltage_v, 6.0 + r);
        near(field.joule_w[1], 2.0 + r * 0.25);
        near(field.joule_w[2], 4.0 + r * 0.75);
        near(field.total_joule_w, field.delivered_w);
        let mut thermal = ThermalSetup::new(faces);
        thermal.power_w = field.joule_w.clone(); // already includes the contact heat
        let heat = solve_energy(&domain, &FluidProperties::dry_air_300k(), &solids,
            &flow, &thermal, &EnergyConfig::default(), &gate).unwrap();
        // Unit conductances and half-width boundary resistances: exact FV solution.
        let q1 = 2.0 + 0.25*r;
        let q2 = 4.0 + 0.75*r;
        let left_flux = (2.5*q1 + 1.5*q2) / 4.0;
        for (got, rise) in heat.temperature.iter().zip([
            0.5*left_flux, 1.5*left_flux, 2.5*left_flux-q1, 3.5*left_flux-2.0*q1-q2,
        ]) { near(*got, 300.0 + rise); }
        near(heat.report.balance.source_w, 6.0+r);
        near(heat.report.balance.boundary_outflow_w, field.delivered_w);
        near(heat.report.balance.residual_w, 0.0);
        temperatures.push(heat.temperature[2]);
        electric.electrodes[1].drive = ElectrodeDrive::Voltage(9.0);
        let voltage = electric.solve_with_contacts(&domain, &[t; 4], &contacts,
            |_, m, _| Ok(if m == 0 {2.0} else {4.0}), &gate).unwrap();
        near(voltage.electrodes[1].current_a, 9.0/(6.0+r));
    }
    assert!(temperatures.windows(2).all(|t| t[1] > t[0]));
}

#[test]
fn electrical_contact_units_cannot_be_replaced_by_bulk_or_thermal_resistance() {
    for wrong in [Dims([3,1,-3,0,-2,0]), Dims([0,-1,3,1,0,0]), Dims::NONE] {
        assert!(resolve_electrical_interface_state(&card(wrong, 1.0), &point(350.0),
            MaterialPropertySelection::SingleClaimOnly).is_err());
    }
}

#[test]
fn source_domain_surface_order_and_processing_state_remain_mandatory() {
    let source = card(CONTACT_DIMS, 1.0);
    for query in [point(299.0), point(401.0), QueryPoint::new().with("T", 350.0).unwrap(),
        QueryPoint::new().with("T", 350.0).unwrap().with("pressure", 3.0).unwrap()]
    {
        assert!(resolve_electrical_interface_state(&source, &query,
            MaterialPropertySelection::SingleClaimOnly).is_err());
    }
    let resolved = resolve_electrical_interface_state(&source, &point(350.0),
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert!(resolved.require_surface_materials(&material("b"), &material("a")).is_err());
    let mut changed = material("b"); changed.revision = 1;
    assert!(resolved.require_surface_materials(&material("a"), &changed).is_err());
    let again = resolve_electrical_interface_state(&source, &point(350.0),
        MaterialPropertySelection::SingleClaimOnly).unwrap();
    assert_eq!(resolved, again);
}

#[test]
fn sourced_zero_is_distinct_from_missing_contact_data() {
    let zero = card(CONTACT_DIMS, 0.0);
    assert_eq!(resolve_electrical_interface_state(&zero, &point(350.0),
        MaterialPropertySelection::SingleClaimOnly).unwrap().resistance_ohm_m2(), 0.0);
    let missing = InterfaceSystemCard::assemble(zero.surface_a().clone(), zero.surface_b().clone(),
        zero.context().clone(), ClaimSet::new(), vec![]).unwrap();
    assert!(resolve_electrical_interface_state(&missing, &point(350.0),
        MaterialPropertySelection::SingleClaimOnly).is_err());
    assert!(resolve_electrical_interface_state(&card(CONTACT_DIMS, -1.0), &point(350.0),
        MaterialPropertySelection::SingleClaimOnly).is_err());
}
