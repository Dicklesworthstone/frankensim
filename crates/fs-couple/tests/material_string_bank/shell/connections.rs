//! Manufactured material-interface fixtures; no measured holding/contact data.
use super::*;
use fs_couple::shear_shell::connections::{MaterialShellConnection, ShellInterfacePoint,
    ShellInterfaceSide, stiffness_per_area_dims, viscosity_per_area_dims};
use fs_material::state_point::{ResolvedInterfaceStatePoint, resolve_interface_state_point};
use fs_matdb::{InterfaceSystemCard, SurfaceSpec, SystemContext};
use fs_plate::shell::mitc::connections::ShellConnectionBudget;

fn fixed_material() -> MaterialStateId {
    MaterialStateId { chemistry: "synthetic support".into(), phase: "solid".into(),
        process: "immobile fixture".into(), revision: 0 }
}
fn interface(k: f64, c: f64, right: MaterialStateId, upper: f64) -> ResolvedInterfaceStatePoint {
    let mut claims = ClaimSet::new();
    let mut requirements = Vec::new();
    for (name, value, dims) in [("contact-stiffness", k, stiffness_per_area_dims()),
        ("contact-viscosity", c, viscosity_per_area_dims())] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(name, dims), value: PropertyValue::Scalar { value, dims },
            validity: ValidityDomain::unconstrained().with("T", 280.0, 310.0).with("omega", 0.0, upper),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::ConstantWithinValidity,
            provenance: Provenance { source: "manufactured shell support fixture".into(),
                license: "CC0-1.0".into(), artifact: None }, observations: vec![],
        }).unwrap();
        requirements.push(ScalarPropertyRequirement::try_new(name, dims, ScalarAdmissibility::Finite).unwrap());
    }
    let card = InterfaceSystemCard::assemble(
        SurfaceSpec { material: shell_material(700.0, 1000.0).material().clone(), texture_frame: "shell foot".into() },
        SurfaceSpec { material: right, texture_frame: "support foot".into() },
        SystemContext { medium: "synthetic bilateral layer".into(), third_body: None,
            environment: "manufactured".into(), history: "unloaded".into() }, claims, vec![]).unwrap();
    resolve_interface_state_point(&card,
        &QueryPoint::new().with("T", 293.15).unwrap().with("omega", 1.0).unwrap(),
        &requirements, MaterialPropertySelection::SingleClaimOnly).unwrap()
}
fn connection(k: f64, c: f64) -> MaterialShellConnection {
    MaterialShellConnection {
        name: "holding".into(), interface: interface(k, c, fixed_material(), 100_000.0),
        stiffness_property: "contact-stiffness".into(), viscosity_property: "contact-viscosity".into(),
        area_m2: 1e-4, direction: [0.0, 0.0, 1.0],
        left: ShellInterfaceSide::Region { name: "wood".into(),
            points: vec![ShellInterfacePoint { node: 5, weight: 1.0, offset_m: [0.01, 0.0, 0.0] }] },
        right: ShellInterfaceSide::Fixed { material: fixed_material() },
    }
}
fn allowance() -> ShellConnectionBudget {
    ShellConnectionBudget { max_connections: 4, max_terms: 8, max_work: 1_000_000 }
}
fn build_with_body(cx: &Cx<'_>, gate: &CancelGate, body: SharedMaterialShellBody)
    -> Result<ReciprocalStringBank, fs_couple::bowed_string::runtime::reciprocal::ReciprocalStringError>
{
    let mut o = options(0); o.network.max_connections = 16; o.network.max_setup_terms = 131_072;
    ReciprocalStringBank::from_material_shear_shell(cx, vec![member(7800.0).source, layered().source],
        gas(101_325.0), body, None, law(), o, gate)
}
fn build_interface(cx: &Cx<'_>, gate: &CancelGate, k: f64, c: f64) -> ReciprocalStringBank {
    let mut body = shell(700.0, 0.02, 0.17, gate);
    body.specimen = body.specimen.with_interfaces(vec![connection(k, c)], allowance(), gate).unwrap();
    build_with_body(cx, gate, body).unwrap()
}

#[test]
fn material_support_changes_body_modes_loss_and_actual_microphone_pressure() {
    with_cx(|cx, gate| {
        let mut soft = build_interface(cx, gate, 1e8, 2e4);
        let mut stiff = build_interface(cx, gate, 4e8, 2e4);
        let mut viscous = build_interface(cx, gate, 1e8, 8e4);
        let frequency = |bank: &ReciprocalStringBank|
            bank.system().components()[2].modes()[0].angular_frequency_rad_s;
        assert!(frequency(&stiff) > frequency(&soft) * (1.0 + 1e-5));
        assert_eq!(frequency(&soft).to_bits(), frequency(&viscous).to_bits());
        assert_eq!(soft.cards()[0].fundamental_hz(), stiff.cards()[0].fundamental_hz());
        let image = soft.source().material_shell_image().unwrap();
        assert_eq!(image.source().interfaces()[0].input().interface.identity(),
            interface(1e8, 2e4, fixed_material(), 100_000.0).identity());
        assert_eq!(image.source().interfaces()[0].law().stiffness_n_m, 1e4);
        assert_eq!(image.source().interfaces()[0].law().damping_n_s_m, 2.0);
        let row = &image.interface_columns()[0];
        for (i, mode) in image.physical_modes().iter().enumerate() {
            assert!((row[i] - (mode[6 * 5 + 2] - 0.01 * mode[6 * 5 + 4])).abs() < 1e-12);
        }
        let channel = soft.source().body_connections().last().unwrap();
        assert_eq!(channel.left.shapes, *row);
        assert_eq!(channel.stiffness_n_m, 0.0);
        assert_eq!(channel.damping_n_s_m, 2.0);
        let zero = vec![0.0; soft.system().mode_count()];
        let mut stiffness_delta = 0.0_f64;
        let mut viscosity_delta = 0.0_f64;
        let mut peak = 0.0_f64;
        for _ in 0..96 {
            let a = *soft.step(&zero, drive()).unwrap();
            let b = *stiff.step(&zero, drive()).unwrap();
            let c = *viscous.step(&zero, drive()).unwrap();
            for frame in [a, b, c] { assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j); }
            peak = peak.max(a.observer_pressure_pa.abs());
            stiffness_delta = stiffness_delta.max((a.observer_pressure_pa - b.observer_pressure_pa).abs());
            viscosity_delta = viscosity_delta.max((a.observer_pressure_pa - c.observer_pressure_pa).abs());
        }
        assert!(peak > 1e-12 && stiffness_delta > 1e-15 && viscosity_delta > 1e-15);
    });
}

#[test]
fn source_rebind_retains_offsets_and_interface_only_semidefinite_loss_is_live() {
    with_cx(|cx, gate| {
        let mut body = shell(700.0, 0.02, 0.17, gate);
        body.specimen = body.specimen.with_region_material("wood", &shell_material(700.0, 0.0), gate).unwrap();
        body.specimen = body.specimen.with_interfaces(vec![connection(1e8, 2e4)], allowance(), gate).unwrap();
        let next = interface(2e8, 4e4, fixed_material(), 100_000.0);
        body.specimen = body.specimen.with_interface_state("holding", &next, gate).unwrap();
        assert_eq!(body.specimen.interfaces()[0].law().stiffness_n_m, 2e4);
        assert_eq!(body.specimen.interfaces()[0].input().interface.identity(), next.identity());
        let mut bank = build_with_body(cx, gate, body).unwrap();
        let image = bank.source().material_shell_image().unwrap();
        assert_eq!(image.damping().channel_count(), 0);
        assert!(image.damping().matrix_per_s().iter().all(|v| *v == 0.0));
        assert_eq!(bank.source().body_connections().len(), 1);
        assert_eq!(bank.source().body_connections()[0].left.shapes, image.interface_columns()[0]);
        let zero = vec![0.0; bank.system().mode_count()];
        for _ in 0..32 { bank.step(&zero, drive()).unwrap(); }
        assert!(bank.system().total_energy_j().unwrap() > 0.0);
    });
}

#[test]
fn original_interface_roles_geometry_band_and_complete_allowances_are_enforced() {
    with_cx(|cx, gate| {
        let body = shell(700.0, 0.02, 0.17, gate);
        for case in 0..5 {
            let mut input = connection(1e8, 2e4);
            let mut budget = allowance();
            match case {
                0 => input.stiffness_property = "missing".into(),
                1 => if let ShellInterfaceSide::Region { points, .. } = &mut input.left { points[0].weight = 0.5; },
                2 => budget.max_terms = 0,
                3 => input.right = ShellInterfaceSide::Fixed { material: shell_material(700.0, 1000.0).material().clone() },
                _ => {
                    input.interface = interface(1e8, 2e4, shell_material(700.0, 1000.0).material().clone(), 100_000.0);
                    input.right = ShellInterfaceSide::Region { name: "wood".into(),
                        points: vec![ShellInterfacePoint { node: 10, weight: 1.0, offset_m: [0.0; 3] }] };
                }
            }
            assert!(body.specimen.with_interfaces(vec![input], budget, gate).is_err());
        }
        let mut input = connection(1e8, 2e4);
        input.interface = interface(1e8, 2e4, fixed_material(), 10.0);
        let mut narrow = body.clone();
        narrow.specimen = body.specimen.with_interfaces(vec![input], allowance(), gate).unwrap();
        assert!(build_with_body(cx, gate, narrow).is_err());
        let mut detached = body.clone();
        let mut input = connection(1e8, 2e4);
        if let ShellInterfaceSide::Region { points, .. } = &mut input.left { points[0].node = 0; }
        detached.specimen = body.specimen.with_interfaces(vec![input], allowance(), gate).unwrap();
        assert!(build_with_body(cx, gate, detached).is_err());
        gate.request();
        assert!(body.specimen.with_interfaces(vec![connection(1e8, 2e4)], allowance(), gate).is_err());
    });
}
