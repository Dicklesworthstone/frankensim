//! G1/G3: integrated expansion reaches native frame compatibility, prestress,
//! layer stresses, modal stiffness and the unchanged reference-mass operator.
//! Synthetic curves are analytical fixtures, not empirical material validation.

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_evidence::ValidityDomain;
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_matdb::{ClaimSet, InterpolationPolicy, MaterialCard, MaterialStateId,
    PropertyClaim, PropertyKey, PropertyValue, Provenance, QueryPoint, UncertaintyModel};
use fs_material::section::{ConcentricElasticSection, ConcentricMaterialLayer};
use fs_material::section::thermal::ConcentricThermoelasticSection;
use fs_material::state_point::{DENSITY_PROPERTY, INVERSE_TEMPERATURE_DIMS,
    LINEAR_THERMAL_EXPANSION_COEFFICIENT_PROPERTY as ALPHA, MaterialPropertySelection,
    ScalarAdmissibility, ScalarPropertyRequirement, YOUNG_MODULUS_PROPERTY,
    integrate_isotropic_thermal_expansion, resolve_material_state_point};
use fs_qty::{Density, Pressure, Temperature};
use fs_solid::frame2d::FrameMassFormulation;
use fs_solid::{Frame2d, Geometry, Rayleigh, SolidError};

fn card(name: &str, young: f64, alpha: f64, softening: f64) -> MaterialCard {
    let mut claims = ClaimSet::new();
    for (key, dims, value, slope) in [
        (DENSITY_PROPERTY, Density::DIMS, 4000.0, -8.0),
        (YOUNG_MODULUS_PROPERTY, Pressure::DIMS, young, -young * softening / 100.0),
        (ALPHA, INVERSE_TEMPERATURE_DIMS, alpha, alpha / 100.0),
    ] {
        claims.insert_claim(PropertyClaim {
            key: PropertyKey::new(key, dims),
            value: PropertyValue::Curve {
                abscissa: "T".into(), abscissa_dims: Temperature::DIMS,
                knots: vec![(260.0, value - 40.0 * slope), (400.0, value + 100.0 * slope)], dims,
            },
            validity: ValidityDomain::unconstrained().with("T", 260.0, 400.0),
            uncertainty: UncertaintyModel::Unstated,
            interpolation: InterpolationPolicy::LinearInside, observations: vec![],
            provenance: Provenance { source: "synthetic thermal frame".into(),
                license: "test-only".into(), artifact: None },
        }).unwrap();
    }
    MaterialCard::assemble(MaterialStateId { chemistry: name.into(), phase: "solid".into(),
        process: "synthetic".into(), revision: 0 }, claims, vec![]).unwrap()
}
fn point(t: f64) -> QueryPoint {
    QueryPoint::new().with("T", t).unwrap()
}
fn thermal(cards: &[(&MaterialCard, f64)], t: f64, reference_t: f64) -> ConcentricThermoelasticSection {
    let mut layers = Vec::new();
    let mut expansions = Vec::new();
    for &(card, radius) in cards {
        let material = resolve_material_state_point(card, &point(t), &[
            ScalarPropertyRequirement::try_new(DENSITY_PROPERTY, Density::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
            ScalarPropertyRequirement::try_new(YOUNG_MODULUS_PROPERTY, Pressure::DIMS, ScalarAdmissibility::StrictlyPositive).unwrap(),
        ], MaterialPropertySelection::SingleClaimOnly).unwrap();
        layers.push(ConcentricMaterialLayer { outer_radius_m: radius, material });
        expansions.push(integrate_isotropic_thermal_expansion(card, &point(reference_t),
            &point(t), MaterialPropertySelection::SingleClaimOnly).unwrap());
    }
    ConcentricThermoelasticSection::try_new(
        ConcentricElasticSection::try_new(0.0, layers).unwrap(), expansions,
    ).unwrap()
}
fn frame(reference: ConcentricElasticSection, angle: f64) -> Frame2d {
    let mut f = Frame2d::with_mass_formulation(Geometry::Linear, FrameMassFormulation::ConsistentReference);
    let a = f.add_node(0.0, 0.0);
    let b = f.add_node(angle.cos(), angle.sin());
    f.fix(a, [true; 3]);
    f.add_material_element(a, b, 3, reference).unwrap();
    f
}
fn with_cx(f: impl FnOnce(&Cx<'_>)) {
    let gate = CancelGate::new_clock_free();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(&Cx::new(&gate, arena,
            StreamKey { seed: 7, kernel_id: 7, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic));
    });
}
fn close(a: f64, b: f64, scale: f64) {
    assert!((a-b).abs() <= 2e-9 * scale.max(a.abs()).max(b.abs()), "{a:e} != {b:e}");
}
fn free_strain(alpha: f64, t: f64) -> f64 {
    let dt = t - 300.0;
    alpha * (dt + dt * dt / 200.0)
}

#[test]
fn free_heating_and_current_modulus_share_original_mass_and_modal_pencil() {
    let card = card("softening", 2e9, 1e-5, 0.5);
    let reference = thermal(&[(&card, 0.01)], 300.0, 300.0).section().clone();
    let original = frame(reference.clone(), 0.0);
    let old_modes = original.natural_frequencies().unwrap();
    let acceleration = vec![0.1, -9.81, 0.7, 0.3, -9.81, -0.2];
    let original_mass_action = original.inertial_forces(&acceleration).unwrap();
    with_cx(|cx| {
        let section = thermal(&[(&card, 0.01)], 350.0, 300.0);
        let mut hot = original.with_thermal_material_sections(0, vec![section.clone(); 3], 0.02, cx).unwrap();
        let mut force = vec![0.0; 6];
        force[3] = 1.0;
        force[4] = 0.01;
        hot.static_load(&force, 1).unwrap();
        let ea = 1.5e9 * std::f64::consts::PI * 0.01_f64.powi(2);
        let ei = 1.5e9 * std::f64::consts::PI * 0.01_f64.powi(4) / 4.0;
        close(hot.displacements()[3], free_strain(1e-5, 350.0) + 1.0 / ea, 0.01);
        close(hot.displacements()[4], 0.01 / (3.0 * ei), 0.01);
        assert_eq!(hot.inertial_forces(&acceleration).unwrap(), original_mass_action);
        assert_eq!(hot.material_section(0), Some(&reference));
        assert!(hot.thermal_material_section(0, 0).unwrap().source().section().line_density_kg_m()
            < reference.line_density_kg_m());
        for (cold, current) in old_modes.iter().zip(hot.natural_frequencies().unwrap()) {
            close(current, 0.75_f64.sqrt() * cold, 1.0);
        }
        let response = hot.thermal_section_response(0, 0).unwrap();
        close(response.axial.force_n(), 1.0, 1.0);
        close(response.bending_moment_n_m, 0.01, 1.0);
        let axial = 1.0 / ea;
        close(response.strain_energy_j_m, 0.5 * ea * axial * axial + 0.5 * 0.01 * 0.01 / ei, 1.0);
        let spread = 1.5e9 * 0.01 * 0.01 / ei;
        close(response.layer_stress_range_pa[0][0], 1.5e9 * axial - spread, 1e5);
        close(response.layer_stress_range_pa[0][1], 1.5e9 * axial + spread, 1e5);
    });
}

#[test]
fn restrained_spatial_temperature_preserves_compatibility_and_compressive_prestress() {
    let card = card("spatial", 2e9, 1e-5, 0.0);
    let reference = thermal(&[(&card, 0.01)], 300.0, 300.0).section().clone();
    let original = frame(reference, 0.0);
    with_cx(|cx| {
        let stations = [300.0, 330.0, 360.0].map(|t| thermal(&[(&card, 0.01)], t, 300.0));
        let average = 1e-5 * (30.0 + 3600.0 / 600.0);
        let ea = 2e9 * std::f64::consts::PI * 0.01_f64.powi(2);
        let mut free = original.with_thermal_material_sections(0, stations.to_vec(), 0.02, cx).unwrap();
        free.static_load(&[0.0; 6], 1).unwrap();
        close(free.displacements()[3], average, 0.01);
        let mut fixed = original.clone();
        fixed.fix(1, [true; 3]);
        let mut fixed = fixed.with_thermal_material_sections(0, stations.to_vec(), 0.02, cx).unwrap();
        fixed.static_load(&[0.0; 6], 1).unwrap();
        assert_eq!(fixed.displacements(), &[0.0; 6]);
        close(fixed.element_forces(0)[0], -ea * average, 1.0);
        let mut compatible = 0.0;
        for (station, (t, weight)) in [300.0, 330.0, 360.0].into_iter().zip([1.0/6.0, 2.0/3.0, 1.0/6.0]).enumerate() {
            let response = fixed.thermal_section_response(0, station).unwrap();
            close(response.axial.total_strain(), free_strain(1e-5, t) - average, 0.01);
            close(response.axial.layers()[0].stress_pa(), -2e9 * average, 1.0);
            compatible += weight * response.axial.total_strain();
        }
        close(compatible, 0.0, 0.01);
    });
}

#[test]
fn free_layered_member_retains_opposing_thermal_stresses_and_storage() {
    let a = card("core", 2e9, 1e-5, 0.0);
    let b = card("skin", 1e9, 3e-5, 0.0);
    let layers = [(&a, 0.01), (&b, 0.02)];
    let original = frame(thermal(&layers, 300.0, 300.0).section().clone(), 0.0);
    with_cx(|cx| {
        let mut hot = original.with_thermal_material_sections(0, vec![thermal(&layers, 350.0, 300.0); 3], 0.02, cx).unwrap();
        hot.static_load(&[0.0; 6], 1).unwrap();
        let e = (2.0 * free_strain(1e-5, 350.0) + 3.0 * free_strain(3e-5, 350.0)) / 5.0;
        close(hot.displacements()[3], e, 0.01);
        let response = hot.thermal_section_response(0, 1).unwrap();
        close(response.axial.force_n(), 0.0, 1.0);
        close(response.axial.layers()[0].stress_pa(), 2e9 * (e - free_strain(1e-5, 350.0)), 1.0);
        close(response.axial.layers()[1].stress_pa(), 1e9 * (e - free_strain(3e-5, 350.0)), 1.0);
        assert!(response.axial.layers()[0].stress_pa() > 0.0);
        assert!(response.axial.layers()[1].stress_pa() < 0.0);
        assert!(response.strain_energy_j_m > 0.0);
    });
}

#[test]
fn prescribed_heating_drives_native_inertia_without_creating_mechanical_energy() {
    let card = card("thermal impulse", 2e9, 1e-5, 0.0);
    let mut original = frame(
        thermal(&[(&card, 0.01)], 300.0, 300.0).section().clone(),
        0.0,
    );
    // One axial free coordinate: native consistent mass is rho*A*L/3.
    original.fix(1, [false, true, true]);
    with_cx(|cx| {
        let mut hot = original.with_thermal_material_sections(
            0,
            vec![thermal(&[(&card, 0.01)], 350.0, 300.0); 3],
            0.02,
            cx,
        ).unwrap();
        let initial_energy = hot.thermal_section_response(0, 0).unwrap().strain_energy_j_m;
        let dt = 1e-5;
        let history = hot.newmark(
            &[0.0; 51], dt, (1.0, 0.0), Rayleigh { a0: 0.0, a1: 0.0 },
        ).unwrap();
        let omega = (3.0_f64 * 2e9 / 4000.0).sqrt();
        let phase = 2.0 * (0.5 * omega * dt).atan();
        let free = free_strain(1e-5, 350.0);
        for (step, displacement) in history.u.iter().enumerate() {
            close(displacement[3], free * (1.0 - (step as f64 * phase).cos()), 0.01);
        }
        let final_storage = hot.thermal_section_response(0, 0).unwrap().strain_energy_j_m;
        close(final_storage + history.kinetic_energy.last().unwrap(), initial_energy, 1.0);
        close(*history.internal_work.last().unwrap(), final_storage - initial_energy, 1.0);
        assert!(hot.velocities()[3] > 0.0);
    });
}

#[test]
fn cooling_negative_expansion_and_rotation_use_original_reference_without_drift() {
    with_cx(|cx| {
        for alpha in [1e-5, -1e-5] {
            let card = card("signed expansion", 2e9, alpha, 0.0);
            let angle = 0.73_f64;
            let original = frame(thermal(&[(&card, 0.01)], 300.0, 300.0).section().clone(), angle);
            let mut current = original.clone();
            for t in [350.0, 280.0, 370.0, 300.0] {
                current = current.with_thermal_material_sections(0,
                    vec![thermal(&[(&card, 0.01)], t, 300.0); 3], 0.02, cx).unwrap();
                current.static_load(&[0.0; 6], 1).unwrap();
                let e = free_strain(alpha, t);
                close(current.displacements()[3], e * angle.cos(), 0.01);
                close(current.displacements()[4], e * angle.sin(), 0.01);
                close(current.thermal_section_response(0, 1).unwrap().strain_energy_j_m, 0.0, 1.0);
            }
            assert_eq!(original.displacements(), &[0.0; 6]);
        }
    });
}

#[test]
fn invalid_source_reference_bounds_and_cancellation_leave_original_frame_intact() {
    let a = card("original", 2e9, 1e-5, 0.0);
    let b = card("foreign", 2e9, 1e-5, 0.0);
    let original = frame(thermal(&[(&a, 0.01)], 300.0, 300.0).section().clone(), 0.0);
    with_cx(|cx| {
        for source in [
            thermal(&[(&b, 0.01)], 350.0, 300.0),
            thermal(&[(&a, 0.01)], 350.0, 310.0),
            thermal(&[(&a, 0.02)], 350.0, 300.0),
        ] {
            assert!(original.with_thermal_material_sections(0, vec![source; 3], 0.02, cx).is_err());
        }
        let source = thermal(&[(&a, 0.01)], 350.0, 300.0);
        assert!(original.with_thermal_material_sections(0, vec![source.clone(); 2], 0.02, cx).is_err());
        assert!(original.with_thermal_material_sections(0, vec![source; 3], 1e-6, cx).is_err());
        assert_eq!(original.displacements(), &[0.0; 6]);
        assert_eq!(original.element_forces(0), [0.0; 3]);
        assert!(original.thermal_material_section(0, 0).is_none());
    });
    let gate = CancelGate::new_clock_free();
    gate.request();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        let cx = Cx::new(&gate, arena,
            StreamKey { seed: 7, kernel_id: 7, tile: 0, iteration: 0 },
            Budget::INFINITE, ExecMode::Deterministic);
        assert!(matches!(original.with_thermal_material_sections(0,
            vec![thermal(&[(&a, 0.01)], 350.0, 300.0); 3], 0.02, &cx), Err(SolidError::Cancelled)));
    });
}
