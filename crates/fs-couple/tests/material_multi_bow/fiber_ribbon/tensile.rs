//! G0/G1/G3: source-bound tensile admission in the actual fiber contact owner.
//! Manufactured elastic data, not measured fiber rupture validation.
use super::*;
use fs_couple::render::schedule::force::coupled::contact::friction::nonlinear::ModalPotential;
use fs_material::section::relaxation::LayerRelaxationProperties;
use fs_phs::Storage;

#[path = "tensile/hereditary.rs"]
mod hereditary;
#[path = "tensile/viscous.rs"]
mod viscous;

const LIMIT: &str = "fiber_tensile_limit";
const RADIUS: f64 = 40e-6;
const YOUNG: f64 = 3e9;

fn ribbon(cells: usize, limit: f64, gate: &CancelGate) -> TensionedFiberRibbon {
    ribbon_with_viscosities(cells, limit, 0.0, 0.0, gate)
}

fn ribbon_with_viscosities(
    cells: usize,
    limit: f64,
    axial_eta: f64,
    bending_eta: f64,
    gate: &CancelGate,
) -> TensionedFiberRibbon {
    let section = ConcentricElasticSection::try_new(
        0.0,
        vec![ConcentricMaterialLayer {
            outer_radius_m: RADIUS,
            material: fiber_material_with_limit(1300.0, YOUNG, axial_eta, bending_eta, Some(limit)),
        }],
    )
    .unwrap();
    ribbon_with_section(cells, section, 1, 1, gate)
}

fn ribbon_with_section(
    cells: usize,
    section: ConcentricElasticSection,
    axial_modes: usize,
    transverse_modes: usize,
    gate: &CancelGate,
) -> TensionedFiberRibbon {
    TensionedFiberRibbon::try_new(
        FiberRibbonGeometry {
            length_m: 0.6,
            width_m: 0.01,
            fiber_count: 180,
            tension_n: 40.0,
        },
        FiberRibbonResolution {
            cells,
            axial_modes,
            transverse_modes,
        },
        section
            .clone()
            .with_axial_viscosity(YOUNG_VISCOSITY_PROPERTY)
            .unwrap(),
        section.with_kelvin_voigt_bending(ETA).unwrap(),
        48000,
        ModalAcousticTimeBudget::audible_reference(),
        gate,
    )
    .unwrap()
    .with_axial_transverse_strain(strain_limits(), strain_config())
    .unwrap()
}

fn near(a: f64, b: f64) {
    assert!((a - b).abs() < 3e-13 * b.abs().max(1e-8), "{a:e} != {b:e}");
}

#[test]
fn tensile_fiber_stress_uses_original_radius_and_is_width_refinement_invariant() {
    let gate = CancelGate::new_clock_free();
    let mut previous = None;
    for cells in [1, 3, 9] {
        let original = ribbon(cells, 100e6, &gate);
        let guarded = original.clone().with_tensile_limit(LIMIT).unwrap();
        let law = guarded.strain().unwrap();
        let scale = (0.5 * guarded.cell_line_density_kg_m() * guarded.geometry().length_m).sqrt();
        let axial = 0.0001;
        let normal = 0.002;
        let q = [axial * scale, normal * scale];
        let assessment = law.tensile_assessment(&q).unwrap().unwrap();
        let wave = core::f64::consts::PI / 0.6;
        let area = core::f64::consts::PI * RADIUS * RADIUS;
        let strain = 40.0 / (180.0 * YOUNG * area) + wave * axial + 0.25 * (wave * normal).powi(2);
        let curvature = wave * wave * normal;
        let stress = YOUNG * (strain + RADIUS * curvature);
        near(assessment.per_fiber_tension_bound_n, YOUNG * area * strain);
        near(assessment.curvature_bound_per_m, curvature);
        near(assessment.governing_tensile_stress_bound_pa, stress);
        near(assessment.utilization, stress / 100e6);
        assert_eq!(
            law.tensile_capacity().unwrap().section(),
            original.axial_material().elastic()
        );
        assert_eq!(
            law.hamiltonian(&q),
            original.strain().unwrap().hamiltonian(&q)
        );
        if let Some(old) = previous {
            near(assessment.governing_tensile_stress_bound_pa, old);
        }
        previous = Some(assessment.governing_tensile_stress_bound_pa);
        // Axial demand passes this source value; actual per-fiber bending does not.
        let bending_limit = YOUNG * strain + 0.5 * YOUNG * RADIUS * curvature;
        let too_weak = ribbon(cells, bending_limit, &gate)
            .with_tensile_limit(LIMIT)
            .unwrap();
        assert!(
            too_weak
                .strain()
                .unwrap()
                .tensile_assessment(&q)
                .unwrap()
                .unwrap()
                .exceeds_limit()
        );
        assert!(too_weak.strain().unwrap().validate_coordinates(&q).is_err());
    }
}

#[test]
fn tensile_fiber_preload_sources_refuse_and_memory_preserves_its_original_owner() {
    let gate = CancelGate::new_clock_free();
    let preload = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
    assert!(
        ribbon(3, 0.99 * preload, &gate)
            .with_tensile_limit(LIMIT)
            .is_err()
    );
    assert!(
        fibers(3, 1, 1300.0, YOUNG, 0.0, 40.0, 48000, &gate)
            .with_axial_transverse_strain(strain_limits(), strain_config())
            .unwrap()
            .with_tensile_limit(LIMIT)
            .is_err()
    );
    let original = ribbon(3, 100e6, &gate);
    let spectrum = original
        .axial_material()
        .elastic()
        .clone()
        .with_prony_bending(&[LayerRelaxationProperties {
            equilibrium_modulus: "young_modulus".into(),
            branches: vec![("delta".into(), "tau".into())],
        }])
        .unwrap();
    let guarded = original.clone().with_tensile_limit(LIMIT).unwrap();
    assert!(guarded.clone().with_tensile_limit(LIMIT).is_err());
    let forward = guarded
        .with_material_relaxation(spectrum.clone(), 64, &gate)
        .unwrap();
    let reverse = original
        .with_material_relaxation(spectrum, 64, &gate)
        .unwrap()
        .with_tensile_limit(LIMIT)
        .unwrap();
    for source in [forward, reverse] {
        assert!(source.strain().unwrap().tensile_capacity().is_some());
        // A public strain clone cannot silently assess total hereditary stress
        // without the original two memory owners. Atomic bank binding supplies them.
        assert!(source.strain().unwrap().requires_velocity_validation());
        assert!(
            source
                .strain()
                .unwrap()
                .validate_state(&[0.0; 2], &[0.0; 2])
                .is_err()
        );
    }
}

#[test]
fn tensile_fiber_viscous_sources_activate_native_endpoint_validation() {
    let gate = CancelGate::new_clock_free();
    // Either source direction needs actual endpoint rates, even when its
    // coefficient is tiny or loss is owned by the original linear propagator.
    for eta in [50.0, 1e-12] {
        for (axial_eta, bending_eta) in [(eta, 0.0), (0.0, eta)] {
            let original = ribbon_with_viscosities(3, 100e6, axial_eta, bending_eta, &gate);
            assert_eq!(
                original.axial_material().viscous_rigidity_n_s() > 0.0,
                axial_eta > 0.0
            );
            assert_eq!(
                original.bending_material().viscous_stiffness_n_m2_s() > 0.0,
                bending_eta > 0.0
            );
            let guarded = original.with_tensile_limit(LIMIT).unwrap();
            let strain = guarded.strain().unwrap();
            assert!(strain.requires_velocity_validation());
            strain.validate_state(&[0.0; 2], &[0.0; 2]).unwrap();
        }
    }
    // Explicitly lossless laws retain their original coordinate-only path.
    let elastic = ribbon(3, 100e6, &gate).with_tensile_limit(LIMIT).unwrap();
    assert!(!elastic.strain().unwrap().requires_velocity_validation());
}

#[test]
fn tensile_fiber_crossing_rolls_back_the_original_contact_bank() {
    with_cx(|cx, gate| {
        let preload = 40.0 / (180.0 * core::f64::consts::PI * RADIUS * RADIUS);
        let source = ribbon(3, preload + 100.0, gate)
            .with_tensile_limit(LIMIT)
            .unwrap();
        let mut bank = base(cx, gate, 7800.0)
            .with_tensioned_fiber_ribbon(source, 0.4, ribbon_limits(), gate)
            .unwrap();
        let first = bank.source().fiber_ribbon().unwrap().component_start();
        let start: usize = bank.system().components()[..first]
            .iter()
            .map(|c| c.modes().len())
            .sum();
        let mut force = vec![0.0; bank.system().mode_count()];
        force[start] = 1.0;
        let before = states(&bank);
        let energy = bank.system().total_energy_j().unwrap();
        let error = bank.step(&force, &vec![drive(0.0, 0.0); 3]).unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(states(&bank), before);
        assert_eq!(bank.system().total_energy_j().unwrap(), energy);
        assert_eq!(bank.system().samples_rendered(), 0);
        assert!(bank.system().last_frame().is_none());
        // A refused candidate must not poison the accepted state for a retry.
        force.fill(0.0);
        bank.step(&force, &vec![drive(0.0, 0.0); 3]).unwrap();
        assert_eq!(bank.system().samples_rendered(), 1);
    });
}
