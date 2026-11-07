//! G0/G1/G3: local axial waves, actual moving ends and original layer strengths.
//! Manufactured elastic specimens, not measured string fracture validation.

use super::*;
use fs_couple::bowed_string::runtime::reciprocal::bank::nonlinear::{
    NonlinearMaterialString, longitudinal::LongitudinalStringResolution,
    supported::longitudinal::MaterialLongitudinalSupport,
};
use fs_couple::render::plate::impact::{
    ImpactError,
    relaxation::InitialMemory,
    supported::{
        StringEndMotion,
        longitudinal::{StringSupportGeometry, relaxation::SupportedRelaxationArm},
    },
};
use fs_exec::{CancelGate, Cx};

const RATE: u32 = 1_536_000;
#[path = "supported/hereditary.rs"]
mod hereditary;

fn rail(position: f64) -> TranslatingSupport {
    TranslatingSupport {
        mass_kg: 0.05,
        stiffness_n_m: 1000.0,
        damping_n_s_m: 0.0,
        initial_position_m: position,
        initial_velocity_m_s: 0.0,
        maximum_travel_m: 0.02,
        maximum_slope: 0.2,
    }
}

fn polar(
    cx: &Cx<'_>,
    gate: &CancelGate,
    source: MaterialStringSource,
    amplitude: [f64; 2],
    bridge: [f64; 2],
) -> NonlinearMaterialString {
    let mut state = initial(6, amplitude[0]);
    state[3].displacement_m_sqrt_kg = initial(1, amplitude[1])[0].displacement_m_sqrt_kg;
    source
        .with_transverse_polarizations(TransverseDirections::try_new([1.0, 0.0], bridge).unwrap())
        .unwrap()
        .into_nonlinear(cx, gas(), RATE, 0.2, state, gate)
        .unwrap()
}

fn attachment(
    modes: &[usize],
    axial: Vec<ModalAcousticState>,
    bridge: [f64; 2],
) -> MaterialLongitudinalSupport {
    MaterialLongitudinalSupport {
        moving_end: StringEndMotion::Second,
        geometry: StringSupportGeometry::try_new(
            [[0.0; 3], [LENGTH, 0.0, 0.0]],
            [0.0, 1.0, 0.0],
            [0.6, 0.8 * bridge[0], 0.8 * bridge[1]],
        )
        .unwrap(),
        resolution: LongitudinalStringResolution {
            axial_modes: modes.to_vec(),
            maximum_slope: 0.2,
            maximum_strain: 0.01,
        },
        initial_axial: axial,
    }
}

fn image(
    source: NonlinearMaterialString,
    attachment: MaterialLongitudinalSupport,
    support: TranslatingSupport,
    gate: &CancelGate,
) -> Result<SupportedMaterialStrings, ImpactError> {
    SupportedMaterialStrings::in_stationary_gas_with_longitudinal_waves(
        vec![source],
        support,
        vec![attachment],
        300.0,
        gate,
    )
}

fn initial_coordinates(image: &SupportedMaterialStrings) -> Vec<f64> {
    image.bodies()[0]
        .initial
        .iter()
        .map(|s| s.displacement_m_sqrt_kg)
        .collect()
}

#[test]
fn g1_sparse_axial_peak_and_actual_endpoint_extension_govern_tensile_admission() {
    with_cx(|cx, gate| {
        let displacement = 2e-5;
        let amplitudes = [1e-5, 2e-6];
        let axial: Vec<_> = amplitudes.iter().map(|&a| initial(1, a)[0]).collect();
        let source = polar(cx, gate, tensile_source(1e9), [0.0; 2], [1.0, 0.0])
            .with_tensile_limit(LIMIT)
            .unwrap();
        let ea = source.axial_rigidity_n();
        let loaded = image(
            source,
            attachment(&[1, 3], axial.clone(), [1.0, 0.0]),
            rail(displacement),
            gate,
        )
        .unwrap();
        let q = initial_coordinates(&loaded);
        let assessment = loaded
            .supported()
            .tensile_assessment(0, &q)
            .unwrap()
            .unwrap();

        // Independent derivative of the actual affine endpoint displacement
        // plus two solid-mass-normalized sine waves. At x=0 both cosines are 1.
        let wave = core::f64::consts::PI / LENGTH;
        let mean = 0.6 * displacement / LENGTH + 0.5 * (0.8 * displacement / LENGTH).powi(2);
        let spread = wave * (amplitudes[0] + 3.0 * amplitudes[1]);
        let tension = 20.0 + ea * (mean + spread);
        let area = core::f64::consts::PI * RADIUS * RADIUS;
        near(assessment.tension_n, tension);
        near(assessment.governing_tensile_stress_bound_pa, tension / area);
        assert_eq!(assessment.curvature_bound_per_m, 0.0);
        let mean_stress = (20.0 + ea * mean) / area;
        let limit = 0.5 * (mean_stress + tension / area);
        assert!(mean_stress < limit);

        // The original fixed-end state and the mean stress both pass. Only the
        // retained local axial wave makes this initial supported image fail.
        let source = polar(cx, gate, tensile_source(limit), [0.0; 2], [1.0, 0.0])
            .with_tensile_limit(LIMIT)
            .unwrap();
        let error = image(
            source.clone(),
            attachment(&[1, 3], axial, [1.0, 0.0]),
            rail(displacement),
            gate,
        )
        .unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        let admitted = image(
            source,
            attachment(&[1, 3], vec![ModalAcousticState::default(); 2], [1.0, 0.0]),
            rail(displacement),
            gate,
        )
        .unwrap();
        // A diagnostic accepts an overlimit trial without mutating or advancing
        // the guarded image, and uses the same mass map as actual admission.
        assert!(
            admitted
                .supported()
                .tensile_assessment(0, &q)
                .unwrap()
                .unwrap()
                .exceeds_limit()
        );
    });
}

#[test]
fn g3_layered_two_plane_tensile_envelope_keeps_the_endpoint_and_frame_physics() {
    with_cx(|cx, gate| {
        let section = ConcentricElasticSection::try_new(
            0.0,
            vec![
                ConcentricMaterialLayer {
                    outer_radius_m: 0.0002,
                    material: strength_material(YOUNG, 300e6, false),
                },
                ConcentricMaterialLayer {
                    outer_radius_m: RADIUS,
                    material: strength_material(50e9, 45e6, false),
                },
            ],
        )
        .unwrap();
        let ea = section.axial_stiffness_n();
        let source = MaterialStringSource::Concentric {
            template: template(3),
            section: section.clone().with_kelvin_voigt_bending(ETA).unwrap(),
            prestress: StringPrestress::FixedTension(20.0),
        };
        let (amplitude, axial, displacement) = ([3e-4, 4e-4], [2e-6, 1e-6], 2e-5);
        let wave = core::f64::consts::PI / LENGTH;
        let mut results = Vec::new();
        for (bridge, transverse) in [
            ([1.0, 0.0], amplitude),
            (
                [0.6, 0.8],
                [
                    0.6 * amplitude[0] - 0.8 * amplitude[1],
                    0.8 * amplitude[0] + 0.6 * amplitude[1],
                ],
            ),
        ] {
            let source = polar(cx, gate, source.clone(), transverse, bridge)
                .with_tensile_limit(LIMIT)
                .unwrap();
            let loaded = image(
                source,
                attachment(
                    &[1, 2],
                    axial.iter().map(|&a| initial(1, a)[0]).collect(),
                    bridge,
                ),
                rail(displacement),
                gate,
            )
            .unwrap();
            assert_eq!(
                loaded.sources()[0].tensile_capacity().unwrap().section(),
                &section
            );
            let assessment = loaded
                .supported()
                .tensile_assessment(0, &initial_coordinates(&loaded))
                .unwrap()
                .unwrap();

            // Squaring h + a*cos(k*x) gives a mean and exactly the first two
            // cosine coefficients; this independent derivation includes the
            // actual transverse chord and its cross term with bending.
            let chord = [
                0.8 * bridge[0] * displacement / LENGTH,
                0.8 * bridge[1] * displacement / LENGTH,
            ];
            let slopes = transverse.map(|a| wave * a);
            let mean = 0.6 * displacement / LENGTH
                + 0.5 * (chord[0] * chord[0] + chord[1] * chord[1])
                + 0.25 * (slopes[0] * slopes[0] + slopes[1] * slopes[1]);
            let first = wave * axial[0] + chord[0] * slopes[0] + chord[1] * slopes[1];
            let second =
                2.0 * wave * axial[1] + 0.25 * (slopes[0] * slopes[0] + slopes[1] * slopes[1]);
            let curvature = wave * wave * transverse[0].hypot(transverse[1]);
            let peak_strain = 20.0 / ea + mean + first.abs() + second.abs();
            let stress = 50e9 * (peak_strain + RADIUS * curvature);
            assert_eq!(assessment.governing_layer, 1);
            near(assessment.tension_n, ea * peak_strain);
            near(assessment.curvature_bound_per_m, curvature);
            near(assessment.governing_tensile_stress_bound_pa, stress);
            near(assessment.utilization, stress / 45e6);
            for station in 0..=128 {
                let phase = core::f64::consts::PI * station as f64 / 128.0;
                let strain = 20.0 / ea + mean + first * phase.cos() + second * (2.0 * phase).cos();
                let actual = 50e9 * (strain + RADIUS * curvature * phase.sin().abs()).max(0.0);
                assert!(actual <= stress * (1.0 + 1e-14));
            }
            results.push(assessment);
        }
        near(results[0].tension_n, results[1].tension_n);
        near(
            results[0].curvature_bound_per_m,
            results[1].curvature_bound_per_m,
        );
        near(results[0].utilization, results[1].utilization);
    });
}

#[test]
fn g0_axial_wave_tensile_crossing_rolls_back_native_and_prepared_support_steps() {
    with_cx(|cx, gate| {
        let limit = 20.0 / (core::f64::consts::PI * RADIUS * RADIUS) + 5e6;
        let original = polar(cx, gate, tensile_source(limit), [0.0; 2], [1.0, 0.0]);
        let axial = vec![ModalAcousticState {
            displacement_m_sqrt_kg: 0.0,
            velocity_m_sqrt_kg_per_s: 0.5,
        }];
        let selection = attachment(&[1], axial, [1.0, 0.0]);
        let guarded = image(
            original.clone().with_tensile_limit(LIMIT).unwrap(),
            selection.clone(),
            rail(0.0),
            gate,
        )
        .unwrap();
        let plain = image(original, selection, rail(0.0), gate).unwrap();
        let force = vec![0.0; guarded.mode_count()];

        // The same unguarded native mechanics accepts the wave. Its actual
        // local tension independently proves a source-limit crossing, while
        // its original slope, positive-tension and energy domains still pass.
        let mut reference = ImpactSystem::new(
            plain.bodies().to_vec(),
            vec![],
            vec![],
            vec![],
            config_at(RATE),
        )
        .unwrap();
        reference.step(&force, gate).unwrap();
        let q: Vec<_> = reference.state().chunks_exact(2).map(|s| s[0]).collect();
        let axial = plain.supported().longitudinal_observation(0, &q).unwrap();
        assert!(axial.maximum_tension_n / (core::f64::consts::PI * RADIUS * RADIUS) > limit);

        let mut native = ImpactSystem::new(
            guarded.bodies().to_vec(),
            vec![],
            vec![],
            vec![],
            config_at(RATE),
        )
        .unwrap();
        let before = native.state().to_vec();
        let energy = native.stored_energy_j();
        let error = native.step(&force, gate).unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(native.state(), before);
        assert_eq!(native.stored_energy_j(), energy);
        assert_eq!(native.samples(), 0);

        let mut prepared = guarded.prepare(config_at(RATE), gate).unwrap();
        let assessment = prepared.tensile_assessment(0).unwrap().unwrap();
        let error = prepared.step(&force, gate).unwrap_err();
        assert!(format!("{error:?}").contains("tensile"), "{error:?}");
        assert_eq!(prepared.system().state(), before);
        assert_eq!(prepared.system().stored_energy_j(), energy);
        assert_eq!(prepared.system().samples(), 0);
        assert_eq!(prepared.tensile_assessment(0).unwrap().unwrap(), assessment);
    });
}

#[test]
fn g0_supported_tensile_keeps_accepted_dynamics_and_refuses_unassessed_history() {
    with_cx(|cx, gate| {
        let source = polar(cx, gate, tensile_source(1e9), [1e-4, -2e-4], [1.0, 0.0]);
        let selection = attachment(&[1, 2], vec![initial(1, 1e-6)[0]; 2], [1.0, 0.0]);
        let plain = image(source.clone(), selection.clone(), rail(1e-5), gate).unwrap();
        let guarded = image(
            source.with_tensile_limit(LIMIT).unwrap(),
            selection,
            rail(1e-5),
            gate,
        )
        .unwrap();
        let q = initial_coordinates(&guarded);
        assert!(
            plain
                .supported()
                .tensile_assessment(0, &q)
                .unwrap()
                .is_none()
        );
        assert!(guarded.supported().tensile_assessment(1, &q).is_err());
        assert!(guarded.supported().tensile_assessment(0, &q[1..]).is_err());
        assert!(
            guarded
                .supported()
                .tensile_assessment(0, &vec![f64::NAN; q.len()])
                .is_err()
        );
        let arms = vec![SupportedRelaxationArm {
            string: 0,
            axial_rigidity_n: 100.0,
            relaxation_time_s: 0.001,
        }];
        assert!(
            plain
                .supported()
                .clone()
                .with_longitudinal_relaxation(arms.clone(), InitialMemory::Relaxed, 3)
                .is_ok()
        );
        let error = guarded
            .supported()
            .clone()
            .with_longitudinal_relaxation(arms, InitialMemory::Relaxed, 3)
            .unwrap_err();
        assert!(format!("{error:?}").contains("layer memory stresses"));

        let force = vec![0.0; guarded.mode_count()];
        let mut a = plain.prepare(config_at(RATE), gate).unwrap();
        let mut b = guarded.prepare(config_at(RATE), gate).unwrap();
        for _ in 0..16 {
            assert_eq!(
                a.step(&force, gate).unwrap().stored_energy_j,
                b.step(&force, gate).unwrap().stored_energy_j
            );
            assert_eq!(a.system().state(), b.system().state());
            assert!(!b.tensile_assessment(0).unwrap().unwrap().exceeds_limit());
        }
    });
}
