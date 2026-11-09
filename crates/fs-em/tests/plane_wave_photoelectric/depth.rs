use super::*;
use fs_em::plane_wave::layers::photoelectric::depth::PreparedLayerDepthPhotoelectric;
use fs_material::photoelectric::source::depth::{
    DepthPhotoelectricMaterial, PhotoelectricDepthBin,
};

const D: f64 = 80e-9;

fn knots(near_face: bool) -> Vec<PhotoelectricKnot> {
    let names = if near_face {
        [
            (200.0, "near_q_short", "near_f_short"),
            (2000.0, "near_q_long", "near_f_long"),
        ]
    } else {
        [
            (200.0, "far_q_short", "far_f_short"),
            (2000.0, "far_q_long", "far_f_long"),
        ]
    };
    names
        .into_iter()
        .map(|(wavelength_nm, q, f)| PhotoelectricKnot {
            wavelength_nm,
            absorbed_quantum_yield_property: PropertyKey::new(q, Dims::NONE),
            mean_excess_fraction_property: PropertyKey::new(f, Dims::NONE),
        })
        .collect()
}

fn profile(
    original: &ResolvedMaterialStatePoint,
    bins: &[(f64, bool)],
) -> DepthPhotoelectricMaterial {
    let bins: Vec<_> = bins
        .iter()
        .map(|&(end_fraction, near_face)| PhotoelectricDepthBin {
            end_fraction,
            knots: knots(near_face),
        })
        .collect();
    DepthPhotoelectricMaterial::try_from_resolved(
        original.clone(),
        PropertyKey::new("removal_work", ENERGY),
        PropertyKey::new("emission_film_thickness", Length::DIMS),
        &bins,
        1800.0,
    )
    .unwrap()
}

fn fixture() -> (
    ResolvedConductorState,
    PhotoelectricMaterial,
    DepthPhotoelectricMaterial,
) {
    let basis = ElectricalPropertyBasis::Conductivity;
    let original = source(
        &[
            (EPS_NAME, EPS, 2.0 * EPSILON),
            (MU_NAME, MU, PERMEABILITY),
            (basis.property_name(), basis.dimensions(), 1e5),
            ("removal_work", ENERGY, WORK),
            ("emission_film_thickness", Length::DIMS, D),
            ("near_q_short", Dims::NONE, 0.8),
            ("near_q_long", Dims::NONE, 0.8),
            ("near_f_short", Dims::NONE, 0.2),
            ("near_f_long", Dims::NONE, 0.2),
            ("far_q_short", Dims::NONE, 0.1),
            ("far_q_long", Dims::NONE, 0.1),
            ("far_f_short", Dims::NONE, 0.7),
            ("far_f_long", Dims::NONE, 0.7),
        ],
        0,
    );
    let conductor =
        ResolvedConductorState::try_from_material_with_basis(original.clone(), basis).unwrap();
    let whole = PhotoelectricMaterial::try_from_resolved(
        original.clone(),
        PropertyKey::new("removal_work", ENERGY),
        PropertyKey::new("emission_film_thickness", Length::DIMS),
        &knots(true),
        1800.0,
    )
    .unwrap();
    let depth = profile(&original, &[(0.5, true), (1.0, false)]);
    (conductor, whole, depth)
}

#[test]
fn g1_depth_escape_integrates_native_coherent_electric_work_not_volume_average() {
    let exterior = vacuum();
    let (metal, _, profile) = fixture();
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: D,
        }],
    )
    .unwrap();
    for polarization in [Polarization::Te, Polarization::Tm] {
        let illumination = stack.prepare_illumination(0.8, polarization).unwrap();
        let owner = PreparedLayerDepthPhotoelectric::try_new(
            &illumination,
            0,
            &profile,
            EscapeFace::EntranceVacuum,
            0.002,
        )
        .unwrap();
        let mut efficiencies = Vec::new();
        for sign in [-1.0, 1.0] {
            let input = [C64::from_re(100.0), C64::new(40.0 * sign, 20.0 * sign)];
            let result = owner.evaluate(input, || false).unwrap();
            let mut reference = [0.0; 2];
            let points = 512;
            let dz = D / points as f64;
            for i in 0..points {
                let field = illumination
                    .prepare_point(0, (i as f64 + 0.5) * dz)
                    .unwrap()
                    .evaluate(input)
                    .unwrap();
                // Independent local Ohmic work, including all three E components.
                let loss = 0.5 * 1e5 * field.electric_v_m.iter().map(|e| e.norm_sq()).sum::<f64>();
                reference[usize::from(i >= points / 2)] += loss * dz * 0.002;
            }
            let photon = 6.626_070_15e-34 * 299_792_458.0 / 250e-9;
            let expected_rate = (0.8 * reference[0] + 0.1 * reference[1]) / photon;
            assert!(
                (result.totals().emitted_electron_rate_s - expected_rate).abs()
                    < 2e-5 * expected_rate
            );
            for (bin, expected) in result.bins().zip(reference) {
                assert!((bin.photoelectric().absorbed_power_w - expected).abs() < 2e-5 * expected);
                near(
                    bin.photoelectric().absorbed_power_w,
                    bin.photoelectric().heat_power_w
                        + bin.photoelectric().removal_power_w
                        + bin.photoelectric().kinetic_power_w,
                );
            }
            efficiencies.push(
                result.totals().emitted_electron_rate_s / result.totals().absorbed_photon_rate_s,
            );
        }
        assert!(
            (efficiencies[0] - efficiencies[1]).abs() > 0.01,
            "phase-dependent absorption must change the profile-weighted escape probability"
        );
    }
}

#[test]
fn g3_uniform_profile_refinement_recovers_the_original_whole_film_law() {
    let exterior = vacuum();
    let (metal, whole_source, _) = fixture();
    let uniform = profile(
        whole_source.resolved(),
        &[(0.25, true), (0.75, true), (1.0, true)],
    );
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: D,
        }],
    )
    .unwrap();
    let illumination = stack.prepare_illumination(0.67, Polarization::Tm).unwrap();
    for face in [EscapeFace::EntranceVacuum, EscapeFace::ExitVacuum] {
        let whole =
            PreparedLayerPhotoelectric::try_new(&illumination, 0, &whole_source, face, 0.01)
                .unwrap();
        let depth =
            PreparedLayerDepthPhotoelectric::try_new(&illumination, 0, &uniform, face, 0.01)
                .unwrap();
        let input = [C64::new(30.0, 70.0), C64::new(20.0, -25.0)];
        let a = whole.evaluate(input).unwrap();
        let b = depth.evaluate(input, || false).unwrap();
        assert_eq!(b.bin_count(), 3);
        assert_eq!(b.bins().count(), 3);
        assert!(b.bin(3).is_none());
        for (actual, expected) in [
            (
                b.totals().absorbed_power_w,
                a.photoelectric().absorbed_power_w,
            ),
            (
                b.totals().emitted_electron_rate_s,
                a.photoelectric().emitted_electron_rate_s,
            ),
            (
                b.totals().surface_charging_current_a,
                a.photoelectric().surface_charging_current_a,
            ),
            (
                b.totals().kinetic_power_w,
                a.photoelectric().kinetic_power_w,
            ),
            (
                b.totals().removal_power_w,
                a.photoelectric().removal_power_w,
            ),
            (b.totals().heat_power_w, a.photoelectric().heat_power_w),
        ] {
            near(actual, expected);
        }
        assert!(std::ptr::eq(b.material(), &uniform));
        assert!(std::ptr::eq(b.illumination(), &illumination));
    }
}

#[test]
fn g3_escape_face_reverses_physical_depths_and_preserves_absorbed_power() {
    let exterior = vacuum();
    let (metal, _, profile) = fixture();
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: D,
        }],
    )
    .unwrap();
    let illumination = stack.prepare_illumination(1.0, Polarization::Te).unwrap();
    let entrance = PreparedLayerDepthPhotoelectric::try_new(
        &illumination,
        0,
        &profile,
        EscapeFace::EntranceVacuum,
        1.0,
    )
    .unwrap();
    let exit = PreparedLayerDepthPhotoelectric::try_new(
        &illumination,
        0,
        &profile,
        EscapeFace::ExitVacuum,
        1.0,
    )
    .unwrap();
    assert_eq!(entrance.bin_depth_interval_m(0), Some([0.0, D / 2.0]));
    assert_eq!(exit.bin_depth_interval_m(0), Some([D / 2.0, D]));
    let left = [C64::from_re(100.0), C64::ZERO];
    let a = entrance.evaluate(left, || false).unwrap();
    let b = exit.evaluate(left, || false).unwrap();
    let reversed = exit.evaluate([left[1], left[0]], || false).unwrap();
    near(a.totals().absorbed_power_w, b.totals().absorbed_power_w);
    near(
        a.totals().emitted_electron_rate_s,
        reversed.totals().emitted_electron_rate_s,
    );
    assert!(a.totals().emitted_electron_rate_s > 1.2 * b.totals().emitted_electron_rate_s);
    let original = illumination
        .prepare_depth_interval(0, 0.0, D)
        .unwrap()
        .evaluate(left)
        .unwrap();
    near(a.totals().absorbed_power_w, original);
    near(
        a.bins().map(|bin| bin.absorbed_power_w_m2()).sum(),
        original,
    );
    assert_eq!(
        a.totals().surface_charging_current_a,
        -a.totals().emitted_electron_current_a
    );
}

#[test]
fn g0_total_bin_mode_budget_cancellation_and_collapsed_support_refuse() {
    let exterior = vacuum();
    let (metal, source, profile_source) = fixture();
    let stack = PlanarDielectricLayers::try_new(
        &exterior,
        DielectricResponse::Lossless(&exterior),
        omega(250e-9),
        &[DielectricLayer {
            response: DielectricResponse::Ohmic(&metal),
            thickness_m: D,
        }],
    )
    .unwrap();
    let illumination = stack.prepare_illumination(1.0, Polarization::Te).unwrap();
    let owner = PreparedLayerDepthPhotoelectric::try_new(
        &illumination,
        0,
        &profile_source,
        EscapeFace::EntranceVacuum,
        1.0,
    )
    .unwrap();
    let factors = [[C64::ONE, C64::ZERO], [C64::ZERO, C64::ONE]];
    let exact = IlluminationBudget {
        maximum_modes: 2,
        maximum_layer_mode_pairs: 4,
    };
    assert!(
        owner
            .evaluate_coherence_modes(
                &factors,
                IlluminationBudget {
                    maximum_layer_mode_pairs: 3,
                    ..exact
                },
                || false
            )
            .is_err()
    );
    let before = owner
        .evaluate_coherence_modes(&factors, exact, || false)
        .unwrap();
    let independent = owner
        .evaluate_independent_sources([C64::ONE; 2], || false)
        .unwrap();
    assert_eq!(before.totals(), independent.totals());
    let mut checkpoints = 0;
    assert!(
        owner
            .evaluate_coherence_modes(&factors, exact, || {
                checkpoints += 1;
                checkpoints == 5
            })
            .is_err()
    );
    let after = owner
        .evaluate_coherence_modes(&factors, exact, || false)
        .unwrap();
    assert_eq!(before.totals(), after.totals());
    let dark = owner.evaluate([C64::ZERO; 2], || false).unwrap();
    assert_eq!(dark.totals().absorbed_power_w, 0.0);
    assert_eq!(dark.totals().emitted_electron_rate_s, 0.0);
    let collapsed = profile(source.resolved(), &[(f64::from_bits(1), true), (1.0, true)]);
    for face in [EscapeFace::EntranceVacuum, EscapeFace::ExitVacuum] {
        assert!(
            PreparedLayerDepthPhotoelectric::try_new(&illumination, 0, &collapsed, face, 1.0)
                .is_err()
        );
    }
}
