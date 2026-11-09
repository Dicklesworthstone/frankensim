//! Independent native field controls for analytic arbitrary-axis boundary jets.
use super::super::{Incidence, solve, solve_with_witness};
use super::*;

const MU: f64 = 1.256_637_062_12e-6;
const OMEGA: f64 = core::f64::consts::TAU * C0 / 550e-9;
const IMPEDANCE: f64 = MU * C0;

fn native() -> NativeInput {
    NativeInput {
        omega: OMEGA,
        ni: 1.2,
        no: 1.5,
        ne: 1.9,
        mui: MU,
        mut_: 1.1 * MU,
        axis: [1.0, 2.0, 3.0],
        incidence: Incidence::TangentialIndex(-0.7),
        electric: [C64::new(0.6, 0.2), C64::new(-0.3, 0.4)],
    }
}

// Each original n(s)=n0*(1+a*(s-1)+b*(s-1)^2/2), s=omega/omega0.
// The relative epsilon jets and principal K jets are independently analytic.
fn source_jets(data: NativeInput, laws: [[f64; 2]; 3]) -> ([[f64; 2]; 3], [[f64; 2]; 2]) {
    let relative = laws.map(|[a, b]| [2.0 * a, 2.0 * (a * a + b)]);
    let principal = core::array::from_fn(|j| {
        let n = [data.no, data.ne][j];
        let [a, b] = laws[j + 1];
        [n * (1.0 + a) / C0, n * (2.0 * a + b) / C0]
    });
    (relative, principal)
}

fn evaluated(
    data: NativeInput,
    laws: [[f64; 2]; 3],
    constraint: IncidenceFrequencyConstraint,
) -> (OrientedCrystalResponse, BoundaryDerivatives) {
    let (boundary, witness) = solve_with_witness(data).unwrap();
    let (relative, principal) = source_jets(data, laws);
    let jet = derivatives(Input {
        original: data,
        boundary: &boundary,
        witness: &witness,
        relative_permittivity_derivatives: relative,
        principal_wavenumber_derivatives: principal,
        constraint,
    })
    .unwrap();
    (boundary, jet)
}

fn vectors(boundary: OrientedCrystalResponse) -> [Vector; 8] {
    [
        boundary.incident_electric_v_m,
        boundary.incident_magnetic_a_m,
        boundary.reflected_electric_v_m,
        boundary.reflected_magnetic_a_m,
        boundary.transmitted[0].electric_v_m,
        boundary.transmitted[0].magnetic_a_m,
        boundary.transmitted[1].electric_v_m,
        boundary.transmitted[1].magnetic_a_m,
    ]
}
fn jets(jet: BoundaryDerivatives) -> [[Vector; 2]; 8] {
    [
        jet.incident.electric,
        jet.incident.magnetic,
        jet.reflected.electric,
        jet.reflected.magnetic,
        jet.transmitted[0].electric,
        jet.transmitted[0].magnetic,
        jet.transmitted[1].electric,
        jet.transmitted[1].magnetic,
    ]
}
fn near(actual: C64, expected: C64, tolerance: f64) {
    assert!(
        (actual - expected).abs() <= tolerance,
        "{actual:?} != {expected:?}, tolerance {tolerance:e}"
    );
}
fn energy_norm(vector: Vector, magnetic: bool) -> f64 {
    vector
        .into_iter()
        .map(|z| z.abs() * if magnetic { IMPEDANCE } else { 1.0 })
        .fold(0.0, f64::max)
}

fn neighboring(
    data: NativeInput,
    laws: [[f64; 2]; 3],
    constraint: IncidenceFrequencyConstraint,
    offset: f64,
) -> OrientedCrystalResponse {
    let s = 1.0 + offset;
    let indices = [data.ni, data.no, data.ne];
    let changed = core::array::from_fn::<_, 3, _>(|j| {
        let [a, b] = laws[j];
        indices[j] * (1.0 + a * offset + 0.5 * b * offset * offset)
    });
    let incidence = match (constraint, data.incidence) {
        (IncidenceFrequencyConstraint::FixedAngle, Incidence::Cosine(c)) => Incidence::Cosine(c),
        (IncidenceFrequencyConstraint::FixedAngle, Incidence::TangentialIndex(p)) => {
            Incidence::TangentialIndex(p * changed[0] / data.ni)
        }
        (IncidenceFrequencyConstraint::FixedTransverseWavenumber, Incidence::Cosine(c)) => {
            Incidence::TangentialIndex(data.ni * ((1.0 - c) * (1.0 + c)).sqrt() / s)
        }
        (
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            Incidence::TangentialIndex(p),
        ) => Incidence::TangentialIndex(p / s),
    };
    solve(NativeInput {
        omega: data.omega * s,
        ni: changed[0],
        no: changed[1],
        ne: changed[2],
        incidence,
        ..data
    })
    .unwrap()
}

fn frequency_control(
    data: NativeInput,
    laws: [[f64; 2]; 3],
    constraint: IncidenceFrequencyConstraint,
) {
    let (_, jet) = evaluated(data, laws, constraint);
    let actual = jets(jet);
    let h = 5e-4;
    let samples = [-2.0 * h, -h, 0.0, h, 2.0 * h]
        .map(|delta| vectors(neighboring(data, laws, constraint, delta)));
    for field in 0..8 {
        let unit = if field % 2 == 0 { 1.0 } else { IMPEDANCE };
        for component in 0..3 {
            let [a, b, c, d, e] = samples.map(|field_values| field_values[field][component]);
            let first = (a - b.scale(8.0) + d.scale(8.0) - e).scale(1.0 / (12.0 * h));
            let second = (-a + b.scale(16.0) - c.scale(30.0) + d.scale(16.0) - e)
                .scale(1.0 / (12.0 * h * h));
            near(
                actual[field][0][component].scale(unit),
                first.scale(unit),
                3e-8,
            );
            near(
                actual[field][1][component].scale(unit),
                second.scale(unit),
                8e-7,
            );
        }
    }
    // Differentiated tangential Maxwell continuity in physical field units,
    // independent of the eigenbasis normalization and coefficient convention.
    for order in 0..2 {
        for j in 0..2 {
            near(
                jet.incident.electric[order][j] + jet.reflected.electric[order][j],
                jet.transmitted[0].electric[order][j] + jet.transmitted[1].electric[order][j],
                2e-13,
            );
            near(
                (jet.incident.magnetic[order][j] + jet.reflected.magnetic[order][j])
                    .scale(IMPEDANCE),
                (jet.transmitted[0].magnetic[order][j] + jet.transmitted[1].magnetic[order][j])
                    .scale(IMPEDANCE),
                2e-13,
            );
        }
    }
}

#[test]
fn nondispersive_fixed_angle_has_exactly_zero_vector_jets_including_evanescent_modes() {
    for data in [
        native(),
        NativeInput {
            ni: 1.9,
            no: 1.2,
            ne: 1.4,
            incidence: Incidence::TangentialIndex(1.7),
            ..native()
        },
    ] {
        let (_, jet) = evaluated(
            data,
            [[0.0; 2]; 3],
            IncidenceFrequencyConstraint::FixedAngle,
        );
        assert_eq!(jets(jet), [[[C64::ZERO; 3]; 2]; 8]);
    }
}

#[test]
fn common_tiny_index_slopes_keep_physical_magnetic_jets_when_bulk_group_slopes_round_away() {
    let alpha = 1e-20;
    let beta = -3e-20;
    let data = native();
    let laws = [[alpha, beta]; 3];
    let (_, principal) = source_jets(data, laws);
    assert_eq!(principal[0][0], data.no / C0);
    let (boundary, jet) = evaluated(data, laws, IncidenceFrequencyConstraint::FixedAngle);
    let original = vectors(boundary);
    for (field, orders) in jets(jet).into_iter().enumerate() {
        for (order, factor) in [alpha, beta].into_iter().enumerate() {
            for j in 0..3 {
                if field % 2 == 0 {
                    assert_eq!(orders[order][j], C64::ZERO);
                } else {
                    let expected = original[field][j].scale(factor);
                    near(orders[order][j], expected, expected.abs() * 3e-15);
                }
            }
        }
    }
}

#[test]
fn distinct_tiny_ordinary_slope_survives_the_raw_normal_index_seam() {
    // At normal incidence, this polarization is perpendicular to the optic
    // axis's transverse projection. Its independent scalar admittance is no/mu.
    // ne has a different local law, so the common-index identity cannot apply.
    let data = NativeInput {
        ni: 1.0,
        no: 1.5,
        ne: 1.9,
        mui: MU,
        mut_: MU,
        incidence: Incidence::Cosine(1.0),
        electric: [C64::from_re(-2.0), C64::ONE],
        ..native()
    };
    let [alpha, beta] = [1e-20, -3e-20];
    let laws = [[0.0; 2], [alpha, beta], [0.0; 2]];
    let (_, principal) = source_jets(data, laws);
    assert_eq!(principal[0][0], data.no / C0);
    let (_, jet) = evaluated(data, laws, IncidenceFrequencyConstraint::FixedAngle);
    let sum = data.ni + data.no;
    let first = -2.0 * data.ni * data.no * alpha / (sum * sum);
    let second = -2.0 * data.ni * data.no * beta / (sum * sum)
        + 4.0 * data.ni * data.no * data.no * alpha * alpha / (sum * sum * sum);
    for (order, derivative) in [first, second].into_iter().enumerate() {
        for component in 0..2 {
            let expected = data.electric[component].scale(derivative);
            near(
                jet.reflected.electric[order][component],
                expected,
                expected.abs() * 3e-12,
            );
            near(
                jet.transmitted[0].electric[order][component]
                    + jet.transmitted[1].electric[order][component],
                expected,
                expected.abs() * 3e-12,
            );
        }
    }
}

#[test]
fn arbitrary_axis_complex_fields_match_original_frequency_laws_at_both_constraints() {
    let laws = [[0.03, -0.02], [-0.04, 0.06], [0.07, -0.03]];
    for data in [
        native(),
        NativeInput {
            ni: 1.9,
            no: 1.2,
            ne: 1.4,
            incidence: Incidence::TangentialIndex(1.7),
            ..native()
        },
    ] {
        for constraint in [
            IncidenceFrequencyConstraint::FixedAngle,
            IncidenceFrequencyConstraint::FixedTransverseWavenumber,
        ] {
            frequency_control(data, laws, constraint);
        }
    }
}

#[test]
fn stationary_bulk_wavenumbers_do_not_erase_constitutive_vector_derivatives() {
    // A local n proportional to 1/omega has K'=K''=0. At fixed kx all
    // directions and electric fields stay fixed, but physical H scales 1/omega.
    let data = native();
    let laws = [[-1.0, 2.0]; 3];
    assert_eq!(source_jets(data, laws).1, [[0.0; 2]; 2]);
    let (boundary, jet) = evaluated(
        data,
        laws,
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    );
    let original = vectors(boundary);
    for (field, orders) in jets(jet).into_iter().enumerate() {
        for (order, factor) in [-1.0, 2.0].into_iter().enumerate() {
            let unit = if field % 2 == 0 { 1.0 } else { IMPEDANCE };
            for j in 0..3 {
                let expected = if field % 2 == 0 {
                    C64::ZERO
                } else {
                    original[field][j].scale(factor)
                };
                near(orders[order][j].scale(unit), expected.scale(unit), 3e-13);
            }
        }
    }
}

#[test]
fn an_unexcited_native_mode_can_acquire_a_nonzero_coupled_frequency_field() {
    let mut data = native();
    data.incidence = Incidence::TangentialIndex(0.7);
    let (_, witness) = solve_with_witness(data).unwrap();
    let basis = witness.basis;
    let ordinary = basis.modes[0];
    // Choose the incident tangential field that excites just this one original
    // native mode. The source polarization is then held fixed in frequency.
    data.electric = [
        (ordinary.e[0] + ordinary.h[1].scale(basis.qi / (basis.ratio_mu * basis.i * basis.i)))
            .scale(0.5),
        (ordinary.e[1] - ordinary.h[0].scale(1.0 / (basis.ratio_mu * basis.qi))).scale(0.5),
    ];
    let laws = [[0.03, -0.02], [-0.04, 0.06], [0.07, -0.03]];
    let (boundary, jet) = evaluated(data, laws, IncidenceFrequencyConstraint::FixedAngle);
    assert!(energy_norm(boundary.transmitted[1].electric_v_m, false) < 2e-14);
    assert!(energy_norm(jet.transmitted[1].electric[0], false) > 1e-5);
    frequency_control(data, laws, IncidenceFrequencyConstraint::FixedAngle);
}

#[test]
fn exact_critical_or_spectrally_split_degenerate_carriers_refuse_vector_jets() {
    for data in [
        NativeInput {
            ni: 2.0,
            no: 1.0,
            ne: 1.7,
            axis: [0.3, 0.8, 0.6],
            incidence: Incidence::TangentialIndex(1.0),
            ..native()
        },
        NativeInput {
            ni: 6.0,
            no: 1.0,
            ne: 7.0,
            axis: [1.0, 0.0, 1.0],
            incidence: Incidence::TangentialIndex(5.0),
            ..native()
        },
        NativeInput {
            no: 1.5,
            ne: 1.5,
            ..native()
        },
    ] {
        let (boundary, witness) = solve_with_witness(data).unwrap();
        let (relative, principal) = source_jets(data, [[0.01, 0.02], [0.03, 0.04], [0.03, 0.05]]);
        assert!(
            derivatives(Input {
                original: data,
                boundary: &boundary,
                witness: &witness,
                relative_permittivity_derivatives: relative,
                principal_wavenumber_derivatives: principal,
                constraint: IncidenceFrequencyConstraint::FixedTransverseWavenumber
            })
            .is_err()
        );
    }
    // Equal phase indices in the original principal plane remain separate
    // admitted TE/TM spectral branches, including unequal second-order laws.
    let data = NativeInput {
        no: 1.5,
        ne: 1.5,
        axis: [1.0, 0.0, 1.0],
        ..native()
    };
    evaluated(
        data,
        [[0.01, 0.02], [0.03, 0.04], [0.03, 0.05]],
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    );
}

#[test]
fn matched_grazing_scalar_and_principal_te_modes_keep_finite_physical_field_jets() {
    for extraordinary in [1.5, 2.0] {
        for cosine in [0.8, 1e-10, 1e-20, 1e-100] {
            let data = NativeInput {
                ni: 1.5,
                no: 1.5,
                ne: extraordinary,
                mui: MU,
                mut_: MU,
                axis: if extraordinary == 1.5 {
                    [1.0, 2.0, 3.0]
                } else {
                    [1.0, 0.0, 0.0]
                },
                incidence: Incidence::Cosine(cosine),
                electric: [C64::ZERO, C64::ONE],
                ..native()
            };
            // The principal TE mode remains exactly matched even when the
            // independent extraordinary material has different local jets.
            let laws = [
                [0.0; 2],
                [0.0; 2],
                if extraordinary == 1.5 {
                    [0.0; 2]
                } else {
                    [0.07, -0.03]
                },
            ];
            let (boundary, jet) = evaluated(
                data,
                laws,
                IncidenceFrequencyConstraint::FixedTransverseWavenumber,
            );
            assert_eq!(boundary.reflected_electric_v_m, [C64::ZERO; 3]);
            assert_eq!(jet.reflected.electric, [[C64::ZERO; 3]; 2]);
            assert_eq!(jet.reflected.magnetic, [[C64::ZERO; 3]; 2]);
            assert_eq!(jet.transmitted[0].electric, [[C64::ZERO; 3]; 2]);
            assert_eq!(jet.transmitted[1].electric, [[C64::ZERO; 3]; 2]);
            assert_eq!(jet.transmitted[1].magnetic, [[C64::ZERO; 3]; 2]);
            let ratio = (((1.0 - cosine) * (1.0 + cosine)) / cosine) / cosine;
            let original = boundary.incident_magnetic_a_m;
            let first = original[0].scale(ratio);
            // Weight by the native Hx before the second grazing factor.
            // ratio² alone overflows at 1e-100 although physical H'' is finite.
            let second = first.scale(-3.0 - ratio);
            for (order, expected) in [first, second].into_iter().enumerate() {
                assert!(expected.re.is_finite() && expected.im.is_finite());
                near(
                    jet.incident.magnetic[order][0],
                    expected,
                    expected.abs() * 3e-13,
                );
                near(
                    jet.transmitted[0].magnetic[order][0],
                    expected,
                    expected.abs() * 3e-13,
                );
            }
            near(
                jet.incident.magnetic[0][2],
                -original[2],
                original[2].abs() * 2e-14,
            );
            near(
                jet.incident.magnetic[1][2],
                original[2].scale(2.0),
                original[2].abs() * 2e-14,
            );
        }
    }
    // Linearity keeps the same matched TE result when a simultaneous TM
    // illumination excites the distinct extraordinary spectral branch.
    let cosine = 1e-20;
    let data = NativeInput {
        ni: 1.5,
        no: 1.5,
        ne: 2.0,
        mui: MU,
        mut_: MU,
        axis: [1.0, 0.0, 0.0],
        incidence: Incidence::Cosine(cosine),
        electric: [C64::new(0.2, 0.1), C64::ONE],
        ..native()
    };
    let (boundary, jet) = evaluated(
        data,
        [[0.0; 2], [0.0; 2], [0.07, -0.03]],
        IncidenceFrequencyConstraint::FixedTransverseWavenumber,
    );
    assert_eq!(jet.transmitted[0].electric, [[C64::ZERO; 3]; 2]);
    let ratio = (((1.0 - cosine) * (1.0 + cosine)) / cosine) / cosine;
    let first = boundary.incident_magnetic_a_m[0].scale(ratio);
    let second = first.scale(-3.0 - ratio);
    for (order, expected) in [first, second].into_iter().enumerate() {
        assert_eq!(jet.reflected.electric[order][1], C64::ZERO);
        assert_eq!(jet.reflected.magnetic[order][0], C64::ZERO);
        assert_eq!(jet.reflected.magnetic[order][2], C64::ZERO);
        near(
            jet.transmitted[0].magnetic[order][0],
            expected,
            expected.abs() * 3e-13,
        );
    }
    // In this exact axial section q_e=ne*cos(theta_i), so the full TM
    // tangential reflection is (ni-ne)/(ni+ne). Its material derivatives
    // remain finite even while longitudinal field derivatives grow at grazing.
    let sum = data.ni + data.ne;
    let alpha = 0.07;
    let beta = -0.03;
    let first = -2.0 * data.ni * data.ne * alpha / (sum * sum);
    let second = -2.0 * data.ni * data.ne * beta / (sum * sum)
        + 4.0 * data.ni * data.ne * data.ne * alpha * alpha / (sum * sum * sum);
    for (order, factor) in [first, second].into_iter().enumerate() {
        let expected = data.electric[0].scale(factor);
        near(
            jet.reflected.electric[order][0],
            expected,
            expected.abs() * 3e-12,
        );
        near(
            jet.transmitted[1].electric[order][0],
            expected,
            expected.abs() * 3e-12,
        );
    }
}
