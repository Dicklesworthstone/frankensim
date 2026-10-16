//! IRP M1.2-W1a static obstacle course. No recording fit or DKT fallback.
//!
//! References: MacNeal & Harder, FEAD 1 (1985) 3–20, and Belytschko et al.,
//! CMAME 51 (1985) 221–258. Author reproduction: Ko, Lee, Lee & Bathe,
//! Computers & Structures 193 (2017) 187–206, doi:10.1016/j.compstruc.2017.08.003,
//! Figs.9/11/17/19 and Tables 8/10/15/17, available at
//! https://web.mit.edu/kjb/www/Principal_Publications/Performance_of_the_MITC3%2B_and_MITC4%2B_shell_elements_in_widely_used_benchmark_problems.pdf
//!
//! The bead gates use roof .3024 (not original .3086), cylinder 1.8248e-5,
//! cut-out hemisphere .094 (MacNeal–Harder), and full hemisphere .0924.
//! MacNeal–Harder Table 5(b), p.11, doi:10.1016/0168-874X(85)90003-4,
//! declares .0940 as its cut-out normalization and discusses .0924 without
//! the hole; these are reference choices rather than exact shell solutions.
//! Ko et al. 2017 Tables 15/16 (pp.202/203) normalize cut-out results by .093;
//! Table 17 (p.204) normalizes full-hemisphere results by .0924. Both cut-out
//! normalizations remain gated. The .47 ruling of 2026-10-16 replaces error-to-
//! reference monotonicity for hemispheres with contracting sequence changes
//! and a signed Richardson limit, allowing the observed finest-grid reversal.
//! These are pointwise benchmark bands, not rigorous error enclosures.

#[path = "shell_obstacle_course/fixtures.rs"]
mod fixtures;
use fixtures::{Case, Fixture};
use fs_material::sheet::SheetMatrix;
use fs_plate::shell::mitc::{
    MitcBudget, ShearShellSection, ShellMass, assemble_mitc3p, solve_static,
};
use fs_sparse::{Coo, DirectOrdering, LdltOptions, SymbolicLdlt};

fn response(f: &Fixture, alpha: f64) -> f64 {
    let model = assemble_mitc3p(
        &f.mesh,
        &f.directors,
        &f.sections,
        &f.constraints,
        ShellMass::Consistent,
        MitcBudget {
            max_nodes: f.mesh.nodes.len(),
            max_triangles: f.mesh.tris.len(),
        },
        alpha,
        || false,
    )
    .unwrap();
    let u = solve_static(&model, &f.load).unwrap();
    f.sign * u[6 * f.observation + f.component]
}

/// Three-point Richardson/Aitken extrapolation of u_h=L+c*q^j. A negative
/// signed q admits alternating convergence; |q|<1 still requires the last
/// two step sizes to decrease strictly. No assumed continuum order is fitted
/// to a literature normalization.
fn sequence_limit(values: &[f64]) -> Option<(f64, f64)> {
    let [x, y, z]: [f64; 3] = values
        .get(values.len().checked_sub(3)?..)?
        .try_into()
        .ok()?;
    let previous = y - x;
    let latest = z - y;
    let q = latest / previous;
    if !q.is_finite() || q == 0.0 || q.abs() >= 1.0 {
        return None;
    }
    let limit = z + latest * q / (1.0 - q);
    let order = -fs_math::det::ln(q.abs()) / fs_math::det::ln(2.0);
    (limit.is_finite() && limit > 0.0).then_some((limit, order))
}

fn hemisphere_reference_band(case: Case, value: f64, tolerance: f64) -> bool {
    let references: &[f64] = match case {
        Case::Hole => &[0.094, 0.093],
        Case::Full => &[0.0924],
        _ => unreachable!("hemisphere reference gate called for another geometry"),
    };
    value.is_finite()
        && references
            .iter()
            .all(|r| (value / r - 1.0).abs() <= tolerance)
}

#[test]
fn g0_signed_richardson_recovers_monotone_and_alternating_limits() {
    for q in [0.5, -0.5] {
        let mut correction = 0.5;
        let values: Vec<_> = (0..4)
            .map(|_| {
                let u = 2.0 + correction;
                correction *= q;
                u
            })
            .collect();
        let (limit, order) = sequence_limit(&values).unwrap();
        assert!((limit - 2.0).abs() < 1e-14);
        assert!((order - 1.0).abs() < 1e-14);
    }
    for invalid in [
        &[1.0, 1.0, 1.0][..],
        &[1.0, 2.0, 4.0],
        &[1.0, 2.0],
        &[1.0, f64::NAN, 2.0],
    ] {
        assert!(sequence_limit(invalid).is_none());
    }
}

/// Change only D in these explicitly isotropic, rho=1, angle=0 fixtures.
/// A, transverse shear, mass, geometry, loads and supports retain their inputs.
fn plant_bending_error(f: &mut Fixture) {
    let shear = f.e / (2.0 * (1.0 + f.nu));
    for section in &mut f.sections {
        let (a, d) = section.principal_matrices();
        let d = d.coefficients();
        *section = ShearShellSection::new(
            *a,
            SheetMatrix::orthotropic(1.05 * d[0], 1.05 * d[4], 1.05 * d[1], 1.05 * d[8]).unwrap(),
            [shear * f.h, 0.0, shear * f.h],
            5.0 / 6.0,
            f.h,
            f.h * f.h * f.h / 12.0,
            0.0,
        )
        .unwrap();
    }
}

fn ladder(case: Case, distorted: bool, tolerance: f64, test: &str) {
    let started = std::time::Instant::now();
    let mut rungs = case.rungs().to_vec();
    // The recorded distorted-hole N=8,16,32 steps do not contract. Retain all
    // original rungs and resolve that row at N=64 rather than erase its gate.
    if distorted && matches!(case, Case::Hole) {
        rungs.push(64);
    }
    let reference = case.reference();
    let mut values = Vec::new();
    let mut errors = Vec::new();
    let hemisphere = matches!(case, Case::Hole | Case::Full);
    let mut planted_values = Vec::new();
    for &n in &rungs {
        let mut fixture = fixtures::build(case, n, distorted);
        let observed = response(&fixture, 1e-3);
        assert!(observed.is_finite() && observed > 0.0);
        values.push(observed);
        errors.push((observed / reference - 1.0).abs());
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1a\",\"test\":\"{test}\",\"tier\":\"G2\",\"inputs\":{{\"benchmark\":\"{}\",\"element\":\"mitc3+\",\"mass\":\"consistent\",\"n\":{n},\"distorted\":{distorted},\"distortion_key\":{}}},\"observed\":{{\"reference\":{reference},\"displacement\":{observed},\"relative_error\":{}}},\"tolerance\":{tolerance},\"elapsed_s\":{}}}",
            case.name(),
            fixtures::DISTORTION_KEY,
            errors.last().unwrap(),
            started.elapsed().as_secs_f64()
        );
        if hemisphere && !distorted {
            plant_bending_error(&mut fixture);
            planted_values.push(response(&fixture, 1e-3));
        }
    }
    let hs: Vec<_> = rungs.iter().map(|&n| 1.0 / n as f64).collect();
    let fit = fs_mms::fit_order(&fs_mms::RefinementLadder::new(hs, errors.clone()).unwrap());
    let changes: Vec<_> = values.windows(2).map(|p| (p[1] - p[0]).abs()).collect();
    let error_monotone = errors[errors.len() - 3..].windows(2).all(|p| p[1] < p[0]);
    let extrapolated = if hemisphere {
        sequence_limit(&values)
    } else {
        None
    };
    let monotone = if hemisphere {
        changes[changes.len() - 1] < changes[changes.len() - 2] && extrapolated.is_some()
    } else {
        error_monotone
    };
    let best_reference = if matches!(case, Case::Hole) {
        0.093
    } else {
        reference
    };
    let finest_pass = (values.last().unwrap() / best_reference - 1.0).abs() <= tolerance;
    let limit_pass = !hemisphere
        || extrapolated.is_some_and(|(u, _)| hemisphere_reference_band(case, u, tolerance));
    let passed = finest_pass && monotone && limit_pass;
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1a\",\"test\":\"{test}\",\"tier\":\"G2\",\"inputs\":{{\"benchmark\":\"{}\",\"model\":\"{}\",\"element\":\"mitc3+\",\"mass\":\"consistent\",\"distorted\":{distorted},\"distortion_key\":{}}},\"provenance\":\"MacNeal-Harder 1985; Belytschko 1985; Ko et al. 2017\",\"rungs\":{{\"n_per_side\":{rungs:?}}},\"observed\":{{\"reference\":{reference},\"u_h\":{values:?},\"rel_err\":{errors:?},\"reference_error_order\":{},\"reference_error_monotone\":{error_monotone},\"last_three_converge\":{monotone}}},\"alternate_reference\":{},\"tolerance\":{tolerance},\"verdict\":{},\"elapsed_s\":{}}}",
        case.name(),
        if matches!(case, Case::Cylinder) {
            "one-eighth"
        } else {
            "quarter"
        },
        fixtures::DISTORTION_KEY,
        fit.observed,
        match case {
            Case::Roof => 0.3086,
            Case::Hole => 0.093,
            _ => reference,
        },
        passed,
        started.elapsed().as_secs_f64()
    );
    assert!(
        finest_pass,
        "{} displacement {values:?}; reference {reference}; relative errors {errors:?}",
        case.name()
    );
    assert!(
        monotone,
        "{} last-three-rung convergence failed: values {values:?}, changes {changes:?}",
        case.name()
    );
    if hemisphere {
        let (limit, order) = extrapolated.expect("last triple must contract");
        println!(
            "{{\"suite\":\"irp-s1\",\"test\":\"{test}\",\"stage\":\"hemisphere-sequence-limit\",\"benchmark\":\"{}\",\"n_per_side\":{rungs:?},\"changes\":{changes:?},\"signed_richardson_limit\":{limit},\"observed_order\":{order},\"finest_reference\":{best_reference},\"original_reference\":{},\"ko_2017_reference\":{best_reference},\"relative_error_ko\":{},\"relative_error_original_reference\":{},\"tolerance\":{tolerance},\"verdict\":{limit_pass}}}",
            case.name(),
            if matches!(case, Case::Hole) {
                0.094
            } else {
                reference
            },
            (limit / best_reference - 1.0).abs(),
            (limit / reference - 1.0).abs()
        );
        assert!(
            limit_pass,
            "{} extrapolated limit {limit} outside cited reference bands",
            case.name()
        );
        if !distorted {
            let (planted, _) = sequence_limit(&planted_values)
                .expect("planted ladder must still contract so its limit is tested");
            let admitted = hemisphere_reference_band(case, planted, 0.02);
            println!(
                "{{\"suite\":\"irp-s1\",\"test\":\"{test}\",\"stage\":\"planted-bending-falsifier\",\"benchmark\":\"{}\",\"bending_scale\":1.05,\"u_h\":{planted_values:?},\"extrapolated_limit\":{planted},\"tolerance\":0.02,\"admitted\":{admitted}}}",
                case.name()
            );
            assert!(
                !admitted,
                "5% bending-stiffness error survived the extrapolated-limit gate"
            );
        }
    }
}

#[test]
fn g2_scordelis_lo_roof() {
    ladder(Case::Roof, false, 0.02, "g2_scordelis_lo_roof");
}
#[test]
fn g2_pinched_cylinder_with_diaphragms() {
    ladder(
        Case::Cylinder,
        false,
        0.02,
        "g2_pinched_cylinder_with_diaphragms",
    );
}
#[test]
fn g2_pinched_hemisphere_18deg_hole() {
    ladder(Case::Hole, false, 0.02, "g2_pinched_hemisphere_18deg_hole");
}
#[test]
fn g2_pinched_hemisphere_full() {
    ladder(Case::Full, false, 0.02, "g2_pinched_hemisphere_full");
}
#[test]
fn g2_obstacle_course_distorted_meshes() {
    for case in [Case::Roof, Case::Cylinder, Case::Hole] {
        ladder(case, true, 0.05, "g2_obstacle_course_distorted_meshes");
    }
}

#[test]
fn g0_drilling_alpha_insensitivity() {
    for case in [Case::Roof, Case::Cylinder, Case::Hole, Case::Full] {
        let n = *case.rungs().last().unwrap();
        let f = fixtures::build(case, n, false);
        let values: Vec<_> = [1e-4, 1e-3, 1e-2]
            .into_iter()
            .map(|a| response(&f, a))
            .collect();
        let relative = values
            .iter()
            .map(|u| (u / values[1] - 1.0).abs())
            .fold(0.0, f64::max);
        println!(
            "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1a\",\"test\":\"g0_drilling_alpha_insensitivity\",\"benchmark\":\"{}\",\"n\":{n},\"alpha\":[0.0001,0.001,0.01],\"displacements\":{values:?},\"relative_spread\":{relative},\"tolerance\":0.0001}}",
            case.name()
        );
        assert!(
            relative < 1e-4,
            "{} drilling variation {relative:e}",
            case.name()
        );
    }
}

fn dkt_response(f: &Fixture) -> Result<f64, fs_plate::PlateError> {
    use fs_plate::shell::{ShellSupport, assemble_shell};
    let section = fs_plate::PlateSection::isotropic(f.e, f.nu, f.h, 1.0)?;
    let original = assemble_shell(&f.mesh, &section, &[], ShellSupport::Free)?;
    let mut fixed = vec![false; 6 * f.mesh.nodes.len()];
    for c in &f.constraints {
        fixed[6 * c.node..6 * c.node + 6].copy_from_slice(&c.mask);
    }
    let mut count = 0;
    let map: Vec<_> = fixed
        .iter()
        .map(|fixed| {
            if *fixed {
                None
            } else {
                let i = count;
                count += 1;
                Some(i)
            }
        })
        .collect();
    let mut k = Coo::new(count, count);
    let mut force = vec![0.0; count];
    for (i, row) in map.iter().enumerate() {
        if let Some(row) = row {
            force[*row] = f.load[i];
            let (cols, vals) = original.k.row(i);
            for (&col, &value) in cols.iter().zip(vals) {
                if let Some(col) = map[col] {
                    k.push(*row, col, value);
                }
            }
        }
    }
    let k = k.assemble();
    let sym = SymbolicLdlt::analyze(&k, DirectOrdering::Amd)
        .map_err(fs_plate::PlateError::StaticFactor)?;
    let factor = sym
        .factor(&k, &LdltOptions::default())
        .map_err(fs_plate::PlateError::StaticFactor)?;
    if factor.inertia().negative != 0 {
        return Err(fs_plate::PlateError::StaticSolve {
            what: "DKT comparison is indefinite",
        });
    }
    let u = factor.solve(&force);
    Ok(f.sign * u[map[6 * f.observation + f.component].unwrap()])
}

#[test]
fn report_dkt_rung_on_obstacle_course() {
    for case in [Case::Roof, Case::Cylinder, Case::Hole, Case::Full] {
        for n in [4, 8, 16, 32] {
            let f = fixtures::build(case, n, false);
            match dkt_response(&f) {
                Ok(value) => println!(
                    "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1a\",\"test\":\"report_dkt_rung_on_obstacle_course\",\"model_form_rung\":\"CST-DKT\",\"benchmark\":\"{}\",\"n\":{n},\"observed\":{value},\"reference\":{},\"relative_error\":{},\"gated\":false}}",
                    case.name(),
                    case.reference(),
                    (value / case.reference() - 1.0).abs()
                ),
                Err(error) => println!(
                    "{{\"suite\":\"irp-s1\",\"test\":\"report_dkt_rung_on_obstacle_course\",\"benchmark\":\"{}\",\"n\":{n},\"gated\":false,\"verdict\":\"refused\",\"error\":\"{error}\"}}",
                    case.name()
                ),
            }
        }
    }
}
