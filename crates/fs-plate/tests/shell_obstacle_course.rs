//! IRP M1.2-W1a static obstacle course. No recording fit or DKT fallback.
//!
//! References: MacNeal & Harder, FEAD 1 (1985) 3–20, and Belytschko et al.,
//! CMAME 51 (1985) 221–258. Author reproduction: Ko, Lee, Lee & Bathe,
//! Computers & Structures 193 (2017) 187–206, doi:10.1016/j.compstruc.2017.08.003,
//! Figs.9/11/17/19 and Tables 8/10/15/17, available at
//! https://web.mit.edu/kjb/www/Principal_Publications/Performance_of_the_MITC3%2B_and_MITC4%2B_shell_elements_in_widely_used_benchmark_problems.pdf
//!
//! The bead gates use roof .3024 (not original .3086), cylinder 1.8248e-5,
//! cut-out hemisphere .094, and full hemisphere .0924. The 2017 cut-out
//! reference is .093; both are logged. These are pointwise benchmark bands,
//! not rigorous error enclosures or instrument/specimen validation.

#[path = "shell_obstacle_course/fixtures.rs"]
mod fixtures;
use fixtures::{Case, Fixture};
use fs_plate::shell::mitc::{MitcBudget, ShellMass, assemble_mitc3p, solve_static};
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

fn ladder(case: Case, distorted: bool, tolerance: f64, test: &str) {
    let started = std::time::Instant::now();
    let rungs = case.rungs();
    let reference = case.reference();
    let mut values = Vec::new();
    let mut errors = Vec::new();
    for &n in rungs {
        let fixture = fixtures::build(case, n, distorted);
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
    }
    let hs: Vec<_> = rungs.iter().map(|&n| 1.0 / n as f64).collect();
    let fit = fs_mms::fit_order(&fs_mms::RefinementLadder::new(hs, errors.clone()).unwrap());
    let monotone = errors[errors.len() - 3..].windows(2).all(|p| p[1] < p[0]);
    let passed = *errors.last().unwrap() <= tolerance && monotone;
    println!(
        "{{\"suite\":\"irp-s1\",\"bead\":\"M1.2-W1a\",\"test\":\"{test}\",\"tier\":\"G2\",\"inputs\":{{\"benchmark\":\"{}\",\"model\":\"{}\",\"element\":\"mitc3+\",\"mass\":\"consistent\",\"distorted\":{distorted},\"distortion_key\":{}}},\"provenance\":\"MacNeal-Harder 1985; Belytschko 1985; Ko et al. 2017\",\"rungs\":{{\"n_per_side\":{rungs:?}}},\"observed\":{{\"reference\":{reference},\"u_h\":{values:?},\"rel_err\":{errors:?},\"order\":{},\"monotone_last_three\":{monotone}}},\"alternate_reference\":{},\"tolerance\":{tolerance},\"verdict\":{},\"elapsed_s\":{}}}",
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
        *errors.last().unwrap() <= tolerance,
        "{} displacement {values:?}; reference {reference}; relative errors {errors:?}",
        case.name()
    );
    assert!(
        monotone,
        "{} last-three-rung error is not monotone: {errors:?}",
        case.name()
    );
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
