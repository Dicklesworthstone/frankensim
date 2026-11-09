//! Barycentric projection of finite chemical-potential differences.
use fs_math::dd::Dd;
use fs_qty::AmountConcentration;
use fs_thermochem::UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R;
use fs_thermochem::transport::maxwell_stefan::{
    MaxwellStefanErrorV1 as Error, MaxwellStefanEvaluationV1, MaxwellStefanSolveV1,
    MaxwellStefanV1, MaxwellStefanWorkspaceV1,
};

pub(super) struct Evaluation<'a> {
    pub forces: Vec<f64>,
    pub reduced_gradients: Vec<f64>,
    pub diffusion: MaxwellStefanEvaluationV1<'a>,
    pub drag_entropy: f64,
    pub force_entropy: f64,
}

/// Original force/friction algebra, shared by the source-bound face and native
/// physical tests. Concentration and gradients have already been supplied by
/// the actual EOS; this function cannot manufacture ideal-gas state data.
#[allow(clippy::too_many_arguments)]
pub(super) fn evaluate<'a>(
    model: &MaxwellStefanV1,
    x: &[f64],
    gradients: &[f64],
    concentration: f64,
    binary: &[f64],
    controls: MaxwellStefanSolveV1,
    entropy_tolerance: f64,
    work: &'a mut MaxwellStefanWorkspaceV1,
    mut cancelled: impl FnMut() -> bool,
) -> Result<Evaluation<'a>, Error> {
    if cancelled() {
        return Err(Error::Cancelled);
    }
    if !entropy_tolerance.is_finite() || entropy_tolerance <= 0.0 {
        return Err(Error::Invalid(
            "positive real-gas entropy pairing tolerance required",
        ));
    }
    let masses: Vec<_> = model.molar_masses().iter().map(|m| m.value()).collect();
    let mut forces = vec![0.0; x.len()];
    let reduced = project_with_gauge(x, &masses, gradients, &mut forces).map_err(Error::Invalid)?;
    let diffusion = model.evaluate_driving_forces(
        x,
        &forces,
        AmountConcentration::new(concentration),
        binary,
        controls,
        work,
        &mut cancelled,
    )?;
    let drag = diffusion
        .drag_entropy_production_w_m3_k()
        .ok_or(Error::Invalid(
            "positive common support needs resolved real-gas friction entropy",
        ))?;
    let mut power = Dd::ZERO;
    for (&molar, &gradient) in diffusion.molar_fluxes().iter().zip(&reduced) {
        if cancelled() {
            return Err(Error::Cancelled);
        }
        // Cancel the physically irrelevant mass mode BEFORE multiplication by
        // a rounded flux. Its small frame residual is a separate admission.
        power = power - Dd::from_f64(molar) * gradient * Dd::from_f64(R);
    }
    let force_power = power.to_f64();
    let allowance =
        entropy_tolerance + controls.relative_tolerance * drag.abs().max(force_power.abs());
    if ![power.hi, power.lo, force_power, allowance]
        .iter()
        .all(|v| v.is_finite())
        || force_power < -allowance
        || (force_power - drag).abs() > allowance
    {
        return Err(Error::Arithmetic(
            "real-gas chemical work did not close original pair friction",
        ));
    }
    Ok(Evaluation {
        forces,
        reduced_gradients: reduced[..x.len()].iter().map(|g| g.to_f64()).collect(),
        diffusion,
        drag_entropy: drag,
        force_entropy: force_power,
    })
}

/// d_i=x_i g_i-Y_i sum(x_j g_j). Remove one mass-proportional gauge before
/// reduction; subtracting a common scalar force would select the molar frame.
/// The same original positive x supplies both drag and entropy conjugacy.
pub(super) fn project_with_gauge(
    x: &[f64],
    mass: &[f64],
    g: &[f64],
    d: &mut [f64],
) -> Result<[Dd; 32], &'static str> {
    let n = x.len();
    if n == 0
        || n > 32
        || mass.len() != n
        || g.len() != n
        || d.len() != n
        || x.iter().chain(mass).any(|v| !v.is_finite() || *v <= 0.0)
        || g.iter().any(|v| !v.is_finite())
    {
        return Err("positive bounded real-gas force arrays required");
    }
    let pivot = (0..n).max_by(|&i, &j| mass[i].total_cmp(&mass[j])).unwrap();
    let mut total = Dd::ZERO;
    let mut mean_mass = Dd::ZERO;
    let mut sum_x = Dd::ZERO;
    let mut reduced = [Dd::ZERO; 32];
    let mut scaled_mass = [Dd::ZERO; 32];
    for i in 0..n {
        // Ratios and mean stay double-double. Rounding 1/3 before multiplying
        // a large mechanical gauge would create a fictitious diffusion force.
        scaled_mass[i] = Dd::from_f64(mass[i]) / Dd::from_f64(mass[pivot]);
        reduced[i] = Dd::from_f64(g[i]) - scaled_mass[i] * Dd::from_f64(g[pivot]);
        total = total + Dd::from_f64(x[i]) * reduced[i];
        mean_mass = mean_mass + Dd::from_f64(x[i]) * scaled_mass[i];
        sum_x = sum_x + Dd::from_f64(x[i]);
    }
    if !mean_mass.hi.is_finite()
        || !mean_mass.lo.is_finite()
        || mean_mass.to_f64() <= 0.0
        || (sum_x.to_f64() - 1.0).abs() > 256.0 * f64::EPSILON
    {
        return Err("real-gas force mean mass or composition is invalid");
    }
    for i in 0..n {
        let specific = reduced[i] - (scaled_mass[i] / mean_mass) * total;
        d[i] = (Dd::from_f64(x[i]) * specific).to_f64();
        if !d[i].is_finite() || (specific.to_f64() != 0.0 && d[i] == 0.0) {
            return Err("real-gas projected chemical force is unrepresentable");
        }
    }
    Ok(reduced)
}

#[cfg(test)]
fn project(x: &[f64], mass: &[f64], g: &[f64], d: &mut [f64]) -> Result<(), &'static str> {
    project_with_gauge(x, mass, g, d).map(|_| ())
}

#[cfg(test)]
mod tests {
    use super::*;
    #[test]
    fn mass_gauge_and_binary_force_follow_the_constrained_potential() {
        let (x, m, g) = ([0.25, 0.75], [0.02, 0.04], [4.0, -2.0]);
        let mut d = [0.0; 2];
        project(&x, &m, &g, &mut d).unwrap();
        let expected = x[0] * x[1] * (m[1] * g[0] - m[0] * g[1]) / (x[0] * m[0] + x[1] * m[1]);
        assert!((d[0] - expected).abs() < 4e-16);
        assert!((d[0] + d[1]).abs() < 4e-16);
        let shifted = [g[0] + 128.0 * m[0], g[1] + 128.0 * m[1]];
        let mut next = [0.0; 2];
        project(&x, &m, &shifted, &mut next).unwrap();
        for i in 0..2 {
            assert!((next[i] - d[i]).abs() < 8e-16);
        }
        project(&x, &m, &[16.0, 32.0], &mut next).unwrap();
        assert_eq!(next, [0.0; 2]);
    }
    #[test]
    fn non_dyadic_mass_ratios_remove_a_large_exact_mechanical_gauge() {
        use fs_qty::MolarMass;
        use fs_thermochem::SpeciesId;
        let x = [0.2, 0.3, 0.5];
        let masses = [1.0, 2.0, 3.0];
        let gauge = 2.0_f64.powi(40);
        let base = [-3.0, 7.0, -1.0];
        let pure = masses.map(|m| gauge * m);
        let shifted = std::array::from_fn::<_, 3, _>(|i| pure[i] + base[i]);
        let mut expected = [0.0; 3];
        let mut actual = [0.0; 3];
        project(&x, &masses, &base, &mut expected).unwrap();
        project(&x, &masses, &shifted, &mut actual).unwrap();
        for i in 0..3 {
            assert!((actual[i] - expected[i]).abs() < 1e-14);
        }
        project(&x, &masses, &pure, &mut actual).unwrap();
        assert!(actual.iter().all(|d| d.abs() < 1e-18));
        let model = MaxwellStefanV1::new(
            masses
                .iter()
                .enumerate()
                .map(|(i, &mass)| {
                    (
                        SpeciesId::new(format!("G{i}")).unwrap(),
                        MolarMass::new(mass),
                    )
                })
                .collect(),
        )
        .unwrap();
        let binary = [0.0, 1e-5, 2e-5, 1e-5, 0.0, 3e-5, 2e-5, 3e-5, 0.0];
        let mut work = model.workspace().unwrap();
        let first = evaluate(
            &model,
            &x,
            &base,
            40.0,
            &binary,
            MaxwellStefanSolveV1::default(),
            1e-12,
            &mut work,
            || false,
        )
        .unwrap();
        let flux = first.diffusion.molar_fluxes().to_vec();
        let power = first.force_entropy;
        let second = evaluate(
            &model,
            &x,
            &shifted,
            40.0,
            &binary,
            MaxwellStefanSolveV1::default(),
            1e-12,
            &mut work,
            || false,
        )
        .unwrap();
        for i in 0..3 {
            assert!((flux[i] - second.diffusion.molar_fluxes()[i]).abs() < 1e-15);
        }
        assert!((power - second.force_entropy).abs() < 1e-14);
    }
    #[test]
    fn ternary_friction_pairing_and_trace_force_remain_resolved() {
        let x = [0.2, 0.3, 0.5];
        let m = [0.004, 0.028, 0.044];
        let g = [-3.0, 7.0, -1.0];
        let mut d = [0.0; 3];
        project(&x, &m, &g, &mut d).unwrap();
        assert!(d.iter().sum::<f64>().abs() < 1e-15);
        // Equal D has N_i=-cD*d_i+x_i*a; mass closure determines a.
        let c_d = 0.03;
        let mean = x.iter().zip(m).map(|(x, m)| x * m).sum::<f64>();
        let a = c_d * d.iter().zip(m).map(|(d, m)| d * m).sum::<f64>() / mean;
        let n: Vec<_> = (0..3).map(|i| -c_d * d[i] + x[i] * a).collect();
        assert!(n.iter().zip(m).map(|(n, m)| n * m).sum::<f64>().abs() < 1e-18);
        let force = -n.iter().zip(g).map(|(n, g)| n * g).sum::<f64>();
        let mut friction = 0.0;
        for i in 0..3 {
            for j in i + 1..3 {
                let slip = n[i] / x[i] - n[j] / x[j];
                friction += x[i] * x[j] * slip * slip / c_d;
            }
        }
        assert!(force > 0.0);
        assert!((force - friction).abs() < 1e-14);
        let mut trace = [0.0; 2];
        project(&[1e-250, 1.0], &[0.004, 0.028], &[1.0, 0.0], &mut trace).unwrap();
        assert!(trace[0] > 0.0 && trace[1] < 0.0);
        assert!((trace[0] / 1e-250 - 1.0).abs() < 1e-14);
    }
    #[test]
    fn unsupported_force_inputs_refuse() {
        let mut d = [0.0; 2];
        assert!(project(&[0.0, 1.0], &[1.0, 2.0], &[0.0, 1.0], &mut d).is_err());
        assert!(project(&[0.5, 0.5], &[1.0, 2.0], &[f64::INFINITY, 1.0], &mut d).is_err());
        assert!(project(&[0.5, 0.6], &[1.0, 2.0], &[0.0, 1.0], &mut d).is_err());
        assert!(project(&[0.5, 0.5], &[f64::MAX, f64::MAX], &[1.0, 0.0], &mut d).is_ok());
    }

    fn binary() -> MaxwellStefanV1 {
        use fs_qty::MolarMass;
        use fs_thermochem::SpeciesId;
        MaxwellStefanV1::new(vec![
            (SpeciesId::new("A").unwrap(), MolarMass::new(0.028)),
            (SpeciesId::new("B").unwrap(), MolarMass::new(0.004)),
        ])
        .unwrap()
    }

    #[test]
    fn original_ms_binary_matches_the_mass_frame_closed_form() {
        let model = binary();
        let mut work = model.workspace().unwrap();
        let x = [0.55, 0.45];
        let g = [-37.0, 48.0];
        let c = 1500.;
        let d = 3e-8;
        let evaluated = evaluate(
            &model,
            &x,
            &g,
            c,
            &[0.0, d, d, 0.0],
            MaxwellStefanSolveV1::default(),
            1e-12,
            &mut work,
            || false,
        )
        .unwrap();
        let m = [0.028, 0.004];
        let mean = x[0] * m[0] + x[1] * m[1];
        let drive = x[0] * x[1] * (m[1] * g[0] - m[0] * g[1]) / mean;
        let expected = -c * d * drive * m[1] / mean;
        assert!((evaluated.diffusion.molar_fluxes()[0] / expected - 1.0).abs() < 2e-14);
        assert!(
            (evaluated.diffusion.molar_fluxes()[1] / (-m[0] / m[1] * expected) - 1.0).abs() < 2e-14
        );
        assert!((evaluated.force_entropy - evaluated.drag_entropy).abs() < 1e-12);
        assert!(
            evaluate(
                &model,
                &x,
                &g,
                c,
                &[0.0, d, d, 0.0],
                MaxwellStefanSolveV1::default(),
                1e-12,
                &mut work,
                || true
            )
            .is_err()
        );
    }

    #[test]
    fn original_pr_and_ms_match_an_independent_dense_binary_oracle() {
        use fs_qty::{MolarMass, Pressure, Temperature};
        use fs_thermochem::eos::peng_robinson::{
            PengRobinsonComponentV1, PengRobinsonMixtureV1, PengRobinsonPhaseV1,
        };
        let component = |mass, tc, pc, w| PengRobinsonComponentV1 {
            molar_mass: MolarMass::new(mass),
            critical_temperature: Temperature::new(tc),
            critical_pressure: Pressure::new(pc),
            acentric_factor: w,
        };
        let components = vec![
            component(0.028, 190.564, 4.5992e6, 0.01142),
            component(0.004, 150.687, 4.863e6, -0.00219),
        ];
        let left = PengRobinsonMixtureV1::new(
            components.clone(),
            vec![0.75, 0.25],
            vec![0., 0.03, 0.03, 0.],
        )
        .unwrap();
        let right =
            PengRobinsonMixtureV1::new(components, vec![0.35, 0.65], vec![0., 0.03, 0.03, 0.])
                .unwrap();
        let l = left
            .evaluate_pt(Temperature::new(330.), Pressure::new(4e6))
            .unwrap();
        let r = right
            .evaluate_pt(Temperature::new(330.), Pressure::new(4e6))
            .unwrap();
        let (l, r) = (
            l.select(PengRobinsonPhaseV1::Vapor),
            r.select(PengRobinsonPhaseV1::Vapor),
        );
        let x = [0.55, 0.45];
        let mut g = [0.; 2];
        for i in 0..2 {
            g[i] = (fs_math::det::ln(right.mole_fractions()[i] / left.mole_fractions()[i])
                + r.log_fugacity_coefficients()[i]
                - l.log_fugacity_coefficients()[i])
                / 0.02;
        }
        let c = 0.5 / l.molar_volume_m3_per_mol() + 0.5 / r.molar_volume_m3_per_mol();
        let model = binary();
        let mut work = model.workspace().unwrap();
        let value = evaluate(
            &model,
            &x,
            &g,
            c,
            &[0., 3e-8, 3e-8, 0.],
            MaxwellStefanSolveV1::default(),
            1e-12,
            &mut work,
            || false,
        )
        .unwrap();
        // Independent scalar PR cubic bisection and closed-form binary MS,
        // using the declared synthetic coefficients rather than either solver.
        assert!((c / 1514.9804470555855 - 1.0).abs() < 2e-13);
        assert!((value.diffusion.molar_fluxes()[0] / 0.0002253251762677151 - 1.0).abs() < 2e-12);
        assert!((value.force_entropy / 0.6938842813864166 - 1.0).abs() < 2e-12);
        assert!((c / (4e6 / (R * 330.)) - 1.0).abs() > 0.03);
    }
}
