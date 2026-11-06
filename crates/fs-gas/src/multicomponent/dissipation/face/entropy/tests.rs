use super::*;
use crate::multicomponent::CaloricDomain;
use crate::multicomponent::tests::{bank, source};
use fs_qty::DynViscosity;
use fs_thermochem::equilibrium::composition::IdealGasEquilibriumV1;
use fs_thermochem::transport::kinetic::ThermalConductivityQuantityV1;
use fs_thermochem::{MassSpecificEnergyQuantityV1, SpeciesId};

// Each row is [Cp/R constant, Cp/R linear slope, M/R]. The entropy and
// divided-enthalpy oracles below integrate these polynomials independently.
const VARIABLE: [[f64; 3]; 2] = [[3.5, 0.2, 1.0], [4.0, 0.4, 2.0]];
const T0: f64 = 0.001;

fn gas(rows: [[f64; 3]; 2], formation: f64) -> MixtureGas {
    MixtureGas::new(
        bank(vec![
            source("A", rows[0][2] * R, rows[0][0], rows[0][1], formation),
            source(
                "B",
                rows[1][2] * R,
                rows[1][0],
                rows[1][1],
                -2.0 * formation,
            ),
        ]),
        CaloricDomain::new(Temperature::new(T0), Temperature::new(50.0)),
        || false,
    )
    .unwrap()
}
fn domain() -> TransportDomain {
    TransportDomain {
        minimum_temperature_k: T0,
        maximum_temperature_k: 50.0,
        minimum_pressure_pa: 1e-8,
        maximum_pressure_pa: 1e8,
    }
}
fn molecular(gas: &MixtureGas, mu: f64, bulk: f64, k: f64, d: f64) -> ConstantTransport {
    ConstantTransport::new(
        gas.source().clone(),
        "manufactured-conjugate-face",
        1,
        domain(),
        TransportCoefficients {
            dynamic_viscosity: DynViscosity::new(mu),
            bulk_viscosity: DynViscosity::new(bulk),
            thermal_conductivity: ThermalConductivityQuantityV1::new(k),
        },
        vec![0.0, d, d, 0.0],
    )
    .unwrap()
}
fn trace<'a>(gas: &MixtureGas, rho: &'a [f64; 2], u: f64, t: f64) -> Trace<'a> {
    let p = R
        * t
        * rho
            .iter()
            .zip(gas.source().models())
            .map(|(r, m)| r / m.molar_mass().value())
            .sum::<f64>();
    let q = gas.conserved(rho, u, p).unwrap();
    let primitive = gas.decode_cell(q.as_cell(), Some(t), &mut || false).unwrap();
    Trace {
        densities: rho,
        primitive,
        transverse_velocity: 0.0,
        cv_density: heat_capacity_density(gas, rho, primitive).unwrap(),
    }
}
#[allow(clippy::too_many_arguments)]
fn flux<L: TransportLaw>(
    gas: &MixtureGas,
    law: &mut L,
    left: &[f64; 2],
    right: &[f64; 2],
    u: [f64; 2],
    t: [f64; 2],
    widths: [f64; 2],
) -> ([f64; 4], TransportFaceReport) {
    let (diffusion, mass_ratio) = prepare_diffusion(gas).unwrap();
    let mut workspace = diffusion.workspace().unwrap();
    let mut binary = [0.0; 4];
    let mut out = [0.0; 4];
    let report = interior(
        gas,
        law,
        &diffusion,
        &mut workspace,
        trace(gas, left, u[0], t[0]),
        trace(gas, right, u[1], t[1]),
        widths[0],
        widths[1],
        mass_ratio,
        &mut binary,
        &mut out,
        MaxwellStefanSolveV1::default(),
        &mut || false,
    )
    .unwrap();
    (out, report)
}
fn near(actual: f64, expected: f64) {
    assert!(
        (actual - expected).abs() < 3e-12 * expected.abs().max(1.0),
        "actual={actual:e}, expected={expected:e}"
    );
}
fn beta_enthalpy(row: [f64; 3], left_t: f64, right_t: f64) -> f64 {
    let [a, b, mass] = row;
    let t_star = if left_t == right_t {
        left_t
    } else {
        left_t * right_t * ((right_t - left_t) / left_t).ln_1p() / (right_t - left_t)
    };
    (a * t_star + 0.5 * b * left_t * right_t - (a - 1.0) * T0 - 0.5 * b * T0 * T0) / mass
}
fn entropy_variables(rows: [[f64; 3]; 2], rho: [f64; 2], u: f64, t: f64) -> [f64; 4] {
    let mut w = [0.0; 4];
    for i in 0..2 {
        let [a, b, mass] = rows[i];
        let e = ((a - 1.0) * (t - T0) + 0.5 * b * (t * t - T0 * T0)) / mass;
        let s = ((a - 1.0) * t.ln() + b * t - (rho[i] / mass).ln()) / mass;
        w[i] = s - e / t - 1.0 / mass + 0.5 * u * u / t;
    }
    w[2] = -u / t;
    w[3] = 1.0 / t;
    w
}
fn pair_rate(
    rows: [[f64; 3]; 2],
    rho: [[f64; 2]; 2],
    u: [f64; 2],
    t: [f64; 2],
    f: [f64; 4],
) -> f64 {
    let a = entropy_variables(rows, rho[0], u[0], t[0]);
    let b = entropy_variables(rows, rho[1], u[1], t[1]);
    (0..4).map(|k| (b[k] - a[k]) * f[k]).sum()
}

#[test]
fn source_divided_enthalpy_repairs_the_energy_conserving_negative_rate_face() {
    let rows = [[3.5, 0.0, 1.0], [6.5, 0.0, 1.0]];
    let gas = gas(rows, 0.0);
    let rho = [[0.1, 0.9], [0.9, 0.1]];
    let (f, report) = flux(
        &gas,
        &mut molecular(&gas, 0.0, 0.0, 0.0, 0.01),
        &rho[0],
        &rho[1],
        [0.0; 2],
        [1.0, 10.0],
        [0.5; 2],
    );
    let j = -0.01 * 9.0_f64.ln();
    near(f[0], j);
    near(f[1], -j);
    near(
        f[3],
        (beta_enthalpy(rows[0], 1.0, 10.0) - beta_enthalpy(rows[1], 1.0, 10.0)) * j,
    );
    // Independent derivatives of the two actual ideal-mixture states. The
    // original arithmetic-enthalpy/Fick face has negative production.
    let legacy = [-0.016, 0.016, 0.0, (3.0 * 5.5 - 3.0 * T0) * 0.016];
    near(
        0.7 * pair_rate(rows, rho, [0.0; 2], [1.0, 10.0], legacy),
        -0.039735310343068805,
    );
    let production = 0.7 * pair_rate(rows, rho, [0.0; 2], [1.0, 10.0], f);
    near(production, 0.7 * 2.0 * 0.01 * 9.0_f64.ln().powi(2));
    assert!(production > 0.0);
    assert_eq!(report.mass_frame_residual_kg_m2_s, 0.0);
    let (reverse, _) = flux(
        &gas,
        &mut molecular(&gas, 0.0, 0.0, 0.0, 0.01),
        &rho[1],
        &rho[0],
        [0.0; 2],
        [10.0, 1.0],
        [0.5; 2],
    );
    for k in 0..4 {
        near(reverse[k], -f[k]);
    }
}

#[test]
fn variable_original_calorics_keep_the_divided_mean_at_equal_and_adjacent_temperatures() {
    let rho = [[0.8, 0.2], [0.3, 0.7]];
    for t in [[2.0, 7.0], [7.0, 2.0], [2.0, 2.0], [2.0, 2.0_f64.next_up()]] {
        let mut unshifted = None;
        for datum in [0.0, 1e18] {
            let g = gas(VARIABLE, datum);
            let (f, _) = flux(
                &g,
                &mut molecular(&g, 0.0, 0.0, 0.0, 0.03),
                &rho[0],
                &rho[1],
                [0.0; 2],
                t,
                [0.2, 0.8],
            );
            let expected = (beta_enthalpy(VARIABLE[0], t[0], t[1])
                - beta_enthalpy(VARIABLE[1], t[0], t[1]))
                * f[0];
            near(f[3], expected);
            assert!(pair_rate(VARIABLE, rho, [0.0; 2], t, f) > 0.0);
            if let Some(previous) = unshifted {
                assert_eq!(f, previous);
            }
            unshifted = Some(f);
        }
    }
    // Constant partial pressures and no authored Soret heat must not produce
    // diffusion merely because the species have different, variable Cp.
    let g = gas(VARIABLE, 1e18);
    let (f, _) = flux(
        &g,
        &mut molecular(&g, 0.0, 0.0, 0.0, 0.03),
        &[0.4, 0.6],
        &[0.08, 0.12],
        [0.0; 2],
        [2.0, 10.0],
        [0.2, 0.8],
    );
    assert!(f[0].abs() < 1e-15 && f[1].abs() < 1e-15);
    assert!(f[3].abs() < 1e-14);
}

#[test]
fn original_soret_dufour_and_viscous_work_close_the_actual_pair_entropy() {
    let g = gas(VARIABLE, 1e18);
    let rho = [[0.8, 0.2], [0.3, 0.7]];
    let t = [2.0, 7.0];
    let u = [-0.3, 1.2];
    let widths = [0.2, 0.8];
    let d = 0.5;
    let x_left = rho[0][0] / (rho[0][0] + 0.5 * rho[0][1]);
    let x_right = rho[1][0] / (rho[1][0] + 0.5 * rho[1][1]);
    let x = 0.8 * x_left + 0.2 * x_right;
    let tf = 0.8 * t[0] + 0.2 * t[1];
    let pf =
        0.8 * (rho[0][0] + 0.5 * rho[0][1]) * t[0] + 0.2 * (rho[1][0] + 0.5 * rho[1][1]) * t[1];
    let concentration = pf / (R * tf);
    let mut reference = None;
    for offset in [0.0, 1e12] {
        let mut law = SoretDufourTransport::new(
            molecular(&g, 0.07, 0.02, 0.09, 0.03),
            "declared-cross-heat",
            1,
            domain(),
            vec![
                (
                    SpeciesId::new("A").unwrap(),
                    MassSpecificEnergyQuantityV1::new(125.0 + offset),
                ),
                (
                    SpeciesId::new("B").unwrap(),
                    MassSpecificEnergyQuantityV1::new(-75.0 + offset),
                ),
            ],
        )
        .unwrap();
        let (f, report) = flux(&g, &mut law, &rho[0], &rho[1], u, t, widths);
        let relative_velocity =
            f[0] / (R * concentration * x) - f[1] / (2.0 * R * concentration * (1.0 - x));
        let friction = R * d * concentration * x * (1.0 - x) / 0.03 * relative_velocity.powi(2);
        let fourier = 0.09 * (t[1] - t[0]).powi(2) / (t[0] * t[1] * d);
        let viscous =
            (4.0 / 3.0 * 0.07 + 0.02) * (0.5 / t[0] + 0.5 / t[1]) * (u[1] - u[0]).powi(2) / d;
        let actual = pair_rate(VARIABLE, rho, u, t, f);
        near(actual, friction + fourier + viscous);
        assert!(actual > 0.0);
        near(
            report.reduced_heat_flux_w_m2,
            -0.09 * (t[1] - t[0]) / d + 200.0 * f[0],
        );
        near(
            f[3],
            (beta_enthalpy(VARIABLE[0], t[0], t[1]) - beta_enthalpy(VARIABLE[1], t[0], t[1]))
                * f[0]
                + report.reduced_heat_flux_w_m2
                + 0.5 * (u[0] + u[1]) * f[2],
        );
        if let Some(previous) = reference {
            assert_eq!(f, previous);
        }
        reference = Some(f);
    }
}

#[test]
fn divided_energy_retains_original_piecewise_caloric_intervals() {
    use fs_evidence::ValidityDomain;
    use fs_thermochem::{
        MolarEnergyQuantityV1, MolarThermalQuantityV1, Nasa9RegionV1, Nasa9StandardStateModelV1,
    };
    let lower = source("A", R, 3.5, 0.2, 0.0);
    let upper = source("A", R, 5.5, 0.1, 0.0);
    let mut first = lower.regions()[0].card().clone();
    first.validity = ValidityDomain::unconstrained().with("T", T0, 3.0);
    let mut second = upper.regions()[0].card().clone();
    second.validity = ValidityDomain::unconstrained().with("T", 3.0, 50.0);
    let original = Nasa9StandardStateModelV1::new(
        lower.species().clone(),
        lower.molar_mass(),
        lower.convention().clone(),
        vec![
            Nasa9RegionV1::from_card(first).unwrap(),
            Nasa9RegionV1::from_card(second).unwrap(),
        ],
    )
    .unwrap()
    .with_caloric_continuity(
        Temperature::new(2.0),
        MolarEnergyQuantityV1::new(1e6),
        MolarThermalQuantityV1::new(1e6),
    )
    .unwrap();
    let g = MixtureGas::new(
        bank(vec![original, source("B", 2.0 * R, 4.0, 0.4, 0.0)]),
        CaloricDomain::new(Temperature::new(T0), Temperature::new(50.0)),
        || false,
    )
    .unwrap();
    let (f, _) = flux(
        &g,
        &mut molecular(&g, 0.0, 0.0, 0.0, 0.03),
        &[0.8, 0.2],
        &[0.3, 0.7],
        [0.0; 2],
        [2.0, 7.0],
        [0.2, 0.8],
    );
    // Independently integrate h(T) d(1/T), retaining the actual Cp jump and
    // the explicit continuity offset at T=3. This is not h at a mean T.
    let c1 = 2.5 * T0 + 0.1 * T0 * T0;
    let c2 = c1 + (5.5 - 3.5) * 3.0 + 0.5 * (0.1 - 0.2) * 9.0;
    let interval = |a: f64, b: f64, c: f64, lo: f64, hi: f64| {
        -a * (hi / lo).ln() - 0.5 * b * (hi - lo) - c * (1.0 / hi - 1.0 / lo)
    };
    let h_a = (interval(3.5, 0.2, c1, 2.0, 3.0) + interval(5.5, 0.1, c2, 3.0, 7.0))
        / (1.0 / 7.0 - 1.0 / 2.0);
    near(f[3], (h_a - beta_enthalpy(VARIABLE[1], 2.0, 7.0)) * f[0]);
}

struct Counted {
    inner: ConstantTransport,
    visits: usize,
}
impl TransportLaw for Counted {
    fn source(&self) -> &IdealGasEquilibriumV1 {
        self.inner.source()
    }
    fn evaluate(
        &mut self,
        t: Temperature,
        p: Pressure,
        x: &[f64],
        binary: &mut [f64],
        cancelled: &mut dyn FnMut() -> bool,
    ) -> Result<TransportCoefficients, DissipationError> {
        self.visits += 1;
        self.inner.evaluate(t, p, x, binary, cancelled)
    }
}

#[test]
fn common_absence_is_admitted_and_one_sided_zero_or_cancellation_precedes_source_queries() {
    let g = gas(VARIABLE, 0.0);
    let mut law = Counted {
        inner: molecular(&g, 0.0, 0.0, 0.09, 0.03),
        visits: 0,
    };
    let (f, _) = flux(
        &g,
        &mut law,
        &[1.0, 0.0],
        &[0.5, 0.0],
        [0.0; 2],
        [2.0, 4.0],
        [0.5; 2],
    );
    assert_eq!(law.visits, 1);
    assert_eq!(f[0], 0.0);
    assert_eq!(f[1], 0.0);
    near(f[3], -0.09 * 2.0 / 0.5);
    let (diffusion, mass_ratio) = prepare_diffusion(&g).unwrap();
    let mut workspace = diffusion.workspace().unwrap();
    let mut binary = [0.0; 4];
    let left = [1.0, 0.0];
    let right = [0.5, 0.5];
    for cancel in [false, true] {
        let mut out = [137.0; 4];
        let error = interior(
            &g,
            &mut law,
            &diffusion,
            &mut workspace,
            trace(&g, &left, 0.0, 2.0),
            trace(&g, &right, 0.0, 2.0),
            0.5,
            0.5,
            mass_ratio,
            &mut binary,
            &mut out,
            MaxwellStefanSolveV1::default(),
            &mut || cancel,
        )
        .unwrap_err();
        if cancel {
            assert_eq!(error, DissipationError::Cancelled);
        } else {
            assert!(matches!(
                error,
                DissipationError::Invalid(
                    "entropy-conjugate diffusion requires common positive species support"
                )
            ));
        }
        assert_eq!(out, [137.0; 4]);
        assert_eq!(law.visits, 1);
    }
}
