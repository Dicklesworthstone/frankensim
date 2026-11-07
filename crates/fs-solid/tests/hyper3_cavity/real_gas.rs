//! Dense homogeneous gas in the native deformable cavity, with independent
//! scalar Peng--Robinson/Helmholtz references. The source is deliberately
//! synthetic: its molar mass and calorics are not a measured methane card.

use super::*;
use fs_material::thermochemical_gas::energy::GasEnergySolve;
use fs_material::thermochemical_gas::energy::entropy::GasEntropySolve;
use fs_material::thermochemical_gas::nonideal::{
    PengRobinsonSpeciesBinding, ThermochemicalPengRobinsonModel,
    energy::parcel::ClosedRealGasParcel,
};
use fs_qty::{Pressure, Temperature};
use fs_solid::hyper3::pressure::cavity::AdiabaticCavity;
use fs_thermochem::eos::peng_robinson::PengRobinsonPhaseV1;
use fs_thermochem::{SpeciesId, UNIVERSAL_GAS_CONSTANT_J_PER_MOL_K as R};

const MOLAR_MASS: f64 = 0.028;
const TC: f64 = 190.564;
const PC: f64 = 4.5992e6;
const OMEGA: f64 = 0.01142;
const CP_SLOPE: f64 = 0.001;
const INITIAL_PRESSURE: f64 = 5.0e6;

fn controls() -> GasEntropySolve {
    GasEntropySolve {
        minimum_temperature_k: 200.0,
        maximum_temperature_k: 600.0,
        entropy_tolerance_j_kg_k: 1e-10,
        maximum_iterations: 128,
    }
}

fn parcel(volume: f64, cp_slope: f64, datum: f64) -> ClosedRealGasParcel {
    let original = super::gas_fixture::gas(volume, INITIAL_PRESSURE, cp_slope, datum);
    let model = ThermochemicalPengRobinsonModel::new(
        original.model().clone(),
        vec![PengRobinsonSpeciesBinding {
            species: SpeciesId::new("synthetic-cavity-gas").unwrap(),
            critical_temperature: Temperature::new(TC),
            critical_pressure: Pressure::new(PC),
            acentric_factor: OMEGA,
        }],
        vec![],
    )
    .unwrap();
    ClosedRealGasParcel::new(
        model,
        original.point().clone(),
        volume,
        PengRobinsonPhaseV1::Vapor,
    )
    .unwrap()
}

#[derive(Clone, Copy)]
struct PrPoint {
    temperature: f64,
    pressure: f64,
    residual_energy_j_mol: f64,
    attractive_entropy_j_mol_k: f64,
    bulk_modulus: f64,
}

fn covolume() -> f64 {
    0.07780 * R * TC / PC
}

/// Independent pure-fluid formulas in molar volume v. The residual Helmholtz
/// energy is -RT ln(1-b/v) - a(T) I(v), where I is the logarithmic attraction
/// integral. Its T and v derivatives give pressure, Cv and the isentropic bulk
/// modulus; no production EOS, entropy inversion or acoustic kernel is called.
fn pr_point(temperature: f64, volume: f64, cp_slope: f64) -> PrPoint {
    let kappa = 0.37464 + 1.54226 * OMEGA - 0.26992 * OMEGA * OMEGA;
    let a0 = 0.45724 * R * R * TC * TC / PC;
    let alpha_root = 1.0 + kappa * (1.0 - (temperature / TC).sqrt());
    let a = a0 * alpha_root * alpha_root;
    let da = -a0 * kappa * alpha_root / (TC * temperature).sqrt();
    let dda = a0 * kappa * (1.0 + kappa) / (2.0 * TC.sqrt() * temperature.powf(1.5));
    let b = covolume();
    let root_two = 2.0_f64.sqrt();
    let integral = ((volume + (1.0 + root_two) * b) / (volume + (1.0 - root_two) * b)).ln()
        / (2.0 * root_two * b);
    let denominator = volume * volume + 2.0 * b * volume - b * b;
    let pressure = R * temperature / (volume - b) - a / denominator;
    let p_t = R / (volume - b) - da / denominator;
    let p_v =
        -R * temperature / (volume - b).powi(2) + 2.0 * a * (volume + b) / denominator.powi(2);
    let cv = R * (2.5 + cp_slope * temperature) + temperature * dda * integral;
    PrPoint {
        temperature,
        pressure,
        residual_energy_j_mol: (temperature * da - a) * integral,
        attractive_entropy_j_mol_k: da * integral,
        bulk_modulus: -volume * (p_v - temperature * p_t * p_t / cv),
    }
}

fn reference_molar_volume() -> f64 {
    // T=300 K is above this synthetic critical temperature: the admitted
    // homogeneous branch is unique and pressure decreases with molar volume.
    let (mut lo, mut hi) = (1.001 * covolume(), 20.0 * R * 300.0 / INITIAL_PRESSURE);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if pr_point(300.0, mid, CP_SLOPE).pressure > INITIAL_PRESSURE {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    0.5 * (lo + hi)
}

fn entropy_change(t0: f64, v0: f64, t1: f64, v1: f64, cp_slope: f64) -> f64 {
    // Ideal-volume entropy and excluded-volume entropy combine into
    // R ln((v1-b)/(v0-b)); common entropy/energy datums disappear analytically.
    R * (2.5 * (t1 / t0).ln() + cp_slope * (t1 - t0) + ((v1 - covolume()) / (v0 - covolume())).ln())
        + pr_point(t1, v1, cp_slope).attractive_entropy_j_mol_k
        - pr_point(t0, v0, cp_slope).attractive_entropy_j_mol_k
}

fn isentrope(volume_ratio: f64) -> (PrPoint, f64) {
    let v0 = reference_molar_volume();
    let v1 = volume_ratio * v0;
    let (mut lo, mut hi) = (200.0, 600.0);
    assert!(entropy_change(300.0, v0, lo, v1, CP_SLOPE) < 0.0);
    assert!(entropy_change(300.0, v0, hi, v1, CP_SLOPE) > 0.0);
    for _ in 0..80 {
        let mid = 0.5 * (lo + hi);
        if entropy_change(300.0, v0, mid, v1, CP_SLOPE) < 0.0 {
            lo = mid;
        } else {
            hi = mid;
        }
    }
    let point = pr_point(0.5 * (lo + hi), v1, CP_SLOPE);
    let t = point.temperature;
    let energy_change_j_mol = R * (2.5 * (t - 300.0) + 0.5 * CP_SLOPE * (t * t - 300.0 * 300.0))
        + point.residual_energy_j_mol
        - pr_point(300.0, v0, CP_SLOPE).residual_energy_j_mol;
    (point, energy_change_j_mol)
}

#[test]
fn g1_dense_cavity_force_potential_and_tangent_share_the_original_real_gas() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let gas = parcel(cavity.reference_volume_m3(), CP_SLOPE, 1e24);
        let saved = gas.clone();
        let displacement = affine_displacement(&nodes);
        let ambient = 1.2e6;
        let load = cavity
            .adiabatic_gas_load(&displacement, &gas, ambient, controls(), cx)
            .unwrap();
        let (volume, gradient) = tetra_volume_gradient(current_inner(&nodes, &displacement));
        let (oracle, energy_per_mole) = isentrope(volume / cavity.reference_volume_m3());
        let moles = cavity.reference_volume_m3() / reference_molar_volume();
        near(gas.mass_kg(), MOLAR_MASS * moles, 2e-11);
        near(load.volume_m3, volume, 2e-15);
        near(load.gas.state().temperature_k(), oracle.temperature, 3e-9);
        near(load.gas.state().pressure_pa(), oracle.pressure, 3e-4);
        near(
            load.potential_change_j,
            moles * energy_per_mole + ambient * (volume - cavity.reference_volume_m3()),
            5e-6,
        );
        let bulk = load.gas.state().density_kg_m3() * load.gas.state().sound_speed_m_s().powi(2);
        near(bulk, oracle.bulk_modulus, 8e-4);
        assert!(
            (bulk / (load.gas.state().gamma() * oracle.pressure) - 1.0).abs() > 0.02,
            "a dense-fluid tangent must not substitute gamma*p for rho*c^2"
        );
        assert_eq!(&load.force_n[..12], &[0.0; 12]);
        for i in 0..12 {
            near(
                load.force_n[12 + i],
                (oracle.pressure - ambient) * gradient[i],
                1e-4,
            );
            near(
                load.force_tangent.acoustic_rank_one_vector()[12 + i],
                (oracle.bulk_modulus / volume).sqrt() * gradient[i],
                2e-7,
            );
        }
        assert_eq!(load.gas.model(), gas.model());
        assert_eq!(load.gas.mass_kg().to_bits(), gas.mass_kg().to_bits());
        assert_eq!(load.gas.point().axes().get("material-lot"), Some(&47.0));

        let unshifted = parcel(cavity.reference_volume_m3(), CP_SLOPE, 0.0);
        let without_datum = cavity
            .adiabatic_gas_load(&displacement, &unshifted, ambient, controls(), cx)
            .unwrap();
        assert_eq!(load.potential_change_j, without_datum.potential_change_j);
        assert_eq!(load.force_n, without_datum.force_n);

        let direction: Vec<_> = (0..24)
            .map(|i| ((i * 7 % 17) as f64 - 8.0) / 25.0)
            .collect();
        let h = 2e-5;
        let plus: Vec<_> = displacement
            .iter()
            .zip(&direction)
            .map(|(u, d)| u + h * d)
            .collect();
        let minus: Vec<_> = displacement
            .iter()
            .zip(&direction)
            .map(|(u, d)| u - h * d)
            .collect();
        let plus = cavity
            .adiabatic_gas_load(&plus, &gas, ambient, controls(), cx)
            .unwrap();
        let minus = cavity
            .adiabatic_gas_load(&minus, &gas, ambient, controls(), cx)
            .unwrap();
        let work_derivative = pairing(&load.force_n, &direction);
        near(
            (plus.potential_change_j - minus.potential_change_j) / (2.0 * h),
            -work_derivative,
            2e-7 * work_derivative.abs(),
        );
        let action = load.force_tangent.apply(&direction, cx).unwrap();
        let scale = action.iter().map(|v| v.abs()).fold(0.0_f64, f64::max);
        for (i, actual) in action.iter().enumerate() {
            near(
                *actual,
                (plus.force_n[i] - minus.force_n[i]) / (2.0 * h),
                2e-7 * scale,
            );
        }
        assert_eq!(gas, saved);
    });
}

#[test]
fn g1_dense_gas_equilibrates_the_native_solid_with_independent_finite_compression() {
    with_cx(|_, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let stretch = 0.96_f64;
        let (target, energy_per_mole) = isentrope(stretch.powi(3));
        // At F=sI, the compressible Neo-Hookean law has
        // P=mu(s-1/s)+3 lambda ln(s)/s. Pick lambda=1.5 mu and the
        // independently computed real pressure to manufacture exact equilibrium.
        let mu = -target.pressure * stretch.powi(2)
            / (stretch - 1.0 / stretch + 4.5 * stretch.ln() / stretch);
        let lambda = 1.5 * mu;
        let material =
            Hyperelastic::new(HyperelasticModel::NeoHookean { mu, lambda }, 4.0).unwrap();
        let piola = mu * (stretch - 1.0 / stretch) + 3.0 * lambda * stretch.ln() / stretch;
        let ambient = 1.2e6;
        let outer_traction = piola + ambient * stretch.powi(2);
        let (_, outer_gradient) = tetra_volume_gradient([nodes[0], nodes[1], nodes[2], nodes[3]]);
        let fixed = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let dead: Vec<_> = (0..12)
            .filter(|dof| !fixed.iter().any(|&(id, _)| id == *dof))
            .map(|dof| (dof, outer_traction * outer_gradient[dof]))
            .collect();
        let gas = parcel(cavity.reference_volume_m3(), CP_SLOPE, 1e24);
        let saved = gas.clone();
        let problem = HyperTetProblem {
            nodes_m: &nodes,
            tetrahedra: &tets,
            material: HyperelasticLawField::Uniform(&material),
            prescribed_m: &fixed,
            nodal_forces_n: &dead,
            body_force_n_m3: [0.0; 3],
            budget: TetAssemblyBudget::standard(),
            settings: HyperTetSettings {
                load_steps: 8,
                ..HyperTetSettings::default()
            },
        };
        let solution = problem
            .solve_with_adiabatic_cavities(
                &[AdiabaticCavity {
                    cavity: &cavity,
                    reference: &gas,
                    exterior_pressure_pa: ambient,
                    controls: controls(),
                }],
                None,
                Some((cavity.surface(), TetPressure::Uniform(ambient))),
                cx,
            )
            .unwrap();
        let u: Vec<_> = solution
            .solid
            .displacement_m
            .iter()
            .flatten()
            .copied()
            .collect();
        for (&actual, &reference) in u.iter().zip(nodes.iter().flatten()) {
            near(actual, (stretch - 1.0) * reference, 4e-8);
        }
        near(solution.solid.min_det_f, stretch.powi(3), 4e-8);
        let load = &solution.cavities[0];
        near(load.gas.state().temperature_k(), target.temperature, 4e-6);
        near(load.gas.state().pressure_pa(), target.pressure, 0.3);
        near(
            load.volume_m3,
            cavity.reference_volume_m3() * stretch.powi(3),
            3e-9,
        );
        let moles = cavity.reference_volume_m3() / reference_molar_volume();
        near(
            load.potential_change_j,
            moles * energy_per_mole
                + ambient * cavity.reference_volume_m3() * (stretch.powi(3) - 1.0),
            0.03,
        );
        for &(dof, reaction) in &solution.solid.reactions_n {
            near(reaction, outer_traction * outer_gradient[dof], 0.5);
        }
        for axis in 0..3 {
            let reaction: f64 = solution
                .solid
                .reactions_n
                .iter()
                .filter(|(dof, _)| dof % 3 == axis)
                .map(|(_, value)| value)
                .sum();
            let applied: f64 = dead
                .iter()
                .filter(|(dof, _)| dof % 3 == axis)
                .map(|(_, value)| value)
                .sum();
            near(reaction + applied, 0.0, 0.5);
        }
        assert_eq!(load.gas.model(), gas.model());
        assert_eq!(load.gas.mass_kg().to_bits(), gas.mass_kg().to_bits());
        assert_eq!(gas, saved);
    });
}

#[test]
fn g4_real_cavity_source_refusal_and_cancellation_leave_the_inventory_unchanged() {
    with_cx(|gate, cx| {
        let (nodes, tets) = shell();
        let cavity = cavity(&nodes, &tets, cx);
        let gas = parcel(cavity.reference_volume_m3(), CP_SLOPE, 1e24);
        let saved = gas.clone();
        let displacement = affine_displacement(&nodes);
        let narrow = GasEntropySolve {
            maximum_temperature_k: 301.0,
            ..controls()
        };
        assert!(
            cavity
                .adiabatic_gas_load(&displacement, &gas, 0.0, narrow, cx)
                .is_err()
        );
        assert!(
            cavity
                .adiabatic_gas_load(&displacement, &gas, -1.0, controls(), cx)
                .is_err()
        );
        let foreign_volume = parcel(2.0 * cavity.reference_volume_m3(), CP_SLOPE, 0.0);
        assert!(
            cavity
                .adiabatic_gas_load(&displacement, &foreign_volume, 0.0, controls(), cx)
                .is_err()
        );
        assert_eq!(gas, saved);
        let accepted = cavity
            .adiabatic_gas_load(&displacement, &gas, 0.0, controls(), cx)
            .unwrap();
        assert!(accepted.gas.state().temperature_k() > 301.0);
        gate.request();
        assert!(matches!(
            cavity.adiabatic_gas_load(&displacement, &gas, 0.0, controls(), cx),
            Err(HyperTetError::Cancelled)
        ));
        assert_eq!(gas, saved);
    });
}

#[test]
fn g1_heating_one_dense_chamber_moves_the_other_without_merging_sources() {
    with_cx(|_, cx| {
        let (nodes, tets) = super::multiple::connected_shells();
        let cavities = super::multiple::cavities(&nodes, &tets, cx);
        let material = Hyperelastic::new(
            HyperelasticModel::NeoHookean {
                mu: 2e7,
                lambda: 3e7,
            },
            4.0,
        )
        .unwrap();
        let problem = super::multiple::problem(&nodes, &tets, &material, &[]);
        let slopes = [0.001, 0.002];
        let gases: [_; 2] = std::array::from_fn(|i| {
            parcel(
                cavities[i].reference_volume_m3(),
                slopes[i],
                if i == 0 { 0.0 } else { 1e24 },
            )
        });
        let sources: [_; 2] = std::array::from_fn(|i| AdiabaticCavity {
            cavity: &cavities[i],
            reference: &gases[i],
            exterior_pressure_pa: INITIAL_PRESSURE,
            controls: controls(),
        });
        let baseline = problem
            .solve_with_adiabatic_cavities(&sources, None, None, cx)
            .unwrap();
        assert_eq!(baseline.solid.displacement_m, vec![[0.0; 3]; nodes.len()]);

        let mut changed = gases.clone();
        let supplied_heat = 10_000.0 * changed[0].mass_kg();
        changed[0]
            .add_heat(
                supplied_heat,
                GasEnergySolve {
                    minimum_temperature_k: 200.0,
                    maximum_temperature_k: 600.0,
                    energy_tolerance_j_kg: 1e-7,
                    maximum_iterations: 128,
                },
                || false,
            )
            .unwrap();
        assert!(changed[0].state().temperature_k() > 300.0);
        let saved = changed.clone();
        let sources: [_; 2] = std::array::from_fn(|i| AdiabaticCavity {
            cavity: &cavities[i],
            reference: &changed[i],
            exterior_pressure_pa: INITIAL_PRESSURE,
            controls: controls(),
        });
        let response = problem
            .solve_with_adiabatic_cavities(&sources, None, None, cx)
            .unwrap();
        assert!(response.cavities[0].volume_m3 > baseline.cavities[0].volume_m3);
        assert!(
            (response.cavities[1].volume_m3 - baseline.cavities[1].volume_m3).abs() > 1e-9,
            "the second chamber must respond through the actual connected solid"
        );
        assert!((response.cavities[1].gas.state().temperature_k() - 300.0).abs() > 1e-6);
        for i in 0..2 {
            let endpoint = &response.cavities[i];
            let moles = changed[i].mass_kg() / MOLAR_MASS;
            let v0 = changed[i].volume_m3() / moles;
            let v1 = endpoint.volume_m3 / moles;
            let t0 = changed[i].state().temperature_k();
            let t1 = endpoint.gas.state().temperature_k();
            near(entropy_change(t0, v0, t1, v1, slopes[i]), 0.0, 2e-11);
            let initial = pr_point(t0, v0, slopes[i]);
            let oracle = pr_point(t1, v1, slopes[i]);
            near(endpoint.gas.state().pressure_pa(), oracle.pressure, 3e-4);
            let stored = moles
                * (R * (2.5 * (t1 - t0) + 0.5 * slopes[i] * (t1 * t1 - t0 * t0))
                    + oracle.residual_energy_j_mol
                    - initial.residual_energy_j_mol);
            near(
                endpoint.potential_change_j,
                stored + INITIAL_PRESSURE * (endpoint.volume_m3 - changed[i].volume_m3()),
                5e-6,
            );
            assert_eq!(
                endpoint.gas.mass_kg().to_bits(),
                gases[i].mass_kg().to_bits()
            );
            assert_eq!(endpoint.gas.model(), gases[i].model());
        }
        for axis in 0..3 {
            let reaction: f64 = response
                .solid
                .reactions_n
                .iter()
                .filter(|(dof, _)| dof % 3 == axis)
                .map(|(_, force)| force)
                .sum();
            near(reaction, 0.0, 0.003);
        }
        assert_eq!(changed, saved);
    });
}
