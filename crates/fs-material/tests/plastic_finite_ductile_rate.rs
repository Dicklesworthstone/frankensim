//! G1/G3: physical overstress, relaxation and progressive finite ductility.
//! Synthetic coefficients exercise the constitutive model, not an alloy fit.
#![deny(unsafe_code)]

use fs_material::ductile::evolution::{
    DifferentiableDuctileInitiationLaw, DuctileJ2Conditions, DuctileJ2Settings,
    DuctileLocusGradient,
};
use fs_material::ductile::{
    DuctileDomain, DuctileError, DuctileInitiationLaw, DuctileStressState, TabulatedDuctileLocus,
};
use fs_material::plastic::energy::EnergyTolerances;
use fs_material::plastic::finite::ductile::{FiniteDuctileError, FiniteDuctileJ2};
use fs_material::plastic::finite::rate::{FiniteOverstressJ2, FiniteOverstressParameters};
use fs_material::plastic::finite::{FiniteJ2, FiniteJ2Error, FiniteJ2Limits, FiniteJ2Parameters};
use fs_material::plastic::tabulated::{HardeningPoint, TabulatedHardening};

const MU: f64 = 500.0;
const Y: f64 = 3.0;
const H: f64 = 40.0;
const ETA: f64 = 400.0;
const ONSET: f64 = 0.002;

fn native() -> FiniteJ2 {
    FiniteJ2::new(
        FiniteJ2Parameters {
            shear_modulus_pa: MU,
            bulk_modulus_pa: 1000.0,
            yield_stress_pa: Y,
            hardening_modulus_pa: H,
        },
        FiniteJ2Limits {
            minimum_jacobian: 0.1,
            maximum_elastic_log_strain: 0.5,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    ).unwrap()
}
fn backbones() -> [FiniteJ2; 3] {
    let table = TabulatedHardening::new(
        [(0.0, Y), (0.01, 5.0), (0.03, 6.0), (0.1, 6.0), (3.0, 20.0)]
            .into_iter()
            .map(|(plastic_strain, yield_stress_pa)| HardeningPoint { plastic_strain, yield_stress_pa })
            .collect(),
    ).unwrap();
    [
        native(),
        native().with_voce_hardening(15.0, 12.0).unwrap(),
        FiniteJ2::new(
            FiniteJ2Parameters { hardening_modulus_pa: 0.0, ..*native().parameters() },
            *native().limits(),
        ).unwrap().with_tabulated_hardening(table).unwrap(),
    ]
}
fn parameters(exponent: f64) -> FiniteOverstressParameters {
    FiniteOverstressParameters { viscosity_pa_s: ETA, stress_scale_pa: 4.0, exponent }
}
fn domain() -> DuctileDomain {
    DuctileDomain {
        triaxiality: [-50.0, 50.0],
        strain_rate_s_inv: [1e-12, 1000.0],
        temperature_k: [300.0, 300.0],
    }
}
fn settings() -> DuctileJ2Settings {
    DuctileJ2Settings {
        characteristic_length_m: 0.02,
        post_initiation_plastic_work_j_m2: 0.004,
        maximum_damage_increment: 1.0,
        energy_tolerances: EnergyTolerances { energy_j_m3: 1e-10, yield_pa: 1e-8 },
    }
}
fn conditions(duration_s: f64) -> DuctileJ2Conditions {
    DuctileJ2Conditions { duration_s, temperature_k: 300.0 }
}
fn locus(onset: f64) -> TabulatedDuctileLocus {
    TabulatedDuctileLocus::new(vec![(-50.0, onset), (50.0, onset)], domain()).unwrap()
}
fn material(onset: f64) -> FiniteDuctileJ2<TabulatedDuctileLocus> {
    FiniteDuctileJ2::new(native(), locus(onset), settings()).unwrap()
        .with_overstress(parameters(1.0)).unwrap()
}
fn coaxial(e: f64, dilation: f64) -> [f64; 9] {
    [
        (dilation + e).exp(), 0.0, 0.0,
        0.0, (dilation - 0.5 * e).exp(), 0.0,
        0.0, 0.0, (dilation - 0.5 * e).exp(),
    ]
}
fn multiply(a: &[f64; 9], b: &[f64; 9]) -> [f64; 9] {
    core::array::from_fn(|ij| {
        (0..3).map(|k| a[3 * (ij / 3) + k] * b[3 * k + ij % 3]).sum()
    })
}
fn near(actual: f64, expected: f64, absolute: f64, relative: f64) {
    assert!(
        (actual - expected).abs() <= absolute + relative * actual.abs().max(expected.abs()),
        "{actual:.16e} != {expected:.16e}; difference {:.8e}", (actual - expected).abs(),
    );
}
fn equilibrium_work(alpha: f64) -> f64 {
    let post = (alpha - ONSET).max(0.0);
    post * (Y + H * (ONSET + 0.5 * post))
}
fn integrity(alpha: f64) -> f64 {
    (-settings().characteristic_length_m * equilibrium_work(alpha)
        / settings().post_initiation_plastic_work_j_m2).exp()
}

#[test]
fn linear_overstress_matches_independent_finite_damage_storage_and_heat() {
    let model = material(ONSET);
    let old = model.initial_state();
    let e = 0.04;
    let dt = 0.2;
    let f = coaxial(e, 0.0);
    let step = model.trial(&f, &old, conditions(dt)).unwrap();
    let alpha = (3.0 * MU * e - Y) / (3.0 * MU + H + ETA / dt);
    let viscous = ETA * alpha / dt;
    let g = integrity(alpha);
    let native_storage = 1.5 * MU * (e - alpha).powi(2) + 0.5 * H * alpha * alpha;
    near(step.state.native_state().alpha, alpha, 1e-14, 1e-13);
    near(step.state.integrity(), g, 2e-14, 2e-13);
    near(step.state.stored_energy_j_m3(), g * native_storage, 1e-13, 2e-12);
    near(step.energy.plastic_dissipation_j_m3, g * (Y + viscous) * alpha, 1e-13, 2e-12);
    near(step.energy.equivalent_plastic_work_increment_j_m2,
        settings().characteristic_length_m * equilibrium_work(alpha), 1e-15, 2e-12);
    // Viscous heat is physical, but is NOT part of the supplied equilibrium
    // resistance integral that controls post-initiation attenuation.
    assert!(step.energy.plastic_dissipation_j_m3 > 5.0 * g * Y * alpha);
    for (i, factor) in [(0, 2.0), (4, -1.0), (8, -1.0)] {
        let tau = g * factor * MU * (e - alpha);
        near(step.kirchhoff_pa[i], tau, 2e-12, 2e-12);
        near(step.first_piola_pa[i], tau / f[i], 2e-12, 2e-12);
    }
    assert_eq!(step.energy.damage_release_j_m3, 0.0);
    near(step.energy.residual_j_m3, 0.0, 1e-13, 0.0);
    assert_eq!(model.stress(&f, &old, conditions(dt)).unwrap(), step.first_piola_pa);
    assert!(matches!(
        model.native().validate_committed(&f, step.state.native_state(), 1e-8),
        Err(FiniteJ2Error::YieldResidual { .. }),
    ));
    model.validate_committed(&step.state, 1e-8).unwrap();
    assert_eq!(model.overstress_parameters(), Some(&parameters(1.0)));
    assert_eq!(old, model.initial_state());
    assert_eq!(step, model.trial(&f, &old, conditions(dt)).unwrap());
}

#[test]
fn actual_duration_changes_flow_strength_and_damage_at_the_same_deformation() {
    let model = material(ONSET);
    let old = model.initial_state();
    let f = coaxial(0.04, 0.0);
    let fast = model.trial(&f, &old, conditions(0.02)).unwrap();
    let slow = model.trial(&f, &old, conditions(2.0)).unwrap();
    assert!(fast.state.native_state().alpha < slow.state.native_state().alpha);
    assert!(fast.first_piola_pa[0] > slow.first_piola_pa[0]);
    assert!(fast.state.damage() < slow.state.damage());
    assert!(fast.state.damage() > 0.0);
    assert_eq!(old, model.initial_state());
}

#[test]
fn held_specimen_relaxes_and_continues_damage_with_convergent_physical_time() {
    let model = material(ONSET);
    let f = coaxial(0.04, 0.0);
    let loaded = model.trial(&f, &model.initial_state(), conditions(0.02)).unwrap();
    let snapshot = loaded.state.clone();
    model.validate_committed(&loaded.state, 1e-8).unwrap();
    assert_eq!(*loaded.state.first_piola_pa(), loaded.first_piola_pa);
    assert_eq!(*loaded.state.kirchhoff_pa(), loaded.kirchhoff_pa);
    assert_eq!(loaded.state, snapshot);
    let hold = model.trial(&f, &loaded.state, conditions(0.2)).unwrap();
    let before = loaded.state.native_state().alpha;
    let increment = (3.0 * MU * 0.04 - Y - (3.0 * MU + H) * before)
        / (3.0 * MU + H + ETA / 0.2);
    near(hold.state.native_state().alpha, before + increment, 1e-14, 2e-12);
    near(hold.state.integrity(), integrity(before + increment), 1e-14, 2e-12);
    assert!(hold.state.damage() > loaded.state.damage());
    assert!(hold.first_piola_pa[0] < loaded.first_piola_pa[0]);
    assert!(hold.state.stored_energy_j_m3() < loaded.state.stored_energy_j_m3());
    assert!(hold.energy.plastic_dissipation_j_m3 > 0.0);
    assert!(hold.energy.damage_release_j_m3 > 0.0);
    assert_eq!(hold.energy.mechanical_work_j_m3, 0.0);
    near(hold.energy.final_stored_j_m3 - hold.energy.initial_stored_j_m3
        + hold.energy.plastic_dissipation_j_m3 + hold.energy.damage_release_j_m3
        + hold.energy.endpoint_work_defect_j_m3, 0.0, 1e-12, 0.0);
    assert_eq!(model.stress(&f, &loaded.state, conditions(0.2)).unwrap(), hold.first_piola_pa);
    assert_eq!(loaded.state, snapshot);

    let duration = 0.2;
    let equilibrium_alpha = (3.0 * MU * 0.04 - Y) / (3.0 * MU + H);
    let exact = equilibrium_alpha + (before - equilibrium_alpha)
        * (-(3.0 * MU + H) * duration / ETA).exp();
    let run = |steps: usize| {
        let mut state = loaded.state.clone();
        for _ in 0..steps {
            state = model.trial(&f, &state, conditions(duration / steps as f64)).unwrap().state;
        }
        state
    };
    let coarse = run(10);
    let fine = run(100);
    assert!((fine.native_state().alpha - exact).abs() < (coarse.native_state().alpha - exact).abs() / 8.0);
    near(fine.integrity(), integrity(fine.native_state().alpha), 2e-13, 2e-11);
}

#[test]
fn intact_overstress_preserves_original_linear_voce_and_tabulated_returns() {
    for backbone in backbones() {
        for exponent in [1.0, 3.0] {
            let model = FiniteDuctileJ2::new(backbone.clone(), locus(10.0), settings()).unwrap()
                .with_overstress(parameters(exponent)).unwrap();
            let mut state = model.initial_state();
            for (e, dt) in [(0.003, 0.1), (0.017, 0.2), (0.05, 0.4), (0.05, 0.1), (0.04, 0.3)] {
                let f = coaxial(e, 0.001);
                let native = FiniteOverstressJ2::new(backbone.clone(), parameters(exponent), dt).unwrap();
                let original = native.response(&f, state.native_state()).unwrap();
                let step = model.trial(&f, &state, conditions(dt)).unwrap();
                assert_eq!(step.state.native_state(), &original.state);
                assert_eq!(step.first_piola_pa, original.first_piola_pa);
                assert_eq!(step.kirchhoff_pa, original.kirchhoff_pa);
                assert_eq!(step.tangent_pa, original.tangent_pa);
                assert_eq!(step.state.stored_energy_j_m3(), original.stored_energy_j_m3);
                assert_eq!(step.energy.plastic_dissipation_j_m3, original.plastic_dissipation_j_m3);
                assert_eq!(step.state.damage(), 0.0);
                model.validate_committed(&step.state, 1e-8).unwrap();
                state = step.state;
            }
        }
    }
}

#[derive(Clone, Debug, PartialEq)]
struct SmoothLocus(DuctileDomain);
impl DuctileInitiationLaw for SmoothLocus {
    fn domain(&self) -> &DuctileDomain { &self.0 }
    fn fracture_strain(&self, _: f64, _: f64, _: f64) -> Result<f64, DuctileError> {
        Err(DuctileError::Invalid("full stress state required"))
    }
    fn fracture_strain_at(&self, s: DuctileStressState, rate: f64, _: f64) -> Result<f64, DuctileError> {
        Ok(0.003 * (1.0 + 0.06 * s.triaxiality + 0.15 * s.lode_parameter
            + 0.025 * (rate / 0.01).ln()))
    }
}
impl DifferentiableDuctileInitiationLaw for SmoothLocus {
    fn fracture_strain_gradient(&self, _: f64, _: f64, _: f64) -> Result<DuctileLocusGradient, DuctileError> {
        Err(DuctileError::Invalid("full stress state required"))
    }
    fn fracture_strain_gradient_at(&self, _: DuctileStressState, rate: f64, _: f64) -> Result<DuctileLocusGradient, DuctileError> {
        Ok(DuctileLocusGradient {
            triaxiality: 0.003 * 0.06,
            lode_parameter: 0.003 * 0.15,
            strain_rate_s_inv: 0.003 * 0.025 / rate,
        })
    }
}

#[test]
fn all_81_tangent_entries_include_kinetics_moving_onset_and_relaxation_damage() {
    for backbone in backbones() {
        for exponent in [1.0, 3.0] {
            let model = FiniteDuctileJ2::new(backbone.clone(), SmoothLocus(domain()), settings()).unwrap()
                .with_overstress(parameters(exponent)).unwrap();
            let before = model.trial(&coaxial(0.003, 0.0), &model.initial_state(), conditions(0.2))
                .unwrap().state;
            assert!(!before.initiation().initiated());
            let f = [1.04, 0.017, -0.006, 0.004, 0.987, 0.009, 0.001, -0.003, 0.995];
            let crossed = model.trial(&f, &before, conditions(0.2)).unwrap();
            assert!(crossed.state.damage() > 0.0);
            let compressed = multiply(&coaxial(0.012, -0.04), &f);
            for (old, endpoint) in [(&before, f), (&crossed.state, f), (&crossed.state, compressed)] {
                let response = model.trial(&endpoint, old, conditions(0.2)).unwrap();
                assert!(response.state.native_state().alpha > old.native_state().alpha);
                assert!(response.state.damage() > old.damage());
                for j in 0..9 {
                    let h = 2e-7;
                    let mut plus = endpoint;
                    let mut minus = endpoint;
                    plus[j] += h;
                    minus[j] -= h;
                    let rp = model.stress(&plus, old, conditions(0.2)).unwrap();
                    let rm = model.stress(&minus, old, conditions(0.2)).unwrap();
                    for i in 0..9 {
                        near(response.tangent_pa[i][j], (rp[i] - rm[i]) / (2.0 * h), 3e-5, 3e-6);
                    }
                }
            }
        }
    }
}

#[test]
fn kinetic_damage_and_domain_refusals_leave_accepted_history_unchanged() {
    let model = material(ONSET);
    let old = model.initial_state();
    for dt in [0.0, -1.0, f64::NAN] {
        assert!(model.trial(old.deformation(), &old, conditions(dt)).is_err());
    }
    assert!(model.validate_committed(&old, f64::NAN).is_err());
    let other = model.clone();
    assert!(matches!(other.validate_committed(&old, 1e-8), Err(FiniteDuctileError::StateMismatch)));
    assert!(matches!(
        other.trial(old.deformation(), &old, conditions(0.2)),
        Err(FiniteDuctileError::StateMismatch),
    ));
    for p in [
        FiniteOverstressParameters { viscosity_pa_s: 0.0, ..parameters(1.0) },
        parameters(0.5),
    ] {
        assert!(model.clone().with_overstress(p).is_err());
    }
    let limited = FiniteDuctileJ2::new(native(), locus(ONSET), DuctileJ2Settings {
        maximum_damage_increment: 1e-5, ..settings()
    }).unwrap().with_overstress(parameters(1.0)).unwrap();
    let virgin = limited.initial_state();
    let f = coaxial(0.04, 0.0);
    assert!(matches!(limited.trial(&f, &virgin, conditions(0.2)), Err(FiniteDuctileError::DamageIncrement { .. })));
    assert!(matches!(limited.stress(&f, &virgin, conditions(0.2)), Err(FiniteDuctileError::DamageIncrement { .. })));
    assert_eq!(virgin, limited.initial_state());
    assert_eq!(old, model.initial_state());
}
