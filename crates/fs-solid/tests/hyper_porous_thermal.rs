//! G0/G1/G3/G4: retained porous thermal feedback and original reference inertia.

#[path = "hyper_porous_thermal/pressure.rs"]
mod pressure;

use fs_alloc::{ArenaConfig, ArenaPool};
use fs_exec::{Budget, CancelGate, Cx, ExecMode, StreamKey};
use fs_material::plastic::finite::porous::thermal::{
    AffineThermalFinitePorous, FinitePorousTemperatureSlopes, FinitePorousThermalParameters,
};
use fs_material::plastic::finite::porous::{FinitePorousLimits, FinitePorousParameters};
use fs_solid::hyper3::plastic::FinitePlasticStepControls;
use fs_solid::hyper3::plastic::porous::thermal::{
    FiniteThermalPorousTetSolid, FiniteThermoplasticControls,
};
use fs_solid::{HyperTetSettings, TetAssemblyBudget};

const NODES: [[f64; 3]; 4] = [[0.0; 3], [1.0, 0.0, 0.0], [0.0, 1.0, 0.0], [0.0, 0.0, 1.0]];
const TETS: [[usize; 4]; 1] = [[0, 1, 2, 3]];
const V: f64 = 1.0 / 6.0;

fn with_cx(f: impl FnOnce(&Cx<'_>, &CancelGate)) {
    let gate = CancelGate::new();
    ArenaPool::new(ArenaConfig::default()).scope(|arena| {
        f(
            &Cx::new(
                &gate,
                arena,
                StreamKey {
                    seed: 11,
                    kernel_id: 43,
                    tile: 0,
                    iteration: 0,
                },
                Budget::INFINITE,
                ExecMode::Deterministic,
            ),
            &gate,
        )
    });
}
fn family() -> AffineThermalFinitePorous {
    AffineThermalFinitePorous::new(
        FinitePorousThermalParameters {
            reference_temperature_k: 300.0,
            temperature_range_k: [250.0, 550.0],
            reference: FinitePorousParameters {
                shear_modulus_pa: 9000.0,
                bulk_modulus_pa: 20000.0,
                yield_stress_pa: 200.0,
                hardening_modulus_pa: 1200.0,
                pressure_sensitivity: 1.5,
            },
            slopes: FinitePorousTemperatureSlopes {
                shear_modulus_pa_k: -6.0,
                bulk_modulus_pa_k: -10.0,
                yield_stress_pa_k: -0.3,
                hardening_modulus_pa_k: -1.0,
            },
            logarithmic_expansion_per_k: 3e-4,
            logarithmic_expansion_slope_per_k2: 1e-7,
        },
        FinitePorousLimits {
            minimum_jacobian: 0.2,
            maximum_elastic_log_strain: 0.7,
            maximum_plastic_log_strain: 2.0,
            maximum_accumulated_plastic_strain: 3.0,
        },
    )
    .unwrap()
}
fn controls() -> FiniteThermoplasticControls {
    FiniteThermoplasticControls {
        material: FinitePlasticStepControls {
            equilibrium_work_tolerance_j: 1e-6,
            maximum_endpoint_work_defect_j: 200.0,
            yield_tolerance_pa: 1e-5,
        },
        maximum_numerical_kinetic_loss_j: 10.0,
    }
}
fn settings() -> HyperTetSettings {
    HyperTetSettings {
        load_steps: 1,
        relative_tolerance: 1e-11,
        ..HyperTetSettings::default()
    }
}
fn close(a: f64, b: f64, tolerance: f64) {
    assert!(
        a.is_finite() && b.is_finite() && (a - b).abs() <= tolerance,
        "{a:.16e} != {b:.16e}, allowance {tolerance:e}"
    );
}
fn initial(cx: &Cx<'_>) -> FiniteThermalPorousTetSolid<'static> {
    FiniteThermalPorousTetSolid::try_new(
        &NODES,
        &TETS,
        vec![family()],
        &[300.0],
        TetAssemblyBudget::standard(),
        cx,
    )
    .unwrap()
}
fn prescribed(stretch: f64) -> Vec<(usize, f64)> {
    NODES
        .iter()
        .enumerate()
        .flat_map(|(node, x)| (0..3).map(move |axis| (3 * node + axis, (stretch - 1.0) * x[axis])))
        .collect()
}
fn thermal_log(t: f64) -> f64 {
    let d = t - 300.0;
    3e-4 * d + 0.5e-7 * d * d
}

#[test]
fn free_thermal_expansion_retains_plastic_volume_and_conserves_supplied_mass() {
    with_cx(|cx, _| {
        let old = initial(cx);
        let supports = [(0, 0.0), (1, 0.0), (2, 0.0), (4, 0.0), (5, 0.0), (8, 0.0)];
        let step = old
            .advance(
                &[500.0],
                &supports,
                &[],
                [0.0; 3],
                0.2,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let theta = thermal_log(500.0).exp();
        for (i, x) in NODES.iter().enumerate() {
            for axis in 0..3 {
                close(
                    step.solid.displacement_m()[i][axis],
                    (theta - 1.0) * x[axis],
                    1e-9,
                );
            }
        }
        assert!(step.equilibrium.steps[0].linear_iterations > 0);
        close(step.solid.element_jacobians()[0], theta.powi(3), 1e-8);
        close(step.solid.element_plastic_jacobians()[0], 1.0, 1e-12);
        assert_eq!(step.solid.element_states(), old.element_states());
        assert_eq!(step.solid.families(), old.families());
        assert_eq!(step.solid.temperatures_k(), &[500.0]);
        assert_eq!(step.energy.plastic_dissipation_j, 0.0);
        close(step.energy.caloric_heat_j, 0.0, 1e-5);
        close(step.energy.balance_error_j, 0.0, 1e-6);
        close(
            step.solid.current_density_kg_m3(&[2.0], cx).unwrap()[0],
            2.0 / theta.powi(3),
            1e-8,
        );
        assert!(
            FiniteThermalPorousTetSolid::try_new(
                &NODES,
                &TETS,
                vec![family()],
                &[500.0],
                TetAssemblyBudget::standard(),
                cx
            )
            .is_err()
        );
        assert!(old.current_density_kg_m3(&[], cx).is_err());
    });
}

#[derive(Clone, Copy)]
struct Hydro {
    stretch: f64,
    t: f64,
    log_jp: f64,
    alpha: f64,
    g: f64,
    bulk: f64,
    hardening: f64,
    free: f64,
    entropy: f64,
    internal: f64,
    piola: f64,
    heat: f64,
}
impl Hydro {
    // Independent scalar associated hydrostatic return: d ln(Jp)=sign(g)*s*dalpha,
    // s=a/sqrt(1+a²/9). No native constitutive response enters this oracle.
    fn at(log_stretch: f64, t: f64, old: Option<Self>) -> Self {
        let d = t - 300.0;
        let (bulk, yield_stress, hardening) = (20000.0 - 10.0 * d, 200.0 - 0.3 * d, 1200.0 - d);
        let s = 1.5 / 1.25_f64.sqrt();
        let (mut log_jp, mut alpha) = old.map_or((0.0, 0.0), |old| (old.log_jp, old.alpha));
        let gm = 3.0 * (log_stretch - thermal_log(t));
        let trial = gm - log_jp;
        let increment = ((s * bulk * trial.abs() - yield_stress - hardening * alpha)
            / (bulk * s * s + hardening))
            .max(0.0);
        log_jp += trial.signum() * s * increment;
        alpha += increment;
        let g = gm - log_jp;
        let free = 0.5 * bulk * g * g + 0.5 * hardening * alpha * alpha;
        let beta = 3e-4 + 1e-7 * d;
        let entropy = 3.0 * bulk * beta * g + 5.0 * g * g + 0.5 * alpha * alpha;
        Self {
            stretch: log_stretch.exp(),
            t,
            log_jp,
            alpha,
            g,
            bulk,
            hardening,
            free,
            entropy,
            internal: free + t * entropy,
            piola: bulk * g / log_stretch.exp(),
            heat: yield_stress * increment,
        }
    }
    fn at_load(piola: f64, t: f64, old: Self) -> Self {
        let (mut lo, mut hi) = (-0.35, 0.2);
        assert!(Self::at(lo, t, Some(old)).piola < piola);
        assert!(Self::at(hi, t, Some(old)).piola > piola);
        for _ in 0..100 {
            let mid = 0.5 * (lo + hi);
            if Self::at(mid, t, Some(old)).piola < piola {
                lo = mid;
            } else {
                hi = mid;
            }
        }
        Self::at(0.5 * (lo + hi), t, Some(old))
    }
    fn energy(self, old: Self) -> [f64; 4] {
        let work = 3.0 * self.piola * (self.stretch - old.stretch);
        let theta = thermal_log(self.t).exp();
        let theta0 = thermal_log(old.t).exp();
        let native_work =
            3.0 * (self.piola * theta) * (self.stretch / theta - old.stretch / theta0);
        let coefficients = 0.5 * (self.bulk - old.bulk) * old.g * old.g
            + 0.5 * (self.hardening - old.hardening) * old.alpha * old.alpha;
        let ts = self.t * self.entropy - old.t * old.entropy;
        let reversible = work - native_work - coefficients - ts;
        let defect = native_work - (self.free - old.free - coefficients) - self.heat;
        [work, reversible, defect, coefficients]
    }
}

#[test]
fn force_driven_crush_heating_and_unloading_follow_independent_thermal_compaction() {
    with_cx(|cx, _| {
        let mut solid = initial(cx);
        let mut old = Hydro::at(0.0, 300.0, None);
        let target = Hydro::at(-0.025, 300.0, Some(old)).piola;
        let supports: Vec<_> = (0..12)
            .filter(|i| ![3, 7, 11].contains(i))
            .map(|i| (i, 0.0))
            .collect();
        let mut hot_alpha = 0.0;
        for (index, (t, p)) in [(300.0, target), (450.0, target), (450.0, 0.0)]
            .into_iter()
            .enumerate()
        {
            let expected = Hydro::at_load(p, t, old);
            let forces = [(3, V * p), (7, V * p), (11, V * p)];
            let step = solid
                .advance(
                    &[t],
                    &supports,
                    &forces,
                    [0.0; 3],
                    0.2,
                    settings(),
                    controls(),
                    cx,
                )
                .unwrap();
            close(step.solid.element_states()[0].alpha, expected.alpha, 1e-9);
            close(
                step.solid.element_jacobians()[0],
                expected.stretch.powi(3),
                1e-9,
            );
            close(
                step.solid.element_plastic_jacobians()[0],
                expected.log_jp.exp(),
                1e-9,
            );
            for i in [0, 4, 8] {
                close(
                    step.solid.element_first_piola_pa()[0][i],
                    expected.piola,
                    1e-6,
                );
            }
            close(step.solid.free_energy_j(), V * expected.free, 1e-7);
            close(step.solid.internal_energy_j(), V * expected.internal, 1e-6);
            close(step.energy.plastic_dissipation_j, V * expected.heat, 1e-7);
            let [work, reversible, defect, coefficients] = expected.energy(old);
            close(step.energy.dead_load_endpoint_work_j, V * work, 1e-7);
            close(step.energy.reversible_heat_j, V * reversible, 1e-6);
            close(step.energy.endpoint_work_defect_j, V * defect, 1e-7);
            close(
                step.energy.coefficient_storage_change_j,
                V * coefficients,
                1e-7,
            );
            close(step.energy.balance_error_j, 0.0, 1e-6);
            if index == 1 {
                assert!(expected.alpha > old.alpha);
                assert!(expected.log_jp < old.log_jp);
                hot_alpha = expected.alpha;
            }
            if index == 2 {
                close(step.energy.plastic_dissipation_j, 0.0, 1e-10);
                close(step.solid.element_states()[0].alpha, hot_alpha, 1e-9);
            }
            old = expected;
            solid = step.solid;
        }
        assert!(solid.element_plastic_jacobians()[0] < 0.9);
    });
}

#[test]
fn graded_inertial_thermal_crush_preserves_mass_and_rejects_late_loss_atomically() {
    with_cx(|cx, gate| {
        let rho = [2.0, 4.0, 6.0, 8.0];
        let solid = initial(cx)
            .with_nodal_inertia(&rho, &[[0.0; 3]; 4], cx)
            .unwrap();
        let stretch = (-0.02_f64).exp();
        let held = prescribed(stretch);
        let dt = 0.2;
        let step = solid
            .advance(
                &[350.0],
                &held,
                &[],
                [0.0; 3],
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        let sum: f64 = rho.iter().sum();
        let mass: [[f64; 4]; 4] = core::array::from_fn(|i| {
            core::array::from_fn(|j| {
                if i == j {
                    V * (sum + 2.0 * rho[i]) / 60.0
                } else {
                    V * (sum + rho[i] + rho[j]) / 120.0
                }
            })
        });
        let velocities: [[f64; 3]; 4] = NODES.map(|x| x.map(|v| (stretch - 1.0) * v / dt));
        let mut kinetic = 0.0;
        for i in 0..4 {
            for j in 0..4 {
                for axis in 0..3 {
                    kinetic += 0.5 * velocities[i][axis] * mass[i][j] * velocities[j][axis];
                }
            }
        }
        let inertia = step.inertia.as_ref().unwrap();
        close(inertia.kinetic_change_j, kinetic, 1e-12);
        close(inertia.numerical_loss_j, kinetic, 1e-12);
        close(
            step.solid.reference_mass_kg().unwrap(),
            V * sum / 4.0,
            1e-12,
        );
        assert_eq!(
            step.solid.reference_nodal_density_kg_m3(),
            Some(rho.as_slice())
        );
        close(
            step.solid.current_density_kg_m3(&[sum / 4.0], cx).unwrap()[0],
            sum / 4.0 / stretch.powi(3),
            1e-12,
        );
        assert!(step.solid.element_plastic_jacobians()[0] < 1.0);
        close(step.energy.balance_error_j, 0.0, 1e-6);
        let old_states = solid.element_states().to_vec();
        let old_t = solid.temperatures_k().to_vec();
        let strict = FiniteThermoplasticControls {
            maximum_numerical_kinetic_loss_j: 0.0,
            ..controls()
        };
        assert!(
            solid
                .advance(&[350.0], &held, &[], [0.0; 3], dt, settings(), strict, cx)
                .is_err()
        );
        assert_eq!(solid.element_states(), old_states);
        assert_eq!(solid.temperatures_k(), old_t);
        let retry = solid
            .advance(
                &[350.0],
                &held,
                &[],
                [0.0; 3],
                dt,
                settings(),
                controls(),
                cx,
            )
            .unwrap();
        assert_eq!(retry.solid.element_states(), step.solid.element_states());
        assert_eq!(retry.energy, step.energy);
        assert!(
            solid
                .advance(
                    &[551.0],
                    &held,
                    &[],
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.element_states(), old_states);
        gate.request();
        assert!(
            solid
                .advance(
                    &[350.0],
                    &held,
                    &[],
                    [0.0; 3],
                    dt,
                    settings(),
                    controls(),
                    cx
                )
                .is_err()
        );
        assert_eq!(solid.element_states(), old_states);
    });
}
