//! Declared continuum resolution changes the moving instrument before rendering.
use super::*;

const FLOW: &str = "cylinder-gas-loading 0.001 30 12";

fn text() -> String {
    include_str!("../cylinder-flow.material-strings")
        .replace("48000 4096 0.01 128", "48000 256 0.01 128")
}

fn small() -> input::Spec {
    input::Spec::parse(text().as_bytes()).unwrap()
}

#[test]
fn cylinder_data_precedes_surface_preload_and_reaches_both_planes_torsion_cavity_and_pcm() {
    scoped(|cx, gate| {
        let mut direct = compile(cx, gate, small(), 7800.0);
        let source = direct.performance.bank().source().bank();
        let loading = source.cylinder_gas_loading().unwrap();
        assert_eq!(loading.resolution(), small().cylinder_gas_loading.unwrap());
        assert_eq!(loading.strings().len(), 1);
        assert_eq!(loading.strings()[0].flow().gas(), source.gas());
        assert_eq!(loading.strings()[0].flow().radius_m(), 0.0003);
        assert_eq!(loading.strings()[0].flow().branches().len(), 12);
        assert_eq!(source.material_cavity().unwrap().cavity.gas(), source.gas());
        assert_eq!(direct.catalog.len(), 3);
        let MaterialStringSource::WithPolarizations(planes) = &source.strings()[0].source else {
            panic!("original polarization source lost");
        };
        assert_eq!(planes.directions().bow(), [0.6, 0.8]);
        let MaterialStringSource::WithTorsion(twist) = planes.source() else {
            panic!("original shear source lost");
        };
        let MaterialStringSource::Concentric { section, .. } = twist.bending() else {
            panic!("original layered section lost");
        };
        assert_eq!(section.elastic().layers().len(), 2);
        let loaded_mass = direct.performance.bank().cards()[0].linear_density_kg_m;
        assert_eq!(loaded_mass, loading.strings()[0].solid_card().linear_density_kg_m
            + loading.strings()[0].flow().added_mass_per_length_kg_m());

        let components = direct.performance.bank().system().components();
        assert_eq!(components.iter().map(|c| c.modes().len()).collect::<Vec<_>>(), vec![3, 2, 1, 36]);
        assert_eq!(components[0].modes()[0], components[0].modes()[1]);
        for state in &components[0].states()[..2] {
            assert_ne!(state.displacement_m_sqrt_kg, 0.0);
        }
        assert!(components[3].states().iter().all(|s|
            s.displacement_m_sqrt_kg == 0.0 && s.velocity_m_sqrt_kg_per_s == 0.0));
        assert!(direct.performance.projected_string_forces().unwrap()
            .initial_forces()[0][6..].iter().all(|f| *f == 0.0));
        assert!(direct.performance.preload_energy_j().unwrap() > 0.0);

        let expected = pressure(&mut direct, 37);
        assert!(expected.iter().any(|p| p.abs() > 1e-12));
        assert!(direct.performance.bank().system().components()[3].states().iter()
            .any(|s| s.velocity_m_sqrt_kg_per_s != 0.0),
            "the released string never transferred momentum to its fluid");
        let frame = direct.performance.bank().system().last_frame().unwrap();
        assert_eq!(frame.external_work_j, 0.0);
        assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        let mut replay = compile(cx, gate, small(), 7800.0);
        let (wav, summary) = replay.write(Cursor::new(Vec::<u8>::new()), gate).unwrap();
        let encoded = fs_couple::pcm_wav::encode_pcm16_wav(&expected, 48000, 0.01).unwrap();
        assert_eq!(wav.into_inner(), encoded.0);
        assert_eq!(summary.clipped_samples, encoded.1 as u64);
        assert_eq!(summary.samples, 256);
    });
}

#[test]
fn pressure_composition_and_core_material_are_independent_data_changes_in_the_same_state() {
    scoped(|cx, gate| {
        let mut original = compile(cx, gate, small(), 7800.0);
        let lowered = text().replace("mixture 293.15 101325", "mixture 293.15 60000");
        let mixture = text().replace("mole 2 N2 0.8 O2 0.2", "mole 3 N2 0.7 O2 0.2 CO2 0.1");
        let mut pressure_swap = compile(cx, gate, input::Spec::parse(lowered.as_bytes()).unwrap(), 7800.0);
        let mut composition_swap = compile(cx, gate, input::Spec::parse(mixture.as_bytes()).unwrap(), 7800.0);
        let mut material_swap = compile(cx, gate, small(), 15_600.0);
        let a = &original.performance.bank().source().bank().cylinder_gas_loading().unwrap().strings()[0];
        let b = &pressure_swap.performance.bank().source().bank().cylinder_gas_loading().unwrap().strings()[0];
        let c = &composition_swap.performance.bank().source().bank().cylinder_gas_loading().unwrap().strings()[0];
        let d = &material_swap.performance.bank().source().bank().cylinder_gas_loading().unwrap().strings()[0];
        assert_eq!(a.solid_card().linear_density_kg_m, b.solid_card().linear_density_kg_m);
        assert_eq!(a.solid_card().linear_density_kg_m, c.solid_card().linear_density_kg_m);
        assert_eq!(a.flow(), d.flow(), "a solid-material swap must not invent another atmosphere");
        assert!(d.solid_card().linear_density_kg_m > a.solid_card().linear_density_kg_m);
        assert!(b.flow().added_mass_per_length_kg_m() < a.flow().added_mass_per_length_kg_m());
        assert!(b.flow().branches()[0].rate_per_s > a.flow().branches()[0].rate_per_s);
        assert_ne!(c.flow().branches(), a.flow().branches());
        assert!(pressure_swap.performance.bank().cards()[0].fundamental_hz()
            > original.performance.bank().cards()[0].fundamental_hz());
        assert!(material_swap.performance.bank().cards()[0].fundamental_hz()
            < original.performance.bank().cards()[0].fundamental_hz());
        let p = pressure(&mut original, 128);
        for changed in [&mut pressure_swap, &mut composition_swap, &mut material_swap] {
            let q = pressure(changed, 37);
            assert!(p.iter().zip(&q).any(|(p, q)| (p-q).abs() > 1e-15));
            assert_ne!(original.performance.bank().system().components()[1].states(),
                changed.performance.bank().system().components()[1].states(),
                "only an observer scaling changed; actual shared body motion did not");
            let frame = changed.performance.bank().system().last_frame().unwrap();
            assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
        }
    });
}

#[test]
fn causal_gas_loading_keeps_construction_loss_for_homogeneous_bonded_and_wound_strings() {
    scoped(|cx, gate| {
        let row = input::tension_loss_rows().unwrap().into_iter()
            .find(|row| row.id == "lampis-2024-dominant-cello-g").unwrap();
        let bonded = text();
        let homogeneous = bonded.replace(
            "section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003", "section homogeneous");
        let wound = bonded.replace(
            "section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003",
            "section close-wound 0.0002 wrap.fsmatpk 0.00005 0.000101")
            .replace("tension 20 1 1", "tension 20 1 0");
        for input in [homogeneous, bonded, wound] {
            let lossy = input.replace(FLOW, &format!("tension-loss {}\n{FLOW}", row.id));
            let mut a = compile(cx, gate, input::Spec::parse(input.as_bytes()).unwrap(), 7800.0);
            let mut b = compile(cx, gate, input::Spec::parse(lossy.as_bytes()).unwrap(), 7800.0);
            let (plain, loss) = (&a.performance.bank().cards()[0], &b.performance.bank().cards()[0]);
            assert_eq!(plain.linear_density_kg_m, loss.linear_density_kg_m);
            assert_eq!(plain.bending_stiffness_n_m2, loss.bending_stiffness_n_m2);
            assert_eq!(plain.viscous_bending_n_m2_s, loss.viscous_bending_n_m2_s);
            for k in 1..=loss.mode_count {
                let wave = k as f64 * core::f64::consts::PI / loss.length_m;
                let omega = loss.mode_omega_rad_s(k);
                assert_eq!(omega, plain.mode_omega_rad_s(k));
                // eta_F multiplies the tensile stiffness, not the air's inertia
                // or its memory. Both polarizations retain that same loss.
                let expected = 0.5 * row.eta_f * loss.tension_n * wave * wave
                    / (loss.linear_density_kg_m * omega * omega);
                assert!(((loss.zetas[k-1] - plain.zetas[k-1]) / expected - 1.0).abs() < 1e-10);
                let modes = b.performance.bank().system().components()[0].modes();
                assert_eq!(modes[k-1].damping_ratio, loss.zetas[k-1]);
                assert_eq!(modes[k-1+loss.mode_count].damping_ratio, loss.zetas[k-1]);
            }
            let p = pressure(&mut a, 128);
            let q = pressure(&mut b, 37);
            assert!(p.iter().zip(&q).any(|(p, q)| (p-q).abs() > 1e-15),
                "tension loss disappeared before the actual body-pressure output");
            for ready in [&a, &b] {
                let frame = ready.performance.bank().system().last_frame().unwrap();
                assert!(frame.energy_residual_j.abs() <= frame.energy_tolerance_j);
            }
        }
    });
}

#[test]
fn fluid_count_is_a_live_axis_and_complete_dimensions_preserve_explicit_work_limits() {
    scoped(|cx, gate| {
        let input = small();
        assert_eq!(input.network_dimensions(false).unwrap(), (42, 39, 42*40*40));
        assert_eq!(input.network_dimensions(true).unwrap(), (42, 41, 42*42*42));
        let finer = input::Spec::parse(text().replace(FLOW,
            "cylinder-gas-loading 0.001 30 16").as_bytes()).unwrap();
        assert_eq!(finer.network_dimensions(false).unwrap(), (54, 51, 54*52*52));
        let mut coarse = compile(cx, gate, input, 7800.0);
        let mut refined = compile(cx, gate, finer, 7800.0);
        assert_eq!(coarse.performance.bank().cards()[0].fundamental_hz(),
            refined.performance.bank().cards()[0].fundamental_hz());
        assert_eq!(refined.performance.bank().system().components()[3].modes().len(), 48);
        assert!(pressure(&mut coarse, 128).iter().zip(pressure(&mut refined, 37))
            .any(|(p, q)| (p-q).abs() > 1e-17));

        // Count fits the file's 128-mode ceiling. Giving BOTH planes a body
        // projection makes one >42-mode exact island, which the owner refuses.
        let connected = text().replace(FLOW, "cylinder-gas-loading 0.001 30 24")
            .replace("two 0.6 0.8 1 0", "two 0.6 0.8 0.6 0.8");
        let connected = input::Spec::parse(connected.as_bytes()).unwrap();
        assert_eq!(connected.network_dimensions(false).unwrap().0, 78);
        let error = connected.compile(cx, gate, |name| Ok(pack(name, 7800.0, 200_000.0, None)))
            .err().expect("the legacy exact owner must refuse an oversized connected island").to_string();
        assert!(error.contains("dense exponential block capacity"), "{error}");
    });
}

#[test]
fn optional_flow_record_is_strict_ordered_and_never_truncates_declared_memory() {
    let bytes = text();
    assert!(input::Spec::parse(bytes.replace(&format!("{FLOW}\n"), "").as_bytes())
        .unwrap().cylinder_gas_loading.is_none());
    let v2 = bytes.replace("frankensim-material-strings-v3", "frankensim-material-strings-v2")
        .replace("polarizations two 0.6 0.8 1 0\n", "");
    assert_eq!(input::Spec::parse(v2.as_bytes()).unwrap().cylinder_gas_loading,
        small().cylinder_gas_loading);
    let v1 = v2.replace("frankensim-material-strings-v2", "frankensim-material-strings-v1")
        .replace("tension 20", "20")
        .replace("section bonded 0.0002 1\nlayer sleeve.fsmatpk 0.0003\n", "");
    assert_eq!(input::Spec::parse(v1.as_bytes()).unwrap().cylinder_gas_loading,
        small().cylinder_gas_loading);
    for invalid in [
        "cylinder-gas-loading 0 30 12", "cylinder-gas-loading 0.6 30 12",
        "cylinder-gas-loading 0.001 0.0001 12", "cylinder-gas-loading 0.001 3001 12",
        "cylinder-gas-loading 0.001 30 0", "cylinder-gas-loading 0.001 30 4097",
        "cylinder-gas-loading NaN 30 12", "cylinder-gas-loading 0.001 30 12 ignored",
        "cylinder-gas-loading 0.001 30 64",
    ] {
        assert!(input::Spec::parse(bytes.replace(FLOW, invalid).as_bytes()).is_err(), "accepted {invalid}");
    }
    assert!(input::Spec::parse(bytes.replace(FLOW, &format!("{FLOW}\n{FLOW}")).as_bytes()).is_err());
    let reordered = bytes.replace(&format!("{FLOW}\n"), "")
        .replace("friction stribeck", &format!("{FLOW}\nfriction stribeck"));
    assert!(input::Spec::parse(reordered.as_bytes()).is_err());
}
