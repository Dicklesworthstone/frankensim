//! Play keys of the embedded Model D with strings built from a complete string
//! construction and compiled material packs, and report what the materials
//! predict next to what the engine renders (native).
//!
//! Per key: tension, core stress utilization, the inharmonicity coefficient
//! B = EI pi^2/(T L^2), the first partials' frequencies, and each partial's
//! amplitude decay rate from the reference bank's string loss (material
//! bending loss, thermoelastic where declared) plus the gas drag of the same
//! atmosphere. Then the key is struck at 2 m/s and rendered for two seconds at
//! 48 kHz: peak pressure, RMS over 0.5-1.0 s and 1.5-2.0 s, and the energy
//! ledger. The render consumes these same coefficients; nothing is tuned.
//!
//! usage: string_materials_render <construction.csv> <bending-loss property> <keys, e.g. 21,60,96> <pack>...
use fs_live_piano::{Config, Description, PianoAudio, preparation::{self, materials::{StringMaterials, operating_point}},
    refmods::{bank_string_losses, linear::Bank}};
use std::sync::Arc;

fn main() {
    let args: Vec<String> = std::env::args().skip(1).collect();
    assert!(args.len() >= 4, "usage: string_materials_render <construction.csv> <loss property> <keys> <pack>...");
    let construction = std::fs::read_to_string(&args[0]).expect("construction file");
    let keys: Vec<u8> = args[2].split(',').map(|k| k.trim().parse().expect("key")).collect();
    let packs: Vec<Vec<u8>> = args[3..].iter().map(|path| std::fs::read(path).expect("pack file")).collect();
    let slices: Vec<&[u8]> = packs.iter().map(Vec::as_slice).collect();
    let point = operating_point("T", 293.15, "frequency", 1000.0).expect("operating point");
    let materials = match StringMaterials::from_packs(&slices, &construction, &point, &args[1]) {
        Ok(m) => Arc::new(m),
        Err(e) => { println!("{{\"constructionRefused\":{e:?}}}"); std::process::exit(2); }
    };
    let mut desc = Description::embedded().expect("embedded description");
    if let Err(e) = preparation::prepare_with_materials(&mut desc, None, None, Some(&materials)) {
        println!("{{\"applyRefused\":{e:?}}}");
        std::process::exit(3);
    }
    let report = materials.strength_report(&desc.courses).expect("strength report");
    let gas = fs_material::gas::GasState::try_new(&fs_material::gas::GasSpec::dry_air_ussa1976(), 293.15, 101_325.0)
        .expect("gas");
    let cfg = Config::live(48_000);
    // The live kernel's own bank arguments (kernel.rs build): band 0.45 x rate.
    let bank = Bank::new_with_string_damping(&desc.courses, &desc.board, 48_000, 0.45 * 48_000.0,
        cfg.max_partials, true, None, false).expect("bank");
    let (losses, _) = bank_string_losses(&bank);
    println!("{{\"material\":{},\n\"keys\":[", materials.json());
    for (n, &key) in keys.iter().enumerate() {
        let index = desc.courses.iter().position(|c| c.midi == key).expect("key in scale");
        let c = desc.courses[index];
        let stress = report.keys[index];
        let b = c.flexural_rigidity_nm2 * (std::f64::consts::PI / c.length_m).powi(2) / c.tension_n;
        let string = bank.strings.iter().find(|s| s.course == index).expect("speaking string");
        let mut audio = PianoAudio::new_with_materials(&desc, cfg.clone(), 1, gas, materials.clone()).expect("engine");
        let gas_rates = audio.piano().string_gas_drag_rates_per_s(key).map(<[f64]>::to_vec).unwrap_or_default();
        let partials = string.modes.clone().take(8).enumerate().map(|(i, m)| {
            let rate = losses[m].mass_rate + losses[m].strain_rate + gas_rates.get(i).copied().unwrap_or(0.0);
            format!("{{\"hz\":{:.3},\"materialRatePerS\":{:.6e},\"gasRatePerS\":{:.6e},\"amplitudeDecayPerS\":{:.6e}}}",
                bank.modes[m].omega / std::f64::consts::TAU, losses[m].strain_rate + losses[m].mass_rate,
                gas_rates.get(i).copied().unwrap_or(0.0), 0.5 * rate)
        }).collect::<Vec<_>>().join(",");
        audio.note_on(key, 2.0).expect("strike");
        let mut pcm = vec![0.0_f32; 96_000];
        audio.render(&mut pcm).expect("render");
        let rms = |a: usize, b: usize| (pcm[a..b].iter().map(|x| f64::from(*x).powi(2)).sum::<f64>() / (b - a) as f64).sqrt();
        let peak = pcm.iter().fold(0.0_f64, |m, x| m.max(f64::from(x.abs())));
        let finite = pcm.iter().all(|x| x.is_finite());
        let ledger = audio.close_ledger();
        println!("{}{{\"key\":{key},\"tensionN\":{:.3},\"utilization\":{:.4},\"fMaxHz\":{:.2},\"b\":{:.6e},\"f1Hz\":{:.4},\
            \"partials\":[{partials}],\"render\":{{\"finite\":{finite},\"peak\":{peak:.4e},\"rms0p5to1s\":{:.4e},\
            \"rms1p5to2s\":{:.4e},\"inputJ\":{:.4e},\"stringLossJ\":{:.4e},\"residualJ\":{:.3e}}}}}",
            if n == 0 { "" } else { ",\n" }, c.tension_n, stress.utilization, stress.f_max_hz, b,
            c.partial_hz(1, c.tension_n), rms(24_000, 48_000), rms(72_000, 96_000), ledger.input_j,
            ledger.string_loss_j, ledger.residual_j);
    }
    println!("]}}");
}
