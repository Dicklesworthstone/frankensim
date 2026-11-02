//! pexfp.378's byte-identity evidence (MB 12550): a piece without a closed-body table renders byte
//! for byte as at the pin. Cymbals and hats are not re-rendered for the observer; these digests of
//! the first 60 ms of `kit::render` (after each piece's settling) are the evidence. The expected
//! digests were printed by this same file built against the pin's fs-vkplate (origin/main before
//! the observer existed), so a change to the plain render loop, the microphone step or the
//! preparation shows here.
use super::*;

const SECONDS: f64 = 0.06;

/// (piece, BLAKE3 of the little-endian f64 bits of every channel, in order) at the pin.
const PIN: [(&str, &str); 7] = [
    ("hihat-closed", "08f1dff1780d651a24388f0c8bdbe24e3fdfe0e546bb906b04911388c1d3ef2d"),
    ("hihat-open", "de9a36e54cea3fcfef8ec84ba09eb4c37eb1206808e3415dd382a60723f71180"),
    ("splash", "004ec059777e0a68f81c0434fb028f9197fba2e52a7771bdf71c7e493356f130"),
    ("crash", "4183cdf76c08dc38b95ea5f426ba35cd335b791a00c5c773cd1a790c44c9cabb"),
    ("ride", "8a25b23d3bca7e148929b7b194f2610ff175e1827900edc111b7f9f5b7380870"),
    ("tom-12 (no table)", "a40b9136a55d3075de9a7e211c4001f7efe0eeff217a42f3aae6e26ee941744a"),
    ("snare (no table)", "9a8104fae7ac0684e30ea7763a6e054316a11ba62d369fae0ecc92efa2664792"),
];

fn piece(name: &str) -> Result<Piece, String> {
    match name {
        "hihat-closed" => drum_cards::hihat(0.0, 80.0, &[(0.005, 2.5)], false),
        "hihat-open" => drum_cards::open_hihat(&[(0.005, 2.5)], false),
        "splash" => drum_cards::cymbal(&kit::splash(), &[(0.005, 2.0, 0.78, 0.005)], false),
        "crash" => drum_cards::cymbal(&kit::crash(), &[(0.005, 6.0, 0.95, 0.008)], false),
        "ride" => drum_cards::cymbal(&drum_cards::ride(false), &[(0.005, 2.5, 0.6, 0.004)], false),
        "tom-12 (no table)" => drum_cards::drum(&kit::tom(12.0), 1.0, &[(0.005, 3.5)], false, false),
        "snare (no table)" => drum_cards::drum(&kit::snare(), 1.0, &[(0.005, 4.5)], false, false),
        other => Err(format!("no piece {other}")),
    }
}

fn digest(out: &[Vec<f64>]) -> String {
    let mut h = fs_blake3::DomainHasher::new("frankensim.vkplate.pexfp378.byte-identity");
    for channel in out {
        for v in channel {
            h.update(&v.to_bits().to_le_bytes());
        }
    }
    h.finalize().to_hex()
}

#[test]
fn renders_without_a_table_are_byte_identical_to_the_pin() {
    let mut differ = Vec::new();
    for (name, pin) in PIN {
        let mut p = piece(name).unwrap_or_else(|e| panic!("{name}: {e}"));
        let out = kit::render(&mut p, SECONDS, |_| {});
        let samples: usize = out.iter().map(Vec::len).sum();
        let got = digest(&out);
        eprintln!("byte-identity\t{name}\t{samples} samples\t{got}\t{}", if got == pin { "identical" } else { "DIFFERS" });
        if got != pin {
            differ.push(name);
        }
    }
    assert!(differ.is_empty(), "not byte-identical to the pin: {differ:?}");
}
