//! fs-soa battery (wf9.5): derive-macro end-to-end fixtures (plain,
//! nested, generic, Qty-typed), 128-byte alignment assertions,
//! AoS↔SoA round-trip bitwise equality, iterator equivalence vs a
//! scalar reference, chunked access with masked tails, view/layout
//! descriptors (logged for auditability), and a property battery of
//! random op sequences mirrored against a `Vec` reference.

use std::fmt::Write as _;

use fs_qty::{Length, Time};
use fs_rand::StreamKey;
use fs_soa::{
    RawView, SOA_ALIGN, Soa, chunk_count, chunks_with_tail, chunks_with_tail_mut, leaf_layout,
};

fn escape_json_string(value: &str) -> String {
    let mut escaped = String::with_capacity(value.len());
    for character in value.chars() {
        match character {
            '"' => escaped.push_str("\\\""),
            '\\' => escaped.push_str("\\\\"),
            '\u{0008}' => escaped.push_str("\\b"),
            '\u{000c}' => escaped.push_str("\\f"),
            '\n' => escaped.push_str("\\n"),
            '\r' => escaped.push_str("\\r"),
            '\t' => escaped.push_str("\\t"),
            character if character.is_control() => {
                write!(&mut escaped, "\\u{:04x}", u32::from(character))
                    .expect("writing a JSON escape to String cannot fail");
            }
            character => escaped.push(character),
        }
    }
    escaped
}

fn log_line(case: &str, verdict: &str, detail: &str) -> String {
    format!(
        "{{\"suite\":\"fs-soa\",\"case\":\"{}\",\"verdict\":\"{}\",\"detail\":\"{}\"}}",
        escape_json_string(case),
        escape_json_string(verdict),
        escape_json_string(detail),
    )
}

fn log(case: &str, verdict: &str, detail: &str) {
    println!("{}", log_line(case, verdict, detail));
}

#[derive(Soa, Clone, Copy, Debug, PartialEq)]
struct Particle {
    pos: [f64; 3],
    vel: [f64; 3],
    mass: f64,
    id: u32,
}

#[derive(Soa, Clone, Copy, Debug, PartialEq)]
struct Inner {
    a: f64,
    b: f32,
}

#[derive(Soa, Clone, Copy, Debug, PartialEq)]
struct Outer {
    #[soa(nested)]
    inner: Inner,
    flag: u8,
}

#[derive(Soa, Clone, Copy, Debug, PartialEq)]
struct Generic<T: Copy + std::fmt::Debug + PartialEq, const N: usize> {
    data: [T; N],
    weight: f64,
}

#[derive(Soa, Clone, Copy, Debug, PartialEq)]
struct Dosed {
    dist: Length,
    lag: Time,
}

fn mk_particle(s: &mut fs_rand::Stream) -> Particle {
    Particle {
        pos: [s.next_f64(), s.next_f64(), s.next_f64()],
        vel: [s.next_f64(), s.next_f64(), s.next_f64()],
        mass: s.next_f64() + 0.5,
        id: u32::try_from(s.next_below(1 << 30)).expect("bounded"),
    }
}

#[test]
fn leaf_buffers_are_128_byte_aligned() {
    let mut soa = ParticleSoa::new();
    let mut stream = StreamKey {
        seed: 7,
        kernel: 0x50A0,
        tile: 1,
    }
    .stream();
    for _ in 0..1000 {
        soa.push(mk_particle(&mut stream));
    }
    for v in soa.field_views() {
        assert!(v.addr != 0, "field {} unallocated", v.name);
        assert_eq!(
            v.addr % SOA_ALIGN,
            0,
            "field {} addr misaligned: %128 = {}",
            v.name,
            v.addr % SOA_ALIGN
        );
        assert!(
            v.achieved_align >= SOA_ALIGN,
            "field {}: achieved {}",
            v.name,
            v.achieved_align
        );
        assert_eq!(v.len, 1000);
        assert_eq!(v.stride_bytes, v.elem_bytes, "views are dense");
        log("alignment", "pass", &v.descr());
    }
}

#[test]
fn aos_soa_roundtrip_bitwise() {
    let mut soa = ParticleSoa::with_capacity(64);
    let mut reference = Vec::new();
    let mut stream = StreamKey {
        seed: 11,
        kernel: 0x50A0,
        tile: 2,
    }
    .stream();
    for _ in 0..500 {
        let p = mk_particle(&mut stream);
        soa.push(p);
        reference.push(p);
    }
    assert_eq!(soa.len(), reference.len());
    for (i, r) in reference.iter().enumerate() {
        let g = soa.get(i);
        for k in 0..3 {
            assert_eq!(
                g.pos[k].to_bits(),
                r.pos[k].to_bits(),
                "pos[{k}] bits at {i}"
            );
            assert_eq!(
                g.vel[k].to_bits(),
                r.vel[k].to_bits(),
                "vel[{k}] bits at {i}"
            );
        }
        assert_eq!(g.mass.to_bits(), r.mass.to_bits());
        assert_eq!(g.id, r.id);
    }
    // Iterator equivalence vs the scalar reference.
    assert!(
        soa.iter().zip(&reference).all(|(a, b)| a == *b),
        "iter mismatch"
    );
    // Column access agrees with gathered values.
    for (i, m) in soa.mass().iter().enumerate() {
        assert_eq!(m.to_bits(), reference[i].mass.to_bits());
    }
    log(
        "roundtrip",
        "pass",
        &format!("n={} bitwise", reference.len()),
    );
}

#[test]
fn scatter_set_and_clear_reuse() {
    let mut soa = ParticleSoa::new();
    let mut stream = StreamKey {
        seed: 13,
        kernel: 0x50A0,
        tile: 3,
    }
    .stream();
    for _ in 0..100 {
        soa.push(mk_particle(&mut stream));
    }
    let replacement = mk_particle(&mut stream);
    soa.set(41, replacement);
    assert_eq!(soa.get(41), replacement);
    // Column mutation shows through gather.
    soa.mass_mut()[7] = 42.0;
    assert_eq!(soa.get(7).mass.to_bits(), 42.0f64.to_bits());
    let cap_before = soa.capacity();
    soa.clear();
    assert_eq!(soa.len(), 0);
    assert!(soa.is_empty());
    assert_eq!(soa.capacity(), cap_before, "clear keeps allocations");
    soa.push(replacement);
    assert_eq!(soa.get(0), replacement);
    log("scatter-clear", "pass", &format!("cap_kept={cap_before}"));
}

#[test]
fn capacity_hint_prevents_regrowth() {
    let mut soa = ParticleSoa::with_capacity(256);
    let mut stream = StreamKey {
        seed: 17,
        kernel: 0x50A0,
        tile: 4,
    }
    .stream();
    soa.push(mk_particle(&mut stream));
    let cap = soa.capacity();
    assert!(cap >= 256, "hint not honored: {cap}");
    let addr0 = soa.field_views()[0].addr;
    for _ in 1..256 {
        soa.push(mk_particle(&mut stream));
    }
    assert_eq!(
        soa.field_views()[0].addr,
        addr0,
        "no reallocation within hinted capacity"
    );
    log("capacity-hint", "pass", &format!("cap={cap}"));
}

#[test]
fn nested_containers_compose() {
    let mut soa = OuterSoa::new();
    for i in 0..300u32 {
        soa.push(Outer {
            inner: Inner {
                a: f64::from(i) * 1.5,
                b: f32::from(u16::try_from(i).expect("small")),
            },
            flag: u8::try_from(i % 251).expect("bounded"),
        });
    }
    assert_eq!(soa.len(), 300);
    let got = soa.get(123);
    assert_eq!(got.inner.a.to_bits(), (123.0f64 * 1.5).to_bits());
    assert_eq!(got.flag, 123);
    // Nested accessor drills into the inner container's columns.
    assert_eq!(soa.inner().a().len(), 300);
    assert_eq!(soa.inner().a()[10].to_bits(), 15.0f64.to_bits());
    // Dotted view paths.
    let names: Vec<String> = soa.field_views().into_iter().map(|v| v.name).collect();
    assert_eq!(names, vec!["inner.a", "inner.b", "flag"]);
    log("nested", "pass", &format!("views={names:?}"));
}

#[test]
fn generic_struct_with_bounds_and_const_param() {
    let mut soa: GenericSoa<f32, 4> = GenericSoa::new();
    for i in 0..64i32 {
        let base = i * 4;
        soa.push(Generic {
            data: [
                base as f32,
                (base + 1) as f32,
                (base + 2) as f32,
                (base + 3) as f32,
            ],
            weight: f64::from(i),
        });
    }
    let expect = [40.0f32, 41.0, 42.0, 43.0];
    for (a, b) in soa.get(10).data.iter().zip(&expect) {
        assert_eq!(a.to_bits(), b.to_bits());
    }
    assert_eq!(soa.weight()[63].to_bits(), 63.0f64.to_bits());
    log("generic", "pass", "GenericSoa<f32, 4> roundtrip");
}

#[test]
fn qty_typed_fields_stay_typed() {
    let mut soa = DosedSoa::new();
    for i in 0..50u32 {
        soa.push(Dosed {
            dist: Length::new(f64::from(i) * 0.25),
            lag: Time::new(f64::from(i)),
        });
    }
    // The column keeps the dimensional type — no unit erasure.
    let col: &[Length] = soa.dist();
    assert_eq!(col[8].value().to_bits(), 2.0f64.to_bits());
    assert_eq!(soa.get(8).lag.value().to_bits(), 8.0f64.to_bits());
    // Qty is an 8-byte leaf: full 128-byte alignment applies.
    for v in soa.field_views() {
        assert_eq!(v.elem_bytes, 8);
        assert_eq!(v.addr % SOA_ALIGN, 0, "{} misaligned", v.name);
    }
    log("qty", "pass", "dimensional columns preserved");
}

#[test]
fn chunked_access_with_masked_tail() {
    let mut soa = ParticleSoa::new();
    let mut stream = StreamKey {
        seed: 23,
        kernel: 0x50A0,
        tile: 5,
    }
    .stream();
    for _ in 0..1030 {
        soa.push(mk_particle(&mut stream));
    }
    let (chunks, tail) = chunks_with_tail(soa.mass(), 8);
    assert_eq!(chunks.len(), 128);
    assert_eq!(tail.len(), 6);
    let mut total = 0usize;
    for c in chunks {
        assert_eq!(c.len(), 8);
        total += c.len();
    }
    total += tail.len();
    assert_eq!(
        total, 1030,
        "chunks + tail cover every element exactly once"
    );
    // Mutable chunk pass: scale the mass column, tail included.
    let mut it = chunks_with_tail_mut(soa.mass_mut(), 8);
    for c in &mut it {
        for m in c {
            *m *= 2.0;
        }
    }
    for m in it.into_remainder() {
        *m *= 2.0;
    }
    assert!(
        soa.mass().iter().all(|&m| m >= 1.0),
        "tail elements scaled too"
    );
    // Tile-identity hook: quantum grouping.
    assert_eq!(chunk_count(1030, 512), 3);
    assert_eq!(chunk_count(1024, 512), 2);
    assert_eq!(chunk_count(0, 512), 0);
    log("chunks", "pass", "128 full + 6 tail, quantum groups 3");
}

#[test]
fn layout_and_view_descriptions_are_stable() {
    let descr = ParticleSoa::layout_descr();
    println!("{descr}");
    let lines: Vec<&str> = descr.lines().collect();
    assert_eq!(lines.len(), 4);
    assert_eq!(
        lines[2],
        "{\"field\":\"mass\",\"elem_bytes\":8,\"elem_align\":8,\"dtype\":\"f64\"}"
    );
    assert!(lines[0].contains("\"field\":\"pos\"") && lines[0].contains("\"elem_bytes\":24"));
    assert!(lines[3].contains("\"field\":\"id\"") && lines[3].contains("\"elem_bytes\":4"));
    // Nested layout recursion carries dotted paths, instance-free.
    let outer = OuterSoa::layout_descr();
    println!("{outer}");
    assert!(outer.contains("\"field\":\"inner.a\""));
    assert!(outer.contains("\"field\":\"flag\""));
    // View descr() excludes addresses (deterministic logs).
    let soa = {
        let mut s = ParticleSoa::new();
        s.push(Particle {
            pos: [0.0; 3],
            vel: [0.0; 3],
            mass: 1.0,
            id: 0,
        });
        s
    };
    let d = soa.field_views()[2].descr();
    assert_eq!(
        d,
        "{\"field\":\"mass\",\"len\":1,\"elem_bytes\":8,\"stride_bytes\":8,\"dtype\":\"f64\"}"
    );
    log("layout", "pass", "descriptions address-free and exact");
}

#[test]
fn descriptor_json_escapes_every_dynamic_string() {
    const HOSTILE: &str = "field\"\\\u{0008}\u{000c}\n\r\t\u{0000}\u{001f}";
    const ESCAPED: &str = "field\\\"\\\\\\b\\f\\n\\r\\t\\u0000\\u001f";
    let view = RawView {
        name: HOSTILE.to_owned(),
        addr: 0,
        len: 3,
        elem_bytes: 4,
        stride_bytes: 4,
        achieved_align: SOA_ALIGN,
        dtype: HOSTILE,
    };
    let descr = view.descr();
    assert_eq!(
        descr,
        format!(
            "{{\"field\":\"{ESCAPED}\",\"len\":3,\"elem_bytes\":4,\"stride_bytes\":4,\"dtype\":\"{ESCAPED}\"}}"
        )
    );
    assert_eq!(descr.lines().count(), 1, "{descr}");
    assert!(!descr.chars().any(|ch| ch < ' '), "{descr}");

    let layout = leaf_layout::<u8>(HOSTILE);
    assert_eq!(
        layout,
        format!("{{\"field\":\"{ESCAPED}\",\"elem_bytes\":1,\"elem_align\":1,\"dtype\":\"u8\"}}")
    );
    assert_eq!(layout.lines().count(), 1, "{layout}");
    assert!(!layout.chars().any(|ch| ch < ' '), "{layout}");
}

#[test]
fn battery_log_line_escapes_hostile_nested_descriptor() {
    assert_eq!(
        log_line("roundtrip", "pass", "n=500 bitwise"),
        r#"{"suite":"fs-soa","case":"roundtrip","verdict":"pass","detail":"n=500 bitwise"}"#,
        "ordinary battery output must remain byte-for-byte stable"
    );

    const HOSTILE: &str = "field\"\\\u{0008}\u{000c}\n\r\t\u{0000}\u{001f}";
    let nested_descriptor = RawView {
        name: HOSTILE.to_owned(),
        addr: 0,
        len: 3,
        elem_bytes: 4,
        stride_bytes: 4,
        achieved_align: SOA_ALIGN,
        dtype: HOSTILE,
    }
    .descr();
    let line = log_line(HOSTILE, HOSTILE, &nested_descriptor);

    assert_eq!(
        line,
        r#"{"suite":"fs-soa","case":"field\"\\\b\f\n\r\t\u0000\u001f","verdict":"field\"\\\b\f\n\r\t\u0000\u001f","detail":"{\"field\":\"field\\\"\\\\\\b\\f\\n\\r\\t\\u0000\\u001f\",\"len\":3,\"elem_bytes\":4,\"stride_bytes\":4,\"dtype\":\"field\\\"\\\\\\b\\f\\n\\r\\t\\u0000\\u001f\"}"}"#
    );
    assert_eq!(
        line.lines().count(),
        1,
        "one JSON record must occupy one physical line"
    );
    assert!(
        !line.chars().any(char::is_control),
        "the physical JSON line must contain no raw control characters: {line:?}"
    );
}

#[test]
fn property_random_ops_match_vec_reference() {
    // Mirror a random op sequence against Vec<Particle>; any layout or
    // bookkeeping bug diverges the gather.
    let mut soa = ParticleSoa::new();
    let mut reference: Vec<Particle> = Vec::new();
    let mut stream = StreamKey {
        seed: 31,
        kernel: 0x50A0,
        tile: 6,
    }
    .stream();
    for step in 0..20_000u32 {
        match stream.next_below(10) {
            0..=5 => {
                let p = mk_particle(&mut stream);
                soa.push(p);
                reference.push(p);
            }
            6 | 7 if !reference.is_empty() => {
                let i =
                    usize::try_from(stream.next_below(reference.len() as u64)).expect("index fits");
                let p = mk_particle(&mut stream);
                soa.set(i, p);
                reference[i] = p;
            }
            8 if !reference.is_empty() => {
                let i =
                    usize::try_from(stream.next_below(reference.len() as u64)).expect("index fits");
                assert_eq!(soa.get(i), reference[i], "gather diverged at step {step}");
            }
            9 if step % 4001 == 0 => {
                soa.clear();
                reference.clear();
            }
            _ => {}
        }
        assert_eq!(soa.len(), reference.len());
    }
    assert!(
        soa.iter().zip(&reference).all(|(a, b)| a == *b),
        "final sweep diverged"
    );
    log(
        "property",
        "pass",
        &format!("20000 ops, final n={}", reference.len()),
    );
}
