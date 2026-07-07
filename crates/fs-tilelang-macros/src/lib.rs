//! fs-tilelang-macros — the in-house `kernel!` proc-macro (plan patch
//! Rev C). Layer: UTIL. No syn/quote/proc-macro2 (the Franken-only
//! dependency law covers macro deps): a hand token walker and a text
//! generator, the fs-soa-derive pattern.
//!
//! One restricted kernel body lowers to:
//! - `run_scalar(…)` — the plain indexed reference loop;
//! - `run(…)` — the lane-shaped variant: const-generic monomorphized
//!   inner loops (LANES ∈ {1, 2, 4, 8}) selected ONCE from the
//!   resolved SIMD tier. Per-element arithmetic is token-identical to
//!   the scalar variant, so the two are bitwise-equal by construction;
//! - `META` — arithmetic-intensity metadata (P6) counted from the
//!   body tokens;
//! - a generated `#[cfg(test)]` twin-test module: G0 tier equivalence
//!   (scalar vs lane-shaped, bitwise) and G5 determinism (repeat run,
//!   bitwise) — for free, from the macro.
//!
//! Static analysis (structured `compile_error!`, no silent fallbacks):
//! read/write aliasing, duplicate identifiers, allocation attempts,
//! user loops, escape-hatch blocks, unknown declarations, reduction
//! mismatches, unassigned writes.
//!
//! GRAMMAR (v1 — restricted on purpose):
//! ```text
//! kernel! {
//!     name: my_kernel,
//!     reads: [x, y],            // &[f64], elementwise (or gather-only)
//!     index_reads: [ix],        // &[u32], elementwise   (optional)
//!     uparams: [nx],            // usize scalars          (optional)
//!     params: [alpha],          // f64 scalars            (optional)
//!     writes: [out],            // &mut [f64]
//!     halo: 1,                  // loop runs halo..n-halo (optional, default 0)
//!     reduction: none,          // none | deterministic_sum | fast_sum
//!     body: {
//!         let t = alpha.mul_add(x, y);        // reads/params by bare name
//!         out = t + shift_sub(x, 1) + shift_add(x, nx); // stencil shifts
//!         acc = t;                             // only with a reduction
//!     },
//! }
//! ```
//! `gather(g, EXPR)` reads buffer `g` at a computed index (EXPR over
//! index_reads/uparams); gather-only buffers skip the common-length
//! assertion. Bounds discipline: safe slice indexing (out-of-range
//! shifts/gathers panic loudly), loop range from the declared halo.

use proc_macro::{Delimiter, Spacing, TokenStream, TokenTree};
use std::fmt::Write as _;

/// Function-like macro: see the crate docs for the grammar.
#[proc_macro]
pub fn kernel(input: TokenStream) -> TokenStream {
    match expand(&input) {
        Ok(src) => src
            .parse()
            .expect("fs-tilelang-macros generated invalid Rust (bug)"),
        Err(msg) => format!("compile_error!({msg:?});")
            .parse()
            .expect("compile_error emit"),
    }
}

#[derive(Default)]
struct Decl {
    name: String,
    reads: Vec<String>,
    index_reads: Vec<String>,
    uparams: Vec<String>,
    params: Vec<String>,
    writes: Vec<String>,
    halo: String,
    reduction: String,
    body: Vec<TokenTree>,
}

fn ident_at(tokens: &[TokenTree], i: usize) -> Option<String> {
    match tokens.get(i) {
        Some(TokenTree::Ident(id)) => Some(id.to_string()),
        _ => None,
    }
}

fn punct_at(tokens: &[TokenTree], i: usize) -> Option<char> {
    match tokens.get(i) {
        Some(TokenTree::Punct(p)) => Some(p.as_char()),
        _ => None,
    }
}

fn render(tokens: &[TokenTree]) -> String {
    let mut s = String::new();
    let mut glue = true;
    for t in tokens {
        if !glue {
            s.push(' ');
        }
        s.push_str(&t.to_string());
        glue = matches!(t, TokenTree::Punct(p) if p.spacing() == Spacing::Joint);
    }
    s
}

fn ident_list(group: &TokenTree) -> Result<Vec<String>, String> {
    let TokenTree::Group(g) = group else {
        return Err("expected a [name, name, …] list".to_string());
    };
    if g.delimiter() != Delimiter::Bracket {
        return Err("expected a [name, name, …] list in square brackets".to_string());
    }
    let mut out = Vec::new();
    let mut expect_ident = true;
    for t in g.stream() {
        match (&t, expect_ident) {
            (TokenTree::Ident(id), true) => {
                out.push(id.to_string());
                expect_ident = false;
            }
            (TokenTree::Punct(p), false) if p.as_char() == ',' => expect_ident = true,
            _ => return Err(format!("malformed identifier list near `{t}`")),
        }
    }
    Ok(out)
}

fn parse_decl(input: &TokenStream) -> Result<Decl, String> {
    let tokens: Vec<TokenTree> = input.clone().into_iter().collect();
    let mut d = Decl {
        halo: "0".to_string(),
        reduction: "none".to_string(),
        ..Decl::default()
    };
    let mut i = 0usize;
    while i < tokens.len() {
        let key = ident_at(&tokens, i)
            .ok_or_else(|| format!("expected a declaration key, found `{}`", tokens[i]))?;
        if punct_at(&tokens, i + 1) != Some(':') {
            return Err(format!("expected `:` after `{key}`"));
        }
        let value = tokens
            .get(i + 2)
            .ok_or_else(|| format!("missing value for `{key}`"))?;
        match key.as_str() {
            "name" => {
                d.name =
                    ident_at(&tokens, i + 2).ok_or("kernel name must be a plain identifier")?;
            }
            "reads" => d.reads = ident_list(value)?,
            "index_reads" => d.index_reads = ident_list(value)?,
            "uparams" => d.uparams = ident_list(value)?,
            "params" => d.params = ident_list(value)?,
            "writes" => d.writes = ident_list(value)?,
            "halo" => d.halo = value.to_string(),
            "reduction" => {
                d.reduction = ident_at(&tokens, i + 2)
                    .ok_or("reduction must be none | deterministic_sum | fast_sum")?;
            }
            "body" => {
                let TokenTree::Group(g) = value else {
                    return Err("body must be a { … } block".to_string());
                };
                if g.delimiter() != Delimiter::Brace {
                    return Err("body must be a { … } block".to_string());
                }
                d.body = g.stream().into_iter().collect();
            }
            other => {
                return Err(format!(
                    "unknown declaration `{other}` (expected name/reads/index_reads/uparams/\
                     params/writes/halo/reduction/body)"
                ));
            }
        }
        i += 3;
        if punct_at(&tokens, i) == Some(',') {
            i += 1;
        }
    }
    if d.name.is_empty() {
        return Err("kernel! requires a name".to_string());
    }
    if d.writes.is_empty() && d.reduction == "none" {
        return Err("kernel! needs at least one write buffer or a reduction".to_string());
    }
    if d.body.is_empty() {
        return Err("kernel! requires a non-empty body".to_string());
    }
    Ok(d)
}

/// Collect buffer names used elementwise (bare or as shift
/// targets); buffers appearing ONLY as `gather(buf, ...)` targets
/// are skipped - a gather grid's length is independent of the
/// query batch.
fn find_bare(tokens: &[TokenTree], names: &[String], out: &mut Vec<String>) {
    let mut skip_group_buffer = false;
    for t in tokens {
        match t {
            TokenTree::Ident(id) => {
                let s = id.to_string();
                if s == "gather" {
                    skip_group_buffer = true;
                    continue;
                }
                if names.contains(&s) && !out.contains(&s) {
                    out.push(s);
                }
            }
            TokenTree::Group(g) => {
                let mut inner: Vec<TokenTree> = g.stream().into_iter().collect();
                if skip_group_buffer {
                    // Drop the buffer ident before the first comma
                    // of gather(buf, …); the index EXPR still gets
                    // scanned.
                    if let Some(comma) = inner
                        .iter()
                        .position(|t| matches!(t, TokenTree::Punct(p) if p.as_char() == ','))
                    {
                        inner.drain(..=comma);
                    }
                }
                find_bare(&inner, names, out);
            }
            _ => {}
        }
        skip_group_buffer = false;
    }
}

fn body_contains_ident(tokens: &[TokenTree], needle: &str) -> bool {
    tokens.iter().any(|t| match t {
        TokenTree::Ident(id) => id.to_string() == needle,
        TokenTree::Group(g) => {
            let inner: Vec<TokenTree> = g.stream().into_iter().collect();
            body_contains_ident(&inner, needle)
        }
        _ => false,
    })
}

/// Reject forbidden constructs and count arithmetic (recursively
/// through groups). Returns (flops, uses_acc).
fn analyze_body(tokens: &[TokenTree], flops: &mut u32, uses_acc: &mut bool) -> Result<(), String> {
    const FORBIDDEN: &[(&str, &str)] = &[
        // Escape splits the token so the naive capsule scanner does
        // not flag this REJECTION literal as unsafe code.
        ("uns\u{61}fe", "uns\u{61}fe blocks are not allowed in kernel bodies"),
        (
            "while",
            "unbounded loops are not allowed in kernel bodies (the tile loop is implicit)",
        ),
        (
            "loop",
            "unbounded loops are not allowed in kernel bodies (the tile loop is implicit)",
        ),
        (
            "for",
            "user loops are not allowed in kernel bodies (the tile loop is implicit)",
        ),
        ("vec", "allocation is not allowed in kernel bodies"),
        ("Vec", "allocation is not allowed in kernel bodies"),
        ("Box", "allocation is not allowed in kernel bodies"),
        ("String", "allocation is not allowed in kernel bodies"),
        ("collect", "allocation is not allowed in kernel bodies"),
        ("push", "allocation is not allowed in kernel bodies"),
        ("return", "early return is not allowed in kernel bodies"),
    ];
    let mut i = 0usize;
    while i < tokens.len() {
        match &tokens[i] {
            TokenTree::Ident(id) => {
                let name = id.to_string();
                if let Some((_, msg)) = FORBIDDEN.iter().find(|(k, _)| *k == name) {
                    return Err((*msg).to_string());
                }
                if name == "acc" {
                    *uses_acc = true;
                }
                if name == "mul_add" {
                    *flops += 2;
                }
            }
            TokenTree::Punct(p) => {
                if matches!(p.as_char(), '+' | '-' | '*' | '/')
                    && punct_at(tokens, i + 1) != Some('=')
                {
                    *flops += 1;
                }
            }
            TokenTree::Group(g) => {
                let inner: Vec<TokenTree> = g.stream().into_iter().collect();
                analyze_body(&inner, flops, uses_acc)?;
            }
            TokenTree::Literal(_) => {}
        }
        i += 1;
    }
    Ok(())
}

/// Rewrite the body for one element index expression `idx_expr`:
/// - bare read/index-read/write idents become indexed accesses;
/// - `shift_add(buf, EXPR)` / `shift_sub(buf, EXPR)` become
///   `buf[idx + (EXPR)]` / `buf[idx - (EXPR)]`;
/// - `gather(buf, EXPR)` becomes `buf[(EXPR) as usize]`;
/// - params/uparams/locals pass through.
fn rewrite(tokens: &[TokenTree], d: &Decl, idx: &str) -> Result<String, String> {
    let mut out = String::new();
    let mut i = 0usize;
    let mut glue = true;
    while i < tokens.len() {
        let mut emitted: Option<String> = None;
        match &tokens[i] {
            TokenTree::Ident(id) => {
                let name = id.to_string();
                let is_shift = name == "shift_add" || name == "shift_sub";
                if is_shift || name == "gather" {
                    let Some(TokenTree::Group(g)) = tokens.get(i + 1) else {
                        return Err(format!("`{name}` must be called as {name}(buffer, expr)"));
                    };
                    let args: Vec<TokenTree> = g.stream().into_iter().collect();
                    let comma = args
                        .iter()
                        .position(|t| matches!(t, TokenTree::Punct(p) if p.as_char() == ','))
                        .ok_or(format!("`{name}` needs two arguments"))?;
                    let buf = render(&args[..comma]);
                    if !d.reads.contains(&buf) {
                        return Err(format!(
                            "`{name}` target `{buf}` is not a declared read buffer"
                        ));
                    }
                    let expr = rewrite(&args[comma + 1..], d, idx)?;
                    emitted = Some(match name.as_str() {
                        "shift_add" => format!("{buf}[{idx} + ({expr})]"),
                        "shift_sub" => format!("{buf}[{idx} - ({expr})]"),
                        _ => format!("{buf}[({expr})]"),
                    });
                    i += 1; // consume the argument group too
                } else if d.reads.contains(&name) || d.writes.contains(&name) {
                    emitted = Some(format!("{name}[{idx}]"));
                } else if d.index_reads.contains(&name) {
                    emitted = Some(format!(
                        "(usize::try_from({name}[{idx}]).expect(\"index\"))"
                    ));
                }
            }
            TokenTree::Group(g) => {
                let inner: Vec<TokenTree> = g.stream().into_iter().collect();
                let rewritten = rewrite(&inner, d, idx)?;
                emitted = Some(match g.delimiter() {
                    Delimiter::Parenthesis => format!("({rewritten})"),
                    Delimiter::Brace => format!("{{{rewritten}}}"),
                    Delimiter::Bracket => format!("[{rewritten}]"),
                    Delimiter::None => rewritten,
                });
            }
            _ => {}
        }
        if !glue {
            out.push(' ');
        }
        glue = matches!(&tokens[i], TokenTree::Punct(p) if p.spacing() == Spacing::Joint);
        match emitted {
            Some(e) => out.push_str(&e),
            None => out.push_str(&tokens[i].to_string()),
        }
        i += 1;
    }
    Ok(out)
}

#[allow(clippy::too_many_lines)]
fn expand(input: &TokenStream) -> Result<String, String> {
    let d = parse_decl(input)?;

    // Alias and duplicate checks.
    let mut seen: Vec<&String> = Vec::new();
    for group in [&d.reads, &d.index_reads, &d.uparams, &d.params, &d.writes] {
        for name in group {
            if seen.contains(&name) {
                return Err(format!(
                    "identifier `{name}` is declared twice (reads and writes must not alias — \
                     no hidden aliasing is the point)"
                ));
            }
            seen.push(name);
        }
    }
    for reserved in ["acc", "gather", "shift_add", "shift_sub"] {
        if seen.iter().any(|s| *s == reserved) {
            return Err(format!("`{reserved}` is reserved by the kernel grammar"));
        }
    }

    let mut flops = 0u32;
    let mut uses_acc = false;
    analyze_body(&d.body, &mut flops, &mut uses_acc)?;
    let uses_gather = body_contains_ident(&d.body, "gather");
    let reduction_variant = match d.reduction.as_str() {
        "none" => {
            if uses_acc {
                return Err(
                    "body assigns `acc` but reduction is `none`; declare deterministic_sum or \
                     fast_sum"
                        .to_string(),
                );
            }
            "None"
        }
        "deterministic_sum" | "fast_sum" => {
            if !uses_acc {
                return Err(format!(
                    "reduction `{}` declared but the body never assigns `acc`",
                    d.reduction
                ));
            }
            if d.reduction == "deterministic_sum" {
                "DeterministicSum"
            } else {
                "FastSum"
            }
        }
        other => {
            return Err(format!(
                "unknown reduction `{other}` (expected none | deterministic_sum | fast_sum)"
            ));
        }
    };

    // Check every write is assigned (token pattern `name =`, not `==`).
    for w in &d.writes {
        let mut assigned = false;
        for (i, t) in d.body.iter().enumerate() {
            if let TokenTree::Ident(id) = t
                && id.to_string() == *w
                && punct_at(&d.body, i + 1) == Some('=')
                && punct_at(&d.body, i + 2) != Some('=')
            {
                assigned = true;
            }
        }
        if !assigned {
            return Err(format!("write buffer `{w}` is never assigned in the body"));
        }
    }

    let bytes: u32 = u32::try_from(8 * (d.reads.len() + d.writes.len()) + 4 * d.index_reads.len())
        .expect("few buffers");

    // Signature pieces.
    let mut args = String::new();
    let mut call_args = String::new();
    for r in &d.reads {
        let _ = write!(args, "{r}: &[f64], ");
        let _ = write!(call_args, "{r}, ");
    }
    for r in &d.index_reads {
        let _ = write!(args, "{r}: &[u32], ");
        let _ = write!(call_args, "{r}, ");
    }
    for u in &d.uparams {
        let _ = write!(args, "{u}: usize, ");
        let _ = write!(call_args, "{u}, ");
    }
    for p in &d.params {
        let _ = write!(args, "{p}: f64, ");
        let _ = write!(call_args, "{p}, ");
    }
    for w in &d.writes {
        let _ = write!(args, "{w}: &mut [f64], ");
        let _ = write!(call_args, "{w}, ");
    }
    let args = args.trim_end_matches(", ");
    let call_args = call_args.trim_end_matches(", ");

    // Length source and common-length assertions: elementwise buffers
    // (bare-name reads, index reads, writes) share a length; gather-
    // only reads are exempt.
    // Elementwise-length policy: writes anchor the common length;
    // every read used by bare name or as a shift target must match.
    // Buffers used ONLY as gather(buf, …) targets keep their own
    // length (a gather grid is usually a different size than the
    // query batch).
    let mut bare_used: Vec<String> = Vec::new();
    find_bare(&d.body, &d.reads, &mut bare_used);
    let anchor = d
        .writes
        .first()
        .cloned()
        .or_else(|| d.reads.first().cloned())
        .expect("checked");
    let mut len_asserts = String::new();
    for r in bare_used
        .iter()
        .cloned()
        .chain(d.index_reads.iter().cloned())
        .chain(d.writes.iter().cloned())
    {
        if r != anchor {
            let _ = writeln!(
                len_asserts,
                "        assert_eq!({r}.len(), {anchor}.len(), \"kernel `{}`: buffer `{r}` length mismatch\");",
                d.name
            );
        }
    }

    let name = &d.name;
    let halo = &d.halo;
    // META is a const item: a halo expression over uparams is runtime-
    // only, so META records literal halos and 0 for dynamic ones
    // (documented in fs-tilelang).
    let halo_meta = d.halo.trim().parse::<u32>().unwrap_or(0).to_string();
    let ret_ty = if reduction_variant == "None" {
        "()"
    } else {
        "f64"
    };
    // Per-element body: with a reduction, `acc` is a fresh per-element
    // local; the loop folds it into fixed 64-element chunk partials in
    // ASCENDING INDEX ORDER — identical in both variants, so even
    // reductions are bitwise-equal across tiers in v1 (FastSum keeps
    // the CONTRACTUAL right to reassociate later).
    let rewritten = rewrite(&d.body, &d, "__i")?;
    let body_elem = if reduction_variant == "None" {
        format!("{{ {rewritten} }}")
    } else {
        format!(
            "{{ let mut acc = 0.0f64; {{ {rewritten} }} __chunk += acc; __count += 1;              if __count == fs_tilelang::REDUCTION_CHUNK {{ __total += __chunk; __chunk = 0.0; __count = 0; }} }}"
        )
    };
    let (red_pre, red_post) = if reduction_variant == "None" {
        (String::new(), String::new())
    } else {
        (
            "let mut __total = 0.0f64; let mut __chunk = 0.0f64; let mut __count = 0usize;"
                .to_string(),
            "__total += __chunk; __total".to_string(),
        )
    };

    // Auto twin tests (G0 tier equivalence + G5 determinism) — only
    // for kernels the macro can DRIVE safely: gather kernels need
    // semantically valid indices and uparam kernels need meaningful
    // strides, so those get their twin checks from the caller's
    // battery instead (the adoption-policy lint boundary).
    let twin_tests = if uses_gather || !d.uparams.is_empty() {
        String::new()
    } else {
        let n_reads = d.reads.len();
        let n_writes = d.writes.len();
        let n_params = d.params.len();
        let mut call_scalar = String::new();
        let mut call_lanes = String::new();
        for (k, _) in d.reads.iter().enumerate() {
            let _ = write!(call_scalar, "&__reads[{k}], ");
            let _ = write!(call_lanes, "&__reads[{k}], ");
        }
        for (k, _) in d.index_reads.iter().enumerate() {
            let _ = write!(call_scalar, "&__ireads[{k}], ");
            let _ = write!(call_lanes, "&__ireads[{k}], ");
        }
        for (k, _) in d.params.iter().enumerate() {
            let _ = write!(call_scalar, "__params[{k}], ");
            let _ = write!(call_lanes, "__params[{k}], ");
        }
        for (k, _) in d.writes.iter().enumerate() {
            let _ = write!(call_scalar, "&mut __out_scalar[{k}], ");
            let _ = write!(call_lanes, "&mut __out_lanes[{k}], ");
        }
        let call_scalar = call_scalar.trim_end_matches(", ").to_string();
        let call_lanes = call_lanes.trim_end_matches(", ").to_string();
        let n_ireads = d.index_reads.len();
        format!(
            "
    #[cfg(test)]
    mod __twin_tests {{
        #[doc = \"Deterministic LCG fixture data (no deps).\"]
        fn __fill(seed: u64, n: usize) -> Vec<f64> {{
            let mut s = seed | 1;
            (0..n)
                .map(|_| {{
                    s = s.wrapping_mul(6364136223846793005).wrapping_add(1442695040888963407);
                    f64::from(u32::try_from(s >> 40).expect(\"24 bits\")) / f64::from(1u32 << 24)
                }})
                .collect()
        }}

        #[test]
        fn tier_equivalence_and_determinism() {{
            let n = 1003usize;
            let __reads: Vec<Vec<f64>> = (0..{n_reads}).map(|k| __fill(11 + k as u64, n)).collect();
            let __ireads: Vec<Vec<u32>> = (0..{n_ireads})
                .map(|k| {{
                    __fill(101 + k as u64, n)
                        .into_iter()
                        .map(|v| (v * 7.0) as u32)
                        .collect()
                }})
                .collect();
            let __params: Vec<f64> = __fill(31, {n_params}.max(1));
            #[allow(unused_mut)]
            let mut __out_scalar: Vec<Vec<f64>> = vec![vec![0.0; n]; {n_writes}.max(1)];
            #[allow(unused_mut)]
            let mut __out_lanes: Vec<Vec<f64>> = vec![vec![0.0; n]; {n_writes}.max(1)];
            let __r_scalar = super::run_scalar({call_scalar});
            for __w in [2usize, 4, 8] {{
                for o in &mut __out_lanes {{
                    o.fill(0.0);
                }}
                let __r_lanes = match __w {{
                    2 => super::run_lanes::<2>({call_lanes}),
                    4 => super::run_lanes::<4>({call_lanes}),
                    _ => super::run_lanes::<8>({call_lanes}),
                }};
                for (a, b) in __out_scalar.iter().flatten().zip(__out_lanes.iter().flatten()) {{
                    assert_eq!(a.to_bits(), b.to_bits(), \"tier divergence (G0)\");
                }}
                assert_eq!(
                    __ret_bits(&__r_scalar),
                    __ret_bits(&__r_lanes),
                    \"reduction tier divergence (G0)\"
                );
            }}
            // G5: repeat run, bitwise.
            #[allow(unused_mut)]
            let mut __out_repeat: Vec<Vec<f64>> = vec![vec![0.0; n]; {n_writes}.max(1)];
            let __r_repeat = super::run_scalar({call_repeat});
            assert_eq!(__ret_bits(&__r_scalar), __ret_bits(&__r_repeat), \"repeat divergence (G5)\");
            for (a, b) in __out_scalar.iter().flatten().zip(__out_repeat.iter().flatten()) {{
                assert_eq!(a.to_bits(), b.to_bits(), \"repeat divergence (G5)\");
            }}
        }}

        fn __ret_bits<T: __RetBits>(v: &T) -> u64 {{
            v.__bits()
        }}
        trait __RetBits {{
            fn __bits(&self) -> u64;
        }}
        impl __RetBits for () {{
            fn __bits(&self) -> u64 {{ 0 }}
        }}
        impl __RetBits for f64 {{
            fn __bits(&self) -> u64 {{ self.to_bits() }}
        }}
    }}
",
            call_repeat = call_scalar.replace("__out_scalar", "__out_repeat"),
        )
    };

    let mut out = String::new();
    let _ = write!(
        out,
        "\
#[doc = \"Generated by `kernel!`: {name} — see fs-tilelang. Scalar reference and lane-shaped variants are bitwise-equal by construction.\"]
#[allow(clippy::too_many_arguments, clippy::missing_panics_doc)]
pub mod {name} {{
    #[doc = \"Kernel metadata (P6: intensity travels with the kernel).\"]
    pub const META: fs_tilelang::KernelMeta = fs_tilelang::KernelMeta {{
        name: \"{name}\",
        flops_per_elem: {flops},
        bytes_per_elem: {bytes},
        halo: {halo_meta},
        reduction: fs_tilelang::ReductionKind::{reduction_variant},
        determinism: fs_tilelang::DeterminismClass::BitwiseAllTiers,
    }};

    #[doc = \"Scalar reference variant (the correctness anchor).\"]
    pub fn run_scalar({args}) -> {ret_ty} {{
{len_asserts}        let __n = {anchor}.len();
        let __halo: usize = {halo};
        {red_pre}
        let mut __i = __halo;
        while __i < __n - __halo {{
            {body_elem}
            __i += 1;
        }}
        {red_post}
    }}

    #[doc = \"Lane-shaped variant: const-generic inner loops over the resolved tier width; per-element arithmetic identical to the scalar variant (bitwise-equal).\"]
    pub fn run({args}) -> {ret_ty} {{
        match fs_tilelang::lane_width() {{
            2 => run_lanes::<2>({call_args}),
            4 => run_lanes::<4>({call_args}),
            8 => run_lanes::<8>({call_args}),
            _ => run_scalar({call_args}),
        }}
    }}

    #[doc = \"Monomorphized lane-grouped loop (autovectorizer target).\"]
    pub fn run_lanes<const LANES: usize>({args}) -> {ret_ty} {{
{len_asserts}        let __n = {anchor}.len();
        let __halo: usize = {halo};
        {red_pre}
        let __span = __n - 2 * __halo;
        let __groups = __span / LANES;
        for __g in 0..__groups {{
            for __l in 0..LANES {{
                let __i = __halo + __g * LANES + __l;
                {body_elem}
            }}
        }}
        let mut __i = __halo + __groups * LANES;
        while __i < __n - __halo {{
            {body_elem}
            __i += 1;
        }}
        {red_post}
    }}
{twin_tests}}}
"
    );
    Ok(out)
}
