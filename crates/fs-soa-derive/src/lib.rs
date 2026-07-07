//! fs-soa-derive — the in-house `#[derive(Soa)]` proc-macro (plan
//! §5.3). Layer: UTIL.
//!
//! In-house because the Franken-only dependency law covers macro deps
//! too: no syn/quote/proc-macro2 — the parser walks std `proc_macro`
//! tokens by hand and the generator emits source text re-parsed into a
//! `TokenStream`. Deterministic: output depends only on input tokens.
//!
//! For `struct P { a: f64, #[soa(nested)] inner: Q, … }` it generates
//! `PSoa` (one aligned `fs_soa::FieldBuf` per leaf field, nested
//! containers for `#[soa(nested)]` fields), AoS gather/scatter,
//! per-field slice accessors, zip-style value iteration, view
//! descriptors, a layout description, and `SoaAble`/`SoaContainer`
//! impls so containers compose.
//!
//! Diagnostics are compile errors via `compile_error!` — no silent
//! fallbacks (unsupported shapes: tuple/unit structs, enums, unions,
//! lifetimes, generic parameter defaults, zero fields, field names
//! colliding with the generated API).

use proc_macro::{Delimiter, Spacing, TokenStream, TokenTree};
use std::fmt::Write as _;

/// Derive a structure-of-arrays container (`<Name>Soa`) for a named-
/// field struct. Field attribute: `#[soa(nested)]` stores that field
/// in its own generated container (the field type must also derive
/// `Soa`).
#[proc_macro_derive(Soa, attributes(soa))]
pub fn derive_soa(input: TokenStream) -> TokenStream {
    match expand(&input) {
        Ok(src) => src
            .parse()
            .expect("fs-soa-derive generated invalid Rust (bug)"),
        Err(msg) => format!("compile_error!({msg:?});")
            .parse()
            .expect("compile_error emit"),
    }
}

struct Field {
    name: String,
    ty: String,
    nested: bool,
}

/// Method names the generated container claims; struct fields must not
/// collide (each field also becomes an accessor method).
const RESERVED: &[&str] = &[
    "new",
    "with_capacity",
    "len",
    "is_empty",
    "capacity",
    "clear",
    "reserve",
    "push",
    "get",
    "set",
    "iter",
    "field_views",
    "layout_descr",
];

fn expand(input: &TokenStream) -> Result<String, String> {
    let tokens: Vec<TokenTree> = input.clone().into_iter().collect();
    let mut i = 0usize;

    skip_attrs(&tokens, &mut i)?;
    let vis = parse_vis(&tokens, &mut i);

    match ident_at(&tokens, i).as_deref() {
        Some("struct") => i += 1,
        Some(k @ ("enum" | "union")) => {
            return Err(format!(
                "#[derive(Soa)] supports only structs with named fields, not {k}s"
            ));
        }
        _ => return Err("#[derive(Soa)] expects a struct definition".to_string()),
    }

    let name = ident_at(&tokens, i).ok_or("#[derive(Soa)] expects a struct name")?;
    i += 1;

    let (generics, ty_args) = parse_generics(&tokens, &mut i)?;

    // Where clause: everything from `where` up to the body group.
    let mut where_clause = String::new();
    if ident_at(&tokens, i).as_deref() == Some("where") {
        let start = i;
        while i < tokens.len()
            && !matches!(&tokens[i], TokenTree::Group(g) if g.delimiter() == Delimiter::Brace)
        {
            i += 1;
        }
        where_clause = render(&tokens[start..i]);
    }

    let body = match tokens.get(i) {
        Some(TokenTree::Group(g)) if g.delimiter() == Delimiter::Brace => g.stream(),
        Some(TokenTree::Group(g)) if g.delimiter() == Delimiter::Parenthesis => {
            return Err(
                "#[derive(Soa)] supports only named fields; tuple structs have no field names \
                 to become accessors"
                    .to_string(),
            );
        }
        _ => return Err(
            "#[derive(Soa)] supports only structs with named fields (unit structs have no fields)"
                .to_string(),
        ),
    };

    let fields = parse_fields(&body)?;
    if fields.is_empty() {
        return Err("#[derive(Soa)] requires at least one field".to_string());
    }
    for f in &fields {
        if RESERVED.contains(&f.name.as_str()) {
            return Err(format!(
                "field `{}` collides with the generated container API; rename it (reserved: {})",
                f.name,
                RESERVED.join(", ")
            ));
        }
    }

    Ok(generate(
        &vis,
        &name,
        &generics,
        &ty_args,
        &where_clause,
        &fields,
    ))
}

// ------------------------------------------------------------- token walking

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

fn skip_attrs(tokens: &[TokenTree], i: &mut usize) -> Result<(), String> {
    while punct_at(tokens, *i) == Some('#') {
        *i += 1;
        match tokens.get(*i) {
            Some(TokenTree::Group(g)) if g.delimiter() == Delimiter::Bracket => *i += 1,
            _ => return Err("malformed attribute".to_string()),
        }
    }
    Ok(())
}

fn parse_vis(tokens: &[TokenTree], i: &mut usize) -> String {
    if ident_at(tokens, *i).as_deref() == Some("pub") {
        *i += 1;
        if let Some(TokenTree::Group(g)) = tokens.get(*i)
            && g.delimiter() == Delimiter::Parenthesis
        {
            *i += 1;
            return format!("pub({})", g.stream());
        }
        return "pub".to_string();
    }
    String::new()
}

/// Parse `<…>` generics after the struct name. Returns the generics
/// with bounds (verbatim, for `impl<…>` and the type definition) and
/// the bare argument list (`<T, N>`) for type positions.
fn parse_generics(tokens: &[TokenTree], i: &mut usize) -> Result<(String, String), String> {
    if punct_at(tokens, *i) != Some('<') {
        return Ok((String::new(), String::new()));
    }
    let start = *i;
    let mut depth = 0i32;
    let mut prev_joint_minus = false;
    while *i < tokens.len() {
        if let TokenTree::Punct(p) = &tokens[*i] {
            let c = p.as_char();
            if c == '<' && !prev_joint_minus {
                depth += 1;
            } else if c == '>' && !prev_joint_minus {
                depth -= 1;
                if depth == 0 {
                    *i += 1;
                    break;
                }
            }
            prev_joint_minus = c == '-' && p.spacing() == Spacing::Joint;
        } else {
            prev_joint_minus = false;
        }
        *i += 1;
    }
    if depth != 0 {
        return Err("unbalanced generics".to_string());
    }
    let inner = &tokens[start + 1..*i - 1];
    let mut names: Vec<String> = Vec::new();
    for seg in split_top_level(inner) {
        if seg.is_empty() {
            continue;
        }
        if matches!(&seg[0], TokenTree::Punct(p) if p.as_char() == '\'') {
            return Err(
                "#[derive(Soa)] does not support lifetime parameters: SoA containers own their \
                 storage (plain-old-data only)"
                    .to_string(),
            );
        }
        if seg.iter().any(
            |t| matches!(t, TokenTree::Punct(p) if p.as_char() == '=' && p.spacing() == Spacing::Alone),
        ) {
            return Err(
                "#[derive(Soa)] does not support generic parameter defaults; remove the `= …` \
                 default"
                    .to_string(),
            );
        }
        let mut j = 0usize;
        if ident_at(&seg, j).as_deref() == Some("const") {
            j += 1;
        }
        match ident_at(&seg, j) {
            Some(n) => names.push(n),
            None => return Err("unsupported generic parameter shape".to_string()),
        }
    }
    let generics = render(&tokens[start..*i]);
    let ty_args = format!("<{}>", names.join(", "));
    Ok((generics, ty_args))
}

/// Split a token run on top-level commas (`<>` depth 0; groups are
/// already atomic trees).
fn split_top_level(tokens: &[TokenTree]) -> Vec<Vec<TokenTree>> {
    let mut out: Vec<Vec<TokenTree>> = vec![Vec::new()];
    let mut depth = 0i32;
    let mut prev_joint_minus = false;
    for t in tokens {
        if let TokenTree::Punct(p) = t {
            let c = p.as_char();
            if c == '<' && !prev_joint_minus {
                depth += 1;
            } else if c == '>' && !prev_joint_minus {
                depth -= 1;
            } else if c == ',' && depth == 0 {
                out.push(Vec::new());
                prev_joint_minus = false;
                continue;
            }
            prev_joint_minus = c == '-' && p.spacing() == Spacing::Joint;
        } else {
            prev_joint_minus = false;
        }
        out.last_mut().expect("non-empty").push(t.clone());
    }
    out
}

fn parse_fields(body: &TokenStream) -> Result<Vec<Field>, String> {
    let tokens: Vec<TokenTree> = body.clone().into_iter().collect();
    let mut fields = Vec::new();
    let mut i = 0usize;
    while i < tokens.len() {
        // Field attributes: recognize #[soa(…)], skip the rest (docs…).
        let mut nested = false;
        while punct_at(&tokens, i) == Some('#') {
            let Some(TokenTree::Group(g)) = tokens.get(i + 1) else {
                return Err("malformed field attribute".to_string());
            };
            if g.delimiter() != Delimiter::Bracket {
                return Err("malformed field attribute".to_string());
            }
            let attr: Vec<TokenTree> = g.stream().into_iter().collect();
            if ident_at(&attr, 0).as_deref() == Some("soa") {
                let Some(TokenTree::Group(args)) = attr.get(1) else {
                    return Err("expected #[soa(nested)]".to_string());
                };
                let args: Vec<TokenTree> = args.stream().into_iter().collect();
                if args.len() == 1 && ident_at(&args, 0).as_deref() == Some("nested") {
                    nested = true;
                } else {
                    return Err(format!(
                        "unknown #[soa(…)] argument `{}`; supported: nested",
                        render(&args)
                    ));
                }
            }
            i += 2;
        }
        let _field_vis = parse_vis(&tokens, &mut i);
        let Some(name) = ident_at(&tokens, i) else {
            if i >= tokens.len() {
                break;
            }
            return Err("expected a field name".to_string());
        };
        i += 1;
        if punct_at(&tokens, i) != Some(':') {
            return Err(format!("expected `:` after field `{name}`"));
        }
        i += 1;
        // Type: tokens until a top-level comma.
        let start = i;
        let mut depth = 0i32;
        let mut prev_joint_minus = false;
        while i < tokens.len() {
            if let TokenTree::Punct(p) = &tokens[i] {
                let c = p.as_char();
                if c == '<' && !prev_joint_minus {
                    depth += 1;
                } else if c == '>' && !prev_joint_minus {
                    depth -= 1;
                } else if c == ',' && depth == 0 {
                    break;
                }
                prev_joint_minus = c == '-' && p.spacing() == Spacing::Joint;
            } else {
                prev_joint_minus = false;
            }
            i += 1;
        }
        if start == i {
            return Err(format!("field `{name}` has an empty type"));
        }
        let ty = render(&tokens[start..i]);
        i += 1; // past the comma (or the end)
        fields.push(Field { name, ty, nested });
    }
    Ok(fields)
}

fn render(tokens: &[TokenTree]) -> String {
    // Respect Punct spacing: a Joint punct glues to the next token
    // (`::`, `->`); putting a space inside them produces unparsable
    // output (`: :`).
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

// --------------------------------------------------------------- generation

#[allow(clippy::too_many_lines)]
fn generate(
    vis: &str,
    name: &str,
    generics: &str,
    ty_args: &str,
    where_clause: &str,
    fields: &[Field],
) -> String {
    let soa = format!("{name}Soa");
    let value_ty = format!("{name} {ty_args}");
    let soa_ty = format!("{soa} {ty_args}");
    let storage = |f: &Field| {
        if f.nested {
            format!("< {} as ::fs_soa::SoaAble >::Soa", f.ty)
        } else {
            format!("::fs_soa::FieldBuf<{}>", f.ty)
        }
    };

    let mut struct_fields = String::new();
    let mut new_fields = String::new();
    let mut cap_fields = String::new();
    let mut clear_stmts = String::new();
    let mut reserve_stmts = String::new();
    let mut push_stmts = String::new();
    let mut get_fields = String::new();
    let mut set_stmts = String::new();
    let mut capacity_expr = String::from("usize::MAX");
    let mut views_stmts = String::new();
    let mut layout_stmts = String::new();
    let mut accessors = String::new();
    let names: Vec<&str> = fields.iter().map(|f| f.name.as_str()).collect();
    let destructure = names.join(", ");

    for f in fields {
        let (n, st) = (&f.name, storage(f));
        let _ = writeln!(struct_fields, "    {n}: {st},");
        if f.nested {
            // Nested storage is driven through the trait in UFCS form:
            // no `use fs_soa::SoaContainer` demanded of the caller, no
            // inference ambiguity.
            let sc = format!("< {st} as ::fs_soa::SoaContainer<{}> >", f.ty);
            let _ = writeln!(new_fields, "            {n}: {sc}::c_new(),");
            let _ = writeln!(cap_fields, "            {n}: {sc}::c_with_capacity(cap),");
            let _ = writeln!(clear_stmts, "        {sc}::c_clear(&mut self.{n});");
            let _ = writeln!(
                reserve_stmts,
                "        {sc}::c_reserve(&mut self.{n}, additional);"
            );
            let _ = writeln!(push_stmts, "        {sc}::c_push(&mut self.{n}, {n});");
            let _ = writeln!(get_fields, "            {n}: {sc}::c_get(&self.{n}, i),");
            let _ = writeln!(set_stmts, "        {sc}::c_set(&mut self.{n}, i, {n});");
            let _ = writeln!(
                views_stmts,
                "        {sc}::c_views(&self.{n}, &::fs_soa::view_name(prefix, \"{n}\"), out);"
            );
            let _ = writeln!(
                layout_stmts,
                "        {sc}::c_layout(&::fs_soa::view_name(prefix, \"{n}\"), out);"
            );
            let _ = writeln!(
                accessors,
                "    /// Nested SoA container for field `{n}`.\n    #[must_use]\n    \
                 {vis} fn {n}(&self) -> &{st} {{ &self.{n} }}\n\
                     /// Mutable nested SoA container for field `{n}`.\n    \
                 {vis} fn {n}_mut(&mut self) -> &mut {st} {{ &mut self.{n} }}"
            );
        } else {
            let _ = writeln!(new_fields, "            {n}: ::fs_soa::FieldBuf::new(),");
            let _ = writeln!(
                cap_fields,
                "            {n}: ::fs_soa::FieldBuf::with_capacity(cap),"
            );
            let _ = writeln!(clear_stmts, "        self.{n}.clear();");
            let _ = writeln!(reserve_stmts, "        self.{n}.reserve(additional);");
            let _ = writeln!(push_stmts, "        self.{n}.push({n});");
            let _ = writeln!(get_fields, "            {n}: self.{n}.as_slice()[i],");
            let _ = writeln!(set_stmts, "        self.{n}.as_mut_slice()[i] = {n};");
            capacity_expr = format!("{capacity_expr}.min(self.{n}.capacity())");
            let _ = writeln!(
                views_stmts,
                "        out.push(self.{n}.view(&::fs_soa::view_name(prefix, \"{n}\")));"
            );
            let _ = writeln!(
                layout_stmts,
                "        out.push(::fs_soa::leaf_layout::<{}>(&::fs_soa::view_name(prefix, \
                 \"{n}\")));",
                f.ty
            );
            let _ = writeln!(
                accessors,
                "    /// Field `{n}` as a dense aligned slice.\n    #[must_use]\n    \
                 {vis} fn {n}(&self) -> &[{ty}] {{ self.{n}.as_slice() }}\n\
                     /// Field `{n}` as a mutable dense aligned slice.\n    \
                 {vis} fn {n}_mut(&mut self) -> &mut [{ty}] {{ self.{n}.as_mut_slice() }}",
                ty = f.ty
            );
        }
    }
    // Nested containers track their own length; capacity() reports the
    // leaf minimum (usize::MAX for an all-nested container is fine:
    // capacity is a leaf concept and such containers delegate).

    format!(
        "\
#[doc = \"Structure-of-arrays container for [`{name}`], generated by `#[derive(Soa)]`: one \
128-byte-aligned buffer per leaf field, SIMD lanes run across elements.\"]
#[derive(Debug, Clone)]
{vis} struct {soa} {generics} {where_clause} {{
{struct_fields}    len: usize,
}}

#[automatically_derived]
impl {generics} {soa_ty} {where_clause} {{
    #[doc = \"Empty container.\"]
    #[must_use]
    {vis} fn new() -> Self {{
        Self {{
{new_fields}            len: 0,
        }}
    }}

    #[doc = \"Empty container with a per-field capacity hint.\"]
    #[must_use]
    {vis} fn with_capacity(cap: usize) -> Self {{
        Self {{
{cap_fields}            len: 0,
        }}
    }}

    #[doc = \"Elements stored.\"]
    #[must_use]
    {vis} fn len(&self) -> usize {{ self.len }}

    #[doc = \"True when no elements are stored.\"]
    #[must_use]
    {vis} fn is_empty(&self) -> bool {{ self.len == 0 }}

    #[doc = \"Elements storable without leaf reallocation (minimum across leaf fields).\"]
    #[must_use]
    {vis} fn capacity(&self) -> usize {{ {capacity_expr} }}

    #[doc = \"Drop all elements (allocations kept).\"]
    {vis} fn clear(&mut self) {{
{clear_stmts}        self.len = 0;
    }}

    #[doc = \"Ensure room for `additional` more elements in every field.\"]
    {vis} fn reserve(&mut self, additional: usize) {{
{reserve_stmts}    }}

    #[doc = \"Append one value, scattered across the field buffers.\"]
    {vis} fn push(&mut self, value: {value_ty}) {{
        let {name} {{ {destructure} }} = value;
{push_stmts}        self.len += 1;
    }}

    #[doc = \"Gather element `i` back into a value. Panics if out of bounds.\"]
    #[must_use]
    {vis} fn get(&self, i: usize) -> {value_ty} {{
        assert!(i < self.len, \"index {{i}} out of bounds (len {{}})\", self.len);
        {name} {{
{get_fields}        }}
    }}

    #[doc = \"Scatter `value` into slot `i`. Panics if out of bounds.\"]
    {vis} fn set(&mut self, i: usize, value: {value_ty}) {{
        assert!(i < self.len, \"index {{i}} out of bounds (len {{}})\", self.len);
        let {name} {{ {destructure} }} = value;
{set_stmts}    }}

    #[doc = \"Zip-style iteration: gathers one value per element.\"]
    {vis} fn iter(&self) -> impl Iterator<Item = {value_ty}> + '_ {{
        (0..self.len).map(|i| self.get(i))
    }}

    #[doc = \"Zero-copy view descriptors for every leaf field (dotted paths for nested \
containers) — the FrankenNumpy membrane shape.\"]
    #[must_use]
    {vis} fn field_views(&self) -> ::std::vec::Vec<::fs_soa::RawView> {{
        let mut out = ::std::vec::Vec::new();
        <Self as ::fs_soa::SoaContainer<{value_ty}>>::c_views(self, \"\", &mut out);
        out
    }}

    #[doc = \"Address-free layout description (one JSON line per leaf field), for \
auditability logs.\"]
    #[must_use]
    {vis} fn layout_descr() -> ::std::string::String {{
        let mut out = ::std::vec::Vec::new();
        <Self as ::fs_soa::SoaContainer<{value_ty}>>::c_layout(\"\", &mut out);
        out.join(\"\\n\")
    }}

{accessors}}}

#[automatically_derived]
impl {generics} ::std::default::Default for {soa_ty} {where_clause} {{
    fn default() -> Self {{ Self::new() }}
}}

#[automatically_derived]
impl {generics} ::fs_soa::SoaAble for {value_ty} {where_clause} {{
    type Soa = {soa_ty};
}}

#[automatically_derived]
impl {generics} ::fs_soa::SoaContainer<{value_ty}> for {soa_ty} {where_clause} {{
    fn c_new() -> Self {{ Self::new() }}
    fn c_with_capacity(cap: usize) -> Self {{ Self::with_capacity(cap) }}
    fn c_len(&self) -> usize {{ self.len }}
    fn c_push(&mut self, value: {value_ty}) {{ self.push(value); }}
    fn c_get(&self, i: usize) -> {value_ty} {{ self.get(i) }}
    fn c_set(&mut self, i: usize, value: {value_ty}) {{ self.set(i, value); }}
    fn c_clear(&mut self) {{ self.clear(); }}
    fn c_reserve(&mut self, additional: usize) {{ self.reserve(additional); }}
    fn c_views(&self, prefix: &str, out: &mut ::std::vec::Vec<::fs_soa::RawView>) {{
{views_stmts}    }}
    fn c_layout(prefix: &str, out: &mut ::std::vec::Vec<::std::string::String>) {{
{layout_stmts}    }}
}}
"
    )
}
