#!/usr/bin/env python3
"""Triage table of every top-level fs-couple and fs-material module (bead rc-1109 .8.2).

usage: python3 docs/coupling-triage/triage.py [REPO_ROOT [OUT]]

Writes OUT (default docs/coupling-triage/triage.tsv under REPO_ROOT, itself by default two levels above this script)
and prints a summary.
Python 3 standard library only; it reads the files on disk, so run it on a clean checkout of the commit it describes.

One row per module declared at the top of crates/fs-couple/src/lib.rs or crates/fs-material/src/lib.rs. Columns:
  crate, module, path        where the module lives (`#[path]` respected; inline modules point at lib.rs)
  lines, files               every .rs line and file in the module's tree (its file and its directory)
  unit_tests                 `#[test]` attributes inside the tree
  test_files                 files under the crate's tests/ that name the module
  features                   cargo features in the `#[cfg(...)]` on its declaration
  consumers                  who names the module from outside its tree, by kind, with file counts:
                               music_render[:<input file>] (the bin and its per-mode input files), bin:<name>,
                               cli:<file> or cli:bin/<name> (fs-cli src: the verb's file), crate:<name> (another
                               crate's src),
                               app:<name>, module:<crate>/<module> (a sibling module's live src),
                               example:<name>, tests, src-tests (`#[cfg(test)]` code in either crate's src),
                               other-tests:<crate>, other:<dir>; a crate root re-export (`pub use module::Item`)
                               counts as the module when the item is named
  reached_via                the lanes that reach the module through live code: instrument (music_render, the
                               bowed bins, apps/instrument-*), product (cli, other crates' src, other apps); a
                               sibling module passes on what reaches it; examples and tests reach nothing
  class                      K kernel-forming (the list in KERNEL below, each with the reason read from the module's
                               own documentation), I instrument lane (reached_via has instrument, with or without
                               product), P pairing with a live consumer (product only), U unreached (no live
                               consumer; it may still have examples and tests)
  kernel_refs                references to fs_system:: or fs_ports::component in the tree's live code
  owners                     the physics owners the tree's live code calls: other crates declaring layer L2 or L3
                               in [package.metadata.frankensim], without the kernel's fs-ports and fs-system
  owners_with_component      those of them whose src implements fs_ports' Component: the pieces that exist
  nearest_sibling            the module sharing the most normalized code lines with this one (either crate)
  shared_lines               the fraction of this module's distinct normalized code lines found in that sibling
  doc                        the first line of the module's `//!` documentation

Code under `#[cfg(test)]` (inline items and the files `#[cfg(test)] mod name;` declares) is test code: it is listed
as src-tests and reaches nothing. A private module (`mod name;`) nobody names, whose live code adds inherent
`impl Type` blocks to a type that exactly one other module defines, is listed as impl-of:<that module> and reached
with it.

What this does not do: reachability is per module, so a reached module passes it to everything its whole tree names,
whichever function the live consumer calls. Glob imports and macro-generated paths are not resolved, so a consumer
that reaches a module only through `use fs_couple::*` or a macro is missed; `super::` paths are not followed. The
connection laws a migration needs are not derived: the owners columns name the Components to compose, not the laws
between them. Duplication counts identical normalized lines of at least 30 characters, not near copies.
"""
import collections
import os
import re
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(sys.argv[1]) if len(sys.argv) > 1 else os.path.abspath(os.path.join(HERE, "..", ".."))
OUT = sys.argv[2] if len(sys.argv) > 2 else os.path.join(ROOT, "docs", "coupling-triage", "triage.tsv")
CRATES = {"fs-couple": "fs_couple", "fs-material": "fs_material"}
SKIP_DIRS = {"target", "node_modules"}
INSTRUMENT_ROOTS = ("music_render", "bin:", "app:instrument-")
# Kernel-forming modules: they define composition machinery rather than a physical pairing.
KERNEL = {
    ("fs-couple", "constitutive_graph"): "K: balance equations as data over graph ports, closed by the shared Newton-Krylov owner",
    ("fs-couple", "iqn_ils"): "K: interface quasi-Newton (IQN-ILS) acceleration of a partitioned coupling iteration",
    ("fs-couple", "port_adoption"): "K: port and account observations of accepted steps for embedding modules in a graph",
    ("fs-material", "graph"): "K: the executable ConstitutiveGraph and its law-node protocol",
}
# Crates whose declared layer is L2 or L3 are physics owners; the kernel's own contract crates are counted apart.
KERNEL_CRATES = {"fs-ports", "fs-system"}


def rs_files():
    found = []
    for top, dirs, files in os.walk(ROOT):
        dirs[:] = sorted(d for d in dirs if d not in SKIP_DIRS and not d.startswith("."))
        for name in sorted(files):
            if name.endswith(".rs"):
                found.append(os.path.relpath(os.path.join(top, name), ROOT))
    return found


def read(rel):
    with open(os.path.join(ROOT, rel), encoding="utf-8", errors="replace") as handle:
        return handle.read()


def declarations(crate):
    """Top-level modules of a crate's lib.rs: name, features, file, directory, inline lines, root re-exports."""
    lib = f"crates/{crate}/src/lib.rs"
    src = f"crates/{crate}/src"
    modules, reexports = [], {}
    attrs, depth, inline = [], 0, None
    for line in read(lib).splitlines():
        code = line.split("//", 1)[0].strip()
        if inline is not None:
            inline["lines"] += 1
            depth += code.count("{") - code.count("}")
            if depth == 0:
                modules.append(inline)
                inline = None
            continue
        if depth > 0:
            depth += code.count("{") - code.count("}")
            continue
        if not code:
            continue
        if code.startswith("#["):
            attrs.append(code)
            continue
        m = re.match(r"(pub(?:\([^)]*\))?\s+)?mod\s+(\w+)\s*(;|\{)", code)
        if m:
            public, name, kind = m.groups()
            features = sorted({f for a in attrs if a.startswith("#[cfg") for f in re.findall(r'feature\s*=\s*"([^"]+)"', a)})
            path = next((re.search(r'"([^"]+)"', a).group(1) for a in attrs if a.startswith("#[path")), None)
            if kind == "{":
                inline = {"name": name, "features": features, "file": lib, "dir": f"{src}/{name}", "lines": 1,
                          "public": bool(public)}
                depth = code.count("{") - code.count("}")
                if depth == 0:
                    modules.append(inline)
                    inline = None
            else:
                if path:
                    file = f"{src}/{path}"
                elif os.path.exists(os.path.join(ROOT, src, f"{name}.rs")):
                    file = f"{src}/{name}.rs"
                else:
                    file = f"{src}/{name}/mod.rs"
                directory = os.path.dirname(file) if file.endswith("/mod.rs") else file[: -len(".rs")]
                modules.append({"name": name, "features": features, "file": file, "dir": directory, "lines": 0,
                                "public": bool(public)})
            attrs = []
            continue
        m = re.match(r"pub\s+use\s+(\w+)::(\{[^}]*\}|\w+)\s*;", code)
        if m:
            module, items = m.groups()
            for item in re.findall(r"\w+", items):
                reexports[item] = module
        depth += code.count("{") - code.count("}")
        attrs = []
    return modules, reexports


def named_paths(text, prefix):
    """First segments after `prefix::`, including `prefix::{a, b::c, ...}` groups (depth-0 items only)."""
    found = set(re.findall(rf"\b{prefix}::(\w+)", text))
    for m in re.finditer(rf"\b{prefix}::\{{", text):
        depth, start, items = 1, m.end(), []
        i = start
        while i < len(text) and depth:
            c = text[i]
            if c == "{":
                depth += 1
            elif c == "}":
                depth -= 1
            elif c == "," and depth == 1:
                items.append(text[start:i])
                start = i + 1
            i += 1
        items.append(text[start : i - 1])
        for item in items:
            word = re.match(r"\s*(\w+)", item)
            if word:
                found.add(word.group(1))
    return found


def consumer_kind(rel, crate_of_file, owner_of_file):
    parts = rel.split("/")
    if owner_of_file:
        return f"module:{owner_of_file[0]}/{owner_of_file[1]}"
    if parts[0] == "crates" and len(parts) > 3:
        crate, area = parts[1], parts[2]
        if area == "src":
            stem = [p[: -len(".rs")] if p.endswith(".rs") else p for p in parts[3:]]
            if crate in CRATES and stem[0] == "bin":
                if stem[1] == "music_render":
                    return "music_render" if len(stem) == 2 else f"music_render:{stem[2]}"
                return f"bin:{stem[1]}"
            if crate in CRATES:
                return None  # lib.rs and files outside any module tree
            if crate == "fs-cli":
                return f"cli:bin/{stem[1]}" if stem[0] == "bin" and len(stem) > 1 else f"cli:{stem[0]}"
            return f"crate:{crate}"
        if crate == crate_of_file:
            if area == "examples":
                return f"example:{parts[3][:-3] if parts[3].endswith('.rs') else parts[3]}"
            return "tests" if area == "tests" else f"{area}"
        return f"other-tests:{crate}"
    if parts[0] == "apps" and len(parts) > 1:
        return f"app:{parts[1]}"
    return f"other:{parts[0]}"


TEST_ATTR = re.compile(r"#\[cfg\((?:all\([^)]*)?\btest\b[^\]]*\]")


def split_tests(rel, text):
    """The text with its `#[cfg(test)]` items removed, and the files its `#[cfg(test)] mod name;` lines declare."""
    live, test_files, cursor = [], [], 0
    for m in TEST_ATTR.finditer(text):
        if m.start() < cursor:
            continue
        i = m.end()
        path = None
        while True:  # further attributes on the same item
            rest = text[i:].lstrip()
            i = len(text) - len(rest)
            if not rest.startswith("#["):
                break
            close = text.find("]", i)
            attr = text[i : close + 1]
            if attr.startswith("#[path"):
                path = re.search(r'"([^"]+)"', attr).group(1)
            i = close + 1
        decl = re.match(r"(?:pub(?:\([^)]*\))?\s+)?mod\s+(\w+)\s*;", text[i:])
        if decl:
            base = os.path.dirname(rel) if rel.endswith(("/mod.rs", "/lib.rs", "/main.rs")) else rel[: -len(".rs")]
            if path:
                test_files.append(os.path.normpath(os.path.join(os.path.dirname(rel), path)))
            else:
                test_files += [f"{base}/{decl.group(1)}.rs", f"{base}/{decl.group(1)}/mod.rs"]
            end = i + decl.end()
        else:
            semicolon, brace = text.find(";", i), text.find("{", i)
            if brace == -1 or (semicolon != -1 and semicolon < brace):
                end = semicolon + 1 if semicolon != -1 else len(text)
            else:
                depth, end = 0, brace
                while end < len(text):
                    if text[end] == "{":
                        depth += 1
                    elif text[end] == "}":
                        depth -= 1
                        if depth == 0:
                            break
                    end += 1
                end += 1
        live.append(text[cursor : m.start()])
        cursor = end
    live.append(text[cursor:])
    return "".join(live), test_files


def normalized_lines(text):
    lines = set()
    for raw in text.splitlines():
        line = " ".join(raw.split())
        if len(line) >= 30 and not line.startswith(("//", "#[", "use ", "pub use ")):
            lines.add(line)
    return lines


def main():
    files = rs_files()
    texts = {rel: read(rel) for rel in files}
    rows, reexports = [], {}
    for crate in CRATES:
        modules, crate_reexports = declarations(crate)
        reexports[crate] = crate_reexports
        for module in modules:
            module["crate"] = crate
            prefix = module["dir"] + "/"
            tree = [module["file"]] if module["lines"] == 0 else []
            tree += [rel for rel in files if rel.startswith(prefix) and rel != module["file"]]
            module["tree"] = sorted(set(tree))
            rows.append(module)
    owner = {}
    for row in rows:
        for rel in row["tree"]:
            owner[rel] = (row["crate"], row["name"])
    by_key = {(row["crate"], row["name"]): row for row in rows}
    layers = {}
    for name in sorted(os.listdir(os.path.join(ROOT, "crates"))):
        manifest = os.path.join(ROOT, "crates", name, "Cargo.toml")
        if os.path.exists(manifest):
            with open(manifest, encoding="utf-8") as handle:
                found = re.search(r'^layer\s*=\s*"([^"]+)"', handle.read(), flags=re.M)
            layers[name] = found.group(1) if found else None
    component_crates = set()
    for rel, text in texts.items():
        parts = rel.split("/")
        if parts[0] == "crates" and len(parts) > 3 and parts[2] == "src":
            if re.search(r"impl(<[^>]*>)?\s+(fs_ports::component::)?Component(<[^>]*>)?\s+for\b", text):
                component_crates.add(parts[1])
    live, test_only = {}, set()
    for rel, text in texts.items():
        live[rel], declared = split_tests(rel, text)
        test_only.update(declared)
    test_prefixes = tuple(
        (os.path.dirname(f) if f.endswith("/mod.rs") else f[: -len(".rs")]) + "/" for f in test_only
    )

    def used_modules(text, crate, crate_of):
        named = named_paths(text, CRATES[crate])
        if crate_of == crate:
            named |= named_paths(text, "crate")
        used = set()
        for word in named:
            module = word if (crate, word) in by_key else reexports[crate].get(word)
            if module and (crate, module) in by_key:
                used.add(module)
        return used

    consumers = collections.defaultdict(collections.Counter)
    edges = collections.defaultdict(set)  # used module -> modules that use it
    roots = collections.defaultdict(set)  # module -> lanes of live consumers
    test_files = collections.Counter()
    for rel, text in texts.items():
        crate_of = rel.split("/")[1] if rel.startswith("crates/") else None
        here = owner.get(rel)
        in_src = rel.split("/")[2:3] == ["src"]
        test_file = rel in test_only or rel.startswith(test_prefixes)
        for crate in CRATES:
            live_used = set() if test_file else used_modules(live[rel], crate, crate_of)
            test_used = used_modules(text, crate, crate_of) - live_used
            for module in test_used:
                key = (crate, module)
                if here != key and in_src:
                    # Test code inside src/: a test, never a live consumer.
                    consumers[key]["src-tests" if crate_of in CRATES else f"other-tests:{crate_of}"] += 1
            for module in live_used:
                key = (crate, module)
                if here == key:
                    continue
                kind = consumer_kind(rel, crate, here)
                if kind is None:
                    continue
                consumers[key][kind] += 1
                if kind == "tests":
                    test_files[key] += 1
                if here:
                    edges[key].add(here)
                elif kind.startswith(INSTRUMENT_ROOTS):
                    roots[key].add("instrument")
                elif kind.startswith(("cli:", "crate:", "app:")):
                    roots[key].add("product")
    # A module nobody names may only add inherent `impl Type` blocks to a type a sibling module defines: its
    # methods are then that sibling's, and it is reached with it.
    defined = collections.defaultdict(set)
    for key, row in by_key.items():
        for rel in row["tree"]:
            for name in re.findall(r"^\s*(?:pub(?:\([^)]*\))?\s+)?(?:struct|enum)\s+(\w+)", live[rel], flags=re.M):
                defined[(key[0], name)].add(key)
    for key, row in by_key.items():
        if consumers[key] or row["public"]:
            continue
        for rel in row["tree"]:
            for name in re.findall(r"^\s*impl(?:<[^>{]*>)?\s+(\w+)\s*(?:<[^>{]*>)?\s*\{", live[rel], flags=re.M):
                homes = defined[(key[0], name)]
                if len(homes) == 1 and key not in homes:
                    home = next(iter(homes))
                    consumers[key][f"impl-of:{home[0]}/{home[1]}"] += 1
                    edges[key].add(home)
    # A module is reached through each lane that reaches a module using it.
    reached = {key: set(lanes) for key, lanes in roots.items()}
    changed = True
    while changed:
        changed = False
        for used, users in edges.items():
            for user in users:
                for lane in reached.get(user, ()):
                    if lane not in reached.setdefault(used, set()):
                        reached[used].add(lane)
                        changed = True
    # Duplication: distinct normalized code lines shared with the nearest sibling.
    line_sets = {key: set().union(*(normalized_lines(texts[rel]) for rel in by_key[key]["tree"] if rel in texts))
                 if by_key[key]["tree"] else set() for key in by_key}
    index = collections.defaultdict(list)
    for key, lines in line_sets.items():
        for line in lines:
            index[line].append(key)
    out = []
    for key in sorted(by_key, key=lambda k: (k[0], -sum(texts[r].count("\n") for r in by_key[k]["tree"] if r in texts), k[1])):
        row = by_key[key]
        tree_texts = [texts[rel] for rel in row["tree"] if rel in texts]
        lines = sum(text.count("\n") for text in tree_texts) + row["lines"]
        unit_tests = sum(len(re.findall(r"#\[test\]", text)) for text in tree_texts)
        shared = collections.Counter()
        for line in line_sets[key]:
            for other in index[line]:
                if other != key:
                    shared[other] += 1
        # The most shared lines; ties go to the first sibling in (crate, module) order.
        sibling, count = min(shared.items(), key=lambda item: (-item[1], item[0])) if shared else (None, 0)
        fraction = count / len(line_sets[key]) if line_sets[key] else 0.0
        own = row["crate"]
        called = collections.Counter()
        kernel_refs = 0
        for rel in row["tree"]:
            for ident in re.findall(r"\b(fs_[a-z0-9_]+)::", live[rel]):
                name = ident.replace("_", "-")
                if name != own and name not in KERNEL_CRATES and layers.get(name) in ("L2", "L3"):
                    called[name] += 1
            kernel_refs += len(re.findall(r"\bfs_system::|\bfs_ports::component\b", live[rel]))
        doc = ""
        for text in tree_texts[:1] if row["lines"] == 0 else []:
            first = re.search(r"^//!\s?(.*\S)", text, flags=re.M)
            doc = first.group(1) if first else ""
        lanes = sorted(reached.get(key, ()))
        if key in KERNEL:
            klass = "K"
        elif "instrument" in lanes:
            klass = "I"
        elif "product" in lanes:
            klass = "P"
        else:
            klass = "U"
        cells = [
            row["crate"],
            row["name"],
            os.path.relpath(os.path.join(ROOT, row["file"]), os.path.join(ROOT, f"crates/{row['crate']}")),
            str(lines),
            str(len(row["tree"])),
            str(unit_tests),
            str(test_files[key]),
            ",".join(row["features"]) or "-",
            ";".join(f"{kind} {n}" for kind, n in sorted(consumers[key].items())) or "none",
            ",".join(lanes) or "-",
            klass,
            str(kernel_refs),
            ",".join(sorted(called)) or "-",
            ",".join(sorted(c for c in called if c in component_crates)) or "-",
            f"{sibling[0]}/{sibling[1]}" if sibling else "-",
            f"{fraction:.3f}",
            (KERNEL[key] + "; " if key in KERNEL else "") + doc[:160],
        ]
        out.append("\t".join(cell.replace("\t", " ") for cell in cells))
    header = ("crate\tmodule\tpath\tlines\tfiles\tunit_tests\ttest_files\tfeatures\tconsumers\treached_via\tclass\t"
              "kernel_refs\towners\towners_with_component\tnearest_sibling\tshared_lines\tdoc")
    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    with open(OUT, "w", encoding="utf-8") as handle:
        handle.write(header + "\n" + "\n".join(out) + "\n")
    summary = collections.defaultdict(lambda: [0, 0])
    for line in out:
        cells = line.split("\t")
        summary[(cells[0], cells[10])][0] += 1
        summary[(cells[0], cells[10])][1] += int(cells[3])
    print(f"wrote {OUT}: {len(out)} rows; Component crates: {', '.join(sorted(component_crates))}")
    for (crate, klass), (count, lines) in sorted(summary.items()):
        print(f"{crate}\t{klass}\t{count} modules\t{lines} lines")


if __name__ == "__main__":
    main()
