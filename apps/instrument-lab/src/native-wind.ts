import "./native-wind.css";
// Canonical native fixtures, bundled verbatim by Vite; no copied physics defaults.
import reedExampleSource from "../../../crates/fs-couple/examples/reed-two-medium-tone-hole.performance?raw";
import flueExampleSource from "../../../crates/fs-couple/examples/flue-nolle.performance?raw";
import { nativeWindStatus, renderNativeWindBundle, type NativeWindStatus, type NativeWindAttachment, type NativeWindGasSources } from "./native-wind-client.mjs";
import { el } from "./util";
import { createNativeGasControls } from "./native-wind-gases";

const EXAMPLES = [
  {
    name: "reed-two-medium-tone-hole.performance",
    label: "Reed — two-medium tone hole",
    source: reedExampleSource,
    intro: "Starter: two-medium tone hole — a short boundary diagnostic with a reed and three tube sections in different gas regions. Change the gas mixtures, temperature, geometry or gestures in the description.",
    help: "ambient sets the inlet and chimney gas; a segment’s gas entry overrides its local medium; exterior sets the surrounding gas. The interior and exterior must share a static pressure because this boundary model omits mean flow. This short numerical fixture is not a measured clarinet or oboe.",
  },
  {
    name: "flue-nolle.performance",
    label: "Flute — Nolle/Price jet",
    source: flueExampleSource,
    intro: "Starter: flute jet — a short blowing-pressure phrase through a cylindrical tube, using the native Nolle/Price jet transfer. Change the gas, tube, mouth geometry or pressure gesture in the description.",
    help: "gas owns the internal medium; exterior sets the surrounding gas. jet_transfer selects the short-channel Nolle/Price source. The native description owns the mechanical and output clocks, and the jet and acoustic mouth are solved together. This synthetic example is not a measured flute.",
  },
];
const exampleUrl = (name: string): string =>
  `https://github.com/Dicklesworthstone/frankensim/blob/main/crates/fs-couple/examples/${name}`;
const LOCAL_FILE_LIMIT = 1_048_576; // Bound file reading even on a static host without a renderer.
const message = (error: unknown): string => error instanceof Error ? error.message : String(error);

/** A source editor and native render client. Synthesis stays in the Rust renderer. */
export function createNativeWindPanel(): HTMLElement {
  let selectedExample = EXAMPLES[0]!;
  const starter = el("select", { id: "native-wind-starter", class: "btn" },
    ...EXAMPLES.map((example) => el("option", { value: example.name }, example.label)));
  const introduction = el("p", { class: "nw-intro" }, selectedExample.intro);
  const sourceHelp = el("p", { id: "native-wind-source-help" }, selectedExample.help);
  const canonical = el("a", { href: exampleUrl(selectedExample.name), target: "_blank", rel: "noopener" }, "Canonical example");
  const source = el("textarea", { id: "native-wind-source", class: "nw-source", rows: "15",
    spellcheck: "false", autocapitalize: "off", autocomplete: "off", wrap: "off",
    "aria-describedby": "native-wind-source-help" });
  source.value = selectedExample.source;
  // A textarea normalizes CRLF. Keep the loaded UTF-8 text until the user edits it.
  let loadedSource = selectedExample.source;
  let loadedEditorValue = source.value;
  const currentSource = (): string => source.value === loadedEditorValue ? loadedSource : source.value;
  const sourceInfo = el("span", { class: "small muted" });
  const connection = el("span", { class: "nw-connection", role: "status", "aria-live": "polite" }, "Checking local renderer…");
  const result = el("p", { class: "nw-result small", role: "status", "aria-live": "polite", "aria-atomic": "true" });
  const render = el("button", { type: "button", class: "btn primary" }, "Render WAV");
  const cancel = el("button", { type: "button", class: "btn" }, "Cancel render");
  const check = el("button", { type: "button", class: "btn" }, "Check renderer");
  const open = el("button", { type: "button", class: "btn" }, "Open performance…");
  const reset = el("button", { type: "button", class: "btn" }, "Reset example");
  const attach = el("button", { type: "button", class: "btn" }, "Attach source files…");
  const attachmentFile = el("input", { type: "file", multiple: "",
    "aria-label": "Attach material packs or other referenced native source files" });
  attachmentFile.hidden = true;
  const decimate = el("input", { type: "checkbox", id: "native-wind-decimate" });
  const decimation = el("label", { for: "native-wind-decimate", class: "small nw-decimation" }, decimate,
    " Use native decimation for a performance above 48 kHz (mechanics keep their original clock)");
  const attachmentList = el("div", { class: "nw-attachments" });
  const attachmentInfo = el("p", { class: "small muted" });
  const attachmentPanel = el("section", { "aria-label": "Attached source files" },
    attachmentInfo, attachmentList);
  type Attached = NativeWindAttachment & { sha256: string };
  let attachments: Attached[] = [];
  const attachmentIdentity = (): string => JSON.stringify(attachments.map(file => [file.path, file.sha256]));
  const gas = createNativeGasControls(() => edited());
  const file = el("input", { type: "file", accept: ".performance,.txt,text/plain",
    "aria-label": "Open a native wind performance" });
  file.hidden = true;
  const audio = el("audio", { controls: "", preload: "metadata", "aria-label": "Your native wind render" });
  const download = el("a", { class: "btn", download: "wind-render.wav" }, "Download WAV");
  const renderedFrom = el("p", { class: "small muted" });
  const renderedFiles = el("p", { class: "small muted nw-file-identities" });
  const playback = el("div", { class: "nw-playback" }, audio, download, renderedFrom, renderedFiles);
  playback.hidden = true;
  cancel.hidden = true;

  let status: NativeWindStatus | null = null;
  let checking: AbortController | null = null;
  let rendering: AbortController | null = null;
  let loadingFile = false;
  let previewUrl: string | null = null;
  let renderedSource: string | null = null;
  let renderedAttachmentIdentity: string | null = null;
  let renderedDecimate = false;
  let renderedGasIdentity: string | null = null;
  let sourceName = selectedExample.name;

  const controls = (): void => {
    const active = rendering !== null;
    render.disabled = active || loadingFile || checking !== null || !status?.available || status.busy || !source.value.trim();
    cancel.hidden = !active;
    check.disabled = active || checking !== null;
    open.disabled = active || loadingFile;
    reset.disabled = active || loadingFile;
    starter.disabled = active || loadingFile;
    attach.disabled = active || loadingFile;
    decimate.disabled = active || loadingFile;
    gas.setDisabled(active || loadingFile);
    attachmentList.querySelectorAll<HTMLInputElement | HTMLButtonElement>("input,button")
      .forEach(control => { control.disabled = active || loadingFile; });
    source.readOnly = active || loadingFile;
    source.setAttribute("aria-busy", String(active));
    render.textContent = active ? "Rendering…" : "Render WAV";
  };
  const edited = (): void => {
    sourceInfo.textContent = `${sourceName} · ${new TextEncoder().encode(currentSource()).byteLength.toLocaleString()} bytes`;
    if (renderedSource !== null) {
      renderedFrom.textContent = currentSource() === renderedSource && attachmentIdentity() === renderedAttachmentIdentity && decimate.checked === renderedDecimate && gas.identity() === renderedGasIdentity
        ? "This WAV was rendered from the performance and source files currently shown."
        : "This WAV is from your previous performance or source files. Render again to hear your edits.";
    }
    controls();
  };
  const clearPreview = (): void => {
    audio.pause();
    audio.removeAttribute("src");
    audio.load();
    download.removeAttribute("href");
    playback.hidden = true;
    renderedSource = null;
    renderedAttachmentIdentity = null;
    renderedGasIdentity = null;
    renderedFiles.textContent = "";
    if (previewUrl !== null) URL.revokeObjectURL(previewUrl);
    previewUrl = null;
  };
  const showResult = (text: string, refused = false): void => {
    result.textContent = text;
    result.classList.toggle("nw-refused", refused);
  };

  const checkConnection = async (): Promise<void> => {
    checking?.abort();
    const controller = new AbortController();
    checking = controller;
    status = null;
    connection.textContent = "Checking local renderer…";
    controls();
    const timeout = window.setTimeout(() => controller.abort(), 5000);
    try {
      const next = await nativeWindStatus(controller.signal);
      if (checking !== controller) return;
      status = next;
      connection.textContent = !next.available ? "Native renderer unavailable"
        : next.busy ? "Local renderer is busy" : "Local renderer ready";
      if (!next.available) {
        showResult("Open the local Instrument Lab server to render and listen. You can edit the source here.");
      } else if (next.busy) {
        showResult("Another render is running. Check the renderer again when it finishes.");
      } else {
        showResult(`Ready to render. The local render limit is ${Math.ceil(next.timeoutMs / 1000)} seconds.`);
      }
    } catch {
      if (checking !== controller) return;
      connection.textContent = "Native renderer unavailable";
      showResult("This host has no native renderer connection. Open the local Instrument Lab server to render and listen.");
    } finally {
      window.clearTimeout(timeout);
      if (checking === controller) {
        checking = null;
        controls();
      }
    }
  };

  render.addEventListener("click", () => {
    if (!status?.available || status.busy || rendering || checking || loadingFile) return;
    const connected = status;
    const submitted = currentSource();
    const submittedFiles = attachments.map(file => ({ ...file }));
    const submittedAttachmentIdentity = attachmentIdentity();
    const submittedDecimate = decimate.checked;
    let submittedGas: NativeWindGasSources | undefined;
    try { submittedGas = gas.read(); }
    catch (error) { showResult(message(error), true); return; }
    const submittedGasIdentity = gas.identity();
    const controller = new AbortController();
    rendering = controller;
    audio.pause();
    showResult("Rendering your performance from its material, gas and geometry data…");
    controls();
    let timedOut = false;
    const timeout = window.setTimeout(() => {
      timedOut = true;
      controller.abort();
    }, connected.timeoutMs + 5000);
    void (async () => {
      try {
        const wave = await renderNativeWindBundle(submitted, submittedFiles, connected, {decimate:submittedDecimate, gasSources:submittedGas}, controller.signal);
        if (rendering !== controller || controller.signal.aborted) return;
        const nextUrl = URL.createObjectURL(wave);
        clearPreview();
        previewUrl = nextUrl;
        renderedSource = submitted;
        renderedAttachmentIdentity = submittedAttachmentIdentity;
        renderedDecimate = submittedDecimate;
        renderedGasIdentity = submittedGasIdentity;
        renderedFiles.textContent = submittedFiles.length
          ? "Rendered source files: " + submittedFiles.map(file => `${file.path} · ${file.bytes.byteLength.toLocaleString()} bytes · SHA-256 ${file.sha256}`).join("\n")
          : "This render used no attached source files.";
        if (submittedGas) renderedFiles.textContent += "\nGas associations: " + submittedGas.components
          .map(component => `${component.speciesPath} / ${component.nasa9Path} / ${component.transportPath}`).join("; ")
          + ` · phase ${submittedGas.phase} · water ${submittedGas.waterSpecies ?? "none"}`;
        audio.src = nextUrl;
        download.href = nextUrl;
        download.download = `${sourceName.replace(/\.[^.]*$/, "") || "wind"}-render.wav`;
        playback.hidden = false;
        edited();
        showResult("Render complete. Press play to listen, or download the WAV.");
      } catch (error) {
        if (rendering !== controller) return;
        showResult(controller.signal.aborted
          ? timedOut ? "The render timed out. Check the local renderer before trying again." : "Render cancelled."
          : message(error), !controller.signal.aborted || timedOut);
      } finally {
        window.clearTimeout(timeout);
        if (rendering === controller) {
          rendering = null;
          controls();
        }
      }
    })();
  });
  cancel.addEventListener("click", () => {
    rendering?.abort();
    showResult("Cancelling render…");
  });
  check.addEventListener("click", () => void checkConnection());
  source.addEventListener("input", edited);
  open.addEventListener("click", () => file.click());
  const loadExample = (): void => {
    source.value = selectedExample.source;
    loadedSource = selectedExample.source;
    loadedEditorValue = source.value;
    sourceName = selectedExample.name;
    introduction.textContent = selectedExample.intro;
    sourceHelp.textContent = selectedExample.help;
    canonical.href = exampleUrl(selectedExample.name);
    attachments = [];
    decimate.checked = false;
    gas.reset();
    drawAttachments();
    edited();
    showResult("Example loaded. Edit the source or render it.");
  };
  reset.addEventListener("click", loadExample);
  starter.addEventListener("change", () => {
    if (rendering || loadingFile) return;
    const selected = EXAMPLES.find((example) => example.name === starter.value);
    if (!selected) return;
    selectedExample = selected;
    loadExample();
  });
  const drawAttachments = (): void => {
    attachmentInfo.textContent = attachments.length
      ? `${attachments.length} source files · ${attachments.reduce((sum, file) => sum + file.bytes.byteLength, 0).toLocaleString()} bytes. Paths must match the performance exactly.`
      : "Attach the original material packs or artifacts referenced by your performance. Keep their source paths; edit a path below when the description uses a subdirectory.";
    attachmentList.replaceChildren();
    gas.setAttachments(attachments.map(file => file.path));
    for (const attached of attachments) {
      const path = el("input", { type: "text", value: attached.path, spellcheck: "false",
        "aria-label": `Source path for ${attached.path}` });
      path.addEventListener("input", () => {
        attached.path = path.value; gas.setAttachments(attachments.map(file => file.path)); edited();
      });
      const remove = el("button", { type: "button", class: "btn" }, "Remove");
      remove.addEventListener("click", () => {
        attachments = attachments.filter(file => file !== attached);
        drawAttachments(); edited();
      });
      attachmentList.append(el("div", { class: "nw-attachment" }, path,
        el("span", { class: "small muted", title: `SHA-256 ${attached.sha256}` },
          `${attached.bytes.byteLength.toLocaleString()} bytes · ${attached.sha256.slice(0, 12)}…`), remove));
    }
    controls();
  };
  decimate.addEventListener("change", edited);
  attach.addEventListener("click", () => attachmentFile.click());
  attachmentFile.addEventListener("change", () => {
    const selected = Array.from(attachmentFile.files ?? []);
    attachmentFile.value = "";
    if (!selected.length || rendering || loadingFile) return;
    const maxFiles = Math.min(32, status?.maximumAttachments ?? 32);
    const maxFileBytes = Math.min(8 * 1024 * 1024, status?.maximumAttachmentFileBytes ?? 8 * 1024 * 1024);
    const maxBytes = Math.min(32 * 1024 * 1024, status?.maximumAttachmentBytes ?? 32 * 1024 * 1024);
    const total = attachments.reduce((sum, file) => sum + file.bytes.byteLength, 0)
      + selected.reduce((sum, file) => sum + file.size, 0);
    if (attachments.length + selected.length > maxFiles || total > maxBytes
        || selected.some(file => !file.size || file.size > maxFileBytes)) {
      showResult(`Choose at most ${maxFiles} nonempty files, ${maxFileBytes.toLocaleString()} bytes per file and ${maxBytes.toLocaleString()} bytes combined.`, true);
      return;
    }
    loadingFile = true; controls();
    void Promise.all(selected.map(async file => {
      const buffer = await file.arrayBuffer();
      const digest = await crypto.subtle.digest("SHA-256", buffer);
      const sha256 = Array.from(new Uint8Array(digest), byte => byte.toString(16).padStart(2, "0")).join("");
      return {path:file.name, bytes:new Uint8Array(buffer), sha256};
    })).then(added => {
      attachments.push(...added);
      drawAttachments(); edited();
      showResult("Source files attached. Check that each path matches the performance, then render.");
    }).catch(() => {
      showResult("Source files could not be read. Your previous attachments are unchanged.", true);
    }).finally(() => { loadingFile = false; controls(); });
  });
  file.addEventListener("change", () => {
    const selected = file.files?.[0];
    file.value = "";
    if (!selected || rendering || loadingFile) return;
    const limit = Math.min(LOCAL_FILE_LIMIT, status?.maximumSourceBytes ?? LOCAL_FILE_LIMIT);
    if (selected.size > limit) {
      showResult(`Choose a performance file no larger than ${limit.toLocaleString()} bytes.`, true);
      return;
    }
    loadingFile = true;
    controls();
    void selected.arrayBuffer().then((bytes) => {
      const text = new TextDecoder("utf-8", { fatal: true, ignoreBOM: true }).decode(bytes);
      source.value = text;
      loadedSource = text;
      loadedEditorValue = source.value;
      sourceName = selected.name;
      showResult("Performance loaded. Edit the source or render it.");
      edited();
    }).catch(() => {
      showResult("This performance could not be read as UTF-8 text. Your previous source is still in the editor.", true);
    }).finally(() => {
      loadingFile = false;
      controls();
    });
  });
  window.addEventListener("pagehide", () => {
    checking?.abort();
    rendering?.abort();
    clearPreview();
  });
  window.addEventListener("pageshow", (event) => {
    if (event.persisted) void checkConnection();
  });

  const panel = el("section", { class: "panel nw-panel", "aria-labelledby": "native-wind-title" },
    el("div", { class: "panel-head" },
      el("div", {}, el("h3", { id: "native-wind-title" }, "Render a wind performance"),
        el("p", { class: "small muted" }, "Edit the source, run the native physics, and listen to your result.")),
      connection),
    el("div", { class: "nw-source-head" }, el("label", { for: "native-wind-starter" }, "Starter performance"), starter),
    introduction,
    el("details", { class: "nw-help small" },
      el("summary", {}, "About the source and local renderer"),
      sourceHelp,
      el("p", {}, "Reed, flue and material-valve descriptions use the original native loader. Attach the exact material packs and artifacts named by the source, including subdirectory paths. Changing a pack changes its physical properties; missing data and unsupported states remain native refusals. Keep the same declared pressure scale when comparing renders."),
      el("p", {}, "Rendering uses the local music_render executable. The hosted review page can show and edit the source; audio is available when you open it through the local renderer server."),
      el("p", {}, canonical,
        " · ", el("a", { href: "https://github.com/Dicklesworthstone/frankensim/tree/main/apps/instrument-lab#native-wind-renderer", target: "_blank", rel: "noopener" }, "Local renderer setup"))),
    el("div", { class: "nw-source-head" }, el("label", { for: "native-wind-source" }, "Performance source"),
      sourceInfo),
    source, attachmentPanel, gas.element, decimation,
    el("div", { class: "btn-row nw-actions" }, render, cancel, open, attach, reset, check, file, attachmentFile),
    result, playback);
  drawAttachments();
  edited();
  void checkConnection();
  return panel;
}
