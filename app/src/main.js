const { invoke } = window.__TAURI__.core;
const { listen } = window.__TAURI__.event;

const defaultSettings = () => ({
  active: "false",
  originalCHM: "",
  originalCHMBand: "1",
  doSmoothing: "true",
  doTops: "true",
  doCrowns: "true",
  smoothWindowSize: "3",
  smoothSigma: "1.000000",
  smoothedCHM: "",
  treetopsDatabase: "",
  crownsRaster: "",
  crownsDatabase: "",
  crownsDoDatabase: "true",
  crownsUpdateHeights: "true",
  crownsRemoveHoles: "true",
  crownsRemoveDangles: "true",
  crownsKeepSmoothed: "true",
  topsMaxNulls: "0.5",
  smoothedCHMDriver: "GTiff",
  crownsRasterDriver: "GTiff",
  treetopsDatabaseDriver: "Spatialite",
  crownsDatabaseDriver: "Spatialite",
  topsThresholds: [
    { threshold: 3.0, window: 3 },
    { threshold: 5.0, window: 7 },
  ],
  crownsThresholds: [
    { threshold: 5.0, fraction: 0.4, radius: 5.0 },
    { threshold: 8.0, fraction: 0.4, radius: 7.0 },
  ],
});

let running = false;

function $(id) {
  return document.getElementById(id);
}

function dirname(path) {
  const idx = Math.max(path.lastIndexOf("/"), path.lastIndexOf("\\"));
  return idx >= 0 ? path.slice(0, idx) : "";
}

function basename(path) {
  const idx = Math.max(path.lastIndexOf("/"), path.lastIndexOf("\\"));
  return idx >= 0 ? path.slice(idx + 1) : path;
}

function stem(path) {
  const name = basename(path);
  const dot = name.lastIndexOf(".");
  return dot >= 0 ? name.slice(0, dot) : name;
}

function join(dir, file) {
  if (!dir) return file;
  const sep = dir.includes("\\") ? "\\" : "/";
  return dir.endsWith(sep) ? `${dir}${file}` : `${dir}${sep}${file}`;
}

function deriveOutputPaths(chmPath) {
  const dir = dirname(chmPath);
  const base = stem(chmPath);
  $("smoothed-chm").value = join(dir, `${base}_smooth.tif`);
  $("treetops-database").value = join(dir, "tops.sqlite");
  $("crowns-raster").value = join(dir, `${base}_crowns.tif`);
  $("crowns-database").value = join(dir, `${base}_crowns.sqlite`);
}

function boolToString(value) {
  return value ? "true" : "false";
}

function renderTopThresholds(items) {
  const container = $("top-thresholds");
  container.innerHTML = "";
  items.forEach((item, index) => {
    const row = document.createElement("div");
    row.className = "threshold-row";
    row.innerHTML = `
      <label class="field">
        <span>Height (m)</span>
        <input type="number" step="0.1" data-kind="top" data-index="${index}" data-field="threshold" value="${item.threshold}" />
      </label>
      <label class="field">
        <span>Window (px)</span>
        <input type="number" step="1" data-kind="top" data-index="${index}" data-field="window" value="${item.window}" />
      </label>
      <button type="button" class="secondary small" data-remove-top="${index}">Remove</button>
    `;
    container.appendChild(row);
  });
}

function renderCrownThresholds(items) {
  const container = $("crown-thresholds");
  container.innerHTML = "";
  items.forEach((item, index) => {
    const row = document.createElement("div");
    row.className = "threshold-row";
    row.innerHTML = `
      <label class="field">
        <span>Height (m)</span>
        <input type="number" step="0.1" data-kind="crown" data-index="${index}" data-field="threshold" value="${item.threshold}" />
      </label>
      <label class="field">
        <span>Fraction</span>
        <input type="number" step="0.05" min="0" max="1" data-kind="crown" data-index="${index}" data-field="fraction" value="${item.fraction}" />
      </label>
      <label class="field">
        <span>Radius (m)</span>
        <input type="number" step="0.1" data-kind="crown" data-index="${index}" data-field="radius" value="${item.radius}" />
      </label>
      <button type="button" class="secondary small" data-remove-crown="${index}">Remove</button>
    `;
    container.appendChild(row);
  });
}

function readTopThresholds() {
  return [...document.querySelectorAll("[data-kind='top'][data-field]")].reduce((acc, input) => {
    const index = Number(input.dataset.index);
    const field = input.dataset.field;
    acc[index] = acc[index] || {};
    acc[index][field] = Number(input.value);
    return acc;
  }, []).filter(Boolean);
}

function readCrownThresholds() {
  return [...document.querySelectorAll("[data-kind='crown'][data-field]")].reduce((acc, input) => {
    const index = Number(input.dataset.index);
    const field = input.dataset.field;
    acc[index] = acc[index] || {};
    acc[index][field] = Number(input.value);
    return acc;
  }, []).filter(Boolean);
}

function collectSettings() {
  const originalCHM = $("original-chm").value.trim();
  const settings = {
    ...defaultSettings(),
    originalCHM,
    originalCHMBand: $("original-chm-band").value,
    doSmoothing: boolToString($("do-smoothing").checked),
    doTops: boolToString($("do-tops").checked),
    doCrowns: boolToString($("do-crowns").checked),
    smoothWindowSize: $("smooth-window-size").value,
    smoothSigma: Number($("smooth-sigma").value).toFixed(6),
    smoothedCHM: $("smoothed-chm").value.trim(),
    treetopsDatabase: $("treetops-database").value.trim(),
    crownsRaster: $("crowns-raster").value.trim(),
    crownsDatabase: $("crowns-database").value.trim(),
    topsThresholds: readTopThresholds(),
    crownsThresholds: readCrownThresholds(),
  };

  if (originalCHM) {
    settings.topsWindowsRaster = join(dirname(originalCHM), "tops_windows.tif");
    settings.topsIdsRaster = join(dirname(originalCHM), "tops_ids.tif");
  }

  return settings;
}

function applySettings(settings) {
  $("original-chm").value = settings.originalCHM || "";
  $("original-chm-band").value = settings.originalCHMBand || "1";
  $("do-smoothing").checked = settings.doSmoothing !== "false";
  $("do-tops").checked = settings.doTops !== "false";
  $("do-crowns").checked = settings.doCrowns !== "false";
  $("smooth-window-size").value = settings.smoothWindowSize || "3";
  $("smooth-sigma").value = settings.smoothSigma || "1";
  $("smoothed-chm").value = settings.smoothedCHM || "";
  $("treetops-database").value = settings.treetopsDatabase || "";
  $("crowns-raster").value = settings.crownsRaster || "";
  $("crowns-database").value = settings.crownsDatabase || "";
  renderTopThresholds(settings.topsThresholds || defaultSettings().topsThresholds);
  renderCrownThresholds(settings.crownsThresholds || defaultSettings().crownsThresholds);
}

function setStatus({ phase, progress, message }) {
  $("status-phase").textContent = phase || "idle";
  $("status-message").textContent = message || "";
  $("progress-fill").style.width = `${Math.max(0, Math.min(100, progress || 0))}%`;
  $("progress-label").textContent = `${Math.max(0, Math.min(100, progress || 0))}%`;
  document.querySelector(".status-bar").dataset.state = phase || "idle";
}

function appendLog(line, stream = "stdout") {
  const output = $("log-output");
  const prefix = stream === "stderr" ? "[stderr] " : "";
  output.textContent += `${prefix}${line}\n`;
  output.scrollTop = output.scrollHeight;
}

function setRunning(isRunning) {
  running = isRunning;
  $("run-btn").disabled = isRunning;
  $("cancel-btn").disabled = !isRunning;
}

async function browse(targetId, mode) {
  const title = mode === "save" ? "Save settings file" : "Select file";
  const picker = mode === "save" ? "pick_save_file" : "pick_open_file";
  const selected = await invoke(picker, { title });
  if (selected) {
    $(targetId).value = selected;
    if (targetId === "original-chm") {
      deriveOutputPaths(selected);
    }
  }
}

async function loadSettings() {
  const selected = await invoke("pick_open_file", { title: "Load settings file" });
  if (!selected) {
    return;
  }

  $("config-path").value = selected;
  const content = await invoke("read_settings_file", { path: selected });
  applySettings(JSON.parse(content));
  setStatus({ phase: "idle", progress: 0, message: `Loaded ${selected}` });
}

async function saveSettings() {
  const path = $("config-path").value.trim();
  if (!path) {
    setStatus({ phase: "error", progress: 0, message: "Choose a settings file path first." });
    return;
  }
  const content = JSON.stringify(collectSettings(), null, 2);
  await invoke("write_settings_file", { path, content });
  setStatus({ phase: "idle", progress: 0, message: `Saved ${path}` });
}

async function runProcessing() {
  const configPath = $("config-path").value.trim();
  if (!configPath) {
    setStatus({ phase: "error", progress: 0, message: "Set a config file path before running." });
    return;
  }
  if (!$("original-chm").value.trim()) {
    setStatus({ phase: "error", progress: 0, message: "Select an input CHM before running." });
    return;
  }

  $("log-output").textContent = "";
  await saveSettings();
  setRunning(true);
  setStatus({ phase: "starting", progress: 0, message: "Launching treetops-cli..." });

  try {
    await invoke("run_processing", { configPath });
  } catch (error) {
    setRunning(false);
    setStatus({ phase: "error", progress: 0, message: String(error) });
  }
}

async function cancelProcessing() {
  await invoke("cancel_processing");
  setRunning(false);
  setStatus({ phase: "cancelled", progress: 0, message: "Processing cancelled." });
}

window.addEventListener("DOMContentLoaded", async () => {
  applySettings(defaultSettings());
  setStatus({ phase: "idle", progress: 0, message: "Ready" });

  $("load-config-btn").addEventListener("click", () => loadSettings().catch((e) => setStatus({ phase: "error", progress: 0, message: String(e) })));
  $("save-config-btn").addEventListener("click", () => saveSettings().catch((e) => setStatus({ phase: "error", progress: 0, message: String(e) })));
  $("run-btn").addEventListener("click", () => runProcessing());
  $("cancel-btn").addEventListener("click", () => cancelProcessing().catch((e) => setStatus({ phase: "error", progress: 0, message: String(e) })));

  document.querySelectorAll(".browse-btn").forEach((button) => {
    button.addEventListener("click", () => {
      browse(button.dataset.target, button.dataset.mode).catch((e) => setStatus({ phase: "error", progress: 0, message: String(e) }));
    });
  });

  $("original-chm").addEventListener("change", (event) => {
    if (event.target.value.trim()) {
      deriveOutputPaths(event.target.value.trim());
    }
  });

  $("add-top-threshold").addEventListener("click", () => {
    const items = readTopThresholds();
    items.push({ threshold: 5.0, window: 5 });
    renderTopThresholds(items);
  });

  $("add-crown-threshold").addEventListener("click", () => {
    const items = readCrownThresholds();
    items.push({ threshold: 5.0, fraction: 0.4, radius: 5.0 });
    renderCrownThresholds(items);
  });

  document.body.addEventListener("click", (event) => {
    const topIndex = event.target.dataset?.removeTop;
    if (topIndex !== undefined) {
      const items = readTopThresholds();
      items.splice(Number(topIndex), 1);
      renderTopThresholds(items);
    }
    const crownIndex = event.target.dataset?.removeCrown;
    if (crownIndex !== undefined) {
      const items = readCrownThresholds();
      items.splice(Number(crownIndex), 1);
      renderCrownThresholds(items);
    }
  });

  await listen("status", (event) => setStatus(event.payload));
  await listen("log", (event) => appendLog(event.payload.line, event.payload.stream));
  await listen("processing-finished", (event) => {
    setRunning(false);
    if (!event.payload) {
      setStatus({ phase: "error", progress: 0, message: "Processing failed." });
    }
  });
});
