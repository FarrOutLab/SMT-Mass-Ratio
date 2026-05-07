/**
 * SMT Mass-Ratio Widget — Standalone Embeddable Version
 *
 * Interactive widget for exploring mass-ratio predictions from a simple
 * stable mass transfer (SMT) model of binary black hole formation.
 *
 * Based on: van Son et al. (2022), ApJ 940, 184
 * Repository: https://github.com/FarrOutLab/SMT-Mass-Ratio
 *
 * Usage:
 *   <div id="smt-widget"></div>
 *   <script src="smt-widget.js"></script>
 *   <script>SMTWidget.init("smt-widget");</script>
 */

// eslint-disable-next-line no-var
var SMTWidget = (function () {
  "use strict";

  // ════════════════════════════════════════════════════════════════
  //  PHYSICS ENGINE  (ported 1:1 from SMTWidget.tsx / notebook)
  // ════════════════════════════════════════════════════════════════

  var N_Q = 2000;
  var N_BETA = 2000;
  var Q_GRID_MIN = 0.01;
  var Q_GRID_MAX = 8;
  var N_PLOT = 500;

  function linspace(start, stop, n) {
    var arr = new Float64Array(n);
    var step = (stop - start) / (n - 1);
    for (var i = 0; i < n; i++) arr[i] = start + i * step;
    return arr;
  }

  function trapezoid(y, x) {
    var sum = 0;
    for (var i = 1; i < x.length; i++) {
      sum += 0.5 * (y[i] + y[i - 1]) * (x[i] - x[i - 1]);
    }
    return sum;
  }

  function zetaRLSingle(beta, q) {
    var part1 =
      -2 * (1 - beta * q - (1 - beta) * (q + 0.5) * (q / (q + 1)));
    var q13 = Math.cbrt(q);
    var q23 = q13 * q13;
    var A = q13 / 3;
    var B = 2 / q13;
    var C =
      (1.2 * q13 + 1 / (1 + q13)) / (0.6 * q23 + Math.log(1 + q13));
    var part2 = A * (B - C);
    var part3 = 1 + beta * q;
    return part1 + part2 * part3;
  }

  // Pre-compute grids once
  var qsGrid = linspace(Q_GRID_MIN, Q_GRID_MAX, N_Q);
  var betasGrid = linspace(0, 1, N_BETA);
  var qobsPlot = linspace(0.01, 1, N_PLOT);
  var qbhPlot = linspace(0.01, 2, N_PLOT);

  var zetaRLGrid = [];
  for (var bi = 0; bi < betasGrid.length; bi++) {
    var row = new Float64Array(qsGrid.length);
    for (var qj = 0; qj < qsGrid.length; qj++) {
      row[qj] = zetaRLSingle(betasGrid[bi], qsGrid[qj]);
    }
    zetaRLGrid.push(row);
  }

  function qCrit(zetaEff, beta) {
    var bestBetaIdx = 0;
    var bestBetaDist = Math.abs(betasGrid[0] - beta);
    for (var i = 1; i < betasGrid.length; i++) {
      var d = Math.abs(betasGrid[i] - beta);
      if (d < bestBetaDist) {
        bestBetaDist = d;
        bestBetaIdx = i;
      }
    }
    var r = zetaRLGrid[bestBetaIdx];
    var bestQIdx = 0;
    var bestQDist = Math.abs(r[0] - zetaEff);
    for (var j = 1; j < r.length; j++) {
      var dj = Math.abs(r[j] - zetaEff);
      if (dj < bestQDist) {
        bestQDist = dj;
        bestQIdx = j;
      }
    }
    return qsGrid[bestQIdx];
  }

  function qZamsUniformPdf(qzams, qzamsMin, qzamsMax) {
    var qmax = qzamsMax > 1 ? 1 : qzamsMax;
    if (qzams < qzamsMin || qzams > qzamsMax || qzams > 1) return 0;
    var denom = qmax - qzamsMin;
    return denom > 0 ? 1 / denom : 0;
  }

  function qToQzams(q, fsn_a, fsn_b, beta, fcore) {
    return (q * (1 - fsn_a)) / (1 - fsn_b) - beta * (1 - fcore);
  }

  function pQobsFromZams(qArr, beta, fcore, fsn_a, fsn_b, qzamsMin, qzamsMax) {
    var dqNoMR = (1 - fsn_a) / (1 - fsn_b);
    var nPts = 1000;
    var qObsNorm = linspace(0.01, 1, nPts);
    var normVals = new Float64Array(nPts);
    for (var k = 0; k < nPts; k++) {
      var qo = qObsNorm[k];
      var qzNormNoMR = qToQzams(qo, fsn_a, fsn_b, beta, fcore);
      var pdfNormNoMR = qZamsUniformPdf(qzNormNoMR, qzamsMin, qzamsMax) * dqNoMR;
      var qzNormMR = qToQzams(1 / qo, fsn_a, fsn_b, beta, fcore);
      var dqMR_k = dqNoMR / (qo * qo);
      var pdfNormMR = qZamsUniformPdf(qzNormMR, qzamsMin, qzamsMax) * dqMR_k;
      normVals[k] = pdfNormNoMR + pdfNormMR;
    }
    var norm = trapezoid(normVals, qObsNorm);

    var pdf = new Float64Array(qArr.length);
    for (var i = 0; i < qArr.length; i++) {
      var q = qArr[i];
      var qzamsNoMR = qToQzams(q, fsn_a, fsn_b, beta, fcore);
      var pdfNoMR = qZamsUniformPdf(qzamsNoMR, qzamsMin, qzamsMax) * dqNoMR;
      var qzamsMR = qToQzams(1 / q, fsn_a, fsn_b, beta, fcore);
      var dqMR = dqNoMR / (q * q);
      var pdfMR = qZamsUniformPdf(qzamsMR, qzamsMin, qzamsMax) * dqMR;
      pdf[i] = norm > 0 ? (pdfNoMR + pdfMR) / norm : 0;
    }
    return pdf;
  }

  function pQbhFromZams(qArr, beta, fcore, fsn_a, fsn_b, qzamsMin, qzamsMax) {
    var dq = (1 - fsn_a) / (1 - fsn_b);
    var nPts = 1000;
    var qzNorm = linspace(0.01, 1, nPts);
    var qbhNorm = new Float64Array(nPts);
    var normVals = new Float64Array(nPts);
    for (var k = 0; k < nPts; k++) {
      qbhNorm[k] = ((1 - fsn_b) / (1 - fsn_a)) * (qzNorm[k] + beta * (1 - fcore));
      normVals[k] = qZamsUniformPdf(qzNorm[k], qzamsMin, qzamsMax) * dq;
    }
    var norm = trapezoid(normVals, qbhNorm);

    var pdf = new Float64Array(qArr.length);
    for (var i = 0; i < qArr.length; i++) {
      var qzams = qToQzams(qArr[i], fsn_a, fsn_b, beta, fcore);
      var pdfVal = qZamsUniformPdf(qzams, qzamsMin, qzamsMax) * dq;
      pdf[i] = norm > 0 ? pdfVal / norm : 0;
    }
    return pdf;
  }

  // ════════════════════════════════════════════════════════════════
  //  SLIDER DEFINITIONS
  // ════════════════════════════════════════════════════════════════

  var sliders = [
    { id: "fsn-a", label: "f_{\\text{SN,a}}", min: 0, max: 1, step: 0.01, initial: 0.2, decimals: 2 },
    { id: "fsn-b", label: "f_{\\text{SN,b}}", min: 0, max: 1, step: 0.01, initial: 0.2, decimals: 2 },
    { id: "fcore", label: "f_{\\text{core}}", min: 0.01, max: 1, step: 0.01, initial: 0.34, decimals: 2 },
    { id: "beta", label: "\\beta", min: 0, max: 1, step: 0.01, initial: 0.5, decimals: 2 },
    { id: "zeta", label: "\\zeta_{\\text{eff}}", min: 1, max: 10, step: 0.1, initial: 6.5, decimals: 1 },
  ];

  // ════════════════════════════════════════════════════════════════
  //  STYLES (injected once)
  // ════════════════════════════════════════════════════════════════

  var STYLE_ID = "smt-widget-styles";

  function injectStyles() {
    if (document.getElementById(STYLE_ID)) return;
    var style = document.createElement("style");
    style.id = STYLE_ID;
    style.textContent = [
      /* Theme toggle */
      ".smt-theme-toggle{display:flex;justify-content:flex-end;margin-bottom:1rem}",
      ".smt-theme-btn{display:flex;align-items:center;gap:.5rem;font-size:.85rem;color:#a1a1aa;cursor:pointer;user-select:none;background:rgba(255,255,255,.05);border:1px solid rgba(255,255,255,.1);border-radius:9999px;padding:.4rem 1rem;font-family:inherit;transition:background .2s,border-color .2s}",
      ".smt-theme-btn:hover{background:rgba(255,255,255,.1);border-color:rgba(255,255,255,.2)}",

      /* Widget container */
      ".smt-widget{display:flex;gap:2rem;align-items:stretch;min-height:480px;width:100%;transition:background .3s,border-color .3s;font-family:'Inter','Segoe UI',system-ui,-apple-system,sans-serif;box-sizing:border-box}",
      ".smt-widget *,.smt-widget *::before,.smt-widget *::after{box-sizing:border-box}",

      /* Controls panel */
      ".smt-controls{display:flex;flex-direction:column;gap:1.25rem;min-width:240px;max-width:280px;justify-content:center;padding:2rem 1.5rem;border:1px solid rgba(255,100,50,.1);border-radius:1rem;background:rgba(255,100,50,.03);backdrop-filter:blur(12px);transition:background .3s,border-color .3s,color .3s}",

      ".smt-slider-group{display:flex;flex-direction:column;gap:6px}",

      ".smt-label{display:flex;justify-content:space-between;align-items:center;font-size:.95rem;color:#d4d4d8;transition:color .3s}",
      ".smt-value{font-family:'Geist Mono','Fira Code','SF Mono',monospace;color:#60a5fa;min-width:3em;text-align:right;font-size:.9rem;transition:color .3s}",

      /* Range slider */
      ".smt-slider{-webkit-appearance:none;appearance:none;width:100%;height:6px;border-radius:3px;background:linear-gradient(90deg,#27272a,#3f3f46);outline:none;transition:background .2s}",
      ".smt-slider:hover{background:linear-gradient(90deg,#3f3f46,#52525b)}",
      ".smt-slider::-webkit-slider-thumb{-webkit-appearance:none;appearance:none;width:18px;height:18px;border-radius:50%;background:radial-gradient(circle,#60a5fa 60%,#3b82f6 100%);cursor:pointer;border:2px solid rgba(96,165,250,.3);box-shadow:0 0 8px rgba(96,165,250,.4);transition:transform .15s,box-shadow .15s}",
      ".smt-slider::-webkit-slider-thumb:hover{transform:scale(1.15);box-shadow:0 0 14px rgba(96,165,250,.6)}",
      ".smt-slider::-moz-range-thumb{width:18px;height:18px;border-radius:50%;background:radial-gradient(circle,#60a5fa 60%,#3b82f6 100%);cursor:pointer;border:2px solid rgba(96,165,250,.3);box-shadow:0 0 8px rgba(96,165,250,.4)}",

      /* Checkbox toggle */
      ".smt-toggle-group{display:flex;align-items:center;gap:.6rem;font-size:.95rem;color:#d4d4d8;margin-top:.5rem;padding-top:.75rem;border-top:1px solid rgba(255,255,255,.06);transition:color .3s,border-color .3s}",
      ".smt-checkbox{accent-color:#F46A35;width:18px;height:18px;cursor:pointer}",

      /* Plot area */
      ".smt-plot-wrap{flex:1;min-width:0;border:1px solid rgba(255,100,50,.1);border-radius:1rem;background:rgba(255,100,50,.03);padding:1rem;display:flex;align-items:center;transition:background .3s,border-color .3s}",
      ".smt-plot{width:100%;height:100%;min-height:420px}",

      /* ── Light theme ─────────────────────────────────────────── */
      ".smt-widget.smt-light{border-radius:1rem;padding:1.5rem;background:#fff;box-shadow:0 2px 12px rgba(0,0,0,.08)}",
      ".smt-light .smt-controls{background:#f4f4f5;border-color:#e4e4e7}",
      ".smt-light .smt-label{color:#27272a}",
      ".smt-light .smt-value{color:#2563eb}",
      ".smt-light .smt-slider{background:linear-gradient(90deg,#d4d4d8,#e4e4e7)}",
      ".smt-light .smt-slider:hover{background:linear-gradient(90deg,#a1a1aa,#d4d4d8)}",
      ".smt-light .smt-slider::-webkit-slider-thumb{background:radial-gradient(circle,#2563eb 60%,#1d4ed8 100%);border-color:rgba(37,99,235,.3);box-shadow:0 0 6px rgba(37,99,235,.3)}",
      ".smt-light .smt-slider::-moz-range-thumb{background:radial-gradient(circle,#2563eb 60%,#1d4ed8 100%);border-color:rgba(37,99,235,.3);box-shadow:0 0 6px rgba(37,99,235,.3)}",
      ".smt-light .smt-toggle-group{color:#3f3f46;border-top-color:#e4e4e7}",
      ".smt-light .smt-plot-wrap{background:#fafafa;border-color:#e4e4e7}",
      ".smt-light svg{color:#27272a}",
      ".smt-light .smt-theme-btn{color:#52525b;background:rgba(0,0,0,.04);border-color:rgba(0,0,0,.1)}",
      ".smt-light .smt-theme-btn:hover{background:rgba(0,0,0,.08);border-color:rgba(0,0,0,.15)}",

      /* ── Responsive ──────────────────────────────────────────── */
      "@media(max-width:768px){.smt-widget{flex-direction:column}.smt-controls{max-width:100%;min-width:unset}.smt-plot{min-height:320px}}",
    ].join("\n");
    document.head.appendChild(style);
  }

  // ════════════════════════════════════════════════════════════════
  //  CDN LOADER HELPERS
  // ════════════════════════════════════════════════════════════════

  function loadScript(src, onLoad) {
    var s = document.createElement("script");
    s.src = src;
    s.async = true;
    s.onload = onLoad;
    document.head.appendChild(s);
  }

  function ensurePlotly(cb) {
    if (window.Plotly) return cb();
    loadScript("https://cdn.plot.ly/plotly-2.35.2.min.js", cb);
  }

  function ensureMathJax(cb) {
    if (window.MathJax && window.MathJax.typesetPromise) return cb();
    // Configure MathJax before loading
    window.MathJax = {
      tex: { inlineMath: [["$", "$"], ["\\(", "\\)"]] },
      svg: { fontCache: "global" },
      startup: {
        ready: function () {
          window.MathJax.startup.defaultReady();
          cb();
        },
      },
    };
    loadScript("https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-svg.js", function () {
      // MathJax calls our startup.ready callback above
    });
  }

  // ════════════════════════════════════════════════════════════════
  //  DOM BUILDER & RENDERING
  // ════════════════════════════════════════════════════════════════

  function el(tag, attrs, children) {
    var e = document.createElement(tag);
    if (attrs) {
      Object.keys(attrs).forEach(function (k) {
        if (k === "className") e.className = attrs[k];
        else if (k.indexOf("on") === 0) e.addEventListener(k.slice(2).toLowerCase(), attrs[k]);
        else e.setAttribute(k, attrs[k]);
      });
    }
    if (children) {
      if (typeof children === "string") e.textContent = children;
      else if (Array.isArray(children)) children.forEach(function (c) { if (c) e.appendChild(c); });
      else e.appendChild(children);
    }
    return e;
  }

  function init(containerId, options) {
    options = options || {};
    var container = document.getElementById(containerId);
    if (!container) {
      console.error("[SMTWidget] Container #" + containerId + " not found.");
      return;
    }

    injectStyles();

    // State
    var values = {};
    sliders.forEach(function (s) { values[s.id] = s.initial; });
    var showQbh = false;
    var lightMode = !!options.light;

    // Build DOM structure
    var plotDiv = el("div", { className: "smt-plot" });
    var controlsDiv = el("div", { className: "smt-controls" });
    var plotWrap = el("div", { className: "smt-plot-wrap" }, [plotDiv]);
    var widgetDiv = el("div", { className: "smt-widget" + (lightMode ? " smt-light" : "") }, [controlsDiv, plotWrap]);

    // Theme toggle
    var themeIcon = el("span", null, lightMode ? "☀️" : "🌙");
    var themeLabel = el("span", null, lightMode ? "Light" : "Dark");
    var themeBtn = el("button", {
      className: "smt-theme-btn",
      type: "button",
      "aria-label": "Toggle light/dark mode",
    }, [themeIcon, themeLabel]);
    var themeWrap = el("div", { className: "smt-theme-toggle" }, [themeBtn]);

    themeBtn.addEventListener("click", function () {
      lightMode = !lightMode;
      themeIcon.textContent = lightMode ? "☀️" : "🌙";
      themeLabel.textContent = lightMode ? "Light" : "Dark";
      if (lightMode) widgetDiv.classList.add("smt-light");
      else widgetDiv.classList.remove("smt-light");
      updatePlot();
    });

    container.appendChild(themeWrap);
    container.appendChild(widgetDiv);

    // Build sliders
    var valueSpans = {};
    sliders.forEach(function (s) {
      var valSpan = el("span", { className: "smt-value" }, values[s.id].toFixed(s.decimals));
      valueSpans[s.id] = valSpan;

      var labelEl = el("label", { className: "smt-label" }, [
        el("span", null, "$" + s.label + "$"),
        valSpan,
      ]);

      var input = el("input", {
        type: "range",
        className: "smt-slider",
        min: String(s.min),
        max: String(s.max),
        step: String(s.step),
        value: String(s.initial),
      });
      input.addEventListener("input", function () {
        var v = parseFloat(this.value);
        values[s.id] = v;
        valSpan.textContent = v.toFixed(s.decimals);
        updatePlot();
      });

      var group = el("div", { className: "smt-slider-group" }, [labelEl, input]);
      controlsDiv.appendChild(group);
    });

    // Show q_BBH toggle
    var checkboxId = "smt-toggle-qbh-" + containerId;
    var checkbox = el("input", {
      type: "checkbox",
      id: checkboxId,
      className: "smt-checkbox",
    });
    checkbox.addEventListener("change", function () {
      showQbh = this.checked;
      updatePlot();
    });
    var toggleLabel = el("label", { "for": checkboxId });
    toggleLabel.innerHTML = "Show $p(q_{\\text{BBH}})$";
    var toggleGroup = el("div", { className: "smt-toggle-group" }, [checkbox, toggleLabel]);
    controlsDiv.appendChild(toggleGroup);

    // ── Plot update function ─────────────────────────────────
    var plotlyReady = false;

    function updatePlot() {
      if (!plotlyReady) return;

      var fsn_a = values["fsn-a"];
      var fsn_b = values["fsn-b"];
      var fcore = values["fcore"];
      var beta = values["beta"];
      var zeta = values["zeta"];

      var qcrit_beta = qCrit(zeta, beta);
      var qcrit_0 = qCrit(zeta, 0);
      var qzamsMin = 1 / qcrit_beta;
      var qzamsMax = fcore * (1 - fsn_a) * qcrit_0 - beta * (1 - fcore);

      var pQobs = pQobsFromZams(qobsPlot, beta, fcore, fsn_a, fsn_b, qzamsMin, qzamsMax);

      // Theme-aware colors
      var axisColor = lightMode ? "#52525b" : "#71717a";
      var axisTitleColor = lightMode ? "#3f3f46" : "#a1a1aa";
      var gridColor = lightMode ? "rgba(0,0,0,0.08)" : "rgba(255,255,255,0.06)";
      var zerolineColor = lightMode ? "rgba(0,0,0,0.15)" : "rgba(255,255,255,0.1)";
      var fontColor = lightMode ? "#3f3f46" : "#a1a1aa";
      var shapeLine = lightMode ? "rgba(0,0,0,0.2)" : "rgba(255,255,255,0.25)";
      var traceBlue = lightMode ? "#2563eb" : "#60a5fa";
      var traceOrange = lightMode ? "#ea580c" : "#F46A35";

      var traces = [
        {
          x: Array.from(qobsPlot),
          y: Array.from(pQobs),
          mode: "lines",
          name: "p(q<sub>obs</sub>)",
          line: { color: traceBlue, width: 2.5 },
        },
      ];

      if (showQbh) {
        var pQbh = pQbhFromZams(qbhPlot, beta, fcore, fsn_a, fsn_b, qzamsMin, qzamsMax);
        traces.push({
          x: Array.from(qbhPlot),
          y: Array.from(pQbh),
          mode: "lines",
          name: "p(q<sub>BBH</sub>)",
          line: { color: traceOrange, width: 2.5 },
        });
      }

      var shapes = [
        {
          type: "line",
          x0: 1, x1: 1,
          y0: 0, y1: 1,
          yref: "paper",
          line: { color: shapeLine, width: 1.5, dash: "dash" },
        },
      ];

      var layout = {
        xaxis: {
          title: { text: "q", font: { size: 18, color: axisTitleColor } },
          range: [0, showQbh ? 2 : 1],
          color: axisColor,
          gridcolor: gridColor,
          zerolinecolor: zerolineColor,
        },
        yaxis: {
          title: { text: "p(q)", font: { size: 18, color: axisTitleColor } },
          color: axisColor,
          gridcolor: gridColor,
          zerolinecolor: zerolineColor,
          rangemode: "tozero",
        },
        paper_bgcolor: "rgba(0,0,0,0)",
        plot_bgcolor: "rgba(0,0,0,0)",
        font: { color: fontColor, family: "Inter, system-ui, sans-serif" },
        margin: { l: 60, r: 20, t: 10, b: 50 },
        shapes: shapes,
        legend: {
          x: 0.05,
          y: 0.95,
          font: { size: 14 },
          bgcolor: "rgba(0,0,0,0)",
        },
        showlegend: true,
      };

      var config = { responsive: true, displayModeBar: false };

      window.Plotly.react(plotDiv, traces, layout, config);
    }

    // ── Load dependencies & boot ─────────────────────────────
    ensurePlotly(function () {
      plotlyReady = true;
      updatePlot();
    });

    ensureMathJax(function () {
      if (window.MathJax && window.MathJax.typesetPromise) {
        window.MathJax.typesetPromise([controlsDiv]).catch(function () { });
      }
    });
  }

  // ── Public API ─────────────────────────────────────────────
  return { init: init };
})();
