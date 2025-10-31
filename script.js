// ========================================================================
// === Constants & Helpers ================================================
// ========================================================================

const gamma = 1.4;
const g = gamma; // Shorthand
const gp1 = gamma + 1;
const gm1 = gamma - 1;

// Convert degrees to radians
function deg2rad(deg) {
  return (deg * Math.PI) / 180;
}

// Convert radians to degrees
function rad2deg(rad) {
  return (rad * 180) / Math.PI;
}

// Generic numerical solver
function newtonRaphson(func, deriv, guess, target, maxIter = 50, tol = 1e-9) {
  let x = guess;
  for (let i = 0; i < maxIter; i++) {
    const f = func(x, gamma) - target;
    if (Math.abs(f) < tol) return x; // Converged
    const df = deriv(x, gamma);
    if (Math.abs(df) < 1e-15) break; // Avoid division by zero
    x = x - f / df;
    if (x < 0) x = tol; // Prevent negative Mach numbers
  }
  return NaN; // Failed to converge
}

// Helper to update a DOM element's text
function updateText(id, value, precision = 4) {
  const el = document.getElementById(id);
  if (el) {
    el.textContent = isNaN(value) ? "–" : value.toFixed(precision);
  }
}

// Helper to update a DOM element's visibility
function setVisible(id, visible) {
  const el = document.getElementById(id);
  if (el) {
    el.style.display = visible ? "block" : "none";
  }
}

// ========================================================================
// === Gas Dynamics Models (The "Math") ===================================
// ========================================================================

// --- 1. Isentropic Flow ---
const isentropic = {
  // Ratios from M
  T_T0: (M, g = gamma) => (1 + (g - 1) / 2 * M ** 2) ** -1,
  p_p0: (M, g = gamma) => (1 + (g - 1) / 2 * M ** 2) ** (-g / (g - 1)),
  rho_rho0: (M, g = gamma) => (1 + (g - 1) / 2 * M ** 2) ** (-1 / (g - 1)),
  A_Astar: (M, g = gamma) => (1 / M) * (((1 + (g - 1) / 2 * M ** 2)) / ((g + 1) / 2)) ** ((g + 1) / (2 * (g - 1))),

  // M from Ratios
  M_from_T: (T_T0, g = gamma) => Math.sqrt((1 / T_T0 - 1) * 2 / (g - 1)),
  M_from_p: (p_p0, g = gamma) => Math.sqrt(((p_p0) ** (-(g - 1) / g) - 1) * 2 / (g - 1)),
  M_from_rho: (rho_rho0, g = gamma) => Math.sqrt(((rho_rho0) ** -(g - 1) - 1) * 2 / (g - 1)),

  // A/A* inversion (needs solver)
  A_Astar_func: (M, g = gamma) => isentropic.A_Astar(M, g),
  A_Astar_deriv: (M, g = gamma) => {
    const term1 = (1 + (g - 1) / 2 * M ** 2) / ((g + 1) / 2);
    const exponent = (g + 1) / (2 * (g - 1));
    return isentropic.A_Astar(M, g) / M * ((1 - M ** 2) / (1 + (g - 1) / 2 * M ** 2));
  },
  M_from_A: (A_Astar, g = gamma, supersonic = false) => {
    if (A_Astar < 1) return NaN;
    if (A_Astar === 1) return 1;
    const guess = supersonic ? (1 + 2 * A_Astar) : (1 - 0.5 * A_Astar);
    return newtonRaphson(isentropic.A_Astar_func, isentropic.A_Astar_deriv, guess, A_Astar);
  }
};

// --- 2. Normal Shock ---
function normalShockRatios(M1, g = gamma) {
  if (M1 <= 1) return {};
  const M1_sq = M1 ** 2;
  const M2 = Math.sqrt((1 + (g - 1) / 2 * M1_sq) / (g * M1_sq - (g - 1) / 2));
  const p2p1 = 1 + 2 * g / (g + 1) * (M1_sq - 1);
  const rho2rho1 = (g + 1) * M1_sq / ((g - 1) * M1_sq + 2);
  const T2T1 = p2p1 / rho2rho1;
  
  // Stagnation pressure ratio (from p02/p01 = (p02/p2) * (p2/p1) * (p1/p01))
  const p02p01 = isentropic.p_p0(M2, g) * p2p1 * (1 / isentropic.p_p0(M1, g));

  return { M2, p2p1, rho2rho1, T2T1, p02p01 };
}

// --- 3. Oblique Shock (θ–β–M relation) ---
function betaAnalytical(M1, thetaDeg, g = gamma, n = 0) {
  const theta = deg2rad(thetaDeg);
  const mu = Math.asin(1 / M1);
  const c = Math.tan(mu) ** 2;
  const a = ((g - 1) / 2 + ((g + 1) / 2) * c) * Math.tan(theta);
  const b = ((g + 1) / 2 + ((g + 3) / 2) * c) * Math.tan(theta);
  const num = 4 * (1 - 3 * a * b) ** 3;
  const den = (27 * a * a * c + 9 * a * b - 2) ** 2;
  const d = Math.sqrt(num / den - 1);
  const term1 = (b + 9 * a * c) / (2 * (1 - 3 * a * b));
  const term2 = (d * (27 * a * a * c + 9 * a * b - 2)) / (6 * a * (1 - 3 * a * b));
  const angle = (n * Math.PI) / 3 + (1 / 3) * Math.atan(1 / d);
  return Math.atan(term1 - term2 * Math.tan(angle));
}

function obliqueShockRatios(M1, beta, theta, g = gamma) {
  const Mn1 = M1 * Math.sin(beta);
  if (Mn1 <= 1) return {}; // Not a shock
  
  const { M2: Mn2, p2p1, rho2rho1, T2T1, p02p01 } = normalShockRatios(Mn1, g);
  if (isNaN(Mn2)) return {};
  
  const M2 = Mn2 / Math.sin(beta - theta);
  return { Mn1, Mn2, M2, p2p1, rho2rho1, T2T1, p02p01 };
}

// --- 4. Prandtl-Meyer Expansion ---
const prandtlMeyer = {
  nu: (M, g = gamma) => {
    if (M < 1) return NaN;
    const g_ratio = Math.sqrt((g + 1) / (g - 1));
    const M_sq_term = Math.sqrt((g - 1) / (g + 1) * (M ** 2 - 1));
    const M_term = Math.sqrt(M ** 2 - 1);
    return rad2deg(g_ratio * Math.atan(M_sq_term) - Math.atan(M_term));
  },

  nu_rad: (M, g = gamma) => deg2rad(prandtlMeyer.nu(M, g)), // For solver

  nu_deriv: (M, g = gamma) => {
    if (M <= 1) return 0;
    const M_sq = M ** 2;
    return (Math.sqrt(M_sq - 1) / (1 + (g - 1) / 2 * (M_sq - 1))) * (1 / M);
  },

  M_from_nu: (nu_deg, g = gamma) => {
    if (nu_deg < 0) return NaN;
    // Max PM angle
    const nu_max = rad2deg(Math.sqrt((g + 1) / (g - 1)) * Math.PI / 2 - Math.PI / 2);
    if (nu_deg > nu_max) return Infinity;

    const nu_rad_target = deg2rad(nu_deg);
    // Guess M based on nu (simple empirical guess)
    const guess = 1.5 + (nu_deg / 130) * 10;
    return newtonRaphson(prandtlMeyer.nu_rad, prandtlMeyer.nu_deriv, guess, nu_rad_target);
  }
};


// ========================================================================
// === Diagram Updaters ===================================================
// ========================================================================

// SVG Arc Path helper (from original)
function arcPath(cx, cy, r, a1, a2) {
  const x1 = cx + r * Math.cos(a1);
  const y1 = cy - r * Math.sin(a1);
  const x2 = cx + r * Math.cos(a2);
  const y2 = cy - r * Math.sin(a2);
  const largeArcFlag = a2 - a1 <= Math.PI ? "0" : "1";
  return `M ${x1},${y1} A ${r},${r} 0 ${largeArcFlag} 0 ${x2},${y2}`;
}

// --- 1. Isentropic (No diagram) ---

// --- 2. Normal Shock Diagram ---
function updateNormalShockDiagram(M1, M2) {
  updateText("ns-M1-label", M1, 2);
  updateText("ns-M2-label", M2, 2);
}

// --- 3. Oblique Shock Diagram ---
function updateObliqueDiagram(thetaDeg, betaDeg) {
  const x0 = 280, y0 = 320;
  const thetaRad = deg2rad(thetaDeg);

  // Wedge and Theta (always draw)
  const Lbase = 220, Lwedge = 150;
  document.getElementById("ob-wedgeBase").setAttribute("x2", x0 + Lbase);
  document.getElementById("ob-wedgeRamp").setAttribute("x2", x0 + Lwedge * Math.cos(thetaRad));
  document.getElementById("ob-wedgeRamp").setAttribute("y2", y0 - Lwedge * Math.sin(thetaRad));
  document.getElementById("ob-thetaArc").setAttribute("d", arcPath(x0, y0, 40, 0, thetaRad));
  document.getElementById("ob-thetaLabel").setAttribute("x", x0 + 60 * Math.cos(thetaRad / 2));
  document.getElementById("ob-thetaLabel").setAttribute("y", y0 - 60 * Math.sin(thetaRad / 2));
  document.getElementById("ob-thetaLabel").textContent = `θ=${thetaDeg.toFixed(1)}°`;

  // Shock and Beta (only if valid)
  const isValid = !isNaN(betaDeg) && betaDeg > thetaDeg;
  setVisible("ob-shockLine", isValid);
  setVisible("ob-betaArc", isValid);
  setVisible("ob-betaLabel", isValid);

  if (isValid) {
    const Lshock = 220;
    const betaRad = deg2rad(betaDeg);
    document.getElementById("ob-shockLine").setAttribute("x2", x0 + Lshock * Math.cos(betaRad));
    document.getElementById("ob-shockLine").setAttribute("y2", y0 - Lshock * Math.sin(betaRad));
    document.getElementById("ob-betaArc").setAttribute("d", arcPath(x0, y0, 60, 0, betaRad));
    document.getElementById("ob-betaLabel").setAttribute("x", x0 + 85 * Math.cos(betaRad / 2));
    document.getElementById("ob-betaLabel").setAttribute("y", y0 - 85 * Math.sin(betaRad / 2));
    document.getElementById("ob-betaLabel").textContent = `β=${betaDeg.toFixed(1)}°`;
  }
}

// --- 4. Prandtl-Meyer Diagram ---
function updatePMDiagram(M1, thetaDeg) {
  const x0 = 280, y0 = 200;
  const thetaRad = deg2rad(thetaDeg);
  
  // Walls
  const Lwall = 220;
  document.getElementById("pm-wall1").setAttribute("x2", x0 + Lwall);
  document.getElementById("pm-wall2").setAttribute("x2", x0 + Lwall * Math.cos(thetaRad));
  document.getElementById("pm-wall2").setAttribute("y2", y0 + Lwall * Math.sin(thetaRad));

  // Theta Arc
  document.getElementById("pm-thetaArc").setAttribute("d", arcPath(x0, y0, 60, -thetaRad, 0));
  document.getElementById("pm-thetaLabel").setAttribute("x", x0 + 80 * Math.cos(thetaRad / 2));
  document.getElementById("pm-thetaLabel").setAttribute("y", y0 + 80 * Math.sin(thetaRad / 2));
  document.getElementById("pm-thetaLabel").textContent = `θ=${thetaDeg.toFixed(1)}°`;

  // Mach Waves (μ = asin(1/M))
  const mu1Rad = Math.asin(1 / M1);
  const mu1Deg = rad2deg(mu1Rad);
  const mu2Rad = Math.asin(1 / (prandtlMeyer.M_from_nu(prandtlMeyer.nu(M1) + thetaDeg) || M1));

  const Lfan = 240;
  document.getElementById("pm-fan1").setAttribute("x2", x0 + Lfan * Math.cos(mu1Rad));
  document.getElementById("pm-fan1").setAttribute("y2", y0 - Lfan * Math.sin(mu1Rad));
  document.getElementById("pm-fan3").setAttribute("x2", x0 + Lfan * Math.cos(thetaRad + mu2Rad));
  document.getElementById("pm-fan3").setAttribute("y2", y0 + Lfan * Math.sin(thetaRad + mu2Rad));

  // Mu1 Arc
  document.getElementById("pm-mu1Arc").setAttribute("d", arcPath(x0, y0, 40, 0, mu1Rad));
  document.getElementById("pm-mu1Label").setAttribute("x", x0 + 55 * Math.cos(mu1Rad / 2));
  document.getElementById("pm-mu1Label").setAttribute("y", y0 - 55 * Math.sin(mu1Rad / 2));
  document.getElementById("pm-mu1Label").textContent = `μ₁=${mu1Deg.toFixed(1)}°`;
}

// ========================================================================
// === Calculator Controllers =============================================
// ========================================================================

// --- 1. Isentropic Controller ---
function calculateIsentropic() {
  const inputType = document.getElementById("is-inputType").value;
  const inputValue = parseFloat(document.getElementById("is-inputValue").value);
  const warning = document.getElementById("is-warning");

  let M = NaN, T_T0, p_p0, rho_rho0, A_Astar;
  warning.textContent = "";

  try {
    switch (inputType) {
      case "M":
        M = inputValue;
        if (M <= 0) throw new Error("M must be > 0.");
        T_T0 = isentropic.T_T0(M);
        p_p0 = isentropic.p_p0(M);
        rho_rho0 = isentropic.rho_rho0(M);
        A_Astar = isentropic.A_Astar(M);
        break;
      case "T_T0":
        T_T0 = inputValue;
        if (T_T0 <= 0 || T_T0 > 1) throw new Error("T/T₀ must be in (0, 1].");
        M = isentropic.M_from_T(T_T0);
        p_p0 = isentropic.p_p0(M);
        rho_rho0 = isentropic.rho_rho0(M);
        A_Astar = isentropic.A_Astar(M);
        break;
      case "p_p0":
        p_p0 = inputValue;
        if (p_p0 <= 0 || p_p0 > 1) throw new Error("p/p₀ must be in (0, 1].");
        M = isentropic.M_from_p(p_p0);
        T_T0 = isentropic.T_T0(M);
        rho_rho0 = isentropic.rho_rho0(M);
        A_Astar = isentropic.A_Astar(M);
        break;
      case "rho_rho0":
        rho_rho0 = inputValue;
        if (rho_rho0 <= 0 || rho_rho0 > 1) throw new Error("ρ/ρ₀ must be in (0, 1].");
        M = isentropic.M_from_rho(rho_rho0);
        T_T0 = isentropic.T_T0(M);
        p_p0 = isentropic.p_p0(M);
        A_Astar = isentropic.A_Astar(M);
        break;
      case "A_Astar_sub":
      case "A_Astar_super":
        A_Astar = inputValue;
        if (A_Astar < 1) throw new Error("A/A* must be ≥ 1.");
        M = isentropic.M_from_A(A_Astar, gamma, inputType === "A_Astar_super");
        T_T0 = isentropic.T_T0(M);
        p_p0 = isentropic.p_p0(M);
        rho_rho0 = isentropic.rho_rho0(M);
        break;
    }
    if (isNaN(M)) throw new Error("Calculation failed. Check inputs.");
  } catch (err) {
    warning.textContent = err.message;
  }

  updateText("is-M", M);
  updateText("is-p_p0", p_p0);
  updateText("is-rho_rho0", rho_rho0);
  updateText("is-T_T0", T_T0);
  updateText("is-A_Astar", A_Astar);
}

// --- 2. Normal Shock Controller ---
function calculateNormalShock() {
  const M1 = parseFloat(document.getElementById("ns-M1").value);
  const warning = document.getElementById("ns-warning");

  if (isNaN(M1) || M1 <= 1) {
    warning.textContent = "M₁ must be > 1 for a normal shock.";
    ["ns-M2", "ns-p2p1", "ns-rho2rho1", "ns-T2T1", "ns-p02p01"].forEach(id => updateText(id, NaN));
    updateNormalShockDiagram(M1, NaN);
    return;
  }
  warning.textContent = "";
  const { M2, p2p1, rho2rho1, T2T1, p02p01 } = normalShockRatios(M1);

  updateText("ns-M2", M2);
  updateText("ns-p2p1", p2p1);
  updateText("ns-rho2rho1", rho2rho1);
  updateText("ns-T2T1", T2T1);
  updateText("ns-p02p01", p02p01);
  updateNormalShockDiagram(M1, M2);
}

// --- 3. Oblique Shock Controller (from original) ---
function calculateObliqueShock() {
  const M1 = parseFloat(document.getElementById("ob-M1").value);
  const thetaDeg = parseFloat(document.getElementById("ob-theta").value);
  const warning = document.getElementById("ob-warning");

  if (isNaN(M1) || isNaN(thetaDeg) || M1 <= 1) {
    warning.textContent = "M₁ must be > 1 for an attached oblique shock.";
    updateDiagram(thetaDeg, NaN);
    return;
  }

  // Calculate weak shock (n=0)
  const betaRad = betaAnalytical(M1, thetaDeg, gamma, 0);
  const betaDeg = rad2deg(betaRad);
  
  if (isNaN(betaRad) || betaDeg <= thetaDeg || betaDeg >= 90) {
    warning.textContent = "No attached weak oblique shock for this (M₁, θ).";
    ["ob-beta", "ob-M1n", "ob-M2n", "ob-M2", "ob-p2p1", "ob-rho2rho1", "ob-T2T1", "ob-p02p01"].forEach(id => updateText(id, NaN));
    updateObliqueDiagram(thetaDeg, NaN);
    return;
  }
  
  warning.textContent = "";
  const thetaRad = deg2rad(thetaDeg);
  const { Mn1, Mn2, M2, p2p1, rho2rho1, T2T1, p02p01 } = obliqueShockRatios(M1, betaRad, thetaRad);

  updateText("ob-thetaOut", thetaDeg, 1);
  updateText("ob-beta", betaDeg, 2);
  updateText("ob-M1n", Mn1);
  updateText("ob-M2n", Mn2);
  updateText("ob-M2", M2);
  updateText("ob-p2p1", p2p1);
  updateText("ob-rho2rho1", rho2rho1);
  updateText("ob-T2T1", T2T1);
  updateText("ob-p02p01", p02p01);

  updateObliqueDiagram(thetaDeg, betaDeg);
}

// --- 4. Prandtl-Meyer Controller ---
function calculatePrandtlMeyer() {
  const M1 = parseFloat(document.getElementById("pm-M1").value);
  const thetaDeg = parseFloat(document.getElementById("pm-theta").value);
  const warning = document.getElementById("pm-warning");
  
  if (isNaN(M1) || M1 <= 1) {
    warning.textContent = "M₁ must be > 1 for an expansion.";
    return;
  }
  
  const nu1 = prandtlMeyer.nu(M1);
  const nu2 = nu1 + thetaDeg;

  const nu_max = rad2deg(Math.sqrt((g + 1) / (g - 1)) * Math.PI / 2 - Math.PI / 2);
  if (nu2 > nu_max) {
    warning.textContent = `Expansion angle is too large. Max ν is ${nu_max.toFixed(1)}°.`;
    return;
  }
  
  warning.textContent = "";
  const M2 = prandtlMeyer.M_from_nu(nu2);
  
  // Get ratios using isentropic relations (expansion is isentropic)
  const p01_p1 = 1 / isentropic.p_p0(M1);
  const p02_p2 = 1 / isentropic.p_p0(M2);
  const p2p1 = p01_p1 / p02_p2; // p01 = p02

  const T01_T1 = 1 / isentropic.T_T0(M1);
  const T02_T2 = 1 / isentropic.T_T0(M2);
  const T2T1 = T01_T1 / T02_T2; // T01 = T02
  
  const rho01_rho1 = 1 / isentropic.rho_rho0(M1);
  const rho02_rho2 = 1 / isentropic.rho_rho0(M2);
  const rho2rho1 = rho01_rho1 / rho02_rho2; // rho01 = rho02

  updateText("pm-nu1", nu1, 2);
  updateText("pm-nu2", nu2, 2);
  updateText("pm-M2", M2);
  updateText("pm-p2p1", p2p1);
  updateText("pm-rho2rho1", rho2rho1);
  updateText("pm-T2T1", T2T1);
  
  updatePMDiagram(M1, thetaDeg);
}

// ========================================================================
// === Initialization =====================================================
// ========================================================================

document.addEventListener("DOMContentLoaded", () => {
  // --- Main Tab Navigation ---
  const navButtons = document.querySelectorAll(".nav-btn");
  const panels = document.querySelectorAll(".calculator-panel");

  navButtons.forEach(btn => {
    btn.addEventListener("click", () => {
      // Deactivate all
      navButtons.forEach(b => b.classList.remove("active"));
      panels.forEach(p => p.classList.remove("active"));

      // Activate clicked
      btn.classList.add("active");
      const targetPanelId = btn.dataset.calculator;
      document.getElementById(`${targetPanelId}-calculator`).classList.add("active");
    });
  });

  // --- Setup Calculator Event Listeners ---
  document.getElementById("is-computeBtn").addEventListener("click", calculateIsentropic);
  document.getElementById("is-inputType").addEventListener("change", calculateIsentropic);
  
  document.getElementById("ns-computeBtn").addEventListener("click", calculateNormalShock);
  
  document.getElementById("ob-computeBtn").addEventListener("click", calculateObliqueShock);
  
  document.getElementById("pm-computeBtn").addEventListener("click", calculatePrandtlMeyer);

  // --- Run all calculators on initial load ---
  calculateIsentropic();
  calculateNormalShock();
  calculateObliqueShock();
  calculatePrandtlMeyer();
});