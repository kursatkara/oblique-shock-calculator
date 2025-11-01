// ========================================================================
// === Constants & Helpers ================================================
// ========================================================================

function getGamma() {
  return parseFloat(document.getElementById("gamma").value) || 1.4;
}

function deg2rad(deg) { return (deg * Math.PI) / 180; }
function rad2deg(rad) { return (rad * 180) / Math.PI; }

function newtonRaphson(func, deriv, guess, target, g, maxIter = 50, tol = 1e-9) {
  let x = guess;
  for (let i = 0; i < maxIter; i++) {
    const f = func(x, g) - target;
    if (Math.abs(f) < tol) return x;
    const df = deriv(x, g);
    if (Math.abs(df) < 1e-15) break;
    x = x - f / df;
    if (x < 0) x = tol;
  }
  return NaN;
}

function updateText(id, value, precision = 4) {
  const el = document.getElementById(id);
  if (el) {
    const text = isNaN(value) ? "–" : value.toFixed(precision);
    if (el.tagName === "INPUT") {
      el.value = text;
    } else {
      el.textContent = text;
    }
  }
}

function setVisible(id, visible) {
  const el = document.getElementById(id);
  if (el) { el.style.display = visible ? "block" : "none"; }
}

function arcPath(cx, cy, r, a1, a2) {
  const x1 = cx + r * Math.cos(a1);
  const y1 = cy - r * Math.sin(a1);
  const x2 = cx + r * Math.cos(a2);
  const y2 = cy - r * Math.sin(a2);
  return `M ${x1},${y1} A ${r},${r} 0 0 0 ${x2},${y2}`;
}

// Helper for flashing new results
function flashResults(panelId) {
    const results = document.querySelectorAll(`#${panelId}-calculator .results .row`);
    results.forEach(el => {
        el.classList.remove('is-new-result-row');
        void el.offsetWidth; // Trigger reflow
        el.classList.add('is-new-result-row');
        setTimeout(() => {
            el.classList.remove('is-new-result-row');
        }, 800); // Must match animation duration
    });
}

// ========================================================================
// === Gas Dynamics Models (The "Math") ===================================
// ========================================================================

const isentropic = {
  T0_T: (M, g) => (1 + (g - 1) / 2 * M ** 2),
  p0_p: (M, g) => (1 + (g - 1) / 2 * M ** 2) ** (g / (g - 1)),
  rho0_rho: (M, g) => (1 + (g - 1) / 2 * M ** 2) ** (1 / (g - 1)),
  A_Astar: (M, g) => (1 / M) * (((1 + (g - 1) / 2 * M ** 2)) / ((g + 1) / 2)) ** ((g + 1) / (2 * (g - 1))),
  M_from_T: (T0_T, g) => Math.sqrt((T0_T - 1) * 2 / (g - 1)),
  M_from_p: (p0_p, g) => Math.sqrt(((p0_p) ** ((g - 1) / g) - 1) * 2 / (g - 1)),
  M_from_rho: (rho0_rho, g) => Math.sqrt(((rho0_rho) ** (g - 1) - 1) * 2 / (g - 1)),
  A_Astar_func: (M, g) => isentropic.A_Astar(M, g),
  A_Astar_deriv: (M, g) => {
    return isentropic.A_Astar(M, g) * (M ** 2 - 1) / (M * (1 + (g - 1) / 2 * M ** 2));
  },
  M_from_A: (A_Astar, g, supersonic = false) => {
    if (A_Astar < 1) return NaN;
    if (A_Astar === 1) return 1;
    // *** FIX: Improved guess for subsonic branch ***
    const guess = supersonic ? (1 + 0.5 * A_Astar) : (1 / A_Astar);
    return newtonRaphson(isentropic.A_Astar_func, isentropic.A_Astar_deriv, guess, A_Astar, g);
  }
};

function normalShockRatios(M1, g) {
  if (M1 <= 1) return {};
  const M1_sq = M1 ** 2;
  const gm1 = g - 1;
  const gp1 = g + 1;

  const M2 = Math.sqrt((1 + gm1 / 2 * M1_sq) / (g * M1_sq - gm1 / 2));
  const p2p1 = 1 + 2 * g / gp1 * (M1_sq - 1);
  const rho2rho1 = gp1 * M1_sq / (gm1 * M1_sq + 2);
  const T2T1 = p2p1 / rho2rho1;
  const p02p01 = isentropic.p0_p(M2, g) * (p2p1) * (1 / isentropic.p0_p(M1, g));
  const p02p1 = isentropic.p0_p(M2, g) * p2p1;

  return { M1, M2, p2p1, rho2rho1, T2T1, p02p01, p02p1 };
}

const normalShockInverse = {
  M1_from_M2: (M2, g) => {
    if (M2 <= 0) return NaN;
    const M2_sq = M2 ** 2;
    const gm1 = g - 1;
    const M1_sq = (1 + gm1 / 2 * M2_sq) / (g * M2_sq - gm1 / 2);
    return M1_sq < 1 ? NaN : Math.sqrt(M1_sq);
  },
  M1_from_p2p1: (p2p1, g) => {
    if (p2p1 <= 1) return NaN;
    const M1_sq = (p2p1 - 1) * (g + 1) / (2 * g) + 1;
    return Math.sqrt(M1_sq);
  },
  M1_from_rho2rho1: (rho2rho1, g) => {
    const gp1 = g + 1;
    const gm1 = g - 1;
    if (rho2rho1 <= 1 || rho2rho1 > gp1 / gm1) return NaN;
    const M1_sq = (2 * rho2rho1) / (gp1 - gm1 * rho2rho1);
    return Math.sqrt(M1_sq);
  },
  M1_from_T2T1: (T2T1, g) => {
    if (T2T1 <= 1) return NaN;
    const gp1 = g + 1;
    const gm1 = g - 1;
    const a = 2 * g * gm1;
    const b = (4 * g - gm1**2 - gp1**2 * T2T1);
    const c = -2 * gm1;
    const discriminant = b**2 - 4*a*c;
    if (discriminant < 0) return NaN;
    const X1 = (-b + Math.sqrt(discriminant)) / (2*a);
    return X1 < 1 ? NaN : Math.sqrt(X1);
  },
  p02p01_func: (M1, g) => normalShockRatios(M1, g).p02p01,
  p02p01_deriv: (M1, g) => {
    const h = 1e-6;
    const f_plus = normalShockRatios(M1 + h, g).p02p01 || 0;
    const f_minus = normalShockRatios(M1 - h, g).p02p01 || 0;
    return (f_plus - f_minus) / (2 * h);
  },
  M1_from_p02p01: (p02p01, g) => {
    if (p02p01 >= 1 || p02p01 <= 0) return NaN;
    const guess = 1 - 2 * (p02p01 - 1);
    return newtonRaphson(normalShockInverse.p02p01_func, normalShockInverse.p02p01_deriv, guess < 1 ? 1.1 : guess, p02p01, g);
  },
  p02p1_func: (M1, g) => normalShockRatios(M1, g).p02p1,
  p02p1_deriv: (M1, g) => {
    const h = 1e-6;
    const f_plus = normalShockRatios(M1 + h, g).p02p1 || 0;
    const f_minus = normalShockRatios(M1 - h, g).p02p1 || 0;
    return (f_plus - f_minus) / (2 * h);
  },
  M1_from_p02p1: (p02p1, g) => {
    if (p02p1 <= 1) return NaN;
    const guess = 1 + (p02p1 / 5);
    return newtonRaphson(normalShockInverse.p02p1_func, normalShockInverse.p02p1_deriv, guess < 1 ? 1.1 : guess, p02p1, g);
  }
};

function obliqueShockRatios(M1, beta, theta, g) {
  const Mn1 = M1 * Math.sin(beta);
  if (Mn1 <= 1) return { M1, theta: rad2deg(theta), beta: rad2deg(beta) };
  
  const { M2: Mn2, p2p1, rho2rho1, T2T1, p02p01 } = normalShockRatios(Mn1, g);
  if (isNaN(Mn2)) return { M1, theta: rad2deg(theta), beta: rad2deg(beta) };
  
  const M2 = Mn2 / Math.sin(beta - theta);
  return { M1, M2, Mn1, Mn2, p2p1, rho2rho1, T2T1, p02p01, theta: rad2deg(theta), beta: rad2deg(beta) };
}

function ob_beta_from_M1_theta(M1, theta, g) {
  const thetaRad = deg2rad(theta);
  const mu = Math.asin(1 / M1);
  const c = Math.tan(mu) ** 2;
  const a = ((g - 1) / 2 + ((g + 1) / 2) * c) * Math.tan(thetaRad);
  const b = ((g + 1) / 2 + ((g + 3) / 2) * c) * Math.tan(thetaRad);
  const num = 4 * (1 - 3 * a * b) ** 3;
  const den = (27 * a * a * c + 9 * a * b - 2) ** 2;
  const d = Math.sqrt(num / den - 1);
  const term1 = (b + 9 * a * c) / (2 * (1 - 3 * a * b));
  const term2 = (d * (27 * a * a * c + 9 * a * b - 2)) / (6 * a * (1 - 3 * a * b));
  const angle = (1 / 3) * Math.atan(1 / d); // n=0 for weak shock
  return Math.atan(term1 - term2 * Math.tan(angle));
}

function ob_theta_from_M1_beta(M1, beta, g) {
  const betaRad = deg2rad(beta);
  const M1_sq = M1 ** 2;
  const num = (M1_sq * Math.sin(betaRad) ** 2 - 1);
  const den = M1_sq * (g + Math.cos(2 * betaRad)) + 2;
  const thetaRad = Math.atan(2 * (1/Math.tan(betaRad)) * num / den);
  return rad2deg(thetaRad);
}

function ob_M1_from_theta_beta(theta, beta, g) {
  const thetaRad = deg2rad(theta);
  const betaRad = deg2rad(beta);
  if (beta <= theta) return NaN;
  const T = Math.tan(thetaRad);
  const C = 1 / Math.tan(betaRad);
  const S2 = Math.sin(betaRad) ** 2;
  const G = g + Math.cos(2 * betaRad);
  const M1_sq = 2 * (C + T) / (2 * C * S2 - T * G);
  return M1_sq < 1 ? NaN : Math.sqrt(M1_sq);
}

const prandtlMeyer = {
  nu: (M, g) => {
    if (M < 1) return 0;
    if (M === Infinity) {
        const gp1 = g + 1;
        const gm1 = g - 1;
        return rad2deg(Math.sqrt(gp1 / gm1) * Math.PI / 2 - Math.PI / 2);
    }
    const gp1 = g + 1;
    const gm1 = g - 1;
    const g_ratio = Math.sqrt(gp1 / gm1);
    const M_sq = M ** 2;
    const term1 = Math.sqrt(gm1 / gp1 * (M_sq - 1));
    const term2 = Math.sqrt(M_sq - 1);
    return rad2deg(g_ratio * Math.atan(term1) - Math.atan(term2));
  },
  nu_rad: (M, g) => deg2rad(prandtlMeyer.nu(M, g)),
  nu_deriv: (M, g) => {
    if (M <= 1) return 0;
    const M_sq = M ** 2;
    const gm1 = g - 1;
    return (Math.sqrt(M_sq - 1) / (1 + gm1 / 2 * M_sq)) * (1 / M);
  },
  M_from_nu: (nu_deg, g) => {
    if (nu_deg < 0) return 1;
    const gp1 = g + 1;
    const gm1 = g - 1;
    const nu_max = rad2deg(Math.sqrt(gp1 / gm1) * Math.PI / 2 - Math.PI / 2);
    if (nu_deg > nu_max) return Infinity;
    const nu_rad_target = deg2rad(nu_deg);
    const guess = 1.5 + (nu_deg / nu_max) * 10;
    return newtonRaphson(prandtlMeyer.nu_rad, prandtlMeyer.nu_deriv, guess, nu_rad_target, g);
  }
};


// ========================================================================
// === Diagram Updaters ===================================================
// ========================================================================

function updateNormalShockDiagram(M1, M2) {
  updateText("ns-M1-label", M1, 2);
  updateText("ns-M2-label", M2, 2);
}

function updateObliqueDiagram(thetaDeg, betaDeg) {
  const x0 = 280, y0 = 320;
  const Lbase = 220, Lwedge = 150;
  const isValidTheta = !isNaN(thetaDeg) && thetaDeg > 0;
  const thetaRad = isValidTheta ? deg2rad(thetaDeg) : 0;
  
  setVisible("ob-wedgeBase", isValidTheta);
  setVisible("ob-wedgeRamp", isValidTheta);
  setVisible("ob-thetaArc", isValidTheta);
  setVisible("ob-thetaLabel", isValidTheta);

  if (isValidTheta) {
    document.getElementById("ob-wedgeBase").setAttribute("x2", x0 + Lbase);
    document.getElementById("ob-wedgeRamp").setAttribute("x2", x0 + Lwedge * Math.cos(thetaRad));
    document.getElementById("ob-wedgeRamp").setAttribute("y2", y0 - Lwedge * Math.sin(thetaRad));
    document.getElementById("ob-thetaArc").setAttribute("d", arcPath(x0, y0, 40, 0, thetaRad));
    // *** Updated position to be further right ***
    document.getElementById("ob-thetaLabel").setAttribute("x", x0 + 135 * Math.cos(thetaRad / 2));
    document.getElementById("ob-thetaLabel").setAttribute("y", y0 - 135 * Math.sin(thetaRad / 2));
    document.getElementById("ob-thetaLabel").textContent = `θ=${thetaDeg.toFixed(1)}°`;
  }

  const isValidBeta = !isNaN(betaDeg) && betaDeg > 0;
  setVisible("ob-shockLine", isValidBeta);
  setVisible("ob-betaArc", isValidBeta);
  setVisible("ob-betaLabel", isValidBeta);

  if (isValidBeta) {
    const Lshock = 220;
    const betaRad = deg2rad(betaDeg);
    document.getElementById("ob-shockLine").setAttribute("x2", x0 + Lshock * Math.cos(betaRad));
    document.getElementById("ob-shockLine").setAttribute("y2", y0 - Lshock * Math.sin(betaRad));
    document.getElementById("ob-betaArc").setAttribute("d", arcPath(x0, y0, 60, 0, betaRad));
    // *** Updated position to be further right ***
    document.getElementById("ob-betaLabel").setAttribute("x", x0 + 160 * Math.cos(betaRad / 2));
    document.getElementById("ob-betaLabel").setAttribute("y", y0 - 160 * Math.sin(betaRad / 2));
    document.getElementById("ob-betaLabel").textContent = `β=${betaDeg.toFixed(1)}°`;
  }
}

// ========================================================================
// === Calculator Controllers =============================================
// ========================================================================

function clearResults(prefix, fields) {
  fields.forEach(id => updateText(`${prefix}-${id}`, NaN));
}

// --- 1. Isentropic Controller ---
function calculateIsentropic() {
  const g = getGamma();
  const inputType = document.getElementById("is-inputType").value;
  const inputValue = parseFloat(document.getElementById("is-inputValue").value);
  const warning = document.getElementById("is-warning");
  warning.textContent = "";
  let M = NaN, T0_T, p0_p, rho0_rho, A_Astar;
  const resultIDs = ["is-M", "is-p0_p", "is-rho0_rho", "is-T0_T", "is-A_Astar"];

  try {
    switch (inputType) {
      case "M":
        M = inputValue;
        if (M <= 0) throw new Error("M must be > 0.");
        break;
      case "T0_T":
        T0_T = inputValue;
        if (T0_T < 1) throw new Error("T₀/T must be ≥ 1.");
        M = isentropic.M_from_T(T0_T, g);
        break;
      case "p0_p":
        p0_p = inputValue;
        if (p0_p < 1) throw new Error("p₀/p must be ≥ 1.");
        M = isentropic.M_from_p(p0_p, g);
        break;
      case "rho0_rho":
        rho0_rho = inputValue;
        if (rho0_rho < 1) throw new Error("ρ₀/ρ must be ≥ 1.");
        M = isentropic.M_from_rho(rho0_rho, g);
        break;
      case "A_Astar_sub":
      case "A_Astar_super":
        A_Astar = inputValue;
        if (A_Astar < 1) throw new Error("A/A* must be ≥ 1.");
        M = isentropic.M_from_A(A_Astar, g, inputType === "A_Astar_super");
        break;
    }
    if (isNaN(M)) throw new Error("Calculation failed. Check inputs or solver failed.");
    T0_T = isentropic.T0_T(M, g);
    p0_p = isentropic.p0_p(M, g);
    rho0_rho = isentropic.rho0_rho(M, g);
    A_Astar = isentropic.A_Astar(M, g);
    flashResults("isentropic");
  } catch (err) {
    warning.textContent = err.message;
    clearResults("is", resultIDs.map(id => id.split('-')[1]));
  }
  updateText("is-M", M);
  updateText("is-p0_p", p0_p);
  updateText("is-rho0_rho", rho0_rho);
  updateText("is-T0_T", T0_T);
  updateText("is-A_Astar", A_Astar);
}

// --- 2. Normal Shock Controller ---
function calculateNormalShock() {
  const g = getGamma();
  const inputType = document.getElementById("ns-inputType").value;
  const inputValue = parseFloat(document.getElementById("ns-inputValue").value);
  const warning = document.getElementById("ns-warning");
  warning.textContent = "";
  let M1 = NaN;
  const resultIDs = ["M1", "M2", "p2p1", "rho2rho1", "T2T1", "p02p01", "p02p1"];
  
  try {
    switch (inputType) {
      case "M1":
        M1 = inputValue;
        if (M1 <= 1) throw new Error("M₁ must be > 1 for a normal shock.");
        break;
      case "M2":
        M1 = normalShockInverse.M1_from_M2(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid M₂. Must be < 1 and > sqrt((g-1)/2g).");
        break;
      case "p2_p1":
        M1 = normalShockInverse.M1_from_p2p1(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid p₂/p₁. Must be > 1.");
        break;
      case "rho2_rho1":
        M1 = normalShockInverse.M1_from_rho2rho1(inputValue, g);
        if (isNaN(M1)) throw new Error(`Invalid ρ₂/ρ₁. Must be 1 < ratio < ${(g+1)/(g-1)}.`);
        break;
      case "T2_T1":
        M1 = normalShockInverse.M1_from_T2T1(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid T₂/T₁. Must be > 1.");
        break;
      case "p02_p01":
        M1 = normalShockInverse.M1_from_p02p01(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid p₀₂/p₀₁. Must be < 1 and > 0.");
        break;
      case "p02_p1":
        M1 = normalShockInverse.M1_from_p02p1(inputValue, g);
        if (isNaN(M1)) throw new Error("Invalid p₀₂/p₁. Must be > 1.");
        break;
    }
    if (isNaN(M1) || M1 <= 1) throw new Error("Calculation failed. Ensure input is in a valid range.");
    
    const { M2, p2p1, rho2rho1, T2T1, p02p01, p02p1 } = normalShockRatios(M1, g);
    updateText("ns-M1", M1);
    updateText("ns-M2", M2);
    updateText("ns-p2p1", p2p1);
    updateText("ns-rho2rho1", rho2rho1);
    updateText("ns-T2T1", T2T1);
    updateText("ns-p02p01", p02p01);
    updateText("ns-p02p1", p02p1);
    updateNormalShockDiagram(M1, M2);
    flashResults("normal-shock");
  } catch (err) {
    warning.textContent = err.message;
    clearResults("ns", resultIDs);
    updateNormalShockDiagram(NaN, NaN);
  }
}

// --- 3. Oblique Shock Controller ---
function clearObliqueInputs() {
    const inputs = ["ob-M1", "ob-theta", "ob-beta"];
    inputs.forEach(id => {
        document.getElementById(id).value = "";
        document.getElementById(id).classList.remove("is-calculated-input");
    });
    document.getElementById("ob-warning").textContent = "";
}

function calculateObliqueShock() {
  const g = getGamma();
  const M1_in = parseFloat(document.getElementById("ob-M1").value);
  const theta_in = parseFloat(document.getElementById("ob-theta").value);
  const beta_in = parseFloat(document.getElementById("ob-beta").value);
  const warning = document.getElementById("ob-warning");
  warning.textContent = "";
  
  const allInputs = [document.getElementById("ob-M1"), document.getElementById("ob-theta"), document.getElementById("ob-beta")];
  allInputs.forEach(el => el.classList.remove("is-calculated-input"));

  const resultIDs = ["M1n", "M2n", "M2", "p2p1", "rho2rho1", "T2T1", "p02p01"];
  let M1 = M1_in, theta = theta_in, beta = beta_in;
  let betaRad, thetaRad, calculatedField = "";

  try {
    const inputs = [!isNaN(M1), !isNaN(theta), !isNaN(beta)];
    const inputCount = inputs.filter(Boolean).length;
    if (inputCount !== 2) {
      throw new Error("Please provide exactly two inputs (M₁, θ, or β).");
    }

    if (!isNaN(M1) && !isNaN(theta)) { // Case 1: (M1, θ) -> β
      if (M1 <= 1) throw new Error("M₁ must be > 1.");
      if (theta < 0) throw new Error("θ must be ≥ 0.");
      betaRad = ob_beta_from_M1_theta(M1, theta, g);
      beta = rad2deg(betaRad);
      if (isNaN(beta)) throw new Error("No attached shock solution (θ > θ_max).");
      thetaRad = deg2rad(theta);
      calculatedField = "ob-beta";
    } else if (!isNaN(M1) && !isNaN(beta)) { // Case 2: (M1, β) -> θ
      if (M1 <= 1) throw new Error("M₁ must be > 1.");
      if (beta <= 0 || beta > 90) throw new Error("β must be in (0, 90].");
      const mu = rad2deg(Math.asin(1/M1));
      if (beta <= mu) throw new Error(`β (${beta.toFixed(1)}°) must be > μ (${mu.toFixed(1)}°).`);
      theta = ob_theta_from_M1_beta(M1, beta, g);
      if (theta < 0) throw new Error("Detached shock region (strong shock).");
      betaRad = deg2rad(beta);
      thetaRad = deg2rad(theta);
      calculatedField = "ob-theta";
    } else if (!isNaN(theta) && !isNaN(beta)) { // Case 3: (θ, β) -> M1
      if (beta <= theta) throw new Error("β must be > θ.");
      if (theta < 0) throw new Error("θ must be ≥ 0.");
      M1 = ob_M1_from_theta_beta(theta, beta, g);
      if (isNaN(M1)) throw new Error("Invalid (θ, β) pair.");
      const mu = rad2deg(Math.asin(1/M1));
      if (beta <= mu) throw new Error(`Implied Mach wave (β ≤ μ).`);
      betaRad = deg2rad(beta);
      thetaRad = deg2rad(theta);
      calculatedField = "ob-M1";
    }
    
    const { Mn1, Mn2, M2, p2p1, rho2rho1, T2T1, p02p01 } = obliqueShockRatios(M1, betaRad, thetaRad, g);
    if (isNaN(M2)) throw new Error("Calculation failed. Check inputs.");

    updateText("ob-M1", M1, 3);
    updateText("ob-theta", theta, 2);
    updateText("ob-beta", beta, 2);
    if (calculatedField) {
        document.getElementById(calculatedField).classList.add("is-calculated-input");
    }
    
    updateText("ob-M1n", Mn1);
    updateText("ob-M2n", Mn2);
    updateText("ob-M2", M2);
    updateText("ob-p2p1", p2p1);
    updateText("ob-rho2rho1", rho2rho1);
    updateText("ob-T2T1", T2T1);
    updateText("ob-p02p01", p02p01);
    
    updateObliqueDiagram(theta, beta);
    flashResults("oblique-shock");

  } catch (err) {
    warning.textContent = err.message;
    clearResults("ob", resultIDs); 
    if (isNaN(M1_in)) updateText("ob-M1", NaN);
    if (isNaN(theta_in)) updateText("ob-theta", NaN);
    if (isNaN(beta_in)) updateText("ob-beta", NaN);
    updateObliqueDiagram(isNaN(theta) ? theta_in : theta, isNaN(beta) ? beta_in : beta);
  }
}

// --- 4. Prandtl-Meyer Controller ---
function calculatePrandtlMeyer() {
  const g = getGamma();
  const M1 = parseFloat(document.getElementById("pm-M1").value);
  const thetaDeg = parseFloat(document.getElementById("pm-theta").value);
  const warning = document.getElementById("pm-warning");
  warning.textContent = "";
  const resultIDs = ["nu1", "mu1", "nu2", "mu2", "M2", "p2p1", "rho2rho1", "T2T1"];
  
  try {
    if (isNaN(M1) || M1 <= 1) throw new Error("M₁ must be > 1 for an expansion.");
    if (isNaN(thetaDeg) || thetaDeg <= 0) throw new Error("θ must be > 0.");

    const nu1 = prandtlMeyer.nu(M1, g);
    const mu1 = rad2deg(Math.asin(1/M1));
    const nu2 = nu1 + thetaDeg;
    
    const nu_max = prandtlMeyer.nu(Infinity, g);
    if (nu2 > nu_max) {
      throw new Error(`Expansion angle is too large. Max ν is ${nu_max.toFixed(1)}°.`);
    }
    
    const M2 = prandtlMeyer.M_from_nu(nu2, g);
    const mu2 = M2 === Infinity ? 0 : rad2deg(Math.asin(1/M2));

    const p2p1 = isentropic.p0_p(M1, g) / isentropic.p0_p(M2, g);
    const rho2rho1 = isentropic.rho0_rho(M1, g) / isentropic.rho0_rho(M2, g);
    const T2T1 = isentropic.T0_T(M1, g) / isentropic.T0_T(M2, g);

    updateText("pm-nu1", nu1, 2);
    updateText("pm-mu1", mu1, 2);
    updateText("pm-nu2", nu2, 2);
    updateText("pm-mu2", mu2, 2);
    updateText("pm-M2", M2);
    updateText("pm-p2p1", p2p1);
    updateText("pm-rho2rho1", rho2rho1);
    updateText("pm-T2T1", T2T1);
    flashResults("prandtl-meyer");

  } catch (err) {
    warning.textContent = err.message;
    clearResults("pm", resultIDs);
  }
}


// ========================================================================
// === Initialization =====================================================
// ========================================================================

function runAllCalculators() {
  calculateIsentropic();
  calculateNormalShock();
  // Don't run oblique shock on load
  calculatePrandtlMeyer();
}

document.addEventListener("DOMContentLoaded", () => {
  // --- Main Tab Navigation ---
  const navButtons = document.querySelectorAll(".nav-btn");
  const panels = document.querySelectorAll(".calculator-panel");
  navButtons.forEach(btn => {
    btn.addEventListener("click", () => {
      navButtons.forEach(b => b.classList.remove("active"));
      panels.forEach(p => p.classList.remove("active"));
      btn.classList.add("active");
      document.getElementById(`${btn.dataset.calculator}-calculator`).classList.add("active");
    });
  });

  // --- Theme Toggle ---
  const themeToggle = document.getElementById("theme-toggle-checkbox");
  const storedTheme = localStorage.getItem("theme");
  const prefersLight = window.matchMedia("(prefers-color-scheme: light)").matches;
  const applyTheme = (theme) => {
    document.body.dataset.theme = theme;
    themeToggle.checked = (theme === "light");
  };
  if (storedTheme) {
    applyTheme(storedTheme);
  } else {
    applyTheme(prefersLight ? "light" : "dark");
  }
  themeToggle.addEventListener("change", () => {
    const newTheme = themeToggle.checked ? "light" : "dark";
    localStorage.setItem("theme", newTheme);
    applyTheme(newTheme);
  });

  // --- Font Size Control ---
  const fontIncrease = document.getElementById("font-increase");
  const fontDecrease = document.getElementById("font-decrease");
  const fontSizes = ["xs", "sm", "md", "lg", "xl"];
  let currentSizeIndex = fontSizes.indexOf(localStorage.getItem("font-size") || "md");
  if (currentSizeIndex === -1) currentSizeIndex = 2; // Default to 'md'
  
  const setFontSize = (index) => {
    if (index < 0) index = 0;
    if (index >= fontSizes.length) index = fontSizes.length - 1;
    currentSizeIndex = index;
    const size = fontSizes[index];
    document.body.dataset.fontSize = size;
    localStorage.setItem("font-size", size);
  };
  
  fontIncrease.addEventListener("click", () => setFontSize(currentSizeIndex + 1));
  fontDecrease.addEventListener("click", () => setFontSize(currentSizeIndex - 1));
  
  // Set initial font size
  setFontSize(currentSizeIndex);


  // --- Setup Calculator Event Listeners ---
  document.getElementById("is-computeBtn").addEventListener("click", calculateIsentropic);
  document.getElementById("is-inputType").addEventListener("change", calculateIsentropic);
  
  document.getElementById("ns-computeBtn").addEventListener("click", calculateNormalShock);
  document.getElementById("ns-inputType").addEventListener("change", calculateNormalShock);
  
  document.getElementById("ob-computeBtn").addEventListener("click", calculateObliqueShock);
  document.getElementById("ob-clearBtn").addEventListener("click", clearObliqueInputs);
  
  document.getElementById("pm-computeBtn").addEventListener("click", calculatePrandtlMeyer);
  
  document.getElementById("gamma").addEventListener("change", runAllCalculators);

  // --- Run calculators on initial load ---
  runAllCalculators();
});
